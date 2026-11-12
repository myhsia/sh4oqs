#!/usr/bin/env python3
"""Verify a frozen driver through its **normal** command line.

    python3 tools/verify_binary.py <moments|dmrg|sh> <executable> [--keep]

The binary takes its configs as arguments and carries none, so this writes a
temporary YAML with the reference-scenario parameters and then runs the binary
exactly as a user would: a few moments by both engines, a short tDMRG run, or a
seeded MASH ensemble.  The produced ``.npz`` is read back and compared with the
values below.  Nothing is added to the driver for this: the check lives here,
the physics code stays as it was.

The archive is parsed with the standard library (``zipfile`` + ``struct`` on the
``.npy`` members), so no numpy is needed — the same file works in CI, on a Mac,
and on the cluster's ``python3``.
"""
import argparse
import ast
import json
import os
import shutil
import struct
import subprocess
import sys
import tempfile
import time
import zipfile

# The reference scenario: the same parameters as config/A-base_line.yaml (the
# plotting-only keys are irrelevant here).  JSON is valid YAML, so this is read
# by the driver's yaml.safe_load as written.
CONFIG = {
    "Pade_m": 9, "Pade_n": 14, "linear": True, "quadratic": False,
    "Delta": 1.0, "epsilon": 1.0, "beta": 0.6, "lambd": 0.4,
    "theta_B": 1.33, "wB": 1.0, "wC": 4.0, "anharmon": 0.4, "scale": 0.3,
    "Cayley_depth": 16, "Cayley_orders": 16, "xlims": [[0, 12], 4],
}

# Measured on macOS arm64 (Python 3.14.7, numpy 2.5.3/2.5.4) and reproduced on
# the CI runners and the cluster; the tolerances are far wider than the platform
# spread of the mpmath quadrature, the CSR iterate, the TDVP steps and the seeded
# ensemble.
MOMENTS = {  # complex Omega_1..3 of qhm, by engine
    "block_sp": [(0.0, 0.6), (-1.137824048904, 0.0), (0.0, -1.832083288027)],
    "symbolic": [(0.0, 0.6), (-1.137824048904, 0.0), (0.0, -1.832083288027)],
}
DMRG_C00 = 5.047893617910e-01   # Re C_00 at t = 1, N=4 nmax=2 chi=8 dt=0.25
SH_C00 = 0.786909497             # C_00 at t = 0.5, MASH 8000 trajectories

#: dtype -> struct code, by (descriptor, itemsize).
_CODES = {("<f8", 8): "d", (">f8", 8): "d", ("<f4", 4): "f", (">f4", 4): "f",
          ("<c16", 16): "dd", (">c16", 16): "dd",
          ("<c8", 8): "ff", (">c8", 8): "ff"}


# --------------------------------------------------------------------------- #
# .npy reading without numpy
# --------------------------------------------------------------------------- #
def read_npy(buffer):
    """Decode one ``.npy`` buffer into ``(shape, flat values)``."""
    if buffer[:6] != b"\x93NUMPY":
        raise ValueError("not a .npy buffer")
    if buffer[6] == 1:
        header_len = struct.unpack("<H", buffer[8:10])[0]
        start = 10
    else:
        header_len = struct.unpack("<I", buffer[8:12])[0]
        start = 12
    header = ast.literal_eval(buffer[start:start + header_len].decode("latin1"))
    descr = header["descr"]
    shape = tuple(header["shape"])
    raw = buffer[start + header_len:]
    try:
        code = _CODES[(descr, int(descr[2:]))]
    except (KeyError, ValueError):
        raise ValueError(f"unsupported dtype {descr!r}") from None
    order = ">" if descr.startswith(">") else "<"
    count = len(raw) // struct.calcsize(code)
    values = list(struct.unpack(order + code * count, raw))
    if len(code) == 2:                       # complex: interleaved re/im
        values = [complex(values[2 * i], values[2 * i + 1])
                  for i in range(count)]
    return shape, values


def load_npz(path):
    """Read the numeric ``.npy`` members of an ``.npz`` into ``{key: (shape, values)}``.

    The metadata members (``branch``/``component``/``identifier``) are arrays of
    strings, which this reader does not decode and the checks do not need, so a
    non-numeric member is skipped rather than treated as an error.
    """
    arrays = {}
    with zipfile.ZipFile(path) as archive:
        for member in archive.namelist():
            if not member.endswith(".npy"):
                continue
            try:
                arrays[member[:-4]] = read_npy(archive.read(member))
            except ValueError:
                continue
    return arrays


def element(array, *index):
    """Row-major element of a ``(shape, values)`` pair."""
    shape, values = array
    if len(index) != len(shape):
        raise ValueError(f"{len(index)} indices for shape {shape}")
    offset = 0
    for size, i in zip(shape, index):
        offset = offset * size + i
    return values[offset]


def identity_defect(C):
    """``max |C(0) - I|`` over the 4x4 basis block."""
    return max(abs(element(C, 0, i, j) - (1.0 if i == j else 0.0))
               for i in range(4) for j in range(4))


# --------------------------------------------------------------------------- #
# Checks
# --------------------------------------------------------------------------- #
def report(name, ok, detail):
    print(f"verify: {name} {'ok' if ok else 'FAIL'} ({detail})")
    return ok


def run(binary, args, cwd):
    """Run the binary; on failure show the tail of its own output and stop."""
    done = subprocess.run([binary] + args, cwd = cwd, capture_output = True,
                          text = True, timeout = 1800)
    if done.returncode != 0:
        sys.stderr.write(done.stdout[-2000:] + done.stderr[-2000:])
        raise SystemExit(f"error: {os.path.basename(binary)} exited "
                         f"{done.returncode} for {' '.join(args)}")
    return done.stdout


def verify_moments(binary, config, tmp):
    ok = True
    data = os.path.join(tmp, "data")
    for engine, reference in MOMENTS.items():
        run(binary, ["-prefix", config, "-method", engine, "-nmax", "3",
                     "-nproc", "1", "-njobs", "1", "-outdir", data], tmp)
        got = load_npz(os.path.join(data, f"moments.{engine}.scl.npz"))["qhm"]
        dev = 0.0
        for i, (re, im) in enumerate(reference):
            want = complex(re, im)
            dev = max(dev, abs(element(got, i) - want) / max(abs(want), 1e-300))
        ok &= report(f"moments {engine}", dev <= 1e-6, f"{dev:.2e} relative")
    return ok


def verify_dmrg(binary, config, tmp):
    data = os.path.join(tmp, "data")
    run(binary, [config, "-case", "qhm", "-geometry", "chain", "-N", "4",
                 "-nmax", "2", "-chi", "8", "-dt", "0.25", "-dtau", "0.1",
                 "-wmax_factor", "4", "-t1", "1", "-outdir", data], tmp)
    stem = os.path.splitext(os.path.basename(config))[0]
    arrays = load_npz(os.path.join(data, stem, "corrmat_dmrg_qhm.npz"))
    t, C = arrays["t_qhm"], arrays["qhm"]
    last = t[0][0] - 1                     # (n,) -> the final sample
    c00 = element(C, last, 0, 0)
    ok = report("dmrg C(0)", identity_defect(C) <= 1e-8,
                f"|C(0)-I|max {identity_defect(C):.2e}")
    dev = abs(c00.real - DMRG_C00) / abs(DMRG_C00)
    ok &= report("dmrg C00", dev <= 1e-3, f"Re C00 {c00.real:.12e} ({dev:.2e})")
    return ok


def verify_sh(binary, config, tmp):
    data = os.path.join(tmp, "data")
    run(binary, [config, "-method", "mash", "-n_traj", "8000", "-N", "8",
                 "-wmax_factor", "4", "-dt", "0.05", "-seed", "0", "-nproc", "1",
                 "-njobs", "1", "-outdir", data], tmp)
    stem = os.path.splitext(os.path.basename(config))[0]
    arrays = load_npz(os.path.join(data, stem, "corrmat_mash.npz"))
    t, C = arrays["t_qhm"], arrays["qhm"]
    times = t[1]
    index = min(range(len(times)), key = lambda i: abs(times[i] - 0.5))
    c00 = element(C, index, 0, 0)
    ok = report("sh C(0)", identity_defect(C) <= 5e-2,
                f"|C(0)-I|max {identity_defect(C):.2e}")
    dev = abs(c00 - SH_C00) / abs(SH_C00)
    ok &= report("sh C00", dev <= 2e-2,
                 f"t={times[index]:g} {c00.real:.12f} ({dev:.2e})")
    return ok


VERIFIERS = {"moments": verify_moments, "dmrg": verify_dmrg, "sh": verify_sh}


def main():
    parser = argparse.ArgumentParser(
        description = "Verify a frozen driver through its normal command line.")
    parser.add_argument("driver", choices = sorted(VERIFIERS))
    parser.add_argument("binary", help = "path to the frozen executable")
    parser.add_argument("--keep", action = "store_true",
                        help = "leave the temporary run directory in place")
    args = parser.parse_args()

    binary = os.path.abspath(args.binary)
    if not os.access(binary, os.X_OK):
        raise SystemExit(f"error: {binary} is not executable")
    tmp = tempfile.mkdtemp(prefix = f"oqs-verify-{args.driver}-")
    started = time.time()
    try:
        config = os.path.join(tmp, "config", "smoke.yaml")
        os.makedirs(os.path.dirname(config), exist_ok = True)
        with open(config, "w") as f:
            json.dump(CONFIG, f, indent = 2)
        print(f"verify: {os.path.basename(binary)} ({args.driver})")
        ok = VERIFIERS[args.driver](binary, config, tmp)
        print(f"verify: {args.driver} {'passed' if ok else 'FAILED'} "
              f"in {time.time() - started:.1f} s")
        return 0 if ok else 1
    finally:
        if args.keep:
            print(f"verify: kept {tmp}")
        else:
            shutil.rmtree(tmp, ignore_errors = True)


if __name__ == "__main__":
    raise SystemExit(main())
