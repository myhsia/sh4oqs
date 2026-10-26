# %%
"""Accuracy of the surface-hopping references against HEOM.

The metric is the shared correlation-matrix error of ``library/mre.py``: over
the calculation window ``[0, tf]`` of each config, the relative ``L2`` error of
the 32 real series of the basis correlation matrix (real and imaginary parts of
its 16 elements), each weighted equally, against the HEOM reference matrix and
averaged over ``mre.N_CORRMAT = 32``, reported in percent ``(E / 32) * 100``.
Only the window is compared: numerator and denominator both use the same slice.

Every method stores ``corrmat_<method>_<case>.npz`` (merged into
``corrmat_<method>.npz``) with the same basis as the canonical
``mkct4oqs`` ``corrmat_mkct.npz`` / ``corrmat_heom.npz``, so no observable
contraction is involved.

This is the runner; the algorithm is imported from the shared ``mre`` module.
It writes ``output/sh_error.json`` (one value per config/method, each row
carrying the window ``tf``) and prints the table.

Usage::

    python sh_error.py                       # RUNLIST x all methods
    python sh_error.py --methods mash fssh
    python sh_error.py --configs A_base_line --labels "" base all
"""
# %%
# 1. Requirements
import argparse, json, os, sys, warnings
import numpy as np
sys.path.append(os.path.join(os.path.dirname(__file__), "./library"))
sys.pycache_prefix = os.path.abspath("./output/pycache")
import mre as core
from fsutils import load_config, load_group, ensure_dir
from pltools import parse_xlims


# 2. Run selection
RUNLIST = [
    "A-base_line",
    "B-slow_system",
    "C-no_bias",
    "D-strong_tunneling",
    "E-fast_bath",
    "F-slow_bath",
    "G-lowT",
    "H-highT",
    "I-weak_coupling",
    "J-strong_coupling",
    "K-linear_quad",
]
CASES   = ["qhm", "qah", "chm", "cah"]
METHODS = ["mkct", "mash", "fssh", "ehrenfest", "lsc"]
# Suffixes of labelled runs (the parameter sweep of the reference scenario).
LABELS  = ["", "base", "N400", "N800", "wf8", "wf14", "traj160k", "all"]

# The HEOM and MKCT references live in the canonical mkct4oqs checkout.
REF  = os.environ.get("SH_REF_DATA", "../../mkct4oqs/source/output/data")
JSON = "./output/sh_error.json"
MD   = "./output/sh_error.md"


def _corrmat_path(prefix, method, case, label=""):
    """The archive holding one method's basis correlation matrix."""
    if method == "mkct":
        return os.path.join(REF, prefix) + os.sep + "corrmat_mkct.npz"
    d = "./output/data/" + prefix + "/"
    if label:
        # labelled sweep runs stay per case (they are not merged)
        return d + f"corrmat_{method}_{case}_{label}.npz"
    merged = d + f"corrmat_{method}.npz"
    if os.path.exists(merged):
        return merged
    return d + f"corrmat_{method}_{case}.npz"


def _corrmat(prefix, method, case, label=""):
    """``(t, C)`` of one (method, case, label) correlation matrix."""
    path = _corrmat_path(prefix, method, case, label)
    if not os.path.exists(path):
        return None
    return load_group(path)[0][case]


def _ref(prefix):
    """``(t, C)`` of the HEOM reference correlation matrix."""
    path = os.path.join(REF, prefix) + os.sep + "corrmat_heom.npz"
    return load_group(path)[0]["qhm"]


def error(prefix, method, case, label=""):
    """Mean relative error (percent) over the 32 corrmat series; None if absent."""
    got = _corrmat(prefix, method, case, label)
    if got is None:
        return None
    t, C = got
    tf = parse_xlims(load_config(prefix)["xlims"], reorder = True)[0][1]
    win = t <= tf * (1.0 + 1e-9)
    if t[win][-1] < tf * (1.0 - 1e-6):
        warnings.warn(
            f"[{prefix}] {method}_{case} reaches t={t[win][-1]:g} but the "
            f"metric window is [0, {tf:g}]; the overlap is scored.",
            stacklevel = 2)
    t, C = t[win], C[win]
    th, Ch = _ref(prefix)
    hw = th <= tf * (1.0 + 1e-9)
    Cr = core._interp_mat(Ch[hw], th[hw], t)
    return 100.0 * core._corrmat_error(C, Cr) / core.N_CORRMAT


def main():
    ap = argparse.ArgumentParser(description = __doc__.splitlines()[0])
    ap.add_argument("--configs", nargs = "*", default = RUNLIST)
    ap.add_argument("--methods", nargs = "*", default = METHODS)
    ap.add_argument("--labels", nargs = "*", default = LABELS)
    args = ap.parse_args()

    ensure_dir("./output")
    raw = {}
    for prefix in args.configs:
        rows = raw.setdefault(prefix, {})
        for method in args.methods:
            for label in args.labels:
                key = method if not label else f"{method}[{label}]"
                per = {c: error(prefix, method, c, label) for c in CASES}
                got = [v for v in per.values() if v is not None]
                if not got:
                    continue
                rows[key] = {"E": float(np.mean(got)),
                             "cases": {c: v for c, v in per.items()
                                       if v is not None},
                             "tf": float(
                    parse_xlims(load_config(prefix)["xlims"],
                                reorder = True)[0][1])}
    if os.path.exists(JSON):
        with open(JSON) as f:
            old = json.load(f)
        for p, rows in raw.items():          # merge, do not lose other runs
            old.setdefault(p, {}).update(rows)
        raw = old
    with open(JSON, "w") as f:
        json.dump(raw, f, indent = 2, sort_keys = True)
    lines = ["# Surface-hopping accuracy against HEOM", "",
             "Metric `(E/32)*100 %` over `[0, tf]` of the 32 real series of "
             "the basis correlation matrix (16 elements, real and imaginary "
             "parts, each weighted equally).", "",
             "| config | " + " | ".join(args.methods) + " |",
             "|---|" + "---:|" * len(args.methods)]
    for prefix in args.configs:
        cells = []
        for m in args.methods:
            v = raw.get(prefix, {}).get(m)
            cells.append(f"{v['E']:.1f}" if v else "--")
        lines.append(f"| {prefix} | " + " | ".join(cells) + " |")
    with open(MD, "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print("WROTE", JSON, "and", MD, flush = True)


if __name__ == "__main__":
    main()

# %%
