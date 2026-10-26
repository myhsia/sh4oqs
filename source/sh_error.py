# %%
"""Accuracy of the surface-hopping references against HEOM.

The metric is the six-curve variant of the shared relative ``L2`` error
(``library/mre.py``): over the calculation window ``[0, tf]`` of each config,
the relative ``L2`` error of six curves

    {Re P0, Re P1, Re C_xx, Re C_xy, Im C_xx, Im C_xy}

against the HEOM reference, averaged over ``mre.N_OBS = 6`` and reported in
percent, ``(E / 6) * 100``.  Only the window is compared: numerator and
denominator both use the same slice.  The reference observables are rebuilt
from the HEOM correlation matrix (``corrmat_heom.npz``) and the MKCT ones from
``corrmat_mkct.npz``, both in the canonical ``mkct4oqs`` checkout, with
``library/corrmat.py``.

This is the runner; the algorithm is imported from the shared ``mre`` module.
It writes ``output/sh_error.json`` (one value per config/method, each row
carrying the window ``tf``) and prints the table.

Usage::

    python sh_error.py                       # RUNLIST x all methods
    python sh_error.py --methods mash fssh
    python sh_error.py --configs A_base_line --labels "" best maxA
"""
# %%
# 1. Requirements
import argparse, json, os, sys
import numpy as np
sys.path.append(os.path.join(os.path.dirname(__file__), "./library"))
sys.pycache_prefix = os.path.abspath("./output/pycache")
import mre as core
from corrmat import observables
from fsutils import load_config, load_npz, load_group, ensure_dir
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
LABELS  = ["", "base", "wf8", "wf14", "N400", "N800", "dt0025", "traj160k",
           "best", "maxA", "maxB"]

# The HEOM and MKCT references live in the canonical mkct4oqs checkout.
REF  = os.environ.get("SH_REF_DATA", "../../mkct4oqs/source/output/data")
JSON = "./output/sh_error.json"
MD   = "./output/sh_error.md"


def _series(prefix, method, case, quantity, label=""):
    """``(t, y)`` of one observable for one (method, case, label)."""
    d = "./output/data/" + prefix + "/"
    tag = f"{method}_{case}" + (f"_{label}" if label else "")
    merged = d + ("cohe" if quantity == "cohem" else quantity)
    merged += f"_{method}" + (f"_{label}" if label else "") + ".npz"
    if not os.path.exists(merged):
        if not label:
            raise FileNotFoundError(f"{merged}: run sh_merge.py first")
        # labelled sweep runs stay per-case (they are not merged)
        return load_npz(d + f"{quantity}_{tag}.npz")
    data, _ = load_group(merged)
    t, y = data[case]
    if quantity == "cohem":
        y = y[1]                     # branch ("+", "-") on the leading axis
    elif quantity == "cohe":
        y = y[0]
    return t, y


def _ref(prefix):
    """``(t, obs)`` of the HEOM reference, rebuilt from its correlation matrix."""
    path = os.path.join(REF, prefix) + os.sep + "corrmat_heom.npz"
    t, C = load_group(path)[0]["qhm"]
    return t, observables(C)


def _mkct(prefix, case):
    """``(t, obs)`` of one MKCT case, rebuilt from the canonical corrmat."""
    path = os.path.join(REF, prefix) + os.sep + "corrmat_mkct.npz"
    t, C = load_group(path)[0][case]
    return t, observables(C)


def error(prefix, method, case, label=""):
    """Mean relative error (percent) over the six curves; None if absent."""
    tf = parse_xlims(load_config(prefix)["xlims"], reorder = True)[0][1]
    if method == "mkct":
        path = os.path.join(REF, prefix) + os.sep + "corrmat_mkct.npz"
        if not os.path.exists(path):
            return None
        t, obs = _mkct(prefix, case)
        pop, c, x = obs["pop"], obs["corr"], obs["xorr"]
    else:
        d = "./output/data/" + prefix + "/"
        tag = f"{method}_{case}" + (f"_{label}" if label else "")
        if not os.path.exists(d + f"pop_{tag}.npz"):
            return None
        t, pop = _series(prefix, method, case, "pop", label)
        _, c = _series(prefix, method, case, "corr", label)
        _, x = _series(prefix, method, case, "xorr", label)
    win = t <= tf * (1.0 + 1e-9)
    t, pop, c, x = t[win], pop[win], c[win], x[win]

    th, ref = _ref(prefix)
    hw = th <= tf * (1.0 + 1e-9)
    thw = th[hw]
    cmu = core._interp(ref["corr"][hw], thw, t)
    cx = core._interp(ref["xorr"][hw], thw, t)
    # The six curves of the canonical metric, exactly as mre.scan_config.
    e = (core._rel(pop[:, 0].real, np.interp(t, thw, ref["pop"][hw].real[:, 0]))
         + core._rel(pop[:, 1].real, np.interp(t, thw, ref["pop"][hw].real[:, 1]))
         + core._rel(c.real, cmu.real)
         + core._rel(x.real, cx.real)
         + core._rel(c.imag, cmu.imag)
         + core._rel(x.imag, cx.imag))
    return 100.0 * e / core.N_OBS


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
             "Metric `(E/6)*100 %` over `[0, tf]` for "
             "`{Re P0, Re P1, Re C_xx, Re C_xy, Im C_xx, Im C_xy}`.", "",
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
