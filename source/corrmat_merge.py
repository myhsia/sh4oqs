# %%
"""Merge the per-case correlation matrices into one archive per method.

A CI job computes one ``(config, case)`` and writes
``corrmat_<method>_<case>.npz``; the analysis reads ``corrmat_<method>.npz``
with one matrix per case, the same layout as ``corrmat_mkct.npz``.  Only the
four matrices are merged -- the observables are contractions evaluated at plot
time, so there is nothing else to combine.

Usage::

    python corrmat_merge.py                      # every config found
    python corrmat_merge.py --configs A-base_line
"""
import argparse, glob, os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), "./library"))
sys.pycache_prefix = os.path.abspath("./output/pycache")
from fsutils import load_group, save_group
from corrmat import BASIS

CASES = ("qhm", "qah", "chm", "cah")


def merge(directory, method):
    series = {}
    for case in CASES:
        path = os.path.join(directory, f"corrmat_{method}_{case}.npz")
        if not os.path.exists(path):
            return None
        series[case] = load_group(path)[0][case]
    save_group(series, directory, f"corrmat_{method}.npz",
               identifier = r"$\mathcal{C}(t)$", branch = BASIS,
               component = BASIS, compress = True)
    return f"corrmat_{method}.npz"


def main():
    ap = argparse.ArgumentParser(description = __doc__.splitlines()[0])
    ap.add_argument("--configs", nargs = "*", default = None)
    ap.add_argument("--methods", nargs = "*",
                    default = ["fssh", "mash", "ehrenfest", "lsc"])
    ap.add_argument("--root", default = "./output/data")
    args = ap.parse_args()
    dirs = [os.path.join(args.root, c) for c in args.configs] if args.configs \
        else sorted(d for d in glob.glob(os.path.join(args.root, "*")) if os.path.isdir(d))
    n = 0
    for d in dirs:
        for m in args.methods:
            if merge(d, m):
                n += 1
    print(f"merged {n} corrmat archives under {args.root}")


if __name__ == "__main__":
    main()

# %%
