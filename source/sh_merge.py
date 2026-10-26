# %%
"""Group the per-case surface-hopping correlation matrices into one archive.

Each CI job computes one ``(config, case)``, so the raw output is
``corrmat_<method>_<case>.npz`` (one file per case, unique per job: no
concurrent writes; the archive holds ``t_<case>`` and the ``(nt, 4, 4)``
matrix under ``<case>``).  This script groups the four cases into
``corrmat_<method>.npz`` with ``fsutils.save_group``, the same layout as
``corrmat_mkct.npz`` and ``corrmat_heom.npz``.

When several independent seeds were run per case
(``corrmat_<method>_<case>_s0.npz``, ``_s1``, ...), their matrices are
averaged first: every file is already an ensemble average, so the mean of
equal-size chunks is the average over the union.

Labelled parameter-sweep runs (``corrmat_<method>_<case>_<label>.npz``) stay
per case and are read directly by the analysis; they are not merged.

Usage::

    python sh_merge.py                      # every config/method found
    python sh_merge.py --configs A-base_line
    python sh_merge.py --methods mash fssh
"""
# %%
import argparse, glob, os, sys
import numpy as np
sys.path.append(os.path.join(os.path.dirname(__file__), "./library"))
sys.pycache_prefix = os.path.abspath("./output/pycache")
from fsutils import load_group, save_group
from corrmat import BASIS

CASES = ("qhm", "qah", "chm", "cah")
METHODS = ["fssh", "mash", "ehrenfest", "lsc"]


def _seed_files(directory, method, case):
    """All files of one case, one per seed (or a single one)."""
    stem = f"corrmat_{method}_{case}"
    # Only the plain file and explicit seeds ``_s0``/``_s1``/...: the labelled
    # sweep files (``_base``, ``_maxA``, ...) belong to a different setting and
    # are not merged.
    files = []
    plain = os.path.join(directory, stem + ".npz")
    if os.path.exists(plain):
        files.append(plain)
    files += sorted(glob.glob(os.path.join(directory, stem + "_s[0-9]*.npz")))
    return files


def _load_avg(files, case):
    """Average the ``(t, C)`` matrices of several seed files."""
    pairs = [load_group(f)[0][case] for f in files]
    t = pairs[0][0]
    C = sum(C for _, C in pairs) / len(pairs)
    return t, C


def merge_config(prefix, methods, labels, root, clean=False):
    directory = os.path.join(root, prefix)
    made = []
    for method in methods:
        series = {}
        for case in CASES:
            files = _seed_files(directory, method, case)
            if files:
                series[case] = _load_avg(files, case)
        if len(series) != len(CASES):
            continue
        name = f"corrmat_{method}.npz"
        save_group(series, directory, name, identifier = r"$\mathcal{C}(t)$",
                   branch = BASIS, component = BASIS, compress = True)
        made.append(name)
        if clean:
            for case in CASES:
                for f in _seed_files(directory, method, case):
                    os.remove(f)
    return made


def main():
    ap = argparse.ArgumentParser(description = __doc__.splitlines()[0])
    ap.add_argument("--configs", nargs = "*", default = None)
    ap.add_argument("--methods", nargs = "*", default = METHODS)
    ap.add_argument("--labels", nargs = "*", default = [""],
                    help = "kept for compatibility; labelled runs are not merged")
    ap.add_argument("--root", default = "./output/data")
    ap.add_argument("--clean", action = "store_true",
                    help = "delete the per-case files once merged")
    args = ap.parse_args()
    configs = args.configs or sorted(
        d for d in os.listdir(args.root)
        if os.path.isdir(os.path.join(args.root, d)))
    total = 0
    for prefix in configs:
        made = merge_config(prefix, args.methods, args.labels, args.root,
                            clean = args.clean)
        if made:
            print(f"[{prefix}] merged {len(made)} archives", flush = True)
            total += len(made)
    print(f"merged {total} archives under {args.root}")


if __name__ == "__main__":
    main()

# %%
