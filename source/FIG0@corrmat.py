# %%
"""Basis correlation matrix C_ij(t) of one method, one page per case.

Every method stores the same 4x4 matrix in the Liouville basis
``|0><0|, |0><1|, |1><0|, |1><1|`` (with ``C(0) = I``):
``corrmat_mkct.npz`` (matrix MKCT), ``corrmat_fssh.npz`` and
``corrmat_mash.npz`` (merged surface-hopping archives), one matrix per case
under the case name.  This runner draws the real (blue) and imaginary (red)
parts of every element as a 4x4 grid of panels, one page per case, in the
style of the tensorial-MKCT supporting information (the basis autocorrelation
function), and writes ``corrmat-<method>@<prefix>.pdf``.

Usage (from ``./source``, with the ``mkct`` environment active):
    python FIG0@corrmat.py                 # every method and prefix
    python FIG0@corrmat.py -method mash
    python FIG0@corrmat.py -prefix A-base_line
When ``corrmat_heom.npz`` exists (from ``python reference.py``), the ``qhm``
page of every method overlays the HEOM/extDEOM reference ``C(t)`` in dashed
black (the shared harmonic reference); the other cases stay method-only.
"""
import os, sys, argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
sys.path.append(os.path.join(os.path.dirname(__file__), "./library"))
sys.pycache_prefix = os.path.abspath("./output/pycache")
from pltools import GridPlt, parse_xlims
from fsutils import ensure_dir, load_config
from corrmat import available_methods, load_corrmat
from corrmat_page import page


# ``mkct.py`` basis_list: |0><0|, |0><1|, |1><0|, |1><1| (both axes of C_ij).
BASIS   = ("00", "01", "10", "11")
CASES   = ("qhm", "qah", "chm", "cah")
# The SI draws Re blue and Im red, both solid; the reference is dashed black.
RE_COLOR, IM_COLOR = "C0", "C3"
REF_COLOR = "k"
# Full-window ticks: five ticks (four intervals) from 0 to the output time.
NTICKS = 5
# Common panel y-range (ticks at -1/0/1, limits a hair beyond), as in the
# 4-grid bodies of the dynamics figures.
YLIM = [-1, 1, 2]


def main(prefix: str, ref_dir: str = "./output/data", label: str = "",
         tmax: float = None):
    """Draw ``corrmat-<method>@<prefix>.pdf``: one 4x4 page per case.

    ``ref_dir`` is a data root (the layout of ``./output/data``) whose
    ``<prefix>/corrmat_heom.npz`` is overlaid in dashed black on the ``qhm``
    page when present.  Every method whose matrix exists in the data directory
    is drawn (``corrmat.available_methods``, grouped or per-case archives);
    ``-label`` selects a tagged run.
    """
    d = "./output/data/" + prefix + "/"
    params = load_config(prefix)
    # The full computed window [0, tf] (not the cut of the dynamics figures);
    # ``tmax`` overrides it for runs that stop earlier (or go further).
    tf = (parse_xlims(params["xlims"], reorder = True)[0][1]
          if tmax is None else float(tmax))
    # Shared reference C(t) for the qhm page (harmonic bath only).
    ref, ref_label = None, ""
    ref_series = load_corrmat(os.path.join(ref_dir, prefix), "heom")
    if ref_series:
        ref = ref_series.get("qhm")
        ref_label = r"\ttfamily extDEOM" if params["quadratic"] \
            else r"\ttfamily HEOM"
    ensure_dir("./output/figures")
    for method in available_methods(d, label = label):
        series = load_corrmat(d, method, label = label)
        name = ("corrmat-" + method + "@" + prefix
                + ("_" + label if label else "")
                + (f"_t{tmax:g}" if tmax is not None else ""))
        with PdfPages("./output/figures/" + name + ".pdf") as pdf:
            for case in CASES:
                if case not in series:
                    continue
                t, C = series[case]
                fig = page(case, [(t, C, None, "-", "")], tf,
                           ref = ref if case == "qhm" else None,
                           ref_label = ref_label)
                pdf.savefig(fig)
                plt.close(fig)


# Define the running demo list
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
if __name__ == "__main__":
    ap = argparse.ArgumentParser(
        description = "Basis correlation matrix figures.")
    ap.add_argument("-ref", type = str, default = "./output/data",
                    help = "data root holding the corrmat_heom.npz references "
                           "(default: %(default)s)")
    ap.add_argument("-label", type = str, default = "",
                    help = "tagged run suffix (default: the baseline)")
    ap.add_argument("-prefix", type = str, default = None,
                    help = "single config (default: every prefix of RUNLIST)")
    ap.add_argument("-tmax", type = float, default = None,
                    help = "axis / window upper limit (default: the config tf)")
    args = ap.parse_args()
    for prefix in ([args.prefix] if args.prefix else RUNLIST):
        main(prefix = prefix, ref_dir = args.ref, label = args.label,
             tmax = args.tmax)

# %%
