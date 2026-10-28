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
from corrmat import load_corrmat


# ``mkct.py`` basis_list: |0><0|, |0><1|, |1><0|, |1><1| (both axes of C_ij).
BASIS   = ("00", "01", "10", "11")
CASES   = ("qhm", "qah", "chm", "cah")
METHODS = ("mkct", "fssh", "mash", "dmrg")
# The SI draws Re blue and Im red, both solid; the reference is dashed black.
RE_COLOR, IM_COLOR = "C0", "C3"
REF_COLOR = "k"
# Full-window ticks: five ticks (four intervals) from 0 to the output time.
NTICKS = 5
# Common panel y-range (ticks at -1/0/1, limits a hair beyond), as in the
# 4-grid bodies of the dynamics figures.
YLIM = [-1, 1, 2]


def _page(case: str, t: np.ndarray, C: np.ndarray, tf: float,
          ref = None, ref_label: str = ""):
    """One 4x4 page of ``C[:, i, j]``: real (blue) and imaginary (red).

    ``ref = (t_ref, C_ref)`` optionally overlays the reference matrix in
    dashed black on the same panels (only passed for the ``qhm`` page).
    """
    n = len(BASIS)
    if C.ndim != 3 or C.shape[1] != n or C.shape[2] != n:
        raise ValueError(
            "corrmat_mkct.npz must hold (Nt, %d, %d) arrays; got %s."
            % (n, n, C.shape))
    # Square panels: the height follows the cell width, with the same tight
    # gaps as the combined bodies (wspace = hspace = 0.1).  The shared y-axis
    # (identical -1/0/1 ticks everywhere) keeps tick labels on the left column
    # only, so the narrow horizontal gap does not collide with the numbers.
    left, right, top, bottom, gap = 0.09, 0.93, 0.93, 0.10, 0.1
    cell = GridPlt.W * (right - left) / (n + (n - 1) * gap)
    height = cell * (n + (n - 1) * gap) / (top - bottom)
    fig, axs = plt.subplots(
        n, n, sharex = True, sharey = True,
        figsize = (GridPlt.W, height),
        gridspec_kw = {'left': left, 'right': right, 'top': top,
                       'bottom': bottom, 'wspace': gap, 'hspace': gap})
    axs = np.atleast_2d(axs)
    xticks = np.linspace(0.0, tf, NTICKS)
    badges = []
    for i in range(n):
        for j in range(n):
            ax = axs[i, j]
            ax.set_box_aspect(1)
            ax.plot(t, C[:, i, j].real, color = RE_COLOR, lw = 1.2)
            ax.plot(t, C[:, i, j].imag, color = IM_COLOR, lw = 1.2)
            if ref is not None:
                t_ref, C_ref = ref
                ax.plot(t_ref, C_ref[:, i, j].real, color = REF_COLOR,
                        ls = "--", lw = 1.2)
                ax.plot(t_ref, C_ref[:, i, j].imag, color = REF_COLOR,
                        ls = "--", lw = 1.2)
            # Same y-formatting as GridPlt._style_axes, but without its y=0
            # guide line (the matrix elements already oscillate about zero).
            ax.ticklabel_format(style = 'sci', scilimits = (-2, 2),
                                axis = 'y', useMathText = True)
            GridPlt._normalize(ax, YLIM)
            ax.set_xlim(0.0, tf)
            ax.set_xticks(xticks)
            GridPlt._anchor_end_xticks(ax, 0.0, tf)
            # Element tag in the south-east corner, in the Re/Im rounded box.
            txt = ax.text(0.96, 0.04, rf"$\mathcal{{C}}_{{{i},{j}}}(t)$",
                          transform = ax.transAxes, ha = "right", va = "bottom",
                          fontsize = 8, bbox = GridPlt.TAG_BBOX)
            badges.append((txt, ax, "SE"))
            ax.tick_params(labelsize = 8)
            if j > 0:
                # Only the left column carries y ticks (the scale is shared).
                ax.tick_params(axis = "y", left = False, labelleft = False)
            if i < n - 1:
                ax.tick_params(labelbottom = False, bottom = False)
    GridPlt._align_corner_badges(fig, badges)
    # Page tag (case) on the left, the Re/Im key centred on the same line.
    fig.text(0.09, 0.965, r"\ttfamily " + case.upper(),
             fontsize = 12, ha = "left", va = "center")
    handles = [Line2D([], [], color = RE_COLOR, lw = 1.2),
               Line2D([], [], color = IM_COLOR, lw = 1.2)]
    labels = [r"$\Re$", r"$\Im$"]
    if ref is not None:
        handles.append(Line2D([], [], color = REF_COLOR, ls = "--", lw = 1.2))
        labels.append(ref_label)
    fig.legend(handles, labels, loc = "center",
               bbox_to_anchor = (0.5, 0.965), ncol = len(labels),
               frameon = False, fontsize = 10, handlelength = 1.6,
               columnspacing = 1.2)
    fig.supxlabel("Time", y = 0.01, fontsize = 12)
    return fig


def main(prefix: str, methods = METHODS, ref_dir: str = ".",
         label: str = ""):
    """Draw ``corrmat-<method>@<prefix>.pdf``: one 4x4 page per case.

    ``ref_dir`` holds the shared reference (``corrmat_heom.npz``), overlaid in
    dashed black on the ``qhm`` page when present; a method without its matrix
    is skipped.  Grouped and per-case archives are both accepted
    (``corrmat.load_corrmat``).
    """
    d = "./output/data/" + prefix + "/"
    params = load_config(prefix)
    # The full computed window [0, tf] (not the cut of the dynamics figures).
    tf = parse_xlims(params["xlims"], reorder = True)[0][1]
    # Shared reference C(t) for the qhm page (harmonic bath only).
    ref, ref_label = None, ""
    ref_series = load_corrmat(os.path.join(ref_dir, prefix), "heom")
    if ref_series:
        ref = ref_series.get("qhm")
        ref_label = r"\ttfamily extDEOM" if params["quadratic"] \
            else r"\ttfamily HEOM"
    ensure_dir("./output/figures")
    for method in methods:
        series = load_corrmat(d, method, label = label)
        if not series:
            print("skip " + prefix + ": no corrmat for " + method)
            continue
        name = "corrmat-" + method + "@" + prefix + ("_" + label if label else "")
        with PdfPages("./output/figures/" + name + ".pdf") as pdf:
            for case in CASES:
                if case not in series:
                    continue
                t, C = series[case]
                fig = _page(case, t, C, tf,
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
    ap.add_argument("-method", nargs = "+", default = ["all"],
                    choices = ["all"] + list(METHODS),
                    help = "methods to draw (default: all)")
    ap.add_argument("-ref", type = str, default = ".",
                    help = "directory holding the corrmat_heom.npz references")
    ap.add_argument("-label", type = str, default = "",
                    help = "tagged run suffix (default: the baseline)")
    ap.add_argument("-prefix", type = str, default = None,
                    help = "single config (default: every prefix of RUNLIST)")
    args = ap.parse_args()
    methods = METHODS if "all" in args.method else args.method
    for prefix in ([args.prefix] if args.prefix else RUNLIST):
        main(prefix = prefix, methods = methods, ref_dir = args.ref,
             label = args.label)

# %%
