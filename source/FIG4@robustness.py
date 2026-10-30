# %%
from matplotlib.backends.backend_pdf import PdfPages
import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), "./library"))
sys.pycache_prefix = os.path.abspath("./output/pycache")
from pltools import GridPlt, parse_xlims
from corrmat import available_methods, observables, load_corrmat
from fsutils import load_config, ensure_dir


# Styling keyed by the four cases + the reference (shared by all groups).  The
# reference is the standard HEOM for linear coupling and the extended DEOM for
# polynomial coupling, so its label depends on the config (see ``labels``).
VARIANTS = {
    "qhm":  ("C3", "-", 1.5, 1.0),
    "qah":  ("C0", ":", 1.5, 1.0),
    "chm":  ("C4", "-", 1.8, 0.5),
    "cah":  ("C2", ":", 1.8, 0.5),
    "heom": ("k",  "-", 1.2, 1.0),
}

# Groups are the methods whose correlation matrix the data directory holds
# (``corrmat_<group>[_<label>].npz``, or the per-case archives); each repo draws
# only its own methods, and ``main`` is called once per group found.

# Fixed y-limits for the observables (the same for the matrix and scalar rows).
YLIMS = {"corr":  [-1, 1, 2],   "xorr":  [-1, 1, 2],
         "cohe+": [-.5, .5, 2], "cohe-": [-.5, .5, 2]}


def _load(d: str, ref_d: str, group: str, names: list,
          label: str = "") -> dict:
    """Load the observables that exist for one group, per case.

    Returns ``{name: {quantity: (t, y, ylim)}}`` (``pop`` as ``(t, y)``).
    Every group stores the same basis correlation matrix
    (``corrmat_<group>[_<label>].npz``); the reference ``corrmat_heom.npz`` is
    read from ``ref_d`` and drawn with every group.  ``corrmat.observables``
    contracts each matrix to the five curves.
    """
    out = {n: {} for n in names}

    def _fill(name, t, C):
        obs = observables(C)
        out[name]["xorr"] = (t, obs["xorr"], YLIMS["xorr"])
        out[name]["cohe+"] = (t, obs["cohe+"], YLIMS["cohe+"])
        out[name]["cohe-"] = (t, obs["cohe-"], YLIMS["cohe-"])
        out[name]["corr"] = (t, obs["corr"], YLIMS["corr"])
        out[name]["pop"] = (t, obs["pop"].real)

    if group != "heom":
        for n, (t, C) in load_corrmat(d, group, label = label).items():
            _fill(n, t, C)
    if "heom" in names:
        for n, (t, C) in load_corrmat(ref_d, "heom").items():
            _fill(n, t, C)
    return out


def main(prefix: str, group: str = "mkct", ref_dir: str = ".",
         label: str = ""):
    names  = list(VARIANTS)
    cases  = [n for n in names if n != "heom"]
    obs    = _load("./output/data/" + prefix + "/",
                   os.path.join(ref_dir, prefix) + "/", group, names,
                   label = label)

    # A quantity is drawn when every case has it (HEOM is optional, added when
    # present); each row is built only from the cases that carry the data.
    def have(key):
        return all(key in obs[n] for n in cases)
    if not any(have(q) for q in ("corr", "xorr", "cohe+", "cohe-")):
        return

    params = load_config(prefix)
    xlims  = parse_xlims(params["xlims"], reorder = True)
    # A two-window config (K) keeps the cut only for the population page; its
    # correlation/coherence pages use the standard single-window 4-grid over
    # [0, 24] with ticks every 6 (xlims[0][1]/xlims[1]*2 in config order).
    two_windows = len(xlims) == 4
    std_xlims   = [[0, 24], 4] if two_windows else xlims

    # Linear coupling is compared against the standard HEOM; the quadratic case
    # against the extended DEOM (dissipaton) hierarchy.
    ref_label = r"\ttfamily extDEOM" if params["quadratic"] else r"\ttfamily HEOM"
    labels = {n: (ref_label if n == "heom" else r"\ttfamily " + n.upper())
              for n in names}

    def row(key):
        out = []
        for n in names:
            if key not in obs[n]:
                continue
            t, y, ylim = obs[n][key]
            spec = (t, y, *VARIANTS[n])
            spec += (ylim,) if ylim is not None else ()
            out.append(spec + (labels[n],))
        return out

    # Population lines (P0 labelled, P1 unlabelled)
    pop_lines = None
    if have("pop"):
        pop_lines = []
        for n in names:
            if "pop" not in obs[n]:
                continue
            t, P = obs[n]["pop"]
            if n == "heom":
                color, ls, lw, alpha = VARIANTS[n]
                pop_lines.append((t, P[:, 0], color, ls, lw, alpha,
                                  labels[n]))
                pop_lines.append((t, P[:, 1], color, ls, lw, alpha / 2, ""))
                continue
            p, q = cases[cases.index(n) % 2], cases[2 + cases.index(n) % 2]
            color, ls = VARIANTS[n][0], VARIANTS[n][1]
            pop_lines.append((t, P[:, 0], color, ls,
                              VARIANTS[p][2], VARIANTS[p][3],
                              labels[n]))
            pop_lines.append((t, P[:, 1], color, ls,
                              VARIANTS[q][2], VARIANTS[q][3], ""))

    rows = [row(q) for q in ("corr", "xorr") if have(q)]
    biy  = [r"$\expval*{\sigma_x(t) \sigma_x(0)}$",
            r"$\expval*{\sigma_x(t) \sigma_y(0)}$"][:len(rows)]

    ensure_dir("./output/figures")
    name = group + "@" + prefix + ("_" + label if label else "")
    with PdfPages("./output/figures/" + name + ".pdf") as pdf:
        # Figures 1-3. Population (two-window cut) + autocorrelation &
        # cross-correlation (standard single-window 4-grid).
        if two_windows and pop_lines is not None:
            GridPlt.reim4gridplt(
                biylabel  = biy,
                supxlabel = "Time",
                xlims     = xlims,
                pop_lines = pop_lines,
                plot_data = [],
                pdf = pdf, close_fig = True
            )
        if rows:
            GridPlt.reim4gridplt(
                biylabel  = biy,
                supxlabel = "Time",
                xlims     = std_xlims,
                pop_lines = None if two_windows else pop_lines,
                plot_data = rows,
                pdf = pdf, close_fig = True
            )
        # Figure 4 & 5. Coherence ±
        if have("cohe+") and have("cohe-"):
            GridPlt.reim4gridplt(
                biylabel  = [r"$\expval*{\sigma_+(t)}$",
                             r"$\expval*{\sigma_-(t)}$"],
                supxlabel = "Time",
                xlims     = std_xlims,
                plot_data = [row("cohe+"), row("cohe-")],
                pdf = pdf, close_fig = True
            )


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
    import argparse
    ap = argparse.ArgumentParser(description = "Per-method comparison figures.")
    ap.add_argument("-ref", type = str, default = ".",
                    help = "directory holding the corrmat_heom.npz references")
    ap.add_argument("-label", type = str, default = "",
                    help = "tagged run suffix (default: the baseline)")
    ap.add_argument("-prefix", type = str, default = None,
                    help = "single config (default: every prefix of RUNLIST)")
    args = ap.parse_args()
    for prefix in ([args.prefix] if args.prefix else RUNLIST):
        d = "./output/data/" + prefix + "/"
        for group in available_methods(d, label = args.label):
            main(prefix = prefix, group = group, ref_dir = args.ref,
                 label = args.label)

# %%
