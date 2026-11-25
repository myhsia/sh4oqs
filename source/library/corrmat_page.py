# %%
"""The 4x4 basis-correlation page: one figure specification, two callers.

``FIG0@corrmat.py`` draws one method's ``C_ij(t)`` grid, and
``FIG9@bary_depth.py`` draws the *same* grid with one curve per moment depth.  The
layout is therefore defined once, here, and both figures import it -- a copy would
let the two drift apart at the first style change.

The page is the ``mkct.py`` basis grid, ``C_ij(t)`` for the four Liouville basis
elements, with the real part blue and the imaginary part red (the SI convention),
the reference matrix dashed black, and the element tag in the south-east corner of
every panel.

:func:`page` takes a **list** of ``(t, C, color, dash, label)`` series, so one
method is a one-element list and four depths is a four-element one.  A caller that
draws only the real part passes ``parts = ("re",)``.
"""
# %%
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from pltools import GridPlt

#: The Liouville basis of ``corrmat``: |0><0|, |0><1|, |1><0|, |1><1|.
BASIS = ("00", "01", "10", "11")
#: The SI draws Re blue and Im red, both solid; the reference is dashed black.
RE_COLOR, IM_COLOR = "C0", "C3"
REF_COLOR = "k"
#: Full-window ticks: five ticks (four intervals) from 0 to the output time.
NTICKS = 5
#: Common panel y-range (ticks at -1/0/1, limits a hair beyond).
YLIM = [-1, 1, 2]
#: Opacity of the imaginary part when both parts share one page.
IM_ALPHA = 0.45


def page(case: str, series: list, tf: float, parts = ("re", "im"),
         ref = None, ref_label: str = "", tag: str = None, note: str = "",
         legend: bool = True, fade_imag: bool = False):
    """One 4x4 page of ``C[:, i, j]`` over any number of series.

    ``series`` is a list of ``(t, C, color, linestyle, label)``; ``C`` is
    ``(nt, 4, 4)`` on the basis grid and ``t`` its own axis (the caller passes one
    grid per pair of axes).  ``parts`` selects which of the real and imaginary
    parts are drawn -- the per-depth figure passes ``("re",)`` and ``("im",)`` on
    two pages, because four series of both parts is eight lines per panel.

    ``ref = (t_ref, C_ref)`` overlays the reference matrix in dashed black.
    ``tag`` replaces the page's case label (default: the upper-cased case).
    ``note`` is prepended to the key as a handle-less entry: it is where the
    *fixed* variable of a controlled comparison is stated, inside the figure.
    ``legend = False`` draws no key at all, for a figure whose key is its own
    legend-bar page (``pltools.GridPlt._legend_bar``), and then ``tag`` is the only
    thing the page carries on that line.

    ``fade_imag`` draws both parts on one page: the real part solid and the
    imaginary part in the same colour at ``IM_ALPHA``, so a controlled comparison
    needs one body page instead of two.  The convention (solid Re, faded Im) is not
    a label on the page -- it belongs in the caption.
    """
    n = len(BASIS)
    for item in series:
        C = item[1]
        if C.ndim != 3 or C.shape[1] != n or C.shape[2] != n:
            raise ValueError(
                "the correlation-matrix archive must hold (Nt, %d, %d) arrays; "
                "got %s." % (n, n, C.shape))
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
        gridspec_kw = {"left": left, "right": right, "top": top,
                       "bottom": bottom, "wspace": gap, "hspace": gap})
    axs = np.atleast_2d(axs)
    xticks = np.linspace(0.0, tf, NTICKS)
    badges = []
    for i in range(n):
        for j in range(n):
            ax = axs[i, j]
            ax.set_box_aspect(1)
            si = len(series) == 1 and len(parts) == 2
            im_alpha = IM_ALPHA if fade_imag else 1.0
            for (t_s, C_s, color, ls, _label) in series:
                if "re" in parts:
                    ax.plot(t_s, C_s[:, i, j].real, ls = ls, lw = 1.2,
                            color = RE_COLOR if si else color)
                if "im" in parts:
                    ax.plot(t_s, C_s[:, i, j].imag, ls = ls, lw = 1.2,
                            alpha = im_alpha,
                            color = IM_COLOR if si else color)
            if ref is not None:
                t_ref, C_ref = ref
                if "re" in parts:
                    ax.plot(t_ref, C_ref[:, i, j].real, color = REF_COLOR,
                            ls = "--", lw = 1.2)
                if "im" in parts:
                    ax.plot(t_ref, C_ref[:, i, j].imag, color = REF_COLOR,
                            ls = "--", lw = 1.2,
                            alpha = IM_ALPHA if fade_imag else 1.0)
            # Same y-formatting as GridPlt._style_axes, but without its y=0
            # guide line (the matrix elements already oscillate about zero).
            ax.ticklabel_format(style = "sci", scilimits = (-2, 2),
                                axis = "y", useMathText = True)
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
    # Page tag (case) on the left, the key centred on the same line.  When the
    # caller states a control in ``note``, the two would collide at this width, so
    # the tag joins the key as its first entry instead.
    if tag is not None or not legend:
        fig.text(0.09, 0.965, r"\ttfamily " + (tag or case.upper()),
                 fontsize = 12, ha = "left", va = "center")
    if not legend:
        fig.supxlabel("Time", y = 0.01, fontsize = 12)
        return fig
    handles, labels = [], []
    if len(series) == 1 and len(parts) == 2 and "re" in parts and "im" in parts:
        # One method: keep the SI identity (blue Re, red Im) rather than a key
        # per series, exactly as the single-method page has always read.
        handles += [Line2D([], [], color = RE_COLOR, lw = 1.2),
                    Line2D([], [], color = IM_COLOR, lw = 1.2)]
        labels += [r"$\Re$", r"$\Im$"]
    else:
        for (_t, _C, color, ls, label) in series:
            handles.append(Line2D([], [], color = color, ls = ls, lw = 1.2))
            labels.append(label)

    if ref is not None:
        handles.append(Line2D([], [], color = REF_COLOR, ls = "--", lw = 1.2))
        labels.append(ref_label)
    if note:
        handles.insert(0, Line2D([], [], color = "none", lw = 0))
        labels.insert(0, (tag + " \u00b7 " if tag else "") + note)
    fig.legend(handles, labels, loc = "center",
               bbox_to_anchor = (0.5, 0.965), ncol = len(labels),
               frameon = False, fontsize = 9, handlelength = 1.4,
               columnspacing = 1.2)
    fig.supxlabel("Time", y = 0.01, fontsize = 12)
    return fig

# %%
