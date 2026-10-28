# %%
"""The basis correlation matrix ``C(t)`` and the observables it carries.

The matrix MKCT (``mkct.py``) and the reference solver (``reference.py``) store
only the correlation matrix in the Liouville basis
``|0><0|, |0><1|, |1><0|, |1><1|``, with ``C(0) = I``, under the keys
``corrmat_mkct.npz`` / ``corrmat_heom.npz``.  Every reported observable is a
linear contraction of its elements with the operators below (via
``qutils.get_expval_ABrho0``), so the archive alone is enough and the plotting
and comparison scripts rebuild the five curves on the fly.

Both axes of ``C`` are labelled by ``BASIS`` and indexed by ``basis_list``
(``phi_i = |i><j|``).  The measured operators ``A`` and the time-zero states
``B rho(0)`` give ``<A(t) B(0)>``; the autocorrelation is normalised to one at
``t = 0``.

``load_corrmat`` reads one method's matrices: the grouped archive
``corrmat_<method>[_<label>].npz`` (one matrix per case) is canonical, and when
it is absent the per-case archives ``corrmat_<method>_<case>[_<label>].npz``
(CI jobs run one case each, segmented tDMRG runs are per case) are assembled.
``mkct`` and ``heom`` are always grouped.
"""
import os

import numpy as np

from qutils import Pauli, basis, qubit, get_expval_ABrho0
from fsutils import load_group


# Labels of basis_list, both axes of C_ij.
BASIS      = ("00", "01", "10", "11")
# Report cases, in the archive order.
CASES      = ("qhm", "qah", "chm", "cah")
# Basis: |0><0|, |0><1|, |1><0|, |1><1|
basis_list = np.array([basis(2, (0, 0)), basis(2, (0, 1)),
                       basis(2, (1, 0)), basis(2, (1, 1))])
X, Y       = Pauli("x"), Pauli("y")
P_list     = [basis(2, 0, kind = "ketbra"),
              basis(2, 1, kind = "ketbra")]
## Initial system density matrices
rho0       = {"0": basis(2, 0, kind = "ketbra"),
              "x": qubit(theta = 2 * np.pi/3, kind = "ketbra")}
## Measured operator A and time-zero state Brho0  ->  <A(t) B(0)>
A          = {"pop": P_list, "+": Pauli("+"), "-": Pauli("-"),
              "corr": X,     "xorr": X}
Brho0      = {"pop": rho0["0"], "+": rho0["0"], "-": rho0["0"],
              "corr": X @ rho0["0"], "xorr": Y @ rho0["x"]}


def observables(C: np.ndarray, normalize: bool = True) -> dict:
    """The five reported series of a basis correlation matrix.

    ``C`` is ``(Nt, 4, 4)``; returns ``{"pop": (Nt, 2), "cohe+": (Nt,),
    "cohe-": (Nt,), "corr": (Nt,), "xorr": (Nt,)}``, all complex.  With
    ``normalize`` (the default) the autocorrelation is scaled to
    ``corr(0) = 1``; a segmented propagation passes ``normalize=False`` and
    rescales once for the whole series.
    """
    C = np.asarray(C)
    out = {
        "pop": np.stack(
            [get_expval_ABrho0(C, A["pop"][i], Brho0["pop"], basis_list)
             for i in range(2)], axis = 1),
        "cohe+": get_expval_ABrho0(C, A["+"], Brho0["+"], basis_list),
        "cohe-": get_expval_ABrho0(C, A["-"], Brho0["-"], basis_list),
        "corr": get_expval_ABrho0(C, A["corr"], Brho0["corr"], basis_list),
        "xorr": get_expval_ABrho0(C, A["xorr"], Brho0["xorr"], basis_list),
    }
    if normalize and out["corr"][0] != 0:
        out["corr"] = out["corr"] / out["corr"][0]
    return out


def load_corrmat(directory, method: str, case: str = None,
                 label: str = "") -> dict:
    """Load one method's basis correlation matrices as ``{case: (t, C)}``.

    The grouped ``corrmat_<method>[_<label>].npz`` is canonical; when absent,
    the per-case ``corrmat_<method>_<case>[_<label>].npz`` archives are
    assembled (per-case CI jobs and segmented tDMRG runs).  ``case`` restricts
    the result to one key; an absent archive gives an empty dict.
    """
    d = str(directory).rstrip("/") + "/"
    tag = method + ("_" + label if label else "")
    grouped = d + "corrmat_" + tag + ".npz"
    if os.path.exists(grouped):
        series, _ = load_group(grouped)
        out = dict(series)
    else:
        out = {}
        for c in CASES:
            candidates = (d + "corrmat_" + tag + "_" + c + ".npz",
                          d + "corrmat_" + method + "_" + c +
                          ("_" + label if label else "") + ".npz")
            path = next((p for p in candidates if os.path.exists(p)), None)
            if path is None:
                continue
            series, _ = load_group(path)
            if c in series:
                out[c] = series[c]
    if case is not None:
        return {case: out[case]} if case in out else {}
    return out

# %%
