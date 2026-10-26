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
"""
import numpy as np

from qutils import Pauli, basis, qubit, get_expval_ABrho0


# Labels of basis_list, both axes of C_ij.
BASIS      = ("00", "01", "10", "11")
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


def observables(C: np.ndarray) -> dict:
    """The five reported series of a basis correlation matrix.

    ``C`` is ``(Nt, 4, 4)``; returns ``{"pop": (Nt, 2), "cohe+": (Nt,),
    "cohe-": (Nt,), "corr": (Nt,), "xorr": (Nt,)}``, all complex, with the
    autocorrelation normalised to ``corr(0) = 1``.
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
    out["corr"] = out["corr"] / out["corr"][0]
    return out

# %%
