# %%
"""Reference solver: the basis correlation matrix of the harmonic reference.

Linear coupling is propagated with qutip's ``HEOMSolver`` and polynomial
(quadratic) coupling with the extended DEOM; both supply the quantum-harmonic
reference used by the MKCT, DMRG and surface-hopping comparisons.  Only the
correlation matrix is stored (``corrmat_heom.npz``, one ``(t, C)`` pair under
``qhm``, ``C(0) = I``); its contraction to the five reported observables is
done by the analysis scripts from ``library/corrmat.py``.

Usage (from ``./source``)::

    python reference.py               # every prefix of RUNLIST
    python reference.py -prefix A-base_line
"""
# 1. Requirements
import mpmath as mp
mp.mp.dps = 31
import numpy as np
import os, sys, argparse
sys.path.append(os.path.join(os.path.dirname(__file__), "./library"))
from qutip import expect, Qobj
from qutip.core.environment import OhmicEnvironment
from qutip.solver.heom import HEOMSolver
sys.pycache_prefix = os.path.abspath("./output/pycache")
from extended_deom import ExtendedDEOMSolver, build_ohmic_bath
from fsutils import load_config, save_group
from qutils import spin_boson, parse_coeffs
from corrmat import BASIS, basis_list
from pltools import parse_xlims


# 2. Common Configurations
def corrmat_ref(RefMats, tlist):
    """Basis correlation matrix ``C_ij(t)`` of the quantum harmonic case.

    ``C[:, j, k] = Tr[phi_j rho_k(t)]`` with the basis-unit initial states
    ``rho_k(0) = phi_k^dagger`` (so ``C(0) = I``), matching the ``(Nt, 4, 4)``
    matrix of ``mkct.py``.  The four basis units are propagated one at a time
    and all four operators measured on each run.
    """
    n = len(basis_list)
    C = np.zeros((len(tlist), n, n), dtype = np.complex128)
    A_ops = [Qobj(phi) for phi in basis_list]
    for k, phi in enumerate(basis_list):
        res = RefMats.run(Qobj(phi.conj().T), tlist)
        for j, A_j in enumerate(A_ops):
            C[:, j, k] = expect(A_j, res.states)
        del res
    return C


def main(prefix: str, corrmat: bool = False):
    # 3. Parameter
    ## 3.1 Loading parameters from YAML
    Lmax          = 6
    params        = load_config(prefix)
    linear        = params["linear"]
    quadratic     = params["quadratic"]
    rescale       = params["scale"]
    lambd         = params["lambd"]
    theta_B       = params["theta_B"]
    wB            = params["wB"]
    wC            = params["wC"]
    Delta         = params["Delta"]
    epsilon       = params["epsilon"]
    beta          = params["beta"]
    anharmon      = float(params["anharmon"] * mp.sqrt(7/6))
    xlims         = parse_xlims(params["xlims"], reorder = True)
    ## 4.2 Rescale the parameters to avoid numerical issues
    beta /= rescale
    wC, Delta, epsilon, lambd, wB, anharmon = [x * rescale \
                                for x in [wC, Delta, epsilon, lambd, wB, anharmon]]
    tf     = xlims[0][1]
    dt     = 1e-3
    tscale = 1/rescale
    tlist  = np.linspace(0, tf,  int(tf/dt))
    wlist  = np.linspace(0, 100, int(1e3))
    ## 3.2 System Hamiltonian, Interaction Operator
    Hs, V = spin_boson(
        model   = "QSBM" if quadratic else "LSBM",
        epsilon = epsilon, Delta = Delta, nohalf = True)
    # 4. Calculation
    ## 4.1 Setup
    #### 4.1.1 Environment & Decompose
    if quadratic:
        poly_coeffs = parse_coeffs(
            [0.5 * wB * (theta_B ** 2 - 1.0),
             -(2.0 * lambd * wB) ** 0.5 * theta_B ** 2,
             lambd * theta_B ** 2],
            linear = linear, quadratic = quadratic)
        a2, a1, a0 = float(poly_coeffs[0]), float(poly_coeffs[1]), \
                     float(poly_coeffs[2])
        tmax_solver = tf * tscale
        fit_tmax    = min(tmax_solver, 30.0)
        alpha       = np.pi * lambd / wC / rescale ** 2
        eta, gamma, q2, env, approx_env, fit_info = build_ohmic_bath(
            alpha = alpha,
            wc    = wC,
            beta  = beta,
            tmax  = fit_tmax,
            Nr    = 6,
        )
    else:
        env = OhmicEnvironment(
            T     = 1 / beta,
            alpha = np.pi * lambd / wC,
            wc    = wC,
            s     = 1.0,
        )
        approx_env, fit_info = env.approximate(
            method   = "espira-I",
            tlist    = wlist,
            Nr       = 4,
            Ni       = 4,
            separate = True,
        )
    #### 4.1.2 Set up the reference solver
    ####       - linear coupling: qutip's standard HEOM
    ####       - polynomial coupling: the extended DEOM (dissipaton) hierarchy
    if quadratic:
        default_options = {
            "rtol": 1e-9,
            "atol": 1e-11,
            "chunk": 200,
        }
        RefMats = ExtendedDEOMSolver(
            Qobj(Hs), Qobj(V), a0, a1, a2, eta, gamma, q2,
            max_depth = Lmax, options = default_options)
    else:
        default_options = {
            "nsteps": 1500,
            "store_states": True,
            "rtol": 1e-12,
            "atol": 1e-12,
            "method": "vern9",
            "progress_bar": "enhanced",
        }
        RefMats = HEOMSolver(
            Qobj(Hs), (approx_env, Qobj(V)),
            max_depth = Lmax, options = default_options)
    ## 4.2 Propagate the four basis-unit initial states phi_k^dagger
    C = corrmat_ref(RefMats, tlist * tscale)
    save_group({"qhm": (tlist, C)}, "./output/data/" + prefix,
               "corrmat_heom.npz", identifier = r"$\mathcal{C}(t)$",
               branch = BASIS, component = BASIS, compress = True)


RUNLIST  = [
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
    ap = argparse.ArgumentParser(description = "HEOM/extDEOM reference solver.")
    ap.add_argument("-prefix", type = str, default = None,
                    help = "single config (default: every prefix of RUNLIST)")
    args = ap.parse_args()
    for prefix in ([args.prefix] if args.prefix else RUNLIST):
        main(prefix = prefix)

# %%
