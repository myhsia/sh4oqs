# %%
# 1. Requirements
import os, sys, argparse, multiprocessing
from concurrent.futures import ProcessPoolExecutor
# Single-threaded BLAS in the process-parallel trajectory workers.
for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
             "VECLIB_MAXIMUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_var, "1")
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "./library"))
import paths
paths.install_cache_env()
if not paths.is_frozen():
    sys.pycache_prefix = os.path.abspath("./output/pycache")
from surface_hopping import mash, fssh
from qutils import ohmic_modes, parse_coeffs
from fsutils import load_config, save_group, resolve_configs, parse_xlims
from corrmat import BASIS


# Run selection — comment out entries to skip.
METHODS = ["fssh", "mash"]
CASES   = ["qhm", "qah", "chm", "cah"]

# Per-method trajectory budgets (used when ``-n_traj`` is not given).
N_TRAJ = {"mash": 40000, "fssh": 20000, "lsc": 40000, "ehrenfest": 20000}

# Mirror of the MKCT ``cases`` matrix: case -> (quantum sampling, anharmonic).
_CASE_MATRIX = {
    "qhm": (True,  False),
    "qah": (True,  True),
    "chm": (False, False),
    "cah": (False, True),
}


def _case_name(quantum: bool, anharm: bool) -> str:
    for name, (q, ah) in _CASE_MATRIX.items():
        if q == quantum and ah == anharm:
            return name
    raise ValueError(f"No case for quantum={quantum}, anharm={anharm}")


def _parse_args():
    parser = argparse.ArgumentParser(
        description = "Surface-hopping reference dynamics for the spin-boson model.")
    parser.add_argument("configs", nargs = "*", default = None,
                        help = "Configs: names, .yaml paths or glob patterns "
                               "such as 'config/*.yaml' (default: run RUNLIST)")
    parser.add_argument("-config-dir", type = str, default = None,
                        help = "where bare config names are looked up (default: "
                               "./config, else the config/ beside a frozen binary)")
    parser.add_argument("-method",
                        choices = ["mash", "fssh", "lsc", "ehrenfest", "all"],
                        default = "all", help = "Trajectory method")
    parser.add_argument("-n_traj", type = int, default = None,
                        help = "Number of trajectories (default: per-method budget)")
    parser.add_argument("-N", type = int, default = 200,
                        help = "Number of bath modes")
    parser.add_argument("-wmax_factor", type = float, default = 10.0,
                        help = "Quadrature cutoff, in units of wC")
    parser.add_argument("-dt", type = float, default = 0.01,
                        help = "Output time step (default: 0.01)")
    parser.add_argument("-seed", type = int, default = 0,
                        help = "RNG seed")
    parser.add_argument("-sampling", choices = ["auto", "quantum", "classical"],
                        default = "auto",
                        help = "Bath sampling (default: by case)")
    parser.add_argument("-label", type = str, default = "",
                        help = "Suffix appended to the output tag (sweeps)")
    parser.add_argument("-nproc", type = int, default = 4,
                        help = "Number of cases to run in parallel (default: %(default)s)")
    parser.add_argument("-njobs", type = int, default = -1,
                        help = "Inner workers per case (default: -1, share the cores)")
    parser.add_argument("-outdir", type = str, default = "./output/data",
                        help = "Root of the data tree (default: %(default)s)")
    parser.add_argument("-version", "--version", action = "version",
                        version = paths.version_string("sh"),
                        help = "print the build identity and exit")
    return parser.parse_args()


def _compute(source: str, method: str = "mash", quantum: bool = True,
             anharm: bool = False, n_traj: int = 40000, N: int = 200,
             wmax_factor: float = 10.0, dt: float = 0.01, seed: int = 0,
             n_jobs: int = 1, label: str = "", config_dir: str = None):
    """Compute one ``(config, method, case)``; return ``(t, tag, case, C)``."""
    # 3. Parameters (physical values; ``scale`` is a numerical device only)
    name      = os.path.splitext(os.path.basename(source))[0]
    params    = load_config(source, config_dir)
    epsilon   = params["epsilon"]
    Delta     = params["Delta"]
    lambd     = params["lambd"]
    wC        = params["wC"]
    beta      = params["beta"]
    quadratic = params["quadratic"]
    model     = None
    if quadratic:
        # QSBM: H_sys = 2 epsilon |1><1| + Delta sigma_x, V = |1><1|,
        # U(q) = a0 + a1 Q + a2 Q^2, so Vbar = epsilon + U/2 and A = -Vbar.
        theta_B = params["theta_B"]
        wB      = params["wB"]
        alp0    = lambd * theta_B ** 2
        alp1    = -(2.0 * lambd * wB) ** 0.5 * theta_B ** 2
        alp2    = 0.5 * wB * (theta_B ** 2 - 1.0)
        a2, a1, a0 = parse_coeffs([alp2, alp1, alp0],
                                  linear = params["linear"], quadratic = True)
        pv0, pv1, pv2 = epsilon + a0 / 2.0, a1 / 2.0, a2 / 2.0
        model = (pv0, pv1, pv2, -pv0, -pv1, -pv2)
    a = float(params["anharmon"]) * np.sqrt(7.0 / 6.0) if anharm else 0.0
    xlims = parse_xlims(params["xlims"], reorder = True)
    tf    = xlims[0][1]
    omega, c = ohmic_modes(lambd, wC, N = N, wmax_factor = wmax_factor)
    # The quartic stiffening raises the fastest effective frequency, so an
    # anharmonic run needs a shorter step than a harmonic one.
    if a > 0.0:
        dt = min(dt, 0.4 / float(np.max(omega)))
    t     = np.linspace(0.0, tf, int(round(tf / dt)) + 1)

    case = _case_name(quantum, anharm)
    tag = f"{method}_{case}"
    if label:
        tag += "_" + label
    print(f"[{name}] {tag}: N={N} n_traj={n_traj} tf={tf:g} dt={dt:g} "
          f"beta={beta} wC={wC} lambd={lambd} a={a:.4f} n_jobs={n_jobs}",
          flush = True)

    # 4. Dynamics: the basis correlation matrix (the observables are its
    # contractions, rebuilt by the analysis scripts).
    if method in ("mash", "lsc"):
        # "lsc" is the same mapping dynamics with the mean-field (continuous Sz)
        # force instead of deterministic hops.
        C = mash.run_corrmat(epsilon, Delta, omega, c, beta, t, quantum = quantum,
                             anharmon = a, n_traj = n_traj, seed = seed,
                             model = model, n_jobs = n_jobs,
                             meanfield = (method == "lsc"))
    else:
        # "ehrenfest" is FSSH with the mean-field force and no hops.
        C = fssh.run_corrmat(epsilon, Delta, omega, c, beta, t,
                             quantum = quantum, anharmon = a,
                             n_traj = n_traj, seed = seed,
                             model = model, n_jobs = n_jobs,
                             meanfield = (method == "ehrenfest"))
    return t, tag, case, C


def _save(name: str, method: str, results: dict, label: str = "",
          outdir: str = "./output/data"):
    """Write one archive per method, with one matrix per case.

    Only the parent writes, and all four cases are in hand, so the file has the
    same layout as ``corrmat_mkct.npz`` (``t_<case>`` plus ``<case>``) and no
    post-processing merge is needed.
    """
    tag = method + (f"_{label}" if label else "")
    save_group({c: results[c] for c in results}, os.path.join(outdir, name),
               "corrmat_" + tag + ".npz", identifier = r"$\mathcal{C}(t)$",
               branch = BASIS, component = BASIS, compress = True)


def _case_job(spec):
    """Worker entry point: one ``(prefix, method, case)`` specification."""
    t, tag, case, C = _compute(*spec)
    return t, tag, case, C


RUNLIST = [
    "A-base_line",
    # "B-slow_system",
    # "C-no_bias",
    # "D-strong_tunneling",
    # "E-fast_bath",
    # "F-slow_bath",
    # "G-lowT",
    # "H-highT",
    # "I-weak_coupling",
    # "J-strong_coupling",
    # "K-linear_quad",
]
def run(source: str, methods = METHODS, cases = CASES, n_traj = None,
        N: int = 200, wmax_factor: float = 10.0, dt: float = 0.01,
        seed: int = 0, nproc: int = 4, njobs: int = -1,
        sampling: str = "auto", label: str = "",
        outdir: str = "./output/data", name: str = None,
        config_dir: str = None):
    """Run one or more methods across the selected cases for a single config.

    ``n_traj`` defaults to the per-method ``N_TRAJ`` budget.  The cases run as
    ``nproc`` parallel processes (at most one per case) and, when there are
    spare cores, each case splits its trajectories over ``njobs`` inner
    workers (``-1``: share the cores over the cases).  The parent process
    writes every archive under ``outdir``.
    """
    name    = name or os.path.splitext(os.path.basename(source))[0]
    nproc   = max(1, int(nproc))
    workers = max(1, min(nproc, len(cases)))
    if njobs is not None and int(njobs) > 0:
        n_inner = int(njobs)
    elif workers >= len(cases) and len(cases) > 1:
        n_inner = 1          # nested pools only add overhead here
    else:
        n_inner = max(1, paths.n_cpus() // workers)
    for method in methods:
        nt = N_TRAJ[method] if n_traj is None else n_traj
        specs = []
        for case in cases:
            quantum, anharm = _CASE_MATRIX[case]
            if sampling != "auto":
                quantum = sampling == "quantum"
            specs.append((source, method, quantum, anharm, nt, N, wmax_factor,
                          dt, seed, n_inner, label, config_dir))
        # Collect every case in the parent, then write one archive per method.
        results = {}
        if workers == 1 or len(specs) == 1:
            for spec in specs:
                t, _tag, case, C = _compute(*spec)
                results[case] = (t, C)
        else:
            with ProcessPoolExecutor(max_workers = workers,
                    mp_context = multiprocessing.get_context("fork")) as pool:
                for t, _tag, case, C in pool.map(_case_job, specs):
                    results[case] = (t, C)
        _save(name, method, results, label, outdir)


if __name__ == "__main__":
    args = _parse_args()
    methods = METHODS if args.method == "all" else [args.method]
    try:
        configs = resolve_configs(args.configs or RUNLIST, args.config_dir)
    except FileNotFoundError as exc:
        raise SystemExit(f"error: {exc}")
    for name, path in configs:
        run(path, name = name, config_dir = args.config_dir,
            methods = methods, cases = CASES, n_traj = args.n_traj,
            N = args.N, wmax_factor = args.wmax_factor, dt = args.dt,
            seed = args.seed, nproc = args.nproc, njobs = args.njobs,
            sampling = args.sampling, label = args.label, outdir = args.outdir)

# %%
