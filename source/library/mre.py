"""Mean relative error of MKCT against HEOM, and the Padé-order scan core.

Shared tool (parallel to ``pltools.py``): the definition of the error, the
per-order scan, and the markdown/JSON reports.  The figure helpers
(``matrix_panel``, ``log_range``) live in ``pltools.GridPlt``.

For each config, scan all Padé orders (m, n) with m, n >= 2 and m + n <= N - 2
(N = number of moments), reconstruct the basis correlation matrix ``C(t)`` with
``MKCT_solver.solve_pade``, and score it against the HEOM reference matrix with
the global relative Frobenius error

    E(m, n) = ||C - C'||_F / ||C'||_F,

the Frobenius norm running over the full ``[0, tf]`` grid and all 16 complex
elements (so there is no per-element division by a vanishing series), reported
in percent.  Every pair is ranked, and the script writes

    output/pade_scan_results.json   all pairs (infeasible pairs as ``null``),
                                    each row carrying its metric window ``tf``
    output/pade_scan.md             the worst-to-best ranking per config, above
                                    ``MAX_ERR_PCT`` (the JSON keeps all pairs)

The driver (``source/pade_scan.py``) owns the run list, the solver settings and
the output paths::

    python pade_scan.py            # scan (or reuse cache)
    python pade_scan.py --rescan   # force a fresh scan
"""

import contextlib, io, json, os, time, warnings
import numpy as np

from fsutils import (load_config, load_group, ensure_dir, moments_filename,
                     load_moments)
from pltools import parse_xlims


MAX_ERR_PCT = 1e3    # default cutoff for the markdown tables (%)


def _rnorm(y):
    return float(np.sqrt(np.sum(np.abs(y) ** 2)))


def _interp(yh, th, t):
    """Resample the HEOM reference ``yh(th)`` onto the MKCT times ``t``."""
    return np.interp(t, th, yh.real) + 1j * np.interp(t, th, yh.imag)


def _interp_mat(C, th, t):
    """Resample a ``(Nt, d, d)`` complex series element by element."""
    flat = np.asarray(C).reshape(len(th), -1)
    out = np.empty((len(t), flat.shape[1]), dtype = np.complex128)
    for k in range(flat.shape[1]):
        out[:, k] = _interp(flat[:, k], th, t)
    return out.reshape((len(t),) + np.shape(C)[1:])


def _corrmat_error(C, C_ref):
    """Global relative Frobenius error ``||C - C_ref||_F / ||C_ref||_F``."""
    den = _rnorm(C_ref)
    return _rnorm(C - C_ref) / den if den > 0 else 0.0


def scan_config(prefix, dt=1e-2, tf_factor=20.0, solve=None,
                data_dir="./output/data", method="block_sp"):
    """Return ``[{'m','n','Nmax','E','tf'}]`` over all admissible Padé orders.

    ``E`` is the global relative Frobenius error in percent
    (``_corrmat_error``); ``solve`` holds the ``solve_pade`` keyword
    arguments other than ``tf``/``tfout``/``dt``.
    """
    # Imported lazily: the error helpers also serve repositories that do
    # not ship the MKCT solver (e.g. the surface-hopping reference).
    from mkct_pade.mkct_solver import MKCT_solver
    solve = dict(solve or {})
    params = load_config(prefix)
    rescale = params["scale"]
    tf = parse_xlims(params["xlims"], reorder=True)[0][1]
    d = os.path.join(data_dir, prefix) + os.sep
    mom = load_moments(d + moments_filename(method, "matrix"))["qhm"]
    N = mom.shape[0]
    Mmax = N - 2
    th, Ch = load_group(d + "corrmat_heom.npz")[0]["qhm"]
    if th[-1] > tf * (1.0 + 1e-6):
        warnings.warn(
            f"[{prefix}] HEOM reference reaches t={th[-1]:g} but the metric "
            f"window is [0, {tf:g}] (config xlims); the extra range is ignored.",
            stacklevel=2)
    # Keep only the calculation window: numerator and denominator both use it.
    win = th <= tf * (1.0 + 1e-9)
    thw, Chw = th[win], Ch[win]

    solver = MKCT_solver.init(mom, rescale=rescale)
    res = []
    t0 = time.time()
    for m in range(2, Mmax - 1):
        for n in range(2, Mmax - m + 1):
            try:
                with contextlib.redirect_stdout(io.StringIO()):
                    t, Ct = solver.solve_pade(
                        tf=tf_factor * tf, tfout=tf, dt=dt, kernel_order=1,
                        pade_order=(m, n), **solve)
                E = 100.0 * _corrmat_error(Ct, _interp_mat(Chw, thw, t))
                if not np.isfinite(E):
                    E = np.inf
            except Exception:
                E = np.inf
            res.append({"m": m, "n": n, "Nmax": m + n + 2, "E": E,
                        "tf": float(tf)})
    feas = sorted([r for r in res if np.isfinite(r["E"])], key=lambda r: r["E"])
    print(f"[{prefix}] N={N} pairs={len(res)} feasible={len(feas)} "
          f"time={time.time()-t0:.0f}s argmin=({feas[0]['m']},{feas[0]['n']}) "
          f"N{feas[0]['Nmax']} E={feas[0]['E']:.4e}", flush=True)
    return res


def write_markdown(raw, path, max_pct=MAX_ERR_PCT, runlist=None):
    """Worst-to-best ranking of the feasible pairs, above the error cutoff."""
    runlist = list(runlist if runlist is not None else raw)
    lines = [
        "# Pade-order scan (QHM vs HEOM)", "",
        "Error metric: the global relative Frobenius error "
        "`E = ||C - C'||_F / ||C'||_F` over the full `[0, tf]` (the config's "
        "window from `xlims`, listed per config) and all 16 complex elements "
        "of the 4x4 basis correlation matrix, QHM vs HEOM; the tables report "
        "`100*E` in percent.",
        "`Nmax = m + n + 2`; pairs satisfy `m, n >= 2` and singular/infeasible "
        f"pairs are excluded. Pairs with an error above `{max_pct:g}%` are "
        "omitted from the tables (the JSON keeps them); each config lists the "
        "remaining pairs ordered from **largest (#1)** to **smallest** error.", "",
        "| config | best (m, n) | Nmax | error (%) | tf | feasible pairs |",
        "|---|---|---:|---:|---:|---:|",
    ]
    for prefix in runlist:
        feas = sorted([r for r in raw[prefix] if np.isfinite(r["E"])],
                      key=lambda r: r["E"])
        b = feas[0]
        lines.append(f"| {prefix} | ({b['m']}, {b['n']}) | {b['Nmax']} | "
                     f"{b['E']:.4g} | {b['tf']:g} | {len(feas)} |")
    for prefix in runlist:
        feas = sorted([r for r in raw[prefix] if np.isfinite(r["E"])],
                      key=lambda r: r["E"])
        shown = [r for r in feas if r["E"] <= max_pct]
        lines += ["", f"## {prefix} ({len(shown)} of {len(feas)} shown)", "",
                  "| # | m | n | Nmax | error (%) |",
                  "|---:|---:|---:|---:|---:|"]
        for i, r in enumerate(reversed(shown), start=1):
            lines.append(f"| {i} | {r['m']} | {r['n']} | {r['Nmax']} | "
                         f"{r['E']:.4g} |")
    ensure_dir(os.path.dirname(path))
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print("WROTE", path, flush=True)


def _to_jsonable(raw):
    """Replace non-finite ``E`` (infeasible) with JSON ``null``."""
    out = {}
    for p, rows in raw.items():
        out[p] = []
        for r in rows:
            r = dict(r)
            if not np.isfinite(r["E"]):
                r["E"] = None
            out[p].append(r)
    return out


def _from_jsonable(raw):
    """Restore JSON ``null`` back to ``float('inf')`` for the report logic."""
    for rows in raw.values():
        for r in rows:
            if r["E"] is None:
                r["E"] = float("inf")
    return raw


def run(runlist, json_path, md_path, *,
        max_pct=MAX_ERR_PCT, rescan=False, **scan_kwargs):
    """Scan ``runlist`` (or load the JSON cache), then write both artifacts."""
    ensure_dir(os.path.dirname(json_path))
    if rescan or not os.path.exists(json_path):
        raw = {p: scan_config(p, **scan_kwargs) for p in runlist}
        print("scanned", len(runlist), "configs", flush=True)
    else:
        with open(json_path) as f:
            raw = _from_jsonable(json.load(f))
        print("loaded", json_path, flush=True)
    with open(json_path, "w") as f:
        json.dump(_to_jsonable(raw), f, indent=2)
    print("WROTE", json_path, flush=True)
    write_markdown(raw, md_path, max_pct=max_pct, runlist=runlist)
    return raw
