"""Filesystem and IO helpers: config loading, output writers and loaders."""
import copy, glob, os, yaml
import numpy as np
from pathlib import Path
import paths


def _apply_relative(value, base):
    if not isinstance(value, str) or base is None:
        return value
    s = value.strip()
    if len(s) < 2 or s[0] not in "*/+-":
        return value
    try:
        factor = float(s[1:])
    except ValueError:
        return value
    if s[0] == "*":
        return base * factor
    if s[0] == "/":
        return base / factor
    if s[0] == "+":
        return base + factor
    return base - factor


def _merge(base, override):
    merged = copy.deepcopy(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = _merge(merged[key], value)
        else:
            merged[key] = _apply_relative(value, merged.get(key))
    return merged


def parse_xlims(raw_xlims, reorder: bool = True):
    """Normalize YAML ``xlims`` into ``[[x_lo, x_hi], x_step, ...]``.

    When ``reorder`` is true and two windows are given, the two groups are
    swapped so that ``xlims[0]`` is the larger (later) window, i.e.
    ``xlims[0][1]`` is the maximum time (used for the calculation range).
    """
    xlims = [tuple(x) if isinstance(x, list) else x for x in raw_xlims]
    if reorder and len(xlims) == 4:
        xlims = [xlims[2], xlims[3], xlims[0], xlims[1]]
    return xlims


def _config_path(source, config_dir: str = None) -> str:
    """Resolve one spec: a path (or ``.yaml`` name) stays, a bare name gains one."""
    text = str(source)
    if (text.endswith((".yaml", ".yml")) or os.sep in text
            or (os.altsep is not None and os.altsep in text)):
        return text
    return os.path.join(paths.config_dir(config_dir), text + ".yaml")


def load_config(source: str, config_dir: str = None):
    """Load one YAML config and resolve its ``base:`` chain.

    ``source`` is a bare name (looked up in ``config_dir``, which defaults to
    :func:`paths.config_dir`) or a path to a ``.yaml`` file.  A ``base:`` entry
    resolves **next to the file that names it**, so a config outside the
    default directory keeps its chain intact.

    Nothing is ever loaded from a bundle: a frozen driver reads the YAMLs the
    user points it at, and a missing file raises ``FileNotFoundError`` naming
    the path it tried.
    """
    path = _config_path(source, config_dir)
    with open(path, "r") as f:
        override = yaml.safe_load(f)
    base = override.pop("base", None)
    if base is None:
        return override
    directory = os.path.dirname(path) or "."
    base_name = os.path.splitext(str(base))[0] + ".yaml"
    base_path = (base_name if os.path.isabs(base_name)
                 else os.path.join(directory, base_name))
    with open(base_path, "r") as f:
        params = yaml.safe_load(f)
    return _merge(params, override)


def resolve_configs(specs, config_dir: str = None):
    """Expand config names, paths and glob patterns into ``[(name, path), ...]``.

    The driver expands the pattern itself, so ``'config/*.yaml'`` works in any
    shell and on a batch node that would not expand it.  A run is named by its
    file stem, which is what names ``output/data/<name>/``; the given order is
    kept and each glob's matches are sorted.  A pattern that matches nothing,
    or a name that resolves to no file, raises ``FileNotFoundError`` here
    rather than failing later inside a worker.
    """
    resolved, seen = [], set()
    for spec in specs:
        text = str(spec)
        if glob.has_magic(text):
            matched = sorted(glob.glob(text))
            if not matched:
                raise FileNotFoundError(f"no config matches {text!r}")
        else:
            matched = [_config_path(text, config_dir)]
        for path in matched:
            if not os.path.isfile(path):
                raise FileNotFoundError(f"no config file at {path!r}")
            key = os.path.abspath(path)
            if key in seen:
                continue
            seen.add(key)
            resolved.append((os.path.splitext(os.path.basename(path))[0], path))
    if not resolved:
        raise FileNotFoundError("no config given")
    return resolved


def ensure_dir(directory) -> None:
    """Create ``directory`` (and any missing parents) if it does not exist."""
    Path(directory).mkdir(parents = True, exist_ok = True)


def save_time(t, Y, directory, filename, identifier = "",
              compress: bool = False) -> None:
    filepath = Path(directory) / filename
    ensure_dir(filepath.parent)
    if compress:
        np.savez(filepath, t = np.asarray(t), y = np.asarray(Y),
                 identifier = np.array(identifier))
        return
    if np.iscomplexobj(Y):
        datout = np.column_stack((t, Y.real, Y.imag))
        header_list = (["t", "Re", "Im"] if identifier == "" else
                       ["t", f"Re[{identifier}]", f"Im[{identifier}]"])
    else:
        Y = np.asarray(Y)
        Y = Y[:, None] if Y.ndim == 1 else Y
        datout = np.column_stack((t, Y))
        prefix = identifier if identifier else "Y"
        header_list = ["t"] + [f"{prefix}{i}" for i in range(Y.shape[1])]
    fmt = "%15.6e" * datout.shape[1]
    header = f"{header_list[0]:>13}" + "".join(f"{h:>15}" for h in header_list[1:])
    np.savetxt(filepath, datout, fmt = fmt, header = header)


def load_npz(path):
    data = np.load(path)
    return data["t"], data["y"]


# Metadata keys of a ``save_group`` archive (not per-case series).
GROUP_META = ("identifier", "component", "branch")


def save_group(series, directory, filename, identifier = "",
               component = None, branch = None, compress: bool = True) -> None:
    """Save per-case series with their own time axes in a single ``.npz``.

    ``series`` maps a case name (or any key, e.g. the Padé order ``"2/3"``) to
    its ``(t, y)`` pair; each pair is stored as ``t_<key>`` and ``<key>``, so
    the archive reads as ``np.load(path)["qhm"]`` with ``["t_qhm"]``.
    ``identifier`` records the observable, ``component`` labels the inner axis
    of a two-component series (e.g. ``("P0", "P1")`` for populations) and
    ``branch`` the first axis of stacked branches (e.g. ``("+", "-")`` for the
    two coherences).
    """
    filepath = Path(directory) / filename
    ensure_dir(filepath.parent)
    payload = {}
    for key, (t, y) in series.items():
        payload["t_" + str(key)] = np.asarray(t)
        payload[str(key)]        = np.asarray(y)
    if identifier:
        payload["identifier"] = np.array(identifier)
    if component is not None:
        payload["component"] = np.asarray(component)
    if branch is not None:
        payload["branch"] = np.asarray(branch)
    (np.savez_compressed if compress else np.savez)(filepath, **payload)


def load_group(path):
    """Load a ``save_group`` archive as ``(series, meta)``.

    ``series`` maps each case (or order) key to its ``(t, y)`` pair; ``meta``
    holds the ``identifier``/``component``/``branch`` labels that are present.
    """
    with np.load(path) as f:
        series = {k: (f["t_" + k], f[k]) for k in f.files
                  if k not in GROUP_META and "t_" + k in f.files}
        meta   = {k: f[k] for k in GROUP_META if k in f.files}
    return series, meta


def moments_filename(method: str, projection: str = "matrix") -> str:
    """Canonical moments archive ``moments.<method>.{mat,scl}.npz``.

    ``method`` is ``"block_sp"`` or ``"symbolic"`` and ``projection`` is
    ``"matrix"`` (``mat``) or ``"scalar"`` (``scl``); both methods' archives
    live in the same ``output/data/<prefix>/`` directory, one array per case.
    """
    return f"moments.{method}.{'scl' if projection == 'scalar' else 'mat'}.npz"


def save_moments(path, series, timing: dict | None = None) -> None:
    """Save one moments array per case, each under its case name.

    ``series`` maps a case (``qhm``/``qah``/``chm``/``cah``) to its array; a
    scalar ``(n_moments,)`` or matrix ``(n_moments, d, d)`` complex array keeps
    its native shape (no flattening).  ``timing`` optionally maps a case to its
    ``(n_orders + 1, 3)`` array of ``[order, cpu_s, wall_s]`` (order 0 is the
    setup), stored as ``timing_<case>``.
    """
    filepath = Path(path)
    ensure_dir(filepath.parent)
    data = {str(k): np.asarray(v) for k, v in series.items()}
    for case, arr in (timing or {}).items():
        data[f"timing_{case}"] = np.asarray(arr, dtype = float)
    np.savez(filepath, **data)


def load_moments(path) -> dict:
    """Load a moments archive as ``{case: array}`` (the ``timing_*`` keys are
    skipped; use :func:`load_timing` for those).

    Each array keeps its native shape: ``(n_moments, d, d)`` for the matrix
    projection or ``(n_moments,)`` for the scalar projection.
    """
    with np.load(path) as f:
        return {k: f[k] for k in f.files if not k.startswith("timing_")}


def load_timing(path) -> dict:
    """Load the per-case timing arrays ``{case: (n_orders + 1, 3)}`` from a
    moments archive (columns ``order, cpu_s, wall_s``; row 0 is the setup)."""
    with np.load(path) as f:
        return {k[len("timing_"):]: f[k] for k in f.files
                if k.startswith("timing_")}


def tab_str(s, N = 80, char = '=', border = '!'):
    """Pad ``s`` to a border line of constant width ``N`` for console logging."""
    n = len(s)
    n_left = (N - n - 4) // 2
    n_right = N - n - 4 - n_left
    return border + char * n_left + ' ' + s + ' ' + char * n_right + border
