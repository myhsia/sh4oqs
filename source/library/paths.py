"""Where the resources, caches, CPUs and build identity come from.

The drivers are the same code in two worlds: a checkout (configs in
``./config``, outputs under ``./output``) and a PyInstaller bundle (configs
inside the bundle, no writable install directory).  Only the differences live
here, so a source run keeps its ``./config``/``./output`` convention exactly as
before, and a frozen binary finds the configs it shipped with.  A local
``./config`` always wins: that is how a binary takes an unsaved override on the
cluster without being rebuilt.

The rest is what a batch system needs and a laptop does not: a CPU count that
respects a Slurm allocation (``os.cpu_count()`` reports the whole node, so
``-njobs -1`` on a 64-core node would oversubscribe two 32-core jobs), and the
build identity printed by ``-version`` (git sha, build date and the dependency
versions actually imported), baked in by ``tools/build_binary.sh``.
"""
import json
import os
import platform
import subprocess
import sys
import tempfile

#: Modules whose frozen version ``-version`` reports, when the driver imported them.
#: Only loaded modules are listed: importing the rest just to print a number would
#: cost seconds (tenpy, matplotlib) and would change what the flag tests.
_VERSIONED = ("numpy", "scipy", "numba", "sympy", "mpmath", "tenpy", "joblib",
              "matplotlib")


def is_frozen() -> bool:
    """True inside a PyInstaller bundle."""
    return bool(getattr(sys, "frozen", False))


def bundle_dir() -> str:
    """Directory holding bundled data.

    PyInstaller sets ``sys._MEIPASS`` to the extraction directory (onefile) or
    to ``_internal`` (onedir).  Outside a bundle the executable's directory is
    the best available answer; no caller depends on it there, because a source
    run resolves ``./config`` first.
    """
    return getattr(sys, "_MEIPASS",
                   os.path.dirname(os.path.abspath(sys.executable)))


def config_dir(explicit: str = None) -> str:
    """Directory that holds the YAML configs.

    ``explicit`` wins.  Otherwise ``./config``, then — when frozen — the
    directory beside the bundle (``<exe dir>/../config``), which is the layout
    the release tarball extracts to (``moments/moments`` beside ``config/``).

    **No configuration is frozen into a bundle.**  A packed binary takes its
    YAMLs from the filesystem, so a missing config is a loud error in the
    driver rather than a silent fallback to a copy nobody can see or edit.
    """
    if explicit:
        return explicit
    if os.path.isdir("config"):
        return "config"
    if is_frozen():
        beside = os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(sys.executable))),
            "config")
        if os.path.isdir(beside):
            return beside
    return "config"


def install_cache_env() -> None:
    """Keep numba's cache out of the install directory.

    Only meaningful when frozen: there ``./output`` may be unwritable, and a
    cache written next to the executable would be shared by every run.
    ``OQS_CACHE_DIR`` overrides the location (a cluster run wants it on a fast
    local filesystem rather than on the NFS home).  Matplotlib is not part of
    the compute binary, so it has no cache to place.
    """
    if not is_frozen():
        return
    cache = os.environ.get("OQS_CACHE_DIR",
                           os.path.join(tempfile.gettempdir(), "oqs-cache"))
    os.makedirs(os.path.join(cache, "numba"), exist_ok = True)
    os.environ.setdefault("NUMBA_CACHE_DIR", os.path.join(cache, "numba"))


def n_cpus() -> int:
    """CPUs this process may use; a Slurm allocation wins over the node count.

    Checked in the order Slurm itself fills them in: ``SLURM_CPUS_PER_TASK``
    (one task's share), ``SLURM_JOB_CPUS_PER_NODE`` (the allocation, possibly
    ``"64(x2)"`` over several nodes), then ``SLURM_CPUS_ON_NODE``.
    """
    for var in ("SLURM_CPUS_PER_TASK", "SLURM_JOB_CPUS_PER_NODE",
                "SLURM_CPUS_ON_NODE"):
        raw = os.environ.get(var, "").strip().split("(")[0]
        if not raw:
            continue
        try:
            value = int(raw)
        except ValueError:
            continue
        if value > 0:
            return value
    return os.cpu_count() or 1


def _git_sha() -> str:
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    try:
        out = subprocess.run(["git", "-C", root, "rev-parse", "--short=12", "HEAD"],
                             capture_output = True, text = True, timeout = 5)
    except (OSError, subprocess.SubprocessError):
        return "unknown"
    return out.stdout.strip() or "unknown"


def build_info() -> dict:
    """Identity of this build: sha, date, interpreter, platform, dependencies.

    A frozen binary carries ``_buildinfo.json`` (written by
    ``tools/build_binary.sh``); a source run falls back to the checkout's git
    sha, and to ``unknown`` when git is absent or the tree is not a repository.
    """
    info = {"frozen": is_frozen(), "python": platform.python_version(),
            "platform": f"{sys.platform}-{platform.machine()}",
            "sha": "unknown", "built": "unknown"}
    path = os.path.join(bundle_dir(), "_buildinfo.json")
    if os.path.isfile(path):
        try:
            with open(path) as f:
                info.update(json.load(f))
        except (OSError, ValueError):
            pass
    if info.get("sha") in (None, "", "unknown") and not is_frozen():
        info["sha"] = _git_sha()
    deps = {}
    for mod in _VERSIONED:
        module = sys.modules.get(mod)
        if module is None:
            continue
        deps[mod] = str(getattr(module, "__version__", "?"))
    info["deps"] = deps
    return info


def version_string(name: str) -> str:
    """One parseable line: name, sha, build date, interpreter, deps."""
    info = build_info()
    deps = " ".join(f"{k}={v}" for k, v in info["deps"].items())
    return (f"{name} sha={info['sha']} built={info['built']} "
            f"python={info['python']} platform={info['platform']} "
            f"frozen={str(bool(info['frozen'])).lower()} {deps}")
