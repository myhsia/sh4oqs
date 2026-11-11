#!/usr/bin/env bash
# The repository's session check, declared in `AGENTS.md` as
# `` Session check: `source/tools/session_check.sh` `` and run by the `upkeep`
# plugin at session start with this repository as the working directory.
#
# Scope, deliberately lighter than the synced check in `mkct4oqs`: this
# repository's archives are produced by CI and already validated there (the
# trajectory suite and `sh_error.py` read them), and its data carry no authored
# symlink graph apart from `link_heom.sh`.  So this check covers only what CI
# cannot see -- dangling `output/data` symlinks and a config that no longer
# parses -- and needs no numpy.
#
# Exit codes, per the plugin contract: 0 clean (silent), 1 something to report
# (one line per problem; the plugin keeps the first three), 2 cannot run here.
set -uo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"     # .../<repo>
data="$root/source/output/data"
bad=0

# 1. Dangling symlinks under data/ (a renamed checkout leaves them in bulk).
if [ -d "$data" ]; then
  while IFS= read -r path; do
    [ -n "$path" ] || continue
    printf 'dangling symlink %s -> %s\n' "${path#"$root"/}" "$(readlink "$path")"
    bad=1
  done < <(find "$data" -type l ! -exec test -e {} \; -print 2>/dev/null)
fi

# 2. Every config parses (truncated or mis-indented YAML silently defaults).
#    A missing interpreter is not a failure: the link half above already ran.
py=""
for cand in "${SH_PYTHON:-}" "${CONDA_PREFIX:-}/bin/python" \
            "$(conda info --base 2>/dev/null || true)/envs/sh/bin/python" \
            "$(command -v python3 2>/dev/null || true)" \
            "$(command -v python 2>/dev/null || true)"; do
  if [ -n "$cand" ] && [ -x "$cand" ] && "$cand" -c 'import yaml' 2>/dev/null; then
    py="$cand"
    break
  fi
done

if [ -n "$py" ]; then
  out="$("$py" - "$root" <<'PY' 2>/dev/null
import glob, os, sys
import yaml

root = sys.argv[1]
for path in sorted(glob.glob(os.path.join(root, "source", "config", "*.yaml"))):
    try:
        with open(path) as handle:
            yaml.safe_load(handle)
    except Exception as error:
        print("config %s: %s" % (os.path.relpath(path, root), error))
PY
)"
  if [ -n "$out" ]; then
    printf '%s\n' "$out"
    bad=1
  fi
fi

exit "$bad"
