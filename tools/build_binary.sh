#!/usr/bin/env bash
# Build the frozen binary for one driver, verify it, and pack the release assets.
#
#   bash tools/build_binary.sh <moments|dmrg|sh>
#
# The interpreter is ${PYTHON:-python}; the CI passes the manylinux container's
# own /opt/python/cp314-cp314/bin/python so the Linux binary links only glibc
# 2.28 (the cluster is Rocky 9.6 / glibc 2.34, and ubuntu-22.04's 2.35 is
# already too new).  Everything lands in ./dist:
#
#   <driver>/<driver>                 the onedir bundle (executable + _internal)
#   <driver>-<os>-<arch>.tar.gz       the shipped asset
#   <driver>-<os>-<arch>.tar.gz.sha256
#   <driver>-<os>-<arch>.json         build identity: sha, date, python, pip freeze
#
# The selftest runs before packing: --version and --selftest on the frozen
# executable, in a temporary directory, so the bundled configs and the process
# pool are exercised exactly as a user would exercise them.
set -euo pipefail

driver="${1:?usage: build_binary.sh <moments|dmrg|sh>}"
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root"
python="${PYTHON:-python}"

case "$(uname -s)" in
  Darwin) os=macos ;;
  Linux)  os=linux ;;
  *) echo "unsupported OS: $(uname -s)" >&2; exit 1 ;;
esac
case "$(uname -m)" in
  arm64|aarch64) arch=arm64 ;;
  x86_64|amd64)  arch=x64 ;;
  *) echo "unsupported architecture: $(uname -m)" >&2; exit 1 ;;
esac
platform="$os-$arch"

build="${TMPDIR:-/tmp}/oqs-build-$driver"
rm -rf "$build" "dist/$driver"
mkdir -p "$build" dist

# The manylinux image has no git, so the workflow passes the commit in.
sha="${OQS_SHA:-$(git -C "$root" rev-parse HEAD 2>/dev/null || echo unknown)}"

# The bundle is the compute path only: no plotting, no notebooks, no qutip (the
# HEOM reference is a separate driver).  matplotlib reaches the graph through a
# lazy import under sympy even when nothing plots, and its data plus compiled
# deps were 23 MB per bundle before this exclude.
excluded="qutip matplotlib PIL lxml tkinter IPython pytest"
excludes=()
for _mod in $excluded; do excludes+=(--exclude-module "$_mod"); done

"$python" - "$build/_buildinfo.json" "$driver" "$sha" "$platform" "$excluded" <<'PY'
import json, platform, subprocess, sys
path, name, sha, plat, excluded = sys.argv[1:6]
freeze = subprocess.run([sys.executable, "-m", "pip", "freeze"],
                        capture_output = True, text = True, check = False)
date = subprocess.run(["date", "-u", "+%Y-%m-%dT%H:%M:%SZ"],
                      capture_output = True, text = True, check = False)
with open(path, "w") as f:
    json.dump({"name": name, "sha": sha,
               "built": date.stdout.strip() or "unknown",
               "python": platform.python_version(), "platform": plat,
               "excluded": excluded.split(),
               "packages": sorted(freeze.stdout.split())},
              f, indent = 2, sort_keys = True)
PY

# The moment engines hand `mpmath` to their joblib workers, and mpmath 1.4 cannot
# be pickled ("cannot pickle '_MPMathModule' object"); pin the version the
# parallel path was measured with, independent of the unpinned requirements.txt.
if [ "$driver" = moments ]; then
  "$python" -m pip install --quiet "mpmath==1.3.0"
fi

# Freeze a copy under a distinct module name: dmrg.py's own `from dmrg import
# chain` refers to the library package of that name, and an entry script called
# `dmrg` shadowed it ("partially initialized module 'dmrg'").  The copy keeps
# the driver byte-identical and the package unambiguous.
entry="$build/${driver}_main.py"
cp "$root/source/$driver.py" "$entry"

"$python" -m PyInstaller --noconfirm --clean --onedir --name "$driver" \
  --paths "$root/source/library" \
  --add-data "$build/_buildinfo.json:." \
  "${excludes[@]}" \
  --distpath "$root/dist" --workpath "$build/work" --specpath "$build" \
  "$entry"

echo "== version"
"dist/$driver/$driver" --version | tee "$build/version.txt"
if [ -n "${OQS_SHA:-}" ] && ! grep -q "sha=$OQS_SHA" "$build/version.txt"; then
  echo "build identity mismatch: expected sha=$OQS_SHA" >&2
  exit 1
fi
echo "== verify"
"$python" tools/verify_binary.py "$driver" "dist/$driver/$driver" \
  | tee "$build/verify.txt"

tar -czf "dist/$driver-$platform.tar.gz" -C dist "$driver"
cp "$build/_buildinfo.json" "dist/$driver-$platform.json"
if command -v shasum > /dev/null 2>&1; then
  ( cd dist && shasum -a 256 "$driver-$platform.tar.gz" \
      > "$driver-$platform.tar.gz.sha256" )
else
  ( cd dist && sha256sum "$driver-$platform.tar.gz" \
      > "$driver-$platform.tar.gz.sha256" )
fi
ls -lh "dist/$driver-$platform.tar.gz"
