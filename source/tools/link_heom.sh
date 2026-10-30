#!/usr/bin/env bash
# Link the reference correlation matrix produced in ../dmrg4oqs into ./output/data.
# The reference is one archive per config, corrmat_heom.npz (save_group schema,
# keyed by qhm), generated once (``python reference.py``) and shared by the three
# checkouts: ``dmrg4oqs`` holds the real files, ``mkct4oqs``/``sh4oqs`` link them
# here so that every figure script finds ``output/data/<prefix>/corrmat_heom.npz``
# with its default ``-ref ./output/data``.  The links are local (output/ is
# gitignored).  Missing sources are skipped, so existing data stays usable.
#
# ``HEOM_DATA=<dir>`` overrides the source data root.  When the source resolves
# to this checkout's own data root (i.e. running in ``dmrg4oqs``) the script is a
# no-op: relinking a file onto itself would destroy the real archive.
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"        # .../<repo>/source
dst="$here/output/data"
mkdir -p "$dst"

src="${HEOM_DATA:-}"
if [ -z "$src" ]; then
  # The sibling reference checkout: ../dmrg4oqs in a plain layout,
  # ../../dmrg4oqs inside the oqs workspace.
  for cand in "$here/../dmrg4oqs/source/output/data" \
              "$here/../../dmrg4oqs/source/output/data"; do
    if [ -d "$cand" ]; then src="$cand"; break; fi
  done
fi
[ -n "$src" ] && [ -d "$src" ] || { echo "no reference data (set HEOM_DATA)"; exit 0; }

src="$(cd "$src" && pwd)"
dst="$(cd "$dst" && pwd)"
if [ "$src" = "$dst" ]; then
  echo "reference data are canonical here ($src); nothing to link"
  exit 0
fi

n=0
for dir in "$src"/*/; do
  p="$(basename "$dir")"
  mkdir -p "$dst/$p"
  # prune links whose target is gone (e.g. a retired prefix)
  for f in "$dst/$p"/*.npz; do
    [ -L "$f" ] && [ ! -e "$f" ] && rm -f "$f"
  done
  # drop the obsolete per-quantity reference links of the pre-corrmat layout
  for f in "$dst/$p"/{pop,cohe,cohem,corr,xorr}_heom.npz; do
    [ -L "$f" ] && rm -f "$f"
  done
  f="$dir"corrmat_heom.npz
  [ -e "$f" ] || continue
  # never replace a real archive of this checkout with a link
  [ -e "$dst/$p/corrmat_heom.npz" ] && [ ! -L "$dst/$p/corrmat_heom.npz" ] && continue
  ln -sfn "$(cd "$(dirname "$f")" && pwd)/$(basename "$f")" "$dst/$p/corrmat_heom.npz"
  n=$((n + 1))
done
echo "linked $n corrmat_heom.npz files from $src"
