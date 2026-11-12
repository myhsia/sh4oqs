#!/usr/bin/env bash
# Download, verify and install the frozen binary for one driver.
#
#   bash tools/fetch-binaries.sh <moments|dmrg|sh> [-p PLATFORM]... [git-ref]
#
# All published platforms land in ./dist as <driver>-<platform>.tar.gz with its
# .sha256 and manifest; only the host platform is extracted to dist/<driver>/
# and run.  A foreign asset cannot execute here, so it is verified without
# being run: the checksum, the manifest sha (it must be the requested commit),
# the presence of <driver>/<driver> in the archive, and the executable's
# machine type (x86-64, aarch64, Mach-O arm64).  The runtime proof for those is
# the CI run on a native runner (plus the bare rockylinux:9 container), and for
# linux-x64 the cluster run.
#
# The release tag is derived from the commit, bin-<sha7>, so source and binary
# cannot drift; a commit that built nothing falls back to the newest release and
# prints the sha it received.  The previous host bundle is restored if the new
# one fails its selftest.
set -euo pipefail

driver=""
platforms=()
ref="HEAD"
while [ $# -gt 0 ]; do
  case "$1" in
    -p|--platform)
      IFS=, read -r -a _split <<< "${2:?-p needs a platform}"
      platforms+=("${_split[@]}")
      shift 2 ;;
    -a|--all) shift ;;
    -h|--help) sed -n '2,18p' "$0"; exit 0 ;;
    -*) echo "unknown option: $1" >&2; exit 2 ;;
    *) if [ -z "$driver" ]; then driver="$1"; else ref="$1"; fi; shift ;;
  esac
done
[ -n "$driver" ] || { echo "usage: fetch-binaries.sh <driver> [-p PLATFORM]... [ref]" >&2; exit 2; }

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root"

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
host="$os-$arch"

known=(linux-x64 linux-arm64 macos-arm64)
if [ ${#platforms[@]} -eq 0 ]; then platforms=("${known[@]}"); fi
for platform in "${platforms[@]}"; do
  printf '%s\n' "${known[@]}" | grep -qx "$platform" \
    || { echo "unknown platform: $platform (choose from ${known[*]})" >&2; exit 2; }
done

sha="$(git rev-parse "$ref")"
short="${sha:0:7}"
tag="bin-$short"

proxy_file="$HOME/.config/opencode/proxy"
if [ -r "$proxy_file" ]; then
  proxy="$(cat "$proxy_file")"
  export https_proxy="$proxy" http_proxy="$proxy"
fi

mkdir -p dist
# A commit that changed nothing under source/ builds no binary, so HEAD does not
# always have a release; fall back to the newest one and say which commit it is.
if ! gh release view "$tag" > /dev/null 2>&1; then
  echo "no release $tag; falling back to the newest bin-* release" >&2
  tag="$(gh release list --limit 100 --json tagName,publishedAt \
        --jq '[.[] | select(.tagName | startswith("bin-"))]
              | sort_by(.publishedAt) | last | .tagName')"
  [ -n "$tag" ] && [ "$tag" != "null" ] \
    || { echo "no bin-* release in this repository" >&2; exit 1; }
fi
echo "== $driver $tag (host $host)"

for platform in "${platforms[@]}"; do
  asset="$driver-$platform"
  echo "-- $asset"
  gh release download "$tag" --pattern "$asset.*" --dir dist --clobber

  if command -v shasum > /dev/null 2>&1; then
    ( cd dist && shasum -a 256 -c "$asset.tar.gz.sha256" )
  else
    ( cd dist && sha256sum -c "$asset.tar.gz.sha256" )
  fi

  got="$(sed -n 's/.*"sha": "\([^"]*\)".*/\1/p' "dist/$asset.json" | head -1)"
  if [ "$tag" = "bin-$short" ] && [ "$got" != "$sha" ]; then
    echo "manifest sha $got does not match the requested commit $sha" >&2
    exit 1
  fi
  [ -n "$got" ] || { echo "no sha in dist/$asset.json" >&2; exit 1; }

  member="$driver/$driver"
  scratch="$(mktemp -d)"
  # No pipe into grep/file: both exit early and would SIGPIPE tar under pipefail,
  # which is how the host asset (the largest member) failed the first time.
  tar -tzf "dist/$asset.tar.gz" > "$scratch/list.txt"
  grep -qx "$member" "$scratch/list.txt" \
    || { echo "archive has no $member" >&2; rm -rf "$scratch"; exit 1; }
  tar -xzf "dist/$asset.tar.gz" -C "$scratch" "$member"
  desc="$(file -b "$scratch/$member")"
  rm -rf "$scratch"
  case "$platform" in
    linux-x64)   shape="ELF 64-bit";    machine="x86-64" ;;
    linux-arm64) shape="ELF 64-bit";    machine="aarch64" ;;
    macos-arm64) shape="Mach-O 64-bit"; machine="arm64" ;;
  esac
  case "$desc" in
    *"$shape"*"$machine"*) echo "   $desc" ;;
    *) echo "machine type does not match $platform: $desc" >&2; exit 1 ;;
  esac

  if [ "$platform" != "$host" ]; then
    echo "   checksum, manifest and machine type only (not run here)"
    continue
  fi

  backup="dist/.$driver.previous"
  rm -rf "$backup"
  if [ -d "dist/$driver" ]; then mv "dist/$driver" "$backup"; fi
  tar -xzf "dist/$asset.tar.gz" -C dist
  "dist/$driver/$driver" --version
  # The numeric check drives the binary through its normal flags; it needs a
  # python3 (stdlib only), and says so rather than pretending when there is none.
  py="$(command -v python3 || command -v python || true)"
  if [ -n "$py" ]; then
    if ! "$py" tools/verify_binary.py "$driver" "dist/$driver/$driver"; then
      echo "verification FAILED; restoring the previous bundle" >&2
      rm -rf "dist/$driver"
      [ -d "$backup" ] && mv "$backup" "dist/$driver"
      exit 1
    fi
  else
    echo "   no python3: numeric verification skipped"
  fi
  rm -rf "$backup"
done

echo "== dist"
ls -1 "dist" | grep -E "^$driver-" | sed 's/^/   /'
