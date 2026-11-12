#!/usr/bin/env bash
# Push the published linux-x64 bundle and the configs to the HPC cluster.
#
#   bash tools/push-binaries.sh <moments|dmrg|sh> [git-ref]
#
# Fetches the linux-x64 asset if ./dist does not hold it, uploads it with
# source/config (the bundles carry no configuration), extracts it in
# ~/oqs-bin and runs --version + --selftest there.  Extraction is retried up to
# three times: the login node has killed tar mid-extraction, and the result of
# that is a half-written tree, not a bad archive (the checksum matches on both
# sides).  Manual only — it spends one OATH code and is never run by CI.
set -euo pipefail

driver="${1:?usage: push-binaries.sh <driver> [git-ref]}"
ref="${2:-HEAD}"
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root"

asset="$driver-linux-x64"
if [ ! -f "dist/$asset.tar.gz" ]; then
  bash tools/fetch-binaries.sh "$driver" -p linux-x64 "$ref"
fi
sha="$(git rev-parse "$ref")"
echo "== $asset @ ${sha:0:7}"

# Prepare the remote tree and upload the tarball plus the configs.  The config
# directory is removed first so a repeated run does not nest config/config.
oathpass 'bash -lc "mkdir -p ~/oqs-bin && rm -rf ~/oqs-bin/config"' > /dev/null
oathpass scp "dist/$asset.tar.gz" :oqs-bin/ > /dev/null
oathpass scp -r source/config :oqs-bin/config > /dev/null
oathpass scp tools/verify_binary.py :oqs-bin/verify_binary.py > /dev/null
echo "   uploaded (bundle, configs, verifier)"

remote="cd ~/oqs-bin && rm -rf $driver && \
for i in 1 2 3; do tar -xzf $asset.tar.gz && break; echo extract-retry; sleep 5; done; \
test -x ./$driver/$driver && ./$driver/$driver --version | head -1 && \
python3 ~/oqs-bin/verify_binary.py $driver ~/oqs-bin/$driver/$driver"
oathpass "bash -lc \"$remote\""
oathpass close

echo "== $driver on the cluster: ~/oqs-bin/$driver/$driver (configs in ~/oqs-bin/config)"
