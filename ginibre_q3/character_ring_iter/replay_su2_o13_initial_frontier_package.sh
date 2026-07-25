#!/usr/bin/env bash
set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
cert=$(cd "$here/../certificates" && pwd)
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

cat "$cert"/su2_o13_initial_frontier_package.b64.part* \
  | base64 -d > "$tmp/package.tar.gz"

echo "a407bf8d0cfb683e54422a84c36533401f5359b751908c3bb36053c78f6ad8b2  $tmp/package.tar.gz" \
  | sha256sum -c -

tar -xzf "$tmp/package.tar.gz" -C "$tmp"
cd "$tmp/o13_initial_frontier_package"
./replay.sh
