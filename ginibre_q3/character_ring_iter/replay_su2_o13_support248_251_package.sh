#!/usr/bin/env bash
set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
cert=$(cd "$here/../certificates" && pwd)
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

cat "$cert"/su2_o13_support248_251_package.b64.part* \
  | base64 -d > "$tmp/package.tar.gz"

echo "0422d1fd71b20e3229730abfcd07f4b67d17497a1566e6f9c5b053733f37fdaf  $tmp/package.tar.gz" \
  | sha256sum -c -

tar -xzf "$tmp/package.tar.gz" -C "$tmp"
cd "$tmp/o13_support248_251_package"
./replay.sh
