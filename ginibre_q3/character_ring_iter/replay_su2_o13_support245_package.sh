#!/usr/bin/env bash
set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
cert=$(cd "$here/../certificates" && pwd)
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

# The package is small enough to ship unsplit.  Accept either layout so that
# re-chunking for GitHub size limits does not break this wrapper.
shopt -s nullglob
parts=("$cert"/su2_o13_support245_package.b64.part*)
shopt -u nullglob
if [ ${#parts[@]} -eq 0 ]; then
  parts=("$cert/su2_o13_support245_package.b64")
fi
for part in "${parts[@]}"; do
  if [ ! -f "$part" ]; then
    echo "missing package chunk: $part" >&2
    exit 1
  fi
done

cat "${parts[@]}" | base64 -d > "$tmp/package.tar.gz"

echo "da31783b2297c692d7dbc4b10e9f41bc7a03b07b664b43c137403ad06995c76f  $tmp/package.tar.gz" \
  | sha256sum -c -

tar -xzf "$tmp/package.tar.gz" -C "$tmp"
cd "$tmp/o13_support245_package"
./replay.sh
