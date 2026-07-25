#!/usr/bin/env bash
set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
cert=$(cd "$here/../certificates" && pwd)
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

# Accept either the single-file or the chunked layout.  Chunking is not
# required here -- the payload is ~17 kB against GitHub's 50 MB warning -- and
# every past corruption of these packages was a truncated chunk write.
shopt -s nullglob
parts=("$cert"/su2_o13_initial_frontier_package.b64.part*)
shopt -u nullglob
if [ ${#parts[@]} -eq 0 ]; then
  parts=("$cert/su2_o13_initial_frontier_package.b64")
fi
for part in "${parts[@]}"; do
  if [ ! -f "$part" ]; then
    echo "missing package chunk: $part" >&2
    exit 1
  fi
done

cat "${parts[@]}" | base64 -d > "$tmp/package.tar.gz"

echo "a407bf8d0cfb683e54422a84c36533401f5359b751908c3bb36053c78f6ad8b2  $tmp/package.tar.gz" \
  | sha256sum -c -

tar -xzf "$tmp/package.tar.gz" -C "$tmp"
cd "$tmp/o13_initial_frontier_package"
./replay.sh
