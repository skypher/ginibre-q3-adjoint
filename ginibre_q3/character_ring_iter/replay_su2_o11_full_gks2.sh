#!/bin/sh
set -eu
HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
GQ3=$(CDPATH= cd -- "$HERE/.." && pwd)
PKG="$GQ3/certificates/su2_o11_full_gks2_package"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT HUP INT TERM
cat "$PKG"/part*.b64 | tr -d '\n' | base64 -d > "$TMP/certificate.tar.gz"
printf '%s  %s\n' 'c2d50ce9fbb48cbe3caec61aff18f1d768b1c067a630eeb33524e5a6c73d1091' "$TMP/certificate.tar.gz" | sha256sum -c -
tar -xzf "$TMP/certificate.tar.gz" -C "$TMP"
cd "$TMP"
./replay.sh "$HERE/verify_su2_o11_first_chamber_exact.cpp"
