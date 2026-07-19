#!/usr/bin/env bash
set -u -o pipefail

cd "$(dirname "$0")"

threads="${OMP_NUM_THREADS:-2}"
echo "__RUN_ID=full_q3_extension_acyclic_source_optimus"
echo "__HOST=$(hostname)"
echo "__START_UTC=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "__CPU_MODEL=$(sed -n 's/^model name[[:space:]]*: //p' /proc/cpuinfo | head -1)"
echo "__REPLAY_THREADS=$threads"
sha256sum \
  Makefile \
  full_q3_extension.tex \
  verify_full_q3_document.py \
  run_fullq3complete0002_optimus.sh \
  certificates/full_q3/full_q3_source_manifest.sha256

OMP_NUM_THREADS="$threads" nice -n 19 make full-q3-extension
status=$?
echo "__EXIT_STATUS=$status"
echo "__END_UTC=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
exit "$status"
