#!/usr/bin/env bash
set -u

run_root=/root/q3_full_match_clean_20260715
project="$run_root/ginibre_q3"
log="$run_root/clean_room_machine_c.log"
cd "$project" || exit 1

{
    echo "__RUN_ID=full_q3_clean_current_source_machine_c"
    echo "__HOST=$(hostname)"
    echo "__START_UTC=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "__CPU_MODEL=$(lscpu | awk -F: '/Model name/{sub(/^[[:space:]]+/, "", $2); print $2; exit}')"
    echo "__LOGICAL_CPUS=$(nproc)"
    echo "__REPLAY_THREADS=32"
    sha256sum run_full_q3_clean_replay_machine_c.sh \
        clean_room_replay.py replay_sources.sha256 \
        certificates/full_q3/full_q3_source_manifest.sha256
    python3 clean_room_replay.py --threads 32
    status=$?
    echo "__EXIT_STATUS=$status"
    echo "__END_UTC=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    exit "$status"
} >"$log" 2>&1
