#!/usr/bin/env bash
set -u

run_dir=/root/q3_full_bounded_littlewood_20260715
log="$run_dir/fullq3bcd0010_machine_c.log"
cd "$run_dir" || exit 1

{
    echo "__RUN_ID=fullq3bcd0010_fail_closed_final"
    echo "__HOST=$(hostname)"
    echo "__START_UTC=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "__OMP_NUM_THREADS=24"
    sha256sum verify_full_q3_bcd_bounded_littlewood_gmp.cpp \
        full_q3_bcd_remaining_data.hpp
    g++ --version | head -n 1
    time env OMP_NUM_THREADS=24 \
        ./verify_full_q3_bcd_bounded_littlewood_gmp --progress
    status=$?
    echo "__EXIT_STATUS=$status"
    echo "__END_UTC=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    exit "$status"
} >"$log" 2>&1
