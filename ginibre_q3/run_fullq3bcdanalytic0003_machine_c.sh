#!/usr/bin/env bash
set -u

run_dir=/root/q3_full_analytic_20260715
log="$run_dir/fullq3bcdanalytic0003_machine_c.log"
cd "$run_dir" || exit 1

{
    echo "__RUN_ID=fullq3bcdanalytic0003_centered_trace_and_monotonicity"
    echo "__HOST=$(hostname)"
    echo "__START_UTC=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "__CPU_MODEL=$(lscpu | awk -F: '/Model name/{sub(/^[[:space:]]+/, "", $2); print $2; exit}')"
    echo "__LOGICAL_CPUS=$(nproc)"
    sha256sum trace_square_cutoff_mpfr.cpp \
        verify_full_q3_bcd_high_mpfr.cpp \
        verify_full_q3_bcd_mid_mpfr.cpp \
        full_q3_bcd_low_tail_data.hpp
    g++ --version | head -n 1

    g++ -O2 -std=c++17 -Wall -Wextra -Werror \
        trace_square_cutoff_mpfr.cpp -l:libmpfr.so.6 -lgmp \
        -o trace_square_cutoff_mpfr
    compile_trace=$?
    echo "__COMPILE_TRACE_STATUS=$compile_trace"

    g++ -O2 -std=c++17 -Wall -Wextra -Werror \
        verify_full_q3_bcd_high_mpfr.cpp -l:libmpfr.so.6 -lgmp \
        -o verify_full_q3_bcd_high_mpfr
    compile_high=$?
    echo "__COMPILE_HIGH_STATUS=$compile_high"

    g++ -O2 -std=c++17 -Wall -Wextra -Werror \
        verify_full_q3_bcd_mid_mpfr.cpp -l:libmpfr.so.6 -lgmp \
        -o verify_full_q3_bcd_mid_mpfr
    compile_mid=$?
    echo "__COMPILE_MID_STATUS=$compile_mid"

    trace_status=125
    high_status=125
    mid_status=125
    if [[ $compile_trace -eq 0 ]]; then
        ./trace_square_cutoff_mpfr
        trace_status=$?
    fi
    echo "__TRACE_EXIT_STATUS=$trace_status"
    if [[ $compile_high -eq 0 ]]; then
        ./verify_full_q3_bcd_high_mpfr
        high_status=$?
    fi
    echo "__HIGH_EXIT_STATUS=$high_status"
    if [[ $compile_mid -eq 0 ]]; then
        ./verify_full_q3_bcd_mid_mpfr
        mid_status=$?
    fi
    echo "__MID_EXIT_STATUS=$mid_status"

    status=0
    for component_status in \
        "$compile_trace" "$compile_high" "$compile_mid" \
        "$trace_status" "$high_status" "$mid_status"; do
        if [[ $component_status -ne 0 ]]; then
            status=1
        fi
    done
    echo "__EXIT_STATUS=$status"
    echo "__END_UTC=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    exit "$status"
} >"$log" 2>&1
