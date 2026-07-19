#!/usr/bin/env bash
set -Eeuo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
output_dir=${FULL_Q3_INDEPENDENT_LOG_DIR:-}
threads=${FULL_Q3_MODULAR_THREADS:-$(nproc)}

if [[ -z "$output_dir" || "$output_dir" != /* ]]; then
    echo "FULL_Q3_BCD_INDEPENDENT: FAIL: set FULL_Q3_INDEPENDENT_LOG_DIR to an absolute new directory" >&2
    exit 2
fi
if [[ -e "$output_dir" ]]; then
    echo "FULL_Q3_BCD_INDEPENDENT: FAIL: refusing to reuse $output_dir" >&2
    exit 2
fi
mkdir -p -- "$output_dir"

cd -- "$script_dir"
g++ -O3 -DNDEBUG -std=c++20 -fopenmp -Wall -Wextra -Wpedantic \
    -Wconversion -Wsign-conversion -Werror \
    character_ring_iter/verify_full_q3_bcd_modular_moment_checker.cpp \
    -lgmpxx -lgmp -o "$output_dir/modular_moment_checker"

printf 'FULL_Q3_BCD_INDEPENDENT stage=compressed-modular event=start utc=%s threads=%s\n' \
    "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$threads" \
    | tee "$output_dir/driver.log"
OMP_NUM_THREADS="$threads" stdbuf -oL -eL \
    "$output_dir/modular_moment_checker" \
    --moments-log certificates/full_q3/fullq3bcd0010_machine_c.log \
    --reverse-log certificates/full_q3/fullq3bcd0006_machine_c_prefix40.log \
    --progress | tee "$output_dir/modular_moment_checker.log"
printf 'FULL_Q3_BCD_INDEPENDENT stage=compressed-modular event=pass utc=%s\n' \
    "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$output_dir/driver.log"
echo "FULL_Q3_BCD_INDEPENDENT: ALL PASS" | tee -a "$output_dir/driver.log"
sha256sum character_ring_iter/verify_full_q3_bcd_modular_moment_checker.cpp \
    character_ring_iter/full_q3_bcd_remaining_data.hpp \
    certificates/full_q3/fullq3bcd0010_machine_c.log \
    certificates/full_q3/fullq3bcd0006_machine_c_prefix40.log \
    "$output_dir/modular_moment_checker.log" \
    "$output_dir/driver.log" > "$output_dir/logs.sha256"
