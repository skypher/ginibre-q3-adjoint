#!/usr/bin/env bash
set -Eeuo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
log_path=${FULL_Q3_MODULAR_LOG:-}
threads=${FULL_Q3_MODULAR_THREADS:-$(nproc)}

if [[ -z "$log_path" || "$log_path" != /* ]]; then
    echo "FULL_Q3_MODULAR_CURRENT: FAIL: set FULL_Q3_MODULAR_LOG to an absolute new path" >&2
    exit 2
fi
if [[ -e "$log_path" ]]; then
    echo "FULL_Q3_MODULAR_CURRENT: FAIL: refusing to overwrite $log_path" >&2
    exit 2
fi

build_dir=$(mktemp -d)
trap 'rm -rf -- "$build_dir"' EXIT
mkdir -p -- "$(dirname -- "$log_path")"

source_file="$script_dir/character_ring_iter/verify_full_q3_bcd_modular_moment_checker.cpp"
data_file="$script_dir/character_ring_iter/full_q3_bcd_remaining_data.hpp"
moments_log="$script_dir/certificates/full_q3/fullq3bcd0010_machine_c.log"
reverse_log="$script_dir/certificates/full_q3/fullq3bcd0006_machine_c_prefix40.log"
binary="$build_dir/verify_full_q3_bcd_modular_moment_checker"

exec > >(tee -a -- "$log_path") 2>&1

echo "__RUN_ID=fullq3bcdmodularfinal0001_current_source"
echo "__START_UTC=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "__HOST=$(uname -a)"
echo "__OMP_NUM_THREADS=$threads"
echo "__GIT_COMMIT=$(git -C "$script_dir" rev-parse HEAD)"
sha256sum "$source_file" "$data_file" "$moments_log" "$reverse_log" "$script_dir/Makefile"
g++ --version | head -n 1

g++ -O3 -DNDEBUG -std=c++20 -fopenmp -Wall -Wextra -Wpedantic \
    -Wconversion -Wsign-conversion -Werror \
    "$source_file" -lgmpxx -lgmp -o "$binary"

set +e
OMP_NUM_THREADS="$threads" stdbuf -oL -eL "$binary" \
    --moments-log "$moments_log" \
    --reverse-log "$reverse_log" \
    --progress
status=$?
set -e

echo "__EXIT_STATUS=$status"
echo "__END_UTC=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
if [[ $status -eq 0 ]]; then
    echo "FULL_Q3_MODULAR_CURRENT: ALL PASS"
else
    echo "FULL_Q3_MODULAR_CURRENT: FAIL"
fi
exit "$status"
