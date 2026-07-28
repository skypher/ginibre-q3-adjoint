#!/usr/bin/env bash
set -u

run_dir=${1:?usage: run_su2_corrected_308_machine_c_when_free.sh RUN_DIR}
old_pattern='[v]erify_su2_seven_shallow_z3 --bounded'

while pgrep -f "$old_pattern" >/dev/null; do
    printf 'WAIT_OLD_BOUNDED %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    sleep 60
done

expected_shallow=5c4f65cd45ca6d5237f43d0641b97cf288edc6a2f251405f605e9b2325ed86aa
expected_residual=d71e4634955def0fed2a1a8078cdae3544cb3d5b18c8cb68b3eabbc99889660c
expected_binary=efb0b72fb59940e56466fd33eef7a020700f05c138dc08fe49197de7ab6e1e1c

printf '%s  %s\n' "$expected_shallow" \
    "$run_dir/verify_su2_seven_shallow_z3.cpp" \
    | sha256sum -c -
printf '%s  %s\n' "$expected_residual" \
    "$run_dir/verify_su2_seven_residual_z3.cpp" \
    | sha256sum -c -
printf '%s  %s\n' "$expected_binary" \
    "$run_dir/verify_su2_seven_shallow_z3" \
    | sha256sum -c -

cores=$(nproc)
runnable=$(ps -eLo stat= | awk '$1 ~ /^R/ {count++} END {print count + 0}')
cpu_workers=$((cores - runnable))
available_kib=$(awk '/^MemAvailable:/ {print $2}' /proc/meminfo)
reserve_kib=$((8 * 1024 * 1024))
per_worker_kib=$((1536 * 1024))
if ((available_kib <= reserve_kib)); then
    memory_workers=0
else
    memory_workers=$(((available_kib - reserve_kib) / per_worker_kib))
fi

workers=$cores
if ((cpu_workers < workers)); then
    workers=$cpu_workers
fi
if ((memory_workers < workers)); then
    workers=$memory_workers
fi
while ((workers < 1)); do
    printf 'WAIT_CAPACITY %s runnable=%d available_kib=%d\n' \
        "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$runnable" "$available_kib"
    sleep 60
    runnable=$(ps -eLo stat= | awk '$1 ~ /^R/ {count++} END {print count + 0}')
    cpu_workers=$((cores - runnable))
    available_kib=$(awk '/^MemAvailable:/ {print $2}' /proc/meminfo)
    if ((available_kib <= reserve_kib)); then
        memory_workers=0
    else
        memory_workers=$(((available_kib - reserve_kib) / per_worker_kib))
    fi
    workers=$cores
    if ((cpu_workers < workers)); then
        workers=$cpu_workers
    fi
    if ((memory_workers < workers)); then
        workers=$memory_workers
    fi
done

log_path="$run_dir/su2_corrected_308_machine_c.log"
printf 'LAUNCH %s workers=%d runnable=%d available_kib=%d\n' \
    "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$workers" "$runnable" "$available_kib"

LD_LIBRARY_PATH="$run_dir" Q3_MAX_THREADS="$workers" \
    nice -n 10 "$run_dir/verify_su2_seven_shallow_z3" --small \
    >"$log_path" 2>&1
status=$?

printf 'EXIT %s status=%d\n' \
    "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$status"
sha256sum "$log_path"
exit "$status"
