#!/usr/bin/env bash
set -euo pipefail

root=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
collector_path=$(realpath "$0")
collector_start_sha=$(sha256sum "$collector_path" | cut -d' ' -f1)
archive="$root/certificates/full_q3/distributed0001"
a_root=/tmp/fullq3cleanfleet0001_a
b_host=root@38.49.219.164
b_root=/root/fullq3_frontier_b_20260715
c_host=root@149.50.101.50
c_root=/root/fullq3_frontier_c_20260715_w1

ledger_sha=28b51756e053cb160e9cf5dde4f97a27ae416f084f7ef23ffedd02409098464a
wrapper_sha=784c47f797a61e805d509c8e6262ccf5abd7f3bf7757a3b8fb4dad6d84b80adf
verifier_sha=4097d664d30b50bb379aa1fdf1aef1ae930556c19214330283a75f5ae68f4aa4
audit_wrapper_sha=208062defb8d3686cb254950b88b48b3a5380e0e5d5ddecb33bd01d39ca22a24
machine_a_monitor_sha=59290ccecea9d51f401a7c7637587656bd80c69f57711a5f977107ecd0edd867
resource_monitor_sha=49ff73fe1eebbaf0fef5344f0d8884103fa28d57c1dfd00e6723f20a09143040
c_min_available_kib=12582912
c_max_swap_used_kib=2097152

require_file() {
  test -f "$1" || { echo "missing file: $1" >&2; exit 1; }
}

require_sha() {
  local path=$1
  local expected=$2
  require_file "$path"
  local actual
  actual=$(sha256sum "$path" | cut -d' ' -f1)
  test "$actual" = "$expected" || {
    echo "identity mismatch for $path: expected=$expected actual=$actual" >&2
    exit 1
  }
}

require_exact_line() {
  local path=$1
  local line=$2
  test "$(grep -Fxc -- "$line" "$path")" -eq 1 || {
    echo "missing or duplicate terminal line in $path: $line" >&2
    exit 1
  }
}

verify_resource_watch() {
  local path=$1
  local samples=0
  local terminals=0
  local line
  while IFS= read -r line; do
    if [[ $line =~ ^RESOURCE_ENVELOPE\ utc=[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}Z\ available_kib=([0-9]+)\ swap_used_kib=([0-9]+)$ ]]; then
      test "${BASH_REMATCH[1]}" -ge "$c_min_available_kib" || {
        echo "resource-watch available-memory breach in $path: $line" >&2
        exit 1
      }
      test "${BASH_REMATCH[2]}" -le "$c_max_swap_used_kib" || {
        echo "resource-watch swap breach in $path: $line" >&2
        exit 1
      }
      samples=$((samples + 1))
    elif test "$line" = "RESOURCE_ENVELOPE VERIFICATION: PASS"; then
      terminals=$((terminals + 1))
    else
      echo "malformed resource-watch line in $path: $line" >&2
      exit 1
    fi
  done < "$path"
  test "$samples" -gt 0 || {
    echo "resource-watch transcript has no samples: $path" >&2
    exit 1
  }
  test "$terminals" -eq 1 || {
    echo "resource-watch transcript has $terminals terminal markers: $path" >&2
    exit 1
  }
  test "$(tail -n 1 "$path")" = "RESOURCE_ENVELOPE VERIFICATION: PASS" || {
    echo "resource-watch terminal marker is not last: $path" >&2
    exit 1
  }
}

declare -A snapshot_manifest_visited=()

copy_manifest_snapshot() {
  local manifest_rel=$1
  local destination_root=$2
  local manifest="$root/$manifest_rel"
  local manifest_dir
  local destination_manifest_dir
  local line
  local value
  local recorded
  local source
  local destination_rel
  local -a nested_manifests=()

  if [[ ${snapshot_manifest_visited[$manifest_rel]+present} ]]; then
    return
  fi
  snapshot_manifest_visited[$manifest_rel]=1

  require_file "$manifest"
  manifest_dir=$(dirname "$manifest")
  destination_manifest_dir="$destination_root/$(dirname "$manifest_rel")"
  mkdir -p "$destination_manifest_dir"
  cp -p "$manifest" "$destination_manifest_dir/$(basename "$manifest_rel")"

  while IFS= read -r line; do
    [[ $line =~ ^([0-9a-f]{64})\ \ (.+)$ ]] || {
      echo "malformed snapshot manifest row in $manifest: $line" >&2
      exit 1
    }
    value=${BASH_REMATCH[1]}
    recorded=${BASH_REMATCH[2]}
    if [[ $recorded = /* ]]; then
      echo "absolute path in snapshot manifest $manifest: $recorded" >&2
      exit 1
    elif [[ $recorded == ginibre_q3/* ]]; then
      source=$(realpath -m "$root/../$recorded")
    elif test -e "$manifest_dir/$recorded"; then
      source=$(realpath "$manifest_dir/$recorded")
    elif test -e "$root/$recorded"; then
      source=$(realpath "$root/$recorded")
    else
      echo "missing artifact named by snapshot manifest $manifest: $recorded" >&2
      exit 1
    fi
    case "$source" in
      "$root"/*) ;;
      *)
        echo "snapshot manifest path escapes the source root: $recorded" >&2
        exit 1
        ;;
    esac
    require_sha "$source" "$value"
    destination_rel=${source#"$root/"}
    mkdir -p "$destination_root/$(dirname "$destination_rel")"
    cp -p "$source" "$destination_root/$destination_rel"
    if [[ $destination_rel == *.sha256 ]]; then
      nested_manifests+=("$destination_rel")
    fi
  done < "$manifest"

  for manifest_rel in "${nested_manifests[@]}"; do
    copy_manifest_snapshot "$manifest_rel" "$destination_root"
  done
}

remote_require_exact_line() {
  local host=$1
  local path=$2
  local line=$3
  ssh -o BatchMode=yes -o ConnectTimeout=10 "$host" \
    "test -f '$path' && test \"\$(grep -Fxc -- '$line' '$path')\" -eq 1" || {
      echo "missing remote terminal line on $host in $path: $line" >&2
      exit 1
    }
}

copy_local_program() {
  local destination=$1
  local name=$2
  cp -p "$a_root/source/character_ring_iter/$name" "$destination/$name"
}

copy_remote_program() {
  local host=$1
  local run_root=$2
  local destination=$3
  local name=$4
  scp -q -p -o BatchMode=yes -o ConnectTimeout=10 \
    "$host:$run_root/source/character_ring_iter/$name" "$destination/$name"
}

require_sha "$root/verify_full_q3_distributed_replay.py" "$verifier_sha"
require_sha "$root/run_full_q3_distributed_audit.sh" "$audit_wrapper_sha"
require_sha "$root/monitor_full_q3_machine_a_low.py" "$machine_a_monitor_sha"
require_sha "$root/monitor_full_q3_resource_envelope.py" "$resource_monitor_sha"

require_file "$a_root/machine_a_low_freeze.log"
require_file "$a_root/results/machine_a_low_logs.sha256"
require_exact_line "$a_root/machine_a_low_freeze.log" \
  "MACHINE_A_LOW FREEZE: ALL PASS"

remote_require_exact_line "$b_host" \
  "$b_root/fullq3fleetbmid0001_b.log" \
  "FULL_Q3_FLEET_PARTITION VERIFICATION: ALL PASS"
remote_require_exact_line "$b_host" \
  "$b_root/fullq3fleetbmid0001_b.log" \
  "__ORCHESTRATOR_EXIT_STATUS=0"
remote_require_exact_line "$c_host" \
  "$c_root/fullq3fleetcheavy0004_c.log" \
  "FULL_Q3_FLEET_PARTITION VERIFICATION: ALL PASS"
remote_require_exact_line "$c_host" \
  "$c_root/fullq3fleetcheavy0004_c.log" \
  "__ORCHESTRATOR_EXIT_STATUS=0"
remote_require_exact_line "$c_host" \
  "$c_root/fullq3fleetcheavy0004_resource_receipt.log" \
  "RESOURCE_ENVELOPE VERIFICATION: PASS"

for destination in \
  fleet_a_low fleet_b_mid fleet_c_heavy execution_source_snapshot
do
  test ! -e "$archive/$destination" || {
    echo "refusing to replace existing archive path: $archive/$destination" >&2
    exit 1
  }
done

staging=$(mktemp -d "$archive/.frontier-collection.XXXXXX")
snapshot_replay_parent=
cleanup() {
  rm -rf "$staging"
  if test -n "$snapshot_replay_parent"; then
    rm -rf "$snapshot_replay_parent"
  fi
}
trap cleanup EXIT

for destination in fleet_a_low fleet_b_mid fleet_c_heavy; do
  mkdir -p "$staging/$destination/results"
  mkdir -p "$staging/$destination/source/character_ring_iter"
done

# Freeze the exact source trees and composition tools consumed by the
# distributed verifier.  This snapshot is intentionally prepared before the
# prepromotion audit and moved into the archive in the same transaction as the
# fleet evidence.  A later publication-only edit can therefore be compared
# with, but can never replace, the source bytes actually replayed here.
snapshot="$staging/execution_source_snapshot/ginibre_q3"
mkdir -p "$snapshot"
require_sha "$collector_path" "$collector_start_sha"
for manifest_rel in \
  replay_sources.sha256 \
  certificates/full_q3/full_q3_source_manifest.sha256 \
  certificates/post_m29/post_m29_accepted_manifest.sha256 \
  certificates/post_m29/bc_interval_tail_code.sha256 \
  certificates/post_m29/post_m29_input_manifest.sha256 \
  certificates/post_m29/b11_delta_logs_bonly.sha256 \
  certificates/post_m29/b12_delta_logs_bonly.sha256 \
  certificates/post_m29/b13_delta_logs_bonly.sha256 \
  certificates/post_m29/c13_delta_logs_conly.sha256 \
  certificates/post_m29/c14_delta_logs_conly.sha256 \
  certificates/post_m29/c17_delta_logs_conly.sha256 \
  certificates/classical_bridge/accepted_log_manifest.sha256 \
  certificates/character_ring_pairing/character_ring_pairing_manifest.sha256 \
  references/references_manifest.sha256
do
  copy_manifest_snapshot "$manifest_rel" "$snapshot"
done

for source_rel in \
  verify_full_q3_distributed_replay.py \
  run_full_q3_distributed_audit.sh \
  collect_full_q3_distributed_frontier.sh \
  run_full_q3_frontier_fleet_queue.py \
  run_full_q3_frontier_partition.py \
  certificates/post_m29/post_m29_artifact_classification.tsv \
  character_ring_iter/verify_frontier_horizontal_strip_kernel.cpp \
  character_ring_iter/gen_header.py \
  character_ring_iter/lie_data.py \
  character_ring_iter/lie_data.h
do
  require_file "$root/$source_rel"
  mkdir -p "$snapshot/$(dirname "$source_rel")"
  cp -p "$root/$source_rel" "$snapshot/$source_rel"
done

# The classical and exceptional post-build audits authenticate the exact
# locally rebuilt executables as well as their manifested sources.  Preserve
# those build artifacts in the execution snapshot so the archived verifier
# does not depend on binaries remaining in the live publication tree.
for binary_rel in \
  character_ring_iter/post_m29_bc_layered_mgf_mpfr \
  character_ring_iter/verify_d10_exact_tail_gmp \
  character_ring_iter/post_m29_d_stable_window_interval_gmp \
  character_ring_iter/verify_character_ring_moment_sources_gmp \
  character_ring_iter/dunkl_exceptional_coefficients_gmp \
  character_ring_iter/direct_chain_rect_g2_f4_certificate \
  character_ring_iter/direct_chain_rect_e6_e7_rank_certificate \
  character_ring_iter/verify_e8_negative_bounds_gmp \
  character_ring_iter/direct_chain_rect_e8_parallel_certificate
do
  require_file "$root/$binary_rel"
  mkdir -p "$snapshot/$(dirname "$binary_rel")"
  cp -p "$root/$binary_rel" "$snapshot/$binary_rel"
done

(
  cd "$snapshot"
  find . -type f ! -name execution_source_snapshot.sha256 -print0 \
    | sort -z | xargs -0 sha256sum > execution_source_snapshot.sha256
  sha256sum -c execution_source_snapshot.sha256 >/dev/null
)
echo "EXECUTION_SOURCE_SNAPSHOT files=$(find "$snapshot" -type f | wc -l) manifests=${#snapshot_manifest_visited[@]} failures=0"

# Machine A: admit only the exact 24-log subset named by the freeze manifest.
cp -p "$a_root/build.log" "$staging/fleet_a_low/build.log"
cp -p "$a_root/fullq3cleanfleet0002_a.log" \
  "$staging/fleet_a_low/fullq3cleanfleet0002_a.log"
cp -p "$a_root/machine_a_low_freeze.log" \
  "$staging/fleet_a_low/machine_a_low_freeze.log"
cp -p "$a_root/machine_a_low_monitor.log" \
  "$staging/fleet_a_low/machine_a_low_monitor.log"
cp -p "$root/monitor_full_q3_machine_a_low.py" \
  "$staging/fleet_a_low/monitor_full_q3_machine_a_low.py"
cp -p "$a_root/run_full_q3_frontier_fleet_queue.py" \
  "$staging/fleet_a_low/run_full_q3_frontier_fleet_queue.py"
cp -p "$a_root/results/machine_a_low_logs.sha256" \
  "$staging/fleet_a_low/results/machine_a_low_logs.sha256"
require_exact_line "$staging/fleet_a_low/machine_a_low_freeze.log" \
  "__STATUS=PASS"
require_exact_line "$staging/fleet_a_low/machine_a_low_freeze.log" \
  "__RUN_ID=fullq3cleanfleet0002"
require_exact_line "$staging/fleet_a_low/machine_a_low_freeze.log" \
  "__LEDGER_SHA256=$ledger_sha"
require_exact_line "$staging/fleet_a_low/machine_a_low_freeze.log" \
  "MACHINE_A_LOW tasks_passed=24 tasks_expected=24 failures=0"
a_receipt_top_sha=$(sed -n 's/^__TOP_LOG_SHA256=//p' \
  "$staging/fleet_a_low/machine_a_low_freeze.log")
test "$(grep -c '^__TOP_LOG_SHA256=' \
  "$staging/fleet_a_low/machine_a_low_freeze.log")" -eq 1
[[ $a_receipt_top_sha =~ ^[0-9a-f]{64}$ ]] || {
  echo "malformed machine-A top-log digest in freeze receipt" >&2
  exit 1
}
require_sha "$staging/fleet_a_low/fullq3cleanfleet0002_a.log" \
  "$a_receipt_top_sha"
while read -r value name; do
  [[ $value =~ ^[0-9a-f]{64}$ && $name =~ ^[a-z0-9_]+\.log$ ]] || {
    echo "malformed machine-A subset manifest row: $value $name" >&2
    exit 1
  }
  cp -p "$a_root/results/$name" "$staging/fleet_a_low/results/$name"
done < "$a_root/results/machine_a_low_logs.sha256"
for stem in \
  post_m29_bc_twentyfifth_frontier_gmp \
  post_m29_bc_twentysixth_frontier_gmp
do
  copy_local_program "$staging/fleet_a_low/source/character_ring_iter" "$stem"
  copy_local_program "$staging/fleet_a_low/source/character_ring_iter" "$stem.cpp"
done

# Machines B and C: their result directories are already exact partitions.
scp -q -p -o BatchMode=yes -o ConnectTimeout=10 \
  "$b_host:$b_root/build.log" "$staging/fleet_b_mid/build.log"
scp -q -p -o BatchMode=yes -o ConnectTimeout=10 \
  "$b_host:$b_root/fullq3fleetbmid0001_b.log" \
  "$staging/fleet_b_mid/fullq3fleetbmid0001_b.log"
scp -q -p -o BatchMode=yes -o ConnectTimeout=10 \
  "$b_host:$b_root/run_full_q3_frontier_fleet_queue.py" \
  "$staging/fleet_b_mid/run_full_q3_frontier_fleet_queue.py"
scp -q -p -o BatchMode=yes -o ConnectTimeout=10 \
  "$b_host:$b_root/run_full_q3_frontier_partition.py" \
  "$staging/fleet_b_mid/run_full_q3_frontier_partition.py"
scp -q -pr -o BatchMode=yes -o ConnectTimeout=10 \
  "$b_host:$b_root/results/." "$staging/fleet_b_mid/results/"
copy_remote_program "$b_host" "$b_root" \
  "$staging/fleet_b_mid/source/character_ring_iter" \
  post_m29_bc_twentyseventh_frontier_gmp
copy_remote_program "$b_host" "$b_root" \
  "$staging/fleet_b_mid/source/character_ring_iter" \
  post_m29_bc_twentyseventh_frontier_gmp.cpp

scp -q -p -o BatchMode=yes -o ConnectTimeout=10 \
  "$c_host:$c_root/build.log" "$staging/fleet_c_heavy/build.log"
scp -q -p -o BatchMode=yes -o ConnectTimeout=10 \
  "$c_host:$c_root/fullq3fleetcheavy0004_c.log" \
  "$staging/fleet_c_heavy/fullq3fleetcheavy0004_c.log"
scp -q -p -o BatchMode=yes -o ConnectTimeout=10 \
  "$c_host:$c_root/fullq3fleetcheavy0004_resource_receipt.log" \
  "$staging/fleet_c_heavy/fullq3fleetcheavy0004_resource_receipt.log"
scp -q -p -o BatchMode=yes -o ConnectTimeout=10 \
  "$c_host:$c_root/fullq3fleetcheavy0004_resource_watch.log" \
  "$staging/fleet_c_heavy/fullq3fleetcheavy0004_resource_watch.log"
scp -q -p -o BatchMode=yes -o ConnectTimeout=10 \
  "$c_host:$c_root/monitor_full_q3_resource_envelope.py" \
  "$staging/fleet_c_heavy/monitor_full_q3_resource_envelope.py"
scp -q -p -o BatchMode=yes -o ConnectTimeout=10 \
  "$c_host:$c_root/run_full_q3_frontier_fleet_queue.py" \
  "$staging/fleet_c_heavy/run_full_q3_frontier_fleet_queue.py"
scp -q -p -o BatchMode=yes -o ConnectTimeout=10 \
  "$c_host:$c_root/run_full_q3_frontier_partition.py" \
  "$staging/fleet_c_heavy/run_full_q3_frontier_partition.py"
scp -q -pr -o BatchMode=yes -o ConnectTimeout=10 \
  "$c_host:$c_root/results/." "$staging/fleet_c_heavy/results/"
for stem in \
  post_m29_bc_twentyseventh_frontier_gmp \
  post_m29_bc_twentyeighth_frontier_gmp \
  post_m29_bc_twentyeighth_b_absorption_gmp \
  post_m29_bc_twentyninth_c_frontier_gmp \
  post_m29_bc_twentyninth_b_absorption_gmp
do
  copy_remote_program "$c_host" "$c_root" \
    "$staging/fleet_c_heavy/source/character_ring_iter" "$stem"
  copy_remote_program "$c_host" "$c_root" \
    "$staging/fleet_c_heavy/source/character_ring_iter" "$stem.cpp"
done

test "$(find "$staging/fleet_a_low/results" -maxdepth 1 -type f -name '*.log' | wc -l)" -eq 24
test "$(find "$staging/fleet_b_mid/results" -maxdepth 1 -type f -name '*.log' | wc -l)" -eq 7
test "$(find "$staging/fleet_c_heavy/results" -maxdepth 1 -type f -name '*.log' | wc -l)" -eq 40
(cd "$staging/fleet_a_low/results" && sha256sum -c machine_a_low_logs.sha256)
(cd "$staging/fleet_b_mid/results" && sha256sum -c fullq3fleetbmid0001_logs.sha256)
(cd "$staging/fleet_c_heavy/results" && sha256sum -c fullq3fleetcheavy0004_logs.sha256)

test "$(sha256sum "$staging/fleet_a_low/run_full_q3_frontier_fleet_queue.py" | cut -d' ' -f1)" = "$ledger_sha"
require_sha "$staging/fleet_a_low/monitor_full_q3_machine_a_low.py" \
  "$machine_a_monitor_sha"
require_exact_line "$staging/fleet_a_low/machine_a_low_monitor.log" \
  "MACHINE_A_LOW FREEZE: ALL PASS"
for destination in fleet_b_mid fleet_c_heavy; do
  test "$(sha256sum "$staging/$destination/run_full_q3_frontier_fleet_queue.py" | cut -d' ' -f1)" = "$ledger_sha"
  test "$(sha256sum "$staging/$destination/run_full_q3_frontier_partition.py" | cut -d' ' -f1)" = "$wrapper_sha"
done
require_sha "$staging/fleet_c_heavy/monitor_full_q3_resource_envelope.py" \
  "$resource_monitor_sha"
require_exact_line "$staging/fleet_c_heavy/fullq3fleetcheavy0004_resource_watch.log" \
  "RESOURCE_ENVELOPE VERIFICATION: PASS"
require_exact_line "$staging/fleet_c_heavy/fullq3fleetcheavy0004_resource_receipt.log" \
  "__STATUS=PROCESS_GROUP_EXITED_WITHIN_ENVELOPE"
require_exact_line "$staging/fleet_c_heavy/fullq3fleetcheavy0004_resource_receipt.log" \
  "__RUN_TOKEN=fullq3fleetcheavy0004"
require_exact_line "$staging/fleet_c_heavy/fullq3fleetcheavy0004_resource_receipt.log" \
  "RESOURCE_ENVELOPE VERIFICATION: PASS"
verify_resource_watch \
  "$staging/fleet_c_heavy/fullq3fleetcheavy0004_resource_watch.log"

for destination in fleet_a_low fleet_b_mid fleet_c_heavy; do
  for source in "$staging/$destination/source/character_ring_iter/"*.cpp; do
    cmp -s "$source" "$root/character_ring_iter/$(basename "$source")" || {
      echo "staged/current source mismatch: $source" >&2
      exit 1
    }
  done
done

# Exercise the complete verifier against the staged paths before any archive
# destination is made visible.  The final wrapper reruns the same audit after
# atomic relocation.
preaudit="$staging/prepromotion_composition.log"
python3 "$root/verify_full_q3_distributed_replay.py" \
  --root "$root" \
  --current-prefix-log "$archive/current_prefix/current_source_machine_c.log" \
  --current-prefix-stages "$archive/current_prefix/replay-logs" \
  --legacy-root "$archive/legacy_c_v6" \
  --legacy-log "$archive/legacy_c_v6/full_replay_v6.log" \
  --legacy-stages "$archive/legacy_c_v6/replay-logs" \
  --fleet-ledger-script "$root/run_full_q3_frontier_fleet_queue.py" \
  --fleet-partition-script "$root/run_full_q3_frontier_partition.py" \
  --fleet-a-build-log "$staging/fleet_a_low/build.log" \
  --fleet-a-log "$staging/fleet_a_low/fullq3cleanfleet0002_a.log" \
  --fleet-a-results "$staging/fleet_a_low/results" \
  --fleet-a-manifest "$staging/fleet_a_low/results/machine_a_low_logs.sha256" \
  --fleet-b-build-log "$staging/fleet_b_mid/build.log" \
  --fleet-b-log "$staging/fleet_b_mid/fullq3fleetbmid0001_b.log" \
  --fleet-b-results "$staging/fleet_b_mid/results" \
  --fleet-b-manifest "$staging/fleet_b_mid/results/fullq3fleetbmid0001_logs.sha256" \
  --fleet-c-build-log "$staging/fleet_c_heavy/build.log" \
  --fleet-c-log "$staging/fleet_c_heavy/fullq3fleetcheavy0004_c.log" \
  --fleet-c-results "$staging/fleet_c_heavy/results" \
  --fleet-c-manifest "$staging/fleet_c_heavy/results/fullq3fleetcheavy0004_logs.sha256" \
  --classical-build-log "$archive/current_classical/current_classical_build.log" \
  --classical-log "$archive/current_classical/current_classical_replay.log" \
  --classical-stages "$archive/current_classical/stages" \
  --exceptional-build-log "$archive/current_exceptional/current_exceptional_build.log" \
  --exceptional-log "$archive/current_exceptional/current_exceptional_replay.log" \
  --exceptional-stages "$archive/current_exceptional/stages" \
  > "$preaudit" 2>&1
require_exact_line "$preaudit" "FULL_Q3_DISTRIBUTED VERIFICATION: ALL PASS"
require_exact_line "$preaudit" "__EXIT_STATUS=0"
rm -f "$preaudit"

mv "$staging/fleet_a_low" "$archive/fleet_a_low"
mv "$staging/fleet_b_mid" "$archive/fleet_b_mid"
mv "$staging/fleet_c_heavy" "$archive/fleet_c_heavy"
mv "$staging/execution_source_snapshot" "$archive/execution_source_snapshot"

# Prove that the archived source snapshot is itself sufficient to replay the
# complete composition.  This is distinct from merely hashing the snapshot:
# the archived verifier is executed with the archived source tree as --root
# and with the final, relocated fleet paths.  The transcript remains in the
# top-level archive and is covered by the final archive manifest.
snapshot_root="$archive/execution_source_snapshot/ginibre_q3"
snapshot_replay_parent=$(mktemp -d /tmp/fullq3-execution-source-replay.XXXXXX)
snapshot_replay_root="$snapshot_replay_parent/ginibre_q3"
cp -a "$snapshot_root" "$snapshot_replay_root"
(
  cd "$snapshot_replay_root"
  sha256sum -c execution_source_snapshot.sha256 >/dev/null
)
snapshot_audit_tmp="$archive/fullq3distributed0001_execution_snapshot.log.tmp"
if ! python3 "$snapshot_replay_root/verify_full_q3_distributed_replay.py" \
  --root "$snapshot_replay_root" \
  --current-prefix-log "$archive/current_prefix/current_source_machine_c.log" \
  --current-prefix-stages "$archive/current_prefix/replay-logs" \
  --legacy-root "$archive/legacy_c_v6" \
  --legacy-log "$archive/legacy_c_v6/full_replay_v6.log" \
  --legacy-stages "$archive/legacy_c_v6/replay-logs" \
  --fleet-ledger-script "$snapshot_replay_root/run_full_q3_frontier_fleet_queue.py" \
  --fleet-partition-script "$snapshot_replay_root/run_full_q3_frontier_partition.py" \
  --fleet-a-build-log "$archive/fleet_a_low/build.log" \
  --fleet-a-log "$archive/fleet_a_low/fullq3cleanfleet0002_a.log" \
  --fleet-a-results "$archive/fleet_a_low/results" \
  --fleet-a-manifest "$archive/fleet_a_low/results/machine_a_low_logs.sha256" \
  --fleet-b-build-log "$archive/fleet_b_mid/build.log" \
  --fleet-b-log "$archive/fleet_b_mid/fullq3fleetbmid0001_b.log" \
  --fleet-b-results "$archive/fleet_b_mid/results" \
  --fleet-b-manifest "$archive/fleet_b_mid/results/fullq3fleetbmid0001_logs.sha256" \
  --fleet-c-build-log "$archive/fleet_c_heavy/build.log" \
  --fleet-c-log "$archive/fleet_c_heavy/fullq3fleetcheavy0004_c.log" \
  --fleet-c-results "$archive/fleet_c_heavy/results" \
  --fleet-c-manifest "$archive/fleet_c_heavy/results/fullq3fleetcheavy0004_logs.sha256" \
  --classical-build-log "$archive/current_classical/current_classical_build.log" \
  --classical-log "$archive/current_classical/current_classical_replay.log" \
  --classical-stages "$archive/current_classical/stages" \
  --exceptional-build-log "$archive/current_exceptional/current_exceptional_build.log" \
  --exceptional-log "$archive/current_exceptional/current_exceptional_replay.log" \
  --exceptional-stages "$archive/current_exceptional/stages" \
  > "$snapshot_audit_tmp" 2>&1
then
  mv "$snapshot_audit_tmp" \
    "$archive/fullq3distributed0001_execution_snapshot_failed.log"
  echo "archived-source composition failed; evidence retained for diagnosis" >&2
  exit 1
fi
mv "$snapshot_audit_tmp" \
  "$archive/fullq3distributed0001_execution_snapshot.log"
rm -rf "$snapshot_replay_parent"
snapshot_replay_parent=
require_exact_line "$archive/fullq3distributed0001_execution_snapshot.log" \
  "FULL_Q3_DISTRIBUTED VERIFICATION: ALL PASS"
require_exact_line "$archive/fullq3distributed0001_execution_snapshot.log" \
  "__EXIT_STATUS=0"

audit_tmp="$archive/fullq3distributed0001.log.tmp"
if ! "$root/run_full_q3_distributed_audit.sh" > "$audit_tmp" 2>&1; then
  mv "$audit_tmp" "$archive/fullq3distributed0001_failed.log"
  echo "distributed composition failed; evidence retained for diagnosis" >&2
  exit 1
fi
mv "$audit_tmp" "$archive/fullq3distributed0001.log"
require_exact_line "$archive/fullq3distributed0001.log" \
  "FULL_Q3_DISTRIBUTED VERIFICATION: ALL PASS"
require_exact_line "$archive/fullq3distributed0001.log" "__EXIT_STATUS=0"
(cd "$archive" && sha256sum fullq3distributed0001.log \
  > fullq3distributed0001.log.sha256)

{
  echo "__UTC=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  require_sha "$collector_path" "$collector_start_sha"
  echo "__COLLECTOR_SHA256=$collector_start_sha"
  echo "__LEDGER_SHA256=$ledger_sha"
  echo "__PARTITION_WRAPPER_SHA256=$wrapper_sha"
  echo "__VERIFIER_SHA256=$verifier_sha"
  echo "__AUDIT_WRAPPER_SHA256=$audit_wrapper_sha"
  echo "__MACHINE_A_MONITOR_SHA256=$machine_a_monitor_sha"
  echo "__RESOURCE_MONITOR_SHA256=$resource_monitor_sha"
  echo "FULL_Q3_DISTRIBUTED_COLLECTION fleet_a=24 fleet_b=7 fleet_c=40 source_snapshot_replay=PASS failures=0"
  echo "FULL_Q3_DISTRIBUTED COLLECTION: ALL PASS"
} > "$archive/fullq3distributed0001_collection.log"

manifest_tmp=$(mktemp /tmp/fullq3distributed0001_manifest.XXXXXX)
(cd "$archive" && find . -type f \
  ! -name 'distributed0001_manifest.sha256' -print0 \
  | sort -z | xargs -0 sha256sum) > "$manifest_tmp"
mv "$manifest_tmp" "$archive/distributed0001_manifest.sha256"

echo "FULL_Q3_DISTRIBUTED COLLECTION AND VERIFICATION: ALL PASS"
