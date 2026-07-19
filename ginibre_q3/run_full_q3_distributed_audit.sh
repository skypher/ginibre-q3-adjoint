#!/usr/bin/env bash
set -eu

root=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
archive="$root/certificates/full_q3/distributed0001"

python3 "$root/verify_full_q3_distributed_replay.py" \
  --root "$root" \
  --current-prefix-log "$archive/current_prefix/current_source_machine_c.log" \
  --current-prefix-stages "$archive/current_prefix/replay-logs" \
  --legacy-root "$archive/legacy_c_v6" \
  --legacy-log "$archive/legacy_c_v6/full_replay_v6.log" \
  --legacy-stages "$archive/legacy_c_v6/replay-logs" \
  --fleet-ledger-script "$root/run_full_q3_frontier_fleet_queue.py" \
  --fleet-partition-script "$root/run_full_q3_frontier_partition.py" \
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
  --exceptional-stages "$archive/current_exceptional/stages"
