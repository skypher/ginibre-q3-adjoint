#!/usr/bin/env bash
set -Eeuo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(git -C "$script_dir" rev-parse --show-toplevel)
log_path=${FINAL_REPLAY_LOG:-}

if [[ -z "$log_path" ]]; then
    echo "FINAL_PUBLICATION_REPLAY: FAIL: set FINAL_REPLAY_LOG to an absolute output path" >&2
    exit 2
fi
if [[ "$log_path" != /* ]]; then
    echo "FINAL_PUBLICATION_REPLAY: FAIL: FINAL_REPLAY_LOG must be absolute" >&2
    exit 2
fi
if [[ -n "$(git -C "$repo_root" status --porcelain=v1 --untracked-files=all)" ]]; then
    echo "FINAL_PUBLICATION_REPLAY: FAIL: the final-source worktree is not clean" >&2
    exit 2
fi
if [[ -e "$log_path" ]]; then
    echo "FINAL_PUBLICATION_REPLAY: FAIL: refusing to overwrite $log_path" >&2
    exit 2
fi

mkdir -p -- "$(dirname -- "$log_path")"
exec > >(tee -a -- "$log_path") 2>&1

run_stage() {
    local stage=$1
    shift
    printf 'FINAL_PUBLICATION_REPLAY stage=%s event=start utc=%s\n' \
        "$stage" "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    "$@"
    printf 'FINAL_PUBLICATION_REPLAY stage=%s event=pass utc=%s\n' \
        "$stage" "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
}

printf 'FINAL_PUBLICATION_REPLAY commit=%s\n' "$(git -C "$repo_root" rev-parse HEAD)"
printf 'FINAL_PUBLICATION_REPLAY source_manifest_sha256=%s\n' \
    "$(sha256sum "$script_dir/certificates/full_q3/full_q3_source_manifest.sha256" | awk '{print $1}')"
printf 'FINAL_PUBLICATION_REPLAY source_date_epoch=%s\n' \
    "${PUBLICATION_SOURCE_DATE_EPOCH:-1784419200}"
printf 'FINAL_PUBLICATION_REPLAY host=%s\n' "$(uname -a)"
printf 'FINAL_PUBLICATION_REPLAY compiler=%s\n' "$(g++ --version | head -n 1)"
printf 'FINAL_PUBLICATION_REPLAY python=%s\n' "$(python3 --version 2>&1)"

run_stage publication-preflight \
    make -C "$script_dir" publication-preflight
run_stage parts-i-ii-clean-room \
    make -C "$script_dir" clean-room-replay \
        REPLAY_THREADS="${REPLAY_THREADS:-$(nproc)}"
run_stage part-iii-full-replay \
    make -C "$script_dir" full-q3-extension \
        REPLAY_THREADS="${REPLAY_THREADS:-$(nproc)}" \
        FULL_Q3_BOUNDED_THREADS="${FULL_Q3_BOUNDED_THREADS:-$(nproc)}"

if [[ -n "$(git -C "$repo_root" status --porcelain=v1 --untracked-files=all)" ]]; then
    echo "FINAL_PUBLICATION_REPLAY: FAIL: replay changed tracked or untracked source-tree files"
    git -C "$repo_root" status --short
    exit 1
fi

printf 'FINAL_PUBLICATION_REPLAY: ALL PASS utc=%s\n' \
    "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
