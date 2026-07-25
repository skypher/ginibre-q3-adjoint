#!/usr/bin/env bash
# Activate the repository hooks in this clone.
#
# core.hooksPath is local configuration, so every clone (and every agent
# working directory) must run this once.  Without it the pre-commit manifest
# guard is inert.
set -euo pipefail
cd "$(dirname "$0")/.."
git config core.hooksPath .githooks
chmod +x .githooks/pre-commit
echo "hooks active: core.hooksPath=$(git config core.hooksPath)"
