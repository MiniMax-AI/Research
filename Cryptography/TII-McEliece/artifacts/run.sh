#!/usr/bin/env bash
set -euo pipefail

exec "$(dirname "$0")/goppa_key_recovery/run.sh" "$@"
