#!/bin/bash
set -e
cd "$(dirname "$0")/unified_collision"
exec ./run.sh "$@"
