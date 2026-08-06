#!/bin/bash
# run.sh - thin shim that delegates to the Python wrapper run.
#
# Usage:  bash run.sh [OPTIONS]
# Run    'python run --help'  for the full list of options.
#
# Example (Slurm):
#   sbatch --cpus-per-task=64 --mem=256G run.sh --slurm --cores 64 --samples sample1

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

exec python3 "${SCRIPT_DIR}/run" "$@"
