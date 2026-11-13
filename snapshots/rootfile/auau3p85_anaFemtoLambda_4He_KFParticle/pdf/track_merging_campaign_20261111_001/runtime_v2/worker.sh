#!/bin/bash
set -euo pipefail
[[ $# -eq 4 ]] || { echo 'Usage: worker.sh CACHE_ROOT ORIGINAL_PICO_LIST OUTPUT_DIR POLICY_JSON' >&2; exit 2; }
mc_runtime="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
mc_scratch="${_CONDOR_SCRATCH_DIR:-/tmp}"
exec /cvmfs/star.sdcc.bnl.gov/star-spack/spack/opt/spack/linux-rhel7-x86_64/gcc-4.8.5/python-3.8.13-myixjs7jesskpbhwiowibizezfh52rky/bin/python3 "$mc_runtime/worker.py" --runtime "$mc_runtime" --cache "$1" --pico-list "$2" --output-dir "$3" --policy "$4" --scratch-base "$mc_scratch"
