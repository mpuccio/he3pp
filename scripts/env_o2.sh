#!/usr/bin/env bash
set -euo pipefail

# Source this file: `source scripts/env_o2.sh`
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
  echo "Source this script instead of executing it: source scripts/env_o2.sh"
  exit 1
fi

eval "$(alienv load O2/latest)"

ROOT_LIB_DEFAULT="/Users/mpuccio/alice/sw/osx_arm64/ROOT/v6-36-04-alice9-local3/lib"
export ROOT_LIB="${ROOT_LIB:-$ROOT_LIB_DEFAULT}"
export CPPYY_BACKEND_LIBRARY="${CPPYY_BACKEND_LIBRARY:-$ROOT_LIB/libcppyy_backend.so}"
export DYLD_LIBRARY_PATH="$ROOT_LIB:${DYLD_LIBRARY_PATH:-}"

export NUCLEI_INPUT="${NUCLEI_INPUT:-/Users/mpuccio/alice/run3/nuclei-analysis/input}"
export NUCLEI_OUTPUT="${NUCLEI_OUTPUT:-/Users/mpuccio/alice/run3/nuclei-analysis/output}"

echo "Loaded O2 environment"
echo "ROOT_LIB=$ROOT_LIB"
echo "NUCLEI_INPUT=$NUCLEI_INPUT"
echo "NUCLEI_OUTPUT=$NUCLEI_OUTPUT"
