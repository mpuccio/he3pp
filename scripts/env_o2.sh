#!/usr/bin/env bash
# Intentionally avoid `set -euo pipefail` here: this file is sourced in an
# interactive shell and must not leak strict-mode flags to the caller.

# Source this file: `source scripts/env_o2.sh`
is_sourced=0
if [[ -n "${BASH_VERSION-}" ]]; then
  [[ "${BASH_SOURCE[0]-}" != "$0" ]] && is_sourced=1
elif [[ -n "${ZSH_VERSION-}" ]]; then
  [[ "${ZSH_EVAL_CONTEXT-}" == *:file* ]] && is_sourced=1
fi

if [[ "$is_sourced" -ne 1 ]]; then
  echo "Source this script instead of executing it: source scripts/env_o2.sh"
  exit 1
fi

alienv load O2/latest

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
