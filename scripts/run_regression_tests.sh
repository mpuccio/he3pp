#!/usr/bin/env bash
set -euo pipefail

# Ensure ROOT Python backend can be resolved in environments where only
# alienv module vars are loaded.
ROOT_LIB_DEFAULT="${ROOTSYS:+$ROOTSYS/lib}"
export ROOT_LIB="${ROOT_LIB:-$ROOT_LIB_DEFAULT}"
export CPPYY_BACKEND_LIBRARY="${CPPYY_BACKEND_LIBRARY:-$ROOT_LIB/libcppyy_backend.so}"
export DYLD_LIBRARY_PATH="$ROOT_LIB:${DYLD_LIBRARY_PATH:-}"

python3 -m unittest discover -s tests/regression -p 'test_*.py' -v
