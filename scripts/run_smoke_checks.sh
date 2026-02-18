#!/usr/bin/env bash
set -euo pipefail

if [[ -z "${NUCLEI_INPUT:-}" || -z "${NUCLEI_OUTPUT:-}" ]]; then
  echo "NUCLEI_INPUT/NUCLEI_OUTPUT are not set. Run: source scripts/env_o2.sh"
  exit 1
fi

# Ensure ROOT dynamic libs are available even if caller did not export them.
ROOT_LIB_DEFAULT="/Users/mpuccio/alice/sw/osx_arm64/ROOT/v6-36-04-alice9-local3/lib"
export ROOT_LIB="${ROOT_LIB:-$ROOT_LIB_DEFAULT}"
export CPPYY_BACKEND_LIBRARY="${CPPYY_BACKEND_LIBRARY:-$ROOT_LIB/libcppyy_backend.so}"
export DYLD_LIBRARY_PATH="$ROOT_LIB:${DYLD_LIBRARY_PATH:-}"

mkdir -p /tmp/he3pp_validation

python3 he3_cli.py --config tests/smoke/validate_data.toml
python3 scripts/validate_smoke.py \
  --reference "$NUCLEI_INPUT/data/LHC22/apass4/DataHistosgiovanni.root" \
  --candidate /tmp/he3pp_validation/DataHistosgiovanni_py.root

python3 he3_cli.py --config tests/smoke/validate_mc.toml
python3 scripts/validate_smoke.py \
  --reference "$NUCLEI_INPUT/MC/LHC23j6b/MChistosgiovanni.root" \
  --candidate /tmp/he3pp_validation/MChistosgiovanni_py.root \
  --normalize-weff \
  --ignore-errors

echo "Smoke checks passed"
