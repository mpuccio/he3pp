#!/usr/bin/env bash
set -euo pipefail

if [[ -z "${NUCLEI_INPUT:-}" || -z "${NUCLEI_OUTPUT:-}" ]]; then
  echo "NUCLEI_INPUT/NUCLEI_OUTPUT are not set. Run: source scripts/env_o2.sh"
  exit 1
fi

ROOT_LIB_DEFAULT="/Users/mpuccio/alice/sw/osx_arm64/ROOT/v6-36-04-alice9-local3/lib"
export ROOT_LIB="${ROOT_LIB:-$ROOT_LIB_DEFAULT}"
export CPPYY_BACKEND_LIBRARY="${CPPYY_BACKEND_LIBRARY:-$ROOT_LIB/libcppyy_backend.so}"
export DYLD_LIBRARY_PATH="$ROOT_LIB:${DYLD_LIBRARY_PATH:-}"

mkdir -p /tmp/he3pp_validation

REF_ROOT="$NUCLEI_INPUT/smoke_references/python_single/LHC24_apass1__LHC25b9"
mkdir -p "$REF_ROOT"

DATA_BASE="$NUCLEI_INPUT/data/LHC24/apass1"
for candidate in AO2D.root AO2D_skimmed.root AO2D_sampled.root; do
  if [[ -f "$DATA_BASE/$candidate" ]]; then
    DATA_INPUT_FILE="$DATA_BASE/$candidate"
    break
  fi
done
if [[ -z "${DATA_INPUT_FILE:-}" ]]; then
  echo "Could not find a valid LHC24 apass1 data input under: $DATA_BASE"
  exit 1
fi
export DATA_INPUT_FILE

python3 - <<'PY'
import os
from he3pp import tasks
from he3pp.settings import current_runtime_config

cfg = current_runtime_config()

tasks.analyse_data(
    os.environ["DATA_INPUT_FILE"],
    "/tmp/he3pp_validation/DataHistos_he3_py.root",
    "he3",
    skim=False,
    draw=False,
    runtime_config=cfg,
)
PY
cp /tmp/he3pp_validation/DataHistos_he3_py.root "$REF_ROOT/DataHistos_he3_ref.root"

python3 - <<'PY'
import os
from he3pp import tasks
from he3pp.settings import current_runtime_config

cfg = current_runtime_config()

tasks.analyse_mc(
    os.path.expandvars("$NUCLEI_INPUT/MC/LHC25b9/AO2D_coalescence.root"),
    "/tmp/he3pp_validation/MChistos_he3_py.root",
    "he3",
    enable_trials=True,
    draw=False,
    runtime_config=cfg,
)
PY
cp /tmp/he3pp_validation/MChistos_he3_py.root "$REF_ROOT/MChistos_he3_ref.root"

echo "Updated smoke references under: $REF_ROOT"
