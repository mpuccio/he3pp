# Migration Notes: C++ Macros to Modular Python

## Status

The workflow is now Python-first (PyROOT) and uses a modular package:

- `he3pp/cli.py`
- `he3pp/settings.py`
- `he3pp/root_io.py`
- `he3pp/tasks.py`

Legacy C++ analysis macro entrypoints were removed.

## Implemented tasks in Python

- `merge_trees`
- `analyse_data` (He3/He4)
- `analyse_mc` (He3/He4)
- `signal`
- `systematics`
- `checkpoint`
- `full_chain`

## Configuration model

Configurable from TOML:

- `[common]`: constants, binning, ranges, naming.
- `[selections]`: selection/filter strings for data and MC.
- `[cuts]`: trial cut grids.
- `[run]`: runtime flags/task selection.
- `[paths]`: inputs/outputs.

Reference template: `config.example.toml`.

## Remaining C++ dependency

Python `signal` currently loads custom RooFit components from `src/` via PyROOT for model parity/performance.

## Validation summary

Using real inputs under `/Users/mpuccio/alice/run3/nuclei-analysis`:

- Data output matches reference exactly.
- MC output matches in content (up to floating precision), with differences in:
  - weighted-efficiency naming convention (Python standard: `Weff*`),
  - uncertainty treatment in efficiencies (Python uses binomial errors where configured).

## New quality-of-life additions

- Task-level logging via `[run].log_level` and optional `[paths].log_file`.
- Run metadata JSON output via `[paths].metadata_output`.
- Automated smoke comparator: `scripts/validate_smoke.py`.
- Automated smoke runner: `scripts/run_smoke_checks.sh`.
- O2/ROOT environment helper: `scripts/env_o2.sh`.

## Environment notes (macOS/O2)

```bash
eval "$(alienv load O2/latest)"
export ROOT_LIB=/Users/mpuccio/alice/sw/osx_arm64/ROOT/v6-36-04-alice9-local3/lib
export CPPYY_BACKEND_LIBRARY="$ROOT_LIB/libcppyy_backend.so"
export DYLD_LIBRARY_PATH="$ROOT_LIB:$DYLD_LIBRARY_PATH"
```
