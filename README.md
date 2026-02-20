# He3pp Python CLI

Unified, config-driven PyROOT workflow for anti-He3/anti-He4 analysis.

## Current structure

- `he3_cli.py`: thin entrypoint.
- `he3pp/cli.py`: config loading + task dispatch.
- `he3pp/defaults.toml`: canonical default configuration (common, selections, cuts, run/report, particle profiles).
- `he3pp/settings.py`: runtime settings/state loaded from `he3pp/defaults.toml` + user TOML overrides.
- `he3pp/root_io.py`: ROOT/PyROOT helpers and column definitions.
- `he3pp/tasks.py`: analysis tasks.
- `config.example.toml`: full config template.

Legacy C++ macro entrypoints have been removed.

## Environment setup

Use the helper script (recommended):

```bash
source scripts/env_o2.sh
```

Or manually load O2 in the current shell and export ROOT runtime paths:

```bash
eval "$(alienv load O2/latest)"
export ROOT_LIB=/Users/mpuccio/alice/sw/osx_arm64/ROOT/v6-36-04-alice9-local3/lib
export CPPYY_BACKEND_LIBRARY="$ROOT_LIB/libcppyy_backend.so"
export DYLD_LIBRARY_PATH="$ROOT_LIB:$DYLD_LIBRARY_PATH"
```

Quick check:

```bash
python3 -c 'import ROOT; print(ROOT.gROOT.GetVersion())'
```

## Usage

```bash
python3 he3_cli.py --config config.example.toml
```

Outputs are organized per variant and species under:
`$NUCLEI_OUTPUT/<period>/<reco_pass>/<variant>/<species>/`
including analysis ROOT files, metadata, and report assets.

Show default merged config:

```bash
python3 he3_cli.py --config config.example.toml --dump-default-config
```

## Config schema

- Built-in defaults live in `he3pp/defaults.toml`; your `--config` file only needs to override what you want to change.
- `[common]`: former `Common.h` style constants (period, reco pass, pt bins, ranges, etc.)
  - Input paths are derived from `period`/`reco_pass` (data) and `mc_production` (MC) plus optional basename keys:
  `data_tree_basename`, `data_analysis_results_basename`, `mc_tree_basename`, `mc_analysis_results_basename`
- `[selections.common]`: shared selection snippets (e.g. skim template)
- `[selections.<species>]`: species-specific selections (required for the species selected in `[run].species`)
- `[cuts]`: trial scan grids (`nsigmaDCAz`, `fTPCnCls`, `nITScls`, `nsigmaTPC`)
- `[run]`: task + runtime flags
- `[particle.<species>]`: particle-profile definitions/overrides (mass, PDG, labels, key column names); extra species can be added via `template = "he3"` or `template = "he4"`
- `[paths]`: optional overrides only (all standard IO paths are auto-derived from `[common]` + `run.species`)

Logging and metadata:

- `[run].log_level`: `DEBUG|INFO|WARNING|ERROR`
- `[paths].log_file`: optional log file path
- `[paths].metadata_output`: optional JSON metadata output path (otherwise auto-derived per species)

Report controls:

- `[report].sections`: ordered list of report sections to render
- `[report].fit_alpha`: Pearson threshold for signal-fit `OK/KO` labels
- `[report].fit_tail`: `single` (default) or `two` for p-value computation
- `[report].tpc_signal_model`: TPC-only model used for extraction plots + summary table
- `[run].species`: single processing species key matching one section in `[particle]`

Available report sections include:
`signal_tof`, `signal_tpc`, `tof_tpc_2d`, `efficiency`, `pt_resolution`, `corrected_spectrum`.

Single-particle mode:
- Select one species with `[run].species`
- Provide the matching selection section (`[selections.<species>]`)
- Use task-specific keys under `[paths]` only when you need custom routing (`data_tree`, `mc_tree`, `*_output`, `*_input`, `report_dir`, etc.)

## Tasks

Supported tasks:

`analyse_data | analyse_mc | signal | systematics | checkpoint | report | full_chain`

## Notes

- Uses **PyROOT only** (no uproot/pandas).
- `signal` still uses `src/` RooFit C++ components loaded by PyROOT.
- Weighted efficiency histograms follow the Python naming policy `Weff*`.

## Smoke validation

Smoke (LHC24 data + LHC25b9 MC, He3):

```bash
scripts/update_smoke_references.sh   # one-time or when intentionally refreshing references
scripts/run_smoke_checks.sh
```

Reference location used by smoke checks:
`$NUCLEI_INPUT/smoke_references/python_single/LHC24_apass1__LHC25b9/`

Legacy `tests/smoke/*.toml` files were retired; smoke is driven directly by the shell scripts above.

## Regression tests

Lightweight regression checks (kept separate from smoke):

```bash
scripts/run_regression_tests.sh
```

Coverage includes:
- config/path derivation
- TPC model production matrix mapping/validation
- report fit-status classification (`OK`/`KO`/`UNK`/`MISSING`)
- MC trial variation changing at least one trial histogram vs default (using smoke reference MC ROOT file)
