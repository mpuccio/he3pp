# He3pp Python CLI

Unified, config-driven PyROOT workflow for anti-He3/anti-He4 analysis.

## Current structure

- `he3_cli.py`: thin entrypoint.
- `he3pp/cli.py`: config loading + task dispatch.
- `he3pp/settings.py`: defaults and runtime overrides from TOML.
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

Outputs are organized per variant under:
`$NUCLEI_OUTPUT/<period>/<reco_pass>/<variant>/`
including analysis ROOT files, metadata, and report assets.

Show default merged config:

```bash
python3 he3_cli.py --config config.example.toml --dump-default-config
```

## Config schema

- `[common]`: former `Common.h` style constants (period, reco pass, pt bins, ranges, etc.)
  - Input paths are derived from `period`/`reco_pass` (data) and `mc_production` (MC) plus optional basename keys:
  `data_tree_basename`, `data_analysis_results_basename`, `mc_tree_basename`, `mc_analysis_results_basename`
- `[selections.common]`: shared selection snippets (e.g. skim template)
- `[selections.he3]`: he3-specific selections (required when `he3` is in `run.species`)
- `[selections.he4]`: he4-specific selections (required when `he4` is in `run.species`)
- `[cuts]`: trial scan grids (`nsigmaDCAz`, `fTPCnCls`, `nITScls`, `nsigmaTPC`)
- `[run]`: task + runtime flags
- `[paths]`: optional non-derivable overrides (most IO paths are auto-derived)

Logging and metadata:

- `[run].log_level`: `DEBUG|INFO|WARNING|ERROR`
- `[paths].log_file`: optional log file path
- `[paths].metadata_output`: JSON metadata output for each run

Report controls:

- `[report].sections`: ordered list of report sections to render
- `[report].fit_n_parameters`: free-parameter count used for chi2 NDF estimate
- `[report].fit_alpha`: Pearson threshold for signal-fit `OK/KO` labels
- `[report].fit_tail`: `single` (default) or `two` for p-value computation
- `[report].tpc_signal_model`: TPC-only model used for extraction plots + summary table
- `[run].species`: processing/report species list (e.g. `["he3"]` or `["he3","he4"]`)
- `[species.<name>.paths]`: optional species-specific path overrides (not required)

Available report sections include:
`signal_tof`, `signal_tpc`, `tof_tpc_2d`, `efficiency`, `pt_resolution`, `corrected_spectrum`.

When two species are requested:
- Both selection sections must exist: `[selections.he3]` and `[selections.he4]`
- Generates per-species subreports in `<report_dir>/he3/` and `<report_dir>/he4/`
- Generates a top-level index at `<report_dir>/index.html` linking both pages
- If provided, uses species-specific path overrides from sections like:
`[species.he3.paths]` and `[species.he4.paths]` (`data_input`, `mc_input`, `signal_input`, `systematics_input`, `metadata_output`)

## Tasks

Supported tasks:

`merge_trees | analyse_data | analyse_mc | signal | systematics | checkpoint | report | full_chain`

## Notes

- Uses **PyROOT only** (no uproot/pandas).
- `signal` still uses `src/` RooFit C++ components loaded by PyROOT.
- Weighted efficiency histograms follow the Python naming policy `Weff*`.

## Smoke validation

Dual smoke (LHC24 data + LHC25b9 MC, He3 + He4):

```bash
scripts/update_smoke_references.sh   # one-time or when intentionally refreshing references
scripts/run_smoke_checks.sh
```

Reference location used by smoke checks:
`$NUCLEI_INPUT/smoke_references/python_dual/LHC24_apass1__LHC25b9/`
