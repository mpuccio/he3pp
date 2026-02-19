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
- `[selections]`: all data/MC filter expressions
- `[cuts]`: trial scan grids (`nsigmaDCAz`, `fTPCnCls`, `nITScls`, `nsigmaTPC`)
- `[run]`: task + runtime flags
- `[paths]`: all input/output paths

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
- `[run].report_mode`: `single` (default) or `dual` for `task=report`
- `[run].report_species_dual`: species list for dual report mode (typically `["antihe3","antihe4"]`)

Available report sections include:
`signal_tof`, `signal_tpc`, `tof_tpc_2d`, `efficiency`, `pt_resolution`, `corrected_spectrum`.

Dual report mode:
- Generates per-species subreports in `<report_dir>/he3/` and `<report_dir>/he4/`
- Generates a top-level index at `<report_dir>/index.html` linking both pages
- Uses species-specific input paths when provided:
`data_input_he3/he4`, `mc_input_he3/he4`, `signal_input_he3/he4`, `systematics_input_he3/he4`, `metadata_output_he3/he4`

## Tasks

Supported tasks:

`merge_trees | analyse_data | analyse_mc | signal | systematics | checkpoint | report | full_chain`

## Notes

- Uses **PyROOT only** (no uproot/pandas).
- `signal` still uses `src/` RooFit C++ components loaded by PyROOT.
- Weighted efficiency histograms follow the Python naming policy `Weff*`.

## Smoke validation

Data smoke:

```bash
python3 he3_cli.py --config tests/smoke/validate_data.toml
python3 scripts/validate_smoke.py \
  --reference "$NUCLEI_INPUT/data/LHC22/apass4/DataHistosgiovanni.root" \
  --candidate /tmp/he3pp_validation/DataHistosgiovanni_py.root
```

MC smoke (content-only + naming normalization):

```bash
python3 he3_cli.py --config tests/smoke/validate_mc.toml
python3 scripts/validate_smoke.py \
  --reference "$NUCLEI_INPUT/MC/LHC23j6b/MChistosgiovanni.root" \
  --candidate /tmp/he3pp_validation/MChistosgiovanni_py.root \
  --normalize-weff \
  --ignore-errors
```

One-command smoke run:

```bash
scripts/run_smoke_checks.sh
```
