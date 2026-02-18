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

Use O2 in the current shell and export ROOT runtime paths:

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

## Tasks

Supported tasks:

`merge_trees | analyse_data | analyse_mc | signal | systematics | checkpoint | full_chain`

## Notes

- Uses **PyROOT only** (no uproot/pandas).
- `signal` still uses `src/` RooFit C++ components loaded by PyROOT.
