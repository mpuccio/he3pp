import argparse
import sys
import tomllib

import ROOT

from . import settings as s


def default_config() -> dict:
    return {
        "common": {
            "period": s.PERIOD,
            "reco_pass": s.RECO_PASS,
            "mc_production": s.MC_PRODUCTION,
            "variant": s.VARIANT,
            "base_input_dir": s.BASE_INPUT_DIR,
            "base_output_root": s.BASE_OUTPUT_ROOT,
            "filter_list_name": s.FILTER_LIST_NAME,
            "pt_bins": s.PT_BINS,
            "cent_pt_limits": s.CENT_PT_LIMITS,
            "tpc_max_pt": s.TPC_MAX_PT,
            "tof_min_pt": s.TOF_MIN_PT,
            "pt_range": s.PT_RANGE,
            "tpc_function_names": s.TPC_FUNCTION_NAMES,
        },
        "selections": {
            "base_rec": s.BASE_REC_SELECTIONS,
            "default_rec": s.DEFAULT_REC_SELECTIONS,
            "he4_base_rec": s.HE4_BASE_SELECTION,
            "he4_primary_rec": s.HE4_PRIMARY_SELECTION,
            "secondary_rec": s.SECONDARY_SELECTION,
            "he3_trial_dca": s.HE3_TRIAL_DCA_SELECTION,
            "he3_nsigma_tof": s.HE3_NSIGMA_TOF_CUT,
            "he4_nsigma_tof": s.HE4_NSIGMA_TOF_CUT,
            "skim_template": s.SKIM_SELECTION_TEMPLATE,
            "he3_mc_reco_append": s.HE3_MC_RECO_APPEND,
            "he3_mc_gen": s.HE3_MC_GEN_SELECTION,
            "he4_mc_pid": s.HE4_MC_PID_SELECTION,
            "he4_mc_reco": s.HE4_MC_RECO_SELECTION,
            "he4_mc_gen": s.HE4_MC_GEN_SELECTION,
            "he4_mc_signal_tracking": s.HE4_MC_SIGNAL_TRACKING,
        },
        "cuts": {
            "nsigmaDCAz": s.CUT_NAMES["nsigmaDCAz"],
            "fTPCnCls": s.CUT_NAMES["fTPCnCls"],
            "nITScls": s.CUT_NAMES["nITScls"],
            "nsigmaTPC": s.CUT_NAMES["nsigmaTPC"],
        },
        "run": {
            "task": "analyse_data",
            "particle": "he3",
            "enable_mt": True,
            "nthreads": 0,
            "enable_trials": True,
            "skim": False,
            "draw": False,
            "is_mc": True,
        },
        "paths": {
            "input": s.DATA_TREE_FILENAME,
            "output": s.DATA_FILENAME,
            "data_tree": s.DATA_TREE_FILENAME,
            "data_output": s.DATA_FILENAME,
            "mc_tree": s.MC_TREE_FILENAME,
            "mc_output": s.MC_FILENAME,
            "signal_output": s.SIGNAL_OUTPUT,
            "systematics_output": s.SYSTEMATICS_OUTPUT,
            "data_analysis_results": s.DATA_ANALYSIS_RESULTS,
            "mc_analysis_results": s.MC_ANALYSIS_RESULTS,
        },
    }


def run(cfg: dict) -> None:
    s.apply_runtime_overrides(cfg)
    from . import tasks  # lazy import so tasks bind runtime-overridden settings

    run_cfg = cfg.get("run", {})
    path_cfg = dict(default_config()["paths"])
    path_cfg.update(cfg.get("paths", {}))

    task = str(run_cfg.get("task", "analyse_data")).lower()
    particle = str(run_cfg.get("particle", "he3")).lower()
    enable_mt = bool(run_cfg.get("enable_mt", True))
    nthreads = int(run_cfg.get("nthreads", 0))
    draw = bool(run_cfg.get("draw", False))

    if enable_mt:
        if nthreads > 0:
            ROOT.EnableImplicitMT(nthreads)
        else:
            ROOT.EnableImplicitMT()

    if task == "merge_trees":
        tasks.merge_trees(path_cfg.get("input", s.MC_TREE_FILENAME), path_cfg.get("output", "MergedAO2D.root"), bool(run_cfg.get("is_mc", True)))
    elif task == "analyse_data":
        input_file = path_cfg.get("input", s.DATA_TREE_FILENAME)
        output_file = path_cfg.get("output", s.DATA_FILENAME_HE4 if particle == "he4" else s.DATA_FILENAME)
        tasks.analyse_data(input_file, output_file, particle, bool(run_cfg.get("skim", False)), draw)
    elif task == "analyse_mc":
        input_file = path_cfg.get("input", s.MC_TREE_FILENAME)
        output_file = path_cfg.get("output", s.MC_FILENAME_HE4 if particle == "he4" else s.MC_FILENAME)
        tasks.analyse_mc(input_file, output_file, particle, bool(run_cfg.get("enable_trials", True)), draw)
    elif task == "signal":
        tasks.signal(path_cfg.get("input", s.DATA_FILENAME), path_cfg.get("output", s.SIGNAL_OUTPUT))
    elif task == "systematics":
        tasks.systematics(
            path_cfg.get("signal_input", s.SIGNAL_OUTPUT),
            path_cfg.get("mc_input", s.MC_FILENAME),
            path_cfg.get("data_analysis_results", s.DATA_ANALYSIS_RESULTS),
            path_cfg.get("output", s.SYSTEMATICS_OUTPUT),
        )
    elif task == "checkpoint":
        tasks.checkpoint(
            path_cfg.get("systematics_input", s.SYSTEMATICS_OUTPUT),
            path_cfg.get("data_analysis_results", s.DATA_ANALYSIS_RESULTS),
            path_cfg.get("mc_input", s.MC_FILENAME),
            path_cfg.get("mc_analysis_results", s.MC_ANALYSIS_RESULTS),
            path_cfg.get("signal_input", s.SIGNAL_OUTPUT),
            path_cfg.get("output"),
        )
    elif task == "full_chain":
        tasks.analyse_data(path_cfg.get("data_tree", s.DATA_TREE_FILENAME), path_cfg.get("data_output", s.DATA_FILENAME), particle, bool(run_cfg.get("skim", False)), draw)
        tasks.analyse_mc(path_cfg.get("mc_tree", s.MC_TREE_FILENAME), path_cfg.get("mc_output", s.MC_FILENAME), particle, bool(run_cfg.get("enable_trials", True)), draw)
        tasks.signal(path_cfg.get("data_output", s.DATA_FILENAME), path_cfg.get("signal_output", s.SIGNAL_OUTPUT))
        tasks.systematics(
            path_cfg.get("signal_output", s.SIGNAL_OUTPUT),
            path_cfg.get("mc_output", s.MC_FILENAME),
            path_cfg.get("data_analysis_results", s.DATA_ANALYSIS_RESULTS),
            path_cfg.get("systematics_output", s.SYSTEMATICS_OUTPUT),
        )
        tasks.checkpoint(
            path_cfg.get("systematics_output", s.SYSTEMATICS_OUTPUT),
            path_cfg.get("data_analysis_results", s.DATA_ANALYSIS_RESULTS),
            path_cfg.get("mc_output", s.MC_FILENAME),
            path_cfg.get("mc_analysis_results", s.MC_ANALYSIS_RESULTS),
            path_cfg.get("signal_output", s.SIGNAL_OUTPUT),
            path_cfg.get("checkpoint_output"),
        )
    else:
        raise ValueError(f"Unsupported task: {task}")


def main() -> int:
    parser = argparse.ArgumentParser(description="Unified PyROOT CLI for anti-He3/anti-He4 analysis")
    parser.add_argument("--config", required=True, help="Path to TOML config")
    parser.add_argument("--dump-default-config", action="store_true", help="Print default config and exit")
    args = parser.parse_args()

    if args.dump_default_config:
        print(default_config())
        return 0

    with open(args.config, "rb") as f:
        cfg = tomllib.load(f)

    merged = default_config()
    merged.update(cfg)
    for section in ("common", "selections", "cuts", "run", "paths"):
        merged.setdefault(section, {})
        merged[section].update(cfg.get(section, {}))

    run(merged)
    return 0


if __name__ == "__main__":
    sys.exit(main())
