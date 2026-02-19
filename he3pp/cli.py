import argparse
import copy
from datetime import datetime, timezone
import json
import logging
import subprocess
import sys
import time
import tomllib

import ROOT

from . import settings as s


LOGGER = logging.getLogger("he3pp")


def _setup_logging(level: str = "INFO", log_file: str | None = None) -> None:
    lvl = getattr(logging, level.upper(), logging.INFO)
    logging.basicConfig(
        level=lvl,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=[logging.StreamHandler(sys.stdout)] + ([logging.FileHandler(log_file)] if log_file else []),
        force=True,
    )


def _git_revision() -> str:
    try:
        out = subprocess.run(["git", "rev-parse", "HEAD"], check=True, capture_output=True, text=True)
        return out.stdout.strip()
    except Exception:
        return "unknown"


def _write_metadata(path: str, payload: dict) -> None:
    from .root_io import ensure_parent, expand

    out = expand(path)
    ensure_parent(out)
    with open(out, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, sort_keys=True)


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
            "log_level": "INFO",
            "report_species": "antihe3",
        },
        "report": {
            "sections": ["signal_tof", "signal_tpc", "tof_tpc_2d", "raw_spectrum", "efficiency", "pt_resolution", "corrected_spectrum"],
            "fit_n_parameters": 6,
            "fit_alpha": 0.05,
            "fit_tail": "single",
            "tpc_signal_model": "ExpGaus",
        },
        "paths": {
            "input": s.DATA_TREE_FILENAME,
            "output": s.DATA_FILENAME,
            "data_tree": s.DATA_TREE_FILENAME,
            "data_output": s.DATA_FILENAME,
            "data_input": s.DATA_FILENAME,
            "mc_tree": s.MC_TREE_FILENAME,
            "mc_output": s.MC_FILENAME,
            "signal_output": s.SIGNAL_OUTPUT,
            "systematics_output": s.SYSTEMATICS_OUTPUT,
            "data_analysis_results": s.DATA_ANALYSIS_RESULTS,
            "mc_analysis_results": s.MC_ANALYSIS_RESULTS,
            "metadata_output": f"{s.BASE_VARIANT_OUTPUT_DIR}run_metadata.json",
            "log_file": "",
            "report_dir": f"{s.BASE_VARIANT_OUTPUT_DIR}report",
        },
    }


def run(cfg: dict) -> None:
    run_cfg = cfg.get("run", {})
    report_cfg = cfg.get("report", {})
    path_cfg_for_log = cfg.get("paths", {})
    _setup_logging(str(run_cfg.get("log_level", "INFO")), str(path_cfg_for_log.get("log_file", "") or ""))
    LOGGER.info("Starting run task=%s particle=%s", run_cfg.get("task", "analyse_data"), run_cfg.get("particle", "he3"))

    s.apply_runtime_overrides(cfg)
    from . import tasks  # lazy import so tasks bind runtime-overridden settings

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

    t0 = time.time()
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
    elif task == "report":
        from .reporting import generate_report

        report_index = generate_report(
            report_dir=path_cfg.get("report_dir", f"{s.BASE_VARIANT_OUTPUT_DIR}report"),
            signal_file_path=path_cfg.get("signal_input", path_cfg.get("signal_output", s.SIGNAL_OUTPUT)),
            mc_file_path=path_cfg.get("mc_input", path_cfg.get("mc_output", s.MC_FILENAME)),
            systematics_file_path=path_cfg.get("systematics_input", path_cfg.get("systematics_output", s.SYSTEMATICS_OUTPUT)),
            metadata_path=path_cfg.get("metadata_output", "run_metadata.json"),
            species=str(run_cfg.get("report_species", "antihe3")),
            sections=list(report_cfg.get("sections", [])),
            fit_n_parameters=int(report_cfg.get("fit_n_parameters", 6)),
            fit_alpha=float(report_cfg.get("fit_alpha", 0.05)),
            fit_tail=str(report_cfg.get("fit_tail", "single")),
            tpc_signal_model=str(report_cfg.get("tpc_signal_model", "ExpGaus")),
            data_file_path=path_cfg.get("data_input", path_cfg.get("data_output", s.DATA_FILENAME)),
        )
        LOGGER.info("Report generated: %s", report_index)
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
        from .reporting import generate_report

        report_index = generate_report(
            report_dir=path_cfg.get("report_dir", f"{s.BASE_VARIANT_OUTPUT_DIR}report"),
            signal_file_path=path_cfg.get("signal_output", s.SIGNAL_OUTPUT),
            mc_file_path=path_cfg.get("mc_output", s.MC_FILENAME),
            systematics_file_path=path_cfg.get("systematics_output", s.SYSTEMATICS_OUTPUT),
            metadata_path=path_cfg.get("metadata_output", "run_metadata.json"),
            species=str(run_cfg.get("report_species", "antihe3")),
            sections=list(report_cfg.get("sections", [])),
            fit_n_parameters=int(report_cfg.get("fit_n_parameters", 6)),
            fit_alpha=float(report_cfg.get("fit_alpha", 0.05)),
            fit_tail=str(report_cfg.get("fit_tail", "single")),
            tpc_signal_model=str(report_cfg.get("tpc_signal_model", "ExpGaus")),
            data_file_path=path_cfg.get("data_output", s.DATA_FILENAME),
        )
        LOGGER.info("Report generated: %s", report_index)
    else:
        raise ValueError(f"Unsupported task: {task}")
    LOGGER.info("Finished run task=%s elapsed_sec=%.2f", task, time.time() - t0)


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
    for section in ("common", "selections", "cuts", "run", "report", "paths"):
        merged.setdefault(section, {})
        merged[section].update(cfg.get(section, {}))

    started = datetime.now(timezone.utc)
    status = "success"
    error = ""
    try:
        run(merged)
    except Exception as exc:
        status = "failed"
        error = str(exc)
        raise
    finally:
        ended = datetime.now(timezone.utc)
        metadata = {
            "status": status,
            "error": error,
            "started_utc": started.isoformat(),
            "ended_utc": ended.isoformat(),
            "duration_sec": (ended - started).total_seconds(),
            "git_revision": _git_revision(),
            "config": copy.deepcopy(merged),
        }
        try:
            _write_metadata(merged.get("paths", {}).get("metadata_output", "run_metadata.json"), metadata)
        except Exception as meta_exc:
            LOGGER.error("Failed to write metadata: %s", meta_exc)
    return 0


if __name__ == "__main__":
    sys.exit(main())
