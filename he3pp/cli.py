import argparse
import copy
from datetime import datetime, timezone
import json
import logging
import subprocess
import sys
import time
try:
    import tomllib  # Python 3.11+
except ModuleNotFoundError:  # pragma: no cover - runtime compatibility path
    import tomli as tomllib

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
            "data_tree_basename": s.DATA_TREE_BASENAME,
            "data_analysis_results_basename": s.DATA_ANALYSIS_RESULTS_BASENAME,
            "mc_tree_basename": s.MC_TREE_BASENAME,
            "mc_analysis_results_basename": s.MC_ANALYSIS_RESULTS_BASENAME,
            "pt_bins": s.PT_BINS,
            "cent_pt_limits": s.CENT_PT_LIMITS,
            "tpc_max_pt": s.TPC_MAX_PT,
            "tof_min_pt": s.TOF_MIN_PT,
            "pt_range": s.PT_RANGE,
            "tpc_function_names": s.TPC_FUNCTION_NAMES,
        },
        "selections": {
            "common": {
                "skim_template": s.SKIM_SELECTION_TEMPLATE,
            },
            "he3": {
                "base_rec": s.BASE_REC_SELECTIONS,
                "default_rec": s.DEFAULT_REC_SELECTIONS,
                "secondary_rec": s.SECONDARY_SELECTION,
                "trial_dca": s.HE3_TRIAL_DCA_SELECTION,
                "nsigma_tof": s.HE3_NSIGMA_TOF_CUT,
                "mc_reco_append": s.HE3_MC_RECO_APPEND,
                "mc_gen": s.HE3_MC_GEN_SELECTION,
            },
            "he4": {
                "base_rec": s.HE4_BASE_SELECTION,
                "primary_rec": s.HE4_PRIMARY_SELECTION,
                "secondary_rec": s.SECONDARY_SELECTION,
                "nsigma_tof": s.HE4_NSIGMA_TOF_CUT,
                "mc_pid": s.HE4_MC_PID_SELECTION,
                "mc_reco": s.HE4_MC_RECO_SELECTION,
                "mc_gen": s.HE4_MC_GEN_SELECTION,
                "mc_signal_tracking": s.HE4_MC_SIGNAL_TRACKING,
            },
        },
        "cuts": {
            "nsigmaDCAz": s.CUT_NAMES["nsigmaDCAz"],
            "fTPCnCls": s.CUT_NAMES["fTPCnCls"],
            "nITScls": s.CUT_NAMES["nITScls"],
            "nsigmaTPC": s.CUT_NAMES["nsigmaTPC"],
        },
        "run": {
            "task": "analyse_data",
            "species": ["he3"],
            "enable_mt": True,
            "nthreads": 0,
            "enable_trials": True,
            "skim": False,
            "draw": False,
            "log_level": "INFO",
        },
        "species": {},
        "report": {
            "sections": ["signal_tof", "signal_tpc", "tof_tpc_2d", "efficiency", "pt_resolution", "corrected_spectrum"],
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


def _parse_species(run_cfg: dict) -> list[str]:
    raw = run_cfg.get("species")
    if raw is None:
        raise ValueError("Missing run.species in config.")
    if isinstance(raw, str):
        values = [raw]
    else:
        values = list(raw)
    out: list[str] = []
    for item in values:
        s_item = str(item).lower()
        if s_item not in ("he3", "he4"):
            raise ValueError(f"Unsupported species '{item}'. Allowed: he3, he4.")
        if s_item not in out:
            out.append(s_item)
    if not out:
        raise ValueError("run.species must contain at least one element.")
    return out


def _anti_species_name(species: str) -> str:
    return f"anti{species}"


def _species_paths(cfg: dict, species: str) -> dict:
    node = cfg.get("species", {}).get(species, {})
    if not isinstance(node, dict):
        return {}
    paths_node = node.get("paths", node)
    return paths_node if isinstance(paths_node, dict) else {}


def _species_output_defaults(species: str) -> dict[str, str]:
    base = f"{s.BASE_VARIANT_OUTPUT_DIR}{species}/"
    return {
        "data_output": f"{base}DataHistos.root",
        "data_input": f"{base}DataHistos.root",
        "mc_output": f"{base}MChistos.root",
        "mc_input": f"{base}MChistos.root",
        "signal_output": f"{base}signal.root",
        "signal_input": f"{base}signal.root",
        "systematics_output": f"{base}systematics.root",
        "systematics_input": f"{base}systematics.root",
        "metadata_output": f"{base}run_metadata.json",
    }


def _resolve_species_path(
    cfg: dict,
    path_cfg: dict,
    species: str,
    key: str,
    *,
    legacy_key: str | None = None,
    generic_key: str | None = None,
    default: str | None = None,
) -> str:
    sp_paths = _species_paths(cfg, species)
    if key in sp_paths:
        return sp_paths[key]
    if legacy_key and legacy_key in path_cfg:
        return path_cfg[legacy_key]
    if generic_key and generic_key in path_cfg:
        return path_cfg[generic_key]
    if default is not None:
        return default
    defaults = _species_output_defaults(species)
    if key in defaults:
        return defaults[key]
    raise KeyError(f"Missing path for species={species}: key={key}")


def run(cfg: dict) -> None:
    run_cfg = cfg.get("run", {})
    report_cfg = cfg.get("report", {})
    path_cfg_for_log = cfg.get("paths", {})
    _setup_logging(str(run_cfg.get("log_level", "INFO")), str(path_cfg_for_log.get("log_file", "") or ""))
    species = _parse_species(run_cfg)
    if str(run_cfg.get("task", "analyse_data")).lower() in {"analyse_data", "analyse_mc", "full_chain"}:
        sel_cfg = cfg.get("selections", {})
        missing_sel_sections = [sp for sp in species if not isinstance(sel_cfg.get(sp), dict)]
        if missing_sel_sections:
            raise ValueError(
                "Missing species selection sections for requested run.species: "
                + ", ".join(missing_sel_sections)
                + ". Add [selections.he3] / [selections.he4] as needed."
            )
    LOGGER.info("Starting run task=%s species=%s", run_cfg.get("task", "analyse_data"), species)

    s.apply_runtime_overrides(cfg)
    from . import tasks  # lazy import so tasks bind runtime-overridden settings

    path_cfg = dict(default_config()["paths"])
    path_cfg.update(cfg.get("paths", {}))

    task = str(run_cfg.get("task", "analyse_data")).lower()
    enable_mt = bool(run_cfg.get("enable_mt", True))
    nthreads = int(run_cfg.get("nthreads", 0))
    draw = bool(run_cfg.get("draw", False))

    if enable_mt:
        if nthreads > 0:
            ROOT.EnableImplicitMT(nthreads)
        else:
            ROOT.EnableImplicitMT()

    t0 = time.time()
    if task == "analyse_data":
        input_file = path_cfg.get("data_tree", s.DATA_TREE_FILENAME)
        selected_outputs = {
            sp: _resolve_species_path(
                cfg,
                path_cfg,
                sp,
                "data_output",
                legacy_key=f"data_output_{sp}",
                generic_key="data_output",
                default=None,
            )
            for sp in species
        }
        tasks.analyse_data_multi(input_file, selected_outputs, bool(run_cfg.get("skim", False)), draw)
    elif task == "analyse_mc":
        input_file = path_cfg.get("mc_tree", s.MC_TREE_FILENAME)
        selected_outputs = {
            sp: _resolve_species_path(
                cfg,
                path_cfg,
                sp,
                "mc_output",
                legacy_key=f"mc_output_{sp}",
                generic_key="mc_output",
                default=None,
            )
            for sp in species
        }
        tasks.analyse_mc_multi(input_file, selected_outputs, bool(run_cfg.get("enable_trials", True)), draw)
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
        from .reporting import generate_dual_report_index, generate_report

        root_report_dir = path_cfg.get("report_dir", f"{s.BASE_VARIANT_OUTPUT_DIR}report")
        if len(species) == 1:
            sp = species[0]
            report_index = generate_report(
                report_dir=root_report_dir,
                signal_file_path=_resolve_species_path(cfg, path_cfg, sp, "signal_input", legacy_key=f"signal_input_{sp}", generic_key="signal_input", default=None),
                mc_file_path=_resolve_species_path(cfg, path_cfg, sp, "mc_input", legacy_key=f"mc_input_{sp}", generic_key="mc_input", default=None),
                systematics_file_path=_resolve_species_path(cfg, path_cfg, sp, "systematics_input", legacy_key=f"systematics_input_{sp}", generic_key="systematics_input", default=None),
                metadata_path=_resolve_species_path(cfg, path_cfg, sp, "metadata_output", legacy_key=f"metadata_output_{sp}", generic_key="metadata_output", default="run_metadata.json"),
                species=_anti_species_name(sp),
                sections=list(report_cfg.get("sections", [])),
                fit_alpha=float(report_cfg.get("fit_alpha", 0.05)),
                fit_tail=str(report_cfg.get("fit_tail", "single")),
                tpc_signal_model=str(report_cfg.get("tpc_signal_model", "ExpGaus")),
                data_file_path=_resolve_species_path(cfg, path_cfg, sp, "data_input", legacy_key=f"data_input_{sp}", generic_key="data_input", default=None),
            )
            LOGGER.info("Report generated: %s", report_index)
        else:
            entries = []
            for sp in species:
                sub_report_dir = f"{root_report_dir}/{sp}"
                report_index = generate_report(
                    report_dir=sub_report_dir,
                    signal_file_path=_resolve_species_path(cfg, path_cfg, sp, "signal_input", legacy_key=f"signal_input_{sp}", generic_key="signal_input", default=None),
                    mc_file_path=_resolve_species_path(cfg, path_cfg, sp, "mc_input", legacy_key=f"mc_input_{sp}", generic_key="mc_input", default=None),
                    systematics_file_path=_resolve_species_path(cfg, path_cfg, sp, "systematics_input", legacy_key=f"systematics_input_{sp}", generic_key="systematics_input", default=None),
                    metadata_path=_resolve_species_path(cfg, path_cfg, sp, "metadata_output", legacy_key=f"metadata_output_{sp}", generic_key="metadata_output", default="run_metadata.json"),
                    species=_anti_species_name(sp),
                    sections=list(report_cfg.get("sections", [])),
                    fit_alpha=float(report_cfg.get("fit_alpha", 0.05)),
                    fit_tail=str(report_cfg.get("fit_tail", "single")),
                    tpc_signal_model=str(report_cfg.get("tpc_signal_model", "ExpGaus")),
                    data_file_path=_resolve_species_path(cfg, path_cfg, sp, "data_input", legacy_key=f"data_input_{sp}", generic_key="data_input", default=None),
                )
                entries.append({"label": sp.upper(), "species": _anti_species_name(sp), "href": f"{sp}/index.html", "path": report_index})
                LOGGER.info("Subreport generated [%s]: %s", sp, report_index)
            dual_index = generate_dual_report_index(
                report_dir=root_report_dir,
                entries=entries,
                metadata_path=path_cfg.get("metadata_output", "run_metadata.json"),
            )
            LOGGER.info("Combined report index generated: %s", dual_index)
    elif task == "full_chain":
        data_outputs = {
            sp: _resolve_species_path(
                cfg,
                path_cfg,
                sp,
                "data_output",
                legacy_key=f"data_output_{sp}",
                generic_key="data_output",
                default=None,
            )
            for sp in species
        }
        mc_outputs = {
            sp: _resolve_species_path(
                cfg,
                path_cfg,
                sp,
                "mc_output",
                legacy_key=f"mc_output_{sp}",
                generic_key="mc_output",
                default=None,
            )
            for sp in species
        }
        tasks.analyse_data_multi(path_cfg.get("data_tree", s.DATA_TREE_FILENAME), data_outputs, bool(run_cfg.get("skim", False)), draw)
        tasks.analyse_mc_multi(path_cfg.get("mc_tree", s.MC_TREE_FILENAME), mc_outputs, bool(run_cfg.get("enable_trials", True)), draw)

        from .reporting import generate_dual_report_index, generate_report

        entries = []
        root_report_dir = path_cfg.get("report_dir", f"{s.BASE_VARIANT_OUTPUT_DIR}report")
        for sp in species:
            sig_out = _resolve_species_path(cfg, path_cfg, sp, "signal_output", legacy_key=f"signal_output_{sp}", generic_key="signal_output", default=None)
            syst_out = _resolve_species_path(cfg, path_cfg, sp, "systematics_output", legacy_key=f"systematics_output_{sp}", generic_key="systematics_output", default=None)
            tasks.signal(data_outputs[sp], sig_out)
            tasks.systematics(sig_out, mc_outputs[sp], path_cfg.get("data_analysis_results", s.DATA_ANALYSIS_RESULTS), syst_out)

            if path_cfg.get("checkpoint_output"):
                tasks.checkpoint(
                    syst_out,
                    path_cfg.get("data_analysis_results", s.DATA_ANALYSIS_RESULTS),
                    mc_outputs[sp],
                    path_cfg.get("mc_analysis_results", s.MC_ANALYSIS_RESULTS),
                    sig_out,
                    path_cfg.get("checkpoint_output"),
                )

            sub_report_dir = root_report_dir if len(species) == 1 else f"{root_report_dir}/{sp}"
            report_index = generate_report(
                report_dir=sub_report_dir,
                signal_file_path=sig_out,
                mc_file_path=mc_outputs[sp],
                systematics_file_path=syst_out,
                metadata_path=_resolve_species_path(cfg, path_cfg, sp, "metadata_output", legacy_key=f"metadata_output_{sp}", generic_key="metadata_output", default="run_metadata.json"),
                species=_anti_species_name(sp),
                sections=list(report_cfg.get("sections", [])),
                fit_alpha=float(report_cfg.get("fit_alpha", 0.05)),
                fit_tail=str(report_cfg.get("fit_tail", "single")),
                tpc_signal_model=str(report_cfg.get("tpc_signal_model", "ExpGaus")),
                data_file_path=data_outputs[sp],
            )
            entries.append({"label": sp.upper(), "species": _anti_species_name(sp), "href": f"{sp}/index.html", "path": report_index})
            LOGGER.info("Report generated [%s]: %s", sp, report_index)
        if len(species) > 1:
            dual_index = generate_dual_report_index(
                report_dir=root_report_dir,
                entries=entries,
                metadata_path=path_cfg.get("metadata_output", "run_metadata.json"),
            )
            LOGGER.info("Combined report index generated: %s", dual_index)
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
