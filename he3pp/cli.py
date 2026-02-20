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
    return s.default_config_template()


def _parse_species(run_cfg: dict, allowed_species: set[str]) -> str:
    raw = run_cfg.get("species")
    if raw is None:
        raise ValueError("Missing run.species in config.")
    if isinstance(raw, str):
        value = raw.strip().lower()
    else:
        values = [str(v).strip().lower() for v in list(raw)]
        if len(values) != 1:
            raise ValueError("Single-particle mode: run.species must contain exactly one value.")
        value = values[0]
    if value not in allowed_species:
        raise ValueError(
            f"Unsupported run.species '{value}'. Allowed: {', '.join(sorted(allowed_species))}."
        )
    return value


def run(cfg: dict, user_cfg: dict | None = None) -> None:
    run_cfg = cfg.get("run", {})
    report_cfg = cfg.get("report", {})
    path_cfg_for_log = cfg.get("paths", {})
    _setup_logging(str(run_cfg.get("log_level", "INFO")), str(path_cfg_for_log.get("log_file", "") or ""))
    particle_cfg = cfg.get("particle", {})
    allowed_species = {
        str(key).strip().lower()
        for key, value in (particle_cfg.items() if isinstance(particle_cfg, dict) else [])
        if isinstance(value, dict)
    }
    if not allowed_species:
        allowed_species = {"he3", "he4"}
    species = _parse_species(run_cfg, allowed_species)
    if str(run_cfg.get("task", "analyse_data")).lower() in {"analyse_data", "analyse_mc", "full_chain"}:
        sel_cfg = cfg.get("selections", {})
        if not isinstance(sel_cfg.get(species), dict):
            raise ValueError(
                f"Missing selection section for run.species='{species}'. "
                f"Add [selections.{species}] in config."
            )
    LOGGER.info("Starting run task=%s species=%s", run_cfg.get("task", "analyse_data"), species)

    runtime_cfg = s.current_runtime_config(cfg)
    from . import tasks

    user_paths = (user_cfg or {}).get("paths", {})
    if not isinstance(user_paths, dict):
        user_paths = {}

    path_cfg = dict(cfg.get("paths", {}))
    species_defaults = runtime_cfg.paths.species_stage_paths(species)
    for key, value in species_defaults.items():
        if key not in user_paths:
            path_cfg[key] = value

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
        input_file = path_cfg.get("data_tree", runtime_cfg.paths.data_tree_filename)
        tasks.analyse_data(
            input_file,
            path_cfg.get("data_output", runtime_cfg.paths.data_filename),
            species,
            bool(run_cfg.get("skim", False)),
            draw,
            runtime_config=runtime_cfg,
        )
    elif task == "analyse_mc":
        input_file = path_cfg.get("mc_tree", runtime_cfg.paths.mc_tree_filename)
        tasks.analyse_mc(
            input_file,
            path_cfg.get("mc_output", runtime_cfg.paths.mc_filename),
            species,
            bool(run_cfg.get("enable_trials", True)),
            draw,
            runtime_config=runtime_cfg,
        )
    elif task == "signal":
        tasks.signal(
            path_cfg.get("data_input", runtime_cfg.paths.data_filename),
            path_cfg.get("signal_output", runtime_cfg.paths.signal_output),
            species,
            runtime_config=runtime_cfg,
        )
    elif task == "systematics":
        tasks.systematics(
            path_cfg.get("signal_input", runtime_cfg.paths.signal_output),
            path_cfg.get("mc_input", runtime_cfg.paths.mc_filename),
            path_cfg.get("data_analysis_results", runtime_cfg.paths.data_analysis_results),
            path_cfg.get("systematics_output", runtime_cfg.paths.systematics_output),
            species,
            runtime_config=runtime_cfg,
        )
    elif task == "checkpoint":
        tasks.checkpoint(
            path_cfg.get("systematics_input", runtime_cfg.paths.systematics_output),
            path_cfg.get("data_analysis_results", runtime_cfg.paths.data_analysis_results),
            path_cfg.get("mc_input", runtime_cfg.paths.mc_filename),
            path_cfg.get("mc_analysis_results", runtime_cfg.paths.mc_analysis_results),
            path_cfg.get("signal_input", runtime_cfg.paths.signal_output),
            path_cfg.get("checkpoint_output"),
            species,
            runtime_config=runtime_cfg,
        )
    elif task == "report":
        from .reporting import generate_report

        particle_profile = runtime_cfg.get_particle_config(species)
        report_index = generate_report(
            report_dir=path_cfg.get("report_dir", runtime_cfg.paths.species_stage_paths(species)["report_dir"]),
            signal_file_path=path_cfg.get("signal_input", runtime_cfg.paths.signal_output),
            mc_file_path=path_cfg.get("mc_input", runtime_cfg.paths.mc_filename),
            systematics_file_path=path_cfg.get("systematics_input", runtime_cfg.paths.systematics_output),
            metadata_path=path_cfg.get("metadata_output", "run_metadata.json"),
            species=str(particle_profile["anti_name"]),
            sections=list(report_cfg.get("sections", [])),
            fit_alpha=float(report_cfg.get("fit_alpha", 0.05)),
            fit_tail=str(report_cfg.get("fit_tail", "single")),
            tpc_signal_model=str(report_cfg.get("tpc_signal_model", "ExpGaus")),
            mc_hist_suffix=str(particle_profile.get("mc_hist_suffix", "He3")),
            data_file_path=path_cfg.get("data_input", runtime_cfg.paths.data_filename),
            runtime_config=runtime_cfg,
        )
        LOGGER.info("Report generated: %s", report_index)
    elif task == "full_chain":
        data_out = path_cfg.get("data_output", runtime_cfg.paths.data_filename)
        mc_out = path_cfg.get("mc_output", runtime_cfg.paths.mc_filename)
        sig_out = path_cfg.get("signal_output", runtime_cfg.paths.signal_output)
        syst_out = path_cfg.get("systematics_output", runtime_cfg.paths.systematics_output)

        tasks.analyse_data(
            path_cfg.get("data_tree", runtime_cfg.paths.data_tree_filename),
            data_out,
            species,
            bool(run_cfg.get("skim", False)),
            draw,
            runtime_config=runtime_cfg,
        )
        tasks.analyse_mc(
            path_cfg.get("mc_tree", runtime_cfg.paths.mc_tree_filename),
            mc_out,
            species,
            bool(run_cfg.get("enable_trials", True)),
            draw,
            runtime_config=runtime_cfg,
        )
        tasks.signal(data_out, sig_out, species, runtime_config=runtime_cfg)
        tasks.systematics(
            sig_out,
            mc_out,
            path_cfg.get("data_analysis_results", runtime_cfg.paths.data_analysis_results),
            syst_out,
            species,
            runtime_config=runtime_cfg,
        )

        if "checkpoint_output" in user_paths and path_cfg.get("checkpoint_output"):
            tasks.checkpoint(
                syst_out,
                path_cfg.get("data_analysis_results", runtime_cfg.paths.data_analysis_results),
                mc_out,
                path_cfg.get("mc_analysis_results", runtime_cfg.paths.mc_analysis_results),
                sig_out,
                path_cfg.get("checkpoint_output"),
                species,
                runtime_config=runtime_cfg,
            )

        from .reporting import generate_report

        particle_profile = runtime_cfg.get_particle_config(species)
        report_index = generate_report(
            report_dir=path_cfg.get("report_dir", runtime_cfg.paths.species_stage_paths(species)["report_dir"]),
            signal_file_path=sig_out,
            mc_file_path=mc_out,
            systematics_file_path=syst_out,
            metadata_path=path_cfg.get("metadata_output", "run_metadata.json"),
            species=str(particle_profile["anti_name"]),
            sections=list(report_cfg.get("sections", [])),
            fit_alpha=float(report_cfg.get("fit_alpha", 0.05)),
            fit_tail=str(report_cfg.get("fit_tail", "single")),
            tpc_signal_model=str(report_cfg.get("tpc_signal_model", "ExpGaus")),
            mc_hist_suffix=str(particle_profile.get("mc_hist_suffix", "He3")),
            data_file_path=data_out,
            runtime_config=runtime_cfg,
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

    merged = s.merge_config(cfg)

    started = datetime.now(timezone.utc)
    status = "success"
    error = ""
    try:
        run(merged, cfg)
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
            metadata_path = (cfg.get("paths", {}) if isinstance(cfg, dict) else {}).get("metadata_output")
            if not metadata_path:
                try:
                    runtime_cfg = s.current_runtime_config(merged)
                    allowed_species = set(runtime_cfg.particle_configs.keys()) or {"he3", "he4"}
                    species = _parse_species(merged.get("run", {}), allowed_species)
                    metadata_path = runtime_cfg.paths.species_stage_paths(species)["metadata_output"]
                except Exception:
                    metadata_path = merged.get("paths", {}).get("metadata_output", "run_metadata.json")
            _write_metadata(metadata_path, metadata)
        except Exception as meta_exc:
            LOGGER.error("Failed to write metadata: %s", meta_exc)
    return 0


if __name__ == "__main__":
    sys.exit(main())
