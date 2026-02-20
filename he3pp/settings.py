from array import array
import copy
from dataclasses import dataclass
from pathlib import Path
from typing import Any

try:
    import tomllib  # Python 3.11+
except ModuleNotFoundError:  # pragma: no cover - runtime compatibility path
    import tomli as tomllib


DEFAULTS_DIR = Path(__file__).parent
DEFAULTS_BY_SPECIES = {
    "he3": DEFAULTS_DIR / "defaults_he3.toml",
    "he4": DEFAULTS_DIR / "defaults_he4.toml",
}
CENT_LENGTH = 1
# Weighted efficiency histograms are intentionally named "Weff*".
WEIGHTED_EFF_NAMING_POLICY = "prefix_W"


def _load_toml(path: Path) -> dict[str, Any]:
    with open(path, "rb") as f:
        cfg = tomllib.load(f)
    if not isinstance(cfg, dict):
        raise ValueError(f"Invalid defaults TOML at {path}: top-level table is missing.")
    return cfg


_DEFAULT_CONFIG_CACHE: dict[str, dict[str, Any]] = {}


def _normalize_species_name(species: str | None) -> str:
    key = str(species or "he3").strip().lower()
    if key not in DEFAULTS_BY_SPECIES:
        raise ValueError(
            f"Unsupported species '{key}'. Supported defaults: {', '.join(sorted(DEFAULTS_BY_SPECIES))}."
        )
    return key


def _species_from_run_cfg(run_cfg: dict[str, Any] | None) -> str | None:
    if not isinstance(run_cfg, dict):
        return None
    raw = run_cfg.get("species")
    if raw is None:
        return None
    if isinstance(raw, str):
        return str(raw).strip().lower()
    values = [str(v).strip().lower() for v in list(raw)]
    if len(values) != 1:
        raise ValueError("Single-particle mode: run.species must contain exactly one value.")
    return values[0]


def _species_hint_from_cfg(cfg: dict[str, Any] | None) -> str | None:
    if not isinstance(cfg, dict):
        return None
    return _species_from_run_cfg(cfg.get("run"))


def default_config_template(species: str | None = None) -> dict[str, Any]:
    key = _normalize_species_name(species)
    if key not in _DEFAULT_CONFIG_CACHE:
        path = DEFAULTS_BY_SPECIES[key]
        _DEFAULT_CONFIG_CACHE[key] = _load_toml(path)
    return copy.deepcopy(_DEFAULT_CONFIG_CACHE[key])


def _required_table(table: dict[str, Any], key: str, context: str = "defaults") -> dict[str, Any]:
    value = table.get(key)
    if not isinstance(value, dict):
        raise ValueError(f"Missing or invalid [{key}] table in {context} config")
    return value


def _required_value(table: dict[str, Any], key: str, context: str) -> Any:
    if key not in table:
        raise ValueError(f"Missing required key '{context}.{key}'")
    return table[key]


def _deep_merge_dict(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    out = copy.deepcopy(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(out.get(key), dict):
            out[key] = _deep_merge_dict(out[key], value)
        else:
            out[key] = copy.deepcopy(value)
    return out


def merge_config(cfg: dict[str, Any] | None, species: str | None = None) -> dict[str, Any]:
    hint = species or _species_hint_from_cfg(cfg) or "he3"
    merged = default_config_template(hint)
    if not isinstance(cfg, dict):
        return merged
    return _deep_merge_dict(merged, cfg)


LETTER = ["M", "A"]


@dataclass(frozen=True)
class RuntimePaths:
    base_output_dir: str
    base_variant_output_dir: str
    data_tree_filename: str
    data_filename: str
    data_analysis_results: str
    mc_analysis_results: str
    mc_tree_filename: str
    mc_filename: str
    signal_output: str
    systematics_output: str

    def species_stage_paths(self, species: str) -> dict[str, str]:
        base = f"{self.base_variant_output_dir}{species}/"
        return {
            "data_output": f"{base}DataHistos.root",
            "data_input": f"{base}DataHistos.root",
            "mc_output": f"{base}MChistos.root",
            "mc_input": f"{base}MChistos.root",
            "signal_output": f"{base}signal.root",
            "signal_input": f"{base}signal.root",
            "systematics_output": f"{base}systematics.root",
            "systematics_input": f"{base}systematics.root",
            "checkpoint_output": f"{base}checkpoint.root",
            "metadata_output": f"{base}run_metadata.json",
            "report_dir": f"{base}report",
        }


@dataclass(frozen=True)
class RuntimeConfig:
    variant: str
    letter: list[str]
    filter_list_name: str
    n_pt_bins: int
    pt_bins: list[float]
    pt_bin_array: array
    cent_length: int
    cent_pt_limits: list[float]
    tpc_max_pt: float
    tof_min_pt: float
    pt_range: list[float]
    tpc_function_names: list[str]
    cut_names: dict[str, list[float]]
    skim_selection_template: str
    particle_configs: dict[str, dict[str, Any]]
    paths: RuntimePaths

    def get_particle_config(self, species: str) -> dict[str, Any]:
        if species not in self.particle_configs:
            raise ValueError(
                f"Unsupported species '{species}'. Available: {', '.join(sorted(self.particle_configs))}."
            )
        return copy.deepcopy(self.particle_configs[species])


def _resolve_particle_profile(species: str, particle_table: dict[str, Any], stack: list[str]) -> dict[str, Any]:
    if species in stack:
        cycle = " -> ".join(stack + [species])
        raise ValueError(f"Circular particle template chain: {cycle}")

    profile = particle_table.get(species)
    if not isinstance(profile, dict):
        raise ValueError(f"Missing [particle.{species}] table in config.")

    own = copy.deepcopy(profile)
    template = own.pop("template", None)
    if template is None:
        return own

    template_key = str(template).strip().lower()
    base = _resolve_particle_profile(template_key, particle_table, stack + [species])
    base.update(own)
    return base


def _build_particle_configs(cfg: dict[str, Any], species: str) -> dict[str, dict[str, Any]]:
    selections = _required_table(cfg, "selections", "config")
    sel_species = _required_table(selections, species, "selections")
    particle_table = _required_table(cfg, "particle", "config")

    profile = _resolve_particle_profile(species, particle_table, [])

    base_sel = str(_required_value(sel_species, "base_rec", f"selections.{species}"))
    primary_raw = sel_species.get("primary_rec", sel_species.get("default_rec"))
    if primary_raw is None:
        raise ValueError(
            f"Missing required key 'selections.{species}.primary_rec' (or 'default_rec')."
        )
    primary_sel = str(primary_raw)
    secondary_sel = str(_required_value(sel_species, "secondary_rec", f"selections.{species}"))

    profile["base_sel"] = base_sel
    profile["primary_sel"] = primary_sel
    profile["secondary_sel"] = secondary_sel
    profile["trial_dca_sel"] = str(sel_species.get("trial_dca", profile.get("trial_dca_sel", "")))
    profile["tof_nsigma_cut"] = float(_required_value(sel_species, "nsigma_tof", f"selections.{species}"))

    if "mc_reco_base" in sel_species:
        profile["mc_reco_base_sel"] = str(sel_species["mc_reco_base"])
    elif "mc_reco_append" in sel_species:
        profile["mc_reco_base_sel"] = f"{base_sel}{str(sel_species['mc_reco_append'])}"
    else:
        profile["mc_reco_base_sel"] = str(profile.get("mc_reco_base_sel", ""))

    if "mc_reco" in sel_species:
        profile["mc_reco_sel"] = str(sel_species["mc_reco"])
    else:
        profile["mc_reco_sel"] = str(profile.get("mc_reco_sel", primary_sel))

    if "mc_gen" in sel_species:
        profile["mc_gen_sel"] = str(sel_species["mc_gen"])
    elif "mc_gen_sel" not in profile:
        raise ValueError(f"Missing required key 'selections.{species}.mc_gen' or particle profile mc_gen_sel.")

    profile["mc_signal_tracking"] = str(sel_species.get("mc_signal_tracking", profile.get("mc_signal_tracking", "")))

    return {species: profile}


def _build_runtime_paths(common: dict[str, Any]) -> RuntimePaths:
    period = str(_required_value(common, "period", "common"))
    reco_pass = str(_required_value(common, "reco_pass", "common"))
    mc_production = str(_required_value(common, "mc_production", "common"))
    variant = str(_required_value(common, "variant", "common"))
    base_input_dir = str(_required_value(common, "base_input_dir", "common"))
    base_output_root = str(_required_value(common, "base_output_root", "common"))
    data_tree_basename = str(_required_value(common, "data_tree_basename", "common"))
    data_analysis_results_basename = str(_required_value(common, "data_analysis_results_basename", "common"))
    mc_tree_basename = str(_required_value(common, "mc_tree_basename", "common"))
    mc_analysis_results_basename = str(_required_value(common, "mc_analysis_results_basename", "common"))

    base_output_dir = f"{base_output_root}{period}/{reco_pass}/"
    base_variant_output_dir = f"{base_output_dir}{variant}/"

    return RuntimePaths(
        base_output_dir=base_output_dir,
        base_variant_output_dir=base_variant_output_dir,
        data_tree_filename=f"{base_input_dir}data/{period}/{reco_pass}/{data_tree_basename}",
        data_filename=f"{base_variant_output_dir}DataHistos.root",
        data_analysis_results=f"{base_input_dir}data/{period}/{reco_pass}/{data_analysis_results_basename}",
        mc_analysis_results=f"{base_input_dir}MC/{mc_production}/{mc_analysis_results_basename}",
        mc_tree_filename=f"{base_input_dir}MC/{mc_production}/{mc_tree_basename}",
        mc_filename=f"{base_variant_output_dir}MChistos.root",
        signal_output=f"{base_variant_output_dir}signal.root",
        systematics_output=f"{base_variant_output_dir}systematics.root",
    )


def current_runtime_config(cfg: dict[str, Any] | None = None) -> RuntimeConfig:
    merged = merge_config(cfg)

    run_cfg = _required_table(merged, "run", "config")
    species = _species_from_run_cfg(run_cfg)
    if species is None:
        raise ValueError("Missing run.species in config.")

    common = _required_table(merged, "common", "config")
    selections = _required_table(merged, "selections", "config")
    sel_common = _required_table(selections, "common", "selections")
    cuts = _required_table(merged, "cuts", "config")

    pt_bins = [float(v) for v in list(_required_value(common, "pt_bins", "common"))]
    if len(pt_bins) < 2:
        raise ValueError("common.pt_bins must contain at least 2 edges.")

    cent_pt_limits = [float(v) for v in list(_required_value(common, "cent_pt_limits", "common"))]
    if len(cent_pt_limits) != CENT_LENGTH:
        raise ValueError(f"common.cent_pt_limits must have length {CENT_LENGTH}.")

    pt_range = [float(v) for v in list(_required_value(common, "pt_range", "common"))]
    if len(pt_range) != 2:
        raise ValueError("common.pt_range must have 2 values: [min, max].")

    tpc_function_names = [str(v) for v in list(_required_value(common, "tpc_function_names", "common"))]
    if not tpc_function_names:
        raise ValueError("common.tpc_function_names must contain at least one entry.")

    cut_names = {
        "nsigmaDCAz": [float(v) for v in list(_required_value(cuts, "nsigmaDCAz", "cuts"))],
        "fTPCnCls": [float(v) for v in list(_required_value(cuts, "fTPCnCls", "cuts"))],
        "nITScls": [float(v) for v in list(_required_value(cuts, "nITScls", "cuts"))],
        "nsigmaTPC": [float(v) for v in list(_required_value(cuts, "nsigmaTPC", "cuts"))],
    }

    paths = _build_runtime_paths(common)
    particle_configs = _build_particle_configs(merged, species)

    return RuntimeConfig(
        variant=str(_required_value(common, "variant", "common")),
        letter=copy.deepcopy(LETTER),
        filter_list_name=str(_required_value(common, "filter_list_name", "common")),
        n_pt_bins=len(pt_bins) - 1,
        pt_bins=pt_bins,
        pt_bin_array=array("d", pt_bins),
        cent_length=CENT_LENGTH,
        cent_pt_limits=cent_pt_limits,
        tpc_max_pt=float(_required_value(common, "tpc_max_pt", "common")),
        tof_min_pt=float(_required_value(common, "tof_min_pt", "common")),
        pt_range=pt_range,
        tpc_function_names=tpc_function_names,
        cut_names=cut_names,
        skim_selection_template=str(_required_value(sel_common, "skim_template", "selections.common")),
        particle_configs=particle_configs,
        paths=paths,
    )


def get_particle_config(species: str, cfg: dict[str, Any] | None = None) -> dict[str, Any]:
    return current_runtime_config(cfg).get_particle_config(species)
