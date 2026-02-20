from array import array
import copy
from dataclasses import dataclass
from pathlib import Path
from typing import Any

try:
    import tomllib  # Python 3.11+
except ModuleNotFoundError:  # pragma: no cover - runtime compatibility path
    import tomli as tomllib


DEFAULTS_TOML = Path(__file__).with_name("defaults.toml")
CENT_LENGTH = 1
# Weighted efficiency histograms are intentionally named "Weff*".
WEIGHTED_EFF_NAMING_POLICY = "prefix_W"


def _load_default_config_from_toml() -> dict[str, Any]:
    with open(DEFAULTS_TOML, "rb") as f:
        cfg = tomllib.load(f)
    if not isinstance(cfg, dict):
        raise ValueError(f"Invalid defaults TOML at {DEFAULTS_TOML}: top-level table is missing.")
    return cfg


DEFAULT_CONFIG = _load_default_config_from_toml()


def default_config_template() -> dict[str, Any]:
    return copy.deepcopy(DEFAULT_CONFIG)


def _required_table(table: dict[str, Any], key: str, context: str = "defaults") -> dict[str, Any]:
    value = table.get(key)
    if not isinstance(value, dict):
        raise ValueError(f"Missing or invalid [{key}] table in {context} TOML: {DEFAULTS_TOML}")
    return value


def _required_value(table: dict[str, Any], key: str, context: str) -> Any:
    if key not in table:
        raise ValueError(f"Missing required key '{context}.{key}' in {DEFAULTS_TOML}")
    return table[key]


def _deep_merge_dict(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    out = copy.deepcopy(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(out.get(key), dict):
            out[key] = _deep_merge_dict(out[key], value)
        else:
            out[key] = copy.deepcopy(value)
    return out


def merge_config(cfg: dict[str, Any] | None) -> dict[str, Any]:
    merged = default_config_template()
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
    data_filename_he4: str
    data_analysis_results: str
    mc_analysis_results: str
    mc_tree_filename: str
    mc_filename: str
    mc_filename_he4: str
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


def _resolve_particle_templates(particle_profiles: dict[str, dict[str, Any]]) -> dict[str, dict[str, Any]]:
    resolved: dict[str, dict[str, Any]] = {}

    def resolve(species: str, stack: list[str]) -> dict[str, Any]:
        if species in resolved:
            return copy.deepcopy(resolved[species])
        if species in stack:
            cycle = " -> ".join(stack + [species])
            raise ValueError(f"Circular particle template chain: {cycle}")
        if species not in particle_profiles:
            raise ValueError(f"Unknown particle template target '{species}'.")

        profile = copy.deepcopy(particle_profiles[species])
        template = profile.get("template")
        if template is None:
            resolved_profile = profile
        else:
            template_key = str(template)
            base = resolve(template_key, stack + [species])
            for key, value in profile.items():
                if key == "template":
                    continue
                base[key] = value
            resolved_profile = base

        resolved[species] = resolved_profile
        return copy.deepcopy(resolved_profile)

    for species_key in particle_profiles:
        resolve(species_key, [])
    return resolved


def _build_particle_configs(cfg: dict[str, Any]) -> dict[str, dict[str, Any]]:
    selections = _required_table(cfg, "selections", "config")
    sel_he3 = _required_table(selections, "he3", "selections")
    sel_he4 = _required_table(selections, "he4", "selections")

    particle_cfg = _required_table(cfg, "particle", "config")
    base_profiles = {
        str(key): copy.deepcopy(value)
        for key, value in particle_cfg.items()
        if isinstance(value, dict)
    }

    if "he3" not in base_profiles or "he4" not in base_profiles:
        raise ValueError(
            f"config must define [particle.he3] and [particle.he4]. Found: {', '.join(sorted(base_profiles))}."
        )

    # Keep known species defaults centralized in configuration-derived values.
    base_profiles["he3"]["tof_nsigma_cut"] = float(_required_value(sel_he3, "nsigma_tof", "selections.he3"))
    base_profiles["he3"]["trial_dca_sel"] = str(_required_value(sel_he3, "trial_dca", "selections.he3"))
    base_profiles["he3"]["base_sel"] = str(_required_value(sel_he3, "base_rec", "selections.he3"))
    base_profiles["he3"]["primary_sel"] = str(_required_value(sel_he3, "default_rec", "selections.he3"))
    base_profiles["he3"]["secondary_sel"] = str(_required_value(sel_he3, "secondary_rec", "selections.he3"))
    base_profiles["he3"]["mc_reco_base_sel"] = (
        f"{str(_required_value(sel_he3, 'base_rec', 'selections.he3'))}"
        f"{str(_required_value(sel_he3, 'mc_reco_append', 'selections.he3'))}"
    )
    base_profiles["he3"]["mc_reco_sel"] = str(_required_value(sel_he3, "default_rec", "selections.he3"))
    base_profiles["he3"]["mc_gen_sel"] = str(_required_value(sel_he3, "mc_gen", "selections.he3"))
    base_profiles["he3"]["mc_signal_tracking"] = ""

    base_profiles["he4"]["tof_nsigma_cut"] = float(_required_value(sel_he4, "nsigma_tof", "selections.he4"))
    base_profiles["he4"]["trial_dca_sel"] = ""
    base_profiles["he4"]["base_sel"] = str(_required_value(sel_he4, "base_rec", "selections.he4"))
    base_profiles["he4"]["primary_sel"] = str(_required_value(sel_he4, "primary_rec", "selections.he4"))
    base_profiles["he4"]["secondary_sel"] = str(_required_value(sel_he4, "secondary_rec", "selections.he4"))
    base_profiles["he4"]["mc_reco_base_sel"] = ""
    base_profiles["he4"]["mc_reco_sel"] = str(_required_value(sel_he4, "mc_reco", "selections.he4"))
    base_profiles["he4"]["mc_gen_sel"] = str(_required_value(sel_he4, "mc_gen", "selections.he4"))
    base_profiles["he4"]["mc_signal_tracking"] = str(_required_value(sel_he4, "mc_signal_tracking", "selections.he4"))

    resolved_profiles = _resolve_particle_templates(base_profiles)

    # Apply per-species selection overrides for all configured particles.
    for species_key, species_cfg in resolved_profiles.items():
        sel_species = selections.get(species_key)
        if not isinstance(sel_species, dict):
            continue

        if "base_rec" in sel_species:
            species_cfg["base_sel"] = str(sel_species["base_rec"])
        if "primary_rec" in sel_species or "default_rec" in sel_species:
            species_cfg["primary_sel"] = str(sel_species.get("primary_rec", sel_species.get("default_rec")))
        if "secondary_rec" in sel_species:
            species_cfg["secondary_sel"] = str(sel_species["secondary_rec"])
        if "trial_dca" in sel_species:
            species_cfg["trial_dca_sel"] = str(sel_species["trial_dca"])
        if "nsigma_tof" in sel_species:
            species_cfg["tof_nsigma_cut"] = float(sel_species["nsigma_tof"])
        if "mc_reco_base" in sel_species:
            species_cfg["mc_reco_base_sel"] = str(sel_species["mc_reco_base"])
        elif "mc_reco_append" in sel_species:
            species_cfg["mc_reco_base_sel"] = f"{species_cfg.get('base_sel', '')}{str(sel_species['mc_reco_append'])}"
        if "mc_reco" in sel_species:
            species_cfg["mc_reco_sel"] = str(sel_species["mc_reco"])
        if "mc_gen" in sel_species:
            species_cfg["mc_gen_sel"] = str(sel_species["mc_gen"])
        if "mc_signal_tracking" in sel_species:
            species_cfg["mc_signal_tracking"] = str(sel_species["mc_signal_tracking"])

    return resolved_profiles


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
        data_filename_he4=f"{base_variant_output_dir}DataHistosHe4.root",
        data_analysis_results=f"{base_input_dir}data/{period}/{reco_pass}/{data_analysis_results_basename}",
        mc_analysis_results=f"{base_input_dir}MC/{mc_production}/{mc_analysis_results_basename}",
        mc_tree_filename=f"{base_input_dir}MC/{mc_production}/{mc_tree_basename}",
        mc_filename=f"{base_variant_output_dir}MChistos.root",
        mc_filename_he4=f"{base_variant_output_dir}MChistosHe4.root",
        signal_output=f"{base_variant_output_dir}signal.root",
        systematics_output=f"{base_variant_output_dir}systematics.root",
    )


def current_runtime_config(cfg: dict[str, Any] | None = None) -> RuntimeConfig:
    merged = merge_config(cfg)

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
    particle_configs = _build_particle_configs(merged)

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
