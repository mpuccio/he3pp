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


LETTER = ["M", "A"]
NAMES = ["he3", "antihe3"]
LABELS = ["^{3}He", "^{3}#bar{He}"]

_common_defaults = _required_table(DEFAULT_CONFIG, "common")
_sel_defaults = _required_table(DEFAULT_CONFIG, "selections")
_sel_common_defaults = _required_table(_sel_defaults, "common", "selections")
_sel_he3_defaults = _required_table(_sel_defaults, "he3", "selections")
_sel_he4_defaults = _required_table(_sel_defaults, "he4", "selections")
_cut_defaults = _required_table(DEFAULT_CONFIG, "cuts")
_particle_defaults = _required_table(DEFAULT_CONFIG, "particle")

MC_PRODUCTION = str(_required_value(_common_defaults, "mc_production", "common"))
RECO_PASS = str(_required_value(_common_defaults, "reco_pass", "common"))
PERIOD = str(_required_value(_common_defaults, "period", "common"))
VARIANT = str(_required_value(_common_defaults, "variant", "common"))
BASE_INPUT_DIR = str(_required_value(_common_defaults, "base_input_dir", "common"))
BASE_OUTPUT_ROOT = str(_required_value(_common_defaults, "base_output_root", "common"))
FILTER_LIST_NAME = str(_required_value(_common_defaults, "filter_list_name", "common"))
DATA_TREE_BASENAME = str(_required_value(_common_defaults, "data_tree_basename", "common"))
DATA_ANALYSIS_RESULTS_BASENAME = str(_required_value(_common_defaults, "data_analysis_results_basename", "common"))
MC_TREE_BASENAME = str(_required_value(_common_defaults, "mc_tree_basename", "common"))
MC_ANALYSIS_RESULTS_BASENAME = str(_required_value(_common_defaults, "mc_analysis_results_basename", "common"))

BASE_REC_SELECTIONS = str(_required_value(_sel_he3_defaults, "base_rec", "selections.he3"))
DEFAULT_REC_SELECTIONS = str(_required_value(_sel_he3_defaults, "default_rec", "selections.he3"))
HE4_BASE_SELECTION = str(_required_value(_sel_he4_defaults, "base_rec", "selections.he4"))
HE4_PRIMARY_SELECTION = str(_required_value(_sel_he4_defaults, "primary_rec", "selections.he4"))
SECONDARY_SELECTION = str(_required_value(_sel_he3_defaults, "secondary_rec", "selections.he3"))
HE3_TRIAL_DCA_SELECTION = str(_required_value(_sel_he3_defaults, "trial_dca", "selections.he3"))
HE3_NSIGMA_TOF_CUT = float(_required_value(_sel_he3_defaults, "nsigma_tof", "selections.he3"))
HE4_NSIGMA_TOF_CUT = float(_required_value(_sel_he4_defaults, "nsigma_tof", "selections.he4"))
SKIM_SELECTION_TEMPLATE = str(_required_value(_sel_common_defaults, "skim_template", "selections.common"))
HE3_MC_RECO_APPEND = str(_required_value(_sel_he3_defaults, "mc_reco_append", "selections.he3"))
HE3_MC_GEN_SELECTION = str(_required_value(_sel_he3_defaults, "mc_gen", "selections.he3"))
HE4_MC_PID_SELECTION = str(_required_value(_sel_he4_defaults, "mc_pid", "selections.he4"))
HE4_MC_RECO_SELECTION = str(_required_value(_sel_he4_defaults, "mc_reco", "selections.he4"))
HE4_MC_GEN_SELECTION = str(_required_value(_sel_he4_defaults, "mc_gen", "selections.he4"))
HE4_MC_SIGNAL_TRACKING = str(_required_value(_sel_he4_defaults, "mc_signal_tracking", "selections.he4"))

PT_BINS = [float(v) for v in _required_value(_common_defaults, "pt_bins", "common")]
if len(PT_BINS) < 2:
    raise ValueError(f"common.pt_bins in {DEFAULTS_TOML} must contain at least 2 edges.")
N_PT_BINS = len(PT_BINS) - 1
PT_BIN_ARRAY = array("d", PT_BINS)
CENT_LENGTH = 1
CENT_PT_LIMITS = [float(v) for v in _required_value(_common_defaults, "cent_pt_limits", "common")]
if len(CENT_PT_LIMITS) != CENT_LENGTH:
    raise ValueError(f"common.cent_pt_limits in {DEFAULTS_TOML} must have length {CENT_LENGTH}.")
TPC_MAX_PT = float(_required_value(_common_defaults, "tpc_max_pt", "common"))
TOF_MIN_PT = float(_required_value(_common_defaults, "tof_min_pt", "common"))
PT_RANGE = [float(v) for v in _required_value(_common_defaults, "pt_range", "common")]
TPC_FUNCTION_NAMES = [str(v) for v in _required_value(_common_defaults, "tpc_function_names", "common")]
if not TPC_FUNCTION_NAMES:
    raise ValueError(f"common.tpc_function_names in {DEFAULTS_TOML} must not be empty.")
# Weighted efficiency histograms are intentionally named "Weff*".
WEIGHTED_EFF_NAMING_POLICY = "prefix_W"

CUT_NAMES = {
    "nsigmaDCAz": [float(v) for v in _required_value(_cut_defaults, "nsigmaDCAz", "cuts")],
    "fTPCnCls": [float(v) for v in _required_value(_cut_defaults, "fTPCnCls", "cuts")],
    "nITScls": [float(v) for v in _required_value(_cut_defaults, "nITScls", "cuts")],
    "nsigmaTPC": [float(v) for v in _required_value(_cut_defaults, "nsigmaTPC", "cuts")],
}

BASE_OUTPUT_DIR = ""
BASE_VARIANT_OUTPUT_DIR = ""
DATA_TREE_FILENAME = ""
DATA_FILENAME = ""
DATA_FILENAME_HE4 = ""
DATA_ANALYSIS_RESULTS = ""
MC_ANALYSIS_RESULTS = ""
MC_TREE_FILENAME = ""
MC_FILENAME = ""
MC_FILENAME_HE4 = ""
SIGNAL_OUTPUT = ""
SYSTEMATICS_OUTPUT = ""

DEFAULT_PARTICLE_CONFIGS = {
    str(k): copy.deepcopy(v)
    for k, v in (_particle_defaults.items() if isinstance(_particle_defaults, dict) else [])
    if isinstance(v, dict)
}
if "he3" not in DEFAULT_PARTICLE_CONFIGS or "he4" not in DEFAULT_PARTICLE_CONFIGS:
    raise ValueError(
        f"defaults.toml must define [particle.he3] and [particle.he4]. Found: {', '.join(sorted(DEFAULT_PARTICLE_CONFIGS))}."
    )

PARTICLE_CONFIGS = copy.deepcopy(DEFAULT_PARTICLE_CONFIGS)


def _refresh_particle_configs() -> None:
    global PARTICLE_CONFIGS
    PARTICLE_CONFIGS = copy.deepcopy(DEFAULT_PARTICLE_CONFIGS)

    if "he3" in PARTICLE_CONFIGS:
        PARTICLE_CONFIGS["he3"]["tof_nsigma_cut"] = HE3_NSIGMA_TOF_CUT
        PARTICLE_CONFIGS["he3"]["trial_dca_sel"] = HE3_TRIAL_DCA_SELECTION
        PARTICLE_CONFIGS["he3"]["base_sel"] = BASE_REC_SELECTIONS
        PARTICLE_CONFIGS["he3"]["primary_sel"] = DEFAULT_REC_SELECTIONS
        PARTICLE_CONFIGS["he3"]["secondary_sel"] = SECONDARY_SELECTION
        PARTICLE_CONFIGS["he3"]["mc_reco_base_sel"] = BASE_REC_SELECTIONS + HE3_MC_RECO_APPEND
        PARTICLE_CONFIGS["he3"]["mc_reco_sel"] = DEFAULT_REC_SELECTIONS
        PARTICLE_CONFIGS["he3"]["mc_gen_sel"] = HE3_MC_GEN_SELECTION
        PARTICLE_CONFIGS["he3"]["mc_signal_tracking"] = ""

    if "he4" in PARTICLE_CONFIGS:
        PARTICLE_CONFIGS["he4"]["tof_nsigma_cut"] = HE4_NSIGMA_TOF_CUT
        PARTICLE_CONFIGS["he4"]["trial_dca_sel"] = ""
        PARTICLE_CONFIGS["he4"]["base_sel"] = HE4_BASE_SELECTION
        PARTICLE_CONFIGS["he4"]["primary_sel"] = HE4_PRIMARY_SELECTION
        PARTICLE_CONFIGS["he4"]["secondary_sel"] = SECONDARY_SELECTION
        PARTICLE_CONFIGS["he4"]["mc_reco_base_sel"] = ""
        PARTICLE_CONFIGS["he4"]["mc_reco_sel"] = HE4_MC_RECO_SELECTION
        PARTICLE_CONFIGS["he4"]["mc_gen_sel"] = HE4_MC_GEN_SELECTION
        PARTICLE_CONFIGS["he4"]["mc_signal_tracking"] = HE4_MC_SIGNAL_TRACKING


def recompute_derived_globals() -> None:
    global BASE_OUTPUT_DIR
    global BASE_VARIANT_OUTPUT_DIR
    global DATA_TREE_FILENAME, DATA_FILENAME, DATA_FILENAME_HE4, DATA_ANALYSIS_RESULTS
    global MC_ANALYSIS_RESULTS, MC_TREE_FILENAME, MC_FILENAME, MC_FILENAME_HE4
    global SIGNAL_OUTPUT, SYSTEMATICS_OUTPUT

    BASE_OUTPUT_DIR = f"{BASE_OUTPUT_ROOT}{PERIOD}/{RECO_PASS}/"
    BASE_VARIANT_OUTPUT_DIR = f"{BASE_OUTPUT_DIR}{VARIANT}/"
    DATA_TREE_FILENAME = f"{BASE_INPUT_DIR}data/{PERIOD}/{RECO_PASS}/{DATA_TREE_BASENAME}"
    DATA_FILENAME = f"{BASE_VARIANT_OUTPUT_DIR}DataHistos.root"
    DATA_FILENAME_HE4 = f"{BASE_VARIANT_OUTPUT_DIR}DataHistosHe4.root"
    DATA_ANALYSIS_RESULTS = f"{BASE_INPUT_DIR}data/{PERIOD}/{RECO_PASS}/{DATA_ANALYSIS_RESULTS_BASENAME}"
    MC_ANALYSIS_RESULTS = f"{BASE_INPUT_DIR}MC/{MC_PRODUCTION}/{MC_ANALYSIS_RESULTS_BASENAME}"
    MC_TREE_FILENAME = f"{BASE_INPUT_DIR}MC/{MC_PRODUCTION}/{MC_TREE_BASENAME}"
    MC_FILENAME = f"{BASE_VARIANT_OUTPUT_DIR}MChistos.root"
    MC_FILENAME_HE4 = f"{BASE_VARIANT_OUTPUT_DIR}MChistosHe4.root"
    SIGNAL_OUTPUT = f"{BASE_VARIANT_OUTPUT_DIR}signal.root"
    SYSTEMATICS_OUTPUT = f"{BASE_VARIANT_OUTPUT_DIR}systematics.root"


def apply_runtime_overrides(cfg: dict[str, Any]) -> None:
    global MC_PRODUCTION, RECO_PASS, PERIOD, VARIANT, BASE_INPUT_DIR, BASE_OUTPUT_ROOT, FILTER_LIST_NAME
    global DATA_TREE_BASENAME, DATA_ANALYSIS_RESULTS_BASENAME, MC_TREE_BASENAME, MC_ANALYSIS_RESULTS_BASENAME
    global BASE_REC_SELECTIONS, DEFAULT_REC_SELECTIONS, HE4_BASE_SELECTION, HE4_PRIMARY_SELECTION
    global SECONDARY_SELECTION, HE3_TRIAL_DCA_SELECTION, HE3_NSIGMA_TOF_CUT, HE4_NSIGMA_TOF_CUT
    global SKIM_SELECTION_TEMPLATE
    global HE3_MC_RECO_APPEND, HE3_MC_GEN_SELECTION, HE4_MC_PID_SELECTION, HE4_MC_RECO_SELECTION
    global HE4_MC_GEN_SELECTION, HE4_MC_SIGNAL_TRACKING
    global N_PT_BINS, PT_BINS, PT_BIN_ARRAY, CENT_PT_LIMITS, TPC_MAX_PT, TOF_MIN_PT, PT_RANGE
    global TPC_FUNCTION_NAMES, CUT_NAMES, PARTICLE_CONFIGS

    common = cfg.get("common", {})
    selections = cfg.get("selections", {})
    sel_common = selections.get("common", {}) if isinstance(selections, dict) else {}
    sel_he3 = selections.get("he3", {}) if isinstance(selections, dict) else {}
    sel_he4 = selections.get("he4", {}) if isinstance(selections, dict) else {}
    cuts = cfg.get("cuts", {})

    MC_PRODUCTION = str(common.get("mc_production", MC_PRODUCTION))
    RECO_PASS = str(common.get("reco_pass", RECO_PASS))
    PERIOD = str(common.get("period", PERIOD))
    VARIANT = str(common.get("variant", VARIANT))
    BASE_INPUT_DIR = str(common.get("base_input_dir", BASE_INPUT_DIR))
    BASE_OUTPUT_ROOT = str(common.get("base_output_root", BASE_OUTPUT_ROOT))
    FILTER_LIST_NAME = str(common.get("filter_list_name", FILTER_LIST_NAME))
    DATA_TREE_BASENAME = str(common.get("data_tree_basename", DATA_TREE_BASENAME))
    DATA_ANALYSIS_RESULTS_BASENAME = str(common.get("data_analysis_results_basename", DATA_ANALYSIS_RESULTS_BASENAME))
    MC_TREE_BASENAME = str(common.get("mc_tree_basename", MC_TREE_BASENAME))
    MC_ANALYSIS_RESULTS_BASENAME = str(common.get("mc_analysis_results_basename", MC_ANALYSIS_RESULTS_BASENAME))

    pt_bins = list(common.get("pt_bins", PT_BINS))
    if len(pt_bins) < 2:
        raise ValueError("common.pt_bins must contain at least 2 edges.")
    PT_BINS = [float(v) for v in pt_bins]
    N_PT_BINS = len(PT_BINS) - 1
    PT_BIN_ARRAY = array("d", PT_BINS)

    cent_pt_limits = [float(v) for v in common.get("cent_pt_limits", CENT_PT_LIMITS)]
    if len(cent_pt_limits) != CENT_LENGTH:
        raise ValueError(f"common.cent_pt_limits must have length {CENT_LENGTH}.")
    CENT_PT_LIMITS = cent_pt_limits

    TPC_MAX_PT = float(common.get("tpc_max_pt", TPC_MAX_PT))
    TOF_MIN_PT = float(common.get("tof_min_pt", TOF_MIN_PT))
    PT_RANGE = [float(v) for v in common.get("pt_range", PT_RANGE)]
    if len(PT_RANGE) != 2:
        raise ValueError("common.pt_range must have 2 values: [min, max].")

    tpc_names = list(common.get("tpc_function_names", TPC_FUNCTION_NAMES))
    if not tpc_names:
        raise ValueError("common.tpc_function_names must contain at least one entry.")
    TPC_FUNCTION_NAMES = [str(v) for v in tpc_names]

    BASE_REC_SELECTIONS = str(sel_he3.get("base_rec", BASE_REC_SELECTIONS))
    DEFAULT_REC_SELECTIONS = str(sel_he3.get("default_rec", DEFAULT_REC_SELECTIONS))
    HE4_BASE_SELECTION = str(sel_he4.get("base_rec", HE4_BASE_SELECTION))
    HE4_PRIMARY_SELECTION = str(sel_he4.get("primary_rec", HE4_PRIMARY_SELECTION))
    SECONDARY_SELECTION = str(sel_he3.get("secondary_rec", sel_he4.get("secondary_rec", SECONDARY_SELECTION)))
    HE3_TRIAL_DCA_SELECTION = str(sel_he3.get("trial_dca", HE3_TRIAL_DCA_SELECTION))
    HE3_NSIGMA_TOF_CUT = float(sel_he3.get("nsigma_tof", HE3_NSIGMA_TOF_CUT))
    HE4_NSIGMA_TOF_CUT = float(sel_he4.get("nsigma_tof", HE4_NSIGMA_TOF_CUT))
    SKIM_SELECTION_TEMPLATE = str(sel_common.get("skim_template", SKIM_SELECTION_TEMPLATE))
    HE3_MC_RECO_APPEND = str(sel_he3.get("mc_reco_append", HE3_MC_RECO_APPEND))
    HE3_MC_GEN_SELECTION = str(sel_he3.get("mc_gen", HE3_MC_GEN_SELECTION))
    HE4_MC_PID_SELECTION = str(sel_he4.get("mc_pid", HE4_MC_PID_SELECTION))
    HE4_MC_RECO_SELECTION = str(sel_he4.get("mc_reco", HE4_MC_RECO_SELECTION))
    HE4_MC_GEN_SELECTION = str(sel_he4.get("mc_gen", HE4_MC_GEN_SELECTION))
    HE4_MC_SIGNAL_TRACKING = str(sel_he4.get("mc_signal_tracking", HE4_MC_SIGNAL_TRACKING))

    CUT_NAMES["nsigmaDCAz"] = [float(v) for v in cuts.get("nsigmaDCAz", CUT_NAMES["nsigmaDCAz"])]
    CUT_NAMES["fTPCnCls"] = [float(v) for v in cuts.get("fTPCnCls", CUT_NAMES["fTPCnCls"])]
    CUT_NAMES["nITScls"] = [float(v) for v in cuts.get("nITScls", CUT_NAMES["nITScls"])]
    CUT_NAMES["nsigmaTPC"] = [float(v) for v in cuts.get("nsigmaTPC", CUT_NAMES["nsigmaTPC"])]

    _refresh_particle_configs()

    particle_cfg = cfg.get("particle", {})
    if isinstance(particle_cfg, dict):
        for key, override in particle_cfg.items():
            if not isinstance(override, dict):
                continue
            if key not in PARTICLE_CONFIGS:
                template = str(override.get("template", "he3"))
                if template not in PARTICLE_CONFIGS:
                    raise ValueError(
                        f"Unknown particle template '{template}' for particle '{key}'. "
                        f"Available templates: {', '.join(sorted(PARTICLE_CONFIGS))}."
                    )
                PARTICLE_CONFIGS[key] = copy.deepcopy(PARTICLE_CONFIGS[template])
            for pkey, pval in override.items():
                if pkey == "template":
                    continue
                PARTICLE_CONFIGS[key][pkey] = pval

    # Apply per-species selections for all configured particles.
    if isinstance(selections, dict):
        for species_key, species_cfg in PARTICLE_CONFIGS.items():
            sel_species = selections.get(species_key, {})
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

    recompute_derived_globals()


def get_particle_config(species: str) -> dict[str, Any]:
    if species not in PARTICLE_CONFIGS:
        raise ValueError(f"Unsupported species '{species}'. Available: {', '.join(sorted(PARTICLE_CONFIGS))}.")
    return copy.deepcopy(PARTICLE_CONFIGS[species])


@dataclass
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

    def get_particle_config(self, species: str) -> dict[str, Any]:
        if species not in self.particle_configs:
            raise ValueError(
                f"Unsupported species '{species}'. Available: {', '.join(sorted(self.particle_configs))}."
            )
        return copy.deepcopy(self.particle_configs[species])


def current_runtime_config() -> RuntimeConfig:
    return RuntimeConfig(
        variant=str(VARIANT),
        letter=copy.deepcopy(LETTER),
        filter_list_name=str(FILTER_LIST_NAME),
        n_pt_bins=int(N_PT_BINS),
        pt_bins=copy.deepcopy(PT_BINS),
        pt_bin_array=array("d", PT_BIN_ARRAY),
        cent_length=int(CENT_LENGTH),
        cent_pt_limits=copy.deepcopy(CENT_PT_LIMITS),
        tpc_max_pt=float(TPC_MAX_PT),
        tof_min_pt=float(TOF_MIN_PT),
        pt_range=copy.deepcopy(PT_RANGE),
        tpc_function_names=copy.deepcopy(TPC_FUNCTION_NAMES),
        cut_names=copy.deepcopy(CUT_NAMES),
        skim_selection_template=str(SKIM_SELECTION_TEMPLATE),
        particle_configs=copy.deepcopy(PARTICLE_CONFIGS),
    )


recompute_derived_globals()
_refresh_particle_configs()
