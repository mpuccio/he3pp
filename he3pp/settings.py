from array import array
from typing import Any

LETTER = ["M", "A"]
NAMES = ["he3", "antihe3"]
LABELS = ["^{3}He", "^{3}#bar{He}"]

MC_PRODUCTION = "LHC23j6b"
RECO_PASS = "apass4"
PERIOD = "LHC22"
VARIANT = "giovanni"
BASE_INPUT_DIR = "$NUCLEI_INPUT/"
BASE_OUTPUT_ROOT = "$NUCLEI_OUTPUT/"
FILTER_LIST_NAME = "nuclei"

BASE_REC_SELECTIONS = "fTPCnCls >= 110 && nITScls >= 5 && std::abs(fEta) < 0.9 && std::abs(fDCAxy) < 0.7 && pt > 0.8 && pt < 9.0"
DEFAULT_REC_SELECTIONS = "fTPCnCls > 120 && nITScls >= 6 && std::abs(nsigmaDCAz) < 7 && std::abs(fDCAxy) < 0.2"
HE4_BASE_SELECTION = "fTPCnCls >= 110 && nITScls >= 5 && abs(fEta) < 0.9 && std::abs(fDCAxy) < 0.7 && ptHe4 > 0.5 && ptHe4 < 9.0"
HE4_PRIMARY_SELECTION = "fTPCnCls > 120 && nITScls >= 6 && std::abs(nsigmaDCAz) < 7 && std::abs(fDCAxy) < 0.2"
SECONDARY_SELECTION = "fTPCnCls > 120 && nITScls >= 6 && std::abs(nsigmaDCAz) > 7 && std::abs(fDCAxy) < 0.2"
HE3_TRIAL_DCA_SELECTION = "std::abs(fDCAxy) < 0.2"
HE3_NSIGMA_TOF_CUT = 3.5
HE4_NSIGMA_TOF_CUT = 3.0
SKIM_SELECTION_TEMPLATE = "std::abs(nsigmaDCAz) < 8 && std::abs(fDCAxy) < 0.2 && std::abs({nsigma}) < 5"
HE3_MC_RECO_APPEND = "&& isPrimary"
HE3_MC_GEN_SELECTION = "isPrimary && std::abs(yMC) < 0.5"
HE4_MC_PID_SELECTION = "isHe4 && isPrimary"
HE4_MC_RECO_SELECTION = "nITScls > 4 && fTPCnCls > 110 && std::abs(fEta) < 0.9 && std::abs(rapidity) < 0.5"
HE4_MC_GEN_SELECTION = "std::abs(yMC) < 0.5"
HE4_MC_SIGNAL_TRACKING = "fTPCnCls > 120 && nITScls >= 6 && std::abs(fDCAz) < 0.7"

N_PT_BINS = 12
PT_BINS = [1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5.0]
PT_BIN_ARRAY = array("d", PT_BINS)
CENT_LENGTH = 1
CENT_PT_LIMITS = [7.0]
TPC_MAX_PT = 7.0
TOF_MIN_PT = 1.0
PT_RANGE = [1.4, 7.0]
NTPC_FUNCTIONS = 3
TPC_FUNCTION_NAMES = ["GausGaus", "ExpGaus", "ExpTailGaus", "LognormalLognormal"]
# Weighted efficiency histograms are intentionally named "Weff*".
WEIGHTED_EFF_NAMING_POLICY = "prefix_W"

CUT_NAMES = {
    "nsigmaDCAz": [6.0, 7.0, 8.0],
    "fTPCnCls": [110.0, 120.0, 130.0],
    "nITScls": [5.0, 6.0, 7.0],
    "nsigmaTPC": [3.0, 4.0, 5.0],
}

BASE_OUTPUT_DIR = ""
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


def recompute_derived_globals() -> None:
    global BASE_OUTPUT_DIR
    global DATA_TREE_FILENAME, DATA_FILENAME, DATA_FILENAME_HE4, DATA_ANALYSIS_RESULTS
    global MC_ANALYSIS_RESULTS, MC_TREE_FILENAME, MC_FILENAME, MC_FILENAME_HE4
    global SIGNAL_OUTPUT, SYSTEMATICS_OUTPUT

    BASE_OUTPUT_DIR = f"{BASE_OUTPUT_ROOT}{PERIOD}/{RECO_PASS}/"
    DATA_TREE_FILENAME = f"{BASE_INPUT_DIR}data/{PERIOD}/{RECO_PASS}/MergedAO2D.root"
    DATA_FILENAME = f"{BASE_INPUT_DIR}data/{PERIOD}/{RECO_PASS}/DataHistos{VARIANT}.root"
    DATA_FILENAME_HE4 = f"{BASE_INPUT_DIR}data/{PERIOD}/{RECO_PASS}/DataHistosHe4{VARIANT}.root"
    DATA_ANALYSIS_RESULTS = f"{BASE_INPUT_DIR}data/{PERIOD}/{RECO_PASS}/AnalysisResults.root"
    MC_ANALYSIS_RESULTS = f"{BASE_INPUT_DIR}MC/{MC_PRODUCTION}/AnalysisResults.root"
    MC_TREE_FILENAME = f"{BASE_INPUT_DIR}MC/{MC_PRODUCTION}/MergedAO2D.root"
    MC_FILENAME = f"{BASE_INPUT_DIR}MC/{MC_PRODUCTION}/MChistos{VARIANT}.root"
    MC_FILENAME_HE4 = f"{BASE_INPUT_DIR}MC/{MC_PRODUCTION}/MChistosHe4{VARIANT}.root"
    SIGNAL_OUTPUT = f"{BASE_OUTPUT_DIR}signal{VARIANT}.root"
    SYSTEMATICS_OUTPUT = f"{BASE_OUTPUT_DIR}systematics{VARIANT}.root"


def apply_runtime_overrides(cfg: dict[str, Any]) -> None:
    global MC_PRODUCTION, RECO_PASS, PERIOD, VARIANT, BASE_INPUT_DIR, BASE_OUTPUT_ROOT, FILTER_LIST_NAME
    global BASE_REC_SELECTIONS, DEFAULT_REC_SELECTIONS, HE4_BASE_SELECTION, HE4_PRIMARY_SELECTION
    global SECONDARY_SELECTION, HE3_TRIAL_DCA_SELECTION, HE3_NSIGMA_TOF_CUT, HE4_NSIGMA_TOF_CUT
    global SKIM_SELECTION_TEMPLATE
    global HE3_MC_RECO_APPEND, HE3_MC_GEN_SELECTION, HE4_MC_PID_SELECTION, HE4_MC_RECO_SELECTION
    global HE4_MC_GEN_SELECTION, HE4_MC_SIGNAL_TRACKING
    global N_PT_BINS, PT_BINS, PT_BIN_ARRAY, CENT_PT_LIMITS, TPC_MAX_PT, TOF_MIN_PT, PT_RANGE
    global TPC_FUNCTION_NAMES, CUT_NAMES

    common = cfg.get("common", {})
    selections = cfg.get("selections", {})
    cuts = cfg.get("cuts", {})

    MC_PRODUCTION = str(common.get("mc_production", MC_PRODUCTION))
    RECO_PASS = str(common.get("reco_pass", RECO_PASS))
    PERIOD = str(common.get("period", PERIOD))
    VARIANT = str(common.get("variant", VARIANT))
    BASE_INPUT_DIR = str(common.get("base_input_dir", BASE_INPUT_DIR))
    BASE_OUTPUT_ROOT = str(common.get("base_output_root", BASE_OUTPUT_ROOT))
    FILTER_LIST_NAME = str(common.get("filter_list_name", FILTER_LIST_NAME))

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
    if len(tpc_names) < NTPC_FUNCTIONS:
        raise ValueError(f"common.tpc_function_names must have at least {NTPC_FUNCTIONS} entries.")
    TPC_FUNCTION_NAMES = [str(v) for v in tpc_names]

    BASE_REC_SELECTIONS = str(selections.get("base_rec", BASE_REC_SELECTIONS))
    DEFAULT_REC_SELECTIONS = str(selections.get("default_rec", DEFAULT_REC_SELECTIONS))
    HE4_BASE_SELECTION = str(selections.get("he4_base_rec", HE4_BASE_SELECTION))
    HE4_PRIMARY_SELECTION = str(selections.get("he4_primary_rec", HE4_PRIMARY_SELECTION))
    SECONDARY_SELECTION = str(selections.get("secondary_rec", SECONDARY_SELECTION))
    HE3_TRIAL_DCA_SELECTION = str(selections.get("he3_trial_dca", HE3_TRIAL_DCA_SELECTION))
    HE3_NSIGMA_TOF_CUT = float(selections.get("he3_nsigma_tof", HE3_NSIGMA_TOF_CUT))
    HE4_NSIGMA_TOF_CUT = float(selections.get("he4_nsigma_tof", HE4_NSIGMA_TOF_CUT))
    SKIM_SELECTION_TEMPLATE = str(selections.get("skim_template", SKIM_SELECTION_TEMPLATE))
    HE3_MC_RECO_APPEND = str(selections.get("he3_mc_reco_append", HE3_MC_RECO_APPEND))
    HE3_MC_GEN_SELECTION = str(selections.get("he3_mc_gen", HE3_MC_GEN_SELECTION))
    HE4_MC_PID_SELECTION = str(selections.get("he4_mc_pid", HE4_MC_PID_SELECTION))
    HE4_MC_RECO_SELECTION = str(selections.get("he4_mc_reco", HE4_MC_RECO_SELECTION))
    HE4_MC_GEN_SELECTION = str(selections.get("he4_mc_gen", HE4_MC_GEN_SELECTION))
    HE4_MC_SIGNAL_TRACKING = str(selections.get("he4_mc_signal_tracking", HE4_MC_SIGNAL_TRACKING))

    CUT_NAMES["nsigmaDCAz"] = [float(v) for v in cuts.get("nsigmaDCAz", CUT_NAMES["nsigmaDCAz"])]
    CUT_NAMES["fTPCnCls"] = [float(v) for v in cuts.get("fTPCnCls", CUT_NAMES["fTPCnCls"])]
    CUT_NAMES["nITScls"] = [float(v) for v in cuts.get("nITScls", CUT_NAMES["nITScls"])]
    CUT_NAMES["nsigmaTPC"] = [float(v) for v in cuts.get("nsigmaTPC", CUT_NAMES["nsigmaTPC"])]

    recompute_derived_globals()


recompute_derived_globals()
