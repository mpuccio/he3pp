from array import array
import logging

import ROOT

from .root_io import ensure_parent, expand
from .settings import RuntimeConfig
from .tasks_common import default_tpc_model_name, find_object_by_class


LOGGER = logging.getLogger("he3pp.tasks")


def systematics(
    signal_file: str,
    mc_file: str,
    data_analysis_results: str,
    output_file: str,
    particle: str = "he3",
    *,
    runtime_config: RuntimeConfig,
) -> None:
    cfg = runtime_config
    LOGGER.info("systematics start particle=%s signal=%s mc=%s output=%s", particle, signal_file, mc_file, output_file)
    p = cfg.get_particle_config(particle)
    species_names = [str(p["name"]), str(p["anti_name"])]
    f_data = ROOT.TFile(expand(signal_file))
    f_mc = ROOT.TFile(expand(mc_file))
    an_results = ROOT.TFile(expand(data_analysis_results))

    norm = None
    zorro_summary = find_object_by_class(an_results, "ZorroSummary")
    if zorro_summary:
        for method_name in ("getNormalisationFactor", "GetNormalisationFactor"):
            method = getattr(zorro_summary, method_name, None)
            if callable(method):
                norm = float(method())
                LOGGER.info("systematics normalization from ZorroSummary.%s(): %.6g", method_name, norm)
                break

    h_ntvx = an_results.Get("bc-selection-task/hCounterTVX")
    if not h_ntvx:
        h_ntvx = an_results.Get("bc-selection-task/hCounterTVXafterBCcuts")
    h_nvtx = an_results.Get("nuclei-spectra/spectra/hRecVtxZData")
    if not h_nvtx:
        raise RuntimeError("Missing normalization histogram: nuclei-spectra/spectra/hRecVtxZData")
    if h_ntvx and hasattr(h_ntvx, "GetEntries"):
        _ = h_ntvx.GetEntries()
    if norm is None:
        norm = h_nvtx.GetEntries()
        LOGGER.warning("systematics fallback normalization from hRecVtxZData entries: %.6g", norm)
    if norm <= 0:
        raise RuntimeError(f"Invalid normalization factor: {norm}")

    syst_tpc = [ROOT.TH2D(f"systTPC{species_names[i]}", ";#it{p}_{T} (GeV/#it{c});Relative systematics TPC", cfg.n_pt_bins, cfg.pt_bin_array, 50, -0.5, 0.5) for i in range(2)]
    syst_tof = [ROOT.TH2D(f"systTOF{species_names[i]}", ";#it{p}_{T} (GeV/#it{c});Relative systematics TOF", cfg.n_pt_bins, cfg.pt_bin_array, 50, -0.5, 0.5) for i in range(2)]

    default_eff_tpc = [None, None]
    default_eff_tof = [None, None]
    default_tpc = [None, None]
    default_tof = [None, None]
    default_tpc_uncorr = [None, None]
    default_tof_uncorr = [None, None]

    for key in f_data.GetListOfKeys():
        key_name = key.GetName()
        if cfg.filter_list_name not in key_name:
            continue

        list_data = f_data.Get(key_name)
        list_mc = f_mc.Get(key_name)

        h_data_tof = [[[None for _ in range(3)] for _ in range(2)] for _ in range(2)]
        tpc_model_names = [str(name) for name in cfg.tpc_function_names]
        chosen_tpc_model_name = default_tpc_model_name(cfg)
        h_data_tpc = [[None for _ in range(len(tpc_model_names))] for _ in range(2)]
        h_eff_tpc = [None, None]
        h_eff_tof = [None, None]

        tof_names = ["hRawCounts", "hRawCountsBinCounting"]
        tof_presel = ["", "_loose", "_tight"]

        for i_s in range(2):
            data_dir = list_data.Get(species_names[i_s])
            h_eff_tpc[i_s] = list_mc.Get(f"effTPC{cfg.letter[i_s]}")
            h_eff_tof[i_s] = list_mc.Get(f"effTOF{cfg.letter[i_s]}")

            if default_eff_tpc[i_s] is None:
                default_eff_tpc[i_s] = h_eff_tpc[i_s].Clone(f"defaultEffTPC{species_names[i_s]}")
            if default_eff_tof[i_s] is None:
                default_eff_tof[i_s] = h_eff_tof[i_s].Clone(f"defaultEffTOF{species_names[i_s]}")

            for i_tof in range(2):
                name = tof_names[i_tof] + cfg.letter[i_s]
                h_data_tof[i_s][i_tof][0] = data_dir.Get(f"GausExp/{name}0{tof_presel[0]}")
                if default_tof_uncorr[i_s] is None:
                    default_tof_uncorr[i_s] = h_data_tof[i_s][i_tof][0].Clone(f"defaultTOFuncorr{species_names[i_s]}")
                h_data_tof[i_s][i_tof][0].Divide(h_eff_tof[i_s])
                if default_tof[i_s] is None:
                    default_tof[i_s] = h_data_tof[i_s][i_tof][0].Clone(f"defaultTOF{species_names[i_s]}")

            for i_tpc, name in enumerate(tpc_model_names):
                h_data_tpc[i_s][i_tpc] = data_dir.Get(f"TPConly/hTPConly{cfg.letter[i_s]}0_{name}")
                if h_data_tpc[i_s][i_tpc] is None:
                    raise RuntimeError(f"Missing TPC signal histogram for model '{name}' in {species_names[i_s]}.")
                if default_tpc_uncorr[i_s] is None and name == chosen_tpc_model_name:
                    default_tpc_uncorr[i_s] = h_data_tpc[i_s][i_tpc].Clone(f"defaultTPCuncorr{species_names[i_s]}")
                h_data_tpc[i_s][i_tpc].Divide(h_eff_tpc[i_s])
                if default_tpc[i_s] is None and name == chosen_tpc_model_name:
                    default_tpc[i_s] = h_data_tpc[i_s][i_tpc].Clone(f"defaultTPC{species_names[i_s]}")
            if default_tpc_uncorr[i_s] is None or default_tpc[i_s] is None:
                raise RuntimeError(
                    f"Default TPC model '{chosen_tpc_model_name}' not found for {species_names[i_s]}. "
                    f"Available models: {', '.join(tpc_model_names)}."
                )

        for i_s in range(2):
            for i_b in range(1, cfg.n_pt_bins + 1):
                pt = h_data_tpc[i_s][0].GetBinCenter(i_b)
                default_val_tpc = default_tpc[i_s].GetBinContent(i_b)
                default_val_tof = default_tof[i_s].GetBinContent(i_b)
                if default_val_tpc != 0:
                    for i_tpc in range(len(tpc_model_names)):
                        value = h_data_tpc[i_s][i_tpc].GetBinContent(i_b)
                        syst_tpc[i_s].Fill(pt, (value - default_val_tpc) / default_val_tpc)
                if default_val_tof != 0:
                    for i_tof in range(2):
                        value = h_data_tof[i_s][i_tof][0].GetBinContent(i_b)
                        syst_tof[i_s].Fill(pt, (value - default_val_tof) / default_val_tof)

    h_syst_tpc = [ROOT.TH1D(f"hSystTPC{cfg.letter[i]}", ";#it{p}_{T} (GeV/#it{c});Relative systematics TPC", cfg.n_pt_bins, cfg.pt_bin_array) for i in range(2)]
    h_syst_tof = [ROOT.TH1D(f"hSystTOF{cfg.letter[i]}", ";#it{p}_{T} (GeV/#it{c});Relative systematics TOF", cfg.n_pt_bins, cfg.pt_bin_array) for i in range(2)]

    for i_s in range(2):
        for i_b in range(1, cfg.n_pt_bins + 1):
            h_syst_tpc[i_s].SetBinContent(i_b, syst_tpc[i_s].ProjectionY("", i_b, i_b).GetRMS())
            h_syst_tof[i_s].SetBinContent(i_b, syst_tof[i_s].ProjectionY("", i_b, i_b).GetRMS())

    out_name = expand(output_file)
    ensure_parent(out_name)
    out = ROOT.TFile(out_name, "recreate")

    tof_matching_m = default_tof_uncorr[0].Clone(f"TOFmatching{species_names[0]}")
    tof_matching_a = default_tof_uncorr[1].Clone(f"TOFmatching{species_names[1]}")
    tof_matching_m.Divide(default_tpc_uncorr[0])
    tof_matching_a.Divide(default_tpc_uncorr[1])
    tof_matching_m.Write()
    tof_matching_a.Write()

    for h in [syst_tpc[0], syst_tpc[1], syst_tof[0], syst_tof[1], h_syst_tpc[0], h_syst_tpc[1], h_syst_tof[0], h_syst_tof[1]]:
        h.Write()

    pub_bins = array("d", [1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0])
    y = [1.2241e-07, 8.4801e-08, 5.0085e-08, 3.2333e-08, 1.7168e-08, 4.8137e-09]
    ey = [1.769e-08, 7.5127e-09, 6.0035e-09, 4.8788e-09, 2.5057e-09, 1.3356e-09]
    sy = [1.3346e-08, 1.0763e-08, 3.2452e-09, 2.1084e-09, 1.1316e-09, 3.1345e-10]

    h_pub = ROOT.TH1D("hPub", ";#it{p}_{T} (GeV/#it{c});#frac{1}{N_{ev}} #frac{d^{2}N}{dyd#it{p}_{T}}", 6, pub_bins)
    h_pub_syst = ROOT.TH1D("hPubSyst", ";#it{p}_{T} (GeV/#it{c});#frac{1}{N_{ev}} #frac{d^{2}N}{dyd#it{p}_{T}}", 6, pub_bins)
    for i in range(1, 7):
        h_pub.SetBinContent(i, y[i - 1])
        h_pub.SetBinError(i, ey[i - 1])
        h_pub_syst.SetBinContent(i, y[i - 1])
        h_pub_syst.SetBinError(i, sy[i - 1])
    h_pub.Write("pubStat")
    h_pub_syst.Write("pubSyst")

    for i_s in range(2):
        f_stat_tpc = default_tpc_uncorr[i_s].Clone(f"fStatTPC{cfg.letter[i_s]}")
        f_syst_tpc = default_tpc_uncorr[i_s].Clone(f"fSystTPC{cfg.letter[i_s]}")
        f_stat_tof = default_tof_uncorr[i_s].Clone(f"fStatTOF{cfg.letter[i_s]}")
        f_syst_tof = default_tof_uncorr[i_s].Clone(f"fSystTOF{cfg.letter[i_s]}")

        for i_bin in range(1, cfg.n_pt_bins + 1):
            yield_tpc = default_tpc_uncorr[i_s].GetBinContent(i_bin)
            yield_tof = default_tof_uncorr[i_s].GetBinContent(i_bin)
            eff_tpc = default_eff_tpc[i_s].GetBinContent(i_bin)
            eff_tof = default_eff_tof[i_s].GetBinContent(i_bin)

            if eff_tpc >= 1e-2:
                f_stat_tpc.SetBinContent(i_bin, yield_tpc / eff_tpc)
                f_stat_tpc.SetBinError(i_bin, default_tpc_uncorr[i_s].GetBinError(i_bin) / eff_tpc)
                f_syst_tpc.SetBinContent(i_bin, yield_tpc / eff_tpc)
                f_syst_tpc.SetBinError(i_bin, h_syst_tpc[i_s].GetBinContent(i_bin) * yield_tpc / eff_tpc)
            else:
                f_stat_tpc.SetBinContent(i_bin, 0.0)
                f_stat_tpc.SetBinError(i_bin, 0.0)
                f_syst_tpc.SetBinContent(i_bin, 0.0)
                f_syst_tpc.SetBinError(i_bin, 0.0)

            if eff_tof >= 1e-2:
                f_stat_tof.SetBinContent(i_bin, yield_tof / eff_tof)
                f_stat_tof.SetBinError(i_bin, default_tof_uncorr[i_s].GetBinError(i_bin) / eff_tof)
                f_syst_tof.SetBinContent(i_bin, yield_tof / eff_tof)
                f_syst_tof.SetBinError(i_bin, h_syst_tof[i_s].GetBinContent(i_bin) * yield_tof / eff_tof)
            else:
                f_stat_tof.SetBinContent(i_bin, 0.0)
                f_stat_tof.SetBinError(i_bin, 0.0)
                f_syst_tof.SetBinContent(i_bin, 0.0)
                f_syst_tof.SetBinError(i_bin, 0.0)

        f_stat_tpc.Scale(1.0 / norm, "width")
        f_syst_tpc.Scale(1.0 / norm, "width")
        f_stat_tof.Scale(1.0 / norm, "width")
        f_syst_tof.Scale(1.0 / norm, "width")

        f_stat_tpc.Write()
        f_syst_tpc.Write()
        f_stat_tof.Write()
        f_syst_tof.Write()

    out.Close()
    LOGGER.info("systematics done particle=%s output=%s", particle, output_file)
