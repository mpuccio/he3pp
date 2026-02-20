from array import array
from datetime import datetime
import logging
from typing import Any

import ROOT

from .settings import *  # Loaded after runtime overrides in cli.run
from .root_io import build_rdf_from_ao2d, define_columns_for_data, ensure_parent, expand, h1_model, h2_model, load_fit_modules, ptr, write_hist

LOGGER = logging.getLogger("he3pp.tasks")


def _weighted_eff_name(base: str) -> str:
    # Weighted efficiencies intentionally use "Weff*" naming.
    return f"W{base}"


def _collect_rresult_ptrs(obj: Any) -> list[Any]:
    out: list[Any] = []
    if isinstance(obj, dict):
        for v in obj.values():
            out.extend(_collect_rresult_ptrs(v))
        return out
    if isinstance(obj, (list, tuple)):
        for v in obj:
            out.extend(_collect_rresult_ptrs(v))
        return out
    if hasattr(obj, "GetValue") and hasattr(obj, "GetPtr"):
        out.append(obj)
    return out


def _run_graphs(actions: list[Any]) -> None:
    if not actions:
        return
    run_graphs = getattr(getattr(ROOT, "RDF", None), "RunGraphs", None)
    if run_graphs:
        run_graphs(actions)
        return
    # Fallback for ROOT builds without RunGraphs.
    for action in actions:
        action.GetValue()


def _find_object_by_class(root_dir: Any, class_name: str) -> Any | None:
    """Depth-first search for the first object with the requested ROOT class."""
    if not root_dir or not hasattr(root_dir, "GetListOfKeys"):
        return None
    for key in root_dir.GetListOfKeys():
        obj = key.ReadObj()
        if obj and hasattr(obj, "ClassName") and obj.ClassName() == class_name:
            return obj
        if obj and hasattr(obj, "InheritsFrom") and obj.InheritsFrom("TDirectory"):
            found = _find_object_by_class(obj, class_name)
            if found:
                return found
    return None


def _book_data_species(df_all: Any, particle: str, skim: bool = False, tag: str = "") -> dict[str, Any]:
    if particle == "he4":
        df_base = df_all.Filter(HE4_BASE_SELECTION)
        df_primary = df_base.Filter(HE4_PRIMARY_SELECTION)
        df_secondary = df_base.Filter(SECONDARY_SELECTION)
        nsigma = "nsigmaHe4"
        mass_cut = "hasGoodTOFmassHe4"
        pt_name = "ptHe4"
        dmass = "deltaMassHe4"
        nominal_mass = 3.72738
        tpc_anti_model = ROOT.RDF.TH2DModel(f"fATPCcounts{tag}", ";#it{p}_{T}^{rec} (GeV/#it{c});^{4}#bar{He} n#sigma_{TPC};Counts", 160, 0.5, 4.5, 100, -5, 5)
        tpc_mat_model = ROOT.RDF.TH2DModel(f"fMTPCcounts{tag}", ";#it{p}_{T}^{rec} (GeV/#it{c});^{4}He n#sigma_{TPC};Counts", N_PT_BINS, PT_BIN_ARRAY, 100, -5, 5)
    else:
        df_base = df_all.Filter(BASE_REC_SELECTIONS)
        df_primary = df_base.Filter(DEFAULT_REC_SELECTIONS)
        df_secondary = df_base.Filter(SECONDARY_SELECTION)
        nsigma = "nsigmaHe3"
        mass_cut = "hasGoodTOFmassHe3"
        pt_name = "pt"
        dmass = "deltaMassHe3"
        nominal_mass = 2.80839
        tpc_anti_model = ROOT.RDF.TH2DModel(f"fATPCcounts{tag}", ";#it{p}_{T}^{rec} (GeV/#it{c});^{3}#bar{He} n#sigma_{TPC};Counts", N_PT_BINS, PT_BIN_ARRAY, 100, -5, 5)
        tpc_mat_model = ROOT.RDF.TH2DModel(f"fMTPCcounts{tag}", ";#it{p}_{T}^{rec} (GeV/#it{c});^{3}He n#sigma_{TPC};Counts", N_PT_BINS, PT_BIN_ARRAY, 100, -5, 5)

    if skim:
        df_base.Filter(SKIM_SELECTION_TEMPLATE.format(nsigma=nsigma)).Snapshot("nucleiTree", "data/skimmed.root")

    matter_map = {"A": "!matter", "M": "matter"}
    h_dca_xy: dict[str, list[Any]] = {}
    h_dca_z: dict[str, list[Any]] = {}
    h_dca_xy_secondary: dict[str, list[Any]] = {}
    h_tpc: dict[str, list[Any]] = {}
    h_tof: dict[str, list[Any]] = {}

    for label, matter_sel in matter_map.items():
        h_dca_xy[label] = [df_primary.Filter(f"{matter_sel} && {nsigma} > -0.5 && {nsigma} < 3 && {mass_cut}").Histo2D(ROOT.RDF.TH2DModel(f"hDCAxy{label}He{tag}", ";#it{p}_{T}^{rec} (GeV/#it{c});DCA_{xy} (cm);Counts", N_PT_BINS, PT_BIN_ARRAY, 100, -0.2, 0.2), pt_name, "fDCAxy")]
        h_dca_z[label] = [df_primary.Filter(f"{matter_sel} && {nsigma} > -0.5 && {nsigma} < 3 && {mass_cut}").Histo2D(ROOT.RDF.TH2DModel(f"hDCAz{label}He{tag}", ";#it{p}_{T}^{rec} (GeV/#it{c});DCA_{z} (cm);Counts", N_PT_BINS, PT_BIN_ARRAY, 100, -0.2, 0.2), pt_name, "fDCAz")]
        h_dca_xy_secondary[label] = [df_secondary.Filter(f"{matter_sel} && {nsigma} > -0.5 && {nsigma} < 3 && {mass_cut}").Histo2D(ROOT.RDF.TH2DModel(f"hDCAxySecondary{label}He{tag}", ";#it{p}_{T}^{rec} (GeV/#it{c});DCA_{xy} (cm);Counts", N_PT_BINS, PT_BIN_ARRAY, 100, -0.2, 0.2), pt_name, "fDCAxy")]

        tpc_sel = f"{matter_sel} && {mass_cut}" if particle == "he3" else matter_sel
        tpc_model = tpc_anti_model if label == "A" else tpc_mat_model
        tpc_pt_name = "ptUncorr" if (particle == "he4" and label == "A") else pt_name
        h_tpc[label] = [df_primary.Filter(tpc_sel).Histo2D(tpc_model, tpc_pt_name, nsigma)]

    nsigma_tof_cut = HE4_NSIGMA_TOF_CUT if particle == "he4" else HE3_NSIGMA_TOF_CUT
    for label, matter_sel in matter_map.items():
        isotope = "4" if particle == "he4" else "3"
        species_label = f"^{{{isotope}}}#bar{{He}}" if label == "A" else f"^{{{isotope}}}He"
        h_tof[label] = [df_primary.Filter(f"{matter_sel} && std::abs({nsigma}) < {nsigma_tof_cut}").Histo2D(ROOT.RDF.TH2DModel(f"f{label}TOFsignal{tag}", f";#it{{p}}_{{T}}^{{rec}} (GeV/#it{{c}});m_{{TOF}}-m_{{{species_label}}};Counts", N_PT_BINS, PT_BIN_ARRAY, 100, -0.9, 1.1), pt_name, dmass)]

    df_primary = df_primary.Define(f"tofMassDeltaPOI{tag}", f"tofMass - {nominal_mass}")
    tof_delta_mass_range = (-0.6, 0.6) if particle == "he4" else (-1.0, 1.0)
    h_tof_mass_vs_tpc_nsigma: dict[str, list[Any]] = {"A": [], "M": []}
    for i_pt in range(N_PT_BINS):
        pt_low = PT_BINS[i_pt]
        pt_high = PT_BINS[i_pt + 1]
        for label, matter_sel in matter_map.items():
            h_tof_mass_vs_tpc_nsigma[label].append(
                df_primary.Filter(f"{matter_sel} && hasTOF && std::abs({nsigma}) < 6 && {pt_name} >= {pt_low} && {pt_name} < {pt_high}").Histo2D(
                    ROOT.RDF.TH2DModel(
                        f"hTOFMassVsTPCnsigma{label}_pt{i_pt:02d}{tag}",
                        ";n#sigma_{TPC};m_{TOF}-m_{0} (GeV/#it{c}^{2});Counts",
                        120,
                        -6.0,
                        6.0,
                        120,
                        tof_delta_mass_range[0],
                        tof_delta_mass_range[1],
                    ),
                    nsigma,
                    f"tofMassDeltaPOI{tag}",
                )
            )

    i_trial = 0
    if particle == "he3":
        for cut_dca_z in CUT_NAMES["nsigmaDCAz"]:
            df_dca_z = df_base.Filter(f"std::abs(nsigmaDCAz) < {cut_dca_z}")
            for cut_tpc in CUT_NAMES["fTPCnCls"]:
                df_tpc = df_dca_z.Filter(f"fTPCnCls > {cut_tpc}")
                for cut_its in CUT_NAMES["nITScls"]:
                    df_its = df_tpc.Filter(f"nITScls >= {cut_its}")
                    for label, matter_sel in matter_map.items():
                        species_label = "^{3}#bar{He}" if label == "A" else "^{3}He"
                        h_dca_xy[label].append(df_its.Filter(f"{matter_sel} && {nsigma} > -0.5 && {nsigma} < 3 && {mass_cut}").Histo2D(ROOT.RDF.TH2DModel(f"hDCAxy{label}He3{i_trial}{tag}", ";#it{p}_{T}^{rec} (GeV/#it{c});DCA_{xy} (cm);Counts", N_PT_BINS, PT_BIN_ARRAY, 560, -0.7, 0.7), pt_name, "fDCAxy"))
                        h_dca_z[label].append(df_its.Filter(f"{matter_sel} && {nsigma} > -0.5 && {nsigma} < 3 && {mass_cut}").Histo2D(ROOT.RDF.TH2DModel(f"hDCAz{label}He3{i_trial}{tag}", ";#it{p}_{T}^{rec} (GeV/#it{c});DCA_{z} (cm);Counts", N_PT_BINS, PT_BIN_ARRAY, 560, -0.7, 0.7), pt_name, "fDCAz"))
                        h_tpc[label].append(df_its.Filter(f"{matter_sel} && {HE3_TRIAL_DCA_SELECTION} && {mass_cut}").Histo2D(ROOT.RDF.TH2DModel(f"f{label}TPCcounts{i_trial}{tag}", f";#it{{p}}_{{T}}^{{rec}} (GeV/#it{{c}});{species_label} n#sigma_{{TPC}};Counts", N_PT_BINS, PT_BIN_ARRAY, 100, -5, 5), pt_name, nsigma))
                        h_tof[label].append(df_its.Filter(f"{matter_sel} && {HE3_TRIAL_DCA_SELECTION} && std::abs({nsigma}) < {HE3_NSIGMA_TOF_CUT}").Histo2D(ROOT.RDF.TH2DModel(f"f{label}TOFsignal{i_trial}{tag}", f";#it{{p}}_{{T}}^{{rec}} (GeV/#it{{c}});m_{{TOF}}-m_{{{species_label}}};Counts", N_PT_BINS, PT_BIN_ARRAY, 100, -0.9, 1.1), pt_name, dmass))
                    i_trial += 1

    return {
        "particle": particle,
        "h_tpc": h_tpc,
        "h_tof": h_tof,
        "h_dca_xy": h_dca_xy,
        "h_dca_z": h_dca_z,
        "h_dca_xy_secondary": h_dca_xy_secondary,
        "h_tof_mass_vs_tpc_nsigma": h_tof_mass_vs_tpc_nsigma,
        "i_trial": i_trial,
    }


def _write_data_bundle(bundle: dict[str, Any], output_file: str) -> None:
    output_file = expand(output_file)
    ensure_parent(output_file)
    out = ROOT.TFile(output_file, "recreate")
    d0 = out.mkdir("nuclei")
    d0.cd()
    for label in ("A", "M"):
        write_hist(bundle["h_tpc"][label][0], f"f{label}TPCcounts")
        write_hist(bundle["h_tof"][label][0], f"f{label}TOFsignal")
        write_hist(bundle["h_dca_xy"][label][0], f"hDCAxy{label}He3")
        write_hist(bundle["h_dca_z"][label][0], f"hDCAz{label}He3")
    for label in ("M", "A"):
        write_hist(bundle["h_dca_xy_secondary"][label][0], f"hDCAxySecondary{label}He3")
    for i_pt in range(N_PT_BINS):
        for label in ("A", "M"):
            write_hist(bundle["h_tof_mass_vs_tpc_nsigma"][label][i_pt], f"hTOFMassVsTPCnsigma{label}_pt{i_pt:02d}")

    for idx in range(bundle["i_trial"]):
        dtrial = out.mkdir(f"nuclei{idx}")
        dtrial.cd()
        for label in ("A", "M"):
            write_hist(bundle["h_tpc"][label][idx + 1], f"f{label}TPCcounts")
            write_hist(bundle["h_tof"][label][idx + 1], f"f{label}TOFsignal")
            write_hist(bundle["h_dca_xy"][label][idx + 1], f"hDCAxy{label}He3")
            write_hist(bundle["h_dca_z"][label][idx + 1], f"hDCAz{label}He3")
    out.Close()


def analyse_data(input_file: str, output_file: str, particle: str, skim: bool = False, draw: bool = False) -> None:
    LOGGER.info("analyse_data start particle=%s input=%s output=%s", particle, input_file, output_file)
    ROOT.gStyle.SetOptStat(0)
    rdf = build_rdf_from_ao2d("O2nucleitable", input_file)
    df_all = define_columns_for_data(rdf)
    bundle = _book_data_species(df_all, particle, skim=skim, tag=f"_{particle}")
    _run_graphs(
        _collect_rresult_ptrs(bundle["h_tpc"])
        + _collect_rresult_ptrs(bundle["h_tof"])
        + _collect_rresult_ptrs(bundle["h_dca_xy"])
        + _collect_rresult_ptrs(bundle["h_dca_z"])
        + _collect_rresult_ptrs(bundle["h_dca_xy_secondary"])
        + _collect_rresult_ptrs(bundle["h_tof_mass_vs_tpc_nsigma"])
    )
    if draw:
        draw_order = [
            (bundle["h_tpc"], ("A", "M")),
            (bundle["h_tof"], ("A", "M")),
            (bundle["h_dca_xy"], ("A", "M")),
            (bundle["h_dca_z"], ("A", "M")),
            (bundle["h_dca_xy_secondary"], ("M", "A")),
        ]
        for hist_map, labels in draw_order:
            for label in labels:
                hist = hist_map[label][0]
                c = ROOT.TCanvas()
                c.SetRightMargin(0.15)
                hist.DrawClone("colz")
    _write_data_bundle(bundle, output_file)
    LOGGER.info("analyse_data done output=%s", output_file)


def analyse_data_multi(input_file: str, outputs: dict[str, str], skim: bool = False, draw: bool = False) -> None:
    particles = [p for p in ("he3", "he4") if p in outputs]
    if not particles:
        raise ValueError("analyse_data_multi requires at least one output among {'he3','he4'}.")
    LOGGER.info("analyse_data_multi start input=%s particles=%s", input_file, particles)
    ROOT.gStyle.SetOptStat(0)
    rdf = build_rdf_from_ao2d("O2nucleitable", input_file)
    df_all = define_columns_for_data(rdf)

    bundles: dict[str, dict[str, Any]] = {}
    all_actions: list[Any] = []
    for particle in particles:
        bundle = _book_data_species(df_all, particle, skim=(skim and particle == "he3"), tag=f"_{particle}")
        bundles[particle] = bundle
        all_actions += (
            _collect_rresult_ptrs(bundle["h_tpc"])
            + _collect_rresult_ptrs(bundle["h_tof"])
            + _collect_rresult_ptrs(bundle["h_dca_xy"])
            + _collect_rresult_ptrs(bundle["h_dca_z"])
            + _collect_rresult_ptrs(bundle["h_dca_xy_secondary"])
            + _collect_rresult_ptrs(bundle["h_tof_mass_vs_tpc_nsigma"])
        )
    _run_graphs(all_actions)

    if draw:
        for particle in particles:
            bundle = bundles[particle]
            for hist_map, labels in (
                (bundle["h_tpc"], ("A", "M")),
                (bundle["h_tof"], ("A", "M")),
            ):
                for label in labels:
                    c = ROOT.TCanvas()
                    c.SetRightMargin(0.15)
                    hist_map[label][0].DrawClone("colz")
    for particle in particles:
        _write_data_bundle(bundles[particle], outputs[particle])
    LOGGER.info("analyse_data_multi done outputs=%s", {k: outputs[k] for k in particles})


def _book_mc_species(df_data: Any, particle: str, enable_trials: bool, tag: str = "") -> dict[str, Any]:
    if particle == "he4":
        df = (
            df_data.Define("gP", "fgPt * std::cosh(fgEta)")
            .Define("gM", "std::abs(fPDGcode) == 1000020030 ? 2.809230089 : (std::abs(fPDGcode) == 1000010030 ? 2.80892 : (std::abs(fPDGcode) == 1000020040 ? 3.72738 : (std::abs(fPDGcode) == 1000010020 ? 1.87561 : 0.1)))")
            .Define("gMt", "std::hypot(gM, fgPt)")
            .Define("yMC", "std::asinh(fgPt / gMt * std::sinh(fgEta))")
            .Define("deltaPtUncorrected", "ptUncorr - fgPt")
            .Define("deltaPt", "ptHe4 - fgPt")
            .Define("ptWeight", "(5.04194/1.3645054) * fgPt * std::exp(-fgPt * 1.35934)")
            .Define("rapidity", "std::asinh(pt / std::hypot(pt, gM) * std::sinh(fEta))")
            .Define("isHe4", "std::abs(fPDGcode) == 1000020040")
            .Filter(HE4_MC_PID_SELECTION)
        )
        df_cut_reco = df.Filter(HE4_MC_RECO_SELECTION)
        df_cut_gen = df.Filter(HE4_MC_GEN_SELECTION)
        h_delta_pt = df_cut_reco.Histo2D(ROOT.RDF.TH2DModel(f"hDeltaPtHe4{tag}", ";#it{p}_{T}^{rec} (GeV/#it{c});#it{p}_{T}^{rec}-#it{p}_{T}^{gen} (GeV/#it{c})", 100, 0, 5, 120, -0.4, 0.2), "ptUncorr", "deltaPtUncorrected")
        h_delta_pt_corr = df_cut_reco.Histo2D(ROOT.RDF.TH2DModel(f"hDeltaPtCorrHe4{tag}", ";#it{p}_{T}^{rec} (GeV/#it{c});#it{p}_{T}^{rec}-#it{p}_{T}^{gen} (GeV/#it{c})", 100, 0, 5, 100, -0.4, 0.2), "pt", "deltaPt")
        h_mom_res = df_cut_reco.Histo2D(ROOT.RDF.TH2DModel(f"hMomResHe4{tag}", ";#it{p}_{T}^{rec} (GeV/#it{c});#it{p}_{T}^{rec}-#it{p}_{T}^{gen} (GeV/#it{c})", 44, 0.9, 5.3, 80, -0.2, 0.2), "pt", "deltaPt")

        h_reco_tpc_a = [df_cut_reco.Filter(f"!matter && {HE4_MC_SIGNAL_TRACKING}").Histo1D(h1_model(f"TPCAHe4{tag}", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts"), "pt")]
        h_reco_tpc_m = [df_cut_reco.Filter(f"matter && {HE4_MC_SIGNAL_TRACKING}").Histo1D(h1_model(f"TPCMHe4{tag}", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts"), "pt")]
        h_reco_tof_a = [df_cut_reco.Filter(f"!matter && {HE4_MC_SIGNAL_TRACKING} && hasTOF").Histo1D(h1_model(f"TOFAHe4{tag}", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts"), "pt")]
        h_reco_tof_m = [df_cut_reco.Filter(f"matter && {HE4_MC_SIGNAL_TRACKING} && hasTOF").Histo1D(h1_model(f"TOFMHe4{tag}", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts"), "pt")]
        h_gen_a = [df_cut_gen.Filter("fPDGcode < 0").Histo1D(h1_model(f"genAHe4{tag}", ";#it{p}_{T}^{gen} (GeV/#it{c});Counts"), "fgPt")]
        h_gen_m = [df_cut_gen.Filter("fPDGcode > 0").Histo1D(h1_model(f"genMHe4{tag}", ";#it{p}_{T}^{gen} (GeV/#it{c});Counts"), "fgPt")]

        h_reco_tpc_a_w = [df_cut_reco.Filter(f"!matter && {HE4_MC_SIGNAL_TRACKING}").Histo1D(h1_model(f"TPCAHe4W{tag}", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts"), "pt", "ptWeight")]
        h_reco_tpc_m_w = [df_cut_reco.Filter(f"matter && {HE4_MC_SIGNAL_TRACKING}").Histo1D(h1_model(f"TPCMHe4W{tag}", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts"), "pt", "ptWeight")]
        h_reco_tof_a_w = [df_cut_reco.Filter(f"!matter && {HE4_MC_SIGNAL_TRACKING} && hasTOF").Histo1D(h1_model(f"TOFAHe4W{tag}", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts"), "pt", "ptWeight")]
        h_reco_tof_m_w = [df_cut_reco.Filter(f"matter && {HE4_MC_SIGNAL_TRACKING} && hasTOF").Histo1D(h1_model(f"TOFMHe4W{tag}", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts"), "pt", "ptWeight")]
        h_gen_a_w = [df_cut_gen.Filter("fPDGcode < 0").Histo1D(h1_model(f"genAHe4W{tag}", ";#it{p}_{T}^{gen} (GeV/#it{c});Counts"), "fgPt", "ptWeight")]
        h_gen_m_w = [df_cut_gen.Filter("fPDGcode > 0").Histo1D(h1_model(f"genMHe4W{tag}", ";#it{p}_{T}^{gen} (GeV/#it{c});Counts"), "fgPt", "ptWeight")]

        reco_df_for_trials = df_cut_reco
        tpc_name = "He4"
        gen_name = "He4"
    else:
        df = (
            df_data.Define("gP", "fgPt * std::cosh(fgEta)")
            .Define("gM", "std::abs(fPDGcode) == 1000020030 ? 2.809230089 : (std::abs(fPDGcode) == 1000010030 ? 2.80892 : (std::abs(fPDGcode) == 1000020040 ? 3.72738 : (std::abs(fPDGcode) == 1000010020 ? 1.87561 : 0.1)))")
            .Define("gMt", "std::hypot(gM, fgPt)")
            .Define("yMC", "std::asinh(fgPt / gMt * std::sinh(fgEta))")
            .Define("deltaPtUncorrected", "ptUncorr - fgPt")
            .Define("deltaPt", "pt - fgPt")
            .Define("ptWeight", "(5.04194/1.3645054) * fgPt * std::exp(-fgPt * 1.35934)")
            .Define("isHe3", "std::abs(fPDGcode) == 1000020030")
            .Filter("isHe3")
        )
        df_cut_reco_base = df.Filter(BASE_REC_SELECTIONS + HE3_MC_RECO_APPEND)
        df_cut_reco = df_cut_reco_base.Filter(DEFAULT_REC_SELECTIONS)
        df_cut_gen = df.Filter(HE3_MC_GEN_SELECTION)
        h_delta_pt = df_cut_reco.Histo2D(ROOT.RDF.TH2DModel(f"hDeltaPtHe3{tag}", ";#it{p}_{T}^{rec} (GeV/#it{c});#it{p}_{T}^{rec}-#it{p}_{T}^{gen} (GeV/#it{c})", 100, 0, 5, 120, -0.4, 0.2), "pt", "deltaPtUncorrected")
        h_delta_pt_corr = df_cut_reco.Histo2D(ROOT.RDF.TH2DModel(f"hDeltaPtCorrHe3{tag}", ";#it{p}_{T}^{rec} (GeV/#it{c});#it{p}_{T}^{rec}-#it{p}_{T}^{gen} (GeV/#it{c})", 100, 0, 5, 100, -0.4, 0.2), "pt", "deltaPt")
        h_mom_res = df_cut_reco.Histo2D(ROOT.RDF.TH2DModel(f"hMomResHe3{tag}", ";#it{p}_{T}^{rec} (GeV/#it{c});#it{p}_{T}^{rec}-#it{p}_{T}^{gen} (GeV/#it{c})", 44, 0.9, 5.3, 80, -0.2, 0.2), "pt", "deltaPt")

        h_reco_tpc_a = [df_cut_reco.Filter("!matter").Histo1D(h1_model(f"TPCAHe3{tag}", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts"), "pt")]
        h_reco_tpc_m = [df_cut_reco.Filter("matter").Histo1D(h1_model(f"TPCMHe3{tag}", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts"), "pt")]
        h_reco_tof_a = [df_cut_reco.Filter("!matter && hasTOF").Histo1D(h1_model(f"TOFAHe3{tag}", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts"), "pt")]
        h_reco_tof_m = [df_cut_reco.Filter("matter && hasTOF").Histo1D(h1_model(f"TOFMHe3{tag}", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts"), "pt")]
        h_gen_a = [df_cut_gen.Filter("fPDGcode < 0").Histo1D(h1_model(f"genAHe3{tag}", ";#it{p}_{T}^{gen} (GeV/#it{c});Counts"), "fgPt")]
        h_gen_m = [df_cut_gen.Filter("fPDGcode > 0").Histo1D(h1_model(f"genMHe3{tag}", ";#it{p}_{T}^{gen} (GeV/#it{c});Counts"), "fgPt")]

        h_reco_tpc_a_w = [df_cut_reco.Filter("!matter").Histo1D(h1_model(f"TPCAHe3W{tag}", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts"), "pt", "ptWeight")]
        h_reco_tpc_m_w = [df_cut_reco.Filter("matter").Histo1D(h1_model(f"TPCMHe3W{tag}", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts"), "pt", "ptWeight")]
        h_reco_tof_a_w = [df_cut_reco.Filter("!matter && hasTOF").Histo1D(h1_model(f"TOFAHe3W{tag}", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts"), "pt", "ptWeight")]
        h_reco_tof_m_w = [df_cut_reco.Filter("matter && hasTOF").Histo1D(h1_model(f"TOFMHe3W{tag}", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts"), "pt", "ptWeight")]
        h_gen_a_w = [df_cut_gen.Filter("fPDGcode < 0").Histo1D(h1_model(f"genAHe3W{tag}", ";#it{p}_{T}^{gen} (GeV/#it{c});Counts"), "fgPt", "ptWeight")]
        h_gen_m_w = [df_cut_gen.Filter("fPDGcode > 0").Histo1D(h1_model(f"genMHe3W{tag}", ";#it{p}_{T}^{gen} (GeV/#it{c});Counts"), "fgPt", "ptWeight")]

        reco_df_for_trials = df_cut_reco_base
        tpc_name = "He3"
        gen_name = "He3"

    matter_map = {"A": "!matter", "M": "matter"}
    h_gen = {"A": h_gen_a, "M": h_gen_m}
    h_reco_tpc = {"A": h_reco_tpc_a, "M": h_reco_tpc_m}
    h_reco_tof = {"A": h_reco_tof_a, "M": h_reco_tof_m}
    h_gen_w = {"A": h_gen_a_w, "M": h_gen_m_w}
    h_reco_tpc_w = {"A": h_reco_tpc_a_w, "M": h_reco_tpc_m_w}
    h_reco_tof_w = {"A": h_reco_tof_a_w, "M": h_reco_tof_m_w}

    n_trials = 0
    if enable_trials:
        for cut_dca_z in CUT_NAMES["nsigmaDCAz"]:
            d_dca = reco_df_for_trials.Filter(f"std::abs(nsigmaDCAz) < {cut_dca_z}")
            for cut_tpc in CUT_NAMES["fTPCnCls"]:
                d_tpc = d_dca.Filter(f"fTPCnCls > {cut_tpc}")
                for cut_its in CUT_NAMES["nITScls"]:
                    _ = d_tpc.Filter(f"nITScls >= {cut_its}")
                    for label, matter_sel in matter_map.items():
                        h_reco_tpc[label].append(reco_df_for_trials.Filter(matter_sel).Histo1D(h1_model(f"TPC{label}{tpc_name}{n_trials}{tag}", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts"), "pt"))
                        h_reco_tof[label].append(reco_df_for_trials.Filter(f"{matter_sel} && hasTOF").Histo1D(h1_model(f"TOF{label}{tpc_name}{n_trials}{tag}", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts"), "pt"))
                        h_reco_tpc_w[label].append(reco_df_for_trials.Filter(matter_sel).Histo1D(h1_model(f"TPC{label}{tpc_name}W{n_trials}{tag}", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts"), "pt", "ptWeight"))
                        h_reco_tof_w[label].append(reco_df_for_trials.Filter(f"{matter_sel} && hasTOF").Histo1D(h1_model(f"TOF{label}{tpc_name}W{n_trials}{tag}", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts"), "pt", "ptWeight"))
                    n_trials += 1

    actions = (
        _collect_rresult_ptrs(h_gen)
        + _collect_rresult_ptrs(h_reco_tpc)
        + _collect_rresult_ptrs(h_reco_tof)
        + _collect_rresult_ptrs(h_gen_w)
        + _collect_rresult_ptrs(h_reco_tpc_w)
        + _collect_rresult_ptrs(h_reco_tof_w)
        + _collect_rresult_ptrs([h_delta_pt, h_delta_pt_corr, h_mom_res])
    )

    return {
        "particle": particle,
        "tpc_name": tpc_name,
        "gen_name": gen_name,
        "h_gen": h_gen,
        "h_reco_tpc": h_reco_tpc,
        "h_reco_tof": h_reco_tof,
        "h_gen_w": h_gen_w,
        "h_reco_tpc_w": h_reco_tpc_w,
        "h_reco_tof_w": h_reco_tof_w,
        "h_delta_pt": h_delta_pt,
        "h_delta_pt_corr": h_delta_pt_corr,
        "h_mom_res": h_mom_res,
        "n_trials": n_trials,
        "actions": actions,
    }


def _draw_mc_eff(bundle: dict[str, Any], canvas_name: str) -> None:
    c = ROOT.TCanvas(canvas_name, canvas_name)
    c.cd()
    h_eff = {}
    for label in ("A", "M"):
        h_eff[label] = bundle["h_reco_tpc"][label][0].GetValue().Clone(f"{canvas_name}_hEff{label}")
        h_eff[label].Divide(bundle["h_reco_tpc"][label][0].GetPtr(), bundle["h_gen"][label][0].GetPtr(), 1.0, 1.0, "B")
    h_eff["A"].SetLineColor(ROOT.kRed)
    h_eff["A"].Draw()
    h_eff["M"].Draw("same")


def _write_mc_bundle(bundle: dict[str, Any], output_file: str) -> None:
    out_name = expand(output_file)
    ensure_parent(out_name)
    out = ROOT.TFile(out_name, "recreate")
    d0 = out.mkdir("nuclei")
    d0.cd()

    tpc_name = bundle["tpc_name"]
    gen_name = bundle["gen_name"]

    for label in ("A", "M"):
        write_hist(bundle["h_gen"][label][0], f"gen{label}{gen_name}")
        write_hist(bundle["h_reco_tpc"][label][0], f"TPC{label}{tpc_name}")
        write_hist(bundle["h_reco_tof"][label][0], f"TOF{label}{tpc_name}")
    write_hist(bundle["h_delta_pt"], f"hDeltaPt{tpc_name}")
    write_hist(bundle["h_delta_pt_corr"], f"hDeltaPtCorr{tpc_name}")
    write_hist(bundle["h_mom_res"], f"hMomRes{tpc_name}")

    for label in ("A", "M"):
        write_hist(bundle["h_gen_w"][label][0], f"gen{label}{gen_name}W")
        write_hist(bundle["h_reco_tpc_w"][label][0], f"TPC{label}{tpc_name}W")
        write_hist(bundle["h_reco_tof_w"][label][0], f"TOF{label}{tpc_name}W")

    def write_eff(prefix: str, reco_tpc: dict[str, Any], reco_tof: dict[str, Any], gen: dict[str, Any], index: int) -> dict[str, Any]:
        eff_tpc = {}
        for label in ("A", "M"):
            eff_tpc[label] = reco_tpc[label][index].GetValue().Clone(f"{prefix}TPC{label}")
            eff_tof = reco_tof[label][index].GetValue().Clone(f"{prefix}TOF{label}")
            eff_tpc[label].Divide(reco_tpc[label][index].GetPtr(), gen[label][0].GetPtr(), 1.0, 1.0, "B")
            eff_tof.Divide(reco_tof[label][index].GetPtr(), gen[label][0].GetPtr(), 1.0, 1.0, "B")
            eff_tpc[label].GetYaxis().SetTitle("Efficiency #times Acceptance")
            eff_tof.GetYaxis().SetTitle("Efficiency #times Acceptance")
            eff_tpc[label].Write(f"{prefix}effTPC{label}")
            eff_tof.Write(f"{prefix}effTOF{label}")
        return eff_tpc

    write_eff("", bundle["h_reco_tpc"], bundle["h_reco_tof"], bundle["h_gen"], 0)
    write_eff(_weighted_eff_name(""), bundle["h_reco_tpc_w"], bundle["h_reco_tof_w"], bundle["h_gen_w"], 0)

    for i in range(bundle["n_trials"]):
        dtrial = out.mkdir(f"nuclei{i}")
        dtrial.cd()
        for label in ("A", "M"):
            write_hist(bundle["h_gen"][label][0], f"gen{label}{gen_name}")
            write_hist(bundle["h_reco_tpc"][label][i + 1], f"TPC{label}{tpc_name}")
            write_hist(bundle["h_reco_tof"][label][i + 1], f"TOF{label}{tpc_name}")
            write_hist(bundle["h_gen_w"][label][0], f"gen{label}{gen_name}W")
            write_hist(bundle["h_reco_tpc_w"][label][i + 1], f"TPC{label}{tpc_name}W")
            write_hist(bundle["h_reco_tof_w"][label][i + 1], f"TOF{label}{tpc_name}W")

        eff_tpc = write_eff("", bundle["h_reco_tpc"], bundle["h_reco_tof"], bundle["h_gen"], i + 1)
        for label in ("A", "M"):
            matching_tof = bundle["h_reco_tof"][label][i + 1].GetValue().Clone(f"matchingTOF{label}{i}")
            matching_tof.Divide(eff_tpc[label])
            matching_tof.Write()

        eff_w_tpc = write_eff(_weighted_eff_name(""), bundle["h_reco_tpc_w"], bundle["h_reco_tof_w"], bundle["h_gen_w"], i + 1)
        for label in ("A", "M"):
            matching_w_tof = bundle["h_reco_tof_w"][label][i + 1].GetValue().Clone(f"matchingWTOF{label}{i}")
            matching_w_tof.Divide(eff_w_tpc[label])
            matching_w_tof.Write()
    out.Close()


def analyse_mc(input_file: str, output_file: str, particle: str, enable_trials: bool, draw: bool = False) -> None:
    LOGGER.info("analyse_mc start particle=%s input=%s output=%s", particle, input_file, output_file)
    ROOT.gStyle.SetOptStat(0)
    rdf = build_rdf_from_ao2d("O2nucleitablemc", input_file)
    df_data = define_columns_for_data(rdf)
    bundle = _book_mc_species(df_data, particle, enable_trials, tag=f"_{particle}")
    _run_graphs(bundle["actions"])
    if draw:
        _draw_mc_eff(bundle, f"effMatterAntiMatter_{particle}")
    _write_mc_bundle(bundle, output_file)
    LOGGER.info("analyse_mc done output=%s", output_file)


def analyse_mc_multi(input_file: str, outputs: dict[str, str], enable_trials: bool, draw: bool = False) -> None:
    particles = [p for p in ("he3", "he4") if p in outputs]
    if not particles:
        raise ValueError("analyse_mc_multi requires at least one output among {'he3','he4'}.")
    LOGGER.info("analyse_mc_multi start input=%s particles=%s", input_file, particles)
    ROOT.gStyle.SetOptStat(0)
    rdf = build_rdf_from_ao2d("O2nucleitablemc", input_file)
    df_data = define_columns_for_data(rdf)

    bundles: dict[str, dict[str, Any]] = {}
    all_actions: list[Any] = []
    for particle in particles:
        bundle = _book_mc_species(df_data, particle, enable_trials, tag=f"_{particle}")
        bundles[particle] = bundle
        all_actions += bundle["actions"]
    _run_graphs(all_actions)

    if draw:
        for particle in particles:
            _draw_mc_eff(bundles[particle], f"effMatterAntiMatter_{particle}")
    for particle in particles:
        _write_mc_bundle(bundles[particle], outputs[particle])
    LOGGER.info("analyse_mc_multi done outputs=%s", {k: outputs[k] for k in particles})


def signal(input_file: str, output_file: str) -> None:
    LOGGER.info("signal start input=%s output=%s", input_file, output_file)
    load_fit_modules()
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)
    ROOT.RooMsgService.instance().setSilentMode(True)
    ROOT.gErrorIgnoreLevel = ROOT.kError

    in_file = ROOT.TFile(expand(input_file))
    out_name = expand(output_file)
    ensure_parent(out_name)
    out_file = ROOT.TFile(out_name, "recreate")

    m = ROOT.RooRealVar("dm2", "m - m_{^{3}He}", -1.2, 1.5, "GeV/#it{c}^{2}")
    m.setBins(1000, "cache")
    m.setRange("Full", -1.2, 1.5)

    f_sig = ROOT.FitExpTailGaus(m)
    ptr(f_sig.mMu).setRange(-1.0, 1.0)
    ptr(f_sig.mMu).setVal(0.1)
    ptr(f_sig.mSigma).setRange(0.05, 0.40)
    ptr(f_sig.mSigma).setVal(0.1)
    ptr(f_sig.mAlpha0).setRange(0.8, 3.0)
    ptr(f_sig.mAlpha0).setVal(1.2)
    ptr(f_sig.mSigCounts).setRange(0.0, 5000.0)

    m_bis = ROOT.RooRealVar("dm2_bis", "m - m_{^{3}He}", -1.2, 1.5, "GeV/#it{c}^{2}")
    m_bis.setBins(1000, "cache")
    m_bis.setRange("Full", -1.2, 1.5)
    f_bkg = ROOT.FitExpExpTailGaus(m_bis)
    f_bkg.UseSignal(False)

    ns = ROOT.RooRealVar("ns", "n#sigma_{^{3}He}", -5.0, 5.0, "a. u.")
    ns.setBins(1000, "cache")
    ns.setRange("Full", -5.0, 5.0)
    ns.setRange("Special", -4.0, 5.0)

    f_gg = ROOT.FitGausGaus(ns)
    ptr(f_gg.mSigma).setRange(0.2, 1.2)
    ptr(f_gg.mSigma).setVal(1.0)
    ptr(f_gg.mMu).setRange(-0.5, 0.5)
    ptr(f_gg.mMuBkg).setRange(-10.0, -4.0)
    ptr(f_gg.mMuBkg).setVal(-7.0)
    ptr(f_gg.mSigmaBkg).setRange(0.2, 6.0)

    f_eg = ROOT.FitExpGaus(ns)
    ptr(f_eg.mSigma).setRange(0.2, 1.2)
    ptr(f_eg.mSigma).setVal(1.0)
    ptr(f_eg.mMu).setRange(-0.5, 0.5)

    f_etg = ROOT.FitExpTailGaus(ns)
    ptr(f_etg.mSigma).setRange(0.2, 1.2)
    ptr(f_etg.mSigma).setVal(1.0)
    ptr(f_etg.mMu).setRange(-0.5, 0.5)

    f_lnl = ROOT.FitLogNormalLogNormal(ns)
    ptr(f_lnl.mSigma).setRange(1.01, 20.0)
    ptr(f_lnl.mSigma).setVal(ROOT.TMath.Exp(1.0))
    ptr(f_lnl.mMu).setRange(-0.5, 0.5)

    tpc_functions = [f_gg, f_eg, f_etg, f_lnl]

    for key in in_file.GetListOfKeys():
        key_name = key.GetName()
        if FILTER_LIST_NAME not in key_name:
            continue
        d_in = in_file.Get(key_name)
        d_out = out_file.mkdir(key_name)
        out_file.cd(key_name)

        f_a_tof = d_in.Get("fATOFsignal")
        f_m_tof = d_in.Get("fMTOFsignal")
        f_a_tpc = d_in.Get("fATPCcounts")
        f_m_tpc = d_in.Get("fMTPCcounts")

        n_pt_bins = f_a_tof.GetNbinsX()
        pt_axis = f_a_tof.GetXaxis()
        pt_labels = pt_axis.GetXbins()

        tof_h = [f_m_tof, f_a_tof]
        tpc_h = [f_m_tpc, f_a_tpc]

        h_raw = [[None] * CENT_LENGTH for _ in range(2)]
        h_raw_bc = [[None] * CENT_LENGTH for _ in range(2)]
        h_sig = [[None] * CENT_LENGTH for _ in range(2)]
        h_signif = [[None] * CENT_LENGTH for _ in range(2)]
        h_chi = [[None] * CENT_LENGTH for _ in range(2)]
        h_chi_tpc = [[None] * CENT_LENGTH for _ in range(2)]
        h_tpconly = [[[None] * NTPC_FUNCTIONS for _ in range(CENT_LENGTH)] for _ in range(2)]
        h_widen = [[None] * CENT_LENGTH for _ in range(2)]
        h_shift = [[None] * CENT_LENGTH for _ in range(2)]
        h_widen_tpc = [[None] * CENT_LENGTH for _ in range(2)]
        h_shift_tpc = [[None] * CENT_LENGTH for _ in range(2)]

        n_sigma_vec = [3.0]
        v_shift = []

        for i_s in range(2):
            d_species = d_out.mkdir(NAMES[i_s])
            d_species.cd()
            d_sig = d_species.mkdir("GausExp")
            for i_c in range(CENT_LENGTH):
                d_sig.mkdir(f"C_{i_c}")
            d_side = d_species.mkdir("Sidebands")
            for i_c in range(CENT_LENGTH):
                d_side.mkdir(f"C_{i_c}")
            d_species.mkdir("Significance")
            d_species.mkdir("Systematic")
            d_species.mkdir("TPConly")
            d_species.mkdir("ChiSquare")

        for i_s in range(2):
            i_c = 0
            for i_t in range(NTPC_FUNCTIONS):
                h_tpconly[i_s][i_c][i_t] = ROOT.TH1D(f"hTPConly{LETTER[i_s]}{i_c}_{TPC_FUNCTION_NAMES[i_t]}", ";p_{T} GeV/c; TPC raw counts", n_pt_bins, pt_labels.GetArray())
            h_signif[i_s][i_c] = ROOT.TH1D(f"hSignificance{LETTER[i_s]}{i_c}", "; p_{T}(GeV/c); #frac{S}{#sqrt{S+B}}", n_pt_bins, pt_labels.GetArray())
            h_chi[i_s][i_c] = ROOT.TH1D(f"hChiSquare{LETTER[i_s]}{i_c}", "; p_{T}(GeV/c); #chi^{2}/NDF", n_pt_bins, pt_labels.GetArray())
            h_chi_tpc[i_s][i_c] = ROOT.TH1D(f"hChiSquareTPC{LETTER[i_s]}{i_c}", "; p_{T}(GeV/c); #chi^{2}/NDF", n_pt_bins, pt_labels.GetArray())
            h_raw[i_s][i_c] = ROOT.TH1D(f"hRawCounts{LETTER[i_s]}{i_c}", "; p_{T}(GeV/c); RawCounts", n_pt_bins, pt_labels.GetArray())
            h_raw_bc[i_s][i_c] = ROOT.TH1D(f"hRawCountsBinCounting{LETTER[i_s]}{i_c}", "; p_{T}(GeV/c); RawCounts", n_pt_bins, pt_labels.GetArray())
            h_sig[i_s][i_c] = ROOT.TH1D(f"hSignalGausExpGaus{LETTER[i_s]}{i_c}", "; p_{T}(GeV/c); RawCounts", n_pt_bins, pt_labels.GetArray())
            h_widen[i_s][i_c] = ROOT.TH1D(f"hWidenRangeSyst{LETTER[i_s]}{i_c}", "; p_{T}(GeV/c); RMS", n_pt_bins, pt_labels.GetArray())
            h_shift[i_s][i_c] = ROOT.TH1D(f"hShiftRangeSyst{LETTER[i_s]}{i_c}", "; p_{T}(GeV/c); RMS", n_pt_bins, pt_labels.GetArray())
            h_widen_tpc[i_s][i_c] = ROOT.TH1D(f"hWidenRangeSystTPC{LETTER[i_s]}{i_c}", "; p_{T}(GeV/c); RMS", n_pt_bins, pt_labels.GetArray())
            h_shift_tpc[i_s][i_c] = ROOT.TH1D(f"hShiftRangeSystTPC{LETTER[i_s]}{i_c}", "; p_{T}(GeV/c); RMS", n_pt_bins, pt_labels.GetArray())

        for i_b in range(n_pt_bins):
            center_pt = pt_axis.GetBinCenter(i_b + 1)
            if center_pt < PT_RANGE[0] or center_pt > PT_RANGE[1]:
                continue

            for i_s in range(2):
                i_c = 0
                if center_pt > CENT_PT_LIMITS[i_c]:
                    continue

                i_title = f"{pt_labels[i_b]:1.1f} #leq #it{{p}}_{{T}} < {pt_labels[i_b + 1]:1.1f} GeV/#it{{c}}"
                i_name = f"d{i_c}_{i_b}"
                dat = tof_h[i_s].ProjectionY(f"data{i_c}_{i_b}", i_b + 1, i_b + 1)

                ptr(f_sig.mTau0).setVal(-0.3)
                ptr(f_sig.mTau0).setVal(0.5)

                d_out.cd(f"{NAMES[i_s]}/GausExp/C_{i_c}")
                fit_plot = f_sig.FitData(dat, i_name, i_title, "Full", "Full", False, -1.2, 1.5)
                ptr(f_sig.mSigma).setConstant(False)
                if center_pt > TOF_MIN_PT:
                    fit_plot.Write()
                h_sig[i_s][i_c].SetBinContent(i_b + 1, ptr(f_sig.mSigCounts).getVal())
                h_sig[i_s][i_c].SetBinError(i_b + 1, ptr(f_sig.mSigCounts).getError())
                h_raw[i_s][i_c].SetBinContent(i_b + 1, ptr(f_sig.mSigCounts).getVal())
                h_raw[i_s][i_c].SetBinError(i_b + 1, ptr(f_sig.mSigCounts).getError())

                residual_vec = []
                for i_sigma, n_sigma in enumerate(n_sigma_vec):
                    left_sigma = ptr(f_sig.mMu).getVal() - n_sigma * ptr(f_sig.mSigma).getVal()
                    right_sigma = ptr(f_sig.mMu).getVal() + (n_sigma + 2.0) * ptr(f_sig.mSigma).getVal()
                    left_edge_bin = dat.FindBin(left_sigma)
                    left_edge_float = dat.GetBinLowEdge(left_edge_bin)
                    right_edge_bin = dat.FindBin(right_sigma)
                    right_edge_float = dat.GetBinLowEdge(right_edge_bin + 1)
                    f_bkg.mX.setRange("signal", left_edge_float, right_edge_float)

                    if i_sigma == 0:
                        f_bkg.mX.setRange("left", dat.GetXaxis().GetXmin(), left_edge_float)
                        f_bkg.mX.setRange("right", right_edge_float, dat.GetXaxis().GetXmax())
                        bkg_plot = f_bkg.FitData(dat, f"{i_name}_sideband", i_title, "left,right", "Full")
                        d_out.cd(f"{NAMES[i_s]}/Sidebands/C_{i_c}")
                        bkg_plot.Write()

                    bkg_integral = 0.0
                    if i_b > 8:
                        bkg_integral = f_bkg.mBackground.createIntegral(m_bis, ROOT.RooFit.NormSet(m_bis), ROOT.RooFit.Range("signal")).getVal() * ptr(f_bkg.mBkgCounts).getVal()
                        h_chi[i_s][i_c].SetBinContent(i_b + 1, f_bkg.mChi2)

                    tot_integral = dat.Integral(left_edge_bin, right_edge_bin)
                    sig_integral = tot_integral - bkg_integral
                    sig_err = ROOT.TMath.Sqrt(tot_integral + bkg_integral)

                    if i_sigma == 0:
                        h_raw_bc[i_s][i_c].SetBinContent(i_b + 1, sig_integral)
                        h_raw_bc[i_s][i_c].SetBinError(i_b + 1, sig_err)
                        if tot_integral > 0:
                            h_signif[i_s][i_c].SetBinContent(i_b + 1, sig_integral / ROOT.TMath.Sqrt(tot_integral))

                    residual_vec.append(sig_integral)

                raw_val = h_raw[i_s][i_c].GetBinContent(i_b + 1)
                if raw_val > 0:
                    h_widen[i_s][i_c].SetBinContent(i_b + 1, ROOT.TMath.RMS(len(residual_vec), array("d", residual_vec)) / raw_val)

                shift_vec = []
                for shift in v_shift:
                    left_sigma = ptr(f_sig.mMu).getVal() - 3.0 * ptr(f_sig.mSigma).getVal() - shift
                    right_sigma = ptr(f_sig.mMu).getVal() + 5.0 * ptr(f_sig.mSigma).getVal() - shift
                    left_edge_bin = dat.FindBin(left_sigma)
                    right_edge_bin = dat.FindBin(right_sigma)
                    left_edge_float = dat.GetBinLowEdge(left_edge_bin)
                    right_edge_float = dat.GetBinLowEdge(right_edge_bin + 1)
                    f_bkg.mX.setRange("signal", left_edge_float, right_edge_float)
                    bkg_integral = 0.0
                    if i_b > 7:
                        bkg_integral = f_bkg.mBackground.createIntegral(m_bis, ROOT.RooFit.NormSet(m_bis), ROOT.RooFit.Range("signal")).getVal() * ptr(f_bkg.mBkgCounts).getVal()
                    tot_integral = dat.Integral(left_edge_bin, right_edge_bin)
                    shift_vec.append(tot_integral - bkg_integral)
                if shift_vec and raw_val > 0:
                    h_shift[i_s][i_c].SetBinContent(i_b + 1, ROOT.TMath.RMS(len(shift_vec), array("d", shift_vec)) / raw_val)

                if center_pt < TPC_MAX_PT:
                    d_out.cd(f"{NAMES[i_s]}/TPConly")
                    tpc_dat = tpc_h[i_s].ProjectionY(f"tpc_data{i_c}_{i_b}", i_b + 1, i_b + 1)
                    for i_t in range(NTPC_FUNCTIONS):
                        fit_range = "Special" if (i_t and center_pt < 1.8) else "Full"
                        tpc_plot = tpc_functions[i_t].FitData(tpc_dat, f"TPC_d_{i_c}_{i_b}_{TPC_FUNCTION_NAMES[i_t]}", i_title, fit_range, fit_range)
                        ptr(tpc_functions[i_t].mSigma).setConstant(False)
                        tpc_plot.Write()
                        h_tpconly[i_s][i_c][i_t].SetBinContent(i_b + 1, ptr(tpc_functions[i_t].mSigCounts).getVal())
                        h_tpconly[i_s][i_c][i_t].SetBinError(i_b + 1, ptr(tpc_functions[i_t].mSigCounts).getError())

        for i_s in range(2):
            i_c = 0
            d_out.cd(f"{NAMES[i_s]}/GausExp")
            h_raw[i_s][i_c].Write()
            h_raw_bc[i_s][i_c].Write()
            h_sig[i_s][i_c].Write()

            d_out.cd(f"{NAMES[i_s]}/Systematic")
            h_shift[i_s][i_c].Write()
            h_widen[i_s][i_c].Write()
            h_widen_tpc[i_s][i_c].Write()
            h_shift_tpc[i_s][i_c].Write()

            d_out.cd(f"{NAMES[i_s]}/Significance")
            h_signif[i_s][i_c].Write()

            d_out.cd(f"{NAMES[i_s]}/TPConly")
            for i_t in range(NTPC_FUNCTIONS):
                h_tpconly[i_s][i_c][i_t].Write()

            d_out.cd(f"{NAMES[i_s]}/ChiSquare")
            h_chi[i_s][i_c].Write()
            h_chi_tpc[i_s][i_c].Write()

    out_file.Close()
    LOGGER.info("signal done output=%s", output_file)


def systematics(signal_file: str, mc_file: str, data_analysis_results: str, output_file: str) -> None:
    LOGGER.info("systematics start signal=%s mc=%s output=%s", signal_file, mc_file, output_file)
    f_data = ROOT.TFile(expand(signal_file))
    f_mc = ROOT.TFile(expand(mc_file))
    an_results = ROOT.TFile(expand(data_analysis_results))

    # Skimmed productions: prefer ZorroSummary normalization factor when available.
    norm = None
    zorro_summary = _find_object_by_class(an_results, "ZorroSummary")
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

    syst_tpc = [ROOT.TH2D(f"systTPC{NAMES[i]}", ";#it{p}_{T} (GeV/#it{c});Relative systematics TPC", N_PT_BINS, PT_BIN_ARRAY, 50, -0.5, 0.5) for i in range(2)]
    syst_tof = [ROOT.TH2D(f"systTOF{NAMES[i]}", ";#it{p}_{T} (GeV/#it{c});Relative systematics TOF", N_PT_BINS, PT_BIN_ARRAY, 50, -0.5, 0.5) for i in range(2)]

    default_eff_tpc = [None, None]
    default_eff_tof = [None, None]
    default_tpc = [None, None]
    default_tof = [None, None]
    default_tpc_uncorr = [None, None]
    default_tof_uncorr = [None, None]

    for key in f_data.GetListOfKeys():
        key_name = key.GetName()
        if FILTER_LIST_NAME not in key_name:
            continue

        list_data = f_data.Get(key_name)
        list_mc = f_mc.Get(key_name)

        h_data_tof = [[[None for _ in range(3)] for _ in range(2)] for _ in range(2)]
        h_data_tpc = [[None for _ in range(3)] for _ in range(2)]
        h_eff_tpc = [None, None]
        h_eff_tof = [None, None]

        tof_names = ["hRawCounts", "hRawCountsBinCounting"]
        tof_presel = ["", "_loose", "_tight"]

        for i_s in range(2):
            data_dir = list_data.Get(NAMES[i_s])
            h_eff_tpc[i_s] = list_mc.Get(f"effTPC{LETTER[i_s]}")
            h_eff_tof[i_s] = list_mc.Get(f"effTOF{LETTER[i_s]}")

            if default_eff_tpc[i_s] is None:
                default_eff_tpc[i_s] = h_eff_tpc[i_s].Clone(f"defaultEffTPC{NAMES[i_s]}")
            if default_eff_tof[i_s] is None:
                default_eff_tof[i_s] = h_eff_tof[i_s].Clone(f"defaultEffTOF{NAMES[i_s]}")

            for i_tof in range(2):
                name = tof_names[i_tof] + LETTER[i_s]
                h_data_tof[i_s][i_tof][0] = data_dir.Get(f"GausExp/{name}0{tof_presel[0]}")
                if default_tof_uncorr[i_s] is None:
                    default_tof_uncorr[i_s] = h_data_tof[i_s][i_tof][0].Clone(f"defaultTOFuncorr{NAMES[i_s]}")
                h_data_tof[i_s][i_tof][0].Divide(h_eff_tof[i_s])
                if default_tof[i_s] is None:
                    default_tof[i_s] = h_data_tof[i_s][i_tof][0].Clone(f"defaultTOF{NAMES[i_s]}")

            for i_tpc in range(3):
                name = TPC_FUNCTION_NAMES[i_tpc]
                h_data_tpc[i_s][i_tpc] = data_dir.Get(f"TPConly/hTPConly{LETTER[i_s]}0_{name}")
                if default_tpc_uncorr[i_s] is None and i_tpc == 1:
                    default_tpc_uncorr[i_s] = h_data_tpc[i_s][i_tpc].Clone(f"defaultTPCuncorr{NAMES[i_s]}")
                h_data_tpc[i_s][i_tpc].Divide(h_eff_tpc[i_s])
                if default_tpc[i_s] is None and i_tpc == 1:
                    default_tpc[i_s] = h_data_tpc[i_s][i_tpc].Clone(f"defaultTPC{NAMES[i_s]}")

        for i_s in range(2):
            for i_b in range(1, N_PT_BINS + 1):
                pt = h_data_tpc[i_s][0].GetBinCenter(i_b)
                default_val_tpc = default_tpc[i_s].GetBinContent(i_b)
                default_val_tof = default_tof[i_s].GetBinContent(i_b)
                if default_val_tpc != 0:
                    for i_tpc in range(3):
                        value = h_data_tpc[i_s][i_tpc].GetBinContent(i_b)
                        syst_tpc[i_s].Fill(pt, (value - default_val_tpc) / default_val_tpc)
                if default_val_tof != 0:
                    for i_tof in range(2):
                        value = h_data_tof[i_s][i_tof][0].GetBinContent(i_b)
                        syst_tof[i_s].Fill(pt, (value - default_val_tof) / default_val_tof)

    h_syst_tpc = [ROOT.TH1D(f"hSystTPC{LETTER[i]}", ";#it{p}_{T} (GeV/#it{c});Relative systematics TPC", N_PT_BINS, PT_BIN_ARRAY) for i in range(2)]
    h_syst_tof = [ROOT.TH1D(f"hSystTOF{LETTER[i]}", ";#it{p}_{T} (GeV/#it{c});Relative systematics TOF", N_PT_BINS, PT_BIN_ARRAY) for i in range(2)]

    for i_s in range(2):
        for i_b in range(1, N_PT_BINS + 1):
            h_syst_tpc[i_s].SetBinContent(i_b, syst_tpc[i_s].ProjectionY("", i_b, i_b).GetRMS())
            h_syst_tof[i_s].SetBinContent(i_b, syst_tof[i_s].ProjectionY("", i_b, i_b).GetRMS())

    out_name = expand(output_file)
    ensure_parent(out_name)
    out = ROOT.TFile(out_name, "recreate")

    tof_matching_m = default_tof_uncorr[0].Clone(f"TOFmatching{NAMES[0]}")
    tof_matching_a = default_tof_uncorr[1].Clone(f"TOFmatching{NAMES[1]}")
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
        f_stat_tpc = default_tpc_uncorr[i_s].Clone(f"fStatTPC{LETTER[i_s]}")
        f_syst_tpc = default_tpc_uncorr[i_s].Clone(f"fSystTPC{LETTER[i_s]}")
        f_stat_tof = default_tof_uncorr[i_s].Clone(f"fStatTOF{LETTER[i_s]}")
        f_syst_tof = default_tof_uncorr[i_s].Clone(f"fSystTOF{LETTER[i_s]}")

        for i_bin in range(1, N_PT_BINS + 1):
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
    LOGGER.info("systematics done output=%s", output_file)


def checkpoint(systematics_file: str, data_ar_file: str, mc_file: str, mc_ar_file: str, signal_file: str, output_file: str | None = None) -> None:
    LOGGER.info("checkpoint start output=%s", output_file if output_file else "auto")
    syst = ROOT.TFile(expand(systematics_file))
    data_ar = ROOT.TFile(expand(data_ar_file))
    mc = ROOT.TFile(expand(mc_file))
    mc_ar = ROOT.TFile(expand(mc_ar_file))
    sig = ROOT.TFile(expand(signal_file))

    if output_file is None:
        output_file = f"checkpoint-{datetime.now().strftime('%d%m%y')}.root"

    out_name = expand(output_file)
    ensure_parent(out_name)
    out = ROOT.TFile(out_name, "RECREATE")
    out.cd()

    syst.Get("pubStat").Clone("published_stat").Write()
    syst.Get("pubSyst").Clone("published_syst").Write()
    syst.Get("fStatTPCA").Clone("tpc_spectrum_stat").Write()
    syst.Get("fSystTPCA").Clone("tpc_spectrum_syst").Write()
    syst.Get("fStatTOFA").Clone("tof_spectrum_stat").Write()
    syst.Get("fSystTOFA").Clone("tof_spectrum_syst").Write()
    mc.Get("nuclei/effTPCA").Clone("tpc_efficiency").Write()
    mc.Get("nuclei/effTOFA").Clone("tof_efficiency").Write()

    out.mkdir("MC")
    out.cd("MC")
    mc.Get("nuclei/genAHe3").Clone("generated").Write()
    mc.Get("nuclei/TPCAHe3").Clone("tpc_reconstructed").Write()
    mc.Get("nuclei/TOFAHe3").Clone("tof_reconstructed").Write()
    mc_ar.Get("nuclei-spectra/spectra/hRecVtxZData").Clone("events_reconstructed").Write()

    out.mkdir("Data")
    out.cd("Data")
    data_ar.Get("nuclei-spectra/spectra/hRecVtxZData").Clone("events_reconstructed").Write()
    sig.Get("nuclei/antihe3/TPConly/hTPConlyA0_ExpGaus").Clone("tpc_rawcounts").Write()
    sig.Get("nuclei/antihe3/GausExp/hRawCountsA0").Clone("tof_rawcounts").Write()

    out.Close()
    LOGGER.info("checkpoint done output=%s", output_file if output_file else "auto")

