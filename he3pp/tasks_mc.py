import logging
from typing import Any

import ROOT

from .root_io import build_rdf_from_ao2d, define_columns_for_data, ensure_parent, expand, h1_model, write_hist
from .settings import *  # Loaded after runtime overrides in cli.run
from .tasks_common import collect_rresult_ptrs, run_graphs, weighted_eff_name


LOGGER = logging.getLogger("he3pp.tasks")


def _book_mc_species(df_data: Any, particle: str, enable_trials: bool, tag: str = "") -> dict[str, Any]:
    p = get_particle_config(particle)
    pdg_abs = int(p["pdg_abs"])
    mass = float(p.get("mc_mass", p["mass"]))
    tpc_name = str(p["mc_hist_suffix"])
    gen_name = str(p["mc_hist_suffix"])
    mc_signal_tracking = str(p.get("mc_signal_tracking", ""))

    df = (
        df_data.Define("gMt", f"std::hypot({mass}, fgPt)")
        .Define("yMC", "std::asinh(fgPt / gMt * std::sinh(fgEta))")
        .Define("deltaPtUncorrected", "ptUncorr - fgPt")
        .Define("deltaPt", f"{p['mc_pt_corr_col']} - fgPt")
        .Define("ptWeight", "(5.04194/1.3645054) * fgPt * std::exp(-fgPt * 1.35934)")
        .Define("isTarget", f"std::abs(fPDGcode) == {pdg_abs}")
        .Filter("isTarget")
    )

    reco_base = str(p.get("mc_reco_base_sel", "")).strip()
    df_cut_reco_base = df.Filter(reco_base) if reco_base else df
    df_cut_reco = df_cut_reco_base.Filter(str(p["mc_reco_sel"]))
    df_cut_gen = df.Filter(str(p["mc_gen_sel"]))
    h_delta_pt = df_cut_reco.Histo2D(
        ROOT.RDF.TH2DModel(f"hDeltaPt{tpc_name}{tag}", ";#it{p}_{T}^{rec} (GeV/#it{c});#it{p}_{T}^{rec}-#it{p}_{T}^{gen} (GeV/#it{c})", 100, 0, 5, 120, -0.4, 0.2),
        str(p["mc_delta_x_col"]),
        "deltaPtUncorrected",
    )
    h_delta_pt_corr = df_cut_reco.Histo2D(ROOT.RDF.TH2DModel(f"hDeltaPtCorr{tpc_name}{tag}", ";#it{p}_{T}^{rec} (GeV/#it{c});#it{p}_{T}^{rec}-#it{p}_{T}^{gen} (GeV/#it{c})", 100, 0, 5, 100, -0.4, 0.2), "pt", "deltaPt")
    h_mom_res = df_cut_reco.Histo2D(ROOT.RDF.TH2DModel(f"hMomRes{tpc_name}{tag}", ";#it{p}_{T}^{rec} (GeV/#it{c});#it{p}_{T}^{rec}-#it{p}_{T}^{gen} (GeV/#it{c})", 44, 0.9, 5.3, 80, -0.2, 0.2), "pt", "deltaPt")

    def _reco_filter(matter_expr: str, with_tof: bool = False) -> str:
        expr = matter_expr
        if mc_signal_tracking:
            expr = f"{expr} && {mc_signal_tracking}"
        if with_tof:
            expr = f"{expr} && hasTOF"
        return expr

    matter_defs = {
        "A": {"matter_sel": "!matter", "pdg_sel": "fPDGcode < 0"},
        "M": {"matter_sel": "matter", "pdg_sel": "fPDGcode > 0"},
    }
    h_reco_tpc: dict[str, list[Any]] = {}
    h_reco_tof: dict[str, list[Any]] = {}
    h_gen: dict[str, list[Any]] = {}
    h_reco_tpc_w: dict[str, list[Any]] = {}
    h_reco_tof_w: dict[str, list[Any]] = {}
    h_gen_w: dict[str, list[Any]] = {}
    for label, defs in matter_defs.items():
        matter_sel = defs["matter_sel"]
        pdg_sel = defs["pdg_sel"]
        h_reco_tpc[label] = [df_cut_reco.Filter(_reco_filter(matter_sel)).Histo1D(h1_model(f"TPC{label}{tpc_name}{tag}", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts"), "pt")]
        h_reco_tof[label] = [df_cut_reco.Filter(_reco_filter(matter_sel, with_tof=True)).Histo1D(h1_model(f"TOF{label}{tpc_name}{tag}", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts"), "pt")]
        h_gen[label] = [df_cut_gen.Filter(pdg_sel).Histo1D(h1_model(f"gen{label}{gen_name}{tag}", ";#it{p}_{T}^{gen} (GeV/#it{c});Counts"), "fgPt")]
        h_reco_tpc_w[label] = [df_cut_reco.Filter(_reco_filter(matter_sel)).Histo1D(h1_model(f"TPC{label}{tpc_name}W{tag}", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts"), "pt", "ptWeight")]
        h_reco_tof_w[label] = [df_cut_reco.Filter(_reco_filter(matter_sel, with_tof=True)).Histo1D(h1_model(f"TOF{label}{tpc_name}W{tag}", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts"), "pt", "ptWeight")]
        h_gen_w[label] = [df_cut_gen.Filter(pdg_sel).Histo1D(h1_model(f"gen{label}{gen_name}W{tag}", ";#it{p}_{T}^{gen} (GeV/#it{c});Counts"), "fgPt", "ptWeight")]

    reco_df_for_trials = df_cut_reco_base
    matter_map = {"A": "!matter", "M": "matter"}

    n_trials = 0
    if enable_trials and bool(p.get("trial_enabled", False)):
        for cut_dca_z in CUT_NAMES["nsigmaDCAz"]:
            d_dca = reco_df_for_trials.Filter(f"std::abs(nsigmaDCAz) < {cut_dca_z}")
            for cut_tpc in CUT_NAMES["fTPCnCls"]:
                d_tpc = d_dca.Filter(f"fTPCnCls > {cut_tpc}")
                for cut_its in CUT_NAMES["nITScls"]:
                    d_its = d_tpc.Filter(f"nITScls >= {cut_its}")
                    for label, matter_sel in matter_map.items():
                        h_reco_tpc[label].append(d_its.Filter(matter_sel).Histo1D(h1_model(f"TPC{label}{tpc_name}{n_trials}{tag}", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts"), "pt"))
                        h_reco_tof[label].append(d_its.Filter(f"{matter_sel} && hasTOF").Histo1D(h1_model(f"TOF{label}{tpc_name}{n_trials}{tag}", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts"), "pt"))
                        h_reco_tpc_w[label].append(d_its.Filter(matter_sel).Histo1D(h1_model(f"TPC{label}{tpc_name}W{n_trials}{tag}", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts"), "pt", "ptWeight"))
                        h_reco_tof_w[label].append(d_its.Filter(f"{matter_sel} && hasTOF").Histo1D(h1_model(f"TOF{label}{tpc_name}W{n_trials}{tag}", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts"), "pt", "ptWeight"))
                    n_trials += 1

    actions = (
        collect_rresult_ptrs(h_gen)
        + collect_rresult_ptrs(h_reco_tpc)
        + collect_rresult_ptrs(h_reco_tof)
        + collect_rresult_ptrs(h_gen_w)
        + collect_rresult_ptrs(h_reco_tpc_w)
        + collect_rresult_ptrs(h_reco_tof_w)
        + collect_rresult_ptrs([h_delta_pt, h_delta_pt_corr, h_mom_res])
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
    write_eff(weighted_eff_name(""), bundle["h_reco_tpc_w"], bundle["h_reco_tof_w"], bundle["h_gen_w"], 0)

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

        eff_w_tpc = write_eff(weighted_eff_name(""), bundle["h_reco_tpc_w"], bundle["h_reco_tof_w"], bundle["h_gen_w"], i + 1)
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
    run_graphs(bundle["actions"])
    if draw:
        _draw_mc_eff(bundle, f"effMatterAntiMatter_{particle}")
    _write_mc_bundle(bundle, output_file)
    LOGGER.info("analyse_mc done output=%s", output_file)
