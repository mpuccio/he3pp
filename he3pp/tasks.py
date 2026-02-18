from array import array
from datetime import datetime
import logging
from typing import Any

import ROOT

from .settings import *  # Loaded after runtime overrides in cli.run
from .root_io import define_columns_for_data, ensure_parent, expand, h1_model, h2_model, load_fit_modules, ptr, write_hist

LOGGER = logging.getLogger("he3pp.tasks")


def _weighted_eff_name(base: str) -> str:
    # Weighted efficiencies intentionally use "Weff*" naming.
    return f"W{base}"


def analyse_data(input_file: str, output_file: str, particle: str, skim: bool = False, draw: bool = False) -> None:
    LOGGER.info("analyse_data start particle=%s input=%s output=%s", particle, input_file, output_file)
    ROOT.gStyle.SetOptStat(0)
    rdf = ROOT.RDataFrame("O2nucleitable", expand(input_file))
    df_base = define_columns_for_data(rdf)

    if particle == "he4":
        df_base = df_base.Filter(HE4_BASE_SELECTION)
        df_primary = df_base.Filter(HE4_PRIMARY_SELECTION)
        df_secondary = df_base.Filter(SECONDARY_SELECTION)
        nsigma = "nsigmaHe4"
        mass_cut = "hasGoodTOFmassHe4"
        pt_name = "ptHe4"
        dmass = "deltaMassHe4"
        tpc_anti_model = ROOT.RDF.TH2DModel("fATPCcounts", ";#it{p}_{T}^{rec} (GeV/#it{c});^{4}#bar{He} n#sigma_{TPC};Counts", 160, 0.5, 4.5, 100, -5, 5)
        tpc_mat_model = h2_model("fMTPCcounts", ";#it{p}_{T}^{rec} (GeV/#it{c});^{4}He n#sigma_{TPC};Counts", 100, -5, 5)
    else:
        df_base = df_base.Filter(BASE_REC_SELECTIONS)
        df_primary = df_base.Filter(DEFAULT_REC_SELECTIONS)
        df_secondary = df_base.Filter(SECONDARY_SELECTION)
        nsigma = "nsigmaHe3"
        mass_cut = "hasGoodTOFmassHe3"
        pt_name = "pt"
        dmass = "deltaMassHe3"
        tpc_anti_model = h2_model("fATPCcounts", ";#it{p}_{T}^{rec} (GeV/#it{c});^{3}#bar{He} n#sigma_{TPC};Counts", 100, -5, 5)
        tpc_mat_model = h2_model("fMTPCcounts", ";#it{p}_{T}^{rec} (GeV/#it{c});^{3}He n#sigma_{TPC};Counts", 100, -5, 5)

    if skim:
        df_base.Filter(SKIM_SELECTION_TEMPLATE.format(nsigma=nsigma)).Snapshot("nucleiTree", "data/skimmed.root")

    h_dca_xy_a = [df_primary.Filter(f"!matter && {nsigma} > -0.5 && {nsigma} < 3 && {mass_cut}").Histo2D(h2_model("hDCAxyAHe", ";#it{p}_{T}^{rec} (GeV/#it{c});DCA_{xy} (cm);Counts", 100, -0.2, 0.2), pt_name, "fDCAxy")]
    h_dca_xy_m = [df_primary.Filter(f"matter && {nsigma} > -0.5 && {nsigma} < 3 && {mass_cut}").Histo2D(h2_model("hDCAxyMHe", ";#it{p}_{T}^{rec} (GeV/#it{c});DCA_{xy} (cm);Counts", 100, -0.2, 0.2), pt_name, "fDCAxy")]
    h_dca_z_a = [df_primary.Filter(f"!matter && {nsigma} > -0.5 && {nsigma} < 3 && {mass_cut}").Histo2D(h2_model("hDCAzAHe", ";#it{p}_{T}^{rec} (GeV/#it{c});DCA_{z} (cm);Counts", 100, -0.2, 0.2), pt_name, "fDCAz")]
    h_dca_z_m = [df_primary.Filter(f"matter && {nsigma} > -0.5 && {nsigma} < 3 && {mass_cut}").Histo2D(h2_model("hDCAzMHe", ";#it{p}_{T}^{rec} (GeV/#it{c});DCA_{z} (cm);Counts", 100, -0.2, 0.2), pt_name, "fDCAz")]
    h_dca_xy_sec_m = [df_secondary.Filter(f"matter && {nsigma} > -0.5 && {nsigma} < 3 && {mass_cut}").Histo2D(h2_model("hDCAxySecondaryMHe", ";#it{p}_{T}^{rec} (GeV/#it{c});DCA_{xy} (cm);Counts", 100, -0.2, 0.2), pt_name, "fDCAxy")]
    h_dca_xy_sec_a = [df_secondary.Filter(f"!matter && {nsigma} > -0.5 && {nsigma} < 3 && {mass_cut}").Histo2D(h2_model("hDCAxySecondaryAHe", ";#it{p}_{T}^{rec} (GeV/#it{c});DCA_{xy} (cm);Counts", 100, -0.2, 0.2), pt_name, "fDCAxy")]

    h_tpc_a = [df_primary.Filter(f"!matter && {mass_cut}" if particle == "he3" else "!matter").Histo2D(tpc_anti_model, "ptUncorr" if particle == "he4" else pt_name, nsigma)]
    h_tpc_m = [df_primary.Filter(f"matter && {mass_cut}" if particle == "he3" else "matter").Histo2D(tpc_mat_model, pt_name, nsigma)]

    nsigma_tof_cut = HE4_NSIGMA_TOF_CUT if particle == "he4" else HE3_NSIGMA_TOF_CUT
    h_tof_a = [df_primary.Filter(f"!matter && std::abs({nsigma}) < {nsigma_tof_cut}").Histo2D(h2_model("fATOFsignal", f";#it{{p}}_{{T}}^{{rec}} (GeV/#it{{c}});m_{{TOF}}-m_{{^{{{'4' if particle == 'he4' else '3'}}}#bar{{He}}}};Counts", 100, -0.9, 1.1), pt_name, dmass)]
    h_tof_m = [df_primary.Filter(f"matter && std::abs({nsigma}) < {nsigma_tof_cut}").Histo2D(h2_model("fMTOFsignal", f";#it{{p}}_{{T}}^{{rec}} (GeV/#it{{c}});m_{{TOF}}-m_{{^{{{'4' if particle == 'he4' else '3'}}}He}};Counts", 100, -0.9, 1.1), pt_name, dmass)]

    i_trial = 0
    if particle == "he3":
        for cut_dca_z in CUT_NAMES["nsigmaDCAz"]:
            df_dca_z = df_base.Filter(f"std::abs(nsigmaDCAz) < {cut_dca_z}")
            for cut_tpc in CUT_NAMES["fTPCnCls"]:
                df_tpc = df_dca_z.Filter(f"fTPCnCls > {cut_tpc}")
                for cut_its in CUT_NAMES["nITScls"]:
                    df_its = df_tpc.Filter(f"nITScls >= {cut_its}")
                    h_dca_xy_a.append(df_its.Filter(f"!matter && {nsigma} > -0.5 && {nsigma} < 3 && {mass_cut}").Histo2D(h2_model(f"hDCAxyAHe3{i_trial}", ";#it{p}_{T}^{rec} (GeV/#it{c});DCA_{xy} (cm);Counts", 560, -0.7, 0.7), pt_name, "fDCAxy"))
                    h_dca_xy_m.append(df_its.Filter(f"matter && {nsigma} > -0.5 && {nsigma} < 3 && {mass_cut}").Histo2D(h2_model(f"hDCAxyMHe3{i_trial}", ";#it{p}_{T}^{rec} (GeV/#it{c});DCA_{xy} (cm);Counts", 560, -0.7, 0.7), pt_name, "fDCAxy"))
                    h_dca_z_a.append(df_its.Filter(f"!matter && {nsigma} > -0.5 && {nsigma} < 3 && {mass_cut}").Histo2D(h2_model(f"hDCAzAHe3{i_trial}", ";#it{p}_{T}^{rec} (GeV/#it{c});DCA_{z} (cm);Counts", 560, -0.7, 0.7), pt_name, "fDCAz"))
                    h_dca_z_m.append(df_its.Filter(f"matter && {nsigma} > -0.5 && {nsigma} < 3 && {mass_cut}").Histo2D(h2_model(f"hDCAzMHe3{i_trial}", ";#it{p}_{T}^{rec} (GeV/#it{c});DCA_{z} (cm);Counts", 560, -0.7, 0.7), pt_name, "fDCAz"))

                    h_tpc_a.append(df_its.Filter(f"!matter && {HE3_TRIAL_DCA_SELECTION} && {mass_cut}").Histo2D(h2_model(f"fATPCcounts{i_trial}", ";#it{p}_{T}^{rec} (GeV/#it{c});^{3}#bar{He} n#sigma_{TPC};Counts", 100, -5, 5), pt_name, nsigma))
                    h_tpc_m.append(df_its.Filter(f"matter && {HE3_TRIAL_DCA_SELECTION} && {mass_cut}").Histo2D(h2_model(f"fMTPCcounts{i_trial}", ";#it{p}_{T}^{rec} (GeV/#it{c});^{3}He n#sigma_{TPC};Counts", 100, -5, 5), pt_name, nsigma))
                    h_tof_a.append(df_its.Filter(f"!matter && {HE3_TRIAL_DCA_SELECTION} && std::abs({nsigma}) < {HE3_NSIGMA_TOF_CUT}").Histo2D(h2_model(f"fATOFsignal{i_trial}", ";#it{p}_{T}^{rec} (GeV/#it{c});m_{TOF}-m_{^{3}#bar{He}};Counts", 100, -0.9, 1.1), pt_name, dmass))
                    h_tof_m.append(df_its.Filter(f"matter && {HE3_TRIAL_DCA_SELECTION} && std::abs({nsigma}) < {HE3_NSIGMA_TOF_CUT}").Histo2D(h2_model(f"fMTOFsignal{i_trial}", ";#it{p}_{T}^{rec} (GeV/#it{c});m_{TOF}-m_{^{3}He};Counts", 100, -0.9, 1.1), pt_name, dmass))
                    i_trial += 1

    if draw:
        for hist in [h_tpc_a[0], h_tpc_m[0], h_tof_a[0], h_tof_m[0], h_dca_xy_a[0], h_dca_xy_m[0], h_dca_z_a[0], h_dca_z_m[0], h_dca_xy_sec_m[0], h_dca_xy_sec_a[0]]:
            c = ROOT.TCanvas()
            c.SetRightMargin(0.15)
            hist.DrawClone("colz")

    output_file = expand(output_file)
    ensure_parent(output_file)
    out = ROOT.TFile(output_file, "recreate")
    d0 = out.mkdir("nuclei")
    d0.cd()
    write_hist(h_tpc_a[0], "fATPCcounts")
    write_hist(h_tpc_m[0], "fMTPCcounts")
    write_hist(h_tof_a[0], "fATOFsignal")
    write_hist(h_tof_m[0], "fMTOFsignal")
    write_hist(h_dca_xy_a[0], "hDCAxyAHe3")
    write_hist(h_dca_xy_m[0], "hDCAxyMHe3")
    write_hist(h_dca_z_a[0], "hDCAzAHe3")
    write_hist(h_dca_z_m[0], "hDCAzMHe3")
    write_hist(h_dca_xy_sec_m[0], "hDCAxySecondaryMHe3")
    write_hist(h_dca_xy_sec_a[0], "hDCAxySecondaryAHe3")

    for idx in range(i_trial):
        dtrial = out.mkdir(f"nuclei{idx}")
        dtrial.cd()
        write_hist(h_tpc_a[idx + 1], "fATPCcounts")
        write_hist(h_tpc_m[idx + 1], "fMTPCcounts")
        write_hist(h_tof_a[idx + 1], "fATOFsignal")
        write_hist(h_tof_m[idx + 1], "fMTOFsignal")
        write_hist(h_dca_xy_a[idx + 1], "hDCAxyAHe3")
        write_hist(h_dca_xy_m[idx + 1], "hDCAxyMHe3")
        write_hist(h_dca_z_a[idx + 1], "hDCAzAHe3")
        write_hist(h_dca_z_m[idx + 1], "hDCAzMHe3")
    out.Close()
    LOGGER.info("analyse_data done output=%s", output_file)


def analyse_mc(input_file: str, output_file: str, particle: str, enable_trials: bool, draw: bool = False) -> None:
    LOGGER.info("analyse_mc start particle=%s input=%s output=%s", particle, input_file, output_file)
    ROOT.gStyle.SetOptStat(0)
    rdf = ROOT.RDataFrame("O2nucleitablemc", expand(input_file))
    df_data = define_columns_for_data(rdf)

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

        h_reco_tpc_a = [df_cut_reco.Filter(f"!matter && {HE4_MC_SIGNAL_TRACKING}").Histo1D(h1_model("TPCAHe4", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts"), "pt")]
        h_reco_tpc_m = [df_cut_reco.Filter(f"matter && {HE4_MC_SIGNAL_TRACKING}").Histo1D(h1_model("TPCMHe4", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts"), "pt")]
        h_reco_tof_a = [df_cut_reco.Filter(f"!matter && {HE4_MC_SIGNAL_TRACKING} && hasTOF").Histo1D(h1_model("TOFAHe4", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts"), "pt")]
        h_reco_tof_m = [df_cut_reco.Filter(f"matter && {HE4_MC_SIGNAL_TRACKING} && hasTOF").Histo1D(h1_model("TOFMHe4", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts"), "pt")]
        h_gen_a = [df_cut_gen.Filter("fPDGcode < 0").Histo1D(h1_model("genAHe4", ";#it{p}_{T}^{gen} (GeV/#it{c});Counts"), "fgPt")]
        h_gen_m = [df_cut_gen.Filter("fPDGcode > 0").Histo1D(h1_model("genMHe4", ";#it{p}_{T}^{gen} (GeV/#it{c});Counts"), "fgPt")]

        h_reco_tpc_a_w = [df_cut_reco.Filter(f"!matter && {HE4_MC_SIGNAL_TRACKING}").Histo1D(h1_model("TPCAHe4W", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts"), "pt", "ptWeight")]
        h_reco_tpc_m_w = [df_cut_reco.Filter(f"matter && {HE4_MC_SIGNAL_TRACKING}").Histo1D(h1_model("TPCMHe4W", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts"), "pt", "ptWeight")]
        h_reco_tof_a_w = [df_cut_reco.Filter(f"!matter && {HE4_MC_SIGNAL_TRACKING} && hasTOF").Histo1D(h1_model("TOFAHe4W", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts"), "pt", "ptWeight")]
        h_reco_tof_m_w = [df_cut_reco.Filter(f"matter && {HE4_MC_SIGNAL_TRACKING} && hasTOF").Histo1D(h1_model("TOFMHe4W", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts"), "pt", "ptWeight")]
        h_gen_a_w = [df_cut_gen.Filter("fPDGcode < 0").Histo1D(h1_model("genAHe4W", ";#it{p}_{T}^{gen} (GeV/#it{c});Counts"), "fgPt", "ptWeight")]
        h_gen_m_w = [df_cut_gen.Filter("fPDGcode > 0").Histo1D(h1_model("genMHe4W", ";#it{p}_{T}^{gen} (GeV/#it{c});Counts"), "fgPt", "ptWeight")]

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

        h_reco_tpc_a = [df_cut_reco.Filter("!matter").Histo1D(h1_model("TPCAHe3", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts"), "pt")]
        h_reco_tpc_m = [df_cut_reco.Filter("matter").Histo1D(h1_model("TPCMHe3", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts"), "pt")]
        h_reco_tof_a = [df_cut_reco.Filter("!matter && hasTOF").Histo1D(h1_model("TOFAHe3", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts"), "pt")]
        h_reco_tof_m = [df_cut_reco.Filter("matter && hasTOF").Histo1D(h1_model("TOFMHe3", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts"), "pt")]
        h_gen_a = [df_cut_gen.Filter("fPDGcode < 0").Histo1D(h1_model("genAHe3", ";#it{p}_{T}^{gen} (GeV/#it{c});Counts"), "fgPt")]
        h_gen_m = [df_cut_gen.Filter("fPDGcode > 0").Histo1D(h1_model("genMHe3", ";#it{p}_{T}^{gen} (GeV/#it{c});Counts"), "fgPt")]

        h_reco_tpc_a_w = [df_cut_reco.Filter("!matter").Histo1D(h1_model("TPCAHe3W", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts"), "pt", "ptWeight")]
        h_reco_tpc_m_w = [df_cut_reco.Filter("matter").Histo1D(h1_model("TPCMHe3W", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts"), "pt", "ptWeight")]
        h_reco_tof_a_w = [df_cut_reco.Filter("!matter && hasTOF").Histo1D(h1_model("TOFAHe3W", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts"), "pt", "ptWeight")]
        h_reco_tof_m_w = [df_cut_reco.Filter("matter && hasTOF").Histo1D(h1_model("TOFMHe3W", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts"), "pt", "ptWeight")]
        h_gen_a_w = [df_cut_gen.Filter("fPDGcode < 0").Histo1D(h1_model("genAHe3W", ";#it{p}_{T}^{gen} (GeV/#it{c});Counts"), "fgPt", "ptWeight")]
        h_gen_m_w = [df_cut_gen.Filter("fPDGcode > 0").Histo1D(h1_model("genMHe3W", ";#it{p}_{T}^{gen} (GeV/#it{c});Counts"), "fgPt", "ptWeight")]

        reco_df_for_trials = df_cut_reco_base
        tpc_name = "He3"
        gen_name = "He3"

    n_trials = 0
    if enable_trials:
        for cut_dca_z in CUT_NAMES["nsigmaDCAz"]:
            d_dca = reco_df_for_trials.Filter(f"std::abs(nsigmaDCAz) < {cut_dca_z}")
            for cut_tpc in CUT_NAMES["fTPCnCls"]:
                d_tpc = d_dca.Filter(f"fTPCnCls > {cut_tpc}")
                for cut_its in CUT_NAMES["nITScls"]:
                    _ = d_tpc.Filter(f"nITScls >= {cut_its}")
                    h_reco_tpc_a.append(reco_df_for_trials.Filter("!matter").Histo1D(h1_model(f"TPCA{tpc_name}{n_trials}", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts"), "pt"))
                    h_reco_tpc_m.append(reco_df_for_trials.Filter("matter").Histo1D(h1_model(f"TPCM{tpc_name}{n_trials}", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts"), "pt"))
                    h_reco_tof_a.append(reco_df_for_trials.Filter("!matter && hasTOF").Histo1D(h1_model(f"TOFA{tpc_name}{n_trials}", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts"), "pt"))
                    h_reco_tof_m.append(reco_df_for_trials.Filter("matter && hasTOF").Histo1D(h1_model(f"TOFM{tpc_name}{n_trials}", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts"), "pt"))

                    h_reco_tpc_a_w.append(reco_df_for_trials.Filter("!matter").Histo1D(h1_model(f"TPCA{tpc_name}W{n_trials}", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts"), "pt", "ptWeight"))
                    h_reco_tpc_m_w.append(reco_df_for_trials.Filter("matter").Histo1D(h1_model(f"TPCM{tpc_name}W{n_trials}", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts"), "pt", "ptWeight"))
                    h_reco_tof_a_w.append(reco_df_for_trials.Filter("!matter && hasTOF").Histo1D(h1_model(f"TOFA{tpc_name}W{n_trials}", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts"), "pt", "ptWeight"))
                    h_reco_tof_m_w.append(reco_df_for_trials.Filter("matter && hasTOF").Histo1D(h1_model(f"TOFM{tpc_name}W{n_trials}", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts"), "pt", "ptWeight"))
                    n_trials += 1

    if draw:
        c = ROOT.TCanvas("effMatterAntiMatter", "effMatterAntiMatter")
        c.cd()
        h_eff_a = h_reco_tpc_a[0].GetValue().Clone("hEffA")
        h_eff_a.Divide(h_reco_tpc_a[0].GetPtr(), h_gen_a[0].GetPtr(), 1.0, 1.0, "B")
        h_eff_a.SetLineColor(ROOT.kRed)
        h_eff_a.Draw()
        h_eff_m = h_reco_tpc_m[0].GetValue().Clone("hEffM")
        h_eff_m.Divide(h_reco_tpc_m[0].GetPtr(), h_gen_m[0].GetPtr(), 1.0, 1.0, "B")
        h_eff_m.Draw("same")

    out_name = expand(output_file)
    ensure_parent(out_name)
    out = ROOT.TFile(out_name, "recreate")
    d0 = out.mkdir("nuclei")
    d0.cd()

    write_hist(h_gen_a[0], f"genA{gen_name}")
    write_hist(h_gen_m[0], f"genM{gen_name}")
    write_hist(h_reco_tpc_a[0], f"TPCA{tpc_name}")
    write_hist(h_reco_tpc_m[0], f"TPCM{tpc_name}")
    write_hist(h_reco_tof_a[0], f"TOFA{tpc_name}")
    write_hist(h_reco_tof_m[0], f"TOFM{tpc_name}")

    write_hist(h_gen_a_w[0], f"genA{gen_name}W")
    write_hist(h_gen_m_w[0], f"genM{gen_name}W")
    write_hist(h_reco_tpc_a_w[0], f"TPCA{tpc_name}W")
    write_hist(h_reco_tpc_m_w[0], f"TPCM{tpc_name}W")
    write_hist(h_reco_tof_a_w[0], f"TOFA{tpc_name}W")
    write_hist(h_reco_tof_m_w[0], f"TOFM{tpc_name}W")

    def write_eff(prefix: str, reco_a: Any, reco_m: Any, reco_tof_a: Any, reco_tof_m: Any, gen_a: Any, gen_m: Any) -> tuple[Any, Any]:
        eff_tpc_a = reco_a.GetValue().Clone(f"{prefix}TPCA")
        eff_tpc_m = reco_m.GetValue().Clone(f"{prefix}TPCM")
        eff_tof_a = reco_tof_a.GetValue().Clone(f"{prefix}TOFA")
        eff_tof_m = reco_tof_m.GetValue().Clone(f"{prefix}TOFM")
        eff_tpc_a.Divide(reco_a.GetPtr(), gen_a.GetPtr(), 1.0, 1.0, "B")
        eff_tpc_m.Divide(reco_m.GetPtr(), gen_m.GetPtr(), 1.0, 1.0, "B")
        eff_tof_a.Divide(reco_tof_a.GetPtr(), gen_a.GetPtr(), 1.0, 1.0, "B")
        eff_tof_m.Divide(reco_tof_m.GetPtr(), gen_m.GetPtr(), 1.0, 1.0, "B")
        for h in (eff_tpc_a, eff_tpc_m, eff_tof_a, eff_tof_m):
            h.GetYaxis().SetTitle("Efficiency #times Acceptance")
        eff_tpc_a.Write(f"{prefix}effTPCA")
        eff_tpc_m.Write(f"{prefix}effTPCM")
        eff_tof_a.Write(f"{prefix}effTOFA")
        eff_tof_m.Write(f"{prefix}effTOFM")
        return eff_tpc_a, eff_tpc_m

    write_eff("", h_reco_tpc_a[0], h_reco_tpc_m[0], h_reco_tof_a[0], h_reco_tof_m[0], h_gen_a[0], h_gen_m[0])
    write_eff(_weighted_eff_name(""), h_reco_tpc_a_w[0], h_reco_tpc_m_w[0], h_reco_tof_a_w[0], h_reco_tof_m_w[0], h_gen_a_w[0], h_gen_m_w[0])

    for i in range(n_trials):
        dtrial = out.mkdir(f"nuclei{i}")
        dtrial.cd()
        write_hist(h_gen_a[0], f"genA{gen_name}")
        write_hist(h_gen_m[0], f"genM{gen_name}")
        write_hist(h_reco_tpc_a[i + 1], f"TPCA{tpc_name}")
        write_hist(h_reco_tpc_m[i + 1], f"TPCM{tpc_name}")
        write_hist(h_reco_tof_a[i + 1], f"TOFA{tpc_name}")
        write_hist(h_reco_tof_m[i + 1], f"TOFM{tpc_name}")

        write_hist(h_gen_a_w[0], f"genA{gen_name}W")
        write_hist(h_gen_m_w[0], f"genM{gen_name}W")
        write_hist(h_reco_tpc_a_w[i + 1], f"TPCA{tpc_name}W")
        write_hist(h_reco_tpc_m_w[i + 1], f"TPCM{tpc_name}W")
        write_hist(h_reco_tof_a_w[i + 1], f"TOFA{tpc_name}W")
        write_hist(h_reco_tof_m_w[i + 1], f"TOFM{tpc_name}W")

        eff_tpc_a, eff_tpc_m = write_eff("", h_reco_tpc_a[i + 1], h_reco_tpc_m[i + 1], h_reco_tof_a[i + 1], h_reco_tof_m[i + 1], h_gen_a[0], h_gen_m[0])
        matching_tof_a = h_reco_tof_a[i + 1].GetValue().Clone(f"matchingTOFA{i}")
        matching_tof_m = h_reco_tof_m[i + 1].GetValue().Clone(f"matchingTOFM{i}")
        matching_tof_a.Divide(eff_tpc_a)
        matching_tof_m.Divide(eff_tpc_m)
        matching_tof_a.Write()
        matching_tof_m.Write()

        eff_w_tpc_a, eff_w_tpc_m = write_eff(_weighted_eff_name(""), h_reco_tpc_a_w[i + 1], h_reco_tpc_m_w[i + 1], h_reco_tof_a_w[i + 1], h_reco_tof_m_w[i + 1], h_gen_a_w[0], h_gen_m_w[0])
        matching_w_tof_a = h_reco_tof_a_w[i + 1].GetValue().Clone(f"matchingWTOFA{i}")
        matching_w_tof_m = h_reco_tof_m_w[i + 1].GetValue().Clone(f"matchingWTOFM{i}")
        matching_w_tof_a.Divide(eff_w_tpc_a)
        matching_w_tof_m.Divide(eff_w_tpc_m)
        matching_w_tof_a.Write()
        matching_w_tof_m.Write()

    out.Close()
    LOGGER.info("analyse_mc done output=%s", output_file)


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

    h_ntvx = an_results.Get("bc-selection-task/hCounterTVX")
    h_nvtx = an_results.Get("nuclei-spectra/spectra/hRecVtxZData")
    _ = h_ntvx.GetEntries()
    norm = h_nvtx.GetEntries()

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


def merge_trees(input_file: str, output_file: str, is_mc: bool = True) -> None:
    LOGGER.info("merge_trees start input=%s output=%s is_mc=%s", input_file, output_file, is_mc)
    in_f = ROOT.TFile(expand(input_file))
    keys = ROOT.gDirectory.GetListOfKeys()
    out_list = ROOT.TList()
    tree_name = "O2nucleitablemc" if is_mc else "O2nucleitable"

    for i in range(keys.GetEntries()):
        key = keys.At(i)
        key_name = key.GetName()
        if key.IsFolder() and key_name.startswith("DF"):
            d = in_f.Get(key_name)
            tree = d.Get(tree_name)
            if tree:
                out_list.Add(tree)

    out_name = expand(output_file)
    ensure_parent(out_name)
    out_f = ROOT.TFile(out_name, "RECREATE")
    if out_list.GetSize() == 0:
        raise RuntimeError("No trees found to merge.")
    new_tree = ROOT.TTree.MergeTrees(out_list)
    new_tree.SetName(tree_name)
    new_tree.Write()
    out_f.Close()
    LOGGER.info("merge_trees done output=%s", output_file)

