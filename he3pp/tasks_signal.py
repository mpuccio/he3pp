from array import array
import logging

import ROOT

from .root_io import ensure_parent, expand, load_fit_modules, ptr
from .settings import *  # Loaded after runtime overrides in cli.run


LOGGER = logging.getLogger("he3pp.tasks")


def signal(input_file: str, output_file: str, particle: str = "he3") -> None:
    LOGGER.info("signal start particle=%s input=%s output=%s", particle, input_file, output_file)
    p = get_particle_config(particle)
    species_names = [str(p["name"]), str(p["anti_name"])]
    matter_label = str(p["label_matter"])
    load_fit_modules()
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)
    ROOT.RooMsgService.instance().setSilentMode(True)
    ROOT.gErrorIgnoreLevel = ROOT.kError

    in_file = ROOT.TFile(expand(input_file))
    out_name = expand(output_file)
    ensure_parent(out_name)
    out_file = ROOT.TFile(out_name, "recreate")

    m = ROOT.RooRealVar("dm2", f"m - m_{{{matter_label}}}", -1.2, 1.5, "GeV/#it{c}^{2}")
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

    m_bis = ROOT.RooRealVar("dm2_bis", f"m - m_{{{matter_label}}}", -1.2, 1.5, "GeV/#it{c}^{2}")
    m_bis.setBins(1000, "cache")
    m_bis.setRange("Full", -1.2, 1.5)
    f_bkg = ROOT.FitExpExpTailGaus(m_bis)
    f_bkg.UseSignal(False)

    ns = ROOT.RooRealVar("ns", f"n#sigma_{{{matter_label}}}", -5.0, 5.0, "a. u.")
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
    tpc_model_names = [str(name) for name in TPC_FUNCTION_NAMES]
    if len(tpc_model_names) > len(tpc_functions):
        raise RuntimeError(
            f"Configured {len(tpc_model_names)} TPC function names but only {len(tpc_functions)} fit functions are available."
        )
    n_tpc_functions = len(tpc_model_names)

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
        h_npar = [[None] * CENT_LENGTH for _ in range(2)]
        h_npar_tpc = [[[None] * n_tpc_functions for _ in range(CENT_LENGTH)] for _ in range(2)]
        h_tpconly = [[[None] * n_tpc_functions for _ in range(CENT_LENGTH)] for _ in range(2)]
        h_widen = [[None] * CENT_LENGTH for _ in range(2)]
        h_shift = [[None] * CENT_LENGTH for _ in range(2)]
        h_widen_tpc = [[None] * CENT_LENGTH for _ in range(2)]
        h_shift_tpc = [[None] * CENT_LENGTH for _ in range(2)]

        n_sigma_vec = [3.0]
        v_shift = []

        for i_s in range(2):
            d_species = d_out.mkdir(species_names[i_s])
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
            for i_t, model_name in enumerate(tpc_model_names):
                h_tpconly[i_s][i_c][i_t] = ROOT.TH1D(f"hTPConly{LETTER[i_s]}{i_c}_{model_name}", ";p_{T} GeV/c; TPC raw counts", n_pt_bins, pt_labels.GetArray())
                h_npar_tpc[i_s][i_c][i_t] = ROOT.TH1D(f"hNFloatParsTPC{LETTER[i_s]}{i_c}_{model_name}", "; p_{T}(GeV/c); N_{float} (TPC fit)", n_pt_bins, pt_labels.GetArray())
            h_signif[i_s][i_c] = ROOT.TH1D(f"hSignificance{LETTER[i_s]}{i_c}", "; p_{T}(GeV/c); #frac{S}{#sqrt{S+B}}", n_pt_bins, pt_labels.GetArray())
            h_chi[i_s][i_c] = ROOT.TH1D(f"hChiSquare{LETTER[i_s]}{i_c}", "; p_{T}(GeV/c); #chi^{2}/NDF", n_pt_bins, pt_labels.GetArray())
            h_chi_tpc[i_s][i_c] = ROOT.TH1D(f"hChiSquareTPC{LETTER[i_s]}{i_c}", "; p_{T}(GeV/c); #chi^{2}/NDF", n_pt_bins, pt_labels.GetArray())
            h_npar[i_s][i_c] = ROOT.TH1D(f"hNFloatPars{LETTER[i_s]}{i_c}", "; p_{T}(GeV/c); N_{float} (TOF fit)", n_pt_bins, pt_labels.GetArray())
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

                d_out.cd(f"{species_names[i_s]}/GausExp/C_{i_c}")
                fit_plot = f_sig.FitData(dat, i_name, i_title, "Full", "Full", False, -1.2, 1.5)
                ptr(f_sig.mSigma).setConstant(False)
                if center_pt > TOF_MIN_PT:
                    fit_plot.Write()
                h_npar[i_s][i_c].SetBinContent(i_b + 1, float(getattr(f_sig, "mNFloatPars", 0)))
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
                        d_out.cd(f"{species_names[i_s]}/Sidebands/C_{i_c}")
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
                    d_out.cd(f"{species_names[i_s]}/TPConly")
                    tpc_dat = tpc_h[i_s].ProjectionY(f"tpc_data{i_c}_{i_b}", i_b + 1, i_b + 1)
                    for i_t, model_name in enumerate(tpc_model_names):
                        fit_range = "Special" if (i_t and center_pt < 1.8) else "Full"
                        tpc_plot = tpc_functions[i_t].FitData(tpc_dat, f"TPC_d_{i_c}_{i_b}_{model_name}", i_title, fit_range, fit_range)
                        ptr(tpc_functions[i_t].mSigma).setConstant(False)
                        tpc_plot.Write()
                        h_npar_tpc[i_s][i_c][i_t].SetBinContent(i_b + 1, float(getattr(tpc_functions[i_t], "mNFloatPars", 0)))
                        h_tpconly[i_s][i_c][i_t].SetBinContent(i_b + 1, ptr(tpc_functions[i_t].mSigCounts).getVal())
                        h_tpconly[i_s][i_c][i_t].SetBinError(i_b + 1, ptr(tpc_functions[i_t].mSigCounts).getError())

        for i_s in range(2):
            i_c = 0
            d_out.cd(f"{species_names[i_s]}/GausExp")
            h_raw[i_s][i_c].Write()
            h_raw_bc[i_s][i_c].Write()
            h_sig[i_s][i_c].Write()

            d_out.cd(f"{species_names[i_s]}/Systematic")
            h_shift[i_s][i_c].Write()
            h_widen[i_s][i_c].Write()
            h_widen_tpc[i_s][i_c].Write()
            h_shift_tpc[i_s][i_c].Write()

            d_out.cd(f"{species_names[i_s]}/Significance")
            h_signif[i_s][i_c].Write()

            d_out.cd(f"{species_names[i_s]}/TPConly")
            for i_t in range(n_tpc_functions):
                h_tpconly[i_s][i_c][i_t].Write()
                h_npar_tpc[i_s][i_c][i_t].Write()

            d_out.cd(f"{species_names[i_s]}/ChiSquare")
            h_chi[i_s][i_c].Write()
            h_chi_tpc[i_s][i_c].Write()
            h_npar[i_s][i_c].Write()

    out_file.Close()
    LOGGER.info("signal done particle=%s output=%s", particle, output_file)
