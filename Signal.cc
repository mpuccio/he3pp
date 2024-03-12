#include "src/Common.h"
#include "src/FitModules.h"
#include "src/Utils.h"
using namespace utils;

#include <memory>
#include <functional>
using std::function;
#include <utility>
using std::pair;
#include <vector>
using std::vector;

#include <TAxis.h>
#include <TDirectory.h>
#include <TEfficiency.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TList.h>
#include <TMath.h>
#include <TColor.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TError.h>

#include <RooArgList.h>
#include <RooFitResult.h>
#include <RooDataHist.h>
#include <RooMsgService.h>
#include <RooExponential.h>
#include <RooGaussian.h>
#include <RooAddPdf.h>
#include <RooPlot.h>
#include <RooRealVar.h>

using namespace RooFit;

void Signal()
{

  /// Suppressing the output
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  RooMsgService::instance().setSilentMode(true);
  gErrorIgnoreLevel = kError; // Suppressing warning outputs

  /// Taking all the histograms from the MC file
  TFile input_file(gSystem->ExpandPathName(kDataFilename.data()));
  TFile output_file(gSystem->ExpandPathName(kSignalOutput.data()), "recreate");

  /// Setting up the fitting environment for TOF analysis
  RooRealVar m("dm2", "m - m_{^{3}He}", -1.2, 1.5, "GeV/#it{c}^{2}");
  m.setBins(1000, "cache");
  m.setRange("Full", -1.2, 1.5);

  FitExpTailGaus fExpExpTailGaus(&m);
  fExpExpTailGaus.mMu->setRange(-1, 1.);
  fExpExpTailGaus.mMu->setVal(0.1);
  fExpExpTailGaus.mMu->setUnit("GeV/#it{c}^{2}");
  fExpExpTailGaus.mSigma->setRange(0.05, 0.40);
  fExpExpTailGaus.mSigma->setVal(0.1);
  fExpExpTailGaus.mSigma->setUnit("GeV/#it{c}^{2}");
  fExpExpTailGaus.mAlpha0->setRange(0.8, 3.);
  fExpExpTailGaus.mAlpha0->setVal(1.2);
  fExpExpTailGaus.mAlpha0->setUnit("GeV/#it{c}^{2}");
  fExpExpTailGaus.mSigCounts->setRange(0., 5000.);
  fExpExpTailGaus.mTau0->setUnit("GeV#it{c}^{2}");

  // Background
  RooRealVar m_bis("dm2_bis", "m - m_{^{3}He}", -1.2, 1.5, "GeV/#it{c}^{2}");
  m_bis.setBins(1000, "cache");
  m_bis.setRange("Full", -1.2, 1.5);
  FitExpExpTailGaus fBkg(&m_bis);
  fBkg.UseSignal(false);
  fBkg.mTau0->setUnit("GeV#it{c}^{2}");
  fBkg.mTau1->setUnit("GeV#it{c}^{2}");

  // Setting up the fitting environment for the TPC analysis
  RooRealVar ns("ns", "n#sigma_{^{3}He}", -5., 5, "a. u.");
  ns.setBins(1000, "cache");
  ns.setRange("Full", -5., 5.);
  ns.setRange("Special", -4, 5.);

  // TPC analysis
  FitGausGaus fGausGaus(&ns);
  fGausGaus.mSigma->setRange(0.2, 1.2);
  fGausGaus.mSigma->setVal(1.);
  fGausGaus.mSigma->setUnit("a. u.");
  fGausGaus.mMu->setRange(-0.5, 0.5);
  fGausGaus.mMu->setUnit("a. u.");
  fGausGaus.mMuBkg->setRange(-10., -4.);
  fGausGaus.mMuBkg->setVal(-7);
  fGausGaus.mMuBkg->setUnit("a. u.");
  fGausGaus.mSigmaBkg->setRange(0.2, 6.);
  fGausGaus.mSigmaBkg->setUnit("a. u.");

  FitExpGaus fExpGausTPC(&ns);
  fExpGausTPC.mSigma->setRange(0.2, 1.2);
  fExpGausTPC.mSigma->setVal(1.);
  fExpGausTPC.mSigma->setUnit("a. u.");
  fExpGausTPC.mMu->setRange(-0.5, 0.5);
  fExpGausTPC.mMu->setUnit("a. u.");


  FitExpTailGaus fExpTailGausTPC(&ns);
  fExpTailGausTPC.mSigma->setRange(0.2, 1.2);
  fExpTailGausTPC.mSigma->setVal(1.);
  fExpTailGausTPC.mSigma->setUnit("a. u.");
  fExpTailGausTPC.mMu->setRange(-0.5, 0.5);
  fExpTailGausTPC.mMu->setUnit("a. u.");

  FitLogNormalLogNormal fLogNormalLogNormalTPC(&ns);
  fLogNormalLogNormalTPC.mSigma->setRange(1.01, 20.);
  fLogNormalLogNormalTPC.mSigma->setVal(TMath::Exp(1.));
  fLogNormalLogNormalTPC.mSigma->setUnit("a. u.");
  fLogNormalLogNormalTPC.mMu->setRange(-0.5, 0.5);
  fLogNormalLogNormalTPC.mMu->setUnit("a. u.");

  FitModule *tpcFunctions[] = {&fGausGaus, &fExpGausTPC, &fExpTailGausTPC, &fLogNormalLogNormalTPC};

  for (auto list_key : *input_file.GetListOfKeys())
  {
    /// Preliminary operation to read the list and create an output dir
    if (string(list_key->GetName()).find(kFilterListNames.data()) == string::npos)
      continue;

    TDirectoryFile *list = (TDirectoryFile *)input_file.Get(list_key->GetName());
    TDirectory *base_dir = output_file.mkdir(list_key->GetName());
    output_file.cd(list_key->GetName());
    std::cout << "Analysing directory " << list_key->GetName() << std::endl;

    /// Taking all the necessary histogram to perform the analysis
    TH2 *fATOFsignal = (TH2 *)list->Get("fATOFsignal");
    TH2 *fMTOFsignal = (TH2 *)list->Get("fMTOFsignal");
    TH2 *fATPCcounts = (TH2 *)list->Get("fATPCcounts");
    TH2 *fMTPCcounts = (TH2 *)list->Get("fMTPCcounts");

    /// Taking information about centrality bins
    /// Taking information about \f$p_{\mathrm{T}}\f$ bins
    const int n_pt_bins = fATOFsignal->GetNbinsX();
    auto pt_axis = fATOFsignal->GetXaxis();
    auto pt_labels = *(pt_axis->GetXbins());

    TH2 *tof_histo[2] = {fMTOFsignal, fATOFsignal};
    TH2 *tpc_histo[2] = {fMTPCcounts, fATPCcounts};

    TH1D *hRawCounts[2][kCentLength];
    TH1D *hRawCountsBinCounting[2][kCentLength];
    TH1D *hSignalGausExpGaus[2][kCentLength];
    TH1D *hSystFit[2][kCentLength];
    TH1D *hSignificance[2][kCentLength];
    TH1D *hChiSquare[2][kCentLength];
    TH1D *hChiSquareTPC[2][kCentLength];
    TH1D *hTPConly[2][kCentLength][kNTPCfunctions];

    vector<float> n_sigma_vec{3.}; // = {3.0, 3.1, 3.2, 3.3, 3.4, 3.5};
    vector<float> v_shift;         // = {-0.1, 0.05, 0., 0.05, 0.1};
    int n_shifts = v_shift.size();
    int kNewGreen = kGreen + 3;
    int color_vector[] = {kBlack, kBlue, kNewGreen, kOrange, kRed};
    TH1D *hWidenRangeSyst[2][kCentLength];
    TH1D *hShiftRangeSyst[2][kCentLength];
    TH1D *hWidenRangeSystTPC[2][kCentLength];
    TH1D *hShiftRangeSystTPC[2][kCentLength];

    /// Creating the directories to be used to store the results
    for (int iS = 0; iS < 2; ++iS)
    {
      TDirectory *dir = base_dir->mkdir(kNames[iS].data());
      dir->cd();
      // dir->mkdir("Fits");
      TDirectory *sig_dir = dir->mkdir("GausExp");
      sig_dir->cd();
      for (int iC = 0; iC < kCentLength; ++iC)
        sig_dir->mkdir(Form("C_%d", iC));
      dir->cd();
      TDirectory *side_dir = dir->mkdir("Sidebands");
      side_dir->cd();
      for (int iC = 0; iC < kCentLength; ++iC)
        side_dir->mkdir(Form("C_%d", iC));
      dir->cd();
      dir->mkdir("Significance");
      dir->mkdir("Systematic");
      dir->mkdir("TPConly");
      dir->mkdir("ChiSquare");
    }

    for (int iS = 0; iS < 2; ++iS)
    {
      int iC{0};
      for (int iT{0}; iT < kNTPCfunctions; ++iT)
      {
        hTPConly[iS][iC][iT] = new TH1D(Form("hTPConly%c%i_%s", kLetter[iS], iC, kTPCfunctName[iT].data()), ";p_{T} GeV/c; TPC raw counts", n_pt_bins, pt_labels.GetArray());
      }
      hSignificance[iS][iC] = new TH1D(Form("hSignificance%c%i", kLetter[iS], iC), "; p_{T}(GeV/c); #frac{S}{#sqrt{S+B}}", n_pt_bins, pt_labels.GetArray());
      hChiSquare[iS][iC] = new TH1D(Form("hChiSquare%c%i", kLetter[iS], iC), "; p_{T}(GeV/c); #chi^{2}/NDF", n_pt_bins, pt_labels.GetArray());
      hChiSquareTPC[iS][iC] = new TH1D(Form("hChiSquareTPC%c%i", kLetter[iS], iC), "; p_{T}(GeV/c); #chi^{2}/NDF", n_pt_bins, pt_labels.GetArray());
      hRawCounts[iS][iC] = new TH1D(Form("hRawCounts%c%i", kLetter[iS], iC), "; p_{T}(GeV/c); RawCounts", n_pt_bins, pt_labels.GetArray());
      hRawCountsBinCounting[iS][iC] = new TH1D(Form("hRawCountsBinCounting%c%i", kLetter[iS], iC), "; p_{T}(GeV/c); RawCounts", n_pt_bins, pt_labels.GetArray());
      hSignalGausExpGaus[iS][iC] = new TH1D(Form("hSignalGausExpGaus%c%i", kLetter[iS], iC), "; p_{T}(GeV/c); RawCounts", n_pt_bins, pt_labels.GetArray());
      hWidenRangeSyst[iS][iC] = new TH1D(Form("hWidenRangeSyst%c%i", kLetter[iS], iC), "; p_{T}(GeV/c); RMS", n_pt_bins, pt_labels.GetArray());
      hShiftRangeSyst[iS][iC] = new TH1D(Form("hShiftRangeSyst%c%i", kLetter[iS], iC), "; p_{T}(GeV/c); RMS", n_pt_bins, pt_labels.GetArray());
      hWidenRangeSystTPC[iS][iC] = new TH1D(Form("hWidenRangeSystTPC%c%i", kLetter[iS], iC), "; p_{T}(GeV/c); RMS", n_pt_bins, pt_labels.GetArray());
      hShiftRangeSystTPC[iS][iC] = new TH1D(Form("hShiftRangeSystTPC%c%i", kLetter[iS], iC), "; p_{T}(GeV/c); RMS", n_pt_bins, pt_labels.GetArray());
    }

    float width_range_syst = 0.;
    float width_range_syst_tpc = 0.;
    float pos_range_syst = 0.;
    float pos_range_syst_tpc = 0.;

    for (int iB = 0; iB < n_pt_bins; ++iB)
    {
      if (pt_axis->GetBinCenter(iB + 1) < kPtRange[0] || pt_axis->GetBinCenter(iB + 1) > kPtRange[1])
        continue;
      float sigma_deut[kCentLength];
      float sigma_deut_tpc[kCentLength];
      for (int iS = 0; iS < 2; ++iS)
      {
        int iC = 0;
        // TOF analysis
        if (pt_axis->GetBinCenter(iB + 1) > kCentPtLimits[iC])
          continue;
        TString iTitle = Form("%1.1f #leq #it{p}_{T} < %1.1f GeV/#it{c}", pt_labels[iB], pt_labels[iB + 1]);
        TString iName = Form("d%i_%i", iC, iB);
        TH1D *dat = tof_histo[iS]->ProjectionY(Form("data%i_%i", iC, iB), iB + 1, iB + 1);

        fExpExpTailGaus.mTau0->setVal(-0.3);
        fExpExpTailGaus.mTau0->setVal(0.5);
        RooDataHist data("data", "data", RooArgList(m), Import(*dat));

        /// GausExp
        base_dir->cd(Form("%s/GausExp/C_%d", kNames[iS].data(), iC));

        RooPlot *expExpGausExpGausPlot = fExpExpTailGaus.FitData(dat, iName, iTitle, "Full", "Full", false, -1.2, 1.5);
        fExpExpTailGaus.mSigma->setConstant(false);
        if (iS == 0)
          sigma_deut[iC] = fExpExpTailGaus.mSigma->getVal();
        if (pt_axis->GetBinCenter(iB + 1) > kTOFminPt)
          expExpGausExpGausPlot->Write();
        hSignalGausExpGaus[iS][iC]->SetBinContent(iB + 1, fExpExpTailGaus.mSigCounts->getVal());
        hSignalGausExpGaus[iS][iC]->SetBinError(iB + 1, fExpExpTailGaus.mSigCounts->getError());
        hRawCounts[iS][iC]->SetBinContent(iB + 1, fExpExpTailGaus.mSigCounts->getVal());
        hRawCounts[iS][iC]->SetBinError(iB + 1, fExpExpTailGaus.mSigCounts->getError());

        /// Bin counting TOF
        float residual_vector[n_sigma_vec.size()];
        for (size_t iSigma = 0; iSigma < n_sigma_vec.size(); iSigma++)
        {
          float left_sigma = fExpExpTailGaus.mMu->getVal() - n_sigma_vec[iSigma] * fExpExpTailGaus.mSigma->getVal();
          float right_sigma = fExpExpTailGaus.mMu->getVal() + (float(n_sigma_vec[iSigma]) + 2.) * fExpExpTailGaus.mSigma->getVal();
          int left_edge_bin = dat->FindBin(left_sigma);
          float left_edge_float = dat->GetBinLowEdge(left_edge_bin);
          int right_edge_bin = dat->FindBin(right_sigma);
          float right_edge_float = dat->GetBinLowEdge(right_edge_bin + 1);
          fBkg.mX->setRange("signal", left_edge_float, right_edge_float);
          if (iSigma == 0)
          {
            fBkg.mX->setRange("left", dat->GetXaxis()->GetXmin(), left_edge_float);
            fBkg.mX->setRange("right", right_edge_float, dat->GetXaxis()->GetXmax());
            RooPlot *bkgPlot = fBkg.FitData(dat, Form("%s_sideband", iName.Data()), iTitle, "left,right", "Full");
            base_dir->cd(Form("%s/Sidebands/C_%d", kNames[iS].data(), iC));
            bkgPlot->Write();
          }
          float bkg_integral = (iB > 8) ? fBkg.mBackground->createIntegral(m_bis, NormSet(m_bis), Range("signal"))->getVal() * fBkg.mBkgCounts->getVal() : 0;
          if (iB > 8)
          {
            hChiSquare[iS][iC]->SetBinContent(iB + 1, fBkg.mChi2);
            hChiSquare[iS][iC]->SetBinError(iB + 1, 0);
          }
          float tot_integral = dat->Integral(left_edge_bin, right_edge_bin);
          float sig_integral = tot_integral - bkg_integral;
          float sig_err = TMath::Sqrt(tot_integral + bkg_integral);
          if (iSigma == 0)
          {
            hRawCountsBinCounting[iS][iC]->SetBinContent(iB + 1, sig_integral);
            hRawCountsBinCounting[iS][iC]->SetBinError(iB + 1, sig_err);
            hSignificance[iS][iC]->SetBinContent(iB + 1, sig_integral / TMath::Sqrt(tot_integral));
          }
          residual_vector[iSigma] = sig_integral;
        }
        width_range_syst = TMath::RMS(n_sigma_vec.size(), residual_vector);
        width_range_syst /= hRawCounts[iS][iC]->GetBinContent(iB + 1);
        hWidenRangeSyst[iS][iC]->SetBinContent(iB + 1, width_range_syst);
        // Moving the counting range
        float shift_vector[n_shifts];
        for (int iShift = 0; iShift < n_shifts; iShift++)
        {
          float left_sigma = fExpExpTailGaus.mMu->getVal() - 3. * fExpExpTailGaus.mSigma->getVal() - v_shift[iShift];
          float right_sigma = fExpExpTailGaus.mMu->getVal() + 5. * fExpExpTailGaus.mSigma->getVal() - v_shift[iShift];
          int left_edge_bin = dat->FindBin(left_sigma);
          float left_edge_float = dat->GetBinLowEdge(left_edge_bin);
          int right_edge_bin = dat->FindBin(right_sigma);
          float right_edge_float = dat->GetBinLowEdge(right_edge_bin + 1);
          fBkg.mX->setRange("signal", left_edge_float, right_edge_float);
          float bkg_integral = (iB > 7) ? fBkg.mBackground->createIntegral(m_bis, NormSet(m_bis), Range("signal"))->getVal() * fBkg.mBkgCounts->getVal() : 0;
          float tot_integral = dat->Integral(left_edge_bin, right_edge_bin);
          float sig_integral = tot_integral - bkg_integral;
          float sig_err = TMath::Sqrt(tot_integral + bkg_integral);
          shift_vector[iShift] = sig_integral;
        }
        pos_range_syst = TMath::RMS(n_shifts, shift_vector);
        pos_range_syst /= hRawCounts[iS][iC]->GetBinContent(iB + 1);
        hShiftRangeSyst[iS][iC]->SetBinContent(iB + 1, pos_range_syst);

        /// TPC analysis
        float currentPt = pt_axis->GetBinCenter(iB + 1);
        if (currentPt < kTPCmaxPt)
        {
          base_dir->cd(Form("%s/TPConly", kNames[iS].data()));
          TH1D *tpc_dat = tpc_histo[iS]->ProjectionY(Form("tpc_data%i_%i", iC, iB), iB + 1, iB + 1);
          RooDataHist tpc_data("tpc_data", "tpc_data", RooArgList(ns), Import(*tpc_dat));

          // if (pt_axis->GetBinCenter(iB + 1) >  2.6)
          //   fGausGaus.UseBackground(false);
          // else
          //   fGausGaus.UseBackground(true);
          for (int iT{0}; iT < kNTPCfunctions; ++iT)
          {
            const char* range = iT && currentPt < 1.8 ? "Special" : "Full";
            RooPlot *gausGausPlot = tpcFunctions[iT]->FitData(tpc_dat, Form("TPC_d_%i_%i_%s", iC, iB, kTPCfunctName[iT].data()), iTitle, range, range);
            tpcFunctions[iT]->mSigma->setConstant(false);
            gausGausPlot->Write();

            hTPConly[iS][iC][iT]->SetBinContent(iB + 1, tpcFunctions[iT]->mSigCounts->getVal());
            hTPConly[iS][iC][iT]->SetBinError(iB + 1, tpcFunctions[iT]->mSigCounts->getError());
          }
        }
      }
    }
    for (int iS = 0; iS < 2; ++iS)
    {
      int iC = 0;
      base_dir->cd(Form("%s/GausExp", kNames[iS].data()));
      hRawCounts[iS][iC]->Write();
      hRawCountsBinCounting[iS][iC]->Write();
      hSignalGausExpGaus[iS][iC]->Write();
      base_dir->cd(Form("%s/Systematic", kNames[iS].data()));
      hShiftRangeSyst[iS][iC]->Write();
      hWidenRangeSyst[iS][iC]->Write();
      hWidenRangeSystTPC[iS][iC]->Write();
      hShiftRangeSystTPC[iS][iC]->Write();
      base_dir->cd(Form("%s/Significance", kNames[iS].data()));
      hSignificance[iS][iC]->Write();
      base_dir->cd(Form("%s/TPConly", kNames[iS].data()));
      for (int iT{0}; iT < kNTPCfunctions; ++iT)
      {
        hTPConly[iS][iC][iT]->Write();
      }
      base_dir->cd(Form("%s/ChiSquare", kNames[iS].data()));
      hChiSquare[iS][iC]->Write();
      hChiSquareTPC[iS][iC]->Write();
    }
    base_dir->Close();
  }
  output_file.Close();
}
