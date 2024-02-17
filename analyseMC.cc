#include "ROOT/RDataFrame.hxx"
#include "TF1.h"
#include "TBox.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TStyle.h"
#include <cmath>
#include <iostream>
#include "src/Common.h"

// enum {
//   kDeuteron = BIT(0),
//   kTriton = BIT(1),
//   kHe3 = BIT(2),
//   kHe4 = BIT(3),
//   kHasTOF = BIT(4),
//   kIsReconstructed = BIT(5),
//   kIsAmbiguous = BIT(6), /// just a placeholder now
//   kPositive = BIT(7),
//   kIsPhysicalPrimary = BIT(8), /// MC flags starting from the second half of the short
//   kIsSecondaryFromMaterial = BIT(9),
//   kIsSecondaryFromWeakDecay = BIT(10)
// };

float weight(float pt)
{
  return (5.04194/1.3645054) * pt * std::exp(-pt * 1.35934);
}

void analyseMC(std::string inputFileName = "MC/AO2D_merged_LHC23j6.root", std::string outputFileName = kMCfilename)
{
  gStyle->SetOptStat(0);
  ROOT::EnableImplicitMT();
  ROOT::RDataFrame d("O2nucleitablemc", inputFileName);

  // ROOT::RDataFrame d("DF_1/O2nucleitablemc", inputFileName);
  auto df = d.Define("p", "std::abs(fPt) * sinh(fEta)").Define("pt", "2 * std::abs(fPt)").Define("tofMass", "fBeta < 1.e-3 ? 1.e9 : fBeta >= 1. ? 0 : fTPCInnerParam * 2 * sqrt(1.f / (fBeta * fBeta) - 1.f)").Define("matter", "fPt > 0").Define("nsigmaHe3", nsigmaHe3, {"fTPCInnerParam", "fTPCsignal"}).Define("nsigmaH3", nsigmaH3, {"fTPCInnerParam", "fTPCsignal"}).Define("nsigmaHe4", nsigmaHe4, {"fTPCInnerParam", "fTPCsignal"}).Define("nITSclsIB", "int(0) + bool(fITSclsMap & 1) + bool(fITSclsMap & 2) + bool(fITSclsMap & 4)").Define("nITScls", "nITSclsIB + (fITSclsMap & 8) + (fITSclsMap & 16) + (fITSclsMap & 32) + (fITSclsMap & 64)").Define("hasTOF", "fFlags & (1 << 5)").Define("isReconstructed", "fFlags & (1 << 6)").Define("isPrimary", "fFlags & (1 << 9)").Define("isSecondaryFromMaterial", "fFlags & (1 << 10)").Define("isSecondaryFromWeakDecay", "fFlags & (1 << 11)").Define("gP", "fgPt * std::cosh(fgEta)").Define("gM", "std::abs(fPDGcode) == 1000020030 ? 2.809230089 : (std::abs(fPDGcode) == 1000010030 ? 2.80892 : (std::abs(fPDGcode) == 1000020040 ? 3.72738 : (std::abs(fPDGcode) == 1000010020 ? 1.87561 : 0.1)))").Define("gE", "std::hypot(gM, gP)").Define("gMt", "std::hypot(gM, fgPt)").Define("yMC", "std::asinh(fgPt / gMt * std::sinh(fgEta))").Define("deltaPt", "pt - fgPt").Define("corrPt", "pt + 4.02294e-02 + 1.55116 * std::exp(-1.00391 * pt)").Define("ptWeight", "(5.04194/1.3645054) * fgPt * std::exp(-fgPt * 1.35934)").Define("rapidity", "std::asinh(pt / std::hypot(pt, gM) * std::sinh(fEta))").Define("isHe3", "std::abs(fPDGcode) == 1000020030").Define("isHe4", "std::abs(fPDGcode) == 1000020040"); //
  // double bins[] = {0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.8, 2.2, 2.6, 3.0, 3.5, 4.0};
  // auto hHe3PID = df.Histo2D({"he3PID", ";TPC inner wall #it{p};n#sigma;", 15, bins, 100, -5, 5}, "fTPCInnerParam", "nsigmaH3");
  // hHe3PID->DrawClone("colz");
  // hHe3PID->GetYaxis()->SetRangeUser(-0.05, 0.05);
  // hHe3PID->FitSlicesY();
  // gDirectory->Get("he3PID_1")->DrawClone("same");
  // gDirectory->Get("he3PID_2")->DrawClone("same");
  TBox box;
  box.SetFillStyle(0);
  box.SetLineStyle(kDashed);

  auto dfCutReco = df.Filter("isReconstructed && isPrimary && std::abs(fPDGcode) == 1000020030 && std::abs(fEta) < 0.9 && std::abs(rapidity) < 0.5");
  auto dfCutGen = df.Filter("std::abs(fPDGcode) == 1000020030 && isPrimary && std::abs(yMC) < 0.5");

  auto hDeltaPtHe3 = dfCutReco.Histo2D({"hDeltaPtHe3", ";#it{p}_{T}^{rec} (GeV/#it{c});#it{p}_{T}^{rec}-#it{p}_{T}^{gen} (GeV/#it{c})", 100, 0, 5, 120, -0.4, 0.2}, "pt", "deltaPt");
  auto hDeltaPtCorrHe3 = dfCutReco.Define("deltaPtCorr", "corrPt - fgPt").Histo2D({"hDeltaPtCorrHe3", ";#it{p}_{T}^{rec} (GeV/#it{c});#it{p}_{T}^{rec}-#it{p}_{T}^{gen} (GeV/#it{c})", 100, 0, 5, 100, -0.4, 0.2}, "pt", "deltaPtCorr");
  auto hMomResHe3 = dfCutReco.Define("deltaPtCorr", "corrPt - fgPt").Histo2D({"hMomResHe3", ";#it{p}_{T}^{rec} (GeV/#it{c});#it{p}_{T}^{rec}-#it{p}_{T}^{gen} (GeV/#it{c})", 44, 0.9, 5.3, 80, -0.2, 0.2}, "pt", "deltaPtCorr");

  std::vector<ROOT::RDF::RResultPtr<TH2D>> hDCAxyAHe3, hDCAxyMHe3;
  std::vector<ROOT::RDF::RResultPtr<TH1D>> hRecoTPCAHe3, hRecoTPCMHe3, hRecoTOFAHe3, hRecoTOFMHe3, hGenAHe3, hGenMHe3;
  std::vector<ROOT::RDF::RResultPtr<TH1D>> hRecoTPCAHe3W, hRecoTPCMHe3W, hRecoTOFAHe3W, hRecoTOFMHe3W, hGenAHe3W, hGenMHe3W;
  hRecoTPCAHe3.push_back(dfCutReco.Filter("!matter && fTPCnCls > 120 && nITScls >= 6 && std::abs(fDCAz) < 0.7").Histo1D({"TPCAHe3", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts", kNPtBins, kPtBins}, "corrPt"));
  hRecoTPCMHe3.push_back(dfCutReco.Filter("matter && fTPCnCls > 120 && nITScls >= 6 && std::abs(fDCAz) < 0.7").Histo1D({"TPCMHe3", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts", kNPtBins, kPtBins}, "corrPt"));
  hRecoTOFAHe3.push_back(dfCutReco.Filter("!matter && fTPCnCls > 120 && nITScls >= 6 && std::abs(fDCAz) < 0.7 && hasTOF").Histo1D({"TOFAHe3", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts", kNPtBins, kPtBins}, "corrPt"));
  hRecoTOFMHe3.push_back(dfCutReco.Filter("matter && fTPCnCls > 120 && nITScls >= 6 && std::abs(fDCAz) < 0.7 && hasTOF").Histo1D({"TOFMHe3", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts", kNPtBins, kPtBins}, "corrPt"));
  hGenAHe3.push_back(dfCutGen.Filter("fPDGcode < 0").Histo1D({"genAHe3", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts", kNPtBins, kPtBins}, "fgPt"));
  hGenMHe3.push_back(dfCutGen.Filter("fPDGcode > 0").Histo1D({"genMHe3", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts", kNPtBins, kPtBins}, "fgPt"));

  hRecoTPCAHe3W.push_back(dfCutReco.Filter("!matter && fTPCnCls > 120 && nITScls >= 6 && std::abs(fDCAz) < 0.7").Histo1D({"TPCAHe3", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts", kNPtBins, kPtBins}, "corrPt", "ptWeight"));
  hRecoTPCMHe3W.push_back(dfCutReco.Filter("matter && fTPCnCls > 120 && nITScls >= 6 && std::abs(fDCAz) < 0.7").Histo1D({"TPCMHe3", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts", kNPtBins, kPtBins}, "corrPt", "ptWeight"));
  hRecoTOFAHe3W.push_back(dfCutReco.Filter("!matter && fTPCnCls > 120 && nITScls >= 6 && std::abs(fDCAz) < 0.7 && hasTOF").Histo1D({"TOFAHe3", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts", kNPtBins, kPtBins}, "corrPt", "ptWeight"));
  hRecoTOFMHe3W.push_back(dfCutReco.Filter("matter && fTPCnCls > 120 && nITScls >= 6 && std::abs(fDCAz) < 0.7 && hasTOF").Histo1D({"TOFMHe3", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts", kNPtBins, kPtBins}, "corrPt", "ptWeight"));
  hGenAHe3W.push_back(dfCutGen.Filter("fPDGcode < 0").Histo1D({"genAHe3", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts", kNPtBins, kPtBins}, "fgPt", "ptWeight"));
  hGenMHe3W.push_back(dfCutGen.Filter("fPDGcode > 0").Histo1D({"genMHe3", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts", kNPtBins, kPtBins}, "fgPt", "ptWeight"));


  hDeltaPtHe3->DrawClone("col");
  auto hDeltaPtHe3prof = hDeltaPtHe3->ProfileX();
  hDeltaPtHe3prof->SetLineColor(kRed);
  hDeltaPtHe3prof->DrawClone("same");
  TF1 f("f", "[0] + [2] * TMath::Exp([1] * x)", 0, 5);
  f.SetParameters(-2.98019e-02, -1.31641e+00, -7.66100e-01);
  hDeltaPtHe3prof->Fit(&f, "NR", "", 1, 5);
  f.DrawClone("same");
  new TCanvas;
  hDeltaPtCorrHe3->DrawClone("col");
  auto hDeltaPtCorrHe3prof = hDeltaPtCorrHe3->ProfileX();
  hDeltaPtCorrHe3prof->SetLineColor(kRed);
  hDeltaPtCorrHe3prof->DrawClone("same");

  new TCanvas;
  TObjArray aSlices;
  hMomResHe3->FitSlicesY(nullptr, 0, -1, 0, "QLNRG3", &aSlices);
  TH1 *_hMomRes = static_cast<TH1 *>(aSlices[2]);
  _hMomRes->GetYaxis()->SetTitle("#sigma_{p_{T}} (GeV/#it{c})");
  aSlices[2]->DrawClone();

  new TCanvas;
  auto hEffA = (TH1 *)(hRecoTPCAHe3[0]->Clone("hEffA"));
  hEffA->Divide(hRecoTPCAHe3[0].GetPtr(), hGenAHe3[0].GetPtr(), 1., 1., "B");
  hEffA->Draw();

  new TCanvas;
  auto hEffM = (TH1 *)(hRecoTPCMHe3[0]->Clone("hEffM"));
  hEffM->Divide(hRecoTPCMHe3[0].GetPtr(), hGenMHe3[0].GetPtr(), 1., 1., "B");
  hEffM->Draw();

  new TCanvas;
  auto hGenRap = dfCutGen.Filter("fPDGcode < 0").Histo1D({"hGenRap", ";y;Counts", 40, -1, 1}, "yMC");
  auto hGenEta = dfCutGen.Filter("fPDGcode < 0").Histo1D({"hGenEta", ";#eta;Counts", 40, -1, 1}, "fgEta");
  hGenRap->DrawClone();
  hGenEta->SetLineColor(kRed);
  hGenEta->DrawClone("same");

  int iTrial{0};
  size_t nTrials{0};
  //size_t nTrials{kCutNames["fDCAz"].size() * kCutNames["fTPCnCls"].size() * kCutNames["nITScls"].size()};
  /*
  for (size_t iDCAz{0}; iDCAz < kCutNames["fDCAz"].size(); ++iDCAz)
  {
    auto dfDCAz = dfCutReco.Filter("std::abs(fDCAz) < " + std::to_string(kCutNames["fDCAz"][iDCAz]));
    for (size_t iTPCcls{0}; iTPCcls < kCutNames["fTPCnCls"].size(); ++iTPCcls)
    {
      auto dfTPCcls = dfDCAz.Filter("fTPCnCls > " + std::to_string(kCutNames["fTPCnCls"][iTPCcls]));
      for (size_t iITScls{0}; iITScls < kCutNames["nITScls"].size(); ++iITScls)
      {
        auto dfITScls = dfTPCcls.Filter("nITScls >= " + std::to_string(kCutNames["nITScls"][iITScls]));
        hRecoTPCAHe3.push_back(dfCutReco.Filter("!matter").Histo1D({Form("TPCAHe3%i", iTrial), ";#it{p}_{T}^{rec} (GeV/#it{c});Counts", kNPtBins, kPtBins}, "corrPt"));
        hRecoTPCMHe3.push_back(dfCutReco.Filter("matter").Histo1D({Form("TPCMHe3%i", iTrial), ";#it{p}_{T}^{rec} (GeV/#it{c});Counts", kNPtBins, kPtBins}, "corrPt"));
        hRecoTOFAHe3.push_back(dfCutReco.Filter("!matter && hasTOF").Histo1D({Form("TOFAHe3%i", iTrial), ";#it{p}_{T}^{rec} (GeV/#it{c});Counts", kNPtBins, kPtBins}, "corrPt"));
        hRecoTOFMHe3.push_back(dfCutReco.Filter("matter && hasTOF").Histo1D({Form("TOFMHe3%i", iTrial), ";#it{p}_{T}^{rec} (GeV/#it{c});Counts", kNPtBins, kPtBins}, "corrPt"));

        hRecoTPCAHe3W.push_back(dfCutReco.Filter("!matter").Histo1D({Form("TPCAHe3%i", iTrial), ";#it{p}_{T}^{rec} (GeV/#it{c});Counts", kNPtBins, kPtBins}, "corrPt", "ptWeight"));
        hRecoTPCMHe3W.push_back(dfCutReco.Filter("matter").Histo1D({Form("TPCMHe3%i", iTrial), ";#it{p}_{T}^{rec} (GeV/#it{c});Counts", kNPtBins, kPtBins}, "corrPt", "ptWeight"));
        hRecoTOFAHe3W.push_back(dfCutReco.Filter("!matter && hasTOF").Histo1D({Form("TOFAHe3%i", iTrial), ";#it{p}_{T}^{rec} (GeV/#it{c});Counts", kNPtBins, kPtBins}, "corrPt", "ptWeight"));
        hRecoTOFMHe3W.push_back(dfCutReco.Filter("matter && hasTOF").Histo1D({Form("TOFMHe3%i", iTrial), ";#it{p}_{T}^{rec} (GeV/#it{c});Counts", kNPtBins, kPtBins}, "corrPt", "ptWeight"));
        iTrial++;
      }
    }
  }*/

  TFile outputFile(outputFileName.data(), "recreate");
  auto dir = outputFile.mkdir("nuclei");
  dir->cd();
  hGenAHe3[0]->Write("genAHe3");
  hGenMHe3[0]->Write("genMHe3");
  hRecoTPCAHe3[0]->Write("TPCAHe3");
  hRecoTPCMHe3[0]->Write("TPCMHe3");
  hRecoTOFAHe3[0]->Write("TOFAHe3");
  hRecoTOFMHe3[0]->Write("TOFMHe3");

  hGenAHe3W[0]->Write("genAHe3W");
  hGenMHe3W[0]->Write("genMHe3W");
  hRecoTPCAHe3W[0]->Write("TPCAHe3W");
  hRecoTPCMHe3W[0]->Write("TPCMHe3W");
  hRecoTOFAHe3W[0]->Write("TOFAHe3W");
  hRecoTOFMHe3W[0]->Write("TOFMHe3W");

  auto effTPCA = (TH1 *)hRecoTPCAHe3[0]->Clone(Form("effTPCA"));
  auto effTPCM = (TH1 *)hRecoTPCMHe3[0]->Clone(Form("effTPCM"));
  auto effTOFA = (TH1 *)hRecoTOFAHe3[0]->Clone(Form("effTOFA"));
  auto effTOFM = (TH1 *)hRecoTOFMHe3[0]->Clone(Form("effTOFM"));
  effTPCA->Divide(hGenAHe3[0].GetPtr());
  effTPCM->Divide(hGenMHe3[0].GetPtr());
  effTOFA->Divide(hGenAHe3[0].GetPtr());
  effTOFM->Divide(hGenMHe3[0].GetPtr());
  effTPCA->Write("effTPCA");
  effTPCM->Write("effTPCM");
  effTOFA->Write("effTOFA");
  effTOFM->Write("effTOFM");

  auto effWTPCA = (TH1 *)hRecoTPCAHe3W[0]->Clone(Form("effWTPCA"));
  auto effWTPCM = (TH1 *)hRecoTPCMHe3W[0]->Clone(Form("effWTPCM"));
  auto effWTOFA = (TH1 *)hRecoTOFAHe3W[0]->Clone(Form("effWTOFA"));
  auto effWTOFM = (TH1 *)hRecoTOFMHe3W[0]->Clone(Form("effWTOFM"));
  effWTPCA->Divide(hGenAHe3W[0].GetPtr());
  effWTPCM->Divide(hGenMHe3W[0].GetPtr());
  effWTOFA->Divide(hGenAHe3W[0].GetPtr());
  effWTOFM->Divide(hGenMHe3W[0].GetPtr());
  effWTPCA->Write("effWTPCA");
  effWTPCM->Write("effWTPCM");
  effWTOFA->Write("effWTOFA");
  effWTOFM->Write("effWTOFM");

  for (size_t iT{0}; iT < nTrials; ++iT)
  {
    auto dir = outputFile.mkdir(Form("nuclei%zu", iT));
    dir->cd();
    hGenAHe3[0]->Write("genAHe3");
    hGenMHe3[0]->Write("genMHe3");
    hRecoTPCAHe3[iT + 1]->Write("TPCAHe3");
    hRecoTPCMHe3[iT + 1]->Write("TPCMHe3");
    hRecoTOFAHe3[iT + 1]->Write("TOFAHe3");
    hRecoTOFMHe3[iT + 1]->Write("TOFMHe3");

    hGenAHe3W[0]->Write("genAHe3W");
    hGenMHe3W[0]->Write("genMHe3W");
    hRecoTPCAHe3W[iT + 1]->Write("TPCAHe3W");
    hRecoTPCMHe3W[iT + 1]->Write("TPCMHe3W");
    hRecoTOFAHe3W[iT + 1]->Write("TOFAHe3W");
    hRecoTOFMHe3W[iT + 1]->Write("TOFMHe3W");

    auto effTPCA = (TH1 *)hRecoTPCAHe3[iT + 1]->Clone(Form("effTPCA%i", iT));
    auto effTPCM = (TH1 *)hRecoTPCMHe3[iT + 1]->Clone(Form("effTPCM%i", iT));
    auto effTOFA = (TH1 *)hRecoTOFAHe3[iT + 1]->Clone(Form("effTOFA%i", iT));
    auto effTOFM = (TH1 *)hRecoTOFMHe3[iT + 1]->Clone(Form("effTOFM%i", iT));
    effTPCA->Divide(hGenAHe3[0].GetPtr());
    effTPCM->Divide(hGenMHe3[0].GetPtr());
    effTPCA->SetLineColor(kRed);
    effTPCM->SetLineColor(kRed);
    effTOFA->Divide(hGenAHe3[0].GetPtr());
    effTOFM->SetLineColor(kBlue);
    effTOFA->SetLineColor(kBlue);
    effTOFM->Divide(hGenMHe3[0].GetPtr());
    effTPCA->GetYaxis()->SetTitle("Efficiency #times Acceptance");
    effTPCM->GetYaxis()->SetTitle("Efficiency #times Acceptance");
    effTOFA->GetYaxis()->SetTitle("Efficiency #times Acceptance");
    effTOFM->GetYaxis()->SetTitle("Efficiency #times Acceptance");
    effTPCA->Write("effTPCA");
    effTPCM->Write("effTPCM");
    effTOFA->Write("effTOFA");
    effTOFM->Write("effTOFM");
    auto matchingTOFA = (TH1 *)hRecoTOFAHe3[iT + 1]->Clone(Form("matchingTOFA%i", iT));
    auto matchingTOFM = (TH1 *)hRecoTOFMHe3[iT + 1]->Clone(Form("matchingTOFM%i", iT));
    matchingTOFA->Divide(effTPCA);
    matchingTOFM->Divide(effTPCM);
    matchingTOFA->Write();
    matchingTOFM->Write();

    auto effWTPCA = (TH1 *)hRecoTPCAHe3W[iT + 1]->Clone(Form("effWTPCA%i", iT));
    auto effWTPCM = (TH1 *)hRecoTPCMHe3W[iT + 1]->Clone(Form("effWTPCM%i", iT));
    auto effWTOFA = (TH1 *)hRecoTOFAHe3W[iT + 1]->Clone(Form("effWTOFA%i", iT));
    auto effWTOFM = (TH1 *)hRecoTOFMHe3W[iT + 1]->Clone(Form("effWTOFM%i", iT));
    effWTPCA->Divide(hGenAHe3W[0].GetPtr());
    effWTPCM->Divide(hGenMHe3W[0].GetPtr());
    effWTPCA->SetLineColor(kRed);
    effWTPCM->SetLineColor(kRed);
    effWTOFA->Divide(hGenAHe3W[0].GetPtr());
    effWTOFM->SetLineColor(kBlue);
    effWTOFA->SetLineColor(kBlue);
    effWTOFM->Divide(hGenMHe3W[0].GetPtr());
    effWTPCA->GetYaxis()->SetTitle("Efficiency #times Acceptance");
    effWTPCM->GetYaxis()->SetTitle("Efficiency #times Acceptance");
    effWTOFA->GetYaxis()->SetTitle("Efficiency #times Acceptance");
    effWTOFM->GetYaxis()->SetTitle("Efficiency #times Acceptance");
    effWTPCA->Write("effWTPCA");
    effWTPCM->Write("effWTPCM");
    effWTOFA->Write("effWTOFA");
    effWTOFM->Write("effWTOFM");
    auto matchingWTOFA = (TH1 *)hRecoTOFAHe3W[iT + 1]->Clone(Form("matchingWTOFA%i", iT));
    auto matchingWTOFM = (TH1 *)hRecoTOFMHe3W[iT + 1]->Clone(Form("matchingWTOFM%i", iT));
    matchingWTOFA->Divide(effWTPCA);
    matchingWTOFM->Divide(effWTPCM);
    matchingWTOFA->Write();
    matchingWTOFM->Write();


  }
}
