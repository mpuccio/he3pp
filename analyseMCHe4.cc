#include "ROOT/RDataFrame.hxx"
#include "TF1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TStyle.h"
#include <cmath>
#include <iostream>
#include "src/Common.h"


float weight(float pt)
{
  return (5.04194/1.3645054) * pt * std::exp(-pt * 1.35934);
}

void analyseMCHe4(std::string inputFileName = kMCtreeFilename, std::string outputFileName = kMCfilenameHe4, bool enableTrials = false)
{
  inputFileName = gSystem->ExpandPathName(inputFileName.data());
  outputFileName = gSystem->ExpandPathName(outputFileName.data());
  gStyle->SetOptStat(0);
  ROOT::EnableImplicitMT();
  ROOT::RDataFrame d("O2nucleitablemc", inputFileName);
  auto dfData = defineColumnsForData(d);
  auto df = dfData.Define("gP", "fgPt * std::cosh(fgEta)")
                  .Define("gM", "std::abs(fPDGcode) == 1000020030 ? 2.809230089 : (std::abs(fPDGcode) == 1000010030 ? 2.80892 : (std::abs(fPDGcode) == 1000020040 ? 3.72738 : (std::abs(fPDGcode) == 1000010020 ? 1.87561 : 0.1)))")
                  .Define("gE", "std::hypot(gM, gP)")
                  .Define("gMt", "std::hypot(gM, fgPt)")
                  .Define("yMC", "std::asinh(fgPt / gMt * std::sinh(fgEta))")
                  .Define("deltaPtUncorrected", "ptUncorr - fgPt")
                  .Define("deltaPt", "ptHe4 - fgPt")
                  .Define("ptWeight", "(5.04194/1.3645054) * fgPt * std::exp(-fgPt * 1.35934)")
                  .Define("rapidity", "std::asinh(pt / std::hypot(pt, gM) * std::sinh(fEta))")
                  .Define("isHe3", "std::abs(fPDGcode) == 1000020030")
                  .Define("isHe4", "std::abs(fPDGcode) == 1000020040")
                  .Filter("isHe4 && isPrimary"); //

  auto dfCutReco = df.Filter("nITScls > 4 && fTPCnCls > 110 && std::abs(fEta) < 0.9 && std::abs(rapidity) < 0.5");
  auto dfCutGen = df.Filter("std::abs(yMC) < 0.5");

  auto hDeltaPtHe4 = dfCutReco.Histo2D({"hDeltaPtHe4", ";#it{p}_{T}^{rec} (GeV/#it{c});#it{p}_{T}^{rec}-#it{p}_{T}^{gen} (GeV/#it{c})", 100, 0, 5, 120, -0.4, 0.2}, "pt", "deltaPtUncorrected");
  auto hDeltaPtCorrHe4 = dfCutReco.Histo2D({"hDeltaPtCorrHe4", ";#it{p}_{T}^{rec} (GeV/#it{c});#it{p}_{T}^{rec}-#it{p}_{T}^{gen} (GeV/#it{c})", 100, 0, 5, 100, -0.4, 0.2}, "pt", "deltaPt");
  auto hMomResHe4 = dfCutReco.Histo2D({"hMomResHe4", ";#it{p}_{T}^{rec} (GeV/#it{c});#it{p}_{T}^{rec}-#it{p}_{T}^{gen} (GeV/#it{c})", 44, 0.9, 5.3, 80, -0.2, 0.2}, "pt", "deltaPt");

  std::vector<ROOT::RDF::RResultPtr<TH2D>> hDCAxyAHe4, hDCAxyMHe4;
  std::vector<ROOT::RDF::RResultPtr<TH1D>> hRecoTPCAHe4, hRecoTPCMHe4, hRecoTOFAHe4, hRecoTOFMHe4, hGenAHe4, hGenMHe4;
  std::vector<ROOT::RDF::RResultPtr<TH1D>> hRecoTPCAHe4W, hRecoTPCMHe4W, hRecoTOFAHe4W, hRecoTOFMHe4W, hGenAHe4W, hGenMHe4W;
  hRecoTPCAHe4.push_back(dfCutReco.Filter("!matter && fTPCnCls > 120 && nITScls >= 6 && std::abs(fDCAz) < 0.7").Histo1D({"TPCAHe4", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts", kNPtBins, kPtBins}, "pt"));
  hRecoTPCMHe4.push_back(dfCutReco.Filter("matter && fTPCnCls > 120 && nITScls >= 6 && std::abs(fDCAz) < 0.7").Histo1D({"TPCMHe4", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts", kNPtBins, kPtBins}, "pt"));
  hRecoTOFAHe4.push_back(dfCutReco.Filter("!matter && fTPCnCls > 120 && nITScls >= 6 && std::abs(fDCAz) < 0.7 && hasTOF").Histo1D({"TOFAHe4", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts", kNPtBins, kPtBins}, "pt"));
  hRecoTOFMHe4.push_back(dfCutReco.Filter("matter && fTPCnCls > 120 && nITScls >= 6 && std::abs(fDCAz) < 0.7 && hasTOF").Histo1D({"TOFMHe4", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts", kNPtBins, kPtBins}, "pt"));
  hGenAHe4.push_back(dfCutGen.Filter("fPDGcode < 0").Histo1D({"genAHe4", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts", kNPtBins, kPtBins}, "fgPt"));
  hGenMHe4.push_back(dfCutGen.Filter("fPDGcode > 0").Histo1D({"genMHe4", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts", kNPtBins, kPtBins}, "fgPt"));

  hRecoTPCAHe4W.push_back(dfCutReco.Filter("!matter && fTPCnCls > 120 && nITScls >= 6 && std::abs(fDCAz) < 0.7").Histo1D({"TPCAHe4", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts", kNPtBins, kPtBins}, "pt", "ptWeight"));
  hRecoTPCMHe4W.push_back(dfCutReco.Filter("matter && fTPCnCls > 120 && nITScls >= 6 && std::abs(fDCAz) < 0.7").Histo1D({"TPCMHe4", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts", kNPtBins, kPtBins}, "pt", "ptWeight"));
  hRecoTOFAHe4W.push_back(dfCutReco.Filter("!matter && fTPCnCls > 120 && nITScls >= 6 && std::abs(fDCAz) < 0.7 && hasTOF").Histo1D({"TOFAHe4", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts", kNPtBins, kPtBins}, "pt", "ptWeight"));
  hRecoTOFMHe4W.push_back(dfCutReco.Filter("matter && fTPCnCls > 120 && nITScls >= 6 && std::abs(fDCAz) < 0.7 && hasTOF").Histo1D({"TOFMHe4", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts", kNPtBins, kPtBins}, "pt", "ptWeight"));
  hGenAHe4W.push_back(dfCutGen.Filter("fPDGcode < 0").Histo1D({"genAHe4", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts", kNPtBins, kPtBins}, "fgPt", "ptWeight"));
  hGenMHe4W.push_back(dfCutGen.Filter("fPDGcode > 0").Histo1D({"genMHe4", ";#it{p}_{T}^{rec} (GeV/#it{c});Counts", kNPtBins, kPtBins}, "fgPt", "ptWeight"));


  hDeltaPtHe4->DrawClone("col");
  auto hDeltaPtHe4prof = hDeltaPtHe4->ProfileX();
  hDeltaPtHe4prof->SetLineColor(kRed);
  hDeltaPtHe4prof->DrawClone("same");
  TF1 f("f", "[0] + [2] * TMath::Exp([1] * x)", 0, 5);
  f.SetParameters(-0.0419608, -1.4019, -1.75861);
  hDeltaPtHe4prof->Fit(&f, "NR", "", 1, 5);
  f.DrawClone("same");

  new TCanvas("hDeltaPtCorrHe4", "hDeltaPtCorrHe4");
  hDeltaPtCorrHe4->DrawClone("col");
  auto hDeltaPtCorrHe4prof = hDeltaPtCorrHe4->ProfileX();
  hDeltaPtCorrHe4prof->SetLineColor(kRed);
  hDeltaPtCorrHe4prof->DrawClone("same");

  new TCanvas("hMomResHe4", "hMomResHe4");
  TObjArray aSlices;
  hMomResHe4->FitSlicesY(nullptr, 0, -1, 0, "QLNRG3", &aSlices);
  TH1 *_hMomRes = static_cast<TH1 *>(aSlices[2]);
  _hMomRes->GetYaxis()->SetTitle("#sigma_{p_{T}} (GeV/#it{c})");
  aSlices[2]->DrawClone();

  new TCanvas("effMatterAntiMatter", "effMatterAntiMatter");
  auto hEffA = (TH1 *)(hRecoTPCAHe4[0]->Clone("hEffA"));
  hEffA->Divide(hRecoTPCAHe4[0].GetPtr(), hGenAHe4[0].GetPtr(), 1., 1., "B");
  hEffA->SetLineColor(kRed);
  hEffA->Draw();
  auto hEffM = (TH1 *)(hRecoTPCMHe4[0]->Clone("hEffM"));
  hEffM->Divide(hRecoTPCMHe4[0].GetPtr(), hGenMHe4[0].GetPtr(), 1., 1., "B");
  hEffM->Draw("same");

  new TCanvas;
  auto hGenRap = dfCutGen.Filter("fPDGcode < 0").Histo1D({"hGenRap", ";y;Counts", 40, -1, 1}, "yMC");
  auto hGenEta = dfCutGen.Filter("fPDGcode < 0").Histo1D({"hGenEta", ";#eta;Counts", 40, -1, 1}, "fgEta");
  hGenRap->DrawClone();
  hGenEta->SetLineColor(kRed);
  hGenEta->DrawClone("same");

  int iTrial{0};
  size_t nTrials{enableTrials ? kCutNames["nsigmaDCAz"].size() * kCutNames["fTPCnCls"].size() * kCutNames["nITScls"].size() : 0};
  for (size_t iDCAz{0}; enableTrials && iDCAz < kCutNames["nsigmaDCAz"].size(); ++iDCAz)
  {
    auto dnsigmaDCAz = dfCutReco.Filter("std::abs(nsigmaDCAz) < " + std::to_string(kCutNames["nsigmaDCAz"][iDCAz]));
    for (size_t iTPCcls{0}; iTPCcls < kCutNames["fTPCnCls"].size(); ++iTPCcls)
    {
      auto dfTPCcls = dnsigmaDCAz.Filter("fTPCnCls > " + std::to_string(kCutNames["fTPCnCls"][iTPCcls]));
      for (size_t iITScls{0}; iITScls < kCutNames["nITScls"].size(); ++iITScls)
      {
        auto dfITScls = dfTPCcls.Filter("nITScls >= " + std::to_string(kCutNames["nITScls"][iITScls]));
        hRecoTPCAHe4.push_back(dfCutReco.Filter("!matter").Histo1D({Form("TPCAHe4%i", iTrial), ";#it{p}_{T}^{rec} (GeV/#it{c});Counts", kNPtBins, kPtBins}, "pt"));
        hRecoTPCMHe4.push_back(dfCutReco.Filter("matter").Histo1D({Form("TPCMHe4%i", iTrial), ";#it{p}_{T}^{rec} (GeV/#it{c});Counts", kNPtBins, kPtBins}, "pt"));
        hRecoTOFAHe4.push_back(dfCutReco.Filter("!matter && hasTOF").Histo1D({Form("TOFAHe4%i", iTrial), ";#it{p}_{T}^{rec} (GeV/#it{c});Counts", kNPtBins, kPtBins}, "pt"));
        hRecoTOFMHe4.push_back(dfCutReco.Filter("matter && hasTOF").Histo1D({Form("TOFMHe4%i", iTrial), ";#it{p}_{T}^{rec} (GeV/#it{c});Counts", kNPtBins, kPtBins}, "pt"));

        hRecoTPCAHe4W.push_back(dfCutReco.Filter("!matter").Histo1D({Form("TPCAHe4%i", iTrial), ";#it{p}_{T}^{rec} (GeV/#it{c});Counts", kNPtBins, kPtBins}, "pt", "ptWeight"));
        hRecoTPCMHe4W.push_back(dfCutReco.Filter("matter").Histo1D({Form("TPCMHe4%i", iTrial), ";#it{p}_{T}^{rec} (GeV/#it{c});Counts", kNPtBins, kPtBins}, "pt", "ptWeight"));
        hRecoTOFAHe4W.push_back(dfCutReco.Filter("!matter && hasTOF").Histo1D({Form("TOFAHe4%i", iTrial), ";#it{p}_{T}^{rec} (GeV/#it{c});Counts", kNPtBins, kPtBins}, "pt", "ptWeight"));
        hRecoTOFMHe4W.push_back(dfCutReco.Filter("matter && hasTOF").Histo1D({Form("TOFMHe4%i", iTrial), ";#it{p}_{T}^{rec} (GeV/#it{c});Counts", kNPtBins, kPtBins}, "pt", "ptWeight"));
        iTrial++;
      }
    }
  }

  TFile outputFile(outputFileName.data(), "recreate");
  auto dir = outputFile.mkdir("nuclei");
  dir->cd();
  hGenAHe4[0]->Write("genAHe4");
  hGenMHe4[0]->Write("genMHe4");
  hRecoTPCAHe4[0]->Write("TPCAHe4");
  hRecoTPCMHe4[0]->Write("TPCMHe4");
  hRecoTOFAHe4[0]->Write("TOFAHe4");
  hRecoTOFMHe4[0]->Write("TOFMHe4");

  hGenAHe4W[0]->Write("genAHe4W");
  hGenMHe4W[0]->Write("genMHe4W");
  hRecoTPCAHe4W[0]->Write("TPCAHe4W");
  hRecoTPCMHe4W[0]->Write("TPCMHe4W");
  hRecoTOFAHe4W[0]->Write("TOFAHe4W");
  hRecoTOFMHe4W[0]->Write("TOFMHe4W");

  auto effTPCA = (TH1 *)hRecoTPCAHe4[0]->Clone(Form("effTPCA"));
  auto effTPCM = (TH1 *)hRecoTPCMHe4[0]->Clone(Form("effTPCM"));
  auto effTOFA = (TH1 *)hRecoTOFAHe4[0]->Clone(Form("effTOFA"));
  auto effTOFM = (TH1 *)hRecoTOFMHe4[0]->Clone(Form("effTOFM"));
  effTPCA->Divide(hGenAHe4[0].GetPtr());
  effTPCM->Divide(hGenMHe4[0].GetPtr());
  effTOFA->Divide(hGenAHe4[0].GetPtr());
  effTOFM->Divide(hGenMHe4[0].GetPtr());
  effTPCA->Write("effTPCA");
  effTPCM->Write("effTPCM");
  effTOFA->Write("effTOFA");
  effTOFM->Write("effTOFM");

  auto effWTPCA = (TH1 *)hRecoTPCAHe4W[0]->Clone(Form("effWTPCA"));
  auto effWTPCM = (TH1 *)hRecoTPCMHe4W[0]->Clone(Form("effWTPCM"));
  auto effWTOFA = (TH1 *)hRecoTOFAHe4W[0]->Clone(Form("effWTOFA"));
  auto effWTOFM = (TH1 *)hRecoTOFMHe4W[0]->Clone(Form("effWTOFM"));
  effWTPCA->Divide(hGenAHe4W[0].GetPtr());
  effWTPCM->Divide(hGenMHe4W[0].GetPtr());
  effWTOFA->Divide(hGenAHe4W[0].GetPtr());
  effWTOFM->Divide(hGenMHe4W[0].GetPtr());
  effWTPCA->Write("effWTPCA");
  effWTPCM->Write("effWTPCM");
  effWTOFA->Write("effWTOFA");
  effWTOFM->Write("effWTOFM");

  for (int iT{0}; iT < static_cast<int>(nTrials); ++iT)
  {
    auto dir = outputFile.mkdir(Form("nuclei%i", iT));
    dir->cd();
    hGenAHe4[0]->Write("genAHe4");
    hGenMHe4[0]->Write("genMHe4");
    hRecoTPCAHe4[iT + 1]->Write("TPCAHe4");
    hRecoTPCMHe4[iT + 1]->Write("TPCMHe4");
    hRecoTOFAHe4[iT + 1]->Write("TOFAHe4");
    hRecoTOFMHe4[iT + 1]->Write("TOFMHe4");

    hGenAHe4W[0]->Write("genAHe4W");
    hGenMHe4W[0]->Write("genMHe4W");
    hRecoTPCAHe4W[iT + 1]->Write("TPCAHe4W");
    hRecoTPCMHe4W[iT + 1]->Write("TPCMHe4W");
    hRecoTOFAHe4W[iT + 1]->Write("TOFAHe4W");
    hRecoTOFMHe4W[iT + 1]->Write("TOFMHe4W");

    auto effTPCA = (TH1 *)hRecoTPCAHe4[iT + 1]->Clone(Form("effTPCA%i", iT));
    auto effTPCM = (TH1 *)hRecoTPCMHe4[iT + 1]->Clone(Form("effTPCM%i", iT));
    auto effTOFA = (TH1 *)hRecoTOFAHe4[iT + 1]->Clone(Form("effTOFA%i", iT));
    auto effTOFM = (TH1 *)hRecoTOFMHe4[iT + 1]->Clone(Form("effTOFM%i", iT));
    effTPCA->Divide(hGenAHe4[0].GetPtr());
    effTPCM->Divide(hGenMHe4[0].GetPtr());
    effTPCA->SetLineColor(kRed);
    effTPCM->SetLineColor(kRed);
    effTOFA->Divide(hGenAHe4[0].GetPtr());
    effTOFM->SetLineColor(kBlue);
    effTOFA->SetLineColor(kBlue);
    effTOFM->Divide(hGenMHe4[0].GetPtr());
    effTPCA->GetYaxis()->SetTitle("Efficiency #times Acceptance");
    effTPCM->GetYaxis()->SetTitle("Efficiency #times Acceptance");
    effTOFA->GetYaxis()->SetTitle("Efficiency #times Acceptance");
    effTOFM->GetYaxis()->SetTitle("Efficiency #times Acceptance");
    effTPCA->Write("effTPCA");
    effTPCM->Write("effTPCM");
    effTOFA->Write("effTOFA");
    effTOFM->Write("effTOFM");
    auto matchingTOFA = (TH1 *)hRecoTOFAHe4[iT + 1]->Clone(Form("matchingTOFA%i", iT));
    auto matchingTOFM = (TH1 *)hRecoTOFMHe4[iT + 1]->Clone(Form("matchingTOFM%i", iT));
    matchingTOFA->Divide(effTPCA);
    matchingTOFM->Divide(effTPCM);
    matchingTOFA->Write();
    matchingTOFM->Write();

    auto effWTPCA = (TH1 *)hRecoTPCAHe4W[iT + 1]->Clone(Form("effWTPCA%i", iT));
    auto effWTPCM = (TH1 *)hRecoTPCMHe4W[iT + 1]->Clone(Form("effWTPCM%i", iT));
    auto effWTOFA = (TH1 *)hRecoTOFAHe4W[iT + 1]->Clone(Form("effWTOFA%i", iT));
    auto effWTOFM = (TH1 *)hRecoTOFMHe4W[iT + 1]->Clone(Form("effWTOFM%i", iT));
    effWTPCA->Divide(hGenAHe4W[0].GetPtr());
    effWTPCM->Divide(hGenMHe4W[0].GetPtr());
    effWTPCA->SetLineColor(kRed);
    effWTPCM->SetLineColor(kRed);
    effWTOFA->Divide(hGenAHe4W[0].GetPtr());
    effWTOFM->SetLineColor(kBlue);
    effWTOFA->SetLineColor(kBlue);
    effWTOFM->Divide(hGenMHe4W[0].GetPtr());
    effWTPCA->GetYaxis()->SetTitle("Efficiency #times Acceptance");
    effWTPCM->GetYaxis()->SetTitle("Efficiency #times Acceptance");
    effWTOFA->GetYaxis()->SetTitle("Efficiency #times Acceptance");
    effWTOFM->GetYaxis()->SetTitle("Efficiency #times Acceptance");
    effWTPCA->Write("effWTPCA");
    effWTPCM->Write("effWTPCM");
    effWTOFA->Write("effWTOFA");
    effWTOFM->Write("effWTOFM");
    auto matchingWTOFA = (TH1 *)hRecoTOFAHe4W[iT + 1]->Clone(Form("matchingWTOFA%i", iT));
    auto matchingWTOFM = (TH1 *)hRecoTOFMHe4W[iT + 1]->Clone(Form("matchingWTOFM%i", iT));
    matchingWTOFA->Divide(effWTPCA);
    matchingWTOFM->Divide(effWTPCM);
    matchingWTOFA->Write();
    matchingWTOFM->Write();


  }
}
