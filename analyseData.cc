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

constexpr int nDCAbins = 37;
constexpr double DCAbinning[]{-0.2,-0.15,-0.125,-0.1,-0.075,-0.05,-0.025,-0.023,-0.021,-0.019,-0.017,-0.015,-0.013,-0.011,-0.009,-0.007,-0.005,-0.003,-0.001,0.001,0.003,0.005,0.007,0.009,0.011,0.013,0.015,0.017,0.019,0.021,0.023, 0.025, 0.05,0.075,0.1,0.125,0.15,0.2};

void analyseData(std::string inputFileName = "data/MergedTrees.root", std::string outputFileName = "data/dataSmall.root", bool skim = true)
{
  gStyle->SetOptStat(0);
  ROOT::EnableImplicitMT();
  ROOT::RDataFrame d("O2nucleitable", inputFileName);
  auto df = d.Define("ptUncorr", "2 * std::abs(fPt)").Define("pt", "ptUncorr + 2.98019e-02 + 7.66100e-01 * std::exp(-1.31641e+00 * ptUncorr)").Define("p", "pt * cosh(fEta)").Define("tofMass", "fBeta < 1.e-3 ? 1.e9 : fBeta >= 1. ? 0 : fTPCInnerParam * 2 * sqrt(1.f / (fBeta * fBeta) - 1.f)").Define("matter", "fPt > 0").Define("nsigmaHe3", nsigmaHe3, {"fTPCInnerParam", "fTPCsignal"}).Define("nsigmaH3", nsigmaH3, {"fTPCInnerParam", "fTPCsignal"}).Define("nsigmaHe4", nsigmaHe4, {"fTPCInnerParam", "fTPCsignal"}).Define("nITSclsIB", "int(0) + bool(fITSclsMap & 1) + bool(fITSclsMap & 2) + bool(fITSclsMap & 4)").Define("nITScls", "nITSclsIB + bool(fITSclsMap & 8) + bool(fITSclsMap & 16) + bool(fITSclsMap & 32) + bool(fITSclsMap & 64)").Define("hasTOF", "fFlags & (1 << 4)").Filter("abs(fEta) < 0.9").Define("isReconstructed", "fFlags & (1 << 5)").Define("isPrimary", "fFlags & (1 << 8)").Define("isSecondaryFromMaterial", "fFlags & (1 << 9)").Define("isSecondaryFromWeakDecay", "fFlags & (1 << 10)").Define("deltaMass", "tofMass - 2.80839").Define("nsigmaDCAxy", nSigmaDCAxy, {"pt", "fDCAxy"}).Define("nsigmaDCAz", nSigmaDCAz, {"pt", "fDCAz"});

  auto dfBase = df.Filter("isReconstructed && std::abs(fEta) < 0.9");

  auto dfPrimary = dfBase.Filter("fTPCnCls > 120 && nITScls >= 6 && std::abs(nsigmaDCAz) < 7 && std::abs(fDCAxy) < 0.2");
  auto dfSecondary = dfBase.Filter("fTPCnCls > 120 && nITScls >= 6 && std::abs(nsigmaDCAz) > 7 && std::abs(fDCAxy) < 0.2");

  if (skim) {
    dfBase.Filter("fTPCnCls >= 110 && nITScls >= 5 && std::abs(nsigmaDCAz) < 8 && std::abs(fDCAxy) < 0.2 && std::abs(nsigmaHe3) < 5 && pt > 0.8 && pt < 7.0").Snapshot("nucleiTree", "data/skimmed.root");
  }

  std::vector<ROOT::RDF::RResultPtr<TH2D>> hDCAxyAHe3, hDCAxyMHe3, hDCAxySecondaryMHe3, hDCAxySecondaryAHe3, hDCAzAHe3, hDCAzMHe3, hTPCAHe3, hTPCMHe3, hTOFAHe3, hTOFMHe3;

  hDCAxyAHe3.push_back(dfPrimary.Filter("!matter && nsigmaHe3 > -0.5 && nsigmaHe3 < 3 && (!hasTOF || std::abs(deltaMass) < 0.4)").Histo2D({"hDCAxyAHe3", ";#it{p}_{T}^{rec} (GeV/#it{c});DCA_{xy} (cm);Counts", kNPtBins, kPtBins, 100, -0.2, 0.2}, "pt", "fDCAxy"));
  hDCAxyMHe3.push_back(dfPrimary.Filter("matter && nsigmaHe3 > -0.5 && nsigmaHe3 < 3 && (!hasTOF || std::abs(deltaMass) < 0.4)").Histo2D({"hDCAxyMHe3", ";#it{p}_{T}^{rec} (GeV/#it{c});DCA_{xy} (cm);Counts", kNPtBins, kPtBins, 100, -0.2, 0.2}, "pt", "fDCAxy"));
  hDCAzAHe3.push_back(dfPrimary.Filter("!matter && nsigmaHe3 > -0.5 && nsigmaHe3 < 3 && (!hasTOF || std::abs(deltaMass) < 0.4)").Histo2D({"hDCAzAHe3", ";#it{p}_{T}^{rec} (GeV/#it{c});DCA_{z} (cm);Counts", kNPtBins, kPtBins, 100, -0.2, 0.2}, "pt", "fDCAz"));
  hDCAzMHe3.push_back(dfPrimary.Filter("matter && nsigmaHe3 > -0.5 && nsigmaHe3 < 3 && (!hasTOF || std::abs(deltaMass) < 0.4)").Histo2D({"hDCAzMHe3", ";#it{p}_{T}^{rec} (GeV/#it{c});DCA_{z} (cm);Counts", kNPtBins, kPtBins, 100, -0.2, 0.2}, "pt", "fDCAz"));
  hDCAxySecondaryMHe3.push_back(dfSecondary.Filter("matter && nsigmaHe3 > -0.5 && nsigmaHe3 < 3 && (!hasTOF || std::abs(deltaMass) < 0.4)").Histo2D({"hDCAxySecondaryMHe3", ";#it{p}_{T}^{rec} (GeV/#it{c});DCA_{xy} (cm);Counts", kNPtBins, kPtBins, 100, -0.2, 0.2}, "pt", "fDCAxy"));
  hDCAxySecondaryAHe3.push_back(dfSecondary.Filter("!matter && nsigmaHe3 > -0.5 && nsigmaHe3 < 3 && (!hasTOF || std::abs(deltaMass) < 0.4)").Histo2D({"hDCAxySecondaryAHe3", ";#it{p}_{T}^{rec} (GeV/#it{c});DCA_{xy} (cm);Counts", kNPtBins, kPtBins, 100, -0.2, 0.2}, "pt", "fDCAxy"));


  hTPCAHe3.push_back(dfPrimary.Filter("!matter && (!hasTOF || std::abs(deltaMass) < 0.4) && std::abs(fDCAxy) < 0.7").Histo2D({"fATPCcounts", ";#it{p}_{T}^{rec} (GeV/#it{c});^{3}#bar{He} n#sigma_{TPC};Counts", kNPtBins, kPtBins, 100, -5, 5}, "pt", "nsigmaHe3"));
  hTPCMHe3.push_back(dfPrimary.Filter("matter && (!hasTOF || std::abs(deltaMass) < 0.4) && std::abs(fDCAxy) < 0.7").Histo2D({"fMTPCcounts", ";#it{p}_{T}^{rec} (GeV/#it{c});^{3}He n#sigma_{TPC};Counts", kNPtBins, kPtBins, 100, -5, 5}, "pt", "nsigmaHe3"));

  hTOFAHe3.push_back(dfPrimary.Filter("!matter && std::abs(nsigmaHe3) < 3.5").Histo2D({"fATOFsignal", ";#it{p}_{T}^{rec} (GeV/#it{c});m_{TOF}-m_{^{3}#bar{He}};Counts", kNPtBins, kPtBins, 100, -0.9, 1.1}, "pt", "deltaMass"));
  hTOFMHe3.push_back(dfPrimary.Filter("matter && std::abs(nsigmaHe3) < 3.5").Histo2D({"fMTOFsignal", ";#it{p}_{T}^{rec} (GeV/#it{c});m_{TOF}-m_{^{3}He};Counts", kNPtBins, kPtBins, 100, -0.9, 1.1}, "pt", "deltaMass"));

  // const map<string,vector<float> > kCutNames {{"fDCAxy", {0.3, 0.5, 0.8, 1.0}},{"fDCAz", {0.3, 0.5, 0.8, 1.0}},{"fTPCnCls", {110, 120, 130}},{"nITScls", {5, 6, 7}}};
  // int iTrial{0};
  // size_t nTrials{kCutNames["fDCAz"].size() * kCutNames["fTPCnCls"].size() * kCutNames["nITScls"].size() * kCutNames["fDCAxy"].size()};
  // for (size_t iDCAz{0}; iDCAz < kCutNames["fDCAz"].size(); ++iDCAz)
  // {
  //   auto dfDCAz = dfBase.Filter("std::abs(fDCAz) < " + std::to_string(kCutNames["fDCAz"][iDCAz]));
  //   for (size_t iTPCcls{0}; iTPCcls < kCutNames["fTPCnCls"].size(); ++iTPCcls)
  //   {
  //     auto dfTPCcls = dfDCAz.Filter("fTPCnCls > " + std::to_string(kCutNames["fTPCnCls"][iTPCcls]));
  //     for (size_t iITScls{0}; iITScls < kCutNames["nITScls"].size(); ++iITScls)
  //     {
  //       auto dfITScls = dfTPCcls.Filter("nITScls >= " + std::to_string(kCutNames["nITScls"][iITScls]));
  //       for (size_t iDCAxy{0}; iDCAxy < kCutNames["fDCAxy"].size(); ++iDCAxy)
  //       {
  //         auto dfDCAxy = dfITScls.Filter("std::abs(fDCAxy) < " + std::to_string(kCutNames["fDCAxy"][iDCAxy]));
  //         hDCAxyAHe3.push_back(dfITScls.Filter("!matter && nsigmaHe3 > -0.5 && nsigmaHe3 < 3 && (!hasTOF || std::abs(deltaMass) < 0.4)").Histo2D({Form("hDCAxyAHe3%i", iTrial), ";#it{p}_{T}^{rec} (GeV/#it{c});DCA_{xy} (cm);Counts", kNPtBins, kPtBins, 560, -0.7, 0.7}, "pt", "fDCAxy"));
  //         hDCAxyMHe3.push_back(dfITScls.Filter("matter && nsigmaHe3 > -0.5 && nsigmaHe3 < 3 && (!hasTOF || std::abs(deltaMass) < 0.4)").Histo2D({Form("hDCAxyMHe3%i", iTrial), ";#it{p}_{T}^{rec} (GeV/#it{c});DCA_{xy} (cm);Counts", kNPtBins, kPtBins, 560, -0.7, 0.7}, "pt", "fDCAxy"));
  //         hDCAzAHe3.push_back(dfITScls.Filter("!matter && nsigmaHe3 > -0.5 && nsigmaHe3 < 3 && (!hasTOF || std::abs(deltaMass) < 0.4)").Histo2D({Form("hDCAxyAHe3%i", iTrial), ";#it{p}_{T}^{rec} (GeV/#it{c});DCA_{xy} (cm);Counts", kNPtBins, kPtBins, 560, -0.7, 0.7}, "pt", "fDCAz"));
  //         hDCAzMHe3.push_back(dfITScls.Filter("matter && nsigmaHe3 > -0.5 && nsigmaHe3 < 3 && (!hasTOF || std::abs(deltaMass) < 0.4)").Histo2D({Form("hDCAxyMHe3%i", iTrial), ";#it{p}_{T}^{rec} (GeV/#it{c});DCA_{xy} (cm);Counts", kNPtBins, kPtBins, 560, -0.7, 0.7}, "pt", "fDCAz"));

  //         hTPCAHe3.push_back(dfDCAxy.Filter("!matter && (!hasTOF || std::abs(deltaMass) < 0.4)").Histo2D({Form("fATPCcounts%i", iTrial), ";#it{p}_{T}^{rec} (GeV/#it{c});^{3}#bar{He} n#sigma_{TPC};Counts", kNPtBins, kPtBins, 100, -5, 5}, "pt", "nsigmaHe3"));
  //         hTPCMHe3.push_back(dfDCAxy.Filter("matter && (!hasTOF || std::abs(deltaMass) < 0.4)").Histo2D({Form("fMTPCcounts%i", iTrial), ";#it{p}_{T}^{rec} (GeV/#it{c});^{3}He n#sigma_{TPC};Counts", kNPtBins, kPtBins, 100, -5, 5}, "pt", "nsigmaHe3"));

  //         hTOFAHe3.push_back(dfDCAxy.Filter("!matter && std::abs(nsigmaHe3) < 3.5").Histo2D({Form("fATOFsignal%i", iTrial), ";#it{p}_{T}^{rec} (GeV/#it{c});m_{TOF}-m_{^{3}#bar{He}};Counts", kNPtBins, kPtBins, 100, -0.9, 1.1}, "pt", "deltaMass"));
  //         hTOFMHe3.push_back(dfDCAxy.Filter("matter && std::abs(nsigmaHe3) < 3.5").Histo2D({Form("fMTOFsignal%i", iTrial), ";#it{p}_{T}^{rec} (GeV/#it{c});m_{TOF}-m_{^{3}He};Counts", kNPtBins, kPtBins, 100, -0.9, 1.1}, "pt", "deltaMass"));
  //         iTrial++;
  //       }
  //     }
  //   }
  // }

  new TCanvas;
  hTPCAHe3[0]->DrawClone("col");
  new TCanvas;
  hTPCMHe3[0]->DrawClone("col");
  new TCanvas;
  hTOFAHe3[0]->DrawClone("col");
  new TCanvas;
  hTOFMHe3[0]->DrawClone("col");
  new TCanvas;
  hDCAxyAHe3[0]->DrawClone("col");
  new TCanvas;
  hDCAxyMHe3[0]->DrawClone("col");
  new TCanvas;
  hDCAzAHe3[0]->DrawClone("col");
  new TCanvas;
  hDCAzMHe3[0]->DrawClone("col");
  auto cv = new TCanvas();
  cv->SetRightMargin(0.15);
  hDCAxySecondaryMHe3[0]->DrawClone("colz");
  cv = new TCanvas;
  cv->SetRightMargin(0.15);
  hDCAxySecondaryAHe3[0]->GetZaxis()->SetRangeUser(hDCAxySecondaryMHe3[0]->GetMinimum(), hDCAxySecondaryMHe3[0]->GetMaximum());
  hDCAxySecondaryAHe3[0]->DrawClone("colz");

  TFile outputFile(outputFileName.data(), "recreate");
  auto dir = outputFile.mkdir("nuclei");
  dir->cd();
  hTPCAHe3[0]->Write("fATPCcounts");
  hTPCMHe3[0]->Write("fMTPCcounts");
  hTOFAHe3[0]->Write("fATOFsignal");
  hTOFMHe3[0]->Write("fMTOFsignal");
  hDCAxyAHe3[0]->Write();
  hDCAxyMHe3[0]->Write();
  hDCAzAHe3[0]->Write();
  hDCAzMHe3[0]->Write();
  hDCAxySecondaryMHe3[0]->Write();
  hDCAxySecondaryAHe3[0]->Write();

  // for (size_t iT{0}; iT < nTrials; ++iT)
  // {
  //   auto dir = outputFile.mkdir(Form("nuclei%i", iT));
  //   dir->cd();
  //   hTPCAHe3[iT + 1]->Write("fATPCcounts");
  //   hTPCMHe3[iT + 1]->Write("fMTPCcounts");
  //   hTOFAHe3[iT + 1]->Write("fATOFsignal");
  //   hTOFMHe3[iT + 1]->Write("fMTOFsignal");
  //   hDCAxyAHe3[iT + 1]->Write("hDCAxyAHe3");
  //   hDCAxyMHe3[iT + 1]->Write("hDCAxyMHe3");
  //   hDCAzAHe3[iT + 1]->Write("hDCAzAHe3");
  //   hDCAzMHe3[iT + 1]->Write("hDCAzMHe3");
  // }
}
