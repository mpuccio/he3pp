#include "ROOT/RDataFrame.hxx"
#include "TF1.h"
#include "TBox.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TStyle.h"
#include <cmath>
#include <iostream>
#include "src/Common.h"

void analyseData(std::string inputFileName = kDataTreeFilename, std::string outputFileName = kDataFilename, bool skim = false)
{
  gStyle->SetOptStat(0);
  ROOT::EnableImplicitMT();
  ROOT::RDataFrame d("O2nucleitable", inputFileName);
  auto dfBase = defineColumnsForData(d).Filter("fTPCnCls >= 110 && nITScls >= 5 && abs(fEta) < 0.9 && std::abs(fDCAxy) < 0.7 && pt > 0.8 && pt < 9.0");
  auto dfPrimary = dfBase.Filter("fTPCnCls > 120 && nITScls >= 6 && std::abs(nsigmaDCAz) < 7 && std::abs(fDCAxy) < 0.2");
  auto dfSecondary = dfBase.Filter("fTPCnCls > 120 && nITScls >= 6 && std::abs(nsigmaDCAz) > 7 && std::abs(fDCAxy) < 0.2");

  if (skim)
  {
    dfBase.Filter("std::abs(nsigmaDCAz) < 8 && std::abs(fDCAxy) < 0.2 && std::abs(nsigmaHe3) < 5").Snapshot("nucleiTree", "data/skimmed.root");
  }

  std::vector<ROOT::RDF::RResultPtr<TH2D>> hDCAxyAHe3, hDCAxyMHe3, hDCAxySecondaryMHe3, hDCAxySecondaryAHe3, hDCAzAHe3, hDCAzMHe3, hTPCAHe3, hTPCMHe3, hTOFAHe3, hTOFMHe3;

  hDCAxyAHe3.push_back(dfPrimary.Filter("!matter && nsigmaHe3 > -0.5 && nsigmaHe3 < 3 && hasGoodTOFmassHe3").Histo2D({"hDCAxyAHe3", ";#it{p}_{T}^{rec} (GeV/#it{c});DCA_{xy} (cm);Counts", kNPtBins, kPtBins, 100, -0.2, 0.2}, "pt", "fDCAxy"));
  hDCAxyMHe3.push_back(dfPrimary.Filter("matter && nsigmaHe3 > -0.5 && nsigmaHe3 < 3 && hasGoodTOFmassHe3").Histo2D({"hDCAxyMHe3", ";#it{p}_{T}^{rec} (GeV/#it{c});DCA_{xy} (cm);Counts", kNPtBins, kPtBins, 100, -0.2, 0.2}, "pt", "fDCAxy"));
  hDCAzAHe3.push_back(dfPrimary.Filter("!matter && nsigmaHe3 > -0.5 && nsigmaHe3 < 3 && hasGoodTOFmassHe3").Histo2D({"hDCAzAHe3", ";#it{p}_{T}^{rec} (GeV/#it{c});DCA_{z} (cm);Counts", kNPtBins, kPtBins, 100, -0.2, 0.2}, "pt", "fDCAz"));
  hDCAzMHe3.push_back(dfPrimary.Filter("matter && nsigmaHe3 > -0.5 && nsigmaHe3 < 3 && hasGoodTOFmassHe3").Histo2D({"hDCAzMHe3", ";#it{p}_{T}^{rec} (GeV/#it{c});DCA_{z} (cm);Counts", kNPtBins, kPtBins, 100, -0.2, 0.2}, "pt", "fDCAz"));
  hDCAxySecondaryMHe3.push_back(dfSecondary.Filter("matter && nsigmaHe3 > -0.5 && nsigmaHe3 < 3 && hasGoodTOFmassHe3").Histo2D({"hDCAxySecondaryMHe3", ";#it{p}_{T}^{rec} (GeV/#it{c});DCA_{xy} (cm);Counts", kNPtBins, kPtBins, 100, -0.2, 0.2}, "pt", "fDCAxy"));
  hDCAxySecondaryAHe3.push_back(dfSecondary.Filter("!matter && nsigmaHe3 > -0.5 && nsigmaHe3 < 3 && hasGoodTOFmassHe3").Histo2D({"hDCAxySecondaryAHe3", ";#it{p}_{T}^{rec} (GeV/#it{c});DCA_{xy} (cm);Counts", kNPtBins, kPtBins, 100, -0.2, 0.2}, "pt", "fDCAxy"));

  hTPCAHe3.push_back(dfPrimary.Filter("!matter && hasGoodTOFmassHe3").Histo2D({"fATPCcounts", ";#it{p}_{T}^{rec} (GeV/#it{c});^{3}#bar{He} n#sigma_{TPC};Counts", kNPtBins, kPtBins, 100, -5, 5}, "pt", "nsigmaHe3"));
  hTPCMHe3.push_back(dfPrimary.Filter("matter && hasGoodTOFmassHe3").Histo2D({"fMTPCcounts", ";#it{p}_{T}^{rec} (GeV/#it{c});^{3}He n#sigma_{TPC};Counts", kNPtBins, kPtBins, 100, -5, 5}, "pt", "nsigmaHe3"));

  hTOFAHe3.push_back(dfPrimary.Filter("!matter && std::abs(nsigmaHe3) < 3.5").Histo2D({"fATOFsignal", ";#it{p}_{T}^{rec} (GeV/#it{c});m_{TOF}-m_{^{3}#bar{He}};Counts", kNPtBins, kPtBins, 100, -0.9, 1.1}, "pt", "deltaMassHe3"));
  hTOFMHe3.push_back(dfPrimary.Filter("matter && std::abs(nsigmaHe3) < 3.5").Histo2D({"fMTOFsignal", ";#it{p}_{T}^{rec} (GeV/#it{c});m_{TOF}-m_{^{3}He};Counts", kNPtBins, kPtBins, 100, -0.9, 1.1}, "pt", "deltaMassHe3"));

  int iTrial{0};
  for (size_t iDCAz{0}; iDCAz < kCutNames["nsigmaDCAz"].size(); ++iDCAz)
  {
    auto dfDCAz = dfBase.Filter("std::abs(nsigmaDCAz) < " + std::to_string(kCutNames["nsigmaDCAz"][iDCAz]));
    for (size_t iTPCcls{0}; iTPCcls < kCutNames["fTPCnCls"].size(); ++iTPCcls)
    {
      auto dfTPCcls = dfDCAz.Filter("fTPCnCls > " + std::to_string(kCutNames["fTPCnCls"][iTPCcls]));
      for (size_t iITScls{0}; iITScls < kCutNames["nITScls"].size(); ++iITScls)
      {
        auto dfITScls = dfTPCcls.Filter("nITScls >= " + std::to_string(kCutNames["nITScls"][iITScls]));
        hDCAxyAHe3.push_back(dfITScls.Filter("!matter && nsigmaHe3 > -0.5 && nsigmaHe3 < 3 && hasGoodTOFmassHe3").Histo2D({Form("hDCAxyAHe3%i", iTrial), ";#it{p}_{T}^{rec} (GeV/#it{c});DCA_{xy} (cm);Counts", kNPtBins, kPtBins, 560, -0.7, 0.7}, "pt", "fDCAxy"));
        hDCAxyMHe3.push_back(dfITScls.Filter("matter && nsigmaHe3 > -0.5 && nsigmaHe3 < 3 && hasGoodTOFmassHe3").Histo2D({Form("hDCAxyMHe3%i", iTrial), ";#it{p}_{T}^{rec} (GeV/#it{c});DCA_{xy} (cm);Counts", kNPtBins, kPtBins, 560, -0.7, 0.7}, "pt", "fDCAxy"));
        hDCAzAHe3.push_back(dfITScls.Filter("!matter && nsigmaHe3 > -0.5 && nsigmaHe3 < 3 && hasGoodTOFmassHe3").Histo2D({Form("hDCzAHe3%i", iTrial), ";#it{p}_{T}^{rec} (GeV/#it{c});DCA_{xy} (cm);Counts", kNPtBins, kPtBins, 560, -0.7, 0.7}, "pt", "fDCAz"));
        hDCAzMHe3.push_back(dfITScls.Filter("matter && nsigmaHe3 > -0.5 && nsigmaHe3 < 3 && hasGoodTOFmassHe3").Histo2D({Form("hDCAzMHe3%i", iTrial), ";#it{p}_{T}^{rec} (GeV/#it{c});DCA_{xy} (cm);Counts", kNPtBins, kPtBins, 560, -0.7, 0.7}, "pt", "fDCAz"));

        hTPCAHe3.push_back(dfITScls.Filter("!matter && std::abs(fDCAxy) < 0.2 && hasGoodTOFmassHe3").Histo2D({Form("fATPCcounts%i", iTrial), ";#it{p}_{T}^{rec} (GeV/#it{c});^{3}#bar{He} n#sigma_{TPC};Counts", kNPtBins, kPtBins, 100, -5, 5}, "pt", "nsigmaHe3"));
        hTPCMHe3.push_back(dfITScls.Filter("matter && std::abs(fDCAxy) < 0.2 && hasGoodTOFmassHe3").Histo2D({Form("fMTPCcounts%i", iTrial), ";#it{p}_{T}^{rec} (GeV/#it{c});^{3}He n#sigma_{TPC};Counts", kNPtBins, kPtBins, 100, -5, 5}, "pt", "nsigmaHe3"));

        hTOFAHe3.push_back(dfITScls.Filter("!matter && std::abs(fDCAxy) < 0.2 && std::abs(nsigmaHe3) < 3.5").Histo2D({Form("fATOFsignal%i", iTrial), ";#it{p}_{T}^{rec} (GeV/#it{c});m_{TOF}-m_{^{3}#bar{He}};Counts", kNPtBins, kPtBins, 100, -0.9, 1.1}, "pt", "deltaMassHe3"));
        hTOFMHe3.push_back(dfITScls.Filter("matter && std::abs(fDCAxy) < 0.2 && std::abs(nsigmaHe3) < 3.5").Histo2D({Form("fMTOFsignal%i", iTrial), ";#it{p}_{T}^{rec} (GeV/#it{c});m_{TOF}-m_{^{3}He};Counts", kNPtBins, kPtBins, 100, -0.9, 1.1}, "pt", "deltaMassHe3"));
        iTrial++;
      }
    }
  }

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

  for (size_t iT{0}; iT < iTrial; ++iT)
  {
    auto dir = outputFile.mkdir(Form("nuclei%i", iT));
    dir->cd();
    hTPCAHe3[iT + 1]->Write("fATPCcounts");
    hTPCMHe3[iT + 1]->Write("fMTPCcounts");
    hTOFAHe3[iT + 1]->Write("fATOFsignal");
    hTOFMHe3[iT + 1]->Write("fMTOFsignal");
    hDCAxyAHe3[iT + 1]->Write("hDCAxyAHe3");
    hDCAxyMHe3[iT + 1]->Write("hDCAxyMHe3");
    hDCAzAHe3[iT + 1]->Write("hDCAzAHe3");
    hDCAzMHe3[iT + 1]->Write("hDCAzMHe3");
  }
}
