#include "ROOT/RDataFrame.hxx"
#include "TF1.h"
#include "TBox.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TStyle.h"
#include <cmath>
#include <iostream>
#include "src/Common.h"
#include "src/Utils.h"

void getTriggerEfficiency()
{
  gStyle->SetOptStat(0);
  ROOT::EnableImplicitMT();
  TChain chainSampled("O2nucleitable"), chainSkimmed("O2nucleitable");

  utils::createChain(chainSampled, kBaseInputDir + "data/" + kPeriod + "/" + kRecoPass + "/AO2D_sampled.root", "O2nucleitable");
  utils::createChain(chainSkimmed, kBaseInputDir + "data/" + kPeriod + "/" + kRecoPass + "/AO2D_skimmed.root", "O2nucleitable");
  ROOT::RDataFrame dfSampled(chainSampled), dfSkimmed(chainSkimmed);

  auto dfSampledFiltered = defineColumnsForData(dfSampled).Filter(kBaseRecSelections.data()).Filter(kDefaultRecSelections.data());
  auto dfSkimmedFiltered = defineColumnsForData(dfSkimmed).Filter(kBaseRecSelections.data()).Filter(kDefaultRecSelections.data());

  std::vector<ROOT::RDF::RResultPtr<TH2D>> hTPCAHe3, hTPCMHe3, hTOFAHe3, hTOFMHe3, hITSAHe3;std::vector<ROOT::RDF::RResultPtr<TH1D>> hDiffPtDist;

  // Sampled data
  hTPCAHe3.push_back(dfSampledFiltered.Filter("!matter && hasGoodTOFmassHe3").Histo2D({"fATPCcounts_sampled", ";#it{p}_{T}^{rec} (GeV/#it{c});^{3}#bar{He} n#sigma_{TPC};Counts", kNPtBins, kPtBins, 100, -5, 5}, "pt", "nsigmaHe3"));
  hTPCMHe3.push_back(dfSampledFiltered.Filter("matter && hasGoodTOFmassHe3").Histo2D({"fMTPCcounts_sampled", ";#it{p}_{T}^{rec} (GeV/#it{c});^{3}He n#sigma_{TPC};Counts", kNPtBins, kPtBins, 100, -5, 5}, "pt", "nsigmaHe3"));
  hTOFAHe3.push_back(dfSampledFiltered.Filter("!matter && std::abs(nsigmaHe3) < 3.5").Histo2D({"fATOFsignal_sampled", ";#it{p}_{T}^{rec} (GeV/#it{c});m_{TOF}-m_{^{3}#bar{He}};Counts", kNPtBins, kPtBins, 100, -0.9, 1.1}, "pt", "deltaMassHe3"));
  hTOFMHe3.push_back(dfSampledFiltered.Filter("matter && std::abs(nsigmaHe3) < 3.5").Histo2D({"fMTOFsignal_sampled", ";#it{p}_{T}^{rec} (GeV/#it{c});m_{TOF}-m_{^{3}He};Counts", kNPtBins, kPtBins, 100, -0.9, 1.1}, "pt", "deltaMassHe3"));
  hITSAHe3.push_back(dfSampledFiltered.Filter("!matter && std::abs(nsigmaHe3) < 3.5").Histo2D({"fAITScounts_sampled", ";#it{p}_{T}^{rec} (GeV/#it{c});^{3}#bar{He} n#sigma_{ITS};Counts", kNPtBins, kPtBins, 100, -5, 5}, "pt", "nsigmaITS"));

  // Skimmed data
  hTPCAHe3.push_back(dfSkimmedFiltered.Filter("!matter && hasGoodTOFmassHe3").Histo2D({"fATPCcounts_skimmed", ";#it{p}_{T}^{rec} (GeV/#it{c});^{3}#bar{He} n#sigma_{TPC};Counts", kNPtBins, kPtBins, 100, -5, 5}, "pt", "nsigmaHe3"));
  hTPCMHe3.push_back(dfSkimmedFiltered.Filter("matter && hasGoodTOFmassHe3").Histo2D({"fMTPCcounts_skimmed", ";#it{p}_{T}^{rec} (GeV/#it{c});^{3}He n#sigma_{TPC};Counts", kNPtBins, kPtBins, 100, -5, 5}, "pt", "nsigmaHe3"));
  hTOFAHe3.push_back(dfSkimmedFiltered.Filter("!matter && std::abs(nsigmaHe3) < 3.5").Histo2D({"fATOFsignal_skimmed", ";#it{p}_{T}^{rec} (GeV/#it{c});m_{TOF}-m_{^{3}#bar{He}};Counts", kNPtBins, kPtBins, 100, -0.9, 1.1}, "pt", "deltaMassHe3"));
  hTOFMHe3.push_back(dfSkimmedFiltered.Filter("matter && std::abs(nsigmaHe3) < 3.5").Histo2D({"fMTOFsignal_skimmed", ";#it{p}_{T}^{rec} (GeV/#it{c});m_{TOF}-m_{^{3}He};Counts", kNPtBins, kPtBins, 100, -0.9, 1.1}, "pt", "deltaMassHe3"));
  hITSAHe3.push_back(dfSkimmedFiltered.Filter("!matter && std::abs(nsigmaHe3) < 3.5").Histo2D({"fAITScounts_skimmed", ";#it{p}_{T}^{rec} (GeV/#it{c});^{3}#bar{He} n#sigma_{ITS};Counts", kNPtBins, kPtBins, 100, -5, 5}, "pt", "nsigmaITS"));


  hDiffPtDist.push_back(dfSampledFiltered.Filter("!matter && nsigmaHe3 > -2. && nsigmaHe3 < 3").Histo1D({"hPtDist_sampled", ";#it{p}_{T}^{gen} (GeV/#it{c});Counts", kNPtBins, kPtBins}, "pt"));
  hDiffPtDist.push_back(dfSkimmedFiltered.Filter("!matter && nsigmaHe3 > -2. && nsigmaHe3 < 3").Histo1D({"hPtDist_skimmed", ";#it{p}_{T}^{gen} (GeV/#it{c});Counts", kNPtBins, kPtBins}, "pt"));

  TFile anResultSkimmed((kBaseInputDir + "data/" + kPeriod + "/" + kRecoPass + "/AnalysisResults_skimmed.root").data(), "read");
  ZorroSummary* zorroSummarySkimmed = (ZorroSummary*)anResultSkimmed.Get("nuclei-spectra/zorroSummary");
  hDiffPtDist[1]->Scale(1.0 / zorroSummarySkimmed->getNormalisationFactor(0));

  TFile anResultSampled((kBaseInputDir + "data/" + kPeriod + "/" + kRecoPass + "/AnalysisResults_sampled.root").data(), "read");
  TH1* hCounterTVX = (TH1*)anResultSampled.Get("eventselection-run3/luminosity/hCounterTVX");
  hDiffPtDist[0]->Scale(1.0 / hCounterTVX->GetEntries());
  ZorroSummary* zorroSummarySampled = (ZorroSummary*)anResultSampled.Get("nuclei-spectra/zorroSummary");
  std::cout << "Normalisation factor (sampled): " << zorroSummarySampled->getNormalisationFactor(0) << "\t" << hCounterTVX->GetEntries() << ";\t Zorro / TVX: " << zorroSummarySampled->getNormalisationFactor(0) / hCounterTVX->GetEntries() << std::endl;

  TFile outputFile("triggerEfficiency.root", "recreate");
  hTPCAHe3[0]->Write("fATPCcounts_sampled");
  hTPCMHe3[0]->Write("fMTPCcounts_sampled");
  hTOFAHe3[0]->Write("fATOFsignal_sampled");
  hTOFMHe3[0]->Write("fMTOFsignal_sampled");
  hTPCAHe3[1]->Write("fATPCcounts_skimmed");
  hTPCMHe3[1]->Write("fMTPCcounts_skimmed");
  hTOFAHe3[1]->Write("fATOFsignal_skimmed");
  hTOFMHe3[1]->Write("fMTOFsignal_skimmed");
  hITSAHe3[0]->Write("fAITScounts_sampled");
  hITSAHe3[1]->Write("fAITScounts_skimmed");
  hDiffPtDist[0]->Write("hPtDist_sampled");
  hDiffPtDist[1]->Write("hPtDist_skimmed");

  // Create ratio histogram with proper error propagation for correlated histograms
  TH1D* hTriggerEfficiency = (TH1D*)hDiffPtDist[1]->Clone("hTriggerEfficiency");
  hTriggerEfficiency->SetTitle("Trigger Efficiency;#it{p}_{T} (GeV/#it{c});Efficiency");
  hTriggerEfficiency->Reset();

  for (int i = 1; i <= hTriggerEfficiency->GetNbinsX(); ++i) {
    double skimmed = hDiffPtDist[1]->GetBinContent(i);
    double sampled = hDiffPtDist[0]->GetBinContent(i);
    double errorSkimmed = hDiffPtDist[1]->GetBinError(i);
    double errorSampled = hDiffPtDist[0]->GetBinError(i);

    if (skimmed > 0) {
      double ratio = sampled / skimmed;
      // For correlated histograms (sampled contains skimmed):
      // sigma(a/b)^2 = (1/b^2)[sigma_a^2 + (a/b)^2*sigma_b^2 - 2*(a/b)*cov(a,b)]
      // Assuming full correlation for the overlap: cov(a,b) ≈ sigma_b^2
      double error = std::sqrt(std::max(0.0,
                               (errorSampled * errorSampled) / (skimmed * skimmed) +
                               (ratio * ratio * errorSkimmed * errorSkimmed) / (skimmed * skimmed) -
                               2.0 * ratio * errorSkimmed * errorSkimmed / (skimmed * skimmed)));
      hTriggerEfficiency->SetBinContent(i, ratio);
      hTriggerEfficiency->SetBinError(i, error);
    }
  }

  hTriggerEfficiency->Write("hTriggerEfficiency");

}
