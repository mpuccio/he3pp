#include <TH1.h>
#include <TH2D.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include "src/Common.h"
#include "EventFiltering/ZorroSummary.h"

void Systematics() {
  TFile fData(kSignalOutput.data());
  TFile fMC(kMCfilename.data());

  TFile anResults(kDataAnalysisResults.data());
  TH1* hNtvx = (TH1*)anResults.Get("bc-selection-task/hCounterTVX");
  TH1* hNcol = (TH1*)anResults.Get("event-selection-task/hColCounterAcc");
  TH1* hNvtx = (TH1*)anResults.Get("nuclei-spectra/spectra/hRecVtxZData");
  ZorroSummary* zorroSummary = (ZorroSummary*)anResults.Get("nuclei-spectra/zorroSummary");
  double nVerticesNoCut{hNcol->GetEntries()};
  double nVertices{hNvtx->GetEntries()};
  double nTVX{zorroSummary ? zorroSummary->getNormalisationFactor(0) : hNtvx->GetEntries()};
  double effTVX{0.756};
  double vertexingEff{0.921};
  // double norm{nTVX * vertexingEff / effTVX};
  // double norm{nVertices};
  double norm{nTVX / effTVX};

  TFile corr("checkRapidity.root");
  TH1* hRapRatio = (TH1*)corr.Get("rapRatio");

  TH2D *systTPC[2];
  TH2D *systTOF[2];
  for (int iC{0}; iC < 2; ++iC) {
    systTPC[iC] = new TH2D(Form("systTPC%s", kNames[iC].data()), ";#it{p}_{T} (GeV/#it{c});Relative systematics TPC", kNPtBins, kPtBins, 50, -0.5, 0.5);
    systTOF[iC] = new TH2D(Form("systTOF%s", kNames[iC].data()), ";#it{p}_{T} (GeV/#it{c});Relative systematics TOF", kNPtBins, kPtBins, 50, -0.5, 0.5);
  }
  TH1* defaultEffTPC[2]{nullptr};
  TH1* defaultEffTOF[2]{nullptr};
  TH1* defaultTPC[2]{nullptr};
  TH1* defaultTOF[2]{nullptr};
  TH1* defaultTPCuncorr[2]{nullptr};
  TH1* defaultTOFuncorr[2]{nullptr};
  for (auto list_key : *fData.GetListOfKeys())
  {
    /// Preliminary operation to read the list and create an output dir
    if (string(list_key->GetName()).find(kFilterListNames.data()) == string::npos)
      continue;

    TDirectoryFile *listData = (TDirectoryFile *)fData.Get(list_key->GetName());
    TDirectoryFile *listMC = (TDirectoryFile *)fMC.Get(list_key->GetName());

    listData->ls();
    std::cout << "Reading " << list_key->GetName() << "\t" << listData << "\t" << listMC << std::endl;

    std::string TOFnames[2]{"hRawCounts", "hRawCountsBinCounting"};
    std::string TOFpresel[3]{"", "_loose", "_tight"};
    TH1* hDataTOF[2][2][3];//{{nullptr, nullptr}, {nullptr, nullptr}};
    TH1* hDataTPC[2][3]{{nullptr, nullptr}, {nullptr, nullptr}};
    TH1* hEffTPC[2]{nullptr, nullptr};
    TH1* hEffTOF[2]{nullptr, nullptr};
    for (int iS{0}; iS < 2; ++iS) {
      auto dataDir = (TDirectoryFile*)listData->Get(Form("%s", kNames[iS].data()));
      hEffTPC[iS] = (TH1*)listMC->Get(Form("effTPC%c", kLetter[iS]));
      hEffTOF[iS] = (TH1*)listMC->Get(Form("effTOF%c", kLetter[iS]));
      if (!hEffTOF[iS]) {
        std::cout << "Missing " << Form("effTOF%c", kLetter[iS]) << std::endl;
      }
      if (!hEffTPC[iS]) {
        std::cout << "Missing " << Form("effTPC%c", kLetter[iS]) << std::endl;
      }
      if (!defaultEffTPC[iS])
        defaultEffTPC[iS] = (TH1*)hEffTPC[iS]->Clone(Form("defaultEffTPC%s", kNames[iS].data()));
      if (!defaultEffTOF[iS])
        defaultEffTOF[iS] = (TH1*)hEffTOF[iS]->Clone(Form("defaultEffTOF%s", kNames[iS].data()));

      for (int iTOF{0}; iTOF < 2; ++iTOF) {
        std::string name{TOFnames[iTOF] + kLetter[iS]};
        hDataTOF[iS][iTOF][0]=(TH1*)dataDir->Get(Form("GausExp/%s0%s", name.data(), TOFpresel[0].data()));
        if (!hDataTOF[iS][iTOF][0])
          std::cout << "Missing " << Form("%s/GausExp/%s%s0", kNames[iS].data(), TOFnames[iTOF].data(), kLetter[iS]) << std::flush << std::endl;
        if (!defaultTOFuncorr[iS])
          defaultTOFuncorr[iS] = (TH1*)hDataTOF[iS][iTOF][0]->Clone(Form("defaultTOFuncorr%s", kNames[iS].data()));
        hDataTOF[iS][iTOF][0]->Divide(hEffTOF[iS]);
        if (!defaultTOF[iS])
          defaultTOF[iS] = (TH1*)hDataTOF[iS][iTOF][0]->Clone(Form("defaultTOF%s", kNames[iS].data()));
      }
      for (int iTPC{0}; iTPC < 3; ++iTPC) {
        std::string name{kTPCfunctName[iTPC]};

        hDataTPC[iS][iTPC]=(TH1*)dataDir->Get(Form("TPConly/hTPConly%c0_%s", kLetter[iS], name.data()));
        if (!hDataTPC[iS][iTPC])
          std::cout << "Missing " << Form("TPConly/%s0", name.data()) << std::endl;
        if (!defaultTPCuncorr[iS] && iTPC == 1)
          defaultTPCuncorr[iS] = (TH1*)hDataTPC[iS][iTPC]->Clone(Form("defaultTPCuncorr%s", kNames[iS].data()));
        hDataTPC[iS][iTPC]->Divide(hEffTPC[iS]);
        if (!defaultTPC[iS] && iTPC == 1)
          defaultTPC[iS] = (TH1*)hDataTPC[iS][iTPC]->Clone(Form("defaultTPC%s", kNames[iS].data()));
      }
    }

    for (int iS{0}; iS < 2; ++iS) {
      for (int iB{1}; iB <= kNPtBins; ++iB) {
        float pt = hDataTPC[iS][0]->GetBinCenter(iB);
        float defaultValueTPC = defaultTPC[iS]->GetBinContent(iB);
        float defaultValueTOF = defaultTOF[iS]->GetBinContent(iB);
        for (int iTPC{0}; iTPC < 3; ++iTPC) {
          float value = hDataTPC[iS][iTPC]->GetBinContent(iB);
          systTPC[iS]->Fill(pt, (value - defaultValueTPC) / defaultValueTPC);
        }
        for (int iTOF{0}; iTOF < 2; ++iTOF) {
          float value = hDataTOF[iS][iTOF][0]->GetBinContent(iB);
          systTOF[iS]->Fill(pt, (value - defaultValueTOF) / defaultValueTOF);
        }
      }
    }
  }

  TH1D* hSystTPC[2];
  TH1D* hSystTOF[2];
  for (int iS{0}; iS < 2; ++iS) {
    hSystTPC[iS] = new TH1D(Form("hSystTPC%c",kLetter[iS]), ";#it{p}_{T} (GeV/#it{c});Relative systematics TPC", kNPtBins, kPtBins);
    hSystTOF[iS] = new TH1D(Form("hSystTOF%c",kLetter[iS]), ";#it{p}_{T} (GeV/#it{c});Relative systematics TOF", kNPtBins, kPtBins);
    for (int iB{1}; iB <= kNPtBins; ++iB) {
      hSystTPC[iS]->SetBinContent(iB, systTPC[iS]->ProjectionY("", iB, iB)->GetRMS());
      hSystTOF[iS]->SetBinContent(iB, systTOF[iS]->ProjectionY("", iB, iB)->GetRMS());
    }
  }

  TFile syst(kSystematicsOutput.data(), "recreate");

  TH1* tofMatchingM = (TH1*)defaultTOFuncorr[0]->Clone(Form("TOFmatching%s", kNames[0].data()));
  TH1* tofMatchingA = (TH1*)defaultTOFuncorr[1]->Clone(Form("TOFmatching%s", kNames[1].data()));
  tofMatchingM->Divide(defaultTPCuncorr[0]);
  tofMatchingA->Divide(defaultTPCuncorr[1]);
  tofMatchingM->Write();
  tofMatchingA->Write();

  systTPC[0]->Write();
  systTPC[1]->Write();
  systTOF[0]->Write();
  systTOF[1]->Write();

  systTPC[0]->ProjectionY("", 6, 6)->Write();
  systTPC[1]->ProjectionY("", 6, 6)->Write();
  systTOF[0]->ProjectionY("", 6, 6)->Write();
  systTOF[1]->ProjectionY("", 6, 6)->Write();

  hSystTPC[0]->Write();
  hSystTPC[1]->Write();
  hSystTOF[0]->Write();
  hSystTOF[1]->Write();
  new TCanvas;
  systTPC[0]->DrawClone("col");
    new TCanvas;
  systTPC[1]->DrawClone("col");
    new TCanvas;
  systTOF[0]->DrawClone("col");
    new TCanvas;
  systTOF[1]->DrawClone("col");

  constexpr int pubBinning{6};
  constexpr double pubBins[pubBinning + 1]{1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0};
  double y[]{1.2241e-07, 8.4801e-08, 5.0085e-08, 3.2333e-08, 1.7168e-08, 4.8137e-09};
  double ey[]{1.769e-08, 7.5127e-09, 6.0035e-09, 4.8788e-09, 2.5057e-09, 1.3356e-09};
  double sy[]{ 1.3346e-08, 1.0763e-08, 3.2452e-09, 2.1084e-09, 1.1316e-09, 3.1345e-10};

  TH1D* hPub = new TH1D("hPub", ";#it{p}_{T} (GeV/#it{c});#frac{1}{N_{ev}} #frac{d^{2}N}{dyd#it{p}_{T}}", pubBinning, pubBins);
  TH1D* hPubSyst = new TH1D("hPubSyst", ";#it{p}_{T} (GeV/#it{c});#frac{1}{N_{ev}} #frac{d^{2}N}{dyd#it{p}_{T}}", pubBinning, pubBins);
  hPub->SetLineColor(kBlack);
  hPubSyst->SetLineColor(kBlack);
  hPub->SetMarkerColor(kBlack);
  hPubSyst->SetMarkerColor(kBlack);
  hPub->SetMarkerStyle(5);
  hPub->SetMarkerSize(1);
  hPubSyst->SetMarkerStyle(5);
  hPubSyst->SetMarkerSize(1);
  hPubSyst->SetFillStyle(0);
  for (int i{1}; i <= 6; ++i) {
    hPub->SetBinContent(i, y[i-1]);
    hPubSyst->SetBinContent(i, y[i-1]);
    hPub->SetBinError(i, ey[i-1]);
    hPubSyst->SetBinError(i, sy[i-1]);
  }
  hPub->Multiply(hRapRatio);
  hPubSyst->Multiply(hRapRatio);
  hPub->Write("pubStat");
  hPubSyst->Write("pubSyst");

  for (int iS{0}; iS < 2; ++iS) {
    TH1D* fStatTPC = (TH1D*)defaultTPCuncorr[iS]->Clone(Form("fStatTPC%c", kLetter[iS]));
    TH1D* fSystTPC = (TH1D*)defaultTPCuncorr[iS]->Clone(Form("fSystTPC%c", kLetter[iS]));
    TH1D* fStatTOF = (TH1D*)defaultTOFuncorr[iS]->Clone(Form("fStatTOF%c", kLetter[iS]));
    TH1D* fSystTOF = (TH1D*)defaultTOFuncorr[iS]->Clone(Form("fSystTOF%c", kLetter[iS]));
    fStatTPC->SetTitle(";#it{p}_{T} (GeV/#it{c});#frac{1}{N_{ev}} #frac{d^{2}N}{dyd#it{p}_{T}}");
    fSystTPC->SetTitle(";#it{p}_{T} (GeV/#it{c});#frac{1}{N_{ev}} #frac{d^{2}N}{dyd#it{p}_{T}}");
    fStatTPC->SetMarkerStyle(20);
    fSystTPC->SetMarkerStyle(20);
    fStatTOF->SetMarkerStyle(21);
    fSystTOF->SetMarkerStyle(21);
    fStatTPC->SetMarkerColor(kRed);
    fSystTPC->SetMarkerColor(kRed);
    fStatTOF->SetMarkerColor(kBlue);
    fSystTOF->SetMarkerColor(kBlue);
    fStatTPC->SetLineColor(kRed);
    fSystTPC->SetLineColor(kRed);
    fStatTOF->SetLineColor(kBlue);
    fSystTOF->SetLineColor(kBlue);
    fSystTPC->SetFillStyle(0);
    fSystTOF->SetFillStyle(0);
    for (int iBin{1}; iBin <= kNPtBins; ++iBin) {
      double binW = defaultTPCuncorr[iS]->GetBinWidth(iBin);
      double yieldTPC = defaultTPCuncorr[iS]->GetBinContent(iBin);
      double yieldTOF = defaultTOFuncorr[iS]->GetBinContent(iBin);
      double statErrorTPC = defaultTPCuncorr[iS]->GetBinError(iBin);
      double statErrorTOF = defaultTOFuncorr[iS]->GetBinError(iBin);
      double effTPC = defaultEffTPC[iS]->GetBinContent(iBin);
      double effTOF = defaultEffTOF[iS]->GetBinContent(iBin);
      std::cout << "Bin " << iBin << "\t" << yieldTPC << "\t" << statErrorTPC << "\t" << effTPC << "\t" << yieldTOF << "\t" << statErrorTOF << "\t" << effTOF << std::endl;
      if (effTPC >= 1.e-2) {
        fStatTPC->SetBinContent(iBin, yieldTPC / effTPC);
        fStatTPC->SetBinError(iBin, defaultTPCuncorr[iS]->GetBinError(iBin) / effTPC);
        fSystTPC->SetBinContent(iBin, yieldTPC / effTPC);
          fSystTPC->SetBinError(iBin, hSystTPC[iS]->GetBinContent(iBin) * yieldTPC / effTPC);
      } else {
        fStatTPC->SetBinContent(iBin, 0);
        fStatTPC->SetBinError(iBin, 0);
        fSystTPC->SetBinContent(iBin, 0);
        fSystTPC->SetBinError(iBin, 0);
      }
      if (effTOF >= 1.e-2) {
        fStatTOF->SetBinContent(iBin, yieldTOF / effTOF);
        fStatTOF->SetBinError(iBin, defaultTOFuncorr[iS]->GetBinError(iBin) / effTOF);
        fSystTOF->SetBinContent(iBin, yieldTOF / effTOF);
        fSystTOF->SetBinError(iBin, hSystTOF[iS]->GetBinContent(iBin) * yieldTOF / effTOF);
      } else {
        fStatTOF->SetBinContent(iBin, 0);
        fStatTOF->SetBinError(iBin, 0);
        fSystTOF->SetBinContent(iBin, 0);
        fSystTOF->SetBinError(iBin, 0);
      }
    }

    fStatTPC->Scale(1./norm, "width");
    fSystTPC->Scale(1./norm, "width");
    fStatTOF->Scale(1./norm, "width");
    fSystTOF->Scale(1./norm, "width");

    fStatTPC->Write();
    fSystTPC->Write();
    fStatTOF->Write();
    fSystTOF->Write();

    TCanvas *cv = new TCanvas(Form("canvas%c", kLetter[iS]));
    cv->DrawFrame(0, 0.5e-10, 8, 9.e-7, ";#it{p}_{T} (GeV/#it{c});#frac{1}{N_{ev}} #frac{d^{2}N}{dyd#it{p}_{T}}");
    cv->SetLogy();
    fStatTPC->DrawClone("x0esame");
    fSystTPC->DrawClone("e2same");
    fStatTOF->DrawClone("x0esame");
    fSystTOF->DrawClone("e2same");
    hPub->DrawClone("x0esame");
    hPubSyst->DrawClone("e2same");
    TLegend *leg = new TLegend(0.7, 0.75, 0.85, 0.9);
    leg->SetFillStyle(0);
    leg->AddEntry(fSystTPC, Form("TPC %s", kLabels[iS].data()), "lep");
    leg->AddEntry(fSystTOF, Form("TOF %s", kLabels[iS].data()), "lep");
    leg->AddEntry(hPubSyst, Form("Publication %s", kLabels[iS].data()), "lep");
    leg->Draw();
    cv->Write();

    fStatTPC->Divide(hPub);
    fStatTOF->Divide(hPub);
    fStatTPC->Write(Form("ratioPubTPC%c", kLetter[iS]));
    fStatTOF->Write(Form("ratioPubTOF%c", kLetter[iS]));

    if (hRapRatio) {
      std::cout << "Correcting for rapidity" << std::endl;
      // for (int iB{1}; iB <= fStatTPC->GetNbinsX(); ++iB) {
      //   int iR = hRapRatio->FindBin(fStatTPC->GetBinCenter(iB));
      //   fStatTPC->SetBinContent(iB, fStatTPC->GetBinContent(iB) / hRapRatio->GetBinContent(iR));
      //   fStatTOF->SetBinContent(iB, fStatTOF->GetBinContent(iB) / hRapRatio->GetBinContent(iR));
      // }
      fStatTPC->Write(Form("ratioRapTPC%c", kLetter[iS]));
      fStatTOF->Write(Form("ratioRapTOF%c", kLetter[iS]));
    }

  }
}
