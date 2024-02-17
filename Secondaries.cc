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
  TFile input_file(kDataFilename.data());
  TFile output_file(kSignalOutput.data(), "recreate");

  /// Setting up the fitting environment for TOF analysis
  RooRealVar m("dm2", "m^{2} - m^p_{^{3}He}", -1.2, 1.5, "GeV/#it{c}^{2}");
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
  }
}