#include "FitModules.h"
#include "Common.h"
#include "RooGausExp.h"
#include "RooGausDExp.h"
#include "Utils.h"
using namespace utils;

#include <TLegend.h>

#include <RooArgList.h>
#include <RooArgSet.h>
#include <RooCBShape.h>
#include <RooChebychev.h>
#include <RooCurve.h>
#include <RooDataHist.h>
#include <RooFitResult.h>
#include <RooLognormal.h>


RooPlot* FitModule::FitData(TH1* dat,TString name, TString title, TString range, TString plotrange, bool change_range, float low_x, float high_x) {
  if (plotrange == "") plotrange = range;
  RooDataHist data("data","data",RooArgList(*mX),Import(*dat));
  RooPlot* plot = mX->frame();
  plot->SetTitle(title.Data());
  plot->SetName(name.Data());
  plot->GetYaxis()->SetTitle(Form("Counts / (%.2f Gev/#it{c}^{2})",plot->GetXaxis()->GetBinWidth(1)));
  for (int i = 2; i--;) RooFitResult *res = mTemplate->fitTo(data,Extended(),Verbose(kFALSE),PrintEvalErrors(-1),PrintLevel(-1),Range(range));
  if(change_range) plot->GetXaxis()->SetRangeUser(low_x,high_x);
  data.plotOn(plot,Name("data"),DrawOption("pz"));
  mTemplate->plotOn(plot,Name("model"),Range(plotrange),NormRange(range));
  mTemplate->plotOn(plot,Name("bkg"),Components(*mBackground),LineStyle(kDashed),LineColor(kRed),Range(plotrange),NormRange(range));
  mTemplate->plotOn(plot,Name("sig"),Components(*mSignal),LineStyle(kDashed),LineColor(kGreen+3),Range(plotrange),NormRange(range));
  mChi2 = plot->chiSquare("model","data");
  plot->remove("model",false);
  plot->remove("bkg",false);
  plot->remove("sig",false);
  mTemplate->plotOn(plot,Name("model"),Range(plotrange),NormRange(range));
  mTemplate->plotOn(plot,Name("bkg"),Components(*mBackground),LineStyle(kDashed),LineColor(kRed),Range(plotrange),NormRange(range));
  mTemplate->plotOn(plot,Name("sig"),Components(*mSignal),LineStyle(kDashed),LineColor(kGreen+3),Range(plotrange),NormRange(range));
  mTemplate->paramOn(plot,Label(Form("#chi^{2}/NDF = %2.4f",mChi2)),Layout(0.64,0.92,0.86));
  plot->getAttLine()->SetLineWidth(0);
  data.removeSelfFromDir();
  return plot;
}

void FitModule::UseBackground(bool useBkg){
  mBkgCounts->setConstant(!useBkg);
  mBkgCounts->setVal(1000 * int(useBkg));
  auto iter = mBackground->getVariables()->createIterator();
  RooRealVar* var = (RooRealVar*)iter->Next();
  while(var){
    var->setConstant(!useBkg);
    var = (RooRealVar*)iter->Next();
  }
}

void FitModule::UseSignal(bool useSig){
  mSigCounts->setConstant(!useSig);
  mSigCounts->setVal(1000 * int(useSig));
  auto iter = mSignal->getVariables()->createIterator();
  RooRealVar* var = (RooRealVar*)iter->Next();
  while(var){
    var->setConstant(!useSig);
    var = (RooRealVar*)iter->Next();
  }
}

FitGausGaus::FitGausGaus(RooRealVar *x) : FitModule(x){
  mSigmaBkg = std::make_unique<RooRealVar>("#sigma_{bkg}","#sigma_{Bkg}",0.01,2.);
  mMuBkg = std::make_unique<RooRealVar>("#mu_{bkg}","#mu_{Bkg}",-5.,-3.);
  mBackground = std::make_unique<RooGaussian>("mBackground","Background",*mX,*mMuBkg,*mSigmaBkg);
  mSignal = std::make_unique<RooGaussian>("mSignal","Signal",*mX,*mMu,*mSigma);
  mTemplate = std::make_unique<RooAddPdf>("mTemplate","Template",RooArgList(*mSignal,*mBackground),
                                     RooArgList(*mSigCounts,*mBkgCounts));
}

FitLogNormalLogNormal::FitLogNormalLogNormal(RooRealVar *x) : FitModule(x){
  mSigma->setRange(1.0001, 20);
  mSigmaBkg = std::make_unique<RooRealVar>("#sigma_{bkg}","#sigma_{Bkg}",1.0001, 20);
  mMuBkg = std::make_unique<RooRealVar>("#mu_{bkg}","#mu_{Bkg}",-6.,-3.);
  mBackground = std::make_unique<RooLognormal>("mBackground","Background",*mX,*mMuBkg,*mSigmaBkg);
  mSignal = std::make_unique<RooLognormal>("mSignal","Signal",*mX,*mMu,*mSigma);
  mTemplate = std::make_unique<RooAddPdf>("mTemplate","Template",RooArgList(*mSignal,*mBackground),
                                     RooArgList(*mSigCounts,*mBkgCounts));
}

FitExpTailGaus::FitExpTailGaus(RooRealVar *x) : FitModule(x){
  mTau0 = std::make_unique<RooRealVar>("#tau_{0}","tau",-10.,-0.00001);
  mBackground = std::make_unique<RooExponential>("mBackground","Background",*mX,*mTau0);
  mAlpha0 = std::make_unique<RooRealVar>("#alpha_{0}","Alpha0",1.6,3.);
  mSignal = std::make_unique<RooGausExp>("mSignal","Signal",*mX,*mMu,*mSigma,*mAlpha0);
  mTemplate = std::make_unique<RooAddPdf>("mTemplate","Template",RooArgList(*mSignal,*mBackground),
                                     RooArgList(*mSigCounts,*mBkgCounts));
}

FitExpGaus::FitExpGaus(RooRealVar *x) : FitModule(x) {
  mTau = std::make_unique<RooRealVar>("mTau","tau bkg",-5.,0.);
  mBackground = std::make_unique<RooExponential>("mBackground","Background",*mX,*mTau);
  mSignal = std::make_unique<RooGaussian>("mSignal","Signal",*mX,*mMu,*mSigma);
  mTemplate = std::make_unique<RooAddPdf>("mTemplate","Template",RooArgList(*mSignal,*mBackground),
                                     RooArgList(*mSigCounts,*mBkgCounts));
}

FitExpCB::FitExpCB(RooRealVar *x) : FitModule(x) {
  mTau = std::make_unique<RooRealVar>("mTau","tau bkg",-5.,0.);
  mBackground = std::make_unique<RooExponential>("mBackground","Background",*mX,*mTau);
  mAlpha = std::make_unique<RooRealVar>("mAlpha","Alpha",-4.,-1.75);
  mN = std::make_unique<RooRealVar>("mN","n",3.,10.);
  mSignal = std::make_unique<RooCBShape>("mSignal","Signal",*mX,*mMu,*mSigma,*mAlpha,*mN);
  mTemplate = std::make_unique<RooAddPdf>("mTemplate","Template",RooArgList(*mSignal,*mBackground),
                                     RooArgList(*mSigCounts,*mBkgCounts));

}

FitExpExpCB::FitExpExpCB(RooRealVar *x) : FitModule(x) {
  mTau0 = std::make_unique<RooRealVar>("#tau_{0}","#tau_{0}",-6.5,-2.);
  mTau1 = std::make_unique<RooRealVar>("#tau_{1}","#tau_{1}",-2.,-0.2);
  mKbkg = std::make_unique<RooRealVar>("K_{bkg}","K_{bkg}",0.,0.,1.);
  mBkg0 = std::make_unique<RooExponential>("mBkg0","background1",*mX,*mTau0);
  mBkg1 = std::make_unique<RooExponential>("mBkg1","background2",*mX,*mTau1);
  mBackground = std::make_unique<RooAddPdf>("mBackground","Background",RooArgList(*mBkg0,*mBkg1),
                                       RooArgList(*mKbkg));
  mAlpha = std::make_unique<RooRealVar>("mAlpha","Alpha",-3.,-0.8);
  mN = std::make_unique<RooRealVar>("mN","n",3.5,40.);
  mSignal = std::make_unique<RooCBShape>("mSignal","Signal",*mX,*mMu,*mSigma,*mAlpha,*mN);
  mTemplate = std::make_unique<RooAddPdf>("mTemplate","Template",RooArgList(*mSignal,*mBackground),
                                     RooArgList(*mSigCounts,*mBkgCounts));
}

FitExpExpGaus::FitExpExpGaus(RooRealVar *x) : FitModule(x) {
  mTau0 = std::make_unique<RooRealVar>("#tau_{0}","#tau_{0}",-6.5,-2.);
  mTau1 = std::make_unique<RooRealVar>("#tau_{1}","#tau_{1}",-2.,-0.2);
  mKbkg = std::make_unique<RooRealVar>("K_{bkg}","K_{bkg}",0.,0.,1.);
  mBkg0 = std::make_unique<RooExponential>("mBkg0","background1",*mX,*mTau0);
  mBkg1 = std::make_unique<RooExponential>("mBkg1","background2",*mX,*mTau1);
  mBackground = std::make_unique<RooAddPdf>("mBackground","Background",RooArgList(*mBkg0,*mBkg1),
                                       RooArgList(*mKbkg));
  mSignal = std::make_unique<RooGaussian>("mSignal","Signal",*mX,*mMu,*mSigma);
  mTemplate = std::make_unique<RooAddPdf>("mTemplate","Template",RooArgList(*mSignal,*mBackground),
                                     RooArgList(*mSigCounts,*mBkgCounts));
}

FitExpExpTailGaus::FitExpExpTailGaus(RooRealVar *x) : FitModule(x) {
  mTau0 = std::make_unique<RooRealVar>("#tau_{0}","#tau_{0}",-10.,-0.5);
  mTau1 = std::make_unique<RooRealVar>("#tau_{1}","#tau_{1}",-0.5,-0.01);
  mKbkg = std::make_unique<RooRealVar>("K_{bkg}","K_{bkg}",0.,0.,1.);
  mBkg0 = std::make_unique<RooExponential>("mBkg0","background1",*mX,*mTau0);
  mBkg1 = std::make_unique<RooExponential>("mBkg1","background2",*mX,*mTau1);
  mBackground = std::make_unique<RooAddPdf>("mBackground","Background",RooArgList(*mBkg0,*mBkg1),
                                       RooArgList(*mKbkg));
  mAlpha0 = std::make_unique<RooRealVar>("#alpha","Alpha0",1.6,3.); // tight range based on low pT
  mSignal = std::make_unique<RooGausExp>("mSignal","Signal",*mX,*mMu,*mSigma,*mAlpha0);
  mTemplate = std::make_unique<RooAddPdf>("mTemplate","Template",RooArgList(*mSignal,*mBackground),
                                     RooArgList(*mSigCounts,*mBkgCounts));
}

FitExpExpTailTailGaus::FitExpExpTailTailGaus(RooRealVar *x) : FitModule(x) {
  mTau0 = std::make_unique<RooRealVar>("#tau_{0}","#tau_{0}",-6.5,-2.);
  mTau1 = std::make_unique<RooRealVar>("#tau_{1}","#tau_{1}",-2.,-0.001);
  mKbkg = std::make_unique<RooRealVar>("K_{bkg}","K_{bkg}",0.,0.,1.);
  mBkg0 = std::make_unique<RooExponential>("mBkg0","background1",*mX,*mTau0);
  mBkg1 = std::make_unique<RooExponential>("mBkg1","background2",*mX,*mTau1);
  mBackground = std::make_unique<RooAddPdf>("mBackground","Background",RooArgList(*mBkg0,*mBkg1),
                                       RooArgList(*mKbkg));
  mAlpha0 = std::make_unique<RooRealVar>("#alpha_{0}","Alpha0",-3.,1.); // tight range based on low pT
  mAlpha1 = std::make_unique<RooRealVar>("#alpha_{1}","Alpha1",1.,3.); // tight range based on low pT
  mSignal = std::make_unique<RooGausDExp>("mSignal","Signal",*mX,*mMu,*mSigma,*mAlpha0,*mAlpha1);
  mTemplate = std::make_unique<RooAddPdf>("mTemplate","Template",RooArgList(*mSignal,*mBackground),
                                     RooArgList(*mSigCounts,*mBkgCounts));
}

FitGausExpTailGaus::FitGausExpTailGaus(RooRealVar *x) : FitModule(x) {
  mTau0 = std::make_unique<RooRealVar>("#tau_{0}","#tau_{0}",-6.5,-0.1);
  mMuBkg = std::make_unique<RooRealVar>("mMuBkg","#a_{0}",-4,-7.,-4.5);
  mSigmaBkg = std::make_unique<RooRealVar>("mSigmaBkg","#a_{1}",0.01,2.5);
  mKbkg = std::make_unique<RooRealVar>("K_{bkg}","K_{bkg}",0.,0.,1.);
  mBkg0 = std::make_unique<RooExponential>("mBkg0","background1",*mX,*mTau0);
  mBkg1 = std::make_unique<RooGaussian>("mBkg1","background2",*mX,*mMuBkg,*mSigmaBkg);
  mBackground = std::make_unique<RooAddPdf>("mBackground","Background",RooArgList(*mBkg0,*mBkg1),
                                       RooArgList(*mKbkg));
  mAlpha0 = std::make_unique<RooRealVar>("#alpha_{0}","Alpha0",1.,3.); // tight range based on low pT
  mSignal = std::make_unique<RooGausExp>("mSignal","Signal",*mX,*mMu,*mSigma,*mAlpha0);
  mTemplate = std::make_unique<RooAddPdf>("mTemplate","Template",RooArgList(*mSignal,*mBackground),
                                     RooArgList(*mSigCounts,*mBkgCounts));
}

FitExpPolTailGaus::FitExpPolTailGaus(RooRealVar *x) : FitModule(x) {
  mTau0 = std::make_unique<RooRealVar>("#tau_{0}","#tau_{0}",-6.5,-2.);
  mA0 = std::make_unique<RooRealVar>("mA0","#a_{0}",-0.15,-2.,-0.1);
  mA1 = std::make_unique<RooRealVar>("mA1","#a_{1}",-2.,-0.1);
  mKbkg = std::make_unique<RooRealVar>("K_{bkg}","K_{bkg}",0.,0.,1.);
  mBkg0 = std::make_unique<RooExponential>("mBkg0","background1",*mX,*mTau0);
  mBkg1 = std::make_unique<RooChebychev>("mBkg1","background2",*mX,RooArgList(*mA0));
  mBackground = std::make_unique<RooAddPdf>("mBackground","Background",RooArgList(*mBkg0,*mBkg1),
                                       RooArgList(*mKbkg));
  mAlpha0 = std::make_unique<RooRealVar>("#alpha_{0}","Alpha0",1.6,3.); // tight range based on low pT
  mSignal = std::make_unique<RooGausExp>("mSignal","Signal",*mX,*mMu,*mSigma,*mAlpha0);
  mTemplate = std::make_unique<RooAddPdf>("mTemplate","Template",RooArgList(*mSignal,*mBackground),
                                     RooArgList(*mSigCounts,*mBkgCounts));
}
