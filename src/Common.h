#ifndef COMMON_H
#define COMMON_H

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <TList.h>
#include <TF1.h>
#include <TObject.h>
#include <TMath.h>

using std::map;
using std::string;
using std::vector;

const char   kLetter[2] = {'M','A'}; // M -> Matter, A -> Anti-matter
const string kNames[2] = {"he3","antihe3"};
const string kLabels[2] = {"^{3}He","^{3}#bar{He}"};

const string kBaseOutputDir = "/Users/mpuccio/alice/Helium-pp/output/";
const string kBaseInputDir = "/Users/mpuccio/alice/Helium-pp/";

const string kDataFilename = kBaseInputDir + "data/dataMedium.root";
const string kMCfilename = kBaseInputDir + "MC/MCpub_newMC.root";

const string kFilterListNames = "nuclei";
const string kNormalisationList = "mpuccio_NucleiPIDqa";

const string kEfficiencyOutput = kBaseOutputDir + "efficiency.root";
const string kSignalOutput = kBaseOutputDir + "signalMedium.root";
const string kSecondariesOutput = kBaseOutputDir + "secondaries.root";
const string kSecondariesOutputRooFit = kBaseOutputDir + "RooSec.root";
const string kSecondariesTPCoutput = kBaseOutputDir + "secondaries_TPC.root";
const string kSecondariesTPCoutputRooFit = kBaseOutputDir + "RooSecTPC.root";
const string kMaterialOutput = kBaseOutputDir + "materialbudget.root";
const string kSpectraOutput = kBaseOutputDir + "spectra.root";
const string kFitSystematicsOutput = kBaseOutputDir + "fitsystematics.root";
const string kSystematicsOutput = kBaseOutputDir + "systematics.root";
const string kSystematicsOutputTPC = kBaseOutputDir + "systematics_TPC.root";
const string kFinalOutput = kBaseOutputDir + "final.root";

const bool   kPrintFigures{true};
const string kFiguresFolder = "/Users/lbariogl/cernbox/Deuterons13TeV/macros/results/images/";
const string kMacrosFolder = "/Users/lbariogl/cernbox/Deuterons13TeV/macros/results/images/macros/";

// constexpr int    kNPtBins = 11;
constexpr int    kNPtBins = 6;
// constexpr double  kPtBins[kNPtBins+1] = {1.4,1.8,2.2,2.6,3.0,3.4,3.8,4.4,5.2,6.0,7.0};
// constexpr double  kPtBins[kNPtBins+1] = {0.8,1.2,1.6,2.,2.4,2.8,3.2,3.6,4.,4.8,6.0,7.0};
constexpr double  kPtBins[kNPtBins+1] = {1.,1.5,2.0,2.5,3.0,4.0,5.0};
const float  kCentralityBins[2] = {0.f,100.f};
const int    kNCentralityBins = 1;

const int    kCentLength = 1;
const int    kCentBinsArray[kCentLength][2] = {{2,2}};
const float  kCentPtLimits[kCentLength] = {7};
const float  kCentLabels[kCentLength][2] = {{0.,100.}};

const float  kTPCmaxPt = 7.0f;
const float  kTOFminPt = 1.f;
const float  kPtRange[2] = {1.4,7};
const float  kPtRangeMatCorrection[2] = {1.05,1.55};
const float  kPtRangeMatCorrectionTPC[2] = {0.65,1.35};

const bool   kUseBarlow{true};
const bool   kSmoothSystematics{true};
const float  kAbsSyst[2] = {0.08,0.1f};

constexpr int kNTPCfunctions = 3;
const std::string kTPCfunctName[4]{"GausGaus", "ExpGaus", "ExpTailGaus", "LognormalLognormal"};


std::map<string,vector<float> > kCutNames {{"fDCAz", {6, 7, 8}},{"fTPCnCls", {110, 120, 130}},{"nITScls", {5, 6, 7}}, {"nsigmaTPC", {3, 4, 5}}};

double bb(double bg, double kp1, double kp2, double kp3, double kp4, double kp5)
{
  double beta = bg / std::sqrt(1. + bg * bg);
  double aa = std::pow(beta, kp4);
  double bb = std::pow(1. / bg, kp5);
  bb = std::log(kp3 + bb);
  return (kp2 - aa - bb) * kp1 / aa;
}

float bbHe3(float mom) { return bb(mom / 2.80839, -321.34, 0.6539, 1.591, 0.8225, 2.363); }
float nsigmaHe3(float mom, float sig) { return (sig / bbHe3(mom * 2) - 1. + 2.20376e-02) / 0.055; }

float bbH3(float mom) { return bb(mom / 2.80892f, -136.71, 0.441, 0.2269, 1.347, 0.8035); }
float nsigmaH3(float mom, float sig) { return (sig / bbH3(mom) - 1.) / 0.07; }

float bbHe4(float mom) { return bb(mom / 3.72738f, -321.34, 0.6539, 1.591, 0.8225, 2.363); }
float nsigmaHe4(float mom, float sig) { return (((sig / bbHe4(mom * 2) - 1.) / 0.07) + 0.5) / 1.4; }

float DCAxyCut(float pt, float nsigma)
{
  float invPt = 1.f / pt;
  return (7.62783e-04 + 4.59326e-03 * invPt + 6.89163e-03 * invPt * invPt) * nsigma;
}
float nSigmaDCAxy(double pt, float dcaxy) {
  return dcaxy / DCAxyCut(pt, 1);
}

float DCAzCut(float pt, float nsigma)
{
  float invPt = 1.f / pt;
  return (5.00000e-04 + 8.73690e-03 * invPt + 9.62329e-04 * invPt * invPt) * nsigma;
}

float nSigmaDCAz(double pt, float dcaz) {
  return dcaz / DCAzCut(pt, 1);
}



#endif
