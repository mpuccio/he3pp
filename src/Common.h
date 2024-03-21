#ifndef COMMON_H
#define COMMON_H

#include "ROOT/RDataFrame.hxx"

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <TList.h>
#include <TF1.h>
#include <TObject.h>
#include <TMath.h>
#include <TSystem.h>

using std::map;
using std::string;
using std::vector;

const char   kLetter[2] = {'M','A'}; // M -> Matter, A -> Anti-matter
const string kNames[2] = {"he3","antihe3"};
const string kLabels[2] = {"^{3}He","^{3}#bar{He}"};

const string kMCproduction = "LHC23j6b";
const string kRecoPass = "apass4";
const string kPeriod = "LHC22";
const string kVariant = "_alpha";

const string kBaseOutputDir = "$NUCLEI_OUTPUT/" + kPeriod + "/" + kRecoPass + "/";
const string kBaseInputDir = "$NUCLEI_INPUT/";

const string kDataTreeFilename = kBaseInputDir + "data/" + kPeriod + "/" + kRecoPass + "/MergedAO2D.root";
const string kDataFilename = kBaseInputDir + "data/" + kPeriod + "/" + kRecoPass + "/DataHistos" + kVariant + ".root";
const string kDataFilenameHe4 = kBaseInputDir + "data/" + kPeriod + "/" + kRecoPass + "/DataHistosHe4" + kVariant + ".root";
const string kDataAnalysisResults = kBaseInputDir + "data/" + kPeriod + "/" + kRecoPass + "/AnalysisResults.root";
const string kMCtreeFilename = kBaseInputDir + "MC/" + kMCproduction + "/MergedAO2D.root";
const string kMCfilename = kBaseInputDir + "MC/" + kMCproduction + "/MChistos"  + kVariant + ".root";
const string kMCfilenameHe4 = kBaseInputDir + "MC/" + kMCproduction + "/MChistosHe4"  + kVariant + ".root";

const string kFilterListNames = "nuclei";

const string kSignalOutput = kBaseOutputDir + "signal" + kVariant + ".root";
const string kSystematicsOutput = kBaseOutputDir + "systematics" + kVariant + ".root";

constexpr int    kNPtBins = 4;
constexpr double  kPtBins[kNPtBins+1] = {0.5,1.5,2.0,2.5,3.0};

const int    kCentLength = 1;
const int    kCentBinsArray[kCentLength][2] = {{2,2}};
const float  kCentPtLimits[kCentLength] = {7};
const float  kCentLabels[kCentLength][2] = {{0.,100.}};

const float  kTPCmaxPt = 7.0f;
const float  kTOFminPt = 1.f;
const float  kPtRange[2] = {1.4,7};

const bool   kUseBarlow{true};
const float  kAbsSyst[2] = {0.08,0.1f};

constexpr int kNTPCfunctions = 3;
const std::string kTPCfunctName[4]{"GausGaus", "ExpGaus", "ExpTailGaus", "LognormalLognormal"};


std::map<string,vector<float> > kCutNames {{"nsigmaDCAz", {6, 7, 8}},{"fTPCnCls", {110, 120, 130}},{"nITScls", {5, 6, 7}}, {"nsigmaTPC", {3, 4, 5}}};
size_t nTrials{kCutNames["fDCAz"].size() * kCutNames["fTPCnCls"].size() * kCutNames["nITScls"].size()};

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
float nsigmaHe4(float mom, float sig) { return (sig / bbHe4(mom * 2) - 1.) / 0.07; }

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

auto defineColumnsForData(ROOT::RDataFrame& d) {
  return d.Define("ptUncorr", "2 * std::abs(fPt)")
          .Define("pt", "ptUncorr + 0.0343554 + 0.96161 * std::exp(-1.51286 * ptUncorr)")
          .Define("ptHe4step1", "ptUncorr + 0.0419608 + 1.75861 * std::exp(-1.4019 * ptUncorr)")
          .Define("ptHe4", "ptHe4step1 + 0.00385223 - 0.442353 * std::exp(-1.59049 * ptHe4step1)")
          .Define("p", "pt * cosh(fEta)")
          .Define("tofMass", "fBeta < 1.e-3 ? 1.e9 : fBeta >= 1. ? 0 : fTPCInnerParam * 2 * sqrt(1.f / (fBeta * fBeta) - 1.f)")
          .Define("matter", "fPt > 0")
          .Define("pidForTracking", "fFlags >> 12")
          .Define("nsigmaHe3", nsigmaHe3, {"fTPCInnerParam", "fTPCsignal"})
          .Define("nsigmaH3", nsigmaH3, {"fTPCInnerParam", "fTPCsignal"})
          .Define("nsigmaHe4", nsigmaHe4, {"fTPCInnerParam", "fTPCsignal"})
          .Define("nITSclsIB", "int(0) + bool(fITSclsMap & 1) + bool(fITSclsMap & 2) + bool(fITSclsMap & 4)")
          .Define("nITScls", "nITSclsIB + bool(fITSclsMap & 8) + bool(fITSclsMap & 16) + bool(fITSclsMap & 32) + bool(fITSclsMap & 64)")
          .Define("hasTOF", "fFlags & (1 << 5)")
          .Define("isPrimary", "fFlags & (1 << 9)")
          .Define("isSecondaryFromMaterial", "fFlags & (1 << 10)")
          .Define("isSecondaryFromWeakDecay", "fFlags & (1 << 11)")
          .Define("deltaMassHe3", "tofMass - 2.80839")
          .Define("deltaMassHe4", "tofMass - 3.72738")
          .Define("hasGoodTOFmassHe3", "!hasTOF || std::abs(deltaMassHe3) < 0.6")
          .Define("hasGoodTOFmassHe4", "!hasTOF || std::abs(deltaMassHe4) < 0.3")
          .Define("nsigmaDCAxy", nSigmaDCAxy, {"pt", "fDCAxy"})
          .Define("nsigmaDCAz", nSigmaDCAz, {"pt", "fDCAz"});
}


#endif
