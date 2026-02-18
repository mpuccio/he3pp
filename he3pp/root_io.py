import os
from pathlib import Path
from typing import Any

import ROOT

from . import settings as s

_DECLARED = False
_FIT_LOADED = False


def expand(path: str) -> str:
    return os.path.expandvars(os.path.expanduser(str(path)))


def ensure_parent(path: str) -> None:
    Path(path).expanduser().parent.mkdir(parents=True, exist_ok=True)


def declare_helpers() -> None:
    global _DECLARED
    if _DECLARED:
        return
    ROOT.gInterpreter.Declare(
        r'''
        #include <cmath>
        double he3pp_bb(double bg, double kp1, double kp2, double kp3, double kp4, double kp5) {
          double beta = bg / std::sqrt(1. + bg * bg);
          double aa = std::pow(beta, kp4);
          double bb = std::pow(1. / bg, kp5);
          bb = std::log(kp3 + bb);
          return (kp2 - aa - bb) * kp1 / aa;
        }
        float he3pp_bbHe3(float mom) { return he3pp_bb(mom / 2.80839, -321.34, 0.6539, 1.591, 0.8225, 2.363); }
        float he3pp_nsigmaHe3(float mom, float sig) { return (sig / he3pp_bbHe3(mom * 2) - 1. + 2.20376e-02) / 0.055; }
        float he3pp_bbH3(float mom) { return he3pp_bb(mom / 2.80892f, -136.71, 0.441, 0.2269, 1.347, 0.8035); }
        float he3pp_nsigmaH3(float mom, float sig) { return (sig / he3pp_bbH3(mom) - 1.) / 0.07; }
        float he3pp_bbHe4(float mom) { return he3pp_bb(mom / 3.72738f, -321.34, 0.6539, 1.591, 0.8225, 2.363); }
        float he3pp_nsigmaHe4(float mom, float sig) { return (sig / he3pp_bbHe4(mom * 2) - 1.) / 0.07; }
        float he3pp_DCAxyCut(float pt, float nsigma) {
          float invPt = 1.f / pt;
          return (7.62783e-04 + 4.59326e-03 * invPt + 6.89163e-03 * invPt * invPt) * nsigma;
        }
        float he3pp_nSigmaDCAxy(double pt, float dcaxy) { return dcaxy / he3pp_DCAxyCut(pt, 1); }
        float he3pp_DCAzCut(float pt, float nsigma) {
          float invPt = 1.f / pt;
          return (5.00000e-04 + 8.73690e-03 * invPt + 9.62329e-04 * invPt * invPt) * nsigma;
        }
        float he3pp_nSigmaDCAz(double pt, float dcaz) { return dcaz / he3pp_DCAzCut(pt, 1); }
        '''
    )
    _DECLARED = True


def load_fit_modules() -> None:
    global _FIT_LOADED
    if _FIT_LOADED:
        return
    ROOT.gROOT.ProcessLine(".L src/RooGausExp.cxx+")
    ROOT.gROOT.ProcessLine(".L src/RooGausDExp.cxx+")
    ROOT.gROOT.ProcessLine(".L src/FitModules.cxx+")
    _FIT_LOADED = True


def ptr(value: Any) -> Any:
    return value.get() if hasattr(value, "get") else value


def write_hist(obj: Any, name: str | None = None) -> None:
    hist = obj.GetValue() if hasattr(obj, "GetValue") else obj
    if name:
        hist.Write(name)
    else:
        hist.Write()


def define_columns_for_data(df: Any) -> Any:
    declare_helpers()
    return (
        df.Define("ptUncorr", "2 * std::abs(fPt)")
        .Define("pt", "ptUncorr + 0.0343554 + 0.96161 * std::exp(-1.51286 * ptUncorr)")
        .Define("ptHe4step1", "ptUncorr + 0.0419608 + 1.75861 * std::exp(-1.4019 * ptUncorr)")
        .Define("ptHe4", "ptHe4step1 + 0.00385223 - 0.442353 * std::exp(-1.59049 * ptHe4step1)")
        .Define("p", "pt * cosh(fEta)")
        .Define("tofMass", "fBeta < 1.e-3 ? 1.e9 : fBeta >= 1. ? 0 : fTPCInnerParam * 2 * sqrt(1.f / (fBeta * fBeta) - 1.f)")
        .Define("matter", "fPt > 0")
        .Define("pidForTracking", "fFlags >> 12")
        .Define("nsigmaHe3", "he3pp_nsigmaHe3(fTPCInnerParam, fTPCsignal)")
        .Define("nsigmaH3", "he3pp_nsigmaH3(fTPCInnerParam, fTPCsignal)")
        .Define("nsigmaHe4", "he3pp_nsigmaHe4(fTPCInnerParam, fTPCsignal)")
        .Define("nITSclsIB", "int(0) + bool(fITSclsMap & 1) + bool(fITSclsMap & 2) + bool(fITSclsMap & 4)")
        .Define("nITScls", "nITSclsIB + bool(fITSclsMap & 8) + bool(fITSclsMap & 16) + bool(fITSclsMap & 32) + bool(fITSclsMap & 64)")
        .Define("hasTOF", "fFlags & (1 << 5)")
        .Define("isPrimary", "fFlags & (1 << 9)")
        .Define("isSecondaryFromMaterial", "fFlags & (1 << 10)")
        .Define("isSecondaryFromWeakDecay", "fFlags & (1 << 11)")
        .Define("deltaMassHe3", "tofMass - 2.80839")
        .Define("deltaMassHe4", "tofMass - 3.72738")
        .Define("yHe3", "std::asinh(pt / std::hypot(pt, 2.80839) * std::sinh(fEta))")
        .Define("yHe4", "std::asinh(pt / std::hypot(pt, 3.72738) * std::sinh(fEta))")
        .Define("hasGoodTOFmassHe3", "!hasTOF || std::abs(deltaMassHe3) < 0.6")
        .Define("hasGoodTOFmassHe4", "!hasTOF || std::abs(deltaMassHe4) < 0.3")
        .Define("nsigmaDCAxy", "he3pp_nSigmaDCAxy(pt, fDCAxy)")
        .Define("nsigmaDCAz", "he3pp_nSigmaDCAz(pt, fDCAz)")
    )


def h2_model(name: str, title: str, y_bins: int, y_min: float, y_max: float) -> Any:
    return ROOT.RDF.TH2DModel(name, title, s.N_PT_BINS, s.PT_BIN_ARRAY, y_bins, y_min, y_max)


def h1_model(name: str, title: str) -> Any:
    return ROOT.RDF.TH1DModel(name, title, s.N_PT_BINS, s.PT_BIN_ARRAY)
