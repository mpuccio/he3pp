#include "src/Common.h"
#include <chrono>
#include <ctime>
#include <iomanip>
#include <sstream>

void CheckpointCreator()
{
  TFile systematics(kSystematicsOutput.data());
  TFile dataAR(kDataAnalysisResults.data());
  TFile mc(kMCfilename.data());
  TFile mcAR(kMCAnalysisResults.data());
  TFile signal(kSignalOutput.data());

  // Get the current time
  auto now = std::chrono::system_clock::now();
  std::time_t time = std::chrono::system_clock::to_time_t(now);
  std::tm *tm = std::localtime(&time);
  std::ostringstream oss;
  oss << "checkpoint-" << std::put_time(tm, "%d%m%y") << ".root";
  TFile checkpoint(oss.str().data(), "RECREATE");
  checkpoint.cd();
  systematics.Get("pubStat")->Clone("published_stat")->Write();
  systematics.Get("pubSyst")->Clone("published_syst")->Write();
  systematics.Get("fStatTPCA")->Clone("tpc_spectrum_stat")->Write();
  systematics.Get("fSystTPCA")->Clone("tpc_spectrum_syst")->Write();
  systematics.Get("fStatTOFA")->Clone("tof_spectrum_stat")->Write();
  systematics.Get("fSystTOFA")->Clone("tof_spectrum_syst")->Write();
  mc.Get("nuclei/effTPCA")->Clone("tpc_efficiency")->Write();
  mc.Get("nuclei/effTOFA")->Clone("tof_efficiency")->Write();
  std::cout << "Main dir done" << std::endl;
  checkpoint.mkdir("MC");
  checkpoint.cd("MC");
  mc.Get("nuclei/genAHe3")->Clone("generated")->Write();
  mc.Get("nuclei/TPCAHe3")->Clone("tpc_reconstructed")->Write();
  mc.Get("nuclei/TOFAHe3")->Clone("tpc_reconstructed")->Write();
  mcAR.Get("nuclei-spectra/spectra/hRecVtxZData")->Clone("events_reconstructed")->Write();
  std::cout << "MC dir done" << std::endl;
  checkpoint.mkdir("Data");
  checkpoint.cd("Data");
  dataAR.Get("nuclei-spectra/spectra/hRecVtxZData")->Clone("events_reconstructed")->Write();
  signal.Get("nuclei/antihe3/TPConly/hTPConlyA0_ExpGaus")->Clone("tpc_rawcounts")->Write();
  signal.Get("nuclei/antihe3/GausExp/hRawCountsA0")->Clone("tof_rawcounts")->Write();
  std::cout << "Data dir done" << std::endl;
}