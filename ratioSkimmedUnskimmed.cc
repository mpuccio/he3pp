#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>

void ratioSkimmedUnskimmed() {
  TFile *_file0 = TFile::Open("../output/LHC23/apass4_skimmed/signal_skimmed.root");
  TFile *_file1 = TFile::Open("../output/LHC22/apass4/signal_unskimmedcomp.root");
  TH1* skimmed = (TH1*)_file0->Get("nuclei/antihe3/TPConly/hTPConlyA0_ExpGaus");
  TH1* unskimmed = (TH1*)_file1->Get("nuclei/antihe3/TPConly/hTPConlyA0_ExpGaus");
  TH1* ratioH = (TH1*)skimmed->Clone("ratio");
  double nVTXskimmed = 1.4e11;
  double nVTXunskimmed = 5.4364867e11;
  ratioH->Divide(unskimmed);
  ratioH->SetTitle(";#it{p}_{T} (GeV/#it{c});Skimmed 2023 / Unskimmed 2022;");
  ratioH->SetMarkerStyle(20);
  ratioH->SetMarkerSize(0.8);
  ratioH->SetMarkerColor(kBlack);
  ratioH->SetLineColor(kBlack);
  ratioH->SetStats(false);
  ratioH->GetXaxis()->SetRange(2, ratioH->GetXaxis()->GetNbins());

  skimmed->Scale(1. / nVTXskimmed, "width");
  skimmed->SetTitle("2023 skimmed;#it{p}_{T} (GeV/#it{c});1/#it{N}_{TVX} #times d#it{N}_{raw}/d#it{p}_{T} (#it{c}/GeV);");
  unskimmed->Scale(1. / nVTXunskimmed, "width");
  unskimmed->SetTitle("2022 unskimmed;#it{p}_{T} (GeV/#it{c});1/#it{N}_{TVX} #times d#it{N}_{raw}/d#it{p}_{T} (#it{c}/GeV);");
  skimmed->SetMarkerStyle(20);
  skimmed->SetMarkerSize(0.8);
  skimmed->SetMarkerColor(kRed);
  skimmed->SetLineColor(kRed);
  unskimmed->SetMarkerStyle(20);
  unskimmed->SetMarkerSize(0.8);
  unskimmed->SetMarkerColor(kBlue);
  unskimmed->SetLineColor(kBlue);

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
  skimmed->Draw("E1");
  unskimmed->Draw("E1 SAME");
  c1->BuildLegend();

  TCanvas *c = new TCanvas("c", "c", 800, 600);
  ratioH->SetDirectory(0);
  ratioH->Draw("E1");
  TLine *line = new TLine(ratioH->GetXaxis()->GetBinLowEdge(2), nVTXskimmed / nVTXunskimmed, ratioH->GetXaxis()->GetBinUpEdge(ratioH->GetXaxis()->GetNbins()), nVTXskimmed / nVTXunskimmed);
  line->SetLineColor(kRed);
  line->SetLineStyle(2);
  line->Draw();
  TLegend *leg = new TLegend(0.426065, 0.786087, 0.938596, 0.841739);
  leg->AddEntry(line, "Inspected TVX 2023 / 2022", "l");
  leg->Draw();
}