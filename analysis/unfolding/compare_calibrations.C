// Compare unfold_data_wEff outputs. No RooUnfold library is needed.
// ROOT: .x compare_calibrations.C+("out_OrigCalib","out_HSCalib","calibration_comparison")
// Optional fourth argument: histogram key (default: final invariant spectrum).
// Ratio errors assume independent inputs; cross-calibration covariance is unavailable.
#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TSystem.h"
#include "TString.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>

namespace CalibrationComparison {
bool sameBinning(const TH1* a, const TH1* b) {
  if (a->GetNbinsX()!=b->GetNbinsX()) return false;
  for (int i=1;i<=a->GetNbinsX()+1;++i) {
    const double x=a->GetXaxis()->GetBinLowEdge(i), y=b->GetXaxis()->GetBinLowEdge(i);
    if (std::fabs(x-y)>1e-9*std::max(1.,std::max(std::fabs(x),std::fabs(y)))) return false;
  }
  return true;
}
void style(TH1* h, int color, int marker) {
  h->SetStats(false); h->SetTitle("");
  h->SetLineColor(color); h->SetMarkerColor(color);
  h->SetMarkerStyle(marker); h->SetMarkerSize(1.0); h->SetLineWidth(2);
}
}

void compare_calibrations(const char* origFolder, const char* hsFolder,
                          const char* outputFolder="calibration_comparison",
                          const char* histKey="hUnfoldedTruthBins_finalInvariant") {
  using namespace CalibrationComparison;
  const std::string origPath=std::string(origFolder)+"/unfolded_data.root";
  const std::string hsPath=std::string(hsFolder)+"/unfolded_data.root";
  std::unique_ptr<TFile> orig(TFile::Open(origPath.c_str(),"READ"));
  std::unique_ptr<TFile> hs(TFile::Open(hsPath.c_str(),"READ"));
  if (!orig.get() || orig->IsZombie() || !hs.get() || hs->IsZombie()) {
    std::cerr << "Cannot open inputs:\n" << origPath << "\n" << hsPath << std::endl; return;
  }
  if (gSystem->mkdir(outputFolder,true)!=0 && gSystem->AccessPathName(outputFolder)) {
    std::cerr << "Cannot create " << outputFolder << std::endl; return;
  }
  const std::string outputPath=std::string(outputFolder)+"/calibration_comparison.root";
  if (outputPath==origPath || outputPath==hsPath) return;
  std::unique_ptr<TFile> output(TFile::Open(outputPath.c_str(),"RECREATE"));
  if (!output.get() || output->IsZombie()) {
    std::cerr << "Cannot create " << outputPath << std::endl; return;
  }
  const double radii[]={0.2,0.3,0.4};
  const char* cents[]={"CENT_0_10","MID_20_40","PERI_60_80"};
  const char* labels[]={"0-10%","20-40%","60-80%"};
  const int leads[]={0,5,7,9};
  int made=0, skipped=0;
  for (int ir=0;ir<3;++ir) for (int ic=0;ic<3;++ic) for (int il=0;il<4;++il) {
    const TString tag=Form("R%.1f_%s_ptlead%d",radii[ir],cents[ic],leads[il]);
    const TString path=tag+"/"+histKey;
    TH1* a=dynamic_cast<TH1*>(orig->Get(path));
    TH1* b=dynamic_cast<TH1*>(hs->Get(path));
    if (!a || !b || a->GetDimension()!=1 || b->GetDimension()!=1 || !sameBinning(a,b)) {
      std::cerr << "Skip " << path << ": missing histogram or incompatible binning.\n";
      ++skipped; continue;
    }
    std::unique_ptr<TH1> ho(static_cast<TH1*>(a->Clone("OrigCalib_"+tag)));
    std::unique_ptr<TH1> hh(static_cast<TH1*>(b->Clone("HSCalib_"+tag)));
    ho->SetDirectory(0); hh->SetDirectory(0);
    style(ho.get(),kRed+1,24); style(hh.get(),kBlue+1,20);
    TGraphErrors go, gh, ratio;
    go.SetName("spectrum_OrigCalib"); gh.SetName("spectrum_HSCalib"); ratio.SetName("ratio_OrigCalib_over_HSCalib");
    double ymin=1e300,ymax=0., rmin=1.,rmax=1.;
    int invalid=0;
    for (int ib=1;ib<=ho->GetNbinsX();++ib) {
      const double x=ho->GetBinCenter(ib);
      const double va=ho->GetBinContent(ib), vb=hh->GetBinContent(ib);
      const double ea=ho->GetBinError(ib), eb=hh->GetBinError(ib);
      if (std::isfinite(va) && std::isfinite(ea) && va>0.) {
        int n=go.GetN(); go.SetPoint(n,x,va); go.SetPointError(n,0.,ea);
        ymin=std::min(ymin,va); ymax=std::max(ymax,va+ea);
      }
      if (std::isfinite(vb) && std::isfinite(eb) && vb>0.) {
        int n=gh.GetN(); gh.SetPoint(n,x,vb); gh.SetPointError(n,0.,eb);
        ymin=std::min(ymin,vb); ymax=std::max(ymax,vb+eb);
      }
      if (!std::isfinite(va) || !std::isfinite(vb) || !std::isfinite(ea) || !std::isfinite(eb) || vb<=0.) { ++invalid; continue; }
      const double v=va/vb;
      const double e=std::hypot(ea/vb,v*eb/vb);
      if (!std::isfinite(v) || !std::isfinite(e)) { ++invalid; continue; }
      int n=ratio.GetN(); ratio.SetPoint(n,x,v); ratio.SetPointError(n,0.,e);
      rmin=std::min(rmin,v-e); rmax=std::max(rmax,v+e);
    }
    if (ymax<=0.) { std::cerr << "Skip " << tag << ": no positive spectrum bins.\n"; ++skipped; continue; }
    go.SetLineColor(kRed+1); go.SetMarkerColor(kRed+1); go.SetMarkerStyle(24);
    gh.SetLineColor(kBlue+1); gh.SetMarkerColor(kBlue+1); gh.SetMarkerStyle(20);
    ratio.SetLineColor(kBlue+1); ratio.SetMarkerColor(kBlue+1); ratio.SetMarkerStyle(20);
    TCanvas canvas("comparison_"+tag,"Calibration comparison",850,850);
    TPad top("top","",0.,0.30,1.,1.), bottom("bottom","",0.,0.,1.,0.30);
    top.SetLeftMargin(.17); top.SetRightMargin(.04); top.SetTopMargin(.05); top.SetBottomMargin(.025);
    bottom.SetLeftMargin(.17); bottom.SetRightMargin(.04); bottom.SetTopMargin(.025); bottom.SetBottomMargin(.30);
    top.Draw(); bottom.Draw(); top.cd(); top.SetLogy();
    std::unique_ptr<TH1> frame(static_cast<TH1*>(ho->Clone("spectrum_frame_"+tag)));
    frame->SetDirectory(0); frame->Reset(); frame->SetMinimum(ymin*.3); frame->SetMaximum(ymax*30.);
    frame->GetXaxis()->SetLabelSize(0.); frame->GetXaxis()->SetTitle("");
    frame->GetYaxis()->SetTitle(std::string(histKey)=="hUnfoldedTruthBins_finalInvariant" ?
      "(1/N_{evt})(1/2#pi p_{T}) d^{2}N/(dp_{T}d#eta) [(GeV/c)^{-2}]" : a->GetYaxis()->GetTitle());
    frame->GetYaxis()->SetTitleSize(.040); frame->GetYaxis()->SetLabelSize(.040); frame->GetYaxis()->SetTitleOffset(1.8);
    frame->Draw("AXIS"); go.Draw("P SAME"); gh.Draw("P SAME");
    TLegend legend(.56,.77,.93,.91); legend.SetBorderSize(0); legend.SetFillStyle(0); legend.SetTextSize(.038);
    legend.AddEntry(&go,"Original calibration","pe"); legend.AddEntry(&gh,"Hanseul calibration (baseline)","pe"); legend.Draw();
    TLatex text; text.SetNDC(); text.SetTextSize(.038);
    text.DrawLatex(.21,.91,Form("Au+Au %s, #sqrt{s_{NN}} = 200 GeV",labels[ic]));
    text.DrawLatex(.21,.85,Form("anti-#it{k}_{T}, R = %.1f",radii[ir]));
    text.DrawLatex(.21,.79,Form("#it{p}_{T,min}^{lead} = %d GeV/#it{c}",leads[il]));
    bottom.cd();
    std::unique_ptr<TH1> rf(static_cast<TH1*>(ho->Clone("ratio_frame_"+tag)));
    rf->SetDirectory(0); rf->Reset();
    const double span=std::max(.2,rmax-rmin);
    rf->SetMinimum(std::min(.8,rmin-.15*span)); rf->SetMaximum(std::max(1.2,rmax+.15*span));
    rf->GetXaxis()->SetTitle("#it{p}_{T,jet} [GeV/#it{c}]");
    rf->GetYaxis()->SetTitle("OrigCalib / HSCalib");
    rf->GetXaxis()->SetTitleSize(.105); rf->GetXaxis()->SetLabelSize(.09); rf->GetXaxis()->SetTitleOffset(1.1);
    rf->GetYaxis()->SetTitleSize(.09); rf->GetYaxis()->SetLabelSize(.085); rf->GetYaxis()->SetTitleOffset(.85); rf->GetYaxis()->SetNdivisions(505);
    rf->Draw("AXIS");
    TLine unity(rf->GetXaxis()->GetXmin(),1.,rf->GetXaxis()->GetXmax(),1.); unity.SetLineStyle(2); unity.Draw(); ratio.Draw("P SAME");
    canvas.cd(); canvas.SaveAs((std::string(outputFolder)+"/compare_"+tag.Data()+".png").c_str());
    output->cd(); TDirectory* dir=output->mkdir(tag); dir->cd();
    ho->Write("OrigCalib"); hh->Write("HSCalib"); ratio.Write(); canvas.Write();
    ++made;
    if (invalid) std::cout << tag << ": omitted " << invalid << " undefined ratio bins.\n";
  }
  output->Close();
  std::cout << "Saved " << made << " comparisons; skipped " << skipped << ". Output: " << outputFolder << std::endl;
}
