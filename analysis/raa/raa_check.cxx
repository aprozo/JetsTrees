// Preliminary inclusive-jet R_AA from invariant Au+Au yields and pp cross sections.
// ROOT usage: .x raa_check.cxx+("unfolded_data.root","pp_reference","out_raa",5.0)
#include "TFile.h"
#include "TDirectory.h"
#include "TH1.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TParameter.h"
#include "TString.h"
#include "TMath.h"
#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include <algorithm>

namespace {
struct Cent {
  const char* tag;
  const char* label;
  double taa;       // <T_AA> in mb^{-1}
  double taaErr;
};

const Cent cents[] = {
  {"CENT_0_10",  "0-10%",  22.8, 1.6},
  {"PERI_60_80", "60-80%", 0.49, 0.14}
};
const double radii[] = {0.2, 0.3, 0.4};
void ensure(const TString& s) { gSystem->mkdir(s.Data(), true); }
bool IsFinite(double x) { return TMath::Finite(x); }
TH1D* copyAsTH1D(const TH1* src, const TString& name) {
  if (!src || src->GetDimension()!=1) return 0;
  const int n=src->GetNbinsX();
  std::vector<double> edges(n+1);
  for (int i=1;i<=n;++i) edges[i-1]=src->GetXaxis()->GetBinLowEdge(i);
  edges[n]=src->GetXaxis()->GetBinUpEdge(n);
  TH1D* out=new TH1D(name,"",n,&edges[0]);
  out->SetDirectory(0);
  for (int i=1;i<=n;++i) {
    out->SetBinContent(i,src->GetBinContent(i));
    out->SetBinError(i,src->GetBinError(i));
  }
  return out;
}
int matchingBin(const TH1* h, double low, double high) {
  for(int i=1;i<=h->GetNbinsX();++i) {
    const double lo=h->GetXaxis()->GetBinLowEdge(i);
    const double hi=h->GetXaxis()->GetBinUpEdge(i);
    if(std::fabs(lo-low)<1.e-5 && std::fabs(hi-high)<1.e-5) return i;
  }
  return -1;
}
}

void raa_check(const char* unfoldedFile,
               const char* ppDir="analysis/pp_reference",
               const char* outDir="out_raa",
               double minPlotPt=5.0)
{               
  TH1::AddDirectory(kFALSE);
  gStyle->SetOptStat(0);
  ensure(outDir); ensure(TString(outDir)+"/png");
  TFile* fAu=TFile::Open(unfoldedFile,"READ");
  if(!fAu || fAu->IsZombie()) { std::cerr<<"[ERROR] cannot open AuAu "<<unfoldedFile<<"\n"; return; }
  TString output=TString(outDir)+"/raa.root";
  TFile* fOut=TFile::Open(output,"RECREATE");
  if(!fOut || fOut->IsZombie()) { std::cerr<<"[ERROR] cannot write "<<output<<"\n"; fAu->Close(); return; }
  TParameter<double>("ptlead_cut_GeV",0.0).Write();
  int nSuccess=0;
  const int nCents = sizeof(cents)/sizeof(cents[0]);

  for (int iR=0; iR<3; ++iR) {
    const double R=radii[iR];
    TString rtag=Form("R%.1f",R);
    TString ppPath=Form("%s/xsec_R%.1f_noUE_mb_jp1.root",ppDir,R);
    TFile* fPP=TFile::Open(ppPath,"READ");
    if(!fPP || fPP->IsZombie()) {
      std::cerr<<"[WARN] cannot open pp reference "<<ppPath<<"\n";
      if(fPP) {fPP->Close();delete fPP;}
      continue;
    }
    TH1* ppStat=dynamic_cast<TH1*>(fPP->Get("xsec_stat"));
    TH1* ppSyst=dynamic_cast<TH1*>(fPP->Get("xsec_syst"));
    if(!ppStat || ppStat->GetDimension()!=1) {
      std::cerr<<"[WARN] missing 1D xsec_stat in "<<ppPath<<"\n";
      fPP->Close();delete fPP;continue;
    }
    for (int iC=0; iC<nCents; ++iC) {
      const Cent& cent=cents[iC];
      TString tag=Form("%s_%s_ptlead0",rtag.Data(),cent.tag);
      TString auPath=tag+"/hUnfoldedTruthBins_finalInvariant";
      TH1* au=dynamic_cast<TH1*>(fAu->Get(auPath));
      if(!au || au->GetDimension()!=1) {
        std::cerr<<"[WARN] missing 1D AuAu spectrum "<<auPath<<"\n";
        continue;
      }
      TH1D* hRaa=copyAsTH1D(au,"hRaa");
      TH1D* hAu=copyAsTH1D(au,"hAuAuInvariant");
      TH1D* hPP=copyAsTH1D(ppStat,"hPpCrossSection");
      hRaa->Reset("ICES");
      int good=0, unmatched=0;
      for(int i=1;i<=au->GetNbinsX();++i) {
        double lo=au->GetXaxis()->GetBinLowEdge(i);
        double hi=au->GetXaxis()->GetBinUpEdge(i);
        if(lo<minPlotPt-1.e-6) continue;
        int j=matchingBin(ppStat,lo,hi);
        if(j<0) { ++unmatched; std::cerr<<"[WARN] unmatched pp bin: "<<tag<<" ["<<lo<<","<<hi<<"]\n"; continue; }
        double a=au->GetBinContent(i),ea=au->GetBinError(i);
        double p=ppStat->GetBinContent(j),ep=ppStat->GetBinError(j);
        double pt=au->GetXaxis()->GetBinCenter(i);
        if(!IsFinite(a)||!IsFinite(ea)||!IsFinite(p)||!IsFinite(ep)||p<=0||a<0||ea<0||ep<0) {
          std::cerr<<"[WARN] invalid bin "<<tag<<" ["<<lo<<","<<hi<<"]\n";continue;
        }
        const double taaPbInv = cent.taa * 1.e-9;  // mb^-1 -> pb^-1

        const double scale =
            (2.0 * TMath::Pi() * pt) / taaPbInv;

        const double value = scale * a / p;

        const double error =
            scale * std::sqrt(
                (ea/p)*(ea/p)
                +
                (a*ep/(p*p))*(a*ep/(p*p))
            );
        hRaa->SetBinContent(i,value);
        hRaa->SetBinError(i,error);
        ++good;
      }
      if(!good) { std::cerr<<"[WARN] no comparable bins: "<<tag<<"\n";delete hRaa;delete hAu;delete hPP;continue; }
      fOut->cd();
      TDirectory* d=fOut->mkdir(tag);
      d->cd();
      hRaa->Write("hRaa");hAu->Write("hAuAuInvariant");hPP->Write("hPpCrossSection");
      if(ppSyst && ppSyst->GetDimension()==1) {
        TH1D* hs=copyAsTH1D(ppSyst,"hPpSystInput");hs->Write();delete hs;
      }
      TParameter<double>("TAA_mbInv",cent.taa).Write();
      TParameter<double>("TAAErr_mbInv",cent.taaErr).Write();
      TCanvas* c=new TCanvas(Form("canvas_%s",tag.Data()),"",850,650);
      c->SetLeftMargin(.13);c->SetBottomMargin(.12);c->SetRightMargin(.05);c->SetTopMargin(.07);
      TH1D* frame = (TH1D*)hRaa->Clone(Form("plotFrame_%s",tag.Data()));
      frame->SetDirectory(0);
      frame->Reset("ICES");
      frame->SetTitle("");
      frame->GetXaxis()->SetTitle("#it{p}_{T,jet} (GeV/#it{c})");
      frame->GetYaxis()->SetTitle("#it{R}_{AA}");
      frame->GetYaxis()->SetRangeUser(0.,2.0);
      frame->GetXaxis()->SetRangeUser(minPlotPt,60.);
      frame->Draw("AXIS");
      TLine line(minPlotPt,1.,60.,1.);line.SetLineStyle(2);line.Draw("SAME");
      // Draw only validated bins; unmatched bins must not appear as zero-valued points.
      TGraphErrors* graph=new TGraphErrors();
      for(int i=1;i<=hRaa->GetNbinsX();++i) {
        double lo=hRaa->GetXaxis()->GetBinLowEdge(i);
        double hi=hRaa->GetXaxis()->GetBinUpEdge(i);
        if(lo<minPlotPt-1.e-6 || matchingBin(ppStat,lo,hi)<0) continue;
        if(ppStat->GetBinContent(matchingBin(ppStat,lo,hi))<=0) continue;
        double y=hRaa->GetBinContent(i), ey=hRaa->GetBinError(i);
        if(!IsFinite(y)||!IsFinite(ey)) continue;
        int k=graph->GetN();graph->SetPoint(k,(lo+hi)/2.,y);graph->SetPointError(k,0.,ey);
      }
      graph->SetMarkerStyle(20);graph->SetMarkerSize(1.1);graph->SetMarkerColor(kBlue+1);graph->SetLineColor(kBlue+1);
      graph->Draw("P SAME");
      TLatex lat;lat.SetNDC();lat.SetTextSize(.035);
      lat.DrawLatex(.16,.91,Form("Au+Au  #sqrt{s_{NN}} = 200 GeV, %s",cent.label));
      lat.DrawLatex(.16,.85,Form("anti-#it{k}_{T}, #it{R} = %.1f, no leading cut",R));
      lat.DrawLatex(.16,.79,Form("#LT#it{T}_{AA}#GT = %.2f #pm %.2f mb^{-1}",cent.taa,cent.taaErr));      
      c->SaveAs(Form("%s/png/RAA_%s.png",outDir,tag.Data()));
      delete c;
      delete graph;
      delete frame;
      std::cout<<"[OK] "<<tag<<": "<<good<<" matched bins, "<<unmatched<<" unmatched\n";
      ++nSuccess;
      delete hRaa;delete hAu;delete hPP;
    }
    fPP->Close();delete fPP;
  }
  fOut->Close();delete fOut;fAu->Close();delete fAu;
  std::cout<<"[DONE] "<<nSuccess<<" spectra written to "<<output<<"\n";
}
