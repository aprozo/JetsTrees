#include "TFile.h"
#include "TDirectory.h"
#include "TDirectoryFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TKey.h"
#include "TROOT.h"
#include "TString.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TPad.h"
#include "TLine.h"
#include "TGraphErrors.h"
#include "TMath.h"

#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>

// =========================================================
// Jet-quality cuts: keep synchronized with the tree maker
// =========================================================
const double CUT_AREA_02 = 0.07;
const double CUT_AREA_03 = 0.20;
const double CUT_AREA_04 = 0.40;
const double CUT_NEUTRAL_FRACTION = 0.95;
const double RECO_DUMMY_CUT = -500.0;

// Set to 2 to omit pThat bins 0 and 1.
static const int kFirstPtHatBinToUse = 0;

// Reco-leading-particle thresholds.
static const int N_LEAD = 4;
static const double PTLEAD_THR[N_LEAD] = {0.0, 5.0, 7.0, 9.0};

// pThat cross sections and generated-event counts.
static const int kNPthatBins = 11;
static const double kXsecWeights[kNPthatBins] = {
  1.616e+0, 1.355e-01, 2.288e-02, 5.524e-03, 2.203e-03,
  3.437e-04, 4.681e-05, 8.532e-06, 2.178e-06, 1.198e-07, 6.939e-09
};

static const double kNgenEvents[kNPthatBins] = {
  1020062, 1529646, 1275275, 1019532, 1019730,
  1020088, 1019739, 765165, 509510, 305922, 101971
};

// Same significance requirement as in make_hists.C:
// C / sigma_C > sqrt(10), evaluated separately in each pThat bin.
static const double kMinSignif = TMath::Sqrt(10.0);

// Reco binning used to construct the significance mask.
static const int nbins_meas = 24;
static const double bin_meas_edges[nbins_meas+1] = {
  -100,-80,-60,-40,-20,-10,-5,-2.5,0,2.5,5,7.5,10,12.5,15,17.5,
  20,22.5,25,27.5,30,35,40,50,60
};

// Common binning for the MC/reco comparison and truth bins for JES/JER.
static const int nbins_common = 10;
static const double bin_common_edges[nbins_common+1] = {
  0,5,10,15,20,25,30,35,40,50,60
};

double areaMinForR(double R)
{
  if (R < 0.25) return CUT_AREA_02;
  if (R < 0.35) return CUT_AREA_03;
  return CUT_AREA_04;
}

int FindPtHatBin(double xsecW)
{
  const double relTol = 1e-6;
  for (int i = 0; i < kNPthatBins; ++i) {
    const double ref = kXsecWeights[i];
    if (std::fabs(xsecW-ref) <= relTol*std::fabs(ref)) return i;
  }
  return -1;
}

int FindVariableBin(double x, int nbins, const double* edges)
{
  if (x < edges[0] || x >= edges[nbins]) return -1;
  for (int ib = 0; ib < nbins; ++ib) {
    if (x >= edges[ib] && x < edges[ib+1]) return ib;
  }
  return -1;
}

TString CentralityLabel(const std::string& cname)
{
  TString label = cname.c_str();
  label.ReplaceAll("CENT_", "");
  label.ReplaceAll("MID_",  "");
  label.ReplaceAll("PERI_", "");
  label.ReplaceAll("_", "-");
  label += " %";
  return label;
}

struct WeightedMoments {
  double sumW;
  double sumW2;
  double sumWX;
  double sumWX2;

  WeightedMoments()
  {
    sumW   = 0.0;
    sumW2  = 0.0;
    sumWX  = 0.0;
    sumWX2 = 0.0;
  }

  void Fill(double x, double w)
  {
    sumW   += w;
    sumW2  += w*w;
    sumWX  += w*x;
    sumWX2 += w*x*x;
  }

  bool Valid() const
  {
    if (sumW <= 0.0)  return false;
    if (sumW2 <= 0.0) return false;
    return true;
  }

  double Mean() const
  {
    if (!Valid()) return 0.0;
    return sumWX/sumW;
  }

  double Variance() const
  {
    if (!Valid()) return 0.0;

    const double mean = sumWX/sumW;
    const double variance = sumWX2/sumW - mean*mean;

    if (variance > 0.0) return variance;
    return 0.0;
  }

  double Sigma() const
  {
    return TMath::Sqrt(Variance());
  }

  double Neff() const
  {
    if (!Valid()) return 0.0;
    return sumW*sumW/sumW2;
  }

  double MeanError() const
  {
    const double neff = Neff();

    if (neff <= 0.0) return 0.0;

    return Sigma()/TMath::Sqrt(neff);
  }

  double SigmaError() const
  {
    const double neff = Neff();

    if (neff <= 1.0) return 0.0;

    return Sigma()/TMath::Sqrt(2.0*(neff-1.0));
  }
};

void StyleHistogram(TH1D* h, int marker, int color)
{
  h->SetMarkerStyle(marker);
  h->SetMarkerSize(1.0);
  h->SetMarkerColor(color);
  h->SetLineColor(color);
  h->SetLineWidth(2);
}

void make_jes_jer(const char* infile  = "embedding_merged_MCReco1p5.root",
                  const char* outfile = "jes_jer_hists.root",
                  const char* plotTop = "jes_jer_plots")
{
  TH1::SetDefaultSumw2(kTRUE);
  gStyle->SetOptStat(0);
  gROOT->SetBatch(kTRUE);

  gSystem->mkdir(plotTop, kTRUE);

  TFile* fin = TFile::Open(infile, "READ");
  if (!fin || fin->IsZombie()) {
    std::cerr << "Error: cannot open input file " << infile << std::endl;
    return;
  }

  TFile* fout = TFile::Open(outfile, "RECREATE");
  if (!fout || fout->IsZombie()) {
    std::cerr << "Error: cannot create output file " << outfile << std::endl;
    fin->Close();
    return;
  }

  TIter nextR(fin->GetListOfKeys());
  TKey* keyR = 0;

  while ((keyR = (TKey*)nextR())) {
    if (std::string(keyR->GetClassName()) != "TDirectoryFile") continue;

    const std::string rname = keyR->GetName();
    if (rname.empty() || rname[0] != 'R') continue;

    double R = 0.0;
    if (sscanf(rname.c_str(), "R%lf", &R) != 1) continue;

    const double dRmax = 0.6*R;
    const double areaMin = areaMinForR(R);

    TDirectoryFile* dirR = (TDirectoryFile*)keyR->ReadObj();
    if (!dirR) continue;

    fout->cd();
    TDirectory* outR = fout->mkdir(rname.c_str());
    if (!outR) outR = fout->GetDirectory(rname.c_str());

    TString plotRDir = Form("%s/%s", plotTop, rname.c_str());
    gSystem->mkdir(plotRDir.Data(), kTRUE);

    TIter nextC(dirR->GetListOfKeys());
    TKey* keyC = 0;

    while ((keyC = (TKey*)nextC())) {
      if (std::string(keyC->GetClassName()) != "TDirectoryFile") continue;

      const std::string cname = keyC->GetName();
      TDirectoryFile* dirC = (TDirectoryFile*)keyC->ReadObj();
      if (!dirC) continue;

      TTree* tree = (TTree*)dirC->Get("JetTree");
      if (!tree) continue;

      const char* requiredBranches[] = {
        "xsecWeight", "centralityWeight", "mc_pt", "deltaR",
        "reco_pt", "reco_pt_corr", "reco_area",
        "reco_neutral_fraction", "reco_trigger_match", "reco_pt_lead"
      };
      bool missing = false;
      for (unsigned int iRequiredBranch = 0;
          iRequiredBranch <
            sizeof(requiredBranches)/sizeof(requiredBranches[0]);
          ++iRequiredBranch) {

        if (!tree->GetBranch(requiredBranches[iRequiredBranch])) {
          std::cerr << "Missing branch "
                    << requiredBranches[iRequiredBranch]
                    << " in " << rname << "/" << cname
                    << std::endl;
          missing = true;
        }
      }
      if (missing) continue;

      Float_t xsecWeight = 0.f;
      Float_t centralityWeight = 1.f;
      Float_t mc_pt = -999.f;
      Float_t deltaR = -1.f;
      Float_t reco_pt = -999.f;
      Float_t reco_pt_corr = -999.f;
      Float_t reco_area = -999.f;
      Float_t reco_neutral_fraction = 999.f;
      Float_t reco_pt_lead = -999.f;
      Bool_t reco_trigger_match = kFALSE;

      tree->SetBranchAddress("xsecWeight", &xsecWeight);
      tree->SetBranchAddress("centralityWeight", &centralityWeight);
      tree->SetBranchAddress("mc_pt", &mc_pt);
      tree->SetBranchAddress("deltaR", &deltaR);
      tree->SetBranchAddress("reco_pt", &reco_pt);
      tree->SetBranchAddress("reco_pt_corr", &reco_pt_corr);
      tree->SetBranchAddress("reco_area", &reco_area);
      tree->SetBranchAddress("reco_neutral_fraction", &reco_neutral_fraction);
      tree->SetBranchAddress("reco_pt_lead", &reco_pt_lead);
      tree->SetBranchAddress("reco_trigger_match", &reco_trigger_match);

      // =====================================================
      // Pass 1: construct the same pThat-by-reco-bin mask as
      // make_hists.C, using the inclusive pTlead >= 0 sample.
      // =====================================================
      TH1D* hMaskSource[kNPthatBins];
      for (int ip = 0; ip < kNPthatBins; ++ip) {
        hMaskSource[ip] = new TH1D(
          Form("hMaskSource_pthat%d_%s_%s", ip, rname.c_str(), cname.c_str()),
          ";#it{p}_{T,reco}^{corr} [GeV/#it{c}];centrality-weighted counts",
          nbins_meas, bin_meas_edges
        );
      }

      const Long64_t nentries = tree->GetEntries();
      for (Long64_t ie = 0; ie < nentries; ++ie) {
        tree->GetEntry(ie);

        const int ip = FindPtHatBin((double)xsecWeight);
        if (ip < kFirstPtHatBinToUse || ip < 0) continue;

        const bool haveMC = (mc_pt > 0.0f);
        const bool haveReco = (reco_pt > (Float_t)RECO_DUMMY_CUT);
        if (!haveMC || !haveReco) continue;
        if (!(deltaR > 0.0f && deltaR < (Float_t)dRmax)) continue;
        if (!(reco_area >= (Float_t)areaMin)) continue;
        if (!(reco_neutral_fraction <= (Float_t)CUT_NEUTRAL_FRACTION)) continue;
        if (reco_trigger_match != kTRUE) continue;

        hMaskSource[ip]->Fill((double)reco_pt_corr,
                              (double)centralityWeight);
      }

      bool keepRecoBin[kNPthatBins][nbins_meas];
      for (int ip = 0; ip < kNPthatBins; ++ip) {
        for (int ix = 0; ix < nbins_meas; ++ix) {
          const double c = hMaskSource[ip]->GetBinContent(ix+1);
          const double e = hMaskSource[ip]->GetBinError(ix+1);
          const double signif = (e > 0.0) ? c/e : 0.0;
          keepRecoBin[ip][ix] = (signif > kMinSignif);
        }
      }

      // =====================================================
      // Histograms and residual moments for each pTlead cut.
      // Both spectra use identical selected matched entries.
      // =====================================================
      TH1D* hMC[N_LEAD];
      TH1D* hReco[N_LEAD];
      TH1D* hResidual[N_LEAD][nbins_common];
      WeightedMoments moments[N_LEAD][nbins_common];

      for (int it = 0; it < N_LEAD; ++it) {
        hMC[it] = new TH1D(
          Form("hMcPt_ptlead%.0f_%s_%s",
               PTLEAD_THR[it], rname.c_str(), cname.c_str()),
          ";#it{p}_{T,jet} [GeV/#it{c}];weighted counts / GeV",
          nbins_common, bin_common_edges
        );
        hReco[it] = new TH1D(
          Form("hRecoPtCorr_ptlead%.0f_%s_%s",
               PTLEAD_THR[it], rname.c_str(), cname.c_str()),
          ";#it{p}_{T,jet} [GeV/#it{c}];weighted counts / GeV",
          nbins_common, bin_common_edges
        );
          for (int iTruthBin = 0;
          iTruthBin < nbins_common;
          ++iTruthBin) {

          hResidual[it][iTruthBin] = new TH1D(
            Form("hResidual_ptlead%.0f_truthbin%d_%s_%s",
                PTLEAD_THR[it],
                iTruthBin,
                rname.c_str(),
                cname.c_str()),
            Form(
              "%.0f #leq #it{p}_{T,true} < %.0f GeV/#it{c};"
              "(#it{p}_{T,corr}-#it{p}_{T,true})/#it{p}_{T,true};"
              "weighted counts",
              bin_common_edges[iTruthBin],
              bin_common_edges[iTruthBin+1]
            ),
            80, -1.0, 1.0
          );
        }
      }

      // =====================================================
      // Pass 2: apply the mask event by event and fill outputs.
      // =====================================================
      for (Long64_t ie = 0; ie < nentries; ++ie) {
        tree->GetEntry(ie);

        const int ip = FindPtHatBin((double)xsecWeight);
        if (ip < kFirstPtHatBinToUse || ip < 0) continue;

        const bool haveMC = (mc_pt > 0.0f);
        const bool haveReco = (reco_pt > (Float_t)RECO_DUMMY_CUT);
        if (!haveMC || !haveReco) continue;
        if (!(deltaR > 0.0f && deltaR < (Float_t)dRmax)) continue;
        if (!(reco_area >= (Float_t)areaMin)) continue;
        if (!(reco_neutral_fraction <= (Float_t)CUT_NEUTRAL_FRACTION)) continue;
        if (reco_trigger_match != kTRUE) continue;

        const int recoMaskBin =
          FindVariableBin((double)reco_pt_corr, nbins_meas, bin_meas_edges);
        if (recoMaskBin < 0) continue;
        if (!keepRecoBin[ip][recoMaskBin]) continue;

        const double w =
          kXsecWeights[ip]/kNgenEvents[ip] * (double)centralityWeight;

        const int truthBin =
          FindVariableBin((double)mc_pt, nbins_common, bin_common_edges);

        for (int it = 0; it < N_LEAD; ++it) {
          if (reco_pt_lead < (Float_t)PTLEAD_THR[it]) continue;

          hMC[it]->Fill((double)mc_pt, w);
          hReco[it]->Fill((double)reco_pt_corr, w);

          if (truthBin >= 0 && mc_pt > 0.0f) {
            const double residual =
              ((double)reco_pt_corr - (double)mc_pt)/(double)mc_pt;

            moments[it][truthBin].Fill(residual, w);
            hResidual[it][truthBin]->Fill(residual, w);
          }
        }
      }

      // Convert spectra to densities.
      for (int it = 0; it < N_LEAD; ++it) {
        for (int iDensityBin = 1;
            iDensityBin <= nbins_common;
            ++iDensityBin) {

          const double bw = hMC[it]->GetBinWidth(iDensityBin);
          if (bw <= 0.0) continue;

          hMC[it]->SetBinContent(
            iDensityBin,
            hMC[it]->GetBinContent(iDensityBin)/bw
          );
          hMC[it]->SetBinError(
            iDensityBin,
            hMC[it]->GetBinError(iDensityBin)/bw
          );

          hReco[it]->SetBinContent(
            iDensityBin,
            hReco[it]->GetBinContent(iDensityBin)/bw
          );
          hReco[it]->SetBinError(
            iDensityBin,
            hReco[it]->GetBinError(iDensityBin)/bw
          );
        }
      }

      outR->cd();
      TDirectory* outC = outR->mkdir(cname.c_str());
      if (!outC) outC = outR->GetDirectory(cname.c_str());
      outC->cd();

      TString plotCDir = Form("%s/%s", plotRDir.Data(), cname.c_str());
      gSystem->mkdir(plotCDir.Data(), kTRUE);

      const TString centLabel = CentralityLabel(cname);

      for (int it = 0; it < N_LEAD; ++it) {
        // Build JES and JER graphs from the weighted residual moments.
        TGraphErrors* gJES = new TGraphErrors();
        TGraphErrors* gJER = new TGraphErrors();
        gJES->SetName(Form("gJES_ptlead%.0f_%s_%s",
                           PTLEAD_THR[it], rname.c_str(), cname.c_str()));
        gJER->SetName(Form("gJER_ptlead%.0f_%s_%s",
                           PTLEAD_THR[it], rname.c_str(), cname.c_str()));

        int np = 0;
        for (int iMomentBin = 0;
            iMomentBin < nbins_common;
            ++iMomentBin) {

          if (!moments[it][iMomentBin].Valid()) continue;
          if (moments[it][iMomentBin].Neff() <= 1.0) continue;

          const double x =
            0.5 * (
              bin_common_edges[iMomentBin] +
              bin_common_edges[iMomentBin + 1]
            );

          const double ex =
            0.5 * (
              bin_common_edges[iMomentBin + 1] -
              bin_common_edges[iMomentBin]
            );

          gJES->SetPoint(
            np,
            x,
            moments[it][iMomentBin].Mean()
          );
          gJES->SetPointError(
            np,
            ex,
            moments[it][iMomentBin].MeanError()
          );

          gJER->SetPoint(
            np,
            x,
            moments[it][iMomentBin].Sigma()
          );
          gJER->SetPointError(
            np,
            ex,
            moments[it][iMomentBin].SigmaError()
          );

          ++np;
        }

        hMC[it]->Write();
        hReco[it]->Write();

        for (int iTruthBin = 0;
            iTruthBin < nbins_common;
            ++iTruthBin) {
          hResidual[it][iTruthBin]->Write();
        }

        gJES->Write();
        gJER->Write();

        // -------------------------------------------------
        // MC versus reco spectrum: upper pad + reco/MC ratio
        // -------------------------------------------------
        TCanvas* cSpec = new TCanvas(
          Form("cMcReco_ptlead%.0f_%s_%s",
               PTLEAD_THR[it], rname.c_str(), cname.c_str()),
          "", 850, 850
        );
        TPad* pTop = new TPad("pTop","",0.0,0.30,1.0,1.0);
        TPad* pBot = new TPad("pBot","",0.0,0.00,1.0,0.30);
        pTop->SetBottomMargin(0.02);
        pTop->SetLogy();
        pBot->SetTopMargin(0.03);
        pBot->SetBottomMargin(0.32);
        pTop->Draw();
        pBot->Draw();

        StyleHistogram(hReco[it], 20, kBlue+1);
        StyleHistogram(hMC[it],   21, kGreen+2);

        pTop->cd();
        hReco[it]->SetTitle("");
        hReco[it]->GetXaxis()->SetLabelSize(0);
        hReco[it]->GetXaxis()->SetTitleSize(0);
        hReco[it]->GetYaxis()->SetTitle("weighted counts / GeV");

        double ymax = hReco[it]->GetMaximum();
        if (hMC[it]->GetMaximum() > ymax) {
          ymax = hMC[it]->GetMaximum();
        }
        if (ymax > 0.0) {
          hReco[it]->SetMaximum(5.0*ymax);
          double ymin = 1e-5*ymax;
            if (ymin < 1e-14) ymin = 1e-14;
            hReco[it]->SetMinimum(ymin);
        }

        hReco[it]->Draw("E1");
        hMC[it]->Draw("E1 SAME");

        TLegend* leg = new TLegend(0.16,0.16,0.43,0.29);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->AddEntry(hReco[it], "Reco #it{p}_{T}^{corr}", "lep");
        leg->AddEntry(hMC[it], "MC #it{p}_{T}^{true}", "lep");
        leg->Draw();

        TLatex tex;
        tex.SetNDC();
        tex.SetTextFont(42);
        tex.SetTextSize(0.040);
        tex.SetTextAlign(31);
        tex.DrawLatex(0.92,0.89,"STAR embedding Au+Au  #sqrt{#it{s}_{NN}} = 200 GeV");
        tex.DrawLatex(0.92,0.83,Form("R = %.1f, %s",R,centLabel.Data()));
        tex.DrawLatex(0.92,0.77,
                      Form("#it{p}_{T}^{lead} #geq %.0f GeV/#it{c}",
                           PTLEAD_THR[it]));
        tex.DrawLatex(0.92,0.71,"THIS THESIS");

        pBot->cd();
        TH1D* hRatio = (TH1D*)hReco[it]->Clone(
          Form("hRecoOverMc_ptlead%.0f_%s_%s",
               PTLEAD_THR[it],rname.c_str(),cname.c_str())
        );
        hRatio->Divide(hMC[it]);
        hRatio->SetTitle("");
        hRatio->GetXaxis()->SetTitle("#it{p}_{T,jet} [GeV/#it{c}]");
        hRatio->GetYaxis()->SetTitle("Reco / MC");
        hRatio->GetYaxis()->SetRangeUser(0.0,2.0);
        hRatio->GetXaxis()->SetTitleSize(0.11);
        hRatio->GetXaxis()->SetLabelSize(0.09);
        hRatio->GetYaxis()->SetTitleSize(0.09);
        hRatio->GetYaxis()->SetLabelSize(0.08);
        hRatio->GetYaxis()->SetTitleOffset(0.55);
        hRatio->GetYaxis()->SetNdivisions(505);
        hRatio->Draw("E1");

        TLine* unity = new TLine(0.0,1.0,60.0,1.0);
        unity->SetLineStyle(2);
        unity->Draw("SAME");

        TString specBase = Form(
          "%s/mc_reco_ptlead%.0f_%s_%s",
          plotCDir.Data(), PTLEAD_THR[it], rname.c_str(), cname.c_str()
        );
        cSpec->SaveAs(specBase+".png");

        hRatio->Write();

        delete unity;
        delete hRatio;
        delete leg;
        delete pTop;
        delete pBot;
        delete cSpec;

        // --------------------------------------------
        // JES and JER summary; no residual-fit panels.
        // --------------------------------------------
        TCanvas* cJJ = new TCanvas(
          Form("cJesJer_ptlead%.0f_%s_%s",
               PTLEAD_THR[it], rname.c_str(), cname.c_str()),
          "", 1300, 600
        );
        cJJ->Divide(2,1);

        cJJ->cd(1);
        gPad->SetLeftMargin(0.13);
        gPad->SetBottomMargin(0.13);
        gJES->SetTitle("");
        gJES->SetMarkerStyle(20);
        gJES->SetMarkerColor(kBlue+1);
        gJES->SetLineColor(kBlue+1);
        gJES->GetXaxis()->SetTitle("#it{p}_{T}^{true} [GeV/#it{c}]");
        gJES->GetYaxis()->SetTitle("JES = #LT(#it{p}_{T}^{corr}-#it{p}_{T}^{true})/#it{p}_{T}^{true}#GT");
        gJES->GetXaxis()->SetLimits(0.0,60.0);
        gJES->SetMinimum(-0.6);
        gJES->SetMaximum(0.3);
        gJES->Draw("AP");

        TLine* zero = new TLine(0.0,0.0,60.0,0.0);
        zero->SetLineStyle(2);
        zero->SetLineColor(kRed+1);
        zero->Draw("SAME");

        cJJ->cd(2);
        gPad->SetLeftMargin(0.13);
        gPad->SetBottomMargin(0.13);
        gJER->SetTitle("");
        gJER->SetMarkerStyle(20);
        gJER->SetMarkerColor(kBlue+1);
        gJER->SetLineColor(kBlue+1);
        gJER->GetXaxis()->SetTitle("#it{p}_{T}^{true} [GeV/#it{c}]");
        gJER->GetYaxis()->SetTitle("JER = #sigma[(#it{p}_{T}^{corr}-#it{p}_{T}^{true})/#it{p}_{T}^{true}]");
        gJER->GetXaxis()->SetLimits(0.0,60.0);
        gJER->SetMinimum(0.0);
        gJER->SetMaximum(0.8);
        gJER->Draw("AP");

        cJJ->cd(0);
        TLatex head;
        head.SetNDC();
        head.SetTextFont(42);
        head.SetTextAlign(22);
        head.SetTextSize(0.030);
        head.DrawLatex(0.50,0.965,
          Form("STAR embedding Au+Au, #sqrt{#it{s}_{NN}} = 200 GeV, R = %.1f, %s, #it{p}_{T}^{lead} #geq %.0f GeV/#it{c}",
               R,centLabel.Data(),PTLEAD_THR[it]));

        TString jjBase = Form(
          "%s/jes_jer_ptlead%.0f_%s_%s",
          plotCDir.Data(), PTLEAD_THR[it], rname.c_str(), cname.c_str()
        );
        cJJ->SaveAs(jjBase+".png");

        delete zero;
        delete cJJ;
        delete gJES;
        delete gJER;
      }

      for (int it = 0; it < N_LEAD; ++it) {
        delete hMC[it];
        delete hReco[it];

        for (int iTruthBin = 0;
            iTruthBin < nbins_common;
            ++iTruthBin) {
          delete hResidual[it][iTruthBin];
        }
      }
      for (int ip = 0; ip < kNPthatBins; ++ip) {
        delete hMaskSource[ip];
      }
    }
  }

  fout->Write();
  fout->Close();
  fin->Close();

  std::cout << "Done. ROOT output: " << outfile << std::endl;
  std::cout << "Plots written under: " << plotTop << std::endl;
}