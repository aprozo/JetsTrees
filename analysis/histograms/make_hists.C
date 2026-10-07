#include "TFile.h"
#include "TDirectory.h"
#include "TDirectoryFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
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

#include <iostream>
#include <string>
#include <cmath>

// =========================================================
//  Jet-quality cuts (must match maker!)
// =========================================================
const double CUT_AREA_02 = 0.07;  // R = 0.2
const double CUT_AREA_03 = 0.20;  // R = 0.3
const double CUT_AREA_04 = 0.40;  // R = 0.4
const double CUT_NEUTRAL_FRACTION = 0.95;
const double RECO_DUMMY_CUT = -500.0;  // dummy threshold

static const int kFirstPtHatBinToUse = 2;  // if == 2; ignore pThat bins 0 and 1  

// pTlead thresholds
const int    N_LEAD = 4;
const double PTLEAD_THR[N_LEAD] = {0.0, 5.0, 7.0, 9.0};

// helper: choose area cut based on R
double areaMinForR(double R) {
  if (R < 0.25) return CUT_AREA_02;
  if (R < 0.35) return CUT_AREA_03;
  return CUT_AREA_04;
}

static const int kNPthatBins = 11;
static const double kXsecWeights[kNPthatBins] = {
  1.616e+0, 1.355e-01, 2.288e-02, 5.524e-03, 2.203e-03,
  3.437e-04, 4.681e-05, 8.532e-06, 2.178e-06, 1.198e-07, 6.939e-09
};

static const double kNgenEvents[kNPthatBins] = {
  1020062, 1529646, 1275275, 1019532, 1019730,
  1020088, 1019739, 765165, 509510, 305922, 101971
};

static const double kMinSignif = std::sqrt(10.0);

// measured & truth binning used in unfolding
static const int nbins_meas = 24;
static const double bin_meas_edges[nbins_meas+1] = {
  -100,-80,-60,-40,-20,-10,-5,-2.5,0,2.5,5,7.5,10,12.5,15,17.5,
  20,22.5,25,27.5,30,35,40,50,60
};

static const int nbins_truth = 10;
static const double bin_truth_edges[nbins_truth+1] = {
  0,5,10,15,20,25,30,35,40,50,60
};

static int FindPtHatBin(double xsecW)
{
  const double relTol = 1e-6;
  for (int i = 0; i < kNPthatBins; ++i) {
    double ref = kXsecWeights[i];
    double diff = std::fabs(xsecW - ref);
    if (diff <= relTol * std::fabs(ref)) return i;
  }
  return -1;
}

void make_hists(const char *infile  = "embedding_merged.root",
                const char *outfile = "hists.root")
{
  TH1::SetDefaultSumw2(kTRUE);
  gStyle->SetOptStat(0);
  gROOT->SetBatch(kTRUE);

  // where to write PDFs
  const char* pdfTop = "pdf";
  gSystem->mkdir(pdfTop, kTRUE);

  TFile *fin = TFile::Open(infile, "READ");
  if (!fin || fin->IsZombie()) {
    std::cerr << "Error: cannot open input file " << infile << std::endl;
    return;
  }

  TFile *fout = TFile::Open(outfile, "RECREATE");
  if (!fout || fout->IsZombie()) {
    std::cerr << "Error: cannot create output file " << outfile << std::endl;
    fin->Close();
    return;
  }

  std::cout << "Input : " << infile  << std::endl;
  std::cout << "Output: " << outfile << std::endl;

  // Loop over R directories: R0.2, R0.3, R0.4, ...
  TIter nextR(fin->GetListOfKeys());
  TKey *keyR = 0;

  while ((keyR = (TKey*) nextR())) {
    if (std::string(keyR->GetClassName()) != "TDirectoryFile") continue;

    std::string rname = keyR->GetName();
    if (rname.empty() || rname[0] != 'R') continue;

    double R = 0.0;
    if (sscanf(rname.c_str(), "R%lf", &R) != 1) continue;

    const double dRmax = 0.6 * R;
    const double areaMin = areaMinForR(R);

    TDirectoryFile *dirR = (TDirectoryFile*) keyR->ReadObj();
    if (!dirR) continue;

    std::cout << "Radius directory: " << rname
              << " (R=" << R << ", dRmax=" << dRmax
              << ", areaMin=" << areaMin << ")\n";

    // output dir R
    fout->cd();
    TDirectory *outR = fout->mkdir(rname.c_str());
    if (!outR) outR = fout->GetDirectory(rname.c_str());

    // pdf dir for this R
    TString pdfRDir; pdfRDir.Form("%s/%s", pdfTop, rname.c_str());
    gSystem->mkdir(pdfRDir.Data(), kTRUE);

    // Loop over centrality subdirs
    TIter nextC(dirR->GetListOfKeys());
    TKey *keyC = 0;

    while ((keyC = (TKey*) nextC())) {
      if (std::string(keyC->GetClassName()) != "TDirectoryFile") continue;

      std::string cname = keyC->GetName();
      TDirectoryFile *dirC = (TDirectoryFile*) keyC->ReadObj();
      if (!dirC) continue;

      TTree *tree = (TTree*) dirC->Get("JetTree");
      if (!tree) continue;

      Long64_t nentries = tree->GetEntries();
      std::cout << "  Centrality: " << cname << " (entries: " << nentries << ")\n";

      // require embedding-like branches
      if (!tree->GetBranch("mc_pt") ||
          !tree->GetBranch("xsecWeight") ||
          !tree->GetBranch("deltaR") ||
          !tree->GetBranch("centralityWeight") ||
          !tree->GetBranch("reco_pt") ||
          !tree->GetBranch("reco_pt_corr") ||
          !tree->GetBranch("reco_area") ||
          !tree->GetBranch("reco_neutral_fraction") ||
          !tree->GetBranch("reco_trigger_match") ||
          !tree->GetBranch("reco_pt_lead")) {
        std::cerr << "    Missing branches in " << rname << "/" << cname << " -> skipping.\n";
        continue;
      }

      // output dirs
      outR->cd();
      TDirectory *outC = outR->mkdir(cname.c_str());
      if (!outC) outC = outR->GetDirectory(cname.c_str());
      outC->cd();

      // pdf dir for this centrality
      TString pdfCDir; pdfCDir.Form("%s/%s", pdfRDir.Data(), cname.c_str());
      gSystem->mkdir(pdfCDir.Data(), kTRUE);

      // ---------------- branches ----------------
      Float_t xsecWeight = 0.f;
      Float_t centralityWeight = 1.f;

      Float_t mc_pt = -999.f;
      Float_t deltaR = -1.f;

      Float_t reco_pt = -999.f;
      Float_t reco_pt_corr = -999.f;
      Float_t reco_area = -999.f;
      Float_t reco_neutral_fraction = 999.f;
      Float_t reco_pt_lead = -999.f;
      Bool_t  reco_trigger_match = kFALSE;

      tree->SetBranchAddress("xsecWeight",        &xsecWeight);
      tree->SetBranchAddress("centralityWeight",  &centralityWeight);
      tree->SetBranchAddress("mc_pt",             &mc_pt);
      tree->SetBranchAddress("deltaR",            &deltaR);
      tree->SetBranchAddress("reco_pt",           &reco_pt);
      tree->SetBranchAddress("reco_pt_corr",      &reco_pt_corr);
      tree->SetBranchAddress("reco_area",         &reco_area);
      tree->SetBranchAddress("reco_neutral_fraction", &reco_neutral_fraction);
      tree->SetBranchAddress("reco_pt_lead",      &reco_pt_lead);
      tree->SetBranchAddress("reco_trigger_match",&reco_trigger_match);

      // ---------------- histogram booking ----------------
      const int nbins_corr = 100;
      const double xmin_corr = -40.0;
      const double xmax_corr =  60.0;

      const int nbins_pt = 60;
      const double xmin_pt = 0.0;
      const double xmax_pt = 60.0;

      // 2D: reco_pt_corr vs mc_pt (for the selected sample; inclusive ptlead>=0)
    TH2D* h2_recoCorr_vs_mc =
  new TH2D(Form("h2_recoPtCorr_vs_mcPt_%s_%s", rname.c_str(), cname.c_str()),
           "Reco p_{T}^{corr} vs MC p_{T};p_{T}^{reco,corr} [GeV];p_{T}^{MC} [GeV]",
           nbins_meas, bin_meas_edges,
           nbins_truth, bin_truth_edges);

      // spectra per ptlead threshold (reco_pt_corr)
      TH1D* hSpec[N_LEAD];
      for (int it = 0; it < N_LEAD; ++it) {
        hSpec[it] = new TH1D(Form("hSpec_recoPtCorr_ptlead%.0f_%s_%s", PTLEAD_THR[it], rname.c_str(), cname.c_str()),
                            Form("Reco p_{T}^{corr} spectrum (ptlead>=%.0f);p_{T}^{reco,corr} [GeV];weighted counts",
                                  PTLEAD_THR[it]),
                            nbins_meas, bin_meas_edges);
      }

      // Fine/raw spectra before analysis rebinning and insignificant-bin removal
      TH1D* hSpecFine[N_LEAD];
      for (int it = 0; it < N_LEAD; ++it) {
        hSpecFine[it] = new TH1D(
          Form("hSpecFine_recoPtCorr_ptlead%.0f_%s_%s", PTLEAD_THR[it], rname.c_str(), cname.c_str()),
          Form("Fine reco p_{T}^{corr} spectrum before masking (ptlead>=%.0f);p_{T}^{reco,corr} [GeV];weighted counts",
              PTLEAD_THR[it]),
          nbins_corr, xmin_corr, xmax_corr
        );
      }

      TH2D* h2_recoCorr_vs_mc_pthat[kNPthatBins];
      TH1D* hSpec_pthat[N_LEAD][kNPthatBins];

      for (int ip = 0; ip < kNPthatBins; ++ip) {
        h2_recoCorr_vs_mc_pthat[ip] =
          new TH2D(Form("h2_recoPtCorr_vs_mcPt_pthat%d_%s_%s", ip, rname.c_str(), cname.c_str()),
                  "Reco p_{T}^{corr} vs MC p_{T};p_{T}^{reco,corr} [GeV];p_{T}^{MC} [GeV]",
                  nbins_meas, bin_meas_edges,
                  nbins_truth, bin_truth_edges);

        for (int it = 0; it < N_LEAD; ++it) {
          hSpec_pthat[it][ip] =
            new TH1D(Form("hSpec_recoPtCorr_ptlead%.0f_pthat%d_%s_%s",
                          PTLEAD_THR[it], ip, rname.c_str(), cname.c_str()),
                    Form("Reco p_{T}^{corr} spectrum (ptlead>=%.0f);p_{T}^{reco,corr} [GeV];centrality-weighted counts",
                          PTLEAD_THR[it]),
                    nbins_meas, bin_meas_edges);
        }
      }

      // ---------------- loop entries ----------------
      for (Long64_t i = 0; i < nentries; ++i) {
        tree->GetEntry(i);

        // base weights
        const int ip = FindPtHatBin((double)xsecWeight);
        if (ip < 0) continue;

        // skip first two pThat bins
        if (ip < kFirstPtHatBinToUse) continue;

        double wCent = (double)centralityWeight;

        // matched condition (embedding matched tree entry)
        bool haveMC   = (mc_pt > 0.0f);
        bool haveReco = (reco_pt > (Float_t)RECO_DUMMY_CUT);

        if (!haveMC) continue;
        if (!haveReco) continue;
        if (!(deltaR > 0.0f && deltaR < (Float_t)dRmax)) continue;

        // reco quality cuts
        if (!(reco_area >= (Float_t)areaMin)) continue;
        if (!(reco_neutral_fraction <= (Float_t)CUT_NEUTRAL_FRACTION)) continue;

        // trigger match
        if (!(reco_trigger_match == kTRUE)) continue;


        // ---------------------------------------------------------
        // Fine/raw spectra BEFORE analysis rebinning and bin removal
        // Here we apply xsecWeight immediately because we are not
        // merging per-pThat bins later for these histograms.
        // ---------------------------------------------------------
        const double wFull = (double)xsecWeight / kNgenEvents[ip] * (double)centralityWeight;

        for (int it = 0; it < N_LEAD; ++it) {
          if (reco_pt_lead >= (Float_t)PTLEAD_THR[it]) {
            hSpecFine[it]->Fill((double)reco_pt_corr, wFull);
          }
        }


        // fill 2D for inclusive (ptlead>=0)
        h2_recoCorr_vs_mc_pthat[ip]->Fill((double)reco_pt_corr, (double)mc_pt, wCent);

        for (int it = 0; it < N_LEAD; ++it) {
          if (reco_pt_lead >= (Float_t)PTLEAD_THR[it]) {
            hSpec_pthat[it][ip]->Fill((double)reco_pt_corr, wCent);
          }
        }
      }

      h2_recoCorr_vs_mc->Reset();
      for (int it = 0; it < N_LEAD; ++it) {
        hSpec[it]->Reset();
      }

      for (int ip = 0; ip < kNPthatBins; ++ip) {
        const double xw = kXsecWeights[ip] / kNgenEvents[ip];

        for (int ix = 1; ix <= nbins_meas; ++ix) {
          // build mask from inclusive reco spectrum (ptlead >= 0)
          const double binC = hSpec_pthat[0][ip]->GetBinContent(ix);
          const double binE = hSpec_pthat[0][ip]->GetBinError(ix);
          const double signif = (binE > 0.0 ? binC / binE : 0.0);
          const bool keep = (signif > kMinSignif);

          if (!keep) continue;

          // ---- merge 1D spectra ----
          for (int it = 0; it < N_LEAD; ++it) {
            const double oldC = hSpec[it]->GetBinContent(ix);
            const double oldE = hSpec[it]->GetBinError(ix);
            const double addC = xw * hSpec_pthat[it][ip]->GetBinContent(ix);
            const double addE = xw * hSpec_pthat[it][ip]->GetBinError(ix);

            hSpec[it]->SetBinContent(ix, oldC + addC);
            hSpec[it]->SetBinError(ix, std::sqrt(oldE*oldE + addE*addE));
          }

          // ---- merge 2D histogram using same reco-bin mask ----
          for (int iy = 1; iy <= nbins_truth; ++iy) {
            const double oldC = h2_recoCorr_vs_mc->GetBinContent(ix, iy);
            const double oldE = h2_recoCorr_vs_mc->GetBinError(ix, iy);
            const double addC = xw * h2_recoCorr_vs_mc_pthat[ip]->GetBinContent(ix, iy);
            const double addE = xw * h2_recoCorr_vs_mc_pthat[ip]->GetBinError(ix, iy);

            h2_recoCorr_vs_mc->SetBinContent(ix, iy, oldC + addC);
            h2_recoCorr_vs_mc->SetBinError(ix, iy, std::sqrt(oldE*oldE + addE*addE));
          }
        }
      }

      // Convert rebinned spectra from weighted counts per bin
      // to weighted counts per GeV.
      for (int it = 0; it < N_LEAD; ++it) {
        for (int ix = 1; ix <= hSpec[it]->GetNbinsX(); ++ix) {
          const double bw = hSpec[it]->GetBinWidth(ix);
          if (bw <= 0.0) continue;

          hSpec[it]->SetBinContent(ix, hSpec[it]->GetBinContent(ix) / bw);
          hSpec[it]->SetBinError(ix,   hSpec[it]->GetBinError(ix)   / bw);
        }

        hSpec[it]->GetYaxis()->SetTitle("Weighted counts / GeV");
      }

      // ---------------- write hists ----------------
      outC->cd();
      h2_recoCorr_vs_mc->Write();

      for (int it = 0; it < N_LEAD; ++it) {
        hSpec[it]->Write();
        hSpecFine[it]->Write();
      }

      double Rval = 0.0;
      sscanf(rname.c_str(), "R%lf", &Rval);

      TString centLabel = cname.c_str();
      centLabel.ReplaceAll("CENT_", "");
      centLabel.ReplaceAll("MID_",  "");
      centLabel.ReplaceAll("PERI_", "");
      centLabel.ReplaceAll("_", "-");
      centLabel += " %";

      // ---------------- make PDF: overlay of all ptlead + ratios to pTlead 0 ----------------
      {
        TCanvas* cSpec = new TCanvas("cSpec","cSpec",800,800);

        TPad* padTop = new TPad("padTop","padTop",0.0,0.30,1.0,1.0);
        TPad* padBot = new TPad("padBot","padBot",0.0,0.00,1.0,0.30);

        padTop->SetBottomMargin(0.02);
        padTop->SetLogy();

        padBot->SetTopMargin(0.03);
        padBot->SetBottomMargin(0.32);

        padTop->Draw();
        padBot->Draw();

        // styles
        hSpec[0]->SetMarkerStyle(20);
        hSpec[0]->SetLineColor(kBlack);
        hSpec[0]->SetMarkerColor(kBlack);

        hSpec[1]->SetMarkerStyle(21);
        hSpec[1]->SetLineColor(kRed+1);
        hSpec[1]->SetMarkerColor(kRed+1);

        hSpec[2]->SetMarkerStyle(22);
        hSpec[2]->SetLineColor(kBlue+1);
        hSpec[2]->SetMarkerColor(kBlue+1);

        hSpec[3]->SetMarkerStyle(23);
        hSpec[3]->SetLineColor(kGreen+2);
        hSpec[3]->SetMarkerColor(kGreen+2);

        // upper pad: spectra
        padTop->cd();

        hSpec[0]->SetTitle("");
        hSpec[0]->GetXaxis()->SetLabelSize(0);
        hSpec[0]->GetXaxis()->SetTitleSize(0);
        hSpec[0]->GetXaxis()->SetRangeUser(-10, 60);
        hSpec[0]->GetYaxis()->SetTitle("Weighted counts / GeV");
        hSpec[0]->Draw("E1");

        hSpec[1]->Draw("E1 SAME");
        hSpec[2]->Draw("E1 SAME");
        hSpec[3]->Draw("E1 SAME");

        TLegend* leg = new TLegend(0.60, 0.68, 0.88, 0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.04);

        for (int it=0; it<N_LEAD; ++it) {
          leg->AddEntry(hSpec[it],
                        Form("#it{p}_{T}^{lead} #geq %.0f GeV", PTLEAD_THR[it]),
                        "lep");
        }
        leg->Draw();

        TString line1;
        line1.Form("R = %.1f, %s", Rval, centLabel.Data());

        TLatex tex;
        tex.SetNDC();
        tex.SetTextFont(42);
        tex.SetTextSize(0.045);

        tex.DrawLatex(0.32, 0.25, line1.Data());
        tex.DrawLatex(0.32, 0.19, "Au+Au  #sqrt{#it{s}_{NN}} = 200 GeV");
        tex.DrawLatex(0.32, 0.13, "THIS THESIS");

        // lower pad: ratios to pTlead >= 0
        padBot->cd();

        TH1D* hRatio[N_LEAD];
        for (int it = 0; it < N_LEAD; ++it) {
          hRatio[it] = (TH1D*)hSpec[it]->Clone(
            Form("hRatio_ptlead%.0f_over_ptlead0_%s_%s",
                 PTLEAD_THR[it], rname.c_str(), cname.c_str())
          );
          hRatio[it]->Divide(hSpec[0]);
        }

        hRatio[1]->SetTitle("");
        hRatio[1]->GetXaxis()->SetTitle("#it{p}_{T, reco}^{corr} (GeV/#it{c})");
        hRatio[1]->GetYaxis()->SetTitle("ratio to #it{p}_{T}^{lead} #geq 0");
        hRatio[1]->GetXaxis()->SetRangeUser(-10, 60);

        hRatio[1]->GetYaxis()->SetRangeUser(0.0, 1.25);

        hRatio[1]->GetXaxis()->SetTitleSize(0.11);
        hRatio[1]->GetXaxis()->SetLabelSize(0.09);
        hRatio[1]->GetYaxis()->SetTitleSize(0.09);
        hRatio[1]->GetYaxis()->SetLabelSize(0.08);
        hRatio[1]->GetYaxis()->SetTitleOffset(0.55);
        hRatio[1]->GetYaxis()->SetNdivisions(505);

        hRatio[1]->Draw("E1");
        hRatio[2]->Draw("E1 SAME");
        hRatio[3]->Draw("E1 SAME");

        TLine* line = new TLine(-10.0, 1.0, 60.0, 1.0);
        line->SetLineStyle(2);
        line->Draw("SAME");

        TString pdfName;
        pdfName.Form("%s/spectra_ptlead_embed_%s_%s.pdf",
                     pdfCDir.Data(), rname.c_str(), cname.c_str());

        cSpec->SaveAs(pdfName.Data());

        TString pngName = pdfName;
        pngName.ReplaceAll(".pdf", ".png");
        cSpec->SaveAs(pngName.Data());

        delete line;
        for (int it = 0; it < N_LEAD; ++it) delete hRatio[it];
        delete leg;
        delete padTop;
        delete padBot;
        delete cSpec;
      }

           // ---------------- make PDF: fine spectra before masking ----------------
      {
        TCanvas* cSpecFine = new TCanvas("cSpecFine","cSpecFine",800,800);

        TPad* padTop = new TPad("padTopFine","padTopFine",0.0,0.30,1.0,1.0);
        TPad* padBot = new TPad("padBotFine","padBotFine",0.0,0.00,1.0,0.30);

        padTop->SetBottomMargin(0.02);
        padTop->SetLogy();

        padBot->SetTopMargin(0.03);
        padBot->SetBottomMargin(0.32);

        padTop->Draw();
        padBot->Draw();

        // ==========================
        // Upper pad: spectra
        // ==========================
        padTop->cd();

        hSpecFine[0]->SetTitle("");
        hSpecFine[0]->GetXaxis()->SetLabelSize(0);
        hSpecFine[0]->GetXaxis()->SetTitleSize(0);
        hSpecFine[0]->GetXaxis()->SetRangeUser(-10,60);

        hSpecFine[0]->GetYaxis()->SetTitle("Weighted counts");

        hSpecFine[0]->SetMarkerStyle(20);
        hSpecFine[0]->SetLineColor(kBlack);
        hSpecFine[0]->SetMarkerColor(kBlack);

        hSpecFine[1]->SetMarkerStyle(21);
        hSpecFine[1]->SetLineColor(kRed+1);
        hSpecFine[1]->SetMarkerColor(kRed+1);

        hSpecFine[2]->SetMarkerStyle(22);
        hSpecFine[2]->SetLineColor(kBlue+1);
        hSpecFine[2]->SetMarkerColor(kBlue+1);

        hSpecFine[3]->SetMarkerStyle(23);
        hSpecFine[3]->SetLineColor(kGreen+2);
        hSpecFine[3]->SetMarkerColor(kGreen+2);

        hSpecFine[0]->Draw("E1");
        hSpecFine[1]->Draw("E1 SAME");
        hSpecFine[2]->Draw("E1 SAME");
        hSpecFine[3]->Draw("E1 SAME");

        TLegend* legFine = new TLegend(0.60,0.68,0.88,0.88);
        legFine->SetBorderSize(0);
        legFine->SetFillStyle(0);
        legFine->SetTextSize(0.04);

        for (int it=0; it<N_LEAD; ++it) {
          legFine->AddEntry(
            hSpecFine[it],
            Form("#it{p}_{T}^{lead} #geq %.0f GeV",PTLEAD_THR[it]),
            "lep"
          );
        }

        legFine->Draw();

        TString line1;
        line1.Form("R = %.1f, %s", Rval, centLabel.Data());

        TLatex tex;
        tex.SetNDC();
        tex.SetTextFont(42);
        tex.SetTextSize(0.045);

        tex.DrawLatex(0.32,0.25,line1.Data());
        tex.DrawLatex(0.32,0.19,"Au+Au  #sqrt{#it{s}_{NN}} = 200 GeV");
        tex.DrawLatex(0.32,0.13,"Before analysis rebinning / bin removal");

        // ==========================
        // Lower pad: ratios
        // ==========================
        padBot->cd();

        TH1D* hRatioFine[N_LEAD];

        for (int it = 0; it < N_LEAD; ++it) {
          hRatioFine[it] =
            (TH1D*)hSpecFine[it]->Clone(
              Form("hRatioFine_ptlead%.0f_over_ptlead0_%s_%s",
                   PTLEAD_THR[it],
                   rname.c_str(),
                   cname.c_str())
            );

          hRatioFine[it]->Divide(hSpecFine[0]);
        }

        hRatioFine[1]->SetTitle("");

        hRatioFine[1]->GetXaxis()->SetTitle("#it{p}_{T,reco}^{corr} (GeV/#it{c})");
        hRatioFine[1]->GetYaxis()->SetTitle("ratio to p_{T}^{lead}>0");

        hRatioFine[1]->GetXaxis()->SetRangeUser(-10,60);
        hRatioFine[1]->GetYaxis()->SetRangeUser(0.0,1.25);

        hRatioFine[1]->GetXaxis()->SetTitleSize(0.11);
        hRatioFine[1]->GetXaxis()->SetLabelSize(0.09);

        hRatioFine[1]->GetYaxis()->SetTitleSize(0.09);
        hRatioFine[1]->GetYaxis()->SetLabelSize(0.08);
        hRatioFine[1]->GetYaxis()->SetTitleOffset(0.55);
        hRatioFine[1]->GetYaxis()->SetNdivisions(505);

        hRatioFine[1]->Draw("E1");
        hRatioFine[2]->Draw("E1 SAME");
        hRatioFine[3]->Draw("E1 SAME");

        TLine* lineFine = new TLine(-10.0,1.0,60.0,1.0);
        lineFine->SetLineStyle(2);
        lineFine->Draw("SAME");

        TString pdfNameFine;
        pdfNameFine.Form("%s/spectra_ptlead_embed_FINE_beforeMask_%s_%s.pdf",
                         pdfCDir.Data(),
                         rname.c_str(),
                         cname.c_str());

        cSpecFine->SaveAs(pdfNameFine.Data());

        TString pngNameFine = pdfNameFine;
        pngNameFine.ReplaceAll(".pdf",".png");
        cSpecFine->SaveAs(pngNameFine.Data());

        delete lineFine;

        for (int it = 0; it < N_LEAD; ++it)
          delete hRatioFine[it];

        delete legFine;
        delete padTop;
        delete padBot;
        delete cSpecFine;
      }

      delete h2_recoCorr_vs_mc;

      for (int it = 0; it < N_LEAD; ++it) {
        delete hSpec[it];
        delete hSpecFine[it];
      }

      for (int ip = 0; ip < kNPthatBins; ++ip) {
        delete h2_recoCorr_vs_mc_pthat[ip];
        for (int it = 0; it < N_LEAD; ++it) {
          delete hSpec_pthat[it][ip];
        }
      }
    }
  }

  fout->Write();
  fout->Close();
  fin->Close();

  std::cout << "Done. ROOT output: " << outfile << "\n";
  std::cout << "PDFs written under: " << pdfTop << "/R*/CENT_*/\n";
}
