#include "RooUnfoldBayes.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldSvd.h"
#include "TSVDUnfold.h"

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TSystem.h"

#include <vector>
#include <string>
#include <iostream>
#include <cmath>

using std::string;
using std::vector;
using std::cout;
using std::endl;

// ------------------------- config -----------------------------
static const double kTestFrac = 0.50;  // 50/50 split
static const UInt_t kSeed     = 12345;    // deterministic split

static const double kPtLeadCuts[] = {0.0, 5.0, 7.0, 9.0};
static const int    kNPtLeadCuts  = sizeof(kPtLeadCuts)/sizeof(kPtLeadCuts[0]);

static const double kMinSignif = std::sqrt(10.0);  // content/error > sqrt(10)
static const bool   kSavePtHatDebug = true;

static const int kFirstPtHatBinToUse = 0;  // if == 2; ignore pThat bins 0 and 1  

// measured & truth binning
static const int nbins_meas = 24;
static const double bin_meas_edges[nbins_meas+1] = {
  -100,-80,-60,-40,-20,-10,-5,-2.5,0,2.5,5,7.5,10,12.5,15,17.5,
  20,22.5,25,27.5,30,35,40,50,60
};

static const int nbins_truth = 10;
static const double bin_truth_edges[nbins_truth+1] = {
  0,5,10,15,20,25,30,35,40,50,60
};


static const vector<string> kCentralities =
  {"CENT_0_10", "MID_20_40", "PERI_60_80"};

static TString NiceCentLabel(const std::string& centToken)
{
  int a = -1, b = -1;
  if (sscanf(centToken.c_str(), "CENT_%d_%d", &a, &b) == 2 ||
      sscanf(centToken.c_str(), "MID_%d_%d", &a, &b) == 2 ||
      sscanf(centToken.c_str(), "PERI_%d_%d", &a, &b) == 2) {
    return Form("%d#font[52]{#minus}%d %%", a, b);
  }
  return TString(centToken.c_str());
}

static const vector<string> kRadii =
  {"R0.2", "R0.3", "R0.4"};

// ---- pThat bins (upper edges) and xsec weights (same order) ----
static const int kNPthatBins = 11;

// upper edges for each bin: 3_5 -> 5, ..., 40_50 -> 50, 50_-1 -> -1 (no upper bound)
static const double kPtHatMax[kNPthatBins] =
  {5, 7, 9, 11, 15, 20, 25, 30, 40, 50, -1};

// xsecWeight values that tag each bin (from your code)
static const double kXsecWeights[kNPthatBins] =
  {1.616e+0,  1.355e-01, 2.288e-02, 5.524e-03, 2.203e-03,
   3.437e-04, 4.681e-05, 8.532e-06, 2.178e-06, 1.198e-07, 6.939e-09};

static const double kNgenEvents[kNPthatBins] =
  {1020062, 1529646, 1275275, 1019532, 1019730,
   1020088, 1019739, 765165, 509510, 305922, 101971};

// reco dummy sentinel (keep real negative jets, reject dummy ~ -999)
static const double RECO_PTCORR_DUMMY_CUT = -500.0;

// match xsecWeight -> bin index (tolerant compare)
static int FindPtHatBin(double xsecW)
{
  // relative tolerance, because floats are a joy
  const double relTol = 1e-6;
  for (int i = 0; i < kNPthatBins; ++i) {
    double ref = kXsecWeights[i];
    double diff = fabs(xsecW - ref);
    if (diff <= relTol * fabs(ref)) return i;
  }
  return -1; // unknown
}

// ------------------ jet-quality cuts -------------------------

// Area cuts per jet radius
static const double CUT_AREA_02 = 0.07;  // R = 0.2
static const double CUT_AREA_03 = 0.20;  // R = 0.3
static const double CUT_AREA_04 = 0.40;  // R = 0.4

// Max neutral energy fraction
static const double CUT_NEUTRAL_FRACTION = 0.95;

// --------------------------------------------------------------
static bool DrawBayesIter(int iter)
{
  return (iter == 1 || iter == 4 || iter == 5 || iter == 8);
}

static void EnsureDir(const string& path){
  if (gSystem->AccessPathName(path.c_str()))
    gSystem->mkdir(path.c_str(), /*recursive=*/true);
}

// Main macro: builds closure + full-statistics responses
void unfold_embedding_inclusive(const char* inputFile,
                      const char* outDir, const char* method = "BAYES")
{
  std::string m(method);
  std::cout << ">>> Unfolding method = " << m << std::endl;

  gStyle->SetOptStat(0);
  EnsureDir(outDir);
  TH1::SetDefaultSumw2(kTRUE);

  TFile* fin = TFile::Open(inputFile, "READ");
  if (!fin || fin->IsZombie()) {
    cout << "[error] Cannot open input file: " << inputFile << endl;
    return;
  }

  const string outRootPath = string(outDir) + "/responses_embedding.root";
  TFile* fout = TFile::Open(outRootPath.c_str(), "RECREATE");
  if (!fout || fout->IsZombie()) {
    cout << "[error] Cannot create output file: " << outRootPath << endl;
    fin->Close(); delete fin;
    return;
  }

  cout << "Input embedding : " << inputFile  << endl;
  cout << "Output responses: " << outRootPath << endl;

  // loop R, centrality, ptlead
  for (size_t iR = 0; iR < kRadii.size(); ++iR) {
    const string R = kRadii[iR];

    // parse numeric R to choose area cut
    double Rval = 0.0;
    if (sscanf(R.c_str(), "R%lf", &Rval) != 1) {
      cout << "[warn] Could not parse radius from " << R << ", skipping.\n";
      continue;
    }
    double areaMin = 0.0;
    if      (Rval < 0.25) areaMin = CUT_AREA_02; // ~0.2
    else if (Rval < 0.35) areaMin = CUT_AREA_03; // ~0.3
    else                  areaMin = CUT_AREA_04; // ~0.4

    for (size_t iC = 0; iC < kCentralities.size(); ++iC) {
      const string C = kCentralities[iC];

      const string treePath = R + "/" + C + "/JetTree";
      TTree* tr = dynamic_cast<TTree*>(fin->Get(treePath.c_str()));
      if (!tr) {
        cout << "[note] missing tree: " << treePath << " (skip)\n";
        continue;
      }

      // speed up I/O: enable only needed branches
      tr->SetBranchStatus("*", 0);
      tr->SetBranchStatus("mc_pt", 1);
      tr->SetBranchStatus("mc_pt_lead", 1);
      tr->SetBranchStatus("reco_pt_corr", 1);
      tr->SetBranchStatus("reco_pt_lead", 1);
      tr->SetBranchStatus("reco_trigger_match", 1);
      tr->SetBranchStatus("centralityWeight", 1);
      tr->SetBranchStatus("xsecWeight", 1);
      tr->SetBranchStatus("reco_area", 1);
      tr->SetBranchStatus("reco_neutral_fraction", 1);

      // branch addresses
      float mc_pt=0, mc_pt_lead=0;
      float reco_pt_corr=0, reco_pt_lead=0;
      bool  reco_trigger_match=false;
      float centralityWeight=1.0f, xsecWeight=1.0f;
      float reco_area=0.0f, reco_neutral_fraction=0.0f;

      tr->SetBranchAddress("mc_pt", &mc_pt);
      tr->SetBranchAddress("mc_pt_lead", &mc_pt_lead);
      tr->SetBranchAddress("reco_pt_corr", &reco_pt_corr);
      tr->SetBranchAddress("reco_pt_lead", &reco_pt_lead);
      tr->SetBranchAddress("reco_trigger_match", &reco_trigger_match);
      tr->SetBranchAddress("centralityWeight", &centralityWeight);
      tr->SetBranchAddress("xsecWeight", &xsecWeight);
      tr->SetBranchAddress("reco_area", &reco_area);
      tr->SetBranchAddress("reco_neutral_fraction", &reco_neutral_fraction);

      const Long64_t n = tr->GetEntries();
      if (n <= 0) {
        cout << "[note] empty tree: " << treePath << " (skip)\n";
        continue;
      }

      cout << "\n=== " << R << "  " << C << " (entries: " << n << ") ===\n";

      for (int ic = 0; ic < kNPtLeadCuts; ++ic) {
        const double cut = kPtLeadCuts[ic];
        const string tag     = R + "_" + C + Form("_ptlead%.0f", cut);

        cout << "  >> pTlead >= " << cut << " GeV  (tag " << tag << ")\n";

        // --- histos: train/test for closure ---
        TH1D* hMeasTrain = new TH1D(("hMeasTrain_"+tag).c_str(),
            ";reco p_{T}^{corr} [GeV];dN/dp_{T}",
            nbins_meas, bin_meas_edges);
        TH1D* hTrueTrain = new TH1D(("hTrueTrain_"+tag).c_str(),
            ";mc p_{T} [GeV];dN/dp_{T}",
            nbins_truth, bin_truth_edges);
        TH1D* hMeasTest  = (TH1D*)hMeasTrain->Clone(("hMeasTest_"+tag).c_str());
        TH1D* hTrueTest  = (TH1D*)hTrueTrain->Clone(("hTrueTest_"+tag).c_str());

        hMeasTrain->SetDirectory(0);
        hTrueTrain->SetDirectory(0);
        hMeasTest ->SetDirectory(0);
        hTrueTest ->SetDirectory(0);

        // --- full-stat histos (no split) ---
        TH1D* hMeasFull = new TH1D(("hMeasFull_"+tag).c_str(),
            ";reco p_{T}^{corr} [GeV];dN/dp_{T}",
            nbins_meas, bin_meas_edges);
        TH1D* hTrueFull = new TH1D(("hTrueFull_"+tag).c_str(),
            ";mc p_{T} [GeV];dN/dp_{T}",
            nbins_truth, bin_truth_edges);
        hMeasFull->SetDirectory(0);
        hTrueFull->SetDirectory(0);

        // --- response matrices ---
        // For closure (train only)
        TH2D* hRespTrain = new TH2D(("hRespTrain_"+tag).c_str(),
            ";p_{T}^{reco,corr};p_{T}^{mc}",
            nbins_meas, bin_meas_edges,
            nbins_truth, bin_truth_edges);
        hRespTrain->SetDirectory(0);

        // For full-statistics response (train+test together)
        TH2D* hRespFull = new TH2D(("hRespFull_"+tag).c_str(),
            ";p_{T}^{reco,corr};p_{T}^{mc}",
            nbins_meas, bin_meas_edges,
            nbins_truth, bin_truth_edges);
        hRespFull->SetDirectory(0);

        // --- prior (truth-only, using *all* jets passing truth-side cuts) ---
        TH1D* hPrior = new TH1D(("hPrior_"+tag).c_str(),
            ";mc p_{T} [GeV];prior",
            nbins_truth, bin_truth_edges);
        hPrior->SetDirectory(0);

        TH1D* hTrueFull_ptHat[kNPthatBins];
        TH1D* hMeasFull_ptHat[kNPthatBins];
        TH2D* hRespFull_ptHat[kNPthatBins];

        TH1D* hTrueTrain_ptHat[kNPthatBins];
        TH1D* hMeasTrain_ptHat[kNPthatBins];
        TH2D* hRespTrain_ptHat[kNPthatBins];

        TH1D* hTrueTest_ptHat[kNPthatBins];
        TH1D* hMeasTest_ptHat[kNPthatBins];

        TH2D* hRespTest_ptHat[kNPthatBins];

        TH1D* hPrior_ptHat[kNPthatBins];

        TH1D* hTruthMask_ptHat[kNPthatBins];
        TH1D* hTruthSignif_ptHat[kNPthatBins];

        for (int ip = 0; ip < kNPthatBins; ++ip) {
          const TString s = Form("_pthat%d_%s", ip, tag.c_str());

          hTrueFull_ptHat[ip] = new TH1D(("hTrueFull" + s).Data(),
              ";mc p_{T} [GeV];centrality-weighted yield",
              nbins_truth, bin_truth_edges);

          hMeasFull_ptHat[ip] = new TH1D(("hMeasFull" + s).Data(),
              ";reco p_{T}^{corr} [GeV];centrality-weighted yield",
              nbins_meas, bin_meas_edges);

          hRespFull_ptHat[ip] = new TH2D(("hRespFull" + s).Data(),
              ";p_{T}^{reco,corr};p_{T}^{mc}",
              nbins_meas, bin_meas_edges,
              nbins_truth, bin_truth_edges);

          hTrueTrain_ptHat[ip] = new TH1D(("hTrueTrain" + s).Data(),
              ";mc p_{T} [GeV];centrality-weighted yield",
              nbins_truth, bin_truth_edges);

          hMeasTrain_ptHat[ip] = new TH1D(("hMeasTrain" + s).Data(),
              ";reco p_{T}^{corr} [GeV];centrality-weighted yield",
              nbins_meas, bin_meas_edges);

          hRespTrain_ptHat[ip] = new TH2D(("hRespTrain" + s).Data(),
              ";p_{T}^{reco,corr};p_{T}^{mc}",
              nbins_meas, bin_meas_edges,
              nbins_truth, bin_truth_edges);

          hTrueTest_ptHat[ip] = new TH1D(("hTrueTest" + s).Data(),
              ";mc p_{T} [GeV];centrality-weighted yield",
              nbins_truth, bin_truth_edges);

          hMeasTest_ptHat[ip] = new TH1D(("hMeasTest" + s).Data(),
              ";reco p_{T}^{corr} [GeV];centrality-weighted yield",
              nbins_meas, bin_meas_edges);

          hPrior_ptHat[ip] = new TH1D(("hPrior" + s).Data(),
              ";mc p_{T} [GeV];centrality-weighted yield",
              nbins_truth, bin_truth_edges);

          hTruthMask_ptHat[ip] = new TH1D(("hTruthMask" + s).Data(),
              ";mc p_{T} [GeV];mask (0/1)",
              nbins_truth, bin_truth_edges);

          hTruthSignif_ptHat[ip] = new TH1D(("hTruthSignif" + s).Data(),
              ";mc p_{T} [GeV];content/error",
              nbins_truth, bin_truth_edges);

          hRespTest_ptHat[ip] = new TH2D(("hRespTest" + s).Data(),
            ";p_{T}^{reco,corr};p_{T}^{mc}",
            nbins_meas, bin_meas_edges,
            nbins_truth, bin_truth_edges);


          hTrueFull_ptHat[ip]->SetDirectory(0);
          hMeasFull_ptHat[ip]->SetDirectory(0);
          hRespFull_ptHat[ip]->SetDirectory(0);
          hTrueTrain_ptHat[ip]->SetDirectory(0);
          hMeasTrain_ptHat[ip]->SetDirectory(0);
          hRespTrain_ptHat[ip]->SetDirectory(0);
          hTrueTest_ptHat[ip]->SetDirectory(0);
          hMeasTest_ptHat[ip]->SetDirectory(0);
          hPrior_ptHat[ip]->SetDirectory(0);
          hTruthMask_ptHat[ip]->SetDirectory(0);
          hTruthSignif_ptHat[ip]->SetDirectory(0);
          hRespTest_ptHat[ip]->SetDirectory(0);
        }


        TRandom3 rng(kSeed);

        // event loop
        for (Long64_t i = 0; i < n; ++i) {
          if ((i % 200000) == 0)
            cout << "    ["<< tag <<"] " << i << "/" << n << "\r" << std::flush;

          tr->GetEntry(i);

          const int ip = FindPtHatBin((double)xsecWeight);
          if (ip < 0) continue;

          // Ignore first two pThat bins: indices 0 and 1
          if (ip < kFirstPtHatBinToUse) continue;

          const double wCent = (double)centralityWeight;

          const bool haveMC = (mc_pt > 0.0);
          const bool haveReco = (reco_pt_corr > RECO_PTCORR_DUMMY_CUT);

          // ----- fill prior: only MC-side cuts -----
          if (haveMC) {
              hPrior_ptHat[ip]->Fill(mc_pt, wCent);
          }

          // ----- reco-side cuts for response & closure -----
          if (!haveReco) continue;
          if (reco_neutral_fraction > CUT_NEUTRAL_FRACTION) continue;

          // reco quality cuts
          // no trigger requirement here: trigger correction will be external
          // if (!reco_trigger_match) continue;
          if (reco_area < areaMin) continue;

          // dual ptlead cut (both reco & MC)
      //    if (!(reco_pt_lead >= cut && mc_pt_lead >= cut)) continue;
          if (reco_pt_lead < cut) continue;   

          // ===== full-statistics response (temporary, per pThat, centrality-weighted only) =====
          hRespFull_ptHat[ip]->Fill(reco_pt_corr, mc_pt, wCent);
          hMeasFull_ptHat[ip]->Fill(reco_pt_corr, wCent);
          hTrueFull_ptHat[ip]->Fill(mc_pt,        wCent);

          // ===== closure split: train vs test =====
          const bool train = (rng.Uniform() > kTestFrac);
          if (train) {
            hRespTrain_ptHat[ip]->Fill(reco_pt_corr, mc_pt, wCent);
            hMeasTrain_ptHat[ip]->Fill(reco_pt_corr, wCent);
            hTrueTrain_ptHat[ip]->Fill(mc_pt,        wCent);
            } else {
              hRespTest_ptHat[ip]->Fill(reco_pt_corr, mc_pt, wCent);
              hMeasTest_ptHat[ip]->Fill(reco_pt_corr, wCent);
              hTrueTest_ptHat[ip]->Fill(mc_pt,        wCent);
            }
        }
        cout << endl;

        hRespFull->Reset();
        hMeasFull->Reset();
        hTrueFull->Reset();

        hRespTrain->Reset();
        hMeasTrain->Reset();
        hTrueTrain->Reset();

        hMeasTest->Reset();
        hTrueTest->Reset();

        hPrior->Reset();

        bool validTruthBin[kNPthatBins][nbins_truth + 1];

        for (int ip = 0; ip < kNPthatBins; ++ip) {
          cout << "    pThat bin " << ip
              << "  xsecWeight = " << kXsecWeights[ip] << endl;

          for (int jb = 1; jb <= nbins_truth; ++jb) {
            const double c = hTrueFull_ptHat[ip]->GetBinContent(jb);
            const double e = hTrueFull_ptHat[ip]->GetBinError(jb);
            const double signif = (e > 0.0 ? c / e : 0.0);
            const bool keep = (signif > kMinSignif);

            validTruthBin[ip][jb] = keep;
            hTruthMask_ptHat[ip]->SetBinContent(jb, keep ? 1.0 : 0.0);
            hTruthSignif_ptHat[ip]->SetBinContent(jb, signif);

            cout << Form("      truth bin %2d [%.1f, %.1f): c=%10.4e e=%10.4e signif=%7.3f keep=%d",
                        jb,
                        bin_truth_edges[jb-1], bin_truth_edges[jb],
                        c, e, signif, keep ? 1 : 0)
                << endl;
          }
        }


        for (int ip = 0; ip < kNPthatBins; ++ip) {
        const double xw = kXsecWeights[ip] / kNgenEvents[ip];

        // ---- prior ----
        for (int jb = 1; jb <= nbins_truth; ++jb) {
          if (!validTruthBin[ip][jb]) continue;
          const double oldC = hPrior->GetBinContent(jb);
          const double oldE = hPrior->GetBinError(jb);
          const double addC = xw * hPrior_ptHat[ip]->GetBinContent(jb);
          const double addE = xw * hPrior_ptHat[ip]->GetBinError(jb);

          hPrior->SetBinContent(jb, oldC + addC);
          hPrior->SetBinError(jb, std::sqrt(oldE*oldE + addE*addE));
        }
        // ---- truth full/train/test ----
        for (int jb = 1; jb <= nbins_truth; ++jb) {
          if (!validTruthBin[ip][jb]) continue;

          {
            const double oldC = hTrueFull->GetBinContent(jb);
            const double oldE = hTrueFull->GetBinError(jb);
            const double addC = xw * hTrueFull_ptHat[ip]->GetBinContent(jb);
            const double addE = xw * hTrueFull_ptHat[ip]->GetBinError(jb);
            hTrueFull->SetBinContent(jb, oldC + addC);
            hTrueFull->SetBinError(jb, std::sqrt(oldE*oldE + addE*addE));
          }

          {
            const double oldC = hTrueTrain->GetBinContent(jb);
            const double oldE = hTrueTrain->GetBinError(jb);
            const double addC = xw * hTrueTrain_ptHat[ip]->GetBinContent(jb);
            const double addE = xw * hTrueTrain_ptHat[ip]->GetBinError(jb);
            hTrueTrain->SetBinContent(jb, oldC + addC);
            hTrueTrain->SetBinError(jb, std::sqrt(oldE*oldE + addE*addE));
          }

          {
            const double oldC = hTrueTest->GetBinContent(jb);
            const double oldE = hTrueTest->GetBinError(jb);
            const double addC = xw * hTrueTest_ptHat[ip]->GetBinContent(jb);
            const double addE = xw * hTrueTest_ptHat[ip]->GetBinError(jb);
            hTrueTest->SetBinContent(jb, oldC + addC);
            hTrueTest->SetBinError(jb, std::sqrt(oldE*oldE + addE*addE));
          }
        }

        // ---- response full/train (keep only accepted truth columns) ----
        for (int ix = 1; ix <= nbins_meas; ++ix) {
          for (int jy = 1; jy <= nbins_truth; ++jy) {
            if (!validTruthBin[ip][jy]) continue;

            {
              const double oldC = hRespFull->GetBinContent(ix, jy);
              const double oldE = hRespFull->GetBinError(ix, jy);
              const double addC = xw * hRespFull_ptHat[ip]->GetBinContent(ix, jy);
              const double addE = xw * hRespFull_ptHat[ip]->GetBinError(ix, jy);
              hRespFull->SetBinContent(ix, jy, oldC + addC);
              hRespFull->SetBinError(ix, jy, std::sqrt(oldE*oldE + addE*addE));
            }

            {
              const double oldC = hRespTrain->GetBinContent(ix, jy);
              const double oldE = hRespTrain->GetBinError(ix, jy);
              const double addC = xw * hRespTrain_ptHat[ip]->GetBinContent(ix, jy);
              const double addE = xw * hRespTrain_ptHat[ip]->GetBinError(ix, jy);
              hRespTrain->SetBinContent(ix, jy, oldC + addC);
              hRespTrain->SetBinError(ix, jy, std::sqrt(oldE*oldE + addE*addE));
            }
          }
        }

        // ---- measured full/train/test from accepted truth columns only ----
        for (int ix = 1; ix <= nbins_meas; ++ix) {

          double addFullC = 0.0;
          double addFullE2 = 0.0;

          double addTrainC = 0.0;
          double addTrainE2 = 0.0;

          double addTestC = 0.0;
          double addTestE2 = 0.0;

          for (int jy = 1; jy <= nbins_truth; ++jy) {
            if (!validTruthBin[ip][jy]) continue;

            const double cFull = xw * hRespFull_ptHat[ip]->GetBinContent(ix, jy);
            const double eFull = xw * hRespFull_ptHat[ip]->GetBinError(ix, jy);

            addFullC  += cFull;
            addFullE2 += eFull * eFull;

            const double cTrain = xw * hRespTrain_ptHat[ip]->GetBinContent(ix, jy);
            const double eTrain = xw * hRespTrain_ptHat[ip]->GetBinError(ix, jy);

            addTrainC  += cTrain;
            addTrainE2 += eTrain * eTrain;

            const double cTest = xw * hRespTest_ptHat[ip]->GetBinContent(ix, jy);
            const double eTest = xw * hRespTest_ptHat[ip]->GetBinError(ix, jy);

            addTestC  += cTest;
            addTestE2 += eTest * eTest;
          }

          {
            const double oldC = hMeasFull->GetBinContent(ix);
            const double oldE = hMeasFull->GetBinError(ix);
            hMeasFull->SetBinContent(ix, oldC + addFullC);
            hMeasFull->SetBinError(ix, std::sqrt(oldE*oldE + addFullE2));
          }

          {
            const double oldC = hMeasTrain->GetBinContent(ix);
            const double oldE = hMeasTrain->GetBinError(ix);
            hMeasTrain->SetBinContent(ix, oldC + addTrainC);
            hMeasTrain->SetBinError(ix, std::sqrt(oldE*oldE + addTrainE2));
          }

          {
            const double oldC = hMeasTest->GetBinContent(ix);
            const double oldE = hMeasTest->GetBinError(ix);
            hMeasTest->SetBinContent(ix, oldC + addTestC);
            hMeasTest->SetBinError(ix, std::sqrt(oldE*oldE + addTestE2));
          }
        }
      }



        // normalize prior to training truth (as before)
        double intPrior = hPrior->Integral();
        double intTrue  = hTrueTrain->Integral();
        if (intPrior > 0.0 && intTrue > 0.0) {
          hPrior->Scale(intTrue / intPrior);
        }

        cout << "    Integral(Truth train) = " << hTrueTrain->Integral(0, -1) << endl;
        cout << "    Integral(Meas  train) = " << hMeasTrain->Integral(0, -1) << endl;
        cout << "    Integral(Truth test ) = " << hTrueTest ->Integral(0, -1) << endl;
        cout << "    Integral(Meas  test ) = " << hMeasTest ->Integral(0, -1) << endl;
        cout << "    Integral(Prior      ) = " << hPrior    ->Integral(0, -1) << endl;
        cout << "    Integral(Truth full)  = " << hTrueFull ->Integral(0, -1) << endl;
        cout << "    Integral(Meas  full)  = " << hMeasFull ->Integral(0, -1) << endl;

        // --- build responses ---

        // Full-statistics response: this one you will use for unfolding DATA
        RooUnfoldResponse response_full(hMeasFull, hTrueFull, hRespFull);
        response_full.SetName(("response_full_"+tag).c_str());

        // Closure response (train-only): used only for unfolding test sample
        RooUnfoldResponse response_closure(hMeasTrain, hTrueTrain, hRespTrain);
        response_closure.SetName(("response_closure_"+tag).c_str());


        // --------------------------------------------------
        // create output directory once for this tag
        // --------------------------------------------------
        TDirectory* d = fout->mkdir(tag.c_str());
        if (!d) d = fout->GetDirectory(tag.c_str());
        d->cd();

        // common outputs: merged masked response + components
        hRespFull->Write("hRespRecoVsTruth_full");
        response_full.Write("response");

        hMeasFull->Write("hMeasFull");
        hTrueFull->Write("hTrueFull");

        hRespTrain->Write("hRespRecoVsTruth_train");
        response_closure.Write("response_closure");

        hMeasTrain->Write("hMeasTrain");
        hTrueTrain->Write("hTrueTrain");
        hMeasTest ->Write("hMeasTest");
        hTrueTest ->Write("hTrueTest");
        hPrior    ->Write("hPrior");

        // save per-pThat debug objects
        if (kSavePtHatDebug) {
          TDirectory* dbg = d->mkdir("debug_pthat");
          if (!dbg) dbg = d->GetDirectory("debug_pthat");
          dbg->cd();

          for (int ip = 0; ip < kNPthatBins; ++ip) {
            TDirectory* dp = dbg->mkdir(Form("pthat_%02d", ip));
            if (!dp) dp = dbg->GetDirectory(Form("pthat_%02d", ip));
            dp->cd();

            hTrueFull_ptHat[ip]->Write("hTrueFull_centWeight");
            hMeasFull_ptHat[ip]->Write("hMeasFull_centWeight");
            hRespFull_ptHat[ip]->Write("hRespFull_centWeight");

            hTrueTrain_ptHat[ip]->Write("hTrueTrain_centWeight");
            hMeasTrain_ptHat[ip]->Write("hMeasTrain_centWeight");
            hRespTrain_ptHat[ip]->Write("hRespTrain_centWeight");

            hRespTest_ptHat[ip]->Write("hRespTest_centWeight");

            hTrueTest_ptHat[ip]->Write("hTrueTest_centWeight");
            hMeasTest_ptHat[ip]->Write("hMeasTest_centWeight");

            hPrior_ptHat[ip]->Write("hPrior_centWeight");
            hTruthMask_ptHat[ip]->Write("hTruthMask");
            hTruthSignif_ptHat[ip]->Write("hTruthSignif");
          }

          d->cd();
        }

      // --- BAYESIAN unfolding ---

        if (m == "BAYES") {
          cout << "Doing Bayesian unfolding..." << endl;
          static const int kBayesIters[] = {1, 2, 3, 4, 5, 6, 7, 8};
          static const int kNBayesIters = sizeof(kBayesIters)/sizeof(kBayesIters[0]);


          // --- unfolding (closure) with explicit prior ---
          vector<TH1D*> unfolded(kNBayesIters, 0);
          for (int ib = 0; ib < kNBayesIters; ++ib) {
            RooUnfoldBayes u(&response_closure, hMeasTest,
                            kBayesIters[ib], false, hPrior);
            TH1D* hunf = (TH1D*)u.Hunfold();
            hunf->SetDirectory(0);
            hunf->SetName(Form("Unfolded_%s_iter%d", tag.c_str(), kBayesIters[ib]));
            unfolded[ib] = hunf;
          }

          // Rebin unfolded to truth binning (for plotting & ratio)
          vector<TH1D*> unfoldedTruth(kNBayesIters, 0);
          for (int ib = 0; ib < kNBayesIters; ++ib) {
            if (!unfolded[ib]) continue;
            unfoldedTruth[ib] = (TH1D*)unfolded[ib]->Rebin(
                nbins_truth,
                Form("UnfTruthBins_%d_%s", kBayesIters[ib], tag.c_str()),
                bin_truth_edges);
            unfoldedTruth[ib]->SetDirectory(0);
          }

         
          // ========================= PLOTTING (closure) =========================

          // Distinct (non-blending) styles for Bayes iterations
          static const int kUnfCols[8] = {kBlack,kRed+1,kBlue+1,kGreen+2,kCyan+2,kOrange+7,kMagenta+2,kViolet+1};
          static const int kUnfMarks[8]  = { 20, 21, 22, 33, 24, 25, 26, 27 };

          TCanvas* c = new TCanvas(("c_"+tag).c_str(), "", 800, 1000);

          // --- manual pads (instead of Divide) so we can control spacing ---
          TPad* pTop = new TPad(("pTop_"+tag).c_str(), "", 0.0, 0.30, 1.0, 1.0);
          TPad* pBot = new TPad(("pBot_"+tag).c_str(), "", 0.0, 0.00, 1.0, 0.30);

          // tighter margins for slides
          pTop->SetLeftMargin(0.12);
          pTop->SetRightMargin(0.03);
          pTop->SetTopMargin(0.05);
          pTop->SetBottomMargin(0.02);   // tiny: bottom pad will carry x labels

          pBot->SetLeftMargin(0.12);
          pBot->SetRightMargin(0.03);
          pBot->SetTopMargin(0.02);
          pBot->SetBottomMargin(0.30);   // room for x-axis title/labels

          pTop->Draw();
          pBot->Draw();

          // ---------- top pad: shapes in TRUTH binning ----------
          pTop->cd();
          gPad->SetLogy();

          TH1D* hT_plot = (TH1D*)hTrueTest->Clone(("hTrueW_"+tag).c_str());
          hT_plot->SetMarkerStyle(20);
          hT_plot->SetMarkerColor(kBlack);
          hT_plot->SetLineColor(kBlack);

          TH1D* hM_truth = (TH1D*)hMeasTest->Rebin(
            nbins_truth, ("hMeasTruthBins_"+tag).c_str(), bin_truth_edges
          );
          hM_truth->SetMarkerStyle(24);
          hM_truth->SetMarkerColor(kBlue+2);
          hM_truth->SetLineColor(kBlue+2);

          // Top pad: hide x labels to save space (bottom pad will have them)
          hT_plot->GetXaxis()->SetLabelSize(0);
          hT_plot->GetXaxis()->SetTitleSize(0);

          hT_plot->GetYaxis()->SetTitleOffset(1.2);
          hT_plot->Draw("E1");
          hM_truth->Draw("E1 SAME");

          TLegend* leg = new TLegend(0.55,0.60,0.90,0.90);
          leg->SetBorderSize(0);
          leg->SetFillStyle(0);
          leg->SetTextSize(0.035);

          leg->AddEntry(hT_plot, "Truth (test)", "lp");
          leg->AddEntry(hM_truth, "Measured (test)", "lp");

          for (int ib = 0; ib < kNBayesIters; ++ib) {
            if (!DrawBayesIter(kBayesIters[ib])) continue;

            TH1D* w = unfoldedTruth[ib];
            if (!w) continue;

            w->SetMarkerStyle(kUnfMarks[ib]);
            w->SetMarkerColor(kUnfCols[ib]);
            w->SetLineColor(kUnfCols[ib]);

            w->Draw("E1 SAME");
            leg->AddEntry(w, Form("Bayes %d it.", kBayesIters[ib]), "lp");
          }
          leg->Draw();

          {
            TLatex lat;
            lat.SetNDC(true);
            lat.SetTextFont(42);
            lat.SetTextSize(0.040);
            lat.DrawLatex(0.16, 0.28, "Au+Au  #sqrt{#it{s}_{NN}} = 200 GeV");
            lat.DrawLatex(0.16, 0.22,
                          Form("#it{R} = %.1f, %s", Rval, NiceCentLabel(C).Data()));
            lat.DrawLatex(0.16, 0.16, Form("#it{p}_{T}^{lead} #geq %.0f GeV/#it{c}", cut));
            lat.DrawLatex(0.16, 0.10, "Unfolding method: Bayesian");
          }

          // ---------- bottom pad: ratio unfolded / truth ----------
          pBot->cd();
          //pBot->SetGridy();

          TH1D* firstRatio = 0;
          for (int ib = 0; ib < kNBayesIters; ++ib) {
            if (!DrawBayesIter(kBayesIters[ib])) continue;

            TH1D* hunf_reb = unfoldedTruth[ib];
            if (!hunf_reb) continue;

            TH1D* r = (TH1D*)hunf_reb->Clone(Form("ratio_%d_%s", kBayesIters[ib], tag.c_str()));
            r->SetDirectory(0);
            r->Divide(hTrueTest);

            // match styles to top pad
            r->SetMarkerStyle(kUnfMarks[ib]);
            r->SetMarkerColor(kUnfCols[ib]);
            r->SetLineColor(kUnfCols[ib]);

            if (!firstRatio) {
              firstRatio = r;
              firstRatio->SetTitle("");
              firstRatio->GetYaxis()->SetTitle("Unfolded / Truth");
              firstRatio->GetYaxis()->SetRangeUser(0.4, 1.6);

              // bottom pad must carry x-axis title + labels
              firstRatio->GetXaxis()->SetTitle("#it{p}_{T} (GeV)");
              firstRatio->GetXaxis()->SetTitleSize(0.11);
              firstRatio->GetXaxis()->SetLabelSize(0.09);
              firstRatio->GetXaxis()->SetTitleOffset(1.05);

              firstRatio->GetYaxis()->SetTitleSize(0.09);
              firstRatio->GetYaxis()->SetLabelSize(0.08);
              firstRatio->GetYaxis()->SetTitleOffset(0.65);

              firstRatio->Draw("E1");
            } else {
              r->Draw("E1 SAME");
            }
          }
          
          if (firstRatio) {
            double xmin = firstRatio->GetXaxis()->GetXmin();
            double xmax = firstRatio->GetXaxis()->GetXmax();
            TLine* ln1 = new TLine(xmin, 1.05, xmax, 1.05);
            ln1->SetLineStyle(2);
            ln1->SetLineColor(kGray);
            ln1->Draw();
            TLine* ln2 = new TLine(xmin, 0.95, xmax, 0.95);
            ln2->SetLineStyle(2);
            ln2->SetLineColor(kGray);
            ln2->Draw();
            TLine* ln3 = new TLine(xmin, 1.00, xmax, 1.00);
            ln3->SetLineStyle(1);
            ln3->SetLineColor(kBlack);
            ln3->Draw();
            TLine* ln4 = new TLine(xmin, 0.90, xmax, 0.90);
            ln4->SetLineStyle(2);
            ln4->SetLineColor(kGray+2);
            ln4->Draw();
            TLine* ln5 = new TLine(xmin, 1.1, xmax, 1.1);
            ln5->SetLineStyle(2);
            ln5->SetLineColor(kGray+2);
            ln5->Draw();
          }

          // ---------- save ----------
          const string pdfPath = string(outDir) + "/BAYES_closure_" + tag + ".pdf";

          d->cd();
            // optionally save unfolded spectra for closure checks
            for (int ib = 0; ib < kNBayesIters; ++ib) {
              if (unfolded[ib]) {
                unfolded[ib]->Write(
                  Form("Unfolded_iter%d", kBayesIters[ib])
                );
              }
            }
            c->SaveAs(pdfPath.c_str());
            {
              string pngPath = pdfPath;
              const size_t pos = pngPath.rfind(".pdf");
              if (pos != string::npos) pngPath.replace(pos, 4, ".png");
              c->SaveAs(pngPath.c_str());
            }


            // ========================= ITERATION STABILITY PLOT =========================

            TCanvas* cstab = new TCanvas(("c_stability_"+tag).c_str(), "", 800, 1000);

            TPad* pTopS = new TPad(("pTop_stability_"+tag).c_str(), "", 0.0, 0.30, 1.0, 1.0);
            TPad* pBotS = new TPad(("pBot_stability_"+tag).c_str(), "", 0.0, 0.00, 1.0, 0.30);

            pTopS->SetLeftMargin(0.12);
            pTopS->SetRightMargin(0.03);
            pTopS->SetTopMargin(0.05);
            pTopS->SetBottomMargin(0.02);

            pBotS->SetLeftMargin(0.12);
            pBotS->SetRightMargin(0.03);
            pBotS->SetTopMargin(0.02);
            pBotS->SetBottomMargin(0.30);

            pTopS->Draw();
            pBotS->Draw();

            // ---------- top: unfolded spectra ----------
            pTopS->cd();
            gPad->SetLogy();

            TH1D* frameStab = 0;

            TLegend* legStab = new TLegend(0.55, 0.55, 0.90, 0.90);
            legStab->SetBorderSize(0);
            legStab->SetFillStyle(0);
            legStab->SetTextSize(0.030);

            static const int kStabCols[] = {
              kBlack, kRed+1, kAzure+2, kGreen+2, kOrange+7,
              kMagenta+2, kCyan+2, kViolet+1, kGray+2, kBlue+1
            };

            static const int kStabMarks[] = {
              20, 21, 22, 23, 33,
              34, 29, 47, 24, 25
            };

            for (int ib = 0; ib < kNBayesIters; ++ib) {
              TH1D* h = unfoldedTruth[ib];
              if (!h) continue;

              const int col = kStabCols[ib % 10];
              const int mar = kStabMarks[ib % 10];

              h->SetMarkerStyle(mar);
              h->SetMarkerColor(col);
              h->SetLineColor(col);

              h->GetXaxis()->SetLabelSize(0);
              h->GetXaxis()->SetTitleSize(0);
              h->GetYaxis()->SetTitle("Unfolded yield");
              h->GetYaxis()->SetTitleOffset(1.2);

              if (!frameStab) {
                frameStab = h;
                frameStab->SetTitle("");
                frameStab->Draw("E1");
              } else {
                h->Draw("E1 SAME");
              }

              legStab->AddEntry(h, Form("Bayes %d it.", kBayesIters[ib]), "lp");
            }

            legStab->Draw();

            {
              TLatex lat;
              lat.SetNDC(true);
              lat.SetTextFont(42);
              lat.SetTextSize(0.040);
              lat.DrawLatex(0.16, 0.28, "Au+Au  #sqrt{#it{s}_{NN}} = 200 GeV");
              lat.DrawLatex(0.16, 0.22,
                            Form("#it{R} = %.1f, %s", Rval, NiceCentLabel(C).Data()));
              lat.DrawLatex(0.16, 0.16,
                            Form("#it{p}_{T}^{lead} #geq %.0f GeV/#it{c}", cut));
              lat.DrawLatex(0.16, 0.10, "Bayesian unfolding stability");
            }

            // ---------- bottom: consecutive iteration ratios ----------
            pBotS->cd();

            TH1D* firstStabRatio = 0;

            for (int ib = 1; ib < kNBayesIters; ++ib) {
              TH1D* hNum = unfoldedTruth[ib];
              TH1D* hDen = unfoldedTruth[ib-1];

              if (!hNum || !hDen) continue;

              TH1D* r = (TH1D*)hNum->Clone(
                  Form("ratio_iter%d_over_iter%d_%s",
                      kBayesIters[ib], kBayesIters[ib-1], tag.c_str())
              );

              r->SetDirectory(0);
              r->Divide(hDen);

              const int col = kStabCols[ib % 10];
              const int mar = kStabMarks[ib % 10];

              r->SetMarkerStyle(mar);
              r->SetMarkerColor(col);
              r->SetLineColor(col);

              r->SetTitle("");
              r->GetYaxis()->SetTitle("Iter. ratio");
              r->GetYaxis()->SetRangeUser(0.8, 1.2);

              r->GetXaxis()->SetTitle("#it{p}_{T}^{truth} (GeV/#it{c})");
              r->GetXaxis()->SetTitleSize(0.11);
              r->GetXaxis()->SetLabelSize(0.09);
              r->GetXaxis()->SetTitleOffset(1.05);

              r->GetYaxis()->SetTitleSize(0.09);
              r->GetYaxis()->SetLabelSize(0.08);
              r->GetYaxis()->SetTitleOffset(0.65);
              r->GetYaxis()->SetNdivisions(505);

              if (!firstStabRatio) {
                firstStabRatio = r;
                firstStabRatio->Draw("E1");
              } else {
                r->Draw("E1 SAME");
              }

              r->Write(
                Form("BayesStabilityRatio_iter%d_over_iter%d",
                    kBayesIters[ib], kBayesIters[ib-1])
              );
            }

            if (firstStabRatio) {
              const double xmin = firstStabRatio->GetXaxis()->GetXmin();
              const double xmax = firstStabRatio->GetXaxis()->GetXmax();

              TLine* l1 = new TLine(xmin, 1.00, xmax, 1.00);
              l1->SetLineColor(kBlack);
              l1->SetLineStyle(1);
              l1->Draw();

              TLine* l2 = new TLine(xmin, 1.05, xmax, 1.05);
              l2->SetLineColor(kGray+1);
              l2->SetLineStyle(2);
              l2->Draw();

              TLine* l3 = new TLine(xmin, 0.95, xmax, 0.95);
              l3->SetLineColor(kGray+1);
              l3->SetLineStyle(2);
              l3->Draw();

              TLine* l4 = new TLine(xmin, 1.10, xmax, 1.10);
              l4->SetLineColor(kGray+2);
              l4->SetLineStyle(2);
              l4->Draw();

              TLine* l5 = new TLine(xmin, 0.90, xmax, 0.90);
              l5->SetLineColor(kGray+2);
              l5->SetLineStyle(2);
              l5->Draw();
            }

            // ---------- save stability plot ----------
            const string pdfPathStab = string(outDir) + "/BAYES_stability_" + tag + ".pdf";

            d->cd();
            cstab->SaveAs(pdfPathStab.c_str());

            {
              string pngPathStab = pdfPathStab;
              const size_t pos = pngPathStab.rfind(".pdf");
              if (pos != string::npos) pngPathStab.replace(pos, 4, ".png");
              cstab->SaveAs(pngPathStab.c_str());
            }

            delete cstab;
            delete legStab;

            // ---------- cleanup ----------
            delete c;
            delete hT_plot;
            delete hM_truth;
            delete leg;
            for (size_t ib = 0; ib < unfolded.size(); ++ib) {
              if (unfolded[ib]) delete unfolded[ib];
            }
            for (size_t ib = 0; ib < unfoldedTruth.size(); ++ib) {
              if (unfoldedTruth[ib]) delete unfoldedTruth[ib];
            }
        } // end of BAYESIAN unfolding



      // --- SVD unfolding ---
       else if (m == "SVD") {
          cout << "Doing SVD unfolding..." << endl;
          //const int kRegValues[] = {2, 3, 4, 5, 6, 7, 8, 9};
          //const int kRegValues[] = {2, 3, 4, 5, 6};
          const int kRegValues[] = {3, 4, 5};
          const int nSVD = sizeof(kRegValues)/sizeof(kRegValues[0]);

          vector<TH1D*> unfoldedSVD(nSVD, nullptr);
          for (int ir = 0; ir < nSVD; ++ir) {
            const int reg = kRegValues[ir];
            RooUnfoldSvd u(&response_closure, hMeasTest, reg);
            // optionally switch on regularization via SVD settings:
            // u.SetRegParam(reg); // not necessary with ctor but shown for clarity
            
            TH1D* hunf = dynamic_cast<TH1D*>(u.Hunfold());
            if (!hunf) {
              cout << "[error] SVD unfolding failed for reg=" << reg << endl;
              continue;
            }
            hunf->SetDirectory(nullptr);
            hunf->SetName(Form("Unfolded_SVD_%s_reg%d", tag.c_str(), reg));
            unfoldedSVD[ir] = hunf;
              // Now access diagnostics — Hunfold() must have been called first
               auto* svdImpl = u.Impl();
                if (svdImpl) {
                    TVectorD sv = svdImpl->GetSV();
                    cout << "=== Singular values (reg=" << reg << ") ===" << endl;
                    for (int i = 0; i < sv.GetNrows(); i++)
                        cout << Form("  Mode %2d : SV = %10.6e", i, sv[i]) << endl;

                    TH1* d = svdImpl->GetD();
                    if (d) {
                        cout << "=== d-vector (reg=" << reg << ") ===" << endl;
                        for (int i = 1; i <= d->GetNbinsX(); i++)
                            cout << Form("  Mode %2d : |d_i| = %10.6e", i, fabs(d->GetBinContent(i))) << endl;
                    } else {
                        cout << "[note] d still null after Hunfold() for reg=" << reg << endl;
                    }
                }

              // // Access the internal SVD implementation
              // auto* impl = u.Impl();
              // if (!impl) {
              //   cout << "[warning] Could not access SVD implementation" << endl;
              //   continue;
              // }

              // cout << "\n=== SVD diagnostics for k_reg = " << reg << " ===" << endl;
              
              // // Get the singular values and d-vector from implementation
              // TVectorD sv = impl->GetSV();
              // TH1* d = impl->GetD();
              
               

              //   cout << "=== Singular values ===" << endl;
              //   for (int i = 0; i < sv.GetNrows(); i++) {
              //       cout << Form("Mode %2d : SV = %10.6e", i, sv[i]) << endl;
              //   }

              //   for (int i = 0; i < d->GetNbinsX(); i++){
              //     cout << Form("Mode %2d : |d_i| = %10.6e", i+1, fabs(d->GetBinContent(i+1))) << endl;
              //   }

                 cout << "Getting impl..." << endl;
                auto* impl = u.Impl();
                if (!impl) { cout << "[warning] impl is null" << endl; continue; }

                cout << "Getting SV..." << endl;
                TVectorD sv = impl->GetSV();
                cout << "SV rows: " << sv.GetNrows() << endl;

                cout << "=== Singular values ===" << endl;
                for (int i = 0; i < sv.GetNrows(); i++)
                  cout << Form("Mode %2d : SV = %10.6e", i, sv[i]) << endl;

                cout << "Getting D..." << endl;
                TH1* d = impl->GetD();
                cout << "D pointer: " << (void*)d << endl;
                if (!d) { cout << "[warning] d is null" << endl; continue; }

                cout << "D nbins: " << d->GetNbinsX() << endl;
                for (int i = 1; i <= d->GetNbinsX(); i++)
                  cout << Form("Mode %2d : |d_i| = %10.6e", i, fabs(d->GetBinContent(i))) << endl;


          }

          // Rebin unfolded to truth binning (for plotting & ratio)
          vector<TH1D*> unfoldedTruth(nSVD, nullptr);
          for (int iSVD = 0; iSVD < nSVD; ++iSVD) {
            if (!unfoldedSVD[iSVD]) continue;
            unfoldedTruth[iSVD] = (TH1D*)unfoldedSVD[iSVD]->Rebin(
                nbins_truth,
                Form("UnfTruthBins_%d_%s", kRegValues[iSVD], tag.c_str()),
                bin_truth_edges);
            unfoldedTruth[iSVD]->SetDirectory(nullptr);
        }

          // ========================= PLOTTING (closure) =========================

          // Distinct (non-blending) styles for SVD regularization values
          static const int kUnfCols[5]   = { kRed+1, kAzure+2, kGreen+2, kOrange+7, kMagenta+2,};
          static const int kUnfMarks[5]  = { 20, 21, 22, 33, 34 };

          TCanvas* c = new TCanvas(("c_"+tag).c_str(), "", 800, 1000);

          // --- manual pads (instead of Divide) so we can control spacing ---
          TPad* pTop = new TPad(("pTop_"+tag).c_str(), "", 0.0, 0.30, 1.0, 1.0);
          TPad* pBot = new TPad(("pBot_"+tag).c_str(), "", 0.0, 0.00, 1.0, 0.30);

          // tighter margins for slides
          pTop->SetLeftMargin(0.12);
          pTop->SetRightMargin(0.03);
          pTop->SetTopMargin(0.02);
          pTop->SetBottomMargin(0.02);   // tiny: bottom pad will carry x labels

          pBot->SetLeftMargin(0.12);
          pBot->SetRightMargin(0.03);
          pBot->SetTopMargin(0.02);
          pBot->SetBottomMargin(0.30);   // room for x-axis title/labels

          pTop->Draw();
          pBot->Draw();

          // ---------- top pad: shapes in TRUTH binning ----------
          pTop->cd();
          gPad->SetLogy();

          TH1D* hT_plot_svd = (TH1D*)hTrueTest->Clone(("hTrueW_"+tag).c_str());
          hT_plot_svd->SetMarkerStyle(20);
          hT_plot_svd->SetMarkerColor(kBlack);
          hT_plot_svd->SetLineColor(kBlack);

          TH1D* hM_truth_svd = (TH1D*)hMeasTest->Rebin(
            nbins_truth, ("hMeasTruthBins_"+tag).c_str(), bin_truth_edges
          );
          hM_truth_svd->SetMarkerStyle(24);
          hM_truth_svd->SetMarkerColor(kBlue+2);
          hM_truth_svd->SetLineColor(kBlue+2);

          // Top pad: hide x labels to save space (bottom pad will have them)
          hT_plot_svd->GetXaxis()->SetLabelSize(0);
          hT_plot_svd->GetXaxis()->SetTitleSize(0);

          hT_plot_svd->GetYaxis()->SetTitleOffset(1.2);
          hT_plot_svd->Draw("E1");
          hM_truth_svd->Draw("E1 SAME");

          TLegend* leg_svd = new TLegend(0.55,0.60,0.90,0.90);
          leg_svd->SetBorderSize(0);
          leg_svd->SetFillStyle(0);
          leg_svd->SetTextSize(0.035);

          leg_svd->AddEntry(hT_plot_svd, "Truth (test)", "lp");
          leg_svd->AddEntry(hM_truth_svd, "Measured (test)", "lp");

          for (int iSVD = 0; iSVD < nSVD; ++iSVD) {
            TH1D* w = unfoldedTruth[iSVD];
            if (!w) continue;

            w->SetMarkerStyle(kUnfMarks[iSVD]);
            w->SetMarkerColor(kUnfCols[iSVD]);
            w->SetLineColor(kUnfCols[iSVD]);

            w->Draw("E1 SAME");
            leg_svd->AddEntry(w, Form("k_reg = %d ", kRegValues[iSVD]), "lp");
          }
          leg_svd->Draw();

          {
            TLatex lat;
            lat.SetNDC(true);
            lat.SetTextFont(42);
            lat.SetTextSize(0.040);
            lat.DrawLatex(0.16, 0.28, "Au+Au  #sqrt{#it{s}_{NN}} = 200 GeV");
            lat.DrawLatex(0.16, 0.22,
                          Form("#it{R} = %.1f, %s", Rval, NiceCentLabel(C).Data()));
            lat.DrawLatex(0.16, 0.16, Form("#it{p}_{T}^{lead} #geq %.0f GeV/#it{c}", cut));
            lat.DrawLatex(0.16, 0.10, "Unfolding method: SVD");
          }

          // ---------- bottom pad: ratio unfolded / truth ----------
          pBot->cd();
          //pBot->SetGridy();

          TH1D* firstRatio = 0;
          for (int iSVD = 0; iSVD < nSVD; ++iSVD) {
            TH1D* hunf_reb = unfoldedTruth[iSVD];
            if (!hunf_reb) continue;

            TH1D* r = (TH1D*)hunf_reb->Clone(Form("ratio_%d_%s", kRegValues[iSVD], tag.c_str()));
            r->SetDirectory(0);
            r->Divide(hTrueTest);

            // match styles to top pad
            r->SetMarkerStyle(kUnfMarks[iSVD]);
            r->SetMarkerColor(kUnfCols[iSVD]);
            r->SetLineColor(kUnfCols[iSVD]);

            if (!firstRatio) {
              firstRatio = r;
              firstRatio->SetTitle("");
              firstRatio->GetYaxis()->SetTitle("Unfolded / Truth");
              firstRatio->GetYaxis()->SetRangeUser(0.4, 1.6);

              // bottom pad must carry x-axis title + labels
              firstRatio->GetXaxis()->SetTitle("#it{p}_{T} (GeV)");
              firstRatio->GetXaxis()->SetTitleSize(0.11);
              firstRatio->GetXaxis()->SetLabelSize(0.09);
              firstRatio->GetXaxis()->SetTitleOffset(1.05);

              firstRatio->GetYaxis()->SetTitleSize(0.09);
              firstRatio->GetYaxis()->SetLabelSize(0.08);
              firstRatio->GetYaxis()->SetTitleOffset(0.65);

              firstRatio->Draw("E1");
            } else {
              r->Draw("E1 SAME");
            }
          }

          if (firstRatio) {
            double xmin = firstRatio->GetXaxis()->GetXmin();
            double xmax = firstRatio->GetXaxis()->GetXmax();
            TLine* ln1 = new TLine(xmin, 1.05, xmax, 1.05);
            ln1->SetLineStyle(2);
            ln1->SetLineColor(kGray);
            ln1->Draw();
            TLine* ln2 = new TLine(xmin, 0.95, xmax, 0.95);
            ln2->SetLineStyle(2);
            ln2->SetLineColor(kGray);
            ln2->Draw();
            TLine* ln3 = new TLine(xmin, 1.00, xmax, 1.00);
            ln3->SetLineStyle(1);
            ln3->SetLineColor(kBlack);
            ln3->Draw();
            TLine* ln4 = new TLine(xmin, 0.90, xmax, 0.90);
            ln4->SetLineStyle(2);
            ln4->SetLineColor(kGray+2);
            ln4->Draw();
            TLine* ln5 = new TLine(xmin, 1.1, xmax, 1.1);
            ln5->SetLineStyle(2);
            ln5->SetLineColor(kGray+2);
            ln5->Draw();
          }

          // ---------- save ----------
          const string pdfPath_svd = string(outDir) + "/SVD_closure_" + tag + ".pdf";

          d->cd();

        // optionally save unfolded spectra for closure checks
        for (int iSVD = 0; iSVD < nSVD; ++iSVD) {
          if (unfoldedSVD[iSVD]) {
            unfoldedSVD[iSVD]->Write(
              Form("Unfolded_SVD_iter%d", kRegValues[iSVD])
            );
          }
        }

        // save closure plot as before
        c->SaveAs(pdfPath_svd.c_str());
          {
            string pngPath_svd = pdfPath_svd;
            const size_t pos = pngPath_svd.rfind(".pdf");
            if (pos != string::npos) pngPath_svd.replace(pos, 4, ".png");
            c->SaveAs(pngPath_svd.c_str());
          }


          // ---------- cleanup ----------
        delete c;
        delete hT_plot_svd;
        delete hM_truth_svd;
        delete leg_svd;
        for (auto* u  : unfoldedSVD)      delete u;
        for (auto* ut : unfoldedTruth) delete ut;

       } //--- end of SVD unfolding ---

        delete hRespTrain;
        delete hRespFull;
        delete hMeasTrain;
        delete hTrueTrain;
        delete hMeasTest;
        delete hTrueTest;
        delete hMeasFull;
        delete hTrueFull;
        delete hPrior;

        for (int ip = 0; ip < kNPthatBins; ++ip) {
          delete hTrueFull_ptHat[ip];
          delete hMeasFull_ptHat[ip];
          delete hRespFull_ptHat[ip];

          delete hTrueTrain_ptHat[ip];
          delete hMeasTrain_ptHat[ip];
          delete hRespTrain_ptHat[ip];
          delete hRespTest_ptHat[ip];

          delete hTrueTest_ptHat[ip];
          delete hMeasTest_ptHat[ip];

          delete hPrior_ptHat[ip];
          delete hTruthMask_ptHat[ip];
          delete hTruthSignif_ptHat[ip];
        }

      } // ptlead cuts
    } // centralities
  } // radii

  fout->Write();
  fout->Close();
  delete fout;

  fin->Close();
  delete fin;

  cout << "\nDone. Responses (closure + full) written to: "
       << outRootPath << endl;
}
