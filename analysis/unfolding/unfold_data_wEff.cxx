// unfold_data.C
#include "RooUnfoldBayes.h"
#include "RooUnfoldResponse.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TMatrixD.h"
#include "TDirectory.h"
#include "TParameter.h"
#include "TAxis.h"
#include "TList.h"

#include <vector>
#include <string>
#include <iostream>
#include <cmath>

using std::string;
using std::vector;
using std::cout;
using std::endl;

// ------------------------- config -----------------------------

// Same iterations as closure
static const int kBayesIterDefault = 4;

// measured & truth binning (must be IDENTICAL to what you used for response)
static const int nbins_meas = 24;
static const double bin_meas_edges[nbins_meas+1] = {
  -100,-80,-60,-40,-20,-10,-5,-2.5,0,2.5,5,7.5,10,12.5,15,17.5,
  20,22.5,25,27.5,30,35,40,50,60
};

static const int nbins_truth = 10;
static const double bin_truth_edges[nbins_truth+1] = {
  0,5,10,15,20,25,30,35,40,50,60
};

// same centrality & R labels as before
static const vector<string> kCentralities =
  {"CENT_0_10", "MID_20_40", "PERI_60_80"};
static const vector<string> kRadii =
  {"R0.2", "R0.3", "R0.4"};

// pt_lead cuts (same as in closure)
static const double kPtLeadCuts[] = {0.0, 5.0, 7.0, 9.0};
static const int    kNPtLeadCuts  = sizeof(kPtLeadCuts)/sizeof(kPtLeadCuts[0]);

// ------------------ jet-quality cuts -------------------------

static const double CUT_AREA_02 = 0.07;  // R = 0.2
static const double CUT_AREA_03 = 0.20;  // R = 0.3
static const double CUT_AREA_04 = 0.40;  // R = 0.4
static const double CUT_NEUTRAL_FRACTION = 0.95;

// --------------------------------------------------------------

static void EnsureDir(const string& path){
  if (gSystem->AccessPathName(path.c_str()))
    gSystem->mkdir(path.c_str(), /*recursive=*/true);
}

static bool SameBinning1D(const TH1* a, const TH1* b, double tol = 1e-6);

// NOTE: the response now already contains trigger/matching efficiency
// and fake/purity information through hMeasFull, hTrueFull and hRespFull.
// Therefore unfold_data.C must not apply h_trig_eff_reco, h_pur_eff_reco,
// or h_match_eff_truth externally anymore.

static TH1D* MakeFinalInvariantSpectrum(const TH1D* src,
                                         const char* outName,
                                         double nEqMB,
                                         double jetR,
                                         bool verbose = true);

static TH1* GetHistAnyDir(TFile* f, const std::string& name)
{
  if (!f) return nullptr;

  TH1* h = dynamic_cast<TH1*>(f->Get(name.c_str()));
  if (h) return h;

  h = dynamic_cast<TH1*>(f->Get(("QA_histograms/" + name).c_str()));
  if (h) return h;

  h = dynamic_cast<TH1*>(f->Get(("eventQA/" + name).c_str()));
  if (h) return h;

  TList* qa = dynamic_cast<TList*>(f->Get("QA_histograms"));
  if (qa) {
    h = dynamic_cast<TH1*>(qa->FindObject(name.c_str()));
    if (h) return h;
  }

  TList* makerList = dynamic_cast<TList*>(f->Get("StPicoHFJetMaker"));
  if (makerList) {
    h = dynamic_cast<TH1*>(makerList->FindObject(name.c_str()));
    if (h) return h;
  }

  return nullptr;
}

struct EqMbNormInfo {
  double nEqNoR;
  double nEqWithR;
  double nHTAcc;

  EqMbNormInfo()
      : nEqNoR(0.0), nEqWithR(0.0), nHTAcc(0.0) {}
};

static double CalcEqMbWithR(TFile* fHT,
                            TFile* fMB,
                            const std::string& centTag,
                            EqMbNormInfo* normInfo = nullptr,
                            TDirectory* outDir = nullptr)
{
  TH1* hHTAll = GetHistAnyDir(fHT, "hrunId_trig");
  TH1* hMBAll = GetHistAnyDir(fMB, "hrunId_trig");

  TH1* hHTAcc = GetHistAnyDir(fHT, Form("hrunId_acc_%s", centTag.c_str()));
  TH1* hMBAcc = GetHistAnyDir(fMB, Form("hrunId_acc_weighted_%s", centTag.c_str()));

  TH1* hEqNoR = GetHistAnyDir(fHT, Form("hrunId_eqMb_%s", centTag.c_str()));

  if (!hHTAll || !hMBAll || !hHTAcc || !hMBAcc || !hEqNoR) {
    cout << "[error] Missing normalization histograms for " << centTag << endl;
    cout << "        Need hrunId, hrunId_acc_" << centTag
         << ", hrunId_eqMb_" << centTag << " in HT file and hrunId, hrunId_acc_"
         << centTag << " in MB file." << endl;
    return -1.0;
  }

  TH1D* hRHT = (TH1D*)hEqNoR->Clone(Form("hR_HT_%s", centTag.c_str()));
  TH1D* hRMB = (TH1D*)hEqNoR->Clone(Form("hR_MB_%s", centTag.c_str()));
  TH1D* hRatio = (TH1D*)hEqNoR->Clone(Form("hRratio_MB_over_HT_%s", centTag.c_str()));
  TH1D* hEqWithR = (TH1D*)hEqNoR->Clone(Form("hrunId_eqMb_withR_%s", centTag.c_str()));

  hRHT->Reset();
  hRMB->Reset();
  hRatio->Reset();
  hEqWithR->Reset();

  double neqNoR = 0.0;
  double neqWithR = 0.0;
  double nHTAcc = 0.0;

  for (int ib = 1; ib <= hEqNoR->GetNbinsX(); ++ib) {
    const double eqNoR = hEqNoR->GetBinContent(ib);

    const double htAll = hHTAll->GetBinContent(ib);
    const double htAcc = hHTAcc->GetBinContent(ib);

    const double mbAll = hMBAll->GetBinContent(ib);
    const double mbAcc = hMBAcc->GetBinContent(ib);

    nHTAcc += htAcc;

    if (eqNoR <= 0.0) continue;

    neqNoR += eqNoR;

    if (htAll <= 0.0 || mbAll <= 0.0 || htAcc <= 0.0) continue;

    const double RHT = htAcc / htAll;
    const double RMB = mbAcc / mbAll;

    if (RHT <= 0.0) continue;

    const double ratio = RMB / RHT;
    const double eqWithR = eqNoR * ratio;

    hRHT->SetBinContent(ib, RHT);
    hRMB->SetBinContent(ib, RMB);
    hRatio->SetBinContent(ib, ratio);
    hEqWithR->SetBinContent(ib, eqWithR);

    neqWithR += eqWithR;
  }

  if (normInfo) {
    normInfo->nEqNoR   = neqNoR;
    normInfo->nEqWithR = neqWithR;
    normInfo->nHTAcc   = nHTAcc;
  }

  if (outDir) {
    outDir->cd();
    hRHT->Write();
    hRMB->Write();
    hRatio->Write();
    hEqWithR->Write();
  }

  delete hRHT;
  delete hRMB;
  delete hRatio;
  delete hEqWithR;

  return neqWithR;
}

vector<EqMbNormInfo> normSummary(kCentralities.size());
vector<bool> normSummaryValid(kCentralities.size(), false);

void unfold_data_wEff(const char* dataFileHT,
                 const char* dataFileMB,
                 const char* responseFile,
                 const char* outDir,
                 int nIter = kBayesIterDefault)
{
  gStyle->SetOptStat(0);
  EnsureDir(outDir);
  TH1::SetDefaultSumw2(kTRUE);

  TFile* fData = TFile::Open(dataFileHT, "READ");
  if (!fData || fData->IsZombie()) {
    cout << "[error] Cannot open HT data file: " << dataFileHT << endl;
    return;
  }

  TFile* fDataMB = TFile::Open(dataFileMB, "READ");
  if (!fDataMB || fDataMB->IsZombie()) {
    cout << "[error] Cannot open MB data file: " << dataFileMB << endl;
    fData->Close(); delete fData;
    return;
  }

  TFile* fRespAll = TFile::Open(responseFile, "READ");
  if (!fRespAll || fRespAll->IsZombie()) {
    cout << "[error] Cannot open response file: " << responseFile << endl;
    fDataMB->Close(); delete fDataMB;
    fData->Close(); delete fData;
    return;
  }

  const string outRootPath = string(outDir) + "/unfolded_data.root";
  TFile* fOutAll = TFile::Open(outRootPath.c_str(), "RECREATE");
  if (!fOutAll || fOutAll->IsZombie()) {
    cout << "[error] Cannot create output file: " << outRootPath << endl;
    fRespAll->Close(); delete fRespAll;
    fDataMB->Close();  delete fDataMB;
    fData->Close();    delete fData;
    return;
  }

  cout << "\n=== Unfolding real data ===\n";
  cout << "  HT data file  : " << dataFileHT << "\n";
  cout << "  MB data file  : " << dataFileMB << "\n";
  cout << "  Response file : " << responseFile << "\n";
  cout << "  Output file   : " << outRootPath << "\n";
  cout << "  Bayes iters   : " << nIter << "\n\n";

  // loop R, centrality, ptlead
  for (size_t iR = 0; iR < kRadii.size(); ++iR) {
    const string R = kRadii[iR];

    // choose area cut from R
    double Rval = 0.0;
    if (sscanf(R.c_str(), "R%lf", &Rval) != 1) {
      cout << "[warn] Could not parse radius from " << R << ", skipping.\n";
      continue;
    }
    double areaMin = 0.0;
    if      (Rval < 0.25) areaMin = CUT_AREA_02;
    else if (Rval < 0.35) areaMin = CUT_AREA_03;
    else                  areaMin = CUT_AREA_04;

    for (size_t iC = 0; iC < kCentralities.size(); ++iC) {
      const string C = kCentralities[iC];

      const string treePath = R + "/" + C + "/JetTree";
      TTree* tr = dynamic_cast<TTree*>(fData->Get(treePath.c_str()));
      if (!tr) {
        cout << "[note] missing tree in data: " << treePath << " (skip)\n";
        continue;
      }

      cout << "\n>>> Tag: " << R << "  " << C << "\n";

      // enable only needed branches (RECO only!)
      tr->SetBranchStatus("*", 0);
      tr->SetBranchStatus("reco_pt_corr", 1);
      tr->SetBranchStatus("reco_pt_lead", 1);
      tr->SetBranchStatus("reco_trigger_match", 1);
      tr->SetBranchStatus("centralityWeight", 1);
      tr->SetBranchStatus("reco_area", 1);
      tr->SetBranchStatus("reco_neutral_fraction", 1);

      float reco_pt_corr=0, reco_pt_lead=0;
      Bool_t reco_trigger_match=kFALSE;
      float centralityWeight=1.0f;
      float reco_area=0.0f, reco_neutral_fraction=0.0f;

      tr->SetBranchAddress("reco_pt_corr", &reco_pt_corr);
      tr->SetBranchAddress("reco_pt_lead", &reco_pt_lead);
      tr->SetBranchAddress("reco_trigger_match", &reco_trigger_match);
      tr->SetBranchAddress("centralityWeight", &centralityWeight);
      tr->SetBranchAddress("reco_area", &reco_area);
      tr->SetBranchAddress("reco_neutral_fraction", &reco_neutral_fraction);

      const Long64_t n = tr->GetEntries();
      if (n <= 0) {
        cout << "  [note] empty tree in data, skip.\n";
        continue;
      }

      TDirectory* normDir = fOutAll->GetDirectory("normalization");
      if (!normDir) normDir = fOutAll->mkdir("normalization");

      double nEqMB_withR = -1.0;

      if (!normSummaryValid[iC]) {
        EqMbNormInfo tmpNorm;

        nEqMB_withR =
            CalcEqMbWithR(fData, fDataMB, C, &tmpNorm, normDir);

        if (nEqMB_withR > 0.0) {
          normSummary[iC] = tmpNorm;
          normSummaryValid[iC] = true;
        }
      } else {
        nEqMB_withR = normSummary[iC].nEqWithR;
      }

      if (nEqMB_withR <= 0.0) {
        cout << "  [error] Bad equivalent MB normalization for " << C
            << ", skipping centrality." << endl;
        continue;
      }

      // Store scalar normalization too. The run-by-run QA histograms are written
      // inside CalcEqMbWithR().
      normDir->cd();
      TParameter<double>(Form("N_MB_equiv_withR_%s", C.c_str()),
                         nEqMB_withR).Write("", TObject::kOverwrite);

      for (int ic = 0; ic < kNPtLeadCuts; ++ic) {
        const double cut = kPtLeadCuts[ic];
        const string tag      = R + "_" + C + Form("_ptlead%.0f", cut);
        const string tagfile  = tag; // same as used when saving responses

        cout << "  - pTlead cut >= " << cut << " GeV (tag " << tag << ")\n";

        // ---------------- measured spectrum in DATA ----------------
        TH1D* hMeasData = new TH1D(
            ("hMeasData_"+tag).c_str(),
            ";p_{T}^{reco,corr} [GeV];counts",
            nbins_meas, bin_meas_edges);
        hMeasData->SetDirectory(0);

        // event loop
        for (Long64_t i = 0; i < n; ++i) {
          if ((i % 200000) == 0)
            cout << "    filling data: " << i << "/" << n << "\r" << std::flush;

          tr->GetEntry(i);

          const double w = (double)centralityWeight;

          if (!reco_trigger_match)               continue;
          if (reco_area < areaMin)               continue;
          if (reco_neutral_fraction > CUT_NEUTRAL_FRACTION) continue;
          if (reco_pt_lead < cut)                continue;

          hMeasData->Fill(reco_pt_corr, w);
        }
        cout << "\n    data integral = " << hMeasData->Integral(0,-1) << endl;

        TH1D* hMeasData_raw = (TH1D*)hMeasData->Clone(Form("hMeasData_raw_%s", tag.c_str()));
        hMeasData_raw->SetDirectory(0);

        // The response object already contains:
        //   - trigger/reco-side inefficiency: hRespFull projectionY < hTrueFull,
        //   - fake/purity information: hMeasFull can exceed hRespFull projectionX,
        //   - matching efficiency: same truth denominator hTrueFull used in RooUnfoldResponse.
        // Do NOT apply h_trig_eff_reco, h_pur_eff_reco, or h_match_eff_truth here,
        // otherwise the data are corrected twice.

        // ---------------- load response & prior ----------------
        TDirectory* dResp =
          dynamic_cast<TDirectory*>(fRespAll->Get(tagfile.c_str()));
        if (!dResp) {
          cout << "    [ERROR] no directory '" << tagfile << "' in response file --> skipping this tag\n";
          delete hMeasData;
          delete hMeasData_raw;
          continue;
        }

        //   response      -> "response"
        //   prior (TH1D)  -> "hPrior"
        RooUnfoldResponse* response = 0;
        dResp->GetObject("response", response);
        if (!response) {
          cout << "    [ERROR] object 'response' not found in dir " << tagfile
              << " --> skipping\n";

          delete hMeasData;
          delete hMeasData_raw;
          continue;
        }

        TH1D* hPrior = 0;
        dResp->GetObject("hPrior", hPrior);
        if (!hPrior) {
          cout << "    [WARN] 'hPrior' not found in dir " << tagfile
               << " (will use default prior from response)\n";
        }

        // ---------------- run unfolding ----------------
        cout << "    running RooUnfoldBayes with " << nIter << " iterations\n";
        RooUnfoldBayes unfold(response, hMeasData, nIter,
                              /*smoothit=*/false,
                              hPrior); // may be 0, RooUnfold handles it

        TH1D* hUnfoldedRecoBins = (TH1D*)unfold.Hunfold(); // reco-binning
        hUnfoldedRecoBins->SetName(Form("UnfoldedRecoBins_%s", tag.c_str()));
        hUnfoldedRecoBins->SetDirectory(0);

        // Rebin to truth binning (this is what you'll use later)
        TH1D* hUnfoldedTruthBins_raw = (TH1D*)hUnfoldedRecoBins->Rebin(
            nbins_truth,
            Form("hUnfoldedTruthBins_raw_%s", tag.c_str()),
            bin_truth_edges);
        hUnfoldedTruthBins_raw->SetDirectory(0);

        // No external matching-efficiency correction here.
        // hUnfoldedTruthBins_raw is already corrected to the truth denominator
        // encoded in the response. It is the physics-level unfolded spectrum
        // before final invariant-yield normalization.

        // ---------------- final physics normalization ----------------
        //
        // The raw unfolded spectra are counts.  These final spectra are
        //
        //   (1/N_evt^MB) * (1/(2*pi*pT)) * d^2N_jet/(dpT deta)
        //
        // using Delta eta = 2*(1-R), because |eta_jet| < 1-R.
        //
        // Keep both unfolded counts and final invariant spectrum in the ROOT output.
        TH1D* hUnfoldedTruthBins_raw_finalInvariant =
            MakeFinalInvariantSpectrum(hUnfoldedTruthBins_raw,
                                        Form("hUnfoldedTruthBins_raw_finalInvariant_%s", tag.c_str()),
                                        nEqMB_withR, Rval, true);

        // ---------------- save to output file ----------------
        TDirectory* dOut = fOutAll->mkdir(tagfile.c_str());
        if (!dOut) dOut = fOutAll->GetDirectory(tagfile.c_str());
        dOut->cd();

        hMeasData_raw->Write("hMeasData_raw", TObject::kOverwrite);
        hMeasData->Write("hMeasData", TObject::kOverwrite);
        hUnfoldedRecoBins->Write("hUnfoldedRecoBins", TObject::kOverwrite);
        hUnfoldedTruthBins_raw->Write("hUnfoldedTruthBins", TObject::kOverwrite);
        hUnfoldedTruthBins_raw_finalInvariant->Write("hUnfoldedTruthBins_finalInvariant", TObject::kOverwrite);

        TParameter<double>("N_MB_equiv_withR", nEqMB_withR).Write("", TObject::kOverwrite);
        TParameter<double>("jetR", Rval).Write("", TObject::kOverwrite);
        TParameter<double>("deltaEta", 2.0*(1.0-Rval)).Write("", TObject::kOverwrite);

        TMatrixD cov = unfold.Eunfold(RooUnfold::kCovariance);
        cov.Write("covariance", TObject::kOverwrite);
        
        // cleanup
        delete hMeasData;
        delete hMeasData_raw;
        delete hUnfoldedRecoBins;
        delete hUnfoldedTruthBins_raw;
        if (hUnfoldedTruthBins_raw_finalInvariant) delete hUnfoldedTruthBins_raw_finalInvariant;
      } // ptlead cuts
    } // centralities
  } // radii

  cout << "\n\n";
cout << "============================================================\n";
cout << "             EQUIVALENT MINIMUM-BIAS SUMMARY\n";
cout << "============================================================\n";

for (size_t iC = 0; iC < kCentralities.size(); ++iC) {

  const string& C = kCentralities[iC];

  cout << "\nCentrality: " << C << "\n";

  if (!normSummaryValid[iC]) {
    cout << "  [no valid normalization]\n";
    continue;
  }

  const EqMbNormInfo& info = normSummary[iC];

  cout << "  HT accepted events        = "
       << info.nHTAcc << "\n";

  cout << "  N_MB_equiv no R           = "
       << info.nEqNoR << "\n";

  cout << "  N_MB_equiv with R         = "
       << info.nEqWithR << "\n";

  if (info.nHTAcc > 0.0) {
    cout << "  <HT -> MB weight>         = "
         << info.nEqNoR / info.nHTAcc << "\n";

    cout << "  total EqMB / HT accepted  = "
         << info.nEqWithR / info.nHTAcc << "\n";
  }

  if (info.nEqNoR > 0.0) {
    cout << "  acceptance correction     = "
         << info.nEqWithR / info.nEqNoR << "\n";
  }
}

cout << "\n============================================================\n";

  fOutAll->Write();
  fOutAll->Close();
  delete fOutAll;

  fRespAll->Close();
  delete fRespAll;

  fData->Close();
  delete fData;

  fDataMB->Close(); 
  delete fDataMB;

  cout << "\nAll done. Unfolded data written to: "
       << outRootPath << "\n";
}


// Backward-compatible wrapper for old steering scripts.
// The effFile argument is intentionally ignored because the response already
// contains the efficiency/purity information.
void unfold_data_wEff(const char* dataFileHT,
                 const char* dataFileMB,
                 const char* responseFile,
                 const char* effFile,
                 const char* outDir,
                 int nIter = kBayesIterDefault)
{
  cout << "[note] effFile argument is ignored: " << effFile << endl;
  cout << "       efficiencies/purity are already encoded in the response.\n";
  unfold_data_wEff(dataFileHT, dataFileMB, responseFile, outDir, nIter);
}

// ==================================================
// Helper implementations
// ==================================================

static bool SameBinning1D(const TH1* a, const TH1* b, double tol)
{
  if (!a || !b) return false;
  if (a->GetNbinsX() != b->GetNbinsX()) return false;

  const TAxis* ax = a->GetXaxis();
  const TAxis* bx = b->GetXaxis();

  for (int ib = 1; ib <= a->GetNbinsX(); ++ib) {
    if (std::fabs(ax->GetBinLowEdge(ib) - bx->GetBinLowEdge(ib)) > tol) return false;
    if (std::fabs(ax->GetBinUpEdge(ib)  - bx->GetBinUpEdge(ib))  > tol) return false;
  }
  return true;
}

static TH1D* MakeFinalInvariantSpectrum(const TH1D* src,
                                         const char* outName,
                                         double nEqMB,
                                         double jetR,
                                         bool verbose)
{
  if (!src) return 0;

  TH1D* h = (TH1D*)src->Clone(outName);
  h->SetDirectory(0);
  h->SetTitle(Form("%s;#it{p}_{T,jet} (GeV/#it{c});(1/N_{evt})(1/2#pi p_{T}) d^{2}N_{jet}/d#it{p}_{T}d#eta",
                   src->GetTitle()));

  const double deltaEta = 2.0 * (1.0 - jetR);

  if (nEqMB <= 0.0 || deltaEta <= 0.0) {
    if (verbose) {
      cout << "    [ERROR] cannot make final invariant spectrum for "
           << src->GetName() << ": N_MB_equiv=" << nEqMB
           << ", deltaEta=" << deltaEta << endl;
    }
    h->Reset();
    return h;
  }

  for (int ib = 1; ib <= h->GetNbinsX(); ++ib) {
    const double y  = h->GetBinContent(ib);
    const double ey = h->GetBinError(ib);
    const double pt = h->GetXaxis()->GetBinCenter(ib);
    const double bw = h->GetXaxis()->GetBinWidth(ib);

    if (pt <= 0.0 || bw <= 0.0) {
      h->SetBinContent(ib, 0.0);
      h->SetBinError(ib, 0.0);
      continue;
    }

    const double norm = nEqMB * bw * deltaEta * (2.0 * M_PI * pt);

    h->SetBinContent(ib, y  / norm);
    h->SetBinError(ib, ey / norm);
  }

  if (verbose) {
    cout << "    made final invariant spectrum: " << h->GetName() << endl;
    cout << "      N_MB_equiv = " << nEqMB
         << ", DeltaEta = " << deltaEta
         << ", factor includes DeltaPt and 2*pi*pT bin-by-bin" << endl;
  }

  return h;
}
