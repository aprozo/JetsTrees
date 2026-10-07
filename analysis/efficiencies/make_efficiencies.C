#include "TFile.h"
#include "TDirectory.h"
#include "TDirectoryFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TKey.h"
#include "TCollection.h"
#include "TROOT.h"
#include "TString.h"
#include "TSystem.h"

#include <iostream>
#include <string>
#include <vector>
#include <cstdio>
#include <cmath>

// =========================================================
//  Jet-quality cuts (must match analysis / maker)
// =========================================================
const double CUT_AREA_02 = 0.07;  // R = 0.2
const double CUT_AREA_03 = 0.20;  // R = 0.3
const double CUT_AREA_04 = 0.40;  // R = 0.4

// Max neutral energy fraction
const double CUT_NEUTRAL_FRACTION = 0.95;

// pT_lead thresholds (on reco-level leading particle)
const int    N_LEAD_THR = 4;
const double PTLEAD_THR[N_LEAD_THR] = {0.0, 5.0, 7.0, 9.0};

// =========================================================
//  Analysis binning (same as unfolding / embedding)
// =========================================================

// reco (measured) binning
static const int nbins_meas = 24;
static const double bin_meas_edges[nbins_meas+1] = {
  -100,-80,-60,-40,-20,-10,-5,-2.5,0,2.5,5,7.5,10,12.5,15,17.5,
  20,22.5,25,27.5,30,35,40,50,60
};

// truth (MC) binning
static const int nbins_truth = 10;
static const double bin_truth_edges[nbins_truth+1] = {
  0,5,10,15,20,25,30,35,40,50,60
};

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

void make_efficiencies(const char *infile  = "embedding_merged.root",
                       const char *outfile = "efficiencies.root")
{
    TH1::SetDefaultSumw2(kTRUE);

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

    const char* qaDir = "QA_plots";
    gSystem->mkdir(qaDir, kTRUE);

    // Loop over top-level keys (radii: R0.2, R0.3, R0.4, ...)
    TIter nextR(fin->GetListOfKeys());
    TKey *keyR = 0;

    while ((keyR = (TKey *) nextR())) {
        if (std::string(keyR->GetClassName()) != "TDirectoryFile")
            continue;

        std::string rname = keyR->GetName();
        if (rname.size() < 2 || rname[0] != 'R') continue;

        double R = 0.0;
        if (sscanf(rname.c_str(), "R%lf", &R) != 1) {
            std::cerr << "  Warning: could not parse radius from '"
                      << rname << "'. Skipping this directory." << std::endl;
            continue;
        }

        const double dRmax = 0.6 * R;

        double areaMin = 0.0;
        if      (R < 0.25) areaMin = CUT_AREA_02;
        else if (R < 0.35) areaMin = CUT_AREA_03;
        else               areaMin = CUT_AREA_04;

        TDirectoryFile *dirR = (TDirectoryFile *) keyR->ReadObj();
        if (!dirR) continue;

        std::cout << "Radius directory: " << rname
                  << " (R = " << R << ", dRmax = " << dRmax
                  << ", areaMin = " << areaMin << ")" << std::endl;

        // Loop over centrality subdirs
        TIter nextC(dirR->GetListOfKeys());
        TKey *keyC = 0;

        while ((keyC = (TKey *) nextC())) {
            if (std::string(keyC->GetClassName()) != "TDirectoryFile")
                continue;

            std::string cname = keyC->GetName();
            TDirectoryFile *dirC = (TDirectoryFile *) keyC->ReadObj();
            if (!dirC) continue;

            TTree *tree = (TTree *) dirC->Get("JetTree");
            if (!tree) {
                std::cerr << "  Warning: JetTree not found in " << rname
                          << "/" << cname << std::endl;
                continue;
            }

            Long64_t nentries = tree->GetEntries();
            std::cout << "  Centrality: " << cname << " (entries: "
                      << nentries << ")" << std::endl;

            // --- Check that required branches exist ---
            if (!tree->GetBranch("mc_pt") ||
                !tree->GetBranch("mc_pt_lead") ||
                !tree->GetBranch("deltaR") ||
                !tree->GetBranch("reco_pt") ||
                !tree->GetBranch("reco_pt_corr") ||
                !tree->GetBranch("reco_pt_lead") ||
                !tree->GetBranch("reco_area") ||
                !tree->GetBranch("reco_neutral_fraction") ||
                !tree->GetBranch("reco_trigger_match") ||
                !tree->GetBranch("centralityWeight") ||
                !tree->GetBranch("xsecWeight")) {
                std::cerr << "    Missing some required branches in "
                          << rname << "/" << cname
                          << " – skipping this centrality." << std::endl;
                continue;
            }

            // --- Set branch addresses ---
            Float_t mc_pt          = 0.f;
            Float_t mc_pt_lead     = 0.f;
            Float_t deltaR         = 0.f;
            Float_t reco_pt        = 0.f;
            Float_t reco_pt_corr   = 0.f;
            Float_t reco_pt_lead   = 0.f;
            Float_t reco_area      = 0.f;
            Float_t reco_neutral_fraction = 0.f;
            Bool_t  reco_trigger_match    = kFALSE;
            Float_t xsecWeight      = 1.f;
            Float_t centralityWeight = 1.f;

            tree->SetBranchAddress("mc_pt",        &mc_pt);
            tree->SetBranchAddress("mc_pt_lead",   &mc_pt_lead);
            tree->SetBranchAddress("deltaR",       &deltaR);
            tree->SetBranchAddress("reco_pt",      &reco_pt);
            tree->SetBranchAddress("reco_pt_corr", &reco_pt_corr);
            tree->SetBranchAddress("reco_pt_lead", &reco_pt_lead);
            tree->SetBranchAddress("reco_area",    &reco_area);
            tree->SetBranchAddress("reco_neutral_fraction", &reco_neutral_fraction);
            tree->SetBranchAddress("reco_trigger_match",    &reco_trigger_match);
            tree->SetBranchAddress("xsecWeight",       &xsecWeight);
            tree->SetBranchAddress("centralityWeight", &centralityWeight);

            // Arrays per pT_lead threshold
            TH1D *h_match_mc_den [N_LEAD_THR];
            TH1D *h_match_mc_num [N_LEAD_THR];
            TH1D *h_match_mc_eff [N_LEAD_THR];

            TH1D *h_trig_den [N_LEAD_THR];
            TH1D *h_trig_num [N_LEAD_THR];
            TH1D *h_trig_eff [N_LEAD_THR];

            TH1D *h_pur_den  [N_LEAD_THR];
            TH1D *h_pur_num  [N_LEAD_THR];
            TH1D *h_pur_eff  [N_LEAD_THR];

            TH1D* h_trig_den_pthat[N_LEAD_THR][kNPthatBins];
            TH1D* h_trig_num_pthat[N_LEAD_THR][kNPthatBins];
            TH1D* h_pur_den_pthat [N_LEAD_THR][kNPthatBins];
            TH1D* h_pur_num_pthat [N_LEAD_THR][kNPthatBins];

            TDirectory* outTag[N_LEAD_THR] = {0};

            // --- Book histograms for each pT_lead threshold ---
            for (int it = 0; it < N_LEAD_THR; ++it) {
                double thr = PTLEAD_THR[it];

                // Build tag name consistent with unfolding:
                // e.g. "R0.2_CENT_0_10_ptlead5"
                TString tagDirName;
                tagDirName.Form("%s_%s_ptlead%.0f", rname.c_str(), cname.c_str(), thr);

                // Create / get directory for this tag at the TOP level of the file
                fout->cd();
                outTag[it] = (TDirectory*)fout->mkdir(tagDirName);
                if (!outTag[it]) outTag[it] = fout->GetDirectory(tagDirName);
                outTag[it]->cd();

                // Histogram base names can now be simple, they live in separate dirs
                // Matching efficiency vs mc_pt
                h_match_mc_den[it] = new TH1D("h_match_mc_den", "MC jets (denominator)",
                                            nbins_truth, bin_truth_edges);
                h_match_mc_num[it] = new TH1D("h_match_mc_num", "Matched MC jets (numerator)",
                                            nbins_truth, bin_truth_edges);
                h_match_mc_eff[it] = (TH1D*)h_match_mc_num[it]->Clone("h_match_mc_eff");
                h_match_mc_eff[it]->SetTitle("");

                // Trigger efficiency vs reco_pt_corr
                h_trig_den[it] = new TH1D("h_trig_den", "Reco jets (denominator)",
                                        nbins_meas, bin_meas_edges);
                h_trig_num[it] = new TH1D("h_trig_num", "Triggered reco jets (numerator)",
                                        nbins_meas, bin_meas_edges);
                h_trig_eff[it] = (TH1D*)h_trig_num[it]->Clone("h_trig_eff");
                h_trig_eff[it]->SetTitle("");

                // Purity vs reco_pt_corr
                h_pur_den[it] = new TH1D("h_pur_den", "All reco jets (denominator)",
                                        nbins_meas, bin_meas_edges);
                h_pur_num[it] = new TH1D("h_pur_num", "Matched reco jets (numerator)",
                                        nbins_meas, bin_meas_edges);
                h_pur_eff[it] = (TH1D*)h_pur_num[it]->Clone("h_pur_eff");
                h_pur_eff[it]->SetTitle("");

                // Per-pThat trigger/purity histograms (centrality-weight only)
                for (int ip = 0; ip < kNPthatBins; ++ip) {
                    h_trig_den_pthat[it][ip] = new TH1D(
                        Form("h_trig_den_pthat%d", ip),
                        Form("Reco jets denominator, pThat bin %d", ip),
                        nbins_meas, bin_meas_edges);
                    h_trig_den_pthat[it][ip]->SetDirectory(0);

                    h_trig_num_pthat[it][ip] = new TH1D(
                        Form("h_trig_num_pthat%d", ip),
                        Form("Triggered reco jets numerator, pThat bin %d", ip),
                        nbins_meas, bin_meas_edges);
                    h_trig_num_pthat[it][ip]->SetDirectory(0);

                    h_pur_den_pthat[it][ip] = new TH1D(
                        Form("h_pur_den_pthat%d", ip),
                        Form("Reco jets purity denominator, pThat bin %d", ip),
                        nbins_meas, bin_meas_edges);
                    h_pur_den_pthat[it][ip]->SetDirectory(0);

                    h_pur_num_pthat[it][ip] = new TH1D(
                        Form("h_pur_num_pthat%d", ip),
                        Form("Matched reco jets purity numerator, pThat bin %d", ip),
                        nbins_meas, bin_meas_edges);
                    h_pur_num_pthat[it][ip]->SetDirectory(0);
                }

            }

            // --- Loop over tree entries and fill histos ---
            for (Long64_t i = 0; i < nentries; ++i) {
                tree->GetEntry(i);

                const int ip = FindPtHatBin((double)xsecWeight);
                if (ip < 0) continue;

                double w = (double)xsecWeight / kNgenEvents[ip] * (double)centralityWeight;
                double wCent = (double)centralityWeight;

                // Reco quality cuts
                bool haveReco = (reco_pt > -500.0);
                bool haveMC = (mc_pt > 0.0);

                bool passRecoCuts = false;
                if (haveReco) {
                    passRecoCuts = (reco_area >= areaMin &&
                                    reco_neutral_fraction <= CUT_NEUTRAL_FRACTION);
                }


                // -- Loop over pT_lead thresholds (on reco level, nested) --
                for (int it = 0; it < N_LEAD_THR; ++it) {
                    double thr = PTLEAD_THR[it];

                    // ---------------- Matching efficiency (MC-based) ----------------
                    // Denominator: all MC jets
                    if (haveMC) {
                        h_match_mc_den[it]->Fill(mc_pt, w);

                        // Numerator: MC jets that have a matched reco jet
                        // passing area/NEF cuts and with reco_pt_lead >= thr
                        bool isMatched = haveReco &&
                                         (deltaR > 0.0) &&
                                         (deltaR < dRmax) &&
                                         passRecoCuts &&
                                         (reco_pt_lead >= thr);

                        if (isMatched) {
                            h_match_mc_num[it]->Fill(mc_pt, w);
                        }
                    }

                    // ---------------- Trigger efficiency & purity (reco-based) -----
                    if (haveReco && passRecoCuts && reco_pt_lead >= thr) {
                        // Trigger efficiency: fill per-pThat with centrality weight only
                        h_trig_den_pthat[it][ip]->Fill(reco_pt_corr, wCent);
                        if (reco_trigger_match) {
                            h_trig_num_pthat[it][ip]->Fill(reco_pt_corr, wCent);
                        }

                        // Purity: fill per-pThat with centrality weight only
                        bool isMatchedReco = haveMC &&
                                             (deltaR > 0.0) &&
                                             (deltaR < dRmax);

                        h_pur_den_pthat[it][ip]->Fill(reco_pt_corr, wCent);
                        if (isMatchedReco) {
                            h_pur_num_pthat[it][ip]->Fill(reco_pt_corr, wCent);
                        }
                    }

                } // end loop over thresholds
            } // end loop over entries


            // ------------------------------------------------------------
            // Merge trigger/purity from per-pThat histograms with reco-bin mask
            // ------------------------------------------------------------
            for (int it = 0; it < N_LEAD_THR; ++it) {
                h_trig_den[it]->Reset();
                h_trig_num[it]->Reset();
                h_pur_den[it]->Reset();
                h_pur_num[it]->Reset();

                for (int ip = 0; ip < kNPthatBins; ++ip) {
                    const double xw = kXsecWeights[ip] / kNgenEvents[ip];

                    for (int ix = 1; ix <= nbins_meas; ++ix) {
                        const double c = h_trig_den_pthat[it][ip]->GetBinContent(ix);
                        const double e = h_trig_den_pthat[it][ip]->GetBinError(ix);
                        const double signif = (e > 0.0 ? c / e : 0.0);
                        const bool keep = (signif > kMinSignif);

                        if (!keep) continue;

                        // trigger denominator
                        {
                            const double oldC = h_trig_den[it]->GetBinContent(ix);
                            const double oldE = h_trig_den[it]->GetBinError(ix);
                            const double addC = xw * h_trig_den_pthat[it][ip]->GetBinContent(ix);
                            const double addE = xw * h_trig_den_pthat[it][ip]->GetBinError(ix);

                            h_trig_den[it]->SetBinContent(ix, oldC + addC);
                            h_trig_den[it]->SetBinError(ix, std::sqrt(oldE*oldE + addE*addE));
                        }

                        // trigger numerator
                        {
                            const double oldC = h_trig_num[it]->GetBinContent(ix);
                            const double oldE = h_trig_num[it]->GetBinError(ix);
                            const double addC = xw * h_trig_num_pthat[it][ip]->GetBinContent(ix);
                            const double addE = xw * h_trig_num_pthat[it][ip]->GetBinError(ix);

                            h_trig_num[it]->SetBinContent(ix, oldC + addC);
                            h_trig_num[it]->SetBinError(ix, std::sqrt(oldE*oldE + addE*addE));
                        }

                        // purity denominator
                        {
                            const double oldC = h_pur_den[it]->GetBinContent(ix);
                            const double oldE = h_pur_den[it]->GetBinError(ix);
                            const double addC = xw * h_pur_den_pthat[it][ip]->GetBinContent(ix);
                            const double addE = xw * h_pur_den_pthat[it][ip]->GetBinError(ix);

                            h_pur_den[it]->SetBinContent(ix, oldC + addC);
                            h_pur_den[it]->SetBinError(ix, std::sqrt(oldE*oldE + addE*addE));
                        }

                        // purity numerator
                        {
                            const double oldC = h_pur_num[it]->GetBinContent(ix);
                            const double oldE = h_pur_num[it]->GetBinError(ix);
                            const double addC = xw * h_pur_num_pthat[it][ip]->GetBinContent(ix);
                            const double addE = xw * h_pur_num_pthat[it][ip]->GetBinError(ix);

                            h_pur_num[it]->SetBinContent(ix, oldC + addC);
                            h_pur_num[it]->SetBinError(ix, std::sqrt(oldE*oldE + addE*addE));
                        }
                    }
                }
            }

            // --- Build efficiency histograms with binomial errors ---
            for (int it = 0; it < N_LEAD_THR; ++it) {
                if (!outTag[it]) continue;   // safety

                outTag[it]->cd();

                double thr = PTLEAD_THR[it];

                TString tagDirName;
                tagDirName.Form("%s_%s_ptlead%.0f", rname.c_str(), cname.c_str(), thr);

                // ------------------------------------------------------------------
                // 1) (Optional) fine-binned efficiencies for writing to the ROOT file
                // ------------------------------------------------------------------
                h_match_mc_eff[it]->Divide(h_match_mc_num[it],
                                        h_match_mc_den[it],
                                        1.0, 1.0, "b");

                h_trig_eff[it]->Divide(h_trig_num[it],
                                    h_trig_den[it],
                                    1.0, 1.0);

                h_pur_eff[it]->Divide(h_pur_num[it],
                                    h_pur_den[it],
                                    1.0, 1.0);

                // -----------------------------------------------------------
                // 2) REBIN numerator & denominator, then compute efficiency
                // -----------------------------------------------------------

                // --- Matching efficiency in TRUTH binning ---
                TH1D *h_match_den_truth = (TH1D*)h_match_mc_den[it]->Clone("h_match_den_truth");
                h_match_den_truth->SetDirectory(0);

                TH1D *h_match_num_truth = (TH1D*)h_match_mc_num[it]->Clone("h_match_num_truth");
                h_match_num_truth->SetDirectory(0);

                TH1D *h_match_eff_truth = (TH1D*)h_match_num_truth->Clone("h_match_eff_truth");
                h_match_eff_truth->SetDirectory(0);
                h_match_eff_truth->Divide(h_match_num_truth, h_match_den_truth, 1.0, 1.0, "b");

                // now ONLY one copy will be written
                outTag[it]->cd();
                h_match_den_truth->Write("h_match_den_truth", TObject::kOverwrite);
                h_match_num_truth->Write("h_match_num_truth", TObject::kOverwrite);
                h_match_eff_truth->Write("h_match_eff_truth", TObject::kOverwrite);

                // QA plot
                TCanvas* cMatch = new TCanvas(
                    Form("cMatch_%s", tagDirName.Data()),
                    "Matching efficiency (truth binning)", 800, 600);

                h_match_eff_truth->SetLineColor(kBlue+1);
                h_match_eff_truth->SetMarkerColor(kBlue+1);
                h_match_eff_truth->SetMarkerStyle(20);
                h_match_eff_truth->GetXaxis()->SetTitle("#it{p}_{T}^{MC} (GeV)");
                h_match_eff_truth->GetYaxis()->SetTitle("Matching efficiency");
                h_match_eff_truth->SetStats(0);
                h_match_eff_truth->Draw("E1");

                cMatch->SaveAs(Form("QA_plots/QA_match_eff_%s.pdf", tagDirName.Data()));
                cMatch->SaveAs(Form("QA_plots/QA_match_eff_%s.png", tagDirName.Data()));
                delete cMatch;

                // --- Trigger efficiency in MEASURED binning ---
                TH1D *h_trig_den_reco = (TH1D*)h_trig_den[it]->Clone("h_trig_den_reco");
                h_trig_den_reco->SetDirectory(0);

                TH1D *h_trig_num_reco = (TH1D*)h_trig_num[it]->Clone("h_trig_num_reco");
                h_trig_num_reco->SetDirectory(0);

                TH1D *h_trig_eff_reco = (TH1D*)h_trig_num_reco->Clone("h_trig_eff_reco");
                h_trig_eff_reco->SetDirectory(0);
                h_trig_eff_reco->SetTitle("");
                h_trig_eff_reco->Divide(h_trig_num_reco, h_trig_den_reco, 1.0, 1.0);

                // write single copy
                outTag[it]->cd();
                h_trig_den_reco->Write("h_trig_den_reco", TObject::kOverwrite);
                h_trig_num_reco->Write("h_trig_num_reco", TObject::kOverwrite);
                h_trig_eff_reco->Write("h_trig_eff_reco", TObject::kOverwrite);

                // QA plot
                TCanvas* cTrig = new TCanvas(
                    Form("cTrig_%s", tagDirName.Data()),
                    "Trigger efficiency (measured binning)", 800, 600);

                h_trig_eff_reco->SetLineColor(kRed+1);
                h_trig_eff_reco->SetMarkerColor(kRed+1);
                h_trig_eff_reco->SetMarkerStyle(20);
                h_trig_eff_reco->GetXaxis()->SetTitle("#it{p}_{T}^{reco,corr} (GeV)");
                h_trig_eff_reco->GetYaxis()->SetTitle("Trigger efficiency");
                h_trig_eff_reco->GetXaxis()->SetRangeUser(-20, 60);
                h_trig_eff_reco->SetStats(0);
                h_trig_eff_reco->Draw("E1");

                cTrig->SaveAs(Form("QA_plots/QA_trig_eff_%s.pdf", tagDirName.Data()));
                cTrig->SaveAs(Form("QA_plots/QA_trig_eff_%s.png", tagDirName.Data()));
                delete cTrig;

                // --- Purity in MEASURED binning ---
                TH1D *h_pur_den_reco  = (TH1D*)h_pur_den[it]->Clone("h_pur_den_reco");
                h_pur_den_reco->SetDirectory(0);

                TH1D *h_pur_num_reco  = (TH1D*)h_pur_num[it]->Clone("h_pur_num_reco");
                h_pur_num_reco->SetDirectory(0);

                TH1D *h_pur_eff_reco = (TH1D*)h_pur_num_reco->Clone("h_pur_eff_reco");
                h_pur_eff_reco->SetDirectory(0);
                h_pur_eff_reco->SetTitle("");
                h_pur_eff_reco->Divide(h_pur_num_reco, h_pur_den_reco, 1.0, 1.0);

                // write single copy
                outTag[it]->cd();
                h_pur_den_reco->Write("h_pur_den_reco", TObject::kOverwrite);
                h_pur_num_reco->Write("h_pur_num_reco", TObject::kOverwrite);
                h_pur_eff_reco->Write("h_pur_eff_reco", TObject::kOverwrite);
                
                // QA plot
                TCanvas* cPur = new TCanvas(
                    Form("cPur_%s", tagDirName.Data()),
                    "Purity (measured binning)", 800, 600);

                h_pur_eff_reco->SetLineColor(kGreen+2);
                h_pur_eff_reco->SetMarkerColor(kGreen+2);
                h_pur_eff_reco->SetMarkerStyle(20);
                h_pur_eff_reco->GetXaxis()->SetTitle("#it{p}_{T}^{reco,corr} (GeV)");
                h_pur_eff_reco->GetYaxis()->SetTitle("Purity");
                h_pur_eff_reco->GetXaxis()->SetRangeUser(-20, 60);
                h_pur_eff_reco->SetStats(0);
                h_pur_eff_reco->Draw("E1");

                cPur->SaveAs(Form("QA_plots/QA_purity_%s.pdf", tagDirName.Data()));
                cPur->SaveAs(Form("QA_plots/QA_purity_%s.png", tagDirName.Data()));
                delete cPur;
            }
            for (int it = 0; it < N_LEAD_THR; ++it) {
                for (int ip = 0; ip < kNPthatBins; ++ip) {
                    delete h_trig_den_pthat[it][ip];
                    delete h_trig_num_pthat[it][ip];
                    delete h_pur_den_pthat[it][ip];
                    delete h_pur_num_pthat[it][ip];
                }
            }
            for (int it = 0; it < N_LEAD_THR; ++it) {
                delete h_match_mc_den[it];
                delete h_match_mc_num[it];
                delete h_match_mc_eff[it];

                delete h_trig_den[it];
                delete h_trig_num[it];
                delete h_trig_eff[it];

                delete h_pur_den[it];
                delete h_pur_num[it];
                delete h_pur_eff[it];
            }
        } // end loop over centrality dirs
    } // end loop over R dirs

    fout->Write();
    fout->Close();
    fin->Close();

    std::cout << "Done. Efficiencies (matching, trigger, purity, "
              << "with reco pTlead thresholds) written to "
              << outfile << std::endl;
}
