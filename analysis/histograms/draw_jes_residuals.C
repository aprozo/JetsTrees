#include "TFile.h"
#include "TDirectory.h"
#include "TDirectoryFile.h"
#include "TKey.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLatex.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"

#include <iostream>
#include <string>
#include <cstdio>

static const int N_LEAD = 4;
static const double PTLEAD_THR[N_LEAD] = {0.0, 5.0, 7.0, 9.0};

static const int nbins_common = 10;
static const double bin_common_edges[nbins_common+1] = {
  0,5,10,15,20,25,30,35,40,50,60
};

static const int kFirstTruthBinToDraw = 1;
static const double kFitMin = -0.8;
static const double kFitMax =  0.5;

TString CentralityLabelResiduals(const std::string& cname)
{
  TString label = cname.c_str();
  label.ReplaceAll("CENT_", "");
  label.ReplaceAll("MID_",  "");
  label.ReplaceAll("PERI_", "");
  label.ReplaceAll("_", "-");
  label += " %";
  return label;
}

void draw_jes_residuals(const char* infile  = "jes_jer_hists.root",
                        const char* plotTop = "jes_jer_plots")
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TFile* fin = TFile::Open(infile, "READ");
  if (!fin || fin->IsZombie()) {
    std::cerr << "Error: cannot open input file " << infile << std::endl;
    return;
  }

  gSystem->mkdir(plotTop, kTRUE);

  TIter nextR(fin->GetListOfKeys());
  TKey* keyR = 0;

  while ((keyR = (TKey*)nextR())) {
    if (std::string(keyR->GetClassName()) != "TDirectoryFile") continue;

    const std::string rname = keyR->GetName();
    if (rname.empty() || rname[0] != 'R') continue;

    double R = 0.0;
    if (sscanf(rname.c_str(), "R%lf", &R) != 1) continue;

    TDirectoryFile* dirR = (TDirectoryFile*)keyR->ReadObj();
    if (!dirR) continue;

    TString plotRDir = Form("%s/%s", plotTop, rname.c_str());
    gSystem->mkdir(plotRDir.Data(), kTRUE);

    TIter nextC(dirR->GetListOfKeys());
    TKey* keyC = 0;

    while ((keyC = (TKey*)nextC())) {
      if (std::string(keyC->GetClassName()) != "TDirectoryFile") continue;

      const std::string cname = keyC->GetName();
      TDirectoryFile* dirC = (TDirectoryFile*)keyC->ReadObj();
      if (!dirC) continue;

      TString plotCDir = Form("%s/%s", plotRDir.Data(), cname.c_str());
      gSystem->mkdir(plotCDir.Data(), kTRUE);

      const TString centLabel = CentralityLabelResiduals(cname);

      for (int it = 0; it < N_LEAD; ++it) {
        TCanvas* cResidual = new TCanvas(
          Form("cResidual_ptlead%.0f_%s_%s",
               PTLEAD_THR[it], rname.c_str(), cname.c_str()),
          "", 1500, 1350
        );
        cResidual->Divide(3,3,0.001,0.001);

        int iPad = 1;
        int nFound = 0;

        for (int iTruthBin = kFirstTruthBinToDraw;
             iTruthBin < nbins_common;
             ++iTruthBin) {

          TH1D* hRes = (TH1D*)dirC->Get(
            Form("hResidual_ptlead%.0f_truthbin%d_%s_%s",
                 PTLEAD_THR[it],
                 iTruthBin,
                 rname.c_str(),
                 cname.c_str())
          );

          cResidual->cd(iPad);
          gPad->SetLeftMargin(0.15);
          gPad->SetRightMargin(0.04);
          gPad->SetBottomMargin(0.14);
          gPad->SetTopMargin(0.10);
          gPad->SetTicks(1,1);

          if (!hRes) {
            TLatex missing;
            missing.SetNDC();
            missing.SetTextAlign(22);
            missing.SetTextSize(0.050);
            missing.DrawLatex(0.50,0.55,"Histogram not found");
            missing.SetTextSize(0.040);
            missing.DrawLatex(
              0.50,0.45,
              Form("%.0f-%.0f GeV/#it{c}",
                   bin_common_edges[iTruthBin],
                   bin_common_edges[iTruthBin+1])
            );
            ++iPad;
            continue;
          }

          ++nFound;

          hRes->SetDirectory(0);
          hRes->SetTitle(
            Form("#it{p}_{T,true}: [%.0f, %.0f) GeV/#it{c}",
                 bin_common_edges[iTruthBin],
                 bin_common_edges[iTruthBin+1])
          );
          hRes->SetMarkerStyle(20);
          hRes->SetMarkerSize(0.65);
          hRes->SetMarkerColor(kBlack);
          hRes->SetLineColor(kBlack);

          hRes->GetXaxis()->SetTitle(
            "(#it{p}_{T,corr}-#it{p}_{T,true})/#it{p}_{T,true}"
          );
          hRes->GetYaxis()->SetTitle("weighted counts");
          hRes->GetXaxis()->SetRangeUser(-1.0,1.0);
          hRes->GetXaxis()->SetTitleSize(0.045);
          hRes->GetYaxis()->SetTitleSize(0.045);
          hRes->GetXaxis()->SetLabelSize(0.040);
          hRes->GetYaxis()->SetLabelSize(0.040);
          hRes->GetYaxis()->SetTitleOffset(1.45);

          hRes->Draw("E1");

          if (hRes->GetEntries() > 10.0 && hRes->Integral() > 0.0) {
            TF1* fit = new TF1(
              Form("fitResidual_ptlead%.0f_bin%d_%s_%s",
                   PTLEAD_THR[it],
                   iTruthBin,
                   rname.c_str(),
                   cname.c_str()),
              "gaus",
              kFitMin,
              kFitMax
            );

            fit->SetLineColor(kRed);
            fit->SetLineWidth(2);

            const int fitStatus = hRes->Fit(fit, "Q0R");

            if (fitStatus == 0) {
              fit->Draw("SAME");

              TLatex fitText;
              fitText.SetNDC();
              fitText.SetTextFont(42);
              fitText.SetTextSize(0.040);
              fitText.DrawLatex(
                0.18,0.84,
                Form("#mu = %.3f #pm %.3f",
                     fit->GetParameter(1),
                     fit->GetParError(1))
              );
              fitText.DrawLatex(
                0.18,0.78,
                Form("#sigma = %.3f #pm %.3f",
                     TMath::Abs(fit->GetParameter(2)),
                     fit->GetParError(2))
              );
            }
          }

          ++iPad;
        }

        if (nFound > 0) {
          cResidual->cd();
          TLatex head;
          head.SetNDC();
          head.SetTextFont(42);
          head.SetTextAlign(22);
          head.SetTextSize(0.020);
          head.DrawLatex(
            0.50,0.995,
            Form("STAR embedding Au+Au, #sqrt{#it{s}_{NN}} = 200 GeV, R = %.1f, %s, #it{p}_{T}^{lead} #geq %.0f GeV/#it{c}",
                 R, centLabel.Data(), PTLEAD_THR[it])
          );

          TString outputName = Form(
            "%s/residuals_ptlead%.0f_%s_%s.png",
            plotCDir.Data(),
            PTLEAD_THR[it],
            rname.c_str(),
            cname.c_str()
          );

          cResidual->SaveAs(outputName);
          std::cout << "Saved " << outputName << std::endl;
        }

        delete cResidual;
      }
    }
  }

  fin->Close();
  std::cout << "Residual plots written under: " << plotTop << std::endl;
}