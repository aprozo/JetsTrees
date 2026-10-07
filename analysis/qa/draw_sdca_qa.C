#include "TFile.h"
#include "TList.h"
#include "TKey.h"
#include "TObject.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdio>

using std::cout;
using std::endl;


// ============================================================================
// Look for an object either directly in the file or inside the maker TList.
// This makes the macro tolerant of slightly different merged-file structures.
// ============================================================================
TObject* FindQaObject(TFile* file, const char* name)
{
    if (!file) return 0;

    // First try directly in the file
    TObject* obj = file->Get(name);
    if (obj) return obj;

    // Known/likely maker-list names
    const char* listNames[] = {
        "StPicoHFJetMaker",
        "stPicoHFJetMaker",
        "StPicoJetMaker",
        "stPicoJetMaker"
    };

    const int nListNames =
        sizeof(listNames) / sizeof(listNames[0]);

    int i = 0;
    for (i = 0; i < nListNames; ++i) {

        TList* list =
            dynamic_cast<TList*>(file->Get(listNames[i]));

        if (!list) continue;

        obj = list->FindObject(name);

        if (obj) return obj;
    }


    // If the name of the list differs, scan all top-level keys
    TIter next(file->GetListOfKeys());
    TKey* key = 0;

    while ((key = (TKey*)next())) {

        TObject* topObj = file->Get(key->GetName());

        if (!topObj) continue;

        TList* list = dynamic_cast<TList*>(topObj);

        if (!list) continue;

        obj = list->FindObject(name);

        if (obj) return obj;
    }

    return 0;
}


// ============================================================================
// Mirror histogram around x = 0.
// Assumes symmetric binning, as used for the sDCA histograms.
// ============================================================================
TH1D* MakeMirroredHistogram(TH1D* hIn, const char* newName)
{
    if (!hIn) return 0;

    TH1D* hOut =
        (TH1D*)hIn->Clone(newName);

    hOut->Reset();

    const int nBins = hIn->GetNbinsX();

    int ib = 0;

    for (ib = 1; ib <= nBins; ++ib) {

        const int mirrorBin = nBins + 1 - ib;

        hOut->SetBinContent(
            ib,
            hIn->GetBinContent(mirrorBin));

        hOut->SetBinError(
            ib,
            hIn->GetBinError(mirrorBin));
    }

    return hOut;
}


// ============================================================================
// Normalize to unit integral.
// ============================================================================
void NormalizeHistogram(TH1D* h)
{
    if (!h) return;

    const double integral =
        h->Integral(1, h->GetNbinsX());

    if (integral > 0.0)
        h->Scale(1.0 / integral);
}

// ============================================================================
// Main
// ============================================================================
void draw_sdca_qa(
    const char* inputFile = "../../trees/data_merged.root",
    const char* outDir    = "sdca_qa_plots")
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    gSystem->mkdir(outDir, true);


    // ------------------------------------------------------------------------
    // Open input
    // ------------------------------------------------------------------------
    TFile* file =
        TFile::Open(inputFile, "READ");

    if (!file || file->IsZombie()) {

        cout << "ERROR: cannot open input file:" << endl;
        cout << "  " << inputFile << endl;

        return;
    }

    cout << endl;
    cout << "Input file: " << inputFile << endl;


    // ------------------------------------------------------------------------
    // pT bins
    // ------------------------------------------------------------------------
    const int nPtBins = 6;

    const char* ptLabels[nPtBins] = {
        "0p2_1",
        "1_4",
        "4_10",
        "10_15",
        "15_20",
        "20_30"
    };

    const char* ptTitles[nPtBins] = {
        "0.2 < #it{p}_{T} < 1 GeV/#it{c}",
        "1 < #it{p}_{T} < 4 GeV/#it{c}",
        "4 < #it{p}_{T} < 10 GeV/#it{c}",
        "10 < #it{p}_{T} < 15 GeV/#it{c}",
        "15 < #it{p}_{T} < 20 GeV/#it{c}",
        "20 < #it{p}_{T} < 30 GeV/#it{c}"
    };


    // ------------------------------------------------------------------------
    // Retrieve 1D sDCA distributions
    // ------------------------------------------------------------------------
    TH1D* hPos[nPtBins];
    TH1D* hNeg[nPtBins];

    int ipt = 0;

    for (ipt = 0; ipt < nPtBins; ++ipt) {

        hPos[ipt] =
            dynamic_cast<TH1D*>(
                FindQaObject(
                    file,
                    Form("hSDcaPos_%s", ptLabels[ipt])));

        hNeg[ipt] =
            dynamic_cast<TH1D*>(
                FindQaObject(
                    file,
                    Form("hSDcaNeg_%s", ptLabels[ipt])));

        if (!hPos[ipt]) {

            cout << "ERROR: missing histogram "
                 << "hSDcaPos_" << ptLabels[ipt]
                 << endl;

            file->Close();
            return;
        }

        if (!hNeg[ipt]) {

            cout << "ERROR: missing histogram "
                 << "hSDcaNeg_" << ptLabels[ipt]
                 << endl;

            file->Close();
            return;
        }
    }


    // ========================================================================
    //
    // 1. Positive vs mirrored-negative sDCA
    //
    // ========================================================================

    TCanvas* cComparison =
        new TCanvas(
            "cComparison",
            "sDCA charge comparison",
            1200,
            1000);

    cComparison->Divide(2, 3);


    TH1D* hPosNorm[nPtBins];
    TH1D* hNegMirror[nPtBins];

    for (ipt = 0; ipt < nPtBins; ++ipt) {

        cComparison->cd(ipt + 1);

        gPad->SetLogy();
        gPad->SetTicks(1, 1);

        hPosNorm[ipt] =
            (TH1D*)hPos[ipt]->Clone(
                Form("hPosNorm_%d", ipt));

        hNegMirror[ipt] =
            (TH1D*)hNeg[ipt]->Clone(
                Form("hNegMirror_%d", ipt));

        hNegMirror[ipt]->Reset();


        const int nBinsMirror =
            hNeg[ipt]->GetNbinsX();

        int ibMirror = 0;

        for (ibMirror = 1;
            ibMirror <= nBinsMirror;
            ++ibMirror) {

            const double x =
                hNeg[ipt]->GetXaxis()->GetBinCenter(ibMirror);

            const int sourceBin =
                hNeg[ipt]->GetXaxis()->FindBin(-x);


            hNegMirror[ipt]->SetBinContent(
                ibMirror,
                hNeg[ipt]->GetBinContent(sourceBin));

            hNegMirror[ipt]->SetBinError(
                ibMirror,
                hNeg[ipt]->GetBinError(sourceBin));
        }

        printf(
        "MIRROR %-8s : neg=%12.0f  mirror=%12.0f  maxNeg=%g  maxMirror=%g\n",
        ptLabels[ipt],
        hNeg[ipt]->Integral(1, hNeg[ipt]->GetNbinsX()),
        hNegMirror[ipt]->Integral(1, hNegMirror[ipt]->GetNbinsX()),
        hNeg[ipt]->GetMaximum(),
        hNegMirror[ipt]->GetMaximum());

        NormalizeHistogram(hPosNorm[ipt]);
        NormalizeHistogram(hNegMirror[ipt]);


        // Restrict visual range to the region relevant for the proposed cut.
        // Original histograms still retain the full -5...5 cm range.
        hPosNorm[ipt]->GetXaxis()->SetRangeUser(
            -2.0,
            2.0);

        hPosNorm[ipt]->GetXaxis()->SetTitle(
            "sDCA (cm)");

        hPosNorm[ipt]->GetYaxis()->SetTitle(
            "Normalized counts");

        hPosNorm[ipt]->SetLineColor(kBlue + 1);
        hPosNorm[ipt]->SetMarkerColor(kBlue + 1);
        hPosNorm[ipt]->SetMarkerStyle(20);
        hPosNorm[ipt]->SetMarkerSize(0.6);

        hNegMirror[ipt]->SetLineColor(kRed + 1);
        hNegMirror[ipt]->SetMarkerColor(kRed + 1);
        hNegMirror[ipt]->SetMarkerStyle(24);
        hNegMirror[ipt]->SetMarkerSize(0.6);


        // Find a useful minimum for log scale
        hPosNorm[ipt]->SetMinimum(1e-5);

        double maxPos =
            hPosNorm[ipt]->GetMaximum();

        double maxNeg =
            hNegMirror[ipt]->GetMaximum();

        double maxY = maxPos;

        if (maxNeg > maxY) {
            maxY = maxNeg;
        }

        hPosNorm[ipt]->SetMaximum(
            maxY * 3.0);


        hPosNorm[ipt]->Draw("HIST");
        hNegMirror[ipt]->Draw("HIST SAME");


        // Proposed +/- 0.5 cm boundary
        TLine* linePos =
            new TLine(
                0.5,
                1e-5,
                0.5,
                maxY * 3.0);

        linePos->SetLineStyle(2);
        linePos->SetLineWidth(2);
        linePos->Draw();


        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.045);

        latex.DrawLatex(
            0.14,
            0.86,
            ptTitles[ipt]);


        if (ipt == 0) {

            TLegend* leg =
                new TLegend(
                    0.52,
                    0.69,
                    0.88,
                    0.88);

            leg->SetBorderSize(0);
            leg->SetFillStyle(0);

            leg->AddEntry(
                hPosNorm[ipt],
                "positive",
                "lep");

            leg->AddEntry(
                hNegMirror[ipt],
                "negative, mirrored",
                "lep");

            leg->Draw();
        }
    }


    cComparison->SaveAs(
        Form(
            "%s/sdca_charge_comparison.png",
            outDir));


// ========================================================================
//
// 2. q*sDCA overestimation-side vs underestimation-side comparison
//
// ========================================================================

TH1D* hOver[nPtBins];
TH1D* hUnder[nPtBins];

TCanvas* cOverUnder =
    new TCanvas(
        "cOverUnder",
        "sDCA over/under comparison",
        1200,
        1000);

cOverUnder->Divide(2, 3);


for (ipt = 0; ipt < nPtBins; ++ipt) {

    // ------------------------------------------------------------
    // Clone known-good source histogram.
    // Keep original binning exactly.
    // ------------------------------------------------------------

    hOver[ipt] =
        (TH1D*)hPos[ipt]->Clone(
            Form("hOver_%d", ipt));

    hUnder[ipt] =
        (TH1D*)hPos[ipt]->Clone(
            Form("hUnder_%d", ipt));

    hOver[ipt]->Reset();
    hUnder[ipt]->Reset();


    const int nBins =
        hPos[ipt]->GetNbinsX();

    int ib = 0;


    // ------------------------------------------------------------
    // Fill only positive x = |sDCA|
    // ------------------------------------------------------------

    for (ib = 1; ib <= nBins; ++ib) {

        const double x =
            hPos[ipt]->GetXaxis()->GetBinCenter(ib);

        if (x <= 0.0)
            continue;


        // Find the corresponding bin at -x explicitly.
        // This is safer than assuming a bin-number reflection.
        const int mirrorBin =
            hPos[ipt]->GetXaxis()->FindBin(-x);


        // q*sDCA > 0:
        // positive charge at +x
        // negative charge at -x

        const double over =
            hPos[ipt]->GetBinContent(ib) +
            hNeg[ipt]->GetBinContent(mirrorBin);

        const double overErr2 =
            hPos[ipt]->GetBinError(ib) *
            hPos[ipt]->GetBinError(ib) +
            hNeg[ipt]->GetBinError(mirrorBin) *
            hNeg[ipt]->GetBinError(mirrorBin);


        // q*sDCA < 0:
        // positive charge at -x
        // negative charge at +x

        const double under =
            hPos[ipt]->GetBinContent(mirrorBin) +
            hNeg[ipt]->GetBinContent(ib);

        const double underErr2 =
            hPos[ipt]->GetBinError(mirrorBin) *
            hPos[ipt]->GetBinError(mirrorBin) +
            hNeg[ipt]->GetBinError(ib) *
            hNeg[ipt]->GetBinError(ib);


        hOver[ipt]->SetBinContent(
            ib,
            over);

        hOver[ipt]->SetBinError(
            ib,
            sqrt(overErr2));


        hUnder[ipt]->SetBinContent(
            ib,
            under);

        hUnder[ipt]->SetBinError(
            ib,
            sqrt(underErr2));
    }


    // ------------------------------------------------------------
    // Sanity check BEFORE plotting
    // ------------------------------------------------------------

    printf(
        "OVER/UNDER %-8s : nbins=%d  over=%12.0f  under=%12.0f  maxOver=%g  maxUnder=%g\n",
        ptLabels[ipt],
        hOver[ipt]->GetNbinsX(),
        hOver[ipt]->Integral(1, hOver[ipt]->GetNbinsX()),
        hUnder[ipt]->Integral(1, hUnder[ipt]->GetNbinsX()),
        hOver[ipt]->GetMaximum(),
        hUnder[ipt]->GetMaximum());


    // ------------------------------------------------------------
    // Plotting copies
    // ------------------------------------------------------------

    cOverUnder->cd(ipt + 1);

    gPad->SetLogy();
    gPad->SetTicks(1, 1);


    TH1D* hOverPlot =
        (TH1D*)hOver[ipt]->Clone(
            Form("hOverPlot_%d", ipt));

    TH1D* hUnderPlot =
        (TH1D*)hUnder[ipt]->Clone(
            Form("hUnderPlot_%d", ipt));


    // Common normalization preserves relative over/under yield
    double commonNorm =
        hOverPlot->Integral(1, hOverPlot->GetNbinsX()) +
        hUnderPlot->Integral(1, hUnderPlot->GetNbinsX());

    if (commonNorm > 0.0) {

        hOverPlot->Scale(
            1.0 / commonNorm);

        hUnderPlot->Scale(
            1.0 / commonNorm);
    }


    hOverPlot->GetXaxis()->SetRangeUser(
        0.0,
        1.0);

    hOverPlot->GetXaxis()->SetTitle(
        "|sDCA| (cm)");

    hOverPlot->GetYaxis()->SetTitle(
        "Normalized counts");


    hOverPlot->SetLineColor(kRed + 1);
    hOverPlot->SetMarkerColor(kRed + 1);
    hOverPlot->SetMarkerStyle(20);
    hOverPlot->SetMarkerSize(0.6);

    hUnderPlot->SetLineColor(kBlue + 1);
    hUnderPlot->SetMarkerColor(kBlue + 1);
    hUnderPlot->SetMarkerStyle(24);
    hUnderPlot->SetMarkerSize(0.6);


    hOverPlot->SetMinimum(1e-6);


    double maxY =
        hOverPlot->GetMaximum();

    if (hUnderPlot->GetMaximum() > maxY) {
        maxY =
            hUnderPlot->GetMaximum();
    }

    hOverPlot->SetMaximum(
        maxY * 3.0);


    hOverPlot->Draw("E");
    hUnderPlot->Draw("E SAME");


    TLine* cut =
        new TLine(
            0.5,
            1e-6,
            0.5,
            maxY * 3.0);

    cut->SetLineStyle(2);
    cut->SetLineWidth(2);
    cut->Draw();


    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.045);

    latex.DrawLatex(
        0.14,
        0.86,
        ptTitles[ipt]);


    if (ipt == 0) {

        TLegend* leg =
            new TLegend(
                0.50,
                0.68,
                0.88,
                0.88);

        leg->SetBorderSize(0);
        leg->SetFillStyle(0);

        leg->AddEntry(
            hOverPlot,
            "q #times sDCA > 0",
            "lep");

        leg->AddEntry(
            hUnderPlot,
            "q #times sDCA < 0",
            "lep");

        leg->Draw();
    }
}


cOverUnder->SaveAs(
    Form(
        "%s/sdca_over_under.png",
        outDir));

// ========================================================================
//
// 2b. Ratio of overestimation-side / underestimation-side
//
// ========================================================================

TCanvas* cQRatio =
    new TCanvas(
        "cQRatio",
        "sDCA side ratio",
        1200,
        1000);

cQRatio->Divide(2, 3);


for (ipt = 0; ipt < nPtBins; ++ipt) {

    cQRatio->cd(ipt + 1);

    gPad->SetTicks(1, 1);


    if (!hOver[ipt] || !hUnder[ipt])
        continue;


    TH1D* hRatio =
        (TH1D*)hOver[ipt]->Clone(
            Form("hSideRatio_%d", ipt));

    hRatio->Divide(
        hUnder[ipt]);


    hRatio->GetXaxis()->SetRangeUser(
        0.0,
        1.0);

    hRatio->GetXaxis()->SetTitle(
        "|sDCA| (cm)");

    hRatio->GetYaxis()->SetTitle(
        "N(q sDCA > 0) / N(q sDCA < 0)");


    hRatio->SetMinimum(0.0);

    double ratioMax = 3.0;

    if (ipt == 5) {
        ratioMax = 10.0;
    }

    hRatio->SetMaximum(ratioMax);


    hRatio->SetMarkerStyle(20);
    hRatio->SetMarkerSize(0.7);
    hRatio->SetLineColor(kBlack);
    hRatio->SetMarkerColor(kBlack);

    hRatio->Draw("E1");


    TLine* unity =
        new TLine(
            0.0,
            1.0,
            1.0,
            1.0);

    unity->SetLineStyle(2);
    unity->Draw();


    TLine* cut =
        new TLine(
            0.5,
            0.0,
            0.5,
            ratioMax);

    cut->SetLineStyle(3);
    cut->Draw();


    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.045);

    latex.DrawLatex(
        0.14,
        0.86,
        ptTitles[ipt]);
}


cQRatio->SaveAs(
    Form(
        "%s/sdca_over_under_ratio.png",
        outDir));

    // ========================================================================
    //
    // 3. Primary pT vs global pT
    //
    // ========================================================================

    TH2D* hPrimaryVsGlobal =
        dynamic_cast<TH2D*>(
            FindQaObject(
                file,
                "hPrimaryPtVsGlobalPt"));


    if (hPrimaryVsGlobal) {

        TCanvas* cPt =
            new TCanvas(
                "cPt",
                "primary vs global pT",
                850,
                750);

        cPt->SetRightMargin(0.15);
        cPt->SetLogz();
        cPt->SetTicks(1, 1);

        hPrimaryVsGlobal->GetXaxis()->SetTitle(
            "global #it{p}_{T} (GeV/#it{c})");

        hPrimaryVsGlobal->GetYaxis()->SetTitle(
            "primary #it{p}_{T} (GeV/#it{c})");

        hPrimaryVsGlobal->Draw("COLZ");


        TLine* diagonal =
            new TLine(
                0.0,
                0.0,
                40.0,
                40.0);

        diagonal->SetLineStyle(2);
        diagonal->SetLineWidth(2);
        diagonal->Draw();


        cPt->SaveAs(
            Form(
                "%s/primary_vs_global_pt.png",
                outDir));
    }
    else {

        cout << "WARNING: hPrimaryPtVsGlobalPt not found."
             << endl;
    }


    // ========================================================================
    //
    // 4. primary/global pT ratio vs sDCA
    //
    // ========================================================================

    TH2D* hRatioVsSDca =
        dynamic_cast<TH2D*>(
            FindQaObject(
                file,
                "hPtRatioVsSDca"));


    if (hRatioVsSDca) {

        TCanvas* cRatio =
            new TCanvas(
                "cRatio",
                "pT ratio vs sDCA",
                850,
                750);

        cRatio->SetRightMargin(0.15);
        cRatio->SetLogz();
        cRatio->SetTicks(1, 1);

        hRatioVsSDca->GetXaxis()->SetRangeUser(
            -2.0,
            2.0);

        hRatioVsSDca->GetXaxis()->SetTitle(
            "sDCA (cm)");

        hRatioVsSDca->GetYaxis()->SetTitle(
            "#it{p}_{T}^{prim} / #it{p}_{T}^{glob}");

        hRatioVsSDca->Draw("COLZ");


        TLine* ratioOne =
            new TLine(
                -2.0,
                1.0,
                2.0,
                1.0);

        ratioOne->SetLineStyle(2);
        ratioOne->SetLineWidth(2);
        ratioOne->Draw();


        TLine* cutPos =
            new TLine(
                0.5,
                0.0,
                0.5,
                10.0);

        cutPos->SetLineStyle(3);
        cutPos->Draw();


        TLine* cutNeg =
            new TLine(
                -0.5,
                0.0,
                -0.5,
                10.0);

        cutNeg->SetLineStyle(3);
        cutNeg->Draw();


        cRatio->SaveAs(
            Form(
                "%s/pt_ratio_vs_sdca.png",
                outDir));
    }
    else {

        cout << "WARNING: hPtRatioVsSDca not found."
             << endl;
    }


// ========================================================================
//
// 5. q*sDCA tail comparison
//
// q*sDCA > +0.5 : overestimated-pT side
// q*sDCA < -0.5 : opposite side
//
// ========================================================================

printf("\n");
printf("======================================================================\n");
printf("                  q*sDCA TAIL SUMMARY\n");
printf("======================================================================\n");
printf("\n");

printf("%-12s %14s %14s %12s %12s\n",
       "pT bin",
       "qDCA>+0.5",
       "qDCA<-0.5",
       "ratio",
       "excess[%]");

printf("----------------------------------------------------------------------\n");


for (ipt = 0; ipt < nPtBins; ++ipt) {

    if (!hOver[ipt] || !hUnder[ipt])
        continue;


    const int firstTailBin =
        hOver[ipt]->GetXaxis()->FindBin(
            0.5 + 1e-6);

    const int lastBin =
        hOver[ipt]->GetNbinsX();


    const double over =
        hOver[ipt]->Integral(
            firstTailBin,
            lastBin);

    const double under =
        hUnder[ipt]->Integral(
            firstTailBin,
            lastBin);


    double ratio = 0.0;
    double excess = 0.0;

    if (under > 0.0) {

        ratio =
            over / under;

        excess =
            100.0 *
            (over - under) /
            under;
    }


    printf("%-12s %14.0f %14.0f %12.5f %12.3f\n",
           ptLabels[ipt],
           over,
           under,
           ratio,
           excess);
}


printf("\n");
printf("Interpretation:\n");
printf("  ratio ~ 1 : symmetric tails\n");
printf("  ratio > 1 : excess on primary-pT overestimation side\n");
printf("======================================================================\n");


    // ========================================================================
//
// 6. Primary/global pT mismatch summary
//
// ========================================================================

if (hPrimaryVsGlobal) {

    printf("\n");
    printf("======================================================================\n");
    printf("             PRIMARY/GLOBAL pT MISMATCH SUMMARY\n");
    printf("======================================================================\n");
    printf("\n");

    const double primaryThresholds[] = {
        10.0,
        15.0,
        20.0,
        25.0
    };

    const int nThresholds =
        sizeof(primaryThresholds) /
        sizeof(primaryThresholds[0]);


    int ith = 0;

    for (ith = 0; ith < nThresholds; ++ith) {

        const double ptMin =
            primaryThresholds[ith];


        double total = 0.0;
        double ratio15 = 0.0;
        double ratio20 = 0.0;
        double globalBelow5 = 0.0;


        int ix = 0;
        int iy = 0;

        for (ix = 1;
             ix <= hPrimaryVsGlobal->GetNbinsX();
             ++ix) {

            const double globalPt =
                hPrimaryVsGlobal->
                GetXaxis()->
                GetBinCenter(ix);


            for (iy = 1;
                 iy <= hPrimaryVsGlobal->GetNbinsY();
                 ++iy) {

                const double primaryPt =
                    hPrimaryVsGlobal->
                    GetYaxis()->
                    GetBinCenter(iy);


                if (primaryPt <= ptMin)
                    continue;


                const double content =
                    hPrimaryVsGlobal->
                    GetBinContent(ix, iy);


                if (content <= 0.0)
                    continue;


                total += content;


                if (globalPt < 5.0)
                    globalBelow5 += content;


                if (globalPt > 0.0) {

                    const double ratio =
                        primaryPt / globalPt;


                    if (ratio > 1.5)
                        ratio15 += content;


                    if (ratio > 2.0)
                        ratio20 += content;
                }
            }
        }


        double fracBelow5 = 0.0;
        double frac15 = 0.0;
        double frac20 = 0.0;


        if (total > 0.0) {

            fracBelow5 =
                100.0 *
                globalBelow5 /
                total;

            frac15 =
                100.0 *
                ratio15 /
                total;

            frac20 =
                100.0 *
                ratio20 /
                total;
        }


        printf("primary pT > %.1f GeV/c\n",
               ptMin);

        printf("  total tracks                  = %.0f\n",
               total);

        printf("  global pT < 5 GeV/c           = %.0f  (%8.5f %%)\n",
               globalBelow5,
               fracBelow5);

        printf("  pTprim / pTglob > 1.5         = %.0f  (%8.5f %%)\n",
               ratio15,
               frac15);

        printf("  pTprim / pTglob > 2.0         = %.0f  (%8.5f %%)\n",
               ratio20,
               frac20);

        printf("\n");
    }


    printf("======================================================================\n");
}

    file->Close();

    cout << endl;
    cout << "Saved QA plots in:" << endl;
    cout << "  " << outDir << endl;
    cout << endl;
}