#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TBranch.h"

#include <iostream>
#include <iomanip>
#include <cmath>

using std::cout;
using std::endl;

static const int kNPthatBins = 11;

static const double kExpectedWeights[kNPthatBins] = {
    1.616e+0,
    1.355e-1,
    2.288e-2,
    5.524e-3,
    2.203e-3,
    3.437e-4,
    4.681e-5,
    8.532e-6,
    2.178e-6,
    1.198e-7,
    6.939e-9
};

// Adjust these labels if your production uses different pThat intervals.
static const char* kPthatLabels[kNPthatBins] = {
    "bin 0",
    "bin 1",
    "bin 2",
    "bin 3",
    "bin 4",
    "bin 5",
    "bin 6",
    "bin 7",
    "bin 8",
    "bin 9",
    "bin 10"
};

int FindExpectedWeight(double value)
{
    int iWeight = 0;
    double expected = 0.0;
    double relDiff = 0.0;

    for (iWeight = 0; iWeight < kNPthatBins; ++iWeight) {
        expected = kExpectedWeights[iWeight];
        relDiff = std::fabs(value - expected) / std::fabs(expected);

        if (relDiff < 1e-5) {
            return iWeight;
        }
    }

    return -1;
}

void CheckOneTree(TFile* file, const TString& treePath)
{
    TTree* tree = (TTree*)file->Get(treePath);

    cout << "\n==================================================" << endl;
    cout << "Tree: " << treePath << endl;

    if (!tree) {
        cout << "  NOT FOUND" << endl;
        return;
    }

    if (!tree->GetBranch("xsecWeight")) {
        cout << "  Branch xsecWeight NOT FOUND" << endl;
        return;
    }

    Float_t xsecWeight = 0.0f;
    tree->SetBranchAddress("xsecWeight", &xsecWeight);

    Long64_t entryCounts[kNPthatBins];

    int iInit = 0;
    int iCount = 0;
    int iPrint = 0;
    int matchedBin = -1;
    int nRepresentedWeights = 0;

    Long64_t iEntry = 0;
    Long64_t unknownEntries = 0;
    Long64_t nEntries = tree->GetEntries();

    double fraction = 0.0;

    for (iInit = 0; iInit < kNPthatBins; ++iInit) {
        entryCounts[iInit] = 0;
    }

    for (iEntry = 0; iEntry < nEntries; ++iEntry) {
        tree->GetEntry(iEntry);

        matchedBin = FindExpectedWeight((double)xsecWeight);

        if (matchedBin >= 0 && matchedBin < kNPthatBins) {
            entryCounts[matchedBin]++;
        } else {
            unknownEntries++;
        }
    }

    for (iCount = 0; iCount < kNPthatBins; ++iCount) {
        if (entryCounts[iCount] > 0) {
            nRepresentedWeights++;
        }
    }

    cout << "Total jet entries     : " << nEntries << endl;
    cout << "Expected values found : "
         << nRepresentedWeights << " / " << kNPthatBins << endl;
    cout << "Unknown-value entries : " << unknownEntries << endl;
    cout << endl;

    cout << std::setw(4)  << "idx"
         << std::setw(18) << "expected weight"
         << std::setw(18) << "jet entries"
         << std::setw(16) << "fraction [%]"
         << "  status"
         << endl;

    for (iPrint = 0; iPrint < kNPthatBins; ++iPrint) {
        fraction = 0.0;

        if (nEntries > 0) {
            fraction =
                100.0 * (double)entryCounts[iPrint] / (double)nEntries;
        }

        cout << std::setw(4) << iPrint
             << std::scientific
             << std::setprecision(8)
             << std::setw(18) << kExpectedWeights[iPrint]
             << std::fixed
             << std::setprecision(0)
             << std::setw(18) << (double)entryCounts[iPrint]
             << std::fixed
             << std::setprecision(4)
             << std::setw(16) << fraction;

        if (entryCounts[iPrint] > 0) {
            cout << "  OK";
        } else {
            cout << "  MISSING";
        }

        cout << endl;
    }

    cout << endl;

    if (nRepresentedWeights == kNPthatBins && unknownEntries == 0) {
        cout << "RESULT: All expected xsecWeight categories are present."
             << endl;
    } else {
        if (nRepresentedWeights != kNPthatBins) {
            cout << "WARNING: Found "
                 << nRepresentedWeights
                 << " expected weight categories; expected "
                 << kNPthatBins
                 << "."
                 << endl;
        }

        if (unknownEntries > 0) {
            cout << "WARNING: "
                 << unknownEntries
                 << " entries contain an unrecognized xsecWeight."
                 << endl;
        }
    }

    tree->ResetBranchAddresses();
}

void check_xsec_weights(
    const char* inputFile =
        "embedding_merged_MCReco1p5.root")
{
    TFile* file = TFile::Open(inputFile, "READ");

    if (!file || file->IsZombie()) {
        cout << "ERROR: Cannot open file: " << inputFile << endl;
        return;
    }

    const char* radii[] = {
        "R0.2",
        "R0.3",
        "R0.4"
    };

    const char* centralities[] = {
        "CENT_0_10",
        "MID_20_40",
        "PERI_60_80"
    };

    const int nRadii = 3;
    const int nCentralities = 3;

    for (int iR = 0; iR < nRadii; ++iR) {
        for (int iCent = 0; iCent < nCentralities; ++iCent) {
            const TString treePath =
                Form("%s/%s/JetTree",
                     radii[iR],
                     centralities[iCent]);

            CheckOneTree(file, treePath);
        }
    }

    file->Close();
    delete file;
}