#ifndef __CINT__
#include "TChain.h"
#include "TROOT.h"
#include "TString.h"
#include "TSystem.h"

#include "StChain/StChain.h"
#include "StChain/StMaker.h"

#include "StEmcADCtoEMaker/StEmcADCtoEMaker.h"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "St_db_Maker/St_db_Maker.h"

#include "StPicoCuts/StPicoCuts.h"

#include "loadSharedHFLibraries.C"

#include <cstdio>
#include <ctime>
#include <iostream>

#include "StPicoHFJetMaker/StPicoHFJetMaker.h"

#include "StRefMultCorr/CentralityMaker.h"
#include "StRefMultCorr/StRefMultCorr.h"

using namespace std;

#else
class StChain;
#endif

void runPicoHFJetMaker(TString inputFile, TString outputFile = "output",
                       const unsigned int makerMode = 0,
                       TString treeName = "picoDst",
                       bool isEmbedding = true,
                       bool storeOnlyTrigOrMc = false){

#ifdef __CINT__
  gROOT->LoadMacro("loadSharedHFLibraries.C");
  loadSharedHFLibraries();
#endif

  // pThat bins
  std::vector<TString> pt_bins_name; // CINT does not support c++11
  pt_bins_name.push_back("3_5");
  pt_bins_name.push_back("5_7");
  pt_bins_name.push_back("7_9");
  pt_bins_name.push_back("9_11");
  pt_bins_name.push_back("11_15");
  pt_bins_name.push_back("15_20");
  pt_bins_name.push_back("20_25");
  pt_bins_name.push_back("25_30");
  pt_bins_name.push_back("30_40");
  pt_bins_name.push_back("40_50");
  pt_bins_name.push_back("50_-1");

  float pt_bins[] = {3.0,  5.0,  7.0,  9.0,  11.0, 15.0,
                     20.0, 25.0, 30.0, 40.0, 50.0, -1};

  // cross-section weights

  float weights[] = {1.616e+0,  1.355e-01, 2.288e-02, 5.524e-03,
                     2.203e-03, 3.437e-04, 4.681e-05, 8.532e-06,
                     2.178e-06, 1.198e-07, 6.939e-09};

  StChain *chain = new StChain();
  // ========================================================================================
  // makerMode    = StPicoJetMaker::kAnalyze;
  // ========================================================================================

  cout << "Maker Mode    " << makerMode << endl;
  cout << "TreeName      " << treeName << endl;
  cout << "Input file: " << inputFile << endl;
  cout << "Loading bad run list: " << endl;
  TString badRunListFileName = "BadRunList_14.list";
  // check if exists
  if (!gSystem->AccessPathName(badRunListFileName)) {
    cout << "Bad run list file found!" << endl;
  } else {
    cout << "Bad run list file not found! Exiting..." << endl;
    exit(1);
  }

  StMessMgr *msg = StMessMgr::Instance();
  msg->SwitchOff("Could not make BEMC detector");

  StPicoDstMaker *picoDstMaker =
      new StPicoDstMaker(StPicoDstMaker::IoRead, inputFile,
                         "picoDstMaker"); // for local testing only
  St_db_Maker *dbMaker = new St_db_Maker("StarDb", "MySQL:StarDb");
  StEmcADCtoEMaker *adc = new StEmcADCtoEMaker();
  StPicoHFJetMaker *stPicoHFJetMaker =
      new StPicoHFJetMaker("jetTree", picoDstMaker, outputFile);

  stPicoHFJetMaker->setMakerMode(makerMode);
  stPicoHFJetMaker->setTreeName(treeName);
  stPicoHFJetMaker->setMcMode(false);
  stPicoHFJetMaker->setIsEmbedding(isEmbedding);

  StPicoCuts *picoCuts = new StPicoCuts("PicoCuts");
  stPicoHFJetMaker->setPicoCuts(picoCuts);

  // ---------------------------------------------------
  // -- Set Base cuts for HF analysis
  // -- File name of bad run list
  picoCuts->setBadRunListFileName(badRunListFileName);
  // -- ADD USER CUTS HERE ----------------------------
  picoCuts->setCutVzMax(30.);
  picoCuts->setCutVzVpdVzMax(3.);
  // picoCuts->setCutRefMult(396, 100000); //REMEMBER now only central
  picoCuts->setCutRefMult(0, 100000);

  // 2014 MB HFT triggers
  // picoCuts->addTriggerId(450050);    // vpdmb-5-p-nobsmd-hlt
  // picoCuts->addTriggerId(450060);    // vpdmb-5-p-nobsmd-hlt
  // picoCuts->addTriggerId(450005);    // vpdmb-5-p-nobsmd
  // picoCuts->addTriggerId(450015);    // vpdmb-5-p-nobsmd
  // picoCuts->addTriggerId(450025);    // vpdmb-5-p-nobsmd

  // 2014 MB triggers
  if (isEmbedding) {
    picoCuts->addTriggerId(450010); // vpdmb-30
    picoCuts->addTriggerId(450020); // vpdmb-30
    picoCuts->addTriggerId(450008); // vpdmb-5
    picoCuts->addTriggerId(450018); // vpdmb-5
  }
  
  // 2014 BHT2*VPD30 triggers
  if (!isEmbedding) {
    picoCuts->addTriggerId(450202);    // BHT2*VPD30
    picoCuts->addTriggerId(450212);    // BHT2*VPD30
  }

  // 2014 BHT3 triggers
  // picoCuts->addTriggerId(450203);    // BHT3
  // picoCuts->addTriggerId(450213);    // BHT3

  // starprod AuAu 2014 luminosity - split low mid high, without lumi = early
  // runs

  //       //SL16j triggers
  // picoCuts->addTriggerId(520802);    // VPDMB-5-p-hlt
  // picoCuts->addTriggerId(520812);    // VPDMB-5-p-hlt
  // picoCuts->addTriggerId(520822);    // VPDMB-5-p-hlt
  // picoCuts->addTriggerId(520832);    // VPDMB-5-p-hlt
  // picoCuts->addTriggerId(520842);    // VPDMB-5-p-hlt

  // picoCuts->addTriggerId(520001);    // VPDMB-5-p-sst
  // picoCuts->addTriggerId(520011);    // VPDMB-5-p-sst
  // picoCuts->addTriggerId(520021);    // VPDMB-5-p-sst
  // picoCuts->addTriggerId(520031);    // VPDMB-5-p-sst
  // picoCuts->addTriggerId(520041);    // VPDMB-5-p-sst
  // picoCuts->addTriggerId(520051);    // VPDMB-5-p-sst

  picoCuts->setCutNHitsFitMin(14);
  picoCuts->setCutNHitsFitnHitsMax(0.52);
  picoCuts->setCutRequireHFT(false);
  picoCuts->setCutTPCNSigma(30.0);

  picoCuts->setCutDcaMin(1.0);
  picoCuts->setCutEta(1);
  picoCuts->setCutPtRange(0.2, 30.0); // default
  picoCuts->setCutERange(0.2, 30.0);  //

  vector<float> R;
  R.push_back(0.2);
  R.push_back(0.3);
  R.push_back(0.4);

  // TPC setters

  stPicoHFJetMaker->setGhostMaxrap(1.0);
  stPicoHFJetMaker->setR(R);
  stPicoHFJetMaker->setJetPtMin(5); // default
  stPicoHFJetMaker->setCutETmin(0.2);
  stPicoHFJetMaker->setNJetsRemove(1);
  stPicoHFJetMaker->setR_bg(0.3);
  stPicoHFJetMaker->setHadronCorr(1.0);
  stPicoHFJetMaker->setTriggerThreshold(18); //~4.2 GeV is the HT2 threshold = 18 ADC
  stPicoHFJetMaker->setMaxNeutralFraction(0.95); // default
  stPicoHFJetMaker->setMaxDcaZHadronCorr(3.0); // cm, max DCA_z for global tracks used for hadronic correction

  // Systematics setters
  stPicoHFJetMaker->setDoTowErrPlus(false);
  stPicoHFJetMaker->setDoTowErrMinus(false);
  stPicoHFJetMaker->setDoTrackErr(false);

  float pThatmin = -1;
  float pThatmax = -1;
  float xsecWeight = -1;

  // read the first line of the inputfile if it is a list
  // place it into the headFile variable
  // to get the pThat range
  TString headFile;
  if (inputFile.Contains(".list")) {
    ifstream inFile(inputFile.Data());
    if (inFile.is_open()) {
      inFile >> headFile;
      cout << "First line of the input file: " << headFile << endl;
      inFile.close();
    }
  } else {
    headFile = inputFile;
  }

  if (isEmbedding){

    cout << "Running in " << "embedding"<< " mode." << endl;

  for (unsigned int pt_name = 0; pt_name < pt_bins_name.size(); pt_name++) {
    TString binTag = "pt" + pt_bins_name[pt_name] + "_";
    if (outputFile.BeginsWith(binTag) ||
        headFile.Contains(binTag) ||
        inputFile.Contains(binTag)) {
      pThatmin = pt_bins[pt_name];
      pThatmax = pt_bins[pt_name + 1];
      xsecWeight = weights[pt_name];
      cout << "pThat range found: " << pThatmin << " - " << pThatmax
           << " with xsecWeight = " << xsecWeight << endl;
      break;     
    }
  }
  if ( xsecWeight == -1) {
    cout << "No pThat range found for embedding! Exiting..." << endl;
    exit(1);
  }
    stPicoHFJetMaker->setMCparameters(pThatmin, pThatmax, xsecWeight);
  

  }
  else{
        cout << "Running in " << "data" << " mode." << endl;
  }


  // Also add protection for StRefMultCorr
  StRefMultCorr *grefmultCorrUtil;
  try {
    grefmultCorrUtil =
        CentralityMaker::instance()->getRefMultCorr_P18ih_VpdMB30_MidLow();
    if (grefmultCorrUtil) {
      stPicoHFJetMaker->setRefMultCorr(grefmultCorrUtil);
      cout << "RefMultCorr set successfully" << endl;
    } else {
      cout << "WARNING: RefMultCorr is null, skipping..." << endl;
    }
  } catch (...) {
    cout << "ERROR: Exception when getting RefMultCorr" << endl;
  }

  // ========================================================================================

  chain->Init();
  if (!picoDstMaker->chain()) {
    cout << "ERROR: Chain is null. Exiting..." << endl;
    return;
  }

  long int nEvents = 0;

  try {
    nEvents = picoDstMaker->chain()->GetEntries();
    cout << "Number of events in chain: " << nEvents << endl;
  } catch (...) {
    cout << "ERROR: Failed to get number of entries from chain!" << endl;
    return;
  }

  for (Int_t i = 0; i < nEvents; i++) {
    // Inside the event loop
    if (i % 100 == 0)
      cout << "Processing event " << i << "/" << nEvents << endl;

    chain->Clear();
    int iret = chain->Make(i);
    if (iret) {
      cout << "Bad return code! " << iret << " at event " << i << endl;
      break;
    }
  }
  cout << "****************************************** " << endl;
  cout << "Work done... now its time to close up shop!" << endl;
  cout << "****************************************** " << endl;
  chain->Finish();
  // double duration = (double) (clock() - start_t) / (double)CLOCKS_PER_SEC;
  cout << "****************************************** " << endl;
  cout << "total number of events  " << nEvents << endl;
  cout << "****************************************** " << endl;
  // cout << "Time needed " << duration << " s" << endl;
  cout << "****************************************** " << endl;

  delete chain;
}
