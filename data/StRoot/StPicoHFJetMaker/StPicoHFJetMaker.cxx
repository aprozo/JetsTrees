#include "StPicoHFJetMaker.h"
#include "JetInfo.h"

#include "BemcNewCalib.h"
#include "StEmcADCtoEMaker/StBemcData.h"
#include "StEmcADCtoEMaker/StEmcADCtoEMaker.h"
#include "StEmcRawMaker/StBemcRaw.h"
#include "StEmcRawMaker/StBemcTables.h"
#include "StEmcRawMaker/defines.h"

#include "TRandom3.h"
#include "TVector2.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>
#include <algorithm>

#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/Selector.hh"
#include "fastjet/config.h"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"

#include "MyJet.h"

using namespace std;

const char* kCentTag[4] = {
  "",            // index 0 unused
  "CENT_0_10",   // c3 = 1
  "MID_20_40",   // c3 = 2
  "PERI_60_80"   // c3 = 3
};


const double CUT_AREA_02 = 0.07; // R = 0.2
const double CUT_AREA_03 = 0.20; // R = 0.3
const double CUT_AREA_04 = 0.40; // R = 0.4

const double CUT_NEUTRAL_FRACTION = 0.95;

vector<MatchedJetPair> MatchJetsEtaPhi(const vector<MyJet> &McJets,
                                       const vector<MyJet> &RecoJets,
                                       const double &R);

ClassImp(StPicoHFJetMaker)

StPicoHFJetMaker::StPicoHFJetMaker(TString name, StPicoDstMaker *picoMaker,
                                   TString outputBaseFileName)
    : StPicoJetMaker(name, picoMaker, outputBaseFileName)
{
  mRefmultCorrUtil   = NULL;
  fMcSumPt           = 0.0f;
  fEventId = 0;
}


// _________________________________________________________
StPicoHFJetMaker::~StPicoHFJetMaker() {
  // destructor
}

// _________________________________________________________
int StPicoHFJetMaker::InitJets() {
  mADCtoEMaker = dynamic_cast<StEmcADCtoEMaker*>(GetMaker("Eread"));
  assert(mADCtoEMaker);
  mTables = mADCtoEMaker->getBemcData()->getTables();

  TH1::SetDefaultSumw2();

  mOutList->SetName("QA_histograms"); 

  TH1D* hcent9 = new TH1D("hcent9", "centrality9;bin;events", 10, -1, 9);
  mOutList->Add(hcent9);

  TAxis* ax = hcent9->GetXaxis();
  const char* lab9[10] = {
    "undef (-1)",
    "0: 0-5%",
    "1: 5-10%",
    "2: 10-20%",
    "3: 20-30%",
    "4: 30-40%",
    "5: 40-50%",
    "6: 50-60%",
    "7: 60-70%",
    "8: 70-80%"
  };
  int b;
  for (b = 1; b <= 10; ++b) ax->SetBinLabel(b, lab9[b-1]);
  ax->CenterLabels(true);
  mOutList->Add(new TH1D("hTowE", "BEMC tower E;E (GeV);counts", 400, 0.0, 40.0));
  mOutList->Add(new TH2D("hTowEtaPhi", "BEMC towers;#eta;#phi",
                       40, -1.0, 1.0, 120, 0.0, 2.0*TMath::Pi()));

  mOutList->Add(new TH1D("hTrackPt", "Primary tracks;p_{T} (GeV/c);counts", 400, 0.0, 40.0));
  mOutList->Add(new TH2D("hTrackEtaPhi", "Primary tracks;#eta;#phi",
                       40, -1.0, 1.0, 120, 0.0, 2.0*TMath::Pi()));
  mOutList->Add(new TH2D("hRhoVsRefMult",
                       "#rho vs refMult;refMult;#rho (GeV/c per unit area)",
                       800, 0, 800,     // adjust refMult range/binning to your dataset
                       200, 0.0, 200.0  // adjust rho range/binning as needed
                       ));

  if (!mIsEmbedding) {
    // Accepted-event QA split by the three analysis centrality classes.
    for (int c3 = 1; c3 <= 3; ++c3) {
      const char* centTag = kCentTag[c3];

      mOutList->Add(new TH1I(Form("hevents_acc_%s", centTag),
                             Form("accepted events (%s);bin;events", centTag),
                             2, 0, 2));
      mOutList->Add(new TH1D(Form("hevents_acc_weighted_%s", centTag),
                             Form("accepted events weighted (%s);bin;weighted events", centTag),
                             2, 0, 2));

      mOutList->Add(new TH1F(Form("hzvertex_acc_%s", centTag),
                             Form("accepted z-vertex (%s);z [cm];events", centTag),
                             250, -50, 50));
      mOutList->Add(new TH1D(Form("hzvertex_acc_weighted_%s", centTag),
                             Form("accepted z-vertex weighted (%s);z [cm];weighted events", centTag),
                             250, -50, 50));

      mOutList->Add(new TH1I(Form("hrunId_acc_%s", centTag),
                             Form("accepted events runId (%s);runId;events", centTag),
                             90913, 15076101, 15167014));
      mOutList->Add(new TH1D(Form("hrunId_acc_weighted_%s", centTag),
                             Form("accepted events runId weighted (%s);runId;weighted events", centTag),
                             90913, 15076101, 15167014));

      TH1D* hNormSummary = new TH1D(Form("hNormSummary_%s", centTag),
                                    Form("normalization summary (%s);quantity;counts", centTag),
                                    2, 0.5, 2.5);
      hNormSummary->GetXaxis()->SetBinLabel(1, "N_HT_acc");
      hNormSummary->GetXaxis()->SetBinLabel(2, "N_MB_equiv");
      hNormSummary->GetXaxis()->CenterLabels(true);
      mOutList->Add(hNormSummary);

      mOutList->Add(new TH1D(Form("hrunId_eqMb_%s", centTag),
                             Form("equivalent MB by run (%s);runId;N_{MB}^{equiv}", centTag),
                             90913, 15076101, 15167014));
    }
  }

  TDirectory* fileDir = gDirectory;
  fTreeRC.clear();
  fTreeRC.reserve(fR.size());

  for (size_t iR = 0; iR < fR.size(); ++iR) {
    const TString rName = Form("R%.1f", fR[iR]);
    TDirectory* rdir = fileDir->mkdir(rName);
    if (!rdir) rdir = (TDirectory*)fileDir->Get(rName);
    rdir->cd();

    std::vector<TTree*> treesC;  // 3 classes: 1..3 (we'll index 0..2)
    treesC.reserve(3);

    // book three classes (central, midcentral, peripheral)
    for (int c3 = 1; c3 <= 3; ++c3) {
      TDirectory* cdir = rdir->mkdir(kCentTag[c3]);
      if (!cdir) cdir = (TDirectory*)rdir->Get(kCentTag[c3]);
      cdir->cd();
      
      // ---- TTree per (R,class); NO centrality branch
      TTree* jetTree = new TTree("JetTree", "JetTree");
      jetTree->Branch("runId", &fRunNumber, "runId/I");
      jetTree->Branch("eventId", &fEventId,   "eventId/I");
      jetTree->Branch("centralityWeight", &fCentralityWeight, "centralityWeight/F"); 
      if (mIsEmbedding) {
        jetTree->Branch("xsecWeight", &fXsecWeight, "xsecWeight/F");
        jetTree->Branch("deltaR", &fDeltaR, "deltaR/F");
        jetTree->Branch("mc_pt", &fMcJet.pt, "mc_pt/F");
        jetTree->Branch("mc_eta", &fMcJet.eta, "mc_eta/F");
        jetTree->Branch("mc_phi", &fMcJet.phi, "mc_phi/F");
        jetTree->Branch("mc_area", &fMcJet.area, "mc_area/F");
        jetTree->Branch("mc_pt_lead", &fMcJet.pt_lead, "mc_pt_lead/F");
        jetTree->Branch("mc_n_constituents", &fMcJet.n_constituents, "mc_n_constituents/I");
        jetTree->Branch("mc_neutral_fraction", &fMcJet.neutral_fraction, "mc_neutral_fraction/F");
        jetTree->Branch("mc_sum_pt", &fMcSumPt, "mc_sum_pt/F");
        
      }
      jetTree->Branch("reco_pt", &fRecoJet.pt, "reco_pt/F");
      jetTree->Branch("reco_pt_corr", &fRecoJet.pt_corr, "reco_pt_corr/F");
      jetTree->Branch("reco_eta", &fRecoJet.eta, "reco_eta/F");
      jetTree->Branch("reco_phi", &fRecoJet.phi, "reco_phi/F");
      jetTree->Branch("reco_area", &fRecoJet.area, "reco_area/F");
      jetTree->Branch("reco_rho", &fRecoJet.rho, "reco_rho/F");
      jetTree->Branch("reco_pt_lead", &fRecoJet.pt_lead, "reco_pt_lead/F");
      jetTree->Branch("reco_n_constituents", &fRecoJet.n_constituents, "reco_n_constituents/I");
      jetTree->Branch("reco_neutral_fraction", &fRecoJet.neutral_fraction, "reco_neutral_fraction/F");
      jetTree->Branch("reco_trigger_match", &fRecoJet.trigger_match, "reco_trigger_match/O");

      treesC.push_back(jetTree);

      rdir->cd(); // back up for next class
    }

    fTreeRC.push_back(treesC);
    fileDir->cd();
  }

  return kStOK;
}

// _________________________________________________________
void StPicoHFJetMaker::ClearJets(Option_t *opt = "") { return; }

// _________________________________________________________
int StPicoHFJetMaker::FinishJets() {
  TDirectory* fileDir = gDirectory;
  if (!fileDir) return kStOK;

  const size_t nR = fR.size();

for (size_t iR = 0; iR < nR; ++iR) {
  TDirectory* rdir = dynamic_cast<TDirectory*>(fileDir->Get(Form("R%.1f", fR[iR])));
  if (!rdir) continue;

  for (int c3 = 1; c3 <= 3; ++c3) {
    const int ciTree = c3 - 1;  // 0..2

    TDirectory* cdir = dynamic_cast<TDirectory*>(rdir->Get(kCentTag[c3]));
    if (!cdir) continue;
    cdir->cd();

    // --- write the tree
    if (iR < fTreeRC.size() && ciTree >= 0 &&
        ciTree < (int)fTreeRC[iR].size() &&
        fTreeRC[iR][ciTree]) {
      fTreeRC[iR][ciTree]->Write();
    }
  } // c3
}   // iR

  fileDir->cd();
  return kStOK;
}

// _________________________________________________________
int StPicoHFJetMaker::MakeJets() {

    for (int i = 0; i < 4800; ++i) {
      Sump[i] = 0.0;
    }

  fMcSumPt = 0.0f; 
  Bool_t vetoReco = kFALSE; 
  
  TH1D *hcent9 = static_cast<TH1D *>(mOutList->FindObject("hcent9"));

  vector<fastjet::PseudoJet> jetTracks;
  vector<fastjet::PseudoJet> neutraljetTracks; // from bemc towers only
  vector<fastjet::PseudoJet> fullTracks;
  vector<fastjet::PseudoJet> MCjetTracks;

  fRunNumber = mPicoDst->event()->runId();
  fEventId = mPicoDst->event()->eventId(); // eventID
  int refMult = mPicoDst->event()->refMult();
  double vz = mPrimVtx.z();
  mRefmultCorrUtil->setEvent(fRunNumber, refMult, mPicoDst->event()->ZDCx(),
                             vz);
  fCentrality = mRefmultCorrUtil->centrality9(); // 0 = 0-5 %,..., 8 = 70-80
                                                 // %
  if (fCentrality == -1)
    return kStOK; // no fCentrality

  fCentralityWeight = mRefmultCorrUtil->weight();

  if (hcent9) hcent9->Fill(fCentrality, fCentralityWeight);

  TH1D* hTowE = (TH1D*)mOutList->FindObject("hTowE");
  TH2D* hTowEtaPhi = (TH2D*)mOutList->FindObject("hTowEtaPhi");
  TH1D* hTrackPt = (TH1D*)mOutList->FindObject("hTrackPt");
  TH2D* hTrackEtaPhi = (TH2D*)mOutList->FindObject("hTrackEtaPhi");

  // map to 3 classes: 1=0-10%, 2=20-40%, 3=60-80%, else 0 (=skip)
  int c3 = 0;
  if (fCentrality == 0 || fCentrality == 1)       c3 = 1; // 0-10%
  else if (fCentrality == 3 || fCentrality == 4)  c3 = 2; // 20-40%
  else if (fCentrality == 7 || fCentrality == 8)  c3 = 3; // 60-80%
  if (c3 == 0) {
  return kStOK;
  }

  if (!mIsEmbedding) {
    const char* centTag = kCentTag[c3];

    TH1I* hEventsAcc =
        static_cast<TH1I*>(mOutList->FindObject(Form("hevents_acc_%s", centTag)));
    TH1D* hEventsAccWeighted = static_cast<TH1D*>(
        mOutList->FindObject(Form("hevents_acc_weighted_%s", centTag)));
    TH1F* hVzAcc =
        static_cast<TH1F*>(mOutList->FindObject(Form("hzvertex_acc_%s", centTag)));
    TH1D* hVzAccWeighted = static_cast<TH1D*>(
        mOutList->FindObject(Form("hzvertex_acc_weighted_%s", centTag)));
    TH1I* hRunAcc =
        static_cast<TH1I*>(mOutList->FindObject(Form("hrunId_acc_%s", centTag)));
    TH1D* hRunAccWeighted = static_cast<TH1D*>(
        mOutList->FindObject(Form("hrunId_acc_weighted_%s", centTag)));
    TH1D* hNormSummary =
        static_cast<TH1D*>(mOutList->FindObject(Form("hNormSummary_%s", centTag)));
    TH1D* hRunEqMb =
        static_cast<TH1D*>(mOutList->FindObject(Form("hrunId_eqMb_%s", centTag)));
    const double runEqMbWeight = getRunEqMbWeight(fRunNumber);

    if (hEventsAcc) hEventsAcc->Fill(1);
    if (hEventsAccWeighted) hEventsAccWeighted->Fill(1.0, fCentralityWeight);
    if (hVzAcc) hVzAcc->Fill(vz);
    if (hVzAccWeighted) hVzAccWeighted->Fill(vz, fCentralityWeight);
    if (hRunAcc) hRunAcc->Fill(fRunNumber);
    if (hRunAccWeighted) hRunAccWeighted->Fill(fRunNumber, fCentralityWeight);
    if (hNormSummary) {
      hNormSummary->Fill(1.0);
      hNormSummary->Fill(2.0, runEqMbWeight);
    }
    if (hRunEqMb) hRunEqMb->Fill(fRunNumber, runEqMbWeight);
  }

  // MC tracks
  int noMCtracks = mPicoDst->numberOfMcTracks();
  for (int i = 0; i < noMCtracks; i++) {
    StPicoMcTrack *mctrk = (StPicoMcTrack *)mPicoDst->mcTrack(i);
    if (mctrk->idVtxStart() > 1)
      continue; // only primary tracks
    int geantId = mctrk->geantId();
    double mcpt = mctrk->pt();
    double mceta = mctrk->eta();
    if ((geantId > 3 && geantId < 7) || fabs(mceta) > 1.0 || mcpt < 0.2)
      continue;
    fMcSumPt += (Float_t)mcpt;
    TVector3 mcmom = mctrk->p();
    double mcphi = mcmom.Phi();
    if (mcphi < 0.0)
      mcphi += 2.0 * TMath::Pi();
    if (mcphi > 2.0 * TMath::Pi())
      mcphi -= 2.0 * TMath::Pi();

    double mcpx, mcpy, mcpz;
    mcpx = mcmom.x();
    mcpy = mcmom.y();
    mcpz = mcmom.z();

    double mcE = mctrk->energy();

    fastjet::PseudoJet inputMcParticle(mcpx, mcpy, mcpz, mcE);
    if (mctrk->charge() == 0) {
      inputMcParticle.set_user_index(0);
    } else
      inputMcParticle.set_user_index(i);
    MCjetTracks.push_back(inputMcParticle);
  }

if (mIsEmbedding && fpThatmax > 0.0 && !MCjetTracks.empty()) {
  float Rcheck = fR.empty() ? 0.4f : *std::max_element(fR.begin(), fR.end());
  fastjet::JetDefinition mc_jet_def_veto(fastjet::antikt_algorithm, Rcheck);
  fastjet::ClusterSequence mc_cs_veto(MCjetTracks, mc_jet_def_veto);
  std::vector<fastjet::PseudoJet> mcjets_veto =
      sorted_by_pt(mc_cs_veto.inclusive_jets(1.0)); // pT > 1 GeV

  const double ptMaxVeto = 1.5 * fpThatmax;

  for (size_t i = 0; i < mcjets_veto.size(); ++i) {
    if (mcjets_veto[i].perp() > ptMaxVeto) {
      return kStOK;
    }
  }
 
}
  // RC part
  GetCaloTrackMomentum(mPicoDst, mPrimVtx); // fill array Sump with momenta of tracks which are matched to BEMC

  StEmcGeom *mEmcGeom = StEmcGeom::getEmcGeom("bemc");
  StEmcPosition *mEmcPosition = new StEmcPosition();

//  double TOWE = 0;
for (int iTow = 0; iTow < 4800; iTow++) { // get btow info
  StPicoBTowHit *towHit = mPicoDst->btowHit(iTow);
  if (!towHit || towHit->isBad())
    continue;
  int realtowID = towHit->numericIndex2SoftId(iTow);
  if (BadTowerMap[realtowID])
    continue;

  double towE = GetTowerCalibEnergy(iTow + 1);
//  double towE = towHit->energy();

  if (doTowErrPlus == true)  towE = towE + 0.038 * towE;
  if (doTowErrMinus == true) towE = towE - 0.038 * towE;

  towE -= fHadronCorr * Sump[iTow];
  if (towE < 0) towE = 0;

  float Toweta_tmp = 0, Towphi = 0;
  if (mEmcGeom) mEmcGeom->getEtaPhi(realtowID, Toweta_tmp, Towphi);

  StThreeVectorF towerPosition = mEmcPosition->getPosFromVertex(
      StThreeVectorF(mPrimVtx.x(), mPrimVtx.y(), mPrimVtx.z()), realtowID);

  if (Towphi < 0)              Towphi += 2.0 * TMath::Pi();
  if (Towphi >= 2.0*TMath::Pi()) Towphi -= 2.0 * TMath::Pi();

  float Toweta = towerPosition.pseudoRapidity();
  double ET = towE / cosh(Toweta);

  if (hTowE) hTowE->Fill(towE, fCentralityWeight);
  if (hTowEtaPhi) hTowEtaPhi->Fill(Toweta, Towphi, fCentralityWeight);


  if (ET > 30.0) {
    vetoReco = kTRUE;
    break;   // no need to check further towers
  }

  // no clustering if vetoReco; we will check this later
  double px = ET * cos(Towphi);
  double py = ET * sin(Towphi);
  double pz = towE * tanh(Toweta);

  fastjet::PseudoJet inputTower(px, py, pz, towE);
  if (inputTower.perp() > fETmincut) {
    inputTower.set_user_index(0); // neutral
    int ADC = towHit->adc() >> 4;
    if (ADC > fTrgthresh) {
      inputTower.set_user_index(9999); // trigger towers
    }
    neutraljetTracks.push_back(inputTower);
  }
} // end tower loop

delete mEmcPosition;


// loop over primary tracks
for (unsigned int i = 0; i < mIdxPicoParticles.size(); i++) {
  StPicoTrack *trk = mPicoDst->track(mIdxPicoParticles[i]);

  if (doTrackErr) {
    static TRandom3 randGen;
    if (randGen.Rndm() > 0.96) continue;
  }

  const TVector3 p = trk->pMom();
  const double pT  = p.Perp();
  if (!(pT > 0)) continue;

  if (pT > 30.0) {
    vetoReco = kTRUE;
    break;   // no need to look at other tracks
  }

  const float eta = p.PseudoRapidity();
  if (fabs(eta) > 1.0) continue;
  float phi = trk->pMom().Phi();
  if (phi < 0) phi += 2.0*TMath::Pi();/*  */
  if (phi >= 2.0*TMath::Pi()) phi -= 2.0*TMath::Pi();

  if (hTrackPt) hTrackPt->Fill(pT, fCentralityWeight);
  if (hTrackEtaPhi) hTrackEtaPhi->Fill(eta, phi, fCentralityWeight);

  float dca = (mPrimVtx - trk->origin()).Mag();

  float charged = trk->charge();
  (void)dca; (void)charged;

  fastjet::PseudoJet pj(p.x(), p.y(), p.z(), p.Mag());

  if (mIsEmbedding) {
    if (trk->qaTruth() > 95) pj.set_user_index(trk->idTruth() - 1);
    else                     pj.set_user_index(trk->charge() ? 1 : 0);
  } else {
    pj.set_user_index(trk->charge() ? 1 : 0);
  }

  jetTracks.push_back(pj);

} // end loop over primary tracks

fullTracks.clear();

if (!vetoReco) {
  // build fullTracks from towers + tracks
  fullTracks = neutraljetTracks;
  fullTracks.insert(fullTracks.end(), jetTracks.begin(), jetTracks.end());
}

//==================================================================================//
// Jet part
//==================================================================================//
fastjet::AreaDefinition area_def(
    fastjet::active_area_explicit_ghosts,
    fastjet::GhostedAreaSpec(fGhostMaxrap, 1, 0.001));

//====================background estimate=======================//
float rho = 0.0;
TH2D* hRhoVsRefMult = (TH2D*)mOutList->FindObject("hRhoVsRefMult");

fastjet::JetDefinition jet_def_for_rho(fastjet::kt_algorithm, fRBg);
nJetsRemove = (c3 == 1 ? 2 : 1);

fastjet::Selector selector = (!fastjet::SelectorNHardest(nJetsRemove)) *
                             fastjet::SelectorAbsEtaMax(1.0) *
                             fastjet::SelectorPtMin(0.01);

if (!vetoReco && !fullTracks.empty()) {
  fastjet::JetMedianBackgroundEstimator bkgd_estimator(
      selector, jet_def_for_rho, area_def);
  bkgd_estimator.set_particles(fullTracks);
  rho = bkgd_estimator.rho();
}

if (hRhoVsRefMult && !vetoReco && !fullTracks.empty()) {
  hRhoVsRefMult->Fill(refMult, rho, fCentralityWeight);
}

//======================================================================//
const double ptMaxVeto = (mIsEmbedding && fpThatmax > 0.0) ? (2.0 * fpThatmax) : -1.0;
for (unsigned int i = 0; i < fR.size(); i++) {
  fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, fR[i]);
  float maxRapJet = 1 - fR[i];

  //==============================Reco jets===============================//
  std::vector<MyJet> myRecoJets;

  if (!vetoReco && !fullTracks.empty()) {
    fastjet::ClusterSequenceArea reco_cluster_seq(fullTracks, jet_def, area_def);
    std::vector<fastjet::PseudoJet> fjets_all =
        sorted_by_pt(reco_cluster_seq.inclusive_jets(fJetPtMin));

    fastjet::Selector fiducial_cut_selector = fastjet::SelectorAbsEtaMax(maxRapJet);
    std::vector<fastjet::PseudoJet> RecoJets = fiducial_cut_selector(fjets_all);

    myRecoJets.reserve(RecoJets.size());
    for (size_t j = 0; j < RecoJets.size(); ++j) {
      myRecoJets.push_back(MyJet(RecoJets[j], rho));
    }
  }

  if (ptMaxVeto > 0.0 && !fullTracks.empty()) {
  for (size_t jr = 0; jr < myRecoJets.size(); ++jr) {
    if (myRecoJets[jr].pt_corr > ptMaxVeto) {
      myRecoJets.clear();
      return kStOK; // veto whole event
    }
  }
}

  //======================================================================//

  // pick tree for this (R, class)
  TTree* jetTree = 0;
  const int ciTree = c3 - 1;
  if (i < fTreeRC.size() && ciTree >= 0 && ciTree < (int)fTreeRC[i].size())
    jetTree = fTreeRC[i][ciTree];

  //==================== Embedding mode ==================//
  if (mIsEmbedding) {
    //============================== MC jets ===============================//
    fastjet::ClusterSequenceArea mc_cluster_seq(MCjetTracks, jet_def, area_def);
    vector<fastjet::PseudoJet> Mcjets_all =
        sorted_by_pt(mc_cluster_seq.inclusive_jets(1.0));

    fastjet::Selector McFiducial_cut_selector =
        fastjet::SelectorAbsEtaMax(maxRapJet) *
        fastjet::SelectorPtMin(0.01);

    vector<fastjet::PseudoJet> McJets = McFiducial_cut_selector(Mcjets_all);
    vector<MyJet> myMcJets;
    myMcJets.reserve(McJets.size());
    for (size_t j = 0; j < McJets.size(); ++j) {
      myMcJets.push_back(MyJet(McJets[j], 0.0f)); // rho = 0 for truth
    }

    //========================= MC–Reco matching ===========================//
    vector<MatchedJetPair> MatchedJets =
        MatchJetsEtaPhi(myMcJets, myRecoJets, fR[i]);

    for (size_t im = 0; im < MatchedJets.size(); ++im) {
      const MatchedJetPair& mp = MatchedJets[im];

      fMcJet   = mp.first;
      fRecoJet = mp.second;
      fDeltaR  = fMcJet.deltaR(fRecoJet);

      const bool haveReco = (fRecoJet.pt >= 0.0);   // with default -999 this is false if no reco
      const bool haveMC   = (fMcJet.pt   >= 0.0);

      if (jetTree && (haveMC || haveReco)) {
          jetTree->Fill();
      }
    } // end loop over MatchedJets

  } else {
    //==================== Data mode ===========================//
    if (!vetoReco) {
      for (size_t jr = 0; jr < myRecoJets.size(); ++jr) {
        fRecoJet = myRecoJets[jr];
        fMcJet   = MyJet();   // dummy
        fDeltaR  = -1.0;

        if (!jetTree) continue;
          jetTree->Fill(); 
      }     
    } // if (!vetoReco)
  } // end embedding/data

} // end loop over R
  return kStOK;
}


//-----------------------------------------------------------------------------
////Correct tower energy
//-----------------------------------------------------------------------------
Double_t StPicoHFJetMaker::GetTowerCalibEnergy(Int_t TowerId) {
  StPicoBTowHit *tower =
      static_cast<StPicoBTowHit *>(mPicoDst->btowHit(TowerId - 1));
  Float_t pedestal, rms;
  Int_t status;
  mTables->getPedestal(BTOW, TowerId, 0, pedestal, rms);
  mTables->getStatus(BTOW, TowerId, status);
  Double_t *TowerCoeff;
  if (fRunNumber <= 15094020)
    TowerCoeff = CPre;
  else
    TowerCoeff = CLowMidHigh;
  Double_t calibEnergy = TowerCoeff[TowerId - 1] * (tower->adc() - pedestal);
  return calibEnergy;
}

//-----------------------------------------------------------------------------
////Correct tower eta for Vz position //// Not used anymore
//-----------------------------------------------------------------------------
Double_t StPicoHFJetMaker::vertexCorrectedEta(double eta, double vz) {
  double tower_theta = 2.0 * atan(exp(-eta));
  double z = 0.0;
  if (eta != 0.0)
    z = mBarrelRadius / tan(tower_theta);
  double z_diff = z - vz;
  double theta_corr = atan2(mBarrelRadius, z_diff);
  double eta_corr = -log(tan(theta_corr / 2.0));
  return eta_corr;
}

//-----------------------------------------------------------------------------
// Fill array with momentum of BEMC-matched tracks
//-----------------------------------------------------------------------------
Bool_t StPicoHFJetMaker::GetCaloTrackMomentum(StPicoDst *mPicoDst,
                                              TVector3 mPrimVtx) {
  // loop over global tracks  - towers
  UInt_t nTracks = mPicoDst->numberOfTracks();
  for (unsigned int itrack = 0; itrack < nTracks; itrack++) {
    StPicoTrack *trk = mPicoDst->track(itrack);
    TVector3 gMom = trk->gMom();
    // using global tracks
    double pT = gMom.Perp();
    if (pT != pT || pT < 0.2)
      continue;
    float eta = gMom.PseudoRapidity();
    if (fabs(eta) > 1)
      continue;
  //  float phi = gMom.Phi();

    float nHitsFit = trk->nHitsFit();
    float nHitsMax = trk->nHitsMax();
    if (nHitsFit < 15 || nHitsFit / nHitsMax < 0.52)
      continue; // some basic QA cuts
    double Bfield = mPicoDst->event()->bField();

    StPicoPhysicalHelix trkhelix = trk->helix(Bfield);
  //  float vtx_x = mPrimVtx.x();
  //  float vtx_y = mPrimVtx.y();
  //  float vtx_z = mPrimVtx.z();

    float dca_z = abs(trk->gDCAz(mPicoDst->event()->primaryVertex().z()));
    if (fabs(dca_z) > maxdcazhadroncorr)
      continue;
    int TowIndex = -99999;
    TowIndex = trk->bemcTowerIndex();
    float p = 0;
    if (TowIndex >= 0) {
      p = gMom.Mag();
      Sump[TowIndex] += p;
    }
  } // END global track loop
  return true;
}

//-----------------------------------------------------------------------------
// Jet matching function using only eta-phi criteria
//-----------------------------------------------------------------------------
vector<MatchedJetPair> MatchJetsEtaPhi(const vector<MyJet> &McJets,
                                       const vector<MyJet> &RecoJets,
                                       const double &R) {
  const double matchRadius = 0.6*R; 
  vector<char> recoUsed(RecoJets.size(), 0);
  vector<MatchedJetPair> matchedJets;
  matchedJets.reserve(McJets.size() + RecoJets.size());

  // For each MC jet, find the closest unused reco jet within matchRadius
  for (const auto &mcJet : McJets) {
    int bestIdx = -1;
    double bestDr = matchRadius;
    for (size_t j = 0; j < RecoJets.size(); ++j) {
      if (recoUsed[j]) continue;
      double dr = mcJet.deltaR(RecoJets[j]); 
      if (dr < bestDr) {
        bestDr = dr;
        bestIdx = (int)j;
      }
    }
    if (bestIdx >= 0) {
      matchedJets.emplace_back(mcJet, RecoJets[bestIdx]);
      recoUsed[bestIdx] = 1;
    } else {
      matchedJets.emplace_back(mcJet, MyJet()); // unmatched MC -> dummy reco
    }
  }

  // Append unmatched reco jets
  for (size_t j = 0; j < RecoJets.size(); ++j) {
    if (!recoUsed[j]) matchedJets.emplace_back(MyJet(), RecoJets[j]);
  }

  return matchedJets;
}