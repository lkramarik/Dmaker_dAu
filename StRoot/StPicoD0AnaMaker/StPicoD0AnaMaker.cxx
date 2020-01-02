#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoHFMaker/StHFCuts.h"
#include "phys_constants.h"
#include "StPicoD0AnaMaker.h"
#include "StPicoKFVertexFitter/StPicoKFVertexFitter.h"
ClassImp(StPicoD0AnaMaker)

const int nptBins=3;
//float const bdtCuts[nptBins] = {0.36, 0.3, 0.29}; //original CM
float const bdtCuts[nptBins] = {0.46, 0.4, 0.39};
//float const bdtCuts[nptBins] = {-1, -1, -1};
const float momBins[nptBins+1] = {1,2,3,5};
TString ptbin[nptBins] = {"12", "23", "35"};

// _________________________________________________________
StPicoD0AnaMaker::StPicoD0AnaMaker(char const* name, StPicoDstMaker* picoMaker, char const* outputBaseFileName) :
        StPicoHFMaker(name, picoMaker, outputBaseFileName),
        mOutFileBaseName(outputBaseFileName), mSwitchRefit(false){
    // constructor
}

// _________________________________________________________
StPicoD0AnaMaker::~StPicoD0AnaMaker() {
    // destructor
}

// _________________________________________________________
int StPicoD0AnaMaker::InitHF() {
    // -- INITIALIZE USER HISTOGRAMS ETC HERE -------------------
    //    add them to the output list mOutList which is automatically written
    //
    // EXAMPLE //  mOutList->Add(new TH1F(...));
    // EXAMPLE //  TH1F* hist = static_cast<TH1F*>(mOutList->Last());
//    mOutList->Add(new TH2F("h_piTOF","h_piTOF",100,0,10, 250, -1, 1.5));
//    mOutList->Add(new TH2F("h_kTOF","h_kTOF",100,0,10, 250, -1, 1.5));
//    mOutList->Add(new TH2F("h_pTOF","h_pTOF",100,0,10, 250, -1, 1.5));
//
//    mOutList->Add(new TH2F("h_piTOF_20","h_piTOF_20",100,0,10, 300, 0, 1));
//    mOutList->Add(new TH2F("h_kTOF_20","h_kTOF_20",100,0,10, 300, 0, 1));
//    mOutList->Add(new TH2F("h_pTOF_20","h_pTOF_20",100,0,10, 300, 0, 1));
//
//    mOutList->Add(new TH2F("h_piTOF_HFT","h_piTOF_HFT",100,0,10, 300, 0, 1));
//    mOutList->Add(new TH2F("h_kTOF_HFT","h_kTOF_HFT",100,0,10, 300, 0, 1));
//    mOutList->Add(new TH2F("h_pTOF_HFT","h_pTOF_HFT",100,0,10, 300, 0, 1));
//
//    mOutList->Add(new TH2F("h_piTOF_HFT_20","h_piTOF_HFT_20",100,0,10, 300, 0, 1));
//    mOutList->Add(new TH2F("h_kTOF_HFT_20","h_kTOF_HFT_20",100,0,10, 300, 0, 1));
//    mOutList->Add(new TH2F("h_pTOF_HFT_20","h_pTOF_HFT_20",100,0,10, 300, 0, 1));
//
//    mOutList->Add(new TH2F("h_piTOFbeta","h_piTOFbeta",500,0,10, 300, 0, 1));
//    mOutList->Add(new TH2F("h_kTOFbeta","h_kTOFbeta",500,0,10, 300, 0, 1));
//    mOutList->Add(new TH2F("h_pTOFbeta","h_pTOFbeta",500,0,10, 300, 0, 1));
//

    mOutList->Add(new TH2F("hD0VsRemoved","hD0VsRemoved", 100, -0.5, 99.5, 100, -0.5, 99.5));

    mOutList->Add(new TH1F("hNTracksRemoved","hNTracksRemoved", 2000, -0.5, 1999.5));
    mOutList->Add(new TH1F("hNTracksPrimary","hNTracksPrimary", 200, -0.5, 199.5));
    mOutList->Add(new TH1F("hNTracksDiffRemovedPrimary","hNTracksDiffRemovedPrimary", 200, -0.5, 199.5));
    mOutList->Add(new TH1F("hNTracksDiffRemovedGlobal","hNTracksDiffRemovedGlobal", 200, -0.5, 199.5));
    mOutList->Add(new TH1F("hHotSpotDiffRemovedPrimary","hHotSpotDiffRemovedPrimary", 200, -0.5, 199.5));
    mOutList->Add(new TH1F("hNTracksGoodToFit","hNTracksGoodToFit", 200, -0.5, 199.5));

    mOutList->Add(new TH1F("hPionPt","hPionPt", 100, 0, 10));
    mOutList->Add(new TH1F("hKaonPt","hKaonPt", 100, 0, 10));

    mOutList->Add(new TH1F("hRemovedPairMass","hRemovedPairMass", 300, 1.7, 2.0));
    mOutList->Add(new TH1F("hInvMassBDT12","hInvMassBDT12", 300, 1.7, 2.0));
    mOutList->Add(new TH1F("hInvMassBDT23","hInvMassBDT23", 300, 1.7, 2.0));
    mOutList->Add(new TH1F("hInvMassBDT35","hInvMassBDT35", 300, 1.7, 2.0));

    mOutList->Add(new TH1F("hPVDiffX","hPVDiffX", 1001, -0.02002, 0.01998));
    mOutList->Add(new TH1F("hPVDiffY","hPVDiffY", 1001, -0.02002, 0.01998));
    mOutList->Add(new TH1F("hPVDiffZ","hPVDiffZ", 1001, -0.02002, 0.01998));

    mOutList->Add(new TH1F("hPVDiffXRemoved","hPVDiffXRemoved", 1000, -2.0, 2.0));
    mOutList->Add(new TH1F("hPVDiffYRemoved","hPVDiffYRemoved", 1000, -2.0, 2.0));
    mOutList->Add(new TH1F("hPVDiffZRemoved","hPVDiffZRemoved", 1000, -2.0, 2.0));

//    mOutList->Add(new TH2F("h_pnsigma","h_pnsigma",1000,0,10, 99, -5, 5));

//    mOutList->Add(new TH2F("h_dedx","h_dedx", 1000, 0, 10, 1000, 0, 10));
//h_tracktest
    mOutList->Add(new TH1D("h_tracktest","h_tracktest", 6, 0.5, 6.5));
    mOutList->Add(new TH1D("h_tracktest_TOF","h_tracktest_TOF", 6, 0.5, 6.5));

    mOutFileBaseName = mOutFileBaseName.ReplaceAll(".root", "");
    TString ntpVars = "grefMult:refMult:runId:eventId:ZDC:BBC:hotSpot:primary:diffRemovedPrimary:pi1_pt:pi1_dca:pi1_nSigma:pi1_nHitFit:pi1_TOFinvbeta:pi1_betaBase:k_pt:k_dca:k_nSigma:k_nHitFit:k_TOFinvbeta:k_betaBase:dcaDaughters:primVz:primVzVpd:D_theta:cosTheta:D_decayL:dcaD0ToPv:D_cosThetaStar"
                      ":D_pt:D_mass";
//    ntp_kaon = new TNtuple("ntp_kaon", "kaon tree","k_pt:k_phi:k_eta:k_nSigma:k_nHitFit:k_TOFinvbeta:pi_eventId:pi_runId");
//    ntp_pion = new TNtuple("ntp_pion", "pion tree","pi_pt:pi_phi:pi_eta:pi_nSigma:pi_nHitFit:pi_TOFinvbeta:k_eventId:k_runId");
    ntp_DMeson_Signal = new TNtuple("ntp_signal","DMeson TreeSignal", ntpVars);
    ntp_DMeson_Background = new TNtuple("ntp_background","DMeson TreeBackground", ntpVars);

    if (mSwitchRefit) {
        TString dir = "./StRoot/weights/";
        TString prefix = "TMVAClassification";

        for (int pT = 0; pT < nptBins; pT++) {
            reader[pT] = new TMVA::Reader("!Color:!Silent");
            reader[pT]->AddVariable("k_dca", &k_dca[pT]);
            reader[pT]->AddVariable("pi1_dca", &pi1_dca[pT]);
            reader[pT]->AddVariable("dcaDaughters", &dcaDaughters[pT]);
            reader[pT]->AddVariable("cosTheta", &cosTheta[pT]);
            reader[pT]->AddVariable("D_decayL", &D_decayL[pT]);
            reader[pT]->AddVariable("dcaD0ToPv", &dcaD0ToPv[pT]);
//            reader[pT]->AddVariable("D_cosThetaStar", &thetaStar[pT]);

            TString methodName = "BDT method";
            TString weightfile = dir + prefix + TString("_BDT.weights.pt") + ptbin[pT] + TString(".xml");
            reader[pT]->BookMVA(methodName, weightfile);
        }
    }

    return kStOK;
}

// _________________________________________________________
void StPicoD0AnaMaker::ClearHF(Option_t *opt="") {
    return;
}

// _________________________________________________________
int StPicoD0AnaMaker::FinishHF() {
    ntp_DMeson_Signal -> Write(ntp_DMeson_Signal->GetName(), TObject::kOverwrite);
    ntp_DMeson_Background -> Write(ntp_DMeson_Background->GetName(), TObject::kOverwrite);
//    ntp_pion -> Write(ntp_pion->GetName(), TObject::kOverwrite);
//    ntp_kaon -> Write(ntp_kaon->GetName(), TObject::kOverwrite);
    return kStOK;
}
// _________________________________________________________
int StPicoD0AnaMaker::MakeHF() {
    createCandidates();
//    analyzeCandidates();

//    TH2F *h_pTOF = static_cast<TH2F*>(mOutList->FindObject("h_pTOF"));
//
//    TH2F *h_piTOF_20 = static_cast<TH2F*>(mOutList->FindObject("h_piTOF_20"));
//    TH2F *h_kTOF_20 = static_cast<TH2F*>(mOutList->FindObject("h_kTOF_20"));
//    TH2F *h_pTOF_20 = static_cast<TH2F*>(mOutList->FindObject("h_pTOF_20"));
//
//    TH2F *h_piTOF_HFT = static_cast<TH2F*>(mOutList->FindObject("h_piTOF_HFT"));
//    TH2F *h_kTOF_HFT = static_cast<TH2F*>(mOutList->FindObject("h_kTOF_HFT"));
//    TH2F *h_pTOF_HFT = static_cast<TH2F*>(mOutList->FindObject("h_pTOF_HFT"));
//
//    TH2F *h_piTOF_HFT_20 = static_cast<TH2F*>(mOutList->FindObject("h_piTOF_HFT_20"));
//    TH2F *h_kTOF_HFT_20 = static_cast<TH2F*>(mOutList->FindObject("h_kTOF_HFT_20"));
//    TH2F *h_pTOF_HFT_20 = static_cast<TH2F*>(mOutList->FindObject("h_pTOF_HFT_20"));
//
//    TH2F *h_piTOFbeta = static_cast<TH2F*>(mOutList->FindObject("h_piTOFbeta"));
//    TH2F *h_kTOFbeta = static_cast<TH2F*>(mOutList->FindObject("h_kTOFbeta"));
//    TH2F *h_pTOFbeta = static_cast<TH2F*>(mOutList->FindObject("h_pTOFbeta"));
//
//    TH2F *h_pinsigma = static_cast<TH2F*>(mOutList->FindObject("h_pinsigma"));
//    TH2F *h_knsigma = static_cast<TH2F*>(mOutList->FindObject("h_knsigma"));
//    TH2F *h_pnsigma = static_cast<TH2F*>(mOutList->FindObject("h_pnsigma"));
//
//    TH2F *h_dedx = static_cast<TH2F*>(mOutList->FindObject("h_dedx"));
//
//    TVector3 pVtx = mPicoDst->event()->primaryVertex();
//
//    UInt_t nTracks = mPicoDst->numberOfTracks();
//    for (unsigned short iTrack = 0; iTrack < nTracks; ++iTrack){
//        StPicoTrack const* trk = mPicoDst->track(iTrack);
//        if (!trk) continue;
//        StPhysicalHelixD helix = trk->helix(mBField);
//        TVector3 momentum = trk->gMom(pVtx, mPicoDst->event()->bField());
//
//        if (!(trk->nHitsFit()>=15)) continue;
//        if (!(fabs(momentum.pseudoRapidity()) <= 1.0)) continue;
//
//        if (fabs(trk->nSigmaPion())<3.0){
//            if (mHFCuts->isTOFHadron(trk, mHFCuts->getTofBetaBase(trk), StHFCuts::kPion)) {
//                    h_piTOF->Fill(trk->gPt(),mHFCuts->getTofBetaBase(trk));
//                    float oneOverBeta = getOneOverBeta(trk, mHFCuts->getTofBetaBase(trk), StHFCuts::kPion);
//                    h_piTOFbeta->Fill(trk->gPt(),oneOverBeta);
//                    if (trk->nHitsFit()>=20) h_piTOF_20->Fill(trk->gPt(),oneOverBeta);
//                    if (trk->isHFTTrack()) h_piTOF_HFT->Fill(trk->gPt(),oneOverBeta);
//                    if ((trk->isHFTTrack()) && (trk->nHitsFit()>=20)) h_piTOF_HFT_20->Fill(trk->gPt(),oneOverBeta);
//            }
//        }
//
//        if (fabs(trk->nSigmaKaon())<3.0){
//            if (mHFCuts->isTOFHadron(trk, mHFCuts->getTofBetaBase(trk), StHFCuts::kKaon)) {
//                h_kTOF->Fill(trk->gPt(),mHFCuts->getTofBetaBase(trk));
//                float oneOverBeta = getOneOverBeta(trk, mHFCuts->getTofBetaBase(trk), StHFCuts::kKaon);
//                h_kTOFbeta->Fill(trk->gPt(),oneOverBeta);
//                if (trk->nHitsFit()>=20) h_kTOF_20->Fill(trk->gPt(),oneOverBeta);
//                if (trk->isHFTTrack()) h_kTOF_HFT->Fill(trk->gPt(),oneOverBeta);
//                if ((trk->isHFTTrack()) && (trk->nHitsFit()>=20)) h_kTOF_HFT_20->Fill(trk->gPt(),oneOverBeta);
//            }
//        }
//
//        if (fabs(trk->nSigmaProton())<3.0){
//            if (mHFCuts->isTOFHadron(trk, mHFCuts->getTofBetaBase(trk), StHFCuts::kProton)) {
//                h_pTOF->Fill(trk->gPt(),mHFCuts->getTofBetaBase(trk));
//                float oneOverBeta = getOneOverBeta(trk, mHFCuts->getTofBetaBase(trk), StHFCuts::kProton);
//                h_pTOFbeta->Fill(trk->gPt(),oneOverBeta);
//                if (trk->nHitsFit()>=20) h_pTOF_20->Fill(trk->gPt(),oneOverBeta);
//                if (trk->isHFTTrack()) h_pTOF_HFT->Fill(trk->gPt(),oneOverBeta);
//                if ((trk->isHFTTrack()) && (trk->nHitsFit()>=20)) h_pTOF_HFT_20->Fill(trk->gPt(),oneOverBeta);
//            }
//        }
//
//        h_pinsigma->Fill(momentum.Mag(),trk->nSigmaPion());
//        h_knsigma->Fill(momentum.Mag(),trk->nSigmaKaon());
//        h_pnsigma->Fill(momentum.Mag(),trk->nSigmaProton());
//        h_dedx->Fill(momentum.Mag(),trk->dEdx());
//
//    } // .. end tracks loop

    return kStOK;
}

// _________________________________________________________
int StPicoD0AnaMaker::createCandidates() {
    //make it run faster:
    if (!(mPicoEvent->BBCx()<950000)) return kStOK;
    nGoodTracks=0;
//    if (!(abs(mPrimVtx.x())<0.6)) return kStOK;
//    if (!(abs(mPrimVtx.y())<0.6)) return kStOK;

    TH2F *hD0VsRemoved = static_cast<TH2F*>(mOutList->FindObject("hD0VsRemoved"));

    TH1F *hKaonPt = static_cast<TH1F*>(mOutList->FindObject("hKaonPt"));
    TH1F *hPionPt = static_cast<TH1F*>(mOutList->FindObject("hPionPt"));

    TH1F *hNTracksRemoved = static_cast<TH1F*>(mOutList->FindObject("hNTracksRemoved"));
    TH1F *hNTracksPrimary = static_cast<TH1F*>(mOutList->FindObject("hNTracksPrimary"));
    TH1F *hNTracksDiffRemovedPrimary = static_cast<TH1F*>(mOutList->FindObject("hNTracksDiffRemovedPrimary"));
    TH1F *hNTracksDiffRemovedGlobal = static_cast<TH1F*>(mOutList->FindObject("hNTracksDiffRemovedGlobal"));
    TH1F *hHotSpotDiffRemovedPrimary = static_cast<TH1F*>(mOutList->FindObject("hHotSpotDiffRemovedPrimary"));
    TH1F *hNTracksGoodToFit = static_cast<TH1F*>(mOutList->FindObject("hNTracksGoodToFit"));

    UInt_t nTracks = mPicoDst->numberOfTracks();
    Int_t nD0 = 0;
    nPrimary = 0;

    for (unsigned short iTrack = 0; iTrack < nTracks; ++iTrack) {
        StPicoTrack* trk = mPicoDst->track(iTrack);
        if (abs(trk->gMom().PseudoRapidity())>1) continue;

        if (trk->isPrimary()) {
            nPrimary++;
            primaryTracks.push_back(iTrack);

            if (mHFCuts->isGoodPion(trk)) {
                mIdxPicoPions.push_back(iTrack);
                hPionPt->Fill(trk->gPt());
            }
            if (mHFCuts->isGoodKaon(trk)){
                mIdxPicoKaons.push_back(iTrack);
                hKaonPt->Fill(trk->gPt());
            }
        }
    }

    TVector3 useVertex(mPrimVtx.x(), mPrimVtx.y(), mPrimVtx.z());
    if (mSwitchRefit) useVertex=refitVertex(true);
//        if (tracksToRemove.size()==nPrimary) {
//            cout<<useVertex.x()<<" "<<useVertex.y()<<" "<<useVertex.z()<<endl;
//            cout<<mPrimVtx.x()<<" "<<mPrimVtx.y()<<" "<<mPrimVtx.z()<<endl;
//        }

    for (unsigned short idxPion1 = 0; idxPion1 < mIdxPicoPions.size(); ++idxPion1) {
        StPicoTrack *pion1 = mPicoDst->track(mIdxPicoPions[idxPion1]);
        for (unsigned short idxKaon = 0; idxKaon < mIdxPicoKaons.size(); ++idxKaon) {
            StPicoTrack *kaon = mPicoDst->track(mIdxPicoKaons[idxKaon]);
            StHFPair *pair = new StHFPair(pion1, kaon, mHFCuts->getHypotheticalMass(StPicoCutsBase::kPion),mHFCuts->getHypotheticalMass(StPicoCutsBase::kKaon), mIdxPicoPions[idxPion1],mIdxPicoKaons[idxKaon], useVertex, mBField, kTRUE);

            if (!mHFCuts->isGoodSecondaryVertexPair(pair)) continue;

            bool isD0 = false;
            if((kaon->charge() + pion1->charge() == 0) ) isD0=true;

            Float_t primary = 0;
            if (pion1->isPrimary()) primary = 3;
            if (kaon->isPrimary()) primary = 4;
            if (pion1->isPrimary() && kaon->isPrimary()) primary = 2;

            Float_t hotSpot=0;
            if (mHFCuts->checkHotSpot(&mPrimVtx)) hotSpot=1;

            const int nNtVars = ntp_DMeson_Signal->GetNvar();
            float ntVar[nNtVars];
            int ii=0;

            ntVar[ii++] = mPicoEvent->grefMult();
            ntVar[ii++] = mPicoEvent->refMult();
            ntVar[ii++] = mPicoEvent->runId();
            ntVar[ii++] = mPicoEvent->eventId();
            ntVar[ii++] = mPicoEvent->ZDCx();
            ntVar[ii++] = mPicoEvent->BBCx();
            ntVar[ii++] = hotSpot;

            ntVar[ii++] = primary;
            ntVar[ii++] = nGoodTracks;

            ntVar[ii++] = pion1->gPt();
            ntVar[ii++] = pair->particle1Dca();
            ntVar[ii++] = pion1->nSigmaPion();
            ntVar[ii++] = pion1->nHitsFit();
            ntVar[ii++] = mHFCuts->getOneOverBeta(pion1, mHFCuts->getTofBetaBase(pion1), StPicoCutsBase::kPion);
            ntVar[ii++] = mHFCuts->getTofBetaBase(pion1);

            ntVar[ii++] = kaon->gPt();
            ntVar[ii++] = pair->particle2Dca();
            ntVar[ii++] = kaon->nSigmaKaon();
            ntVar[ii++] = kaon->nHitsFit();
            ntVar[ii++] = mHFCuts->getOneOverBeta(kaon, mHFCuts->getTofBetaBase(kaon), StPicoCutsBase::kKaon);
            ntVar[ii++] = mHFCuts->getTofBetaBase(kaon);

            ntVar[ii++] = pair->dcaDaughters();
            ntVar[ii++] = mPrimVtx.z();
            ntVar[ii++] = mPicoEvent->vzVpd();

            ntVar[ii++] = pair->pointingAngle();
            ntVar[ii++] = cos(pair->pointingAngle());
            ntVar[ii++] = pair->decayLength();
            ntVar[ii++] = pair->DcaToPrimaryVertex(); //(pair->decayLength())*sin(pair->pointingAngle());
            ntVar[ii++] = pair->cosThetaStar();

            ntVar[ii++] = pair->pt();
            ntVar[ii++] = pair->m();

            if (isD0) {
                ntp_DMeson_Signal->Fill(ntVar);
                nD0++;
            } else {
                ntp_DMeson_Background->Fill(ntVar);
            }

        }  // for (unsigned short idxKaon = 0; idxKaon < mIdxPicoKaons.size(); ++idxKaon)
    } // for (unsigned short idxPion1 = 0; idxPion1 < mIdxPicoPions.size(); ++idxPion1)

    if (nD0>0) {
        hNTracksRemoved->Fill(tracksToRemove.size());
        hNTracksGoodToFit->Fill(nGoodTracks);
        hNTracksPrimary->Fill(nPrimary);
        hNTracksDiffRemovedPrimary->Fill(nPrimary-tracksToRemove.size());
        hNTracksDiffRemovedGlobal->Fill(nTracks-tracksToRemove.size());
        if(mPrimVtx.x()>-0.25 && mPrimVtx.x()<-0.16 && mPrimVtx.y()>-0.25 && mPrimVtx.y()<-0.16) hHotSpotDiffRemovedPrimary->Fill(nPrimary-tracksToRemove.size());
    }

    hD0VsRemoved->Fill(tracksToRemove.size(), nD0);

    mIdxPicoPions.clear();
    mIdxPicoPions.shrink_to_fit();

    primaryTracks.clear();
    primaryTracks.shrink_to_fit();

    mIdxPicoKaons.clear();
    mIdxPicoKaons.shrink_to_fit();

    tracksToRemove.clear();
    tracksToRemove.shrink_to_fit();

    return kStOK;
}

//____________________________________________________________________________________
TVector3 StPicoD0AnaMaker::refitVertex(bool always){
    bool pairRem=true;
    bool singleTrack=true;
//    bool singleTrack=!pairRem;
    std::vector<int> goodTracksToFit;

    TH1F *hPVDiffX = static_cast<TH1F*>(mOutList->FindObject("hPVDiffX"));
    TH1F *hPVDiffY = static_cast<TH1F*>(mOutList->FindObject("hPVDiffY"));
    TH1F *hPVDiffZ = static_cast<TH1F*>(mOutList->FindObject("hPVDiffZ"));
    TH1F *hPVDiffXRemoved = static_cast<TH1F*>(mOutList->FindObject("hPVDiffXRemoved"));
    TH1F *hPVDiffYRemoved = static_cast<TH1F*>(mOutList->FindObject("hPVDiffYRemoved"));
    TH1F *hPVDiffZRemoved = static_cast<TH1F*>(mOutList->FindObject("hPVDiffZRemoved"));

    TH1F *hRemovedPairMass = static_cast<TH1F*>(mOutList->FindObject("hRemovedPairMass"));
    TH1F *hInvMass[3] = {static_cast<TH1F*>(mOutList->FindObject("hInvMassBDT12")), static_cast<TH1F*>(mOutList->FindObject("hInvMassBDT23")), static_cast<TH1F*>(mOutList->FindObject("hInvMassBDT35"))};

    bool isRemovedtrack=false;

    //removing according single track cuts

    if (singleTrack) {
        float dca;
        for (unsigned short iTrack = 0; iTrack < primaryTracks.size(); ++iTrack) {
            StPicoTrack* trk = mPicoDst->track(primaryTracks[iTrack]);
                dca = (mPrimVtx - trk->origin()).Mag();
                if (dca > 0.009) tracksToRemove.push_back(primaryTracks[iTrack]);
        }
    }


    //removing with pair cuts
    if (pairRem) {
        for (unsigned short idxPion1 = 0; idxPion1 < mIdxPicoPions.size(); ++idxPion1) {
            StPicoTrack const *pion1 = mPicoDst->track(mIdxPicoPions[idxPion1]);
            for (unsigned short idxKaon = 0; idxKaon < mIdxPicoKaons.size(); ++idxKaon) {
                StPicoTrack const *kaon = mPicoDst->track(mIdxPicoKaons[idxKaon]);

                if(kaon->charge()+pion1->charge() != 0) continue;
                StHFPair *pair = new StHFPair(pion1, kaon, mHFCuts->getHypotheticalMass(StPicoCutsBase::kPion), mHFCuts->getHypotheticalMass(StPicoCutsBase::kKaon), mIdxPicoPions[idxPion1], mIdxPicoKaons[idxKaon], mPrimVtx, mBField, kTRUE);
                if(!(pair->m()>1.81 && pair->m()<1.92)) continue;

                //CUTS

                if (cos(pair->pointingAngle()) > 0.9 && pair->dcaDaughters() < 0.007 && pair->m()>1.81 && pair->m()<1.92) {
                    tracksToRemove.push_back(mIdxPicoPions[idxPion1]);
                    tracksToRemove.push_back(mIdxPicoKaons[idxKaon]);
                    hRemovedPairMass->Fill(pair->m());
                    isRemovedtrack = true;
                }


                //BDT
                //find the correct pT bin
                /*
                int pTbin = 0;
                for (int pT = 0; pT < nptBins; pT++) {
                    if(pair->pt() >= momBins[pT] && pair->pt() < momBins[pT+1]) pTbin = pT;
                }

                k_dca[pTbin] = pair->particle2Dca();
                pi1_dca[pTbin] = pair->particle1Dca();
                D_decayL[pTbin] = pair->decayLength();
                cosTheta[pTbin] = cos(pair->pointingAngle());
                dcaD0ToPv[pTbin] = pair->DcaToPrimaryVertex();
                dcaDaughters[pTbin] = pair->dcaDaughters();
//                thetaStar[pTbin] = pair->cosThetaStar();

                //evaluate BDT, continue just pairs that have passed BDT cut
                float valueMVA = reader[pTbin]->EvaluateMVA("BDT method");
                if(valueMVA > bdtCuts[pTbin]) {
                    //filling plots of invariant mass for unlike and like sign pairs
                    hInvMass[pTbin]->Fill(pair->m());
                    //excluding daughter tracks
                    tracksToRemove.push_back(mIdxPicoPions[idxPion1]);
                    tracksToRemove.push_back(mIdxPicoKaons[idxKaon]);

                }*/
            }
        }
    }

    for (int iTrk = 0; iTrk < primaryTracks.size(); ++iTrk) {
        if(std::binary_search(tracksToRemove.begin(), tracksToRemove.end(), primaryTracks[iTrk])) continue;
        goodTracksToFit.push_back(primaryTracks[iTrk]);
    }

    TVector3 newKFVertex(mPrimVtx.x(), mPrimVtx.y(), mPrimVtx.z());
    if (goodTracksToFit.size()>2) {
//    if (isRemovedtrack || always) {
        //Make new vertex and evaluate stuff:
        StPicoKFVertexFitter kfVertexFitter;
//        KFVertex kfVertex = kfVertexFitter.primaryVertexRefit(mPicoDst, tracksToRemove); //when removed tracks need to be checked
        KFVertex kfVertex = kfVertexFitter.primaryVertexRefitUsingTracks(mPicoDst, goodTracksToFit); //when you have array with tracks for refit
        if (kfVertex.GetX()) {
            newKFVertex.SetXYZ(kfVertex.GetX(), kfVertex.GetY(), kfVertex.GetZ());
            hPVDiffX->Fill(mPrimVtx.x() - newKFVertex.x());
            hPVDiffY->Fill(mPrimVtx.y() - newKFVertex.y());
            hPVDiffZ->Fill(mPrimVtx.z() - newKFVertex.z());

            if (tracksToRemove.size() > 0) {
                hPVDiffXRemoved->Fill(mPrimVtx.x() - newKFVertex.x());
                hPVDiffYRemoved->Fill(mPrimVtx.y() - newKFVertex.y());
                hPVDiffZRemoved->Fill(mPrimVtx.z() - newKFVertex.z());
            }
        }
    }
    nGoodTracks=goodTracksToFit.size();

    goodTracksToFit.clear();
    goodTracksToFit.shrink_to_fit();

    return newKFVertex;
}

//____________________________________________________________________________________
int StPicoD0AnaMaker::analyzeCandidates() {
    TH1D *h_tracktest = static_cast<TH1D*>(mOutList->FindObject("h_tracktest"));
    TH1D *h_tracktest_TOF = static_cast<TH1D*>(mOutList->FindObject("h_tracktest_TOF"));

    const char *aNames[]   = {"all", "NHitFit>20", "pT>0.6", "dca<1.5","TPC pion", "TPC kaon"};
    const char *aNamesTOF[]   = {"all TOF match", "NHitFit>20 & TOF", "pT>0.6 & TOF","dca<1.5 & TOF","TPC & TOF pion", "TPC & TOF kaon"};

    for (unsigned int ii = 0; ii < 6; ii++) {
        h_tracktest->GetXaxis()->SetBinLabel(ii+1, aNames[ii]);
        h_tracktest_TOF->GetXaxis()->SetBinLabel(ii+1, aNamesTOF[ii]);
    }

    for(unsigned int i=0;i<mPicoDst->numberOfTracks();i++) {
        StPicoTrack const *t = mPicoDst->track(i);
        if (!t) continue;
        if (!t->isHFTTrack()) continue;
        h_tracktest->Fill(1);
        if (mHFCuts->isGoodTrack(t)) h_tracktest->Fill(2); // NhitsFit

        float pt=t->gPt();
        if ((pt>0.6) && mHFCuts->isGoodTrack(t)) h_tracktest->Fill(3);
        float dca = (mPrimVtx - t->origin()).Mag();
        if (dca<1.5 && (pt>0.6) && mHFCuts->isGoodTrack(t)) h_tracktest->Fill(4);
        bool tpcPion = mHFCuts->isTPCHadron(t, StPicoCutsBase::kPion);
        bool tpcKaon = mHFCuts->isTPCHadron(t, StPicoCutsBase::kKaon);
        if ((pt>0.6) && (dca<1.5 ) && (mHFCuts->isGoodTrack(t))) {
            if (tpcPion) h_tracktest->Fill(5);
            if (tpcKaon) h_tracktest->Fill(6);
        }
        if (mHFCuts->isTOFmatched(t)) {
            h_tracktest_TOF->Fill(1);
            if (mHFCuts->isGoodTrack(t)) h_tracktest_TOF->Fill(2); // NhitsFit
            if (pt>0.6 && mHFCuts->isGoodTrack(t)) h_tracktest_TOF->Fill(3);
            if (dca<1.5 && pt>0.6 && mHFCuts->isGoodTrack(t)) h_tracktest_TOF->Fill(4);
            if ((pt>0.6) && (dca<1.5 ) && (mHFCuts->isGoodTrack(t))) {
//                if (tpcPion && mHFCuts->isTOFHadronPID(t, mHFCuts->getTofBetaBase(t), StPicoCutsBase::kPion))  h_tracktest_TOF->Fill(5);
                if (tpcPion)  h_tracktest_TOF->Fill(5);
                if (tpcKaon && mHFCuts->isTOFHadronPID(t, mHFCuts->getTofBetaBase(t), StPicoCutsBase::kKaon))  h_tracktest_TOF->Fill(6);
            }
        }


    }
//    for (unsigned short idxPion1 = 0; idxPion1 < mIdxPicoPions.size(); ++idxPion1) {
//        StPicoTrack const *t = mPicoDst->track(mIdxPicoPions[idxPion1]);
//        ntp_pion->Fill(t->gPt(), t->gMom().phi(), t->gMom().pseudoRapidity(), t->nSigmaPion(), t->nHitsFit(), getOneOverBeta(t, mHFCuts->getTofBetaBase(t), StPicoCutsBase::kPion), mPicoEvent->eventId(), mPicoEvent->runId());
//    }
//
//    for (unsigned short idxKaon = 0; idxKaon < mIdxPicoKaons.size(); ++idxKaon) {
//        StPicoTrack const *t = mPicoDst->track(mIdxPicoKaons[idxKaon]);
//        ntp_kaon->Fill(t->gPt(), t->gMom().phi(), t->gMom().pseudoRapidity(), t->nSigmaKaon(), t->nHitsFit(), getOneOverBeta(t, mHFCuts->getTofBetaBase(t), StPicoCutsBase::kKaon), mPicoEvent->eventId(), mPicoEvent->runId());
//    }
    return kStOK;
}