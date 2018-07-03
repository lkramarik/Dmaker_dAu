#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoDstMaker/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoHFMaker/StHFCuts.h"
#include "phys_constants.h"
#include "StBTofUtil/tofPathLength.hh"
#include "StPicoD0AnaMaker.h"
ClassImp(StPicoD0AnaMaker)

// _________________________________________________________
StPicoD0AnaMaker::StPicoD0AnaMaker(char const* name, StPicoDstMaker* picoMaker, char const* outputBaseFileName, char const* inputHFListHFtree = "") :
        StPicoHFMaker(name, picoMaker, outputBaseFileName, inputHFListHFtree),
        mOutFileBaseName(outputBaseFileName){
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
//    mOutList->Add(new TH2F("h_pinsigma","h_pinsigma",1000,0,10, 99, -5, 5));
//    mOutList->Add(new TH2F("h_knsigma","h_knsigma",1000,0,10, 99, -5, 5));
//    mOutList->Add(new TH2F("h_pnsigma","h_pnsigma",1000,0,10, 99, -5, 5));
//
//    mOutList->Add(new TH2F("h_dedx","h_dedx", 1000, 0, 10, 1000, 0, 10));
//h_tracktest
    mOutList->Add(new TH1D("h_tracktest","h_tracktest", 6, 0.5, 6.5));
    mOutList->Add(new TH1D("h_tracktest_TOF","h_tracktest_TOF", 6, 0.5, 6.5));

    mOutFileBaseName = mOutFileBaseName.ReplaceAll(".root", "");

//    ntp_kaon = new TNtuple("ntp_kaon", "kaon tree","k_pt:k_phi:k_eta:k_nSigma:k_nHitFit:k_TOFinvbeta:pi_eventId:pi_runId");
//    ntp_pion = new TNtuple("ntp_pion", "pion tree","pi_pt:pi_phi:pi_eta:pi_nSigma:pi_nHitFit:pi_TOFinvbeta:k_eventId:k_runId");
    ntp_DMeson_Signal = new TNtuple("ntp_signal","DMeson TreeSignal","grefMult:runId:eventId:pi1_phi:pi1_eta:pi1_pt:pi1_dca:pi1_dedx:pi1_nSigma:pi1_nHitFit:pi1_nHitdedx:pi1_TOFinvbeta:pi1_betaBase:k_phi:k_eta:k_pt:k_dca:k_dedx:k_nSigma:k_nHitFit:k_nHitdedx:k_TOFinvbeta:k_betaBase:dcaDaughters:flag:primVz:primVzVpd:primVzDiff:D_rapidity:D_theta:cosTheta:D_decayL:dcaD0ToPv:D_phi:D_eta:D_cosThetaStar:D_pt:D_mass");
    ntp_DMeson_Background = new TNtuple("ntp_background","DMeson TreeBackground","grefMult:runId:eventId:pi1_phi:pi1_eta:pi1_pt:pi1_dca:pi1_dedx:pi1_nSigma:pi1_nHitFit:pi1_nHitdedx:pi1_TOFinvbeta:pi1_betaBase:k_phi:k_eta:k_pt:k_dca:k_dedx:k_nSigma:k_nHitFit:k_nHitdedx:k_TOFinvbeta:k_betaBase:dcaDaughters:flag:primVz:primVzVpd:primVzDiff:D_rapidity:D_theta:cosTheta:D_decayL:dcaD0ToPv:D_phi:D_eta:D_cosThetaStar:D_pt:D_mass");

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

//    TH2F *h_piTOF = static_cast<TH2F*>(mOutList->FindObject("h_piTOF"));
//    TH2F *h_kTOF = static_cast<TH2F*>(mOutList->FindObject("h_kTOF"));
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
//    StThreeVectorF pVtx = mPicoDst->event()->primaryVertex();
//
//    UInt_t nTracks = mPicoDst->numberOfTracks();
//    for (unsigned short iTrack = 0; iTrack < nTracks; ++iTrack){
//        StPicoTrack const* trk = mPicoDst->track(iTrack);
//        if (!trk) continue;
//        StPhysicalHelixD helix = trk->helix(mBField);
//        StThreeVectorF momentum = trk->gMom(pVtx, mPicoDst->event()->bField());
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
//        h_pinsigma->Fill(momentum.mag(),trk->nSigmaPion());
//        h_knsigma->Fill(momentum.mag(),trk->nSigmaKaon());
//        h_pnsigma->Fill(momentum.mag(),trk->nSigmaProton());
//        h_dedx->Fill(momentum.mag(),trk->dEdx());
//
//    } // .. end tracks loop

    return kStOK;
}

// _________________________________________________________
int StPicoD0AnaMaker::createCandidates() {
//    for (unsigned short idxPion1 = 0; idxPion1 < mIdxPicoPions.size(); ++idxPion1) {
//        StPicoTrack const *pion1 = mPicoDst->track(mIdxPicoPions[idxPion1]);
//        for (unsigned short idxKaon = 0; idxKaon < mIdxPicoKaons.size(); ++idxKaon) {
//            StPicoTrack const *kaon = mPicoDst->track(mIdxPicoKaons[idxKaon]);
//    StHFPair *pair = new StHFPair(pion1, kaon, mHFCuts->getHypotheticalMass(StPicoCutsBase::kPion),mHFCuts->getHypotheticalMass(StPicoCutsBase::kKaon), mIdxPicoPions[idxPion1],mIdxPicoKaons[idxKaon], mPrimVtx, mBField, kTRUE);

    for(unsigned int i=0;i<mPicoDst->numberOfTracks();i++)  {
        StPicoTrack const* pion1 = mPicoDst->track(i);
        if (!isPion(pion1)) continue;

        for(unsigned  int j=0;j<mPicoDst->numberOfTracks();j++)  {
            StPicoTrack const* kaon = mPicoDst->track(j);
            if (!isKaon(kaon)) continue;

            StHFPair *pair = new StHFPair(pion1, kaon, mHFCuts->getHypotheticalMass(StPicoCutsBase::kPion),mHFCuts->getHypotheticalMass(StPicoCutsBase::kKaon), i, j, mPrimVtx, mBField, kTRUE);
            if (!mHFCuts->isClosePair(pair)) continue;

//            if(pair->pt() < 1) continue;
//            if(pair->pt() > 2) continue;

            float flag = -99.;

            if( kaon->charge()<0 && pion1->charge()>0 ) flag=0.; // -+
            if( kaon->charge()>0 && pion1->charge()<0 ) flag=1.; // +-

            if( kaon->charge()<0 && pion1->charge()<0) flag=4.; // --
            if( kaon->charge()>0 && pion1->charge()>0) flag=5.; // ++

            int ii=0;
            float ntVar[38];
            ntVar[ii++] = mPicoDst->event()->refMult();
            ntVar[ii++] = mPicoEvent->runId();
            ntVar[ii++] = mPicoEvent->eventId();

            ntVar[ii++] = pion1->gMom(mPrimVtx,mBField).phi();
            ntVar[ii++] = pion1->gMom(mPrimVtx,mBField).pseudoRapidity();
            ntVar[ii++] = pion1->gPt();
            ntVar[ii++] = pair->particle1Dca();
            ntVar[ii++] = pion1->dEdx();
            ntVar[ii++] = pion1->nSigmaPion();
            ntVar[ii++] = pion1->nHitsFit();
            ntVar[ii++] = pion1->nHitsDedx();
            ntVar[ii++] = mHFCuts->getOneOverBeta(pion1, mHFCuts->getTofBetaBase(pion1), StPicoCutsBase::kPion);
            ntVar[ii++] = mHFCuts->getTofBetaBase(pion1);

            ntVar[ii++] = kaon->gMom().phi();
            ntVar[ii++] = kaon->gMom().pseudoRapidity();
            ntVar[ii++] = kaon->gPt();
            ntVar[ii++] = pair->particle2Dca();
            ntVar[ii++] = kaon->dEdx();
            ntVar[ii++] = kaon->nSigmaKaon();
            ntVar[ii++] = kaon->nHitsFit();
            ntVar[ii++] = kaon->nHitsDedx();
            ntVar[ii++] = mHFCuts->getOneOverBeta(kaon, mHFCuts->getTofBetaBase(kaon), StPicoCutsBase::kKaon);
            ntVar[ii++] = mHFCuts->getTofBetaBase(kaon);

            ntVar[ii++] = pair->dcaDaughters();
            ntVar[ii++] = flag;
            ntVar[ii++] = mPrimVtx.z();
            ntVar[ii++] = mPicoEvent->vzVpd();
            ntVar[ii++] = fabs(mPicoEvent->primaryVertex().z() - mPicoEvent->vzVpd());

            ntVar[ii++] = pair->rapidity();
            ntVar[ii++] = pair->pointingAngle();
            ntVar[ii++] = cos(pair->pointingAngle());
            ntVar[ii++] = pair->decayLength();
            ntVar[ii++] = pair->DcaToPrimaryVertex(); //(pair->decayLength())*sin(pair->pointingAngle());
            ntVar[ii++] = pair->phi();
            ntVar[ii++] = pair->eta();
            ntVar[ii++] = pair->cosThetaStar();

            ntVar[ii++] = pair->pt();
            ntVar[ii++] = pair->m();

            if ((flag == 0) || (flag == 1)) {
                ntp_DMeson_Signal->Fill(ntVar);
            } else {
                ntp_DMeson_Background->Fill(ntVar);
            }
        }  // for (unsigned short idxKaon = 0; idxKaon < mIdxPicoKaons.size(); ++idxKaon)
    } // for (unsigned short idxPion1 = 0; idxPion1 < mIdxPicoPions.size(); ++idxPion1)

    return kStOK;
}

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
        float dca = (mPrimVtx - t->dcaPoint()).mag();
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

bool StPicoD0AnaMaker::isHadron(StPicoTrack const * const trk, int pidFlag) const {
    return (mHFCuts->isGoodTrack(trk) && mHFCuts->isTPCHadron(trk, pidFlag));
}

bool StPicoD0AnaMaker::isPion(StPicoTrack const * const trk) const {
    if (!mHFCuts->cutMinDcaToPrimVertex(trk, StPicoCutsBase::kPion)) return false;
    if (!mHFCuts->isGoodTrack(trk)) return false; //HFT, NhitsFit, pt range
    float kBeta = mHFCuts->getTofBetaBase(trk);

//    if (!mHFCuts->isTOFmatched(trk)) return false;

    bool tof = false;
    bool tpc = false;
    bool tofAvailable = (kBeta > 0) && (kBeta == kBeta);
    if (mHFCuts->isTOFHadronPID(trk, kBeta, StPicoCutsBase::kPion)) tof = true; // 1/beta diff
    if (mHFCuts->isTPCHadron(trk, StPicoCutsBase::kPion)) tpc = true;

//    hybrid TOF
    bool goodPion = (tofAvailable && tof && tpc) || (!tofAvailable && tpc);
    return goodPion;
    //    if (!mHFCuts->isTOFHadronPID(trk, mHFCuts->getTofBetaBase(trk), StPicoCutsBase::kPion) ) return false; // 1/beta diff
//    if (!mHFCuts->isHybridTOFHadron(trk, mHFCuts->getTofBetaBase(trk), StPicoCutsBase::kPion) ) return false;
}

bool StPicoD0AnaMaker::isKaon(StPicoTrack const * const trk) const {
    if (!mHFCuts->cutMinDcaToPrimVertex(trk, StPicoCutsBase::kKaon)) return false;
    if (!mHFCuts->isGoodTrack(trk)) return false;
    float kBeta = mHFCuts->getTofBetaBase(trk);

//    if (!mHFCuts->isTOFmatched(trk)) return false;

    bool tof = false;
    bool tpc = false;
    bool tofAvailable = (kBeta > 0) && (kBeta == kBeta);
    if (mHFCuts->isTOFHadronPID(trk, kBeta, StPicoCutsBase::kKaon)) tof = true; // 1/beta diff
    if (mHFCuts->isTPCHadron(trk, StPicoCutsBase::kKaon)) tpc = true; //TPC NsigmaKaon

//    hybrid TOF
    bool goodKaon = (tofAvailable && tof && tpc) || (!tofAvailable && tpc);

//    LIANG:
//     bool goodKaon = (tofAvailable && tof) || (!tofAvailable && tpc);

//    LUKAS:
//    bool goodKaon = tof && tpc;
//    bool goodKaon = true;

    return goodKaon;
}

bool StPicoD0AnaMaker::isProton(StPicoTrack const * const trk) const {
    return (mHFCuts->isGoodTrack(trk) && mHFCuts->isTPCHadron(trk, StPicoCutsBase::kProton));
}


