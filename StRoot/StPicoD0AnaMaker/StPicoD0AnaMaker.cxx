#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoDstMaker/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
//#include "StPicoD0EventMaker/StPicoD0Event.h"
//#include "StPicoD0EventMaker/StKaonPion.h"
#include "StPicoHFMaker/StHFCuts.h"

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
    mOutList->Add(new TH2F("h_piTOF","h_piTOF",100,0,10, 250, -1, 1.5));
    mOutList->Add(new TH2F("h_kTOF","h_kTOF",100,0,10, 250, -1, 1.5));
    mOutList->Add(new TH2F("h_pTOF","h_pTOF",100,0,10, 250, -1, 1.5));

    mOutList->Add(new TH2F("h_piTOF_20","h_piTOF_20",100,0,10, 300, 0, 1));
    mOutList->Add(new TH2F("h_kTOF_20","h_kTOF_20",100,0,10, 300, 0, 1));
    mOutList->Add(new TH2F("h_pTOF_20","h_pTOF_20",100,0,10, 300, 0, 1));

    mOutList->Add(new TH2F("h_piTOF_HFT","h_piTOF_HFT",100,0,10, 300, 0, 1));
    mOutList->Add(new TH2F("h_kTOF_HFT","h_kTOF_HFT",100,0,10, 300, 0, 1));
    mOutList->Add(new TH2F("h_pTOF_HFT","h_pTOF_HFT",100,0,10, 300, 0, 1));

    mOutList->Add(new TH2F("h_piTOF_HFT_20","h_piTOF_HFT_20",100,0,10, 300, 0, 1));
    mOutList->Add(new TH2F("h_kTOF_HFT_20","h_kTOF_HFT_20",100,0,10, 300, 0, 1));
    mOutList->Add(new TH2F("h_pTOF_HFT_20","h_pTOF_HFT_20",100,0,10, 300, 0, 1));

    mOutList->Add(new TH2F("h_piTOFbeta","h_piTOFbeta",500,0,10, 300, 0, 1));
    mOutList->Add(new TH2F("h_kTOFbeta","h_kTOFbeta",500,0,10, 300, 0, 1));
    mOutList->Add(new TH2F("h_pTOFbeta","h_pTOFbeta",500,0,10, 300, 0, 1));

    mOutList->Add(new TH2F("h_pinsigma","h_pinsigma",1000,0,10, 99, -5, 5));
    mOutList->Add(new TH2F("h_knsigma","h_knsigma",1000,0,10, 99, -5, 5));
    mOutList->Add(new TH2F("h_pnsigma","h_pnsigma",1000,0,10, 99, -5, 5));

    mOutList->Add(new TH2F("h_dedx","h_dedx", 1000, 0, 10, 1000, 0, 10));

    mOutFileBaseName = mOutFileBaseName.ReplaceAll(".root", "");

    ntp_DMeson_Signal = new TNtuple("ntp_signal","DMeson TreeSignal","grefMult:pi1_runId:pi1_eventId:pi1_phi:pi1_eta:pi1_pt:pi1_dca:pi1_dedx:pi1_nSigma:pi1_nHitFit:pi1_nHitdedx:pi1_TOFinvbeta:pi1_betaBase:k_runId:k_eventId:k_phi:k_eta:k_pt:k_dca:k_dedx:k_nSigma:k_nHitFit:k_nHitdedx:k_TOFinvbeta:k_betaBase:dcaDaughters:flag:primVz:D_theta:cosTheta:D_decayL:dcaD0ToPv:D_phi:D_eta:D_cosThetaStar:D_pt:D_mass:D_mass_LS:D_mass_US");
    ntp_DMeson_Background = new TNtuple("ntp_background","DMeson TreeBackground","grefMult:pi1_runId:pi1_eventId:pi1_phi:pi1_eta:pi1_pt:pi1_dca:pi1_dedx:pi1_nSigma:pi1_nHitFit:pi1_nHitdedx:pi1_TOFinvbeta:pi1_betaBase:k_runId:k_eventId:k_phi:k_eta:k_pt:k_dca:k_dedx:k_nSigma:k_nHitFit:k_nHitdedx:k_TOFinvbeta:k_betaBase:dcaDaughters:flag:primVz:D_theta:cosTheta:D_decayL:dcaD0ToPv:D_phi:D_eta:D_cosThetaStar:D_pt:D_mass:D_mass_LS:D_mass_US");

    return kStOK;
}

// _________________________________________________________
void StPicoD0AnaMaker::ClearHF(Option_t *opt="") {
    return;
}

// _________________________________________________________
int StPicoD0AnaMaker::FinishHF() {
    if( isMakerMode() != StPicoHFMaker::kWrite ){

        ntp_DMeson_Signal -> Write(ntp_DMeson_Signal->GetName(), TObject::kOverwrite);
//        ntp_DMeson_Signal -> Write();
        ntp_DMeson_Background -> Write(ntp_DMeson_Background->GetName(), TObject::kOverwrite);
//        ntp_DMeson_Background -> Write();
    }
//    mOutFile->Close();
    return kStOK;
}
// _________________________________________________________
int StPicoD0AnaMaker::MakeHF() {
//    createCandidates();
//        analyzeCandidates();

    TH2F *h_piTOF = static_cast<TH2F*>(mOutList->FindObject("h_piTOF"));
    TH2F *h_kTOF = static_cast<TH2F*>(mOutList->FindObject("h_kTOF"));
    TH2F *h_pTOF = static_cast<TH2F*>(mOutList->FindObject("h_pTOF"));

    TH2F *h_piTOF_20 = static_cast<TH2F*>(mOutList->FindObject("h_piTOF_20"));
    TH2F *h_kTOF_20 = static_cast<TH2F*>(mOutList->FindObject("h_kTOF_20"));
    TH2F *h_pTOF_20 = static_cast<TH2F*>(mOutList->FindObject("h_pTOF_20"));

    TH2F *h_piTOF_HFT = static_cast<TH2F*>(mOutList->FindObject("h_piTOF_HFT"));
    TH2F *h_kTOF_HFT = static_cast<TH2F*>(mOutList->FindObject("h_kTOF_HFT"));
    TH2F *h_pTOF_HFT = static_cast<TH2F*>(mOutList->FindObject("h_pTOF_HFT"));

    TH2F *h_piTOF_HFT_20 = static_cast<TH2F*>(mOutList->FindObject("h_piTOF_HFT_20"));
    TH2F *h_kTOF_HFT_20 = static_cast<TH2F*>(mOutList->FindObject("h_kTOF_HFT_20"));
    TH2F *h_pTOF_HFT_20 = static_cast<TH2F*>(mOutList->FindObject("h_pTOF_HFT_20"));

    TH2F *h_piTOFbeta = static_cast<TH2F*>(mOutList->FindObject("h_piTOFbeta"));
    TH2F *h_kTOFbeta = static_cast<TH2F*>(mOutList->FindObject("h_kTOFbeta"));
    TH2F *h_pTOFbeta = static_cast<TH2F*>(mOutList->FindObject("h_pTOFbeta"));

    TH2F *h_pinsigma = static_cast<TH2F*>(mOutList->FindObject("h_pinsigma"));
    TH2F *h_knsigma = static_cast<TH2F*>(mOutList->FindObject("h_knsigma"));
    TH2F *h_pnsigma = static_cast<TH2F*>(mOutList->FindObject("h_pnsigma"));

    TH2F *h_dedx = static_cast<TH2F*>(mOutList->FindObject("h_dedx"));

    StThreeVectorF pVtx = mPicoDst->event()->primaryVertex();

    UInt_t nTracks = mPicoDst->numberOfTracks();
    for (unsigned short iTrack = 0; iTrack < nTracks; ++iTrack){
        StPicoTrack const* trk = mPicoDst->track(iTrack);
        if (!trk) continue;
        StPhysicalHelixD helix = trk->helix(mBField);
        StThreeVectorF momentum = trk->gMom(pVtx, mPicoDst->event()->bField());

        if (!(trk->nHitsFit()>=15)) continue;
        if (!(fabs(momentum.pseudoRapidity()) <= 1.0)) continue;

        if (fabs(trk->nSigmaPion())<3.0){
            if (mHFCuts->isTOFHadron(trk, mHFCuts->getTofBetaBase(trk), StHFCuts::kPion)) {
                    h_piTOF->Fill(trk->gPt(),mHFCuts->getTofBetaBase(trk));
                    float oneOverBeta = getOneOverBeta(trk, mHFCuts->getTofBetaBase(trk), StHFCuts::kPion);
                    h_piTOFbeta->Fill(trk->gPt(),oneOverBeta);
                    if (trk->nHitsFit()>=20) h_piTOF_20->Fill(trk->gPt(),oneOverBeta);
                    if (trk->isHFTTrack()) h_piTOF_HFT->Fill(trk->gPt(),oneOverBeta);
                    if ((trk->isHFTTrack()) && (trk->nHitsFit()>=20)) h_piTOF_HFT_20->Fill(trk->gPt(),oneOverBeta);
            }
        }

        if (fabs(trk->nSigmaKaon())<3.0){
            if (mHFCuts->isTOFHadron(trk, mHFCuts->getTofBetaBase(trk), StHFCuts::kKaon)) {
                h_kTOF->Fill(trk->gPt(),mHFCuts->getTofBetaBase(trk));
                float oneOverBeta = getOneOverBeta(trk, mHFCuts->getTofBetaBase(trk), StHFCuts::kKaon);
                h_kTOFbeta->Fill(trk->gPt(),oneOverBeta);
                if (trk->nHitsFit()>=20) h_kTOF_20->Fill(trk->gPt(),oneOverBeta);
                if (trk->isHFTTrack()) h_kTOF_HFT->Fill(trk->gPt(),oneOverBeta);
                if ((trk->isHFTTrack()) && (trk->nHitsFit()>=20)) h_kTOF_HFT_20->Fill(trk->gPt(),oneOverBeta);
            }
        }

        if (fabs(trk->nSigmaProton())<3.0){
            if (mHFCuts->isTOFHadron(trk, mHFCuts->getTofBetaBase(trk), StHFCuts::kProton)) {
                h_pTOF->Fill(trk->gPt(),mHFCuts->getTofBetaBase(trk));
                float oneOverBeta = getOneOverBeta(trk, mHFCuts->getTofBetaBase(trk), StHFCuts::kProton);
                h_pTOFbeta->Fill(trk->gPt(),oneOverBeta);
                if (trk->nHitsFit()>=20) h_pTOF_20->Fill(trk->gPt(),oneOverBeta);
                if (trk->isHFTTrack()) h_pTOF_HFT->Fill(trk->gPt(),oneOverBeta);
                if ((trk->isHFTTrack()) && (trk->nHitsFit()>=20)) h_pTOF_HFT_20->Fill(trk->gPt(),oneOverBeta);
            }
        }

        h_pinsigma->Fill(momentum.mag(),trk->nSigmaPion());
        h_knsigma->Fill(momentum.mag(),trk->nSigmaKaon());
        h_pnsigma->Fill(momentum.mag(),trk->nSigmaProton());
        h_dedx->Fill(momentum.mag(),trk->dEdx());

    } // .. end tracks loop

   return kStOK;
}

// _________________________________________________________
int StPicoD0AnaMaker::createCandidates() {
    for (unsigned short idxPion1 = 0; idxPion1 < mIdxPicoPions.size(); ++idxPion1) {
        StPicoTrack const *pion1 = mPicoDst->track(mIdxPicoPions[idxPion1]);
        for (unsigned short idxKaon = 0; idxKaon < mIdxPicoKaons.size(); ++idxKaon) {
            StPicoTrack const *kaon = mPicoDst->track(mIdxPicoKaons[idxKaon]);
            if (mIdxPicoKaons[idxKaon] == mIdxPicoPions[idxPion1]) continue;
            StHFPair *pair = new StHFPair(pion1, kaon, mHFCuts->getHypotheticalMass(StHFCuts::kPion),mHFCuts->getHypotheticalMass(StHFCuts::kKaon), mIdxPicoPions[idxPion1],mIdxPicoKaons[idxKaon], mPrimVtx, mBField, kTRUE);

            if (!mHFCuts->isClosePair(pair)) continue;
//            Filling ntp
            if(pair->pt() < 1) continue;
            if(pair->pt() > 2) continue;

            float flag = -99.;

            if( kaon->charge()<0 && pion1->charge()>0 ) flag=0.; // -+
            if( kaon->charge()>0 && pion1->charge()<0 ) flag=1.; // +-

            if( kaon->charge()<0 && pion1->charge()<0) flag=4.; // --
            if( kaon->charge()>0 && pion1->charge()>0) flag=5.; // ++

            int ii=0;
            float ntVar[39];
            ntVar[ii++] = mPicoDst->event()->refMult();
            ntVar[ii++] = mPicoHFEvent->runId();
            ntVar[ii++] = mPicoHFEvent->eventId();
            ntVar[ii++] = pion1->gMom(mPrimVtx,mBField).phi();
            ntVar[ii++] = pion1->gMom(mPrimVtx,mBField).pseudoRapidity();
            ntVar[ii++] = pion1->gPt();
            ntVar[ii++] = pair->particle1Dca();
            ntVar[ii++] = pion1->dEdx();
            ntVar[ii++] = pion1->nSigmaPion();
            ntVar[ii++] = pion1->nHitsFit();
            ntVar[ii++] = pion1->nHitsDedx();
            ntVar[ii++] = getOneOverBeta(pion1, mHFCuts->getTofBetaBase(pion1), StHFCuts::kPion);
            ntVar[ii++] = mHFCuts->getTofBetaBase(pion1);

            ntVar[ii++] = mPicoHFEvent->runId();
            ntVar[ii++] = mPicoHFEvent->eventId();
            ntVar[ii++] = kaon->gMom(mPrimVtx,mBField).phi();
            ntVar[ii++] = kaon->gMom(mPrimVtx,mBField).pseudoRapidity();
            ntVar[ii++] = kaon->gPt();
            ntVar[ii++] = pair->particle2Dca();
            ntVar[ii++] = kaon->dEdx();
            ntVar[ii++] = kaon->nSigmaKaon();
            ntVar[ii++] = kaon->nHitsFit();
            ntVar[ii++] = kaon->nHitsDedx();
            ntVar[ii++] = getOneOverBeta(kaon, mHFCuts->getTofBetaBase(kaon), StHFCuts::kKaon);
            ntVar[ii++] = mHFCuts->getTofBetaBase(kaon);

            ntVar[ii++] = pair->dcaDaughters();

            ntVar[ii++] = flag;
            ntVar[ii++] = mPrimVtx.z();
            ntVar[ii++] = pair->pointingAngle();
            ntVar[ii++] = cos(pair->pointingAngle());
            ntVar[ii++] = pair->decayLength();
            ntVar[ii++] = pair->DcaToPrimaryVertex(); //(pair->decayLength())*sin(pair->pointingAngle());
            ntVar[ii++] = pair->phi();
            ntVar[ii++] = pair->eta();
            ntVar[ii++] = pair->cosThetaStar();

            ntVar[ii++] = pair->pt(); //sqrt(pow(pair->px(),2.0)+pow(pair->py(),2.0));
            ntVar[ii++] = pair->m();
            if ((flag == 0) || (flag == 1)) {
                ntVar[ii++] = -5; //D_mass_LS
                ntVar[ii++] = pair->m();//D_mass_US
            } else {
                ntVar[ii++] = pair->m(); //D_mass_LS
                ntVar[ii++] = -5;//D_mass_US
            }

            if ((flag == 0) || (flag == 1)) {
                ntp_DMeson_Signal->Fill(ntVar);
            } else {
                ntp_DMeson_Background->Fill(ntVar);
            }



        }  // for (unsigned short idxKaon = 0; idxKaon < mIdxPicoKaons.size(); ++idxKaon)
    } // for (unsigned short idxPion1 = 0; idxPion1 < mIdxPicoPions.size(); ++idxPion1)

    return kStOK;
}

// no using_________________________________________________________
int StPicoD0AnaMaker::analyzeCandidates() {
    // --- Analyze previously constructed candidates and output to ntuple
    TClonesArray const * aCandidates= mPicoHFEvent->aHFSecondaryVertices();
    if( mPicoHFEvent->nHFSecondaryVertices() > 0 ){
        for (unsigned int idx = 0; idx <  mPicoHFEvent->nHFSecondaryVertices(); ++idx) {
            StHFPair const* pair = static_cast<StHFPair*>(aCandidates->At(idx));
            StPicoTrack const* pion1 = mPicoDst->track(pair->particle1Idx());
            StPicoTrack const* kaon = mPicoDst->track(pair->particle2Idx());

            if(pair->pt() < 1) continue;
            if(pair->pt() > 2) continue;

            float flag = -99.;

            if( kaon->charge()<0 && pion1->charge()>0 ) flag=0.; // -+
            if( kaon->charge()>0 && pion1->charge()<0 ) flag=1.; // +-

            if( kaon->charge()<0 && pion1->charge()<0) flag=4.; // --
            if( kaon->charge()>0 && pion1->charge()>0) flag=5.; // ++

            int ii=0;
            float ntVar[39];
            ntVar[ii++] = mPicoDst->event()->refMult();
            ntVar[ii++] = mPicoHFEvent->runId();
            ntVar[ii++] = mPicoHFEvent->eventId();
            ntVar[ii++] = pion1->gMom(mPrimVtx,mBField).phi();
            ntVar[ii++] = pion1->gMom(mPrimVtx,mBField).pseudoRapidity();
            ntVar[ii++] = pion1->gMom(mPrimVtx,mBField).perp();
            ntVar[ii++] = pair->particle1Dca();
            ntVar[ii++] = pion1->dEdx();
            ntVar[ii++] = pion1->nSigmaPion();
            ntVar[ii++] = pion1->nHitsFit();
            ntVar[ii++] = pion1->nHitsDedx();
            ntVar[ii++] = getOneOverBeta(pion1, mHFCuts->getTofBetaBase(pion1), StHFCuts::kPion);
            ntVar[ii++] = mHFCuts->getTofBetaBase(pion1);

            ntVar[ii++] = mPicoHFEvent->runId();
            ntVar[ii++] = mPicoHFEvent->eventId();
            ntVar[ii++] = kaon->gMom(mPrimVtx,mBField).phi();
            ntVar[ii++] = kaon->gMom(mPrimVtx,mBField).pseudoRapidity();
            ntVar[ii++] = kaon->gMom(mPrimVtx,mBField).perp();
            ntVar[ii++] = pair->particle2Dca();
            ntVar[ii++] = kaon->dEdx();
            ntVar[ii++] = kaon->nSigmaKaon();
            ntVar[ii++] = kaon->nHitsFit();
            ntVar[ii++] = kaon->nHitsDedx();
            ntVar[ii++] = getOneOverBeta(kaon, mHFCuts->getTofBetaBase(kaon), StHFCuts::kKaon);
            ntVar[ii++] = mHFCuts->getTofBetaBase(kaon);

            ntVar[ii++] = pair->dcaDaughters();

            ntVar[ii++] = flag;
            ntVar[ii++] = mPrimVtx.z();
            ntVar[ii++] = pair->pointingAngle();
            ntVar[ii++] = cos(pair->pointingAngle());
            ntVar[ii++] = pair->decayLength();
            ntVar[ii++] = pair->DcaToPrimaryVertex(); //(pair->decayLength())*sin(pair->pointingAngle());
            ntVar[ii++] = pair->phi();
            ntVar[ii++] = pair->eta();
            ntVar[ii++] = pair->cosThetaStar();

            ntVar[ii++] = pair->pt(); //sqrt(pow(pair->px(),2.0)+pow(pair->py(),2.0));
            ntVar[ii++] = pair->m();
            if ((flag == 0) || (flag == 1)) {
                ntVar[ii++] = -5; //D_mass_LS
                ntVar[ii++] = pair->m();//D_mass_US
            } else {
                ntVar[ii++] = pair->m(); //D_mass_LS
                ntVar[ii++] = -5;//D_mass_US
            }

            if ((flag == 0) || (flag == 1)) {
                ntp_DMeson_Signal->Fill(ntVar);
            } else {
                ntp_DMeson_Background->Fill(ntVar);
            }
        } // for (unsigned int idx = 0; idx <  mPicoHFEvent->nHFSecondaryVertices(); ++idx) {
    }
    return kStOK;
}

// _________________________________________________________
bool StPicoD0AnaMaker::isHadron(StPicoTrack const * const trk, int pidFlag) const {
    // -- good hadron
    return (mHFCuts->isGoodTrack(trk) && mHFCuts->isTPCHadron(trk, pidFlag));
}

// _________________________________________________________
bool StPicoD0AnaMaker::isPion(StPicoTrack const * const trk) const {
    // -- good pion
//    StThreeVectorF t = trk->pMom(); //11.12.
    StThreeVectorF t = trk->gMom(mPrimVtx, mBField); //11.12.
    if (fabs(t.pseudoRapidity()) > 1.) return false; //pridano fabs 1212
    if (!mHFCuts->isTOFHadron(trk, mHFCuts->getTofBetaBase(trk), StHFCuts::kPion) ) return false;
//    if (!mHFCuts->isHybridTOFHadron(trk, mHFCuts->getTofBetaBase(trk), StHFCuts::kPion) ) return false;
    if (!mHFCuts->cutMinDcaToPrimVertex(trk, StPicoCutsBase::kPion)) return false;
    return (mHFCuts->isGoodTrack(trk) && mHFCuts->isTPCHadron(trk, StPicoCutsBase::kPion));
}

// _________________________________________________________
bool StPicoD0AnaMaker::isKaon(StPicoTrack const * const trk) const {
    // -- good kaon
//    StThreeVectorF t = trk->pMom(); //11.12.
    StThreeVectorF t = trk->gMom(mPrimVtx, mBField); //11.12.
    if (fabs(t.pseudoRapidity()) > 1.) return false;//pridano fabs 1212
    if (!mHFCuts->isTOFHadron(trk, mHFCuts->getTofBetaBase(trk), StHFCuts::kKaon) ) return false;
//    if (!mHFCuts->isHybridTOFHadron(trk, mHFCuts->getTofBetaBase(trk), StHFCuts::kKaon) ) return false;
    if (!mHFCuts->cutMinDcaToPrimVertex(trk, StPicoCutsBase::kKaon)) return false;
    return (mHFCuts->isGoodTrack(trk) && mHFCuts->isTPCHadron(trk, StPicoCutsBase::kKaon));
}

// _________________________________________________________
bool StPicoD0AnaMaker::isProton(StPicoTrack const * const trk) const {
    // -- good proton
    return (mHFCuts->isGoodTrack(trk) && mHFCuts->isTPCHadron(trk, StPicoCutsBase::kProton));
}

float StPicoD0AnaMaker::getOneOverBeta(StPicoTrack const * const trk,  float const & tofBeta, int pidFlag){
    if (tofBeta <= 0)
        return -5;
    float m2 = mHFCuts->getHypotheticalMass(pidFlag)*mHFCuts->getHypotheticalMass(pidFlag);
    float ptot    = trk->gPtot();
    float betaInv = sqrt(ptot*ptot + m2) / ptot;
    return fabs(1/tofBeta - betaInv);
}