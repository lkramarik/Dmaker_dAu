#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoDstMaker/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
//#include "StPicoD0EventMaker/StPicoD0Event.h"
//#include "StPicoD0EventMaker/StKaonPion.h"
//#include "StPicoD0AnaMaker.h"
#include "StPicoHFMaker/StHFCuts.h"

#include "StPicoD0AnaMaker.h"
ClassImp(StPicoD0AnaMaker)

// _________________________________________________________
StPicoD0AnaMaker::StPicoD0AnaMaker(char const* name, StPicoDstMaker* picoMaker, char const* outputBaseFileName, char const* inputHFListHFtree = "") :
        StPicoHFMaker(name, picoMaker, outputBaseFileName, inputHFListHFtree),
        mOutFileBaseName(outputBaseFileName){ //mDecayChannel(kChannel1), tu bolo

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
    mOutFileBaseName = mOutFileBaseName.ReplaceAll(".root", "");

    ntp_DMeson_Signal = new TNtuple("ntp_signal","DMeson TreeSignal","grefMult:pi1_runId:pi1_eventId:pi1_phi:pi1_eta:pi1_pt:pi1_dca:pi1_dedx:pi1_nSigma:pi1_nHitFit:pi1_nHitdedx:pi1_TOFinvbeta:pi1_betaBase:k_runId:k_eventId:k_phi:k_eta:k_pt:k_dca:k_dedx:k_nSigma:k_nHitFit:k_nHitdedx:k_TOFinvbeta:k_betaBase:dcaDaughters:flag:primVz:D_theta:cosTheta:D_decayL:dcaD0ToPv:D_phi:D_eta:D_cosThetaStar:D_pt:D_mass:D_mass_LS:D_mass_US:centrality:refmultcorr:reweight");
    ntp_DMeson_Background = new TNtuple("ntp_background","DMeson TreeBackground","grefMult:pi1_runId:pi1_eventId:pi1_phi:pi1_eta:pi1_pt:pi1_dca:pi1_dedx:pi1_nSigma:pi1_nHitFit:pi1_nHitdedx:pi1_TOFinvbeta:pi1_betaBase:k_runId:k_eventId:k_phi:k_eta:k_pt:k_dca:k_dedx:k_nSigma:k_nHitFit:k_nHitdedx:k_TOFinvbeta:k_betaBase:dcaDaughters:flag:primVz:D_theta:cosTheta:D_decayL:dcaD0ToPv:D_phi:D_eta:D_cosThetaStar:D_pt:D_mass:D_mass_LS:D_mass_US:centrality:refmultcorr:reweight");

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
    //cout<<"MakeHF start"<<endl;  
    // -- process event
    //    ADD YOUR PROCESSING CODE HERE
    //    ... it is usefull to use the methods below
    //     - createCandidates()
    //     - analyzeCandidates()

    std::clock_t start1 = std::clock();//kvapil
    if (isMakerMode() == StPicoHFMaker::kWrite) {
        createCandidates();
    }
    else if (isMakerMode() == StPicoHFMaker::kRead) {
        // -- the reading back of the perviously written trees happens in the background
        analyzeCandidates();
    }
    else if (isMakerMode() == StPicoHFMaker::kAnalyze) {
        //cout<<"going to create candidates"<<endl;
        createCandidates();
        //cout<<"candidated created, going to analyze them"<<endl;
        analyzeCandidates();
        //cout<<"candidates analysed"<<endl;
        //createQA();
    }

    double duration = (double) (std::clock() - start1) / (double) CLOCKS_PER_SEC;
    return kStOK;
}

// _________________________________________________________
int StPicoD0AnaMaker::createCandidates() {

    for (unsigned short idxPion1 = 0; idxPion1 < mIdxPicoPions.size(); ++idxPion1) {
        StPicoTrack const *pion1 = mPicoDst->track(mIdxPicoPions[idxPion1]);
        //     if( !mHFCuts->isHybridTOFHadron(pion1, mHFCuts->getTofBetaBase(pion1), StHFCuts::kPion) ) continue;
        for (unsigned short idxKaon = 0; idxKaon < mIdxPicoKaons.size(); ++idxKaon) {
            StPicoTrack const *kaon = mPicoDst->track(mIdxPicoKaons[idxKaon]);
            //        if( !mHFCuts->isHybridTOFHadron(kaon, mHFCuts->getTofBetaBase(kaon), StHFCuts::kKaon) ) continue;
            if (mIdxPicoKaons[idxKaon] == mIdxPicoPions[idxPion1]) continue;
            // -- Making pair
            StHFPair pair(pion1, kaon, mHFCuts->getHypotheticalMass(StHFCuts::kPion),mHFCuts->getHypotheticalMass(StHFCuts::kKaon), mIdxPicoPions[idxPion1],mIdxPicoKaons[idxKaon], mPrimVtx, mBField, kTRUE);
            if (!mHFCuts->isClosePair(pair)) continue;
            mPicoHFEvent->addHFSecondaryVertexPair(&pair);
        }  // for (unsigned short idxKaon = 0; idxKaon < mIdxPicoKaons.size(); ++idxKaon)
    } // for (unsigned short idxPion1 = 0; idxPion1 < mIdxPicoPions.size(); ++idxPion1)

    return kStOK;
}

// _________________________________________________________
int StPicoD0AnaMaker::analyzeCandidates() {
    // --- Analyze previously constructed candidates and output to ntuple
    TClonesArray const * aCandidates= mPicoHFEvent->aHFSecondaryVertices();
    if( mPicoHFEvent->nHFSecondaryVertices() > 0 ){
        for (unsigned int idx = 0; idx <  mPicoHFEvent->nHFSecondaryVertices(); ++idx) {
            StHFPair const* pair = static_cast<StHFPair*>(aCandidates->At(idx));
            StPicoTrack const* pion1 = mPicoDst->track(pair->particle1Idx());
            StPicoTrack const* kaon = mPicoDst->track(pair->particle2Idx());

            //TOF ---
            float kaonTOFinvbeta = -999;
            float pion1TOFinvbeta = -999;
            float kaonBetaBase = -999;
            float pion1BetaBase = -999;
            kaonBetaBase = mHFCuts->getTofBetaBase(kaon);
            pion1BetaBase = mHFCuts->getTofBetaBase(pion1);

            // all of the tracks need to have TOF info
            if(isnan(kaonBetaBase) && kaonBetaBase == -999) continue;
            if(isnan(pion1BetaBase) && pion1BetaBase == -999) continue;

            float ptot=9999;
            float betaInv = 9999;


            if(!isnan(kaonBetaBase) && kaonBetaBase > -999){
                ptot = kaon->gPtot();
                betaInv = sqrt(ptot*ptot + M_KAON_PLUS*M_KAON_PLUS) / ptot;
                kaonTOFinvbeta = fabs(1/kaonBetaBase - betaInv);
            }

            if(!isnan(pion1BetaBase) && pion1BetaBase > -999){
                ptot = pion1->gPtot();
                betaInv = sqrt(ptot*ptot + M_PION_PLUS*M_PION_PLUS) / ptot;
                pion1TOFinvbeta = fabs(1/pion1BetaBase - betaInv);
            }

//            if(!isnan(pion1BetaBase) && pion1BetaBase > 0){
//                ptot = pion1->gPtot();
//                pion1TOFinvbeta = fabs(1. / pion1BetaBase - sqrt(1+M_PION_PLUS*M_PION_PLUS/(pion1->gMom(mPrimVtx,mBField).mag()*pion1->gMom(mPrimVtx,mBField).mag())));
//            }

            // -- Flag D0 and background
            float flag = -99.;
            if( kaon->charge()<0 && pion1->charge()>0 ) flag=0.; // -+
            if( kaon->charge()>0 && pion1->charge()<0 ) flag=1.; // +-

            if( kaon->charge()<0 && pion1->charge()<0) flag=4.; // --
            if( kaon->charge()>0 && pion1->charge()>0) flag=5.; // ++

            int const centrality = 1;
            const double reweight = 1;
            const double refmultCor = 1;

            int ii=0;
            float ntVar[42];
            // Saving to NTUPLE
            // float globalTracks = (float)(mPicoHFEvent->numberOfGlobalTracks());
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
            ntVar[ii++] = pion1TOFinvbeta;
            ntVar[ii++] = pion1BetaBase;

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
            ntVar[ii++] = kaonTOFinvbeta;
            ntVar[ii++] = kaonBetaBase;

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
            ntVar[ii++] = centrality;
            ntVar[ii++] = refmultCor;
            ntVar[ii++] = reweight;

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
    StThreeVectorF t = trk->gMom(); //11.12.
    if (fabs(t.pseudoRapidity()) > 1.) return false; //pridano fabs 1212
    if (!mHFCuts->isHybridTOFHadron(trk, mHFCuts->getTofBetaBase(trk), StHFCuts::kPion) ) return false;
    if (!mHFCuts->cutMinDcaToPrimVertex(trk, StPicoCutsBase::kPion)) return false;
    return (mHFCuts->isGoodTrack(trk) && mHFCuts->isTPCHadron(trk, StPicoCutsBase::kPion));
}

// _________________________________________________________
bool StPicoD0AnaMaker::isKaon(StPicoTrack const * const trk) const {
    // -- good kaon
//    StThreeVectorF t = trk->pMom(); //11.12.
    StThreeVectorF t = trk->gMom(); //11.12.
    if (fabs(t.pseudoRapidity()) > 1.) return false;//pridano fabs 1212
    if (!mHFCuts->isHybridTOFHadron(trk, mHFCuts->getTofBetaBase(trk), StHFCuts::kKaon) ) return false;
    if (!mHFCuts->cutMinDcaToPrimVertex(trk, StPicoCutsBase::kKaon)) return false;
    return (mHFCuts->isGoodTrack(trk) && mHFCuts->isTPCHadron(trk, StPicoCutsBase::kKaon));
}

// _________________________________________________________
bool StPicoD0AnaMaker::isProton(StPicoTrack const * const trk) const {
    // -- good proton
    return (mHFCuts->isGoodTrack(trk) && mHFCuts->isTPCHadron(trk, StPicoCutsBase::kProton));
}

bool StPicoD0AnaMaker::isCloseTracks(StPicoTrack const * const trk1, StPicoTrack const * const trk2, StThreeVectorF const & vtx, float bField) const {

    StPhysicalHelixD p1Helix = trk1->helix(mPicoDst->event()->bField());
    StPhysicalHelixD p2Helix = trk2->helix(mPicoDst->event()->bField());
    p1Helix.moveOrigin(p1Helix.pathLength(vtx));
    p2Helix.moveOrigin(p2Helix.pathLength(vtx));
    if( ( p1Helix.origin()-vtx ).mag()>0.2 || ( p2Helix.origin()-vtx ).mag()>0.2 ) return false;

    //Requires loading constants
    StThreeVectorF const p1Mom = p1Helix.momentum(bField * kilogauss);
    StThreeVectorF const p2Mom = p2Helix.momentum(bField * kilogauss);
    StPhysicalHelixD const p1StraightLine(p1Mom, p1Helix.origin(), 0, trk1->charge());
    StPhysicalHelixD const p2StraightLine(p2Mom, p2Helix.origin(), 0, trk2->charge());
    //DCA
    pair<double, double> const ss = p1StraightLine.pathLengths(p2StraightLine);
    StThreeVectorF const p1AtDcaToP2 = p1StraightLine.at(ss.first);
    StThreeVectorF const p2AtDcaToP1 = p2StraightLine.at(ss.second);
    float const dca = (p1AtDcaToP2-p2AtDcaToP1).mag();
    if(dca > 0.50) return false; //kubo : 0.009, in cm
    // -- good pair
    return true;
}
