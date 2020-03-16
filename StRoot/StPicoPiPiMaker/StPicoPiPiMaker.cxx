#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoHFMaker/StHFCuts.h"
#include "phys_constants.h"
#include "StPicoPiPiMaker.h"
ClassImp(StPicoPiPiMaker)

// _________________________________________________________
StPicoPiPiMaker::StPicoPiPiMaker(char const* name, StPicoDstMaker* picoMaker, char const* outputBaseFileName) :
        StPicoHFMaker(name, picoMaker, outputBaseFileName),
        mOutFileBaseName(outputBaseFileName){
    // constructor
}

// _________________________________________________________
StPicoPiPiMaker::~StPicoPiPiMaker() {
    // destructor
}

// _________________________________________________________
int StPicoPiPiMaker::InitHF() {
    // EXAMPLE //  mOutList->Add(new TH1F(...));
    // EXAMPLE //  TH1F* hist = static_cast<TH1F*>(mOutList->Last());

    mOutFileBaseName = mOutFileBaseName.ReplaceAll(".root", "");

    TString ntpVars = "pi1_pt:pi1_dca:pi1_nSigma:pi1_nHitFit:pi1_eta:pi1_phi:pi1_TOFinvbeta:pi2_pt:pi2_dca:pi2_nSigma:pi2_nHitFit:pi2_eta:pi2_phi:pi2_TOFinvbeta:dcaDaughters:hotSpot:primVz:primVzVpd:bbcRate:nTofTracks:nBTOFMatch:nHftTracks:pair_cosTheta:pair_decayL:pair_dcaToPv:pair_pt:pair_mass";

    ntp_signal = new TNtuple("ntp_signal","Ks_Signal", ntpVars);
    ntp_background = new TNtuple("ntp_background","Ks_background",ntpVars);

    return kStOK;
}

// _________________________________________________________
void StPicoPiPiMaker::ClearHF(Option_t *opt="") {
    return;
}

// _________________________________________________________
int StPicoPiPiMaker::FinishHF() {
    ntp_signal -> Write(ntp_signal->GetName(), TObject::kOverwrite);
    ntp_background -> Write(ntp_background->GetName(), TObject::kOverwrite);
    return kStOK;
}
// _________________________________________________________
int StPicoPiPiMaker::MakeHF() {
    createCandidates();
    return kStOK;
}

// _________________________________________________________
int StPicoPiPiMaker::createCandidates() {
    //making array of good pions
    float nTofTracks = 0;
    float nHftTracks = 0;
    for(unsigned int k = 0; k < mPicoDst->numberOfTracks(); ++k) {
        StPicoTrack const *trkTest = mPicoDst->track(k);
        if (mHFCuts->isTOFmatched(trkTest)) nTofTracks += 1;
        if (trkTest->isHFTTrack()) nHftTracks += 1;
        if (abs(trkTest->gMom().PseudoRapidity())>1) continue;
        if (mHFCuts->isGoodPion(trkTest)) mIdxPicoPions.push_back(k);
    }

    for (unsigned short j = 0; j < mIdxPicoPions.size(); ++j) {
        StPicoTrack const *pion1 = mPicoDst->track(mIdxPicoPions[j]);

        for(unsigned int i = j+1; i < mIdxPicoPions.size(); ++i)  {
            StPicoTrack const* pion2 = mPicoDst->track(mIdxPicoPions[i]);

//            if (pion1->id() == pion2->id()) continue;

            StHFPair *pair = new StHFPair(pion1, pion2, mHFCuts->getHypotheticalMass(StPicoCutsBase::kPion),mHFCuts->getHypotheticalMass(StPicoCutsBase::kPion), mIdxPicoPions[j], mIdxPicoPions[i], mPrimVtx, mBField, kTRUE);
            if (!mHFCuts->isGoodSecondaryVertexPair(pair)) continue;

            Float_t hotSpot=0;
            if (mHFCuts->checkHotSpot(&mPrimVtx)) hotSpot=1;

            bool isKs = false;
            if (pion1->charge()+pion2->charge() == 0) isKs = true;

            int ii=0;
            float ntVar[26];

            ntVar[ii++] = pion1->gPt();
            ntVar[ii++] = pair->particle1Dca();
            ntVar[ii++] = pion1->nSigmaPion();
            ntVar[ii++] = pion1->nHitsFit();
            ntVar[ii++] = pion1->gMom().PseudoRapidity();
            ntVar[ii++] = pion1->gMom().Phi();
            ntVar[ii++] = mHFCuts->getOneOverBeta(pion1, mHFCuts->getTofBetaBase(pion1), StPicoCutsBase::kPion);

            ntVar[ii++] = pion2->gPt();
            ntVar[ii++] = pair->particle2Dca();
            ntVar[ii++] = pion2->nSigmaPion();
            ntVar[ii++] = pion2->nHitsFit();
            ntVar[ii++] = pion2->gMom().PseudoRapidity();
            ntVar[ii++] = pion2->gMom().Phi();
            ntVar[ii++] = mHFCuts->getOneOverBeta(pion2, mHFCuts->getTofBetaBase(pion2), StPicoCutsBase::kPion);

            ntVar[ii++] = pair->dcaDaughters();
            ntVar[ii++] = hotSpot;
            ntVar[ii++] = mPrimVtx.z();
            ntVar[ii++] = mPicoEvent->vzVpd();
            ntVar[ii++] = mPicoEvent->BBCx() / 1000.;
            ntVar[ii++] = nTofTracks;
            ntVar[ii++] = (float)mPicoEvent->nBTOFMatch();
            ntVar[ii++] = nHftTracks;

            ntVar[ii++] = cos(pair->pointingAngle());
            ntVar[ii++] = pair->decayLength();

            ntVar[ii++] = pair->DcaToPrimaryVertex();
            ntVar[ii++] = pair->pt();
            ntVar[ii++] = pair->m();

            if (isKs) {
                ntp_signal->Fill(ntVar);
            } else {
                ntp_background->Fill(ntVar);
            }
        }  // for (unsigned short idxKaon = 0; idxKaon < mIdxPicoKaons.size(); ++idxKaon)
    } // for (unsigned short idxPion1 = 0; idxPion1 < mIdxPicoPions.size(); ++idxPion1)
    mIdxPicoPions.clear();
    return kStOK;
}