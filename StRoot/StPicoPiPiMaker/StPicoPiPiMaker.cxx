#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoDstMaker/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoHFMaker/StHFCuts.h"
#include "phys_constants.h"
#include "StBTofUtil/tofPathLength.hh"
#include "StPicoPiPiMaker.h"
ClassImp(StPicoPiPiMaker)

// _________________________________________________________
StPicoPiPiMaker::StPicoPiPiMaker(char const* name, StPicoDstMaker* picoMaker, char const* outputBaseFileName, char const* inputHFListHFtree = "") :
        StPicoHFMaker(name, picoMaker, outputBaseFileName, inputHFListHFtree),
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

    TString ntpVars = "pi1_pt:pi1_dca:pi1_nSigma:pi1_nHitFit:pi1_TOFinvbeta:pi2_pt:pi2_dca:pi2_nSigma:pi2_nHitFit:pi2_TOFinvbeta:dcaDaughters:primVz:primVzVpd:pair_cosTheta:pair_decayL:pair_dcaToPv:pair_cosThetaStar:pair_pt:pair_mass";

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
    for(unsigned int k = 0; k < mPicoDst->numberOfTracks(); ++k) {
        StPicoTrack const *trkTest = mPicoDst->track(k);
        if (mHFCuts->isGoodPion(trkTest) && mHFCuts->isTOFPion(trkTest)) mIdxPicoPions.push_back(k);
//        if (mHFCuts->isGoodPion(trkTest)) mIdxPicoPions.push_back(k);
    }
    for (unsigned short j = 0; j < mIdxPicoPions.size(); ++j) {
        StPicoTrack const *pion1 = mPicoDst->track(mIdxPicoPions[j]);

        for(unsigned int i = j+1; i < mIdxPicoPions.size(); ++i)  {
            StPicoTrack const* pion2 = mPicoDst->track(mIdxPicoPions[i]);

//            if (pion1->id() == pion2->id()) continue;

            StHFPair *pair = new StHFPair(pion1, pion2, mHFCuts->getHypotheticalMass(StPicoCutsBase::kPion),mHFCuts->getHypotheticalMass(StPicoCutsBase::kPion), i, j, mPrimVtx, mBField, kTRUE);
            if (!mHFCuts->isClosePair(pair)) continue;

            bool isKs = false;
            if (pion1->charge()+pion2->charge() == 0) isKs = true;

            int ii=0;
            float ntVar[19];
            ntVar[ii++] = pion1->gPt();
            ntVar[ii++] = pair->particle1Dca();
            ntVar[ii++] = pion1->nSigmaPion();
            ntVar[ii++] = pion1->nHitsFit();
            ntVar[ii++] = mHFCuts->getOneOverBeta(pion1, mHFCuts->getTofBetaBase(pion1), StPicoCutsBase::kPion);

            ntVar[ii++] = pion2->gPt();
            ntVar[ii++] = pair->particle2Dca();
            ntVar[ii++] = pion2->nSigmaPion();
            ntVar[ii++] = pion2->nHitsFit();
            ntVar[ii++] = mHFCuts->getOneOverBeta(pion2, mHFCuts->getTofBetaBase(pion2), StPicoCutsBase::kPion);

            ntVar[ii++] = pair->dcaDaughters();
            ntVar[ii++] = mPrimVtx.z();
            ntVar[ii++] = mPicoEvent->vzVpd();
            ntVar[ii++] = cos(pair->pointingAngle());
            ntVar[ii++] = pair->decayLength();

            ntVar[ii++] = pair->DcaToPrimaryVertex();
            ntVar[ii++] = pair->cosThetaStar();
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