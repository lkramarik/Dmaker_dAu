#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoDstMaker/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoHFMaker/StHFCuts.h"
#include "phys_constants.h"
#include "StBTofUtil/tofPathLength.hh"
#include "StPicoKKMaker.h"
ClassImp(StPicoKKMaker)

// _________________________________________________________
StPicoKKMaker::StPicoKKMaker(char const* name, StPicoDstMaker* picoMaker, char const* outputBaseFileName, char const* inputHFListHFtree = "") :
        StPicoHFMaker(name, picoMaker, outputBaseFileName, inputHFListHFtree),
        mOutFileBaseName(outputBaseFileName){
    // constructor
}

// _________________________________________________________
StPicoKKMaker::~StPicoKKMaker() {
    // destructor
}

// _________________________________________________________
int StPicoKKMaker::InitHF() {
    // EXAMPLE //  mOutList->Add(new TH1F(...));
    // EXAMPLE //  TH1F* hist = static_cast<TH1F*>(mOutList->Last());

    mOutFileBaseName = mOutFileBaseName.ReplaceAll(".root", "");

    TString ntpVars = "pi1_pt:pi1_dca:pi1_nSigma:pi1_nHitFit:pi1_TOFinvbeta:pi2_pt:pi2_dca:pi2_nSigma:pi2_nHitFit:pi2_TOFinvbeta:dcaDaughters:primVz:primVzVpd:pair_cosTheta:pair_decayL:pair_dcaToPv:pair_cosThetaStar:pair_pt:pair_mass";

    ntp_signal = new TNtuple("ntp_signal","phi_Signal", ntpVars);
    ntp_background = new TNtuple("ntp_background","phi_background",ntpVars);

    return kStOK;
}

// _________________________________________________________
void StPicoKKMaker::ClearHF(Option_t *opt="") {
    return;
}

// _________________________________________________________
int StPicoKKMaker::FinishHF() {
    ntp_signal -> Write(ntp_signal->GetName(), TObject::kOverwrite);
    ntp_background -> Write(ntp_background->GetName(), TObject::kOverwrite);
    return kStOK;
}
// _________________________________________________________
int StPicoKKMaker::MakeHF() {
    createCandidates();
    return kStOK;
}

// _________________________________________________________
int StPicoKKMaker::createCandidates() {
    //making array of good kaons
    for(unsigned int k = 0; k < mPicoDst->numberOfTracks(); ++k) {
        StPicoTrack const *trkTest = mPicoDst->track(k);
        if (mHFCuts->isGoodKaon(trkTest)) mIdxPicoKaons.push_back(k);
    }
    for (unsigned short j = 0; j < mIdxPicoKaons.size(); ++j) {
        StPicoTrack const *kaon1 = mPicoDst->track(mIdxPicoKaons[j]);

        for(unsigned int i = j+1; i < mIdxPicoKaons.size(); ++i)  {
            StPicoTrack const* kaon2 = mPicoDst->track(mIdxPicoKaons[i]);

//            if (kaon1->id() == kaon2->id()) continue;

            StHFPair *pair = new StHFPair(kaon1, kaon2, mHFCuts->getHypotheticalMass(StPicoCutsBase::kKaon),mHFCuts->getHypotheticalMass(StPicoCutsBase::kKaon), mIdxPicoKaons[j], mIdxPicoKaons[i], mPrimVtx, mBField, kTRUE);
//            if (!mHFCuts->isClosePair(pair)) continue;
            if (!mHFCuts->isGoodSecondaryVertexPair(pair)) continue;

            bool isPhi = false;
            if (kaon1->charge()+kaon2->charge() == 0) isPhi = true;

            int ii=0;
            float ntVar[19];
            ntVar[ii++] = kaon1->gPt();
            ntVar[ii++] = pair->particle1Dca();
            ntVar[ii++] = kaon1->nSigmaKaon();
            ntVar[ii++] = kaon1->nHitsFit();
            ntVar[ii++] = mHFCuts->getOneOverBeta(kaon1, mHFCuts->getTofBetaBase(kaon1), StPicoCutsBase::kKaon);

            ntVar[ii++] = kaon2->gPt();
            ntVar[ii++] = pair->particle2Dca();
            ntVar[ii++] = kaon2->nSigmaKaon();
            ntVar[ii++] = kaon2->nHitsFit();
            ntVar[ii++] = mHFCuts->getOneOverBeta(kaon2, mHFCuts->getTofBetaBase(kaon2), StPicoCutsBase::kKaon);

            ntVar[ii++] = pair->dcaDaughters();
            ntVar[ii++] = mPrimVtx.z();
            ntVar[ii++] = mPicoEvent->vzVpd();
            ntVar[ii++] = cos(pair->pointingAngle());
            ntVar[ii++] = pair->decayLength();

            ntVar[ii++] = pair->DcaToPrimaryVertex();
            ntVar[ii++] = pair->cosThetaStar();
            ntVar[ii++] = pair->pt();
            ntVar[ii++] = pair->m();

            if (isPhi) {
                ntp_signal->Fill(ntVar);
            } else {
                ntp_background->Fill(ntVar);
            }
        }  // for (unsigned short idxKaon = 0; idxKaon < mIdxPicoKaons.size(); ++idxKaon)
    } // for (unsigned short idxPion1 = 0; idxPion1 < mIdxPicoPions.size(); ++idxPion1)
    mIdxPicoKaons.clear();
    return kStOK;
}