#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoHFMaker/StHFCuts.h"
#include "phys_constants.h"
#include "StPicoKKMaker.h"
ClassImp(StPicoKKMaker)

// _________________________________________________________
StPicoKKMaker::StPicoKKMaker(char const* name, StPicoDstMaker* picoMaker, char const* outputBaseFileName) :
        StPicoHFMaker(name, picoMaker, outputBaseFileName),
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

    TString ntpVars = "pi1_pt:pi1_dca:pi1_nSigma:pi1_nHitFit:pi1_eta:pi1_phi:pi1_TOFinvbeta:pi2_pt:pi2_dca:pi2_nSigma:pi2_nHitFit:pi2_eta:pi2_phi:pi2_TOFinvbeta:dcaDaughters:hotSpot:primVz:primVzVpd:bbcRate:nTofTracks:nBTOFMatch:nHftTracks:pair_cosTheta:pair_decayL:pair_dcaToPv:pair_pt:pair_mass";

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
    float nTofTracks = 0;
    float nHftTracks = 0;
    for(unsigned int k = 0; k < mPicoDst->numberOfTracks(); ++k) {
        StPicoTrack const *trkTest = mPicoDst->track(k);
        if (mHFCuts->isTOFmatched(trkTest)) nTofTracks += 1;
        if (trkTest->isHFTTrack()) nHftTracks += 1;
        if (abs(trkTest->gMom().PseudoRapidity())>1) continue;
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

            Float_t hotSpot=0;
            if (mHFCuts->checkHotSpot(&mPrimVtx)) hotSpot=1;

            bool isPhi = false;
            if (kaon1->charge()+kaon2->charge() == 0) isPhi = true;

            int ii=0;
            float ntVar[26];
            ntVar[ii++] = kaon1->gPt();
            ntVar[ii++] = pair->particle1Dca();
            ntVar[ii++] = kaon1->nSigmaKaon();
            ntVar[ii++] = kaon1->nHitsFit();
            ntVar[ii++] = kaon1->gMom().PseudoRapidity();
            ntVar[ii++] = kaon1->gMom().Phi();
            ntVar[ii++] = mHFCuts->getOneOverBeta(kaon1, mHFCuts->getTofBetaBase(kaon1), StPicoCutsBase::kKaon);

            ntVar[ii++] = kaon2->gPt();
            ntVar[ii++] = pair->particle2Dca();
            ntVar[ii++] = kaon2->nSigmaKaon();
            ntVar[ii++] = kaon2->nHitsFit();
            ntVar[ii++] = kaon2->gMom().PseudoRapidity();
            ntVar[ii++] = kaon2->gMom().Phi();
            ntVar[ii++] = mHFCuts->getOneOverBeta(kaon2, mHFCuts->getTofBetaBase(kaon2), StPicoCutsBase::kKaon);

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