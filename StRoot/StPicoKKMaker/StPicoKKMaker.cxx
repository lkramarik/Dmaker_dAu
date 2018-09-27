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
    cout<<"FinishHF beg"<<endl;
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
    //making array of good pions
    for(unsigned int k = 0; k < mPicoDst->numberOfTracks(); ++k) {
        StPicoTrack const *trkTest = mPicoDst->track(k);
        if (!mHFCuts->isGoodKaon(trkTest)) mIdxPicoKaons.push_back(k);
    }

    for (unsigned short j = 0; j < mIdxPicoKaons.size(); ++j) {
        StPicoTrack const *kaon1 = mPicoDst->track(mIdxPicoKaons[j]);

        for(unsigned int i = j+1; i < mIdxPicoKaons.size(); ++i)  {
            StPicoTrack const* kaon2 = mPicoDst->track(mIdxPicoKaons[i]);

//            if (kaon1->id() == kaon2->id()) continue;

            StHFPair *pair = new StHFPair(kaon1, kaon2, mHFCuts->getHypotheticalMass(StPicoCutsBase::kKaon),mHFCuts->getHypotheticalMass(StPicoCutsBase::kKaon), i, j, mPrimVtx, mBField, kTRUE);
            if (!mHFCuts->isClosePair(pair)) continue;

            bool isPhi = false;
            if (kaon1->charge()+kaon2->charge() == 0) isPhi = true;

            int ii=0;
            float ntVar[19];
            ntVar[ii++] = kaon1->gPt();
            ntVar[ii++] = pair->particle1Dca();
            ntVar[ii++] = kaon1->nSigmaPion();
            ntVar[ii++] = kaon1->nHitsFit();
            ntVar[ii++] = mHFCuts->getOneOverBeta(kaon1, mHFCuts->getTofBetaBase(kaon1), StPicoCutsBase::kKaon);

            ntVar[ii++] = kaon2->gPt();
            ntVar[ii++] = pair->particle2Dca();
            ntVar[ii++] = kaon2->nSigmaPion();
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

    return kStOK;
}

bool StPicoKKMaker::isHadron(StPicoTrack const * const trk, int pidFlag) const {
    return (mHFCuts->isGoodTrack(trk) && mHFCuts->isTPCHadron(trk, pidFlag));
}

bool StPicoKKMaker::isPion(StPicoTrack const * const trk) const {
    if (!mHFCuts->cutMinDcaToPrimVertex(trk, StPicoCutsBase::kPion)) return false;
    if (!mHFCuts->isGoodTrack(trk)) return false; //HFT, NhitsFit, pt range
    if (!mHFCuts->isTPCHadron(trk, StPicoCutsBase::kPion)) return false;

    if (!mHFCuts->isTOFmatched(trk)) return false;

//    hybrid TOF
//    bool goodPion = (tofAvailable && tof && tpc) || (!tofAvailable && tpc);

    return true;

    //    float kBeta = mHFCuts->getTofBetaBase(trk);
    //    bool tofAvailable = (kBeta > 0) && (kBeta == kBeta);
//    if (mHFCuts->isTOFHadronPID(trk, kBeta, StPicoCutsBase::kPion)) tof = true; // 1/beta diff

//    bool goodPion = tpc;

    //    if (!mHFCuts->isTOFHadronPID(trk, mHFCuts->getTofBetaBase(trk), StPicoCutsBase::kPion) ) return false; // 1/beta diff
//    if (!mHFCuts->isHybridTOFHadron(trk, mHFCuts->getTofBetaBase(trk), StPicoCutsBase::kPion) ) return false;
}

bool StPicoKKMaker::isKaon(StPicoTrack const * const trk) const {
    if (!mHFCuts->cutMinDcaToPrimVertex(trk, StPicoCutsBase::kKaon)) return false;
    if (!mHFCuts->isGoodTrack(trk)) return false;
    if (!mHFCuts->isTPCHadron(trk, StPicoCutsBase::kKaon)) return false;

    float kBeta = mHFCuts->getTofBetaBase(trk);

    if (!mHFCuts->isTOFmatched(trk)) return false;
    if (!mHFCuts->isTOFHadronPID(trk, kBeta, StPicoCutsBase::kKaon)) return false; // 1/beta diff

//    hybrid TOF
//    bool goodKaon = (tofAvailable && tof && tpc) || (!tofAvailable && tpc);

//    LIANG:
//     bool goodKaon = (tofAvailable && tof) || (!tofAvailable && tpc);

    return true;
}

bool StPicoKKMaker::isProton(StPicoTrack const * const trk) const {
    return (mHFCuts->isGoodTrack(trk) && mHFCuts->isTPCHadron(trk, StPicoCutsBase::kProton));
}


