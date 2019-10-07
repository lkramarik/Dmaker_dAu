#include <limits>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <string>

#include "StPicoCutsBase.h"

#include "StPicoEvent/StPicoPhysicalHelix.h"
#include "phys_constants.h"
#include "SystemOfUnits.h"

#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"

ClassImp(StPicoCutsBase)

// _________________________________________________________
StPicoCutsBase::StPicoCutsBase() : TNamed("PicoCutsBase", "PicoCutsBase"),
                                   mTOFCorr(new StV0TofCorrection), mPicoDst(NULL), mEventStatMax(6), mTOFResolution(0.013),
                                   mBadRunListFileName("picoList_bad.list"), mVzMax(6.), mVzVpdVzMax(3.),
                                   mNHitsFitMin(15), mRequireHFT(true), mNHitsFitnHitsMax(0.), mPrimaryDCAtoVtxMax(6.0), mPtMin(0.2), mHybridTof(false), mOnlyHotSpot(false) {

    for (Int_t idx = 0; idx < kPicoPIDMax; ++idx) {
        mPtRange[idx][0] = std::numeric_limits<float>::lowest();
        mPtRange[idx][1] = std::numeric_limits<float>::max();
        mDcaMin[idx] = std::numeric_limits<float>::lowest();
        mDcaMinTertiary[idx] = std::numeric_limits<float>::lowest();
        mPtotRangeTOF[idx][0] = std::numeric_limits<float>::lowest();
        mPtotRangeTOF[idx][1] = std::numeric_limits<float>::max();
        mPtotRangeHybridTOF[idx][0] = std::numeric_limits<float>::lowest();
        mPtotRangeHybridTOF[idx][1] = std::numeric_limits<float>::max();
        mTPCNSigmaMax[idx] = std::numeric_limits<float>::max();
        mTOFDeltaOneOverBetaMax[idx] = std::numeric_limits<float>::max();
    }

    mHypotheticalMass[kPion]      = M_PION_PLUS;
    mHypotheticalMass2[kPion]     = M_PION_PLUS*M_PION_PLUS;
    mHypotheticalMass[kKaon]      = M_KAON_PLUS;
    mHypotheticalMass2[kKaon]     = M_KAON_PLUS*M_KAON_PLUS;
    mHypotheticalMass[kProton]    = M_PROTON;
    mHypotheticalMass2[kProton]   = M_PROTON*M_PROTON;
    mHypotheticalMass[kElectron]  = M_ELECTRON;
    mHypotheticalMass2[kElectron] = M_ELECTRON*M_ELECTRON;
    mHypotheticalMass[kMuon]      = M_MUON_PLUS;
    mHypotheticalMass2[kMuon]     = M_MUON_PLUS*M_MUON_PLUS;
    mHypotheticalMass[kK0Short]   = M_KAON_0_SHORT;
    mHypotheticalMass2[kK0Short]  = M_KAON_0_SHORT*M_KAON_0_SHORT;
    mHypotheticalMass[kLambda]    = M_LAMBDA;
    mHypotheticalMass2[kLambda]   = M_LAMBDA*M_LAMBDA;
}

// _________________________________________________________
StPicoCutsBase::StPicoCutsBase(const Char_t *name) : TNamed(name, name),
                                                     mTOFCorr(new StV0TofCorrection), mPicoDst(NULL), mEventStatMax(6), mTOFResolution(0.013),
                                                     mBadRunListFileName("picoList_bad_MB.list"), mVzMax(6.), mVzVpdVzMax(3.),
                                                     mNHitsFitMin(15), mRequireHFT(true), mNHitsFitnHitsMax(0.), mPrimaryDCAtoVtxMax(6.0), mPtMin(0.2), mHybridTof(false), mOnlyHotSpot(false) {
    // -- constructor

    for (Int_t idx = 0; idx < kPicoPIDMax; ++idx) {
        mPtRange[idx][0] = std::numeric_limits<float>::lowest();
        mPtRange[idx][1] = std::numeric_limits<float>::max();
        mDcaMin[idx] = std::numeric_limits<float>::lowest();
        mDcaMinTertiary[idx] = std::numeric_limits<float>::lowest();
        mPtotRangeTOF[idx][0] = std::numeric_limits<float>::lowest();
        mPtotRangeTOF[idx][1] = std::numeric_limits<float>::max();
        mPtotRangeHybridTOF[idx][0] = std::numeric_limits<float>::lowest();
        mPtotRangeHybridTOF[idx][1] = std::numeric_limits<float>::max();
        mTPCNSigmaMax[idx] = std::numeric_limits<float>::max();
        mTOFDeltaOneOverBetaMax[idx] = std::numeric_limits<float>::max();
    }

    mHypotheticalMass[kPion]      = M_PION_PLUS;
    mHypotheticalMass2[kPion]     = M_PION_PLUS*M_PION_PLUS;
    mHypotheticalMass[kKaon]      = M_KAON_PLUS;
    mHypotheticalMass2[kKaon]     = M_KAON_PLUS*M_KAON_PLUS;
    mHypotheticalMass[kProton]    = M_PROTON;
    mHypotheticalMass2[kProton]   = M_PROTON*M_PROTON;
    mHypotheticalMass[kElectron]  = M_ELECTRON;
    mHypotheticalMass2[kElectron] = M_ELECTRON*M_ELECTRON;
    mHypotheticalMass[kMuon]      = M_MUON_PLUS;
    mHypotheticalMass2[kMuon]     = M_MUON_PLUS*M_MUON_PLUS;
    mHypotheticalMass[kK0Short]   = M_KAON_0_SHORT;
    mHypotheticalMass2[kK0Short]  = M_KAON_0_SHORT*M_KAON_0_SHORT;
    mHypotheticalMass[kLambda]    = M_LAMBDA;
    mHypotheticalMass2[kLambda]   = M_LAMBDA*M_LAMBDA;

}
// _________________________________________________________
StPicoCutsBase::~StPicoCutsBase() {
    // destructor

    if (mTOFCorr)
        delete mTOFCorr;
    mTOFCorr = NULL;
}

// _________________________________________________________
void StPicoCutsBase::initBase() {
    // -- init cuts class

    // -- Read in bad run list and fill vector
    // -----------------------------------------

    // -- open list
    ifstream runs;

    // -- open in working dir
    runs.open(mBadRunListFileName.Data());
    if (!runs.is_open()) {
        runs.open(Form("picoLists/%s", mBadRunListFileName.Data()));
        if (!runs.is_open()) {
            cout << "StPicoCutsBase::initBase -- Bad run list NOT found :" << mBadRunListFileName << endl;
            cout << "StPicoCutsBase::initBase -- continue without bad run selection! " << endl;
            //exit(EXIT_FAILURE);
        }
    }

    if (runs.is_open()) {
        Int_t runId = 0;
        while( runs >> runId )
            mVecBadRunList.push_back(runId);

        runs.close();

        // -- sort bad runs vector
        std::sort(mVecBadRunList.begin(), mVecBadRunList.end());
    }
}

// _________________________________________________________
bool StPicoCutsBase::isGoodEvent(StPicoDst const * const picoDst, int *aEventCuts) {
    // -- set current mPicoDst
    mPicoDst = picoDst;

    // -- get picoDst event
    StPicoEvent* picoEvent = mPicoDst->event();

    // -- set current primary vertex
    mPrimVtx = picoEvent->primaryVertex();

//    if(mOnlyHotSpot) cout<<"m hot spor true"<<endl;
    if(mOnlyHotSpot && !(checkHotSpot(&mPrimVtx))) return false;

    // -- quick method without providing stats
    if (!aEventCuts) {
        return (isGoodRun(picoEvent) && isGoodTrigger(picoEvent) &&
                fabs(picoEvent->primaryVertex().z()) < mVzMax &&
                fabs(picoEvent->primaryVertex().z() - picoEvent->vzVpd()) < mVzVpdVzMax);
    }

    // -- reset event cuts
    for (unsigned int ii = 0; ii < mEventStatMax; ++ii)
        aEventCuts[ii] = 0;

    unsigned int iCut = 0;
    // -- 0 - before event cuts
    aEventCuts[iCut] = 0;

    // -- 1 - is bad run
    ++iCut;
    if (!isGoodRun(picoEvent)) aEventCuts[iCut] = 1;

    // -- 2 - No Trigger fired
    ++iCut;
//  trigger - is ok?
    if (!isGoodTrigger(picoEvent)) aEventCuts[iCut] = 1;

    // -- 3 - Vertex z outside cut window
    ++iCut;
    if (fabs(picoEvent->primaryVertex().z()) >= mVzMax) aEventCuts[iCut] = 1;

    // -- 4 Vertex z - vertex_z(vpd) outside cut window
    ++iCut;
    if (fabs(picoEvent->primaryVertex().z() - picoEvent->vzVpd()) >= mVzVpdVzMax) aEventCuts[iCut] = 1;

    ++iCut;

    //if the event is wrong, the array member is 1 (eccept [0])
    // -- is rejected
    bool isGoodEvent = true;
    for (unsigned int ii = 0; ii < mEventStatMax-1; ++ii) {
        if (aEventCuts[ii]) isGoodEvent = false;
    }

    if(!isGoodEvent) aEventCuts[mEventStatMax-1]=1;

    return isGoodEvent;
}

// _________________________________________________________
bool StPicoCutsBase::isGoodRun(StPicoEvent const * const picoEvent) const {
    // -- is good run (not in bad runlist)
    return (!(std::binary_search(mVecBadRunList.begin(), mVecBadRunList.end(), picoEvent->runId())));
}

// _________________________________________________________
bool StPicoCutsBase::isGoodTrigger(StPicoEvent const * const picoEvent) const {
    // -- is good trigger in list of good triggerIds

    for(std::vector<unsigned int>::const_iterator iter = mVecTriggerIdList.begin(); iter != mVecTriggerIdList.end(); ++iter)
        if(picoEvent->isTrigger(*iter))
            return true;

    return false;
}

// _________________________________________________________
bool StPicoCutsBase::isGoodTrack(StPicoTrack const * const trk) const {
    //  int tofIndex = trk->bTofPidTraitsIndex();
//  bool TofMatch = kFALSE;
//  StPicoBTofPidTraits* tofPidTraits;
//  if (tofIndex >= 0)  tofPidTraits = mPicoDst->btofPidTraits(tofIndex);
//  if (tofIndex >= 0 && tofPidTraits && tofPidTraits->btofMatchFlag() > 0)  TofMatch = kTRUE;
    return ((!mRequireHFT || trk->isHFTTrack()) && trk->nHitsFit() >= mNHitsFitMin && cutMaxDcaToPrimVertex(trk) && trk->gPt() > mPtMin);
}

// _________________________________________________________
bool StPicoCutsBase::isGoodPion(StPicoTrack const *const trk) const {
    if (!isGoodTrack(trk)) return false;
    if (!cutMinDcaToPrimVertex(trk, StPicoCutsBase::kPion)) return false;
    if (!isTPCPion(trk)) return false;

    bool tof = false;
    if (mHybridTof) tof = isHybridTOFPion(trk);
    if (!mHybridTof) tof = isTOFmatched(trk);
    return tof;
}

// _________________________________________________________
bool StPicoCutsBase::isGoodKaon(StPicoTrack const *const trk) const {
    if (!isGoodTrack(trk)) return false;
    if (!cutMinDcaToPrimVertex(trk, StPicoCutsBase::kKaon)) return false;
    if (!isTPCKaon(trk)) return false;
    bool tof = false;
    if (mHybridTof) tof = isHybridTOFKaon(trk);
    if (!mHybridTof) tof = isTOFKaon(trk);

    return tof;
}

// _________________________________________________________
bool StPicoCutsBase::isGoodProton(StPicoTrack const *const trk) const {
    if (!isGoodTrack(trk)) return false;
    if (!cutMinDcaToPrimVertex(trk, StPicoCutsBase::kProton)) return false;
    if (!isTPCProton(trk)) return false;
    bool tof = false;
    if (mHybridTof) tof = isHybridTOFProton(trk);
    if (!mHybridTof) tof = isTOFProton(trk);

    return tof;
}

// _________________________________________________________
bool StPicoCutsBase::cutMinDcaToPrimVertex(StPicoTrack const * const trk, int pidFlag) const {
    // -- check on min dca for identified particle
    float dca = (mPrimVtx - trk->origin()).Mag();
    return (dca >= mDcaMin[pidFlag]);
}

// _________________________________________________________
bool StPicoCutsBase::cutMaxDcaToPrimVertex(StPicoTrack const * const trk) const {
    // -- check on max dca for all particles
    float dca = (mPrimVtx - trk->origin()).Mag();
    return (dca <= mPrimaryDCAtoVtxMax);
}

// _________________________________________________________
bool StPicoCutsBase::cutMinDcaToPrimVertexTertiary(StPicoTrack const * const trk, int pidFlag) const {
    // -- check on min dca for identified particle - used for tertiary particles only

    StPicoPhysicalHelix helix = trk->helix(mPicoDst->event()->bField());
    helix.moveOrigin(helix.pathLength(mPrimVtx));
    float dca = (mPrimVtx - helix.origin()).Mag();

    return (dca >= mDcaMinTertiary[pidFlag]);
}

// _________________________________________________________
bool StPicoCutsBase::isTPCHadron(StPicoTrack const * const trk, int pidFlag) const {
    // -- check for good hadron in TPC
    float nSigma = std::numeric_limits<float>::quiet_NaN();

    if (pidFlag == kPion)
        nSigma = fabs(trk->nSigmaPion());
    else if (pidFlag == kKaon)
        nSigma = fabs(trk->nSigmaKaon());
    else if (pidFlag == kProton)
        nSigma = fabs(trk->nSigmaProton());

    return (nSigma < mTPCNSigmaMax[pidFlag] && trk->nHitsFit() >= mNHitsFitMin);
}

// _________________________________________________________
bool StPicoCutsBase::isTOFHadronPID(StPicoTrack const *trk, float const & tofBeta, int pidFlag) const {
    if (tofBeta <= 0) {return false;}
    double ptot    = trk->gPtot();
    float betaInv = ptot / sqrt(ptot*ptot + mHypotheticalMass2[pidFlag]);
    return ( fabs(1/tofBeta - 1/betaInv) < mTOFDeltaOneOverBetaMax[pidFlag] );
}

// _________________________________________________________
bool StPicoCutsBase::isTOFHadron(StPicoTrack const *trk, float const & tofBeta, int pidFlag) const {
    // -- check for good hadron in TOF in ptot range
    //    use for
    //      - primary hadrons
    //      - secondarys from charm decays (as an approximation)
    //    return:
    //      not in ptot range : true

    // -- only apply, if in ptot range


//  float ptot = trk->gPtot();
//  if (ptot < mPtotRangeTOF[pidFlag][0] || ptot >= mPtotRangeTOF[pidFlag][1])
//    return true;

    return isTOFHadronPID(trk, tofBeta, pidFlag);
}

// _________________________________________________________
bool StPicoCutsBase::isTOFmatched(StPicoTrack const *trk) const {
    /*OLD
    int tofIndex = trk->bTofPidTraitsIndex();
    trk->isTofTrack();
    bool TofMatch = kFALSE;
    StPicoBTofPidTraits* tofPidTraits;
    if (tofIndex >= 0)  tofPidTraits = mPicoDst->btofPidTraits(tofIndex);
    if (tofIndex >= 0 && tofPidTraits && tofPidTraits->btofMatchFlag() > 0)  TofMatch = kTRUE;
     */
    return trk->isTofTrack();
}

// _________________________________________________________
bool StPicoCutsBase::isHybridTOFHadron(StPicoTrack const *trk, float const & tofBeta, int pidFlag) const {
    // -- check for good hadron in TOF in ptot range
    //    use for
    //      - primary hadrons
    //      - secondarys from charm decays (as an approximation)
    //    return:
    //      not in ptot range : true
    //      no TOF info       : true

    // -- only apply, if in ptot range
//  float ptot = trk->gPtot();
//  if (ptot < mPtotRangeHybridTOF[pidFlag][0] || ptot >= mPtotRangeHybridTOF[pidFlag][1])
//    return true;

    // -- only apply, if has TOF information
    if (tofBeta <= 0 || tofBeta != tofBeta )
        return true;

    return isTOFHadronPID(trk, tofBeta, pidFlag);
}

// _________________________________________________________
StPicoBTofPidTraits* StPicoCutsBase::hasTofPid(StPicoTrack const * const trk) const {
    // -- check if track has TOF pid information
    //    return NULL otherwise

    int index2tof = trk->bTofPidTraitsIndex();
    return (index2tof >= 0) ? mPicoDst->btofPidTraits(index2tof) : NULL;
}

// _________________________________________________________
float StPicoCutsBase::getTofBetaBase(StPicoTrack const * const trk) const {
    int index2tof = trk->bTofPidTraitsIndex(); //if smaller than 0 => not TOF track
    float beta = std::numeric_limits<float>::quiet_NaN();

    if(index2tof >= 0) {
        StPicoBTofPidTraits *tofPid = mPicoDst->btofPidTraits(index2tof);
        if(tofPid)  beta = tofPid->btofBeta();

        if (beta < 1e-4) {
            TVector3 const btofHitPos = tofPid->btofHitPos();
            StPicoPhysicalHelix helix = trk->helix(mPicoDst->event()->bField());

            float L = tofPathLength(&mPrimVtx, &btofHitPos, helix.curvature());
            float tof = tofPid->btof();
            if (tof > 0) beta = L / (tof * (C_C_LIGHT / 1.e9));
            else beta = std::numeric_limits<float>::quiet_NaN();
        }
    }

    return beta;
}

// _________________________________________________________
float StPicoCutsBase::getOneOverBeta(StPicoTrack const * const trk,  float const & tofBeta, int pidFlag) const {
    if ((tofBeta <= 0) || (tofBeta!=tofBeta))
        return std::numeric_limits<float>::max();
    float m2 = mHypotheticalMass[pidFlag]*mHypotheticalMass[pidFlag];
    float ptot = trk->gPtot();
    float betaInv = ptot / sqrt(ptot*ptot + m2);
    return (1/tofBeta - 1/betaInv);
}

// _________________________________________________________
float StPicoCutsBase::getTofBeta(StPicoTrack const * const trk) const {
    return getTofBetaBase(trk);
}

// _________________________________________________________
float StPicoCutsBase::getTofBeta(StPicoTrack const * const trk,
                                 TVector3 const & secondaryMother, TVector3 const & secondaryVtx) const {
    // -- provide correced beta of TOF for pico track
    //    use for
    //      - secondaries

    float beta = std::numeric_limits<float>::quiet_NaN();

    StPicoBTofPidTraits *tofPid = hasTofPid(trk);
    if (!tofPid)
        return beta;
//
//
//  // -- set waypoints
//  mTOFCorr->setVectors3D(mPrimVtx)(secondaryVtx)(tofHit);
//
//  // -- set mother track
//  mTOFCorr->setMotherTracks(secondaryMother);
//
//  float tof = tofPid->btof();
//  StPicoPhysicalHelix helix = trk->helix(mPicoDst->event()->bField());
//
//  // -- correct beta
//  mTOFCorr->correctBeta(helix, tof, beta);
//
//  // -- clean up
//  mTOFCorr->clearContainers();
//
    return beta;
}

// _________________________________________________________
float StPicoCutsBase::getTofBeta(StPicoTrack const * const trk,
                                 TVector3 const & secondaryMother, TVector3 const & secondaryVtx,
                                 TVector3 const & tertiaryMother,  TVector3 const & tertiaryVtx) const {
    // -- provide correced beta of TOF for pico track
    //    use for
    //      - tertiaries

    float beta = std::numeric_limits<float>::quiet_NaN();

//  StPicoBTofPidTraits *tofPid = hasTofPid(trk);
//  if (!tofPid)
//    return beta;
//
//  StThreeVectorD tofHit = tofPid->btofHitPos();
//
//  // -- set waypoints
//  mTOFCorr->setVectors3D(mPrimVtx)(secondaryVtx)(tertiaryVtx)(tofHit);
//
//  // -- set mother track
//  mTOFCorr->setMotherTracks(secondaryMother)(tertiaryMother);
//
//  float tof = tofPid->btof();
//  StPicoPhysicalHelix helix = trk->helix(mPicoDst->event()->bField());
//
//  // -- correct beta
////  mTOFCorr->correctBeta(helix, tof, beta);
//
//  // -- clean up
//  mTOFCorr->clearContainers();

    return beta;
}

// _________________________________________________________
float StPicoCutsBase::tofPathLength(const TVector3* beginPoint, const TVector3* endPoint, float curvature) const {
    float xdif =  endPoint->x() - beginPoint->x();
    float ydif =  endPoint->y() - beginPoint->y();

    float C = sqrt(xdif*xdif + ydif*ydif);
    float s_perp = C;
    if (curvature){
        float R = 1/curvature;
        s_perp = 2*R * asin(C/(2*R));
    }

    float s_z = fabs(endPoint->z() - beginPoint->z());
    float value = sqrt(s_perp*s_perp + s_z*s_z);

    return value;
}

// _________________________________________________________
bool StPicoCutsBase::checkHotSpot(TVector3* vertex) const {
    if(vertex->x()>-0.25 && vertex->x()<-0.16 && vertex->y()>-0.25 && vertex->y()<-0.16) return true;
    else return false;
}
