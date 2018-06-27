#include <limits>

#include "TNtuple.h"
#include "TTree.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TList.h"

#include "StPicoEventMixer.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoDstMaker/StPicoDst.h"

#include "StPicoMixedEventMaker.h"
#include "StMixerEvent.h"
#include "StMixerHists.h"
#include "StPicoHFMaker/StHFCuts.h"
#include "StPicoHFMaker/StHFPair.h"

//-----------------------------------------------------------
StPicoEventMixer::StPicoEventMixer(char* category):
        mEvents(),
        mHists(NULL),
        mHFCuts(NULL),
        mEventsBuffer(StPicoMixedEventMaker::defaultBufferSize),
        filledBuffer(0),
        mSETupleSig(NULL),
        mSETupleBack(NULL),
        mMETupleSig(NULL),
        mMETupleBack(NULL),
        mSingleParticleList(NULL),
        fillSinglePartHists(false)
{
    mHists = new StMixerHists(category);
}
//-----------------------------------------------------------
StPicoEventMixer::~StPicoEventMixer()
{
    delete mHists;
    for(unsigned int i =0 ; i<mEvents.size() ; i++){
        delete mEvents.at(i);
    }
}
//-----------------------------------------------------------
void StPicoEventMixer::finish() {
    mHists->closeFile();
}
//-----------------------------------------------------------
bool StPicoEventMixer::addPicoEvent(StPicoDst const* const picoDst, float weight)
{
    if( !isGoodEvent(picoDst) )
        return false;
    int nTracks = picoDst->numberOfTracks();
    StThreeVectorF pVertex = picoDst->event()->primaryVertex();
    StMixerEvent* event = new StMixerEvent(pVertex, picoDst->event()->bField());

    event->addPicoEvent(*(picoDst->event()));
    event->setWeight(weight);

    bool isTpcPi = false;
    bool isTofPi = false;
    bool isTpcP  = false;
    bool isTofP  = false;
    bool isTpcK  = false;
    bool isTofK  = false;
    //Event.setNoTracks( nTracks );
    for( int iTrk = 0; iTrk < nTracks; ++iTrk) {
        StPicoTrack const* trk = picoDst->track(iTrk);

        if(!mHFCuts->isGoodTrack(trk))
            continue;

        const float beta = mHFCuts->getTofBeta(trk);

        bool saveTrack = false;
        int pidFlag = StPicoCutsBase::kPion;
        if( isTpcPion(trk) && mHFCuts->isHybridTOFPion(trk, beta, pVertex) && mHFCuts->cutMinDcaToPrimVertex(trk, pidFlag)) {
            isTpcPi = true;
            isTofPi = true;
            saveTrack = true;
            event->addPion(event->getNoTracks());
        }
        pidFlag = StPicoCutsBase::kKaon;
        if(isTpcKaon(trk) && mHFCuts->isTOFKaon(trk, beta, pVertex) && mHFCuts->cutMinDcaToPrimVertex(trk, pidFlag)) {
            isTpcK = true;
            isTofK = true;
            saveTrack = true;
            event->addKaon(event->getNoTracks());
        }
        pidFlag = StPicoCutsBase::kProton;
        if(isTpcProton(trk) && mHFCuts->isTOFProton(trk, beta, pVertex) && mHFCuts->cutMinDcaToPrimVertex(trk, pidFlag)) {
            isTpcP = true;
            isTofP = true;
            saveTrack = true;
            event->addProton(event->getNoTracks());
        }
        if(saveTrack == true){
            StMixerTrack mTrack(pVertex, picoDst->event()->bField(), *trk, isTpcPi, isTofPi, isTpcK, isTofK, isTpcP, isTofP);
            event->addTrack(mTrack);
        }
    }
    // if ( event->getNoPions() > 0 ||  event->getNoKaons() > 0 || event->getNoProtons() > 0) {
    mEvents.push_back(event);
    filledBuffer+=1;
    // }
    // else {
    //   delete event;
    //   return false;
    // }
    //Returns true if need to do mixing, false if buffer has space still
    if ( filledBuffer == mEventsBuffer)
        return true;
    return false;
}
//-----------------------------------------------------------
void StPicoEventMixer::mixEvents() {
    size_t const nEvent = mEvents.size();
    if(StPicoMixedEventMaker::fillSingleTrackHistos)
    {
    }

    int const nTracksEvt1 = mEvents.at(0)->getNoPions();

    // Check if there are kaons in the first evt for saving time (cannot be done if
    // we want to save the single particle ctrl plots)
    if(!StPicoMixedEventMaker::fillSingleTrackHistos && nTracksEvt1 == 0)
    {
        --filledBuffer;
        return;
    }

    // Go through the event buffer
    for( size_t iEvt2 = 0; iEvt2 < nEvent; ++iEvt2) {
        int const nTracksEvt2 = mEvents.at(iEvt2)->getNoKaons();

        if( iEvt2 == 0 )
            mHists->fillSameEvt(mEvents.at(0)->vertex());
        else
            mHists->fillMixedEvt(mEvents.at(0)->vertex());

        // evts trk loops
        for(int iTrk1 = 0; iTrk1 < nTracksEvt1; ++iTrk1) {
            for( int iTrk2 = 0; iTrk2 < nTracksEvt2; ++iTrk2) {
                // check if the tracks are the same
                if(iEvt2 == 0) {
                    if(mEvents.at(0)->pionId(iTrk1) == mEvents.at(iEvt2)->kaonId(iTrk2)) //?
                        continue;
                }

                StPicoTrack const* pion = mEvents.at(0)->pionAt(iTrk1);
//                StMixerTrack const* kaon = mEvents.at(iEvt2)->kaonAt(iTrk2);
                StPicoTrack const* kaon = mEvents.at(iEvt2)->kaonAt(iTrk2);

                StHFPair *pair = new StHFPair(pion, kaon,
                                 mHFCuts->getHypotheticalMass(StHFCuts::kPion),
                                 mHFCuts->getHypotheticalMass(StHFCuts::kKaon),
                                 mEvents.at(0)->vertex(), mEvents.at(iEvt2)->vertex(),
                                 mEvents.at(0)->field() );

                if (!mHFCuts->isClosePair(pair)) continue;

                int ii=0;
                float ntVar[21];
                ntVar[ii++] = pion->gPt();
                ntVar[ii++] = pair->particle1Dca();
                ntVar[ii++] = pion->nSigmaPion();
                ntVar[ii++] = pion->nHitsFit();
                ntVar[ii++] = mHFCuts->getOneOverBeta(pion, mHFCuts->getTofBetaBase(pion), StPicoCutsBase::kPion);
                ntVar[ii++] = kaon->gPt();
                ntVar[ii++] = pair->particle2Dca();
                ntVar[ii++] = kaon->nSigmaKaon();
                ntVar[ii++] = kaon->nHitsFit();
                ntVar[ii++] = mHFCuts->getOneOverBeta(kaon, mHFCuts->getTofBetaBase(kaon), StPicoCutsBase::kKaon);
                ntVar[ii++] = pair->dcaDaughters();
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

                int charge = mEvents.at(0)->pionAt(iTrk1).charge() +  mEvents.at(iEvt2)->kaonAt(iTrk2).charge(); // 0 = signal

                if(iEvt2 == 0)
                    fillNtpSameEvtPair(ntVar, &pair, charge );
                else
                    fillNtpMixedEvtPair(ntVar, &pair, charge);
            }
        } //second track track loop
    } // the first track track loop

--filledBuffer;
delete mEvents.at(0);
mEvents.erase(mEvents.begin());
}

void StPicoEventMixer::fillNtpSameEvtPair(float ntVar[21], StHFPair const & pair, int charge)
{
    if(charge == 0 )
        mSETupleSig -> Fill(ntVar);
    else
        mSETupleBack-> Fill(ntVar);
    return;
}


void StPicoEventMixer::fillNtpMixedEvtPair(float ntVar[21], StHFPair const & pair, int charge)
{
    if(charge == 0 )
        mMETupleSig -> Fill(ntVar);
    else
        mMETupleBack-> Fill(ntVar);
    return;
}

// _________________________________________________________
void StPicoEventMixer::fillTracks(StMixerEvent* evt, bool isSameEvt, int pidFlag)
{
    if (!fillSinglePartHists)
        return;

    // get the corresponting histograms and track vectors
    const std::string evtName = isSameEvt ? "SE" : "ME";
    std::string particleName;
    int nTracks;
    switch (pidFlag)
    {
        case StHFCuts::kProton:
            particleName = "p";
            nTracks = evt->getNoProtons();
            break;
        case StHFCuts::kKaon:
            particleName = "K";
            nTracks = evt->getNoKaons();
            break;
        case StHFCuts::kPion:
            particleName = "pi";
            nTracks = evt->getNoPions();
            break;
        default:
            cerr << "StPicoEventMixer::fillTracks: unknown pidFlag ... exiting" << endl;
            throw;
    }
    TH2D *etaPhiHist = static_cast<TH2D*>(mSingleParticleList->FindObject(Form("%sEtaPhi%s",particleName.data(), evtName.data())));
    TH2D *phiPtHist  = static_cast<TH2D*>(mSingleParticleList->FindObject(Form("%sPhiPt%s",particleName.data(), evtName.data())));
    TH1D *dcaHist = static_cast<TH1D*>(mSingleParticleList->FindObject(Form("%sDCA%s",particleName.data(), evtName.data())));
    TH1D *nTracksHist = static_cast<TH1D*>(mSingleParticleList->FindObject(Form("%stracks%s",particleName.data(), evtName.data())));
    nTracksHist->Fill(nTracks);

    // particle loop
    const float weight = evt->weight();
    for (int i = 0; i < nTracks; ++i)
    {
        StMixerTrack trk;
        switch (pidFlag)
        {
            case StHFCuts::kProton: trk = evt->protonAt(i);
                break;
            case StHFCuts::kKaon:   trk = evt->kaonAt(i);
                break;
            case StHFCuts::kPion:   trk = evt->pionAt(i);
                break;
        }
        const float eta = trk.gMom().pseudoRapidity();
        const float phi = trk.gMom().phi();
        const float pt = trk.gMom().perp();
        const float dca  = (trk.origin() - evt->vertex()).mag();

        etaPhiHist->Fill(phi,eta,weight);
        phiPtHist->Fill(pt,phi, weight);
        dcaHist->Fill(dca, weight);
    }

}
// _________________________________________________________
bool StPicoEventMixer::isMixerPion(StMixerTrack const& track) {
    short info = track.getTrackInfo();
    //TPC pion
    if( (info & (1 << kPionTPCbit)) >> kPionTPCbit != 1) return false;
    //TOF pion
    if( (info & (1 << kPionTOFbit)) >> kPionTOFbit != 1) return false;
    return true;
}
// _________________________________________________________
bool StPicoEventMixer::isMixerKaon(StMixerTrack const& track) {
    short info = track.getTrackInfo();
    //TPC Kaon
    if( (info & (1 << kKaonTPCbit)) >> kKaonTPCbit != 1) return false;
    //TOF Kaon
    if( (info & (1 << kKaonTOFbit)) >> kKaonTOFbit != 1) return false;
    return true;
}
//-----------------------------------------------------------
bool StPicoEventMixer::isMixerProton(StMixerTrack const& track) {
    short info = track.getTrackInfo();
    //TPC Proton
    if( (info & (1 << kProtonTPCbit)) >> kProtonTPCbit != 1) return false;
    //TOF Proton
    if( (info & (1 << kProtonTOFbit)) >> kProtonTOFbit != 1) return false;
    return true;
}
//-----------------------------------------------------------
bool StPicoEventMixer::isGoodEvent(StPicoDst const * const picoDst)
{
    return (mHFCuts->isGoodEvent(picoDst));
}
//-----------------------------------------------------------
bool StPicoEventMixer::isTpcPion(StPicoTrack const * const trk)
{
    return( isTPCHadron(trk, StPicoCutsBase::kPion));
}
//-----------------------------------------------------------
bool StPicoEventMixer::isTpcKaon(StPicoTrack const * const trk)
{
    return( isTPCHadron(trk, StPicoCutsBase::kKaon));
}
//-----------------------------------------------------------
bool StPicoEventMixer::isTpcProton(StPicoTrack const * const trk)
{
    return( isTPCHadron(trk, StPicoCutsBase::kProton));
}
//-----------------------------------------------------------
bool StPicoEventMixer::isTPCHadron(StPicoTrack const * const trk, int pidFlag)
{
    return( mHFCuts->isTPCHadron(trk, pidFlag));
}
//-----------------------------------------------------------
bool StPicoEventMixer::isGoodTrack(StPicoTrack const * const trk)
{
    return (mHFCuts->isGoodTrack(trk));
}
//-----------------------------------------------------------
bool StPicoEventMixer::isGoodTrigger(StPicoEvent const * const mPicoEvent) const
{
    return mHFCuts->isGoodTrigger(mPicoEvent);
}
