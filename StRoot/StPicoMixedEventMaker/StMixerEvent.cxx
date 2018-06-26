#include <limits>

#include "StMixerEvent.h"
#include "StPicoDstMaker/StPicoEvent.h"

StMixerEvent::StMixerEvent() :  mVtx(StThreeVectorF()),
    mBField(std::numeric_limits<float>::quiet_NaN()),
    mWeightFromCentrality(1)
{
}
StMixerEvent::StMixerEvent(StMixerEvent *t) : mVtx(t->mVtx), mBField(t->mBField),
					      mTracks(t->mTracks),
					      mEventKaons(t->mEventKaons), mEventPions(t->mEventPions), mEventProtons(t->mEventProtons),
					      mWeightFromCentrality(1)
{
}
StMixerEvent::StMixerEvent(StThreeVectorF vtx, float b) :  mVtx(StThreeVectorF()),
    mBField(std::numeric_limits<float>::quiet_NaN()),
    mWeightFromCentrality(1)
{
    mVtx = vtx;
    mBField = b;

}
void StMixerEvent::addTrack(StMixerTrack t)
{
  mTracks.push_back(t);
}
void StMixerEvent::addPion(int arrayId)
{
  mEventPions.push_back(arrayId);
}
void StMixerEvent::addKaon(int arrayId)
{
  mEventKaons.push_back(arrayId);
}
void StMixerEvent::addProton(int arrayId)
{
  mEventProtons.push_back(arrayId);
}
void StMixerEvent::addPicoEvent(StPicoEvent const & event)
{
  mEventId = event.eventId();
  mRunId = event.runId();
}
