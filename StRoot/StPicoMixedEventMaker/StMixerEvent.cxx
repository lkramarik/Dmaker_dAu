#include <limits>

#include "StMixerEvent.h"
#include "StPicoEvent/StPicoEvent.h"

StMixerEvent::StMixerEvent() :  mVtx(TVector3()),
    mBField(std::numeric_limits<float>::quiet_NaN())
{
}
StMixerEvent::StMixerEvent(StMixerEvent *t) : mVtx(t->mVtx), mBField(t->mBField),
					      mTracks(t->mTracks),
					      mEventKaons(t->mEventKaons), mEventPions(t->mEventPions), mEventProtons(t->mEventProtons)
{
}
StMixerEvent::StMixerEvent(TVector3 vtx, float b) :  mVtx(TVector3()),
    mBField(std::numeric_limits<float>::quiet_NaN())
{
    mVtx = vtx;
    mBField = b;
}
void StMixerEvent::addTrack(StPicoTrack t)
//void StMixerEvent::addTrack(StMixerTrack t)
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
