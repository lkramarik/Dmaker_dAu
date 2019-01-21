#include <limits>

#include "StMixerTrack.h"
#include "StPicoEvent/StPicoTrack.h"

StMixerTrack::StMixerTrack() : 
  mOrigin(TVector3()),
  mMom(TVector3()),
  mTrackInfo(std::numeric_limits<short>::min())
//  StPicoTrack()
{
}
StMixerTrack::StMixerTrack(TVector3 const & pVtx, float B, StPicoTrack const& picoTrack, bool isTpcPi, bool isTofPi, bool isTpcK, bool isTofK, bool isTpcP, bool isTofP) :
    mOrigin(), 
    mMom(picoTrack.gMom(pVtx,B)), 
    mTrackInfo(0)
//    StPicoTrack()
{
  StPicoPhysicalHelix helix = picoTrack.helix(B);
  helix.moveOrigin(helix.pathLength(pVtx));
  mOrigin = helix.origin();

  if( picoTrack.charge() == 1 ) mTrackInfo = mTrackInfo | 1;
  //Pi
  if( isTpcPi == true ) mTrackInfo = mTrackInfo | (1 << 1);
  if( isTofPi == true ) mTrackInfo = mTrackInfo | (1 << 2);
  //K
  if( isTpcK == true ) mTrackInfo = mTrackInfo | (1 << 3);
  if( isTofK == true ) mTrackInfo = mTrackInfo | (1 << 4);
  //p
  if( isTpcP == true ) mTrackInfo = mTrackInfo | (1 << 5);
  if( isTofP == true ) mTrackInfo = mTrackInfo | (1 << 6);

}
StMixerTrack::StMixerTrack(StMixerTrack const * t) : 
  mOrigin(t->mOrigin), 
  mMom(t->mMom), 
  mTrackInfo(t->mTrackInfo)
{
}

