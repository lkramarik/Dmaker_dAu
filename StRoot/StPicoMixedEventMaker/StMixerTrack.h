#ifndef StMixerTrack_hh
#define StMixerTrack_hh
/* **************************************************
 *
 * Track class used for mixed event buffer, stripped down 
 * to minimum information neede to reconstruct the helix
 * and basic track information. Currently include:
 * 1) charge
 * 2) isTpcPi & isTofPi
 * 3) isTpcKaon & is TofKaon
 *
 * **************************************************
 *
 *  Initial Authors:  
 *            Michael Lomnitz (mrlomnitz@lbl.gov)
 *            Mustafa Mustafa (mmustafa@lbl.gov)
 *         ** Miroslav Simko  (msimko@bnl.gov)
 *
 *  ** Code Maintainer
 *
 * **************************************************
 */

#include "TVector3.h"
#include "TLorentzVector.h"

class StPicoTrack;

class StMixerTrack
        {
 public:
  StMixerTrack();
  StMixerTrack(StMixerTrack const *);
  StMixerTrack(TVector3 const & pVtx, float B,StPicoTrack const&, bool, bool, bool, bool, bool, bool);
  short const getTrackInfo() const;
  int const charge() const ;
  TVector3 const& gMom() const;
  TVector3 const& origin() const;
  ~StMixerTrack(){}
 private:
  TVector3 mOrigin;
  TVector3 mMom;
  short mTrackInfo;
};
inline short const StMixerTrack::getTrackInfo() const { return(mTrackInfo); }
inline TVector3 const & StMixerTrack::gMom() const { return(mMom) ;}
inline TVector3 const & StMixerTrack::origin() const { return(mOrigin) ;}
inline int const StMixerTrack::charge() const { 
  int temp = (mTrackInfo & 1);
  if(temp == 1) return 1;
  else return -1;
}
#endif
