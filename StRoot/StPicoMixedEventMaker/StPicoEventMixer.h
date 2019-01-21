#ifndef StPicoEventMixer_hh
#define StPicoEventMixer_hh

/* **************************************************
 *  
 * Class stores event buffer used in event mixing. Mixing
 * is done automatically once buffer reaches defined maximum.
 * Template provided used for D0 reconstruction, user should personalize 
 * mixEvent() method to cosntruct desired background.
 *
 * **************************************************
 * 
 * Initial Authors:
 *          Michael Lomnitz (mrlomnitz@lbl.gov)
 *          Mustafa Mustafa   (mmustafa@lbl.gov)
 *  Other authors:
 *       ** Miroslav Simko  (msimko@bnl.gov)
 *
 *  ** Code maintainer 
 *
 * **************************************************
 */

#include <vector>
#include "TNtuple.h"
#include "TList.h"

#include "StPicoHFMaker/StHFCuts.h"

class TTree;
class TH2F;
class StPicoEvent;
class StPicoTrack;
class StPicoDst;
class StMixerTrack;
class StMixerEvent;
class StMixerPair;
class StMixerHists;
class StHFCuts;

class StPicoEventMixer {
 public: 
  enum bitNumber { kPionTPCbit, kPionTOFbit, kKaonTPCbit, kKaonTOFbit, kProtonTPCbit, kProtonTOFbit};

  StPicoEventMixer(char* category);
  ~StPicoEventMixer();
  bool addPicoEvent(StPicoDst const* picoDst, float weight = 1.);
  void setEventBuffer(int buffer);
  void setSameEvtNtupleSig(TNtuple *tuple)  { mSETupleSig = tuple; }
  void setSameEvtNtupleBack(TNtuple *tuple)  { mSETupleBack = tuple; }
  void setMixedEvtNtupleSig(TNtuple *tuple) { mMETupleSig = tuple; }
  void setMixedEvtNtupleBack(TNtuple *tuple) { mMETupleBack = tuple; }
  void mixEvents();
  void finish();
  void setHFCuts(StHFCuts * hfCuts) { mHFCuts = hfCuts; }
  StHFCuts * getHFCuts() { return mHFCuts; }

//  void setFillSinglePartHists(bool yesOrNo) { fillSinglePartHists = yesOrNo; }
//  bool isFillingSinglePartHists() { return fillSinglePartHists; }
//  void setSinglePartHistsList(TList *list) { mSingleParticleList = list; }

 private:
  bool isMixerPion(StMixerTrack const&);
  bool isMixerKaon(StMixerTrack const&);
  bool isMixerProton(StMixerTrack const&);
  bool isGoodEvent(StPicoDst const * const picoDst);
  bool isGoodTrigger(StPicoEvent const * const) const;
  bool isGoodTrack(StPicoTrack const * const trk);
  bool isTpcPion(StPicoTrack const * const);
  bool isTpcKaon(StPicoTrack const * const);
  bool isTpcProton(StPicoTrack const * const);
  bool isTPCHadron(StPicoTrack const * const, int pidFlag);
  bool isCloseMixerPair(StMixerPair const * pair) const;

//  void fillTracks(StMixerEvent* evt, bool isSameEvt, int PidFlag);
  void fillNtpSameEvtPair(float ntVar[], int charge);
  void fillNtpMixedEvtPair(float ntVar[], int charge);

  std::vector <StMixerEvent*> mEvents;
//  StMixerHists* mHists;
  StHFCuts *mHFCuts;
  unsigned short int mEventsBuffer; 
  unsigned short int filledBuffer;
  float dca1, dca2, dcaDaughters, theta_hs, decayL_hs;
  float pt_hs, mass_hs, eta_hs, phi_hs;

  TNtuple *mSETupleSig;
  TNtuple *mSETupleBack;
  TNtuple *mMETupleSig;
  TNtuple *mMETupleBack;
};

inline void StPicoEventMixer::setEventBuffer(int buffer){ mEventsBuffer = buffer;}
			    
    
#endif
