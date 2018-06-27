#ifndef StPicoMixedEventMaker_h
#define StPicoMixedEventMaker_h

#include "StMaker.h"

/* **************************************************
 *
 *  Base class for Mixed Event cosntructions
 *  Template implemented for D0 recosntruction. User should use a 
 *  Mixer per category in Event Mixing and define event buffer size (10 by default).
 *  For different decays changes must be made to StPicoEventMixer class
 * 
 *
 * **************************************************
 *
 *  Initial Authors:
 *        **  Michael Lomnitz  (mrlomnitz@lbl.gov)
 *            Mustafa Mustafa  (mmustafa@lbl.gov)
 *
 *  ** Code Maintainer
 * 
 *
 * **************************************************  
 */

class TTree;
class TNtuple;
class TFile;
class TChain;
class TList;

class StPicoDst;
class StPicoDstMaker;
class StPicoEvent;
class StPicoTrack;
class StHFCuts;

class StPicoEventMixer;

class StPicoMixedEventMaker : public StMaker 
{
  public:
  StPicoMixedEventMaker(char const* name, StPicoDstMaker* picoMaker, StHFCuts* hfCuts,
			char const* outputBaseFileName,  
			char const* inputHFListHFtree);
    ~StPicoMixedEventMaker();
    Int_t Init();
    Int_t Make();
    Int_t Finish();
    void  Clear(Option_t* opt="");

    Int_t SetCategories();
    void setBufferSize(int size) { mBufferSize = size; }

    enum  mixerConst { defaultBufferSize = 5, fillSingleTrackHistos = 1}; // enum trick to setting class-speciffic constants
 private:
    int categorize(StPicoDst const*);
    void fillCentralities();
    StPicoDst*      mPicoDst;
    StPicoDstMaker* mPicoDstMaker;      
    StPicoEvent*    mPicoEvent;         
    StHFCuts*	    mHFCuts;

    TNtuple *ntp_signal_SE;
    TNtuple *ntp_signal_ME;
    TNtuple *ntp_background_SE;
    TNtuple *ntp_background_ME;

    StPicoEventMixer* mPicoEventMixer[10][9];

    TString         mOuputFileBaseName; 
    TString         mInputFileName;     

    int             mEventCounter;
    int		    mBufferSize;

    bool loadEventPlaneCorr(int const runId);
                                        
    TNtuple*        mSETuple;
    TNtuple*        mMETuple;
    TFile*          mOutputFileTree; 

    TList* mSingePartHists;

    ClassDef(StPicoMixedEventMaker, 0)
};
#endif
