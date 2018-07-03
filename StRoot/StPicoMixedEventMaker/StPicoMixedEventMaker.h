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

	TNtuple *mSETupleSig;
	TNtuple *mSETupleBack;
	TNtuple *mMETupleSig;
	TNtuple *mMETupleBack;
	static const int m_nmultEdge = 7;
	static float const m_multEdge[m_nmultEdge+1] = {0, 4, 8, 12, 16, 20, 24, 200};
	int getMultIndex(float multiplicity);

private:
    int categorize(StPicoDst const*);
    void fillCentralities();
    StPicoDst*      mPicoDst;
    StPicoDstMaker* mPicoDstMaker;      
    StPicoEvent*    mPicoEvent;         
    StHFCuts*	    mHFCuts;

    StPicoEventMixer* mPicoEventMixer[10][9];

    TString         mOuputFileBaseName; 
    TString         mInputFileName;     

    int             mEventCounter;
    int		    mBufferSize;

    bool loadEventPlaneCorr(int const runId);
                                        
    TFile*          mOutputFileTree;



    ClassDef(StPicoMixedEventMaker, 0)
};
#endif
