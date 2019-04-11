#ifndef StPicoHFMaker_h
#define StPicoHFMaker_h

#include "TVector3.h"

#include "StChain/StMaker.h"

class TTree;
class TFile;
class TChain;

class StPicoDst;
class StPicoDstMaker;
class StPicoEvent;
class StPicoTrack;
class StPicoHFEvent;
class StHFPair;
class StHFCuts;
//class StHFHists;

class StPicoHFMaker : public StMaker 
{
  public:
    StPicoHFMaker(char const* , StPicoDstMaker* , char const*);
    virtual ~StPicoHFMaker();
    
    // -- TO BE IMPLEMENTED BY DAUGHTER CLASS
    virtual Int_t InitHF()                  { return kStOK; }
    virtual Int_t MakeHF()                  { return kStOK; }
    virtual void  ClearHF(Option_t *opt="") { return; }
    virtual Int_t FinishHF()                { return kStOK; }

    void setHFBaseCuts(StHFCuts* cuts);

    // -- TO BE IMPLEMENTED BY DAUGHTER CLASS
    virtual bool  isHadron(StPicoTrack const*, int pidFlag)   const { return true; }
    virtual bool  isPion(StPicoTrack const*)   const { return true; }
    virtual bool  isKaon(StPicoTrack const*)   const { return true; }
    virtual bool  isProton(StPicoTrack const*) const { return true; }

    // -- Inhertited from StMaker 
    //    NOT TO BE OVERWRITTEN by daughter class
    //    daughter class should implement xxxHF()
    //    -> will be declared as "final" when C++11 is used in STAR
    Int_t Init();
    Int_t Make();
    void  Clear(Option_t *opt="");
    Int_t Finish();

  protected:
    StPicoDst      *mPicoDst;
    StHFCuts       *mHFCuts;
    StPicoHFEvent  *mPicoHFEvent;
    StPicoDstMaker* mPicoDstMaker;       // ptr to picoDst maker
    StPicoEvent*    mPicoEvent;          // ptr to picoDstEvent

    float           mBField;
    TVector3        mPrimVtx;
    TList          *mOutList;

private:
    void  resetEvent();
    bool  setupEvent();
    TFile*          mOutputFileList;     // ptr to file saving the list of histograms

    void  initializeEventStats();
    void  fillEventStats(int *aEventStat);

    TString         mOutputFileBaseName; // base name for output files
    ClassDef(StPicoHFMaker, 0)
};

inline void StPicoHFMaker::setHFBaseCuts(StHFCuts* cuts)   { mHFCuts = cuts; }

#endif
