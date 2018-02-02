#ifndef StPicoHFMaker_h
#define StPicoHFMaker_h

#include "StChain/StMaker.h"
#include "StarClassLibrary/StLorentzVectorF.hh"

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
    StPicoHFMaker(char const* name, StPicoDstMaker* picoMaker, char const* outputBaseFileName,  
		       char const* inputHFListHFtree);
    virtual ~StPicoHFMaker();
    
    // -- TO BE IMPLEMENTED BY DAUGHTER CLASS
    virtual Int_t InitHF()                  { return kStOK; }
    virtual Int_t MakeHF()                  { return kStOK; }
    virtual void  ClearHF(Option_t *opt="") { return; }
    virtual Int_t FinishHF()                { return kStOK; }

    void setHFBaseCuts(StHFCuts* cuts);

    void setTreeName(const char* tName);

    void setMakerMode(unsigned short us);
    void setDecayMode(unsigned short us);
    void setMcMode(bool b);

    // -- different modes to use the StPicoHFMaker class
    //    - kAnalyze - don't write candidate trees, just fill histograms
    //    - kWrite   - write candidate trees
    //    - kRead    - read candidate trees and fill histograms
    enum eMakerMode {kAnalyze, kWrite, kRead};

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

    unsigned int isDecayMode() const;
    unsigned int isMakerMode() const;
    bool         isMcMode() const;
    
    // -- protected members ------------------------

    StPicoDst      *mPicoDst;

    StHFCuts       *mHFCuts;
//    StHFHists      *mHFHists; comm. LK 04012018

    StPicoHFEvent  *mPicoHFEvent;
    StPicoEvent*    mPicoEvent;          // ptr to picoDstEvent

    float           mBField;
    StThreeVectorF  mPrimVtx;

    TList          *mOutList;

    std::vector<unsigned short> mIdxPicoPions;
    std::vector<unsigned short> mIdxPicoKaons;
    std::vector<unsigned short> mIdxPicoProtons;

  private:
    void  resetEvent();
    bool  setupEvent();
    
    void  initializeEventStats();
    void  fillEventStats(int *aEventStat);

    // -- private members ------------------------

    unsigned int    mDecayMode;          // use enum of StPicoHFEvent::eHFEventMode
    unsigned int    mMakerMode;          // use enum of StPicoEventMaker::eMakerMode

    bool            mMcMode;             // use MC mode

    TString         mOutputTreeName;     // name for output trees

    TString         mOutputFileBaseName; // base name for output files
                                         //   for tree     -> <mOutputFileBaseName>.<mOutputTreeName>.root
                                         //   for histList -> <mOutputFileBaseName>.GetName().root

    TString         mInputFileName;      // filename of input list of HF trees (needs to be in the 
                                         // same order as the picoDstList

    StPicoDstMaker* mPicoDstMaker;       // ptr to picoDst maker


    TTree*          mTree;               // tree holding "mPicoHFEvent" for writing only

    TChain*         mHFChain;            // chain to read in HF tree
    int             mEventCounter;       // n Processed events in chain

    TFile*          mOutputFileTree;     // ptr to file saving the HFtree
    TFile*          mOutputFileList;     // ptr to file saving the list of histograms
    ClassDef(StPicoHFMaker, 0)
};

inline void StPicoHFMaker::setHFBaseCuts(StHFCuts* cuts)   { mHFCuts = cuts; }
inline void StPicoHFMaker::setTreeName(const char* tName)  { mOutputTreeName = tName; }

inline void StPicoHFMaker::setMakerMode(unsigned short us) { mMakerMode = us; }
inline void StPicoHFMaker::setDecayMode(unsigned short us) { mDecayMode = us; }
inline void StPicoHFMaker::setMcMode(bool b)               { mMcMode = b; }

inline unsigned int StPicoHFMaker::isDecayMode() const     { return mDecayMode; }
inline unsigned int StPicoHFMaker::isMakerMode() const     { return mMakerMode; }
inline bool StPicoHFMaker::isMcMode() const                { return mMcMode; }
#endif
