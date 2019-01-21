#ifndef StPicoPiPiMaker_h
#define StPicoPiPiMaker_h

#include "StPicoHFMaker/StPicoHFMaker.h"
#include "TNtuple.h"
//#include "StRefMultCorr/StRefMultCorr.h"
#include "TH2F.h"
//#include "StPicoD0AnaHists.h"
#include <vector>
#include "TClonesArray.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StPicoHFMaker/StPicoHFEvent.h"
#include "StPicoHFMaker/StHFCuts.h"
#include "StPicoHFMaker/StHFPair.h"
//#include "StPicoHFMaker/StHFTriplet.h"

#include "phys_constants.h"

#include "TH1F.h"
#include "TH3F.h"
#include <ctime>

class StPicoDst;
class StPicoDstMaker;
class StPicoEvent;
class StPicoTrack;
class StPicoHFEvent;

class StHFPair;
class StHFTriplet;
class StHFCuts;

class StPicoPiPiMaker : public StPicoHFMaker
{
public:
    StPicoPiPiMaker(char const*, StPicoDstMaker*, char const*);
    virtual ~StPicoPiPiMaker();

    virtual Int_t InitHF();
    virtual Int_t MakeHF();
    virtual void  ClearHF(Option_t *opt);
    virtual Int_t FinishHF();

protected:
    std::vector<unsigned short> mIdxPicoPions;

private:
    int createCandidates();

    TNtuple *ntp_signal;
    TNtuple *ntp_background;


    int mRunNumber;
    TString mOutFileBaseName;

    TFile* mOutFile;

    ClassDef(StPicoPiPiMaker, 1) //set to 1
};

#endif
