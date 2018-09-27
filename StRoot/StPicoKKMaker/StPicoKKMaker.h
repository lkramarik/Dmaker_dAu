#ifndef StPicoKKMaker_h
#define StPicoKKMaker_h

#include "StPicoHFMaker/StPicoHFMaker.h"
#include "TNtuple.h"
//#include "StRefMultCorr/StRefMultCorr.h"
#include "TH2F.h"
//#include "StPicoD0AnaHists.h"
#include <vector>
#include "TClonesArray.h"
#include "StThreeVectorF.hh"
#include "StLorentzVectorF.hh"
#include "StPicoDstMaker/StPicoDst.h"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StPicoHFMaker/StPicoHFEvent.h"
#include "StPicoHFMaker/StHFCuts.h"
#include "StPicoHFMaker/StHFPair.h"
//#include "StPicoHFMaker/StHFTriplet.h"
#include "StBTofUtil/tofPathLength.hh"

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

class StPicoKKMaker : public StPicoHFMaker
{
public:
    StPicoKKMaker(char const*, StPicoDstMaker*, char const*, char const*);
    virtual ~StPicoKKMaker();

    virtual Int_t InitHF();
    virtual Int_t MakeHF();
    virtual void  ClearHF(Option_t *opt);
    virtual Int_t FinishHF();

protected:
    virtual bool isHadron(StPicoTrack const*, int pidFlag) const;
    virtual bool isPion(StPicoTrack const*) const;
    virtual bool isKaon(StPicoTrack const*) const;
    virtual bool isProton(StPicoTrack const*) const;
    std::vector<unsigned short> mIdxPicoPions;

private:
    int createCandidates();

    TNtuple *ntp_signal;
    TNtuple *ntp_background;


    int mRunNumber;
    TString mOutFileBaseName;

    TFile* mOutFile;

    ClassDef(StPicoKKMaker, 1) //set to 1
};

#endif
