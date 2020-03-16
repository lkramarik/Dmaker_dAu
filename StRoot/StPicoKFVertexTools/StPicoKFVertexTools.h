#ifndef StPicoKFVertexTools_h
#define StPicoKFVertexTools_h

#include "StPicoHFMaker/StPicoHFMaker.h"
#include "TNtuple.h"
#include "TVector3.h"

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
class StPicoHFMaker;

class StPicoKFVertexTools : public StPicoHFMaker
{
public:
    StPicoKFVertexTools(char const*, StPicoDstMaker*, char const*);
    virtual ~StPicoKFVertexTools();

    virtual Int_t InitHF();
    virtual Int_t MakeHF();
    virtual void  ClearHF(Option_t *opt);
    virtual Int_t FinishHF();

protected:

private:
    TNtuple *ntp_vertex;
    TNtuple *ntp_KFReso;
    TString mOutFileBaseName;
    TFile* mOutFile;
    void makeKFReso(std::vector<int>&, int);
    void compareFitters(std::vector<int>&, int, int);

    ClassDef(StPicoKFVertexTools, 1) //set to 1

};

#endif
