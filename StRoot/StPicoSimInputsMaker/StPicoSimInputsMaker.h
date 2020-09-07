#ifndef StPicoSimInputsMaker_h
#define StPicoSimInputsMaker_h

#include "StPicoHFMaker/StPicoHFMaker.h"
#include "TNtuple.h"
#include <vector>
#include "TClonesArray.h"
#include "StLorentzVectorF.hh"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StPicoHFMaker/StPicoHFEvent.h"
#include "StPicoHFMaker/StHFCuts.h"
//#include "StPicoHFMaker/StHFTriplet.h"
#include "StBTofUtil/tofPathLength.hh"
#include "StAnaCuts.h"
#include "phys_constants.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TVector3.h"
#include <ctime>

class StPicoDst;
class StPicoDstMaker;
class StPicoEvent;
class StPicoTrack;
class StPicoHFEvent;
class StHFCuts;

class StPicoSimInputsMaker : public StPicoHFMaker
{
public:
    StPicoSimInputsMaker(char const*, StPicoDstMaker*, char const*);
    virtual ~StPicoSimInputsMaker();

    virtual Int_t InitHF();
    virtual Int_t MakeHF();
    virtual void  ClearHF(Option_t *opt);
    virtual Int_t FinishHF();
    virtual double DCA(StPicoTrack const*, TVector3 const &) const;
    int createQA();

   void histoInit(TString fileBaseName,bool fillQaHists=true);
   void addTpcDenom1(bool IsPion, bool IsKaon, float pt, int centrality, int Eta, int Phi, float Vz, float ZDC);
   void addHFTNumer1(bool IsPion, bool IsKaon, float pt, int centrality, int EtaIndex, int PhiIndex, float Vz, float Zdc);
   void addDcaPtCent(float dca, float dcaXy, float  dcaZ, bool IsPion, bool IsKaon, float pt,  int centrality, int Eta, int Phi, float Vz, int zdc);
   int getMultIndex(float multiplicity);
   int getEtaIndexDca(float Eta) ;
   int getZdcIndex(float zdc);
   int getPhiIndexDca(float Phi) ;
   int getVzIndexDca(float Vz) ;
   int getPtIndexDca(float pt);
   int getEtaIndexRatio(float Eta) ;
   int getPhiIndexRatio(float Phi) ;
   int getVzIndexRatio(float Vz) ;
//   void addCent(const double refmultCor, int centrality, const double reweight, const float vz);
   void closeFile();

private:
    TString mOutFileBaseName;

    bool mFillQaHists;
    TFile* mOutFileRatio;
    TFile* mOutFileDCA;
    TFile* mOutFileTuple;
    TFile* mOutFileTOFRatio;

    TH1F* mh1Cent;
    TH1F* mh1CentWg;
    TH1F* mh1gRefmultCor;
    TH1F* mh1gRefmultCorWg;
    TH2F* mh2CentVz;
    TH2F* mh2CentVzWg;

    //HFT ratio QA
    TH2F* mh2Tpc1PtCent;
    TH2F* mh2Tpc1PhiVz;
    TH2F* mh2HFT1PtCent;
    TH2F* mh2HFT1PhiVz;

    TH1F* mh2Tpc1PtCentPartEtaVzPhi[vars::m_nParticles][vars::m_nEtasRatio][vars::m_nVzsRatio][vars::m_nPhisRatio][vars::m_nmultEdge];
    TH1F* mh2HFT1PtCentPartEtaVzPhi[vars::m_nParticles][vars::m_nEtasRatio][vars::m_nVzsRatio][vars::m_nPhisRatio][vars::m_nmultEdge];

    TH1D* h1Tofmatch[vars::m_nParticles][3];
    TH1D* h1TofmatchTOF[vars::m_nParticles][3];

    //HFT Dca
    TH2F* mh3DcaXyZPtCentPartEtaVzPhi[vars::m_nParticles][vars::m_nEtasDca][vars::m_nVzsDca][vars::m_nmultEdge][vars::m_nPtsDca];

    TH3F* mh3DcaPtCent;
    TH3F* mh3DcaXyPtCent;
    TH3F* mh3DcaZPtCent;

    TH3F* mh3VzZdcMult;

    TH1D* mhVx;
    TH1D* mhVy;
    TH1D* mhVz;
    TH2D* mhVxVy;

    TNtuple *ntp_tracks;

    ClassDef(StPicoSimInputsMaker, 1) //set to 1
};

#endif
