#ifndef StPicoDpmAnaMaker_h
#define StPicoDpmAnaMaker_h

#include "StPicoHFMaker/StPicoHFMaker.h"
#include "TNtuple.h"
//#include "StRefMultCorr/StRefMultCorr.h"
#include "TH2F.h"
//#include "StPicoDpmAnaHists.h"
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

/* **************************************************
 *  Sample class fo HF picoDST analysis
 * --------------------------------------------------
 * 
 *  For more info look also in the .h files in StPicoHFMaker/
 *     StPicoHFMaker/StPicoHFMaker.h      <-- Base Class for analysis
 *     StPicoHFMaker/StPicoHFEvent.h      <-- Holds candidates for one event (written to Tree)
 *     StPicoHFMaker/StHFCuts.h           <-- Cuts, can be set in run macro
 *     StPicoHFMaker/StHFPair.h           <-- Holds a pair candidate of a two body decay
 *     StPicoHFMaker/StHFTriplet.h        <-- Holds a triplet of a three body decay
 *
 *  Usage:
 *   - Implement
 *        InitHF()
 *        MakeHF()
 *        ClearHF()
 *        FinishHF()
 *
 *  - Do not ovewrite Init, Make, Clear, Finish which are inhertited from StPicoHFMaker via StMaker 

 *  - Set StHFCuts class via setHFBaseCuts(...) in run macro
 *
 *  - Set use mode of StPicoHFMaker class  via setMakerMode(...)
 *     use enum of StPicoHFMaker::eMakerMode
 *      StPicoHFMaker::kAnalyze - don't write candidate trees, just fill histograms
 *      StPicoHFMaker::kWrite   - write candidate trees
 *      StPicoHFMaker::kRead    - read candidate trees and fill histograms
 *
 *  - Set decay mode of analysis via setDecayMode(...)
 *     use enum of StPicoHFEvent::eHFEventMode (see there for more info)
 *      StPicoHFEvent::kTwoParticleDecay,
 *      StPicoHFEvent::kThreeParticleDecay
 *      StPicoHFEvent::kTwoAndTwoParticleDecay
 *
 *  - Implement these track selection methods used to fill vectors for 'good' identified particles
 *      (methods from StHFCuts utility class can/should be used)
 *       isPion
 *       isKaon
 *       isProton
 *
 *  --------------------------------------------------
 *  
 *  Initial Authors:  
 *            Xin Dong        (xdong@lbl.gov)
 *            Mustafa Mustafa (mmustafa@lbl.gov)
 *          **Jochen Thaeder  (jmthader@lbl.gov) 
 * 
 *  ** Code Maintainer
 *
 * **************************************************
 */

class StPicoDst;
class StPicoDstMaker;
class StPicoEvent;
class StPicoTrack;
class StPicoHFEvent;

class StHFPair;
class StHFTriplet;
class StHFCuts;

class StRefMultCorr;

class StPicoDpmAnaMaker : public StPicoHFMaker 
{
 public:
  StPicoDpmAnaMaker(char const* name, StPicoDstMaker* picoMaker, char const* outputBaseFileName,  
		       char const* inputHFListHFtree);
  virtual ~StPicoDpmAnaMaker();
  
  virtual Int_t InitHF();
  virtual Int_t MakeHF();
  virtual void  ClearHF(Option_t *opt);
  virtual Int_t FinishHF();
  // -- Lomnitz: Added this cut funtions to to filter iwthout having to make pairs
  virtual bool isCloseTracks(StPicoTrack const*, StPicoTrack const*,StThreeVectorF const & , float) const;
  virtual double DCA(StPicoTrack const*, StThreeVectorF const &) const;
  int createQA();
  
  // -- ADOPT DECAY CHANNELS, if wished ------------------- 
  void setDecayChannel(unsigned int u) { mDecayChannel = u; }

  enum eDecayChannel {kChannel1, kChannel2, kChannel3};

  //void setRefMutCorr(StRefMultCorr* gRefMultCorr) { mRefmultCorrUtil = gRefMultCorr; }
  //StRefMultCorr* getRefMultCorr() { return mRefmultCorrUtil; }

 //  float getTofBeta(StPicoTrack const*,StThreeVectorF const& vtx) const;

   void histoInit(TString fileBaseName,bool fillQaHists=true);
   void addTpcDenom1(bool IsPion, bool IsKaon, bool IsProton, float pt, int centrality, float Eta, float Phi, float Vz);
   void addHFTNumer1(bool IsPion, bool IsKaon, bool IsProton, float pt, int centrality, float Eta, float Phi, float Vz);
   void addDcaPtCent(float dca, float dcaXy, float  dcaZ, bool IsPion, bool IsKaon, bool IsProton, float pt,  int centrality, float Eta, float Phi, float Vz);
   int getEtaIndexDca(float Eta) ;
   int getPhiIndexDca(float Phi) ;
   int getVzIndexDca(float Vz) ;
   int getEtaIndexRatio(float Eta) ;
   int getPhiIndexRatio(float Phi) ;
   int getVzIndexRatio(float Vz) ;
   void addCent(const double refmultCor, int centrality, const double reweight, const float vz);
   void closeFile();

 // virtual float getEta(int index){return m_EtaEdgeDca[index];};


 //  ClassDef(StPicoDpmAnaMaker, 1)

 protected:
  virtual bool isHadron(StPicoTrack const*, int pidFlag) const;
  virtual bool isPion(StPicoTrack const*) const;
  virtual bool isKaon(StPicoTrack const*) const;
  virtual bool isProton(StPicoTrack const*) const;

private:
  int createCandidates();
  int analyzeCandidates();



  // -- private members --------------------------

  unsigned int mDecayChannel;


  // -- ADD USER MEMBERS HERE ------------------- 
   TNtuple *ntp_DMeson_Signal;
   TNtuple *ntp_DMeson_Background;


  // StRefMultCorr* mRefmultCorrUtil;
   int mRunNumber;
       
TString mOutFileBaseName;

  bool mFillQaHists;
   TFile* mOutFile;

   //Cuts----------------------------
   static const int m_nParticles = 3;
   //TString m_ParticleName[m_nParticles];

   static const int m_nEtasDca = 5;
  // float m_EtaEdgeDca[m_nEtasDca+1];
   static const int m_nPhisDca = 11;
   //static float m_PhiEdgeDca[m_nPhisDca + 1];

   static const int m_nVzsDca = 4;
   //static float m_VzEdgeDca[m_nVzsDca + 1];

   static const int m_nCentsDca = 9;
   //static float m_CentEdgeDca[m_nCentsDca + 1];

   static const int m_nPtsDca = 19;
  // static float m_PtEdgeDca[m_nPtsDca + 1];

   static const int m_nEtasRatio = 10;
  // static float m_EtaEdgeRatio[m_nEtasRatio + 1];

   static const int m_nPhisRatio = 11;
   //static float m_PhiEdgeRatio[m_nPhisRatio + 1];

   static const int m_nVzsRatio = 6;
   //static float m_VzEdgeRatio[m_nVzsRatio + 1];

   static const int m_nCentsRatio = 10;
  // static float m_CentEdgeRatio[m_nCentsRatio + 1];

   static const int m_nPtsRatio = 36;
  // static float m_PtEdgeRatio[m_nPtsRatio + 1];

  static const int m_nDcasDca = 144;
 // static float m_DcaEdgeDca[m_nDcasDca + 1];
   //-----------------------------------

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
   TH2F* mh2Tpc1PtCentPartEtaVzPhi[m_nParticles][m_nEtasRatio][m_nVzsRatio][m_nPhisRatio];
   TH2F* mh2HFT1PtCentPartEtaVzPhi[m_nParticles][m_nEtasRatio][m_nVzsRatio][m_nPhisRatio];

   //HFT Dca
   TH3F* mh3DcaXyZPtCentPartEtaVzPhi[m_nParticles][m_nEtasDca][m_nVzsDca][m_nCentsDca];

   TH3F* mh3DcaPtCent;
   TH3F* mh3DcaXyPtCent;
   TH3F* mh3DcaZPtCent;


  // -- ADD USER MEMBERS HERE -------------------

  ClassDef(StPicoDpmAnaMaker, 1) //set to 1
};

#endif
