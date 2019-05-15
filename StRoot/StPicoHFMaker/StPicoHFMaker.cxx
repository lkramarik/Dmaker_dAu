//lukas
#include <vector>

#include "TTree.h"
#include "TFile.h"
#include "TChain.h"

#include "StPicoEvent/StPicoDst.h"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
//#include "StPicoPrescales/StPicoPrescales.h"

#include "StHFCuts.h"
//#include "StHFHists.h"
#include "StPicoHFEvent.h"
#include "StPicoHFMaker.h"
#include "StHFPair.h"

ClassImp(StPicoHFMaker)

// _________________________________________________________
StPicoHFMaker::StPicoHFMaker(char const* name, StPicoDstMaker* picoMaker,
                             char const* outputBaseFileName) :
        StMaker(name), mPicoDst(NULL), mHFCuts(NULL), mPicoHFEvent(NULL), mBField(0.), mOutList(NULL),
        mOutputFileBaseName(outputBaseFileName),
        mPicoDstMaker(picoMaker), mPicoEvent(NULL),
        mOutputFileList(NULL) {
  // -- constructor
}

// _________________________________________________________
StPicoHFMaker::~StPicoHFMaker() {
  // -- destructor

  if (mHFCuts)
    delete mHFCuts;
  mHFCuts = NULL;

  /* mTree is owned by mOutputFile directory, it will be destructed once
   * the file is closed in ::Finish() */
}

// _________________________________________________________
Int_t StPicoHFMaker::Init() {
  // -- Inhertited from StMaker 
  //    NOT TO BE OVERWRITTEN by daughter class
  //    daughter class should implement InitHF()
  // -- check for cut class
  if (!mHFCuts)    mHFCuts = new StHFCuts;
  mHFCuts->init();

  // -- file with outputs
  mOutputFileList = new TFile(Form("%s.%s.root", mOutputFileBaseName.Data(), GetName()), "RECREATE");
  mOutputFileList->SetCompressionLevel(1);

  // -- disable automatic adding of objects to file
  bool oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(false);

  // -- add list which holds all histograms  
  mOutList = new TList();
  mOutList->SetName(GetName());
  mOutList->SetOwner(true);

  // -- create event stat histograms
  initializeEventStats();

  // -- call method of daughter class
  InitHF();
  TH1::AddDirectory(oldStatus);
  // -- reset event to be in a defined state
  resetEvent();

  return kStOK;
}

// _________________________________________________________
Int_t StPicoHFMaker::Finish() {
  mOutputFileList->cd();
  mOutList->Write(mOutList->GetName(),  TObject::kSingleKey); //predtym TObject::kSingleKey

  // -- call method of daughter class
  FinishHF();

  mOutputFileList->Close();

  return kStOK;
}

// _________________________________________________________
void StPicoHFMaker::resetEvent() {
  // -- reset event
//  mIdxPicoPions.clear();
//  mIdxPicoKaons.clear();
//  mIdxPicoProtons.clear();

//  mPicoHFEvent->clear("C");
}

// _________________________________________________________
void StPicoHFMaker::Clear(Option_t *opt) {
  // -- Inhertited from StMaker
  //    NOT TO BE OVERWRITTEN by daughter class
  //    daughter class should implement ClearHF()
  // -- call method of daughter class
  ClearHF();
//  resetEvent();
}

// _________________________________________________________
Int_t StPicoHFMaker::Make() {
// -- Inhertited from StMaker 
  //    NOT TO BE OVERWRITTEN by daughter class
  //    daughter class should implement MakeHF()

  if (!mPicoDstMaker) {
    LOG_WARN << " StPicoHFMaker - No PicoDstMaker! Skip! " << endm;
    return kStWarn;
  }

  mPicoDst = mPicoDstMaker->picoDst();
  if (!mPicoDst) {
    LOG_WARN << " StPicoHFMaker - No PicoDst! Skip! " << endm;
    return kStWarn;
  }

  Int_t iReturn = kStOK;

  if (setupEvent()) {
//    UInt_t nTracks = mPicoDst->numberOfTracks();
//    for (unsigned short iTrack = 0; iTrack < nTracks; ++iTrack) {
//        StPicoTrack* trk = mPicoDst->track(iTrack);
//
//        if (!trk || !mHFCuts->isGoodTrack(trk)) continue;
//        if (isPion(trk))   mIdxPicoPions.push_back(iTrack);   // isPion implemented by daughter class
//        if (isKaon(trk))   mIdxPicoKaons.push_back(iTrack);   // isKaon implemented by daughter class
//        if (isProton(trk)) mIdxPicoProtons.push_back(iTrack); // isProton method to be implemented by daughter class
//    }

    // -- call method of daughter class
    iReturn = MakeHF();
  }

  // -- reset event to be in a defined state
//  resetEvent();
  return (kStOK && iReturn);
}

// _________________________________________________________
bool StPicoHFMaker::setupEvent() {
  // -- fill members from pico event, check for good eventa and fill event statistics

  mPicoEvent = mPicoDst->event();
  mBField = mPicoEvent->bField();
  mPrimVtx = mPicoEvent->primaryVertex();

  int aEventStat[mHFCuts->eventStatMax()];
  bool bResult = mHFCuts->isGoodEvent(mPicoDst, aEventStat);

  fillEventStats(aEventStat);

  return bResult;
}

// _________________________________________________________
void StPicoHFMaker::initializeEventStats() {
  // -- Initialize event statistics histograms

  const char *aEventCutNames[]   = {"all", "good run", "trigger", "#it{v}_{z}", "#it{v}_{z}-#it{v}^{VPD}_{z}", "accepted"};

  mOutList->Add(new TH1F("hEventStat0","Event cut statistics 0;Event Cuts;Events", mHFCuts->eventStatMax(), -0.5, mHFCuts->eventStatMax()-0.5));
  TH1F *hEventStat0 = static_cast<TH1F*>(mOutList->Last());

  mOutList->Add(new TH1F("hEventStat1","Event cut statistics 1;Event Cuts;Events", mHFCuts->eventStatMax(), -0.5, mHFCuts->eventStatMax()-0.5));
  TH1F *hEventStat1 = static_cast<TH1F*>(mOutList->Last());

  for (unsigned int ii = 0; ii < mHFCuts->eventStatMax(); ii++) {
    hEventStat0->GetXaxis()->SetBinLabel(ii+1, aEventCutNames[ii]);
    hEventStat1->GetXaxis()->SetBinLabel(ii+1, aEventCutNames[ii]);
  }
}

//________________________________________________________________________
void StPicoHFMaker::fillEventStats(int *aEventStat) {
  // -- Fill event statistics 

  TH1F *hEventStat0 = static_cast<TH1F*>(mOutList->FindObject("hEventStat0"));
//  TH1F *hEventStat1 = static_cast<TH1F*>(mOutList->FindObject("hEventStat1"));

  for (unsigned int idx = 0; idx < mHFCuts->eventStatMax() ; ++idx) {
    if (!aEventStat[idx])
      hEventStat0->Fill(idx);
  }

//  for (unsigned int idx = 0; idx < mHFCuts->eventStatMax(); ++idx) {
//    if (aEventStat[idx])
//      break;
//    hEventStat1->Fill(idx);
//  }
}

