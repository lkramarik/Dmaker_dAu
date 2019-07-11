#include "TROOT.h"
#include "TSystem.h"
#include "TChain.h"
#include "StMaker.h"
#include "StChain.h"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoHFMaker/StPicoHFEvent.h"
#include "StPicoHFMaker/StHFCuts.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoMixedEventMaker/StPicoMixedEventMaker.h"
#include "macros/loadSharedHFLibraries.C"
#include <iostream>
#include <ctime>
#include <cstdio>
#include "StPicoD0AnaMaker/StPicoD0AnaMaker.h"
using namespace std;

void runPicoMixedEvent(
			const Char_t *inputFile="./picoLists/runs_local_test.list",
			const Char_t *outputFile="outputLocal",
			const Char_t *badRunListFileName = "./picoLists/picoList_bad.list") {

  string SL_version = "SL18f";
  string env_SL = getenv ("STAR");
  if (env_SL.find(SL_version)==string::npos) {
      cout<<"Environment Star Library does not match the requested library in runPicoHFMyAnaMaker.C. Exiting..."<<endl;
      exit(1);
  }
  
  gROOT->LoadMacro("loadSharedHFLibraries.C");
  loadSharedHFLibraries();
  StChain *chain = new StChain();

  TString sInputFile(inputFile);

  if (!sInputFile.Contains(".list") && !sInputFile.Contains("picoDst.root")) {
    cout << "No input list or picoDst root file provided! Exiting..." << endl;
    exit(1);
  }

  StHFCuts* hfCuts = new StHFCuts("hfBaseCuts");

  hfCuts->setBadRunListFileName(badRunListFileName);
  hfCuts->setCutVzMax(6.);
  hfCuts->setCutVzVpdVzMax(6.);
  hfCuts->addTriggerId(530003); //VPD-5

  hfCuts->setCutNHitsFitMin(15);
  hfCuts->setCutRequireHFT(true);
  hfCuts->setHybridTof(true);

  // kaonPion pair cuts
  float dcaDaughtersMax = 0.2;  // maximum
  float decayLengthMin  = 0.000; // minimum
  float decayLengthMax  = 999999; //std::numeric_limits<float>::max();
  float cosThetaMin     = -20.;   // minimum
  float minMass         = 1.7;
  float maxMass         = 2.1;
  float pairDcaMax      = 99.9;

  hfCuts->setCutSecondaryPair(dcaDaughtersMax, decayLengthMin, decayLengthMax, cosThetaMin, minMass, maxMass, pairDcaMax);

  hfCuts->setCutTPCNSigmaPion(3.0);
  hfCuts->setCutTPCNSigmaKaon(2.0);
  hfCuts->setCutTOFDeltaOneOverBetaKaon(0.03);
  hfCuts->setCutTOFDeltaOneOverBetaPion(0.03);
  hfCuts->setCutPtMin(0.15);

  StPicoDstMaker* picoDstMaker = new StPicoDstMaker(static_cast<StPicoDstMaker::PicoIoMode>(StPicoDstMaker::IoRead), inputFile, "picoDstMaker");
  StPicoMixedEventMaker* picoMixedEventMaker = new StPicoMixedEventMaker("picoMixedEventMaker", picoDstMaker, hfCuts, outputFile);
  picoMixedEventMaker->setBufferSize(6);

  clock_t start = clock(); // getting starting time
  chain->Init();
  Int_t nEvents = picoDstMaker->chain()->GetEntries();
  cout << " Total entries = " << nEvents << endl;

  for (Int_t i=0; i<nEvents; ++i) {
//    if(i%10==0)       cout << "Working on eventNumber " << i << endl;
//    cout << "Working on eventNumber " << i << endl;
    chain->Clear();
    int iret = chain->Make(i);
    if (iret) { cout << "Bad return code!" << iret << endl; break;}
  }

  chain->Finish();
  double duration = (double) (clock() - start) / (double) CLOCKS_PER_SEC;
  cout << "****************************************** " << endl;
  cout << "Work done, total number of events  " << nEvents << endl;
  cout << "Time needed " << duration << " s" << endl;
  delete chain;
}

