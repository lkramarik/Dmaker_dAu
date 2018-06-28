#ifndef __CINT__
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

#else
class StChain;
#endif

class StPicoDstMaker;
class StPicoMixedEventMaker;
class StMaker;
StChain *chain;

void runPicoMixedEvent(
			const Char_t *inputFile="/gpfs01/star/pwg/lkramarik/Dmaker_dAu/picoLists/runs_local_test.list",	
			const Char_t *outputFile="outputBaseName",  
			const Char_t *badRunListFileName = "/gpfs01/star/pwg/lkramarik/Dmaker_dAu/picoLists/picoList_bad.list") {

  string SL_version = "SL17d";
  string env_SL = getenv ("STAR");
  if (env_SL.find(SL_version)==string::npos) {
      cout<<"Environment Star Library does not match the requested library in runPicoHFMyAnaMaker.C. Exiting..."<<endl;
      exit(1);
  }
  
  gROOT->LoadMacro("loadSharedHFLibraries.C");
  loadSharedHFLibraries();

  chain = new StChain();

  TString sInputFile(inputFile);
  TString sInputListHF("");

  if (!sInputFile.Contains(".list") && !sInputFile.Contains("picoDst.root")) {
    cout << "No input list or picoDst root file provided! Exiting..." << endl;
    exit(1);
  }

  StHFCuts* hfCuts = new StHFCuts("hfBaseCuts");
  cout<<"event stuff set"<<endl;

  hfCuts->setBadRunListFileName(badRunListFileName);
  hfCuts->setCutVzMax(60.);
  hfCuts->setCutVzVpdVzMax(30.);
  hfCuts->addTriggerId(530003); //VPD-5

  hfCuts->setCutNHitsFitMin(15);
  hfCuts->setCutRequireHFT(true);

   // kaonPion pair cuts
  float dcaDaughtersMax = 0.2;  // maximum
  float decayLengthMin  = 0.000; // minimum
  float decayLengthMax  = 999999; //std::numeric_limits<float>::max();
  float cosThetaMin     = -20.;   // minimum
  float minMass         = 0.6;
  float maxMass         = 2.6;
  hfCuts->setCutSecondaryPair(dcaDaughtersMax, decayLengthMin, decayLengthMax, cosThetaMin, minMass, maxMass);
 
  //Single track pt
  hfCuts->setCutPtRange(0.2,50.0,StHFCuts::kPion); //0.2 , 50.0
  hfCuts->setCutPtRange(0.2,50.0,StHFCuts::kKaon); //0.2, 50.0
  //TPC setters
  hfCuts->setCutTPCNSigmaPion(3.0); //3
  hfCuts->setCutTPCNSigmaKaon(2.0); //3
  //TOF setters, need to set pt range as well
  hfCuts->setCutTOFDeltaOneOverBeta(0.05, StHFCuts::kKaon); // v podstate 5 sigma; nastavene = f * (sigmaTOF), sigma TOF je 0.013 
  hfCuts->setCutPtotRangeHybridTOF(0.2,50.0,StHFCuts::kKaon);
  hfCuts->setCutTOFDeltaOneOverBeta(0.06, StHFCuts::kPion); // v podstate 6 sigma
  hfCuts->setCutPtotRangeHybridTOF(0.2,50.0,StHFCuts::kPion); 

  StPicoDstMaker* picoDstMaker = new StPicoDstMaker(static_cast<StPicoDstMaker::PicoIoMode>(StPicoDstMaker::IoRead), inputFile, "picoDstMaker");
  cout<<"ok, picoDstMaker created"<<endl;
  StPicoMixedEventMaker* picoMixedEventMaker = new StPicoMixedEventMaker("picoMixedEventMaker", picoDstMaker, hfCuts, outputFile, "");

  clock_t start = clock(); // getting starting time
  chain->Init();
  
  cout << " Total entries = " << picoDstMaker->chain()->GetEntries() << endl;

//  for (Int_t i=0; i<nEvents; i++) {
  for (Int_t i=0; i<2000; i++) {
    if(i%10==0)       cout << "Working on eventNumber " << i << endl;
    chain->Clear();
    int iret = chain->Make(i);
    if (iret) { cout << "Bad return code!" << iret << endl; break;}
  }
  
  cout << "****************************************** " << endl;
  cout << "Work done... now its time to close up shop!"<< endl;
  cout << "****************************************** " << endl;
  chain->Finish();
  double duration = (double) (clock() - start) / (double) CLOCKS_PER_SEC;
  cout << "****************************************** " << endl;
  cout << "total number of events  " << picoDstMaker->chain()->GetEntries() << endl;
  cout << "****************************************** " << endl;
  cout << "Time needed " << duration << " s" << endl;
  cout << "****************************************** " << endl;
  
  delete chain;

}

