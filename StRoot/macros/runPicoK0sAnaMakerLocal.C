#include "TROOT.h"
#include "TSystem.h"
#include "TChain.h"
#include "StMaker.h"
#include "StChain.h"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoHFMaker/StPicoHFEvent.h"
#include "StPicoHFMaker/StHFCuts.h"
#include "StPicoEvent/StPicoEvent.h"
#include "macros/loadSharedHFLibraries.C"
#include <iostream>
#include <ctime>
#include <cstdio>
#include "StPicoPiPiMaker/StPicoPiPiMaker.h"

using namespace std;

void runPicoK0sAnaMakerLocal(
			const Char_t *inputFile="./picoLists/runs_local_test.list",
			const Char_t *outputFile="outputLocal",
			const Char_t *badRunListFileName = "./picoLists/picoList_bad.list") {
  string SL_version = "SL18f";
  string env_SL = getenv ("STAR");
  if (env_SL.find(SL_version)==string::npos) {
      cout<<"Environment Star Library does not match the requested library in runPicoHFMyAnaMaker.C. Exiting..."<<endl;
      exit(1);
  }
  
  Int_t nEvents = 1000000;

  StChain *chain = new StChain();

  TString sInputFile(inputFile);

  if (!sInputFile.Contains(".list") && !sInputFile.Contains("picoDst.root")) {
    cout << "No input list or picoDst root file provided! Exiting..." << endl;
    exit(1);
  }

  StHFCuts* hfCuts = new StHFCuts("hfBaseCuts");
  cout<<"event stuff set"<<endl;
  // ---------------------------------------------------
  // -- Set Base cuts for HF analysis

  // -- File name of bad run list
   hfCuts->setBadRunListFileName(badRunListFileName); 

  hfCuts->setCutVzMax(6.);
  hfCuts->setCutVzVpdVzMax(3.);
  hfCuts->addTriggerId(530003); //VPD-5

  hfCuts->setCutNHitsFitMin(15); //default is 20
  hfCuts->setCutRequireHFT(false);
  hfCuts->setHybridTof(true);
  //LK hfCuts->setCutDcaMin(0.009,StHFCuts::kPion);
  //LK  hfCuts->setCutDcaMin(0.007,StHFCuts::kKaon);
  //hfCuts->setCutNHitsFitnHitsMax(0.52);

   // kaonPion pair cuts
  float dcaDaughtersMax = 0.5;  // maximum
  float decayLengthMin  = 0.5; // minimum
  float decayLengthMax  = 999999; //std::numeric_limits<float>::max();
  float cosThetaMin     = 0.85;   // minimum
  float minMass         = 0.4;
  float maxMass         = 0.6;
  float pairDcaMax      = 99.9;

  hfCuts->setCutSecondaryPair(dcaDaughtersMax, decayLengthMin, decayLengthMax, cosThetaMin, minMass, maxMass, pairDcaMax);
 
  //Single track pt
  hfCuts->setCutPtRange(0.15,50.0,StHFCuts::kPion); //0.2 , 50.0
  hfCuts->setCutPtRange(0.15,50.0,StHFCuts::kKaon); //0.2, 50.0
  //TPC setters
  hfCuts->setCutTPCNSigmaPion(3.0); //3
  hfCuts->setCutTPCNSigmaKaon(2.0); //3
  //TOF setters, need to set pt range as well
  hfCuts->setCutTOFDeltaOneOverBeta(0.05, StHFCuts::kKaon); // v podstate 5 sigma; nastavene = f * (sigmaTOF), sigma TOF je 0.013 
  hfCuts->setCutPtotRangeHybridTOF(0.15,50.0,StHFCuts::kKaon);
  hfCuts->setCutTOFDeltaOneOverBeta(0.06, StHFCuts::kPion); // v podstate 6 sigma
  hfCuts->setCutPtotRangeHybridTOF(0.15,50.0,StHFCuts::kPion);

  StPicoDstMaker* picoDstMaker = new StPicoDstMaker(static_cast<StPicoDstMaker::PicoIoMode>(StPicoDstMaker::IoRead), inputFile, "picoDstMaker");

  StPicoPiPiMaker* PicoK0sAnaMaker = new StPicoPiPiMaker("picoK0sAnaMaker", picoDstMaker, outputFile);
  PicoK0sAnaMaker->setHFBaseCuts(hfCuts);

//  StPicoD0AnaMaker* PicoD0AnaMaker = new StPicoD0AnaMaker("picoD0AnaMaker", picoDstMaker, outputFile, sInputListHF);
//  PicoD0AnaMaker->setTreeName(treeName);
//  PicoD0AnaMaker->setDecayMode(StPicoHFEvent::kTwoParticleDecay);
//  PicoD0AnaMaker->setHFBaseCuts(hfCuts);

//  StPicoMixedEventMaker* picoMixedEventMaker = new StPicoMixedEventMaker("picoMixedEventMaker", picoDstMaker, hfCuts, outputFile, inputFile);
//  picoMixedEventMaker->setBufferSize(3);

  clock_t start = clock(); // getting starting time
  chain->Init();
  
  int total = picoDstMaker->chain()->GetEntries();
  cout << " Total entries = " << total << endl;
  if(nEvents>total) nEvents = total;

  for (Int_t i=0; i<nEvents; i++) {
//  for (Int_t i=0; i<2000; i++) {
    if(i%10==0)       cout << "Working on eventNumber " << i << endl;

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

