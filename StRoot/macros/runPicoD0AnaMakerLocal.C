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
#include "StPicoMixedEventMaker/StPicoMixedEventMaker.h"
#include <iostream>
#include <ctime>
#include <cstdio>
#include "StPicoD0AnaMaker/StPicoD0AnaMaker.h"
#include "StPicoD0V2AnaMaker/StPicoD0V2AnaMaker.h"
#include "StPicoQAMaker/StPicoQAMaker.h"
#include "StPicoKFVertexTools/StPicoKFVertexTools.h"

using namespace std;

void runPicoD0AnaMakerLocal(
			const Char_t *inputFile="./picoLists/runs_local_test.list",
			const Char_t *outputFile="outputLocal",
			const Char_t *badRunListFileName = "./picoLists/picoList_bad.list") {
  string SL_version = "SL18f";
  string env_SL = getenv ("STAR");
  if (env_SL.find(SL_version)==string::npos) {
      cout<<"Environment Star Library does not match the requested library. Exiting..."<<endl;
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

  // -- File name of bad run list
   hfCuts->setBadRunListFileName(badRunListFileName); 

  hfCuts->setCutVzMax(6.);
  hfCuts->setCutVzVpdVzMax(6.);
  hfCuts->addTriggerId(530003); //VPD-5

  hfCuts->setCutNHitsFitMin(15); //default is 20
  hfCuts->setCutRequireHFT(true);
  hfCuts->setHybridTof(true);
  hfCuts->setCheckHotSpot(false);


    //LK hfCuts->setCutDcaMin(0.009,StHFCuts::kPion); //federic 1aug2016
  //LK  hfCuts->setCutDcaMin(0.007,StHFCuts::kKaon); //federic 3aug2016
  //hfCuts->setCutNHitsFitnHitsMax(0.52);  kvapil

   // kaonPion pair cuts
  float dcaDaughtersMax = 0.2;  // maximum
  float decayLengthMin  = 0.000; // minimum
  float decayLengthMax  = 999999; //std::numeric_limits<float>::max();
  float cosThetaMin     = -20.;   // minimum
  float minMass         = 0.6;
  float maxMass         = 2.6;
  float pairDcaMax      = 99.9;

  hfCuts->setCutSecondaryPair(dcaDaughtersMax, decayLengthMin, decayLengthMax, cosThetaMin, minMass, maxMass, pairDcaMax);
 
  //Single track pt
  hfCuts->setCutPtRange(0.15,50.0,StHFCuts::kPion); //0.2 , 50.0
  hfCuts->setCutPtRange(0.15,50.0,StHFCuts::kKaon); //0.2, 50.0
  //TPC setters
  hfCuts->setCutTPCNSigmaPion(3.0); //3
  hfCuts->setCutTPCNSigmaKaon(3.0); //2
  //TOF setters, need to set pt range as well
  hfCuts->setCutTOFDeltaOneOverBeta(0.06, StHFCuts::kKaon); // v podstate 5 sigma; nastavene = f * (sigmaTOF), sigma TOF je 0.013
  hfCuts->setCutPtotRangeHybridTOF(0.2,50.0,StHFCuts::kKaon);
  hfCuts->setCutTOFDeltaOneOverBeta(0.06, StHFCuts::kPion); // v podstate 6 sigma
  hfCuts->setCutPtotRangeHybridTOF(0.2,50.0,StHFCuts::kPion);

  //                               ptmin, ptmax, dcaDaughtersMax, decayLengthMin,  cosThetaMin, pairDcaMax, pionDca, kaonDca
  hfCuts->setCutSecondaryPairPtBin(1,      2,              0.007,          0.012,         0.5,      0.005,    0.009, 0.007);
  hfCuts->setCutSecondaryPairPtBin(2,      3,              0.016,          0.003,         0.5,      0.0065,   0.009, 0.01);
  hfCuts->setCutSecondaryPairPtBin(3,      5,              0.015,          0.009,         0.6,      0.0064,   0.0064, 0.0076);


  StPicoDstMaker* picoDstMaker = new StPicoDstMaker(static_cast<StPicoDstMaker::PicoIoMode>(StPicoDstMaker::IoRead), inputFile, "picoDstMaker");
  StPicoD0AnaMaker* PicoD0AnaMaker = new StPicoD0AnaMaker("picoD0AnaMaker", picoDstMaker, outputFile);
    PicoD0AnaMaker->workWithRefit(false);
    PicoD0AnaMaker->setHFBaseCuts(hfCuts);

//  StPicoMixedEventMaker* picoMixedEventMaker = new StPicoMixedEventMaker("picoMixedEventMaker", picoDstMaker, hfCuts, outputFile);
//  picoMixedEventMaker->setBufferSize(3);

//  StPicoD0V2AnaMaker* PicoD0V2AnaMaker = new StPicoD0V2AnaMaker("picoD0V2AnaMaker", picoDstMaker, outputFile);
//  PicoD0V2AnaMaker->setHFBaseCuts(hfCuts);


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

