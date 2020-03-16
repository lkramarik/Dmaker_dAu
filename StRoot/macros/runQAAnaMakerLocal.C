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
#include "StPicoQAMaker/StPicoQAMaker.h"

using namespace std;

class StChain;
class StPicoDstMaker;
class StPicoQAMaker;
class StMaker;

void runQAAnaMakerLocal(
			const Char_t *inputFile="./picoLists/runs_local_test.list",
			const Char_t *outputFile="outputLocal",
			const Char_t *badRunListFileName = "./picoLists/picoList_bad.list") {
  string SL_version = "SL18f";
  string env_SL = getenv ("STAR");
  if (env_SL.find(SL_version)==string::npos) {
      cout<<"Environment Star Library does not match the requested library. Exiting..."<<endl;
      exit(1);
  }
  Int_t nEvents = 2000000;
  StChain *chain = new StChain();
  TString sInputFile(inputFile);

  if (!sInputFile.Contains(".list") && !sInputFile.Contains("picoDst.root")) {
    cout << "No input list or picoDst root file provided! Exiting..." << endl;
    exit(1);
  }
  cout<<"event stuff set"<<endl;


  StHFCuts* hfCuts = new StHFCuts("hfBaseCuts");

  hfCuts->setBadRunListFileName(badRunListFileName);

  hfCuts->setCutVzMax(6.);
  hfCuts->setCutVzVpdVzMax(6.);
  hfCuts->addTriggerId(530003); //VPD-5
  hfCuts->setCheckHotSpot(true);

  hfCuts->setCutNHitsFitMin(15); //default is 20
  hfCuts->setCutRequireHFT(true);
  hfCuts->setHybridTof(false);
  //LK hfCuts->setCutDcaMin(0.009,StHFCuts::kPion);
  //LK  hfCuts->setCutDcaMin(0.007,StHFCuts::kKaon);

   // kaonPion pair cuts
  float dcaDaughtersMax = 0.2;  // maximum
  float decayLengthMin  = 0.000; // minimum
  float decayLengthMax  = 999999; //std::numeric_limits<float>::max();
  float cosThetaMin     = -20.;   // minimum
  float minMass         = 0.6;
  float maxMass         = 2.6;
  float pairDcaMax      = 99.9;

  //Single track pt
  hfCuts->setCutPtRange(0.15,50.0,StHFCuts::kPion); //0.2 , 50.0
  hfCuts->setCutPtRange(0.15,50.0,StHFCuts::kKaon); //0.2, 50.0
  //TPC setters
  hfCuts->setCutTPCNSigmaPion(10); //3
  hfCuts->setCutTPCNSigmaKaon(10); //3
  //TOF setters, need to set pt range as well
  hfCuts->setCutTOFDeltaOneOverBeta(0.1, StHFCuts::kKaon); // v podstate 5 sigma; nastavene = f * (sigmaTOF), sigma TOF je 0.013
  hfCuts->setCutPtotRangeHybridTOF(0.2,50.0,StHFCuts::kKaon);
  hfCuts->setCutTOFDeltaOneOverBeta(0.1, StHFCuts::kPion); // v podstate 6 sigma
  hfCuts->setCutPtotRangeHybridTOF(0.2,50.0,StHFCuts::kPion);

  StPicoDstMaker* picoDstMaker = new StPicoDstMaker(static_cast<StPicoDstMaker::PicoIoMode>(StPicoDstMaker::IoRead), inputFile, "picoDstMaker");
  StPicoQAMaker* PicoQAAnaMaker = new StPicoQAMaker("picoQAAnaMaker", picoDstMaker, outputFile);
  PicoQAAnaMaker->setHFBaseCuts(hfCuts);

  clock_t start = clock(); // getting starting time
  chain->Init();
  
  int total = picoDstMaker->chain()->GetEntries();
  cout << " Total entries = " << total << endl;
  if(nEvents>total) nEvents = total;

  for (Int_t i=0; i<nEvents; i++) {
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
  cout << "total number of events  " << nEvents << endl;
  cout << "****************************************** " << endl;
  cout << "Time needed " << duration << " s" << endl;
  cout << "****************************************** " << endl;
  
  delete chain;

}

