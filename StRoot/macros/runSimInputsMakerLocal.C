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
#include "StPicoSimInputsMaker/StPicoSimInputsMaker.h"

using namespace std;

void runSimInputsMakerLocal(
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
  hfCuts->addTriggerId(530003); //VPD-5
  hfCuts->setCutVzMax(6);
  hfCuts->setCutVzVpdVzMax(6.);

  hfCuts->setCutNHitsFitMin(15);
  hfCuts->setCutRequireHFT(false); //we want to study HFT ratio, thus need non-HFT tracks
  hfCuts->setHybridTof(true);
  hfCuts->setCutPrimaryDCAtoVtxMax(4);

  hfCuts->setCutTPCNSigmaPion(3.0);
  hfCuts->setCutTPCNSigmaKaon(2.0);

  hfCuts->setCutTOFDeltaOneOverBeta(0.03, StHFCuts::kPion); // v podstate 5 sigma; nastavene = f * (sigmaTOF), sigma TOF je 0.013
  hfCuts->setCutTOFDeltaOneOverBeta(0.03, StHFCuts::kKaon); // v podstate 5 sigma; nastavene = f * (sigmaTOF), sigma TOF je 0.013

  hfCuts->setCutPtMin(0.15);

  StPicoDstMaker* picoDstMaker = new StPicoDstMaker(static_cast<StPicoDstMaker::PicoIoMode>(StPicoDstMaker::IoRead), inputFile, "picoDstMaker");
  StPicoSimInputsMaker* picoSimInputs = new StPicoSimInputsMaker("picoSimInputs", picoDstMaker, outputFile);
  picoSimInputs->setHFBaseCuts(hfCuts);

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

