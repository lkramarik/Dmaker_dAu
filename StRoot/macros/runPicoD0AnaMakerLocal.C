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
#include "macros/loadSharedHFLibraries.C"
#include <iostream>
#include <ctime>
#include <cstdio>
#include "StPicoD0AnaMaker/StPicoD0AnaMaker.h"

using namespace std;

#else
class StChain;
#endif
StChain *chain;
void runPicoD0AnaMakerLocal(
			const Char_t *inputFile="/gpfs01/star/pwg/lkramarik/Dmaker_dAu/picoLists/runs_local_test.list",	
			const Char_t *outputFile="outputBaseName",  
            const unsigned int makerMode = 0 ,
			const Char_t *badRunListFileName = "/gpfs01/star/pwg/lkramarik/Dmaker_dAu/picoLists/picoList_bad.list",
            const Char_t *treeName = "picoHFtree",
			const Char_t *productionBasePath = "/gpfs01/star/pwg/lkramarik/Dmaker_dAu/") {
  string SL_version = "SL17d";
  string env_SL = getenv ("STAR");
  if (env_SL.find(SL_version)==string::npos) {
      cout<<"Environment Star Library does not match the requested library in runPicoHFMyAnaMaker.C. Exiting..."<<endl;
      exit(1);
  }
  
  Int_t nEvents = 1000000;

#ifdef __CINT__
  gROOT->LoadMacro("loadSharedHFLibraries.C");
  loadSharedHFLibraries();
#endif

  chain = new StChain();

  TString sInputFile(inputFile);
  TString sInputListHF("");
  TString sProductionBasePath(productionBasePath);
  TString sTreeName(treeName);

  if (!sInputFile.Contains(".list") && !sInputFile.Contains("picoDst.root")) {
    cout << "No input list or picoDst root file provided! Exiting..." << endl;
    exit(1);
  }

  StPicoDstMaker* picoDstMaker = new StPicoDstMaker(static_cast<StPicoDstMaker::PicoIoMode>(StPicoDstMaker::IoRead), inputFile, "picoDstMaker");
  cout<<"ok, picoDstMaker created"<<endl;
  StPicoD0AnaMaker* PicoD0AnaMaker = new StPicoD0AnaMaker("picoD0AnaMaker", picoDstMaker, outputFile, sInputListHF);
  PicoD0AnaMaker->setTreeName(treeName);
  PicoD0AnaMaker->setDecayMode(StPicoHFEvent::kTwoParticleDecay);

  StHFCuts* hfCuts = new StHFCuts("hfBaseCuts");
  PicoD0AnaMaker->setHFBaseCuts(hfCuts);
  cout<<"event stuff set"<<endl;
  // ---------------------------------------------------
  // -- Set Base cuts for HF analysis

  // -- File name of bad run list
   hfCuts->setBadRunListFileName(badRunListFileName); 

  // -- ADD USER CUTS HERE ----------------------------

  hfCuts->setCutVzMax(6.);
  hfCuts->setCutVzVpdVzMax(3.);
  
  hfCuts->addTriggerId(3); //VPD-5
  hfCuts->addTriggerId(6);  //highMult-VPD-5
  hfCuts->addTriggerId(7);  //highMult2-VPD-5
  hfCuts->addTriggerId(15); //BHT1-VPD-10
  hfCuts->addTriggerId(16); //BHT2-VPD-30
  hfCuts->addTriggerId(530003); //VPD-5
  hfCuts->addTriggerId(530101); //highMult-VPD-5
  hfCuts->addTriggerId(530102); //highMult2-VPD-5
  hfCuts->addTriggerId(530201); //BHT1-VPD-10
  hfCuts->addTriggerId(530202); //BHT2-VPD-30 
  hfCuts->addTriggerId(530213); //BHT3
  

  hfCuts->setCutNHitsFitMin(15); //default is 20
  hfCuts->setCutRequireHFT(true);

  //LK hfCuts->setCutDcaMin(0.009,StHFCuts::kPion); //federic 1aug2016
  //LK  hfCuts->setCutDcaMin(0.007,StHFCuts::kKaon); //federic 3aug2016
  //hfCuts->setCutNHitsFitnHitsMax(0.52);  kvapil

  // -- Channel0
  PicoD0AnaMaker->setDecayMode(StPicoHFEvent::kTwoParticleDecay);

  // -- ADD USER CUTS HERE ----------------------------
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
  // set refmultCorr
  //StRefMultCorr* grefmultCorrUtil = CentralityMaker::instance()->getgRefMultCorr_P16id();
  //PicoD0AnaMaker->setRefMutCorr(grefmultCorrUtil);

  clock_t start = clock(); // getting starting time
  chain->Init();
  
  int total = picoDstMaker->chain()->GetEntries();
  cout << " Total entries = " << total << endl;
  if(nEvents>total) nEvents = total;

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
  cout << "total number of events  " << nEvents << endl;
  cout << "****************************************** " << endl;
  cout << "Time needed " << duration << " s" << endl;
  cout << "****************************************** " << endl;
  
  delete chain;

}

