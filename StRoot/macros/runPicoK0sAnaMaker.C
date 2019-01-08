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
#include "StPicoPiPiMaker/StPicoPiPiMaker.h"

using namespace std;

#else
class StChain;
#endif
class StPicoDstMaker;
class StMaker;
class StPicoPiPiMaker;

StChain *chain;
void runPicoK0sAnaMaker(
    const char*  inputFile,
    const Char_t *outputFile,  
    const Char_t *badRunListFileName, const Char_t *treeName,
    const Char_t *productionBasePath) {
    string SL_version = "SL17d";
    string env_SL = getenv ("STAR");
    if (env_SL.find(SL_version)==string::npos) {
        cout<<"Environment Star Library does not match the requested library in run**.C. Exiting..."<<endl;
        exit(1);
    }

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

    StHFCuts* hfCuts = new StHFCuts("hfBaseCuts");

    hfCuts->setBadRunListFileName(badRunListFileName);
    hfCuts->addTriggerId(530003); //VPD-5
    hfCuts->setCutPrimaryDCAtoVtxMax(10);
    hfCuts->setCutVzMax(6);
    hfCuts->setCutVzVpdVzMax(3.);
    hfCuts->setCutNHitsFitMin(15);
    hfCuts->setCutRequireHFT(false);
    hfCuts->setHybridTof(true);

    hfCuts->setCutTPCNSigmaPion(6);
    hfCuts->setCutTOFDeltaOneOverBetaPion(0.08); //0.013*5
    hfCuts->setCutPtMin(0.15);

    hfCuts->setCutDcaMin(0.08,StHFCuts::kPion);

    float dcaDaughtersMax = 0.6;  // maximum toto ide
    float decayLengthMin  = 0.5; // minimum
    float decayLengthMax  = 6;  //
    float cosThetaMin     = 0.7;   // minimum
    float minMass         = 0.42;
    float maxMass         = 0.58;
    float pairDcaMax      = 0.9;

    hfCuts->setCutSecondaryPair(dcaDaughtersMax, decayLengthMin, decayLengthMax, cosThetaMin, minMass, maxMass, pairDcaMax);

    StPicoDstMaker* picoDstMaker = new StPicoDstMaker(static_cast<StPicoDstMaker::PicoIoMode>(StPicoDstMaker::IoRead), inputFile, "picoDstMaker");

    StPicoPiPiMaker* PicoK0sAnaMaker = new StPicoPiPiMaker("picoK0sAnaMaker", picoDstMaker, outputFile, sInputListHF);
    PicoK0sAnaMaker->setTreeName(treeName);
    PicoK0sAnaMaker->setDecayMode(StPicoHFEvent::kTwoParticleDecay);
    PicoK0sAnaMaker->setHFBaseCuts(hfCuts);

    clock_t start = clock(); // getting starting time
    chain->Init();
    Int_t nEvents = picoDstMaker->chain()->GetEntries();
    cout << " Total entries = " << nEvents << endl;

    for (Int_t i=0; i<nEvents; ++i) {
//        if(i%10==0)       cout << "Working on eventNumber " << i << endl;
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
