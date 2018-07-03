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
void runPicoD0AnaMaker(
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
    hfCuts->setCutRequireHFT(true);

    hfCuts->setCutTPCNSigmaPion(3.0);
    hfCuts->setCutTPCNSigmaKaon(2.0);
    hfCuts->setCutTOFDeltaOneOverBetaKaon(0.03);
    hfCuts->setCutTOFDeltaOneOverBetaPion(0.03);
    hfCuts->setCutPtMin(0.15);

    hfCuts->setCutDcaMin(0.00,StHFCuts::kPion);
    hfCuts->setCutDcaMin(0.00,StHFCuts::kKaon);

    float dcaDaughtersMax = 0.04;  // maximum toto ide
    float decayLengthMin  = 0.00; // minimum
    float decayLengthMax  = 999;  //std::numeric_limits<float>::max(); toto ide (cutuje)
    float cosThetaMin     = 0.;   // minimum
    float minMass         = 0.4;
    float maxMass         = 2.4;

//    hfCuts->setCutDcaMin(0.004,StHFCuts::kPion);
//    hfCuts->setCutDcaMin(0.004,StHFCuts::kKaon);
//
//    float dcaDaughtersMax = 0.04;  // maximum toto ide
//    float decayLengthMin  = 0.009; // minimum
//    float decayLengthMax  = 999;  //std::numeric_limits<float>::max(); toto ide (cutuje)
//    float cosThetaMin     = 0.5;   // minimum
//    float minMass         = 0.4;
//    float maxMass         = 2.4;

    hfCuts->setCutSecondaryPair(dcaDaughtersMax, decayLengthMin, decayLengthMax, cosThetaMin, minMass, maxMass);

    StPicoDstMaker* picoDstMaker = new StPicoDstMaker(static_cast<StPicoDstMaker::PicoIoMode>(StPicoDstMaker::IoRead), inputFile, "picoDstMaker");

//    StPicoD0AnaMaker* PicoD0AnaMaker = new StPicoD0AnaMaker("picoD0AnaMaker", picoDstMaker, outputFile, sInputListHF);
//    PicoD0AnaMaker->setTreeName(treeName);
//    PicoD0AnaMaker->setDecayMode(StPicoHFEvent::kTwoParticleDecay);
//    PicoD0AnaMaker->setHFBaseCuts(hfCuts);

    StPicoMixedEventMaker* picoMixedEventMaker = new StPicoMixedEventMaker("picoMixedEventMaker", picoDstMaker, hfCuts, outputFile, inputFile);
    picoMixedEventMaker->setBufferSize(7);
    
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
