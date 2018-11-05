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
#include "StPicoKKMaker/StPicoKKMaker.h"

using namespace std;

#else
class StChain;
#endif
class StPicoDstMaker;
class StPicoKKMaker;
class StMaker;
StChain *chain;
void runPicoPhiAnaMaker(
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

    hfCuts->setCutVzMax(6.);
    hfCuts->setCutVzVpdVzMax(3.);
    hfCuts->addTriggerId(530003); //VPD-5

    hfCuts->setCutNHitsFitMin(15); //default is 20
    hfCuts->setCutRequireHFT(true);
    hfCuts->setHybridTof(false);
    hfCuts->setCutDcaMin(0.00,StHFCuts::kKaon);

    // kaonPion pair cuts
    float dcaDaughtersMax = 0.7;  // maximum
    float decayLengthMin  = 0.000; // minimum
    float decayLengthMax  = 25; //std::numeric_limits<float>::max();
    float cosThetaMin     = 0.2;   // minimum
    float minMass         = 1;
    float maxMass         = 1.05;
    float pairDcaMax      = 0.5;

    hfCuts->setCutPtRange(0.15,50.0,StHFCuts::kKaon);
    hfCuts->setCutTPCNSigmaKaon(999);
    hfCuts->setCutTOFDeltaOneOverBeta(999, StHFCuts::kKaon);
    hfCuts->setCutPtotRangeHybridTOF(0.15,50.0,StHFCuts::kKaon);

    hfCuts->setCutSecondaryPair(dcaDaughtersMax, decayLengthMin, decayLengthMax, cosThetaMin, minMass, maxMass, pairDcaMax);

    StPicoDstMaker* picoDstMaker = new StPicoDstMaker(static_cast<StPicoDstMaker::PicoIoMode>(StPicoDstMaker::IoRead), inputFile, "picoDstMaker");

    StPicoKKMaker* PicoPhiAnaMaker = new StPicoKKMaker("picoPhiAnaMaker", picoDstMaker, outputFile, sInputListHF);
    PicoPhiAnaMaker->setTreeName(treeName);
    PicoPhiAnaMaker->setDecayMode(StPicoHFEvent::kTwoParticleDecay);
    PicoPhiAnaMaker->setHFBaseCuts(hfCuts);

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
