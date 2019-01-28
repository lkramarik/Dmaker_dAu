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

void runPicoPhiAnaMaker(
    const char*  inputFile,
    const Char_t *outputFile,  
    const Char_t *badRunListFileName) {
    string SL_version = "SL18f";
    string env_SL = getenv ("STAR");
    if (env_SL.find(SL_version)==string::npos) {
        cout<<"Environment Star Library does not match the requested library in run**.C. Exiting..."<<endl;
        exit(1);
    }

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

    hfCuts->setCutNHitsFitMin(15); //default is 20
    hfCuts->setCutRequireHFT(false);
    hfCuts->setHybridTof(true);
    hfCuts->setCutDcaMin(0.00,StHFCuts::kKaon);

    // kaonPion pair cuts
    float dcaDaughtersMax = 1;  // maximum
    float decayLengthMin  = 0.000; // minimum
    float decayLengthMax  = 25; //std::numeric_limits<float>::max();
    float cosThetaMin     = 0.85;   // minimum
    float minMass         = 1;
    float maxMass         = 1.05;
    float pairDcaMax      = 1;

    hfCuts->setCutPtRange(0.15,50.0,StHFCuts::kKaon);
    hfCuts->setCutTPCNSigmaKaon(6);
    hfCuts->setCutTOFDeltaOneOverBeta(0.08, StHFCuts::kKaon);
    hfCuts->setCutPtotRangeHybridTOF(0.15,50.0,StHFCuts::kKaon);

    hfCuts->setCutSecondaryPair(dcaDaughtersMax, decayLengthMin, decayLengthMax, cosThetaMin, minMass, maxMass, pairDcaMax);

    StPicoDstMaker* picoDstMaker = new StPicoDstMaker(static_cast<StPicoDstMaker::PicoIoMode>(StPicoDstMaker::IoRead), inputFile, "picoDstMaker");

    StPicoKKMaker* PicoPhiAnaMaker = new StPicoKKMaker("picoPhiAnaMaker", picoDstMaker, outputFile);
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

    chain->Finish();
    double duration = (double) (clock() - start) / (double) CLOCKS_PER_SEC;
    cout << "****************************************** " << endl;
    cout << "Work done, total number of events  " << nEvents << endl;
    cout << "Time needed " << duration << " s" << endl;
    delete chain;
}
