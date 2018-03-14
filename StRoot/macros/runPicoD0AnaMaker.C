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

    StPicoDstMaker* picoDstMaker = new StPicoDstMaker(static_cast<StPicoDstMaker::PicoIoMode>(StPicoDstMaker::IoRead), inputFile, "picoDstMaker");
    cout<<"ok, picoDstMaker created"<<endl;
    StPicoD0AnaMaker* PicoD0AnaMaker = new StPicoD0AnaMaker("picoD0AnaMaker", picoDstMaker, outputFile, sInputListHF);
    PicoD0AnaMaker->setTreeName(treeName);
    PicoD0AnaMaker->setDecayMode(StPicoHFEvent::kTwoParticleDecay);

    StHFCuts* hfCuts = new StHFCuts("hfBaseCuts");
    PicoD0AnaMaker->setHFBaseCuts(hfCuts);

    hfCuts->setBadRunListFileName(badRunListFileName);

    hfCuts->addTriggerId(530003); //VPD-5

    hfCuts->setCutPrimaryDCAtoVtxMax(1.5);
    hfCuts->setCutVzMax(6);
    hfCuts->setCutVzVpdVzMax(1000.);
    hfCuts->setCutNHitsFitMin(15);
    hfCuts->setCutRequireHFT(true);

//    WIDE
    hfCuts->setCutDcaMin(0.003,StHFCuts::kPion); // OK
    hfCuts->setCutDcaMin(0.003,StHFCuts::kKaon); // OK

    float dcaDaughtersMax = 0.05;  // maximum toto ide
    float decayLengthMin  = 0.009; // minimum
    float decayLengthMax  = 999;  //std::numeric_limits<float>::max(); toto ide (cutuje)
    float cosThetaMin     = 0.5;   // minimum
    float minMass         = 0.4;
    float maxMass         = 2.4;

    //    WIDER than pt12 tmva
//    hfCuts->setCutDcaMin(0.007,StHFCuts::kPion); // OK
//   hfCuts->setCutDcaMin(0.007,StHFCuts::kKaon); // OK
    hfCuts->setCutTPCNSigmaPion(3.0);
    hfCuts->setCutTPCNSigmaKaon(2.0);
    hfCuts->setCutTOFDeltaOneOverBetaKaon(0.03);
//    hfCuts->setCutTOFDeltaOneOverBetaPion(999);
    hfCuts->setCutPtMin(0.2);
//    float dcaDaughtersMax = 0.012;  // maximum toto ide
//    float decayLengthMin  = 0.005; // minimum
//    float decayLengthMax  = 4; //std::numeric_limits<float>::max(); toto ide (cutuje)
//    float cosThetaMin     = 0.5;   // minimum

//    tmva ready
//    hfCuts->setCutDcaMin(0.000,StHFCuts::kPion); // OK
//    hfCuts->setCutDcaMin(0.000,StHFCuts::kKaon); // OK
//
//    float dcaDaughtersMax = 1;  // maximum toto ide
//    float decayLengthMin  = 0.001; // minimum
//    float decayLengthMax  = 4; //std::numeric_limits<float>::max(); toto ide (cutuje)
//    float cosThetaMin     = -2;   // minimum
//    float minMass         = 1.761;
//    float maxMass         = 1.971;


//    float dcaToPvMax = 0.003649;
//    hfCuts -> setCutSecondaryPairDcaToPvMax(dcaToPvMax); //getting read for isgooSVpair.

//    float dcaDaughtersMax = 2;  // maximum
//    float decayLengthMin  = 0.000; // minimum
//    float decayLengthMax  = 3; //std::numeric_limits<float>::max();
//    float cosThetaMin     = 0.;   // minimum
//    float minMass         = 0.4;
//    float maxMass         = 2.4;
//
//    hfCuts->setCutSecondaryPair(dcaDaughtersMax, decayLengthMin, decayLengthMax, cosThetaMin, minMass, maxMass); //ok
    
    //Single track pt


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
    
    
