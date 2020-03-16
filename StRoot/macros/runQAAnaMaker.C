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

void runQAAnaMaker(
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
    hfCuts->addTriggerId(530003); //VPD-5
    hfCuts->setCheckHotSpot(true);
    //    hfCuts->addTriggerId(2); //ZDC-VPD-5
//    hfCuts->addTriggerId(3); //VPD-5
//    hfCuts->addTriggerId(6);  //highMult-VPD-5
//    hfCuts->addTriggerId(7);  //highMult2-VPD-5
//    hfCuts->addTriggerId(15); //BHT1-VPD-10
//    hfCuts->addTriggerId(16); //BHT2-VPD-30
//    hfCuts->addTriggerId(17); //BHT3
//    hfCuts->addTriggerId(530002); //ZCD-VPD-5
//    hfCuts->addTriggerId(530101); //highMult-VPD-5
//    hfCuts->addTriggerId(530102); //highMult2-VPD-5
//    hfCuts->addTriggerId(530201); //BHT1-VPD-10
//    hfCuts->addTriggerId(530202); //BHT2-VPD-30
//    hfCuts->addTriggerId(530213); //BHT3

    hfCuts->setCutVzMax(6);
    hfCuts->setCutVzVpdVzMax(6.);

    hfCuts->setCutTOFDeltaOneOverBeta(0.03, StHFCuts::kPion); // v podstate 5 sigma; nastavene = f * (sigmaTOF), sigma TOF je 0.013
    hfCuts->setCutTOFDeltaOneOverBeta(0.03, StHFCuts::kKaon); // v podstate 5 sigma; nastavene = f * (sigmaTOF), sigma TOF je 0.013

    StPicoDstMaker* picoDstMaker = new StPicoDstMaker(static_cast<StPicoDstMaker::PicoIoMode>(StPicoDstMaker::IoRead), inputFile, "picoDstMaker");
    StPicoQAMaker* PicoQAMaker = new StPicoQAMaker("picoQAMaker", picoDstMaker, outputFile);
    PicoQAMaker->setHFBaseCuts(hfCuts);

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
