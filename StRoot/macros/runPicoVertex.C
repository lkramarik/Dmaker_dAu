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
#include "StPicoKFVertexTools/StPicoKFVertexTools.h"

using namespace std;

void runPicoVertex(
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
    hfCuts->setCutVzMax(6);
    hfCuts->setCutVzVpdVzMax(6.);

    hfCuts->setCutPrimaryDCAtoVtxMax(10);
    hfCuts->setCutNHitsFitMin(15);
    hfCuts->setCutRequireHFT(true);
    hfCuts->setHybridTof(true);

    hfCuts->setCutTPCNSigmaPion(3.0);
    hfCuts->setCutTPCNSigmaKaon(2.0);
    hfCuts->setCutTOFDeltaOneOverBetaKaon(0.03);
    hfCuts->setCutTOFDeltaOneOverBetaPion(0.03);
    hfCuts->setCutPtMin(0.15);

    hfCuts->setCutDcaMin(0.001,StHFCuts::kPion);
    hfCuts->setCutDcaMin(0.001,StHFCuts::kKaon);

//   setCutSecondaryPairPtBin(      ptmin,  ptmax,  dcaDaughtersMax,  decayLengthMin,  cosThetaMin,  pairDcaMax, pionDca, kaonDca);
    hfCuts->setCutSecondaryPairPtBin(1,      2,              0.016,          0.012,         0.5,      0.005,    0.009, 0.007);
    hfCuts->setCutSecondaryPairPtBin(2,      3,              0.016,          0.003,         0.5,      0.0065,   0.009, 0.01);
    hfCuts->setCutSecondaryPairPtBin(3,      5,              0.018,          0.009,         0.5,      0.0064,   0.0064, 0.0076);

    StPicoDstMaker* picoDstMaker = new StPicoDstMaker(static_cast<StPicoDstMaker::PicoIoMode>(StPicoDstMaker::IoRead), inputFile, "picoDstMaker");
    StPicoKFVertexTools* PicoVertex = new StPicoKFVertexTools("PicoVertex", picoDstMaker, outputFile);
    PicoVertex->setHFBaseCuts(hfCuts);

    clock_t start = clock(); // getting starting time
    chain->Init();
    Int_t nEvents = picoDstMaker->chain()->GetEntries();
    cout << " Total entries = " << nEvents << endl;

    for (Int_t i=0; i<nEvents; ++i) {
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
