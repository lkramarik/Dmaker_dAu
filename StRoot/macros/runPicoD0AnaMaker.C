/*   root -l -b -q StRoot/macros/loadSharedHFLibraries.C StRoot/macros/runPicoHFMyAnaMaker.C++
 *  - Different modes to use the  class
 *    - StPicoHFMaker::kAnalyze - don't write candidate trees, just fill histograms
 *        inputFile : fileList of PicoDst files or single picoDst file
 *        outputFile: baseName for outfile 
 *    - StPicoHFMaker::kWrite   - write candidate trees
 *        inputFile : path to single picoDist file
 *        outputFile: baseName for outfile 
 *    - StPicoHFMaker::kRead    - read candidate trees and fill histograms
 *        inputFile : fileList of PicoDst files
 *        outputFile: baseName for outfile 
 */

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
//#include "StPicoHFMyAnaMaker/StPicoHFMyAnaMaker.h"
#include "macros/loadSharedHFLibraries.C"
#include <iostream>
#include <ctime>
#include <cstdio>
#include "StPicoD0AnaMaker/StPicoD0AnaMaker.h" //kvapil
//#include "StRefMultCorr/StRefMultCorr.h"
//#include "StRefMultCorr/CentralityMaker.h"

using namespace std;

#else
class StChain;
#endif
StChain *chain;
void runPicoD0AnaMaker(
    //const Char_t *inputFile="/gpfs01/star/pwg/lkramarik/pico_pAu/produced/st_physics_16127031_raw_4000024.picoDst.root",	
    const char*  inputFile, 			 
    const Char_t *outputFile,  
    const unsigned int makerMode,
    const Char_t *badRunListFileName, const Char_t *treeName,
    const Char_t *productionBasePath) {
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
    
    // ========================================================================================
    //makerMode    = StPicoHFMaker::kAnalyze;
    // ========================================================================================
    cout << "Maker Mode    " << makerMode << endl;
    cout << "TreeName      " << treeName << endl; 

    TString sInputFile(inputFile);
    TString sInputListHF("");  
    TString sProductionBasePath(productionBasePath);
    TString sTreeName(treeName);
    
    if (makerMode == StPicoHFMaker::kAnalyze) {
        if (!sInputFile.Contains(".list") && !sInputFile.Contains("picoDst.root")) {
            cout << "No input list or picoDst root file provided! Exiting..." << endl;
            exit(1);
        }
    }
    else if (makerMode == StPicoHFMaker::kWrite) {
        if (!sInputFile.Contains("picoDst.root")) {
            cout << "No input picoDst root file provided! Exiting..." << endl;
            exit(1);
        }
    }
    else if (makerMode == StPicoHFMaker::kRead) {
        if (!sInputFile.Contains(".list")) {
            cout << "No input list provided! Exiting..." << endl;
            exit(1);
        }
        
        // -- prepare filelist for picoDst from trees
        sInputListHF = sInputFile;
        sInputFile = "tmpPico.list";
        
        TString command = "sed 's|" + sTreeName + ".root|picoDst.root|g' " + sInputListHF + " > " + sInputFile;
        cout << "COMMAND : " << command << endl; 
        gSystem->Exec(command.Data());
        
        command = "sed -i 's|^.*" + sTreeName + "|" + sProductionBasePath + "|g' " + sInputFile; // + " > " + sInputFile;
        cout << "COMMAND : " << command << endl; 
        gSystem->Exec(command.Data());
    }
    else {
        cout << "Unknown makerMode! Exiting..." << endl;
        exit(1);
    }
    
    StPicoDstMaker* picoDstMaker = new StPicoDstMaker(static_cast<StPicoDstMaker::PicoIoMode>(StPicoDstMaker::IoRead), inputFile, "picoDstMaker");
    cout<<"ok, picoDstMaker created"<<endl;
    StPicoD0AnaMaker* PicoD0AnaMaker = new StPicoD0AnaMaker("picoD0AnaMaker", picoDstMaker, outputFile, sInputListHF);
    PicoD0AnaMaker->setMakerMode(makerMode);
    PicoD0AnaMaker->setTreeName(treeName);
    PicoD0AnaMaker->setDecayMode(StPicoHFEvent::kTwoParticleDecay);

    StHFCuts* hfCuts = new StHFCuts("hfBaseCuts");
    PicoD0AnaMaker->setHFBaseCuts(hfCuts);

    cout<<"event stuff set"<<endl;

    // -- File name of bad run list
    hfCuts->setBadRunListFileName(badRunListFileName); 

    hfCuts->addTriggerId(2); //ZDC-VPD-5 
    hfCuts->addTriggerId(3); //VPD-5
    hfCuts->addTriggerId(6);  //highMult-VPD-5
    hfCuts->addTriggerId(7);  //highMult2-VPD-5
    hfCuts->addTriggerId(15); //BHT1-VPD-10
    hfCuts->addTriggerId(16); //BHT2-VPD-30
    hfCuts->addTriggerId(17); //BHT3
    hfCuts->addTriggerId(530002); //ZCD-VPD-5
    hfCuts->addTriggerId(530003); //VPD-5
    hfCuts->addTriggerId(530101); //highMult-VPD-5
    hfCuts->addTriggerId(530102); //highMult2-VPD-5
    hfCuts->addTriggerId(530201); //BHT1-VPD-10
    hfCuts->addTriggerId(530202); //BHT2-VPD-30 
    hfCuts->addTriggerId(530213); //BHT3
    
//    hfCuts->setCutDcaMin(0.00,StHFCuts::kPion); // OK
//    hfCuts->setCutDcaMin(0.00,StHFCuts::kKaon); // OK

    hfCuts->setCutPrimaryDCAtoVtxMax(1.); //default is 1.0
    hfCuts->setCutVzMax(6.);
    hfCuts->setCutVzVpdVzMax(3.);
    hfCuts->setCutNHitsFitMin(15); //default is 20
    hfCuts->setCutRequireHFT(true);


//    WIDE
//    hfCuts->setCutDcaMin(0.0005,StHFCuts::kPion); // OK
//    hfCuts->setCutDcaMin(0.0005,StHFCuts::kKaon); // OK
//
//    float dcaDaughtersMax = 1;  // maximum toto ide
//    float decayLengthMin  = 0.001; // minimum
//    float decayLengthMax  = 4; //std::numeric_limits<float>::max(); toto ide (cutuje)
//    float cosThetaMin     = 0.;   // minimum
//    float minMass         = 0.4;
//    float maxMass         = 2.4;

    //    WIDER than pt12 tmva
    hfCuts->setCutDcaMin(0.007,StHFCuts::kPion); // OK
    hfCuts->setCutDcaMin(0.007,StHFCuts::kKaon); // OK

    float dcaDaughtersMax = 0.012;  // maximum toto ide
    float decayLengthMin  = 0.005; // minimum
    float decayLengthMax  = 4; //std::numeric_limits<float>::max(); toto ide (cutuje)
    float cosThetaMin     = 0.5;   // minimum
    float minMass         = 0.4;
    float maxMass         = 2.4;



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
    hfCuts->setCutSecondaryPair(dcaDaughtersMax, decayLengthMin, decayLengthMax, cosThetaMin, minMass, maxMass); //ok
    
    //Single track pt
    hfCuts->setCutPtRange(0.1,50.0,StHFCuts::kPion); //0.2 , 50.0 OK
    hfCuts->setCutPtRange(0.1,50.0,StHFCuts::kKaon); //0.2, 50.0 OK
    //TPC setters
    hfCuts->setCutTPCNSigmaPion(3.0); //3 OK
    hfCuts->setCutTPCNSigmaKaon(3.0); //3 OK
    //TOF setters, need to set pt range as well
    hfCuts->setCutTOFDeltaOneOverBeta(0.04, StHFCuts::kKaon); // nastavene = f * (sigmaTOF), sigma TOF je 0.013 OK
    hfCuts->setCutPtotRangeHybridTOF(0.1,50.0,StHFCuts::kKaon); // OK
    hfCuts->setCutTOFDeltaOneOverBeta(0.1, StHFCuts::kPion); //
    hfCuts->setCutPtotRangeHybridTOF(0.1,50.0,StHFCuts::kPion); // OK
    // set refmultCorr
    //StRefMultCorr* grefmultCorrUtil = CentralityMaker::instance()->getgRefMultCorr_P16id();
    //PicoD0AnaMaker->setRefMutCorr(grefmultCorrUtil);
    
    clock_t start = clock(); // getting starting time
    chain->Init();
    
    int total = picoDstMaker->chain()->GetEntries();
    cout << " Total entries = " << total << endl;
    if(nEvents>total) nEvents = total;
    
    for (Int_t i=0; i<nEvents; i++) {
//        if(i%10000==0)       cout << "Working on eventNumber " << i << endl;
            
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
    
    // -- clean up if in read mode
    if (makerMode == StPicoHFMaker::kRead)
        gSystem->Exec(Form("rm -f %s", sInputFile.Data()));
    }
    
    
