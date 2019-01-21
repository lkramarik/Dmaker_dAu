#include "TSystem.h"

#ifndef __CINT__
#include "StMuDSTMaker/COMMON/macros/loadSharedLibraries.C"
#endif

extern TSystem* gSystem;

void loadSharedHFLibraries() {

  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();

  gSystem->Load("StBTofUtil");
  gSystem->Load("StPicoEvent");
  gSystem->Load("StPicoDstMaker");
  gSystem->Load("StPicoCutsBase");
  gSystem->Load("StPicoHFMaker");
  gSystem->Load("StPicoD0AnaMaker");
  gSystem->Load("StPicoD0V2AnaMaker");
  gSystem->Load("StPicoPiPiMaker");
  gSystem->Load("StPicoKKMaker");
  gSystem->Load("StPicoQAMaker");
  gSystem->Load("StPicoMixedEventMaker");

  cout << " loading of shared HF libraries are done" << endl;
 }
