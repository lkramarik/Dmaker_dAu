#include <string>
#include "TFile.h"
#include "TList.h"
void remove_ntp(TString input){
  TFile *file=new TFile(input,"read");

  //   gDirectory->Delete(object_to_remove.c_str());
//   file->Close();

 TList *list = (TList*)file -> Get("PicoD0AnaMaker;1");
 TFile *fileOut=new TFile("small_"+input,"recreate");
//  fileOut -> cd();
 list -> Clone() -> Write("PicoD0AnaMaker",1,0);

    
}

