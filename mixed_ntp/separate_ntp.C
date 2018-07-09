#include <string>
#include "TFile.h"
#include "TNtuple.h"
#include "TList.h"
void separate_ntp(TString input = "2018-07-03_11-06_983.picoMEtree.root"){
    TFile *file=new TFile(input,"read");
    TString ntpnames[4] = {"ntp_signal_SE", "ntp_signal_ME", "ntp_background_SE", "ntp_background_ME"};


    TList *list = (TList*)file -> Get("picoMixedEventMaker;1");
    for (int i = 0; i < 4; ++i) {
        TNtuple *ntp = (TNtuple*)file -> Get(ntpnames[i]);

        TFile *fileOut=new TFile(ntpnames[i]+"_"+input,"recreate");
        list -> Clone() -> Write("picoMixedEventMaker",1,0);
        ntp -> CloneTree() -> Write();
//        ntp -> Clone() -> Write(ntpnames[i],1,0);
        fileOut -> Close();
//        delete fileOut;

    }

    file -> Close();
//  fileOut -> cd();





}

