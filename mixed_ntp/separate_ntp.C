#include <string>
#include "TFile.h"
#include "TNtuple.h"
#include "TList.h"
//void separate_ntp(TString input = "2018-07-03_11-06_983.picoMEtree.root"){
void separate_ntp(TString input = "/home/lukas/work/dmesons/Dmaker_dAu/res_analyse/ntp/ntp_1307.root"){
    TFile *file=new TFile(input,"read");
    TString ntpnames[2] = {"ntp_signal", "ntp_background"};
//    TString ntpnames[4] = {"ntp_signal_SE", "ntp_signal_ME", "ntp_background_SE", "ntp_background_ME"};


    TList *list = (TList*)file -> Get("picoD0AnaMaker;1");
//    TList *list = (TList*)file -> Get("picoMixedEventMaker;1");
//    for (int i = 0; i < 4; ++i) {
    for (int i = 0; i < 2; ++i) {
        TNtuple *ntp = (TNtuple*)file -> Get(ntpnames[i]);

        TFile *fileOut=new TFile(ntpnames[i]+".all.root","recreate");
//        TFile *fileOut=new TFile(ntpnames[i]+"_"+input,"recreate");
        list -> Clone() -> Write("picoD0AnaMaker",1,0);
//        list -> Clone() -> Write("picoMixedEventMaker",1,0);
        ntp -> CloneTree() -> Write(ntpnames[i],TObject::kOverwrite);
//        ntp -> Clone() -> Write(ntpnames[i],1,0);
        fileOut -> Close();
//        delete fileOut;

    }

    file -> Close();
//  fileOut -> cd();





}

