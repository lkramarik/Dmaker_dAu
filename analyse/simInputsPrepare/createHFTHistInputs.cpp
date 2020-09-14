#include "StAnaCuts.h"


void createHFTHistInputs() {
    TFile fDca1("dca.3007.primary.dca1cm.root");
    TFile *outHistF = new TFile("dcaxy_vs_dcaz_NEW.root", "RECREATE");
    outHistF->SetCompressionSettings(0); //needed to open file in ROOT5







}