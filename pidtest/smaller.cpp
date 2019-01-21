//
// Created by lukas on 6.12.2018.
//
// Divide NTP to signalLike and backgroundLike

#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TString.h"
#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TNtuple.h"

using namespace std;

void smaller(TString input="/media/lukas/376AD6A434B7392F/work/pid/ntp.noHft.picoK0sAnaMaker.root") {
    TFile* data = new TFile(input ,"r");
    TFile *fileOut = new TFile(input+"small.root", "RECREATE");  // output root file

    TString ntpVars = "pi1_pt:pi1_nSigma:pi1_TOFinvbeta:pi1_phi:pi2_pt:pi2_nSigma:pi2_TOFinvbeta:pi2_phi:dcaDaughters:pair_cosTheta:pair_pt:pair_mass";

    TNtuple *ntpO[2] = {new TNtuple("ntp_signal","Ks_Signal", ntpVars), new TNtuple("ntp_background","Ks_background",ntpVars)};
    TString ntpName[2] = {"ntp_signal", "ntp_background"};
    Float_t pi1_pt, pi1_nSigma, pi1_TOFinvbeta, pi2_pt, pi2_nSigma, pi2_TOFinvbeta, dcaDaughters, pair_cosTheta, pair_pt, pair_mass, pi1_phi, pi2_phi;
    for (int i = 0; i < 2; ++i) {
        TNtuple *ntp = (TNtuple*) data->Get(ntpName[i]);

        ntp->SetBranchAddress("pi1_pt", &pi1_pt);
        ntp->SetBranchAddress("pi1_nSigma", &pi1_nSigma);
        ntp->SetBranchAddress("pi1_TOFinvbeta", &pi1_TOFinvbeta);
        ntp->SetBranchAddress("pi1_phi", &pi1_phi);

        ntp->SetBranchAddress("pi2_pt", &pi2_pt);
        ntp->SetBranchAddress("pi2_nSigma", &pi2_nSigma);
        ntp->SetBranchAddress("pi2_TOFinvbeta", &pi2_TOFinvbeta);
        ntp->SetBranchAddress("pi2_phi", &pi2_phi);

        ntp->SetBranchAddress("dcaDaughters", &dcaDaughters);
        ntp->SetBranchAddress("pair_cosTheta", &pair_cosTheta);
        ntp->SetBranchAddress("pair_pt", &pair_pt);
        ntp->SetBranchAddress("pair_mass", &pair_mass);

        cout<<"ok"<<endl;

        for (long int j = 0; j < ntp->GetEntries(); j++) {
            ntp->GetEntry(j);
            if((pair_mass>0.48) && (pair_mass<0.52))
            ntpO[i]->Fill(pi1_pt, pi1_nSigma, pi1_TOFinvbeta, pi1_phi, pi2_pt, pi2_nSigma, pi2_TOFinvbeta, pi2_phi, dcaDaughters, pair_cosTheta, pair_pt, pair_mass);

        }

    }


    fileOut->cd();
    ntpO[0]->Write(ntpO[0]->GetName(), TObject::kOverwrite);
    ntpO[1]->Write(ntpO[1]->GetName(), TObject::kOverwrite);

    fileOut->Close();
    data->Close();
    cout<<"Done."<<endl;
}

