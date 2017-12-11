//D0 mass = 1.864 84
#include <iostream>
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TString.h"
#include "TNtuple.h"

TH1F* hInvMassBackMin = new TH1F("background minus", "background minus", 2000, 0.4, 2.4);
TH1F* hInvMassBackPlus = new TH1F("background plus", "background plus", 2000, 0.4, 2.4);
TH1F* hInvMassSign = new TH1F("signal", "signal", 2000, 0.4, 2.4);
TH1F* hInvMassBack = new TH1F("background", "background", 2000, 0.4, 2.4);

void projectNtp(TFile* data, TFile* dataRes, TString ntpName) {
    Float_t cosDthetaMin=0;
    Float_t D_massMin=0.4;
    Float_t D_massMax=2.4;
    Float_t D_ptMin=2;
    Float_t D_ptMax=3;
    Float_t D_decayLMin=0.02202;
    Float_t pi1_dcaMin=0.008595;
    Float_t k_dcaMin=0.009448;
    Float_t k_nSigmaMax=2;
    Float_t pi1_nSigmaMax=3;
    Float_t pi1_TOFinvbetaMax=999;
    Float_t k_TOFinvbetaMax=0.03;
    Float_t dcaDaughtersMax=0.00525;
    Float_t dca_d0Max=0.003649;

    data->ls();
    TNtuple* ntp = (TNtuple *) data->Get(ntpName);
    Long64_t numberEntr = ntp->GetEntries();
    cout << "Number of entries in Ntuple: " << numberEntr << endl;;

    Float_t flag, D_theta, D_mass, D_pt, D_decayL, k_pt, pi1_pt, pi1_dca, k_dca, k_nSigma, pi1_nSigma, pi1_TOFinvbeta, k_TOFinvbeta, dcaDaughters, pi1_eventId, k_eventId, dca_d0;
    ntp->SetBranchAddress("flag", &flag);
    ntp->SetBranchAddress("D_mass", &D_mass);
    ntp->SetBranchAddress("D_decayL", &D_decayL);
    ntp->SetBranchAddress("D_theta", &D_theta);
    ntp->SetBranchAddress("D_pt", &D_pt);
    ntp->SetBranchAddress("pi1_pt", &pi1_pt);
    ntp->SetBranchAddress("k_pt", &k_pt);
    ntp->SetBranchAddress("pi1_dca", &pi1_dca);
    ntp->SetBranchAddress("k_dca", &k_dca);
    ntp->SetBranchAddress("k_nSigma", &k_nSigma);
    ntp->SetBranchAddress("pi1_nSigma", &pi1_nSigma);
    ntp->SetBranchAddress("pi1_TOFinvbeta", &pi1_TOFinvbeta);
    ntp->SetBranchAddress("k_TOFinvbeta", &k_TOFinvbeta);
    ntp->SetBranchAddress("dcaDaughters", &dcaDaughters);
    ntp->SetBranchAddress("k_eventId", &k_eventId);
    ntp->SetBranchAddress("pi1_eventId", &pi1_eventId);

    TH1F *hpiTOFinvbeta = new TH1F("piTOFinvbeta", "piTOFinvbeta", 2100, -2, 0.1);
    TH1F *hpinSigma = new TH1F("pinSigma", "pinSigma", 800, -4, 4);
    TH1F *hknSigma = new TH1F("knSigma", "knSigma", 800, -4, 4);
    TH1F *hkTOFinvbeta = new TH1F("kTOFinvbeta", "kTOFinvbeta", 2100, -2, 0.1);
    TH1F *hdecayLength = new TH1F("hdecayLength", "hdecayLength", 2000, 0, 0.2);
    TH1F *hpi1_dca = new TH1F("hpi1_dca", "hpi1_dca", 2000, 0, 0.2);
    TH1F *hk_dca = new TH1F("hk_dca", "hk_dca", 2000, 0, 0.2);
    TH1F *hdcaDaughters = new TH1F("hdcaDaughters", "hdcaDaughters", 2000, 0, 0.2);
    TH1F *hcosTheta = new TH1F("hcosTheta", "hcosTheta", 1000, 0, 1);
    TH1F *hdca_d0 = new TH1F("hdca_d0", "hdca_d0", 3000, 0, 0.3);

    for (Long64_t i = 0; i < numberEntr; ++i) {
        if (i % 10000000 == 0) { cout << i << endl; }
        ntp->GetEntry(i);
        if (TMath::Cos(D_theta) > cosDthetaMin) {
            if ((D_mass > D_massMin) && (D_mass < D_massMax)) {
                if ((pi1_dca > pi1_dcaMin) && (k_dca > k_dcaMin) && (dcaDaughters < dcaDaughtersMax) && (D_decayL > D_decayLMin)){
                    if ((fabs(k_nSigma) < k_nSigmaMax) && (fabs(pi1_nSigma) < pi1_nSigmaMax)  ) {
//                        if ((k_TOFinvbeta < 0) || ((k_TOFinvbeta > 0) && (fabs(k_TOFinvbeta) < k_TOFinvbetaMax))){
                        if (((k_TOFinvbeta > 0) && (fabs(k_TOFinvbeta) < k_TOFinvbetaMax)){
                            if ((pi1_TOFinvbeta < 0) || ((pi1_TOFinvbeta > 0) && (fabs(pi1_TOFinvbeta) < pi1_TOFinvbetaMax))){
                                if ((D_pt > D_ptMin) && (D_pt < D_ptMax)) {
                                    dca_d0 = D_decayL * sqrt(1 - TMath::Cos(D_theta) * TMath::Cos(D_theta));
                                    if (dca_d0 < dca_d0Max) {
                                        hpiTOFinvbeta->Fill(pi1_TOFinvbeta);
                                        hkTOFinvbeta->Fill(k_TOFinvbeta);
                                        hpinSigma->Fill(pi1_nSigma);
                                        hknSigma->Fill(k_nSigma);
                                        hdecayLength->Fill(D_decayL);
                                        hpi1_dca->Fill(pi1_dca);
                                        hk_dca->Fill(k_dca);
                                        hdcaDaughters->Fill(dcaDaughters);
                                        hcosTheta->Fill(cos(D_theta));
                                        hdca_d0->Fill(dca_d0);
                                        dca_d0 = 0;

                                        if ((flag == 0) || (flag == 1)) { hInvMassSign->Fill(D_mass); }
                                        if (flag == 4) {
                                            hInvMassBackMin->Fill(D_mass);
                                            hInvMassBack->Fill(D_mass);
                                        }
                                        if (flag == 5) {
                                            hInvMassBackPlus->Fill(D_mass);
                                            hInvMassBack->Fill(D_mass);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    TList *listOut = new TList();

    listOut->Add(hpiTOFinvbeta);
    listOut->Add(hkTOFinvbeta);
    listOut->Add(hpinSigma);
    listOut->Add(hknSigma);
    listOut->Add(hdecayLength);
    listOut->Add(hpi1_dca);
    listOut->Add(hk_dca);
    listOut->Add(hdcaDaughters);
    listOut->Add(hcosTheta);
    listOut->Add(hdca_d0);

    dataRes->cd();
    listOut->Write("res_" + ntpName, 1, 0);

    delete hpiTOFinvbeta;
    delete hpinSigma;
    delete hknSigma;
    delete hkTOFinvbeta;
    delete hdecayLength;
    delete hpi1_dca;
    delete hk_dca;
    delete hdcaDaughters;
    delete hcosTheta;
    delete hdca_d0;
    delete listOut;
    delete ntp;
}

void project(TString input = "testnew.root"){
    TFile* data = new TFile(input ,"r");
    TFile* dataRes = new TFile("res_"+input ,"RECREATE");

    projectNtp(data, dataRes, "ntp_signal");
    projectNtp(data, dataRes, "ntp_background");

    TList* list = (TList*) data -> Get("picoDpmAnaMaker;1");
    TH1F* hStat = (TH1F*) list -> FindObject("hEventStat1");
    dataRes->cd();
    hInvMassSign -> Write();
    hInvMassBack -> Write();
    hInvMassBackPlus -> Write();
    hInvMassBackMin -> Write();
    hStat -> Write();

    delete list;
    data->Close();
    dataRes->Close();
    cout<<"res_"+input<<endl;
    cout<<"done"<<endl;

}




