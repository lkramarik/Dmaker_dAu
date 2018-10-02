#include "TStyle.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TMath.h"
#include "TNtuple.h"
#include "TList.h"
#include "TLatex.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include <iostream>
#include <ctime>
#include <cstdio>
#include <fstream>
#include <unistd.h>
using namespace std;
using namespace TMath;
TH1F* subtractBckg(TString, Float_t, Float_t, Float_t, Float_t, string);
TH1F* hSig = new TH1F();
TFile* dataRes = new TFile();

void pidEff() {
//    bool kaons = true;
//    bool pions = false;
    bool kaons = false;
    bool pions = true;

//    TString input = "ntp.picoPhiAnaMaker.root";
    TString input = "ntp.picoK0sAnaMaker.root";
//    TString input = "outputBaseName.picoK0sAnaMaker.root";

    Float_t ptBins[] = {0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.5, 3, 4};
    int nBins = sizeof(ptBins);

    Float_t massMin, massMax;
    string pair;
    if (kaons) {
        pair = "KK";
        massMin = 1.0;// Phi to KK
        massMax = 1.04;// Phi to KK
    }
    if (pions) {
        pair = "#pi#pi";
        massMin = 0.4;// K to pipi
        massMax = 0.6;// K to pipi
    }

    TString output = "dum.root";
    dataRes = new TFile(output,"RECREATE");

    int i = 0;
//    for (int i = 0; i < 2; ++i) {
//    for (int i = 0; i < nBins-1; ++i) {
        TH1F* signal = (TH1F*) subtractBckg(input, massMin, massMax, ptBins[i], ptBins[i+1], pair);
//        signal -> Draw();
//
}

TH1F* subtractBckg(TString input, Float_t massMin, Float_t massMax, Float_t ptmin, Float_t ptmax, string pair){
    TFile* data = new TFile(input ,"r");
    TNtuple* ntpS = (TNtuple*)data -> Get("ntp_signal");
    TNtuple* ntpB = (TNtuple*)data -> Get("ntp_background");

    TString nameSig = Form("hS_%.3f_%.3f", ptmin, ptmax);
    TString nameBack = Form("hB_%.3f_%.3f", ptmin, ptmax);
    TString nameSubtr = Form("hSsubtr_%.3f_%.3f", ptmin, ptmax);

    TString setCuts = Form("pair_mass>%f && pair_mass<%f", massMin, massMax);
    dataRes->cd();
    TList *listOut = new TList();
    TH1F* hS = new TH1F(nameSig, nameSig, 500, massMin, massMax);
    hS->Sumw2();
    hS->GetYaxis()->SetTitle("Counts");
    hS->GetXaxis()->SetTitle(Form("Mass_{%s} (GeV/c^{2})", pair.c_str()));
    hS->SetMarkerStyle(2);
    hS->SetMarkerColor(1);
    hS->SetLineColor(1);
    hS->SetStats(0);
    hS->SetTitle(Form("p_{T}: %.3f-%.3f GeV/c", ptmin, ptmax));
    hS->GetXaxis()->SetTitleSize(0.045);
    hS->GetYaxis()->SetTitleSize(0.045);
    hS->GetYaxis()->SetTitleOffset(1.1);

    TH1F* hB = (TH1F*)hS->Clone(nameBack);
    hB->SetMarkerColor(46);

    ntpS -> Project(nameSig, "pair_mass", setCuts);
    ntpB -> Project(nameBack, "pair_mass", setCuts);

    hSig = (TH1F*)hS->Clone(nameSubtr);
    hB->SetMarkerColor(9);
    hSig->Add(hB,-1);
    hSig->Draw();
    listOut->Add(hS);
    listOut->Add(hB);
    listOut->Add(hSig);
    listOut->Write(Form("%s_%.3f_%.3f", pair.c_str(), ptmin, ptmax), 1, 0);

    data->Close();
    return hSig;
}
//
//void peakFit(TH1F* hToFit){
//
//
//
//    TF1 *funLS = new TF1("funLS","pol1(0)+gaus(2)", fitRMin,fitRMax);
//    funLS->SetParameters(1.,1.,1.,1.84,0.01);
//    funLS->SetLineColor(2);
//    funLS->SetLineStyle(7);
//    funLS->SetLineStyle(1);
//    funLS->SetParName(2,"height");
//    funLS->SetParName(3,"mean");
//    funLS->SetParName(4,"sigma");
//    funLS->SetParLimits(3,rotwthmin,rotwthmax);
//    funLS->SetLineColor(9);
//
//    mSigma =
//    pairReco -> SetSigma(sigma);
//    pairReco -> SetMean(mean);
//}


