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
TH1F* subtractBckg(TString, Int_t, Float_t, Float_t, Float_t, Float_t, TString, TString, TString, TString);
void peakFit(TH1F*, Float_t, Float_t, Float_t, TString, Float_t, Float_t);
TH1F* hSig = new TH1F();
Float_t mSigma;
Float_t mMean;

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

    Float_t massMin, massMax, mean;
    TString pair;
    if (kaons) {
        pair = "KK";
        massMin = 1.0;// Phi to KK
        massMax = 1.04;// Phi to KK
        mean = 1.1;
    }
    if (pions) {
        pair = "#pi#pi";
        massMin = 0.4;// K to pipi
        massMax = 0.6;// K to pipi
        mean = 0.5;
    }


    int i = 0;
//    for (int i = 0; i < 2; ++i) {
//    for (int i = 0; i < nBins-1; ++i) {
        TH1F* signal = (TH1F*) subtractBckg(input, 500, massMin, massMax, ptBins[i], ptBins[i+1], pair, "", "pair_mass", "Mass_{%s} (GeV/c^{2})");
        peakFit(signal, mean, massMin, massMax, pair, ptBins[i], ptBins[i+1]);
        TH1F* nSigmaSignal = (TH1F*) subtractBckg(input, 50, -5, 5, ptBins[i], ptBins[i+1], pair, Form("pair_mass>%f && pair_mass<%f && ", mMean-mSigma, mMean+mSigma), "pi1_nSigma", "pi1_nSigma");
//        peakFit(hSigmaSignal, 0, )

}

TH1F* subtractBckg(TString input, Int_t nBins, Float_t massMin, Float_t massMax, Float_t ptmin, Float_t ptmax, TString pair, TString setCuts, TString variable, TString varName){
    TFile* data = new TFile(input ,"r");
    TNtuple* ntpS = (TNtuple*)data -> Get("ntp_signal");
    TNtuple* ntpB = (TNtuple*)data -> Get("ntp_background");

    TString nameSig = Form("hS_%.3f_%.3f", ptmin, ptmax);
    TString nameBack = Form("hB_%.3f_%.3f", ptmin, ptmax);
    TString nameSubtr = Form("hSsubtr_%.3f_%.3f", ptmin, ptmax);
    TFile* dataRes = new TFile();
    TString output = "dum.root";
    dataRes = new TFile(output,"RECREATE");

    setCuts += Form("pair_mass>%f && pair_mass<%f", massMin, massMax);
    dataRes->cd();
    TList *listOut = new TList();
    TH1F* hS = new TH1F(nameSig, nameSig, nBins, massMin, massMax);
    hS->Sumw2();
    hS->GetYaxis()->SetTitle("Counts");
    hS->GetXaxis()->SetTitle(Form(varName, pair.Data()));
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
    hB->SetLineColor(46);

    ntpS -> Project(nameSig, variable, setCuts);
    ntpB -> Project(nameBack, variable, setCuts);

    hSig = (TH1F*)hS->Clone(nameSubtr);
    hSig->SetDirectory(0);
    hSig->SetMarkerColor(9);
    hSig->SetLineColor(9);
    hSig->Add(hB,-1);

    listOut->Add(hS);
    listOut->Add(hB);
    listOut->Add(hSig);
    listOut->Write(Form("%s_%s_%.3f_%.3f", pair.Data(), variable.Data(), ptmin, ptmax), 1, 0);

    TCanvas *c = new TCanvas("c","%.3f_%.3f",1000,900);
    TLegend *legend = new TLegend(0.125,0.789, 0.407, 0.884,"","brNDC");
    legend->AddEntry(hB, "LS background", "pl");
    legend->AddEntry(hS, "US signal", "pl");
    hS->Draw();
    hB->Draw("same");
    legend->Draw("same");
    pair.ReplaceAll("#","");
    c->SaveAs(Form("./img/%s/%.3f_%.3f.png", pair.Data(), ptmin, ptmax));
    c->Close();
    data->Close();
    dataRes->Close();
    return hSig;
}

void peakFit(TH1F* hToFit, Float_t mean, Float_t massMin, Float_t massMax, TString pair, Float_t ptmin, Float_t ptmax){
    TF1 *funLS = new TF1("funLS","pol1(0)+gaus(2)", massMin, massMax);
    funLS->SetParameters(1.,1.,1.,mean,0.01);
    funLS->SetLineColor(2);
    funLS->SetLineStyle(7);
    funLS->SetLineStyle(1);
    funLS->SetParName(2,"height");
    funLS->SetParName(3,"mean");
    funLS->SetParName(4,"sigma");
    funLS->SetParLimits(3,mean-0.4*mean,mean+0.4*mean);
    funLS->SetLineColor(9);

    TCanvas *c = new TCanvas("c","%.3f_%.3f",1000,900);
    gStyle->SetOptFit(1);
    hToFit->Fit(funLS, "LRM");
    hToFit->Draw();
    mSigma = funLS->GetParameter(4);
    mMean = funLS->GetParameter(3);
    pair.ReplaceAll("#","");
    c->SaveAs(Form("./img/%s/fit/%.3f_%.3f.png", pair.Data(), ptmin, ptmax));
    c->Close();
}

