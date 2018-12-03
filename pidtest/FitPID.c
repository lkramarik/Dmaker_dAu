//
// Created by lukas on 8.11.2018.
//
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
#include "TCut.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include <iostream>
#include <ctime>
#include <cstdio>
#include <fstream>
#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include "FitPID.h"

using namespace std;
using namespace TMath;
ClassImp(FitPID)

FitPID::FitPID() : TObject(),
                   mSigma(999), mSigmaE(999), mMean(999), mMeanE(999), mHeight(999), mOutputFileName("defaultName.root") {

}

FitPID::~FitPID() {
    // destructor

}
TH1F* FitPID::projectSubtractBckg(TString input, Int_t nBins, Float_t massMin, Float_t massMax, Float_t ptmin, Float_t ptmax, TString pair, TCut setCut, TString variable, TString varName){
    std::vector<TCut> mCuts;
    mCuts.push_back(setCut);
    TFile* data = new TFile(input ,"r");
    TNtuple* ntpS = (TNtuple*)data -> Get("ntp_signal");
    TNtuple* ntpB = (TNtuple*)data -> Get("ntp_background");
    TString pairShort = pair;
    pairShort.ReplaceAll("#","");
    TString nameSig = Form("hS_%.3f_%.3f", ptmin, ptmax);
    TString nameBack = Form("hB_%.3f_%.3f", ptmin, ptmax);
    TString nameSubtr = Form("hSsubtr_%.3f_%.3f", ptmin, ptmax);

    TFile* dataRes = new TFile();
    cout<<mOutputFileName<<endl;
    dataRes = new TFile(mOutputFileName,"UPDATE");

    TCut setCuts = "";
    for(unsigned int k = 0; k < mCuts.size(); ++k) {
        setCuts += mCuts[k];
    }
    const char* cut = setCuts;
    cout<<setCuts<<endl;

    dataRes->cd();
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
    TString imgSave = Form("./img/%s/%s_%.3f_%.3f.png", pairShort.Data(), varName.Data(), ptmin, ptmax);

    TH1F* hSig = subtractBckg(hS,hB,nameSubtr,dataRes,imgSave);

    TList *listOut = new TList();
    listOut->Add(hS);
    listOut->Add(hB);
    listOut->Add(hSig);
    listOut->Write(Form("%s_%s_%.3f_%.3f", pair.Data(), variable.Data(), ptmin, ptmax), 1, 0);

    data->Close();
    dataRes->Close();
    return hSig;
}

TH1F* FitPID::subtractBckg(TH1F* hS, TH1F* hB, TString nameSubtr, TFile* outputF, TString imgSave) {
    outputF->cd();
    hSig = (TH1F*)hS->Clone(nameSubtr);
    hSig->SetDirectory(0);
    hSig->SetMarkerColor(9);
    hSig->SetLineColor(9);
    hSig->Add(hB,-1);

    TCanvas *c = new TCanvas("c","c",1000,900);
    TLegend *legend = new TLegend(0.125,0.789, 0.407, 0.884,"","brNDC");
    legend->AddEntry(hB, "LS background", "pl");
    legend->AddEntry(hS, "US signal", "pl");
    hS->Draw();
    hB->Draw("same");
    legend->Draw("same");
    c->SaveAs(imgSave);
    c->Close();
    return hSig;
}

void FitPID::peakFit(TH1F* hToFit, Float_t mean, Float_t sigma, Float_t massMin, Float_t massMax, TString pair, Float_t ptmin, Float_t ptmax, TString varName){
    TF1 *funLS = new TF1("funLS","pol1(0)+gaus(2)", massMin, massMax);
    funLS->SetParameters(1.,1.,mHeight,mean,sigma);
    funLS->SetLineColor(2);
    funLS->SetLineStyle(7);
    funLS->SetLineStyle(1);
    funLS->SetParName(2,"height");
    funLS->SetParName(3,"mean");
    funLS->SetParName(4,"sigma");
    if (mean!=0)   {
        funLS->SetParLimits(3,0.1*mean,mean+1.1*mean);
    } else {
        funLS->SetParLimits(3,-1.5,1.5);
    }
    funLS->SetParLimits(4,0,4);
    funLS->SetLineColor(9);

    TCanvas *c = new TCanvas("c","%.3f_%.3f",1000,900);
    gStyle->SetOptFit(1);
    hToFit->Fit(funLS, "LRM");
    hToFit->Draw();
    mHeight = funLS->GetParameter(2);
    mSigma = funLS->GetParameter(4);
    mSigmaE = funLS->GetParError(4);
    mMean = funLS->GetParameter(3);
    mMeanE = funLS->GetParError(3);
    pair.ReplaceAll("#","");
    c->SaveAs(Form("./img/%s/fit/%s_%.3f_%.3f.png", pair.Data(), varName.Data(), ptmin, ptmax));
    c->Close();
}

