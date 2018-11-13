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

using namespace std;
using namespace TMath;
TH1F* projectSubtractBckg(TString, Int_t, Float_t, Float_t, Float_t, Float_t, TString, TCut, TString, TString);
TH1F* subtractBckg(TH1F*, TH1F*, TString, TFile*, TString);
void peakFit(TH1F*, Float_t, Float_t, Float_t, Float_t, TString, Float_t, Float_t, TString);
TH1F* hSig = new TH1F();
Float_t mSigma = 999;
Float_t mSigmaE = 999;
Float_t mMean = 999;
Float_t mMeanE = 999;
Float_t mHeight = 10000;
TString mOutputFileName = "defaultName.root";

inline void setOutputFileName(TString name){mOutputFileName=name;}


void pidEff() {
    bool kaons = true;
    bool pions = false;
//    bool kaons = false;
//    bool pions = true;

    TString input = "ntp.picoPhiAnaMaker.small.root";
//    TString input = "ntp.picoK0sAnaMaker.root";
//    TString input = "/gpfs01/star/pwg/lkramarik/Dmaker_dAu/workDir/Phi_large/production/ntp.picoPhiAnaMaker.root";
//    TString input = "outputBaseName.picoK0sAnaMaker.root";

//    Float_t ptBins[] = {2.2, 2.5, 3};
//    Float_t ptBins[] = {0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.5, 3, 4}; //pion
    Float_t ptBins[] = {0.2, 0.35, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.5, 3, 4}; //kaon
    const int nBins = sizeof(ptBins)/sizeof(Float_t);
    cout<<nBins<<endl;
    Float_t means[nBins], binWidth[nBins], xPt[nBins], massMean, massSigma, eff[nBins], effError[nBins], meansE[nBins], sigmasE[nBins];
    Double_t errorClean, errorAna;
    Float_t sigmas[nBins];

    Float_t massMin, massMax, mean, sigma, ptPairMin, ptPairMax;
    TString pair, pairName;
    TCut cut, cutPair;
    if (kaons) {
        pair = "KK";
        pairName = "KK";
        massMin = 1.0;// Phi to KK
        massMax = 1.04;// Phi to KK
        mean = 1.02;
        sigma = 0.04;
        ptPairMin = 0.9;
        ptPairMax = 10;

        setOutputFileName("mass_"+pairName+".root");
        cut=Form("pair_mass>%.3f && pair_mass<%.3f", massMin, massMax);
//    cutPair=Form("pair_pt>%.3f && pair_pt<%.3f && pair_cosTheta>0.6 && dcaDaughters<0.3 && pair_dcaToPv<0.3", ptPairMin, ptPairMax);
        cutPair=Form("pair_pt>%.3f && pair_pt<%.3f && pi1_pt>0.5 && pi1_nSigma<2", ptPairMin, ptPairMax);
        TH1F *signal = (TH1F*) projectSubtractBckg(input,50, massMin, massMax, ptPairMin, ptPairMax, pair, cut+cutPair, "pair_mass", "Mass_{%s} (GeV/c^{2})");
        peakFit(signal, mean, sigma, massMin, massMax, pair, ptPairMin, ptPairMax, "mass");
        massMean = mMean;
        massSigma = mSigma;
    }
    if (pions) {
        pair = "#pi#pi";
        pairName = "pipi";

        massMin = 0.4;// K to pipi
        massMax = 0.6;// K to pipi
        mean = 0.5;
        sigma = 0.004;
        ptPairMin = 0.5;
        ptPairMax = 10;

        setOutputFileName("mass_"+pairName+".root");
        TH1F *signal = (TH1F*) projectSubtractBckg(input,50, massMin, massMax, ptPairMin, ptPairMax, pair, Form("pair_mass>%.3f && pair_mass<%.3f && pair_pt>%.3f && pair_pt<%.3f", massMin, massMax, ptPairMin, ptPairMax), "pair_mass", "Mass_{%s} (GeV/c^{2})");
        peakFit(signal, mean, sigma, massMin, massMax, pair, ptPairMin, ptPairMax, "mass");
        massMean = mMean;
        massSigma = mSigma;
    }

    int analysedBins=0;



//    for (int i = 6; i < 8; ++i) {
    for (int i = 0; i < nBins-1; ++i) {
        //clean pions:
        setOutputFileName("nSigma_"+pairName+"_1.root");
        cut=Form("pair_mass>%f && pair_mass<%f && pi1_pt>%.3f && pi1_pt<%.3f", massMean-2*massSigma, massMean+2*massSigma, ptBins[i], ptBins[i+1]);
        TH1F *hSigmaSignal1 = (TH1F*) projectSubtractBckg(input, 50, -5, 5, ptBins[i], ptBins[i+1], pair, cut+cutPair, "pi1_nSigma", "pi1_nSigma");

        setOutputFileName("nSigma_"+pairName+"_2.root");
        cut=Form("pair_mass>%f && pair_mass<%f && pi2_pt>%.3f && pi2_pt<%.3f", massMean-2*massSigma, massMean+2*massSigma, ptBins[i], ptBins[i+1]);
        TH1F *hSigmaSignal2 = (TH1F*) projectSubtractBckg(input, 50, -5, 5, ptBins[i], ptBins[i+1], pair, cut+cutPair, "pi2_nSigma", "pi2_nSigma");

        hSigmaSignal1->Add(hSigmaSignal2);
        peakFit(hSigmaSignal1, 0, 1,-5, 5, pair, ptBins[i], ptBins[i + 1], "pi_nSigma");

        binWidth[i] = (ptBins[i+1]-ptBins[i])/2;
        xPt[i] = (ptBins[i+1]+ptBins[i]) / 2;
        means[i] = mMean;
        meansE[i] = mMeanE;
        sigmas[i] = mSigma;
        sigmasE[i] = mSigmaE;
        Double_t integralClean = hSigmaSignal1->IntegralAndError(hSigmaSignal1->FindBin(mMean-2*mSigma),hSigmaSignal1->FindBin(mMean+2*mSigma),errorClean,""); //number of it without PID cut
        cout<<integralClean<<endl;

        //tof pions after my PID cut:
        setOutputFileName("nSigma_"+pairName+"_ana_1.root");
        cut=Form("pair_mass>%f && pair_mass<%f && pi1_pt>%f && pi1_pt<%f && pi1_TOFinvbeta<0.03 && pi1_TOFinvbeta>0", massMean-2*massSigma, massMean+2*massSigma, ptBins[i], ptBins[i+1]);
        TH1F *hSigmaSignalAna1 = (TH1F*) projectSubtractBckg(input, 50, -5, 5, ptBins[i], ptBins[i+1], pair, cut+cutPair, "pi1_nSigma", "pi1_nSigma_ana");

        setOutputFileName("nSigma_"+pairName+"_ana_2.root");
        cut=Form("pair_mass>%f && pair_mass<%f && pi2_pt>%f && pi2_pt<%f && pi2_TOFinvbeta<0.03 && pi2_TOFinvbeta>0", massMean-2*massSigma, massMean+2*massSigma, ptBins[i], ptBins[i+1]);
        TH1F *hSigmaSignalAna2 = (TH1F*) projectSubtractBckg(input, 50, -5, 5, ptBins[i], ptBins[i+1], pair, cut+cutPair, "pi2_nSigma", "pi2_nSigma_ana");

        hSigmaSignalAna1->Add(hSigmaSignalAna2);
        peakFit(hSigmaSignalAna1, 0, 1,-5, 5, pair, ptBins[i], ptBins[i+1], "pi_nSigma_ana");

        Double_t integralAna = hSigmaSignalAna1->IntegralAndError(hSigmaSignalAna1->FindBin(mMean-2*mSigma),hSigmaSignalAna1->FindBin(mMean+2*mSigma),errorAna,""); //number of it without PID cut
        cout<<integralAna<<endl;
        eff[i]=(float)integralAna/(float)integralClean;
        effError[i]=sqrt(errorAna*errorAna/(pow(integralClean,2)) + errorClean*errorClean*integralAna*integralAna/(pow(integralClean,4)));
        cout<<eff[i]<<" pm "<<effError[i]<<endl;
        ++analysedBins;
    }

    TGraphErrors *gMean = new TGraphErrors(analysedBins,xPt,means,binWidth, meansE);
    TGraphErrors *gSigmas = new TGraphErrors(analysedBins,xPt,sigmas,binWidth, sigmasE);
    TGraphErrors *gEff = new TGraphErrors(analysedBins,xPt,eff,binWidth, effError);

    TCanvas *c1 = new TCanvas("c1","%.3f_%.3f",1000,900);
    gMean->SetMarkerStyle(20);
    gMean->SetMarkerSize(0.9);
    gMean->SetMarkerColor(kBlack);
    gMean->SetLineColor(kBlack);
    gMean->GetYaxis()->SetTitle("nSigma TPC mean");
    gMean->GetYaxis()->SetTitleOffset(1.1);
    gMean->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    gMean->SetTitle("");
    gMean->Draw("ap");
    c1->SaveAs("mean.pdf");
    c1->Close();

    TCanvas *c2 = new TCanvas("c2","%.3f_%.3f",1000,900);
    gSigmas->SetMarkerStyle(20);
    gSigmas->SetMarkerSize(0.9);
    gSigmas->SetMarkerColor(kBlack);
    gSigmas->SetLineColor(kBlack);
    gSigmas->GetYaxis()->SetTitle("nSigma TPC sigma");
    gSigmas->GetYaxis()->SetTitleOffset(1.1);
    gSigmas->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    gSigmas->SetTitle("");
    gSigmas->Draw("ap");
    c2->SaveAs("sigma.pdf");
    c2->Close();

    TCanvas *c3 = new TCanvas("c3","%.3f_%.3f",1000,900);
    gEff->SetMarkerStyle(20);
    gEff->SetMarkerSize(0.9);
    gEff->SetMarkerColor(kBlack);
    gEff->SetLineColor(kBlack);
    gEff->GetYaxis()->SetTitle("nSigma TOF eff");
    gEff->GetYaxis()->SetTitleOffset(1.1);
    gEff->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    gEff->SetTitle("");
    gEff->Draw("ap");
    c3->SaveAs("effTOF.pdf");
    c3->Close();
}

TH1F* projectSubtractBckg(TString input, Int_t nBins, Float_t massMin, Float_t massMax, Float_t ptmin, Float_t ptmax, TString pair, TCut setCut, TString variable, TString varName){
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

TH1F* subtractBckg(TH1F* hS, TH1F* hB, TString nameSubtr, TFile* outputF, TString imgSave) {
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

void peakFit(TH1F* hToFit, Float_t mean, Float_t sigma, Float_t massMin, Float_t massMax, TString pair, Float_t ptmin, Float_t ptmax, TString varName){
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

