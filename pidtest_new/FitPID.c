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
#include "TEventList.h"
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
                   mSigma(999), mSigmaE(999), mMean(999), mMeanE(999), mHeight(999), mFuncIntegral(999), mFuncIntegralError(999), mOutputFileName("defaultName.root") {

}

//______________________________________________________________________________________________________________________________________________________
FitPID::~FitPID() {
    // destructor

}

//______________________________________________________________________________________________________________________________________________________
TH1F* FitPID::projectSubtractBckg(TString dirName, TString input, Int_t nBins, Float_t massMin, Float_t massMax, Float_t ptmin, Float_t ptmax, TString pair, TCut setCut, TString variable, TString varName, bool save){
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
    hS->SetMarkerStyle(20);
    hS->SetMarkerColor(1);
    hS->SetLineColor(1);
    hS->SetStats(0);
    hS->SetTitle(Form("p_{T}: %.3f-%.3f GeV/c", ptmin, ptmax));
    hS->GetXaxis()->SetTitleSize(0.045);
    hS->GetYaxis()->SetTitleSize(0.045);
    hS->GetYaxis()->SetTitleOffset(1.25);

    TH1F* hB = (TH1F*)hS->Clone(nameBack);
    hB->SetMarkerColor(46);
    hB->SetLineColor(46);

    ntpS -> Project(nameSig, variable, setCuts);
    ntpB -> Project(nameBack, variable, setCuts);
    TString imgSave = Form("./%s/img/%s/%s_%.3f_%.3f.png", dirName.Data(), pairShort.Data(), variable.Data(), ptmin, ptmax);

    TH1F* hSig = subtractBckg(hS,hB,nameSubtr,dataRes,imgSave,save);

    TList *listOut = new TList();
    listOut->Add(hS);
    listOut->Add(hB);
    listOut->Add(hSig);
    listOut->Write(Form("%s_%.3f_%.3f", pair.Data(), ptmin, ptmax), 1, 0);

    data->Close();
    dataRes->Close();
    return hSig;
}

//______________________________________________________________________________________________________________________________________________________
TH1F* FitPID::subtractBckg(TH1F* hS, TH1F* hB, TString nameSubtr, TFile* outputF, TString imgSave, bool save) {
    outputF->cd();
    TH1F* hSig = (TH1F*)hS->Clone(nameSubtr);
    hSig->SetDirectory(0);
    hSig->SetMarkerColor(9);
    hSig->SetMarkerStyle(20);
    hSig->SetLineColor(9);
//    Double_t integralSig = hSig->Integral(hSig->FindBin(1.03), hSig->FindBin(1.035),"");
//    Double_t integralBkg = hB->Integral(hB->FindBin(1.03), hB->FindBin(1.035),"");

//    Double_t integralSig = hSig->Integral(hSig->FindBin(1.005), hSig->FindBin(1.015),"");
//    Double_t integralBkg = hB->Integral(hB->FindBin(1.005), hB->FindBin(1.015),"");

//    hSig->Add(hB,-integralSig/integralBkg);
    hSig->Add(hB,-1);

    TCanvas *c = new TCanvas("c","c",1000,900);
    gPad->SetLeftMargin(0.15);
    TLegend *legend = new TLegend(0.155,0.789, 0.427, 0.884,"","brNDC");
    legend->AddEntry(hB, "LS background", "pl");
    legend->AddEntry(hS, "US signal", "pl");
    hS->Draw();
    hB->Draw("same");
    legend->Draw("same");
    if (save) c->SaveAs(imgSave);
    c->Close();
    return hSig;
}

//______________________________________________________________________________________________________________________________________________________
bool FitPID::peakFit(TString dirName, TH1F* hToFit, Float_t mean, Float_t sigma, Float_t massMin, Float_t massMax, TString pair, Float_t ptmin, Float_t ptmax, TString varName, Float_t nSigmaLine){
    funPID = new TF1("funPID","pol1(0)+gaus(2)", massMin, massMax);
    funPID->SetParameters(1.,1.,mHeight,mean,sigma);
    funPID->SetLineColor(hToFit->GetMarkerColor());
    funPID->SetLineWidth(5);
    funPID->SetLineStyle(1);
    funPID->SetParName(2,"height");
    funPID->SetParName(3,"mean");
    funPID->SetParName(4,"sigma");
    if (mean!=0)   {
        funPID->SetParLimits(3,0.1*mean,mean+0.9*mean);
    } else {
        funPID->SetParLimits(3,-0.9,3);
    }
    funPID->SetParLimits(4,0,4);

    TCanvas *c = new TCanvas("c","%.3f_%.3f",1000,900);
    gPad->SetLeftMargin(0.15);
    gStyle->SetOptFit(1);
    hToFit->GetYaxis()->SetTitleOffset(1.25);
    r = hToFit->Fit(funPID, "LRS");
    Int_t fitStatus = r;
    cout<<"fit status: "<<fitStatus<<endl;
    mHeight = funPID->GetParameter(2);
    mSigma = funPID->GetParameter(4);
    mSigmaE = funPID->GetParError(4);
    mMean = funPID->GetParameter(3);
    mMeanE = funPID->GetParError(3);
    mFuncIntegral = funPID->Integral(massMin, massMax, 1.e-12);
    mFuncIntegralError = funPID->IntegralError(massMin, massMax, r->GetParams(), r->GetCovarianceMatrix().GetMatrixArray());

    TLine *left = new TLine(mMean - nSigmaLine*mSigma, hToFit->GetMaximum(), mMean - nSigmaLine*mSigma, hToFit->GetMinimum());
    left->SetLineColor(46);
    left->SetLineWidth(4);
    left->SetLineStyle(9);
    TLine *right = new TLine(mMean + nSigmaLine*mSigma, hToFit->GetMaximum(), mMean + nSigmaLine*mSigma, hToFit->GetMinimum());
    right->SetLineStyle(9);
    right->SetLineWidth(4);
    right->SetLineColor(46);

    hToFit->Draw();
    left->Draw("same");
    right->Draw("same");

    pair.ReplaceAll("#","");
    c->SaveAs(Form("./%s/img/%s/fit/%s_%.3f_%.3f.png", dirName.Data(), pair.Data(), varName.Data(), ptmin, ptmax));
    c->Close();

    if (fitStatus!=0) return false;
    return true;
}
//______________________________________________________________________________________________________________________________________________________

float FitPID::getFuncIntegral(Double_t min, Double_t max) {
    return funPID->Integral(min, max, 1.e-12);
}
//______________________________________________________________________________________________________________________________________________________
float FitPID::getFuncIntegralError(Double_t min, Double_t max) {
    return funPID->IntegralError(min, max, r->GetParams(), r->GetCovarianceMatrix().GetMatrixArray());
}
//______________________________________________________________________________________________________________________________________________________
bool FitPID::peakFitResSub(TString dirName, TH1F* hToFit, Float_t mean, Float_t sigma, Float_t massMin, Float_t massMax, TString pair, Float_t ptmin, Float_t ptmax, TString varName, Float_t nSigmaLine){
    TF1* funLS = new TF1("funLS","pol1(0)+gaus(2)", massMin, massMax);
    funLS->SetParameters(1.,1.,mHeight,mean,sigma);
    funLS->SetLineColor(2);
    funLS->SetLineStyle(1);
    funLS->SetParName(2,"height");
    funLS->SetParName(3,"mean");
    funLS->SetParName(4,"sigma");
    if (mean!=0)   {
        funLS->SetParLimits(3,0.1*mean,mean+0.9*mean);
    } else {
        funLS->SetParLimits(3,-0.9,3);
    }
    funLS->SetParLimits(4,0,4);
    funLS->SetLineColor(9);

    TCanvas *c = new TCanvas("c","%.3f_%.3f",1000,900);
    gPad->SetLeftMargin(0.15);
    gStyle->SetOptFit(1);
    hToFit->GetYaxis()->SetTitleOffset(1.25);
    hToFit->Fit(funLS, "LRMN");

    TF1 *funLin = new TF1("funLin","pol1(0)", massMin, massMax);
    funLin->SetParameters(funLS->GetParameter(0),funLS->GetParameter(1));
    TH1F* hResSubtr = (TH1F*)hToFit->Clone();

    for (int i = 0; i < hToFit->GetNbinsX(); ++i) {
        hToFit -> SetBinContent(i, hResSubtr->GetBinContent(i) - funLin->Eval(hResSubtr->GetBinCenter(i)));
    }

    hToFit->Sumw2();

    mHeight = funLS->GetParameter(0);
    mSigma = funLS->GetParameter(2);
    mSigmaE = funLS->GetParError(2);
    mMean = funLS->GetParameter(1);
    mMeanE = funLS->GetParError(1);

    funPID = new TF1("funPID","gaus(0)", massMin, massMax);
    funPID->SetLineColor(hToFit->GetLineColor());
//    funPID->SetLineStyle(7);
    funPID->SetParameters(funLS->GetParameter(2),funLS->GetParameter(3),funLS->GetParameter(4));
    r = hToFit->Fit(funPID, "NLRS", "", 0.4*massMin, 0.4*massMax);
    Int_t fitStatus = r;
    cout<<"fit status: "<<fitStatus<<endl;

    TF1* funDraw = new TF1("funDraw","gaus(0)",massMin, massMax);
    funDraw->SetParameters(funPID->GetParameters());
    funDraw->SetLineColor(hToFit->GetMarkerColor());
    funDraw->SetLineWidth(5);

    hToFit->Draw();
    funDraw->Draw("same");
    mHeight = funPID->GetParameter(0);
    mSigma = funPID->GetParameter(2);
    mSigmaE = funPID->GetParError(2);
    mMean = funPID->GetParameter(1);
    mMeanE = funPID->GetParError(1);
    mFuncIntegral = funPID->Integral(massMin, massMax, 1.e-12);
    mFuncIntegralError = funPID->IntegralError(massMin, massMax, r->GetParams(), r->GetCovarianceMatrix().GetMatrixArray());

    TLine *left = new TLine(mMean - nSigmaLine*mSigma, hToFit->GetMaximum(), mMean - nSigmaLine*mSigma, hToFit->GetMinimum());
    left->SetLineColor(46);
    left->SetLineWidth(4);
    left->SetLineStyle(9);
    TLine *right = new TLine(mMean + nSigmaLine*mSigma, hToFit->GetMaximum(), mMean + nSigmaLine*mSigma, hToFit->GetMinimum());
    right->SetLineStyle(9);
    right->SetLineWidth(4);
    right->SetLineColor(46);
    left->Draw("same");
    right->Draw("same");
    pair.ReplaceAll("#","");
    c->SaveAs(Form("./%s/img/%s/fit/%s_%.3f_%.3f.png", dirName.Data(), pair.Data(), varName.Data(), ptmin, ptmax));
    c->Close();

    if (fitStatus!=0) return false;
    return true;
}


//______________________________________________________________________________________________________________________________________________________
void FitPID::peakMassFit(TString dirName, TH1F* hToFit, Float_t mean, Float_t sigma, Float_t massMin, Float_t massMax, TString pair, Float_t ptmin, Float_t ptmax, TString varName){
    TF1 *funLS = new TF1("funLS","pol1(0)+gaus(2)", massMin, massMax);
    funLS->SetParameters(1.,1.,mHeight,mean,sigma);
    funLS->SetLineStyle(1);
    funLS->SetLineWidth(5);
    funLS->SetParName(2,"height");
    funLS->SetParName(3,"mean");
    funLS->SetParName(4,"sigma");
    if (mean!=0)   {
        funLS->SetParLimits(3,0.1*mean,mean+0.9*mean);
    } else {
        funLS->SetParLimits(3,-0.9,3);
    }
//    funLS->SetParLimits(4,0,0.0035);
    funLS->SetLineColor(9);

    TCanvas *c = new TCanvas("c","%.3_%.3f",1000,900);
    gPad->SetLeftMargin(0.15);
    gStyle->SetOptFit(1);
    hToFit->GetYaxis()->SetTitleOffset(1.25);
    hToFit->Fit(funLS, "LR", "", massMin, massMax);
    hToFit->SetTitle("");

    mHeight = funLS->GetParameter(2);
    mSigma = funLS->GetParameter(4);
    mSigmaE = funLS->GetParError(4);
    mMean = funLS->GetParameter(3);
    mMeanE = funLS->GetParError(3);

    TF1* fitDraw = new TF1("fitDraw","pol1(0)+gaus(2)", massMin, massMax);
    fitDraw->SetParameters(funLS->GetParameters());
    fitDraw->SetLineWidth(5);
    fitDraw->SetLineStyle(1);
    fitDraw->SetLineColor(9);


    TLine *left = new TLine(mMean - 2*mSigma, hToFit->GetMaximum(), mMean - 2*mSigma, hToFit->GetMinimum());
    left->SetLineColor(46);
    left->SetLineWidth(4);
    left->SetLineStyle(9);
    TLine *right = new TLine(mMean + 2*mSigma, hToFit->GetMaximum(), mMean + 2*mSigma, hToFit->GetMinimum());
    right->SetLineStyle(9);
    right->SetLineWidth(4);
    right->SetLineColor(46);

    hToFit->Draw();
//    fitDraw->Draw("same");
    left->Draw("same");
    right->Draw("same");
    mHeight = funLS->GetParameter(2);
    mSigma = funLS->GetParameter(4);
    mSigmaE = funLS->GetParError(4);
    mMean = funLS->GetParameter(3);
    mMeanE = funLS->GetParError(3);
    pair.ReplaceAll("#","");
    c->SaveAs(Form("./%s/img/%s/fit/%s_%.3f_%.3f.png", dirName.Data(), pair.Data(), varName.Data(), ptmin, ptmax));
    c->Close();
}
//______________________________________________________________________________________________________________________________________________________
void FitPID::makeTuple(TString input, TCut cuts, bool plot2Part){
    TFile* data = new TFile(input ,"r");
    TNtuple* ntp[2] = {(TNtuple*)data -> Get("ntp_background"), (TNtuple*)data -> Get("ntp_signal")};
    TString outVars = "pi1_pt:pi1_nSigma:pi1_TOFinvbeta:bbcRate:nTofTracks:pi1_isHft";
    TFile *fileOut = new TFile(input+".cutted.root", "RECREATE");
    cout<<"Smaller nTuple: "<<cuts<<endl;
    TNtuple* ntpOut[2] = {new TNtuple("ntp_background","ntp_background",outVars), new TNtuple("ntp_signal","ntp_signal",outVars)};

    Float_t pi1_pt, pi1_nSigma, pi1_TOFinvbeta, pi2_pt, pi2_nSigma, pi2_TOFinvbeta, bbcRate, nTofTracks, pi1_isHft, pi2_isHft;
    Long64_t indexCut;
    for (int j = 0; j < 2; ++j) {
        ntp[j]->Draw(">>elist",cuts);
        TEventList *elist = (TEventList*)gDirectory->Get("elist");
        ntp[j]->SetBranchAddress("pi1_pt", &pi1_pt);
        ntp[j]->SetBranchAddress("pi1_nSigma", &pi1_nSigma);
        ntp[j]->SetBranchAddress("pi1_isHft", &pi1_isHft);
        ntp[j]->SetBranchAddress("pi1_TOFinvbeta", &pi1_TOFinvbeta);

        ntp[j]->SetBranchAddress("pi2_pt", &pi2_pt);
        ntp[j]->SetBranchAddress("pi2_nSigma", &pi2_nSigma);
        ntp[j]->SetBranchAddress("pi2_isHft", &pi2_isHft);
        ntp[j]->SetBranchAddress("pi2_TOFinvbeta", &pi2_TOFinvbeta);

        ntp[j]->SetBranchAddress("bbcRate", &bbcRate);
        ntp[j]->SetBranchAddress("nTofTracks", &nTofTracks);

        ntp[j]->SetEventList(elist);
        cout<<elist->GetN()<<endl;

        for (int i = 0; i < elist->GetN(); ++i) {
            indexCut = elist->GetEntry(i);
            ntp[j]->GetEntry(indexCut);
            ntpOut[j]->Fill(pi1_pt, pi1_nSigma, pi1_TOFinvbeta, bbcRate, nTofTracks, pi1_isHft);
            if(plot2Part) ntpOut[j]->Fill(pi2_pt, pi2_nSigma, pi2_TOFinvbeta, bbcRate, nTofTracks, pi2_isHft);
        }
    }
    fileOut->cd();
    ntpOut[0]->Write(ntpOut[0]->GetName(), TObject::kOverwrite);
    ntpOut[1]->Write(ntpOut[1]->GetName(), TObject::kOverwrite);
    fileOut->Close();
    data->Close();
}
////______________________________________________________________________________________________________________________________________________________
//void FitPID::makeTuple(TString input, TCut cuts, bool plot2Part){
//    TFile* data = new TFile(input ,"r");
//    TNtuple* ntp[2] = {(TNtuple*)data -> Get("ntp_background"), (TNtuple*)data -> Get("ntp_signal")};
//    TString outVars = "pi1_pt:pi1_nSigmaK:pi1_nSigmaPi:pi1_TOFinvbetaK:pi1_TOFinvbetaPi:bbcRate:nTofTracks:pi1_isHft";
//    TFile *fileOut = new TFile(input+".cutted.root", "RECREATE");
//    cout<<"Smaller nTuple: "<<cuts<<endl;
//    TNtuple* ntpOut[2] = {new TNtuple("ntp_background","ntp_background",outVars), new TNtuple("ntp_signal","ntp_signal",outVars)};
//
//    Float_t pi1_pt, pi1_nSigmaK, pi1_nSigmaPi, pi1_TOFinvbetaK, pi1_TOFinvbetaPi;
//    Float_t pi2_pt, pi2_nSigmaK, pi2_nSigmaPi, pi2_TOFinvbetaK, pi2_TOFinvbetaPi, bbcRate, nTofTracks, pi1_isHft, pi2_isHft;
//
//    Long64_t indexCut;
//    for (int j = 0; j < 2; ++j) {
//        ntp[j]->Draw(">>elist",cuts);
//        TEventList *elist = (TEventList*)gDirectory->Get("elist");
//        ntp[j]->SetBranchAddress("pi1_pt", &pi1_pt);
//        ntp[j]->SetBranchAddress("pi1_nSigmaK", &pi1_nSigmaK);
//        ntp[j]->SetBranchAddress("pi1_nSigmaPi", &pi1_nSigmaPi);
//        ntp[j]->SetBranchAddress("pi1_isHft", &pi1_isHft);
//        ntp[j]->SetBranchAddress("pi1_TOFinvbetaK", &pi1_TOFinvbetaK);
//        ntp[j]->SetBranchAddress("pi1_TOFinvbetaPi", &pi1_TOFinvbetaPi);
//
//        ntp[j]->SetBranchAddress("pi2_pt", &pi2_pt);
//        ntp[j]->SetBranchAddress("pi2_nSigmaK", &pi2_nSigmaK);
//        ntp[j]->SetBranchAddress("pi2_nSigmaPi", &pi2_nSigmaPi);
//        ntp[j]->SetBranchAddress("pi2_isHft", &pi2_isHft);
//        ntp[j]->SetBranchAddress("pi2_TOFinvbetaK", &pi2_TOFinvbetaK);
//        ntp[j]->SetBranchAddress("pi2_TOFinvbetaPi", &pi2_TOFinvbetaPi);
//
//        ntp[j]->SetBranchAddress("bbcRate", &bbcRate);
//        ntp[j]->SetBranchAddress("nTofTracks", &nTofTracks);
//
//        ntp[j]->SetEventList(elist);
//        cout<<elist->GetN()<<endl;
//
//        for (int i = 0; i < elist->GetN(); ++i) {
//            indexCut = elist->GetEntry(i);
//            ntp[j]->GetEntry(indexCut);
//            ntpOut[j]->Fill(pi1_pt, pi1_nSigmaK, pi1_nSigmaPi, pi1_TOFinvbetaK, pi1_TOFinvbetaPi, bbcRate, nTofTracks, pi1_isHft);
//            if(plot2Part) ntpOut[j]->Fill(pi2_pt, pi2_nSigmaK, pi2_nSigmaPi, pi2_TOFinvbetaK, pi2_TOFinvbetaPi, bbcRate, nTofTracks, pi2_isHft);
//        }
//    }
//    fileOut->cd();
//    ntpOut[0]->Write(ntpOut[0]->GetName(), TObject::kOverwrite);
//    ntpOut[1]->Write(ntpOut[1]->GetName(), TObject::kOverwrite);
//    fileOut->Close();
//    data->Close();
//}