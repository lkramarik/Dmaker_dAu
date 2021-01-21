//
// Created by lukas on 8.11.2018.
//
#include "FitPID.h"
#include "TStyle.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TMath.h"
#include "TNtuple.h"
#include "TLine.h"
#include "TList.h"
#include "TLatex.h"
#include "TF1.h"
#include "TCut.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include <iostream>
#include <ctime>
#include <cstdio>
#include <fstream>
#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include "TSystem.h"
#include "TROOT.h"
//class FitPID;
void pidEff(bool tofPidEff,  int bbcMin, int bbcMax, int nTof, float nsigma, float tofInvBeta, float ptTrackCut);

void pidEffKaon() {
    gROOT->ProcessLine(".L FitPID.c++");
    pidEff(true, 0, 950, 1, 1., 0.02, 0.5);

    /*
    int bbcMin=0;
    int bbcMax[]={500,800,950};
    int nTof[]={-1,0,1};

    for (int iBBC = 0; iBBC < 3; ++iBBC) {
        for (int iTOF = 0; iTOF < 3; ++iTOF) {
            pidEff(false, bbcMin, bbcMax[iBBC], nTof[iTOF], 3, 0.03, 0.5);
            pidEff(true, bbcMin, bbcMax[iBBC], nTof[iTOF], 3, 0.03, 0.5);
        }
    }

    float nsigma[]={3, 5};
    float tofInvBeta[]={0.03, 0.05};
    float ptTrackCut[]={0.3, 0.5};

    for (int inSigma = 0; inSigma < 2; ++inSigma) {
        for (int iBeta = 0; iBeta < 2; ++iBeta) {
            for (int iPt = 0; iPt < 2; ++iPt) {
                pidEff(false, 0, 950, 1, nsigma[inSigma], tofInvBeta[iBeta], ptTrackCut[iPt]);
                pidEff(false, 0, 950, 0, nsigma[inSigma], tofInvBeta[iBeta], ptTrackCut[iPt]);
                pidEff(true, 0, 950, 1, nsigma[inSigma], tofInvBeta[iBeta], ptTrackCut[iPt]);
                pidEff(true, 0, 950, 0, nsigma[inSigma], tofInvBeta[iBeta], ptTrackCut[iPt]);
            }
        }
    }
     */

}
void pidEff(bool tofPidEff=true,  int bbcMin=0, int bbcMax=950, int nTof=1, float nsigma=3, float tofInvBeta=0.03, float ptTrackCut=0.5) {
    bool plotPart2 = false;

//    TString input = "/media/lukas/376AD6A434B7392F/work/pid/ntp.noHft.picoPhiAnaMaker.2401.root";
//    TString input = "ntp.picoPhiAnaMaker.nonHftPairs.1709.root";
//    TString input = "/media/lukas/376AD6A434B7392F/work/pid/ntp.picoPhiAnaMaker.3001.root";
    TString input = "/home/lukas/work/dmesons/Dmaker_ndAu/Dmaker_dAu/ntp/ntp.phi.2001.hotspot.root";
    TString inputCut = input+".cutted.root";

    Float_t ptBins[] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.2, 1.4, 1.7, 2., 2.5, 3., 3.5, 4.}; //kaon/
//    Float_t ptBins[] = {0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.5, 3, 4}; //pion
//    Float_t ptBins[] = {0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 1.05, 1.2, 1.4, 1.7, 2.3, 3, 4}; //kaon

    const int nBins = sizeof(ptBins) / sizeof(Float_t);
    TH1F* hTotalCount = new TH1F("hTotalCount","hTotalCount",nBins-1,ptBins);
    TH1F* hPassCount = new TH1F("hPassCount","hPassCount",nBins-1,ptBins);
    TH1F* hTotalCountFunc = new TH1F("hTotalCountFunc","hTotalCountFunc",nBins-1,ptBins);
    TH1F* hPassCountFunc = new TH1F("hPassCountFunc","hPassCountFunc",nBins-1,ptBins);

    Float_t means[nBins], binWidth[nBins], xPt[nBins], massMean, massSigma, eff[nBins], effError[nBins], meansE[nBins], sigmasE[nBins];
    Double_t errorClean, errorAna;
    Float_t sigmas[nBins];

    Float_t massMin, massMax, mean, sigma, ptPairMin, ptPairMax;
    Float_t height = 1000;
    TString pair, pairName, effTypeString;
    TCut cut, cutPair, tof1, tpc1;

    pair = "KK";
    pairName = "KK";
    massMin = 1.0;// Phi to KK
    massMax = 1.04;// Phi to KK
    mean = 1.02;
    sigma = 0.002;
    ptPairMin = 0.2;
    ptPairMax = 10;

    FitPID *fitmass = new FitPID();
    TH1F *hSigmaSignal1 = new TH1F();
    TH1F *hSigmaSignal2 = new TH1F();
    TH1F *hSigmaSignalAna1 = new TH1F();
    TH1F *hSigmaSignalAna2 = new TH1F();

    tof1="abs(pi1_TOFinvbeta)<0.03";
    tpc1="abs(pi1_nSigma<2.)";

    cut = Form("pair_mass>%.3f && pair_mass<%.3f", massMin, massMax);

    float nsigmaInFit=1.;
    cutPair=Form("pair_pt>%.3f && pair_pt<%.3f && bbcRate>%i && bbcRate<=%i && nHftTracks>%i && abs(pi2_nSigma)<%.1f && abs(pi2_TOFinvbeta)<%.2f && pi2_pt>%.1f && "
                 "pi1_dca<1. && pi2_dca<0.5 &&  pair_cosTheta>0.85 && pair_dcaToPv<1. && dcaDaughters<0.5 && pair_decayL<10.",
                 ptPairMin, ptPairMax, bbcMin, bbcMax, nTof, nsigma, tofInvBeta, ptTrackCut);

    if (tofPidEff) effTypeString = "tofPidEff";
    else  effTypeString = "tpc";

    TString dirName=Form("%s_bbc%i_%i_nHft%i_nsigma%.1f_tof%.2f_pt%.1f", effTypeString.Data(), bbcMin, bbcMax, nTof, nsigma, tofInvBeta, ptTrackCut);

    gSystem->Exec(Form("mkdir %s", dirName.Data()));
    gSystem->Exec(Form("mkdir %s/rootFiles", dirName.Data()));
    gSystem->Exec(Form("mkdir %s/img", dirName.Data()));
    gSystem->Exec(Form("mkdir %s/img/KK", dirName.Data()));
    gSystem->Exec(Form("mkdir %s/img/KK/fit", dirName.Data()));

    fitmass->setOutputFileName(Form("%s/rootFiles/mass_", dirName.Data()) + pairName + ".root");
    fitmass->setHeight(15000);
    TH1F *signal = (TH1F*) fitmass->projectSubtractBckg(dirName, input, 50, massMin, massMax, ptPairMin, ptPairMax, pair, cut + cutPair, "pair_mass", "Mass_{%s} (GeV/c^{2})", true);
    fitmass->peakMassFit(dirName, signal, mean, sigma, massMin, massMax, pair, ptPairMin, ptPairMax, "mass");
    massMean = fitmass->getMean();
    massSigma = fitmass->getSigma();

    cut = Form("pair_mass>%.3f && pair_mass<%.3f", massMean-1*massSigma, massMean+1*massSigma);
    if (tofPidEff) cut += "abs(pi1_TOFinvbeta)<93.";

    fitmass->makeTuple(input, cut+cutPair, plotPart2);
    Int_t nBinsInFit = 50;
    int analysedBins = 0;
    height = 1600;
    for (int i = 0; i < nBins-1; ++i) {
        //clean pions:
        FitPID *pid1 = new FitPID();
        pid1->setOutputFileName(Form("%s/rootFiles/nSigma_", dirName.Data()) + pairName + "_1.root");
        cut = Form("pi1_pt>%.3f && pi1_pt<%.3f", ptBins[i], ptBins[i+1]);
        if (tofPidEff) hSigmaSignal1 = (TH1F *) pid1->projectSubtractBckg(dirName, inputCut, nBinsInFit, -5, 5, ptBins[i], ptBins[i + 1], pair, cut, "pi1_nSigma", "Kaon n#sigma^{TPC}", true);
        else hSigmaSignal1 = (TH1F *) pid1->projectSubtractBckg(dirName, inputCut, nBinsInFit, -0.07, 0.07, ptBins[i], ptBins[i + 1], pair, cut, "pi1_TOFinvbeta", "Kaon #delta1/#beta_{TOF}", true);

        pid1->setHeight(height);
        if (tofPidEff) {
            if (i==1 || i==0) pid1->peakFit(dirName, hSigmaSignal1, 0.8, 1, -1.5, 5, pair, ptBins[i], ptBins[i + 1], "nSigma", nsigmaInFit);
            else      pid1->peakFit(dirName, hSigmaSignal1, 0, 1, -5, 5, pair, ptBins[i], ptBins[i + 1], "nSigma", nsigmaInFit);
        } else {
            pid1->peakFit(dirName, hSigmaSignal1, 0, 0.02, -0.07, 0.07, pair, ptBins[i], ptBins[i + 1], "1overBeta", nsigmaInFit);
        }

        binWidth[i] = (ptBins[i + 1] - ptBins[i]) / 2;
        xPt[i] = (ptBins[i + 1] + ptBins[i]) / 2;
        means[i] = pid1->getMean();
        meansE[i] = pid1->getMeanError();
        sigmas[i] = pid1->getSigma();
        sigmasE[i] = pid1->getSigmaError();
        Double_t integralClean = hSigmaSignal1->IntegralAndError(hSigmaSignal1->FindBin(pid1->getMean() - nsigmaInFit*pid1->getSigma()), hSigmaSignal1->FindBin(pid1->getMean() + nsigmaInFit*pid1->getSigma()), errorClean, ""); //number of it without PID cut

        hTotalCount->SetBinContent(i+1, integralClean);
        hTotalCount->SetBinError(i+1, errorClean);
        hTotalCountFunc->SetBinContent(i+1, pid1->getFuncIntegral());
        hTotalCountFunc->SetBinError(i+1, pid1->getFuncIntegralError());

        //tof pions after my PID cut:
        FitPID *pidAna1 = new FitPID();
        pidAna1->setOutputFileName(Form("%s/rootFiles/nSigma_", dirName.Data())+pairName+"_ana_1.root");
        if (tofPidEff) hSigmaSignalAna1 = (TH1F*) pidAna1->projectSubtractBckg(dirName, inputCut, nBinsInFit, -5, 5, ptBins[i], ptBins[i+1], pair, cut+tof1, "pi1_nSigma", "Kaon n#sigma^{TPC}", false);
        else hSigmaSignalAna1 = (TH1F*) pidAna1->projectSubtractBckg(dirName, inputCut, nBinsInFit, -0.07, 0.07, ptBins[i], ptBins[i+1], pair, cut+tpc1, "pi1_TOFinvbeta", "Kaon #delta1/#beta_{TOF}", false);

        Double_t integralAna = hSigmaSignalAna1->IntegralAndError(hSigmaSignalAna1->FindBin(pid1->getMean() - nsigmaInFit*pid1->getSigma()), hSigmaSignalAna1->FindBin(pid1->getMean() + nsigmaInFit*pid1->getSigma()), errorAna, ""); //number of it without PID cut
        hPassCount->SetBinContent(i+1, integralAna);
        hPassCount->SetBinError(i+1, errorAna);

        pidAna1->setHeight(height);
        if (tofPidEff) {
            if (i==1 || i==0) pidAna1->peakFit(dirName, hSigmaSignalAna1, 0.8, 1, -1.5, 5, pair, ptBins[i], ptBins[i + 1], "nSigmaAna", nsigmaInFit);
            else      pidAna1->peakFit(dirName, hSigmaSignalAna1, 0, 1, -5, 5, pair, ptBins[i], ptBins[i + 1], "nSigmaAna", nsigmaInFit);
        } else {
            pidAna1->peakFit(dirName, hSigmaSignalAna1, 0, 0.02, -0.07, 0.07, pair, ptBins[i], ptBins[i + 1], "1overBetaAna", nsigmaInFit);
        }

        hPassCountFunc->SetBinContent(i+1, pidAna1->getFuncIntegral());
        hPassCountFunc->SetBinError(i+1, pidAna1->getFuncIntegralError());


        TLine *left = new TLine(pid1->getMean() - nsigmaInFit*pid1->getSigma(), hSigmaSignalAna1->GetMaximum(), pid1->getMean() - nsigmaInFit*pid1->getSigma(), hSigmaSignalAna1->GetMinimum());
        left->SetLineColor(46);
        TLine *right = new TLine(pid1->getMean() + nsigmaInFit*pid1->getSigma(), hSigmaSignalAna1->GetMaximum(), pid1->getMean() + nsigmaInFit*pid1->getSigma(), hSigmaSignalAna1->GetMinimum());
        right->SetLineColor(46);
        TCanvas *c = new TCanvas("c",Form("%.2f_%.2f", ptBins[i],ptBins[i+1]),1000,900);
        gPad->SetLeftMargin(0.15);
        hSigmaSignalAna1->Draw();
        left->Draw("same");
        right->Draw("same");
        c->SaveAs(Form("%s/img/KK/fit/ana.%.2f_%.2f.png", dirName.Data(), ptBins[i], ptBins[i+1]));
        c->Close();
        left=0x0;
        right=0x0;

        ++analysedBins;
        height=pid1->getHeight();
    }

    TGraphErrors *gMean = new TGraphErrors(analysedBins,xPt,means,binWidth, meansE);
    TGraphErrors *gSigmas = new TGraphErrors(analysedBins,xPt,sigmas,binWidth, sigmasE);

    for (int i = 0; i < nBins; ++i) {
        cout<<hPassCount->GetBinCenter(i+1)<<" "<<hPassCount->GetBinContent(i+1)<<" "<<hTotalCount->GetBinContent(i+1)<<endl;
        if (hPassCount->GetBinContent(i+1)>=hTotalCount->GetBinContent(i+1)) hPassCount->SetBinContent(i+1, hTotalCount->GetBinContent(i+1)-1);
        if (hPassCountFunc->GetBinContent(i+1)>=hTotalCountFunc->GetBinContent(i+1)) hPassCountFunc->SetBinContent(i+1, hTotalCountFunc->GetBinContent(i+1)-1);
    }

    TGraphAsymmErrors *gEff = new TGraphAsymmErrors();
    gEff->Divide(hPassCount,hTotalCount,"cp");

    TGraphAsymmErrors *gEffFunc = new TGraphAsymmErrors();
    gEffFunc->Divide(hPassCountFunc,hTotalCountFunc,"cp");

    TCanvas *c1 = new TCanvas("c1","%.3f_%.3f",1000,900);
    gMean->SetMarkerStyle(20);
    gMean->SetMarkerSize(0.9);
    gMean->SetMarkerColor(kBlack);
    gMean->SetLineColor(kBlack);
    if (tofPidEff) gMean->GetYaxis()->SetTitle("K n#sigma mean [n#sigma]");
    else gMean->GetYaxis()->SetTitle("kaon #delta1/#beta mean [#delta1/#beta]");
    gMean->GetYaxis()->SetTitleOffset(1.1);
    gMean->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    gMean->SetTitle("");
    gMean->Draw("ap");
    c1->SaveAs(Form("%s/rootFiles/meanKaon.root", dirName.Data()));
    c1->Close();

    TCanvas *c2 = new TCanvas("c2","%.3f_%.3f",1000,900);
    gSigmas->SetMarkerStyle(20);
    gSigmas->SetMarkerSize(0.9);
    gSigmas->SetMarkerColor(kBlack);
    gSigmas->SetLineColor(kBlack);
    if (tofPidEff) gSigmas->GetYaxis()->SetTitle("K n#sigma sigma [n#sigma]");
    else gSigmas->GetYaxis()->SetTitle("kaon #delta1/#beta sigma [#delta1/#beta]");
    gSigmas->GetYaxis()->SetTitleOffset(1.1);
    gSigmas->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    gSigmas->SetTitle("");
    gSigmas->Draw("ap");
    c2->SaveAs(Form("%s/rootFiles/sigmaKaon.root", dirName.Data()));
    c2->Close();

    TCanvas *c3 = new TCanvas("c3","%.3f_%.3f",1000,900);
    gEff->SetMarkerStyle(20);
    gEff->SetMarkerSize(0.9);
    gEff->SetMarkerColor(kBlack);
    gEff->SetLineColor(kBlack);
    gEff->GetYaxis()->SetTitle("Efficiency");
    gEff->GetYaxis()->SetTitleOffset(1.1);
    gEff->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    gEff->GetYaxis()->SetRangeUser(0.,1.2);
    gEff->SetTitle("");
    gEff->Draw("ap");
    c3->SaveAs(Form("%s/rootFiles/effTOFKaon.root", dirName.Data()));
    c3->SaveAs(Form("%s/rootFiles/%s_effKaon.png", dirName.Data(), effTypeString.Data()));
    c3->SaveAs(Form("%s/rootFiles/%s_effKaon.eps", dirName.Data(), effTypeString.Data()));
    c3->Close();

    TCanvas *c4 = new TCanvas("c4","%.3f_%.3f",1000,900);
    gEffFunc->SetMarkerStyle(20);
    gEffFunc->SetMarkerSize(0.9);
    gEffFunc->SetMarkerColor(kBlack);
    gEffFunc->SetLineColor(kBlack);
    gEffFunc->GetYaxis()->SetTitle("Efficiency");
    gEffFunc->GetYaxis()->SetTitleOffset(1.1);
    gEffFunc->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    gEffFunc->GetYaxis()->SetRangeUser(0.,1.2);
    gEffFunc->SetTitle("");
    gEffFunc->Draw("ap");
    c4->SaveAs(Form("%s/rootFiles/effTOFKaonFunc.root", dirName.Data()));
    c4->SaveAs(Form("%s/rootFiles/%s_effKaonFunc.png", dirName.Data(), effTypeString.Data()));
    c4->SaveAs(Form("%s/rootFiles/%s_effKaonFunc.eps", dirName.Data(), effTypeString.Data()));
    c4->Close();

    TFile* resOut = new TFile(Form("%s/rootFiles/results_KK.root", dirName.Data()) ,"RECREATE");
    gEff->Write("eff");
    gEffFunc->Write("effFunc");
    gMean->Write("mean");
    gSigmas->Write("sigma");
    resOut->Close();

}
