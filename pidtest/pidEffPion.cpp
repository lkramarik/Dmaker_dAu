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

void pidEffPion() {
    gROOT->ProcessLine(".L FitPID.c++");

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

}

void pidEff(bool tofPidEff=true,  int bbcMin=0, int bbcMax=950, int nTof=1, float nsigma=3, float tofInvBeta=0.03, float ptTrackCut=0.5) {
    bool plotPart2 = false;

//    TString input = "ntp.picoK0sAnaMaker.root";
    TString input = "/home/lukas/work/dmesons/Dmaker_ndAu/Dmaker_dAu/ntp/ntp.K0s.hotspot.1901.root";
//    TString input = "/media/lukas/376AD6A434B7392F/work/pid/ntp.picoK0sAnaMaker.1103.root";
//    TString input = "/media/lukas/376AD6A434B7392F/work/pid/ntp.picoK0sAnaMaker.2901.small.root";
    TString inputCut = input+".cutted.root";

    Float_t ptBins[] = {0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.7, 1.9, 2.2, 2.6, 3, 4}; //pion
    const int nBins = sizeof(ptBins) / sizeof(Float_t);
    TH1F* hTotalCount = new TH1F("hTotalCount","hTotalCount",nBins-1,ptBins);
    TH1F* hPassCount = new TH1F("hPassCount","hPassCount",nBins-1,ptBins);

    Float_t means[nBins], binWidth[nBins], xPt[nBins], massMean, massSigma, eff[nBins], effError[nBins], meansE[nBins], sigmasE[nBins];
    Double_t errorClean, errorAna;
    Float_t sigmas[nBins];

    Float_t massMin, massMax, mean, sigma, ptPairMin, ptPairMax;
    TString pair, pairName, effTypeString;
    TCut cut, cutPair, tof1, tpc1;
    FitPID *fitmass = new FitPID();
    TH1F *hSigmaSignal1 = new TH1F();
    TH1F *hSigmaSignal2 = new TH1F();
    TH1F *hSigmaSignalAna1 = new TH1F();
    TH1F *hSigmaSignalAna2 = new TH1F();

    pair = "#pi#pi";
    pairName = "pipi";

    tof1="abs(pi1_TOFinvbeta)<0.03";
    tpc1="abs(pi1_nSigma<3)";

    massMin = 0.42;// K to pipi
    massMax = 0.58;// K to pipi
    mean = 0.5;
    sigma = 0.004;
    ptPairMin = 0.5;
    ptPairMax = 100;
    Float_t height = 40000;

    cut = Form("pair_mass>%.3f && pair_mass<%.3f", massMin, massMax);

    float nsigmaInFit=2;

//    cutPair=Form("pair_pt>%.3f && pair_pt<%.3f && bbcRate>%i && bbcRate<=%i && nHftTracks>%i && abs(pi2_nSigma)<%.1f && abs(pi2_TOFinvbeta)<%.2f && pi2_pt>%.1f && nHftTracks==0", ptPairMin, ptPairMax, bbcMin, bbcMax, nTof, nsigma, tofInvBeta, ptTrackCut);

    cutPair=Form("pair_pt>%.3f && pair_pt<%.3f && bbcRate<%i && nHftTracks>%i && abs(pi2_nSigma)<%.1f && pi2_pt>%.1f && pi1_dca<1. && pi2_dca<1.", ptPairMin, ptPairMax, bbcMax, nTof, nsigma, ptTrackCut);


    if (tofPidEff) effTypeString = "tofPidEff";
    else  effTypeString = "tpc";

    TString dirName=Form("%s_bbc%i_%i_nHft%i_nsigma%.1f_tof%.2f_pt%.1f", effTypeString.Data(), bbcMin, bbcMax, nTof, nsigma, tofInvBeta, ptTrackCut);


//    cutPair=Form("pair_pt>%.3f && pair_pt<%.3f && bbcRate<%i && nHftTracks>%i", ptPairMin, ptPairMax, bbc, nTof);
//    TString dirName=Form("tpc_bbc%i_nHft%i",bbc, nTof);


    gSystem->Exec(Form("mkdir %s", dirName.Data()));
    gSystem->Exec(Form("mkdir %s/rootFiles", dirName.Data()));
    gSystem->Exec(Form("mkdir %s/img", dirName.Data()));
    gSystem->Exec(Form("mkdir %s/img/pipi", dirName.Data()));
    gSystem->Exec(Form("mkdir %s/img/pipi/fit", dirName.Data()));

    fitmass->setOutputFileName(Form("%s/rootFiles/mass_", dirName.Data()) + pairName + ".root");
//    fitmass->setHeight(1000000);
    fitmass->setHeight(15000);
    TH1F *signal = (TH1F*) fitmass->projectSubtractBckg(dirName, input, 50, massMin, massMax, ptPairMin, ptPairMax, pair, cut + cutPair, "pair_mass", "Mass_{%s} (GeV/c^{2})", true);
    fitmass->peakFit(dirName, signal, mean, sigma, massMin, massMax, pair, ptPairMin, ptPairMax, "mass", 2);
    massMean = fitmass->getMean();
    massSigma = fitmass->getSigma();

    cut = Form("pair_mass>%.3f && pair_mass<%.3f", massMean-2*massSigma, massMean+2*massSigma);
    if (tofPidEff) cut += "abs(pi1_TOFinvbeta)<93";

    fitmass->makeTuple(input, cut+cutPair, plotPart2);

    int analysedBins = 0;
    height=1500;
    for (int i = 0; i < nBins-1; ++i) {
        //clean pions:
        FitPID *pid1 = new FitPID();
        pid1->setOutputFileName(Form("%s/rootFiles/nSigma_", dirName.Data()) + pairName + "_1.root");
        cut = Form("pi1_pt>%.3f && pi1_pt<%.3f", ptBins[i], ptBins[i+1]);
        if (tofPidEff) hSigmaSignal1 = (TH1F *) pid1->projectSubtractBckg(dirName, inputCut, 80, -5, 5, ptBins[i], ptBins[i + 1], pair, cut, "pi1_nSigma", "Kaon n#sigma^{TPC}", true);
        else hSigmaSignal1 = (TH1F *) pid1->projectSubtractBckg(dirName, inputCut, 80, -0.07, 0.07, ptBins[i], ptBins[i + 1], pair, cut, "pi1_TOFinvbeta", "Kaon #delta1/#beta_{TOF}", true);

        pid1->setHeight(height);
        if (tofPidEff) pid1->peakFit(dirName, hSigmaSignal1, 0, 1, -5, 5, pair, ptBins[i], ptBins[i + 1], "nSigma", nsigmaInFit);
        else pid1->peakFit(dirName, hSigmaSignal1, 0, 0.02, -0.07, 0.07, pair, ptBins[i], ptBins[i + 1], "1overBeta", nsigmaInFit);

        binWidth[i] = (ptBins[i + 1] - ptBins[i]) / 2;
        xPt[i] = (ptBins[i + 1] + ptBins[i]) / 2;
        means[i] = pid1->getMean();
        meansE[i] = pid1->getMeanError();
        sigmas[i] = pid1->getSigma();
        sigmasE[i] = pid1->getSigmaError();
        Double_t integralClean = hSigmaSignal1->IntegralAndError(hSigmaSignal1->FindBin(pid1->getMean() - nsigmaInFit*pid1->getSigma()), hSigmaSignal1->FindBin(pid1->getMean() + nsigmaInFit*pid1->getSigma()), errorClean, ""); //number of it without PID cut

        hTotalCount->SetBinContent(i+1, integralClean);
        hTotalCount->SetBinError(i+1, errorClean);

        //tof pions after my PID cut:
        FitPID *pidAna1 = new FitPID();
        pidAna1->setOutputFileName(Form("%s/rootFiles/nSigma_", dirName.Data())+pairName+"_ana_1.root");
        if (tofPidEff) hSigmaSignalAna1 = (TH1F*) pidAna1->projectSubtractBckg(dirName, inputCut, 80, -5, 5, ptBins[i], ptBins[i+1], pair, cut+tof1, "pi1_nSigma", "Pion n#sigma^{TPC}", false);
        else hSigmaSignalAna1 = (TH1F*) pidAna1->projectSubtractBckg(dirName, inputCut, 80, -0.07, 0.07, ptBins[i], ptBins[i+1], pair, cut+tpc1, "pi1_TOFinvbeta", "Pion #delta1/#beta_{TOF}", false);

        Double_t integralAna = hSigmaSignalAna1->IntegralAndError(hSigmaSignalAna1->FindBin(pid1->getMean() - nsigmaInFit*pid1->getSigma()), hSigmaSignalAna1->FindBin(pid1->getMean() + nsigmaInFit*pid1->getSigma()), errorAna, ""); //number of it without PID cut

        hPassCount->SetBinContent(i+1, integralAna);
        hPassCount->SetBinError(i+1, errorAna);

        TLine *left = new TLine(pid1->getMean() - nsigmaInFit*pid1->getSigma(), hSigmaSignalAna1->GetMaximum(), pid1->getMean() - nsigmaInFit*pid1->getSigma(), hSigmaSignalAna1->GetMinimum());
        left->SetLineColor(46);
        TLine *right = new TLine(pid1->getMean() + nsigmaInFit*pid1->getSigma(), hSigmaSignalAna1->GetMaximum(), pid1->getMean() + nsigmaInFit*pid1->getSigma(), hSigmaSignalAna1->GetMinimum());
        right->SetLineColor(46);
        TCanvas *c = new TCanvas("c",Form("%.2f_%.2f", ptBins[i],ptBins[i+1]),1000,900);
        gPad->SetLeftMargin(0.15);
        hSigmaSignalAna1->Draw();
        left->Draw("same");
        right->Draw("same");
        c->SaveAs(Form("%s/img/pipi/fit/ana.%.2f_%.2f.png", dirName.Data(), ptBins[i], ptBins[i+1]));
        c->Close();
        left=0x0;
        right=0x0;

        ++analysedBins;
        height=pid1->getHeight();
    }

    for (int i = 0; i < nBins-1; ++i) {
        cout<<hPassCount->GetBinCenter(i+1)<<" "<<hPassCount->GetBinContent(i+1)<<" "<<hTotalCount->GetBinContent(i+1)<<endl;
        if (hPassCount->GetBinContent(i+1)>=hTotalCount->GetBinContent(i+1)) hPassCount->SetBinContent(i+1, hTotalCount->GetBinContent(i+1)-1);
    }

    TGraphErrors *gMean = new TGraphErrors(analysedBins,xPt,means,binWidth, meansE);
    TGraphErrors *gSigmas = new TGraphErrors(analysedBins,xPt,sigmas,binWidth, sigmasE);
    TGraphAsymmErrors *gEff = new TGraphAsymmErrors();
    gEff->Divide(hPassCount,hTotalCount,"cp");

    TCanvas *c1 = new TCanvas("c1","%.3f_%.3f",1000,900);
    gMean->SetMarkerStyle(20);
    gMean->SetMarkerSize(0.9);
    gMean->SetMarkerColor(kBlack);
    gMean->SetLineColor(kBlack);
    if (tofPidEff) gMean->GetYaxis()->SetTitle("#pi n#sigma mean [n#sigma]");
    else gMean->GetYaxis()->SetTitle("#pi #Delta1/#beta mean [#Delta1/#beta]");
    gMean->GetYaxis()->SetTitleOffset(1.1);
    gMean->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    gMean->SetTitle("");
    gMean->Draw("ap");
    c1->SaveAs(Form("%s/rootFiles/meanPion.root", dirName.Data()));
    c1->Close();

    TCanvas *c2 = new TCanvas("c2","%.3f_%.3f",1000,900);
    gSigmas->SetMarkerStyle(20);
    gSigmas->SetMarkerSize(0.9);
    gSigmas->SetMarkerColor(kBlack);
    gSigmas->SetLineColor(kBlack);
    if (tofPidEff) gSigmas->GetYaxis()->SetTitle("#pi n#sigma sigma [n#sigma]");
    else gSigmas->GetYaxis()->SetTitle("#pi #Delta1/#beta sigma [#Delta1/#beta]");
    gSigmas->GetYaxis()->SetTitleOffset(1.1);
    gSigmas->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    gSigmas->SetTitle("");
    gSigmas->Draw("ap");
    c2->SaveAs(Form("%s/rootFiles/sigmaPion.root", dirName.Data()));
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
    c3->SaveAs(Form("%s/rootFiles/effTOFPion.root", dirName.Data()));
    c3->SaveAs(Form("%s/rootFiles/%s_effPion.png", dirName.Data(), effTypeString.Data()));
    c3->SaveAs(Form("%s/rootFiles/%s_effPion.eps", dirName.Data(), effTypeString.Data()));
    c3->Close();

    TFile* resOut = new TFile(Form("%s/rootFiles/results_pipi.root", dirName.Data()),"RECREATE");
    gEff->Write("eff");
    gMean->Write("mean");
    gSigmas->Write("sigma");
    resOut->Close();
}