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

void pidEffPion() {
//    gSystem->Load("FitPID");
    bool tofPid = true;
    bool plotPart2 = true;
    gROOT->ProcessLine(".L FitPID.c++");

    gSystem->Exec("rm rootFiles/nSigma_pipi_1.root rootFiles/nSigma_pipi_2.root");
    gSystem->Exec("rm rootFiles/nSigma_pipi_ana_1.root rootFiles/nSigma_pipi_ana_2.root");
    gSystem->Exec("rm rootFiles/nSigma_pipi_ana.root rootFiles/nSigma_pipi.root");
    gSystem->Exec("rm rootFiles/mass_pipi.root");

//    TString input = "ntp.picoK0sAnaMaker.root";
    TString inputMass = "/media/lukas/376AD6A434B7392F/work/pid/ntp.picoK0sAnaMaker.2901.root";
    TString input = "/media/lukas/376AD6A434B7392F/work/pid/ntp.picoK0sAnaMaker.2901.small.root";
    TString inputCut = input+".cutted.root";

    Float_t ptBins[] = {0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.5, 3, 4}; //pion
    const int nBins = sizeof(ptBins) / sizeof(Float_t);
    cout << nBins << endl;
    Float_t means[nBins], binWidth[nBins], xPt[nBins], massMean, massSigma, eff[nBins], effError[nBins], meansE[nBins], sigmasE[nBins];
    Double_t errorClean, errorAna;
    Float_t sigmas[nBins];

    Float_t massMin, massMax, mean, sigma, ptPairMin, ptPairMax;
    TString pair, pairName;
    TCut cut, cutPair, tof1, tof2, tpc1, tpc2;
    FitPID *fitmass = new FitPID();
    TH1F *hSigmaSignal1 = new TH1F();
    TH1F *hSigmaSignal2 = new TH1F();
    TH1F *hSigmaSignalAna1 = new TH1F();
    TH1F *hSigmaSignalAna2 = new TH1F();

    pair = "#pi#pi";
    pairName = "pipi";

    bool hybridTof=false;

    if (hybridTof) {
        tof1="pi1_TOFinvbeta<0.03 || pi1_TOFinvbeta>93";
        tof2="pi2_TOFinvbeta<0.03 || pi2_TOFinvbeta>93";
    }
    else {
        tof1="pi1_TOFinvbeta<93";
        tof2="pi2_TOFinvbeta<93";
    }

    tpc1="abs(pi1_nSigma<3)";
    tpc2="abs(pi2_nSigma<3)";

    massMin = 0.42;// K to pipi
    massMax = 0.58;// K to pipi
    mean = 0.5;
    sigma = 0.004;
    ptPairMin = 0.5;
    ptPairMax = 10;
    Float_t height = 40000;
    fitmass->setOutputFileName("rootFiles/mass_" + pairName + ".root");
    fitmass->setHeight(1000000);
    cut = Form("pair_mass>%.3f && pair_mass<%.3f", massMin, massMax);
//    cutPair = Form("pair_pt>%.3f && pair_pt<%.3f && pair_decayL<4", ptPairMin, ptPairMax);
//    cutPair = Form("pair_pt>%.3f && pair_pt<%.3f && nTofTracks>0", ptPairMin, ptPairMax);
//    cutPair=Form("pair_pt>%.3f && pair_pt<%.3f && bbcRate<800 && nTofTracks>0", ptPairMin, ptPairMax);

//    cutPair=Form("pair_pt>%.3f && pair_pt<%.3f && bbcRate<500 && nTofTracks>0 && abs(pi2_nSigma)<2 && abs(pi2_TOFinvbeta)<0.03 && pi2_pt>0.5", ptPairMin, ptPairMax);
    cutPair=Form("pair_pt>%.3f && pair_pt<%.3f && bbcRate<500 && nTofTracks>0 && abs(pi2_nSigma)<1 && abs(pi2_TOFinvbeta)<0.02 && pi2_pt>0.5", ptPairMin, ptPairMax);

    TH1F *signal = (TH1F*) fitmass->projectSubtractBckg(inputMass, 50, massMin, massMax, ptPairMin, ptPairMax, pair, cut + cutPair, "pair_mass", "Mass_{%s} (GeV/c^{2})", true);
    fitmass->peakFit(signal, mean, sigma, massMin, massMax, pair, ptPairMin, ptPairMax, "mass", 2);
    massMean = fitmass->getMean();
    massSigma = fitmass->getSigma();

    cut = Form("pair_mass>%.3f && pair_mass<%.3f", massMean-2*massSigma, massMean+2*massSigma);
    fitmass->makeTuple(input, cut+cutPair, plotPart2);

    int analysedBins = 0;
    height=1500;
    for (int i = 0; i < nBins-1; ++i) {
        //clean pions:
        FitPID *pid1 = new FitPID();
        pid1->setOutputFileName("rootFiles/nSigma_" + pairName + "_1.root");
        cut = Form("pi1_pt>%.3f && pi1_pt<%.3f", ptBins[i], ptBins[i+1]);
        if (tofPid) hSigmaSignal1 = (TH1F *) pid1->projectSubtractBckg(inputCut, 50, -5, 5, ptBins[i], ptBins[i + 1], pair, cut, "pi1_nSigma", "Pion n#sigma^{TPC}", true);
        else hSigmaSignal1 = (TH1F *) pid1->projectSubtractBckg(inputCut, 50, -0.07, 0.07, ptBins[i], ptBins[i + 1], pair, cut + cutPair, "pi1_TOFinvbeta", "Pion #Delta1/#beta{TOF}", true);

        pid1->setHeight(height);
        if (tofPid) pid1->peakFit(hSigmaSignal1, 0, 1, -5, 5, pair, ptBins[i], ptBins[i + 1], "nSigma", 1);
        else pid1->peakFit(hSigmaSignal1, 0, 0.02, -0.07, 0.07, pair, ptBins[i], ptBins[i + 1], "1overBeta", 1);

        binWidth[i] = (ptBins[i + 1] - ptBins[i]) / 2;
        xPt[i] = (ptBins[i + 1] + ptBins[i]) / 2;
        means[i] = pid1->getMean();
        meansE[i] = pid1->getMeanError();
        sigmas[i] = pid1->getSigma();
        sigmasE[i] = pid1->getSigmaError();
        Double_t integralClean = hSigmaSignal1->IntegralAndError(hSigmaSignal1->FindBin(pid1->getMean() - 1*pid1->getSigma()), hSigmaSignal1->FindBin(pid1->getMean() + 1*pid1->getSigma()), errorClean, ""); //number of it without PID cut

        //tof pions after my PID cut:
        FitPID *pidAna1 = new FitPID();
        pidAna1->setOutputFileName("rootFiles/nSigma_"+pairName+"_ana_1.root");
        if (tofPid) hSigmaSignalAna1 = (TH1F*) pidAna1->projectSubtractBckg(inputCut, 50, -5, 5, ptBins[i], ptBins[i+1], pair, cut+tof1, "pi1_nSigma", "Kaon n#sigma^{TPC}", false);
        else TH1F *hSigmaSignalAna1 = (TH1F*) pidAna1->projectSubtractBckg(input, 50, -0.07, 0.07, ptBins[i], ptBins[i+1], pair, cut+cutPair+tpc1, "pi1_TOFinvbeta", "Kaon #delta1/#beta_{TOF}", false);

        Double_t integralAna = hSigmaSignalAna1->IntegralAndError(hSigmaSignalAna1->FindBin(pid1->getMean() - 1*pid1->getSigma()), hSigmaSignalAna1->FindBin(pid1->getMean() + 1*pid1->getSigma()), errorAna, ""); //number of it without PID cut
        TLine *left = new TLine(pid1->getMean() - 1*pid1->getSigma(), hSigmaSignalAna1->GetMaximum(), pid1->getMean() - 1*pid1->getSigma(), hSigmaSignalAna1->GetMinimum());
        left->SetLineColor(46);
        TLine *right = new TLine(pid1->getMean() + 1*pid1->getSigma(), hSigmaSignalAna1->GetMaximum(), pid1->getMean() + 1*pid1->getSigma(), hSigmaSignalAna1->GetMinimum());
        right->SetLineColor(46);
        TCanvas *c = new TCanvas("c",Form("%.2f_%.2f", ptBins[i],ptBins[i+1]),1000,900);
        hSigmaSignalAna1->Draw();
        left->Draw("same");
        right->Draw("same");
        c->SaveAs(Form("img/pipi/fit/ana.%.2f_%.2f.png", ptBins[i], ptBins[i+1]));
        c->Close();
        left=0x0;
        right=0x0;

        cout<<integralAna<<endl;
        eff[i]=(float)integralAna/(float)integralClean;
        effError[i]=sqrt(errorAna*errorAna/(pow(integralClean,2)) + errorClean*errorClean*integralAna*integralAna/(pow(integralClean,4)));
        cout<<eff[i]<<" pm "<<effError[i]<<endl;
        ++analysedBins;
        height=pid1->getHeight();
    }

    TGraphErrors *gMean = new TGraphErrors(analysedBins,xPt,means,binWidth, meansE);
    TGraphErrors *gSigmas = new TGraphErrors(analysedBins,xPt,sigmas,binWidth, sigmasE);
    TGraphErrors *gEff = new TGraphErrors(analysedBins,xPt,eff,binWidth, effError);

    TCanvas *c1 = new TCanvas("c1","%.3f_%.3f",1000,900);
    gMean->SetMarkerStyle(20);
    gMean->SetMarkerSize(0.9);
    gMean->SetMarkerColor(kBlack);
    gMean->SetLineColor(kBlack);
    if (tofPid) gMean->GetYaxis()->SetTitle("#pi n#sigma mean [n#sigma]");
    else gMean->GetYaxis()->SetTitle("#pi #Delta1/#beta mean [#Delta1/#beta]");
    gMean->GetYaxis()->SetTitleOffset(1.1);
    gMean->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    gMean->SetTitle("");
    gMean->Draw("ap");
    c1->SaveAs("rootFiles/meanPion.root");
    c1->Close();

    TCanvas *c2 = new TCanvas("c2","%.3f_%.3f",1000,900);
    gSigmas->SetMarkerStyle(20);
    gSigmas->SetMarkerSize(0.9);
    gSigmas->SetMarkerColor(kBlack);
    gSigmas->SetLineColor(kBlack);
    if (tofPid) gSigmas->GetYaxis()->SetTitle("#pi n#sigma sigma [n#sigma]");
    else gSigmas->GetYaxis()->SetTitle("#pi #Delta1/#beta sigma [#Delta1/#beta]");
    gSigmas->GetYaxis()->SetTitleOffset(1.1);
    gSigmas->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    gSigmas->SetTitle("");
    gSigmas->Draw("ap");
    c2->SaveAs("rootFiles/sigmaPion.root");
    c2->Close();

    TCanvas *c3 = new TCanvas("c3","%.3f_%.3f",1000,900);
    gEff->SetMarkerStyle(20);
    gEff->SetMarkerSize(0.9);
    gEff->SetMarkerColor(kBlack);
    gEff->SetLineColor(kBlack);
    gEff->GetYaxis()->SetTitle("Efficiency");
    gEff->GetYaxis()->SetTitleOffset(1.1);
    gEff->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    gEff->SetTitle("");
    gEff->Draw("ap");
    c3->SaveAs("rootFiles/effTOFPion.root");
    c3->Close();

    TFile* resOut = new TFile("rootFiles/results_pipi.root" ,"RECREATE");
    gEff->Write("eff");
    gMean->Write("mean");
    gSigmas->Write("sigma");
    resOut->Close();

//    gSystem->Exec("hadd -k -f rootFiles/nSigma_pipi.root rootFiles/nSigma_pipi_1.root rootFiles/nSigma_pipi_2.root");
//    gSystem->Exec("hadd -k -f rootFiles/nSigma_pipi_ana.root rootFiles/nSigma_pipi_ana_1.root rootFiles/nSigma_pipi_ana_2.root");
//    gSystem->Exec("rm rootFiles/nSigma_pipi_1.root rootFiles/nSigma_pipi_2.root");
//    gSystem->Exec("rm rootFiles/nSigma_pipi_ana_1.root rootFiles/nSigma_pipi_ana_2.root");


}