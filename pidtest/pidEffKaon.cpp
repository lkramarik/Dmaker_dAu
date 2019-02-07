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

void pidEffKaon() {
//    gSystem->Load("FitPID");
    bool tofPid = true;
    bool plotPart2 = true;
    gROOT->ProcessLine(".L FitPID.c++");
    gSystem->Exec("rm rootFiles/nSigma_KK_1.root rootFiles/nSigma_KK_2.root");
    gSystem->Exec("rm rootFiles/nSigma_KK_ana_1.root rootFiles/nSigma_KK_ana_2.root");
    gSystem->Exec("rm rootFiles/nSigma_KK_ana.root rootFiles/nSigma_KK.root");
    gSystem->Exec("rm rootFiles/mass_KK.root");

//    TString input = "/media/lukas/376AD6A434B7392F/work/pid/ntp.noHft.picoPhiAnaMaker.2401.root";
    TString input = "/media/lukas/376AD6A434B7392F/work/pid/ntp.picoPhiAnaMaker.2901.root";
    TString inputCut = input+".cutted.root";

    Float_t ptBins[] = {0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 1.05, 1.2, 1.4, 1.7, 2.3, 3, 4}; //kaon

    const int nBins = sizeof(ptBins) / sizeof(Float_t);
    cout << nBins << endl;
    Float_t means[nBins], binWidth[nBins], xPt[nBins], massMean, massSigma, eff[nBins], effError[nBins], meansE[nBins], sigmasE[nBins];
    Double_t errorClean, errorAna;
    Float_t sigmas[nBins];

    Float_t massMin, massMax, mean, sigma, ptPairMin, ptPairMax;
    Float_t height = 1000;
    TString pair, pairName;
    TCut cut, cutPair, tof1, tof2, tpc1, tpc2;

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

    fitmass->setOutputFileName("rootFiles/mass_" + pairName + ".root");
    fitmass->setHeight(15000);


    bool hybridTof=false;

    if (hybridTof) {
        tof1="pi1_TOFinvbeta<0.03 || pi1_TOFinvbeta>93";
        tof2="pi2_TOFinvbeta<0.03 || pi2_TOFinvbeta>93";
    }
    else {
        tof1="pi1_TOFinvbeta<93";
        tof2="pi2_TOFinvbeta<93";
    }

    tpc1="abs(pi1_nSigma<2)";
    tpc2="abs(pi2_nSigma<2)";

    cut = Form("pair_mass>%.3f && pair_mass<%.3f", massMin, massMax);
//    cutPair=Form("pair_pt>%.3f && pair_pt<%.3f && pair_cosTheta>0.85 && pi2_pt>0.15 && bbcRate<500 && nTofTracks>0 && abs(pi2_nSigma)<2 && abs(pi2_TOFinvbeta)<0.03", ptPairMin, ptPairMax);
    cutPair=Form("pair_pt>%.3f && pair_pt<%.3f && pair_cosTheta>0.85 && pi2_pt>0.15 && bbcRate<800 && nTofTracks>0", ptPairMin, ptPairMax);

//    cutPair=Form("pair_pt>%.3f && pair_pt<%.3f && pair_cosTheta>0.6 && dcaDaughters<0.3 && pair_dcaToPv<0.3", ptPairMin, ptPairMax);
//    cutPair = Form("pair_pt>%.3f && pair_pt<%.3f && pi2_pt>1 && fabs(pi2_nSigma)<2.", ptPairMin, ptPairMax);
//    cutPair = Form("pair_pt>%.3f && pair_pt<%.3f", ptPairMin, ptPairMax);


    TH1F *signal = (TH1F*) fitmass->projectSubtractBckg(input, 50, massMin, massMax, ptPairMin, ptPairMax, pair, cut + cutPair, "pair_mass", "Mass_{%s} (GeV/c^{2})", true);
    fitmass->peakMassFit(signal, mean, sigma, massMin, massMax, pair, ptPairMin, ptPairMax, "mass");
    massMean = fitmass->getMean();
    massSigma = fitmass->getSigma();

    cut = Form("pair_mass>%.3f && pair_mass<%.3f", massMean-2*massSigma, massMean+2*massSigma);
    fitmass->makeTuple(input, cut+cutPair, plotPart2);

    int analysedBins = 0;
    height = 1600;
    for (int i = 0; i < nBins-1; ++i) {
        //clean pions:
        FitPID *pid1 = new FitPID();
        pid1->setOutputFileName("rootFiles/nSigma_" + pairName + "_1.root");
        cut = Form("pi1_pt>%.3f && pi1_pt<%.3f", ptBins[i], ptBins[i+1]);
        if (tofPid) hSigmaSignal1 = (TH1F *) pid1->projectSubtractBckg(inputCut, 50, -5, 5, ptBins[i], ptBins[i + 1], pair, cut, "pi1_nSigma", "Kaon n#sigma^{TPC}", true);
        else hSigmaSignal1 = (TH1F *) pid1->projectSubtractBckg(inputCut, 50, -0.07, 0.07, ptBins[i], ptBins[i + 1], pair, cut + cutPair, "pi1_TOFinvbeta", "Kaon #delta1/#beta{TOF}", true);

        pid1->setHeight(height);
        if (tofPid) {
            if (i==1) pid1->peakFit(hSigmaSignal1, 0, 1, -1.5, 5, pair, ptBins[i], ptBins[i + 1], "nSigma", 1);
            else      pid1->peakFit(hSigmaSignal1, 0, 1, -5, 5, pair, ptBins[i], ptBins[i + 1], "nSigma", 1);
        } else {
            pid1->peakFit(hSigmaSignal1, 0, 0.02, -0.07, 0.07, pair, ptBins[i], ptBins[i + 1], "1overBeta", 1);
        }

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
        c->SaveAs(Form("img/KK/fit/ana.%.2f_%.2f.png", ptBins[i], ptBins[i+1]));
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
    if (tofPid) gMean->GetYaxis()->SetTitle("K n#sigma mean [n#sigma]");
    else gMean->GetYaxis()->SetTitle("kaon #delta1/#beta mean [#delta1/#beta]");
    gMean->GetYaxis()->SetTitleOffset(1.1);
    gMean->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    gMean->SetTitle("");
    gMean->Draw("ap");
    c1->SaveAs("rootFiles/meanKaon.root");
    c1->Close();

    TCanvas *c2 = new TCanvas("c2","%.3f_%.3f",1000,900);
    gSigmas->SetMarkerStyle(20);
    gSigmas->SetMarkerSize(0.9);
    gSigmas->SetMarkerColor(kBlack);
    gSigmas->SetLineColor(kBlack);
    if (tofPid) gSigmas->GetYaxis()->SetTitle("K n#sigma sigma [n#sigma]");
    else gSigmas->GetYaxis()->SetTitle("kaon #delta1/#beta sigma [#delta1/#beta]");
    gSigmas->GetYaxis()->SetTitleOffset(1.1);
    gSigmas->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    gSigmas->SetTitle("");
    gSigmas->Draw("ap");
    c2->SaveAs("rootFiles/sigmaKaon.root");
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
    c3->SaveAs("rootFiles/effTOFKaon.root");
    c3->Close();

    TFile* resOut = new TFile("rootFiles/results_KK.root" ,"RECREATE");
    gEff->Write("eff");
    gMean->Write("mean");
    gSigmas->Write("sigma");
    resOut->Close();

//    gSystem->Exec("hadd -k -f rootFiles/nSigma_KK.root rootFiles/nSigma_KK_1.root rootFiles/nSigma_KK_2.root");
//    gSystem->Exec("hadd -k -f rootFiles/nSigma_KK_ana.root rootFiles/nSigma_KK_ana_1.root rootFiles/nSigma_KK_ana_2.root");
//    gSystem->Exec("rm rootFiles/nSigma_KK_1.root rootFiles/nSigma_KK_2.root");
//    gSystem->Exec("rm rootFiles/nSigma_KK_ana_1.root rootFiles/nSigma_KK_ana_2.root");

}
