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
    gROOT->ProcessLine(".L FitPID.c++");

//    TString input = "ntp.picoK0sAnaMaker.root";
    TString inputMass = "/media/lukas/376AD6A434B7392F/work/pid/ntp.noHft.picoK0sAnaMaker.root";
    TString input = "/media/lukas/376AD6A434B7392F/work/pid/ntp.noHft.picoK0sAnaMaker.small.root";

//    TString input = "/gpfs01/star/pwg/lkramarik/Dmaker_dAu/workDir/K0s_last/production/ntp.picoK0sAnaMaker.root";
//    TString input = "/gpfs01/star/pwg/lkramarik/Dmaker_dAu/workDir/Phi_large/production/ntp.picoPhiAnaMaker.root";
//    TString input = "outputBaseName.picoK0sAnaMaker.root";

    Float_t ptBins[] = {0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.5, 3, 4}; //pion
    const int nBins = sizeof(ptBins) / sizeof(Float_t);
    cout << nBins << endl;
    Float_t means[nBins], binWidth[nBins], xPt[nBins], massMean, massSigma, eff[nBins], effError[nBins], meansE[nBins], sigmasE[nBins];
    Double_t errorClean, errorAna;
    Float_t sigmas[nBins];

    Float_t massMin, massMax, mean, sigma, ptPairMin, ptPairMax;
    TString pair, pairName;
    TCut cut, cutPair;
    FitPID *fitmass = new FitPID();

    pair = "#pi#pi";
    pairName = "pipi";

    massMin = 0.42;// K to pipi
    massMax = 0.58;// K to pipi
    mean = 0.5;
    sigma = 0.004;
    ptPairMin = 0.5;
    ptPairMax = 10;
    Float_t height = 40000;
    fitmass->setOutputFileName("mass_" + pairName + ".root");
    fitmass->setHeight(1000000);
    cut = Form("pair_mass>%.3f && pair_mass<%.3f", massMin, massMax);
//    cutPair = Form("pair_pt>%.3f && pair_pt<%.3f && pair_decayL<4", ptPairMin, ptPairMax);
    cutPair = Form("pair_pt>%.3f && pair_pt<%.3f", ptPairMin, ptPairMax);
    TH1F *signal = (TH1F*) fitmass->projectSubtractBckg(inputMass, 50, massMin, massMax, ptPairMin, ptPairMax, pair, cut + cutPair, "pair_mass", "Mass_{%s} (GeV/c^{2})");


    fitmass->peakFit(signal, mean, sigma, massMin, massMax, pair, ptPairMin, ptPairMax, "mass");
    massMean = fitmass->getMean();
    massSigma = fitmass->getSigma();

    int analysedBins = 0;

    for (int i = 0; i < nBins - 1; ++i) {
        //clean pions:
        FitPID *pid1 = new FitPID();
        pid1->setOutputFileName("nSigma_" + pairName + "_1.root");
        cut = Form("pair_mass>%f && pair_mass<%f && pi1_pt>%.3f && pi1_pt<%.3f", massMean - 2 * massSigma, massMean + 2 * massSigma, ptBins[i], ptBins[i + 1]);
        TH1F *hSigmaSignal1 = (TH1F *) pid1->projectSubtractBckg(input, 50, -5, 5, ptBins[i], ptBins[i + 1], pair, cut + cutPair, "pi1_nSigma", "Pion n#sigma^{TPC}");

        FitPID *pid2 = new FitPID();
        pid2->setOutputFileName("nSigma_" + pairName + "_2.root");
        cut = Form("pair_mass>%f && pair_mass<%f && pi2_pt>%.3f && pi2_pt<%.3f", massMean - 2 * massSigma, massMean + 2 * massSigma, ptBins[i], ptBins[i + 1]);
        TH1F *hSigmaSignal2 = (TH1F *) pid2->projectSubtractBckg(input, 50, -5, 5, ptBins[i], ptBins[i + 1], pair, cut + cutPair, "pi2_nSigma", "Pion n#sigma^{TPC}");

        hSigmaSignal1->Add(hSigmaSignal2);
        pid1->setHeight(height);
        pid1->peakFit(hSigmaSignal1, 0, 1, -5, 5, pair, ptBins[i], ptBins[i + 1], "nSigma TPC");

        binWidth[i] = (ptBins[i + 1] - ptBins[i]) / 2;
        xPt[i] = (ptBins[i + 1] + ptBins[i]) / 2;
        means[i] = pid1->getMean();
        meansE[i] = pid1->getMeanError();
        sigmas[i] = pid1->getSigma();
        sigmasE[i] = pid1->getSigmaError();
        Double_t integralClean = hSigmaSignal1->IntegralAndError(hSigmaSignal1->FindBin(pid1->getMean() - 1*pid1->getSigma()), hSigmaSignal1->FindBin(pid1->getMean() + 1*pid1->getSigma()), errorClean, ""); //number of it without PID cut
        cout << integralClean << endl;

        //tof pions after my PID cut:
        FitPID *pidAna1 = new FitPID();
        pidAna1->setOutputFileName("nSigma_"+pairName+"_ana_1.root");
        cut=Form("pair_mass>%f && pair_mass<%f && pi1_pt>%f && pi1_pt<%f && pi1_TOFinvbeta<3", massMean-2*massSigma, massMean+2*massSigma, ptBins[i], ptBins[i+1]); //tof match
//        cut=Form("pair_mass>%f && pair_mass<%f && pi1_pt>%f && pi1_pt<%f && pi1_TOFinvbeta<0.03", massMean-2*massSigma, massMean+2*massSigma, ptBins[i], ptBins[i+1]); //tof match
        TH1F *hSigmaSignalAna1 = (TH1F*) pidAna1->projectSubtractBckg(input, 50, -5, 5, ptBins[i], ptBins[i+1], pair, cut+cutPair, "pi1_nSigma", "Pion n#sigma^{TPC}");

        FitPID *pidAna2 = new FitPID();
        pidAna2->setOutputFileName("nSigma_"+pairName+"_ana_2.root");
        cut=Form("pair_mass>%f && pair_mass<%f && pi2_pt>%f && pi2_pt<%f && pi2_TOFinvbeta<3", massMean-2*massSigma, massMean+2*massSigma, ptBins[i], ptBins[i+1]);
//        cut=Form("pair_mass>%f && pair_mass<%f && pi2_pt>%f && pi2_pt<%f && pi2_TOFinvbeta<0.03", massMean-2*massSigma, massMean+2*massSigma, ptBins[i], ptBins[i+1]);
        TH1F *hSigmaSignalAna2 = (TH1F*) pidAna2->projectSubtractBckg(input, 50, -5, 5, ptBins[i], ptBins[i+1], pair, cut+cutPair, "pi2_nSigma", "Pion n#sigma^{TPC}"); //tof match

        hSigmaSignalAna1->Add(hSigmaSignalAna2);

        Double_t integralAna = hSigmaSignalAna1->IntegralAndError(hSigmaSignalAna1->FindBin(pid1->getMean()-1*pid1->getSigma()),hSigmaSignalAna1->FindBin(pid1->getMean()+1*pid1->getSigma()),errorAna,""); //number of it without PID cut
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
    gMean->GetYaxis()->SetTitle("#pi n#sigma mean [n#sigma]");
    gMean->GetYaxis()->SetTitleOffset(1.1);
    gMean->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    gMean->SetTitle("");
    gMean->Draw("ap");
    c1->SaveAs("mean.root");
    c1->Close();

    TCanvas *c2 = new TCanvas("c2","%.3f_%.3f",1000,900);
    gSigmas->SetMarkerStyle(20);
    gSigmas->SetMarkerSize(0.9);
    gSigmas->SetMarkerColor(kBlack);
    gSigmas->SetLineColor(kBlack);
    gSigmas->GetYaxis()->SetTitle("#pi n#sigma sigma [n#sigma]");
    gSigmas->GetYaxis()->SetTitleOffset(1.1);
    gSigmas->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    gSigmas->SetTitle("");
    gSigmas->Draw("ap");
    c2->SaveAs("sigma.root");
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
    c3->SaveAs("effTOF.root");
    c3->Close();

    TFile* resOut = new TFile("results_pipi.root" ,"RECREATE");
    gEff->Write("eff");
    gMean->Write("mean");
    gSigmas->Write("sigma");
    resOut->Close();
    gSystem->Exec(“hadd -k -f nSigma_pipi.root nSigma_pipi_1.root nSigma_pipi_2.root”);
    gSystem->Exec(“rm nSigma_pipi_1.root nSigma_pipi_2.root”);

}