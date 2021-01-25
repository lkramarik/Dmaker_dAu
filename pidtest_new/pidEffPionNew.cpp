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
TString pidEff(bool tofPidEff,  int bbcMin, int bbcMax, int nTof, float nsigma, float tofInvBeta, float ptTrackCut, TCut addiCut);
TString pidEff(bool tofPidEff,  int bbcMin, int bbcMax, int nTof, float nsigma, float tofInvBeta, float ptTrackCut);

void pidEffPionNew() {
    gROOT->ProcessLine(".L FitPID.c++");
//    cout<<pidEff(true, 0, 950, 1, 2., 0.02, 0.,"")<<endl;

    TString dir;
    dir=pidEff(true, 0, 950, 1, 100, 100000, 0., "pi1_isHft>0 && pi2_isHft>0");
    gSystem->Exec(Form("mv %s plotPart2/", dir.Data()));
    dir=pidEff(false, 0, 950, 1, 100, 100000, 0., "pi1_isHft>0 && pi2_isHft>0");
    gSystem->Exec(Form("mv %s plotPart2/", dir.Data()));

    /*
    TCut cutHft[] = {"pi1_isHft>-100", "pi1_isHft>0", "pi1_isHft<1"};
    TString folder[] = {"all", "hft", "nonhft"};

    TString dir;
    for (int i = 0; i < 3; ++i) {
        int bbcMin=0;
        int bbcMax[]={500,800,950};
        int nTof[]={-1,0,1};

//        for (int iBBC = 0; iBBC < 3; ++iBBC) {
//            for (int iTOF = 0; iTOF < 3; ++iTOF) {
//                dir=pidEff(false, bbcMin, bbcMax[iBBC], nTof[iTOF], 3, 0.03, 0.5, cutHft[i]);
//                gSystem->Exec(Form("mv %s %s/", dir.Data(), folder[i].Data()));
//                dir=pidEff(true, bbcMin, bbcMax[iBBC], nTof[iTOF], 3, 0.03, 0.5, cutHft[i]);
//                gSystem->Exec(Form("mv %s %s/", dir.Data(), folder[i].Data()));
//            }
//        }

        float nsigma[]={3, 5, 200.};
        float tofInvBeta[]={0.03, 0.05, 200.};
        float ptTrackCut[]={0., 0.3, 0.5};

        for (int inSigma = 0; inSigma < 3; ++inSigma) {
            for (int iBeta = 0; iBeta < 3; ++iBeta) {
                for (int iPt = 0; iPt < 3; ++iPt) {
                    dir=pidEff(false, 0, 950, 1, nsigma[inSigma], tofInvBeta[iBeta], ptTrackCut[iPt], cutHft[i]);
                    gSystem->Exec(Form("mv %s %s/", dir.Data(), folder[i].Data()));

                    dir=pidEff(false, 0, 950, 0, nsigma[inSigma], tofInvBeta[iBeta], ptTrackCut[iPt], cutHft[i]);
                    gSystem->Exec(Form("mv %s %s/", dir.Data(), folder[i].Data()));

                    dir=pidEff(true, 0, 950, 1, nsigma[inSigma], tofInvBeta[iBeta], ptTrackCut[iPt], cutHft[i]);
                    gSystem->Exec(Form("mv %s %s/", dir.Data(), folder[i].Data()));

                    dir=pidEff(true, 0, 950, 0, nsigma[inSigma], tofInvBeta[iBeta], ptTrackCut[iPt], cutHft[i]);
                    gSystem->Exec(Form("mv %s %s/", dir.Data(), folder[i].Data()));
                }
            }
        }

    }
     */


}

TString pidEff(bool tofPidEff=true,  int bbcMin=0, int bbcMax=950, int nTof=1, float nsigma=3, float tofInvBeta=0.03, float ptTrackCut=0.5) {
    return pidEff(tofPidEff,  bbcMin, bbcMax, nTof, nsigma, tofInvBeta, ptTrackCut, "");
}

TString pidEff(bool tofPidEff=true,  int bbcMin=0, int bbcMax=950, int nTof=1, float nsigma=3, float tofInvBeta=0.03, float ptTrackCut=0.5, TCut additionalCut="") {
    bool plotPart2 = true;

    TString input = "/home/lukas/work/dmesons/Dmaker_ndAu/Dmaker_dAu/ntp/ntp.K0s.hotspot.1901.root";
    TString inputCut = input+".cutted.root";

    Float_t ptBins[] = {0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.7, 1.9, 2.2, 2.6, 3, 4}; //pion
    const int nBins = sizeof(ptBins) / sizeof(Float_t);
    TH1F* hTotalCount = new TH1F("hTotalCount","hTotalCount",nBins-1,ptBins);
    TH1F* hPassCount = new TH1F("hPassCount","hPassCount",nBins-1,ptBins);
    TH1F* hTotalCountFunc = new TH1F("hTotalCountFunc","hTotalCountFunc",nBins-1,ptBins);
    TH1F* hPassCountFunc = new TH1F("hPassCountFunc","hPassCountFunc",nBins-1,ptBins);

    Float_t means[nBins], binWidth[nBins], xPt[nBins], massMean, massSigma, eff[nBins], effError[nBins], meansE[nBins], sigmasE[nBins];
    Double_t errorClean, errorAna, anaCut;
    Float_t sigmas[nBins];

    Float_t massMin, massMax, mean, sigma, ptPairMin, ptPairMax;
    TString pair, pairName, effTypeString;
    TCut cut, cutPair;
    FitPID *fitmass = new FitPID();
    TH1F *hSigmaSignal1 = new TH1F();
    TH1F *hSigmaSignal2 = new TH1F();
    TH1F *hSigmaSignalAna1 = new TH1F();
    TH1F *hSigmaSignalAna2 = new TH1F();

    pair = "#pi#pi";
    pairName = "pipi";
    massMin = 0.465;// K to pipi
    massMax = 0.53;// K to pipi
    mean = 0.495;
    sigma = 0.0047;
    ptPairMin = 0.5;
    ptPairMax = 100;
    Float_t height = 40000;

    cut = Form("pair_mass>%.3f && pair_mass<%.3f", massMin, massMax);

    float nsigmaInFit=1;

    cutPair=Form("pair_pt>%.3f && pair_pt<%.3f && bbcRate<%i && nHftTracks>%i && abs(pi2_nSigma)<%.1f && pi2_pt>%.1f && pi1_dca<1.", ptPairMin, ptPairMax, bbcMax, nTof, nsigma, ptTrackCut);
    cutPair+=additionalCut;

    if (tofPidEff) {
        effTypeString = "tofPidEff";
        anaCut=0.03;
    }
    else  {
        effTypeString = "tpc";
        anaCut=3;
    }

    TString dirName=Form("pi_%s_bbc%i_%i_nHft%i_nsigma%.1f_tof%.2f_pt%.1f", effTypeString.Data(), bbcMin, bbcMax, nTof, nsigma, tofInvBeta, ptTrackCut);

    gSystem->Exec(Form("mkdir %s", dirName.Data()));
    gSystem->Exec(Form("mkdir %s/rootFiles", dirName.Data()));
    gSystem->Exec(Form("mkdir %s/img", dirName.Data()));
    gSystem->Exec(Form("mkdir %s/img/pipi", dirName.Data()));
    gSystem->Exec(Form("mkdir %s/img/pipi/fit", dirName.Data()));

    fitmass->setOutputFileName(Form("%s/rootFiles/mass_", dirName.Data()) + pairName + ".root");
//    fitmass->setHeight(1000000);
    fitmass->setHeight(150000);
    TH1F *signal = (TH1F*) fitmass->projectSubtractBckg(dirName, input, 50, massMin, massMax, ptPairMin, ptPairMax, pair, cut + cutPair, "pair_mass", "Mass_{%s} (GeV/c^{2})", true);
//    fitmass->peakFit(dirName, signal, mean, sigma, massMin, massMax, pair, ptPairMin, ptPairMax, "mass", 2);
    fitmass->peakMassFit(dirName, signal, mean, sigma, massMin, massMax, pair, ptPairMin, ptPairMax, "mass");

    massMean = fitmass->getMean();
    massSigma = fitmass->getSigma();

    cut = Form("pair_mass>%.3f && pair_mass<%.3f", massMean-2*massSigma, massMean+2*massSigma);
    if (tofPidEff) cut += "abs(pi1_TOFinvbeta)<93";

    fitmass->makeTuple(input, cut+cutPair, plotPart2);
    Int_t nBinsInFit = 50;
    int analysedBins = 0;
    height=1500;
    bool fitOk = false;

    for (int i = 0; i < nBins-1; ++i) {
        //clean pions:
        FitPID *pid1 = new FitPID();
        pid1->setOutputFileName(Form("%s/rootFiles/nSigma_", dirName.Data()) + pairName + "_1.root");
        cut = Form("pi1_pt>%.3f && pi1_pt<%.3f", ptBins[i], ptBins[i+1]);
        if (tofPidEff) hSigmaSignal1 = (TH1F *) pid1->projectSubtractBckg(dirName, inputCut, nBinsInFit, -0.06, 0.06, ptBins[i], ptBins[i + 1], pair, cut, "pi1_TOFinvbeta", "Pion #Delta1/#beta_{TOF}", true);
        else hSigmaSignal1 = (TH1F *) pid1->projectSubtractBckg(dirName, inputCut, nBinsInFit, -4, 4, ptBins[i], ptBins[i + 1], pair, cut, "pi1_nSigma", "Pion n#sigma^{TPC}", true);

        pid1->setHeight(height);
        if (tofPidEff) {
            if (i==0) fitOk=pid1->peakFitResSub(dirName, hSigmaSignal1, 0.01, 0.02, -0.06, 0.06, pair, ptBins[i], ptBins[i + 1], "1overBeta", nsigmaInFit);
            fitOk=pid1->peakFitResSub(dirName, hSigmaSignal1, 0, 0.01, -0.06, 0.06, pair, ptBins[i], ptBins[i + 1], "1overBeta", nsigmaInFit);
        } else {
//            if (i==1 || i==0) pid1->peakFitResSub(dirName, hSigmaSignal1, 0.8, 1, -1.5, 5, pair, ptBins[i], ptBins[i + 1], "nSigma", nsigmaInFit);
//            else      pid1->peakFitResSub(dirName, hSigmaSignal1, 0, 1, -5, 5, pair, ptBins[i], ptBins[i + 1], "nSigma", nsigmaInFit);
            fitOk=pid1->peakFitResSub(dirName, hSigmaSignal1, 0, 1, -4, 4, pair, ptBins[i], ptBins[i + 1], "nSigma", nsigmaInFit);
//            pid1->peakFitResSub(dirName, hSigmaSignal1, 0, 1, -5, 5, pair, ptBins[i], ptBins[i + 1], "nSigma", nsigmaInFit);
        }
        binWidth[i] = (ptBins[i + 1] - ptBins[i]) / 2;
        xPt[i] = (ptBins[i + 1] + ptBins[i]) / 2;
        means[i] = pid1->getMean();
        meansE[i] = pid1->getMeanError();
        sigmas[i] = pid1->getSigma();
        sigmasE[i] = pid1->getSigmaError();

        Double_t integralClean = hSigmaSignal1->IntegralAndError(1, nBinsInFit, errorClean, ""); //number of it without PID cut
        Double_t integralAna = hSigmaSignal1->IntegralAndError(hSigmaSignal1->FindBin(-1*anaCut), hSigmaSignal1->FindBin(abs(anaCut)), errorAna, ""); //number of it without PID cut

        hTotalCount->SetBinContent(i+1, integralClean);
        hTotalCount->SetBinError(i+1, errorClean);
        hPassCount->SetBinContent(i+1, integralAna);
        hPassCount->SetBinError(i+1, errorAna);

        hTotalCountFunc->SetBinContent(i+1, pid1->getFuncIntegral());
        hTotalCountFunc->SetBinError(i+1, pid1->getFuncIntegralError());
        hPassCountFunc->SetBinContent(i+1, pid1->getFuncIntegral(-1*anaCut,abs(anaCut)));
        hPassCountFunc->SetBinError(i+1, pid1->getFuncIntegralError(-1*anaCut,abs(anaCut)));

        ++analysedBins;
        height=pid1->getHeight();
    }

    for (int i = 0; i < nBins-1; ++i) {
        cout<<hPassCount->GetBinCenter(i+1)<<" "<<hPassCount->GetBinContent(i+1)<<" "<<hTotalCount->GetBinContent(i+1)<<endl;
        if (hPassCount->GetBinContent(i+1)>=hTotalCount->GetBinContent(i+1)) hPassCount->SetBinContent(i+1, hTotalCount->GetBinContent(i+1)-1);
        if (hPassCountFunc->GetBinContent(i+1)>=hTotalCountFunc->GetBinContent(i+1)) hPassCountFunc->SetBinContent(i+1, hTotalCountFunc->GetBinContent(i+1)-1);
    }

    TGraphErrors *gMean = new TGraphErrors(analysedBins,xPt,means,binWidth, meansE);
    TGraphErrors *gSigmas = new TGraphErrors(analysedBins,xPt,sigmas,binWidth, sigmasE);

    TGraphAsymmErrors *gEff = new TGraphAsymmErrors();
    gEff->Divide(hPassCount,hTotalCount,"n");
    TGraphAsymmErrors *gEffFunc = new TGraphAsymmErrors();
    gEffFunc->Divide(hPassCountFunc,hTotalCountFunc,"n");

    TCanvas *c1 = new TCanvas("c1","%.3f_%.3f",1000,900);
    gMean->SetMarkerStyle(20);
    gMean->SetMarkerSize(0.9);
    gMean->SetMarkerColor(kBlack);
    gMean->SetLineColor(kBlack);
    if (tofPidEff) gMean->GetYaxis()->SetTitle("#pi #Delta1/#beta mean [#Delta1/#beta]");
    else gMean->GetYaxis()->SetTitle("#pi n#sigma mean [n#sigma]");
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
    if (tofPidEff) gSigmas->GetYaxis()->SetTitle("#pi #Delta1/#beta sigma [#Delta1/#beta]");
    else gSigmas->GetYaxis()->SetTitle("#pi n#sigma sigma [n#sigma]");
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
    gEff->GetYaxis()->SetRangeUser(0.6,1.1);
    gEff->SetTitle("");
    gEff->Draw("ap");
    c3->SaveAs(Form("%s/rootFiles/effTOFPion.root", dirName.Data()));
    c3->SaveAs(Form("%s/rootFiles/%s_effPion.png", dirName.Data(), effTypeString.Data()));
    c3->SaveAs(Form("%s/rootFiles/%s_effPion.eps", dirName.Data(), effTypeString.Data()));
    c3->Close();

    TCanvas *c4 = new TCanvas("c4","%.3f_%.3f",1000,900);
    gEffFunc->SetMarkerStyle(20);
    gEffFunc->SetMarkerSize(0.9);
    gEffFunc->SetMarkerColor(kBlack);
    gEffFunc->SetLineColor(kBlack);
    gEffFunc->GetYaxis()->SetTitle("Efficiency");
    gEffFunc->GetYaxis()->SetTitleOffset(1.1);
    gEffFunc->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    gEffFunc->GetYaxis()->SetRangeUser(0.6,1.1);
    gEffFunc->SetTitle("");
    gEffFunc->Draw("ap");
    c4->SaveAs(Form("%s/rootFiles/effTOFPionFunc.root", dirName.Data()));
    c4->SaveAs(Form("%s/rootFiles/%s_effPionFunc.png", dirName.Data(), effTypeString.Data()));
    c4->SaveAs(Form("%s/rootFiles/%s_effPionFunc.eps", dirName.Data(), effTypeString.Data()));
    c4->Close();

    TFile* resOut = new TFile(Form("%s/rootFiles/results_pipi.root", dirName.Data()),"RECREATE");
    gEffFunc->Write("effFunc");
    gEff->Write("eff");
    gMean->Write("mean");
    gSigmas->Write("sigma");
    resOut->Close();
    return dirName;
}