//
// Created by lukas on 03.01.20.
// Test of Phi/K0s signal in and out of the hotspot
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


void makeOneCut(TString hotSpotFolder="", TCut hotspot="") {


    TString pair[] = {"KK", "#pi#pi"};
    TString pairName[] = {"KK", "pipi"};
    TString dirName = "hotspot";
    dirName+=hotSpotFolder;
    TString input[] = {"/home/lukas/work/dmesons/Dmaker_ndAu/Dmaker_dAu/ntp/ntp.KK.0301.hotspot.root", "/home/lukas/work/dmesons/Dmaker_ndAu/Dmaker_dAu/ntp/ntp.PiPi.0301.hotspot.root"};

    gSystem->Exec(Form("mkdir %s", dirName.Data()));
//    gSystem->Exec(Form("mkdir %s%s", dirName.Data(), hotSpotFolder.Data()));
    gSystem->Exec(Form("mkdir %s/rootFiles", dirName.Data()));
    gSystem->Exec(Form("mkdir %s/img", dirName.Data()));
    gSystem->Exec(Form("mkdir %s/img/KK", dirName.Data()));
    gSystem->Exec(Form("mkdir %s/img/KK/fit", dirName.Data()));

    gSystem->Exec(Form("mkdir %s", dirName.Data()));
    gSystem->Exec(Form("mkdir %s/rootFiles", dirName.Data()));
    gSystem->Exec(Form("mkdir %s/img", dirName.Data()));
    gSystem->Exec(Form("mkdir %s/img/pipi", dirName.Data()));
    gSystem->Exec(Form("mkdir %s/img/pipi/fit", dirName.Data()));


    Double_t massMean, massSigma;
    Float_t massMin, massMax, mean, sigma, ptPairMin, ptPairMax;
    TCut cut, cutPair;

    for (int i = 0; i < 2; ++i) {
        FitPID *fitmass = new FitPID();

        fitmass->setOutputFileName(Form("%s/rootFiles/mass_", dirName.Data()) + pairName[i] + ".root");
//        fitmass->setOutputFileName(Form("%s%s/rootFiles/mass_", dirName.Data(), hotSpotFolder.Data()) + pairName[i] + ".root");
        fitmass->setHeight(15000);

        if (i==0) {
            massMin = 1.0;// Phi to KK
            massMax = 1.04;// Phi to KK
            mean = 1.02;
            sigma = 0.002;
            ptPairMin = 0.2;
            ptPairMax = 10;

            cut = Form("pair_mass>%.3f && pair_mass<%.3f", massMin, massMax);
            cutPair=Form("pair_pt>%.3f && pair_pt<%.3f", ptPairMin, ptPairMax);
        }

        if (i==1) {
            massMin = 0.42;// K to pipi
            massMax = 0.58;// K to pipi
            mean = 0.5;
            sigma = 0.004;
            ptPairMin = 0.5;
            ptPairMax = 10;

            cut = Form("pair_mass>%.3f && pair_mass<%.3f", massMin, massMax);
            cutPair=Form("pair_pt>%.3f && pair_pt<%.3f", ptPairMin, ptPairMax);
        }

        TH1F *signal = (TH1F*) fitmass->projectSubtractBckg(dirName, input[i], 50, massMin, massMax, ptPairMin, ptPairMax, pair[i], cut + cutPair + hotspot, "pair_mass", "Mass_{%s} (GeV/c^{2})", true);
        fitmass->peakMassFit(dirName, signal, mean, sigma, massMin, massMax, pair[i], ptPairMin, ptPairMax, "mass");
        massMean = fitmass->getMean();
        massSigma = fitmass->getSigma();

        cout<<massMean<<" +- "<<massSigma<<endl;
        delete fitmass;
    }



}



void pidTest() {
    gROOT->ProcessLine(".L FitPID.c++");

    makeOneCut("", "");
    makeOneCut("/in", "hotSpot==1");
//    makeOneCut("/in", "hotSpot==1");

    TString pairName[] = {"KK", "pipi"};

    for (int i = 0; i < 2; ++i) {
        TFile *file = new TFile(Form("hotspot/rootFiles/mass_%s.root", pairName[i].Data()), "READ");
        TFile *fileIn = new TFile(Form("hotspot/in/rootFiles/mass_%s.root", pairName[i].Data()), "READ");

    }



}

