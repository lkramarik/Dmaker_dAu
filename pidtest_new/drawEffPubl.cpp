//
// Created by lukas on 29.1.2019.
//
#include<iostream>
#include<fstream>
#include "TGraphErrors.h"
#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "projectFiles.cpp"
#include <stdio.h>
// this is the best code for efficiency plotting
using namespace std;

void drawEffPubl(){

    std::vector<TString> inputFiles;
    std::vector<TFile*> inputFilesF;
    std::vector<TGraphErrors*> graphs[4];
    std::vector<TGraphErrors*> graphsEff;
    std::vector<TGraphErrors*> graphsEffFunc;
    std::vector<TGraphErrors*> graphsSigma;
    std::vector<TGraphErrors*> graphsMean;
//    std::vector<TString*> legendStrings;
    std::vector<const char *> legendStrings;
    Int_t colors[] = {1, 46, 8, 9, 38};

//    bool tpc= true;
    bool tpc= false;
//    bool pions= true;
    bool pions= false;
    bool compareNhft = false;


    if (!compareNhft) {
        if (tpc) {
            if (pions) {
                projectFiles("pipi", "plotPart2/all/pi_tpc_bbc0_950_nHft1_nsigma100.0_tof1000000.00_pt0.0/");
                inputFiles.push_back("plotPart2/all/pi_tpc_bbc0_950_nHft1_nsigma100.0_tof1000000.00_pt0.0/rootFiles/results_pipi.root");
                legendStrings.push_back("TPC");

                projectFiles("pipi", "plotPart2/hft/pi_tpc_bbc0_950_nHft1_nsigma100.0_tof1000000.00_pt0.0");
                inputFiles.push_back("plotPart2/hft/pi_tpc_bbc0_950_nHft1_nsigma100.0_tof1000000.00_pt0.0/rootFiles/results_pipi.root");
                legendStrings.push_back("TPC, with HFT");

                projectFiles("pipi", "plotPart2/all/pi_tpc_bbc0_950_nHft1_nsigma100.0_tof1000000.00_pt0.0/");
                inputFiles.push_back("plotPart2/all/pi_tpc_bbc0_950_nHft1_nsigma100.0_tof1000000.00_pt0.0/rootFiles/results_pipi.root");
                legendStrings.push_back("TPC, without HFT");

            } else { //KAONS
                inputFiles.push_back("plotPart2/all/tpc_bbc0_950_nHft1_nsigma100.0_tof1000000.00_pt0.0/rootFiles/results_KK.root");
                legendStrings.push_back("TPC");

                inputFiles.push_back("plotPart2/hft/tpc_bbc0_950_nHft1_nsigma100.0_tof1000000.00_pt0.0/rootFiles/results_KK.root");
                legendStrings.push_back("TPC, with HFT");

                inputFiles.push_back("plotPart2/nonhft/tpc_bbc0_950_nHft1_nsigma100.0_tof1000000.00_pt0.0/rootFiles/results_KK.root");
                legendStrings.push_back("TPC, without HFT");
            }
        } else { //TOF
            if (pions) {
                inputFiles.push_back("plotPart2/all/pi_tofPidEff_bbc0_950_nHft1_nsigma100.0_tof1000000.00_pt0.0/rootFiles/results_pipi.root");
                legendStrings.push_back("TPC");

                inputFiles.push_back("plotPart2/hft/pi_tofPidEff_bbc0_950_nHft1_nsigma100.0_tof1000000.00_pt0.0/rootFiles/results_pipi.root");
                legendStrings.push_back("TPC, with HFT");

                inputFiles.push_back("plotPart2/nonhft/pi_tofPidEff_bbc0_950_nHft1_nsigma100.0_tof1000000.00_pt0.0/rootFiles/results_pipi.root");
                legendStrings.push_back("TPC, without HFT");

            } else { //KAONS
                inputFiles.push_back("plotPart2/all/tofPidEff_bbc0_950_nHft1_nsigma100.0_tof1000000.00_pt0.0/rootFiles/results_KK.root");
                legendStrings.push_back("TPC");

                inputFiles.push_back("plotPart2/hft/tofPidEff_bbc0_950_nHft1_nsigma100.0_tof1000000.00_pt0.0/rootFiles/results_KK.root");
                legendStrings.push_back("TPC, with HFT");

                inputFiles.push_back("plotPart2/nonhft/tofPidEff_bbc0_950_nHft1_nsigma100.0_tof1000000.00_pt0.0/rootFiles/results_KK.root");
                legendStrings.push_back("TPC, without HFT");
            }
        }
    }

    if (compareNhft) {
        if (tpc) {
            if (pions) {
                inputFiles.push_back("plotPart2/all/pi_tpc_bbc0_950_nHft-100_nsigma100.0_tof1000000.00_pt0.0/rootFiles/results_pipi.root");
                legendStrings.push_back("no cut on N_{HFT tracks}");

                inputFiles.push_back("plotPart2/all/pi_tpc_bbc0_950_nHft1_nsigma100.0_tof1000000.00_pt0.0/rootFiles/results_pipi.root");
                legendStrings.push_back("N_{HFT tracks} > 1");

            } else { //KAONS
                inputFiles.push_back("plotPart2/all/tpc_bbc0_950_nHft-100_nsigma100.0_tof1000000.00_pt0.0/rootFiles/results_KK.root");
                legendStrings.push_back("no cut on N_{HFT tracks}");

                inputFiles.push_back("plotPart2/all/tpc_bbc0_950_nHft1_nsigma100.0_tof1000000.00_pt0.0/rootFiles/results_KK.root");
                legendStrings.push_back("N_{HFT tracks} > 1");
            }
        } else { //TOF
            if (pions) {
                inputFiles.push_back("plotPart2/all/pi_tofPidEff_bbc0_950_nHft1_nsigma100.0_tof1000000.00_pt0.0/rootFiles/results_pipi.root");
                legendStrings.push_back("TPC");

                inputFiles.push_back("plotPart2/all/pi_tofPidEff_bbc0_950_nHft-100_nsigma100.0_tof1000000.00_pt0.0/rootFiles/results_pipi.root");
                legendStrings.push_back("TPC, with HFT");

            } else { //KAONS
                inputFiles.push_back("plotPart2/all/tofPidEff_bbc0_950_nHft1_nsigma100.0_tof1000000.00_pt0.0/rootFiles/results_KK.root");
                legendStrings.push_back("TPC");

                inputFiles.push_back("plotPart2/all/tofPidEff_bbc0_950_nHft-100_nsigma100.0_tof1000000.00_pt0.0/rootFiles/results_KK.root");
                legendStrings.push_back("TPC, with HFT");
            }
        }
    }


    TString particleShort;
    if (pions) particleShort="#pi";
    else particleShort="K";

    TString grNames[] = {"eff", "effFunc", "mean", "sigma"};
    TString title[4];
    if (tpc) {
        title[0] = Form("#varepsilon_{PID}, TPC n#sigma_{%s}", particleShort.Data());
        title[1] = Form("#varepsilon_{PID}, TPC n#sigma_{%s}", particleShort.Data());
        title[2] = Form("n#sigma_{%s} mean [n#sigma_{%s}]", particleShort.Data(), particleShort.Data());
        title[3] = Form("n#sigma_{%s} sigma [n#sigma_{%s}]", particleShort.Data(), particleShort.Data());
    }
    else{
        title[0] = Form("#varepsilon_{PID}, TOF #delta1/#beta_{%s}", particleShort.Data());
        title[1] = Form("#varepsilon_{PID}, TOF #delta1/#beta_{%s}", particleShort.Data());
        title[2] = Form("#delta1/#beta_{%s} mean [#delta1/#beta_{%s}]", particleShort.Data(), particleShort.Data());
        title[3] = Form("#delta1/#beta_{%s} sigma [#delta1/#beta_{%s}]", particleShort.Data(), particleShort.Data());
    }

    TCanvas *out[4];
    TLegend *legend[4];
    Float_t min[4];
    Float_t max[4];
    if (tpc){
        min[0] = 0.7; min[1] = 0.7; min[2]=-1; min[3]=0.;
        max[0] = 1.1; max[1] = 1.1; max[2]=1; max[3]=1.5;
    } else {
        min[0] = 0.7; min[1] = 0.7; min[2]=-0.04; min[3]=0.;
        max[0] = 1.1; max[1] = 1.1; max[2]=0.04; max[3]=0.02;
    }
    TFile *in[5];


    for (unsigned short i = 0; i < inputFiles.size(); ++i) {
        in[i] = TFile::Open(inputFiles[i]);
        inputFilesF.push_back(in[i]);
        TGraphErrors *grEff = (TGraphErrors*) inputFilesF[i]->Get("eff");
        graphs[0].push_back(grEff);

        TGraphErrors *grMean = (TGraphErrors*) inputFilesF[i]->Get("effFunc");
        graphs[1].push_back(grMean);

        TGraphErrors *grSigma = (TGraphErrors*) inputFilesF[i]->Get("mean");
        graphs[2].push_back(grSigma);

        TGraphErrors *gr = (TGraphErrors*) inputFilesF[i]->Get("sigma");
        graphs[3].push_back(gr);
    }

    for (int i = 0; i < 4; ++i) {
        out[i] = new TCanvas(grNames[i], grNames[i], 1200, 1000);
        out[i]->SetGrid(1, 1);
        out[i]->SetLeftMargin(0.15);

        legend[i] = new TLegend(0.2, 0.133, 0.34, 0.311);
        legend[i] -> SetFillStyle(0);
        legend[i] -> SetLineColor(0);
        legend[i] -> SetTextSize(0.03);
        for (unsigned short j = 0; j < graphs[i].size(); ++j) {
            graphs[i][j]->SetMarkerColor(colors[j]);
            graphs[i][j]->SetLineColor(colors[j]);
            graphs[i][j]->SetName(Form("%i%i",j,i));
            graphs[i][j]->GetYaxis()->SetRangeUser(min[i],max[i]);
            graphs[i][j]->GetYaxis()->SetTitleOffset(1.7);
            graphs[i][j]->GetYaxis()->SetTitle(title[i]);
            graphs[i][j]->SetTitle("Kaons");
            if (pions) graphs[i][j]->SetTitle("Pions");
            if (j==0) graphs[i][j]->Draw("ap");
            graphs[i][j]->Draw("samep");
            legend[i] -> AddEntry(graphs[i][j], legendStrings[j], "pl");
        }

        legend[i]->Draw("same");
        if (tpc) {
            if (pions) out[i]->SaveAs("results/tpcPions" + grNames[i] + ".eps");
            else out[i]->SaveAs("results/tpcKaons" + grNames[i] + ".eps");
            if (pions) out[i]->SaveAs("results/tpcPions" + grNames[i] + ".png");
            else out[i]->SaveAs("results/tpcKaons" + grNames[i] + ".png");
        }
        else {
            if (pions) out[i]->SaveAs("results/tofPions" + grNames[i] + ".eps");
            else out[i]->SaveAs("results/tofKaons" + grNames[i] + ".eps");
            if (pions) out[i]->SaveAs("results/tofPions" + grNames[i] + ".png");
            else out[i]->SaveAs("results/tofKaons" + grNames[i] + ".png");
        }
    }

    // eff bin counts vs eff from fct fit.


}