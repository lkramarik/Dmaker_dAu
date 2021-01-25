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
#include <stdio.h>

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

    bool tpc= true;
    bool pions= true;

    if (tpc){
        if (pions){
            inputFiles.push_back("plotPart2/all/pi_tpc_bbc0_950_nHft1_nsigma100.0_tof100.00_pt0.0/rootFiles/results_pipi.root");
            legendStrings.push_back("part2");

            inputFiles.push_back("plotPart2/hft/pi_tpc_bbc0_950_nHft1_nsigma100.0_tof100000.00_pt0.0/rootFiles/results_pipi.root");
            legendStrings.push_back("part2 hft");
        } else { //KAONS
            inputFiles.push_back("plotPart2/all/tpc_bbc0_950_nHft1_nsigma100.0_tof100.00_pt0.0/rootFiles/results_KK.root");
            legendStrings.push_back("part2");

            inputFiles.push_back("plotPart2/hft/tpc_bbc0_950_nHft1_nsigma100.0_tof100.00_pt0.0/rootFiles/results_KK.root");
            legendStrings.push_back("part2 hft");
        }
    } else { //TOF
        if (pions){
            inputFiles.push_back("plotPart2/all/pi_tofPidEff_bbc0_950_nHft1_nsigma100.0_tof100.00_pt0.0/rootFiles/results_pipi.root");
            legendStrings.push_back("part2");

            inputFiles.push_back("plotPart2/hft/pi_tofPidEff_bbc0_950_nHft1_nsigma100.0_tof100000.00_pt0.0/rootFiles/results_pipi.root");
            legendStrings.push_back("part2 hft");
        } else { //KAONS
            inputFiles.push_back("plotPart2/all/tofPidEff_bbc0_950_nHft1_nsigma100.0_tof100.00_pt0.0/rootFiles/results_KK.root");
            legendStrings.push_back("part2");

            inputFiles.push_back("plotPart2/hft/tofPidEff_bbc0_950_nHft1_nsigma100.0_tof100.00_pt0.0/rootFiles/results_KK.root");
            legendStrings.push_back("part2 hft");
        }
    }
//    inputFiles.push_back("/home/lukas/work/dmesons/Dmaker_ndAu/Dmaker_dAu/pidtest/matchTof_vzDiff6_bbc950_nTof0/rootFiles/results_pipi.root");
//    inputFiles.push_back("/home/lukas/work/dmesons/Dmaker_ndAu/Dmaker_dAu/pidtest/matchTOF_VzDiff6_bbc950_ntof0/rootFiles/results_KK.root");
//    inputFiles.push_back("/home/lukas/work/dmesons/Dmaker_ndAu/Dmaker_dAu/pidtest/hybridTOF_vzDiff6_bbc950_nTof0/rootFiles/results_KK.root");
//    legendStrings.push_back("BBC rate < 950 kHz");

//    inputFiles.push_back("/home/lukas/work/dmesons/Dmaker_ndAu/Dmaker_dAu/pidtest/matchTOF_vzDiff6_bbc800_nTof0/rootFiles/results_pipi.root");
//    inputFiles.push_back("/home/lukas/work/dmesons/Dmaker_ndAu/Dmaker_dAu/pidtest/matchTOF_VzDiff6_bbc800_ntof0/rootFiles/results_KK.root");
//    inputFiles.push_back("/home/lukas/work/dmesons/Dmaker_ndAu/Dmaker_dAu/pidtest/hybridTOF_vzDiff6_bbc800_nTof0/rootFiles/results_KK.root");
//    legendStrings.push_back("BBC rate < 800 kHz");

//    inputFiles.push_back("/home/lukas/work/dmesons/Dmaker_ndAu/Dmaker_dAu/pidtest/bbc950_nHft0_nsigma1.0_tof0.03_pt0.5/rootFiles/results_KK.root");
//    inputFiles.push_back("/home/lukas/work/dmesons/Dmaker_ndAu/Dmaker_dAu/pidtest/bbc950_nHft0_nsigma1.0_tof0.03_pt0.5/rootFiles/results_pipi.root");
//    inputFiles.push_back("/home/lukas/work/dmesons/Dmaker_ndAu/Dmaker_dAu/pidtest/tpc_bbc950_nHft0_nsigma3.0_tof0.03_pt0.5/rootFiles/results_pipi.root");
//    legendStrings.push_back("n#sigma_{TPC}<3 & #Delta1/#beta<0.03");
//
//    inputFiles.push_back("/home/lukas/work/dmesons/Dmaker_ndAu/Dmaker_dAu/pidtest/tpc_bbc950_nHft0_nsigma3.0_tof0.03_pt0.5/rootFiles/results_KK.root");
//    legendStrings.push_back("no cuts on partner");



//    inputFiles.push_back("/home/lukas/work/dmesons/Dmaker_ndAu/Dmaker_dAu/pidtest/bbc950_nHft0_nsigma2.0_tof0.02_pt0.5/rootFiles/results_pipi.root");
//    inputFiles.push_back("/home/lukas/work/dmesons/Dmaker_ndAu/Dmaker_dAu/pidtest/bbc950_nHft0_nsigma2.0_tof0.02_pt0.5/rootFiles/results_KK.root");
//    legendStrings.push_back("n#sigma_{TPC}<2 & #Delta1/#beta<0.02");

//    inputFiles.push_back("hft/tpc_bbc0_950_nHft0_nsigma3.0_tof0.03_pt0.3/rootFiles/results_KK.root");
//    legendStrings.push_back("hft tpc_bbc0_950_nHft0_nsigma3.0_tof0.03_pt0.3");
////    legendStrings.push_back("BBC: 0-750 kHz & n#sigma_{TPC}<3 & #Delta1/#beta<0.03");
//
//    inputFiles.push_back("all/tpc_bbc0_950_nHft0_nsigma3.0_tof0.03_pt0.3/rootFiles/results_KK.root");
//    legendStrings.push_back("tpc_bbc0_950_nHft0_nsigma3.0_tof0.03_pt0.3");
////    legendStrings.push_back("BBC: 750-950 kHz & n#sigma_{TPC}<3 & #Delta1/#beta<0.03");
//
//    inputFiles.push_back("all/tpc_bbc0_950_nHft0_nsigma5.0_tof0.03_pt0.3/rootFiles/results_KK.root");
////    legendStrings.push_back("BBC: 0-950 kHz & n#sigma_{TPC}<3 & #Delta1/#beta<0.03");
//    legendStrings.push_back("tpc_bbc0_950_nHft0_nsigma5.0_tof0.03_pt0.3");
//
//    inputFiles.push_back("all/tpc_bbc0_950_nHft0_nsigma3.0_tof0.03_pt0.0/rootFiles/results_KK.root");
////    legendStrings.push_back("BBC: 0-950 kHz & n#sigma_{TPC}<3 & #Delta1/#beta<0.03");
//    legendStrings.push_back("tpc_bbc0_950_nHft0_nsigma3.0_tof0.03_pt0.0");
//
//    inputFiles.push_back("plotPart2/tpc_bbc0_950_nHft1_nsigma100.0_tof100.00_pt0.0/rootFiles/results_KK.root");
//////    legendStrings.push_back("BBC: 0-950 kHz & n#sigma_{TPC}<3 & #Delta1/#beta<0.03");
//    legendStrings.push_back("part2");



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
        min[0] = 0.5; min[1] = 0.5; min[2]=-0.8; min[3]=0.;
        max[0] = 1.2; max[1] = 1.2; max[2]=0.8; max[3]=1.5;
    } else {
        min[0] = 0.5; min[1] = 0.5; min[2]=-0.06; min[3]=0.;
        max[0] = 1.2; max[1] = 1.2; max[2]=0.06; max[3]=0.06;
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