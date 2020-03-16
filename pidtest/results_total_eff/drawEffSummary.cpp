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

void drawEffSummary(){

    std::vector<TString> inputFiles;
    std::vector<TFile*> inputFilesF;
    std::vector<TGraphErrors*> graphs;
    std::vector<TF1*> fTGr;
    std::vector<TF1*> fTS;
//    std::vector<TString*> legendStrings;
    std::vector<const char *> legendStrings;
    Int_t colors[] = {1, 46, 8, 9};

    TString particle="K";

    inputFiles.push_back("bbc0_750_nHft0_nsigma3.0_tof0.03_pt0.0/totalEff_"+particle+".root");
    legendStrings.push_back("BBC: 0-750 kHz & n#sigma_{TPC}<3 & #Delta1/#beta<0.03");

    inputFiles.push_back("bbc750_950_nHft0_nsigma3.0_tof0.03_pt0.0/totalEff_"+particle+".root");
    legendStrings.push_back("BBC: 750-950 kHz & n#sigma_{TPC}<3 & #Delta1/#beta<0.03");

    inputFiles.push_back("bbc0_950_nHft0_nsigma3.0_tof0.03_pt0.0/totalEff_"+particle+".root");
    legendStrings.push_back("BBC: 0-950 kHz & n#sigma_{TPC}<3 & #Delta1/#beta<0.03");

    //    inputFiles.push_back("/home/lukas/work/dmesons/Dmaker_ndAu/Dmaker_dAu/pidtest/matchTOF_VzDiff6_bbc950_ntof0/rootFiles/results_KK.root");
//    inputFiles.push_back("/home/lukas/work/dmesons/Dmaker_ndAu/Dmaker_dAu/pidtest/hybridTOF_vzDiff6_bbc950_nTof0/rootFiles/results_KK.root");





    for (unsigned short i = 0; i < inputFiles.size(); ++i) {
        TFile *in = TFile::Open(inputFiles[i]);
        inputFilesF.push_back(in);
        TGraphErrors *gr = (TGraphErrors*) inputFilesF[i]->Get("grTotalGraphEffPid_"+particle);
        TF1 *fT = (TF1*) inputFilesF[i]->Get("fTotalGraphEffPid_"+particle);
        TF1 *fT1 = (TF1*) inputFilesF[i]->Get("fTotalEffPid_"+particle);
        graphs.push_back(gr);
        fTGr.push_back(fT);
        fTS.push_back(fT1);
    }
    TCanvas *out = new TCanvas("out", "out", 1200, 1000);
    TCanvas *outFTGR = new TCanvas("outFTGR", "outFTGR", 1200, 1000);
    TCanvas *outFT = new TCanvas("outFT", "outFT", 1200, 1000);

    TMultiGraph *mg = new TMultiGraph();
    mg->GetYaxis()->SetTitle("#varepsilon_{PID}");
    mg->GetXaxis()->SetTitle("p_{T} [GeV/c]");

    if (particle=="pi") mg->SetTitle("Pions");
    if (particle=="K") mg->SetTitle("Kaons");
    TLegend *legend = new TLegend(0.134, 0.16, 0.27, 0.34);
    legend -> SetFillStyle(0);
    legend -> SetLineColor(0);
    legend -> SetTextSize(0.035);
    cout<<graphs.size()<<endl;
    for (unsigned short j = 0; j < graphs.size(); ++j) {
        graphs[j]->SetMarkerColor(colors[j]);
        graphs[j]->SetLineColor(colors[j]);
        fTGr[j]->SetLineColor(colors[j]);
        fTS[j]->SetLineColor(colors[j]);
        graphs[j]->GetYaxis()->SetRangeUser(0,1.2);
        graphs[j]->SetName(Form("%i",j));
        legend -> AddEntry(graphs[j], legendStrings[j], "pl");
        mg->Add(graphs[j]);

        if (particle=="pi"){
            fTGr[j]->SetTitle("Pions");
            fTS[j]->SetTitle("Pions");
        }
        if (particle=="K"){
            fTGr[j]->SetTitle("Kaons");
            fTS[j]->SetTitle("Kaons");
        }

        outFTGR->cd();
        fTGr[j]->GetYaxis()->SetTitle("#varepsilon_{PID} fit total");
        fTGr[j]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
        if (j==0) fTGr[j]->Draw();
        else fTGr[j]->Draw("same");

        outFT->cd();
        fTS[j]->GetYaxis()->SetTitle("#varepsilon_{PID} fit separately");
        fTS[j]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
        if (j==0) fTS[j]->Draw();
        else fTS[j]->Draw("same");

    }
    outFTGR->cd();
    legend->Draw("same");
//    outFTGR->SaveAs

    outFT->cd();
    legend->Draw("same");

//    mg->SetMinimum(0);
    out->cd();
    mg->Draw("ap");

//    mg->GetYaxis()->SetRangeUser(0.,20.);

    legend->Draw("same");
}