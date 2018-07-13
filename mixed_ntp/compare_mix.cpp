//
// Created by lukas on 9.7.2018.
//
#include <string>
#include "TFile.h"
#include "TNtuple.h"
#include "TList.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TLegend.h"
#include<iostream>
using namespace std;

void drawVar(TString var, TString ntpnames[4], Float_t min, Float_t max);

void compare_mix() {
    TString vars[9] = {     "D_mass",       "D_pt",         "k_dca",    "pi1_dca",  "dcaDaughters",         "D_decayL",          "D_theta",             "k_pt",         "pi1_pt"};
    TString titles[9] = {   "pair mass",    "pair p_{T}",   "kaon DCA", "pion DCA", "DCA of pion and kaon", "pair decay lenght", "pair pointing angle", "kaon p_{T}",   "pion p_{T}" };
    TString units[9] = {    " [GeV/c^{2}]", " [GeV/c]",     " [cm]",    " [cm]",    " [cm]",                " [cm]",             " [-]",                " [GeV/c]",     " [GeV/c]"};
    Float_t limsMin[9] = {  0.5,            0.,             0.004,      0.004,      0.0,                    0.009,                 0,                     0,              0};
    Float_t limsMax[9] = {  3.0,            10,             0.04,       0.04,       0.016,                  0.06,                1.7,                   3.5,            3.5};

    TH1F* hVar[6][9] = {{new TH1F(),new TH1F(),new TH1F(),new TH1F(),new TH1F(),new TH1F(),new TH1F(),new TH1F(),new TH1F()},
                        {new TH1F(),new TH1F(),new TH1F(),new TH1F(),new TH1F(),new TH1F(),new TH1F(),new TH1F(),new TH1F()},
                        {new TH1F(),new TH1F(),new TH1F(),new TH1F(),new TH1F(),new TH1F(),new TH1F(),new TH1F(),new TH1F()},
                        {new TH1F(),new TH1F(),new TH1F(),new TH1F(),new TH1F(),new TH1F(),new TH1F(),new TH1F(),new TH1F()},
                        {new TH1F(),new TH1F(),new TH1F(),new TH1F(),new TH1F(),new TH1F(),new TH1F(),new TH1F(),new TH1F()},
                        {new TH1F(),new TH1F(),new TH1F(),new TH1F(),new TH1F(),new TH1F(),new TH1F(),new TH1F(),new TH1F()}};

//    TRatioPlot* rp[6][9] ={
//            {new TRatioPlot(), new TRatioPlot(), new TRatioPlot(), new TRatioPlot(), new TRatioPlot(), new TRatioPlot(), new TRatioPlot(), new TRatioPlot(), new TRatioPlot()},
//            {new TRatioPlot(), new TRatioPlot(), new TRatioPlot(), new TRatioPlot(), new TRatioPlot(), new TRatioPlot(), new TRatioPlot(), new TRatioPlot(), new TRatioPlot()},
//            {new TRatioPlot(), new TRatioPlot(), new TRatioPlot(), new TRatioPlot(), new TRatioPlot(), new TRatioPlot(), new TRatioPlot(), new TRatioPlot(), new TRatioPlot()},
//            {new TRatioPlot(), new TRatioPlot(), new TRatioPlot(), new TRatioPlot(), new TRatioPlot(), new TRatioPlot(), new TRatioPlot(), new TRatioPlot(), new TRatioPlot()},
//            {new TRatioPlot(), new TRatioPlot(), new TRatioPlot(), new TRatioPlot(), new TRatioPlot(), new TRatioPlot(), new TRatioPlot(), new TRatioPlot(), new TRatioPlot()},
//            {new TRatioPlot(), new TRatioPlot(), new TRatioPlot(), new TRatioPlot(), new TRatioPlot(), new TRatioPlot(), new TRatioPlot(), new TRatioPlot(), new TRatioPlot()}};

    TCanvas *c[9] = {new TCanvas("c1","c1",900,1200),
                     new TCanvas("c2","c2",900,1200),
                     new TCanvas("c3","c3",900,1200),
                     new TCanvas("c4","c4",900,1200),
                     new TCanvas("c5","c5",900,1200),
                     new TCanvas("c6","c6",900,1200),
                     new TCanvas("c7","c7",900,1200),
                     new TCanvas("c8","c8",900,1200),
                     new TCanvas("c9","c9",900,1200)};

    Int_t colors[6] = {9, 46, 8, 1, 28, 39};
    TString ntpnames[6] = {"background", "signal", "signal_SE", "signal_ME", "background_SE", "background_ME"};
    TNtuple *ntp[6];
    TFile* file[6];
    for (int i = 0; i < 6; ++i) {
        file[i] = new TFile("ntp_"+ntpnames[i] + ".all.root","r");
        ntp[i] = (TNtuple*)file[i] -> Get("ntp_"+ntpnames[i]);
    }

    //  published
    float decayL[] = {      0.012,  0.003,  0.009,  0.003};
    float dcaDaughters[] = {0.007,  0.016,  0.015,  0.013};
    float dcaD0[] = {       0.005,  0.0065, 0.0064, 0.0076};
    float cos[] = {         0.5,    0.5,   0.6,    0.5};
    float kdca[] = {        0.007,  0.01, 0.0076, 0.007};
    float pidca[] =         {0.009, 0.009, 0.0064, 0.0079}; //{0.009,0.009, 0.0079, 0.0079};

    TCut* detLu = new TCut(Form("D_mass > %1.2f && D_mass < %1.2f &&"
                                "D_pt > %1.2f && D_pt < %1.2f && "
                                "k_dca > %1.2f && k_dca < %1.2f && pi1_dca > %1.2f &&  pi1_dca < %1.2f && "
                           "dcaDaughters > %1.2f && dcaDaughters < %1.2f &&"
                           "D_decayL > %1.2f && D_decayL < %1.2f && "
                           "k_pt > 0.15 && pi1_pt > 0.15 && cos(D_theta) > 0.5",
                           limsMin[0], limsMax[0],
                           limsMin[1], limsMax[1],
                           limsMin[2], limsMax[2], limsMin[3], limsMax[3],
                                limsMin[4], limsMax[4],
                                limsMin[5], limsMax[5]));

    TFile *fOut = new TFile("res_compare.root","recreate");
    Double_t nentr, max;
    for (Int_t k = 0; k < 9 ; ++k) {
	cout<<vars[k]<<endl;
        c[k] -> cd();
        gPad->SetLeftMargin(0.15);
        TLegend *legend = new TLegend(0.6, 0.71, 0.75, 0.89);
//        TLegend *legend = new TLegend(0.6, 0.76, 0.75, 0.89); //for 4 stuff in there
        legend -> SetFillStyle(0);
        legend -> SetLineColor(0);
        legend -> SetTextSize(0.035);

        max = 0;
        for (int i = 0; i < 6; ++i) {
            const int j = k;
//            hVar[i][k] = new TH1F();
//            hVar[i][k] -> SetName(vars[k]+"_"+ntpnames[i]);
//            hVar[i][k] -> SetTitle(titles[k]);
            hVar[i][k] = new TH1F(vars[k]+"_"+ntpnames[i], titles[k], 200, limsMin[k], limsMax[k]);
            hVar[i][k] -> SetMarkerSize(0.4);
            hVar[i][k] -> SetLineColor(colors[i]);
            hVar[i][k] -> SetMarkerColor(46);
            hVar[i][k] -> SetLineWidth(3);
            hVar[i][k] -> GetXaxis()->SetTitle(vars[k]+units[k]);
            hVar[i][k] -> GetYaxis()->SetTitle("1/N_{entries}");
            hVar[i][k] -> GetYaxis() -> SetTitleOffset(2);
            ntp[i] -> Project(hVar[i][k]->GetName(), vars[k], *detLu);
            nentr = hVar[i][k] -> GetEntries();
            hVar[i][k]->Scale(1/nentr);
            hVar[i][k]->SetStats(0);

//            max = (hVar[i][k]->GetMaximum() > max) ? hVar[i][k]->GetMaximum() : max;
//            hVar[0][k]->SetAxisRange(0., 1.12*max,"Y");

//            if (i==0)   hVar[i][k] -> Draw();
//            else   hVar[i][k]->Draw("same");

            if (i==0) continue;
            if (i==1) {
                hVar[i][k] -> Divide(hVar[0][k]);
                hVar[i][k] -> Draw();
            }
            else {
                hVar[i][k] -> Divide(hVar[0][k]);
                hVar[i][k] -> Draw("same");
            }
            legend -> AddEntry(hVar[i][k], ntpnames[i], "pl");

            fOut -> cd();
            hVar[i][k] -> Write();

        }
        legend->Draw();
        c[k] -> Modified();
        c[k] -> Update();
//                     c[k]->SetTicks(0, 1);
//                     rp->Draw();
//                     rp->GetLowYaxis()->SetNdivisions(505);
//                     c[k]->Update();
        c[k]->SaveAs(vars[k]+".png");
    }

    for (int l = 0; l < 6; ++l) {
        file[l] -> Close();

    }
    fOut -> Close();
}

void drawVar(TString var, TString ntpName[4], Float_t min, Float_t max) {
    TFile *f = new TFile("res_compare.root","read");
    TCanvas *c = new TCanvas("c","c",900,1200);
    gPad->SetLeftMargin(0.15);
    TLegend *legend = new TLegend(0.6, 0.76, 0.75, 0.89);
    legend -> SetFillStyle(0);
    legend -> SetLineColor(0);
    legend -> SetTextSize(0.035);
    TH1F* hVar[4] = {new TH1F(),new TH1F(),new TH1F(),new TH1F()};
    for (int i = 0; i < 4; ++i) {
        hVar[i] = (TH1F*)f -> Get(var+"_"+ntpName[i]);
        hVar[i] -> SetAxisRange(min, max,"X");
        if (i==0) hVar[i] -> Draw();
        else hVar[i] -> Draw("same");
        legend -> AddEntry(hVar[i], ntpName[i], "pl");

    }
    legend->Draw();
    c -> Modified();
    c -> Update();
    c -> SaveAs(var+".png");

}

//
//if (i==1) {
//rp[i][k] = new TRatioPlot(hVar[i][k], hVar[0][k]);
//rp[i][k] -> Draw();
//}
//else {
//rp[i][k] = new TRatioPlot(hVar[i][k], hVar[0][k]);
//rp[i][k] -> Draw("same");
//