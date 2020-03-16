//
// Created by lukas on 17.12.19.
//
#include <iostream>
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TString.h"
#include "TROOT.h"
#include "TFile.h"
#include "TPaveText.h"
#include "TNtuple.h"
#include "TBox.h"


void vzVpdProject() {
    TString input = "ntp.PV.1712.root";
    TFile* data = new TFile(input ,"r");
    TNtuple* ntp = (TNtuple*)data -> Get("ntp_vertex");

    Double_t vzVpdCuts[] = {0, 1, 2, 3, 4, 5, 6};
    TH2F* hRanges[6];
    TH2F* hCut[7];

    TString name;
    Double_t axisMin=-1;
    Double_t axisMax=1;
    Double_t nBins = (abs(axisMin)+abs(axisMax))*200;

    for (int j = 0; j < 7; ++j) {
        name = Form("VzTPC_VPD_cut_%.0f", vzVpdCuts[j] );
        cout<<name<<endl;
        hCut[j] = new TH2F(name, "", nBins, axisMin, axisMax, nBins, axisMin, axisMax);
        hCut[j]->GetXaxis()->SetTitle("X [cm]");
        hCut[j]->GetYaxis()->SetTitle("Y [cm]");
//        hCut[j]->GetYaxis()->SetTitleOffset(0.9);
        hCut[j]->GetZaxis()->SetRangeUser(0., 3500000);
        hCut[j]->SetStats(0);
    }

    for (int k = 0; k < 6; ++k) {
        name = Form("VzTPC_VPD_range_%.0f_%.0f", vzVpdCuts[k], vzVpdCuts[k+1]);
        cout<<name<<endl;
        hRanges[k] = new TH2F(name, "", nBins, axisMin, axisMax, nBins, axisMin, axisMax);
        hRanges[k]->GetXaxis()->SetTitle("X [cm]");
        hRanges[k]->GetYaxis()->SetTitle("Y [cm]");
//        hRanges[k]->GetYaxis()->SetTitleOffset(0.9);
        hRanges[k]->GetZaxis()->SetRangeUser(0., 3500000);
        hRanges[k]->SetStats(0);

    }

    Float_t picoDstVz, VzVpd, picoDstVx, picoDstVy;
    ntp -> SetBranchAddress("picoDstVz",&picoDstVz);
    ntp -> SetBranchAddress("VzVpd", &VzVpd);
    ntp -> SetBranchAddress("picoDstVx", &picoDstVx);
    ntp -> SetBranchAddress("picoDstVy", &picoDstVy);

    cout<<ntp->GetEntries()<<endl;
    for (int i = 0; i < ntp->GetEntries(); ++i) {
//    for (int i = 0; i < 1000; ++i) {
        if (i%30000000==0) {cout<<i<<endl;}

        ntp -> GetEntry(i);
        for (int j = 0; j < 7; ++j) {
            if (abs(picoDstVz-VzVpd)<vzVpdCuts[j]) hCut[j]->Fill(picoDstVx,picoDstVy);
        }

        for (int k = 0; k < 6; ++k) {
            if (abs(picoDstVz-VzVpd)>=vzVpdCuts[k] && abs(picoDstVz-VzVpd)<vzVpdCuts[k+1]) hRanges[k]->Fill(picoDstVx,picoDstVy);
        }
    }

    TCanvas *c = new TCanvas("c","c",1000,900);
    gPad->SetLogz();
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.15);

    TBox *box = new TBox(-0.1, -0.1, -0.3, -0.3);
    box->SetFillColor(2);
    box->SetFillStyle(0);
    box->SetLineColor(2);
    box->SetLineWidth(2.);

    for (int l = 0; l < 7; ++l) {
        name = Form("VzTPC_VPD_cut_%.0f", vzVpdCuts[l] );
        hCut[l]->Draw("colz");

        TPaveText *text5 = new TPaveText(0.796, 0.90, 0.8667, 0.9395, "brNDC");
        text5->SetTextSize(0.02);
        text5->SetLineColor(0);
        text5->SetShadowColor(0);
        text5->SetFillColor(0);
        text5->SetTextFont(42);
        text5->AddText(Form("abs(V_{z,TPC}-V_{z,VPD}) < %.0f cm", vzVpdCuts[l]));
        text5->Draw("same");
        box->Draw("same");

        c->SaveAs("vzVpdStudy/"+name+".png");
        c->SaveAs("vzVpdStudy/"+name+".pdf");
    }

    for (int m = 0; m < 6; ++m) {
        name = Form("VzTPC_VPD_range_%.0f_%.0f", vzVpdCuts[m], vzVpdCuts[m+1] );
        hRanges[m]->Draw("colz");
        TPaveText *text5 = new TPaveText(0.796, 0.90, 0.8667, 0.9395, "brNDC");
        text5->SetTextSize(0.02);
        text5->SetLineColor(0);
        text5->SetShadowColor(0);
        text5->SetFillColor(0);
        text5->SetTextFont(42);
        text5->AddText(Form("%.0f cm < abs(V_{z,TPC}-V_{z,VPD}) < %.0f cm", vzVpdCuts[m], vzVpdCuts[m+1]));
        text5->Draw("same");
        box->Draw("same");

        c->SaveAs("vzVpdStudy/"+name+".png");
        c->SaveAs("vzVpdStudy/"+name+".pdf");


    }



}