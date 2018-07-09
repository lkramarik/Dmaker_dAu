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
void compare_mix() {
    TString vars[9] = {"D_mass", "D_pt", "k_dca", "pi1_dca", "dcaDaughters", "D_decayL", "D_theta","k_pt","pi1_pt"};
    TString titles[9] = {"pair mass", "pair p_{T}", "kaon DCA", "pion DCA", "DCA of pion and kaon", "pair decay lenght", "pair pointing angle","kaon p_{T}","pion p_{T}" };
    TString units[9] = {" [GeV/c^{2}]", " [GeV/c]", " [cm]", " [cm]", " [cm]", " [cm]", " [-]", " [GeV/c]", " [GeV/c]"};
    Float_t limsMin[9] = {0.5, 0.9, 0.00, 0.00, 0.0, 0.0, 0,0,0};
    Float_t limsMax[9] = {4.0, 2.1, 0.04, 0.04, 0.016, 0.06, 1.7, 3.5, 3.5};
    TH1F* hVar[4][9] = {{new TH1F(),new TH1F(),new TH1F(),new TH1F(),new TH1F(),new TH1F(),new TH1F(),new TH1F(),new TH1F()},
                        {new TH1F(),new TH1F(),new TH1F(),new TH1F(),new TH1F(),new TH1F(),new TH1F(),new TH1F(),new TH1F()},
                        {new TH1F(),new TH1F(),new TH1F(),new TH1F(),new TH1F(),new TH1F(),new TH1F(),new TH1F(),new TH1F()},
                        {new TH1F(),new TH1F(),new TH1F(),new TH1F(),new TH1F(),new TH1F(),new TH1F(),new TH1F(),new TH1F()}};

    TCanvas *c[9] = {new TCanvas("c1","c1",900,1200),
                     new TCanvas("c2","c2",900,1200),
                     new TCanvas("c3","c3",900,1200),
                     new TCanvas("c4","c4",900,1200),
                     new TCanvas("c5","c5",900,1200),
                     new TCanvas("c6","c6",900,1200),
                     new TCanvas("c7","c7",900,1200),
                     new TCanvas("c8","c8",900,1200),
                     new TCanvas("c9","c9",900,1200)};

    Int_t colors[4] = {9, 46, 8, 1};
    TString ntpnames[4] = {"signal_SE", "signal_ME", "background_SE", "background_ME"};
    TNtuple *ntp[4];
    TFile* file[4];
    for (int i = 0; i < 4; ++i) {
        file[i] = new TFile("ntp_"+ntpnames[i] + ".all.root","r");
        ntp[i] = (TNtuple*)file[i] -> Get("ntp_"+ntpnames[i]);
    }

    TCut* detLu = new TCut("k_dca > 0.00 && pi1_dca > 0.00 && dcaDaughters < 0.02 && D_decayL > 0.001 && k_pt > 0.15 && pi1_pt > 0.15 && cos(D_theta) > 0");

    TFile *fOut = new TFile("res_compare.root","recreate");
    Double_t nentr, max;
    for (Int_t k = 0; k < 9 ; ++k) {
        c[k] -> cd();
        gPad->SetLeftMargin(0.15);
        TLegend *legend = new TLegend(0.6, 0.76, 0.75, 0.89);
        legend -> SetFillStyle(0);
        legend -> SetLineColor(0);
        legend -> SetTextSize(0.035);

        max = 0;
        for (int i = 0; i < 4; ++i) {
            const int j = k;
            hVar[i][k] = new TH1F(vars[k]+"_"+ntpnames[i], titles[k], 100, limsMin[k], limsMax[k]);
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

            max = (hVar[i][k]->GetMaximum() > max) ? hVar[i][k]->GetMaximum() : max;
            hVar[0][k]->SetAxisRange(0., 1.12*max,"Y");
            if (i==0) hVar[i][k] -> Draw();
            else hVar[i][k] -> Draw("same");
            legend -> AddEntry(hVar[i][k], ntpnames[i], "pl");
            fOut -> cd();
            hVar[i][k] -> Write();

        }
        legend->Draw();
        c[k] -> Modified();
        c[k] -> Update();

        c[k]->SaveAs(vars[k]+".png");
    }

    for (int l = 0; l < 4; ++l) {
        file[l] -> Close();

    }
    fOut -> Close();
}







