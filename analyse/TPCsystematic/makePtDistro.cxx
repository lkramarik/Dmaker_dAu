#include "TStyle.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TMath.h"
#include "TNtuple.h"
#include "TLine.h"
#include "TList.h"
#include "TLatex.h"
#include "TF1.h"
#include "TCanvas.h"
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
#include "TRandom.h"
#include "TRandom2.h"
#include "TRandom3.h"
#include "TROOT.h"
#include "toyMcClass.h"
#include "TLorentzVector.h"

double getPtWeight(double pt);
void setHistoStyle(TH1F* histo, Int_t color, Int_t marker, TString titleX, TString titleY);

//____________________________________________________________________________________________________________________________
void setHistoStyle(TH1F* histo, Int_t color, Int_t marker, TString titleX, TString titleY) {
    histo->Sumw2();

    histo->SetMarkerStyle(marker);
    histo->SetLineColor(color);
    histo->SetLineWidth(4);
    histo->SetMarkerColor(color);
    histo->SetMarkerSize(2);

    histo->SetStats(0);
    histo->GetYaxis()->SetTitle(titleY);
    histo->GetYaxis()->SetTitleOffset(1.5);
    histo->GetYaxis()->SetLabelSize(0.05);
    histo->GetYaxis()->SetTitleSize(0.05);
    histo->GetYaxis()->SetLabelFont(42);
    histo->GetYaxis()->SetTitleFont(42);
    histo->GetYaxis()->CenterTitle(kTRUE);
//    histo->GetYaxis()->SetMaxDigits(2);

    histo->GetXaxis()->SetTitle(titleX);
    histo->GetXaxis()->SetLabelFont(42);
    histo->GetXaxis()->SetTitleFont(42);
    histo->GetXaxis()->SetTitleSize(0.05);
    histo->GetXaxis()->SetLabelSize(0.05);
    histo->GetXaxis()->CenterTitle(kTRUE);
    histo->GetYaxis()->CenterTitle(kTRUE);
}

//__________________________________________________________________________________________________________
void makePtDistro(){
    gROOT->ProcessLine(".L toyMcClass.C++");

    TFile* sysFileOut = new TFile("pt_distro_daughters.root","RECREATE");
    Int_t col[]={kMagenta+2, kCyan+2, kGreen+2};

    double_t ptMin[] = {1,2,3};
    double_t ptMax[] = {2,3,5};
    const int nBins = sizeof(ptMin) / sizeof(double_t);

    TH1F* hPions[3];
    TH1F* hKaons[3];
    TH2F* hCorr[3];

    for (int i = 0; i < nBins; ++i) {
        hKaons[i] = new TH1F(Form("hKaons_pt_%.2f_%.2f", ptMin[i], ptMax[i]), Form("%.1f < D^{0} p_{T} < %.1f GeV/c", ptMin[i], ptMax[i]), 50, 0., 5.);
        setHistoStyle(hKaons[i],kRed+2,20,"p_{T} [GeV/c]","Scaled counts");

        hPions[i] = new TH1F(Form("hPions_pt_%.2f_%.2f", ptMin[i], ptMax[i]), Form("%.1f < D^{0} p_{T} < %.1f GeV/c", ptMin[i], ptMax[i]), 50, 0., 5.);
        setHistoStyle(hPions[i],kBlue+2,20,"p_{T} [GeV/c]","Scaled counts");

        hCorr[i] = new TH2F(Form("pions_vs_kaons_pt_%.2f_%.2f", ptMin[i], ptMax[i]), Form("%.1f < D^{0} p_{T} < %.1f GeV/c",ptMin[i], ptMax[i]),
                            80, 0., 5.,
                            80, 0., 5.);
        hCorr[i]->GetXaxis()->SetTitle("kaon p_{T} [GeV/c]");
        hCorr[i]->GetYaxis()->SetTitle("pion p_{T} [GeV/c]");
        hCorr[i]->SetStats(0);
        hCorr[i]->GetXaxis()->CenterTitle(kTRUE);
        hCorr[i]->GetYaxis()->CenterTitle(kTRUE);
        hCorr[i]->GetXaxis()->SetTitleSize(0.05);
        hCorr[i]->GetYaxis()->SetTitleSize(0.05);
        hCorr[i]->GetXaxis()->SetLabelSize(0.05);
        hCorr[i]->GetYaxis()->SetLabelSize(0.05);
    }

    TFile* weightFile = new TFile("HIJING_D0_pt_y.root","READ");
    TH1F* hWeight = (TH1F*)weightFile->Get("hMcD0Pt");
    hWeight->Scale(1/hWeight->GetEntries());

    TFile* fMC = new TFile("/home/lukas/work/tmva_d0/sim/D0.toyMc.global.1.root","read");
    TNtuple* tMC = (TNtuple*) fMC->Get("nt");
    toyMcClass* t = new toyMcClass(tMC);

    long int nD0s = tMC->GetEntries();
    cout<<nD0s<<endl;
    Float_t ptTmpPi, ptTmpK,yPi,yK;

//    for (int iD0 = 0; iD0 < nD0s/20; ++iD0) {
//    for (int iD0 = 0; iD0 < nD0s/10; ++iD0) {
    for (int iD0 = 0; iD0 < nD0s; ++iD0) {
        t->GetEntry(iD0);
        double ptD0 = t->rPt;

        Int_t bin = hWeight->FindBin(ptD0);
        double weightPt = hWeight->GetBinContent(bin);

        double piPtD = t->pRPt;
        double kPtD = t->kRPt;

        for (int i = 0; i < nBins; ++i) {
            if (ptD0 >= ptMin[i] && ptD0 < ptMax[i]) {
                hKaons[i]->Fill(kPtD, weightPt);
                hPions[i]->Fill(piPtD, weightPt);
                hCorr[i]->Fill(kPtD, piPtD);
            }
        }
    }

    TCanvas* canD0peaks = new TCanvas("D0_daughters_pT", "D0_daughters_pT", 1800, 900);
    canD0peaks->Divide(nBins,1, 1E-11, 1E-11);

    TCanvas* canD0Corr = new TCanvas("D0_daughters_pT_Corr", "D0_daughters_pT_Corr", 1800, 900);
    canD0Corr->Divide(nBins,1, 1E-11, 1E-11);

    sysFileOut->cd();
    double max=0;
    for (int i = 0; i < nBins; ++i) {
        hKaons[i]->Scale(1/hKaons[i]->GetMaximum());
        hPions[i]->Scale(1/hPions[i]->GetMaximum());

        if (max < hKaons[i]->GetMaximum()) max = hKaons[i]->GetMaximum();
        if (max < hPions[i]->GetMaximum()) max = hPions[i]->GetMaximum();

        hKaons[i]->Write();
        hPions[i]->Write();
    }
    max*=1.1;

    for (int i = 0; i < nBins; ++i) {
        canD0peaks->cd(i+1);
        gPad->SetTopMargin(0.07);
        gPad->SetBottomMargin(0.13);
        gPad->SetGrid();

        hKaons[i]->GetYaxis()->SetRangeUser(0.,max);
        hKaons[i]->GetXaxis()->SetRangeUser(0.,4.9);

        if (i==0){
            gPad->SetLeftMargin(0.15);
            gPad->SetRightMargin(0.0);
        }

        if (i>0){
            hKaons[i]->GetXaxis()->SetRangeUser(0.1,4.9);
            gPad->SetLeftMargin(0.0);
            gPad->SetRightMargin(0.0);

        }

        if (i==nBins-1){
            gPad->SetRightMargin(0.05);
        }

        hKaons[i]->Draw("hist C");
        hPions[i]->Draw("same hist C");

        ///////////----------------
        canD0Corr->cd(i+1);
        gStyle->SetPalette(53);
        TColor::InvertPalette();

        gPad->SetTopMargin(0.07);
        gPad->SetBottomMargin(0.13);
        hCorr[i]->GetXaxis()->SetRangeUser(0.,4.9);
        hCorr[i]->Scale(1./hCorr[i]->GetMaximum());
        if (i==0){
            gPad->SetLeftMargin(0.15);
            gPad->SetRightMargin(0.0);
        }

        if (i>0){
            hCorr[i]->GetXaxis()->SetRangeUser(0.1,4.9);
            gPad->SetLeftMargin(0.0);
            gPad->SetRightMargin(0.0);

        }

        if (i==nBins-1){
            gPad->SetRightMargin(0.12);
        }
        hCorr[i]->Draw("colz");

    }

    canD0peaks->cd();
    TLegend *legend=new TLegend(0.082,0.21, 0.50, 0.31,"","brNDC");
    legend->SetFillStyle(0);
    legend->SetLineColor(0);
    legend->SetTextSize(0.05);
    legend->AddEntry(hKaons[0], "kaons", "l");
    legend->AddEntry(hPions[0], "pions", "l");
    legend->Draw("same");

    canD0peaks->SaveAs(".png");
    canD0peaks->SaveAs(".eps");
    canD0peaks->SaveAs(".pdf");

    canD0Corr->SaveAs(".png");
    canD0Corr->SaveAs(".eps");
    canD0Corr->SaveAs(".pdf");
}


//__________________________________________
double getPtWeight(double pt){
    TFile* weightFile = new TFile("HIJING_D0_pt_y.root","READ");
    TH1F* hWeight = (TH1F*)weightFile->Get("hMcD0Pt");
    hWeight->Scale(1/hWeight->GetEntries());
    Int_t bin = hWeight->FindBin(pt);

    return hWeight->GetBinContent(bin);
}
