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
    histo->SetLineWidth(3);
    histo->SetMarkerColor(color);
    histo->SetMarkerSize(2);

    histo->SetStats(0);
    histo->GetYaxis()->SetTitle(titleY);
    histo->GetYaxis()->SetTitleOffset(1.5);
    histo->GetYaxis()->SetLabelSize(0.045);
    histo->GetYaxis()->SetTitleSize(0.045);
    histo->GetYaxis()->SetLabelFont(42);
    histo->GetYaxis()->SetTitleFont(42);
    histo->GetYaxis()->CenterTitle(kTRUE);
//    histo->GetYaxis()->SetMaxDigits(2);

    histo->GetXaxis()->SetTitle(titleX);
    histo->GetXaxis()->SetLabelFont(42);
    histo->GetXaxis()->SetTitleFont(42);
    histo->GetXaxis()->SetTitleSize(0.045);
    histo->GetXaxis()->SetLabelSize(0.045);
    histo->GetXaxis()->CenterTitle(kTRUE);
}

//__________________________________________________________________________________________________________
void makeToD0(){
    gROOT->ProcessLine(".L toyMcClass.C++");

    TFile* sysFile = new TFile("tpc_sys.root","READ");
    TFile* sysFileOut = new TFile("tpc_sys_res.root","RECREATE");
    TString grNames[] = {"gr1_fun_dca", "gr1_fun_dca_hft", "gr1_fun_dca_all"};
    TString grTitle[] = {"TPC tracks w/o HFT hits", "HFT tracks", "TPC tracks"};
    Int_t col[]={kMagenta+2, kCyan+2, kGreen+2};
    Int_t lineStyle[]={9, 2, 1};

    TGraph* grSys[3];

    Float_t  ptMaxGraph=2.3;

    double_t ptBins[] = {1,2,3,5};
    const int nBins = sizeof(ptBins) / sizeof(double_t);

    TH1F* hAllD0 = new TH1F("hAllD0", "All D0", 50, 0., 5.);
    TH1F* hAllD0Pion = new TH1F("hAllD0Pion", "All D0 Pion", nBins-1, ptBins);
    TH1F* hAllD0Kaon = new TH1F("hAllD0Kaon", "All D0 Kaon", nBins-1, ptBins);

    TH1F* hAllD0Ana = new TH1F("hAllD0Ana", "Ana All D0", nBins-1, ptBins);
    TH1F* hAllPions = new TH1F("hAllPions", "All Pions", 50, 0.1, 5.);
    TH1F* hAllKaons = new TH1F("hAllKaons", "All Kaons", 50, 0.1, 5.);

    TH1F* hPassD0[3];
    TH1F* hPassD0Pion[3];
    TH1F* hPassD0Kaon[3];

    TH1F* hPassD0Ana[3];
    TH1F* hPassPions[3];
    TH1F* hPassKaons[3];

    for (int i = 0; i < 3; ++i) {
        grSys[i]=(TGraph*)sysFile->Get(grNames[i]);
        hPassD0Ana[i] = new TH1F("hAna"+grNames[i], "hAna"+grNames[i], nBins-1, ptBins);
        setHistoStyle(hPassD0Ana[i],col[i],20,"p_{T} [GeV/c]","TPC DCA embedding sys.");
        hPassD0Ana[i]->SetLineStyle(lineStyle[i]);

        hPassD0[i] = new TH1F("h"+grNames[i], "h"+grNames[i], 50, 0., 5.);
        setHistoStyle(hPassD0[i],1,20,"p_{T} [GeV/c]","");

        hPassD0Pion[i] = new TH1F("h_D0_pion"+grNames[i], "h"+grNames[i], nBins-1, ptBins);
        setHistoStyle(hPassD0Pion[i],kCyan+2,24,"p_{T} [GeV/c]","");
        hPassD0Pion[i]->SetLineStyle(2);

        hPassD0Kaon[i] = new TH1F("h_D0_kaon"+grNames[i], "h"+grNames[i], nBins-1, ptBins);
        setHistoStyle(hPassD0Kaon[i],kMagenta+2,25,"p_{T} [GeV/c]","");
        hPassD0Kaon[i]->SetLineStyle(2);

        hPassKaons[i] = new TH1F("hKaons"+grNames[i], "hKaons"+grNames[i], 50, 0.1, 5.);
        setHistoStyle(hPassKaons[i],kRed+2,20,"p_{T} [GeV/c]","");

        hPassPions[i] = new TH1F("hPions"+grNames[i], "hPions"+grNames[i], 50, 0.1, 5.);
        setHistoStyle(hPassPions[i],kBlue+2,20,"p_{T} [GeV/c]","");

    }


    ///////--------------------------------------------------------
    TFile* weightFile = new TFile("HIJING_D0_pt_y.root","READ");
    TH1F* hWeight = (TH1F*)weightFile->Get("hMcD0Pt");
    hWeight->Scale(1/hWeight->GetEntries());

    TFile* fMC = new TFile("/home/lukas/work/tmva_d0/sim/D0.toyMc.global.1.root","read");
    TNtuple* tMC = (TNtuple*) fMC->Get("nt");
    toyMcClass* t = new toyMcClass(tMC);


    TRandom3 *rndm1 = new TRandom3();
    TRandom3 *rndm2 = new TRandom3();

    long int nD0s = tMC->GetEntries();
    cout<<nD0s<<endl;
    Float_t ptTmpPi, ptTmpK,yPi,yK;
    for (int iD0 = 0; iD0 < nD0s/10; ++iD0) {
//    for (int iD0 = 0; iD0 < nD0s; ++iD0) {
        t->GetEntry(iD0);

        double ptD0 = t->rPt;

        Int_t bin = hWeight->FindBin(ptD0);
        double weightPt = hWeight->GetBinContent(bin);

        //pion:
        double piPtD=t->pRPt;
        if (piPtD < 0.15) continue;
        hAllPions->Fill(piPtD);
        ptTmpPi=piPtD;
        if (piPtD > ptMaxGraph)  ptTmpPi=ptMaxGraph;

        //kaon:
        double kPtD=t->kRPt;
        if (kPtD < 0.15) continue;
        hAllKaons->Fill(kPtD);
        ptTmpK=kPtD;
        if (kPtD > ptMaxGraph)  ptTmpK=ptMaxGraph;

        hAllD0->Fill(ptD0);
        hAllD0Ana->Fill(ptD0);
        hAllD0Pion->Fill(ptD0);
        hAllD0Kaon->Fill(ptD0);

        Double_t  piRan = rndm1->Uniform(2.)/2.;
        Double_t  kRan = rndm2->Uniform(2.)/2.;

        for (int i = 0; i < 3; ++i) {
            bool passPion = false;
            bool passKaon = false;

            yPi = grSys[i]->Eval(ptTmpPi, nullptr, "");
            if (yPi<piRan) {
                hPassPions[i]->Fill(piPtD);
                hPassD0Pion[i]->Fill(ptD0);
                passPion= true;
            }

            yK = grSys[i]->Eval(ptTmpK, nullptr, "");
            if (yK<kRan) {
                hPassKaons[i]->Fill(kPtD);
                hPassD0Kaon[i]->Fill(ptD0);
                passKaon= true;
            }

            if ((passKaon) || (passKaon)) {
//            if ((passPion)  (passPion)) {
                hPassD0[i]->Fill(ptD0);
                hPassD0Ana[i]->Fill(ptD0);
            }
        }
    }

    TCanvas* canD0peaks = new TCanvas("canD0peaks", "canD0peaks", 1500, 800);
    canD0peaks->Divide(3,1,1E-11, 1E-11);

    sysFileOut->cd();
    hAllD0->Write();
    hAllKaons->Write();
    hAllPions->Write();

    for (int i = 0; i < 3; ++i) {
        sysFileOut->cd();
        hPassD0[i]->Write();
        hPassPions[i]->Write();
        hPassKaons[i]->Write();

        hPassD0[i]->Divide(hAllD0);
        hPassPions[i]->Divide(hAllPions);
        hPassKaons[i]->Divide(hAllKaons);
        hPassD0Kaon[i]->Divide(hAllD0Kaon);
        hPassD0Pion[i]->Divide(hAllD0Pion);

        canD0peaks->cd(i+1);
        gPad->SetGrid();
        gPad->SetTopMargin(0.07);
        gPad->SetBottomMargin(0.13);

        hPassD0[i]->SetTitle(grTitle[i]);
        hPassD0[i]->GetXaxis()->CenterTitle();
        hPassD0[i]->GetYaxis()->CenterTitle();
        hPassD0[i]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
        hPassD0[i]->GetXaxis()->SetTitleSize(0.045);
        hPassD0[i]->GetYaxis()->SetTitleSize(0.045);
        hPassD0[i]->GetYaxis()->SetTitleOffset(1.7);
        hPassD0[i]->GetYaxis()->SetLabelSize(0.045);

        if (i==0){
            gPad->SetLeftMargin(0.15);
            gPad->SetRightMargin(0.0);
            hPassD0[i]->GetYaxis()->SetTitle("TPC DCA embedding sys.");
        }
        if (i>0){
            hPassD0[i]->GetXaxis()->SetRangeUser(0.1,4.9);
            gPad->SetLeftMargin(0.0);
            gPad->SetRightMargin(0.0);

        }
        if (i==2){
            gPad->SetLeftMargin(0.0);
            gPad->SetRightMargin(0.01);
        }

        hPassD0[i]->GetYaxis()->SetRangeUser(0.7,1.11);
        hPassD0[i]->GetXaxis()->SetRangeUser(0.,4.9);

        if (i>0){
            gPad->SetRightMargin(0.05);
            hPassD0[i]->GetXaxis()->SetRangeUser(0.1,4.9);
        }

        hPassD0[i]->Draw("hist C");
        hPassPions[i]->Draw("hist C same");
        hPassKaons[i]->Draw("hist C same");
        hPassD0Pion[i]->Draw("hist same");
        hPassD0Kaon[i]->Draw("hist same");

        grSys[i]->Draw("plsame");
    }

    TLegend *legend=new TLegend(0.57,0.13, 0.995, 0.33,"","brNDC");
    legend->SetFillStyle(0);
    legend->SetLineColor(0);
    legend->SetTextSize(0.05);
    legend->AddEntry(hPassD0[0], "D^{0}", "l");
    legend->AddEntry(hPassPions[0], "daughter #pi", "l");
    legend->AddEntry(hPassKaons[0], "daughter K", "l");
    legend->Draw("same");

    //////////////////////////////////////////////////////////////////////
    TCanvas* canD0Ana = new TCanvas("canD0Ana", "canD0Ana", 1000, 1000);
    gPad->SetGrid();
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.01);

    TLegend* legend2=new TLegend(0.57,0.13, 0.995, 0.33,"","brNDC");
    legend2->SetFillStyle(0);
    legend2->SetLineColor(0);
    legend2->SetTextSize(0.05);
    for (int i = 0; i < 3; ++i) {
        for (int j = 1; j < hPassD0Ana[i]->GetNbinsX()+1; ++j) {
            hPassD0Ana[i]->SetBinContent(j, hPassD0Kaon[i]->GetBinContent(j)+hPassD0Pion[i]->GetBinContent(j) );
            hPassD0Ana[i]->SetBinError(j,0.);

        }
        hPassD0Ana[i]->GetYaxis()->SetRangeUser(0,0.5);
        hPassD0Ana[i]->SetTitle("");

        if (i==0) hPassD0Ana[i]->Draw();
        else hPassD0Ana[i]->Draw("same");
        hPassD0Kaon[i]->Draw("same");
        hPassD0Pion[i]->Draw("same");

        sysFileOut->cd();
        hPassD0Ana[i]->Write();
    }
    legend2->Draw("same");

    //////////////////////////////////////////////////////////////////////
    TCanvas* canD0 = new TCanvas("canD0", "canD0", 1000, 1000);
    gPad->SetGrid();
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.01);


    for (int i = 0; i < 3; ++i) {
        hPassD0[i]->SetTitle("");
        hPassD0[i]->GetYaxis()->SetRangeUser(0,0.5);

        if (i==0) hPassD0[i]->Draw();
        else hPassD0[i]->Draw("same");
        hPassD0Kaon[i]->Draw("same");
        hPassD0Pion[i]->Draw("same");
    }
    legend2->Draw("same");


    Int_t colors[]={1, kGreen+2, kMagenta+2, kRed+2, kOrange+10, kYellow+1, kCyan+2, kBlue-4, 42, kRed-6};
    Int_t styles[]={20,21,22,23,24,25,26,27,28,29,30};
}


//__________________________________________
double getPtWeight(double pt){
    TFile* weightFile = new TFile("HIJING_D0_pt_y.root","READ");
    TH1F* hWeight = (TH1F*)weightFile->Get("hMcD0Pt");
    hWeight->Scale(1/hWeight->GetEntries());
    Int_t bin = hWeight->FindBin(pt);

    return hWeight->GetBinContent(bin);
}

//__________________________________________
//void make1minus(TH1)