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
TGraphErrors* histoToGraph(TH1F* h);

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
    histo->GetYaxis()->CenterTitle(kTRUE);
}

//__________________________________________________________________________________________________________
void makeToD0New() {
    gROOT->ProcessLine(".L toyMcClass.C++");

    TFile *sysFile = new TFile("tpc_sys.root", "READ");
    TFile *sysFileOut = new TFile("tpc_sys_res.root", "RECREATE");

    TString particles[] = {"kaon", "pion"};
    int iPart = 0;

    TString grNames[] = {"", "_hft", "_all"};

    TString var[] = {"dca", "nHitsFit"};
    Double_t maxEff[] = {0.3, 0.01};
    TGraphErrors *grD0Total[3][2];

    TString grTitle[] = {"TPC tracks w/o HFT hits", "HFT tracks", "TPC tracks"};
    Int_t col[] = {kMagenta + 2, kCyan + 2, kGreen + 2};
    Int_t lineStyle[] = {9, 2, 1};
    Int_t markerStyle[] = {20, 21, 22};

    TGraph *grSysKaon[3];
    TGraph *grSysPion[3];

    Float_t ptMaxGraph = 2.3;

    double_t ptBins[] = {1, 2, 3, 5};
    const int nBins = sizeof(ptBins) / sizeof(double_t);

    TH1F *hAllD0 = new TH1F("hAllD0", "All D0", 50, 0., 5.);
    TH1F *hAllD0Pion = new TH1F("hAllD0Pion", "All D0 Pion", nBins - 1, ptBins);
    TH1F *hAllD0Kaon = new TH1F("hAllD0Kaon", "All D0 Kaon", nBins - 1, ptBins);

    TH1F *hAllD0Ana = new TH1F("hAllD0Ana", "Ana All D0", nBins - 1, ptBins);
    TH1F *hAllPions = new TH1F("hAllPions", "All Pions", 50, 0.1, 5.);
    TH1F *hAllKaons = new TH1F("hAllKaons", "All Kaons", 50, 0.1, 5.);

    TH1F *hPassD0[3];
    TH1F *hPassD0Pion[3];
    TH1F *hPassD0Kaon[3];

    TH1F *hPassD0Ana[3];


    for (int iVar = 0; iVar < 2; ++iVar) {
        for (int i = 0; i < 3; ++i) {
            grSysKaon[i] = (TGraph *) sysFile->Get("kaon_gr1_fun_" + var[iVar] + grNames[i]);
            grSysPion[i] = (TGraph *) sysFile->Get("pion_gr1_fun_" + var[iVar] + grNames[i]);

            hPassD0Ana[i] = new TH1F("hAna" + grNames[i], "hAna" + grNames[i], nBins - 1, ptBins);
            setHistoStyle(hPassD0Ana[i], 1, 20, "D^{0} p_{T} [GeV/c]", "TPC DCA embedding sys.");
            hPassD0Ana[i]->SetLineStyle(lineStyle[i]);

            hPassD0[i] = new TH1F("h" + grNames[i], "h" + grNames[i], 50, 0., 5.);
            setHistoStyle(hPassD0[i], 1, 20, "D^{0} p_{T} [GeV/c]", "");

            hPassD0Pion[i] = new TH1F("h_D0_pion" + grNames[i], "h" + grNames[i], nBins - 1, ptBins);
            setHistoStyle(hPassD0Pion[i], kRed + 2, 21, "D^{0} p_{T} [GeV/c]", "");
            hPassD0Pion[i]->SetLineStyle(2);

            hPassD0Kaon[i] = new TH1F("h_D0_kaon" + grNames[i], "h" + grNames[i], nBins - 1, ptBins);
            setHistoStyle(hPassD0Kaon[i], kBlue + 2, 22, "D^{0} p_{T} [GeV/c]", "");
            hPassD0Kaon[i]->SetLineStyle(2);

        }


        ///////--------------------------------------------------------
        TFile *weightFile = new TFile("HIJING_D0_pt_y.root", "READ");
        TH1F *hWeight = (TH1F *) weightFile->Get("hMcD0Pt");
        hWeight->Scale(1 / hWeight->GetEntries());

        TFile *fMC = new TFile("/home/lukas/work/tmva_d0/sim/D0.toyMc.global.1.root", "read");
        TNtuple *tMC = (TNtuple *) fMC->Get("nt");
        toyMcClass *t = new toyMcClass(tMC);


        TRandom3 *rndm1 = new TRandom3();
        TRandom3 *rndm2 = new TRandom3();

        long int nD0s = tMC->GetEntries();
        cout << nD0s << endl;
        Float_t ptTmpPi, ptTmpK, yPi, yK;
        for (int iD0 = 0; iD0 < nD0s / 20; ++iD0) {
//    for (int iD0 = 0; iD0 < nD0s; ++iD0) {
            t->GetEntry(iD0);

            double ptD0 = t->rPt;

            Int_t bin = hWeight->FindBin(ptD0);
            double weightPt = hWeight->GetBinContent(bin);

            //pion:
            double piPtD = t->pRPt;
            if (piPtD < 0.15) continue;
            hAllPions->Fill(piPtD);
            ptTmpPi = piPtD;
            if (piPtD > ptMaxGraph) ptTmpPi = ptMaxGraph;

            //kaon:
            double kPtD = t->kRPt;
            if (kPtD < 0.15) continue;
            hAllKaons->Fill(kPtD);
            ptTmpK = kPtD;
            if (kPtD > ptMaxGraph) ptTmpK = ptMaxGraph;

            hAllD0->Fill(ptD0);
            hAllD0Ana->Fill(ptD0);
            hAllD0Pion->Fill(ptD0);
            hAllD0Kaon->Fill(ptD0);

            Double_t piRan = rndm1->Uniform(2.) / 2.;
            Double_t kRan = rndm2->Uniform(2.) / 2.;

            for (int i = 0; i < 3; ++i) {
                bool passPion = false;
                bool passKaon = false;

                yPi = grSysPion[i]->Eval(ptTmpPi, nullptr, "");
                if (yPi < piRan) {
                    hPassD0Pion[i]->Fill(ptD0);
                    passPion = true;
                }

                yK = grSysKaon[i]->Eval(ptTmpK, nullptr, "");
                if (yK < kRan) {
                    hPassD0Kaon[i]->Fill(ptD0);
                    passKaon = true;
                }
            }
        }


        sysFileOut->cd();
        hAllD0->Write();
        hAllKaons->Write();
        hAllPions->Write();

        TCanvas *canD0peaks = new TCanvas("canD0peaks", "canD0peaks", 1500, 800);
        canD0peaks->Divide(3, 1, 1E-11, 1E-11);

        for (int i = 0; i < 3; ++i) {

            //uncertainty calculation for D0 ( = pion + kaon)
            hPassD0Kaon[i]->Divide(hAllD0Kaon);
            hPassD0Pion[i]->Divide(hAllD0Pion);
            for (int j = 1; j < hPassD0Ana[i]->GetNbinsX() + 1; ++j) {
                hPassD0Ana[i]->SetBinContent(j, hPassD0Kaon[i]->GetBinContent(j) + hPassD0Pion[i]->GetBinContent(j));
                hPassD0Ana[i]->SetBinError(j, 0.);
            }
            grD0Total[i][iVar] = histoToGraph(hPassD0Ana[i]);

            canD0peaks->cd(i + 1);
            gPad->SetGrid();
            gPad->SetTopMargin(0.07);
            gPad->SetBottomMargin(0.13);

            grD0Total[i][iVar]->SetTitle(grTitle[i]);
            grD0Total[i][iVar]->GetXaxis()->SetRangeUser(1., 4.9);
//            grD0Total[i][iVar]->GetYaxis()->SetMaxDigits(2);

            if (i == 0) {
                gPad->SetLeftMargin(0.2);
                gPad->SetRightMargin(0.0);
                grD0Total[i][iVar]->GetYaxis()->SetTitle("TPC embedding sys. uncert.");
            }
            if (i > 0) {
                gPad->SetLeftMargin(0.0);
                gPad->SetRightMargin(0.0);

            }
            if (i == 2) {
                gPad->SetRightMargin(0.1);
            }

            grD0Total[i][iVar]->GetYaxis()->SetRangeUser(0., maxEff[iVar]);

            grD0Total[i][iVar]->Draw("ap");
            hPassD0Pion[i]->Draw("same");
            hPassD0Kaon[i]->Draw("same");
        }

        TLegend *legend = new TLegend(0.05, 0.70, 0.48, 0.9, "", "brNDC");
        legend->SetFillStyle(0);
        legend->SetLineColor(0);
        legend->SetTextSize(0.05);
        legend->AddEntry(grD0Total[0][0], "D^{0} = #pi + K uncert. ", "pl");
        legend->AddEntry(hPassD0Pion[0], " daughter #pi", "pl");
        legend->AddEntry(hPassD0Kaon[0], "daughter K", "pl");
        legend->Draw("same");


        canD0peaks->SaveAs("img/uncertainty_" + var[iVar] + ".png");

    }

    //////////////////////////////////////////////////////////////////////
    TGraphErrors *grTotalTotal[3];
    for (int i = 0; i < 3; ++i) {
        grTotalTotal[i] = new TGraphErrors();
        grTotalTotal[i]->SetMarkerColor(col[i]);
        grTotalTotal[i]->SetMarkerStyle(markerStyle[i]);
        grTotalTotal[i]->SetLineColor(col[i]);
        grTotalTotal[i]->SetMarkerSize(2);
        grTotalTotal[i]->SetLineWidth(3);

        grTotalTotal[i]->GetXaxis()->SetTitle("D^{0} p_{T} [GeV/c]");
        grTotalTotal[i]->GetXaxis()->SetTitleSize(0.04);
        grTotalTotal[i]->GetYaxis()->SetTitleSize(0.04);
        grTotalTotal[i]->GetXaxis()->SetTitleOffset(1.7);
        grTotalTotal[i]->GetYaxis()->SetTitle("TPC embedding sys. uncert.");

        grTotalTotal[i]->GetXaxis()->CenterTitle();
        grTotalTotal[i]->GetYaxis()->CenterTitle();
    }

    for (int k = 0; k < grD0Total[0][0]->GetN(); ++k) {
        for (int i = 0; i < 3; ++i) { //types of uncertainty
            double x1 = grD0Total[i][0]->GetY()[k];
            double x2 = grD0Total[i][1]->GetY()[k];

            x1=sqrt(x1*x1 + x2*x2);
            grTotalTotal[i]->SetPoint(k, grD0Total[i][0]->GetX()[k], x1);
            grTotalTotal[i]->SetPointError(k, grD0Total[i][0]->GetEX()[k], 0);
        }
    }

    TLegend* legend2=new TLegend(0.18,0.69, 0.60, 0.89,"","brNDC");
    legend2->SetFillStyle(0);
    legend2->SetLineColor(0);
    legend2->SetTextSize(0.035);

    TCanvas *canD0Ana = new TCanvas("canD0Ana", "canD0Ana", 1000, 1000);
    gPad->SetGrid();
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.01);

    for (int i = 0; i < 3; ++i) {
        grTotalTotal[i]->GetYaxis()->SetRangeUser(0,0.4);
        sysFileOut->cd();
        grTotalTotal[i]->Write("tpc_sys"+grNames[i]);

        canD0Ana->cd();
        if (i==0) grTotalTotal[i]->Draw("ap");
        else grTotalTotal[i]->Draw("same p");
        legend2->AddEntry(grTotalTotal[i], grTitle[i], "pl");

    }
    legend2->Draw("same");

    return;


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

//-------------------------------------------------------------------------
TGraphErrors* histoToGraph(TH1F* h) {
    const int nPoints = h->GetNbinsX();

    TGraphErrors* gr=new TGraphErrors();

    double content=0;
    double binCenter=-1;
    Double_t x, y, xE, yE;
    int nAddedPoints = 0;
    for (int i = 0; i < nPoints; ++i) {
        binCenter=h->GetBinCenter(i+1);
        content=h->GetBinContent(i+1);

        gr->SetPoint(nAddedPoints, binCenter, content);
        gr->SetPointError(nAddedPoints, h->GetBinWidth(i+1)/2., h->GetBinError(i+1));

//        gr->SetPointError(nAddedPoints, 0, 0);
        nAddedPoints++;
    }
    TString name=h->GetName();
    name="gr"+name;
    gr->SetName(name);

    gr->SetMarkerSize(h->GetMarkerSize());
    gr->SetMarkerColor(h->GetMarkerColor());
    gr->SetMarkerStyle(h->GetMarkerStyle());
    gr->SetLineColor(h->GetLineColor());
    gr->SetLineWidth(3);
    gr->SetTitle(h->GetTitle());
    gr->GetXaxis()->SetTitle(h->GetXaxis()->GetTitle());
    gr->GetXaxis()->CenterTitle();
    gr->GetYaxis()->CenterTitle();

    gr->GetXaxis()->SetLabelSize(0.05);
    gr->GetYaxis()->SetLabelSize(0.05);
    gr->GetXaxis()->SetTitleSize(0.05);
    gr->GetYaxis()->SetTitleSize(0.05);

    return gr;
}