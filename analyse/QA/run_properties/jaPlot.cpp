//
// Created on 11.2.2018.
//
#include<iostream>
#include<fstream>
#include"TFile.h"
#include"TChain.h"
#include"TH1.h"
#include"TH2.h"
#include"TMath.h"
#include"TStyle.h"
#include"TPad.h"
#include"TCanvas.h"
#include"TObjArray.h"
#include"TLegend.h"
#include"TPaveText.h"
#include"TList.h"
#include"TAxis.h"
#include"TLatex.h"
#include"TF1.h"
#include"TLegend.h"
#include"TPaveText.h"

using namespace std;
using namespace TMath;

std::vector<int> RunNumberVector;
std::vector<int> BadRunNumberVector;

TFile *fileData = new TFile("qa_av.root", "READ");
TString textData = "dAu 200 GeV, Run16, trigger VPD-5, |V_{z}|< 6 cm, |V_{z,TPC}-V_{z,VPD}|< 6 cm, SL18f";

//TFile *fileData = new TFile("res_track_ntrkhftmin_event.root", "READ");
void make(TString inHist_name, TString histTitle, double nhitfit, double eta, double dca, double ptmin, double ptmax, TString yAxis, double rangeLow, double rangeUp, bool save, bool useForQA, bool makeAv, bool plotLeg, bool dcaetacut);

void jaPlot() {
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetTextFont(42) ;
    gStyle->SetLabelFont(42,"X");
    gStyle->SetLabelFont(42,"Y");
    gStyle->SetTitleFont(42,"X");
    gStyle->SetTitleFont(42,"Y");
    ifstream RunList("./ListOfRuns/Runnumber.list");
    if (RunList.is_open()) {
        int Run;
        while (RunList >> Run) {
            RunNumberVector.push_back(Run);
        }
    } else {
        cout << "Failed to open file!" << endl;
        return;
    }

    ofstream BadRunList("./BadRunList/BadRunList_MB_final.list");
    if(!BadRunList.is_open()) {
        cout<<"Failed to open file BadRunList_MB.list!"<<endl;
    }

    bool save = true;

//    z tuples
//    make("hZDCoverBBC_pr", "", 1, 0, 0.2, -1, "<ZDC/BBC rate>", 0., 0.23, save, true, true, false, false);
//    make("hNTOFtracks_pr", "", 1, 0, 0.2, -1, "<# TOF+TPC tracks>", 0., 12, save, true, true, true, true);
    make("h_mh1gRefmultCor", "", -1, 1, 0, 0.15, -1, "<Ref. multiplicity in #eta < |0.5|>", 0, 22, save, false, true, false,false);
    make("h_QA_nTracks", "", 15, 1, 1.5, 0.15, -1, "<# tracks>", 0., 35, save, true, true, true, true);
    make("h_QA_nTracks_HFT", "", 15, 1, 1.5, 0.15, -1, "<# HFT tracks>", 0., 20, save, false, false, true, true);
    make("h_QA_nTracks_TOF", "", 15, 1, 1.5, 0.15, -1, "<# TOF matched tracks>", 0., 20, save, true, true, true, true);
    make("h_QA_nTracks_HFT_TOF", "", 15, 1, 1.5, 0.15, -1, "<# HFT + TOF matched tracks>", 0., 20, save, false, false, true, true);


    make("h_QA_pT", "", 15, 1, 1.5, 0.15, -1, "TPC <#it{p_{T}}> [GeV/#it{c}]", 0.1, 0.65, save, true, true, true,true);
    make("h_QA_pT_HFT", "", 15, 1, 1.5, 0.15, -1, "HFT+TPC <#it{p_{T}}> [GeV/#it{c}]", 0.1, 0.65, save, true, true, true,true);
    make("h_QA_pT_HFT_TOF", "", 15, 1, 1.5, 0.15, -1, "HFT+TOF+TPC <#it{p_{T}}> [GeV/#it{c}]", 0.1, 0.65, save, true, true, true,true);
    make("h_QA_pT_TOF", "", 15, 1, 1.5, 0.15, -1, "TOF+TPC <#it{p_{T}}> [GeV/#it{c}]", 0.1, 0.65, save, true, true, true,true);

    make("h_QA_HFT_ratio", "", 15, 1, 1.5, 0.15, -1, "<HFT+TPC/TPC>", 0., 1, save, false, false, true,true);
    make("h_QA_HFT_ratio_0406", "", 15, 1, 1.5, 0.4, 0.6, "<HFT+TPC/TPC>", 0., 1, save, false, false, true,true);
    make("h_QA_HFT_ratio_1012", "", 15, 1, 1.5, 1, 1.2, "<HFT+TPC/TPC>", 0., 1, save, false, false, true,true);
    make("h_QA_HFT_ratio_3040", "", 15, 1, 1.5, 3, 4, "<HFT+TPC/TPC>", 0., 1, save, false, false, true,true);

    make("h_QA_ZDC_rate", "", -1, 1, 1.5, 0.15, -1, "<ZDCx> [kHz]", 0., 240, save, false, false, false, false);
    make("h_QA_BBC_rate", "", -1, 1, 0, 0.15, -1, "<BBCx> [kHz]", 0., 1400, save, false, false, false, false);
    make("h_QA_ZDC_over_BBC", "", -1, 1, 0, 0.15, -1, "<ZDCx/BBCx>", 0., 0.35, save, false, false, false, false);

    make("h_QA_Vz", "", -1, 1, 0, 0.15, -1, "V_{z} [cm]", -6, 6., save, false, true, false, false);
    make("h_QA_VzmVzVPD", "", -1, 1, 0, 0.15, -1, "|V_{z,VPD} - V_{z,TPC}| [cm]", 0., 6., save, false, true, false, false);

    make("h_QA_nHitsFit", "", -1, 1, 0, 0.15, -1, "<nHitsFit>", 0, 30, save, true, true, true, false);
    make("h_QA_nHitsFitTOF", "", 15, 1, 0, 0.15, -1, "TOF+TPC <N hits fit.>", 0., 45, save, true, true, true, false);
    make("h_QA_nHitsDedx", "", -1, 1, 0, 0.15, -1, "<nHitsDedx>", 0., 23, save, true, true, true, false);
    make("h_QA_nHFTHits", "", 15, 1, 1.5, 0.15, -1, "HFT+TPC # HFT hits", -0.1, 3, save, true, true, true, true);
    make("h_QA_nHFTHitsTOF", "", 15, 1, 1.5, 0.15, -1, "HFT+TOF+TPC # HFT hits", -0.1, 3, save, true, true, true, true);

    make("h_QA_DCA_z_zoom_HFT", "", 15, 1, 0.1, 0.15, -1, "HFT+TPC <DCA_{z}> [cm]", -0.003, 0.003, save, true, true, true, true);
    make("h_QA_DCA_xy_zoom_HFT", "", 15, 1, 0.1, 0.15, -1, "HFT+TPC <DCA_{xy}> [cm]", -0.003, 0.003, save, true, true, true, true);

    make("h_QA_nSigmaKaon", "", -1, 1, 0, 0.15, -1, "Kaon TPC <n#sigma> [-]", -2.4, -0.2, save, true, true, true, false);
    make("h_QA_nSigmaPion", "", -1, 1, 0, 0.15, -1, "Pion TPC <n#sigma> [-]", -1.3, 1.3, save, true, true, true, false);
    make("h_QA_OneOverBetaDiffKaon", "", 15, 1, 0, 0.2, -1, "Kaon <|1/#beta-1/#beta_{teo}|>", -0.3, 0.0, save, true, true, true, false);
    make("h_QA_OneOverBetaDiffPion", "", 15, 1, 0, 0.2, -1, "Pion <|1/#beta-1/#beta_{teo}|>", 0, 0.06, save, true, true, true, false);
    make("h_QA_Beta", "", 15, 1, 0, 0.15, -1, "Track TOF #beta [-]", 0.25, 1, save, true, true, true, false);
//
////    make("h_QA_DCA_TOF", "", 1, 1.5, 0.2, -1, "TOF+TPC <DCA> [cm]", -1, 1, save, true, true, true, true);
////    make("h_QA_DCA_TPC", "", 1, 1.5, 0.2, -1, "TPC <DCA> [cm]", -1, 1, save, true, true, true, true);
//    make("h_QA_DCA_z_HFT", "", 1, 1.5, 0.2, -1, "HFT+TPC <DCA_{z}> [cm]", -0.05, 0.01, save, true, true, true, true);
//    make("h_QA_DCA_xy_HFT", "", 1, 1.5, 0.2, -1, "HFT+TPC <DCA_{xy}> [cm]", -0.03, 0.005, save, true, true, true, true);
//    make("h_QA_DCA_z_zoom_HFT_TOF", "", 1.5, 0.1, 0.2, -1, "HFT+TOF+TPC <DCA> [cm]", -0.003, 0.003, save, true, true, true, true);
//    make("h_QA_DCA_xy_zoom_HFT_TOF", "", 1.5, 0.1, 0.2, -1, "HFT+TOF+TPC <DCA> [cm]", -0.003, 0.003, save, true, true, true, true);

    std::sort(BadRunNumberVector.begin(), BadRunNumberVector.end());
    for(unsigned int i = 0; i<BadRunNumberVector.size(); i++) {
        BadRunList<<BadRunNumberVector.at(i)<<endl;
    }
    cout<<BadRunNumberVector.size()<<endl;
}

void make(TString inHist_name, TString histTitle, double nHitsFit, double eta, double dca, double ptmin, double ptmax, TString yAxis, double rangeLow, double rangeUp, bool save, bool useForQA, bool makeAv, bool plotLeg, bool dcaetacut){
    TString canvasName = "c_"+inHist_name;
    TCanvas *out = new TCanvas(canvasName, canvasName, 1500, 750);
    TH1D *inHist = (TH1D*) fileData->Get(inHist_name);
//    inHist->SetTitle(histTitle);
    TH1D *inHist_copy = (TH1D*) inHist->Clone("inHist_copy");
    inHist_copy->SetMarkerStyle(2);
    inHist_copy->SetMarkerColor(kBlack);
    inHist_copy->GetXaxis()->SetTitle("RunIndex");
    inHist_copy->GetYaxis()->SetTitle(yAxis);
    inHist_copy->GetXaxis()->SetRange(2, RunNumberVector.size()+1);
    inHist_copy->GetYaxis()->SetRangeUser(rangeLow, rangeUp);
    inHist_copy->GetXaxis()->SetLabelFont(42);
    inHist_copy->GetXaxis()->SetTitleFont(42);

    TH1D *inHist_copy_draw = (TH1D*) inHist_copy->Clone("inHist_copy_draw");
//    Int_t centers[6] = {54,55,103,114,115,134};
    Double_t centers[6] = {54.5,55.5,103.5,114.5,115.5,134.5};
//
//    for (int m = 0; m < inHist_copy_draw->GetNbinsX(); ++m) {
//       for (int l = 0; l < 6; ++l) {
////            if  (m!=centers[l]) inHist_copy_draw->SetBinContent(m,-99);
////            if  (inHist_copy_draw->GetBinCenter(m)!=centers[l]) inHist_copy_draw->SetBinContent(m,-99);
//            if  (inHist_copy_draw->GetBinCenter(m)==centers[l]) inHist_copy->SetBinContent(m,-99); //ide
////            cout<<inHist_copy_draw->GetBinCenter(m)<<endl;
//        }
//    }
    inHist_copy_draw->SetMarkerStyle(2);
    inHist_copy_draw->SetMarkerColor(46);
    inHist_copy_draw->Draw("colz p e");

    inHist_copy->Draw("same");

    if (makeAv) {
        int nOut = 1;
        double av, err;

        while (nOut != 0) {
            int nFullBins = 0; //reset value to 0
            for (unsigned int i = 1; i < RunNumberVector.size()+1; i++) {//start at bin 2
                if (inHist->GetBinContent(i + 1) != 0) {
                    nFullBins += 1;
                }
            }

            av = inHist->Integral() / nFullBins;
            err = 0;
            for (unsigned int j = 1; j < RunNumberVector.size()+1; j++) {
                if (inHist->GetBinContent(j + 1) != 0) {
                    err += (inHist->GetBinContent(j + 1) - av) * (inHist->GetBinContent(j + 1) - av); //sum of squares of individual sub-errors
                }
            }

            err = TMath::Sqrt(err / ((nFullBins - 1)));
            nOut = 0; //reset nOut to zero
            if (useForQA) {
                for (unsigned int k = 1; k < RunNumberVector.size() + 1; k++) {

                    if ((TMath::Abs(inHist->GetBinContent(k + 1) - av) > 4 * err) && (inHist->GetBinContent(k + 1) != 0)) {
//                    if (((TMath::Abs(inHist->GetBinContent(k) - av) > 4 * err) && (inHist->GetBinContent(k) != 0)) || (inHist->GetBinContent(k))>950) {
                        nOut += 1;
//                        cout<<inHist_name<<"::"<<RunNumberVector.at(k - 1)<<"::"<< inHist->GetBinContent(k + 1)<<"::"<<k-1<<"::"<<endl;
                        cout<<inHist_name<<"::"<<RunNumberVector.at(k - 1)<<"::"<< inHist->GetBinContent(k)<<"::"<<k-1<<"::"<<endl;
                        inHist->SetBinContent(k + 1, 0);
                        std::sort(BadRunNumberVector.begin(), BadRunNumberVector.end());
                        if (!(std::binary_search(BadRunNumberVector.begin(), BadRunNumberVector.end(), RunNumberVector.at(k - 1)))) {
                            BadRunNumberVector.push_back(RunNumberVector.at(k - 1));
                        }
                    }
                }
            }
        }
        TF1 *av_func = new TF1("av_func", Form("%.8f", av), 0, RunNumberVector.at(RunNumberVector.size() - 1));
        av_func->SetLineColor(2);
        av_func->SetLineWidth(1);
        av_func->Draw("same");

        TF1 *err_low = new TF1("err_low", Form("%.8f", av - 4 * err), 0, RunNumberVector.at(RunNumberVector.size() - 1));
        err_low->SetLineColor(2);
        err_low->SetLineWidth(1);
        err_low->SetLineStyle(2);
        err_low->Draw("same");

        TF1 *err_up = new TF1("err_up", Form("%.8f", av + 4 * err), 0, RunNumberVector.at(RunNumberVector.size() - 1));
        err_up->SetLineColor(2);
        err_up->SetLineWidth(1);
        err_up->SetLineStyle(2);
        err_up->Draw("same");
        TLegend *legend = new TLegend(0.423, 0.15, 0.656, 0.302);
        legend->SetFillStyle(0);
        legend->AddEntry(av_func, "Average over runs");
        legend->AddEntry(err_low, "#pm4#sigma");
        legend->SetBorderSize(0);
        legend->Draw("same");
    }
    if(plotLeg) {
        TPaveText *text = new TPaveText(0.134269, 0.135881, 0.3837, 0.2887, "NDC");
        if (dcaetacut) text->AddText(Form("|#it{#eta}| < %.1f, |DCA| < %.1f cm", eta, dca));
        if (nHitsFit > 0) text->AddText(Form("NHitsFit > %.0f", nHitsFit));
        if (ptmax == -1) { text->AddText(Form("#it{p_{T}} > %.2f GeV/#it{c}", ptmin)); }
        else { text->AddText(Form("%.1f<#it{p_{T}}<%.1f GeV/#it{c}", ptmin, ptmax)); }
        text->SetFillColor(0);
        text->SetTextSize(0.04);
        text->Draw("same");
    }

    TPaveText *text1 = new TPaveText(0.4438, 0.91, 0.9166, 0.953, "NDC");
    text1->AddText(textData);
    text1->SetFillColor(0);
    text1->Draw("same");
    if (save) {
        TString folder = "./myFinal/";
        TString outName = folder + canvasName + ".png";
        out->SaveAs(outName);
    }

}