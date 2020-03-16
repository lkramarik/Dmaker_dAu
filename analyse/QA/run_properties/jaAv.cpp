#include<iostream>
#include"TFile.h"
#include"TProfile.h"
#include"TChain.h"
#include"TH1.h"
#include"TH2.h"
#include"TMath.h"
#include"TStyle.h"
#include"TPad.h"
#include"TCanvas.h"
#include"TFractionFitter.h"
#include"TObjArray.h"
#include"TLegend.h"
#include"TPaveText.h"
#include"TList.h"
#include"TAxis.h"
#include <fstream>

using namespace std;
using namespace TMath;
void makeProfAv(TString histoIn, TString outName, TString yaxis);
void makeRatio(TString histoIn, TString outName, TString yaxis);

//TFile *fileData = new TFile("../Production/qa_honza_hitfit20_dca15.root", "READ");
//TFile *fileDataI = new TFile("../qa.2201.root", "READ");
TFile *fileDataI = new TFile("../qa.3001.root", "READ");
TFile *outFile = new TFile("qa_av.root", "RECREATE");

void jaAv(){

//makeProfAv("h_QA_pT", "pTperTrack",  "TPC <#it{p_{T}}> [GeV/#it{c}]");
//    makeProfAv("h_QA_pT", "pTperTrack",  "TPC <#it{p_{T}}> [GeV/#it{c}]");
//makeProfAv("pTperTrack",  "TPC <#it{p_{T}}> [GeV/#it{c}]");
    cout<<"ok, lets start with averaging"<<endl;

    makeProfAv("h_QA_ZDC_rate","h_QA_ZDC_rate","");
    makeProfAv("h_QA_BBC_rate","h_QA_BBC_rate","");
    makeProfAv("h_QA_ZDC_over_BBC","h_QA_ZDC_over_BBC","");
    makeProfAv("h_mh1gRefmultCor","h_mh1gRefmultCor","");

    makeProfAv("h_QA_nTracks_0406","h_QA_nTracks_0406","");
    makeProfAv("h_QA_nTracks_1012","h_QA_nTracks_1012","");
    makeProfAv("h_QA_nTracks_3040","h_QA_nTracks_3040","");
    makeProfAv("h_QA_nTracks","h_QA_nTracks","");

    makeProfAv("h_QA_nTracks_TPC","h_QA_nTracks_TPC","");
    makeProfAv("h_QA_nTracks_TOF","h_QA_nTracks_TOF","");

    makeProfAv("h_QA_nTracks_HFT_0406","h_QA_nTracks_HFT_0406","");
    makeProfAv("h_QA_nTracks_HFT_1012","h_QA_nTracks_HFT_1012","");
    makeProfAv("h_QA_nTracks_HFT_3040","h_QA_nTracks_HFT_3040","");
    makeProfAv("h_QA_nTracks_HFT","h_QA_nTracks_HFT","");

    makeProfAv("h_QA_nTracks_HFT_TOF","h_QA_nTracks_HFT_TOF","");

    makeProfAv("h_QA_Vz","h_QA_Vz","");
    makeProfAv("h_QA_VzmVzVPD","h_QA_VzmVzVPD","");

    makeProfAv("h_QA_pT","h_QA_pT","");
    makeProfAv("h_QA_pT_HFT","h_QA_pT_HFT","");
    makeProfAv("h_QA_pT_HFT_TOF","h_QA_pT_HFT_TOF","");
    makeProfAv("h_QA_pT_TOF","h_QA_pT_TOF","");

    makeProfAv("h_QA_nTracks_TOF","h_QA_nTracks_TOF","");
    makeProfAv("h_QA_nTracks","h_QA_nTracks_TPC","");
    makeProfAv("h_QA_nHitsFit","h_QA_nHitsFit","");
    makeProfAv("h_QA_nSigmaKaon","h_QA_nSigmaKaon","");
    makeProfAv("h_QA_nSigmaPion","h_QA_nSigmaPion","");
    makeProfAv("h_QA_nHitsDedx","h_QA_nHitsDedx","");
    makeProfAv("h_QA_OneOverBetaDiffKaon","h_QA_OneOverBetaDiffKaon","");
    makeProfAv("h_QA_OneOverBetaDiffPion","h_QA_OneOverBetaDiffPion","");
    makeProfAv("h_QA_Beta","h_QA_Beta","");
    makeProfAv("h_QA_nHitsFitTOF","h_QA_nHitsFitTOF","");
    makeProfAv("h_QA_nHFTHits","h_QA_nHFTHits","");
    makeProfAv("h_QA_nHFTHitsTOF","h_QA_nHFTHitsTOF","");
//    makeProfAv("h_QA_DCA_TOF","h_QA_DCA_TOF","");
//    makeProfAv("h_QA_DCA_HFT","h_QA_DCA_HFT","");
//    makeProfAv("h_QA_DCA_HFT_TOF","h_QA_DCA_HFT_TOF","");
//    makeProfAv("h_QA_DCA_TPC","h_QA_DCA_TPC","");
//
//
//
//    makeProfAv("h_QA_DCA_xy_HFT_TOF","h_QA_DCA_xy_HFT_TOF","");
//    makeProfAv("h_QA_DCA_xy_HFT","h_QA_DCA_xy_HFT","");
    makeProfAv("h_QA_DCA_xy_zoom_HFT","h_QA_DCA_xy_zoom_HFT","");
//    makeProfAv("h_QA_DCA_xy_zoom_HFT_TOF","h_QA_DCA_xy_HFT_TOF","");
//
//    makeProfAv("h_QA_DCA_xy_TPC","h_QA_DCA_xy_TPC","");
//    makeProfAv("h_QA_DCA_xy_zoom_TPC","h_QA_DCA_xy_zoom_TPC","");
//
//
//    makeProfAv("h_QA_DCA_z_TPC","h_QA_DCA_z_TPC","");
//    makeProfAv("h_QA_DCA_z_zoom_TPC","h_QA_DCA_z_zoom_TPC","");
//    makeProfAv("h_QA_DCA_z_zoom_HFT_TOF","h_QA_DCA_z_zoom_HFT_TOF","");
    makeProfAv("h_QA_DCA_z_zoom_HFT","h_QA_DCA_z_zoom_HFT","");
//    makeProfAv("h_QA_DCA_z_HFT_TOF","h_QA_DCA_z_HFT_TOF","");
//    makeProfAv("h_QA_DCA_z_HFT","h_QA_DCA_z_HFT","");
//    makeProfAv("h_QA_DCA_xy_TOF","h_QA_DCA_xy_TOF","");
//    makeProfAv("h_QA_DCA_z_TOF","h_QA_DCA_z_TOF","");
    cout<<"averaging done"<<endl;

    makeRatio("h_QA_nTracks","h_QA_nTracks_HFT","h_QA_HFT_ratio");
    makeRatio("h_QA_nTracks_0406","h_QA_nTracks_HFT_0406","h_QA_HFT_ratio_0406");
    makeRatio("h_QA_nTracks_1012","h_QA_nTracks_HFT_1012","h_QA_HFT_ratio_1012");
    makeRatio("h_QA_nTracks_3040","h_QA_nTracks_HFT_3040","h_QA_HFT_ratio_3040");
    outFile->Close();

}

void makeProfAv(TString histoIn, TString outName, TString yAxis){
//void makeProfAv(TString histoIn, TString outName, TString yAxis, Double_t ptmin, Double_t ptmax){
    TList *histoList = (TList*) fileDataI->Get("picoQAMaker");
    TH2F *h = static_cast<TH2F*>(histoList->FindObject(histoIn));
    h->Sumw2();
    TProfile *pxy=new TProfile();
    pxy = h->ProfileY("pxy");
    pxy->SetMarkerStyle(2);
    pxy->SetMarkerColor(kBlack);
    pxy->GetXaxis()->SetTitle("RunIndex");
    pxy->GetYaxis()->SetTitle(yAxis);
    pxy -> SetTitle("");

    pxy->GetXaxis()->SetLabelFont(42);
    pxy->GetXaxis()->SetTitleFont(42);
//    pxy->SetTitle(title[l]);
//    tempoName=title[l]+"_mean";

    outFile->cd();
    pxy->Write(outName);
}

void makeRatio(TString histoIn1, TString histoIn2, TString outName){
    TProfile *p1 = (TProfile*) outFile->Get(histoIn1);
    TProfile *p2 = (TProfile*) outFile->Get(histoIn2);
    TH1D *h1 = p1->ProjectionX("prX1");
    TH1D *h2 = p2->ProjectionX("prX2");
//    TProfile *h1 = static_cast<TProfile*>(histoList->FindObject(histoIn1));
//    TProfile *h2 = static_cast<TProfile*>(histoList->FindObject(histoIn2));
    h2->Divide(h1);
    outFile->cd();
    h2->Write(outName);




}