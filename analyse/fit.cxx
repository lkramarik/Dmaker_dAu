//D0 mass = 1.864 84

#include <iostream>
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TString.h"
#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TNtuple.h"

void fit(TString input){    
    TFile* data = new TFile(input ,"r");
    ntp = (TNtuple*)data -> Get("ntp;1");
    Float_t flag, D_mass, D_pt, k_pt, pi1_pt, pi1_dca, k_dca;
    ntp -> SetBranchAddress("flag",&flag);
    ntp -> SetBranchAddress("D_mass", &D_mass);
    //ntp -> SetBranchAddress("", &);
    ntp -> SetBranchAddress("D_pt", &D_pt);
    ntp -> SetBranchAddress("pi1_pt", &pi1_pt);
    ntp -> SetBranchAddress("k_pt", &k_pt);
    ntp -> SetBranchAddress("pi1_dca", &pi1_dca);
    ntp -> SetBranchAddress("k_dca", &k_dca);
    
    TH1F* hInvMassBackMin = new TH1F("background minus", "background minus", 2000, 0.6, 2.6);
    TH1F* hInvMassBackPlus = new TH1F("background plus", "background plus", 2000, 0.6, 2.6);
    TH1F* hInvMassSign = new TH1F("signal", "signal", 2000, 0.6, 2.6);    
    TH1F* hInvMassBack = new TH1F("background", "background", 2000, 0.6, 2.6); 
    
    hInvMassBackMin -> Sumw2();
    hInvMassBackPlus -> Sumw2();
    hInvMassSign -> Sumw2();
    
    Long64_t numberEntr = ntp -> GetEntries();
    cout<<"Number of entries in Ntuple: "<<numberEntr<<endl;;
    for (Long64_t i = 0; i < numberEntr; i++) {
        if (i%10000000==0) {cout<<i<<endl;}
        ntp -> GetEntry(i);        
        //         if ((pi1_dca > 0.008) && (k_dca > 0.007)){
        if ((pi1_dca > 0.008)){
            if ((D_pt > 1) && (D_pt < 6)) {            
                if ((flag >= 0 ) && (flag < 2)) {hInvMassSign -> Fill(D_mass); }
                else if (flag == 4) {hInvMassBackMin -> Fill(D_mass); }
                else {hInvMassBackPlus -> Fill(D_mass); }                
            }
        }
    }

    hInvMassBack -> Add(hInvMassBackPlus,1);

    Float_t value, error, valueM, errorM, valueP, errorP;
    for (Int_t j = 0; j < Nentr; j++) {
        valueP = hInvMassBack -> GetBinContent(j);  
        errorP = hInvMassBack -> GetBinError(j);
        valueM = hInvMassBackMin -> GetBinContent(j);  
        errorM = hInvMassBackMin -> GetBinError(j);
        error = 0.5*sqrt(valueM*errorP*errorP/valueP + valueP*errorM*errorM/valueM);
        hInvMassBack -> SetBinContent(j, sqrt(valueP*valueM) );
        hInvMassBack -> SetBinError(j, error);
    }
    TFile* dataRes = new TFile("res_"+input ,"RECREATE");
    hInvMassSign -> Write();
    hInvMassBack -> Write(); 
    cout<<"res_"+input<<endl;
    cout<<"done"<<endl;
}

void fitting() {
    Float_t fitRMin = 1.65;
    Float_t fitRMax = 2.35;
    int rebin = 10;
    TString input = "res_ntp_all.root";
    TFile* data = new TFile(input ,"r");
    hInvMassBack = (TH1F*)data -> Get("background");
    hInvMassSign = (TH1F*)data -> Get("signal");
   
    hInvMassSign -> Rebin(rebin);
    hInvMassBack -> Rebin(rebin);
    
    hInvMassSign -> SetMarkerColor(46);
    hInvMassSign -> SetLineColor(46);
    
    Double_t hBackIntegral = hInvMassBack -> Integral(hInvMassBack -> FindBin(1.60), hInvMassBack -> FindBin(1.7), "") + hInvMassBack -> Integral(hInvMassBack -> FindBin(1.95), hInvMassBack -> FindBin(2.5), "");
    Double_t hSignIntegral = hInvMassSign -> Integral(hInvMassSign -> FindBin(1.60), hInvMassSign -> FindBin(1.7), "") + hInvMassSign -> Integral(hInvMassSign -> FindBin(1.95), hInvMassSign -> FindBin(2.5), "");
    
    //     Double_t hBackIntegral = hInvMassBack -> Integral(hInvMassBack -> FindBin(2.06), hInvMassBack -> FindBin(2.2), "");// + hInvMassBack -> Integral(hInvMassBack -> FindBin(1.9), hInvMassBack -> FindBin(2), "");
    //     Double_t hSignIntegral = hInvMassSign -> Integral(hInvMassSign -> FindBin(2.06), hInvMassSign -> FindBin(2.2), "");// + hInvMassSign -> Integral(hInvMassSign -> FindBin(1.9), hInvMassSign -> FindBin(2), "");
    
    cout<<hInvMassBack -> GetNbinsX()<<endl;
    cout<<hInvMassSign -> GetNbinsX() << endl;

    cout<<"Backgroung Integral: "<<hBackIntegral<<endl;
    cout<<"Signal Integral :"<<hSignIntegral<<endl;
    
    const int Nbins = 1000/rebin;
    float err[Nbins], errS[Nbins];
    for (int j=0; j<Nbins; j++) {
        err[j] = hInvMassBack -> GetBinError(j);
        errS[j] = hInvMassSign -> GetBinError(j);
    }
    hInvMassBack -> Scale(hSignIntegral/hBackIntegral);
    
    for (j=0; j<Nbins; j++) {
        err[j] = err[j]*hSignIntegral/hBackIntegral;
        hInvMassBack -> SetBinError(j, err[j]);        
    }   //scaling error
       
    hInvMassSign -> Add(hInvMassBack,-1);
    for (j=0; j<Nbins; j++) {        
        hInvMassSign -> SetBinError(j, sqrt(err[j]*err[j] + errS[j]*errS[j]));
    }
    
    TCanvas *c3 = new TCanvas("c3","c3",1200,800);
    TF1 *fun0 = new TF1("fun0","pol1(0)+gaus(2)",fitRMin,fitRMax);
    fun0->SetParameters(1.,1.,1.,1.865,0.015);
    fun0->SetLineColor(2);
    fun0->SetLineStyle(7);
    
    const float MKSize    = 1.;
    const float rotwthmin = 1.80; // peak mean fitting range
    const float rotwthmax = 1.98; //peak mean fitting range
    double mm[Nbins],ym[Nbins],yme[Nbins],ym1[Nbins];
    
    for (int ib=0; ib<Nbins; ib++) {
        mm[ib]  = hInvMassSign -> GetBinCenter(ib);
        ym[ib]  = hInvMassSign -> GetBinContent(ib);
        yme[ib] = hInvMassSign -> GetBinError(ib);
    }
    
    TGraphErrors *gm = new TGraphErrors(Nbins,mm,ym,0,yme);
    gm->SetMarkerStyle(20);
    gm->SetMarkerSize(0.02);
    gm->SetMarkerColor(4);
    gm->SetLineColor(4);
    //gm->Draw("psame");
    //gm->Draw("apsame");
    //gm->Draw("AP");
    fun0->SetParLimits(3,rotwthmin,rotwthmax);
    //fun0->SetParLimits(4,0.010,0.014);
    //fun0->SetParLimits(4,0.007,0.020);
    //     fun0->SetParLimits(4,0.007,0.015);
    gm->Fit(fun0,"NOR");
    
    TF1 *resfunm = new TF1("resfunm","pol1",fitRMin,fitRMax);
    resfunm->SetParameters(fun0->GetParameter(0),fun0->GetParameter(1));
    resfunm->SetLineStyle(7);
    resfunm->SetLineWidth(1);
    //resfunm->Draw("same");
    //fun0->Draw("same");
    
    TF1 *resfunm1 = new TF1("resfunm1","pol1(0)+gaus(2)",fitRMin,fitRMax);
    resfunm1->SetParameters(0.,0.,fun0->GetParameter(2),fun0->GetParameter(3),fun0->GetParameter(4));
    resfunm1->SetLineColor(4);
    resfunm1->SetLineStyle(7);
    //resfunm1->Draw("same");
    
    for(j=0; j<Nbins; j++) ym1[j] = ym[j] - resfunm->Eval(mm[j]);
    
    TGraphErrors *gm1 = new TGraphErrors(Nbins,mm,ym1,0,yme);
    gm1->SetMarkerStyle(24);
    gm1->SetMarkerSize(MKSize-0.1);
    gm1->SetMarkerColor(2);
    gm1->SetLineColor(2);
    //gm1->Draw("psame");
    
    TMultiGraph *mg = new TMultiGraph();
    mg->Add(gm);
    //     mg->Add(gm1);
    mg->Draw("ap");
    //mg->GetYaxis()->SetTitle("Raw Counts (#times 10^{3})");
    mg->GetYaxis()->SetTitle("Raw Counts");
    //mg->GetYaxis()->SetTitleSize(titsize);
    mg->GetYaxis()->SetTitleOffset(1.5);
    //mg->GetYaxis()->SetLabelOffset(0.03);
    //mg->GetYaxis()->SetLabelSize(0.047);
    //mg->GetXaxis()->SetNdivisions(208);
    //  mg->GetXaxis()->CenterTitle();
    mg->GetXaxis()->SetTitle("Mass_{K#pi} (GeV/c^{2})");
    
    resfunm->Draw("same");
    fun0->Draw("same");
    resfunm1->Draw("same");
    
    
}


