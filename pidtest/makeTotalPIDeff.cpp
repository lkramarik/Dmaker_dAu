//
// Created by lukas on 19. 9. 2019.
//
#include "FitPID.h"
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
#include "TCut.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include <iostream>
#include <ctime>
#include <cstdio>
#include <fstream>
#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include "TSystem.h"
#include "TROOT.h"

TF1* FitEff(TGraphErrors* gr, TString folder);
void makeTotalPIDeff();

TF1* FitEff(TGraphErrors* gr, TString folder){
    gSystem->Exec("mkdir results_total_eff/"+folder);


    TF1 *fitF1 = new TF1("eff_fit", "[0]+[1]/x+[2]*x+[3]*x*x+[4]/x/x", 0.15, 3); //momentum resolution fit
    gr->Fit(fitF1, "R");

    TCanvas *out = new TCanvas("out", "out", 900, 1100);
    gr->Draw("ap");
    TString filename = "results_total_eff/"+folder+(TString)gr->GetName()+".png";
    out->SaveAs(filename);
    filename = "results_total_eff/"+folder+(TString)gr->GetName()+".pdf";
    out->SaveAs(filename);
    out->Close();

    cout<<"chi2/ndf is "<<fitF1->GetChisquare()/fitF1->GetNDF()<<endl;

    TFile *file = new TFile("results_total_eff/"+folder+(TString)gr->GetName()+".root", "RECREATE");
    fitF1->Write(gr->GetName(), TObject::kOverwrite);
    file->Close();

    fitF1->SetName(gr->GetName());
    return fitF1;
}

void makeTotalPIDeff(){
    int bbc=800;
    int nTof=0;
    float nsigma=3;
    float tofInvBeta=0.03;
    float ptTrackCut=0.5;
    TString particle = "K";

    TString cutComb=Form("bbc%i_nHft%i_nsigma%.1f_tof%.2f_pt%.1f/",bbc, nTof, nsigma, tofInvBeta, ptTrackCut);
    TString pair = particle+particle;

    TFile *fTpcPID = new TFile("tpc_bbc950_nHft0_nsigma3.0_tof0.03_pt0.5/rootFiles/results_"+pair+".root", "READ");
//    TFile *fTpcPID = new TFile("tpc_"+cutComb+"rootFiles/results_"+pair+".root", "READ");
    TFile *fTofPID = new TFile("tofPidEff_"+cutComb+"rootFiles/results_"+pair+".root", "READ");
    TFile *fTofMatch = new TFile(cutComb+"rootFiles/results_"+pair+".root", "READ");
//    TFile *fTofHybrid = new TFile("hybrid_"+cutComb+"rootFiles/results_"+pair+".root", "READ");
    TFile *fileTofHybrid = new TFile("hybrid_bbc950_nHft0/rootFiles/results_"+pair+".root", "READ");


    TGraphErrors *gTOF = (TGraphErrors*) fTofPID->Get("eff");
    gTOF->SetNameTitle("TOF_PID_"+particle, "TOF PID "+particle);

    TGraphErrors *gTPC = (TGraphErrors*) fTpcPID->Get("eff");
    gTPC->SetNameTitle("TPC_PID_"+particle, "TPC PID "+particle);

    TGraphErrors *gTOFMatch = (TGraphErrors*) fTofMatch->Get("eff");
    gTOFMatch->SetNameTitle("TOF_match_"+particle, "TOF matching "+particle);

    TGraphErrors *gTOFHybrid = (TGraphErrors*) fileTofHybrid->Get("eff");
    gTOFHybrid->SetNameTitle("TOF_hybrid_"+particle, "TOF hybrid "+particle);

    TF1* fTOF = FitEff(gTOF, cutComb);
    fTOF->SetName("fTOF");
    TF1* fTPC = FitEff(gTPC, cutComb);
    fTPC->SetName("fTPC");
    TF1* fTOFMatch = FitEff(gTOFMatch, cutComb);
    fTOFMatch->SetName("fTOFMatch");
    TF1* fTOFHybrid = FitEff(gTOFHybrid, cutComb);
    fTOFHybrid->SetName("fTOFHybrid");

    TF1* fTotal = new TF1("fTotal", "(fTPC(x))*(fTOFMatch(x))*(fTOF(x))+(1-(fTOFMatch(x)))*(fTPC(x))", 0.15, 4);
    TF1* fTotalHybrid = new TF1("fTotalHybrid", "(fTPC(x))*(fTOFHybrid(x))", 0.15, 4);

    TCanvas *out = new TCanvas("out", "out", 900, 1100);
    gPad->SetLeftMargin(0.15);
    fTotal->Draw();
    fTotalHybrid->Draw("same");

    fTotal->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fTotal->GetYaxis()->SetTitle("#varepsilon_{PID}");
    fTotal->SetTitle("");

    out->SaveAs("fTotalEffPid_"+particle+".png");
    out->SaveAs("fTotalEffPid_"+particle+".pdf");

    TFile *file = new TFile("results_total_eff/"+cutComb+"totalEff_"+particle+".root", "RECREATE");
    fTotal->Write("fTotalEffPid_"+particle);
    file->Close();

//    TCanvas *out1 = new TCanvas("out1", "out1", 1500, 1000);
//    hTOF->Draw();

    fTofPID->Close();
    fTofMatch->Close();
    fTpcPID->Close();
    fileTofHybrid->Close();


}

