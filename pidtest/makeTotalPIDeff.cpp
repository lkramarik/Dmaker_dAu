//
// Created by lukas on 19. 9. 2019.
//
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
#include "TVirtualFitter.h"
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
TF1* FitEffRange(TGraphErrors* gr, TString folder, double fitM, double fitMax);
TF1* FitEffRangeLinear(TGraphErrors* gr, TString folder, double fitM, double fitMax);
void makeTotalPIDeff();
double totalEffFct(double *x, double *par);
double ratioFct(double *x, double *par);
TGraphErrors* totalGraph(TGraphErrors*, TGraphErrors*, TGraphErrors*, TString, TString);
TF1* fTOF;
TF1* fTPC;
TF1* fTOFMatch;
TF1* fTotal;
TF1* fTotalGraph;

//---------------------------------------------------------------------------------------
TF1* FitEff(TGraphErrors* gr, TString folder){
    return FitEffRange(gr, folder, 0.15, 1.5);
}

//---------------------------------------------------------------------------------------
TF1* FitEffRange(TGraphErrors* gr, TString folder, double fitRangeMin=0.15, double fitRangeMax=1.5){
    gSystem->Exec("mkdir results_total_eff/"+folder);

//    TF1 *fitF1 = new TF1("eff_fit", "[0]+[1]/x+[2]*x+[3]*x*x", 0.15, 5); //momentum resolution fit
//    TF1 *fitF1 = new TF1("eff_fit", "0.5*([0]+[1]*x+[2]*x*x+[3]*x*x*x)", 0.1, 5);
//    fitF1->SetParameters(1, -0.06, -0.1, 0.02);

//    TF1 *fitF1 = new TF1("eff_fit", "[0]+[1]/x+[2]/x/x", 0.15, 4); //momentum resolution fit
//    TF1 *fitF1 = new TF1("eff_fit", "[0]+[1]/x+[2]*x+[3]*x*x+[4]/x/x", 0.15, 4); //momentum resolution fit

//    TF1 *fitF1 = new TF1("eff_fit", "[0]+[1]/x+[2]*x+[3]/x/x", 0.15, 5); //momentum resolution fit
//    fitF1->SetParameters(1, -0.06, -0.1, 0.02, 0.006);

//    TF1 *fitF1 = new TF1("eff_fit", "[0]+[1]/x+[2]/x/x", 0.15, 5); //momentum resolution fit
//    fitF1->SetParameters(1, -0.06, -0.1, 0.006);

    TF1 *fitF1 = new TF1("eff_fit", "[0]+[1]/x", 0.15, 5); //momentum resolution fit
    fitF1->SetParameters(1, -0.06, -0.1);

//    fitF1->SetParLimits(0, 0, 1);
//    fitF1->SetParLimits(1, 0, 1);
//    fitF1->SetParLimits(2, 0, 1);
//    fitF1->SetParLimits(3, -1, 0);
//    fitF1->SetParLimits(4, -1, 0);

    gr->Fit(fitF1, "RW", "", fitRangeMin, fitRangeMax);
    gr->GetYaxis()->SetRangeUser(0,1.2);

    TCanvas *out = new TCanvas("out", "out", 900, 1100);
    out->SetGrid();
    gr->Draw("ap");

    TH1D *hint = new TH1D("hint", "Fitted gaussian with .95 conf.band", 100, 0.15, 3);
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint);
    //Now the "hint" histogram has the fitted function values as the
    //bin contents and the confidence intervals as bin errors
    hint->GetYaxis()->SetRangeUser(0,1.2);

    hint->SetStats(kFALSE);
    hint->SetFillColor(0);
    hint->SetLineColor(17);
//    hint->Draw("e4 same");

    TString filename = "results_total_eff/"+folder+(TString)gr->GetName()+".png";
    out->SaveAs(filename);
    filename = "results_total_eff/"+folder+(TString)gr->GetName()+".pdf";
    out->SaveAs(filename);
//    out->Close();

    cout<<"chi2/ndf is "<<fitF1->GetChisquare()/fitF1->GetNDF()<<endl;

    TFile *file = new TFile("results_total_eff/"+folder+(TString)gr->GetName()+".root", "RECREATE");
    TString name = gr->GetName();
    name = "f_"+name;
    fitF1->Write(name, TObject::kOverwrite);
    gr->Write();
    file->Close();

    fitF1->SetName(gr->GetName());
    return fitF1;
}

//---------------------------------------------------------------------------------------
TF1* FitEffRangeLinear(TGraphErrors* gr, TString folder, double fitRangeMin=0.15, double fitRangeMax=1.5){
    gSystem->Exec("mkdir results_total_eff/"+folder);

    TF1 *fitF1 = new TF1("eff_fit", "[0]", 0.15, 5); //momentum resolution fit
//    fitF1->SetParameters(1, -0.06, -0.1);

//    fitF1->SetParLimits(0, 0, 1);
//    fitF1->SetParLimits(1, 0, 1);
//    fitF1->SetParLimits(2, 0, 1);
//    fitF1->SetParLimits(3, -1, 0);
//    fitF1->SetParLimits(4, -1, 0);

    gr->Fit(fitF1, "RW", "", fitRangeMin, fitRangeMax);
    gr->GetYaxis()->SetRangeUser(0,1.2);

    TCanvas *out = new TCanvas("out", "out", 900, 1100);
    out->SetGrid();
    gr->Draw("ap");

    TH1D *hint = new TH1D("hint", "Fitted gaussian with .95 conf.band", 100, 0.15, 3);
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint);
    //Now the "hint" histogram has the fitted function values as the
    //bin contents and the confidence intervals as bin errors
    hint->GetYaxis()->SetRangeUser(0,1.2);

    hint->SetStats(kFALSE);
    hint->SetFillColor(0);
    hint->SetLineColor(17);
//    hint->Draw("e4 same");

    TString filename = "results_total_eff/"+folder+(TString)gr->GetName()+".png";
    out->SaveAs(filename);
    filename = "results_total_eff/"+folder+(TString)gr->GetName()+".pdf";
    out->SaveAs(filename);
//    out->Close();

    cout<<"chi2/ndf is "<<fitF1->GetChisquare()/fitF1->GetNDF()<<endl;

    TFile *file = new TFile("results_total_eff/"+folder+(TString)gr->GetName()+".root", "RECREATE");
    TString name = gr->GetName();
    name = "f_"+name;
    fitF1->Write(name, TObject::kOverwrite);
    gr->Write();
    file->Close();

    fitF1->SetName(gr->GetName());
    return fitF1;
}

//---------------------------------------------------------------------------------------
double totalEffFct(double *x, double *par){
    return fTPC->Eval(x[0])*fTOFMatch->Eval(x[0])*fTOF->Eval(x[0])+(1.0-fTOFMatch->Eval(x[0]))*fTPC->Eval(x[0]);
}

//---------------------------------------------------------------------------------------
double ratioFct(double *x, double *par){
    return fTotal->Eval(x[0])/fTotalGraph->Eval(x[0]);
}

//---------------------------------------------------------------------------------------
TGraphErrors* totalGraph(TGraphErrors* gTOF, TGraphErrors* gTPC, TGraphErrors* gTOFMatch, TString folder, TString particle) {
    if ( (gTOF->GetN()!=gTPC->GetN()) || (gTOFMatch->GetN()!=gTPC->GetN()) || (gTOF->GetN()!=gTOFMatch->GetN())) cout<<"Your graphs are not compatible."<<endl;

    const int nBins=gTOF->GetN();
    Double_t tof, tpc, tofMatch, tofE, tpcE, tofMatchE;
    Double_t pT[nBins], ptWidths[nBins], eff[nBins], effE[nBins];

    for (int i = 0; i < gTOF->GetN(); ++i) {
        ptWidths[i] = gTOF->GetErrorX(i);

        tofE=gTOF->GetErrorY(i);
        tofMatchE=gTOFMatch->GetErrorY(i);
        tpcE=gTPC->GetErrorY(i);

        gTOF->GetPoint(i, pT[i], tof);
        gTPC->GetPoint(i, pT[i], tpc);
        gTOFMatch->GetPoint(i, pT[i], tofMatch);

        effE[i]=pow((tof*tpc-tpc)*tofMatchE,2)+pow((tofMatch*tof+1-tofMatch)*tpcE,2)+pow((tofMatch*tpc)*tofE,2);
        effE[i]=sqrt(effE[i]);

        eff[i]=tofMatch*tof*tpc+(1-tofMatch)*tpc;
    }

    TCanvas *out = new TCanvas("outTotal", "outTotal", 900, 1100);

    TGraphErrors *gTotal = new TGraphErrors(nBins, pT, eff, ptWidths, effE);
    gTotal->SetMarkerStyle(20);
    gTotal->SetMarkerSize(0.9);
    gTotal->SetMarkerColor(kBlack);
    gTotal->SetLineColor(kBlack);
    gTotal->GetYaxis()->SetTitle("Efficiency");
    gTotal->GetYaxis()->SetTitleOffset(1.1);
    gTotal->GetYaxis()->SetRangeUser(0,1.2);
    gTotal->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    gTotal->SetTitle("");
    gTotal->Draw("ap");
    out->SaveAs("results_total_eff/"+folder+"graphTotalEff_"+particle+".pdf");
    out->SaveAs("results_total_eff/"+folder+"graphTotalEff_"+particle+".png");
    out->Close();
    return gTotal;
}

//---------------------------------------------------------------------------------------
void makeTpcPidEff() {
    int bbcMin=0;
    int bbcMax=950;
    int nTof=-1;
    float nsigma=300;
    float tofInvBeta=999;
    float ptTrackCut=0.;
    TString particle = "pi";
//    TString particle = "K";
//    TString prefix = "tpc_";
    TString prefix = "tofPidEff_";

    TString cutComb=Form("bbc%i_%i_nHft%i_nsigma%.1f_tof%.2f_pt%.1f/",bbcMin, bbcMax, nTof, nsigma, tofInvBeta, ptTrackCut);
//    TString cutComb1=Form("bbc%i_%i_nHft%i_nsigma%.1f_tof%.2f_pt%.1f/",bbcMin, bbcMax, nTof, nsigma, tofInvBeta, 0.);
    TString pair = particle+particle;

    TFile *fTpcPID = new TFile(prefix+cutComb+"rootFiles/results_"+pair+".root", "READ");

    TGraphErrors *gTPC = (TGraphErrors*) fTpcPID->Get("eff");
    gTPC->SetNameTitle(prefix+particle, prefix+particle);

//    fTPC = FitEffRange(gTPC, cutComb, 0.3, 2);
    fTPC = FitEffRangeLinear(gTPC, cutComb, 0.3, 2);

//    fTpcPID->Close();

}

//---------------------------------------------------------------------------------------
void makeTotalPIDeff(){
    int bbcMin=0;
    int bbcMax=950;
    int nTof=-1;
    float nsigma=300;
    float tofInvBeta=999;
    float ptTrackCut=0.;
    TString particle = "pi";
//    TString particle = "K";

    TString cutComb=Form("bbc%i_%i_nHft%i_nsigma%.1f_tof%.2f_pt%.1f/",bbcMin, bbcMax, nTof, nsigma, tofInvBeta, ptTrackCut);
//    TString cutComb1=Form("bbc%i_%i_nHft%i_nsigma%.1f_tof%.2f_pt%.1f/",bbcMin, bbcMax, nTof, nsigma, tofInvBeta, 0.);
    TString pair = particle+particle;

//    TFile *fTpcPID = new TFile("tpc_bbc950_nHft0_nsigma3.0_tof0.03_pt0.5/rootFiles/results_"+pair+".root", "READ");
    TFile *fTpcPID = new TFile("tpc_"+cutComb+"rootFiles/results_"+pair+".root", "READ");
    TFile *fTofPID = new TFile("tofPidEff_"+cutComb+"rootFiles/results_"+pair+".root", "READ");
    TFile *fTofMatch = new TFile(cutComb+"rootFiles/results_"+pair+".root", "READ");
//    TFile *fTofHybrid = new TFile("hybrid_"+cutComb+"rootFiles/results_"+pair+".root", "READ");
//    TFile *fileTofHybrid = new TFile("hybrid_bbc950_nHft0/rootFiles/results_"+pair+".root", "READ");

    TGraphErrors *gTOF = (TGraphErrors*) fTofPID->Get("eff");
    gTOF->SetNameTitle("TOF_PID_"+particle, "TOF PID "+particle);

    TGraphErrors *gTPC = (TGraphErrors*) fTpcPID->Get("eff");
    gTPC->SetNameTitle("TPC_PID_"+particle, "TPC PID "+particle);

    TGraphErrors *gTOFMatch = (TGraphErrors*) fTofMatch->Get("eff");
    gTOFMatch->SetNameTitle("TOF_match_"+particle, "TOF matching "+particle);

//    TGraphErrors *gTOFHybrid = (TGraphErrors*) fileTofHybrid->Get("eff");
//    gTOFHybrid->SetNameTitle("TOF_hybrid_"+particle, "TOF hybrid "+particle);

    TGraphErrors *gTotal = totalGraph(gTOF, gTPC, gTOFMatch, cutComb, particle);
    gTotal->SetNameTitle("total_eff_"+particle, "total eff. "+particle);
    fTotalGraph = FitEff(gTotal, cutComb);
    fTotalGraph->SetName("fTotalGraph");
    fTotalGraph->SetLineColor(1);
//    fTotalGraph->SetLineStyle(1);


    fTOF = FitEff(gTOF, cutComb);
    fTOF->SetName("fTOF");
    fTPC = FitEff(gTPC, cutComb);
    fTPC->SetName("fTPC");
    fTOFMatch = FitEff(gTOFMatch, cutComb);
    fTOFMatch->SetName("fTOFMatch");
//    TF1* fTOFHybrid = FitEff(gTOFHybrid, cutComb);
//    fTOFHybrid->SetName("fTOFHybrid");

//    TF1* fTotal = new TF1("fTotal", "fTPC*fTOFMatch*fTOF+(1-fTOFMatch)*fTPC", 0.15, 4);
//    TF1* fTotal = new TF1("fTotal", "(fTPC(x))*(fTOFMatch(x))*(fTOF(x))+(1.0-(fTOFMatch(x)))*(fTPC(x))", 0.15, 4);
    fTotal = new TF1("fTotal", totalEffFct, 0.15, 4);

//    TF1* fTotal1 = new TF1("fTotal1", "(fTPC(x))*(fTOFMatch(x))*(fTOF(x))", 0.15, 4);
//    TF1* fTotal2 = new TF1("fTotal2", "(1.0-(fTOFMatch(x)))*(fTPC(x))", 0.15, 4);
//    TF1* fTotal = new TF1("fTotal", "(fTotal1(x))+(fTotal2(x))", 0.15, 4);

    //    TF1* fTotal = new TF1("fTotal", "(fTPC(0))*(fTOFMatch(5))*(fTOF(10))+(1.0-(fTOFMatch(15)))*(fTPC(20))", 0.15, 4);
    cout<<fTotal->Eval(2)<<endl;
    cout<<fTPC->Eval(2)*fTOFMatch->Eval(2)*fTOF->Eval(2)+(1-fTOFMatch->Eval(2))*fTPC->Eval(2)<<endl;

//    TF1* fTotalHybrid = new TF1("fTotalHybrid", "(fTPC(x))*(fTOFHybrid(x))", 0.15, 4);

    TCanvas *out = new TCanvas("out", "out", 900, 1100);
    out->SetGrid();
    gPad->SetLeftMargin(0.15);
    fTotal->Draw();
//    fTotalHybrid->Draw("same");
    fTotalGraph->Draw("same");


    TLegend *legend = new TLegend(0.58, 0.13, 0.924, 0.283);
    legend->SetFillStyle(0);
    legend->AddEntry(fTotal, "Separate fits");
    legend->AddEntry(fTotalGraph, "Fit of total graph");
    legend->SetBorderSize(0);
    legend->Draw("same");

    fTotal->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fTotal->GetYaxis()->SetTitle("#varepsilon_{PID} for "+particle);
    fTotal->GetYaxis()->SetRangeUser(0, 1.2);
    fTotal->SetTitle("");

    out->SaveAs("results_total_eff/"+cutComb+"fTotalEffPid_"+particle+".png");
    out->SaveAs("results_total_eff/"+cutComb+"fTotalEffPid_"+particle+".pdf");

    TF1* fRatio = new TF1("fRatio", ratioFct, 0.15, 4);
    cout<<fTotal->Eval(3.5)<<endl;
    cout<<fTotalGraph->Eval(3.5)<<endl;
    cout<<fTotalGraph->Eval(3.5)+fTotal->Eval(3.5)<<endl;
    cout<<fRatio->Eval(3.5)<<endl;

//    fRatio->Draw();

    TFile *file = new TFile("results_total_eff/"+cutComb+"totalEff_"+particle+".root", "RECREATE");
    file->SetCompressionAlgorithm(1);
    fTotal->Write("fTotalEffPid_"+particle);
    fTotalGraph->Write("fTotalGraphEffPid_"+particle);
    gTotal->Write("grTotalGraphEffPid_"+particle);
    file->Close();

//    TCanvas *out1 = new TCanvas("out1", "out1", 1500, 1000);
//    hTOF->Draw();

    fTofPID->Close();
    fTofMatch->Close();
    fTpcPID->Close();
//    fileTofHybrid->Close();

}

