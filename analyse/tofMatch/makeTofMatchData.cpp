#include "TH3D.h"
#include "TH3F.h"
#include "TFile.h"
#include "TH2D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TString.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TGraphErrors.h"
#include "TSpline.h"


#include <iostream>
using namespace std;
Int_t lineColor=9;
TString xAxis="p_{T} (GeV/c)";
TString yAxis="TOF matching efficiency";
Float_t maxYaxis=0.7;
Int_t color[]={1,9,8};
Int_t selectedSigma=0;

void setHisto(TH1D* h){
    h->Rebin(8);
    h->SetStats(0);
    h->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h->SetTitle("");
    h->GetXaxis()->SetRangeUser(0.,5.);
    h->GetYaxis()->SetRangeUser(0.,1.);
    h->GetYaxis()->SetRangeUser(0.,1.);
}

//-------------------------------------------------------------------------
TGraphErrors* histoToGraph(TH1D* h, Double_t excludeMin, Double_t excludeMax) {
    const int nPoints = h->GetNbinsX();

    TGraphErrors* gr=new TGraphErrors();

    double content=0;
    double binCenter=-1;
    Double_t x, y, xE, yE;
    int nAddedPoints = 0;
    for (int i = 0; i < nPoints; ++i) {
        binCenter=h->GetBinCenter(i+1);
        if ((binCenter>excludeMin && binCenter<excludeMax) || binCenter>5.) {
            continue;
        }
        content=h->GetBinContent(i+1);

        gr->SetPoint(nAddedPoints, binCenter, content);
        gr->SetPointError(nAddedPoints, 0, h->GetBinError(i+1));
//        gr->SetPointError(nAddedPoints, 0, 0);
        nAddedPoints++;
    }
    TString name=h->GetName();
    name="gr"+name;
    gr->SetName(name);
    return gr;
}

//-------------------------------------------------------------------------
TGraphErrors* histoToGraphNotFit(TH1D* h, Double_t excludeMin, Double_t excludeMax) {
    const int nPoints = h->GetNbinsX();

    TGraphErrors* gr=new TGraphErrors();

    double content=0;
    double binCenter=-1;
    Double_t x, y, xE, yE;
    int nAddedPoints = 0;
    for (int i = 0; i < nPoints; ++i) {
        binCenter=h->GetBinCenter(i+1);
        if ((binCenter<excludeMin || binCenter>excludeMax) || binCenter>5.) {
            continue;
        }
        content=h->GetBinContent(i+1);

        gr->SetPoint(nAddedPoints, binCenter, content);
        gr->SetPointError(nAddedPoints, 0, h->GetBinError(i+1));
//        gr->SetPointError(nAddedPoints, 0, 0);
        nAddedPoints++;
    }
    TString name=h->GetName();
    name="gr"+name;
    gr->SetName(name);
    return gr;
}

//-------------------------------------------------------------------------
TSpline3* histoToSpline(TH1D* h, Double_t excludeMin, Double_t excludeMax) {
    const int nPoints = h->GetNbinsX();

    TGraph* gr=new TGraph();

    double content=0;
    double binCenter=-1;
    Double_t x, y, xE, yE;
    int nAddedPoints = 0;
    for (int i = 0; i < nPoints; ++i) {
        binCenter=h->GetBinCenter(i+1);
        if ((binCenter>excludeMin && binCenter<excludeMax) || binCenter>5.) {
            continue;
        }
        content=h->GetBinContent(i+1);

//        splineSpectra = new TSpline3("Spline", pTspline, y, nPoints, "b2e2", 0, 0);
//        splineSpectra->SetLineColor(9);

        gr->SetPoint(nAddedPoints, binCenter, content);
//        spline->SetPointError(nAddedPoints, 0, h->GetBinError(i+1));
//        gr->SetPointError(nAddedPoints, 0, 0);
        nAddedPoints++;
    }

    TSpline3* spline=new TSpline3("apl", gr, "b2e2", 0, 0);
    spline->SetLineColor(9);

    return spline;
}

//-------------------------------------------------------------------------
//void makeTofMatchData(TString inputFileName="tofMatch.primaries.hotspot.root") {
void makeTofMatchData(TString inputFileName="tofratio.all.1509.root") {
    TString inputFileNameOrig=inputFileName;
    inputFileNameOrig.ReplaceAll(".root","");
    gSystem->Exec(Form("mkdir %s ", inputFileNameOrig.Data()));

    const int m_nParticlesCharged = 4;
    const TString m_ParticleChargedName[m_nParticlesCharged] = {"PionPlus", "PionMinus", "KaonPlus", "KaonMinus"};
    const TString m_ParticleChargedName_plot[m_nParticlesCharged] = {"#pi^{+}", "#pi^{-}", "K^{+}", "K^{-}"};

    TFile* inputFile = new TFile(inputFileName, "READ");

    TH1D* h1Tofmatch[m_nParticlesCharged][3];
    TH1D* h1TofmatchTOF[m_nParticlesCharged][3];
    TGraphErrors* grTOF[m_nParticlesCharged][3];
    TGraphErrors* grTOFTotal[m_nParticlesCharged][3];
    TGraphErrors* grTOFNotFit[m_nParticlesCharged][3];
    TSpline3* splineTOF[m_nParticlesCharged][3];

    TString hisName;
    TCanvas* cNsigma=new TCanvas("cNsigma","cNsigma", 1200,1500);
    cNsigma->Divide(2,2);
    TCanvas* cNsigmaGR=new TCanvas("cNsigmaGR","cNsigmaGR", 1200,1500);
    cNsigmaGR->Divide(2,2);
    TCanvas* cTmp=new TCanvas("cTmp","cTmp", 1200,1000);

    TF1 *fitF1[m_nParticlesCharged][3];
    TF1 *fitF1Line[m_nParticlesCharged][3];
    TF1 *fitTotal[m_nParticlesCharged][3];
//    TF1 *fitF1 = new TF1("eff_fit", "[0]+[1]/x+[2]/x/x+[3]/x/x/x+[4]/x/x/x", 0.15, 2.2); //momentum resolution fit
//    TF1 *fitF1 = new TF1("eff_fit", "[0]+[1]/x+[2]*x+[3]*x*x", 0.15, 1); //momentum resolution fit
//    TF1 *fitF1 = new TF1("fitF1","[0]*log(-[1]*x)",0.,3);

//    TF1 *fitF1 = new TF1("eff_fit", "[0]+[1]/x+[2]*x", 0.15, 3); //momentum resolution fit
//    TF1 *fitF1 = new TF1("eff_fit", "[0]+[1]/x", 0.15, 1); //momentum resolution fit
//    fitF1->SetParameters(1, -0.06, -0.1, 0.02, 0.006);
//    fitF1->SetParLimits(0, 0, 1);
//    fitF1->SetParLimits(1, 0, 1);
//    fitF1->SetParLimits(2, 0, 1);
//    fitF1->SetParLimits(3, -1, 0);

    Double_t excludeMinKaon=0.36;
    Double_t excludeMaxKaon=2;

    Double_t excludeMinPion=0.5;
    Double_t excludeMaxPion=2.2;

    inputFileName="out_"+inputFileName;
    TFile* outputFile = new TFile(inputFileName, "RECREATE");

    TLegend *legend = new TLegend(0.48, 0.14, 0.634, 0.32);
    legend -> SetFillStyle(0);
    legend -> SetLineColor(0);
    legend -> SetTextSize(0.05);

    TPaveText *text1 = new TPaveText(0.66, 0.14, 0.87, 0.188, "NDC");
    TString textData=Form("TPC |n#sigma| < %i",selectedSigma+1);
    text1->AddText(textData);
    text1->SetFillColor(0);

    for (int iParticle = 0; iParticle < m_nParticlesCharged; ++iParticle) {
        for (int nsigma = 0; nsigma < 3; ++nsigma) {

            hisName=Form("h1_Tofmatch_tpc_%s_nsigma%d", m_ParticleChargedName[iParticle].Data(), nsigma+1);
            h1Tofmatch[iParticle][nsigma] = (TH1D*)inputFile->Get(hisName);
            hisName=Form("h1_Tofmatch_tof_%s_nsigma%d", m_ParticleChargedName[iParticle].Data(), nsigma+1);
            h1TofmatchTOF[iParticle][nsigma] = (TH1D*)inputFile->Get(hisName);

            setHisto(h1Tofmatch[iParticle][nsigma]);
            setHisto(h1TofmatchTOF[iParticle][nsigma]);

            h1TofmatchTOF[iParticle][nsigma]->SetTitle(m_ParticleChargedName_plot[iParticle]) ;

            h1TofmatchTOF[iParticle][nsigma]->Divide(h1Tofmatch[iParticle][nsigma]);



            if (iParticle<2){ //pion

                splineTOF[iParticle][nsigma]=histoToSpline(h1TofmatchTOF[iParticle][nsigma], excludeMinPion, excludeMaxPion);

                grTOF[iParticle][nsigma]=histoToGraph(h1TofmatchTOF[iParticle][nsigma], excludeMinPion, excludeMaxPion);
                grTOFTotal[iParticle][nsigma]=histoToGraph(h1TofmatchTOF[iParticle][nsigma], 0., 0.);
                grTOFNotFit[iParticle][nsigma]=histoToGraphNotFit(h1TofmatchTOF[iParticle][nsigma], excludeMinPion, excludeMaxPion);

                grTOF[iParticle][nsigma]->SetTitle(m_ParticleChargedName_plot[iParticle]);
                grTOFTotal[iParticle][nsigma]->SetTitle(m_ParticleChargedName_plot[iParticle]);
                fitF1[iParticle][nsigma] = new TF1("eff_fit", "[0]+[1]/x+[3]/x/x", 0.165, 5); //go
                fitF1[iParticle][nsigma]->SetLineColor(lineColor);

//                fitF1[iParticle][nsigma] = new TF1("eff_fit", "[0]+[1]/x/x", 0.15, excludeMinPion);
                grTOF[iParticle][nsigma]->Fit(fitF1[iParticle][nsigma], "REX0W");

//                fitF1[iParticle][nsigma] = new TF1("eff_fit", "gaus(0)", 0.15, excludeMinPion+0.5);
//                fitF1[iParticle][nsigma]->SetLineColor(1);
//                grTOF[iParticle][nsigma]->Fit(fitF1[iParticle][nsigma], "R");
//
//                fitF1Line[iParticle][nsigma] = new TF1("eff_fit_line", "[0]", excludeMaxPion, 5);
//                fitF1Line[iParticle][nsigma]->SetLineColor(1);
//                grTOF[iParticle][nsigma]->Fit(fitF1Line[iParticle][nsigma], "R+");

//                Double_t par[4];
//                fitTotal[iParticle][nsigma]=new TF1("total", "gaus(0)+[3]",0.15,4.5);
//                fitF1[iParticle][nsigma]->GetParameters(&par[0]);
//                fitF1Line[iParticle][nsigma]->GetParameters(&par[3]);
//
//                fitTotal[iParticle][nsigma]->SetParameters(par);
//                fitTotal[iParticle][nsigma]->SetLineColor(8);
//
//                grTOF[iParticle][nsigma]->Fit(fitTotal[iParticle][nsigma],"R+");

                outputFile->cd();
                hisName=Form("gr_tofMatch_%s_nsigma%i", m_ParticleChargedName[iParticle].Data(), nsigma);
                grTOF[iParticle][nsigma]->Write(hisName);
                hisName=Form("f_tofMatch_%s_nsigma%i", m_ParticleChargedName[iParticle].Data(), nsigma);
                fitF1[iParticle][nsigma]->Write(hisName);
            }

            if (iParticle>1){ //kaon
                splineTOF[iParticle][nsigma]=histoToSpline(h1TofmatchTOF[iParticle][nsigma], excludeMinKaon, excludeMaxKaon);

                grTOF[iParticle][nsigma]=histoToGraph(h1TofmatchTOF[iParticle][nsigma], excludeMinKaon, excludeMaxKaon);
                grTOFTotal[iParticle][nsigma]=histoToGraph(h1TofmatchTOF[iParticle][nsigma], 0., 0.);
                grTOFNotFit[iParticle][nsigma]=histoToGraphNotFit(h1TofmatchTOF[iParticle][nsigma], excludeMinKaon, excludeMaxKaon);

                grTOF[iParticle][nsigma]->SetTitle(m_ParticleChargedName_plot[iParticle]);
                grTOFTotal[iParticle][nsigma]->SetTitle(m_ParticleChargedName_plot[iParticle]);

                fitF1[iParticle][nsigma] = new TF1("eff_fit", "[0]+[1]/x+[2]*x+[3]/x/x", 0.15, 5); //good
//                fitF1[iParticle][nsigma] = new TF1("eff_fit", "[0]+[1]/x+[2]*x+[3]/x/x", 0.15, excludeMinKaon);
                fitF1[iParticle][nsigma]->SetLineColor(lineColor);
                grTOF[iParticle][nsigma]->Fit(fitF1[iParticle][nsigma], "REX0W");

                outputFile->cd();
                hisName=Form("gr_tofMatch_%s_nsigma%i", m_ParticleChargedName[iParticle].Data(), nsigma);
                grTOF[iParticle][nsigma]->Write(hisName);
                hisName=Form("f_tofMatch_%s_nsigma%i", m_ParticleChargedName[iParticle].Data(), nsigma);
                fitF1[iParticle][nsigma]->Write(hisName);

            }
//            h1TofmatchTOF[iParticle][nsigma]->Fit("eff_fit", "", "", 0.15, 5);



            grTOF[iParticle][nsigma]->GetXaxis()->SetTitle(xAxis);
            grTOF[iParticle][nsigma]->GetYaxis()->SetTitle(yAxis);
            grTOF[iParticle][nsigma]->GetYaxis()->SetRangeUser(0.,maxYaxis);
            grTOF[iParticle][nsigma]->GetXaxis()->SetRangeUser(0.,5.);
            grTOF[iParticle][nsigma]->SetMarkerStyle(20);
            grTOF[iParticle][nsigma]->SetMarkerSize(0.8);
            grTOF[iParticle][nsigma]->SetMarkerColor(color[nsigma]);
            grTOF[iParticle][nsigma]->SetLineColor(color[nsigma]);

            grTOFTotal[iParticle][nsigma]->GetXaxis()->SetTitle(xAxis);
            grTOFTotal[iParticle][nsigma]->GetYaxis()->SetTitle(yAxis);
            grTOFTotal[iParticle][nsigma]->GetYaxis()->SetRangeUser(0.,maxYaxis);
            grTOFTotal[iParticle][nsigma]->GetXaxis()->SetRangeUser(0.,5.);
            grTOFTotal[iParticle][nsigma]->SetMarkerStyle(20);
            grTOFTotal[iParticle][nsigma]->SetMarkerSize(0.8);
            grTOFTotal[iParticle][nsigma]->SetMarkerColor(color[nsigma]);
            grTOFTotal[iParticle][nsigma]->SetLineColor(color[nsigma]);

            TString legendText=Form("TPC |n#sigma| < %i",nsigma+1);
            if (iParticle==0) legend -> AddEntry(grTOFTotal[iParticle][nsigma], legendText, "pl");


            grTOFNotFit[iParticle][nsigma]->SetMarkerStyle(24);
            grTOFNotFit[iParticle][nsigma]->SetMarkerSize(0.8);
            grTOFNotFit[iParticle][nsigma]->SetMarkerColor(46);

            cNsigma->cd(iParticle+1);
            cNsigma->SetLeftMargin(0.2);
            if (nsigma==0) grTOFTotal[iParticle][nsigma]->Draw("apl");
            else grTOFTotal[iParticle][nsigma]->Draw("same pl");


            if (nsigma==selectedSigma) {
                cNsigmaGR->cd(iParticle + 1);
                cNsigmaGR->SetLeftMargin(0.2);

                if (nsigma == 0) grTOF[iParticle][nsigma]->Draw("ap");
                else grTOF[iParticle][nsigma]->Draw("same ap");
                grTOFNotFit[iParticle][nsigma]->Draw("same p");
                fitF1[iParticle][nsigma]->Draw("same");
                text1->Draw("same");

            }
//            if (iParticle<2) //pion
//                fitTotal[iParticle][nsigma]->Draw("same");
        }

    }
    cNsigma->cd(4);
    legend->Draw("same");
    cTmp->Close();

    cNsigmaGR->cd(4);



    cNsigma->SaveAs(Form("%s/tofMatchAll.png", inputFileNameOrig.Data()));
    cNsigma->SaveAs(Form("%s/tofMatchAll.eps", inputFileNameOrig.Data()));
    cNsigmaGR->SaveAs(Form("%s/tofMatchFit.png", inputFileNameOrig.Data()));
    cNsigmaGR->SaveAs(Form("%s/tofMatchFit.eps", inputFileNameOrig.Data()));

//    for (int i = 0; i < ; ++i) {
//
//    }



}