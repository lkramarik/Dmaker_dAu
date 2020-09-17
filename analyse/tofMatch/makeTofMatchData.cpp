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
#include "TGraphErrors.h"
#include "TSpline.h"


#include <iostream>
using namespace std;

void setHisto(TH1D* h){
//    h->Rebin(4);
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
void makeTofMatchData(TString inputFileName="tofratio.all.1509.root") {
    const int m_nParticlesCharged = 4;
    const TString m_ParticleChargedName[m_nParticlesCharged] = {"PionPlus", "PionMinus", "KaonPlus", "KaonMinus"};

    TFile* inputFile = new TFile(inputFileName, "READ");

    TH1D* h1Tofmatch[m_nParticlesCharged][3];
    TH1D* h1TofmatchTOF[m_nParticlesCharged][3];
    TGraphErrors* grTOF[m_nParticlesCharged][3];
    TSpline3* splineTOF[m_nParticlesCharged][3];

    TString hisName;
    TCanvas* cNsigma=new TCanvas("cNsigma","cNsigma", 1200,1200);
    cNsigma->Divide(2,2);
    TCanvas* cNsigmaGR=new TCanvas("cNsigmaGR","cNsigmaGR", 1200,1200);
    cNsigmaGR->Divide(2,2);
    TCanvas* cTmp=new TCanvas("cTmp","cTmp", 1200,1000);

    TF1 *fitF1[m_nParticlesCharged][3];
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

    Double_t excludeMinKaon=0.3;
    Double_t excludeMaxKaon=2;

    Double_t excludeMinPion=0.5;
    Double_t excludeMaxPion=2.2;

    TFile* outputFile = new TFile("tofMatchEfficiency.root", "RECREATE");

    for (int iParticle = 0; iParticle < m_nParticlesCharged; ++iParticle) {
        for (int nsigma = 1; nsigma < 2; ++nsigma) {

            hisName=Form("h1_Tofmatch_tpc_%s_nsigma%d", m_ParticleChargedName[iParticle].Data(), nsigma+1);
            h1Tofmatch[iParticle][nsigma] = (TH1D*)inputFile->Get(hisName);
            hisName=Form("h1_Tofmatch_tof_%s_nsigma%d", m_ParticleChargedName[iParticle].Data(), nsigma+1);
            h1TofmatchTOF[iParticle][nsigma] = (TH1D*)inputFile->Get(hisName);

            setHisto(h1Tofmatch[iParticle][nsigma]);
            setHisto(h1TofmatchTOF[iParticle][nsigma]);

            h1TofmatchTOF[iParticle][nsigma]->Divide(h1Tofmatch[iParticle][nsigma]);


            if (iParticle<2){ //pion
                splineTOF[iParticle][nsigma]=histoToSpline(h1TofmatchTOF[iParticle][nsigma], excludeMinPion, excludeMaxPion);

                grTOF[iParticle][nsigma]=histoToGraph(h1TofmatchTOF[iParticle][nsigma], excludeMinPion, excludeMaxPion);

//                fitF1[iParticle][nsigma] = new TF1("eff_fit", "[0]+[1]/x+[3]/x/x", 0.165, 4); //good
                fitF1[iParticle][nsigma] = new TF1("eff_fit", "[0]+[1]/x/x", 0.165, 4);
                grTOF[iParticle][nsigma]->Fit(fitF1[iParticle][nsigma], "REX0W");

                outputFile->cd();
                hisName=Form("gr_tofMatch_%s_nsigma%i", m_ParticleChargedName[iParticle].Data(), nsigma);
                grTOF[iParticle][nsigma]->Write(hisName);
                hisName=Form("f_tofMatch_%s_nsigma%i", m_ParticleChargedName[iParticle].Data(), nsigma);
                fitF1[iParticle][nsigma]->Write(hisName);
            }

            if (iParticle>1){ //kaon
                splineTOF[iParticle][nsigma]=histoToSpline(h1TofmatchTOF[iParticle][nsigma], excludeMinKaon, excludeMaxKaon);

                grTOF[iParticle][nsigma]=histoToGraph(h1TofmatchTOF[iParticle][nsigma], excludeMinKaon, excludeMaxKaon);

//                fitF1[iParticle][nsigma] = new TF1("eff_fit", "[0]+[1]/x+[2]*x+[3]/x/x", 0.2, 3); //good
                fitF1[iParticle][nsigma] = new TF1("eff_fit", "[0]+[1]/x+[2]*x+[3]/x/x", 0.2, 3);
                grTOF[iParticle][nsigma]->Fit(fitF1[iParticle][nsigma], "REX0W");

                outputFile->cd();
                hisName=Form("gr_tofMatch_%s_nsigma%i", m_ParticleChargedName[iParticle].Data(), nsigma);
                grTOF[iParticle][nsigma]->Write(hisName);
                hisName=Form("f_tofMatch_%s_nsigma%i", m_ParticleChargedName[iParticle].Data(), nsigma);
                fitF1[iParticle][nsigma]->Write(hisName);

            }
//            h1TofmatchTOF[iParticle][nsigma]->Fit("eff_fit", "", "", 0.15, 5);

            cNsigma->cd(iParticle+1);
            if (nsigma==0) h1TofmatchTOF[iParticle][nsigma]->Draw();
            else h1TofmatchTOF[iParticle][nsigma]->Draw("same");
            splineTOF[iParticle][nsigma]->Draw("same");
            fitF1[iParticle][nsigma]->Draw("same");

            cNsigmaGR->cd(iParticle+1);
            if (nsigma==0) grTOF[iParticle][nsigma]->Draw("ap");
            else grTOF[iParticle][nsigma]->Draw("same ap");
        }
    }

    cTmp->Close();

//    for (int i = 0; i < ; ++i) {
//
//    }



}