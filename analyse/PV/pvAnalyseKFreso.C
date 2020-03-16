#include "TCut.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TNtuple.h"
#include "TPaveText.h"
#include "TFile.h"
#include "TF1.h"
#include "TLatex.h"
#include <iostream>
#include "TSystem.h"

using namespace std;

void pvAnalyseKFreso() {
//    TString input = "ntp.KFVertex.global.2105.root";
    TString input = "ntp.PVrefit.global.root";
//    TString input = "ntp.KFVertexReso.1605.root";
    TFile *data = new TFile(input, "r");
    TNtuple *ntp = (TNtuple *) data->Get("ntp_KFReso");
    TFile *fOut = new TFile("res.KFreso.prim30." + input, "recreate");

    TH1F *hVar = new TH1F();

    bool draw = false;
    Int_t col[] = {46, 8};

    TString var[] = {"", "x", "y", "z"};

    Float_t limsMin[] = {0, -0.6, -1, -1};
    Float_t limsMax[] = {1, 0.6, 1, 1};

    int nPrimMin = 10;
    int nHftMin = 2;
//    TString detCuts = "nPrimTracks>0 && nHftTracks>1.5 && abs(picoDstVx)<0.5 && abs(picoDstVy)<0.5";
    TString detCuts = Form("nPrimTracks>%i && nHftTracks>%i", nPrimMin, nHftMin);
//    TString detCuts = Form("nPrimTracks>%i && nHftTracks>%i", nPrimMin, nHftMin);
    TString folder = Form("nPrimMin%i_nHftMin%i_%s/", nPrimMin, nHftMin, input.Data());
    gSystem->Exec(Form("mkdir %s", folder.Data()));

//    TString detCuts = Form("nPrimTracks>%i", nPrimMin);

    TString hisName, varName;
    TCut cut1, cut2;
    TCanvas *c = new TCanvas("c1", "c1", 900, 1200);
    c->SetGrid();

//    Double_t fwhm = 0, sigma = 0;

    for (int j = 0; j < 4; ++j) { //variable to project
        TPaveText *text4 = new TPaveText(0.35, 0.848, 0.229, 0.873, "brNDC");
        text4->SetTextSize(0.03);
        text4->SetLineColor(0);
        text4->SetShadowColor(0);
        text4->SetFillColor(0);

        varName = Form("KFdiff%s", var[j].Data());
        hVar = new TH1F(var[j], var[j], 1500, limsMin[j], limsMax[j]);
        ntp->Project(var[j], varName, detCuts);
        hVar->Rebin(2);
        hVar->SetStats(0);
        hVar->SetTitle(detCuts);
        hVar->Scale(1/hVar->GetEntries());
        hVar->SetFillColor(46);
        hVar->SetFillStyle(3004);
        hVar->SetMarkerStyle(20);
        hVar->SetLineColor(46);
        hVar->SetMarkerColor(46);
        fOut->cd();
        hVar->GetXaxis()->SetTitle(Form("#Delta%s(KF1,KF2) [cm]", var[j].Data()));
        hVar->GetYaxis()->SetTitleOffset(1.3);
        hVar->GetYaxis()->SetTitle("1/N");
        TF1 *fun0 = new TF1("fun0", "[constant]/( pi*[gamma]*(1+(x-[mean])*(x-[mean])/([gamma]*[gamma])))+[jump]", -0.05, 0.05);
        fun0->SetParameter(1, 0.005);
        fun0->SetParLimits(1, 0.001, 1); //gamma
        fun0->SetParLimits(2, -0.01, 0.01);
        fun0->SetParameter(2, 0.001);
        fun0->SetParLimits(3, -0.01, 0.01);
        fun0->SetLineColor(46);
        if (j!=0) {
            hVar->Fit("fun0", "M", "R", -0.07, 0.07);
            Double_t gamma = fun0->GetParameter(1);
            cout<<gamma<<endl;
            Double_t fwhm = 2*gamma;
            cout<<fwhm<<endl;
            text4->AddText(Form("FWHM: %.0f #mum", fwhm*10000));
        }

        c->SetLogy();
        hVar->Draw("HIST");

        c->SaveAs(folder + var[j] + ".png");
        hVar->Write();
        c->Clear();
        hVar->GetXaxis()->SetRangeUser(-0.07, 0.07);

//        TCanvas *c = new TCanvas("c1", "c1", 900, 1200);
        c->SetLogy();
        hVar->Draw();
        if (j!=0) text4->Draw("same");

        c->SaveAs(folder + var[j] + ".zoom.png");
        c->Clear();

        hVar->Reset();
    }

    c->SetLogy(0);
    hVar = new TH1F("gRefMult", "gRefMult", 60, 0, 60);
    ntp->Project("gRefMult", "grefMult", detCuts);
    hVar->GetXaxis()->SetTitle("gRefMult");
    hVar->GetYaxis()->SetTitle("N_{evts}");
    hVar->GetYaxis()->SetTitleOffset(0.8);
    hVar->SetFillColor(46);
    hVar->SetLineColor(46);
    hVar->SetFillStyle(3004);
    hVar->SetStats(0);
    hVar->SetTitle("");
    hVar->Draw("HIST");

    c->SaveAs(folder + "grefMult.png");

    fOut->Close();
    data->Close();
}
//____________________________________________________________________________________________________________________
void comp(){
    TFile* data1 = new TFile("res.ntp.PicoVertex.chi2Cut35.nhft0.root" ,"r");
    TFile* data2 = new TFile("res.ntp.PicoVertex.global.chi2cut35.nhft0.root" ,"r"); //this one looks better for x_diff

    TH1F *h1 = new TH1F();
    TH1F *h2 = new TH1F();

    TString var[] = {"x_diff", "y_diff", "z_diff", "ErrX_diff", "ErrY_diff", "ErrZ_diff"};

    for (int k = 0; k < 6; k++) {
        h1 = static_cast<TH1F*>(data1->Get(var[k]));
        h1->SetLineColor(46);
        h1->SetMarkerColor(46);
//        h1->Scale(1/h1->GetEntries());

        h2 = static_cast<TH1F*>(data2->Get(var[k]));
//        h2->Scale(1/h2->GetEntries());

        TCanvas *c = new TCanvas("c1", "c1", 900, 1200);
        c->SetLogy();
        h1->Draw();
        h2->Draw("same");
        Double_t intE;
        cout<<"h1 integral: "<<h1->IntegralAndError(h1->FindBin(-0.0005), h1->GetXaxis()->GetLast(), intE, "")<<endl;
        cout<<"h1 integral: "<<h1->Integral(0, 1000,"")<<endl;
        cout<<"h2 integral: "<<h2->Integral(h2->FindBin(-0.001), h2->GetXaxis()->GetLast(),"")<<endl;
//        cout<<"h2 integral: "<<h2->Integral(h2->FindBin(-0.001), h2->GetMaximumBin(),"")<<endl;

        TLegend *legend = new TLegend(0.6, 0.71, 0.75, 0.89);
        legend -> SetFillStyle(0);
        legend -> SetLineColor(0);
        legend -> SetTextSize(0.035);
        legend -> AddEntry(h1, "chi2Cut35.nhft0", "pl");
        legend -> AddEntry(h2, "global.chi2cut35.nhft0", "pl");
        legend -> Draw("same");
        c->SaveAs(var[k]+"_compare.png");
//        c->Close();
    }
}

//_______________________________________________________________________________________________________________________
void projectBins() {
//void projectBins(TH2F *hVar) {
//    TString input = "ntp.PicoVertex.chi2Cut35.nhft0.root";
    TString input = "ntp.PV.2003.staneling.root";
//    TString input = "ntp.PV.D0rem.global.2603.root";
    TFile* data = new TFile("res."+input ,"r");
    TFile* outData = new TFile("1d.err."+input ,"RECREATE");

    TH2F *h2d = new TH2F();
    TH1D *h1d = new TH1D();
    TH1D *h1dNeg = new TH1D();

    TString vertexType[] = {"picoDstV", "KFV"};
    TString varName[6];
    int ii = 0;
    for (int k = 0; k < 2; ++k) {
        varName[ii++] = Form("%sx:%sErrX", vertexType[k].Data(), vertexType[k].Data());
        varName[ii++] = Form("%sy:%sErrY", vertexType[k].Data(), vertexType[k].Data());
        varName[ii++] = Form("%sz:%sErrZ", vertexType[k].Data(), vertexType[k].Data());
    }

    Float_t positionBins[] = {0, 0.4, 0.8, 1.2, 1.6, 2, 2.4, 2.8, 3.2, 3.6, 4, 4.4, 4.8, 5.4, 6}; //this will be projected later from TH2
    int nBins = sizeof(positionBins)/ sizeof(Float_t);

    for (int j = 0; j < 6; ++j) {
        h2d = static_cast<TH2F*>(data->Get(varName[j]));
        cout<<varName[j]<<endl;

        for (int i = 0; i < nBins-1; ++i) {
            TString hisName = Form("%s_%.1f_%.1f", varName[j].Data(), positionBins[i], positionBins[i+1]);
            h1d = h2d->ProjectionX(hisName, h2d->GetYaxis()->FindBin(positionBins[i]), h2d->GetYaxis()->FindBin(positionBins[i+1]), "e"); //y-ova os ostava
            h1dNeg = h2d->ProjectionX("b", h2d->GetYaxis()->FindBin(-1*positionBins[i+1]), h2d->GetYaxis()->FindBin(-1*positionBins[i]), "e"); //y-ova os ostava
            h1d->Add(h1dNeg);
            h1d->SetNameTitle(hisName, hisName);
            outData->cd();
            h1d->Write();
            h1d->Reset();
        }
        h2d->Reset();
    }
    outData->Close();

    TFile* inData = new TFile("1d.err."+input ,"r");

    TH1D *h1dPico = new TH1D();
    TH1D *h1dKF = new TH1D();

    TString hisNamePico[3] = {"picoDstVx:picoDstVErrX", "picoDstVy:picoDstVErrY", "picoDstVz:picoDstVErrZ"};
    TString hisNameKF[3] = {"KFVx:KFVErrX", "KFVy:KFVErrY", "KFVz:KFVErrZ"};
    TString variable[3] = {"x", "y", "z"};

    for (int l = 0; l < 3; ++l) {
        TCanvas *cBig = new TCanvas("cBig", "cBig", 1800, 1600);
        cBig->Divide(5,3,0.005,0.005,0);

        for (int i = 0; i < nBins-1; ++i) {
            TString hisName = hisNamePico[l] + Form("_%.1f_%.1f", positionBins[i], positionBins[i + 1]);
            cout<<hisName<<endl;
            h1dPico = static_cast<TH1D*>(inData->Get(hisName));
            h1dPico->Scale(1/h1dPico->GetEntries());
            h1dPico->SetStats(0);
            h1dPico->Rebin(2);
            h1dPico->GetXaxis()->SetRangeUser(-0.03,0.6);
            h1dPico->SetFillColor(9);
            h1dPico->SetFillStyle(3005);
            if (l==2)  h1dPico->GetXaxis()->SetRangeUser(-0.03,0.3);
            h1dPico->GetXaxis()->SetTitle("#Delta [cm]");
            h1dPico->GetYaxis()->SetTitle("1/N");
            h1dPico->GetYaxis()->SetTitleOffset(0.7);
            h1dPico->GetXaxis()->SetTitleSize(0.05);
            h1dPico->GetYaxis()->SetTitleSize(0.05);
            h1dPico->GetXaxis()->SetLabelSize(0.05);
            h1dPico->GetYaxis()->SetLabelSize(0.05);

            h1dPico->SetTitle(Form("%s_%.1f_%.1f", variable[l].Data(),positionBins[i], positionBins[i + 1]));

            hisName = hisNameKF[l] + Form("_%.1f_%.1f", positionBins[i], positionBins[i + 1]);
            h1dKF = static_cast<TH1D*>(inData->Get(hisName));
            h1dKF->Scale(1/h1dKF->GetEntries());
            h1dKF->SetStats(0);
            h1dKF->Rebin(2);
            h1dKF->SetFillColor(46);
            h1dKF->SetFillStyle(3004);
            h1dKF->SetLineColor(46);

            if (h1dKF->GetMaximum()>h1dPico->GetMaximum()) h1dPico->GetYaxis()->SetRangeUser(0.001,1.2*h1dKF->GetMaximum());

//            TCanvas *c = new TCanvas("c1", "c1", 900, 1200);
            cBig->cd(i+1);
            cBig->cd(i+1)->SetLeftMargin(0.1);

            cBig->cd(i+1)->SetLogy();

            h1dPico->DrawCopy("HIST");
            h1dKF->DrawCopy("HIST same");

            TLegend *legend = new TLegend(0.6, 0.71, 0.75, 0.89);
            legend -> SetFillStyle(0);
            legend -> SetLineColor(0);
            legend -> SetTextSize(0.065);
            legend -> AddEntry(h1dPico, "picoDst", "pl");
            legend -> AddEntry(h1dKF, "KF", "pl");
            legend->Draw("same");
//            h1dPico->Reset();
//            h1dKF->Reset();

//            c->SaveAs(Form("%s_%.1f_%.1f.png", variable[l].Data(), positionBins[i], positionBins[i + 1]));
//            c->Close();
        }
//        cBig->SaveAs("%s.png");
        cBig->SaveAs(Form("%s.png", variable[l].Data()));
        cBig->Close();

    }
    inData->Close();
}
