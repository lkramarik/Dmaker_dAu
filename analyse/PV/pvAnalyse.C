#include "TCut.h"
#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TStyle.h"
#include "TPaveStats.h"
#include <iostream>

using namespace std;

void pvAnalyse(){
//    TString input = "ntp.PicoVertex.chi2Cut35.nhft0.root";
//    TString input = "ntp.PicoVertex.global.chi2cut35.nhft0.root";
//    TString input = "ntp.PV.2003.staneling.root";
    TString input = "ntp.PV.1712.root";
//    TString input = "ntp.PV.0310.root";

//    TString input = "ntp.PV.2003.staneling.root";
//    TString input = "ntp.PV.D0rem.global.2603.root";
    TString folder = "";
//    TString folder = "/media/lukas/376AD6A434B7392F/work/";

    TFile* data = new TFile(folder+input ,"r");
    TNtuple* ntp = (TNtuple*)data -> Get("ntp_vertex");
    TH1F* hVarL[] = {new TH1F(), new TH1F()};
    TH2F* hVar = new TH2F();

    TH1F* hDiff = new TH1F();
    TFile *fOut = new TFile("res."+input,"recreate");

    bool draw= false;
    Int_t col[] = {46, 8};

    TString vertexType[] = {"picoDstV", "KFV"};
    TString var[] = {"x", "y", "z", "ErrX", "ErrY", "ErrZ"};

    Float_t limsMin[] = {-6, -6, -6, -0.01, -0.01, -0.01};
    Float_t limsMax[] = {6, 6, 6, 1, 1, 1};


//    TCut* detCuts = new TCut("nGlobTracks>0 && abs(picoDstVz)<5");
//    TCut* detCuts = new TCut("nGlobTracks>0 && nHftTracks>1.5 && abs(picoDstVx)<0.5 && abs(picoDstVy)<0.5");
//    TString detCuts = "nGlobTracks>0 && nHftTracks>1.5 && abs(picoDstVx)<0.5 && abs(picoDstVy)<0.5 && sqrt(picoDstVx*picoDstVx+picoDstVy*picoDstVy)<0.5";
//    TString detCuts = "nHftTracks>1.5 && sqrt(picoDstVx*picoDstVx+picoDstVy*picoDstVy)<0.5 && picoDstVy>-0.25 && picoDstVy<-0.16 && picoDstVx>-0.25 && picoDstVx<-0.16"; //hotspot
    TString detCuts = "nHftTracks>0 && picoDstVy>-0.25 && picoDstVy<-0.16 && picoDstVx>-0.25 && picoDstVx<-0.16"; //hotspot
//    TString detCuts = "nHftTracks>0"; //all
//    TString detCuts = "nHftTracks>0 && (picoDstVy<-0.25 || picoDstVy>-0.16 || picoDstVx<-0.25 || picoDstVx>-0.16)"; //out of hotspot
//    TString detCuts = "nGlobTracks>0 && nHftTracks>1.5";
    TString hisName, varName;
    TCut cut1, cut2;


    for (int j = 0; j < 6; ++j) { //variable to project
        for (int k = 0; k < 2; ++k) { //pico and KF
            varName = Form("%s%s", vertexType[k].Data(), var[j].Data());
            hisName = varName;
            hVarL[k] = new TH1F(hisName, hisName, 1000, limsMin[j], limsMax[j]);
            ntp->Project(hisName, varName, detCuts);
        }


        hisName = Form("%s_diff", var[j].Data());
        varName = Form("%s%s-%s%s", vertexType[0].Data(), var[j].Data(), vertexType[1].Data(), var[j].Data());
        hDiff = new TH1F(hisName, hisName, 1000, -2, 2);
        ntp->Project(hisName, varName, detCuts);

        varName = Form("%s%s:%s%s", vertexType[0].Data(), var[j].Data(), vertexType[1].Data(), var[j].Data());
        hisName = varName;
        TH2F *hCorrelation = new TH2F(hisName, varName, 2000, limsMin[j], limsMax[j], 2000, limsMin[j], limsMax[j]);
        ntp->Project(hisName, varName, detCuts);

        fOut->cd();

        hVarL[0]->Write();
        hVarL[1]->Write();
        hDiff->Write();
        hCorrelation->Write();

        hVarL[0]->Reset();
        hVarL[1]->Reset();
        hDiff->Reset();
        hCorrelation->Reset();
    }

    varName = "sqrt(picoDstVx*picoDstVx+picoDstVy*picoDstVy)";
    hVarL[0] = new TH1F("picoDstr", "picoDstr", 1000, -0.01, 2);
    ntp->Project("picoDstr", varName, detCuts);
    fOut->cd();
    hVarL[0]->Write();

    varName = "sqrt(KFVx*KFVx+KFVy*KFVy)";
    hVarL[0] = new TH1F("KFr", "KFr", 1000, -0.01, 2);
    ntp->Project("KFr", varName, detCuts);
    fOut->cd();
    hVarL[0]->Write();

    varName = "abs(sqrt(KFVx*KFVx+KFVy*KFVy)-sqrt(picoDstVx*picoDstVx+picoDstVy*picoDstVy))";
    hVarL[0] = new TH1F("r_diff", "r_diff", 1000, -0.001, 2);
    ntp->Project("r_diff", varName, detCuts);
    fOut->cd();
    hVarL[0]->Write();

    varName = "abs(sqrt(KFVx*KFVx+KFVy*KFVy+KFVz*KFVz)-sqrt(picoDstVx*picoDstVx+picoDstVy*picoDstVy+picoDstVz*picoDstVz))";
    hVarL[0] = new TH1F("r3D_diff", "r3D_diff", 1000, -0.001, 2);
    ntp->Project("r3D_diff", varName, detCuts);
    fOut->cd();
    hVarL[0]->Write();



    TH2F* hHotSpots = new TH2F("hHotSpots", "hHotSpots", 2000, -6, 6, 2000, -6, 6);
    ntp->Project("hHotSpots", "picoDstVx:picoDstVy", detCuts);
    hHotSpots->Write();

    TH2F* hPVxVsZ = new TH2F("hPVxVsZ", "hPVxVsZ", 2000, -6, 6, 2000, -6, 6);
    ntp->Project("hPVxVsZ", "picoDstVx:picoDstVz", detCuts);
    hPVxVsZ->Write();

    TH2F* hPVyVsZ = new TH2F("hPVyVsZ", "hPVyVsZ", 2000, -6, 6, 2000, -6, 6);
    ntp->Project("hPVyVsZ", "picoDstVy:picoDstVz", detCuts);
    hPVyVsZ->Write();

    TH2F* hPVerrXVsRefMult = new TH2F("picoDstVErrX_vs_refMult", "picoDstVErrX_vs_refMult", 200, 0, 200, 1000, 0, 1);
    ntp->Project("picoDstVErrX_vs_refMult", "picoDstVErrX:refMult", detCuts);
    hPVerrXVsRefMult->Write();

    TH2F* hPVerrZVsRefMult = new TH2F("picoDstVErrZ_vs_refMult", "picoDstVErrZ_vs_refMult", 200, 0, 200, 1000, 0, 1);
    ntp->Project("picoDstVErrZ_vs_refMult", "picoDstVErrZ:refMult", detCuts);
    hPVerrZVsRefMult->Write();

    varName = "sqrt(picoDstVx**2+picoDstVy**2):sqrt(KFVx**2+KFVy**2)";
    TH2F* hdR = new TH2F("hdRCorrelation", varName, 2000, 0, 6, 2000, 0, 6);
    ntp->Project("hdRCorrelation", varName, detCuts);
    hdR->Write();




    //err vs value
//    for (int l = 0; l < 2; ++l) {
//        varName = Form("%sx:%sErrX", vertexType[l].Data(), vertexType[l].Data());
//        hVar = new TH2F(varName, varName, 2000, 0, 6, 2000, -6, 6);
//        ntp->Project(varName, varName, detCuts);
//        hVar->Write();
//        hVar->Reset();
//
//        varName = Form("%sy:%sErrY", vertexType[l].Data(), vertexType[l].Data());
//        hVar = new TH2F(varName, varName, 2000, 0, 6, 2000, -6, 6);
//        ntp->Project(varName, varName, detCuts);
//        hVar->Write();
//        hVar->Reset();
//
//        varName = Form("%sz:%sErrZ", vertexType[l].Data(), vertexType[l].Data());
//        hVar = new TH2F(varName, varName, 2000, 0, 6, 2000, -6, 6);
//        ntp->Project(varName, varName, detCuts);
//        hVar->Write();
//        hVar->Reset();
//    }

//    TH1F* hRuns = new TH1F("hRuns", "hRuns", 10000, 17100000.5, 17200000.5);
//    ntp->Project("hRuns", "runId", detCuts);

//    TCut* runCut = new TCut("(picoDstVx>-0.05 && picoDstVx<0.1)");
//    TH1F* hRunsCut = new TH1F("hRunsCut", "hRunsCut", 10000, 17100000.5, 17200000.5);
//    ntp->Project("hRunsCut", "runId", *runCut);
//
//    hRunsCut->Divide(hRuns);
//    fOut->cd();
//
//    hRuns->Write();
//    hRunsCut->Write();
//
//    TH1F* hRunOne = new TH1F("hRunOne", "hRunOne", 1000, -5, 5);
//    ntp->Project("hRunOne", "picoDstVx", "runId<17134980");
//    hRunOne->Write();

    fOut->Close();
    data->Close();
}

//____________________________________________________________________________________________________________________
void plotStuff() {
    TString input = "ntp.PV.0310.root";
    TFile* data = new TFile("res."+input ,"r");

    TH1F *hKF = new TH1F();
    TH1F *hPico = new TH1F();

    TH2F *hKF2D = new TH2F();
    TH2F *hPico2D = new TH2F();

    TString vertexType[] = {"picoDst", "KF", "picoDstVErr"};

    TString hisName[] = {"ErrX", "ErrY", "ErrZ", "x", "y", "z", "r"};
//    TString axisName[] =     {};
    TString variableKF[] =   {"KFVErrX",      "KFVErrY",      "KFVErrZ",      "KFVx",      "KFVy",      "KFVz",     "KFr"};
    TString variablePico[] = {"picoDstVErrX", "picoDstVErrY", "picoDstVErrZ", "picoDstVx", "picoDstVy", "picoDstVz","picoDstr"};

    TCanvas *c = new TCanvas("c1", "c1", 900, 1200);
    c->SetLogy();
    c->SetLogx();

    c->SetGrid();
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.15);

    for (int j = 0; j < 6; ++j) {
        TString hisNameGet = variablePico[j]+":"+variableKF[j];
        hKF2D = static_cast<TH2F*>(data->Get(hisNameGet));
        hKF2D->SetStats(0);
//        hKF2D->Rebin2D(2,2,"");
        if (j<3) {
            hKF2D->GetYaxis()->SetRangeUser(0.001,0.2);
            hKF2D->GetXaxis()->SetRangeUser(0.001,0.2);
        }
        hKF2D->GetYaxis()->SetTitle(vertexType[0]+" "+hisName[j]+" [cm]");
        hKF2D->GetXaxis()->SetTitle(vertexType[1]+" "+hisName[j]+" [cm]");
        hKF2D->GetXaxis()->SetTitleOffset(1.1);
        hKF2D->Draw("colz");
        c->SaveAs(hisNameGet+".png");
        c->Clear();
    }
    c->SetLogx(0);

    gPad->SetRightMargin(0.1);

    for (int k = 0; k < 7; ++k) {
        hKF = static_cast<TH1F*>(data->Get(variableKF[k]));
        hPico = static_cast<TH1F*>(data->Get(variablePico[k]));
        hPico->Scale(1/hPico->GetEntries());
        hPico->SetStats(0);
        hPico->Rebin(4);
        hPico->SetFillColor(9);
        hPico->SetFillStyle(3005);
        hPico->SetLineColor(9);
        hPico->SetTitle("");


        hKF->Scale(1/hPico->GetEntries());
        hKF->SetStats(0);
        hKF->Rebin(4);
        hKF->SetFillColor(46);
        hKF->SetFillStyle(3004);
        hKF->SetLineColor(46);
        hKF->GetYaxis()->SetTitle("1/N");
        hKF->GetYaxis()->SetTitleOffset(0.8);
        hKF->GetXaxis()->SetTitle(hisName[k]+" [cm]");
        hKF->SetTitle("");

        hKF->Draw("HIST");
        hPico->Draw("HIST same");
        TLegend *legend = new TLegend(0.6, 0.71, 0.75, 0.89);
        legend -> SetFillStyle(0);
        legend -> SetLineColor(0);
        legend -> SetTextSize(0.05);
        legend -> AddEntry(hPico, "picoDst", "f");
        legend -> AddEntry(hKF, "KF", "f");
        legend->Draw("same");
        c->SaveAs(hisName[k]+".png");
        c->Clear();
    }

    TString var[] = {"x_diff", "y_diff", "z_diff", "r_diff"};
    TString varTitle[] = {"x", "y", "z", "r"};
    for (int i = 0; i < 4 ; ++i) {
        TPaveText *text4 = new TPaveText(0.35, 0.848, 0.229, 0.873, "brNDC");
        text4->SetTextSize(0.03);
        text4->SetLineColor(0);
        text4->SetShadowColor(0);
        text4->SetFillColor(0);
//        gStyle->SetOptStat("M");
        hPico = static_cast<TH1F*>(data->Get(var[i]));
        hPico->Scale(1/hPico->GetEntries());
//        hPico->SetStats(0);
//        hPico->Rebin(2);
        hPico->SetFillColor(9);
        hPico->SetFillStyle(3005);
        hPico->SetMarkerColor(9);
        hPico->SetMarkerStyle(20);
        hPico->SetLineColor(9);
        hPico->SetTitle("");
        hPico->SetStats(0);
        hPico->GetYaxis()->SetTitle("1/N");
        hPico->GetXaxis()->SetRangeUser(-0.07,0.07);
        hPico->GetYaxis()->SetTitleOffset(0.8);
        hPico->GetXaxis()->SetTitle("#Delta"+varTitle[i]+"(Pico,KF) [cm]");

        if (i!=3) {
            TF1 *fun0 = new TF1("fun0", "[constant]/( pi*[gamma]*(1+(x-[mean])*(x-[mean])/([gamma]*[gamma])))+[jump]", -0.05, 0.05);
            fun0->SetLineColor(9);
            fun0->SetParameter(1, 0.005);
            fun0->SetParLimits(1, 0.0001, 1); //gamma
            fun0->SetParLimits(2, -0.01, 0.01);
            fun0->SetParameter(2, 0.001);
            fun0->SetParLimits(3, -0.01, 0.01);
            hPico->Fit("fun0", "M", "R", -0.07, 0.07);
            Double_t gamma = fun0->GetParameter(1);
            cout<<gamma<<endl;
            Double_t fwhm = 2*gamma;
            cout<<fwhm<<endl;
            text4->AddText(Form("FWHM: %.0f #mum", fwhm*10000));
        }

        hPico->Draw();
        if (i!=3) text4->Draw("same");
//        hPico->Draw("HIST");
//        gPad->Update();
//        TPaveStats *st = (TPaveStats*)hPico->FindObject("stats");
//        st->SetX1NDC(0.5); //new x start position
//        st->SetX2NDC(0.9); //new x end position
//        st->SetY1NDC(0.8); //new x start position
//        st->SetY2NDC(0.9); //new x end position
        c->SaveAs(var[i]+".png");
        c->Clear();
    }


    c->Close();
}


//____________________________________________________________________________________________________________________
void comp(){
    TFile* data1 = new TFile("res.ntp.PV.2003.staneling.root" ,"r");
    TFile* data2 = new TFile("res.ntp.PV.2003.staneling.outHotSpot.root" ,"r"); //this one looks better for x_diff

    TH1F *h1 = new TH1F();
    TH1F *h2 = new TH1F();

    TString var[] = {"picoDstVErrX", "picoDstVErrY", "picoDstVErrZ", "ErrX_diff", "ErrY_diff", "ErrZ_diff"};

    for (int k = 0; k < 6; k++) {
        h1 = static_cast<TH1F*>(data1->Get(var[k]));
        h1->SetLineColor(46);
        h1->SetMarkerColor(46);
        h1->Rebin(4);
        h1->Scale(1/h1->GetEntries());

        h2 = static_cast<TH1F*>(data2->Get(var[k]));
        h2->Rebin(4);
        h2->Scale(1/h2->GetEntries());

        TCanvas *c = new TCanvas("c1", "c1", 900, 1200);
        c->SetLogy();
        h1->Draw("HIST");
        h2->Draw("HIST same");
        Double_t intE;
        cout<<"h1 integral: "<<h1->IntegralAndError(h1->FindBin(-0.0005), h1->GetXaxis()->GetLast(), intE, "")<<endl;
        cout<<"h1 integral: "<<h1->Integral(0, 1000,"")<<endl;
        cout<<"h2 integral: "<<h2->Integral(h2->FindBin(-0.001), h2->GetXaxis()->GetLast(),"")<<endl;
//        cout<<"h2 integral: "<<h2->Integral(h2->FindBin(-0.001), h2->GetMaximumBin(),"")<<endl;

        TLegend *legend = new TLegend(0.6, 0.71, 0.75, 0.89);
        legend -> SetFillStyle(0);
        legend -> SetLineColor(0);
        legend -> SetTextSize(0.035);
        legend -> AddEntry(h1, "hotspot", "pl");
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
//    TFile* data = new TFile("res."+input ,"r");
    TFile* data = new TFile("res.ntp.PV.2003.staneling.root","r");
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

    Float_t positionBins[] = {0, 0.5001, 1.2, 1.6, 2, 2.4, 2.8, 3.2, 3.6, 4, 4.4, 4.8, 5.4, 6}; //this will be projected later from TH2
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
            h1dPico->GetXaxis()->SetRangeUser(-0.03,0.2);
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
