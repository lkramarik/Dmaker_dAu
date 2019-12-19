//
// Created by lukas on 3. 10. 2019.
//
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

void hotSpotProject(){
    TString input = "ntp.PV.1712.root";
    TString folder = "";
//    TString folder = "/media/lukas/376AD6A434B7392F/work/";
    TFile* outFile = new TFile("resH."+input ,"recreate");

    TFile* data = new TFile(folder+input ,"r");
    TNtuple* ntp = (TNtuple*)data -> Get("ntp_vertex");

    TH1F *hOut;// = new TH1F();
    TH1F *hIn;// = new TH1F();

//    TString variable[] = {"BBC", "ZDC", "refMult", "nPrimTracks", "nHftTrack", "picoDstVErrX", "picoDstVErrY", "picoDstVErrZ", ""};
    TString variable[] = {"BBC", "ZDC", "refMult", "nPrimTracks", "nHftTrack", "picoDstVErrX", "picoDstVErrY", "picoDstVErrZ", "runId"};
    Float_t limsMin[] = {0, 0, 0, 0, 0, 0, 0, 0, 17130000};
    Float_t limsMax[] = {1200, 200, 100, 100, 100, 0.5, 0.5, 0.5, 17150000};
    Float_t nBins[] = {1200, 200, 100, 100, 100, 100, 100, 100, 20000};
    TString in = "picoDstVy>-0.25 && picoDstVy<-0.16 && picoDstVx>-0.25 && picoDstVx<-0.16 && (nFlag0!=nPrimTracks)"; //hotspot
    TString out = "(picoDstVy<-0.25 || picoDstVy>-0.16 || picoDstVx<-0.25 || picoDstVx>-0.16) && (nFlag0!=nPrimTracks)"; //out of hotspot
    TString inSecondHot = "(picoDstVy<0.2 && picoDstVy>-0.1 && picoDstVx<0.2 && picoDstVx>-0.1)"; //out of hotspot
//    TString inSecondHot = "(picoDstVy<0.15 && picoDstVy>-0.1 && picoDstVx<0.15 && picoDstVx>-0.1)"; //out of hotspot

    //    TString axisName[] =     {};

    TCanvas *c = new TCanvas("c1", "c1", 900, 1200);
    c->SetGrid();

    cout<<"Number of stuff in the second spot: "<<ntp->GetEntries(inSecondHot)<<endl;

    for (int k = 0; k < 9; ++k) {
        hOut = new TH1F("hOut", variable[k], nBins[k], limsMin[k], limsMax[k]);
        hIn = new TH1F("hIn", variable[k], nBins[k], limsMin[k], limsMax[k]);
        cout<<variable[k]<<endl;
//        ntp->Project("hOut", variable[k], out);
        ntp->Project("hOut", variable[k], inSecondHot);
        ntp->Project("hIn", variable[k], in);

        hIn->Scale(1/hIn->GetEntries());
        hIn->SetStats(0);
//        hIn->Rebin(4);
        hIn->SetFillColor(9);
        hIn->SetFillStyle(3005);
        hIn->SetLineColor(9);
        hIn->SetTitle("");

        hOut->Scale(1/hOut->GetEntries());
        hOut->SetStats(0);
//        hOut->Rebin(4);
        hOut->SetFillColor(46);
        hOut->SetFillStyle(3004);
        hOut->SetLineColor(46);
        hOut->GetYaxis()->SetTitle("1/N");
        hOut->GetYaxis()->SetTitleOffset(0.8);
        hOut->GetXaxis()->SetTitle(variable[k]);
        hOut->SetTitle("");

        hOut->Draw("HIST");
        hIn->Draw("HIST same");
        TLegend *legend = new TLegend(0.6, 0.71, 0.75, 0.89);
        legend -> SetFillStyle(0);
        legend -> SetLineColor(0);
        legend -> SetTextSize(0.05);
        legend -> AddEntry(hIn, "in hotspot", "f");
        legend -> AddEntry(hOut, "out of hotspot", "f");
        legend->Draw("same");
        c->SaveAs("hotspot/"+variable[k]+".png");
        outFile->cd();
        hOut->Write(variable[k]+"_out");
        hIn->Write(variable[k]+"_in");
        c->Write();
        c->Clear();
        delete hOut;
        delete hIn;
    }
}