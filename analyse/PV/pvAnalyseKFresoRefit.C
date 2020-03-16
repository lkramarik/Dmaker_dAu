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

void pvAnalyseKFresoRefit() {
//    TString input = "ntp.KFVertex.global.2105.root";
//    TString input = "ntp.PVrefit.global.root";
    TString input = "ntp.KF.D0.meeting.hotspot.root";
//    TString input = "ntp.KFVertexReso.1605.root";
    TFile *data = new TFile(input, "r");
    TList* list = (TList*) data -> Get("picoD0AnaMaker");

    TFile *fOut = new TFile("res.refit." + input, "recreate");

    TH1F *hVar = new TH1F();

    TString hisName, varName;
    TCanvas *c = new TCanvas("c1", "c1", 900, 1200);
    c->SetGrid();

    TString var[] = {"X", "Y", "Z"};
//    Double_t fwhm = 0, sigma = 0;

    for (int j = 0; j < 3; ++j) { //variable to project
        TPaveText *text4 = new TPaveText(0.35, 0.848, 0.229, 0.873, "brNDC");
        text4->SetTextSize(0.05);
        text4->SetLineColor(0);
        text4->SetShadowColor(0);
        text4->SetFillColor(0);

        hVar = (TH1F*) list -> FindObject(Form("hPVDiff%s", var[j].Data()));
        hVar->Rebin(12);
        hVar->SetStats(0);
        hVar->Scale(1/hVar->GetEntries());
        hVar->SetFillColor(46);
        hVar->SetFillStyle(3004);
        hVar->SetMarkerStyle(20);
        hVar->SetLineColor(46);
        hVar->SetMarkerColor(46);
        fOut->cd();
        hVar->GetXaxis()->SetTitle(Form("#Delta%s [cm]", var[j].Data()));
        hVar->GetYaxis()->SetTitleOffset(1.3);
        hVar->GetYaxis()->SetTitle("1/N");
        TF1 *fun0 = new TF1("fun0", "[constant]/( pi*[gamma]*(1+(x-[mean])*(x-[mean])/([gamma]*[gamma])))+[jump]", -0.006, 0.006);
        fun0->SetParameter(0, 10000);
        fun0->SetParameter(1, 0.001);
        fun0->SetParameter(2, 5000);
        fun0->SetParameter(3, 0.0005);

        fun0->SetParLimits(1, -0.1, 0.1); //gamma
//        fun0->SetParLimits(2, -10000, 10000);
//        fun0->SetParameter(2, 0.001);
//        fun0->SetParLimits(3, -0.01, 0.01);
        fun0->SetLineColor(46);

//            hVar->Fit("gaus", "M", "R", -0.02, 0.02);
            hVar->Fit("fun0", "M", "R", -0.007, 0.007);
            Double_t gamma = fun0->GetParameter(1);
            cout<<gamma<<endl;
            Double_t fwhm = 2*gamma;
            cout<<fwhm<<endl;
            text4->AddText(Form("FWHM: %.0f #mum", fwhm*10000));
        c->SetLogy();
        hVar->Draw("HIST");

        c->SaveAs(var[j] + ".refit.png");
        hVar->Write();
        c->Clear();
        hVar->GetXaxis()->SetRangeUser(-0.005, 0.005);

//        TCanvas *c = new TCanvas("c1", "c1", 900, 1200);
        c->SetLogy();
        hVar->Draw();
        text4->Draw("same");
        c->SaveAs(var[j] + ".zoom.refit.png");

//        hVar->Reset();
    }

    fOut->Close();
    data->Close();
}