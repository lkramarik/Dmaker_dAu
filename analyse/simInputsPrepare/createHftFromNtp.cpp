#include<iostream>
#include<fstream>
#include<vector>
#include"TH1.h"
#include"TH2.h"
#include"TF1.h"
#include"TCut.h"
#include"TMath.h"
#include"TCanvas.h"
#include"TFile.h"
#include"TLatex.h"
#include"TStyle.h"
#include"TNtuple.h"
#include"TPad.h"
#include"TLegend.h"
#include"TPaveText.h"
#include"TSystem.h"
#include"TAxis.h"
#include"TTree.h"
#include"TFitResultPtr.h"
#include"TFitResult.h"
#include"TString.h"

using namespace std;

void draw(TH1D* histo, TString xAxis, TString yAxis, Int_t color) {
    histo->GetXaxis()->SetTitle(xAxis);
    histo->GetYaxis()->SetTitle(yAxis);
    histo->SetLineColor(color);
    histo->SetMarkerColor(color);
    histo->SetStats(0);
}
//_______________________________________________________________________
void createHftFromNtp(){
    Int_t colors[] = {1,46,8,4,5,6,7,9,41,42,44,45,16};
//    TString fileNameData = "outputLocal.hists.root";
//    TString fileNameData = "ntp.track.2403.root";
//    TString fileNameData = "ntp.2403.sample.root";
    TString fileNameData = "/home/lukas/work/dmesons/Dmaker_ndAu/Dmaker_dAu/analyse/simInputsPrepare/toHadd/ntp.ratio.1404.half.root";
//    TString fileNameData = "/home/lukas/work/dmesons/Dmaker_ndAu/Dmaker_dAu/analyse/simInputsPrepare/toHadd/ntp.ratio.1404.half.rootprimary.root";
    auto* dataF = new TFile(fileNameData,"r");
    auto* dataOut = new TFile("ratio.root","recreate");
    auto* ntpData = (TNtuple*)dataF -> Get("ntp_tracks");

    Float_t dca, isHFT, zdc, isPrimary, nHitsFit, refMult, nHftTracks;

    double maxPt = 6;

    TH1D *pt_hist_HFT = new TH1D("pT_HFT", "pT_HFT", 100, 0.15, maxPt);
    TH1D *pt_hist = new TH1D("pt_hist", "pt_hist", 100, 0.15, maxPt);

//    ntp_tracks = new TNtuple("ntp_tracks","ntp_tracks", "pt:dca:isHFT:zdc:isPrimary:nHitsFit:refMult:nHftTracks");


    const char* cuts[] = {"dca<=0.1", "dca<=1.5", "dca<=1.0", "dca<=0.7", "dca<=0.5", "dca<=0.3"};
//    const char* cuts[] = {"dca<=1.5 && nHftTracks>1", "dca<=1.0 && nHftTracks>1", "dca<=0.7 && nHftTracks>1", "dca<=0.5 && nHftTracks>1", "dca<=0.3 && nHftTracks>1"};
//    const char* cuts[] = {"dca<=1.5 && nHftTracks>0", "dca<=1.0 && nHftTracks>0", "dca<=0.7 && nHftTracks>0", "dca<=0.5 && nHftTracks>0", "dca<=0.3 && nHftTracks>0",
//                          "dca<=1.5", "dca<=1.0", "dca<=0.7", "dca<=0.5", "dca<=0.3"};

    TString namesSuffix[] = {"dca0.1_TofMatched", "dca1.5_TofMatched", "dca1.0_TofMatched", "dca0.7_TofMatched", "dca0.5_TofMatched", "dca0.3_TofMatched"};
//    TString namesSuffix[] = {"dca1.5_TofMatched", "dca1.0_nHft1_TofMatched", "dca0.7_nHft1_TofMatched", "dca0.5_nHft1_TofMatched", "dca0.3_nHft1_TofMatched"};
//    TString namesSuffix[] = {"dca1.5_nHft1", "dca1.0_nHft1", "dca0.7_nHft1", "dca0.5_nHft1", "dca0.3_nHft1"};

//    TString namesSuffix[] = {"dca1.5_nHft0", "dca1.0_nHft0", "dca0.7_nHft0", "dca0.5_nHft0", "dca0.3_nHft0",
//                             "dca1.5", "dca1.0", "dca0.7", "dca0.5", "dca0.3"};

//    TCut additionalCut = "isTOF>0 && isPrimary>0";
    TCut additionalCut = "isPrimary>0 && nHftTracks>1";

    const int nBins = sizeof(namesSuffix) / sizeof(TString);

    TH1D *pt_hist_clone[nBins];

    TCanvas *c1 = new TCanvas("c1","c1",1200,900);
    TCanvas *cwork = new TCanvas("cwork","cwork",1200,900);
    TLegend *legend = new TLegend(0.155,0.789, 0.427, 0.884,"","brNDC");
    for (int i = 0; i < 6; ++i) {
        cout<<"nbin "<<i<<endl;
        TCut cut = cuts[i];
        cut+="pt>0.15";
        cut+="pt<6";
        cut+=additionalCut;
        cout<<cut<<endl;
        cwork->cd();
        ntpData -> Project("pt_hist", "pt", cut);
        cut+="isHFT>0";
        ntpData -> Project("pT_HFT", "pt", cut);

        pt_hist_clone[i] = (TH1D*)pt_hist_HFT->Clone(Form("HFT_matching_ratio_%i", i));
        pt_hist_clone[i]->SetNameTitle(Form("HFT_matching_ratio_%s", namesSuffix[i].Data()), cuts[i]);
        pt_hist_clone[i]->Sumw2();
        pt_hist_clone[i]->Divide(pt_hist);
        pt_hist_clone[i]->Write();
        pt_hist_HFT->Write(Form("pt_HFT_%i_%s", i, namesSuffix[i].Data()));
        pt_hist->Write(Form("pT_%i_%s", i, namesSuffix[i].Data()));
        draw(pt_hist_clone[i], "p_{T} [GeV/c]", "HFT ratio from data", colors[i]);
    }

    cout<<"draw"<<endl;

    for (int j = 0; j < 1; ++j) {
        c1->cd();
        if (j==0 )pt_hist_clone[j]->DrawCopy(); else pt_hist_clone[j]->DrawCopy("same");
        legend->AddEntry(pt_hist_clone[j], cuts[j], "pl");
    }

    legend->Draw("same");
    c1->SaveAs("hft/ratios.png");
}


void divide_ntp() {
    TString input="/home/lukas/work/dmesons/Dmaker_ndAu/Dmaker_dAu/analyse/simInputsPrepare/toHadd/ntp.ratio.1404.half.root";
    TFile* data = new TFile(input ,"r");
    TFile *fileOut = new TFile(input+"primary.root", "RECREATE");  // output root file

    TNtuple *ntpOut = new TNtuple("ntp_tracks","ntp_tracks", "pt:dca:isHFT:isTOF:zdc:nHitsFit:refMult:nHftTracks:nTofTracks:particleId");

    TNtuple *ntp = (TNtuple*) data->Get("ntp_tracks");
    Float_t  pt, dca, isHFT, isTOF, zdc, nHitsFit, refMult, nHftTracks, nTofTracks, particleId, isPrimary;

    ntp->SetBranchAddress("pt", &pt);
    ntp->SetBranchAddress("dca", &dca);
    ntp->SetBranchAddress("isHFT", &isHFT);
    ntp->SetBranchAddress("isTOF", &isTOF);
    ntp->SetBranchAddress("zdc", &zdc);
    ntp->SetBranchAddress("nHitsFit", &nHitsFit);
    ntp->SetBranchAddress("refMult", &refMult);
    ntp->SetBranchAddress("isPrimary", &isPrimary);
    ntp->SetBranchAddress("nHftTracks", &nHftTracks);
    ntp->SetBranchAddress("nTofTracks", &nTofTracks);
    ntp->SetBranchAddress("particleId", &particleId);


    const int nNtVars = ntpOut->GetNvar();
    float ntVar[nNtVars];

    cout<<"lets do this"<<endl;
    for (long int i = 0; i < ntp->GetEntries(); i++) {
        ntp->GetEntry(i);
        if (isPrimary>0) {
            ntpOut->Fill(pt,dca,isHFT,isTOF,zdc,nHitsFit,refMult,nHftTracks,nTofTracks,particleId);
        }
    }
    fileOut->cd();
    ntpOut->Write(ntpOut->GetName(), TObject::kOverwrite);
    fileOut->Close();

    data->Close();
    cout<<"Done."<<endl;
}