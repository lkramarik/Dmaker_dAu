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
#include"TProfile.h"
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
//    TString fileNameData = "/home/lukas/work/dmesons/Dmaker_ndAu/Dmaker_dAu/analyse/simInputsPrepare/toHadd/ntp.ratio.1404.half.root";
    TString fileNameData = "/media/lukas/1E4183350724EC9C/0.hists.root";
//    TString fileNameData = "/media/lukas/1E4183350724EC9C/tracks.sim.vzvpd3.hotspot.root";
//    TString fileNameData = "/home/lukas/work/dmesons/Dmaker_ndAu/Dmaker_dAu/analyse/simInputsPrepare/toHadd/ntp.ratio.1404.half.rootprimary.root";
    auto* dataF = new TFile(fileNameData,"r");
    auto* dataOut = new TFile("ratio.root","recreate");
    auto* ntpData = (TNtuple*)dataF -> Get("ntp_tracks");

    Float_t  pt, dca, isHFT, isTOF, zdc, nHitsFit, refMult, nHftTracks, nTofTracks, particleId, isPrimary, dcaxy, dcaz, runId, eventId;
    ntpData->SetBranchAddress("runId", &runId);
    ntpData->SetBranchAddress("eventId", &eventId);
    ntpData->SetBranchAddress("pt", &pt);
    ntpData->SetBranchAddress("dca", &dca);
    ntpData->SetBranchAddress("dcaXy", &dcaxy);
    ntpData->SetBranchAddress("dcaZ", &dcaz);
    ntpData->SetBranchAddress("isHFT", &isHFT);
    ntpData->SetBranchAddress("isTOF", &isTOF);
    ntpData->SetBranchAddress("zdc", &zdc);
    ntpData->SetBranchAddress("nHitsFit", &nHitsFit);
    ntpData->SetBranchAddress("refMult", &refMult);
    ntpData->SetBranchAddress("isPrimary", &isPrimary);
    ntpData->SetBranchAddress("nHftTracks", &nHftTracks);
    ntpData->SetBranchAddress("nTofTracks", &nTofTracks);
    ntpData->SetBranchAddress("particleId", &particleId);

//    Float_t dcaCuts[]={0.1, 0.3, 0.5, 0.7, 1, 1.5};
//    const int nDCAcuts = sizeof(dcaCuts) / sizeof(Float_t);

    const int nDCAcuts=6;
    const float dcaCuts[nDCAcuts+1]={0., 0.3, 0.5, 0.7, 1.0, 1.2, 1.5};

    TString trackKinds[]=    {"","primary",     "TOFmatched", "primary_TOFmatched", "_nHftTracks1", "_nHftTracks2","_nTofTracks2"};
    TString trackKindsCuts[]={"","isPrimary>0", "isTOF>0",    "isPrimary>0 && isTOF>0"};
    const int nKinds = sizeof(trackKinds) / sizeof(TString);

    TProfile *prRatio[nKinds][nDCAcuts];
//    TH2F *prRatio[nKinds][nDCAcuts];

    TH1D* hPt[2][nKinds][nDCAcuts];
    TH1D* hPtHft[2][nKinds][nDCAcuts];

    TString nameHisto;
    for (int i = 0; i < nKinds; ++i) {
        for (int j = 0; j < nDCAcuts; ++j) {
            for (int k = 0; k < 2; ++k) {
                nameHisto = Form("hPt_dca%.1f_%.1f_%s_%i", dcaCuts[j], dcaCuts[j + 1], trackKinds[i].Data(), k);
                hPt[k][i][j] = new TH1D(nameHisto, nameHisto, 100, 0, 10);
                hPt[k][i][j]->Sumw2();
                nameHisto = Form("hPtHft_dca%.1f_%.1f_%s_%i", dcaCuts[j], dcaCuts[j + 1], trackKinds[i].Data(), k);
                hPtHft[k][i][j] = new TH1D(nameHisto, nameHisto, 100, 0, 10);
                hPtHft[k][i][j]->Sumw2();
            }
            nameHisto = Form("hPtHft_dca%.1f_%.1f_%s", dcaCuts[j], dcaCuts[j + 1], trackKinds[i].Data());
            prRatio[i][j] = new TProfile(nameHisto,nameHisto,100,0,10,0,20);
//            prRatio[i][j] = new TH2F(nameHisto,nameHisto,100,0.,10.,100,0.,2.);
            prRatio[i][j]->Sumw2();
        }
    }

    Float_t previousEvent=-999;
    Float_t previousRun=-999;
    int nEvts=0;
    cout<<"total entr "<<ntpData->GetEntries()<<endl;

    for (int k = 0; k < ntpData->GetEntries(); ++k) {
//    for (int k = 0; k < ntpData->GetEntries()/100; ++k) {
        ntpData->GetEntry(k);
        if (particleId!=1) continue;
        if (previousEvent!=eventId || previousRun!=runId) { //new event
            previousEvent = eventId;
            previousRun = runId;
            nEvts++;
            for (int i = 0; i < nKinds; ++i) {
                for (int j = 0; j < nDCAcuts; ++j) {
                    hPtHft[1][i][j]->Divide(hPt[1][i][j]);
                    for (int l = 1; l < hPt[1][i][j]->GetNbinsX()+1; ++l) {
                        prRatio[i][j]->Fill(hPtHft[1][i][j]->GetBinCenter(l), hPtHft[1][i][j]->GetBinContent(l));
                    }
                    hPt[1][i][j]->Reset();
                    hPtHft[1][i][j]->Reset();
                }
            }
        }

        for (int i = 0; i < nDCAcuts; ++i) {
            if (dca>dcaCuts[i] && dca<dcaCuts[i+1]) {
                for (int j = 0; j < 2; ++j) { //total vs. evt averaging
                    hPt[j][0][i]->Fill(pt);
                    if (isPrimary > 0) hPt[j][1][i]->Fill(pt);
                    if (isTOF > 0) hPt[j][2][i]->Fill(pt);
                    if (isPrimary > 0 && isTOF > 0) hPt[j][3][i]->Fill(pt);
                    if (isPrimary > 0 && isTOF > 0 && nHftTracks > 1) hPt[j][4][i]->Fill(pt);
                    if (isPrimary > 0 && isTOF > 0 && nHftTracks > 2) hPt[j][5][i]->Fill(pt);
                    if (isPrimary > 0 && isTOF > 0 && nTofTracks > 2) hPt[j][6][i]->Fill(pt);

                    if (isHFT > 0) {
                        hPtHft[j][0][i]->Fill(pt);
                        if (isPrimary > 0) hPtHft[j][1][i]->Fill(pt);
                        if (isTOF > 0) hPtHft[j][2][i]->Fill(pt);
                        if (isPrimary > 0 && isTOF > 0) hPtHft[j][3][i]->Fill(pt);
                        if (isPrimary > 0 && isTOF > 0 && nHftTracks > 1) hPtHft[j][4][i]->Fill(pt);
                        if (isPrimary > 0 && isTOF > 0 && nHftTracks > 2) hPtHft[j][5][i]->Fill(pt);
                        if (isPrimary > 0 && isTOF > 0 && nTofTracks > 2) hPtHft[j][6][i]->Fill(pt);
                    }
                }
            }
        }
    }

    cout<<"nevents "<<nEvts<<endl;
    for (int i = 0; i < nKinds; ++i) {
        for (int j = 0; j < nDCAcuts-2; ++j) {
            hPtHft[0][i][j]->Divide(hPt[0][i][j]);
            cout<<hPtHft[0][i][j]->Integral()/100<<endl;
            dataOut->cd();
            hPt[0][i][j]->Write();
            hPtHft[0][i][j]->Write();
            prRatio[i][j]->Write();
        }
    }
    dataOut->Close();

    /*
    const char* cuts[] = {"dca<=0.1", "dca<=1.5", "dca<=1.0", "dca<=0.7", "dca<=0.5", "dca<=0.3"};
//    const char* cuts[] = {"dca<=1.5 && nHftTracks>1", "dca<=1.0 && nHftTracks>1", "dca<=0.7 && nHftTracks>1", "dca<=0.5 && nHftTracks>1", "dca<=0.3 && nHftTracks>1"};
//    const char* cuts[] = {"dca<=1.5 && nHftTracks>0", "dca<=1.0 && nHftTracks>0", "dca<=0.7 && nHftTracks>0", "dca<=0.5 && nHftTracks>0", "dca<=0.3 && nHftTracks>0",
//                          "dca<=1.5", "dca<=1.0", "dca<=0.7", "dca<=0.5", "dca<=0.3"};

    TString namesSuffix[] = {"dca0.1_TofMatched", "dca1.5_TofMatched", "dca1.0_TofMatched", "dca0.7_TofMatched", "dca0.5_TofMatched", "dca0.3_TofMatched"};
//    TString namesSuffix[] = {"dca1.5_TofMatched", "dca1.0_nHft1_TofMatched", "dca0.7_nHft1_TofMatched", "dca0.5_nHft1_TofMatched", "dca0.3_nHft1_TofMatched"};
//    TString namesSuffix[] = {"dca1.5_nHft1", "dca1.0_nHft1", "dca0.7_nHft1", "dca0.5_nHft1", "dca0.3_nHft1"};

//    TString namesSuffix[] = {"dca1.5_nHft0", "dca1.0_nHft0", "dca0.7_nHft0", "dca0.5_nHft0", "dca0.3_nHft0",
//                             "dca1.5", "dca1.0", "dca0.7", "dca0.5", "dca0.3"};

    TCut additionalCut = "isTOF>0 && isPrimary>0";
//    TCut additionalCut = "isPrimary>0 && nHftTracks>1";

    const int nBins = sizeof(namesSuffix) / sizeof(TString);

    TH1D *pt_hist_clone[nBins];

    TCanvas *c1 = new TCanvas("c1","c1",1200,900);
    TCanvas *cwork = new TCanvas("cwork","cwork",1200,900);
    TLegend *legend = new TLegend(0.155,0.789, 0.427, 0.884,"","brNDC");
    for (int i = 0; i < 6; ++i) {
        cout<<"nbin "<<i<<endl;
        TCut cut = cuts[i];
//        cut+="pt>0.15";
//        cut+="pt<6";
        cut+=additionalCut;
        cout<<cut<<endl;
        cwork->cd();
        ntpData -> Project("pt_hist", "pt", cut);
        cut+="isHFT>0";
        ntpData -> Project("pT_HFT", "pt", cut);

        pt_hist_clone[i] = (TH1D*)pt_hist_HFT->Clone();
        pt_hist_clone[i]->SetNameTitle(Form("hPt_%s", namesSuffix[i].Data()), cuts[i]);
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
    */
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