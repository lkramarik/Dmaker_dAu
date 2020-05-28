#include"TSystem.h"
#include"TAxis.h"
#include"TTree.h"
#include"TCanvas.h"
#include"TFile.h"
#include"TNtuple.h"
#include"TH1.h"
#include"TLegend.h"
#include"TPad.h"
#include"TFitResultPtr.h"
#include"TFitResult.h"
#include"TString.h"
#include<iostream>
#include<fstream>
#include<vector>
using namespace std;

const int nDcaRatio=6;
const float dcaCut[nDcaRatio]={0.3, 0.5, 0.7, 1.0, 1.2, 1.5};

const int nPtDca=6;
const float ptCut[nPtDca+1]={0.15, 0.3, 0.7, 1.2, 2, 3, 5};

//_____________________________________________________________________________________________________________
void draw(TH1D* histo, TString xAxis, TString yAxis, Int_t color) {
    histo->GetXaxis()->SetTitle(xAxis);
    histo->GetYaxis()->SetTitle(yAxis);
    histo->SetLineColor(color);
    histo->SetMarkerColor(color);
    histo->SetMarkerStyle(20);
    histo->SetStats(0);
}

//_____________________________________________________________________________________________________________
void createQATracks() {
    Int_t colors[] = {1,46,8,4,5,6,7,9,41,42,44,45,16};
//    TString fileNameData = "outputLocal.hists.root";
//    TString fileNameData = "ntp.track.2403.root";
//    TString fileNameData = "ntp.2403.sample.root";
    TString fileNameData = "/home/lukas/work/dmesons/Dmaker_ndAu/Dmaker_dAu/analyse/simInputsPrepare/toHadd/ntp.ratio.1404.half.rootprimary.root";
//    TString fileNameData = "/home/lukas/work/dmesons/Dmaker_ndAu/Dmaker_dAu/analyse/simInputsPrepare/toHadd/ntp.ratio.1404.half.root";
    auto* dataF = new TFile(fileNameData,"r");
    auto* dataOut = new TFile("tracksQA.data.root","recreate");
    auto* ntpData = (TNtuple*)dataF -> Get("ntp_tracks");

    double maxPt = 6;

    TString variables[] = {"dca", "nHitsFit"};
    Double_t limMin[] = {0, 0.5};
    Double_t limMax[] = {1.5, 50.5};
    Int_t nBins[] = {100, 51};
    TH1D *his;
    TH1D *hisHft;

    Float_t  pt, dca, isHFT, isTOF, zdc, nHitsFit, refMult, nHftTracks, nTofTracks, particleId, isPrimary;
    ntpData->SetBranchAddress("pt", &pt);
    ntpData->SetBranchAddress("dca", &dca);
    ntpData->SetBranchAddress("isHFT", &isHFT);
    ntpData->SetBranchAddress("isTOF", &isTOF);
    ntpData->SetBranchAddress("zdc", &zdc);
    ntpData->SetBranchAddress("nHitsFit", &nHitsFit);
    ntpData->SetBranchAddress("refMult", &refMult);
//    ntpData->SetBranchAddress("isPrimary", &isPrimary);
    ntpData->SetBranchAddress("nHftTracks", &nHftTracks);
    ntpData->SetBranchAddress("nTofTracks", &nTofTracks);
    ntpData->SetBranchAddress("particleId", &particleId);

    Double_t maxDca=0.2;

    //init histos
    TString histoname;
    TH1D* hnHitsFit = new TH1D("hnHitsFit", "hnHitsFit", 51, 0.5, 50.5);
    TH1D* hnHitsFitHft = new TH1D("hnHitsFitHft", "hnHitsFitHft", 51, 0.5, 50.5);
    TH1D* hDcaAll = new TH1D("hDca_all", "hDca_all", 100, 0, maxDca);
    TH1D* hDcaAllHft = new TH1D("hDca_all_hft", "hDca_all_hft", 100, 0, maxDca);


    TH1D* hPtRecoDCA[10];
    TH1D* hPtRecoHftDCA[10];
    TH1D* hDCAPtBins[10];
    TH1D* hDCAPtBinsHft[10];

    for (int j = 0; j < nDcaRatio; ++j) {
        histoname=Form("hPt_dca%.1f", dcaCut[j]);
        hPtRecoDCA[j]=new TH1D(histoname, histoname, 100,0,10);

        histoname=Form("hPt_hft_dca%.1f", dcaCut[j]);
        hPtRecoHftDCA[j]=new TH1D(histoname, histoname, 100,0,10);
    }

    for (int j = 0; j < nPtDca; ++j) {
        histoname=Form("hDca_pt_%.2f_%.2f", ptCut[j], ptCut[j+1]);
        hDCAPtBins[j]=new TH1D(histoname, histoname, 100, 0, maxDca);

        histoname=Form("hDca_hft_pt_%.2f_%.2f", ptCut[j], ptCut[j+1]);
        hDCAPtBinsHft[j]=new TH1D(histoname, histoname, 100, 0, maxDca);
    }

    const int nVars = sizeof(variables) / sizeof(TString);
    for (int k = 0; k < ntpData->GetEntries(); ++k) {
        ntpData->GetEntry(k);
        if (isTOF==0) continue;
        if (dca>1.5) continue;

        hnHitsFit->Fill(nHitsFit);
        hDcaAll->Fill(dca);

        for (int j = 0; j < nPtDca; ++j) {
            if (pt>ptCut[j] && pt<ptCut[j+1]) hDCAPtBins[j]->Fill(dca);
        }
        for (int i = 0; i < nDcaRatio; ++i) {
            if (dca<dcaCut[i]) hPtRecoDCA[i]->Fill(pt);
        }

        if (isHFT>0) {
            hnHitsFitHft->Fill(nHitsFit);
            hDcaAllHft->Fill(dca);
            for (int j = 0; j < nPtDca; ++j) {
                if (pt>ptCut[j] && pt<ptCut[j+1]) hDCAPtBinsHft[j]->Fill(dca);
            }
            for (int i = 0; i < nDcaRatio; ++i) {
                if (dca<dcaCut[i]) hPtRecoHftDCA[i]->Fill(pt);
            }
        }

//        for (int i = 0; i < nVars; ++i) {
//            cout<<"working on "<<variables[i]<<endl;
//            his = new TH1D("hist", variables[i], nBins[i], limMin[i], limMax[i]);
//            hisHft = new TH1D("hist_HFT", variables[i], nBins[i], limMin[i], limMax[i]);
//            ntpData -> Project("hist", variables[i]);
//            ntpData -> Project("hist_HFT", variables[i], "isHFT>0");
//
//            his->Scale(1/his->GetEntries());
//            hisHft->Scale(1/hisHft->GetEntries());
//
//            dataOut->cd();
//            hisname = Form("%s", variables[i].Data());
//            his->Write(hisname);
//            hisname = Form("%s_hft", variables[i].Data());
//            hisHft->Write(hisname);
//        }
    }
    TList *listOutQA = new TList();
    listOutQA->Add(hnHitsFit);
    listOutQA->Add(hnHitsFitHft);
    listOutQA->Add(hDcaAll);
    listOutQA->Add(hDcaAllHft);
    for (int l = 0; l < nPtDca; ++l) {
        listOutQA->Add(hDCAPtBins[l]);
        listOutQA->Add(hDCAPtBinsHft[l]);
    }

    for (int j = 0; j < nDcaRatio; ++j) {
        listOutQA->Add(hPtRecoDCA[j]);
        listOutQA->Add(hPtRecoHftDCA[j]);
    }

    dataOut->cd();
    listOutQA->Write("hists_QA", 1, 0);

    dataF->Close();
    dataOut->Close();
}

//________________________________________________________________________________________________________________________
void compareEvtEmb() {
//    TString fileNameSim = "/home/lukas/work/D0-fullEvent/analyse/test.out.root";
    TString fileNameSim = "/home/lukas/work/D0-fullEvent/analyse/out_test.7.1905.big.root";
    auto* dataSimF = new TFile(fileNameSim,"r");
    TList* list = (TList*)dataSimF -> Get("hists_event_QA;1");
    TH1F* vzSim=static_cast<TH1F*>(list->FindObject("hRcVzAccepVtx"));

    ///////
//    TH1D* vzdata = new TH1D("vzdata","vzdata",100,-6,6);
//    TH1D* vydata = new TH1D("vydata","vydata",100,-1,1);
//    TH1D* vxdata = new TH1D("vxdata","vxdata",100,-1,1);
//    TFile* dataF = new TFile("/home/lukas/work/dmesons/Dmaker_ndAu/Dmaker_dAu/analyse/PV/ntp.PV.1712.root","r");
//    TNtuple* ntpData = (TNtuple*)dataF -> Get("ntp_vertex;1");
////    ntpData->Scan();
//    cout<<ntpData->GetEntries()<<endl;
//    Float_t picoDstVz, nHftTracks, BBC, picoDstVx, picoDstVy;
//    ntpData->SetBranchAddress("picoDstVx", &picoDstVx);
//    ntpData->SetBranchAddress("picoDstVy", &picoDstVy);
//    ntpData->SetBranchAddress("picoDstVz", &picoDstVz);
//    ntpData->SetBranchAddress("BBC", &BBC);
//    ntpData->SetBranchAddress("nHftTracks", &nHftTracks);
//
//    for (int i = 0; i < ntpData->GetEntries(); ++i) {
//        ntpData->GetEntry(i);
////        if(BBC<900 && nHftTracks>2 && picoDstVy>-0.25 && picoDstVy<-0.16 && picoDstVx>-0.25 && picoDstVx<-0.16) vxdata->Fill(picoDstVx);
////        if(BBC<900 && nHftTracks>2 && picoDstVy>-0.25 && picoDstVy<-0.16 && picoDstVx>-0.25 && picoDstVx<-0.16) vydata->Fill(picoDstVy);
////        if(BBC<900 && nHftTracks>2 && picoDstVy>-0.25 && picoDstVy<-0.16 && picoDstVx>-0.25 && picoDstVx<-0.16) vzdata->Fill(picoDstVz);
//
////        if(BBC<900 && nHftTracks>2) vxdata->Fill(picoDstVx);
////        if(BBC<900 && nHftTracks>2) vydata->Fill(picoDstVy);
////        if(BBC<900 && nHftTracks>2) vzdata->Fill(picoDstVz);
//
//        vxdata->Fill(picoDstVx);
//        vydata->Fill(picoDstVy);
//        vzdata->Fill(picoDstVz);
//
//    }
//    TFile* dataout = new TFile("vzData.root","recreate");
//    vxdata->Write();
//    vydata->Write();
//    vzdata->Write();
/////////////////////

    TFile* dataF = new TFile("vzData.root","r");
    TH1F* vxdata = (TH1F*)dataF -> Get("vxdata");
    TH1F* vydata = (TH1F*)dataF -> Get("vydata");
    TH1F* vzdata = (TH1F*)dataF -> Get("vzdata");

    Double_t nBinsData = vzdata->Integral(vzdata->FindBin(-6), vzdata->FindBin(6));
    Double_t nBinsSim = vzSim->Integral(vzSim->FindBin(-6), vzSim->FindBin(6));

//    cout<<nBinsData<<" "<<nBinsSim<<endl;
//    cout<<nBinsData<<" "<<nBinsSim<<endl;
    cout<<vzdata->GetNbinsX()<<endl;
    cout<<vzSim->GetNbinsX()<<endl;
//    vzdata->Rebin(2);

    vzSim->Scale(1/vzSim->GetEntries());
    vxdata->Scale(1/vxdata->GetEntries());
    vydata->Scale(1/vydata->GetEntries());
    vzdata->Scale(1/vzdata->GetEntries());

    cout<<"vzsim"<<endl;
    vzSim->Fit("gaus");
    cout<<"vx"<<endl;
    vxdata->Fit("gaus");
    cout<<"vy"<<endl;
    vydata->Fit("gaus");
    cout<<"vzDat"<<endl;
    vzdata->Fit("gaus");
//    vzSim->Scale(1/nBinsSim);
//    vzdata->Scale(1/nBinsData);

//    vzdata->Divide(vzSim);
    TCanvas* c = new TCanvas("c","c", 800,900);
    vzSim->Draw();
    vzdata->Draw("same");

    TCanvas* cX = new TCanvas("cX","cX", 800,900);
    vxdata->Draw();

    TCanvas* cY = new TCanvas("cY","cY", 800,900);
    vydata->Draw();

//    dataF->Close();

}

//________________________________________________________________________________________________________________________
void compareEmb() {
    TString fileNameData = "/home/lukas/work/dmesons/Dmaker_ndAu/Dmaker_dAu/analyse/simInputsPrepare/tracksQA.data.root";
//    TString fileNameSim = "/home/lukas/work/D0-fullEvent/analyse/test.out.root";
    TString fileNameSim = "/home/lukas/work/D0-fullEvent/analyse/out_test.big.cc.setup6.root";
    TString fileNameSim1 = "/home/lukas/work/D0-fullEvent/analyse/out_test.big.MB.2705.root";

    auto* dataF = new TFile(fileNameData,"r");
    auto* dataSimF = new TFile(fileNameSim,"r");
    auto* dataSimF1 = new TFile(fileNameSim1,"r");


    TList *list[3];
    list[0] = (TList*)dataF -> Get("hists_QA;1");
    list[1] = (TList*)dataSimF -> Get("hists_QA;1");
    list[2] = (TList*)dataSimF1 -> Get("hists_QA;1");

    TString names[] = {"data", "sim", "sim1"};

    const int nFiles=3;

    Int_t color[]={46,9,8};

    TString histoname;

    TH1D* hnHitsFit[nFiles];
    TH1D* hnHitsFitHft[nFiles];
    TH1D* hDcaAll[nFiles];
    TH1D* hDcaAllHft[nFiles];
    TH1D* hPtRecoDCA[nFiles][10];
    TH1D* hPtRecoHftDCA[nFiles][10];
    TH1D* hDCAPtBins[nFiles][10];
    TH1D* hDCAPtBinsHft[nFiles][10];

    for (int k = 0; k < nFiles; ++k) {
        cout<<names[k]<<endl;
        hDcaAll[k]=new TH1D();
        hDcaAll[k]=static_cast<TH1D*>(list[k]->FindObject("hDca_all"));
        hDcaAll[k]->SetTitle("");
        draw(hDcaAll[k], "DCA [cm]", "1/N_{entries}", color[k]);

        hDcaAllHft[k]=new TH1D();
        hDcaAllHft[k]=static_cast<TH1D*>(list[k]->FindObject("hDca_all_hft"));
        hDcaAllHft[k]->SetTitle("HFT");
        draw(hDcaAllHft[k], "DCA [cm]", "1/N_{entries}", color[k]);

        hnHitsFit[k]=new TH1D();
        hnHitsFit[k]=static_cast<TH1D*>(list[k]->FindObject("hnHitsFit"));
        hnHitsFit[k]->SetTitle("");
        draw(hnHitsFit[k], "nHitsFit", "1/N_{entries}", color[k]);

        hnHitsFitHft[k]=new TH1D();
        hnHitsFitHft[k]=static_cast<TH1D*>(list[k]->FindObject("hnHitsFitHft"));
        hnHitsFitHft[k]->SetTitle("HFT");
        draw(hnHitsFitHft[k], "nHitsFit", "1/N_{entries}", color[k]);

        for (int j = 0; j < nDcaRatio; ++j) {
            histoname=Form("hPt_dca%.1f", dcaCut[j]);
            hPtRecoDCA[k][j]=new TH1D();
            hPtRecoDCA[k][j]=static_cast<TH1D*>(list[k]->FindObject(histoname));
            histoname=Form("DCA < %.1f cm", dcaCut[j]);
            hPtRecoDCA[k][j]->SetTitle(histoname);
            draw(hPtRecoDCA[k][j],"p_{T} [GeV/c]", "1/N_{entries}", color[k]);

            histoname=Form("hPt_hft_dca%.1f", dcaCut[j]);
            hPtRecoHftDCA[k][j]=new TH1D();
            hPtRecoHftDCA[k][j]=static_cast<TH1D*>(list[k]->FindObject(histoname));
            histoname=Form("HFT, DCA < %.1f cm", dcaCut[j]);
            hPtRecoHftDCA[k][j]->SetTitle(histoname);
            draw(hPtRecoHftDCA[k][j], "p_{T} [GeV/c]", "1/N_{entries}", color[k]);

        }

        for (int j = 0; j < nPtDca; ++j) {
            histoname=Form("hDca_pt_%.2f_%.2f", ptCut[j], ptCut[j+1]);
            hDCAPtBins[k][j]=new TH1D();
            hDCAPtBins[k][j]=static_cast<TH1D*>(list[k]->FindObject(histoname));
            histoname=Form("%.2f < p_{T} < %.2f GeV/c", ptCut[j], ptCut[j+1]);
            hDCAPtBins[k][j]->SetTitle(histoname);
            draw(hDCAPtBins[k][j], "DCA [cm]", "1/N_{entries}", color[k]);

            histoname=Form("hDca_hft_pt_%.2f_%.2f", ptCut[j], ptCut[j+1]);
            hDCAPtBinsHft[k][j]=new TH1D();
            hDCAPtBinsHft[k][j]=static_cast<TH1D*>(list[k]->FindObject(histoname));
            histoname=Form("HFT, %.2f < p_{T} < %.2f GeV/c", ptCut[j], ptCut[j+1]);
            hDCAPtBinsHft[k][j]->SetTitle(histoname);
            draw(hDCAPtBinsHft[k][j], "DCA [cm]", "1/N_{entries}", color[k]);

        }

    }

    TLegend *legend = new TLegend(0.767, 0.78, 0.92, 0.89);
    legend -> SetFillStyle(0);
    legend -> SetLineColor(0);
    legend -> SetTextSize(0.04);
    legend -> AddEntry(hDCAPtBins[0][0], names[0], "pl");
    legend -> AddEntry(hDCAPtBins[1][0], names[1], "pl");
    legend -> AddEntry(hDCAPtBins[2][0], names[2], "pl");

    const int nCanvas = nDcaRatio+nPtDca+2;
    cout<<nCanvas<<endl;

    TCanvas* c[nCanvas];
    for (int i = 0; i < nCanvas; ++i) {
        histoname=Form("c%i", i);
        c[i]=new TCanvas(histoname, histoname, 1500,800);
        c[i]->Divide(2,1);
    }

    int canId = -1;

    canId++;
    for (int k = 0; k < nFiles; ++k) {
        c[canId]->cd(1);
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.05);

        hnHitsFit[k]->Scale(1/hnHitsFit[k]->GetEntries());
        if (k==0) hnHitsFit[k]->Draw();
        hnHitsFit[k]->Draw("same");
    }

    for (int k = 0; k < nFiles; ++k) {
        c[canId]->cd(2);
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.05);

        hnHitsFitHft[k]->Scale(1/hnHitsFitHft[k]->GetEntries());
        if (k==0) hnHitsFitHft[k]->Draw();
        hnHitsFitHft[k]->Draw("same");
        legend->Draw("same");
    }
    //______________________________________________________
    canId++;
    for (int k = 0; k < nFiles; ++k) {
        c[canId]->cd(1);
        gPad->SetLogy();
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.05);

        hDcaAll[k]->Scale(1/hDcaAll[k]->GetEntries());
        if (k==0) hDcaAll[k]->Draw();
        hDcaAll[k]->Draw("same");
    }

    for (int k = 0; k < nFiles; ++k) {
        c[canId]->cd(2);
        gPad->SetLogy();
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.05);

        hDcaAllHft[k]->Scale(1/hDcaAllHft[k]->GetEntries());
        if (k==0) hDcaAllHft[k]->Draw();
        hDcaAllHft[k]->Draw("same");
        legend->Draw("same");
    }
    //______________________________________________________
    for (int j = 0; j < nDcaRatio; ++j) {
        canId++;
        for (int k = 0; k < nFiles; ++k) {
            c[canId]->cd(1);
            gPad->SetLogy();
            gPad->SetLeftMargin(0.15);
            gPad->SetRightMargin(0.05);

            hPtRecoDCA[k][j]->Scale(1/hPtRecoDCA[k][j]->GetEntries());
            if (k==0) hPtRecoDCA[k][j]->Draw();
            hPtRecoDCA[k][j]->Draw("same");
        }

        for (int k = 0; k < nFiles; ++k) {
            c[canId]->cd(2);
            gPad->SetLogy();
            gPad->SetLeftMargin(0.15);
            gPad->SetRightMargin(0.05);

            hPtRecoHftDCA[k][j]->Scale(1/hPtRecoHftDCA[k][j]->GetEntries());
            if (k==0) hPtRecoHftDCA[k][j]->Draw();
            hPtRecoHftDCA[k][j]->Draw("same");
            legend->Draw("same");

        }
    }

    for (int j = 0; j < nPtDca; ++j) {
        canId++;
        for (int k = 0; k < nFiles; ++k) {
            c[canId]->cd(1);
            gPad->SetLogy();
            gPad->SetLeftMargin(0.15);
            gPad->SetRightMargin(0.05);

            hDCAPtBins[k][j]->Scale(1/hDCAPtBins[k][j]->GetEntries());
            if (k==0) hDCAPtBins[k][j]->Draw();
            hDCAPtBins[k][j]->Draw("same");
        }

        for (int k = 0; k < nFiles; ++k) {
            c[canId]->cd(2);
            gPad->SetLogy();
            gPad->SetLeftMargin(0.15);
            gPad->SetRightMargin(0.05);

            hDCAPtBinsHft[k][j]->Scale(1/hDCAPtBinsHft[k][j]->GetEntries());
            if (k==0) hDCAPtBinsHft[k][j]->Draw();
            hDCAPtBinsHft[k][j]->Draw("same");
            legend->Draw("same");

        }

    }

    c[0]->SaveAs("trackQA.pdf(");
    for (int i = 1; i < nCanvas-2; ++i) {
        c[i]->SaveAs("trackQA.pdf");
    }
    c[nCanvas-1]->SaveAs("trackQA.pdf)");

    for (int i = 0; i < nCanvas; ++i) {
        c[i]->SaveAs(Form("img/trackQA_%i.png", i));
        c[i]->Close();
    }

}

//________________________________________________________________________________________________________________________
void compareVertex() {
//    TString fileNameSim = "/home/lukas/work/D0-fullEvent/analyse/test.out.root";
    TString fileNameSim = "/home/lukas/work/D0-fullEvent/analyse/out_test.7.1905.big.root";
    auto* dataSimF = new TFile(fileNameSim,"r");
    TList* list = (TList*)dataSimF -> Get("hists_event_QA;1");
    TH1F* vzSim=static_cast<TH1F*>(list->FindObject("hRcVzAccepVtx"));

    ///////
//    TH1D* vzdata = new TH1D("vzdata","vzdata",100,-6,6);
//    TH1D* vydata = new TH1D("vydata","vydata",100,-1,1);
//    TH1D* vxdata = new TH1D("vxdata","vxdata",100,-1,1);
//    TFile* dataF = new TFile("/home/lukas/work/dmesons/Dmaker_ndAu/Dmaker_dAu/analyse/PV/ntp.PV.1712.root","r");
//    TNtuple* ntpData = (TNtuple*)dataF -> Get("ntp_vertex;1");
////    ntpData->Scan();
//    cout<<ntpData->GetEntries()<<endl;
//    Float_t picoDstVz, nHftTracks, BBC, picoDstVx, picoDstVy;
//    ntpData->SetBranchAddress("picoDstVx", &picoDstVx);
//    ntpData->SetBranchAddress("picoDstVy", &picoDstVy);
//    ntpData->SetBranchAddress("picoDstVz", &picoDstVz);
//    ntpData->SetBranchAddress("BBC", &BBC);
//    ntpData->SetBranchAddress("nHftTracks", &nHftTracks);
//
//    for (int i = 0; i < ntpData->GetEntries(); ++i) {
//        ntpData->GetEntry(i);
////        if(BBC<900 && nHftTracks>2 && picoDstVy>-0.25 && picoDstVy<-0.16 && picoDstVx>-0.25 && picoDstVx<-0.16) vxdata->Fill(picoDstVx);
////        if(BBC<900 && nHftTracks>2 && picoDstVy>-0.25 && picoDstVy<-0.16 && picoDstVx>-0.25 && picoDstVx<-0.16) vydata->Fill(picoDstVy);
////        if(BBC<900 && nHftTracks>2 && picoDstVy>-0.25 && picoDstVy<-0.16 && picoDstVx>-0.25 && picoDstVx<-0.16) vzdata->Fill(picoDstVz);
//
////        if(BBC<900 && nHftTracks>2) vxdata->Fill(picoDstVx);
////        if(BBC<900 && nHftTracks>2) vydata->Fill(picoDstVy);
////        if(BBC<900 && nHftTracks>2) vzdata->Fill(picoDstVz);
//
//        vxdata->Fill(picoDstVx);
//        vydata->Fill(picoDstVy);
//        vzdata->Fill(picoDstVz);
//
//    }
//    TFile* dataout = new TFile("vzData.root","recreate");
//    vxdata->Write();
//    vydata->Write();
//    vzdata->Write();
/////////////////////

    TFile* dataF = new TFile("vzData.root","r");
    TH1F* vxdata = (TH1F*)dataF -> Get("vxdata");
    TH1F* vydata = (TH1F*)dataF -> Get("vydata");
    TH1F* vzdata = (TH1F*)dataF -> Get("vzdata");

    Double_t nBinsData = vzdata->Integral(vzdata->FindBin(-6), vzdata->FindBin(6));
    Double_t nBinsSim = vzSim->Integral(vzSim->FindBin(-6), vzSim->FindBin(6));

//    cout<<nBinsData<<" "<<nBinsSim<<endl;
//    cout<<nBinsData<<" "<<nBinsSim<<endl;
    cout<<vzdata->GetNbinsX()<<endl;
    cout<<vzSim->GetNbinsX()<<endl;
//    vzdata->Rebin(2);

    vzSim->Scale(1/vzSim->GetEntries());
    vxdata->Scale(1/vxdata->GetEntries());
    vydata->Scale(1/vydata->GetEntries());
    vzdata->Scale(1/vzdata->GetEntries());

    cout<<"vzsim"<<endl;
    vzSim->Fit("gaus");
    cout<<"vx"<<endl;
    vxdata->Fit("gaus");
    cout<<"vy"<<endl;
    vydata->Fit("gaus");
    cout<<"vzDat"<<endl;
    vzdata->Fit("gaus");
//    vzSim->Scale(1/nBinsSim);
//    vzdata->Scale(1/nBinsData);

//    vzdata->Divide(vzSim);
    TCanvas* c = new TCanvas("c","c", 800,900);
    vzSim->Draw();
    vzdata->Draw("same");

    TCanvas* cX = new TCanvas("cX","cX", 800,900);
    vxdata->Draw();

    TCanvas* cY = new TCanvas("cY","cY", 800,900);
    vydata->Draw();

//    dataF->Close();

}