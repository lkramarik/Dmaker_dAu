#include"TSystem.h"
#include"TAxis.h"
#include"TTree.h"
#include"TCanvas.h"
#include"TFile.h"
#include"TNtuple.h"
#include"TH1.h"
#include"TLegend.h"
#include"TPad.h"
#include"TRatioPlot.h"
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

//const int nPtDca=9;
//const float ptCut[nPtDca+1]={0.15, 0.3, 0.5, 0.7, 1.0, 1.2, 1.5, 2, 3, 5};
//_____________________________________________________________________________________________________________
void draw(TH1D* histo, TString xAxis, TString yAxis, Int_t color) {
    histo->GetXaxis()->SetTitle(xAxis);
    histo->GetYaxis()->SetTitle(yAxis);
    histo->SetLineColor(color);
    histo->SetMarkerColor(color);
    histo->SetMarkerStyle(20);
    histo->SetStats(0);

    Double_t dcaMax=0.2;
    if(xAxis=="p_{T} [GeV/c]") histo->GetXaxis()->SetRangeUser(0,4);
    if(xAxis=="DCA_{xy} [cm]") histo->GetXaxis()->SetRangeUser(-1*dcaMax,dcaMax);
    if(xAxis=="DCA_{z} [cm]") histo->GetXaxis()->SetRangeUser(-1*dcaMax,dcaMax);
    if(xAxis=="DCA [cm]") histo->GetXaxis()->SetRangeUser(0,dcaMax);

//    if(xAxis=="DCA_{xy} [cm]") histo->GetXaxis()->SetRangeUser(-0.1,0.1);
//    if(xAxis=="DCA_{z} [cm]") histo->GetXaxis()->SetRangeUser(-0.1,0.1);
//    if(xAxis=="DCA [cm]") histo->GetXaxis()->SetRangeUser(0,0.1);
}

//_____________________________________________________________________________________________________________
void createQATracks() {
    Int_t colors[] = {1,46,8,4,5,6,7,9,41,42,44,45,16};
//    TString fileNameData = "outputLocal.hists.root";
//    TString fileNameData = "ntp.track.2403.root";
//    TString fileNameData = "ntp.2403.sample.root";
//    TString fileNameData = "/home/lukas/work/dmesons/Dmaker_ndAu/Dmaker_dAu/analyse/simInputsPrepare/toHadd/ntp.ratio.1404.half.root";
//    TString fileNameData = "/home/lukas/work/dmesons/Dmaker_ndAu/Dmaker_dAu/analyse/simInputsPrepare/toHadd/ntp.ratio.1404.half.rootprimary.root";
    TString fileNameData = "/media/lukas/1E4183350724EC9C/tracks.sim.vzvpd3.hotspot.root";
    auto* dataF = new TFile(fileNameData,"r");
    auto* dataOut = new TFile("tracksQA.data.new.root","recreate");
    auto* ntpData = (TNtuple*)dataF -> Get("ntp_tracks");

    double maxPt = 6;

    TString variables[] = {"dca", "nHitsFit"};
    Double_t limMin[] = {0, 0.5};
    Double_t limMax[] = {1.5, 50.5};
    Int_t nBins[] = {100, 51};
    TH1D *his;
    TH1D *hisHft;

    Float_t  pt, dca, isHFT, isTOF, zdc, nHitsFit, refMult, nHftTracks, nTofTracks, particleId, isPrimary, dcaxy, dcaz;
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

    Double_t maxDca=0.2;

    //init histos
    TString histoname;
    TH1D* hnHitsFit = new TH1D("hnHitsFit", "hnHitsFit", 51, 0.5, 50.5);
    TH1D* hnHitsFitHft = new TH1D("hnHitsFitHft", "hnHitsFitHft", 51, 0.5, 50.5);
    TH1D* hDcaAll = new TH1D("hDca_all", "hDca_all", 100, 0, maxDca);
    TH1D* hDcaAllHft = new TH1D("hDca_all_hft", "hDca_all_hft", 100, 0, maxDca);
    TH1D* hDcaxyAll = new TH1D("hDcaxy_all", "hDcaxy_all", 100, -1*maxDca, maxDca);
    TH1D* hDcaxyAllHft = new TH1D("hDcaxy_all_hft", "hDcaxy_all_hft", 100, -1*maxDca, maxDca);
    TH1D* hDcazAll = new TH1D("hDcaz_all", "hDcaz_all", 100, -1*maxDca, maxDca);
    TH1D* hDcazAllHft = new TH1D("hDcaz_all_hft", "hDcaz_all_hft", 100, -1*maxDca, maxDca);

    TH1D* hPtRecoDCA[10];
    TH1D* hPtRecoHftDCA[10];
    TH1D* hDCAPtBins[10];
    TH1D* hDCAPtBinsHft[10];
    TH1D* hDCAxyPtBins[10];
    TH1D* hDCAxyPtBinsHft[10];
    TH1D* hDCAzPtBins[10];
    TH1D* hDCAzPtBinsHft[10];

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

        histoname=Form("hDcaxy_pt_%.2f_%.2f", ptCut[j], ptCut[j+1]);
        hDCAxyPtBins[j]=new TH1D(histoname, histoname, 100, -1*maxDca, maxDca);
        histoname=Form("hDcaxy_hft_pt_%.2f_%.2f", ptCut[j], ptCut[j+1]);
        hDCAxyPtBinsHft[j]=new TH1D(histoname, histoname, 100, -1*maxDca, maxDca);

        histoname=Form("hDcaz_pt_%.2f_%.2f", ptCut[j], ptCut[j+1]);
        hDCAzPtBins[j]=new TH1D(histoname, histoname, 100, -1*maxDca, maxDca);
        histoname=Form("hDcaz_hft_pt_%.2f_%.2f", ptCut[j], ptCut[j+1]);
        hDCAzPtBinsHft[j]=new TH1D(histoname, histoname, 100, -1*maxDca, maxDca);
    }

    const int nVars = sizeof(variables) / sizeof(TString);
    for (int k = 0; k < ntpData->GetEntries()/3; ++k) {
        ntpData->GetEntry(k);
        if (isTOF==0) continue;
        if (isPrimary==0) continue;
        if (particleId!=0) continue;
//        if (nHitsFit<20) continue;
//        if (nTofTracks<1) continue;
        if (dca>1.5) continue;

        hnHitsFit->Fill(nHitsFit);
        hDcaAll->Fill(dca);
        hDcaxyAll->Fill(dcaxy);
        hDcazAll->Fill(dcaz);

        for (int j = 0; j < nPtDca; ++j) {
            if (pt>ptCut[j] && pt<ptCut[j+1]) {
                hDCAPtBins[j]->Fill(dca);
                hDCAxyPtBins[j]->Fill(dcaxy);
                hDCAzPtBins[j]->Fill(dcaz);
            }
        }
        for (int i = 0; i < nDcaRatio; ++i) {
            if (dca<dcaCut[i]) hPtRecoDCA[i]->Fill(pt);
        }

        if (isHFT>0) {
            hnHitsFitHft->Fill(nHitsFit);
            hDcaAllHft->Fill(dca);
            hDcaxyAllHft->Fill(dcaxy);
            hDcazAllHft->Fill(dcaz);
            for (int j = 0; j < nPtDca; ++j) {
                if (pt>ptCut[j] && pt<ptCut[j+1]) {
                    hDCAPtBinsHft[j]->Fill(dca);
                    hDCAxyPtBinsHft[j]->Fill(dcaxy);
                    hDCAzPtBinsHft[j]->Fill(dcaz);
                }
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
    listOutQA->Add(hDcaxyAll);
    listOutQA->Add(hDcaxyAllHft);
    listOutQA->Add(hDcazAll);
    listOutQA->Add(hDcazAllHft);
    for (int l = 0; l < nPtDca; ++l) {
        listOutQA->Add(hDCAPtBins[l]);
        listOutQA->Add(hDCAPtBinsHft[l]);
        listOutQA->Add(hDCAxyPtBins[l]);
        listOutQA->Add(hDCAxyPtBinsHft[l]);
        listOutQA->Add(hDCAzPtBins[l]);
        listOutQA->Add(hDCAzPtBinsHft[l]);
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
    std::vector<TList*> list;
    std::vector<TString> names;


//    auto* dataF41 = new TFile("/home/lukas/work/D0-fullEvent/analyse/out_production.vtx.3M.1308.root","r");
//    list.push_back((TList*)dataF41 -> Get("hists_QA;1"));
//    names.push_back("sim 1308");

//    auto* dataF41a = new TFile("/home/lukas/work/dmesons/Dmaker_ndAu/Dmaker_dAu/analyse/simInputsPrepare/tracksQA.data.new.root","r");
//    list.push_back((TList*)dataF41a -> Get("hists_QA;1"));
//    names.push_back("dat");
//
//    auto* dataF41b = new TFile("/home/lukas/work/D0-fullEvent/analyse/out_production.vtx.3M.geant.primaryTag.1307.root","r");
//    list.push_back((TList*)dataF41b -> Get("hists_QA;1"));
//    names.push_back("sim");
////
//
//    auto* dataF41 = new TFile("/home/lukas/work/D0-fullEvent/analyse/out_production.vtx.3M.geant.primaryTag.1307.pions.goodEv.hotspot.root","r");
//    list.push_back((TList*)dataF41 -> Get("hists_QA;1"));
//    names.push_back("sim goode");
//
    auto* dataF4 = new TFile("/home/lukas/work/D0-fullEvent/analyse/out_production.vtx.3M.geant.primaryTag.1307.nonP.root","r");
    list.push_back((TList*)dataF4 -> Get("hists_QA;1"));
    names.push_back("sim nonP");
//
//
    auto* dataF2 = new TFile("/home/lukas/work/dmesons/Dmaker_ndAu/Dmaker_dAu/analyse/simInputsPrepare/tracksQA.data.new.primary.root","r");
    list.push_back((TList*)dataF2 -> Get("hists_QA;1"));
    names.push_back("primary");

    auto* dataFad = new TFile("/home/lukas/work/dmesons/Dmaker_ndAu/Dmaker_dAu/analyse/simInputsPrepare/tracksQA.data.new.TOFmatched.root","r");
    list.push_back((TList*)dataFad -> Get("hists_QA;1"));
    names.push_back("TOF match");


//
//    auto* dataF1 = new TFile("/home/lukas/work/dmesons/Dmaker_ndAu/Dmaker_dAu/analyse/simInputsPrepare/tracksQA.data.new.pions.primary.tofmatched.root","r");
//    list.push_back((TList*)dataF1 -> Get("hists_QA;1"));
//    names.push_back("pions p TOF");
//
    auto* dataF3 = new TFile("/home/lukas/work/dmesons/Dmaker_ndAu/Dmaker_dAu/analyse/simInputsPrepare/tracksQA.data.new.all.root","r");
    list.push_back((TList*)dataF3 -> Get("hists_QA;1"));
    names.push_back("all");//

//    auto* dataF3q = new TFile("/home/lukas/work/dmesons/Dmaker_ndAu/Dmaker_dAu/analyse/simInputsPrepare/tracksQA.data.nTof0.root","r");
//    list.push_back((TList*)dataF3q -> Get("hists_QA;1"));
//    names.push_back("all ntof0");




//    auto* dataF2 = new TFile("/home/lukas/work/dmesons/Dmaker_ndAu/Dmaker_dAu/analyse/simInputsPrepare/tracksQA.data.new.pions.primary.root","r");
//    list.push_back((TList*)dataF2 -> Get("hists_QA;1"));
//    names.push_back("pions p");


//
//    auto* dataF30 = new TFile("/home/lukas/work/dmesons/Dmaker_ndAu/Dmaker_dAu/analyse/simInputsPrepare/tracksQA.data.new.nhitsfit20.allpions.root","r");
//    list.push_back((TList*)dataF30 -> Get("hists_QA;1"));
//    names.push_back("all pions nhitsfit20");
//

//
//    auto* dataF = new TFile("/home/lukas/work/dmesons/Dmaker_ndAu/Dmaker_dAu/analyse/simInputsPrepare/tracksQA.data.new.all.root","r");
//    list.push_back((TList*)dataF -> Get("hists_QA;1"));
//    names.push_back("all");
//
//
//
//
//

//
//
//
////
////    auto* dataF5 = new TFile("/home/lukas/work/D0-fullEvent/analyse/out_production.vtx.3M.geant.primaryTag.1307.pions.all.root","r");
////    list.push_back((TList*)dataF5 -> Get("hists_QA;1"));
////    names.push_back("sim all pions");
////
//    auto* dataF6 = new TFile("/home/lukas/work/D0-fullEvent/analyse/out_production.vtx.3M.geant.primaryTag.1307.pions.primary.nhits20.root","r");
//    list.push_back((TList*)dataF6 -> Get("hists_QA;1"));
//    names.push_back("sim P n20 pions");

    const int nFiles=list.size();

    Int_t color[]={46,9,8,7,1,42,40,5};
    Double_t scale[]={ 10000000, 957., 744.};

    TLegend *legend = new TLegend(0.576, 0.68, 0.73, 0.89);
    legend -> SetFillStyle(0);
    legend -> SetLineColor(0);
    legend -> SetTextSize(0.04);

    TString histoKinds[]={"hDca_all", "hDcaxy_all", "hDcaz_all", "hDca_all_hft", "hDcaxy_all_hft", "hDcaz_all_hft", "hnHitsFit", "hnHitsFitHft"};
//    TString xAxisKinds[]

    TString histKinds[]={"hDca_pt", "hDca_hft_pt", "hDcaxy_pt", "hDcaxy_hft_pt", "hDcaz_pt", "hDcaz_hft_pt"};
    TString titleKinds[]={"", "HFT tracks ","", "HFT tracks ","", "HFT tracks "};
    TString xAxisKinds[]={"DCA [cm]", "DCA [cm]","DCA_{xy} [cm]", "DCA_{xy} [cm]","DCA_{z} [cm]", "DCA_{z} [cm]"};
    const int nHistKinds = sizeof(histKinds) / sizeof(TString);

    TString nameHisto, titleHisto;
    TH1D* hDCApTs[nFiles][nHistKinds][nPtDca];
    TRatioPlot* ratios[nFiles][nHistKinds][nPtDca];
    TCanvas* canPtBins[2][nHistKinds/2][nPtDca];
    TCanvas* canRatiosBins[nFiles][nHistKinds][nPtDca];
    TGraph *g[nFiles][nHistKinds][nPtDca];
    for (int k = 0; k < nFiles; ++k) {
        for (int i = 0; i < nHistKinds; ++i) {
            for (int j = 0; j < nPtDca; ++j) {
                if(k==0 && i<nHistKinds/2){
                    for (int canI = 0; canI < 2; ++canI) {
                        nameHisto = Form("c%i%i%i", i, j, canI);
                        canPtBins[canI][i][j] = new TCanvas(nameHisto, nameHisto, 1500, 800);
                        canPtBins[canI][i][j]->Divide(2, 1);
                        canPtBins[canI][i][j]->cd(1);
                        gPad->SetLogy();
                        gPad->SetLeftMargin(0.15);
                        gPad->SetRightMargin(0.05);
                        canPtBins[canI][i][j]->cd(2);
                        gPad->SetLogy();
                        gPad->SetLeftMargin(0.15);
                        gPad->SetRightMargin(0.05);
                    }
                }
                nameHisto = Form("cRatio%i%i%i", i, j, k);
                canRatiosBins[k][i][j]=new TCanvas(nameHisto, nameHisto, 1500, 800);
                nameHisto=Form("%s_%.2f_%.2f", histKinds[i].Data(), ptCut[j], ptCut[j+1]);
                titleHisto=Form("%s %.2f < p_{T} < %.2f GeV/c", titleKinds[i].Data(), ptCut[j], ptCut[j+1]);
                hDCApTs[k][i][j]=static_cast<TH1D*>(list[k]->FindObject(nameHisto));
                hDCApTs[k][i][j]->SetTitle(titleHisto);
                hDCApTs[k][i][j]->Scale(10/hDCApTs[k][i][j]->Integral());
                draw(hDCApTs[k][i][j], xAxisKinds[i], "Normalized Yield (a.u.)", color[k]);
                ratios[k][i][j]=new TRatioPlot(hDCApTs[k][i][j],hDCApTs[0][i][j]);
                canRatiosBins[k][i][j]->cd();
                ratios[k][i][j]->Draw();
                ratios[k][i][j]->GetLowerRefGraph() -> SetMinimum(0.001); //if this is 0.0, problems with X axis
                ratios[k][i][j]->GetLowerRefGraph() -> SetMaximum(2);
            }
        }
    }

    for (int k = 0; k < nFiles; ++k) {
        legend -> AddEntry(hDCApTs[k][0][0], names[k], "pl");
    }

    for (int j = 0; j < nPtDca; ++j) {
        for (int i = 0; i < nHistKinds; ++i) {
            for (int k = 0; k < nFiles; ++k) {
                if (i%2==1) canPtBins[0][(i-1)/2][j]->cd(2);
                if (i%2==0) canPtBins[0][i/2][j]->cd(1);
                if (k==0) hDCApTs[k][i][j]->Draw();
                else hDCApTs[k][i][j]->Draw("same");

                if (i%2==1) canPtBins[1][(i-1)/2][j]->cd(2);
                if (i%2==0) canPtBins[1][i/2][j]->cd(1);
//                if (k==1) ratios[k][i][j]->Draw();
                if (k==0) ratios[k][i][j]->Draw(); else {
                    g[k][i][j] = ratios[k][i][j]->GetLowerRefGraph();
                    TH1F *hUp = (TH1F *) ratios[k][i][j]->GetUpperRefObject();
                    ratios[0][i][j]->GetLowerPad()->cd();
                    g[k][i][j]->SetMarkerColor(hUp->GetMarkerColor());
                    g[k][i][j]->SetLineColor(hUp->GetLineColor());
                    g[k][i][j]->Draw("same,p");
                    ratios[0][i][j]->GetUpperPad()->cd();
                    hUp->DrawCopy("same");
                }
            }
            canPtBins[0][(i)/2][j]->cd(2);
            legend->Draw("same");
            canPtBins[1][(i)/2][j]->cd(2);
            legend->Draw("same");
        }
    }


    for (int k = 0; k < nFiles; ++k) {
        for (int i = 0; i < nHistKinds; ++i) {
            for (int j = 0; j < nPtDca; ++j) {
                canRatiosBins[k][i][j]->Close();
                if (k == 0 && i < nHistKinds / 2) {
                    for (int canI = 0; canI < 2; ++canI) {
                        nameHisto = Form("img/new/%i_%s_%i.png", canI, histoKinds[i].Data(), j);
                        canPtBins[canI][i][j]->SaveAs(nameHisto);
                        canPtBins[canI][i][j]->Close();
                    }
                }
            }
        }
    }





    /*
    for (int k = 0; k < nFiles; ++k) {
        cout<<names[k]<<endl;
        hDcaAll[k]=new TH1D();
        hDcaAll[k]=static_cast<TH1D*>(list[k]->FindObject("hDca_all"));
        hDcaAll[k]->SetTitle("");
        draw(hDcaAll[k], "DCA [cm]", "1/N_{entries}", color[k]);

        hDcaxyAll[k]=new TH1D();
        hDcaxyAll[k]=static_cast<TH1D*>(list[k]->FindObject("hDcaxy_all"));
        hDcaxyAll[k]->SetTitle("");
        draw(hDcaxyAll[k], "DCA_{xy} [cm]", "1/N_{entries}", color[k]);

        hDcazAll[k]=new TH1D();
        hDcazAll[k]=static_cast<TH1D*>(list[k]->FindObject("hDcaz_all"));
        hDcazAll[k]->SetTitle("");
        draw(hDcazAll[k], "DCA_{z} [cm]", "1/N_{entries}", color[k]);

        hDcaAllHft[k]=new TH1D();
        hDcaAllHft[k]=static_cast<TH1D*>(list[k]->FindObject("hDca_all_hft"));
        hDcaAllHft[k]->SetTitle("HFT tracks");
        draw(hDcaAllHft[k], "DCA [cm]", "1/N_{entries}", color[k]);

        hDcaxyAllHft[k]=new TH1D();
        hDcaxyAllHft[k]=static_cast<TH1D*>(list[k]->FindObject("hDcaxy_all_hft"));
        hDcaxyAllHft[k]->SetTitle("HFT tracks");
        draw(hDcaxyAllHft[k], "DCA_{xy} [cm]", "1/N_{entries}", color[k]);

        hDcazAllHft[k]=new TH1D();
        hDcazAllHft[k]=static_cast<TH1D*>(list[k]->FindObject("hDcaz_all_hft"));
        hDcazAllHft[k]->SetTitle("HFT tracks");
        draw(hDcazAllHft[k], "DCA_{z} [cm]", "1/N_{entries}", color[k]);

        hnHitsFit[k]=new TH1D();
        hnHitsFit[k]=static_cast<TH1D*>(list[k]->FindObject("hnHitsFit"));
        hnHitsFit[k]->SetTitle("");
        draw(hnHitsFit[k], "nHitsFit", "1/N_{entries}", color[k]);

        hnHitsFitHft[k]=new TH1D();
        hnHitsFitHft[k]=static_cast<TH1D*>(list[k]->FindObject("hnHitsFitHft"));
        hnHitsFitHft[k]->SetTitle("HFT tracks");
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
            histoname=Form("HFT tracks, %.2f < p_{T} < %.2f GeV/c", ptCut[j], ptCut[j+1]);
            hDCAPtBinsHft[k][j]->SetTitle(histoname);
            draw(hDCAPtBinsHft[k][j], "DCA [cm]", "1/N_{entries}", color[k]);

            //xy
            histoname=Form("hDcaxy_pt_%.2f_%.2f", ptCut[j], ptCut[j+1]);
            hDCAxyPtBins[k][j]=new TH1D();
            hDCAxyPtBins[k][j]=static_cast<TH1D*>(list[k]->FindObject(histoname));
            histoname=Form("%.2f < p_{T} < %.2f GeV/c", ptCut[j], ptCut[j+1]);
            hDCAxyPtBins[k][j]->SetTitle(histoname);
            draw(hDCAxyPtBins[k][j], "DCA_{xy} [cm]", "1/N_{entries}", color[k]);

            histoname=Form("hDcaxy_hft_pt_%.2f_%.2f", ptCut[j], ptCut[j+1]);
            hDCAxyPtBinsHft[k][j]=new TH1D();
            hDCAxyPtBinsHft[k][j]=static_cast<TH1D*>(list[k]->FindObject(histoname));
            histoname=Form("HFT tracks, %.2f < p_{T} < %.2f GeV/c", ptCut[j], ptCut[j+1]);
            hDCAxyPtBinsHft[k][j]->SetTitle(histoname);
            draw(hDCAxyPtBinsHft[k][j], "DCA_{xy} [cm]", "1/N_{entries}", color[k]);

            //z
            histoname=Form("hDcaz_pt_%.2f_%.2f", ptCut[j], ptCut[j+1]);
            hDCAzPtBins[k][j]=new TH1D();
            hDCAzPtBins[k][j]=static_cast<TH1D*>(list[k]->FindObject(histoname));
            histoname=Form("%.2f < p_{T} < %.2f GeV/c", ptCut[j], ptCut[j+1]);
            hDCAzPtBins[k][j]->SetTitle(histoname);
            draw(hDCAzPtBins[k][j], "DCA_{z} [cm]", "1/N_{entries}", color[k]);

            histoname=Form("hDcaz_hft_pt_%.2f_%.2f", ptCut[j], ptCut[j+1]);
            hDCAzPtBinsHft[k][j]=new TH1D();
            hDCAzPtBinsHft[k][j]=static_cast<TH1D*>(list[k]->FindObject(histoname));
            histoname=Form("HFT tracks, %.2f < p_{T} < %.2f GeV/c", ptCut[j], ptCut[j+1]);
            hDCAzPtBinsHft[k][j]->SetTitle(histoname);
            draw(hDCAzPtBinsHft[k][j], "DCA_{z} [cm]", "1/N_{entries}", color[k]);
        }
        legend -> AddEntry(hDCAPtBins[k][0], names[k], "pl");

    }

    cout<<"loaded"<<endl;

    const int nCanvas = nDcaRatio*3+nPtDca+4;
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

    //XY______________________________________________________
    canId++;
    for (int k = 0; k < nFiles; ++k) {
        c[canId]->cd(1);
        gPad->SetLogy();
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.05);

        hDcaxyAll[k]->Scale(1/hDcaxyAll[k]->GetEntries());
        if (k==0) hDcaxyAll[k]->Draw();
        hDcaxyAll[k]->Draw("same");
    }

    for (int k = 0; k < nFiles; ++k) {
        c[canId]->cd(2);
        gPad->SetLogy();
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.05);

        hDcaxyAllHft[k]->Scale(1/hDcaxyAllHft[k]->GetEntries());
        if (k==0) hDcaxyAllHft[k]->Draw();
        hDcaxyAllHft[k]->Draw("same");
        legend->Draw("same");
    }

    //ZZZ______________________________________________________
    canId++;
    for (int k = 0; k < nFiles; ++k) {
        c[canId]->cd(1);
        gPad->SetLogy();
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.05);

        hDcazAll[k]->Scale(1/hDcazAll[k]->GetEntries());
        if (k==0) hDcazAll[k]->Draw();
        hDcazAll[k]->Draw("same");
    }

    for (int k = 0; k < nFiles; ++k) {
        c[canId]->cd(2);
        gPad->SetLogy();
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.05);

        hDcazAllHft[k]->Scale(1/hDcazAllHft[k]->GetEntries());
        if (k==0) hDcazAllHft[k]->Draw();
        hDcazAllHft[k]->Draw("same");
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
    //______________________________________________________
    for (int j = 0; j < nPtDca; ++j) {
        canId++;
        for (int k = 0; k < nFiles; ++k) {
            c[canId]->cd(1);
            gPad->SetLogy();
            gPad->SetLeftMargin(0.15);
            gPad->SetRightMargin(0.05);
            hDCAPtBins[k][j]->Scale(1./hDCAPtBins[k][j]->Integral());
            if (k==0) hDCAPtBins[k][j]->Draw();
            hDCAPtBins[k][j]->Draw("same");
        }

        for (int k = 0; k < nFiles; ++k) {
            c[canId]->cd(2);
            gPad->SetLogy();
            gPad->SetLeftMargin(0.15);
            gPad->SetRightMargin(0.05);
            hDCAPtBinsHft[k][j]->Scale(1./hDCAPtBinsHft[k][j]->GetEntries());
            if (k==0) hDCAPtBinsHft[k][j]->Draw();
            hDCAPtBinsHft[k][j]->Draw("same");
            legend->Draw("same");
        }
    }
    //______________________________________________________    //______________________________________________________
    for (int j = 0; j < nPtDca; ++j) {
        canId++;
        for (int k = 0; k < nFiles; ++k) {
            c[canId]->cd(1);
            gPad->SetLogy();
            gPad->SetLeftMargin(0.15);
            gPad->SetRightMargin(0.05);
            hDCAxyPtBins[k][j]->Scale(1./hDCAxyPtBins[k][j]->Integral());
            if (k==0) hDCAxyPtBins[k][j]->Draw();
            hDCAxyPtBins[k][j]->Draw("same");
        }

        for (int k = 0; k < nFiles; ++k) {
            c[canId]->cd(2);
            gPad->SetLogy();
            gPad->SetLeftMargin(0.15);
            gPad->SetRightMargin(0.05);
            hDCAxyPtBinsHft[k][j]->Scale(1./hDCAxyPtBinsHft[k][j]->GetEntries());
            if (k==0) hDCAxyPtBinsHft[k][j]->Draw();
            hDCAxyPtBinsHft[k][j]->Draw("same");
            legend->Draw("same");
        }
    }
    //______________________________________________________    //______________________________________________________
    for (int j = 0; j < nPtDca; ++j) {
        canId++;
        for (int k = 0; k < nFiles; ++k) {
            c[canId]->cd(1);
            gPad->SetLogy();
            gPad->SetLeftMargin(0.15);
            gPad->SetRightMargin(0.05);
            hDCAzPtBins[k][j]->Scale(1./hDCAzPtBins[k][j]->Integral());
            if (k==0) hDCAzPtBins[k][j]->Draw();
            hDCAzPtBins[k][j]->Draw("same");
        }

        for (int k = 0; k < nFiles; ++k) {
            c[canId]->cd(2);
            gPad->SetLogy();
            gPad->SetLeftMargin(0.15);
            gPad->SetRightMargin(0.05);
            hDCAzPtBinsHft[k][j]->Scale(1./hDCAzPtBinsHft[k][j]->GetEntries());
            if (k==0) hDCAzPtBinsHft[k][j]->Draw();
            hDCAzPtBinsHft[k][j]->Draw("same");
            legend->Draw("same");
        }
    }
    //______________________________________________________




    c[0]->SaveAs("trackQA.pdf(");
    for (int i = 1; i < nCanvas-2; ++i) {
        c[i]->SaveAs("trackQA.pdf");
    }
    c[nCanvas-1]->SaveAs("trackQA.pdf)");

    for (int i = 0; i < nCanvas; ++i) {
        c[i]->SaveAs(Form("img/trackQA_%i.png", i));
        c[i]->Close();
    }
*/
}

//________________________________________________________________________________________________________________________



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