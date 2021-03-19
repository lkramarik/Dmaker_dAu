#include<iostream>
#include<fstream>
#include<vector>
#include"TH1.h"
#include"TH2.h"
#include"TF1.h"
#include"TMath.h"
#include"TCanvas.h"
#include"TFile.h"
#include"TLatex.h"
#include"TStyle.h"
#include"TPad.h"
#include"TLegend.h"
#include"TPaveText.h"
#include"TAxis.h"
#include"TTree.h"
#include"TRatioPlot.h"
#include"TFitResultPtr.h"
#include"TFitResult.h"
#include"TString.h"

using namespace std;
//const int nmultEdge = 4;
//float const multEdge[nmultEdge+1] = {0, 8, 12, 20, 200};
//
const int nmultEdge = 1;
float const multEdge[nmultEdge+1] = {0, 200};

Float_t ptBins[]={0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.65, 0.75, 0.85, 0.95, 1.2, 1.45, 1.7, 2.0, 2.3, 2.6, 2.9, 3.2, 3.5, 4., 4.5, 5};

//Float_t maxDCA=3.0;
Float_t maxDCA=1.5;
Int_t nBinsDCA=50;

//void effCalc(TString, TString);
int getPtBin(float pt, int nPtBins);
void initHistosPtBins(TH1F** histoArray, int nHistos, TString histoName, Int_t nBins, Double_t axMin, Double_t axMax, TString particle);
TString makeHistosMC(TString inputFile, TString particle);
TString makeHistosData(TString inputFile,TString particle);
TH1F* makeRatioPlot(TH1F** histoArray, TString histoName);
void setHistoStyle(TH1F* histo, Int_t color, Int_t marker, TString titleX, TString titleY);
//void makeSpline(TGraph* gr);
TGraph* makeSpline(TGraph* gr);

//________________________________________________________________________________________
void effSys(){
    bool plotting= true;

    TString particles[] = {"kaon", "pion"};
    TFile* fOutputs=new TFile("tpc_sys.root", "RECREATE");

    TString fileNames[2];

    int iPart=0;

    for (int iPart = 0; iPart < 2; ++iPart) {

        fileNames[0]=makeHistosMC("/home/lukas/work/embedding_dAu/analyse/ntp.2303.all.root", particles[iPart]);
        fileNames[1]=makeHistosData("ntp.trk.0903.smaller.root", particles[iPart]);

        int  color[] = {1, kMagenta+2};


        TString legendNames[]={"embedding", "data"};
        TString legendTitle[]={"w/o HFT hits: DCA<1. / DCA<1.5", "w/o HFT hits: nHitsFit>20 / nHitsFit>15",
                               "HFT required: DCA<1. / DCA<1.5", "HFT required: nHitsFit>20 / nHitsFit>15",
                               "DCA<1. / DCA<1.5", "nHitsFit>20 / nHitsFit>15"};
        TString varNames[]= {"dca", "nHitsFit",
                             "dca_hft", "nHitsFit_hft",
                             "dca_all", "nHitsFit_all"};
        TString varNamesAxis[]= {"DCA [cm]", "nHitsFit",
                                 "DCA [cm]", "nHitsFit",
                                 "DCA [cm]", "nHitsFit"};
        const int nVars = sizeof(varNames) / sizeof(TString);

        int nPoints = sizeof(ptBins) / sizeof(Float_t);

        TH1F* hRatio[nVars][2]; //[dca/nHitsfit] [data/MC]
        TH1F *h[nVars][2][nPoints-1]; //[dca/nHitsfit] [data/MC] [pTbins]

        TLegend *legend[nVars];
        TLegend *legendOne;
        legendOne=new TLegend(0.596,0.114, 0.76, 0.312,"","brNDC");
        legendOne->SetFillStyle(0);
        legendOne->SetLineColor(0);
        legendOne->SetTextSize(0.035);

        TCanvas *canVar[nVars][nPoints-1];

        TString tmpName;
        for (int iVar = 0; iVar < nVars; ++iVar) {
            legend[iVar] = new TLegend(0.187,0.112, 0.35, 0.30,"","brNDC");
            legend[iVar] -> SetFillStyle(0);
            legend[iVar] -> SetLineColor(0);
            legend[iVar] -> SetTextSize(0.035);

            for (int iFile = 0; iFile < 2; ++iFile) { //histo laoading
                TFile* fIn=new TFile(fileNames[iFile],"READ");
                tmpName=varNames[iVar]+"_ratio";
                cout<<tmpName<<endl;
                hRatio[iVar][iFile] = (TH1F*)fIn->Get(tmpName);

                setHistoStyle(hRatio[iVar][iFile], color[iFile], 20, "p_{T} [GeV/c]", "Ratio");
                hRatio[iVar][iFile]->GetXaxis()->SetRangeUser(0.2, 5.);
                hRatio[iVar][iFile]->GetYaxis()->SetRangeUser(0.5, 1.1);
                hRatio[iVar][iFile]->SetTitle("");

                legend[iVar]->AddEntry(hRatio[iVar][iFile], legendNames[iFile], "pl");
                if (iVar==0) legendOne->AddEntry(hRatio[iVar][iFile], legendNames[iFile], "pl");

                for (int iPt = 0; iPt < nPoints-1; ++iPt) {
                    tmpName=Form("%s_%.2f_%.2f", varNames[iVar].Data(), ptBins[iPt], ptBins[iPt+1]);
                    h[iVar][iFile][iPt] = (TH1F*)fIn->Get(tmpName);
                    setHistoStyle(h[iVar][iFile][iPt], color[iFile], 20, varNamesAxis[iVar], "Counts");
                }
            }

            if (!plotting) continue;
            for (int iFile = 0; iFile < 2; ++iFile) {
                for (int iPt = 0; iPt < nPoints-1; ++iPt) {
                    tmpName=Form("%s_%.2f_%.2f", varNames[iVar].Data(), ptBins[iPt], ptBins[iPt+1]);
                    if (iFile==0) canVar[iVar][iPt] = new TCanvas(tmpName,tmpName,900,1000);
                    canVar[iVar][iPt]->cd();
                    gPad->SetLeftMargin(0.15);
                    gPad->SetRightMargin(0.05);
                    if (iVar==0 || iVar==2 || iVar==4) gPad->SetLogy();

                    h[iVar][iFile][iPt]->Scale(1. / h[iVar][iFile][iPt]->Integral(1., h[iVar][iFile][iPt]->GetNbinsX()));

                    if (iFile==0) h[iVar][iFile][iPt]->Draw();
                    else  h[iVar][iFile][iPt]->Draw("same");
                    legendOne->Draw("same");
                }
            }

        }

        if (plotting) {
            for (int iVar = 0; iVar < nVars; ++iVar) {
                for (int iPt = 0; iPt < nPoints - 1; ++iPt) {
                    tmpName = "img/" + particles[iPart] + "/" + (TString) canVar[iVar][iPt]->GetName() + ".png";
                    canVar[iVar][iPt]->SaveAs(tmpName);
                    canVar[iVar][iPt]->Close();
                }
            }
        }

        TCanvas *c1[nVars];
        TRatioPlot *rp[nVars];
        TF1* fun1[nVars];
        TGraph* grLowRef[nVars];
        TGraph* grLowRefSpline[nVars];

        for (int iVar = 0; iVar < nVars; ++iVar) {
            tmpName=particles[iPart]+"_"+varNames[iVar]+"_ratio";

            c1[iVar] = new TCanvas(tmpName, tmpName, 900, 1100);
            gPad->SetLeftMargin(0.15);
            gPad->SetRightMargin(0.05);
            rp[iVar] = new TRatioPlot(hRatio[iVar][1], hRatio[iVar][0]);
            rp[iVar] -> SetH1DrawOpt("E");
            rp[iVar] -> SetH2DrawOpt("E");
            rp[iVar] -> SetGraphDrawOpt("EP");
//        rp[iVar] -> SetTitle("");
            rp[iVar]->Draw();

            //Y:
            rp[iVar] -> GetLowerRefYaxis() -> SetTitle("Data/Embedding");
//        rp[iVar] -> GetLowerRefGraph() -> SetMinimum(0.75);
//        rp[iVar] -> GetLowerRefGraph() -> SetMaximum(1.25);
            rp[iVar] -> GetLowerRefYaxis() -> CenterTitle(kTRUE);
            rp[iVar] -> GetLowerRefYaxis() -> SetLabelSize(0.04);
            rp[iVar] -> GetLowerRefYaxis() -> SetTitleSize(0.045);
            rp[iVar] -> GetLowerRefYaxis() -> SetTitleOffset(1.5);

            //        //X:
            rp[iVar] -> GetLowerRefXaxis() -> CenterTitle(kTRUE);
            rp[iVar] -> GetLowerRefXaxis() -> SetLabelSize(0.04);
            rp[iVar] -> GetLowerRefXaxis() -> SetTitleSize(0.04);
            rp[iVar] -> GetLowerRefXaxis() -> SetTitleOffset(1.17);

            rp[iVar] -> GetUpperPad() -> SetMargin(0.15,0.05,0.13,0.08);
            rp[iVar] -> GetUpperPad() -> cd();
            gPad->SetGrid();

            legend[iVar]->SetHeader(particles[iPart]+" "+legendTitle[iVar]);
            legend[iVar]->Draw("same");

            rp[iVar] -> GetLowerPad() -> SetMargin(0.15,0.05,0.13,0.08);
            rp[iVar] -> GetLowerPad() -> SetGrid();

            rp[iVar] -> SetLowBottomMargin(0.35);
            rp[iVar] -> SetLeftMargin(0.14);
            rp[iVar] -> SetRightMargin(0.04);
            rp[iVar] -> SetUpTopMargin(0.05);
            rp[iVar] -> SetSeparationMargin(0.03);

            grLowRef[iVar] = new TGraph();
            grLowRef[iVar] = rp[iVar] -> GetLowerRefGraph();

            grLowRefSpline[iVar] = new TGraph();

            tmpName=Form("%s_gr_fun_%s", particles[iPart].Data(), varNames[iVar].Data());
            grLowRef[iVar]->SetName(tmpName);

            fOutputs->cd();
            grLowRef[iVar]->Write();

            grLowRefSpline[iVar] = makeSpline(grLowRef[iVar]);

            tmpName=Form("%s_gr1_fun_%s", particles[iPart].Data(), varNames[iVar].Data());
            grLowRefSpline[iVar]->SetMarkerColor(46);
            grLowRefSpline[iVar]->SetLineColor(46);
            grLowRefSpline[iVar]->SetName(tmpName);
            grLowRefSpline[iVar]->Write();

//        grLowRefSpline[iVar]->SetName(tmpName);
            rp[iVar] -> GetLowerPad() -> cd();
            gPad->SetGrid();
            grLowRefSpline[iVar]->Draw("lsame");

            tmpName="img/"+(TString)c1[iVar]->GetTitle();
            c1[iVar]->SaveAs(tmpName+".png");
        }
    }

}

//________________________________________________________________________________________
TString makeHistosMC(TString inputFile, TString particle){
    TString outFileName=particle+"_tpc_sys_eff_embedding_MC.root";

    TFile *outFile = new TFile(outFileName, "recreate");
    outFile->SetCompressionSettings(0);

    TFile *soubor = new TFile(inputFile, "READ"); //open file to read
    TTree *tree = (TTree*)soubor->Get("tracks");
    TTree *treeEvents = (TTree*)soubor->Get("eventCount");
    Long64_t nEntries = tree->GetEntries();
    Long64_t nEvents = treeEvents->GetEntries();

    float geantIdRequiredMin, geantIdRequiredMax, geantId;

    Float_t gPt, gEta, gPhi, nFit, nCom, dca, dcaXY, centrality, nMcTracks, vx, vy, vz, pt, eta, hftTopo, startVtxX, startVtxY, startVtxZ;

    //MC tracks
    tree->SetBranchAddress("pt", &pt);
    tree->SetBranchAddress("eta", &eta);
    //reconstructed tracks
    tree->SetBranchAddress("gPt", &gPt);
    tree->SetBranchAddress("gEta", &gEta);
    tree->SetBranchAddress("gPhi", &gPhi);
    tree->SetBranchAddress("nFit", &nFit);
    tree->SetBranchAddress("nCom", &nCom);
    tree->SetBranchAddress("dca", &dca);
    tree->SetBranchAddress("dcaXY", &dcaXY);
    tree->SetBranchAddress("refMult", &centrality);
    tree->SetBranchAddress("geantId", &geantId);
    tree->SetBranchAddress("hftTopo", &hftTopo);
    tree->SetBranchAddress("startVtxX", &startVtxX);
    tree->SetBranchAddress("startVtxY", &startVtxY);
    tree->SetBranchAddress("startVtxZ", &startVtxZ);
    tree->SetBranchAddress("hftTopo", &hftTopo);
    //event tree
    treeEvents->SetBranchAddress("nMcTracks", &nMcTracks);
    treeEvents->SetBranchAddress("vx", &vx);
    treeEvents->SetBranchAddress("vy", &vy);
    treeEvents->SetBranchAddress("vz", &vz);

    int nPoints = sizeof(ptBins) / sizeof(Float_t);
    cout<<nPoints<<endl;

    TH1F *hDCA[nPoints-1];
    initHistosPtBins(hDCA, nPoints-1, "dca", nBinsDCA, 0., maxDCA, particle);
    TH1F *hNhitsFit[nPoints-1];
    initHistosPtBins(hNhitsFit, nPoints-1, "nHitsFit", 50, 0.5, 50.5, particle);

    TH1F *hDCAHft[nPoints-1];
    initHistosPtBins(hDCAHft, nPoints-1, "dca_hft", nBinsDCA, 0., maxDCA, particle);
    TH1F *hNhitsFitHft[nPoints-1];
    initHistosPtBins(hNhitsFitHft, nPoints-1, "nHitsFit_hft", 50, 0.5, 50.5, particle);

    TH1F *hDCAAll[nPoints-1];
    initHistosPtBins(hDCAAll, nPoints-1, "dca_all", nBinsDCA, 0., maxDCA, particle);
    TH1F *hNhitsFitAll[nPoints-1];
    initHistosPtBins(hNhitsFitAll, nPoints-1, "nHitsFit_all", 50, 0.5, 50.5, particle);

    int iTrack = 0;
    for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
        treeEvents->GetEntry(iEvent);
        iTrack += nMcTracks;

        bool hotSpot= false;
        if(vx>-0.25 && vx<-0.16 && vy>-0.25 && vy<-0.16) hotSpot=true;
//        if (!hotSpot) continue;

        for (Long64_t i = iTrack - nMcTracks; i < iTrack; i++) {
            tree->GetEntry(i);
            if (!(startVtxX==vx && startVtxY==vy && startVtxZ==vz)) continue;

            if (abs(eta) > 1.) continue; //eta of MC track
            if (abs(gEta) > 1.) continue;
            if (!(nCom > 10)) continue;
//            if (abs(dca) < 0.005) continue;

            bool goodPidGeant= false;

            if ((particle=="kaon") && ((geantId==11) || (geantId==12)))  goodPidGeant=true;
            if ((particle=="pion") && ((geantId==8) || (geantId==9)))  goodPidGeant=true;

            if (!goodPidGeant) continue;

            bool pxl1Hit = (UInt_t)hftTopo >> 0 & 0x1;
            bool pxl2Hit = (UInt_t)hftTopo >> 1 & 0x3;
            bool istHit = (UInt_t)hftTopo  >> 3 & 0x3;
            bool sstHit = (UInt_t)hftTopo >> 5 & 0x3;
            bool isHftTrack= false;
            bool isHftHit= false;
            if(pxl1Hit && pxl2Hit && (istHit || sstHit)) isHftTrack= true;
            if(pxl1Hit || pxl2Hit || istHit || sstHit) isHftHit= true;
            const int bin = getPtBin(gPt, nPoints - 1);
            if (bin == -1) continue;

            if (!isHftTrack && !isHftHit){
                if (abs(dca) < 1.)
                    hNhitsFit[bin]->Fill(nFit);
                if (abs(nFit) > 15.)
                    hDCA[bin]->Fill(dca);
            }

            if (isHftTrack){
                if (abs(dca) < 1.)
                    hNhitsFitHft[bin]->Fill(nFit);
                if (abs(nFit) > 15.)
                    hDCAHft[bin]->Fill(dca);
            }

            if (abs(dca) < 1.)
                hNhitsFitAll[bin]->Fill(nFit);
            if (abs(nFit) > 15.)
                hDCAAll[bin]->Fill(dca);
        }
    }

    outFile->cd();

    TH1F* hDcaRatio = makeRatioPlot(hDCA, "dca");
    setHistoStyle(hDcaRatio, 1, 20, "p_{T} [GeV/c]", "Ratio");
    TH1F* hNhitsFitRatio = makeRatioPlot(hNhitsFit, "nHitsFit");
    setHistoStyle(hNhitsFitRatio, 1, 20, "p_{T} [GeV/c]", "Ratio");

    TH1F* hDcaRatioHft = makeRatioPlot(hDCAHft, "dca_hft");
    setHistoStyle(hDcaRatioHft, 1, 20, "p_{T} [GeV/c]", "Ratio");
    TH1F* hNhitsFitRatioHft = makeRatioPlot(hNhitsFitHft, "nHitsFit_hft");
    setHistoStyle(hNhitsFitRatioHft, 1, 20, "p_{T} [GeV/c]", "Ratio");

    TH1F* hDcaRatioAll = makeRatioPlot(hDCAAll, "dca_all");
    setHistoStyle(hDcaRatioAll, 1, 20, "p_{T} [GeV/c]", "Ratio");
    TH1F* hNhitsFitRatioAll = makeRatioPlot(hNhitsFitAll, "nHitsFit_all");
    setHistoStyle(hNhitsFitRatioAll, 1, 20, "p_{T} [GeV/c]", "Ratio");

    hDcaRatio->Write("dca_ratio");
    hNhitsFitRatio->Write("nHitsFit_ratio");
    hDcaRatioHft->Write("dca_hft_ratio");
    hNhitsFitRatioHft->Write("nHitsFit_hft_ratio");
    hDcaRatioAll->Write("dca_all_ratio");
    hNhitsFitRatioAll->Write("nHitsFit_all_ratio");

    for (int i = 0; i < nPoints-1; ++i) {
        hNhitsFit[i]->Write();
        hDCA[i]->Write();
        hNhitsFitHft[i]->Write();
        hDCAHft[i]->Write();
        hNhitsFitAll[i]->Write();
        hDCAAll[i]->Write();
    }

    soubor->Close();
    outFile->Close();

    return outFileName;
}

//________________________________________________________________________________________
TString makeHistosData(TString inputFile, TString particle){
    TString outFileName=particle+"_tpc_sys_eff_embedding_data.root";
    TFile *outFile = new TFile(outFileName, "recreate");
    outFile->SetCompressionSettings(0);

    TFile *soubor = new TFile(inputFile, "READ"); //open file to read
    TTree *tree = (TTree*)soubor->Get("ntp_tracks");
    Long64_t nEntries = tree->GetEntries();

    Float_t gPt, gEta, gPhi, nFit, nCom, dca, dcaXY, centrality, nMcTracks, vx, vy, pt, eta, isHft, nSigmaPion, nSigmaKaon, invBetaPion, invBetaKaon, hasHftHit, isPrimaryTrk;

    tree->SetBranchAddress("pt", &pt);
    tree->SetBranchAddress("eta", &eta);
    tree->SetBranchAddress("nHitsFitTrk", &nFit);
    tree->SetBranchAddress("dca", &dca);
    tree->SetBranchAddress("dcaXy", &dcaXY);
    tree->SetBranchAddress("isHft", &isHft);
    tree->SetBranchAddress("nSigmaPion", &nSigmaPion);
    tree->SetBranchAddress("nSigmaKaon", &nSigmaKaon);
    tree->SetBranchAddress("invBetaPion", &invBetaPion);
    tree->SetBranchAddress("invBetaKaon", &invBetaKaon);
    tree->SetBranchAddress("hasHftHit", &hasHftHit);
    tree->SetBranchAddress("isPrimaryTrk", &isPrimaryTrk);

    int nPoints = sizeof(ptBins)/sizeof(Float_t);
    cout<<nPoints<<endl;

    TH1F *hDCA[nPoints-1];
    initHistosPtBins(hDCA, nPoints-1, "dca", nBinsDCA, 0., maxDCA, particle);
    TH1F *hNhitsFit[nPoints-1];
    initHistosPtBins(hNhitsFit, nPoints-1, "nHitsFit", 50, 0.5, 50.5, particle);

    TH1F *hDCAHft[nPoints-1];
    initHistosPtBins(hDCAHft, nPoints-1, "dca_hft", nBinsDCA, 0., maxDCA, particle);
    TH1F *hNhitsFitHft[nPoints-1];
    initHistosPtBins(hNhitsFitHft, nPoints-1, "nHitsFit_hft", 50, 0.5, 50.5, particle);

    TH1F *hDCAAll[nPoints-1];
    initHistosPtBins(hDCAAll, nPoints-1, "dca_all", nBinsDCA, 0., maxDCA, particle);
    TH1F *hNhitsFitAll[nPoints-1];
    initHistosPtBins(hNhitsFitAll, nPoints-1, "nHitsFit_all", 50, 0.5, 50.5, particle);

    TH1F* hPtDca1 = new TH1F("ptAllDca1", "ptAllDca1", 50, 0., 5.);
    TH1F* hPtDca1Hft = new TH1F("ptAllDca1Hft", "ptAllDca1Hft", 50, 0., 5.);
    TH1F* hPtDca15 = new TH1F("ptAllDca15", "ptAllDca15", 50, 0., 5.);
    TH1F* hPtDca15Hft = new TH1F("ptAllDca15Hft", "ptAllDca15Hft", 50, 0., 5.);

    hPtDca1->Sumw2();
    hPtDca1Hft->Sumw2();
    hPtDca15->Sumw2();
    hPtDca15Hft->Sumw2();

    int iTrack = 0;

//    for (Long64_t i=0; i < nEntries; i++) {
    for (Long64_t i=0; i < nEntries/10; i++) {
        tree->GetEntry(i);
        if (TMath::Abs(eta) > 1) continue;

        if ((particle=="kaon") && (abs(nSigmaKaon)>2)) continue;
        if ((particle=="pion") && (abs(nSigmaPion)>3)) continue;

        if (isPrimaryTrk<0) continue;

        //HFT ratio plots:
        if (abs(nFit) > 15.) {
            if (abs(dca)<1.){
                hPtDca1->Fill(pt);
                if (isHft>0) hPtDca1Hft->Fill(pt);
            }
            if (abs(dca)<1.5){
                hPtDca15->Fill(pt);
                if (isHft>0) hPtDca15Hft->Fill(pt);
            }
        }
        // end HFT ratio plots

        const int bin = getPtBin(pt, nPoints - 1);
        if (bin == -1) continue;

        bool tof=true;
//        bool tof=false;
//        if (abs(nSigmaPion)<3) {
//            if (invBetaPion < 0) tof = true;
//            else if (abs(invBetaPion) < 0.03) tof = true;
//        }
//        if (abs(nSigmaKaon)<2) {
//            if (invBetaKaon < 0) tof = true;
//            else if (abs(invBetaKaon) < 0.03) tof = true;
//        }
//        if (invBetaKaon>0 || invBetaPion>0) tof = true;


        if (!tof) continue;
        // HFT tracks
        if (isHft>0){
            if (abs(dca) < 1.)
                hNhitsFitHft[bin]->Fill(nFit);
            if (abs(nFit) > 15.)
                hDCAHft[bin]->Fill(dca);
        }

        if (abs(dca) < 1.)
            hNhitsFitAll[bin]->Fill(nFit);
        if (abs(nFit) > 15.)
            hDCAAll[bin]->Fill(dca);

        // NON-HFT tracks
        if (isHft>0) continue;
        if (hasHftHit>0) continue;

        if (abs(dca) < 1.)
            hNhitsFit[bin]->Fill(nFit);
        if (abs(nFit) > 15.)
            hDCA[bin]->Fill(dca);
    }

    outFile->cd();
    TH1F* hDcaRatio = makeRatioPlot(hDCA, "dca");
    setHistoStyle(hDcaRatio, 1, 20, "p_{T} [GeV/c]", "Ratio");
    TH1F* hNhitsFitRatio = makeRatioPlot(hNhitsFit, "nHitsFit");
    setHistoStyle(hNhitsFitRatio, 1, 20, "p_{T} [GeV/c]", "Ratio");

    TH1F* hDcaRatioHft = makeRatioPlot(hDCAHft, "dca_hft");
    setHistoStyle(hDcaRatioHft, 1, 20, "p_{T} [GeV/c]", "Ratio");
    TH1F* hNhitsFitRatioHft = makeRatioPlot(hNhitsFitHft, "nHitsFit_hft");
    setHistoStyle(hNhitsFitRatioHft, 1, 20, "p_{T} [GeV/c]", "Ratio");

    TH1F* hDcaRatioAll = makeRatioPlot(hDCAAll, "dca_all");
    setHistoStyle(hDcaRatioAll, 1, 20, "p_{T} [GeV/c]", "Ratio");
    TH1F* hNhitsFitRatioAll = makeRatioPlot(hNhitsFitAll, "nHitsFit_all");
    setHistoStyle(hNhitsFitRatioAll, 1, 20, "p_{T} [GeV/c]", "Ratio");

    hDcaRatio->Write("dca_ratio");
    hNhitsFitRatio->Write("nHitsFit_ratio");
    hDcaRatioHft->Write("dca_hft_ratio");
    hNhitsFitRatioHft->Write("nHitsFit_hft_ratio");
    hDcaRatioAll->Write("dca_all_ratio");
    hNhitsFitRatioAll->Write("nHitsFit_all_ratio");

    for (int i = 0; i < nPoints-1; ++i) {
        hNhitsFit[i]->Write();
        hDCA[i]->Write();
        hNhitsFitHft[i]->Write();
        hDCAHft[i]->Write();
        hNhitsFitAll[i]->Write();
        hDCAAll[i]->Write();
    }


    hPtDca15Hft->Divide(hPtDca15);
    hPtDca1Hft->Divide(hPtDca1);

    ///--------------------------------------------------------------------------------
    hPtDca15Hft->Write();
    hPtDca1Hft->Write();
    setHistoStyle(hPtDca15Hft, 1, 24, "p_{T} [GeV/c]", "HFT ratio");
    setHistoStyle(hPtDca1Hft, kMagenta+2, 20, "p_{T} [GeV/c]", "HFT ratio");


    TCanvas *c1 = new TCanvas("hftRatio", "hftRatio", 900, 1100);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);

    TRatioPlot *rp = new TRatioPlot(hPtDca15Hft, hPtDca1Hft);
    rp -> SetH1DrawOpt("E");
    rp -> SetH2DrawOpt("E");
    rp -> SetGraphDrawOpt("EP");
    rp->Draw();

    //Y:
    rp -> GetLowerRefYaxis() -> SetTitle("DCA<1.5 / DCA<1");
    rp -> GetLowerRefGraph() -> SetMinimum(0.75);
    rp -> GetLowerRefGraph() -> SetMaximum(1.25);
    rp -> GetLowerRefYaxis() -> CenterTitle(kTRUE);
    rp -> GetLowerRefYaxis() -> SetLabelSize(0.04);
    rp -> GetLowerRefYaxis() -> SetTitleSize(0.045);
    rp -> GetLowerRefYaxis() -> SetTitleOffset(1.5);

    rp -> GetLowerRefXaxis() -> CenterTitle(kTRUE);
    rp -> GetLowerRefXaxis() -> SetLabelSize(0.04);
    rp -> GetLowerRefXaxis() -> SetTitleSize(0.04);
    rp -> GetLowerRefXaxis() -> SetTitleOffset(1.17);

    rp -> GetUpperPad() -> SetMargin(0.15,0.05,0.13,0.08);
    rp -> GetUpperPad() -> cd();

    TLegend *legend=new TLegend(0.596,0.114, 0.76, 0.312,"","brNDC");
    legend->SetFillStyle(0);
    legend->SetLineColor(0);
    legend->SetTextSize(0.035);
    legend->AddEntry(hPtDca15Hft, "DCA < 1.5 cm", "pl");
    legend->AddEntry(hPtDca1Hft, "DCA < 1 cm", "pl");
    legend->Draw("same");

    rp -> GetLowerPad() -> SetMargin(0.15,0.05,0.13,0.08);

    c1->SaveAs(".png");


//    soubor->Close();
//    outFile->Close();

    return outFileName;
}

//________________________________________________________________________________________
int getPtBin(float pt, int nPtBins){
    for (int i = 0; i < nPtBins-1; i++){
        if ((pt >= ptBins[i]) && (pt < ptBins[i+1]))
            return i;
    }
    return -1;
}

//________________________________________________________________________________________
void initHistosPtBins(TH1F** histoArray, int nHistos, TString histoName, Int_t nBins, Double_t axMin, Double_t axMax, TString particle){
    for (int i = 0; i < nHistos; ++i) {
        TString tmp = Form("%s_%.2f_%.2f", histoName.Data(), ptBins[i], ptBins[i+1]);

        histoArray[i] = new TH1F(tmp, tmp, nBins, axMin, axMax);
        histoArray[i]->SetName(tmp);

        tmp = Form("%s, %.2f < p_{T} < %.2f GeV/c", particle.Data(), ptBins[i], ptBins[i+1]);
        histoArray[i]->SetTitle(tmp);
        histoArray[i]->Sumw2();
    }
}

//________________________________________________________________________________________
TH1F* makeRatioPlot(TH1F** histoArray, TString histoName){
    int nPoints = sizeof(ptBins) / sizeof(Float_t);
    cout<<nPoints<<endl;

    TH1F* histo = new TH1F(histoName, histoName, nPoints-1, ptBins);
    for (int i = 0; i < nPoints-2; ++i) {
        histoArray[i]->ClearUnderflowAndOverflow();
        histoArray[i]->Scale(1/histoArray[i]->GetEntries());

        if (histoName=="dca" || histoName=="dca_hft" || histoName=="dca_all"){
            double errDca1, errDca15;
            double integralDca1=histoArray[i]->IntegralAndError(1,
                                                                histoArray[i]->FindBin(1.),
                                                                errDca1);
            errDca1=errDca1/integralDca1;

            double integralDca15=histoArray[i]->IntegralAndError(1,
                                                                 histoArray[i]->FindBin(maxDCA),
                                                                 errDca15);
            errDca15=errDca15/integralDca15;
            if (integralDca15>0.) {
                histo->SetBinContent(i+1, integralDca1 / integralDca15);
                histo->SetBinError(i+1, integralDca1 / integralDca15 * sqrt(errDca15 * errDca15 + errDca1 * errDca1));
            }
        }

        if (histoName=="nHitsFit" || histoName=="nHitsFit_hft" || histoName=="nHitsFit_all"){
            double err15, err20;
            double integral15=histoArray[i]->IntegralAndError(histoArray[i]->FindBin(15.),
                                                              histoArray[i]->FindBin(50.5),
                                                              err15);
            err15=err15/integral15;

            double integral20=histoArray[i]->IntegralAndError(histoArray[i]->FindBin(20.),
                                                              histoArray[i]->FindBin(50.5),
                                                              err20);
            err20=err20/integral20;

            histo->SetBinContent(i+1, integral20/integral15);
            histo->SetBinError(i+1, integral20/integral15 * sqrt(err20*err20 + err15*err15));
        }
    }

    return histo;
}

//____________________________________________________________________________________________________________________________
void setHistoStyle(TH1F* histo, Int_t color, Int_t marker, TString titleX, TString titleY) {
    histo->Sumw2();

    histo->SetMarkerStyle(marker);
    histo->SetLineColor(color);
    histo->SetLineWidth(2);
    histo->SetMarkerColor(color);
    histo->SetMarkerSize(2);

    histo->SetStats(0);
    histo->GetYaxis()->SetTitle(titleY);
    histo->GetYaxis()->SetTitleOffset(1.5);
    histo->GetYaxis()->SetLabelSize(0.045);
    histo->GetYaxis()->SetTitleSize(0.045);
    histo->GetYaxis()->SetLabelFont(42);
    histo->GetYaxis()->SetTitleFont(42);
    histo->GetYaxis()->CenterTitle(kTRUE);
//    histo->GetYaxis()->SetMaxDigits(2);

    histo->GetXaxis()->SetTitle(titleX);
    histo->GetXaxis()->SetLabelFont(42);
    histo->GetXaxis()->SetTitleFont(42);
    histo->GetXaxis()->SetTitleSize(0.045);
    histo->GetXaxis()->SetLabelSize(0.045);
    histo->GetXaxis()->CenterTitle(kTRUE);
}

//___________________________________--
//void makeSpline(TGraph* gr) {
TGraph* makeSpline(TGraph* gr) {
    const int nPtPoints=100;
    TGraph *grOut = new TGraph();

    for (int i = 0; i < nPtPoints; ++i) {
        double pt = 0.15+i*(2.5-0.15)/nPtPoints;
        double y = gr->Eval(pt, nullptr, "");
        grOut->SetPoint(i, pt, y);
//        gr->SetPoint(i, pt, y);
    }
    return grOut;
//    return;
}