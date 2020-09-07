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
#include <ctime>

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
    clock_t start = clock(); // getting starting time

    Int_t colors[] = {1,46,8,4,5,6,7,9,41,42,44,45,16};
    TString fileNameData = "/media/lukas/1E4183350724EC9C/1.hists.root";
//    TString fileNameData = "/media/lukas/1E4183350724EC9C/tracks.sim.vzvpd3.hotspot.root";
//    TString fileNameData = "/home/lukas/work/dmesons/Dmaker_ndAu/Dmaker_dAu/analyse/simInputsPrepare/toHadd/ntp.ratio.1404.half.rootprimary.root";
    auto* dataF = new TFile(fileNameData,"r");
    auto* dataOut = new TFile("ratio.root","recreate");
    auto* ntpData = (TNtuple*)dataF -> Get("ntp_tracks");
//    ntpData->LoadBaskets(5000000000);
    Float_t  pt, dca, isHFT, isTOF, zdc, nHitsFit, refMult, nHftTracks, nTofTracks, particleId, isPrimary, dcaxy, dcaz, runId, eventId, nSigmaPion,nSigmaKaon,invBetaPion,invBetaKaon, eta, phi;

    ntpData->SetBranchAddress("runId", &runId);
    ntpData->SetBranchAddress("eventId", &eventId);
    ntpData->SetBranchAddress("pt", &pt);
    ntpData->SetBranchAddress("eta", &eta);
    ntpData->SetBranchAddress("phi", &phi);
    ntpData->SetBranchAddress("dca", &dca);
    ntpData->SetBranchAddress("dcaXy", &dcaxy);
    ntpData->SetBranchAddress("dcaZ", &dcaz);
    ntpData->SetBranchAddress("isHft", &isHFT);
//    ntpData->SetBranchAddress("isTOF", &isTOF);
//    ntpData->SetBranchAddress("zdc", &zdc);
    ntpData->SetBranchAddress("nHitsFitTrk", &nHitsFit);
    ntpData->SetBranchAddress("multiplicity", &refMult);
    ntpData->SetBranchAddress("isPrimaryTrk", &isPrimary);
    ntpData->SetBranchAddress("nHftTracks", &nHftTracks);
    ntpData->SetBranchAddress("nTofTracks", &nTofTracks);
//    ntpData->SetBranchAddress("particleId", &particleId);
    ntpData->SetBranchAddress("invBetaKaon", &invBetaKaon);
    ntpData->SetBranchAddress("invBetaPion", &invBetaPion);
    ntpData->SetBranchAddress("nSigmaKaon", &nSigmaKaon);
    ntpData->SetBranchAddress("nSigmaPion", &nSigmaPion);

//    Float_t dcaCuts[]={0.1, 0.3, 0.5, 0.7, 1, 1.5};
//    const int nDCAcuts = sizeof(dcaCuts) / sizeof(Float_t);

    const int nDCAcuts=6;
    const float dcaCuts[nDCAcuts+1]={0., 0.3, 0.5, 0.7, 1.0, 1.2, 1.5};

    const int nPtcuts=9;
    const float ptCuts[nPtcuts+1]={0.15, 0.3, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0, 3.0, 5.0};

    TString particleName[]={"pion","kaon"};
//    TString trackKinds[]=    {"","primary",     "TOFmatched", "primary_TOFmatched", "_nHftTracks1", "_nHftTracks2","_nTofTracks2"};
    TString trackKinds[]=    {"",  "TPC",   "TPC_TOF", "TPC_hybridTOF", "TPC_TOFPid"};
//    TString trackKindsCuts[]={"","isPrimary>0", "isTOF>0",    "isPrimary>0 && isTOF>0"};
    const int nKinds = sizeof(trackKinds) / sizeof(TString);

    TProfile *prRatio[nKinds][nDCAcuts];

    TH1D* hPtAll[2][nKinds][2];
    TH1D* hPt[2][nKinds][nDCAcuts][2];
    TH1D* hEta[2][nKinds][nPtcuts][2];
    TH1D* hPhi[2][nKinds][nPtcuts][2];

    TH1D* hDCAPtBins[2][nKinds][nPtcuts][2];
    TH1D* hDCAxyPtBins[2][nKinds][nPtcuts][2];
    TH1D* hDCAzPtBins[2][nKinds][nPtcuts][2];

    TString nameHisto;
    TString hftText[]={"","_hft"};

    Double_t maxDca=0.1;
    for (int i = 0; i < nKinds; ++i) {
        for (int k = 0; k < 2; ++k) {
            for (int iPart = 0; iPart < 2; ++iPart) {
                nameHisto = Form("hPt_%s_%i_%s", trackKinds[i].Data(), k, particleName[iPart].Data());
                hPtAll[k][i][iPart] = new TH1D(nameHisto, nameHisto, 100, 0, 10);
                hPtAll[k][i][iPart]->Sumw2();

                for (int j = 0; j < nDCAcuts; ++j) {
                    nameHisto = Form("hPt_dca%.1f_%.1f_%s_%i_%s", dcaCuts[j], dcaCuts[j + 1], trackKinds[i].Data(), k, particleName[iPart].Data());
                    hPt[k][i][j][iPart] = new TH1D(nameHisto, nameHisto, 100, 0, 10);
                    hPt[k][i][j][iPart]->Sumw2();
                }

                for (int j = 0; j < nPtcuts; ++j) {
                    nameHisto = Form("hEta_pt%.1f_%.1f_%s_%i_%s", ptCuts[j], ptCuts[j + 1], trackKinds[i].Data(), k, particleName[iPart].Data());
                    hEta[k][i][j][iPart] = new TH1D(nameHisto, nameHisto, 100, -1, 1);
                    hEta[k][i][j][iPart]->Sumw2();

                    nameHisto = Form("hPhi_pt%.1f_%.1f_%s_%i_%s", ptCuts[j], ptCuts[j + 1], trackKinds[i].Data(), k, particleName[iPart].Data());
                    hPhi[k][i][j][iPart] = new TH1D(nameHisto, nameHisto, 100, -3.2, 3.2);
                    hPhi[k][i][j][iPart]->Sumw2();

                    nameHisto = Form("hDCA_pt%.1f_%.1f_%s_%s%s", ptCuts[j], ptCuts[j + 1], trackKinds[i].Data(), particleName[iPart].Data(), hftText[k].Data());
                    hDCAPtBins[k][i][j][iPart] = new TH1D(nameHisto, nameHisto, 100, 0, maxDca);
                    hDCAPtBins[k][i][j][iPart]->Sumw2();

                    nameHisto = Form("hDCAxy_pt%.1f_%.1f_%s_%s%s", ptCuts[j], ptCuts[j + 1], trackKinds[i].Data(), particleName[iPart].Data(), hftText[k].Data());
                    hDCAxyPtBins[k][i][j][iPart] = new TH1D(nameHisto, nameHisto, 100, -1*maxDca, maxDca);
                    hDCAxyPtBins[k][i][j][iPart]->Sumw2();

                    nameHisto = Form("hDCAz_pt%.1f_%.1f_%s_%s%s", ptCuts[j], ptCuts[j + 1], trackKinds[i].Data(), particleName[iPart].Data(), hftText[k].Data());
                    hDCAzPtBins[k][i][j][iPart] = new TH1D(nameHisto, nameHisto, 100, -1*maxDca, maxDca);
                    hDCAzPtBins[k][i][j][iPart]->Sumw2();
                }
            }

        }
//        nameHisto = Form("hPtHft_dca%.1f_%.1f_%s", dcaCuts[j], dcaCuts[j + 1], trackKinds[i].Data());
//        prRatio[i][j] = new TProfile(nameHisto,nameHisto,100,0,10,0,20);
//            prRatio[i][j] = new TH2F(nameHisto,nameHisto,100,0.,10.,100,0.,2.);
//        prRatio[i][j]->Sumw2();
    }


    Float_t previousEvent=-999;
    Float_t previousRun=-999;
    int nEvts=0;
    cout<<"total entr "<<ntpData->GetEntries()<<endl;

    Float_t isHybridTOF[]={-1,-1};
    Float_t isTOFtrk[]={-1,-1};
    Float_t isTOFPidtrk[]={-1,-1};
    Float_t isTPCtrk[]={-1,-1};

    for (int k = 0; k < ntpData->GetEntries(); ++k) {
//    for (int k = 0; k < ntpData->GetEntries()/10; ++k) {
        ntpData->GetEntry(k);
        if (nHitsFit<20) continue;
        if (isPrimary<1) continue;

        if (abs(nSigmaPion)<1) isTPCtrk[0]=1; else isTPCtrk[0]=-1;
        if (abs(invBetaPion)<20) isTOFtrk[0]=1; else isTOFtrk[0]=-1;
        if (abs(invBetaPion)<0.02) isTOFPidtrk[0]=1; else isTOFPidtrk[0]=-1;
        if ((isTOFtrk[0]==-1) ||  (abs(invBetaPion)<0.03)) isHybridTOF[0]=1; else isHybridTOF[0]=-1;

        if (abs(nSigmaKaon)<1) isTPCtrk[1]=1; else isTPCtrk[1]=-1;
        if (abs(invBetaKaon)<20) isTOFtrk[1]=1; else isTOFtrk[1]=-1;
        if (abs(invBetaKaon)<0.02) isTOFPidtrk[1]=1; else isTOFPidtrk[1]=-1;
        if ((isTOFtrk[1]==-1) || (abs(invBetaKaon)<0.03)) isHybridTOF[1]=1; else isHybridTOF[1]=-1;

//_________TPROFILE_________
//        if (previousEvent!=eventId || previousRun!=runId) { //new event
//            previousEvent = eventId;
//            previousRun = runId;
//            nEvts++;
//            for (int i = 0; i < nKinds; ++i) {
//                for (int j = 0; j < nDCAcuts; ++j) {
//                    hPtHft[1][i][j]->Divide(hPt[1][i][j]);
//                    for (int l = 1; l < hPt[1][i][j]->GetNbinsX()+1; ++l) {
//                        prRatio[i][j]->Fill(hPtHft[1][i][j]->GetBinCenter(l), hPtHft[1][i][j]->GetBinContent(l));
//                    }
//                    hPt[1][i][j]->Reset();
//                    hPtHft[1][i][j]->Reset();
//                }
//            }
//        }
//_________END TPROFILE_________


        for (int j = 0; j < 2; ++j) {
            if (j==1 && isHFT<0) continue;

            for (int iPart = 0; iPart < 2; ++iPart) {
                //_______no Binnning_______
                hPtAll[j][0][iPart]->Fill(pt);
                if (isTPCtrk[iPart]>0){
                    hPtAll[j][1][iPart]->Fill(pt);
                    if (isTOFtrk[iPart]>0) {
                        hPtAll[j][2][iPart]->Fill(pt);
                    }
                    if (isHybridTOF[iPart]>0) {
                        hPtAll[j][3][iPart]->Fill(pt);
                    }
                    if (isTOFPidtrk[iPart]>0) {
                        hPtAll[j][4][iPart]->Fill(pt);
                    }
                }

                //_______DCA binning____________
                for (int i = 0; i < nDCAcuts; ++i) {
                    if (dca>dcaCuts[i] && dca<dcaCuts[i+1]) {
                        hPt[j][0][i][iPart]->Fill(pt);
                        if (isTPCtrk[iPart]>0){
                            hPt[j][1][i][iPart]->Fill(pt);
                            if (isTOFtrk[iPart]>0) {
                                hPt[j][2][i][iPart]->Fill(pt);
                            }
                            if (isHybridTOF[iPart]>0) {
                                hPt[j][3][i][iPart]->Fill(pt);
                            }
                            if (isTOFPidtrk[iPart]>0) {
                                hPt[j][4][i][iPart]->Fill(pt);
                            }
                        }
                    }
                }

                //_______PT Binning____________
                for (int i = 0; i < nPtcuts; ++i) {
                    if (pt>ptCuts[i] && pt<ptCuts[i+1]) {
                        hEta[j][0][i][iPart]->Fill(eta);
                        hPhi[j][0][i][iPart]->Fill(phi);
                        hDCAPtBins[j][0][i][iPart]->Fill(dca);
                        hDCAxyPtBins[j][0][i][iPart]->Fill(dcaxy);
                        hDCAzPtBins[j][0][i][iPart]->Fill(dcaz);

                        if (isTPCtrk[iPart]>0){
                            hEta[j][1][i][iPart]->Fill(eta);
                            hPhi[j][1][i][iPart]->Fill(phi);
                            hDCAPtBins[j][1][i][iPart]->Fill(dca);
                            hDCAxyPtBins[j][1][i][iPart]->Fill(dcaxy);
                            hDCAzPtBins[j][1][i][iPart]->Fill(dcaz);

                            if (isTOFtrk[iPart]>0) {
                                hEta[j][2][i][iPart]->Fill(eta);
                                hPhi[j][2][i][iPart]->Fill(phi);
                                hDCAPtBins[j][2][i][iPart]->Fill(dca);
                                hDCAxyPtBins[j][2][i][iPart]->Fill(dcaxy);
                                hDCAzPtBins[j][2][i][iPart]->Fill(dcaz);
                            }
                            if (isHybridTOF[iPart]>0) {
                                hEta[j][3][i][iPart]->Fill(eta);
                                hPhi[j][3][i][iPart]->Fill(phi);
                                hDCAPtBins[j][3][i][iPart]->Fill(dca);
                                hDCAxyPtBins[j][3][i][iPart]->Fill(dcaxy);
                                hDCAzPtBins[j][3][i][iPart]->Fill(dcaz);
                            }
                            if (isTOFPidtrk[iPart]>0) {
                                hEta[j][4][i][iPart]->Fill(eta);
                                hPhi[j][4][i][iPart]->Fill(phi);
                                hDCAPtBins[j][4][i][iPart]->Fill(dca);
                                hDCAxyPtBins[j][4][i][iPart]->Fill(dcaxy);
                                hDCAzPtBins[j][4][i][iPart]->Fill(dcaz);
                            }
                        }
                    }
                }
            }
        }
    }

    cout<<"nevents "<<nEvts<<endl;
    for (int i = 0; i < nKinds; ++i) {
        for (int iPart = 0; iPart < 2; ++iPart) {
            hPtAll[1][i][iPart]->Divide(hPtAll[0][i][iPart]);
            dataOut->cd();
            hPtAll[1][i][iPart]->Write();

            for (int j = 0; j < nDCAcuts-2; ++j) {
                hPt[1][i][j][iPart]->Divide(hPt[0][i][j][iPart]);
                cout << hPt[1][i][j][iPart]->Integral() / 100 << endl;

                dataOut->cd();
                hPt[1][i][j][iPart]->Write();
            }

            for (int j = 0; j < nPtcuts; ++j) {
                hEta[1][i][j][iPart]->Divide(hEta[0][i][j][iPart]);
                hPhi[1][i][j][iPart]->Divide(hPhi[0][i][j][iPart]);

                dataOut->cd();
                hEta[1][i][j][iPart]->Write();
                hPhi[1][i][j][iPart]->Write();

                for (int kHft = 0; kHft < 2; ++kHft) {
                    hDCAPtBins[kHft][i][j][iPart]->Write();
                    hDCAxyPtBins[kHft][i][j][iPart]->Write();
                    hDCAzPtBins[kHft][i][j][iPart]->Write();
                }
            }
        }
    }

    dataOut->Close();

    double duration = (double) (clock() - start) / (double) CLOCKS_PER_SEC;
    cout << "****************************************** " << endl;
    cout << "Time needed " << duration << " s" << endl;
}

//___________________________________________________________________________________________
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