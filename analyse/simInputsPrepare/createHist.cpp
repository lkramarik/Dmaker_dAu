//
// Created by lukas on 3.4.2018.
//
#include "TH3D.h"
#include "TH3F.h"
#include "TFile.h"
#include "TH2D.h"
#include "TString.h"

#include <iostream>
using namespace std;

const int nmultEdge = 7;

TH1D* h1Vz[nmultEdge];
TH1D* h1ZdcX[nmultEdge];

void createHist() {
    const Int_t nParticles = 2; //ok
    int const nVzs = 6; //ok
    int const nEtas = 10; //ok
//    const Int_t nPtBins = 19;
    const Int_t nPhi = 11; //ok

    // input file and output file
    TFile fDca1("2101.hists.root");
//    TFile fDca1("ratio.hists.0810.root");
    TFile *outRatioPion = new TFile("hftratio_vs_pt_dAu_pion.root", "RECREATE");
    outRatioPion->SetCompressionSettings(0);
    TFile *outRatioKaon = new TFile("hftratio_vs_pt_dAu_kaon.root", "RECREATE");
    outRatioKaon->SetCompressionSettings(0);


    TFile *outHist2d = new TFile("2d.root", "RECREATE");
//    TFile *outHist2dLowStats = new TFile("2d_badstats.root", "RECREATE");
    TFile *outEvent = new TFile("inputs.event.root", "RECREATE");
    outEvent->SetCompressionSettings(0);

    int multEdge[nmultEdge + 1] = {0, 4, 8, 12, 16, 20, 24, 200};

//    const Double_t ptEdge[nPtBins + 1] = {0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.5, 4.0, 6.0, 12.0};
//    const int m_nZdc = 5;
//    float const m_zdcEdge[m_nZdc+1] = {0,50,90,130,170,210};

    const int m_nZdc = 2;
    float const m_zdcEdge[m_nZdc+1] = {0,150,210};

//    for(int iParticle = 0; iParticle < nParticles; ++iParticle){
//        for (int iEta = 0; iEta < nEtas; ++iEta){
//            for (int iVz = 0; iVz < nVzs; ++iVz){
//                for(int iPhi = 0; iPhi < nPhi; ++iPhi){

    //VZ and ZDC test, this is directly in simulation
    TH3F* mh3VzZdcMult = (TH3F*)fDca1.Get("mh3VzZdcMult");
    outEvent->cd();
    mh3VzZdcMult->Write();

    int binVzmin = 1;
    int binVzup = mh3VzZdcMult->GetXaxis()->GetNbins(); //ok
    int binZDCmin = 1;
    int binZDCmax = mh3VzZdcMult->GetYaxis()->GetNbins(); //ok

    for (int ii = 0; ii < nmultEdge; ++ii)   {
        int binMultmin = mh3VzZdcMult->GetZaxis()->FindBin(multEdge[ii]);
        int binMultmax = mh3VzZdcMult->GetZaxis()->FindBin(multEdge[ii+1]);
        cout<<multEdge[ii]<<" "<<multEdge[ii+1]<<endl;
        cout<<binMultmin<<" "<<binMultmax<<endl;
        cout<<mh3VzZdcMult->GetZaxis()->GetNbins()<<endl;
        outEvent->cd();

        h1Vz[ii] = mh3VzZdcMult -> ProjectionX("_px",binZDCmin, binZDCmax, binMultmin, binMultmax, ""); //vz zdc
        h1Vz[ii]->SetDirectory(0);
        h1Vz[ii]->Write(Form("vz_mult_%i_%i", (int)multEdge[ii], (int)multEdge[ii+1]));
        h1ZdcX[ii] = mh3VzZdcMult -> ProjectionY("_py",binVzmin, binVzup, binMultmin, binMultmax, ""); //vz zdc
        h1ZdcX[ii]->SetDirectory(0);
        h1ZdcX[ii]->Scale(1/h1ZdcX[ii]->GetEntries());
        h1ZdcX[ii]->Write(Form("zdc_mult_%i_%i", (int)multEdge[ii], (int)multEdge[ii+1]));
    }

    TH1D* hrefMult = mh3VzZdcMult -> ProjectionZ("_pz",binVzmin, binVzup, binZDCmin, binZDCmax, "");
    hrefMult->Draw();
    hrefMult->Write("hrefMult");


    for (int iParticle = 0; iParticle < 2; ++iParticle) {
        for (int iEta = 0; iEta < nEtas; ++iEta) {
            for (int iVz = 0; iVz < nVzs; ++iVz) {
                for (int iPhi = 0; iPhi < nPhi; ++iPhi) {
                    // Getting 3D histogram
                    TH3F *hist3D = 0;
                    const char *h3dName = Form("h3_hft_mult_pt_p%d_eta%d_vz%d_phi%d", iParticle, iEta, iVz, iPhi);
                    hist3D = (TH3F*)(fDca1.Get(h3dName));
                    if (!hist3D) {
                        std::cout << "histogram \"" << h3dName << "\" not found." << endl;
//                        return;
                        continue;
                    }
                    TH3F *hist3Dtpc = 0;
                    const char *h3dNametpc = Form("h3_tpc_mult_pt_p%d_eta%d_vz%d_phi%d", iParticle, iEta, iVz, iPhi);
                    hist3Dtpc = (TH3F*)(fDca1.Get(h3dNametpc));
                    if (!hist3Dtpc) {
                        std::cout << "histogram \"" << h3dNametpc << "\" not found." << endl;
//                        return;
                        continue;
                    }

                    TH2F *hist2Dtpc = (TH2F*)hist3Dtpc->Project3D("xze"); // result: y = pt, x = ZDC; skipping multiplicity
                    TH2F *hist2D = (TH2F*)hist3D->Project3D("xze");

                    TString name = Form("h2_tpc_zdc_pt_p%d_eta%d_vz%d_phi%d", iParticle, iEta, iVz, iPhi);
                    TString nameHft = Form("h2_hft_zdc_pt_p%d_eta%d_vz%d_phi%d", iParticle, iEta, iVz, iPhi);

                    hist2Dtpc->SetNameTitle(name, name);
                    hist2D->SetNameTitle(nameHft, nameHft);

                    outHist2d->cd();
                    hist2Dtpc->Write();
                    hist2D->Write();
//                    hist2D -> Draw();
                    int binlow, binup;

                    for (int iZDC = 0; iZDC < m_nZdc; ++iZDC) {
//                        for (int iCent = 0; iCent < nmultEdge; ++iCent) {
                            binlow = hist2D->GetXaxis()->FindBin(m_zdcEdge[iZDC]);
                            binup = hist2D->GetXaxis()->FindBin(m_zdcEdge[iZDC + 1]);
//                            cout << m_zdcEdge[iZDC] << " " << m_zdcEdge[iZDC + 1] << endl;
//                            cout << binlow << " " << binup << endl;
                            TH1D *hist1d = hist2D->ProjectionY("_py", binlow, binup, "");
                            TH1D *hist1dtpc = hist2Dtpc->ProjectionY("_py", binlow, binup, "");
                            hist1d->Divide(hist1dtpc);
                            hist1d->SetNameTitle(Form("h_hftratio_p%d_eta%d_vz%d_phi%d_z%d", iParticle, iEta, iVz, iPhi, iZDC), Form("h_hftratio_p%d_eta%d_vz%d_phi%d_z%d", iParticle, iEta, iVz, iPhi, iZDC));
                            if (iParticle == 0) {
                                outRatioPion->cd();
                                hist1d->Write();
                            }

                            if (iParticle == 1) {
                                outRatioKaon->cd();
                                hist1d->Write();
                            }


//                        }
                    }




                    // index of pT bin in 3D histogram
//                    int iPt3d = 1;

//                    // slicing 3D histogram into 2D slices
//                    for (int iCent = 2; iCent < 11; iCent++){
//                        // creating the 2D histogram
//                        const char *h1dName = Form("mh1HFT1PtCentPartEtaVzPhi_%i_%i_%i_%i_%i", iParticle, iEta, iVz, iPhi, iCent-2);
//                        const char *h1dNametpc = Form("mh1TPC1PtCentPartEtaVzPhi_%i_%i_%i_%i_%i", iParticle, iEta, iVz, iPhi, iCent-2);
//
//                        //                        char *cent = (char*)" ";
////                        if(iCent-2 == 8) cent = (char*)"_Cent0-5%";
////                        if(iCent-2 == 7) cent = (char*)"_Cent5-10%";
////                        if(iCent-2 == 6) cent = (char*)"_Cent10-20%";
////                        if(iCent-2 == 5) cent = (char*)"_Cent20-30%";
////                        if(iCent-2 == 4) cent = (char*)"_Cent30-40%";
////                        if(iCent-2 == 3) cent = (char*)"_Cent40-50%";
////                        if(iCent-2 == 2) cent = (char*)"_Cent50-60%";
////                        if(iCent-2 == 1) cent = (char*)"_Cent60-70%";
////                        if(iCent-2 == 0) cent = (char*)"_Cent70-80%";
//
//                        const char *h1dTitle = Form("%s%s", h1dName, cent);
//                        TH1D *hist1D = new TH1D(h1dName, h1dTitle, nPtBins, ptEdge);
//                        hist1D->GetXaxis()->SetTitle(xyTitle);
//
//                        const char *h1dTitletpc = Form("%s%s", h1dNametpc, cent);
//                        TH1D *hist1Dtpc = new TH1D(h1dNametpc, h1dTitletpc, nPtBins, ptEdge);
//
//                        for (int ipT = 1; ipT < nPtBins+1; ipT++){
//                            hist1D->SetBinContent(ipT, hist2D->GetBinContent(ipT,iCent));
//                            hist1D->SetBinError(ipT, hist2D->GetBinError(ipT,iCent));
//                            hist1Dtpc->SetBinContent(ipT, hist2Dtpc->GetBinContent(ipT,iCent));
//                            hist1Dtpc->SetBinError(ipT, hist2Dtpc->GetBinError(ipT,iCent));
//                        }
//                        //hist1D->Sumw2();
//                        //hist1Dtpc->Sumw2();
//
//                        hist1D->Divide(hist1Dtpc);
//
//                        outHistF->cd();
//                        hist1D->Write();
//                    }
//                }
//            }
//            cout<<"Finished writing histograms for eta " << iEta << endl;
//        }
//    }
//                    cout << "Done..." << endl;
                }
            }
        }
    }
}