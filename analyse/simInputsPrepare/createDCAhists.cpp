//
// Created by lukas on 9.4.2018.
//

#include "TH3.h"
#include "TFile.h"
#include "TH2.h"
#include "TString.h"

#include <iostream>
using namespace std;

void createDCAhists() {
    // contstants
    const Int_t nParticles = 2; //ok
    int const nVzs = 4; //ok
    const double epsilon = 0.000001;
    const int nEtas = 3; //ok

    // input file and output file
//    TFile fDca1("dca.hists.0810.root");
    TFile fDca1("dca.2802.root");
    TFile *outHistF = new TFile("dcaxy_vs_dcaz.root", "RECREATE");
    outHistF->SetCompressionSettings(0); //needed to open file in ROOT5
    TFile *outEvent = new TFile("inputs.event.root", "RECREATE");
    outEvent->SetCompressionSettings(0);




//    TFile *outHist2dLowStats = new TFile("2d_badstats.root", "RECREATE");

//    const int nZDC = 5;
//    float const m_zdcEdge[nZDC+1] = {0,50,90,130,170,210};

    const int nZDC = 2;
    float const m_zdcEdge[nZDC+1] = {0,150,210};

    const int m_nmultEdge = 7;
    float const m_multEdge[m_nmultEdge+1] = {0, 4, 8, 12, 16, 20, 24, 200};

//    const int m_nmultEdge = 1; //7
//    float const m_multEdge[m_nmultEdge+1] = {0, 200}; //currently not used in dca

    const int nPtBins = 12;
    float const ptEdge[nPtBins + 1] = {0.15, 0.3, 0.4, 0.5, 0.6, 0.8, 1., 1.25, 1.5, 2., 3., 5., 12.};

    const int m_nDcasDca = 148;
    float const  m_DcaEdgeDca[m_nDcasDca + 1] =   {
            -1.5, -1.2, -1 , -0.96 , -0.92 , -0.88 , -0.84 , -0.8 , -0.76 , -0.72 , -0.68 , -0.64 , -0.6 , -0.56 , -0.52 , -0.48 , -0.44 , -0.4 , -0.36 , -0.32 , -0.28 , -0.24 , -0.2 , -0.16 , -0.12 ,  -0.08,
            -0.078 , -0.075 , -0.072 , -0.069 , -0.066 , -0.063 , -0.06 , -0.057 , -0.054 , -0.051 , -0.048 , -0.045 , -0.042 , -0.039 , -0.036 , -0.033 , -0.03 , -0.027 , -0.024 , -0.021 , -0.018 , -0.015 , -0.012 ,
            -0.01 , -0.0096 , -0.0092 , -0.0088 , -0.0084 , -0.008 , -0.0076 , -0.0072 , -0.0068 , -0.0064 , -0.006 , -0.0056 , -0.0052 , -0.0048 , -0.0044 , -0.004 , -0.0036 , -0.0032 , -0.0028 , -0.0024 , -0.002 , -0.0016 , -0.0012 , -0.0008 , -0.0004 , 0 , 0.0004 , 0.0008 , 0.0012 , 0.0016 , 0.002 , 0.0024 , 0.0028 , 0.0032 , 0.0036 , 0.004 , 0.0044 , 0.0048 , 0.0052 , 0.0056 , 0.006 , 0.0064 , 0.0068 , 0.0072 , 0.0076 , 0.008 , 0.0084 , 0.0088 , 0.0092 , 0.0096 , 0.01 ,
            0.012 , 0.015 , 0.018 , 0.021 , 0.024 , 0.027 , 0.03 , 0.033 , 0.036 , 0.039 , 0.042 , 0.045 , 0.048 , 0.051 , 0.054 , 0.057 , 0.06 , 0.063 , 0.066 , 0.069 , 0.072 , 0.075 , 0.078 ,
            0.08 , 0.12 , 0.16 , 0.2 , 0.24 , 0.28 , 0.32 , 0.36 , 0.4 , 0.44 , 0.48 , 0.52 , 0.56 , 0.6 , 0.64 , 0.68 , 0.72 , 0.76 , 0.8 , 0.84 , 0.88 , 0.92 , 0.96 , 1, 1.2, 1.5
    };

    for (int i = 0; i < m_nDcasDca; ++i) {
        if (m_DcaEdgeDca[i+1]<m_DcaEdgeDca[i]) cout<<"wrong"<<endl;
    }

    TH3F* mh3VzZdcMult = (TH3F*)fDca1.Get("mh3VzZdcMult");
    outEvent->cd();
    mh3VzZdcMult->Write();

    TH1D* h1Vz[m_nmultEdge];
    TH1D* h1ZdcX[m_nmultEdge];

    int binVzmin = 1;
    int binVzup = mh3VzZdcMult->GetXaxis()->GetNbins(); //ok
    int binZDCmin = 1;
    int binZDCmax = mh3VzZdcMult->GetYaxis()->GetNbins(); //ok

    for (int ii = 0; ii < m_nmultEdge; ++ii)   {
        int binMultmin = mh3VzZdcMult->GetZaxis()->FindBin(m_multEdge[ii]);
        int binMultmax = mh3VzZdcMult->GetZaxis()->FindBin(m_multEdge[ii+1]);
        cout<<m_multEdge[ii]<<" "<<m_multEdge[ii+1]<<endl;
        cout<<binMultmin<<" "<<binMultmax<<endl;
        cout<<mh3VzZdcMult->GetZaxis()->GetNbins()<<endl;
        outEvent->cd();

        h1Vz[ii] = mh3VzZdcMult -> ProjectionX("_px",binZDCmin, binZDCmax, binMultmin, binMultmax, ""); //vz zdc
        h1Vz[ii]->SetDirectory(0);
        h1Vz[ii]->Write(Form("vz_mult_%i_%i", (int)m_multEdge[ii], (int)m_multEdge[ii+1]));
        h1ZdcX[ii] = mh3VzZdcMult -> ProjectionY("_py",binVzmin, binVzup, binMultmin, binMultmax, ""); //vz zdc
        h1ZdcX[ii]->SetDirectory(0);
        h1ZdcX[ii]->Scale(1/h1ZdcX[ii]->GetEntries());
        h1ZdcX[ii]->Write(Form("zdc_mult_%i_%i", (int)m_multEdge[ii], (int)m_multEdge[ii+1]));
    }

    TH1D* hrefMult = mh3VzZdcMult -> ProjectionZ("_pz",binVzmin, binVzup, binZDCmin, binZDCmax, "");
    hrefMult->Draw();
    outEvent->cd();
    hrefMult->Write("hrefMult");
    outEvent->Close();

    for(int iParticle = 0; iParticle < nParticles; ++iParticle){
        for (int iEta = 0; iEta < nEtas; ++iEta) {
            for (int iVz = 0; iVz < nVzs; ++iVz) {
                for (int iCent = 0; iCent < m_nmultEdge; ++iCent) {
                    TH3F *hist3D = 0;

                    for (int iZDC = 0; iZDC < nZDC; ++iZDC) {
                        TH3F *hist3D1 = 0;
//                        std::cout << iCent << " " << iEta << " " << iVz << std::endl;
                        // Getting 3D histogram
                        const char *h3dName = Form("mh3DcaXyZPt_p%d_eta%d_vz%d_z%d_m%d", iParticle, iEta, iVz, iZDC, iCent); //ok
//                        hist3D = (TH3F * )(fDca1.Get(h3dName));
//                        if (!hist3D) {
//                            std::cout << "histogram \"" << h3dName << "\" not found." << endl;
//                            return;
//                        }

                        if (iZDC == 0) {
                            hist3D = (TH3F * )(fDca1.Get(h3dName));
                            if (!hist3D) {
                                std::cout << "histogram \"" << h3dName << "\" not found." << endl;
                                return;
                            }
                        } else {
                            hist3D1 = (TH3F * )(fDca1.Get(h3dName));
                            if (!hist3D1) {
                                std::cout << "histogram \"" << h3dName << "\" not found." << endl;
                                return;
                            }
                            hist3D->Add(hist3D1);
                        }
                    }
                    // index of pT bin in 3D histogram
                    int iPt3d = 1;
                    int zdc=0;
                    // slicing 3D histogram into 2D slices
                    for (int iPt = 0; iPt < nPtBins; ++iPt) {
                        // creating the 2D histogram
                        const char *h2dName = Form("mh2DcaPtCentPartEtaVzPhi_p%i_eta%i_vz%i_m%i_pt%i_zdc%i", iParticle, iEta, iVz, iCent, iPt, zdc);
//                            const char *h2dName = Form("mh2DcaPtCentPartEtaVzPhi_p%i_eta%i_vz%i_pt%i_zdc%i", iParticle, iEta, iVz, iPt, iZDC);
                        double midPt = (ptEdge[iPt + 1] + ptEdge[iPt]) * 0.5; // middle of the pT bin
                        const char *h2dTitle = Form("%s_pT-%.2f", hist3D->GetTitle(), midPt);
                        TH2D* hist2D = new TH2D(h2dName, h2dTitle, m_nDcasDca, m_DcaEdgeDca, m_nDcasDca, m_DcaEdgeDca);
                        hist2D->GetXaxis()->SetTitle(hist3D->GetYaxis()->GetTitle());
                        hist2D->GetYaxis()->SetTitle(hist3D->GetZaxis()->GetTitle());

                        while ((hist3D->GetXaxis()->GetBinUpEdge(iPt3d)) - epsilon < ptEdge[iPt + 1]) {
                            // filling the 2D histogram
                            for (int iXy = 0; iXy < m_nDcasDca + 2; ++iXy) {
                                for (int iZ = 0; iZ < m_nDcasDca + 2; ++iZ) {
                                    // copy bin content
                                    Double_t last = hist2D->GetBinContent(iXy, iZ);
                                    hist2D->SetBinContent(iXy, iZ, (Double_t)hist3D->GetBinContent(iPt3d, iXy, iZ) + last);
                                }
                            }
                            // raise index of pT bin
                            ++iPt3d;
                        }

                        // testing if some particles were recorded
                        if ((hist2D->Integral() < 1000) && ptEdge[iPt]!=5) {
                            cout << "Number of particles recorded: " << hist2D->Integral()
                                 << " from pT histogram: " << iPt << " between " << ptEdge[iPt] << " and " << ptEdge[iPt + 1]
                                 << ", mult bin number: " << iCent << " between " << m_multEdge[iCent] << endl;
//                            outHist2dLowStats->cd();
//                            hist2D.Write();
                        }
                        // Saving the 2D histogram
                        outHistF->cd();
                        hist2D->Write(h2dName, 1 , 0);

                    }
                }
            }
        }
    }
//            cout<<"Finished writing histograms for centrality " << iCent << endl;

outHistF->Close();

cout << endl;
cout << "Done..." << endl;
}
