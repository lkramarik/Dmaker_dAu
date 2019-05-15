//
// Created by lukas on 9.4.2018.
//

#include "TH3D.h"
#include "TFile.h"
#include "TH2D.h"
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
    TFile fDca1("dca.hists.root");
    TFile *outHistF = new TFile("dcaxy_vs_dcaz.root", "RECREATE");
    outHistF->SetCompressionSettings(0);

//    const int nZDC = 5;
//    float const m_zdcEdge[nZDC+1] = {0,50,90,130,170,210};

    const int nZDC = 3;
    float const m_zdcEdge[nZDC+1] = {0,90,170,210};

//    const int m_nmultEdge = 7;
//    float const m_multEdge[m_nmultEdge+1] = {0, 4, 8, 12, 16, 20, 24, 200};

    const int m_nmultEdge = 1; //7
    float const m_multEdge[m_nmultEdge+1] = {0, 200}; //currently not used in dca


    const int nPtBins = 12;
    float const ptEdge[nPtBins + 1] = {0.2,0.3,0.4,0.5,0.6,0.8,1.,1.25,1.5,2.,3.0,5.,12.};

    const int m_nDcasDca = 148;
    float const  m_DcaEdgeDca[m_nDcasDca + 1] =   {
            -1.5, -1.2, -1 , -0.96 , -0.92 , -0.88 , -0.84 , -0.8 , -0.76 , -0.72 , -0.68 , -0.64 , -0.6 , -0.56 , -0.52 , -0.48 , -0.44 , -0.4 , -0.36 , -0.32 , -0.28 , -0.24 , -0.2 , -0.16 , -0.12 ,  -0.08,
            -0.078 , -0.075 , -0.072 , -0.069 , -0.066 , -0.063 , -0.06 , -0.057 , -0.054 , -0.051 , -0.048 , -0.045 , -0.042 , -0.039 , -0.036 , -0.033 , -0.03 , -0.027 , -0.024 , -0.021 , -0.018 , -0.015 , -0.012 ,
            -0.01 , -0.0096 , -0.0092 , -0.0088 , -0.0084 , -0.008 , -0.0076 , -0.0072 , -0.0068 , -0.0064 , -0.006 , -0.0056 , -0.0052 , -0.0048 , -0.0044 , -0.004 , -0.0036 , -0.0032 , -0.0028 , -0.0024 , -0.002 , -0.0016 , -0.0012 , -0.0008 , -0.0004 , 0 , 0.0004 , 0.0008 , 0.0012 , 0.0016 , 0.002 , 0.0024 , 0.0028 , 0.0032 , 0.0036 , 0.004 , 0.0044 , 0.0048 , 0.0052 , 0.0056 , 0.006 , 0.0064 , 0.0068 , 0.0072 , 0.0076 , 0.008 , 0.0084 , 0.0088 , 0.0092 , 0.0096 , 0.01 ,
            0.012 , 0.015 , 0.018 , 0.021 , 0.024 , 0.027 , 0.03 , 0.033 , 0.036 , 0.039 , 0.042 , 0.045 , 0.048 , 0.051 , 0.054 , 0.057 , 0.06 , 0.063 , 0.066 , 0.069 , 0.072 , 0.075 , 0.078 ,
            0.08 , 0.12 , 0.16 , 0.2 , 0.24 , 0.28 , 0.32 , 0.36 , 0.4 , 0.44 , 0.48 , 0.52 , 0.56 , 0.6 , 0.64 , 0.68 , 0.72 , 0.76 , 0.8 , 0.84 , 0.88 , 0.92 , 0.96 , 1, 1.2, 1.5
    };


    for(int iParticle = 0; iParticle < nParticles; ++iParticle){
            for (int iEta = 0; iEta < nEtas; ++iEta) {
                for (int iVz = 0; iVz < nVzs; ++iVz) {
                    for (int iZDC = 0; iZDC < nZDC; ++iZDC) {

                        TH3F *hist3D = 0;
                        TH3F *hist3D1 = 0;

                        for (int iCent = 0; iCent < m_nmultEdge; ++iCent) {

                            std::cout << iCent << " " << iEta << " " << iVz << std::endl;
                            // Getting 3D histogram

                            const char *h3dName = Form("mh3DcaXyZPt_p%d_eta%d_vz%d_z%d_m%d", iParticle, iEta, iVz, iZDC, iCent); //ok
                            if (iCent==0) {
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


                        // getting parameters of the 3D histogram (nBins, min and max of the axes, and title)
                        int nXyBins = m_nDcasDca;
                        const char *xyTitle = hist3D->GetYaxis()->GetTitle();
                        int nZBins = m_nDcasDca;
                        const char *zTitle = hist3D->GetZaxis()->GetTitle();
                        const char *hist3dTitle = hist3D->GetTitle();

                        // index of pT bin in 3D histogram
                        int iPt3d = 1;

                        // slicing 3D histogram into 2D slices
                        for (int iPt = 0; iPt < nPtBins; ++iPt) {
                            // creating the 2D histogram
//                            const char *h2dName = Form("mh2DcaPtCentPartEtaVzPhi_p%i_eta%i_vz%i_m%i_pt%i_zdc%i", iParticle, iEta, iVz, iCent, iPt, iZDC);
                            const char *h2dName = Form("mh2DcaPtCentPartEtaVzPhi_p%i_eta%i_vz%i_pt%i_zdc%i", iParticle, iEta, iVz, iPt, iZDC);
                            double midPt = (ptEdge[iPt + 1] + ptEdge[iPt]) * 0.5; // middle of the pT bin
                            const char *h2dTitle = Form("%s_pT-%.2f", hist3dTitle, midPt);
                            TH2F hist2D(h2dName, h2dTitle, nXyBins, m_DcaEdgeDca, nZBins, m_DcaEdgeDca);
//                            TH2D hist2D(h2dName, h2dTitle, nXyBins, m_DcaEdgeDca, nZBins, m_DcaEdgeDca);
                            hist2D.GetXaxis()->SetTitle(xyTitle);
                            hist2D.GetYaxis()->SetTitle(zTitle);

                            while ((hist3D->GetXaxis()->GetBinUpEdge(iPt3d)) - epsilon < ptEdge[iPt + 1]) {
                                // filling the 2D histogram
                                for (int iXy = 0; iXy < nXyBins + 2; ++iXy) {
                                    for (int iZ = 0; iZ < nZBins + 2; ++iZ) {
                                        // copy bin content
                                        double last = hist2D.GetBinContent(iXy, iZ);
                                        hist2D.SetBinContent(iXy, iZ, hist3D->GetBinContent(iPt3d, iXy, iZ) + last);
                                    }
                                }
                                // raise index of pT bin
                                ++iPt3d;
                            }



                            // testing if some particles were recorded
                            if (hist2D.Integral() == 0)
                                cout << "Number of particles recorded: " << hist2D.Integral()
                                     << " from pT histogram: " << iPt << " between " << ptEdge[iPt] << " and " << ptEdge[iPt + 1]
                                     << ", pT bin number: " << iPt3d << " between " << hist3D->GetXaxis()->GetBinLowEdge(iPt3d - 1) << " and " << hist3D->GetXaxis()->GetBinUpEdge(iPt3d - 1) << endl;
                            // Saving the 2D histogram
                            outHistF->cd();
                            hist2D.Write();
                        }
                    }
                }
            }
//            cout<<"Finished writing histograms for centrality " << iCent << endl;
        }
    outHistF->Close();

    cout << endl;
    cout << "Done..." << endl;
}
