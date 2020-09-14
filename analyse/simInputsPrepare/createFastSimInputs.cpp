#include "TH3D.h"
#include "TH3F.h"
#include "TFile.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TString.h"
#include "TSystem.h"
#include "TLegend.h"

#include "StAnaCutsData.h"

#include <iostream>
using namespace std;


//-----------------------------------------
void createFastSimInputs(TString inputFileName="hftratio_all.root") {
    inputFileName="produced/"+inputFileName;
    TFile* inputFile = new TFile(inputFileName, "READ");

    TString shortName = inputFileName.ReplaceAll(".root", 5, "/", 1);
    cout<<shortName<<endl;
    gSystem->Exec("mkdir -p "+shortName);


    TH1F* mh2Tpc1PtCentPartEtaVzPhi[vars::m_nParticles][vars::m_nEtasRatio][vars::m_nVzsRatio][vars::m_nPhisRatio][vars::m_nmultEdge];
    TH1F* mh2HFT1PtCentPartEtaVzPhi[vars::m_nParticles][vars::m_nEtasRatio][vars::m_nVzsRatio][vars::m_nPhisRatio][vars::m_nmultEdge];

    TString hisName;
    //getting pT histos
    for (int iParticle = 0; iParticle < vars::m_nParticles; iParticle++) {
        for (int iEta = 0; iEta < vars::m_nEtasRatio; iEta++) {
            for (int iVz = 0; iVz < vars::m_nVzsRatio; iVz++) {
                for (int iPhi = 0; iPhi < vars::m_nPhisRatio; iPhi++) {
                    for (int iMult = 0; iMult < vars::m_nmultEdge; ++iMult) {
                        hisName = Form("tpc/h_tpc_pt_p%d_eta%d_vz%d_phi%d_m%d", iParticle, iEta, iVz, iPhi, iMult);
                        mh2Tpc1PtCentPartEtaVzPhi[iParticle][iEta][iVz][iPhi][iMult] = (TH1F * )(inputFile->Get(hisName));

                        hisName = Form("hft/h_hft_pt_p%d_eta%d_vz%d_phi%d_m%d", iParticle, iEta, iVz, iPhi, iMult);
                        mh2HFT1PtCentPartEtaVzPhi[iParticle][iEta][iVz][iPhi][iMult] = (TH1F * )(inputFile->Get(hisName));
                    }
                }
            }
        }
    }

    TFile* mOutFileRatio;
    TString fileNames[]={"hftratio_vs_pt_dAu_pion_hijing.root", "hftratio_vs_pt_dAu_kaon_hijing.root"};

    for (int iParticle = 0; iParticle < vars::m_nParticles; iParticle++) {
        mOutFileRatio = new TFile(shortName+fileNames[iParticle], "RECREATE");
        mOutFileRatio->cd();

        for (int iEta = 0; iEta < vars::m_nEtasRatio; iEta++) {
            for (int iVz = 0; iVz < vars::m_nVzsRatio; iVz++) {
                for (int iPhi = 0; iPhi < vars::m_nPhisRatio; iPhi++) {
                    for (int iMult = 0; iMult < vars::m_nmultEdge; ++iMult) {
                        mh2HFT1PtCentPartEtaVzPhi[iParticle][iEta][iVz][iPhi][iMult]->Divide(mh2Tpc1PtCentPartEtaVzPhi[iParticle][iEta][iVz][iPhi][iMult]);
                        hisName = Form("h_hftratio_p%d_eta%d_vz%d_phi%d_m%d", iParticle, iEta, iVz, iPhi, iMult);
                        mh2HFT1PtCentPartEtaVzPhi[iParticle][iEta][iVz][iPhi][iMult]->Clone()->Write(hisName, TObject::kOverwrite);
                    }
                }
            }
        }
        mOutFileRatio->Close();
    }

    TH2F* h2PtCent[] = {(TH2F*)inputFile->Get("mh2HFT1PtCent"), (TH2F*)inputFile->Get("mh2Tpc1PtCent")};
    TString tmpName[] = {"hft", "tpc"};
    TH1D* ptCent[2][vars::m_nmultEdge];// [hft/tpc][cent]

    Int_t binMultMin, binMultMax;
    for (int i = 0; i < 2; ++i) {
        for (int iMult = 0; iMult < vars::m_nmultEdge; ++iMult) {
            binMultMin=h2PtCent[i]->GetYaxis()->FindBin(vars::m_multEdge[iMult]);
            binMultMax=h2PtCent[i]->GetYaxis()->FindBin(vars::m_multEdge[iMult+1]);
            hisName=Form("hft_ratio_%.1f_%.1f_%s", vars::m_multEdge[iMult], vars::m_multEdge[iMult+1], tmpName[i].Data());

            ptCent[i][iMult]=h2PtCent[i]->ProjectionX(hisName,binMultMin, binMultMax, "");
        }
    }

    TLegend *legend = new TLegend(0.6, 0.71, 0.75, 0.89);
    legend -> SetFillStyle(0);
    legend -> SetLineColor(0);
    legend -> SetTextSize(0.035);

    TCanvas* c=new TCanvas("c","c", 1200,1000);
    gStyle->SetPalette(kPastel);
    for (int iMult = 0; iMult < vars::m_nmultEdge; ++iMult) {
        ptCent[0][iMult]->Divide(ptCent[1][iMult]);
        hisName=Form("hft_ratio_%.1f_%.1f", vars::m_multEdge[iMult], vars::m_multEdge[iMult+1]);
        ptCent[0][iMult]->SetName(hisName);
        ptCent[0][iMult]->SetTitle("");
        ptCent[0][iMult]->GetYaxis()->SetRangeUser(0,1);
        ptCent[0][iMult]->GetYaxis()->SetTitle("HFT ratio");
        ptCent[0][iMult]->SetStats(0);
        ptCent[0][iMult]->SetLineWidth(2);
//
        if (iMult==0) ptCent[0][iMult]->Draw("hist PLC PMC");
        else ptCent[0][iMult]->Draw("hist same PLC PMC");

        hisName=Form("%.0f < refMult < %.0f", vars::m_multEdge[iMult], vars::m_multEdge[iMult+1]);
        legend -> AddEntry(ptCent[0][iMult], hisName, "pl");

    }

    legend->Draw("same");
    c->SaveAs(inputFileName+"hft_plot.png");
    c->SaveAs(inputFileName+"hft_plot.pdf");
    c->Close();

    inputFile->Close();
}
