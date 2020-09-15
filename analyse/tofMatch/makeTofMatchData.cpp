#include "TH3D.h"
#include "TH3F.h"
#include "TFile.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TString.h"
#include "TSystem.h"
#include "TLegend.h"


#include <iostream>
using namespace std;

void setHisto(TH1D* h){
    h->SetStats(0);
    h->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h->SetTitle("");
    h->GetXaxis()->SetRangeUser(0.,5.);
    h->GetYaxis()->SetRangeUser(0.,1.);

}


void makeTofMatchData(TString inputFileName="tofratio.all.1509.root") {
    const int m_nParticlesCharged = 4;
    const TString m_ParticleChargedName[m_nParticlesCharged] = {"PionPlus", "PionMinus", "KaonPlus", "KaonMinus"};

    TFile* inputFile = new TFile(inputFileName, "READ");

    TH1D* h1Tofmatch[m_nParticlesCharged][3];
    TH1D* h1TofmatchTOF[m_nParticlesCharged][3];

    TString hisName;
    TCanvas* cNsigma=new TCanvas("cNsigma","cNsigma", 1200,1200);
    cNsigma->Divide(2,2);
    //TCanvas* cNsigma=new TCanvas("cNsigma","c", 1200,1000);

    for (int iParticle = 0; iParticle < m_nParticlesCharged; ++iParticle) {
        for (int nsigma = 0; nsigma < 3; ++nsigma) {
            hisName=Form("h1_Tofmatch_tpc_%s_nsigma%d", m_ParticleChargedName[iParticle].Data(), nsigma+1);
            h1Tofmatch[iParticle][nsigma] = (TH1D*)inputFile->Get(hisName);
            hisName=Form("h1_Tofmatch_tof_%s_nsigma%d", m_ParticleChargedName[iParticle].Data(), nsigma+1);
            h1TofmatchTOF[iParticle][nsigma] = (TH1D*)inputFile->Get(hisName);



            h1TofmatchTOF[iParticle][nsigma]->Divide(h1Tofmatch[iParticle][nsigma]);
            setHisto(h1Tofmatch[iParticle][nsigma]);
            setHisto(h1TofmatchTOF[iParticle][nsigma]);
            cNsigma->cd(iParticle+1);
            if (nsigma==0) h1TofmatchTOF[iParticle][nsigma]->Draw();
            else h1TofmatchTOF[iParticle][nsigma]->Draw("same");
        }
    }

//    for (int i = 0; i < ; ++i) {
//
//    }



}