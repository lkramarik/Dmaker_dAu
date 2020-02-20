#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
//#include "StPicoD0EventMaker/StPicoD0Event.h"
//#include "StPicoD0EventMaker/StKaonPion.h"
#include "StPicoHFMaker/StHFCuts.h"
#include "phys_constants.h"
#include "StBTofUtil/tofPathLength.hh"
#include "StPicoSimInputsMaker.h"
#include "TH2F.h"
#include "TVector3.h"
#include "StAnaCuts.h"
#include <iostream>
using namespace std;

ClassImp(StPicoSimInputsMaker)

// _________________________________________________________
StPicoSimInputsMaker::StPicoSimInputsMaker(char const* name, StPicoDstMaker* picoMaker, char const* outputBaseFileName) :
        StPicoHFMaker(name, picoMaker, outputBaseFileName),
        mOutFileBaseName(outputBaseFileName){
    // constructor
}

// _________________________________________________________
StPicoSimInputsMaker::~StPicoSimInputsMaker() {
    // destructor
}

// _________________________________________________________
int StPicoSimInputsMaker::InitHF() {
    mOutFileBaseName = mOutFileBaseName.ReplaceAll(".root", "");
    unsigned int nTracks = mPicoDst->numberOfTracks();

    histoInit(mOutFileBaseName, true); //for createQA()
    return kStOK;
}

// _________________________________________________________
void StPicoSimInputsMaker::ClearHF(Option_t *opt="") {
    return;
}

// _________________________________________________________
int StPicoSimInputsMaker::FinishHF() {
    closeFile();
    return kStOK;
}
// _________________________________________________________
int StPicoSimInputsMaker::MakeHF() {
    createQA();
    return kStOK;
}

//_________________________________________________________
int StPicoSimInputsMaker::createQA(){
//    cout<<"createQA"<<endl;
    unsigned int nTracks = mPicoDst->numberOfTracks();
    int multiplicity = mPicoDst->event()->refMult();
    int ZdcIndex = getZdcIndex(mPicoDst->event()->ZDCx()/1000.);
    if (ZdcIndex==-1) return 0;
    mh3VzZdcMult -> Fill(mPrimVtx.z(),  mPicoDst->event()->ZDCx()/1000., multiplicity);

    for (unsigned int iTrack = 0; iTrack < nTracks; ++iTrack) {
        StPicoTrack const* trk = mPicoDst->track(iTrack);
        if (!trk) continue;
        StPicoPhysicalHelix helix = trk->helix(mBField);
        TVector3 momentum = trk->gMom(mPrimVtx, mBField);

        if (!(mHFCuts->isGoodTrack(trk))) continue; //nHitsFit, pT, max DCA, hft (must be non-HFT)
        if (!(fabs(momentum.PseudoRapidity()) < 1.0)) continue;

        bool goodPion = false;
        bool goodKaon = false;

        bool tofPion = false;
        bool tofKaon = false;

        bool tpcPion = false;
        bool tpcKaon = false;

        if(mHFCuts->isTPCPion(trk)) tpcPion = true;
        if(mHFCuts->isTPCKaon(trk)) tpcKaon = true;

        if(mHFCuts->isHybridTOFHadron(trk, mHFCuts->getTofBetaBase(trk), StPicoCutsBase::kPion)) tofPion = true;
        if(mHFCuts->isHybridTOFHadron(trk, mHFCuts->getTofBetaBase(trk), StPicoCutsBase::kKaon)) tofKaon = true;

        goodPion = (tofPion && tpcPion);
        goodKaon = (tofKaon && tpcKaon);

        TVector3 dcaPoint = trk->origin();
        float dca = (mPrimVtx - trk->origin()).Mag();
        float dcaZ = dcaPoint.z() - mPrimVtx.z();
        double dcaXy = helix.geometricSignedDistance(mPrimVtx.x(), mPrimVtx.y());

        int EtaIndex = getEtaIndexDca(fabs(momentum.PseudoRapidity()));
        int EtaIndexRatio = getEtaIndexRatio(momentum.PseudoRapidity());
//        cout<<"eta is "<<EtaIndex<<endl;
        int PhiIndex = getPhiIndexRatio(momentum.Phi());
        if ((EtaIndex==-1) || (PhiIndex==-1)) continue;

//        if (trk->isHFTTrack()) {
            for (int i = 0; i < 3; ++i) {
                if (trk->nSigmaPion() < i + 1) {
                    if (tofPion) { h1TofmatchTOF[0][i]->Fill(momentum.Perp()); }
                    h1Tofmatch[0][i]->Fill(momentum.Perp());
                }

                if (trk->nSigmaKaon() < i + 1) {
                    if (tofKaon) { h1TofmatchTOF[1][i]->Fill(momentum.Perp()); }
                    h1Tofmatch[1][i]->Fill(momentum.Perp());
                }
            }
//        }

        if (trk->isHFTTrack() && (goodPion || goodKaon) && vars::dcaHists){
            addDcaPtCent(dca, dcaXy, dcaZ, goodPion, goodKaon, momentum.Perp(), multiplicity, EtaIndex, PhiIndex, mPrimVtx.z(), ZdcIndex);
        }

        if (vars::ratioHists) {
            if ((goodPion || goodKaon)) {
                addTpcDenom1(goodPion, goodKaon, momentum.Perp(), multiplicity, EtaIndexRatio, PhiIndex, mPrimVtx.z(), mPicoDst->event()->ZDCx() / 1000.);
                mh2Tpc1PhiVz->Fill(momentum.Phi(), mPrimVtx.z());
            }

            if (trk->isHFTTrack() && (goodPion || goodKaon)) {
//        if (trk && fabs(dca) < 1.5 && trk->isHFTTrack() && (goodPion || goodKaon || goodProton) && fabs(dcaXy) < 1. && fabs(dcaZ) < 1.){
                addHFTNumer1(goodPion, goodKaon, momentum.Perp(), multiplicity, EtaIndexRatio, PhiIndex, mPrimVtx.z(), mPicoDst->event()->ZDCx() / 1000.);
                mh2HFT1PhiVz->Fill(momentum.Phi(), mPrimVtx.z());
            }
        }
    } // .. end tracks loop
    return 0;
}

// _________________________________________________________
void StPicoSimInputsMaker::histoInit(TString fileBaseName, bool fillQaHists) {
    mFillQaHists = fillQaHists;
    mOutFile = new TFile(fileBaseName + ".hists.root", "RECREATE");

    TH1::SetDefaultSumw2();
    if (!mFillQaHists) return;

    for (int iParticle = 0; iParticle < vars::m_nParticles; ++iParticle) {
        for (int nsigma = 0; nsigma < 3; ++nsigma) {
            h1Tofmatch[iParticle][nsigma] = new TH1D(Form("h1_Tofmatch_tpc1_p%d_nsigma%d", iParticle, nsigma+1), Form("h1_Tofmatch_tpc1_p%d_nsigma%d", iParticle, nsigma+1),vars::m_nPtTOF,vars::m_PtTOFedge);
            h1TofmatchTOF[iParticle][nsigma] = new TH1D(Form("h1_TofmatchTOF_tpc1_p%d_nsigma%d", iParticle, nsigma+1), Form("h1_TofmatchTOF_tpc1_p%d_nsigma%d", iParticle, nsigma+1),vars::m_nPtTOF,vars::m_PtTOFedge);
        }
    }

    if(vars::ratioHists) {
        mh2Tpc1PtCent = new TH2F("mh2Tpc1PtCent", "Tpc tracks;p_{T}(GeV/c);cent",  vars::m_nPtsRatio, vars::m_PtEdgeRatio, vars::m_nmultEdge, vars::m_multEdge);
        mh2HFT1PtCent = new TH2F("mh2HFT1PtCent", "HFT tracks;p_{T}(GeV/c);cent",  vars::m_nPtsRatio, vars::m_PtEdgeRatio, vars::m_nmultEdge, vars::m_multEdge);
        mh2Tpc1PhiVz = new TH2F("mh2Tpc1PhiVz", "Tpc tracks;#Phi;Vz", 100, -3.1415, 3.1415, 20, -10, 10);
        mh2HFT1PhiVz = new TH2F("mh2HFT1PhiVz", "HFT tracks;#Phi;Vz", 100, -3.1415, 3.1415, 20, -10, 10);

        for (int iParticle = 0; iParticle < vars::m_nParticles; iParticle++) {
            for (int iEta = 0; iEta < vars::m_nEtasRatio; iEta++) {
                for (int iVz = 0; iVz < vars::m_nVzsRatio; iVz++) {
                    for (int iPhi = 0; iPhi < vars::m_nPhisRatio; iPhi++) {
                        mh2Tpc1PtCentPartEtaVzPhi[iParticle][iEta][iVz][iPhi] = new TH3F(Form("h3_tpc_mult_pt_p%d_eta%d_vz%d_phi%d", iParticle, iEta, iVz, iPhi), Form("h3_tpc_mult_pt_p%d_eta%d_vz%d_phi%d", iParticle, iEta, iVz, iPhi), vars::m_nPtsRatio, vars::m_PtEdgeRatio, vars::m_nmultEdge,
                                                                                         vars::m_multEdge, vars::m_nZdc, vars::m_zdcEdge);
                        mh2HFT1PtCentPartEtaVzPhi[iParticle][iEta][iVz][iPhi] = new TH3F(Form("h3_hft_mult_pt_p%d_eta%d_vz%d_phi%d", iParticle, iEta, iVz, iPhi), Form("h3_hft_mult_pt_p%d_eta%d_vz%d_phi%d", iParticle, iEta, iVz, iPhi), vars::m_nPtsRatio, vars::m_PtEdgeRatio, vars::m_nmultEdge,
                                                                                         vars::m_multEdge, vars::m_nZdc, vars::m_zdcEdge);
                    }
                }
            }
        }
    }

    if (vars::dcaHists) {
        for (int iParticle = 0; iParticle < vars::m_nParticles; iParticle++) {
            for (int iEta = 0; iEta < vars::m_nEtasDca; iEta++) {
                for (int iVz = 0; iVz < vars::m_nVzsDca; iVz++) {
                    for (int iZdc = 0; iZdc < vars::m_nZdcDCA; iZdc++) {
                        for (int iCent = 0; iCent < vars::m_nmultEdge; iCent++) {
//                            cout<<Form("mh3DcaXyZPt_p%d_eta%d_vz%d_z%d_m%d", iParticle, iEta, iVz, iZdc, iCent)<<endl;
                            mh3DcaXyZPtCentPartEtaVzPhi[iParticle][iEta][iVz][iZdc][iCent] = new TH3F(Form("mh3DcaXyZPt_p%d_eta%d_vz%d_z%d_m%d", iParticle, iEta, iVz, iZdc, iCent), Form("mh3DcaXyZPt_p%d_eta%d_vz%d_z%d_m%d", iParticle, iEta, iVz, iZdc, iCent),
                                                                                                      vars::m_nPtsDca, vars::m_PtEdgeDca, vars::m_nDcasDca, vars::m_DcaEdgeDca, vars::m_nDcasDca, vars::m_DcaEdgeDca);
                        }
                    }
                }
            }
        }
        mh3DcaPtCent = new TH3F("mh3DcaPtCent", "mh3DcaPtCent;p_{T}(GeV/c);cent;Dca(cm)", vars::m_nPtsDca, vars::m_PtEdgeDca, vars::m_nmultEdge, vars::m_multEdge, vars::m_nDcasDca, vars::m_DcaEdgeDca);
        mh3DcaXyPtCent = new TH3F("mh3DcaXyPtCent", "mh3DcaXyPtCent;p_{T}(GeV/c);cent;DcaXy(cm)", vars::m_nPtsDca, vars::m_PtEdgeDca, vars::m_nmultEdge, vars::m_multEdge, vars::m_nDcasDca, vars::m_DcaEdgeDca);
        mh3DcaZPtCent = new TH3F("mh3DcaZPtCent", "mh3DcaZPtCent;p_{T}(GeV/c);cent;DcaZ(cm)", vars::m_nPtsDca, vars::m_PtEdgeDca, vars::m_nmultEdge, vars::m_multEdge, vars::m_nDcasDca, vars::m_DcaEdgeDca);

    }

    mh3VzZdcMult = new TH3F("mh3VzZdcMult", "mh3VzZdcMult", 100, -6, 6, 100, 0, 250,  50, 0, 50);
//    cout<<"endHistoint"<<endl;
}

// _________________________________________________________
void StPicoSimInputsMaker::addTpcDenom1(bool IsPion, bool IsKaon, float pt, int multiplicity, int EtaIndex, int PhiIndex, float Vz, float zdc){
    int VzIndex = getVzIndexRatio(Vz);
    if(VzIndex == -1) return;

    if (IsPion){
        mh2Tpc1PtCentPartEtaVzPhi[0][EtaIndex][VzIndex][PhiIndex]->Fill(pt, multiplicity, zdc);
    }

    if (IsKaon){
        mh2Tpc1PtCentPartEtaVzPhi[1][EtaIndex][VzIndex][PhiIndex]->Fill(pt, multiplicity, zdc);
    }

    mh2Tpc1PtCent->Fill(pt, multiplicity);
//    if (fabs(Eta) < mHFCuts->cutEta()  && pt > mHFCuts->cutPt()) mh2Tpc1PhiVz->Fill(Phi, Vz);
}

// _________________________________________________________
void StPicoSimInputsMaker::addHFTNumer1(bool IsPion, bool IsKaon, float pt, int multiplicity, int EtaIndex, int PhiIndex, float Vz, float zdc){
    int VzIndex = getVzIndexRatio(Vz);
    if(VzIndex == -1) return;
    if (IsPion){
        mh2HFT1PtCentPartEtaVzPhi[0][EtaIndex][VzIndex][PhiIndex]->Fill(pt, multiplicity, zdc);
    }
    if (IsKaon){
        mh2HFT1PtCentPartEtaVzPhi[1][EtaIndex][VzIndex][PhiIndex]->Fill(pt, multiplicity, zdc);
    }

    mh2HFT1PtCent->Fill(pt, multiplicity);
//    if (fabs(Eta) < mHFCuts->cutEta()  && pt > mHFCuts->cutPt()) mh2HFT1PhiVz->Fill(Phi, Vz);
}

// _________________________________________________________
void StPicoSimInputsMaker::addDcaPtCent(float dca, float dcaXy, float dcaZ, bool IsPion, bool IsKaon, float pt,  int multiplicity, int EtaIndex, int Phi, float Vz, int ZdcIndex){
    int VzIndex = getVzIndexDca(Vz);
    int centrality = getMultIndex(multiplicity);
    if(VzIndex == -1) return;
    if (centrality == -1) return;
    if (IsPion){
        mh3DcaXyZPtCentPartEtaVzPhi[0][EtaIndex][VzIndex][ZdcIndex][centrality]->Fill(pt, dcaXy, dcaZ);
    }
    if (IsKaon){
        mh3DcaXyZPtCentPartEtaVzPhi[1][EtaIndex][VzIndex][ZdcIndex][centrality]->Fill(pt, dcaXy, dcaZ);
    }

    mh3DcaPtCent->Fill(pt, multiplicity, dca);
    mh3DcaXyPtCent->Fill(pt, multiplicity, dcaXy);
    mh3DcaZPtCent->Fill(pt, multiplicity, dcaZ);
}

// _________________________________________________________
int StPicoSimInputsMaker::getZdcIndex(float zdc){
    for (int i = 0; i < vars::m_nZdcDCA; i++){
        if ((zdc >= vars::m_zdcEdgeDCA[i]) && (zdc < vars::m_zdcEdgeDCA[i + 1]))
            return i;
    }
    return -1;
}

// _________________________________________________________
int StPicoSimInputsMaker::getMultIndex(float multiplicity){
    for (int i = 0; i < vars::m_nmultEdge; i++){
        if ((multiplicity >= vars::m_multEdge[i]) && (multiplicity < vars::m_multEdge[i + 1]))
            return i;
    }
    return -1;
}

// _________________________________________________________
int StPicoSimInputsMaker::getEtaIndexDca(float Eta){
    for (int i = 0; i < vars::m_nEtasDca; i++){
        if ((Eta >= vars::m_EtaEdgeDca[i]) && (Eta < vars::m_EtaEdgeDca[i + 1]))
            return i;
    }
    return -1;
}

// _________________________________________________________
int StPicoSimInputsMaker::getVzIndexDca(float Vz){
    for (int i = 0; i < vars::m_nVzsDca; i++){
        if ((Vz >= vars::m_VzEdgeDca[i]) && (Vz < vars::m_VzEdgeDca[i + 1]))
            return i;
    }
    return -1;
}

// _________________________________________________________
int StPicoSimInputsMaker::getEtaIndexRatio(float Eta){
    for (int i = 0; i < vars::m_nEtasRatio; i++){
        if ((Eta >= vars::m_EtaEdgeRatio[i]) && (Eta < vars::m_EtaEdgeRatio[i + 1]))
            return i;
    }
    return -1;
}

// _________________________________________________________
int StPicoSimInputsMaker::getPhiIndexRatio(float Phi){
    for (int i = 0; i < vars::m_nPhisRatio; i++){
        if ((Phi >= vars::m_PhiEdgeRatio[i]) && (Phi < vars::m_PhiEdgeRatio[i + 1]))
            return i;
    }
    return -1;
}

// _________________________________________________________
int StPicoSimInputsMaker::getVzIndexRatio(float Vz){
    for (int i = 0; i < vars::m_nVzsRatio; i++) {
        if ((Vz >= vars::m_VzEdgeRatio[i]) && (Vz < vars::m_VzEdgeRatio[i + 1]))
            return i;
    }
    return -1;
}

// _________________________________________________________
void StPicoSimInputsMaker::closeFile()
{
    mOutFile->cd();

//    mh1Cent->Write();
//    mh1CentWg->Write();
//    mh1gRefmultCor->Write();
//    mh1gRefmultCorWg->Write();
//    mh2CentVz->Write();
//    mh2CentVzWg->Write();

    //HFT DCA Ratio
    if (vars::dcaHists) {
        for (int iParticle = 0; iParticle < vars::m_nParticles; iParticle++) {
            for (int iEta = 0; iEta < vars::m_nEtasDca; iEta++) {
                for (int iVz = 0; iVz < vars::m_nVzsDca; iVz++) {
                    for (int iZdc = 0; iZdc < vars::m_nZdcDCA; iZdc++) {
                        for (int iCent = 0; iCent < vars::m_nmultEdge; iCent++) {
                            mh3DcaXyZPtCentPartEtaVzPhi[iParticle][iEta][iVz][iZdc][iCent]->Write();
                        }
                    }
                }
            }
        }
        mh3DcaZPtCent->Write();
        mh3DcaPtCent->Write();
        mh3DcaXyPtCent->Write();
    }

    if (vars::ratioHists) {
        mh2Tpc1PtCent->Write();
        mh2Tpc1PhiVz->Write();
        mh2HFT1PhiVz->Write();
        mh2HFT1PtCent->Write();
        for (int iParticle = 0; iParticle < vars::m_nParticles; iParticle++) {
            for (int iEta = 0; iEta < vars::m_nEtasRatio; iEta++) {
                for (int iVz = 0; iVz < vars::m_nVzsRatio; iVz++) {
                    for (int iPhi = 0; iPhi < vars::m_nPhisRatio; iPhi++) {
                        mh2Tpc1PtCentPartEtaVzPhi[iParticle][iEta][iVz][iPhi]->Write();
                        mh2HFT1PtCentPartEtaVzPhi[iParticle][iEta][iVz][iPhi]->Write();
                    }
                }
            }
        }
    }

    for (int iParticle = 0; iParticle < vars::m_nParticles; ++iParticle) {
        for (int nsigma = 0; nsigma < 3; ++nsigma) {
            h1Tofmatch[iParticle][nsigma] -> Write();
            h1TofmatchTOF[iParticle][nsigma] -> Write();
        }
    }

    mh3VzZdcMult -> Write();
//    mOutFile->Write();
    mOutFile->Close();

    //mOutFile->Delete();
}

// _________________________________________________________
double StPicoSimInputsMaker::DCA(StPicoTrack const * const trk, TVector3 const & vtx) const {
    return ((trk->origin() - vtx).Mag());
}
