#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoDstMaker/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoHFMaker/StHFCuts.h"
#include "phys_constants.h"
#include "StBTofUtil/tofPathLength.hh"
#include "StPicoD0AnaMaker.h"
ClassImp(StPicoD0AnaMaker)

// _________________________________________________________
StPicoD0AnaMaker::StPicoD0AnaMaker(char const* name, StPicoDstMaker* picoMaker, char const* outputBaseFileName, char const* inputHFListHFtree = "") :
        StPicoHFMaker(name, picoMaker, outputBaseFileName, inputHFListHFtree),
        mOutFileBaseName(outputBaseFileName){
    // constructor
}

// _________________________________________________________
StPicoD0AnaMaker::~StPicoD0AnaMaker() {
    // destructor
}

// _________________________________________________________
int StPicoD0AnaMaker::InitHF() {
    // -- INITIALIZE USER HISTOGRAMS ETC HERE -------------------
    //    add them to the output list mOutList which is automatically written
    //
    // EXAMPLE //  mOutList->Add(new TH1F(...));
    // EXAMPLE //  TH1F* hist = static_cast<TH1F*>(mOutList->Last());
//    mOutList->Add(new TH2F("h_piTOF","h_piTOF",100,0,10, 250, -1, 1.5));
//    mOutList->Add(new TH2F("h_kTOF","h_kTOF",100,0,10, 250, -1, 1.5));
//    mOutList->Add(new TH2F("h_pTOF","h_pTOF",100,0,10, 250, -1, 1.5));
//
//    mOutList->Add(new TH2F("h_piTOF_20","h_piTOF_20",100,0,10, 300, 0, 1));
//    mOutList->Add(new TH2F("h_kTOF_20","h_kTOF_20",100,0,10, 300, 0, 1));
//    mOutList->Add(new TH2F("h_pTOF_20","h_pTOF_20",100,0,10, 300, 0, 1));
//
//    mOutList->Add(new TH2F("h_piTOF_HFT","h_piTOF_HFT",100,0,10, 300, 0, 1));
//    mOutList->Add(new TH2F("h_kTOF_HFT","h_kTOF_HFT",100,0,10, 300, 0, 1));
//    mOutList->Add(new TH2F("h_pTOF_HFT","h_pTOF_HFT",100,0,10, 300, 0, 1));
//
//    mOutList->Add(new TH2F("h_piTOF_HFT_20","h_piTOF_HFT_20",100,0,10, 300, 0, 1));
//    mOutList->Add(new TH2F("h_kTOF_HFT_20","h_kTOF_HFT_20",100,0,10, 300, 0, 1));
//    mOutList->Add(new TH2F("h_pTOF_HFT_20","h_pTOF_HFT_20",100,0,10, 300, 0, 1));
//
//    mOutList->Add(new TH2F("h_piTOFbeta","h_piTOFbeta",500,0,10, 300, 0, 1));
//    mOutList->Add(new TH2F("h_kTOFbeta","h_kTOFbeta",500,0,10, 300, 0, 1));
//    mOutList->Add(new TH2F("h_pTOFbeta","h_pTOFbeta",500,0,10, 300, 0, 1));
//
//    mOutList->Add(new TH2F("h_pinsigma","h_pinsigma",1000,0,10, 99, -5, 5));
//    mOutList->Add(new TH2F("h_knsigma","h_knsigma",1000,0,10, 99, -5, 5));
//    mOutList->Add(new TH2F("h_pnsigma","h_pnsigma",1000,0,10, 99, -5, 5));
//
//    mOutList->Add(new TH2F("h_dedx","h_dedx", 1000, 0, 10, 1000, 0, 10));
//h_tracktest
	DeclareHistograms();

    mOutList->Add(new TH1D("h_tracktest","h_tracktest", 6, 0.5, 6.5));
    mOutList->Add(new TH1D("h_tracktest_TOF","h_tracktest_TOF", 6, 0.5, 6.5));

    mOutFileBaseName = mOutFileBaseName.ReplaceAll(".root", "");

//    ntp_kaon = new TNtuple("ntp_kaon", "kaon tree","k_pt:k_phi:k_eta:k_nSigma:k_nHitFit:k_TOFinvbeta:pi_eventId:pi_runId");
//    ntp_pion = new TNtuple("ntp_pion", "pion tree","pi_pt:pi_phi:pi_eta:pi_nSigma:pi_nHitFit:pi_TOFinvbeta:k_eventId:k_runId");
    ntp_DMeson_Signal = new TNtuple("ntp_signal","DMeson TreeSignal",            "grefMult:runId:eventId:pi1_pt:pi1_dca:pi1_nSigma:pi1_nHitFit:pi1_TOFinvbeta:pi1_betaBase:k_pt:k_dca:k_nSigma:k_nHitFit:k_TOFinvbeta:k_betaBase:dcaDaughters:flag:primVz:primVzVpd:primVzDiff:D_theta:cosTheta:D_decayL:dcaD0ToPv:D_cosThetaStar:D_pt:D_mass");
    ntp_DMeson_Background = new TNtuple("ntp_background","DMeson TreeBackground","grefMult:runId:eventId:pi1_pt:pi1_dca:pi1_nSigma:pi1_nHitFit:pi1_TOFinvbeta:pi1_betaBase:k_pt:k_dca:k_nSigma:k_nHitFit:k_TOFinvbeta:k_betaBase:dcaDaughters:flag:primVz:primVzVpd:primVzDiff:D_theta:cosTheta:D_decayL:dcaD0ToPv:D_cosThetaStar:D_pt:D_mass");

    return kStOK;
}

// _________________________________________________________
void StPicoD0AnaMaker::ClearHF(Option_t *opt="") {
    return;
}

// _________________________________________________________
int StPicoD0AnaMaker::FinishHF() {
    ntp_DMeson_Signal -> Write(ntp_DMeson_Signal->GetName(), TObject::kOverwrite);
    ntp_DMeson_Background -> Write(ntp_DMeson_Background->GetName(), TObject::kOverwrite);
    WriteHistograms();
//    ntp_pion -> Write(ntp_pion->GetName(), TObject::kOverwrite);
//    ntp_kaon -> Write(ntp_kaon->GetName(), TObject::kOverwrite);
    return kStOK;
}
// _________________________________________________________
int StPicoD0AnaMaker::MakeHF() {
	getHadronCorV2(1);
    createCandidates();
//    analyzeCandidates();

//    TH2F *h_piTOF = static_cast<TH2F*>(mOutList->FindObject("h_piTOF"));
//    TH2F *h_kTOF = static_cast<TH2F*>(mOutList->FindObject("h_kTOF"));
//    TH2F *h_pTOF = static_cast<TH2F*>(mOutList->FindObject("h_pTOF"));
//
//    TH2F *h_piTOF_20 = static_cast<TH2F*>(mOutList->FindObject("h_piTOF_20"));
//    TH2F *h_kTOF_20 = static_cast<TH2F*>(mOutList->FindObject("h_kTOF_20"));
//    TH2F *h_pTOF_20 = static_cast<TH2F*>(mOutList->FindObject("h_pTOF_20"));
//
//    TH2F *h_piTOF_HFT = static_cast<TH2F*>(mOutList->FindObject("h_piTOF_HFT"));
//    TH2F *h_kTOF_HFT = static_cast<TH2F*>(mOutList->FindObject("h_kTOF_HFT"));
//    TH2F *h_pTOF_HFT = static_cast<TH2F*>(mOutList->FindObject("h_pTOF_HFT"));
//
//    TH2F *h_piTOF_HFT_20 = static_cast<TH2F*>(mOutList->FindObject("h_piTOF_HFT_20"));
//    TH2F *h_kTOF_HFT_20 = static_cast<TH2F*>(mOutList->FindObject("h_kTOF_HFT_20"));
//    TH2F *h_pTOF_HFT_20 = static_cast<TH2F*>(mOutList->FindObject("h_pTOF_HFT_20"));
//
//    TH2F *h_piTOFbeta = static_cast<TH2F*>(mOutList->FindObject("h_piTOFbeta"));
//    TH2F *h_kTOFbeta = static_cast<TH2F*>(mOutList->FindObject("h_kTOFbeta"));
//    TH2F *h_pTOFbeta = static_cast<TH2F*>(mOutList->FindObject("h_pTOFbeta"));
//
//    TH2F *h_pinsigma = static_cast<TH2F*>(mOutList->FindObject("h_pinsigma"));
//    TH2F *h_knsigma = static_cast<TH2F*>(mOutList->FindObject("h_knsigma"));
//    TH2F *h_pnsigma = static_cast<TH2F*>(mOutList->FindObject("h_pnsigma"));
//
//    TH2F *h_dedx = static_cast<TH2F*>(mOutList->FindObject("h_dedx"));
//
//    StThreeVectorF pVtx = mPicoDst->event()->primaryVertex();
//
//    UInt_t nTracks = mPicoDst->numberOfTracks();
//    for (unsigned short iTrack = 0; iTrack < nTracks; ++iTrack){
//        StPicoTrack const* trk = mPicoDst->track(iTrack);
//        if (!trk) continue;
//        StPhysicalHelixD helix = trk->helix(mBField);
//        StThreeVectorF momentum = trk->gMom(pVtx, mPicoDst->event()->bField());
//
//        if (!(trk->nHitsFit()>=15)) continue;
//        if (!(fabs(momentum.pseudoRapidity()) <= 1.0)) continue;
//
//        if (fabs(trk->nSigmaPion())<3.0){
//            if (mHFCuts->isTOFHadron(trk, mHFCuts->getTofBetaBase(trk), StHFCuts::kPion)) {
//                    h_piTOF->Fill(trk->gPt(),mHFCuts->getTofBetaBase(trk));
//                    float oneOverBeta = getOneOverBeta(trk, mHFCuts->getTofBetaBase(trk), StHFCuts::kPion);
//                    h_piTOFbeta->Fill(trk->gPt(),oneOverBeta);
//                    if (trk->nHitsFit()>=20) h_piTOF_20->Fill(trk->gPt(),oneOverBeta);
//                    if (trk->isHFTTrack()) h_piTOF_HFT->Fill(trk->gPt(),oneOverBeta);
//                    if ((trk->isHFTTrack()) && (trk->nHitsFit()>=20)) h_piTOF_HFT_20->Fill(trk->gPt(),oneOverBeta);
//            }
//        }
//
//        if (fabs(trk->nSigmaKaon())<3.0){
//            if (mHFCuts->isTOFHadron(trk, mHFCuts->getTofBetaBase(trk), StHFCuts::kKaon)) {
//                h_kTOF->Fill(trk->gPt(),mHFCuts->getTofBetaBase(trk));
//                float oneOverBeta = getOneOverBeta(trk, mHFCuts->getTofBetaBase(trk), StHFCuts::kKaon);
//                h_kTOFbeta->Fill(trk->gPt(),oneOverBeta);
//                if (trk->nHitsFit()>=20) h_kTOF_20->Fill(trk->gPt(),oneOverBeta);
//                if (trk->isHFTTrack()) h_kTOF_HFT->Fill(trk->gPt(),oneOverBeta);
//                if ((trk->isHFTTrack()) && (trk->nHitsFit()>=20)) h_kTOF_HFT_20->Fill(trk->gPt(),oneOverBeta);
//            }
//        }
//
//        if (fabs(trk->nSigmaProton())<3.0){
//            if (mHFCuts->isTOFHadron(trk, mHFCuts->getTofBetaBase(trk), StHFCuts::kProton)) {
//                h_pTOF->Fill(trk->gPt(),mHFCuts->getTofBetaBase(trk));
//                float oneOverBeta = getOneOverBeta(trk, mHFCuts->getTofBetaBase(trk), StHFCuts::kProton);
//                h_pTOFbeta->Fill(trk->gPt(),oneOverBeta);
//                if (trk->nHitsFit()>=20) h_pTOF_20->Fill(trk->gPt(),oneOverBeta);
//                if (trk->isHFTTrack()) h_pTOF_HFT->Fill(trk->gPt(),oneOverBeta);
//                if ((trk->isHFTTrack()) && (trk->nHitsFit()>=20)) h_pTOF_HFT_20->Fill(trk->gPt(),oneOverBeta);
//            }
//        }
//
//        h_pinsigma->Fill(momentum.mag(),trk->nSigmaPion());
//        h_knsigma->Fill(momentum.mag(),trk->nSigmaKaon());
//        h_pnsigma->Fill(momentum.mag(),trk->nSigmaProton());
//        h_dedx->Fill(momentum.mag(),trk->dEdx());
//
//    } // .. end tracks loop

    return kStOK;
}

// _________________________________________________________
int StPicoD0AnaMaker::createCandidates() {
//    for (unsigned short idxPion1 = 0; idxPion1 < mIdxPicoPions.size(); ++idxPion1) {
//        StPicoTrack const *pion1 = mPicoDst->track(mIdxPicoPions[idxPion1]);
//        for (unsigned short idxKaon = 0; idxKaon < mIdxPicoKaons.size(); ++idxKaon) {
//            StPicoTrack const *kaon = mPicoDst->track(mIdxPicoKaons[idxKaon]);
//    StHFPair *pair = new StHFPair(pion1, kaon, mHFCuts->getHypotheticalMass(StPicoCutsBase::kPion),mHFCuts->getHypotheticalMass(StPicoCutsBase::kKaon), mIdxPicoPions[idxPion1],mIdxPicoKaons[idxKaon], mPrimVtx, mBField, kTRUE);

    for(unsigned int i=0;i<mPicoDst->numberOfTracks();i++)  {
        StPicoTrack const* pion1 = mPicoDst->track(i);
        if (!mHFCuts -> isGoodPion(pion1)) continue;

        for(unsigned  int j=0;j<mPicoDst->numberOfTracks();j++)  {
            StPicoTrack const* kaon = mPicoDst->track(j);
            if (pion1->id() == kaon->id()) continue;

            if (!mHFCuts -> isGoodKaon(kaon)) continue;
            StHFPair *pair = new StHFPair(pion1, kaon, mHFCuts->getHypotheticalMass(StPicoCutsBase::kPion),mHFCuts->getHypotheticalMass(StPicoCutsBase::kKaon), i, j, mPrimVtx, mBField, kTRUE);
            if (!mHFCuts->isGoodSecondaryVertexPair(pair)) continue;

            float flag = -99.;

            if( kaon->charge()<0 && pion1->charge()>0 ) flag=0.; // -+
            if( kaon->charge()>0 && pion1->charge()<0 ) flag=1.; // +-

            if( kaon->charge()<0 && pion1->charge()<0) flag=4.; // --
            if( kaon->charge()>0 && pion1->charge()>0) flag=5.; // ++

            int ii=0;
            float ntVar[27];
            ntVar[ii++] = mPicoDst->event()->refMult();
            ntVar[ii++] = mPicoEvent->runId();
            ntVar[ii++] = mPicoEvent->eventId();

            ntVar[ii++] = pion1->gPt();
            ntVar[ii++] = pair->particle1Dca();
            ntVar[ii++] = pion1->nSigmaPion();
            ntVar[ii++] = pion1->nHitsFit();
            ntVar[ii++] = mHFCuts->getOneOverBeta(pion1, mHFCuts->getTofBetaBase(pion1), StPicoCutsBase::kPion);
            ntVar[ii++] = mHFCuts->getTofBetaBase(pion1);

            ntVar[ii++] = kaon->gPt();
            ntVar[ii++] = pair->particle2Dca();
            ntVar[ii++] = kaon->nSigmaKaon();
            ntVar[ii++] = kaon->nHitsFit();
            ntVar[ii++] = mHFCuts->getOneOverBeta(kaon, mHFCuts->getTofBetaBase(kaon), StPicoCutsBase::kKaon);
            ntVar[ii++] = mHFCuts->getTofBetaBase(kaon);

            ntVar[ii++] = pair->dcaDaughters();
            ntVar[ii++] = flag;
            ntVar[ii++] = mPrimVtx.z();
            ntVar[ii++] = mPicoEvent->vzVpd();
            ntVar[ii++] = fabs(mPicoEvent->primaryVertex().z() - mPicoEvent->vzVpd());

            ntVar[ii++] = pair->pointingAngle();
            ntVar[ii++] = cos(pair->pointingAngle());
            ntVar[ii++] = pair->decayLength();
            ntVar[ii++] = pair->DcaToPrimaryVertex(); //(pair->decayLength())*sin(pair->pointingAngle());
            ntVar[ii++] = pair->cosThetaStar();

            ntVar[ii++] = pair->pt();
            ntVar[ii++] = pair->m();

            if ((flag == 0) || (flag == 1)) {
                ntp_DMeson_Signal->Fill(ntVar);
                if(!getCorV2(pair, 1)) continue;
            } else {
                ntp_DMeson_Background->Fill(ntVar);
                if(!getCorV2(pair, 1)) continue;
            }
        }  // for (unsigned short idxKaon = 0; idxKaon < mIdxPicoKaons.size(); ++idxKaon)
    } // for (unsigned short idxPion1 = 0; idxPion1 < mIdxPicoPions.size(); ++idxPion1)

    return kStOK;
}

int StPicoD0AnaMaker::analyzeCandidates() {
    TH1D *h_tracktest = static_cast<TH1D*>(mOutList->FindObject("h_tracktest"));
    TH1D *h_tracktest_TOF = static_cast<TH1D*>(mOutList->FindObject("h_tracktest_TOF"));

    const char *aNames[]   = {"all", "NHitFit>20", "pT>0.6", "dca<1.5","TPC pion", "TPC kaon"};
    const char *aNamesTOF[]   = {"all TOF match", "NHitFit>20 & TOF", "pT>0.6 & TOF","dca<1.5 & TOF","TPC & TOF pion", "TPC & TOF kaon"};

    for (unsigned int ii = 0; ii < 6; ii++) {
        h_tracktest->GetXaxis()->SetBinLabel(ii+1, aNames[ii]);
        h_tracktest_TOF->GetXaxis()->SetBinLabel(ii+1, aNamesTOF[ii]);
    }

    for(unsigned int i=0;i<mPicoDst->numberOfTracks();i++) {
        StPicoTrack const *t = mPicoDst->track(i);
        if (!t) continue;
        if (!t->isHFTTrack()) continue;
        h_tracktest->Fill(1);
        if (mHFCuts->isGoodTrack(t)) h_tracktest->Fill(2); // NhitsFit

        float pt=t->gPt();
        if ((pt>0.6) && mHFCuts->isGoodTrack(t)) h_tracktest->Fill(3);
        float dca = (mPrimVtx - t->dcaPoint()).mag();
        if (dca<1.5 && (pt>0.6) && mHFCuts->isGoodTrack(t)) h_tracktest->Fill(4);
        bool tpcPion = mHFCuts->isTPCHadron(t, StPicoCutsBase::kPion);
        bool tpcKaon = mHFCuts->isTPCHadron(t, StPicoCutsBase::kKaon);
        if ((pt>0.6) && (dca<1.5 ) && (mHFCuts->isGoodTrack(t))) {
            if (tpcPion) h_tracktest->Fill(5);
            if (tpcKaon) h_tracktest->Fill(6);
        }
        if (mHFCuts->isTOFmatched(t)) {
            h_tracktest_TOF->Fill(1);
            if (mHFCuts->isGoodTrack(t)) h_tracktest_TOF->Fill(2); // NhitsFit
            if (pt>0.6 && mHFCuts->isGoodTrack(t)) h_tracktest_TOF->Fill(3);
            if (dca<1.5 && pt>0.6 && mHFCuts->isGoodTrack(t)) h_tracktest_TOF->Fill(4);
            if ((pt>0.6) && (dca<1.5 ) && (mHFCuts->isGoodTrack(t))) {
//                if (tpcPion && mHFCuts->isTOFHadronPID(t, mHFCuts->getTofBetaBase(t), StPicoCutsBase::kPion))  h_tracktest_TOF->Fill(5);
                if (tpcPion)  h_tracktest_TOF->Fill(5);
                if (tpcKaon && mHFCuts->isTOFHadronPID(t, mHFCuts->getTofBetaBase(t), StPicoCutsBase::kKaon))  h_tracktest_TOF->Fill(6);
            }
        }
    }
//    for (unsigned short idxPion1 = 0; idxPion1 < mIdxPicoPions.size(); ++idxPion1) {
//        StPicoTrack const *t = mPicoDst->track(mIdxPicoPions[idxPion1]);
//        ntp_pion->Fill(t->gPt(), t->gMom().phi(), t->gMom().pseudoRapidity(), t->nSigmaPion(), t->nHitsFit(), getOneOverBeta(t, mHFCuts->getTofBetaBase(t), StPicoCutsBase::kPion), mPicoEvent->eventId(), mPicoEvent->runId());
//    }
//
//    for (unsigned short idxKaon = 0; idxKaon < mIdxPicoKaons.size(); ++idxKaon) {
//        StPicoTrack const *t = mPicoDst->track(mIdxPicoKaons[idxKaon]);
//        ntp_kaon->Fill(t->gPt(), t->gMom().phi(), t->gMom().pseudoRapidity(), t->nSigmaKaon(), t->nHitsFit(), getOneOverBeta(t, mHFCuts->getTofBetaBase(t), StPicoCutsBase::kKaon), mPicoEvent->eventId(), mPicoEvent->runId());
//    }
    return kStOK;
}


void StPicoD0AnaMaker::DeclareHistograms() {
    // for v2 calcualtion
  TString flatten[5];
  flatten[0] = "v2";
  flatten[1] = "cosD";
  flatten[2] = "sinD";
  flatten[3] = "cosHadron";
  flatten[4] = "sinHadron";
  TString sb[8] = {"s1like","s3like","hSBlike","lSBlike","s1unlike","s3unlike","hSBunlike","lSBunlike"};
  // float xbin[7] = {0,1,2,3,4,5,10};
  const int xbinSize=5;
  float xbin[6] = {0,1,2,3,5,10};
  float binMass[2001];
  float binPhi[2001];
  //candPt = new TProfile("candPt","",xbinSize,xbin);
  for(int i=0;i<2001;i++)
    binPhi[i] = 0.005*i-5;
  for(int i=0;i<2001;i++)
    binMass[i] = 0.01*i;
  float xWeight[6] = {0,7,12,16,22,100};
  for(int i=0;i!=8;i++)
  {
    for(int k=0;k!=3;k++)
    {
      for(int j=0;j!=5;j++)
      {
        TString name = sb[i]+flatten[j]+Form("_%i",k);
        profV2[i][j][k] = new TProfile(name.Data(),"",xbinSize,xbin);
        profV2[i][j][k]->Sumw2();
      }
      TString weightName = sb[i]+Form("_%i_weigth",k);
      v2Weight[i][k] = new TH2D(weightName.Data(),"",5,xWeight,xbinSize,xbin);
      v2Weight[i][k]->Sumw2();

      TString namehPhi = "hadronPhi_"+sb[i]+Form("_%i",k);
      TString nameDPhi = "DPhi_"+sb[i]+Form("_%i",k);
      // hPhiHadron[i][k] = new TH2F(namehPhi.Data(),"",2000,binPhi,xbinSize,xbin);
      // hPhiD[i][k]= new TH2F(nameDPhi.Data(),"",2000,binPhi,xbinSize,xbin);
      // hPhiD[i][k]->Sumw2();
      // hPhiHadron[i][k]->Sumw2();
    }
  }
  float ptbin1[12] = {0.225,0.375,0.525,0.675,0.825,0.975,1.12,1.27,1.42,1.58,1.73,1.88};
  float ptbin2[11];
  for(int i=0;i<11;i++)
    ptbin2[i] = 0.5*(ptbin1[i]+ptbin1[i+1]);
  for(int k=0;k<3;k++)
  {
    for(int i=0;i<5;i++)
    {
      hadronV2[i][k] = new TH1D(Form("hadron_%s_%i",flatten[i].Data(),k),"",5,xWeight);
      hadronV2[i][k]->Sumw2();
      hadronV2_sum[i][k] = new TH1D(Form("hadronsum_%s_%i",flatten[i].Data(),k),"",5,xWeight);
      hadronV2_sum[i][k]->Sumw2();
      /*
      for(int j=0;j<9;j++)
      {
        hadronV2_excl[i][j][k] = new TH1D(Form("hadron_%s_cent%i_%i",flatten[i].Data(),j,k),"",10,ptbin2);
        hadronV2_excl[i][j][k]->Sumw2();
        hadronV2_excl_sum[i][j][k] = new TH1D(Form("hadronsum_%s_cent%i_%i",flatten[i].Data(),j,k),"",10,ptbin2);
        hadronV2_excl_sum[i][j][k]->Sumw2();
      }
      */
    }
}

  double fitmean[6] = {1.85921,1.8633,1.86403,1.86475,1.86252,1.86534};
  double fitsigma[6] = {0.018139,0.0139476,0.0158346,0.0169282,0.0199567,0.0189131};
  // ifstream ifs("efficiency.txt");
  // for(int i=0; i<6; i++)
  //   for(int j=0; j<4; j++)
  //     ifs>>efficiency[j][i];
  for(int i=0;i<6;i++)//pt bin
  {
    for(int j=0;j<5;j++)//flatten
    {
      TString massName[2];
      massName[0] = Form("likeMass%i",i)+flatten[j];
      massName[1] = Form("unlikeMass%i",i)+flatten[j];
      /*
      for(int k=0;k<2;k++)
      {
        V2Mass[k][i][j] = new TProfile(massName[k].Data(),"",18,fitmean[i]-9*fitsigma[i],fitmean[i]+9*fitsigma[i]);
        V2Mass[k][i][j]->Sumw2();
      }
      */
    }
  }

      printf("Histograms declared! \n");//
}

void StPicoD0AnaMaker::WriteHistograms() {
   //Saving for v2 calculation
  for(int i=0;i!=8;i++)
  {
    for(int k=0;k!=3;k++)
    {
      for(int j=0;j!=5;j++)
      {
        profV2[i][j][k]->Write();
      }
      v2Weight[i][k]->Write();
    }
  }
  /*
  for(int i=0;i<6;i++)
  {
    for(int j=0;j<5;j++)
    {
      for(int k=0;k<2;k++)
        V2Mass[k][i][j]->Write();
    }
  }
  */
  //massLike->Write();
  //candPt->Write();
  //massUnlike->Write();
  for(int k=0;k<3;k++)
  {
    for(int i=0;i<5;i++)
    {
      hadronV2[i][k]->Write();
      hadronV2_sum[i][k]->Write();
    }
  }

      printf("Histograms written! \n");
}

bool StPicoD0AnaMaker::getHadronCorV2(int idxGap)
{
  double etaGap[3] = {0,0.15,0.05};
  double mEtaGap = etaGap[idxGap];
  float hadronFill[7] = {0};
  const double reweight = 1;//mGRefMultCorrUtil->getWeight();
  // int centrality  = mGRefMultCorrUtil->getCentralityBin9();
  StPicoEvent *event = (StPicoEvent *)mPicoDst->event();
  int mult = event->grefMult();
  for(unsigned int i=0;i<mPicoDst->numberOfTracks();++i)
  {
    StPicoTrack const* hadron = mPicoDst->track(i);
    if(hadron->pMom().perp()<0.2) continue;
    if(!isGoodHadron(hadron)) continue;
    float etaHadron = hadron->pMom().pseudoRapidity();
    float phiHadron = hadron->pMom().phi();
    if(etaHadron<-0.5*mEtaGap)//backward sample 
    {
      hadronFill[0]++;
      hadronFill[1] += sin(2 * phiHadron);
      hadronFill[2] += cos(2 * phiHadron);
    }			
    if(etaHadron>0.5*mEtaGap)//forward sample
    {
      hadronFill[3]++;
      hadronFill[4] += sin(2 * phiHadron);
      hadronFill[5] += cos(2 * phiHadron);
    }		
  }
  hadronFill[6] = mult;
  hadronFill[7] = reweight;
  //mHadronTuple->Fill(hadronFill);
  printf("mult %i \n", mult);
  printf("hadron fill 0 %f       1 %f        2 %f         3 %f          4 %f        5 %f \n", hadronFill[0], hadronFill[1], hadronFill[2], hadronFill[3], hadronFill[4], hadronFill[5]);
  if(hadronFill[0]==0 || hadronFill[3]==0)
    return false; 
  double temp = (hadronFill[1]*hadronFill[4]+hadronFill[2]*hadronFill[5]);
  hadronV2[0][idxGap]->Fill(mult,temp*reweight);
  hadronV2[1][idxGap]->Fill(mult,hadronFill[2]*reweight);
  hadronV2[2][idxGap]->Fill(mult,hadronFill[1]*reweight);
  hadronV2[3][idxGap]->Fill(mult,hadronFill[5]*reweight);
  hadronV2[4][idxGap]->Fill(mult,hadronFill[4]*reweight);
  hadronV2_sum[0][idxGap]->Fill(mult,hadronFill[0]*hadronFill[3]*reweight);
  hadronV2_sum[1][idxGap]->Fill(mult,hadronFill[0]*reweight);
  hadronV2_sum[2][idxGap]->Fill(mult,hadronFill[0]*reweight);
  hadronV2_sum[3][idxGap]->Fill(mult,hadronFill[3]*reweight);
  hadronV2_sum[4][idxGap]->Fill(mult,hadronFill[3]*reweight);

  //    StPicoTrack const* hadron = picoDst->track(i);
  //  hadronV2_excl[0][centrality]->Fill(hadron->pMom().perp(),temp*reweight);
  //  hadronV2_excl[1][centrality]->Fill(hadron->pMom().perp(),hadronFill[2]*reweight);
  //  hadronV2_excl[2][centrality]->Fill(hadron->pMom().perp(),hadronFill[1]*reweight);
  //  hadronV2_excl[3][centrality]->Fill(hadron->pMom().perp(),hadronFill[5]*reweight);
  //  hadronV2_excl[4][centrality]->Fill(hadron->pMom().perp(),hadronFill[4]*reweight);
  //  hadronV2_excl_sum[0][centrality]->Fill(hadron->pMom().perp(),hadronFill[0]*hadronFill[3]*reweight);
  //  hadronV2_excl_sum[1][centrality]->Fill(hadron->pMom().perp(),hadronFill[0]*reweight);
  //  hadronV2_excl_sum[2][centrality]->Fill(hadron->pMom().perp(),hadronFill[0]*reweight);
  //  hadronV2_excl_sum[3][centrality]->Fill(hadron->pMom().perp(),hadronFill[3]*reweight);
  //  hadronV2_excl_sum[4][centrality]->Fill(hadron->pMom().perp(),hadronFill[3]*reweight);
      printf("GetCor done! \n");
  return true;
}

bool StPicoD0AnaMaker::isGoodHadron(StPicoTrack const * const trk) const
{
  //return trk->pMom().perp() > mycuts::hadronPtMin &&trk->pMom().perp() < mycuts::hadronPtMax && trk->nHitsFit() >= mycuts::nHitsFit &&fabs(trk->pMom().pseudoRapidity())<1.&&fabs(trk->nSigmaElectron())>3 && (1.0*trk->nHitsFit()/trk->nHitsMax())>0.52;
  int tofIndex = trk->bTofPidTraitsIndex();
  bool TofMatch = kFALSE;
  StPicoBTofPidTraits* tofPidTraits;
  if (tofIndex >= 0)  tofPidTraits = mPicoDst->btofPidTraits(tofIndex); //GNX
  if (tofIndex >= 0 && tofPidTraits && tofPidTraits->btofMatchFlag() > 0)  TofMatch = kTRUE;
  return TofMatch && trk->pMom().perp() > 0.2 &&trk->pMom().perp() < 2.0 && trk->nHitsFit() >= 15 &&fabs(trk->pMom().pseudoRapidity())<1. && (1.0*trk->nHitsFit()/trk->nHitsMax())>0.52;
}



bool StPicoD0AnaMaker::getCorV2(StHFPair *kp,double weight)
{
  // int centrality  = mGRefMultCorrUtil->getCentralityBin9();
  StPicoEvent *event = (StPicoEvent *)mPicoDst->event();
  int mult = event->grefMult();
  // TClonesArray const * aKaonPion = mPicoD0Event->kaonPionArray();
  // StKaonPion const* kp = (StKaonPion*)aKaonPion->At(idxCand);
  StPicoTrack const* kaon = mPicoDst->track(kp->particle1Idx());
  StPicoTrack const* pion = mPicoDst->track(kp->particle2Idx());
  int charge = kaon->charge() * pion->charge();
  double dMass = kp->m();


  if(kp->pt()>10) return false;
  int ptIdx = 5;
  if(kp->pt()<5)
    ptIdx = static_cast<int>(kp->pt());
  double fitmean[6] = {1.85921,1.8633,1.86403,1.86475,1.86252,1.86534};
  double fitsigma[6] = {0.018139,0.0139476,0.0158346,0.0169282,0.0199567,0.0189131};
  // double fitmean[6] = {1.8620,1.8647,1.8617,1.86475,1.8608,1.8626};
  // double fitsigma[6] = {0.0190,0.0143,0.0143,0.0169282,0.0199567,0.0189131};
  double mean = fitmean[ptIdx];
  double sigma = fitsigma[ptIdx];
  bool fillSB[8];
  printf("charge %i \n", charge);
  printf("mass %f \n", dMass);
  fillSB[0] =  (charge>0)&& (dMass>(mean-1*sigma)) &&  (dMass<(mean+1*sigma));
  fillSB[1] =  (charge>0)&& (dMass>(mean-3*sigma)) &&  (dMass<(mean+3*sigma));
  fillSB[2] =  (charge>0) && (((dMass>(mean+4*sigma)) &&  (dMass<(mean+9*sigma))) ||((dMass>(mean-9*sigma)) &&  (dMass<(mean-4*sigma))));
  fillSB[4] = (charge==-1)&& (((dMass>(mean+5.5*sigma)) &&  (dMass<(mean+7.5*sigma))) ||((dMass>(mean-7.5*sigma)) &&  (dMass<(mean-5.5*sigma))));
  fillSB[5] = (charge==-1)&& (dMass>(mean-3*sigma)) &&  (dMass<(mean+3*sigma));
  fillSB[6] = (charge==-1)&& (((dMass>(mean+4*sigma)) &&  (dMass<(mean+9*sigma))) ||((dMass>(mean-9*sigma)) &&  (dMass<(mean-4*sigma))));
  fillSB[3] = fillSB[1] || fillSB[2];
  fillSB[7] = fillSB[1] || fillSB[2] || fillSB[6];
  double etaGap[3] = {0,0.15,0.05};

  for(int k=0;k<3;k++)
  {
    double corFill[7] = {0};
    corFill[0] = 1 ;
    // corFill[1] = sin(2* kp->phi())/sqrt(hadronv2[centBin]);
    // corFill[2] = cos(2* kp->phi())/sqrt(hadronv2[centBin]);
    corFill[1] = sin(2* kp->phi());
    corFill[2] = cos(2* kp->phi());
    int chargeIdx = charge>0 ? 0:1;
    for(int j=0;j<8;j++)
    {
      if(fillSB[j])
      {
        profV2[j][1][k]->Fill(kp->pt(),corFill[1],weight);
        profV2[j][2][k]->Fill(kp->pt(),corFill[2],weight);
      }
    }
    for(unsigned int i=0; i<mPicoDst->numberOfTracks();i++)
    {
      StPicoTrack const* hadron = mPicoDst->track(i);
      if(hadron->pMom().perp()<0.2) continue;
      if(!isGoodHadron(hadron)) continue;
      if(i==kp->particle1Idx() || i==kp->particle2Idx()) continue;
      float etaHadron = hadron->pMom().pseudoRapidity();
      float phiHadron = hadron->pMom().phi();
      if(!isEtaGap(kp->eta(),etaGap[k],etaHadron))  continue;
      corFill[3]++;
      corFill[4] += sin(2 * phiHadron);
      corFill[5] += cos(2 * phiHadron);
      for(int j=0;j<8;j++)
      {
        if(fillSB[j])
        {
          printf("som tutok! \n");
          profV2[j][3][k]->Fill(kp->pt(),sin(2*phiHadron),weight);s
          profV2[j][4][k]->Fill(kp->pt(),cos(2*phiHadron),weight);
          profV2[j][0][k]->Fill(kp->pt(),cos(2*(phiHadron-kp->phi())),weight);
          //v2Weight[j][k]->Fill(mult,kp->pt(),weight);
        }
      }
    }// Loop over charged tracks
    // if(corFill[3]<=0) return false;
    // double cumulant = (corFill[1]*corFill[4]+corFill[2]*corFill[5])/(corFill[3]); 
    // cout<<"kp phi = "<<kp->phi()<<"\tcumulant = "<<cumulant<<endl;
    // cout<<"number of correlated hadrons = "<<corFill[3]<<endl;
    // for(int j=0;j<8;j++)
    // {
    //   if(fillSB[j])
    //   {
    //     // testV2->Fill(kp->pt(),cumulant);
    //   }
    // }
  }//Loop over different eta gap (k)
  return true;
}


bool StPicoD0AnaMaker::isEtaGap(double dEta,double mGap,double hEta)
{
  if(mGap == 0) return true;
  double range =  2. - mGap*2;
  // if(dEta> (1.-2*mGap))
  //   return hEta<(dEta-mGap) && hEta>(dEta-mGap-range);
  // else if(dEta<(-1.+2*mGap))
  //   return hEta>(dEta+mGap) && hEta<(dEta+mGap+range);
  // else 
  //   return (hEta>(dEta+mGap) || hEta<(dEta-mGap));
  if(dEta>0)
    return hEta<-mGap;
  else
    return hEta>mGap;
}