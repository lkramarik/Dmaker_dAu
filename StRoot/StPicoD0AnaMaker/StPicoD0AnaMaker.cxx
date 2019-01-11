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
	createCandidates();
	getHadronCorV2(1);
//    analyzeCandidates();

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

            //double reweight = mGRefMultCorrUtil->getWeight(); THIS DOESNT WORK
            double reweight = 1;
            float d0Pt = pair->pt();
            double dMass = pair->m();
            if ((flag == 0) || (flag == 1)) {
                ntp_DMeson_Signal->Fill(ntVar);
                if(dMass>1.81&&dMass<1.91)
            		candPt->Fill(d0Pt,d0Pt,reweight);
        		massPt->Fill(dMass,d0Pt,reweight);
            } else {
                ntp_DMeson_Background->Fill(ntVar);
                massPtLike->Fill(dMass,d0Pt,reweight);
            }
            
            

            
            //this is messy!! I will rewrite it to a isD0pair function ASAP! 
            if(pair->m() < 1.804 || pair->m() > 1.924 || pair->pt() < 1 || pair->pt() > 5) continue;
            //mean 1.864, sigma 0.02
            if(pair->pt() > 1 && pair->pt() < 2)
            {
            	if(pair->decayLength() > 0.012 && pair->dcaDaughters() < 0.007 && pair->DcaToPrimaryVertex() < 0.005 && cos(pair->pointingAngle()) > 0.5 && pair->particle2Dca() > 0.007 && pair->particle1Dca() > 0.009);
            	{
            		getCorV2(pair, reweight);
            	}
            }
            if(pair->pt() > 2 && pair->pt() < 3)
            {
				if(pair->decayLength() > 0.003 && pair->dcaDaughters() < 0.016 && pair->DcaToPrimaryVertex() < 0.0065 && cos(pair->pointingAngle()) > 0.5 && pair->particle2Dca() > 0.01 && pair->particle1Dca() > 0.009);
            	{
            		getCorV2(pair, reweight);
            	}
            }
        	if(pair->pt() > 3 && pair->pt() < 5)
            {
            	if(pair->decayLength() > 0.009 && pair->dcaDaughters() < 0.015 && pair->DcaToPrimaryVertex() < 0.0064 && cos(pair->pointingAngle()) > 0.6 && pair->particle2Dca() > 0.0076 && pair->particle1Dca() > 0.0064);
            	{
            		getCorV2(pair, reweight);
            	}
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
    return kStOK;
}


void StPicoD0AnaMaker::DeclareHistograms() {
 
  TString names[4] = {"cos_B", "cos_F", "sin_B", "sin_F"}; //backward and forward samples
  float multBin[6] = {0,7,12,16,22,100};
  for(int m = 0; m < 4; m++)
  {
  	TString aa = names[m];
  	qVec[m] = new TProfile(aa.Data(),"Q vector", 5, multBin);
  	aa += "Pow2";
  	qVecPow2[m] = new TProfile(aa.Data(),"Q vector", 5, multBin);
  	qVec[m]->Sumw2();
  	qVecPow2[m]->Sumw2();
  	aa = names[m];
  	aa += "_no_mult";
  	qVec2[m] = new TProfile(aa.Data(),"Q vector", 1, 0, 100);
  }
  float momBins[7] = {0,1,2,3,4,5,10};
  TString multBinNames[6] = {"0","7","12","16","22","100"};
  for(int m = 0; m < 5; m++)
  {
  	TString aa = "cosD_" + multBinNames[m] + "_" + multBinNames[m+1];
  	corrD[0][m] = new TProfile(aa.Data(),"",6,momBins);
  	aa = "sinD_" + multBinNames[m] + "_" + multBinNames[m+1];
  	corrD[1][m] = new TProfile(aa.Data(),"",6,momBins);
  	aa = "dirFlow_" + multBinNames[m] + "_" + multBinNames[m+1];
  	dirFlow[m] = new TProfile(aa.Data(), "", 6 , momBins);
  }
  corrD2[0] = new TProfile("cosD_no_mult","cos D",6,momBins);	
  corrD2[1] = new TProfile("sinD_no_mult","sin D",6,momBins);
  dirFlow2 = new TProfile("dirFlow_no_mult","dir flow",6,momBins);
  refFlow = new TProfile("refFlow", "", 5, multBin);
  refFlow2 = new TProfile("refFlow_no_mult", "", 1, 0, 100);
  
  
  hadron_phi = new TH1D("hadron_phi", "Hadron phi", 2000, -5, 5);
  D_phi = new TH1D("D_phi", "D phi", 2000, -5, 5);


}

void StPicoD0AnaMaker::WriteHistograms() {
 
  for(int m = 0; m < 4; m++)
  {
  	qVec[m]->Write();
  	qVecPow2[m]->Write();
  	qVec2[m]->Write();
  }
  refFlow->Write();
  
  for(int m = 0; m < 5; m++)
  {
  	corrD[0][m]->Write();
  	corrD[1][m]->Write();
  	dirFlow[m]->Write();
  }
  corrD2[0]->Write();
  corrD2[1]->Write();
  dirFlow2->Write();
  refFlow2->Write();

  hadron_phi->Write();
  D_phi->Write();

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
    if(!mHFCuts->isGoodTrack(hadron)) continue;
    if(!mHFCuts->isGoodProton(hadron) || !mHFCuts->isGoodKaon(hadron) || !mHFCuts->isGoodPion(hadron)) continue;
    float etaHadron = hadron->gMom().pseudoRapidity();
    float phiHadron = hadron->gMom().phi();
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
  if(hadronFill[0]==0 || hadronFill[3]==0)
    return false; 

  //Z code: reference flow creation: average sin/cos phi of a hadron in an event.... (no error!)  
  hadron_phi->Fill(phiHadron);

  if(idxGap==1)
  {
 	qVec[0]->Fill(mult,hadronFill[2]/hadronFill[0],reweight);
  	qVec[1]->Fill(mult,hadronFill[5]/hadronFill[3],reweight);
  	qVec[2]->Fill(mult,hadronFill[1]/hadronFill[0],reweight);
  	qVec[3]->Fill(mult,hadronFill[4]/hadronFill[3],reweight);
  	qVecPow2[0]->Fill(mult,(hadronFill[2]*hadronFill[2])/(hadronFill[0]*hadronFill[0]),reweight);
  	qVecPow2[1]->Fill(mult,(hadronFill[5]*hadronFill[5])/(hadronFill[3]*hadronFill[3]),reweight);
  	qVecPow2[2]->Fill(mult,(hadronFill[1]*hadronFill[1])/(hadronFill[0]*hadronFill[0]),reweight);
  	qVecPow2[3]->Fill(mult,(hadronFill[4]*hadronFill[4])/(hadronFill[3]*hadronFill[3]),reweight);
  	refFlow->Fill(mult,((hadronFill[2]*hadronFill[5])/(hadronFill[0]*hadronFill[3])),reweight);
  	//no mult
  	qVec2[0]->Fill(mult,hadronFill[2]/hadronFill[0],reweight);
  	qVec2[1]->Fill(mult,hadronFill[5]/hadronFill[3],reweight);
  	qVec2[2]->Fill(mult,hadronFill[1]/hadronFill[0],reweight);
  	qVec2[3]->Fill(mult,hadronFill[4]/hadronFill[3],reweight);
  	refFlow2->Fill(mult,((hadronFill[2]*hadronFill[5])/(hadronFill[0]*hadronFill[3])),reweight);
  }
  return true;
}


bool StPicoD0AnaMaker::getCorV2(StHFPair *kp,double weight)
{
  StPicoEvent *event = (StPicoEvent *)mPicoDst->event();
  int mult = event->grefMult();
  StPicoTrack const* kaon = mPicoDst->track(kp->particle1Idx());
  StPicoTrack const* pion = mPicoDst->track(kp->particle2Idx());
  int charge = kaon->charge() * pion->charge();
  double dMass = kp->m();
  double hadronv2=1;
  float multBin[6] = {0,7,12,16,22,100};

  if(kp->pt()>10) return false;
  int ptIdx = 5;
  if(kp->pt()<5)
    ptIdx = static_cast<int>(kp->pt());

  double etaGap[3] = {0,0.15,0.05};

  for(int k=0;k<3;k++)
  {
    double corFill[7] = {0};
    corFill[0] = 1 ;
    corFill[1] = sin(2* kp->phi())/sqrt(hadronv2);
    corFill[2] = cos(2* kp->phi())/sqrt(hadronv2);

    int chargeIdx = charge>0 ? 0:1;

    D_phi->Fill(kp->phi());
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
      corFill[4] += sin(2 * phiHadron)/sqrt(hadronv2);
      corFill[5] += cos(2 * phiHadron)/sqrt(hadronv2);
      
  
      if(k==1)
      {
      	for(int m = 0; m < 5; m++)
      	{
      		if(mult >= multBin[m] && mult < multBin[m+1])
      		{
      			corrD[0][m]->Fill(kp->pt(),corFill[2],weight);
      			corrD[1][m]->Fill(kp->pt(),corFill[1],weight);
      			dirFlow[m]->Fill(kp->pt(),corFill[2]*corFill[5]/corFill[3],weight);
      		} 
      	}
      	corrD2[0]->Fill(kp->pt(),corFill[2],weight);
      	corrD2[1]->Fill(kp->pt(),corFill[1],weight);
      	dirFlow2->Fill(kp->pt(),corFill[2]*corFill[5]/corFill[3],weight);
      }
    }
    
  }
  return true;
}


bool StPicoD0AnaMaker::isEtaGap(double dEta,double mGap,double hEta)
{
  if(mGap == 0) return true;
  //double range =  2. - mGap*2;
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