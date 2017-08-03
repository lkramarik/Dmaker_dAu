#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoDstMaker/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
//#include "StPicoD0EventMaker/StPicoD0Event.h"
//#include "StPicoD0EventMaker/StKaonPion.h"
//#include "StPicoD0AnaMaker.h"
#include "StPicoHFMaker/StHFCuts.h"

#include "StPicoDpmAnaMaker.h"
ClassImp(StPicoDpmAnaMaker)

// _________________________________________________________
StPicoDpmAnaMaker::StPicoDpmAnaMaker(char const* name, StPicoDstMaker* picoMaker, char const* outputBaseFileName,  
				       char const* inputHFListHFtree = "") :
  StPicoHFMaker(name, picoMaker, outputBaseFileName, inputHFListHFtree),
   mOutFileBaseName(outputBaseFileName){ //mDecayChannel(kChannel1), tu bolo
   
  // constructor
}

// _________________________________________________________
StPicoDpmAnaMaker::~StPicoDpmAnaMaker() {
  // destructor
}

// _________________________________________________________
int StPicoDpmAnaMaker::InitHF() {
  // -- INITIALIZE USER HISTOGRAMS ETC HERE -------------------
  //    add them to the output list mOutList which is automatically written
  //
  // EXAMPLE //  mOutList->Add(new TH1F(...));
  // EXAMPLE //  TH1F* hist = static_cast<TH1F*>(mOutList->Last());
  mOutList->Add(new TH1F("h_time_per_event","h_time_per_event", 2000., 0., 20.));
  mOutList->Add(new TH2F("h_piTOF","h_piTOF",100,0,10, 2, 0, 2));
  mOutList->Add(new TH2F("h_kTOF","h_kTOF",100,0,10, 2, 0, 2));
  mOutList->Add(new TH2F("h_pTOF","h_pTOF",100,0,10, 2, 0, 2));

  mOutList->Add(new TH2F("h_piTOF_20","h_piTOF_20",100,0,10, 2, 0, 2));
  mOutList->Add(new TH2F("h_kTOF_20","h_kTOF_20",100,0,10, 2, 0, 2));
  mOutList->Add(new TH2F("h_pTOF_20","h_pTOF_20",100,0,10, 2, 0, 2));

  mOutList->Add(new TH2F("h_piTOF_1sig","h_piTOF_1sig",100,0,10, 2, 0, 2));
  mOutList->Add(new TH2F("h_kTOF_1sig","h_kTOF_1sig",100,0,10, 2, 0, 2));
  mOutList->Add(new TH2F("h_pTOF_1sig","h_pTOF_1sig",100,0,10, 2, 0, 2));

  mOutList->Add(new TH2F("h_piTOF_HFT","h_piTOF_HFT",100,0,10, 2, 0, 2));
  mOutList->Add(new TH2F("h_kTOF_HFT","h_kTOF_HFT",100,0,10, 2, 0, 2));
  mOutList->Add(new TH2F("h_pTOF_HFT","h_pTOF_HFT",100,0,10, 2, 0, 2));

  mOutList->Add(new TH2F("h_piTOF_HFT_1sig","h_piTOF_HFT_1sig",100,0,10, 2, 0, 2));
  mOutList->Add(new TH2F("h_kTOF_HFT_1sig","h_kTOF_HFT_1sig",100,0,10, 2, 0, 2));
  mOutList->Add(new TH2F("h_pTOF_HFT_1sig","h_pTOF_HFT_1sig",100,0,10, 2, 0, 2));

  mOutList->Add(new TH2F("h_piTOF_HFT_20","h_piTOF_HFT_20",100,0,10, 2, 0, 2));
  mOutList->Add(new TH2F("h_kTOF_HFT_20","h_kTOF_HFT_20",100,0,10, 2, 0, 2));
  mOutList->Add(new TH2F("h_pTOF_HFT_20","h_pTOF_HFT_20",100,0,10, 2, 0, 2));

  mOutList->Add(new TH2F("h_piTOF_HFT_1sig_20","h_piTOF_HFT_1sig_20",100,0,10, 2, 0, 2));
  mOutList->Add(new TH2F("h_kTOF_HFT_1sig_20","h_kTOF_HFT_1sig_20",100,0,10, 2, 0, 2));
  mOutList->Add(new TH2F("h_pTOF_HFT_1sig_20","h_pTOF_HFT_1sig_20",100,0,10, 2, 0, 2));

  mOutList->Add(new TH2F("h_piTOFbeta","h_piTOFbeta",500,0,10, 500, 0, 5));
  mOutList->Add(new TH2F("h_kTOFbeta","h_kTOFbeta",500,0,10, 500, 0, 5));
  mOutList->Add(new TH2F("h_pTOFbeta","h_pTOFbeta",500,0,10, 500, 0, 5));

  mOutList->Add(new TH2F("h_pinsigma","h_pinsigma",1000,0,10, 99, -5, 5));
  mOutList->Add(new TH2F("h_knsigma","h_knsigma",1000,0,10, 99, -5, 5));
  mOutList->Add(new TH2F("h_pnsigma","h_pnsigma",1000,0,10, 99, -5, 5));


  mOutList->Add(new TH2F("h_dedx","h_dedx", 1000, 0, 10, 1000, 0, 10));
 /* mOutList->Add(new TH2D("h_pidedx","h_pidedx", 500, 0, 10, 500, -10, 10);
  mOutList->Add(new TH2D("h_kdedx","h_kdedx", 500, 0, 10, 500, -10, 10));
  mOutList->Add(new TH2D("h_pdedx","h_pdedx", 500, 0, 10, 500, -10, 10));*/
 mOutFileBaseName = mOutFileBaseName.ReplaceAll(".root", "");

   mOutList->Add(new TH1F("h_mh1Cent", "EventsVsCentrality;cent;Counts", 10, -1.5, 8.5));
   mOutList->Add(new TH1F("h_mh1CentWg", "EventsVsCentrality;cent;Counts", 10, -1.5, 8.5));
   mOutList->Add(new TH1F("h_mh1gRefmultCor", "gRefmultCor;gRefmult;Counts", 700, 0, 700));
   mOutList->Add(new TH1F("h_mh1gRefmultCorWg", "gRefmultCorWg;gRefmultCorWg;Counts", 700, 0, 700));
   mOutList->Add(new TH2F("h_mh2CentVz", "CentralityVsVz;cent;Vz", 10, -1.5, 8.5, 200, -10, 10));
   mOutList->Add(new TH2F("h_mh2CentVzWg", "CentralityVsVzWg;cent;Vz", 10, -1.5, 8.5, 200, -10, 10));

   mOutList->Add(new TH1F("h_mh1Cent_run", "EventsVsCentrality;cent;Counts", 10, -1.5, 8.5));
   mOutList->Add(new TH1F("h_mh1CentWg_run", "EventsVsCentrality;cent;Counts", 10, -1.5, 8.5));
   mOutList->Add(new TH1F("h_mh1gRefmultCor_run", "gRefmultCor;gRefmult;Counts", 700, 0, 700));
   mOutList->Add(new TH1F("h_mh1gRefmultCorWg_run", "gRefmultCorWg;gRefmultCorWg;Counts", 700, 0, 700));
   mOutList->Add(new TH2F("h_mh2CentVz_run", "CentralityVsVz;cent;Vz", 10, -1.5, 8.5, 200, -10, 10));
   mOutList->Add(new TH2F("h_mh2CentVzWg_run", "CentralityVsVzWg;cent;Vz", 10, -1.5, 8.5, 200, -10, 10));
 
  //histoInit(mOutFileBaseName, true);

  
   // -------------- USER VARIABLES -------------------------

  ntp_DMeson = new TNtuple("ntp","DMeson Tree","grefMult:pi1_runId:pi1_eventId:pi1_phi:pi1_eta:pi1_pt:pi1_dca:pi1_dedx:pi1_nSigma:pi1_nHitFit:pi1_nHitdedx:pi1_TOFinvbeta:pi1_betaBase:k_runId:k_eventId:k_phi:k_eta:k_pt:k_dca:k_dedx:k_nSigma:k_nHitFit:k_nHitdedx:k_TOFinvbeta:k_betaBase:dcaMax:flag:primVz:D_theta:D_decayL:D_phi:D_eta:D_pt:D_mass:D_mass_LS:D_mass_US:D_dV0Max:centrality:refmult:refmultcorr:reweight");
  mRunNumber = 0;
  return kStOK;
}

// _________________________________________________________
void StPicoDpmAnaMaker::ClearHF(Option_t *opt="") {
  return;
}

// _________________________________________________________
int StPicoDpmAnaMaker::FinishHF() {
   if( isMakerMode() != StPicoHFMaker::kWrite )
    ntp_DMeson->Write();
   //closeFile();
  return kStOK;
}

// _________________________________________________________
int StPicoDpmAnaMaker::MakeHF() {
  //cout<<"MakeHF start"<<endl;  
// -- process event
  //    ADD YOUR PROCESSING CODE HERE
  //    ... it is usefull to use the methods below
  //     - createCandidates()
  //     - analyzeCandidates()

  std::clock_t start1 = std::clock();//kvapil
  if (isMakerMode() == StPicoHFMaker::kWrite) {
    createCandidates();
  }
  else if (isMakerMode() == StPicoHFMaker::kRead) {
    // -- the reading back of the perviously written trees happens in the background
    analyzeCandidates();
  }
  else if (isMakerMode() == StPicoHFMaker::kAnalyze) {
    //cout<<"going to create candidates"<<endl;
    createCandidates();
    //cout<<"candidated created, going to analyze them"<<endl;
    analyzeCandidates();
    //cout<<"candidates analysed"<<endl;
    //createQA();
  }


  //fill tof and tpc
  TH2F *h_piTOF = static_cast<TH2F*>(mOutList->FindObject("h_piTOF"));
  TH2F *h_kTOF = static_cast<TH2F*>(mOutList->FindObject("h_kTOF"));
  TH2F *h_pTOF = static_cast<TH2F*>(mOutList->FindObject("h_pTOF"));

  TH2F *h_piTOF_20 = static_cast<TH2F*>(mOutList->FindObject("h_piTOF_20"));
  TH2F *h_kTOF_20 = static_cast<TH2F*>(mOutList->FindObject("h_kTOF_20"));
  TH2F *h_pTOF_20 = static_cast<TH2F*>(mOutList->FindObject("h_pTOF_20"));

  TH2F *h_piTOF_1sig = static_cast<TH2F*>(mOutList->FindObject("h_piTOF_1sig"));
  TH2F *h_kTOF_1sig = static_cast<TH2F*>(mOutList->FindObject("h_kTOF_1sig"));
  TH2F *h_pTOF_1sig = static_cast<TH2F*>(mOutList->FindObject("h_pTOF_1sig"));

  TH2F *h_piTOF_HFT = static_cast<TH2F*>(mOutList->FindObject("h_piTOF_HFT"));
  TH2F *h_kTOF_HFT = static_cast<TH2F*>(mOutList->FindObject("h_kTOF_HFT"));
  TH2F *h_pTOF_HFT = static_cast<TH2F*>(mOutList->FindObject("h_pTOF_HFT"));

  TH2F *h_piTOF_HFT_1sig = static_cast<TH2F*>(mOutList->FindObject("h_piTOF_HFT_1sig"));
  TH2F *h_kTOF_HFT_1sig = static_cast<TH2F*>(mOutList->FindObject("h_kTOF_HFT_1sig"));
  TH2F *h_pTOF_HFT_1sig = static_cast<TH2F*>(mOutList->FindObject("h_pTOF_HFT_1sig"));

  TH2F *h_piTOF_HFT_20 = static_cast<TH2F*>(mOutList->FindObject("h_piTOF_HFT_20"));
  TH2F *h_kTOF_HFT_20 = static_cast<TH2F*>(mOutList->FindObject("h_kTOF_HFT_20"));
  TH2F *h_pTOF_HFT_20 = static_cast<TH2F*>(mOutList->FindObject("h_pTOF_HFT_20"));

  TH2F *h_piTOF_HFT_1sig_20 = static_cast<TH2F*>(mOutList->FindObject("h_piTOF_HFT_1sig_20"));
  TH2F *h_kTOF_HFT_1sig_20 = static_cast<TH2F*>(mOutList->FindObject("h_kTOF_HFT_1sig_20"));
  TH2F *h_pTOF_HFT_1sig_20 = static_cast<TH2F*>(mOutList->FindObject("h_pTOF_HFT_1sig_20"));

  TH2F *h_piTOFbeta = static_cast<TH2F*>(mOutList->FindObject("h_piTOFbeta"));
  TH2F *h_kTOFbeta = static_cast<TH2F*>(mOutList->FindObject("h_kTOFbeta"));
  TH2F *h_pTOFbeta = static_cast<TH2F*>(mOutList->FindObject("h_pTOFbeta"));

  TH2F *h_pinsigma = static_cast<TH2F*>(mOutList->FindObject("h_pinsigma"));
  TH2F *h_knsigma = static_cast<TH2F*>(mOutList->FindObject("h_knsigma"));
  TH2F *h_pnsigma = static_cast<TH2F*>(mOutList->FindObject("h_pnsigma"));

  TH2F *h_dedx = static_cast<TH2F*>(mOutList->FindObject("h_dedx"));

   TH1F *h_mh1Cent = static_cast<TH1F*>(mOutList->FindObject("h_mh1Cent"));
   TH1F *h_mh1CentWg = static_cast<TH1F*>(mOutList->FindObject("h_mh1CentWg"));
   TH1F *h_mh1gRefmultCor = static_cast<TH1F*>(mOutList->FindObject("h_mh1gRefmultCor"));
   TH1F *h_mh1gRefmultCorWg = static_cast<TH1F*>(mOutList->FindObject("h_mh1gRefmultCorWg"));
   TH2F *h_mh2CentVz = static_cast<TH2F*>(mOutList->FindObject("h_mh2CentVz"));
   TH2F *h_mh2CentVzWg = static_cast<TH2F*>(mOutList->FindObject("h_mh2CentVzWg"));

   TH1F *h_mh1Cent_run = static_cast<TH1F*>(mOutList->FindObject("h_mh1Cent_run"));
   TH1F *h_mh1CentWg_run = static_cast<TH1F*>(mOutList->FindObject("h_mh1CentWg_run"));
   TH1F *h_mh1gRefmultCor_run = static_cast<TH1F*>(mOutList->FindObject("h_mh1gRefmultCor_run"));
   TH1F *h_mh1gRefmultCorWg_run = static_cast<TH1F*>(mOutList->FindObject("h_mh1gRefmultCorWg_run"));
   TH2F *h_mh2CentVz_run = static_cast<TH2F*>(mOutList->FindObject("h_mh2CentVz_run"));
   TH2F *h_mh2CentVzWg_run = static_cast<TH2F*>(mOutList->FindObject("h_mh2CentVzWg_run"));

    StThreeVectorF pVtx = mPicoDst->event()->primaryVertex();
//     float multiplicity = mPicoDst->event()->refMult();
//     h_refMult -> Fill(multiplicity);

/*
      mRefmultCorrUtil->init(mPicoDst->event()->runId());
      if (!mRefmultCorrUtil){
         LOG_WARN << " No mGRefMultCorrUtil! Skip! " << endl;
         return kStWarn;
      }
      if (mRefmultCorrUtil->isBadRun(mPicoDst->event()->runId())) return kStOK;      
      mRefmultCorrUtil->setVzForWeight(6, -6.0, 6.0);
      mRefmultCorrUtil->readScaleForWeight("StRoot/StRefMultCorr/macros/weight_grefmult_vpd30_vpd5_Run14_P16id.txt");
      mRefmultCorrUtil->initEvent(mPicoDst->event()->grefMult(), mPrimVtx.z(), mPicoDst->event()->ZDCx()) ;
      int const centrality = mRefmultCorrUtil->getCentralityBin9();
      const double reweight = mRefmultCorrUtil->getWeight();
      const double refmultCor = mRefmultCorrUtil->getRefMultCorr();
      h_mh1gRefmultCor->Fill(refmultCor);
      h_mh1gRefmultCorWg->Fill(refmultCor, reweight);
      h_mh1Cent->Fill(centrality);
      h_mh1CentWg->Fill(centrality, reweight);
      h_mh2CentVz->Fill(centrality, pVtx.z());
      h_mh2CentVzWg->Fill(centrality, pVtx.z(), reweight);

	  if((mPicoDst->event()->runId()) >= 15107008){
		  h_mh1gRefmultCor_run->Fill(refmultCor);
		  h_mh1gRefmultCorWg_run->Fill(refmultCor, reweight);
		  h_mh1Cent_run->Fill(centrality);
		  h_mh1CentWg_run->Fill(centrality, reweight);
		  h_mh2CentVz_run->Fill(centrality, pVtx.z());
		  h_mh2CentVzWg_run->Fill(centrality, pVtx.z(), reweight);
	  }


*/
    UInt_t nTracks = mPicoDst->numberOfTracks();
    for (unsigned short iTrack = 0; iTrack < nTracks; ++iTrack){
        StPicoTrack const* trk = mPicoDst->track(iTrack);
        if (!trk) continue;
        StPhysicalHelixD helix = trk->helix(mPicoDst->event()->bField());
        StThreeVectorF momentum = trk->gMom(mPicoDst->event()->primaryVertex(), mPicoDst->event()->bField());

        if (!(trk->nHitsFit()>=20)) continue;
        if (!(fabs(momentum.pseudoRapidity()) <= 1.0)) continue;      

        if (fabs(trk->nSigmaPion())<3.0){
            float piBeta = mHFCuts->getTofBetaBase(trk);
            Int_t piTofAvailable = 0;
            if(!isnan(piBeta) && piBeta > 0){
		 		piTofAvailable = 1;
		 		float tofPion = fabs(1. / piBeta - sqrt(1+M_PION_PLUS*M_PION_PLUS/(momentum.mag()*momentum.mag())));
		 		h_piTOFbeta->Fill(trk->gPt(),tofPion);
	    	}
            h_piTOF->Fill(trk->gPt(),piTofAvailable);     
			if (trk->nHitsFit()>=20) h_piTOF_20->Fill(trk->gPt(),piTofAvailable);   
            if (trk->isHFTTrack()) h_piTOF_HFT->Fill(trk->gPt(),piTofAvailable); 
			if ((trk->isHFTTrack()) && (fabs(trk->nSigmaPion())<1.0)) h_piTOF_HFT_1sig->Fill(trk->gPt(),piTofAvailable); 
			if ((trk->isHFTTrack()) && (trk->nHitsFit()>=20)) h_piTOF_HFT_20->Fill(trk->gPt(),piTofAvailable); 
			if ((trk->isHFTTrack()) && (fabs(trk->nSigmaPion())<1.0) && (trk->nHitsFit()>=20)) h_piTOF_HFT_1sig_20->Fill(trk->gPt(),piTofAvailable); 
            if (fabs(trk->nSigmaPion())<1.0) h_piTOF_1sig->Fill(trk->gPt(),piTofAvailable); 
        }

        if (fabs(trk->nSigmaKaon())<3.0){
            float kBeta = mHFCuts->getTofBetaBase(trk);
            Int_t kTofAvailable = 0;
            if(!isnan(kBeta) && kBeta > 0){
			 	kTofAvailable = 1;
			 	float tofKaon = fabs(1. / kBeta - sqrt(1+M_KAON_PLUS*M_KAON_PLUS/(momentum.mag()*momentum.mag())));
		 		h_kTOFbeta->Fill(trk->gPt(),tofKaon);
	    	}
            h_kTOF->Fill(trk->gPt(),kTofAvailable);
			if (trk->nHitsFit()>=20) h_kTOF_20->Fill(trk->gPt(),kTofAvailable);
            if (trk->isHFTTrack()) h_kTOF_HFT->Fill(trk->gPt(),kTofAvailable); 
		    if ((trk->isHFTTrack()) && (fabs(trk->nSigmaKaon())<1.0)) h_kTOF_HFT_1sig->Fill(trk->gPt(),kTofAvailable); 
			if ((trk->isHFTTrack()) && (trk->nHitsFit()>=20)) h_kTOF_HFT_20->Fill(trk->gPt(),kTofAvailable); 
			if ((trk->isHFTTrack()) && (fabs(trk->nSigmaKaon())<1.0) && (trk->nHitsFit()>=20)) h_kTOF_HFT_1sig_20->Fill(trk->gPt(),kTofAvailable);
            if (fabs(trk->nSigmaKaon())<1.0) h_kTOF_1sig->Fill(trk->gPt(),kTofAvailable);  
        }

        if (fabs(trk->nSigmaProton())<3.0){
            float pBeta = mHFCuts->getTofBetaBase(trk);
            Int_t pTofAvailable = 0;
            if(!isnan(pBeta) && pBeta > 0){
		 		pTofAvailable = 1;
		 		float tofProton = fabs(1. / pBeta - sqrt(1+M_PROTON*M_PROTON/(momentum.mag()*momentum.mag())));
		 		h_pTOFbeta->Fill(trk->gPt(),tofProton);
	    	}
            h_pTOF->Fill(trk->gPt(),pTofAvailable);
			if (trk->nHitsFit()>=20) h_pTOF_20->Fill(trk->gPt(),pTofAvailable); 
			if (trk->isHFTTrack()) h_pTOF_HFT->Fill(trk->gPt(),pTofAvailable); 
			if ((trk->isHFTTrack()) && (fabs(trk->nSigmaProton())<1.0)) h_pTOF_HFT_1sig->Fill(trk->gPt(),pTofAvailable); 
			if ((trk->isHFTTrack()) && (trk->nHitsFit()>=20)) h_pTOF_HFT_20->Fill(trk->gPt(),pTofAvailable); 
			if ((trk->isHFTTrack()) && (fabs(trk->nSigmaProton())<1.0) && (trk->nHitsFit()>=20)) h_pTOF_HFT_1sig_20->Fill(trk->gPt(),pTofAvailable); 
            if (fabs(trk->nSigmaProton())<1.0) h_pTOF_1sig->Fill(trk->gPt(),pTofAvailable); 
        }

        h_pinsigma->Fill(momentum.mag(),trk->nSigmaPion());
        h_knsigma->Fill(momentum.mag(),trk->nSigmaKaon());
        h_pnsigma->Fill(momentum.mag(),trk->nSigmaProton());

        h_dedx->Fill(momentum.mag(),trk->dEdx());
   
  } // .. end tracks loop


      
    double duration = (double) (std::clock() - start1) / (double) CLOCKS_PER_SEC;
  TH1F *h_time_per_event = static_cast<TH1F*>(mOutList->FindObject("h_time_per_event"));
  h_time_per_event->Fill(duration);


  return kStOK;
}

int StPicoDpmAnaMaker::createQA(){
       //int const currentRun = mPicoHFEvent->runId();
       //if(currentRun != mRunNumber)
     //  {
       //mRunNumber = currentRun;
/*  
    mRefmultCorrUtil->init(mPicoDst->event()->runId());
      if (!mRefmultCorrUtil){
         LOG_WARN << " No mGRefMultCorrUtil! Skip! " << endl;
         return kStWarn;
      }
      if (mRefmultCorrUtil->isBadRun(mPicoDst->event()->runId())) return kStOK;
      
      mRefmultCorrUtil->setVzForWeight(6, -6.0, 6.0);
      mRefmultCorrUtil->readScaleForWeight("StRoot/StRefMultCorr/macros/weight_grefmult_vpd30_vpd5_Run14_P16id.txt");

      mRefmultCorrUtil->initEvent(mPicoDst->event()->grefMult(), mPrimVtx.z(), mPicoDst->event()->ZDCx()) ;

       int const centrality = mRefmultCorrUtil->getCentralityBin9();
       const double reweight = mRefmultCorrUtil->getWeight();
       const double refmultCor = mRefmultCorrUtil->getRefMultCorr();
*/  
     //mHists->addCent(refmultCor, centrality, reweight, pVtx.z());
       UInt_t nTracks = mPicoDst->numberOfTracks();
	int const centrality = 0;
       
	for (unsigned short iTrack = 0; iTrack < nTracks; ++iTrack)
       {
          StPicoTrack const* trk = mPicoDst->track(iTrack);
          if (!trk) continue;
          StPhysicalHelixD helix = trk->helix(mPicoDst->event()->bField());
          float dca = float(helix.geometricSignedDistance(mPrimVtx));
          StThreeVectorF momentum = trk->gMom(mPrimVtx, mPicoDst->event()->bField());

         // if (!isGoodQaTrack(trk, momentum, dca)) continue; pt, nhits, pseudorap
		if (!(trk->gPt()>0.3)) continue;
		if (!(trk->nHitsFit()>20)) continue;
        if (!(fabs(momentum.pseudoRapidity()) < 1.0)) continue;

          StThreeVectorF dcaPoint = helix.at(helix.pathLength(mPrimVtx.x(), mPrimVtx.y()));
          float dcaZ = dcaPoint.z() - mPrimVtx.z();
          double dcaXy = helix.geometricSignedDistance(mPrimVtx.x(), mPrimVtx.y());

          bool tpcPion = false;
          bool tpcKaon = false;
          bool tpcProton = false;
		  if(trk->nSigmaPion()<3.0) tpcPion = true;
		  if(trk->nSigmaKaon()<3.0) tpcKaon = true;
	      if(trk->nSigmaProton()<3.0) tpcProton = true;
          float hBeta = mHFCuts->getTofBetaBase(trk);
          bool hTofAvailable = !isnan(hBeta) && hBeta > 0;

          bool tofPion = false;
          bool tofKaon = false;
          bool tofProton = false;

	      if(fabs(1./hBeta - sqrt(1+M_PION_PLUS*M_PION_PLUS/(momentum.mag()*momentum.mag())))<=0.06) tofPion = true;
          if(fabs(1./hBeta - sqrt(1+M_KAON_PLUS*M_KAON_PLUS/(momentum.mag()*momentum.mag())))<=0.05) tofKaon = true;
		  if(fabs(1./hBeta - sqrt(1+M_PROTON*M_PROTON/(momentum.mag()*momentum.mag())))) tofProton = true;

          bool goodPion = (hTofAvailable && tofPion && tpcPion) || (!hTofAvailable && tpcPion);//Always require TPC
          bool goodKaon = (hTofAvailable && tofKaon && tpcKaon) || (!hTofAvailable && tpcKaon);
          bool goodProton = (hTofAvailable && tofProton && tpcProton) || (!hTofAvailable && tpcProton);

          if (trk  && fabs(dca) < 1.5 && trk->isHFTTrack() && (goodPion || goodKaon || goodProton)){
             addDcaPtCent(dca, dcaXy, dcaZ, goodPion, goodKaon, goodProton, momentum.perp(), centrality, momentum.pseudoRapidity(), momentum.phi(), mPrimVtx.z()); //add Dca distribution
          }
          if (trk  && fabs(dca) < 1.5 && (goodPion || goodKaon || goodProton)){
             //std::cout<<"1: "<<goodPion<<" "<< goodKaon<<" "<<  goodProton<<" "<<  momentum.perp()<<" "<<  centrality<<" "<<  momentum.pseudoRapidity()<<" "<<  momentum.phi()<<" "<<  mPrimVtx.z()<<std::endl;
             addTpcDenom1(goodPion, goodKaon, goodProton, momentum.perp(), centrality, momentum.pseudoRapidity(), momentum.phi(), mPrimVtx.z()); //Dca cut on 1.5cm, add Tpc Denominator
          }
          if (trk && fabs(dca) < 1.5 && trk->isHFTTrack() && (goodPion || goodKaon || goodProton) && fabs(dcaXy) < 1. && fabs(dcaZ) < 1.){
             addHFTNumer1(goodPion, goodKaon, goodProton, momentum.perp(), centrality,  momentum.pseudoRapidity(), momentum.phi(), mPrimVtx.z()); //Dca cut on 1.5cm, add HFT Numerator
          }
       } // .. end tracks loop
   return 0;
}

// _________________________________________________________
int StPicoDpmAnaMaker::createCandidates() {

  for (unsigned short idxPion1 = 0; idxPion1 < mIdxPicoPions.size(); ++idxPion1) {
    StPicoTrack const *pion1 = mPicoDst->track(mIdxPicoPions[idxPion1]);
//     if( !mHFCuts->isHybridTOFHadron(pion1, mHFCuts->getTofBetaBase(pion1), StHFCuts::kPion) ) continue;
      for (unsigned short idxKaon = 0; idxKaon < mIdxPicoKaons.size(); ++idxKaon) {
        StPicoTrack const *kaon = mPicoDst->track(mIdxPicoKaons[idxKaon]);
        // -- Kaon selection
        // -- TOF
//        if( !mHFCuts->isHybridTOFHadron(kaon, mHFCuts->getTofBetaBase(kaon), StHFCuts::kKaon) ) continue;
        if (mIdxPicoKaons[idxKaon] == mIdxPicoPions[idxPion1]) continue;
        // -- Making pair
        StHFPair pair(pion1,kaon,mHFCuts->getHypotheticalMass(StHFCuts::kPion),mHFCuts->getHypotheticalMass(StHFCuts::kKaon), mIdxPicoPions[idxPion1],mIdxPicoKaons[idxKaon], mPrimVtx, mBField, kFALSE);

//       if (!mHFCuts->isGoodSecondaryVertexPair(triplet)) continue; LK 200717
      if (!mHFCuts->isClosePair(pair)) continue;
      
       mPicoHFEvent->addHFSecondaryVertexPair(&pair);

      }  // for (unsigned short idxKaon = 0; idxKaon < mIdxPicoKaons.size(); ++idxKaon)
  } // for (unsigned short idxPion1 = 0; idxPion1 < mIdxPicoPions.size(); ++idxPion1)


 return kStOK;
}

// _________________________________________________________
int StPicoDpmAnaMaker::analyzeCandidates() {
  // --- Analyze previously constructed candidates and output to ntuple
  // -- Decay channel1
  TClonesArray const * aCandidates= mPicoHFEvent->aHFSecondaryVertices();
  if( mPicoHFEvent->nHFSecondaryVertices() > 0 ){
    for (unsigned int idx = 0; idx <  mPicoHFEvent->nHFSecondaryVertices(); ++idx) {
      
      StHFPair const* pair = static_cast<StHFPair*>(aCandidates->At(idx));
      StPicoTrack const* pion1 = mPicoDst->track(pair->particle1Idx());
      StPicoTrack const* kaon = mPicoDst->track(pair->particle2Idx());
      
      // Greates distance between tracks
      float const dcaDaughters = pair->dcaDaughters();
      float dcaMax = dcaDaughters;
      
      //TOF ---
      float kaonBetaBase = -1;
      float pion1BetaBase = -1;
      kaonBetaBase = mHFCuts->getTofBetaBase(kaon);
      pion1BetaBase = mHFCuts->getTofBetaBase(pion1);

      float kaonTOFinvbeta = fabs(1. / mHFCuts->getTofBetaBase(kaon) - sqrt(1+M_KAON_PLUS*M_KAON_PLUS/(kaon->gMom(mPrimVtx,mBField).mag()*kaon->gMom(mPrimVtx,mBField).mag())));
      float pion1TOFinvbeta = fabs(1. / mHFCuts->getTofBetaBase(pion1) - sqrt(1+M_PION_PLUS*M_PION_PLUS/(pion1->gMom(mPrimVtx,mBField).mag()*pion1->gMom(mPrimVtx,mBField).mag())));

      // -- Flag D0 and background
      float flag = -99.;
      if( kaon->charge()<0 && pion1->charge()>0 ) flag=0.; // -+
      if( kaon->charge()>0 && pion1->charge()<0 ) flag=1.; // +-

      if( kaon->charge()<0 && pion1->charge()<0) flag=4.; // --
      if( kaon->charge()>0 && pion1->charge()>0) flag=5.; // ++

      // getting centrality
      /*int const currentRun = mPicoHFEvent->runId();

      if(currentRun != mRunNumber)
      {
          // init a new run
          mRunNumber = currentRun;
          mRefmultCorrUtil->init(mRunNumber);
          mRefmultCorrUtil->setVzForWeight(6, -6.0, 6.0);
          mRefmultCorrUtil->readScaleForWeight("StRoot/StRefMultCorr/macros/weight_grefmult_vpd30_vpd5_Run14.txt");
          for(Int_t i=0;i<6;i++){
            mRefmultCorrUtil->get(i, 0);
          }
      }

      mRefmultCorrUtil->initEvent(mPicoDst->event()->grefMult(), mPrimVtx.z(), mPicoDst->event()->ZDCx()) ;
      int const centrality = mRefmultCorrUtil->getCentralityBin9() ;*/
/*
      mRefmultCorrUtil->init(mPicoDst->event()->runId());
      if (!mRefmultCorrUtil){
         LOG_WARN << " No mGRefMultCorrUtil! Skip! " << endl;
         return kStWarn;
      }
      if (mRefmultCorrUtil->isBadRun(mPicoDst->event()->runId())) return kStOK;
      
      mRefmultCorrUtil->setVzForWeight(6, -6.0, 6.0);
      mRefmultCorrUtil->readScaleForWeight("StRoot/StRefMultCorr/macros/weight_grefmult_vpd30_vpd5_Run14_P16id.txt");

      mRefmultCorrUtil->initEvent(mPicoDst->event()->grefMult(), mPrimVtx.z(), mPicoDst->event()->ZDCx()) ;

       int const centrality = mRefmultCorrUtil->getCentralityBin9();
       const double reweight = mRefmultCorrUtil->getWeight();
       const double refmultCor = mRefmultCorrUtil->getRefMultCorr();
*/
        int const centrality = 1;
       const double reweight = 1;
       const double refmultCor = 1;

      int ii=0;
      float ntVar[41];
      // ---
      // Saving to NTUPLE
     // float ref = (float)(mPicoHFEvent->grefMult()); 
     // float globalTracks = (float)(mPicoHFEvent->numberOfGlobalTracks());
      
      ntVar[ii++] = mPicoDst->event()->refMult();
      ntVar[ii++] = mPicoHFEvent->runId();
      ntVar[ii++] = mPicoHFEvent->eventId();
      ntVar[ii++] = pion1->gMom(mPrimVtx,mBField).phi();
      ntVar[ii++] = pion1->gMom(mPrimVtx,mBField).pseudoRapidity();
      ntVar[ii++] = pion1->gPt();
      ntVar[ii++] = pair->particle1Dca();
      ntVar[ii++] = pion1->dEdx();
      ntVar[ii++] = pion1->nSigmaPion();
      ntVar[ii++] = pion1->nHitsFit();
      ntVar[ii++] = pion1->nHitsDedx();
      ntVar[ii++] = pion1TOFinvbeta;
      ntVar[ii++] = pion1BetaBase;

      ntVar[ii++] = mPicoHFEvent->runId();
      ntVar[ii++] = mPicoHFEvent->eventId();
      ntVar[ii++] = kaon->gMom(mPrimVtx,mBField).phi();
      ntVar[ii++] = kaon->gMom(mPrimVtx,mBField).pseudoRapidity();
      ntVar[ii++] = kaon->gPt();
      ntVar[ii++] = pair->particle2Dca();
      ntVar[ii++] = kaon->dEdx();
      ntVar[ii++] = kaon->nSigmaKaon();
      ntVar[ii++] = kaon->nHitsFit();
      ntVar[ii++] = kaon->nHitsDedx();
      ntVar[ii++] = kaonTOFinvbeta;
      ntVar[ii++] = kaonBetaBase;

      ntVar[ii++] = dcaMax;

      ntVar[ii++] = flag;
      ntVar[ii++] = mPrimVtx.z();
      ntVar[ii++] = pair->pointingAngle();
      ntVar[ii++] = pair->decayLength();
      ntVar[ii++] = pair->phi();
      ntVar[ii++] = pair->eta();
      ntVar[ii++] = sqrt(pow(pair->px(),2.0)+pow(pair->py(),2.0));
      ntVar[ii++] = pair->m();
      if ((flag == 0) || (flag == 1)) {
          ntVar[ii++] = 0; //D_mass_LS
          ntVar[ii++] = pair->m();//D_mass_US
      } else { 
          ntVar[ii++] = pair->m(); //D_mass_LS
          ntVar[ii++] = 0;//D_mass_US
      }
      
      ntVar[ii++] = 0; // triplet->dV0Max();
      ntVar[ii++] = centrality;
      ntVar[ii++] = 0; //mPicoDst->event()->refMult();
      ntVar[ii++] = refmultCor;
      ntVar[ii++] = reweight;
      ntp_DMeson->Fill(ntVar);
    } // for (unsigned int idx = 0; idx <  mPicoHFEvent->nHFSecondaryVertices(); ++idx) {
  }
 return kStOK;
}

// _________________________________________________________
bool StPicoDpmAnaMaker::isHadron(StPicoTrack const * const trk, int pidFlag) const {
  // -- good hadron
  return (mHFCuts->isGoodTrack(trk) && mHFCuts->isTPCHadron(trk, pidFlag));
}

// _________________________________________________________
bool StPicoDpmAnaMaker::isPion(StPicoTrack const * const trk) const {
  // -- good pion
   StThreeVectorF t = trk->pMom();
   if (fabs(t.pseudoRapidity()) > 1.) return false; //pridano fabs 1212
   if (!mHFCuts->isHybridTOFHadron(trk, mHFCuts->getTofBetaBase(trk), StHFCuts::kPion) ) return false;
   if (!mHFCuts->cutMinDcaToPrimVertex(trk, StPicoCutsBase::kPion)) return false;
  return (mHFCuts->isGoodTrack(trk) && mHFCuts->isTPCHadron(trk, StPicoCutsBase::kPion));
}

// _________________________________________________________
bool StPicoDpmAnaMaker::isKaon(StPicoTrack const * const trk) const {
  // -- good kaon
  StThreeVectorF t = trk->pMom();
  if (fabs(t.pseudoRapidity()) > 1.) return false;//pridano fabs 1212
  if (!mHFCuts->isHybridTOFHadron(trk, mHFCuts->getTofBetaBase(trk), StHFCuts::kKaon) ) return false;
  if (!mHFCuts->cutMinDcaToPrimVertex(trk, StPicoCutsBase::kKaon)) return false;
  return (mHFCuts->isGoodTrack(trk) && mHFCuts->isTPCHadron(trk, StPicoCutsBase::kKaon));
} 

// _________________________________________________________
bool StPicoDpmAnaMaker::isProton(StPicoTrack const * const trk) const {
  // -- good proton
  return (mHFCuts->isGoodTrack(trk) && mHFCuts->isTPCHadron(trk, StPicoCutsBase::kProton));
}

double StPicoDpmAnaMaker::DCA(StPicoTrack const * const trk, StThreeVectorF const & vtx) const {
  // -- particle DCA
  StPhysicalHelixD pHelix = trk->helix(mPicoDst->event()->bField());
  pHelix.moveOrigin(pHelix.pathLength(vtx));
  return ((pHelix.origin() - vtx).mag());
}


bool StPicoDpmAnaMaker::isCloseTracks(StPicoTrack const * const trk1, StPicoTrack const * const trk2, StThreeVectorF const & vtx, float bField) const {
  
  StPhysicalHelixD p1Helix = trk1->helix(mPicoDst->event()->bField());
  StPhysicalHelixD p2Helix = trk2->helix(mPicoDst->event()->bField());
  p1Helix.moveOrigin(p1Helix.pathLength(vtx));
  p2Helix.moveOrigin(p2Helix.pathLength(vtx));
  if( ( p1Helix.origin()-vtx ).mag()>0.2 || ( p2Helix.origin()-vtx ).mag()>0.2 ) return false;

  //Requires loading constants
  StThreeVectorF const p1Mom = p1Helix.momentum(bField * kilogauss);
  StThreeVectorF const p2Mom = p2Helix.momentum(bField * kilogauss);
  StPhysicalHelixD const p1StraightLine(p1Mom, p1Helix.origin(), 0, trk1->charge());
  StPhysicalHelixD const p2StraightLine(p2Mom, p2Helix.origin(), 0, trk2->charge());
  //DCA
  pair<double, double> const ss = p1StraightLine.pathLengths(p2StraightLine);
  StThreeVectorF const p1AtDcaToP2 = p1StraightLine.at(ss.first);
  StThreeVectorF const p2AtDcaToP1 = p2StraightLine.at(ss.second);
  float const dca = (p1AtDcaToP2-p2AtDcaToP1).mag();
  if(dca > 0.050) return false; //kubo : 0.009, in cm
// -- good pair
  return true;
}

//-----------------------------------------------------------------------------

void StPicoDpmAnaMaker::histoInit(TString fileBaseName, bool fillQaHists){
  TString m_ParticleName[m_nParticles] = {"Pion", "Kaon", "Proton"};

   float m_EtaEdgeDca[m_nEtasDca+1] = {-1.0, -0.6, -0.2, 0.2, 0.6, 1.0}; //replace bottom!!!
   float m_PhiEdgeDca[m_nPhisDca + 1] = {-3.14159, -2.80359, -2.17527, -1.54696, -0.918637, -0.290319, 0.338, 0.966319, 1.59464, 2.22296, 2.85127, 3.14159};
   float m_VzEdgeDca[m_nVzsDca + 1] = { -6.0, -3.0, 0, 3.0, 6.0};//replace bottom!!!
   float m_CentEdgeDca[m_nCentsDca + 1] = { -0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5};
   float m_PtEdgeDca[m_nPtsDca + 1] = {0.3, 0.4, 0.5, 0.6,  0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.5, 4.0, 6.0, 12.0};
   float m_EtaEdgeRatio[m_nEtasRatio + 1] = { -1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4 , 0.6, 0.8, 1.0}; //replace bottom!!!
   float m_PhiEdgeRatio[m_nPhisRatio + 1] = { -3.14159, -2.80359, -2.17527, -1.54696, -0.918637, -0.290319, 0.338, 0.966319, 1.59464, 2.22296, 2.85127, 3.14159};//replace bottom!!!
   float m_VzEdgeRatio[m_nVzsRatio + 1] = { -6.0, -4.0, -2.0, 0, 2.0, 4.0, 6.0};//replace bottom!!!
   float m_CentEdgeRatio[m_nCentsRatio + 1] = { -1.5, -0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5};
   float m_PtEdgeRatio[m_nPtsRatio + 1] =
   {
      0.3, 0.4, 0.5, 0.6 , 0.7 , 0.8 , 0.9 ,
      1. , 1.1 , 1.2 , 1.3 , 1.4 , 1.5 , 1.6 , 1.7 , 1.8 , 1.9 ,
      2. , 2.2 , 2.4 , 2.6 , 2.8 , 3.0 ,
      3.4 , 3.8 , 4.2 , 4.6 , 5.0 ,  5.5 ,
      6. , 6.5 , 7.0 , 8.0 , 9.0 , 10. , 11,  12.0
   };
  float m_DcaEdgeDca[m_nDcasDca + 1] =
   {
     -1 , -0.96 , -0.92 , -0.88 , -0.84 , -0.8 , -0.76 , -0.72 , -0.68 , -0.64 , -0.6 , -0.56 , -0.52 , -0.48 , -0.44 , -0.4 , -0.36 , -0.32 , -0.28 , -0.24 , -0.2 , -0.16 , -0.12 ,  -0.08,
     -0.078 , -0.075 , -0.072 , -0.069 , -0.066 , -0.063 , -0.06 , -0.057 , -0.054 , -0.051 , -0.048 , -0.045 , -0.042 , -0.039 , -0.036 , -0.033 , -0.03 , -0.027 , -0.024 , -0.021 , -0.018 , -0.015 , -0.012 ,
      -0.01 , -0.0096 , -0.0092 , -0.0088 , -0.0084 , -0.008 , -0.0076 , -0.0072 , -0.0068 , -0.0064 , -0.006 , -0.0056 , -0.0052 , -0.0048 , -0.0044 , -0.004 , -0.0036 , -0.0032 , -0.0028 , -0.0024 , -0.002 , -0.0016 , -0.0012 , -0.0008 , -0.0004 , 0 , 0.0004 , 0.0008 , 0.0012 , 0.0016 , 0.002 , 0.0024 , 0.0028 , 0.0032 , 0.0036 , 0.004 , 0.0044 , 0.0048 , 0.0052 , 0.0056 , 0.006 , 0.0064 , 0.0068 , 0.0072 , 0.0076 , 0.008 , 0.0084 , 0.0088 , 0.0092 , 0.0096 , 0.01 ,
      0.012 , 0.015 , 0.018 , 0.021 , 0.024 , 0.027 , 0.03 , 0.033 , 0.036 , 0.039 , 0.042 , 0.045 , 0.048 , 0.051 , 0.054 , 0.057 , 0.06 , 0.063 , 0.066 , 0.069 , 0.072 , 0.075 , 0.078 ,
      0.08 , 0.12 , 0.16 , 0.2 , 0.24 , 0.28 , 0.32 , 0.36 , 0.4 , 0.44 , 0.48 , 0.52 , 0.56 , 0.6 , 0.64 , 0.68 , 0.72 , 0.76 , 0.8 , 0.84 , 0.88 , 0.92 , 0.96 , 1
   };



   //set in private
   // for(int temp = 0;temp<m_nParticles;temp++) m_ParticleName[temp]=temp_ParticleName[temp];
	//for(int temp = 0;temp<m_nEtasDca+1;temp++) m_EtaEdgeDca[temp]=temp_EtaEdgeDca[temp];
	//for(int temp2 = 0;temp2<m_nPhisDca+1;temp2++) m_PhiEdgeDca[temp2]=temp_PhiEdgeDca[temp2];
	/*for(int temp = 0;temp<m_nVzsDca+1;temp++) m_VzEdgeDca[temp]=temp_VzEdgeDca[temp];
	for(int temp = 0;temp<m_nCentsDca+1;temp++) m_CentEdgeDca[temp]=temp_CentEdgeDca[temp];
	for(int temp = 0;temp<m_nPtsDca+1;temp++) m_PtEdgeDca[temp]=temp_PtEdgeDca[temp];
	for(int temp = 0;temp<m_nEtasRatio+1;temp++) m_EtaEdgeRatio[temp]=temp_EtaEdgeRatio[temp];
	for(int temp = 0;temp<m_nPhisRatio+1;temp++) m_PhiEdgeRatio[temp]=temp_PhiEdgeRatio[temp];
	for(int temp = 0;temp<m_nVzsRatio+1;temp++) m_VzEdgeRatio[temp]=temp_VzEdgeRatio[temp];
	for(int temp = 0;temp<m_nCentsRatio+1;temp++) m_CentEdgeRatio[temp]=temp_CentEdgeRatio[temp];
	for(int temp = 0;temp<m_nPtsRatio+1;temp++) m_PtEdgeRatio[temp]=temp_PtEdgeRatio[temp];
	for(int temp = 0;temp<m_nDcasDca+1;temp++) m_DcaEdgeDca[temp]=temp_DcaEdgeDca[temp];*/

   mFillQaHists = fillQaHists;
   mOutFile = new TFile(fileBaseName+".hists.root", "RECREATE");
   //mOutFile = new TFile(Form("%s.hists.root", fileBaseName.Data()), "RECREATE");
/*
   for (int iParticle = 0; iParticle < m_nParticles; iParticle++){
      for (int iEta = 0; iEta < m_nEtasDca; iEta++){
         for (int iVz = 0; iVz < m_nVzsDca; iVz++){
            for (int iCent = 0; iCent < m_nCentsDca; iCent++){
               mh3DcaXyZPtCentPartEtaVzPhi[iParticle][iEta][iVz][iCent] = NULL;
            }
         }
      }
   }

   for (int iParticle = 0; iParticle < m_nParticles; iParticle++){
      for (int iEta = 0; iEta < m_nEtasRatio; iEta++){
         for (int iVz = 0; iVz < m_nVzsRatio; iVz++){
            for (int iPhi = 0; iPhi < m_nPhisRatio; iPhi++){
               mh2Tpc1PtCentPartEtaVzPhi[iParticle][iEta][iVz][iPhi] = NULL;
               mh2HFT1PtCentPartEtaVzPhi[iParticle][iEta][iVz][iPhi] = NULL;
            }
         }
      }
   }*/

  

   TH1::SetDefaultSumw2();
   if (!mFillQaHists) return;

   mh1Cent         = new TH1F("mh1Cent", "EventsVsCentrality;cent;Counts", 10, -1.5, 8.5);
   mh1CentWg         = new TH1F("mh1CentWg", "EventsVsCentrality;cent;Counts", 10, -1.5, 8.5);
   mh1gRefmultCor  = new TH1F("mh1gRefmultCor", "gRefmultCor;gRefmult;Counts", 700, 0, 700);
   mh1gRefmultCorWg  = new TH1F("mh1gRefmultCorWg", "gRefmultCorWg;gRefmultCorWg;Counts", 700, 0, 700);
   mh2CentVz         = new TH2F("mh2CentVz", "CentralityVsVz;cent;Vz", 10, -1.5, 8.5, 200, -10, 10);
   mh2CentVzWg = new TH2F("mh2CentVzWg", "CentralityVsVzWg;cent;Vz", 10, -1.5, 8.5, 200, -10, 10);

   //Add some HFT ratio plots
   mh2Tpc1PtCent  = new TH2F("mh2Tpc1PtCent", "Tpc tacks;p_{T}(GeV/c);cent", 120, 0, 12, 10, -1.5, 8.5); //Dca 1.5cm
   mh2HFT1PtCent  = new TH2F("mh2HFT1PtCent", "HFT tacks;p_{T}(GeV/c);cent", 120, 0, 12, 10, -1.5, 8.5); //Dca 1.5cm
   mh2Tpc1PhiVz  = new TH2F("mh2Tpc1PhiVz", "Tpc tacks;#Phi;Vz", 100, -3.1415, 3.1415, 20, -10, 10); //Dca 1.5cm
   mh2HFT1PhiVz  = new TH2F("mh2HFT1PhiVz", "HFT tacks;#Phi;Vz", 100, -3.1415, 3.1415, 20, -10, 10); //Dca 1.5cm

   for (int iParticle = 0; iParticle < m_nParticles; iParticle++){
      for (int iEta = 0; iEta < m_nEtasRatio; iEta++){
         for (int iVz = 0; iVz < m_nVzsRatio; iVz++){
            for (int iPhi = 0; iPhi < m_nPhisRatio; iPhi++){
             
               mh2Tpc1PtCentPartEtaVzPhi[iParticle][iEta][iVz][iPhi]  = new TH2F(Form("mh2Tpc1PtCentPartEtaVzPhi_%d_%d_%d_%d", iParticle, iEta, iVz, iPhi), "mh2Tpc1PtCent_"+m_ParticleName[iParticle]+Form("_Eta%2.1f_Vz%2.1f_Phi%2.1f;p_{T}(GeV/c);cent", m_EtaEdgeRatio[iEta], m_VzEdgeRatio[iVz], m_PhiEdgeRatio[iPhi]), m_nPtsRatio, m_PtEdgeRatio, m_nCentsRatio, m_CentEdgeRatio); //Dca 1.cm
               mh2HFT1PtCentPartEtaVzPhi[iParticle][iEta][iVz][iPhi]  = new TH2F(Form("mh2HFT1PtCentPartEtaVzPhi_%d_%d_%d_%d", iParticle, iEta, iVz, iPhi), "mh2HFT1PtCent_"+m_ParticleName[iParticle]+Form("_Eta%2.1f_Vz%2.1f_Phi%2.1f;p_{T}(GeV/c);cent", m_EtaEdgeRatio[iEta], m_VzEdgeRatio[iVz], m_PhiEdgeRatio[iPhi]), m_nPtsRatio, m_PtEdgeRatio, m_nCentsRatio, m_CentEdgeRatio); //Dca 1.cm
            }
         }
      }
   }

   // Add some Dca, resolution
   for (int iParticle = 0; iParticle < m_nParticles; iParticle++){
      for (int iEta = 0; iEta < m_nEtasDca; iEta++){
         for (int iVz = 0; iVz < m_nVzsDca; iVz++){
            for (int iCent = 0; iCent < m_nCentsDca; iCent++){
               
   	    	   mh3DcaXyZPtCentPartEtaVzPhi[iParticle][iEta][iVz][iCent]  = new TH3F(Form("mh3DcaXyZPtCentPartEtaVzPhi_%d_%d_%d_%d", iParticle, iEta, iVz, iCent),"mh3DcaXyZPt_"+m_ParticleName[iParticle]+Form("_Eta%2.1f_Vz%2.1f_Cent%2.1f;p_{T}(GeV/c);DcaXy(cm);DcaZ(cm)", m_EtaEdgeDca[iEta], m_VzEdgeDca[iVz], m_CentEdgeDca[iCent]), m_nPtsDca, m_PtEdgeDca, m_nDcasDca, m_DcaEdgeDca, m_nDcasDca, m_DcaEdgeDca); //Dca 1.cm
            }
         }
      }
   }

   mh3DcaPtCent  = new TH3F("mh3DcaPtCent", "mh3DcaPtCent;p_{T}(GeV/c);cent;Dca(cm)", 120, 0, 12, 10, -1.5, 8.5, 1000, -1, 1); //Dca 1.cm
   mh3DcaXyPtCent  = new TH3F("mh3DcaXyPtCent", "mh3DcaXyPtCent;p_{T}(GeV/c);cent;DcaXy(cm)", 120, 0, 12, 10, -1.5, 8.5, 1000, -1, 1); //Dca 1.cm
   mh3DcaZPtCent  = new TH3F("mh3DcaZPtCent", "mh3DcaZPtCent;p_{T}(GeV/c);cent;DcaZ(cm)", 120, 0, 12, 10, -1.5, 8.5, 1000, -1, 1); //Dca 1.cm

}

//-----------------------------------------------------------------------
void StPicoDpmAnaMaker::addTpcDenom1(bool IsPion, bool IsKaon, bool IsProton, float pt, int centrality, float Eta, float Phi, float Vz){
   int EtaIndex = getEtaIndexRatio(Eta);
   int PhiIndex = getPhiIndexRatio(Phi);
   int VzIndex = getVzIndexRatio(Vz);
      if(EtaIndex == -1) return;
   if(PhiIndex == -1) return;
   if(VzIndex == -1) return;
   //std::cout<<"2: "<<IsPion<<" "<<IsKaon<<" "<<IsProton<<" "<<pt<<" "<<centrality<<" "<<Eta<<" "<<Phi<<" "<<Vz<<" "<<EtaIndex<<" "<<PhiIndex<<" "<<VzIndex<<std::endl;
   
   if (IsPion){
      mh2Tpc1PtCentPartEtaVzPhi[0][EtaIndex][VzIndex][PhiIndex]->Fill(pt, centrality);
      //if(mh2Tpc1PtCentPartEtaVzPhi[0][EtaIndex][VzIndex][PhiIndex]) std::cout<<"true"<<<<std::endl;
      //std::cout<<pt<<" "<<centrality<<std::endl;
   }
   if (IsKaon){
      mh2Tpc1PtCentPartEtaVzPhi[1][EtaIndex][VzIndex][PhiIndex]->Fill(pt, centrality);
   }
   if (IsProton){
      mh2Tpc1PtCentPartEtaVzPhi[2][EtaIndex][VzIndex][PhiIndex]->Fill(pt, centrality);
   }
   mh2Tpc1PtCent->Fill(pt, centrality);
   if (fabs(Eta) < 0.1 && pt > 3.0) mh2Tpc1PhiVz->Fill(Phi, Vz);
}
//-----------------------------------------------------------------------
void StPicoDpmAnaMaker::addHFTNumer1(bool IsPion, bool IsKaon, bool IsProton, float pt, int centrality, float Eta, float Phi, float Vz){
   int EtaIndex = getEtaIndexRatio(Eta);
   int PhiIndex = getPhiIndexRatio(Phi);
   int VzIndex = getVzIndexRatio(Vz);
   if(EtaIndex == -1) return;
   if(PhiIndex == -1) return;
   if(VzIndex == -1) return;
   if (IsPion){
      mh2HFT1PtCentPartEtaVzPhi[0][EtaIndex][VzIndex][PhiIndex]->Fill(pt, centrality);
   }
   if (IsKaon){
      mh2HFT1PtCentPartEtaVzPhi[1][EtaIndex][VzIndex][PhiIndex]->Fill(pt, centrality);
   }
   if (IsProton){
      mh2HFT1PtCentPartEtaVzPhi[2][EtaIndex][VzIndex][PhiIndex]->Fill(pt, centrality);
   }
   mh2HFT1PtCent->Fill(pt, centrality);
   if (fabs(Eta) < 0.1 && pt > 3.0) mh2HFT1PhiVz->Fill(Phi, Vz);
}
//---------------------------------------------------------------------
void StPicoDpmAnaMaker::addDcaPtCent(float dca, float dcaXy, float dcaZ, bool IsPion, bool IsKaon, bool IsProton, float pt,  int centrality, float Eta, float Phi, float Vz){
   int EtaIndex = getEtaIndexDca(Eta);
   int VzIndex = getVzIndexDca(Vz);
   if(EtaIndex == -1) return;
   if(VzIndex == -1) return;

   if (centrality < 0) return; // remove bad centrality, only keep 9 centralities
   if (IsPion){
      mh3DcaXyZPtCentPartEtaVzPhi[0][EtaIndex][VzIndex][centrality]->Fill(pt, dcaXy, dcaZ);
   }
   if (IsKaon){
      mh3DcaXyZPtCentPartEtaVzPhi[1][EtaIndex][VzIndex][centrality]->Fill(pt, dcaXy, dcaZ);
   }
   if (IsProton){
      mh3DcaXyZPtCentPartEtaVzPhi[2][EtaIndex][VzIndex][centrality]->Fill(pt, dcaXy, dcaZ);
   }
   mh3DcaPtCent->Fill(pt, centrality, dca);
   mh3DcaXyPtCent->Fill(pt, centrality, dcaXy);
   mh3DcaZPtCent->Fill(pt, centrality, dcaZ);
}
//---------------------------------------------------------------------
int StPicoDpmAnaMaker::getEtaIndexDca(float Eta){
   float EtaEdgeDca[m_nEtasDca+1] = {-1.0, -0.6, -0.2, 0.2, 0.6, 1.0};
   for (int i = 0; i < m_nEtasDca; i++){
	 if ((Eta >= EtaEdgeDca[i]) && (Eta < EtaEdgeDca[i + 1]))
         return i;
   }
   //std::cout<<"SOMETHING WENT TERRIBRU WONG"<<std::endl;
   //return m_nEtasDca -1;
   return -1;
}

//---------------------------------------------------------------------
int StPicoDpmAnaMaker::getVzIndexDca(float Vz){
  float VzEdgeDca[m_nVzsDca + 1] = { -6.0, -3.0, 0, 3.0, 6.0};
   for (int i = 0; i < m_nVzsDca; i++){
      if ((Vz >= VzEdgeDca[i]) && (Vz < VzEdgeDca[i + 1]))
         return i;
   }
//std::cout<<"SOMETHING WENT TERRIBRU WONG"<<std::endl;
   //return m_nVzsDca - 1;
   return -1;
}
//---------------------------------------------------------------------
int StPicoDpmAnaMaker::getEtaIndexRatio(float Eta){
  float EtaEdgeRatio[m_nEtasRatio + 1] = { -1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4 , 0.6, 0.8, 1.0}; 
   for (int i = 0; i < m_nEtasRatio; i++){
      if ((Eta >= EtaEdgeRatio[i]) && (Eta < EtaEdgeRatio[i + 1]))
         return i;
   }
//std::cout<<"SOMETHING WENT TERRIBRU WONG"<<std::endl;
   //return m_nEtasRatio - 1;
   return -1;
}
//---------------------------------------------------------------------
int StPicoDpmAnaMaker::getPhiIndexRatio(float Phi){
  float PhiEdgeRatio[m_nPhisRatio + 1] = { -3.14159, -2.80359, -2.17527, -1.54696, -0.918637, -0.290319, 0.338, 0.966319, 1.59464, 2.22296, 2.85127, 3.14159};
   for (int i = 0; i < m_nPhisRatio; i++){
      if ((Phi >= PhiEdgeRatio[i]) && (Phi < PhiEdgeRatio[i + 1]))
         return i;
   }
//std::cout<<"SOMETHING WENT TERRIBRU WONG"<<std::endl;
  // return m_nPhisRatio - 1;
   return -1;
}
//---------------------------------------------------------------------
int StPicoDpmAnaMaker::getVzIndexRatio(float Vz){
  float VzEdgeRatio[m_nVzsRatio + 1] = { -6.0, -4.0, -2.0, 0, 2.0, 4.0, 6.0};
   for (int i = 0; i < m_nVzsRatio; i++) {
      if ((Vz >= VzEdgeRatio[i]) && (Vz < VzEdgeRatio[i + 1]))
         return i;
   }
//std::cout<<"SOMETHING WENT TERRIBLE WRONG"<<std::endl;
  // return m_nVzsRatio - 1;
  return -1;
}

void StPicoDpmAnaMaker::addCent(const double refmultCor, int centrality, const double reweight, const float vz)
{
   mh1gRefmultCor->Fill(refmultCor);
   mh1gRefmultCorWg->Fill(refmultCor, reweight);
   mh1Cent->Fill(centrality);
   mh1CentWg->Fill(centrality, reweight);
   mh2CentVz->Fill(centrality, vz);
   mh2CentVzWg->Fill(centrality, vz, reweight);
}

//---------------------------------------------------------------------
void StPicoDpmAnaMaker::closeFile()
{
   mOutFile->cd();

   mh1Cent->Write();
   mh1CentWg->Write();
   mh1gRefmultCor->Write();
   mh1gRefmultCorWg->Write();
   mh2CentVz->Write();
   mh2CentVzWg->Write();

   //HFT ratio QA
   mh2Tpc1PtCent->Write();
   mh2Tpc1PhiVz->Write();
   mh2HFT1PhiVz->Write();
   mh2HFT1PtCent->Write();

   //HFT DCA Ratio
   for (int iParticle = 0; iParticle < m_nParticles; iParticle++)
   {
      for (int iEta = 0; iEta < m_nEtasDca; iEta++)
      {
         for (int iVz = 0; iVz < m_nVzsDca; iVz++)
         {
            for (int iCent = 0; iCent < m_nCentsDca; iCent++)
            {
               mh3DcaXyZPtCentPartEtaVzPhi[iParticle][iEta][iVz][iCent]->Write();
            }
         }
      }
   }
  // std::cout<<"tuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuu"<<m_nParticles<<" "<<m_nEtasRatio<<std::endl;

   for (int iParticle = 0; iParticle < m_nParticles; iParticle++)
   {
      for (int iEta = 0; iEta < m_nEtasRatio; iEta++)
      {
         for (int iVz = 0; iVz < m_nVzsRatio; iVz++)
         {
            for (int iPhi = 0; iPhi < m_nPhisRatio; iPhi++)
            {
               mh2Tpc1PtCentPartEtaVzPhi[iParticle][iEta][iVz][iPhi]->Write();
               mh2HFT1PtCentPartEtaVzPhi[iParticle][iEta][iVz][iPhi]->Write();
            }
         }
      }
   }

   mh3DcaPtCent->Write();
   mh3DcaXyPtCent->Write();
   mh3DcaZPtCent->Write();

   // nt->Write();
   mOutFile->Write();
   mOutFile->Close();
   //mOutFile->Delete();
}


