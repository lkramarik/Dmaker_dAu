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
StPicoDpmAnaMaker::StPicoDpmAnaMaker(char const* name, StPicoDstMaker* picoMaker, char const* outputBaseFileName, char const* inputHFListHFtree = "") :
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
     *  mOutList->Add(new TH2D("h_kdedx","h_kdedx", 500, 0, 10, 500, -10, 10));
     *  mOutList->Add(new TH2D("h_pdedx","h_pdedx", 500, 0, 10, 500, -10, 10));*/
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
    mOutList->Add(new TH1F("h_refMult", "Multiplicity", 120, 0, 120));

    ntp_DMeson_Signal = new TNtuple("ntp_signal","DMeson TreeSignal","grefMult:pi1_runId:pi1_eventId:pi1_phi:pi1_eta:pi1_pt:pi1_dca:pi1_dedx:pi1_nSigma:pi1_nHitFit:pi1_nHitdedx:pi1_TOFinvbeta:pi1_betaBase:k_runId:k_eventId:k_phi:k_eta:k_pt:k_dca:k_dedx:k_nSigma:k_nHitFit:k_nHitdedx:k_TOFinvbeta:k_betaBase:dcaDaughters:flag:primVz:D_theta:cosTheta:D_decayL:dcaD0ToPv:D_phi:D_eta:D_cosThetaStar:D_pt:D_mass:D_mass_LS:D_mass_US:centrality:refmultcorr:reweight");
    ntp_DMeson_Background = new TNtuple("ntp_background","DMeson TreeBackground","grefMult:pi1_runId:pi1_eventId:pi1_phi:pi1_eta:pi1_pt:pi1_dca:pi1_dedx:pi1_nSigma:pi1_nHitFit:pi1_nHitdedx:pi1_TOFinvbeta:pi1_betaBase:k_runId:k_eventId:k_phi:k_eta:k_pt:k_dca:k_dedx:k_nSigma:k_nHitFit:k_nHitdedx:k_TOFinvbeta:k_betaBase:dcaDaughters:flag:primVz:D_theta:cosTheta:D_decayL:dcaD0ToPv:D_phi:D_eta:D_cosThetaStar:D_pt:D_mass:D_mass_LS:D_mass_US:centrality:refmultcorr:reweight");

    mRunNumber = 0;
    return kStOK;
}

// _________________________________________________________
void StPicoDpmAnaMaker::ClearHF(Option_t *opt="") {
    return;
}

// _________________________________________________________
int StPicoDpmAnaMaker::FinishHF() {
    if( isMakerMode() != StPicoHFMaker::kWrite ){
        ntp_DMeson_Signal -> Write(ntp_DMeson_Signal->GetName(), TObject::kOverwrite);
        ntp_DMeson_Background -> Write(ntp_DMeson_Background->GetName(), TObject::kOverwrite);
    }
    mOutFile->Close();

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

    TH1F *h_refMult = static_cast<TH1F*>(mOutList->FindObject("h_refMult"));

    StThreeVectorF pVtx = mPicoDst->event()->primaryVertex();

    float multiplicity = mPicoDst->event()->refMult();
    h_refMult -> Fill(multiplicity);

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

// _________________________________________________________
int StPicoDpmAnaMaker::createCandidates() {

    for (unsigned short idxPion1 = 0; idxPion1 < mIdxPicoPions.size(); idxPion1++) {
        StPicoTrack const *pion1 = mPicoDst->track(mIdxPicoPions[idxPion1]);
        //     if( !mHFCuts->isHybridTOFHadron(pion1, mHFCuts->getTofBetaBase(pion1), StHFCuts::kPion) ) continue;
        for (unsigned short idxKaon = 0; idxKaon < mIdxPicoKaons.size(); idxKaon++) {
            StPicoTrack const *kaon = mPicoDst->track(mIdxPicoKaons[idxKaon]);
            //        if( !mHFCuts->isHybridTOFHadron(kaon, mHFCuts->getTofBetaBase(kaon), StHFCuts::kKaon) ) continue;
            if (mIdxPicoKaons[idxKaon] == mIdxPicoPions[idxPion1]) continue;
            // -- Making pair
            StHFPair pair(pion1, kaon, mHFCuts->getHypotheticalMass(StHFCuts::kPion),mHFCuts->getHypotheticalMass(StHFCuts::kKaon), mIdxPicoPions[idxPion1],mIdxPicoKaons[idxKaon], mPrimVtx, mBField, kTRUE);
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

            //TOF ---
            float kaonTOFinvbeta = -1;
            float pion1TOFinvbeta = -1;
            float kaonBetaBase = -1;
            float pion1BetaBase = -1;
            kaonBetaBase = mHFCuts->getTofBetaBase(kaon);
            pion1BetaBase = mHFCuts->getTofBetaBase(pion1);

            if(!isnan(kaonBetaBase) && kaonBetaBase > 0){
                kaonTOFinvbeta = fabs(1. / kaonBetaBase - sqrt(1+M_KAON_PLUS*M_KAON_PLUS/(kaon->gMom(mPrimVtx,mBField).mag()*kaon->gMom(mPrimVtx,mBField).mag())));
            }

            if(!isnan(pion1BetaBase) && pion1BetaBase > 0){
                pion1TOFinvbeta = fabs(1. / pion1BetaBase - sqrt(1+M_PION_PLUS*M_PION_PLUS/(pion1->gMom(mPrimVtx,mBField).mag()*pion1->gMom(mPrimVtx,mBField).mag())));
            }

            // -- Flag D0 and background
            float flag = -99.;
            if( kaon->charge()<0 && pion1->charge()>0 ) flag=0.; // -+
            if( kaon->charge()>0 && pion1->charge()<0 ) flag=1.; // +-

            if( kaon->charge()<0 && pion1->charge()<0) flag=4.; // --
            if( kaon->charge()>0 && pion1->charge()>0) flag=5.; // ++


            int const centrality = 1;
            const double reweight = 1;
            const double refmultCor = 1;

            int ii=0;
            float ntVar[42];
            // ---
            // Saving to NTUPLE
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

            ntVar[ii++] = pair->dcaDaughters();

            ntVar[ii++] = flag;
            ntVar[ii++] = mPrimVtx.z();
            ntVar[ii++] = pair->pointingAngle();
            ntVar[ii++] = cos(pair->pointingAngle());
            ntVar[ii++] = pair->decayLength();
            ntVar[ii++] = pair->DcaToPrimaryVertex();
            ntVar[ii++] = pair->phi();
            ntVar[ii++] = pair->eta();
            ntVar[ii++] = pair->cosThetaStar();

            ntVar[ii++] = sqrt(pow(pair->px(),2.0)+pow(pair->py(),2.0));
            ntVar[ii++] = pair->m();
            if ((flag == 0) || (flag == 1)) {
                ntVar[ii++] = -5; //D_mass_LS
                ntVar[ii++] = pair->m();//D_mass_US
            } else {
                ntVar[ii++] = pair->m(); //D_mass_LS
                ntVar[ii++] = -5;//D_mass_US
            }
            ntVar[ii++] = centrality;
            ntVar[ii++] = refmultCor;
            ntVar[ii++] = reweight;

            if ((flag == 0) || (flag == 1)) {
                ntp_DMeson_Signal->Fill(ntVar);
            } else {
                ntp_DMeson_Background->Fill(ntVar);
            }
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
//    if (!mHFCuts->cutMinDcaToPrimVertex(trk, StPicoCutsBase::kPion)) return false;
    return (mHFCuts->isGoodTrack(trk) && mHFCuts->isTPCHadron(trk, StPicoCutsBase::kPion));
}

// _________________________________________________________
bool StPicoDpmAnaMaker::isKaon(StPicoTrack const * const trk) const {
    // -- good kaon
    StThreeVectorF t = trk->pMom();
    if (fabs(t.pseudoRapidity()) > 1.) return false;//pridano fabs 1212
    if (!mHFCuts->isHybridTOFHadron(trk, mHFCuts->getTofBetaBase(trk), StHFCuts::kKaon) ) return false;
//    if (!mHFCuts->cutMinDcaToPrimVertex(trk, StPicoCutsBase::kKaon)) return false;
    return (mHFCuts->isGoodTrack(trk) && mHFCuts->isTPCHadron(trk, StPicoCutsBase::kKaon));
}

// _________________________________________________________
bool StPicoDpmAnaMaker::isProton(StPicoTrack const * const trk) const {
    // -- good proton
    return (mHFCuts->isGoodTrack(trk) && mHFCuts->isTPCHadron(trk, StPicoCutsBase::kProton));
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
    if(dca > 0.50) return false; //kubo : 0.009, in cm
    // -- good pair
    return true;
}
