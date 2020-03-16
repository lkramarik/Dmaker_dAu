#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoHFMaker/StHFCuts.h"
#include "phys_constants.h"
#include "StPicoQAMaker.h"
#include "StPicoKFVertexFitter/StPicoKFVertexFitter.h"

ClassImp(StPicoQAMaker)

// _________________________________________________________
StPicoQAMaker::StPicoQAMaker(char const* name, StPicoDstMaker* picoMaker, char const* outputBaseFileName) :
        StPicoHFMaker(name, picoMaker, outputBaseFileName),
        mOutFileBaseName(outputBaseFileName){
    // constructor
}

// _________________________________________________________
StPicoQAMaker::~StPicoQAMaker() {
    // destructor
}

// _________________________________________________________
int StPicoQAMaker::InitHF() {
    ifstream RunList("./StRoot/StPicoQAMaker/runs_numbers.list");

    if (RunList.is_open()) {
        int Run;
        while (RunList >> Run) {
            RunNumberVector.push_back(Run);
        }
    } else {
        cout<<"Failed to open file!"<<endl;
    }

    //___detector, centrality and statistics (No. of events, tracks...) histograms
    mOutList->Add(new TH2F("h_mh1Cent", "EventsVsCentrality;cent;CountsvsRunIndex", 10, -1.5, 8.5, RunNumberVector.size() + 1, -1, RunNumberVector.size()));
    mOutList->Add(new TH2F("h_mh1CentWg", "EventsVsCentrality;cent;CountsvsRunIndex", 10, -1.5, 8.5, RunNumberVector.size() + 1, -1, RunNumberVector.size()));
    mOutList->Add(new TH2F("h_gRefmult", "gRefmult;CountsvsRunIndex", 700, 0, 700, RunNumberVector.size() + 1, -1, RunNumberVector.size()));

    mOutList->Add(new TH2F("h_gRefmult_vs_ZDCx", "gRefmult;ZDCx", 250, 0, 250, 101, -0.5, 100.5));
    mOutList->Add(new TH2F("h_gRefmult_vs_BBCx", "gRefmult;BBCx", 1350, 0, 1350, 101, -0.5, 100.5));

    mOutList->Add(new TH1F("h_gRefmult_Vz_min6_min4", "gRefmult for Vz>-6 && Vz<=-4", 101, -0.5, 100.5));
    mOutList->Add(new TH1F("h_gRefmult_Vz_min4_min2", "gRefmult for Vz>-4 && Vz<=-2", 101, -0.5, 100.5));
    mOutList->Add(new TH1F("h_gRefmult_Vz_min2_0", "gRefmult for Vz>-2 && Vz<=0", 101, -0.5, 100.5));
    mOutList->Add(new TH1F("h_gRefmult_Vz_0_2", "gRefmult for Vz>0 && Vz<=2", 101, -0.5, 100.5));
    mOutList->Add(new TH1F("h_gRefmult_Vz_2_4", "gRefmult for Vz>2 && Vz<=4", 101, -0.5, 100.5));
    mOutList->Add(new TH1F("h_gRefmult_Vz_4_6", "gRefmult for Vz>4 && Vz<=6", 101, -0.5, 100.5));

    mOutList->Add(new TH2F("h_gRefmult_HFT", "gRefmult;CountsvsRunIndex", 700, 0, 700, RunNumberVector.size() + 1, -1, RunNumberVector.size()));
    mOutList->Add(new TH2F("h_gRefmult_vs_ZDCx_HFT", "gRefmult;ZDCx", 250, 0, 250, 101, -0.5, 100.5));
    mOutList->Add(new TH2F("h_gRefmult_vs_BBCx_HFT", "gRefmult;BBCx", 1350, 0, 1350, 101, -0.5, 100.5));

    mOutList->Add(new TH1F("h_gRefmult_Vz_min6_min4_HFT", "gRefmult for Vz>-6 && Vz<=-4", 101, -0.5, 100.5));
    mOutList->Add(new TH1F("h_gRefmult_Vz_min4_min2_HFT", "gRefmult for Vz>-4 && Vz<=-2", 101, -0.5, 100.5));
    mOutList->Add(new TH1F("h_gRefmult_Vz_min2_0_HFT", "gRefmult for Vz>-2 && Vz<=0", 101, -0.5, 100.5));
    mOutList->Add(new TH1F("h_gRefmult_Vz_0_2_HFT", "gRefmult for Vz>0 && Vz<=2", 101, -0.5, 100.5));
    mOutList->Add(new TH1F("h_gRefmult_Vz_2_4_HFT", "gRefmult for Vz>2 && Vz<=4", 101, -0.5, 100.5));
    mOutList->Add(new TH1F("h_gRefmult_Vz_4_6_HFT", "gRefmult for Vz>4 && Vz<=6", 101, -0.5, 100.5));

    mOutList->Add(new TH2F("h_gRefmult_HFT_hybridTOF", "gRefmult;CountsvsRunIndex", 700, 0, 700, RunNumberVector.size() + 1, -1, RunNumberVector.size()));
    mOutList->Add(new TH2F("h_gRefmult_vs_ZDCx_HFT_hybridTOF", "gRefmult;ZDCx", 250, 0, 250, 101, -0.5, 100.5));
    mOutList->Add(new TH2F("h_gRefmult_vs_BBCx_HFT_hybridTOF", "gRefmult;BBCx", 1350, 0, 1350, 101, -0.5, 100.5));

    mOutList->Add(new TH1F("h_gRefmult_Vz_min6_min4_HFT_hybridTOF", "gRefmult for Vz>-6 && Vz<=-4", 101, -0.5, 100.5));
    mOutList->Add(new TH1F("h_gRefmult_Vz_min4_min2_HFT_hybridTOF", "gRefmult for Vz>-4 && Vz<=-2", 101, -0.5, 100.5));
    mOutList->Add(new TH1F("h_gRefmult_Vz_min2_0_HFT_hybridTOF", "gRefmult for Vz>-2 && Vz<=0", 101, -0.5, 100.5));
    mOutList->Add(new TH1F("h_gRefmult_Vz_0_2_HFT_hybridTOF", "gRefmult for Vz>0 && Vz<=2", 101, -0.5, 100.5));
    mOutList->Add(new TH1F("h_gRefmult_Vz_2_4_HFT_hybridTOF", "gRefmult for Vz>2 && Vz<=4", 101, -0.5, 100.5));
    mOutList->Add(new TH1F("h_gRefmult_Vz_4_6_HFT_hybridTOF", "gRefmult for Vz>4 && Vz<=6", 101, -0.5, 100.5));

    mOutList->Add(new TH2F("h_mh1gRefmultCorWg", "gRefmultCorWg;gRefmultCorWg;CountsvsRunIndex", 700, 0, 700, RunNumberVector.size() + 1, -1, RunNumberVector.size()));
    mOutList->Add(new TH3F("h_mh2CentVz", "CentralityVsVz;cent;VzvsRunIndex", 10, -1.5, 8.5, 200, -10, 10, RunNumberVector.size() + 1, -1, RunNumberVector.size()));
    mOutList->Add(new TH3F("h_mh2CentVzWg", "CentralityVsVzWg;cent;VzvsRunIndex", 10, -1.5, 8.5, 200, -10, 10, RunNumberVector.size() + 1, -1, RunNumberVector.size()));

    mOutList->Add(new TH2D("h_QA_Vz", "Vz_vs_RunIndex", 200, -10, 10, RunNumberVector.size() + 1, -1, RunNumberVector.size()));
    mOutList->Add(new TH2D("h_QA_VzmVzVPD", "Vz-VzVPD_vs_RunIndex", 100, -10, 10, RunNumberVector.size() + 1, -1, RunNumberVector.size()));

    mOutList->Add(new TH2F("h_QA_Vz_position", "VzDef_vs_VzKF", 1200, -6, 6, 1200, -6, 6));
    mOutList->Add(new TH2F("h_QA_Vx_position", "VxDef_vs_VxKF", 1200, -6, 6, 1200, -6, 6));
    mOutList->Add(new TH2F("h_QA_Vy_position", "VyDef_vs_VyKF", 1200, -6, 6, 1200, -6, 6));

    mOutList->Add(new TH2D("h_QA_ZDC_rate", "ZDC_rateVsRunIndex", 250, 0, 250, RunNumberVector.size() + 1, -1, RunNumberVector.size())); //check binning and range
    mOutList->Add(new TH2D("h_QA_BBC_rate", "BBC_rateVsRunIndex", 2000, 0, 2000, RunNumberVector.size() + 1, -1, RunNumberVector.size())); //check binning and range
    mOutList->Add(new TH2D("h_QA_ZDC_over_BBC", "ZDC/BBC_rateVsRunIndex", 5000, 0, 5, RunNumberVector.size() + 1, -1, RunNumberVector.size())); //check binning and range

    mOutList->Add(new TH2D("h_QA_ZDC_rate_pileUp", "ZDC_rateVsRunIndex", 250, 0, 250, RunNumberVector.size() + 1, -1, RunNumberVector.size())); //check binning and range
    mOutList->Add(new TH2D("h_QA_BBC_rate_pileUp", "BBC_rateVsRunIndex", 2000, 0, 2000, RunNumberVector.size() + 1, -1, RunNumberVector.size())); //check binning and range

    mOutList->Add(new TH2D("h_QA_ZDC_rate_pileUp_TOF", "ZDC_rate_TOFVsRunIndex", 250, 0, 250, RunNumberVector.size() + 1, -1, RunNumberVector.size())); //check binning and range
    mOutList->Add(new TH2D("h_QA_BBC_rate_pileUp_TOF", "BBC_rate_TOFVsRunIndex", 2000, 0, 2000, RunNumberVector.size() + 1, -1, RunNumberVector.size())); //check binning and range

    mOutList->Add(new TH2D("h_QA_ZDC_rate_TOF_0406", "ZDC_rate_TOF_0406VsRunIndex", 250, 0, 250, RunNumberVector.size() + 1, -1, RunNumberVector.size())); //check binning and range
    mOutList->Add(new TH2D("h_QA_BBC_rate_TOF_0406", "BBC_rate_TOF_0406VsRunIndex", 2000, 0, 2000, RunNumberVector.size() + 1, -1, RunNumberVector.size())); //check binning and range

    mOutList->Add(new TH2D("h_QA_ZDC_rate_TOF_1012", "ZDC_rate_TOF_1012VsRunIndex", 250, 0, 250, RunNumberVector.size() + 1, -1, RunNumberVector.size())); //check binning and range
    mOutList->Add(new TH2D("h_QA_BBC_rate_TOF_1012", "BBC_rate_TOF_1012VsRunIndex", 2000, 0, 2000, RunNumberVector.size() + 1, -1, RunNumberVector.size())); //check binning and range

    mOutList->Add(new TH2D("h_QA_ZDC_rate_TOF_3040", "ZDC_rate_TOF_3040VsRunIndex", 250, 0, 250, RunNumberVector.size() + 1, -1, RunNumberVector.size())); //check binning and range
    mOutList->Add(new TH2D("h_QA_BBC_rate_TOF_3040", "BBC_rate_TOF_3040VsRunIndex", 2000, 0, 2000, RunNumberVector.size() + 1, -1, RunNumberVector.size())); //check binning and range

    mOutList->Add(new TH2D("h_QA_ZDC_rate_TOF_HFT", "ZDC_rate_TOF_HFTVsRunIndex", 250, 0, 250, RunNumberVector.size() + 1, -1, RunNumberVector.size())); //check binning and range
    mOutList->Add(new TH2D("h_QA_BBC_rate_TOF_HFT", "BBC_rate_TOF_HFTVsRunIndex", 2000, 0, 2000, RunNumberVector.size() + 1, -1, RunNumberVector.size())); //check binning and range

    mOutList->Add(new TH2D("h_QA_ZDC_rate_TOF_HFT_1012", "ZDC_rate_TOF_HFT_1012VsRunIndex", 250, 0, 250, RunNumberVector.size() + 1, -1, RunNumberVector.size())); //check binning and range
    mOutList->Add(new TH2D("h_QA_BBC_rate_TOF_HFT_1012", "BBC_rate_TOF_HFT_1012VsRunIndex", 2000, 0, 2000, RunNumberVector.size() + 1, -1, RunNumberVector.size())); //check binning and range

    mOutList->Add(new TH2D("h_QA_ZDC_rate_TOF_HFT_3040", "ZDC_rate_TOF_HFT_3040VsRunIndex", 250, 0, 250, RunNumberVector.size() + 1, -1, RunNumberVector.size())); //check binning and range
    mOutList->Add(new TH2D("h_QA_BBC_rate_TOF_HFT_3040", "BBC_rate_TOF_HFT_3040VsRunIndex", 2000, 0, 2000, RunNumberVector.size() + 1, -1, RunNumberVector.size())); //check binning and range


    mOutList->Add(new TH2D("h_QA_BBC_rate_kaons_matching", "BBC_rate_kaons_matching_Vs_Nkaons", 400, 0, 2000, 600, 0, 600)); //check binning and range
    mOutList->Add(new TH2D("h_QA_BBC_rate_pions_matching", "BBC_rate_pions_matching_Vs_Npions", 400, 0, 2000, 600, 0, 600)); //check binning and range

    mOutList->Add(new TH2D("h_QA_BBC_rate_kaons", "BBC_rate_kaons_Vs_Nkaons", 400, 0, 2000, 600, 0, 600)); //check binning and range
    mOutList->Add(new TH2D("h_QA_BBC_rate_pions", "BBC_rate_pions_Vs_Npions", 400, 0, 2000, 600, 0, 600)); //check binning and range

    mOutList->Add(new TH2D("h_QA_BBC_rate_kaons_HFT", "BBC_rate_kaons_Vs_NHFTkaons_", 400, 0, 2000, 600, 0, 600)); //check binning and range
    mOutList->Add(new TH2D("h_QA_BBC_rate_pions_HFT", "BBC_rate_pions_Vs_NHFTpions", 400, 0, 2000, 600, 0, 600)); //check binning and range

    mOutList->Add(new TH2D("h_QA_BBC_rate_kaons_HFT_TOF", "BBC_rate_kaons_Vs_N_HFT_TOF_kaons_", 400, 0, 2000, 600, 0, 600)); //check binning and range
    mOutList->Add(new TH2D("h_QA_BBC_rate_pions_HFT_TOF", "BBC_rate_pions_Vs_N_HFT_TOF_pions", 400, 0, 2000, 600, 0, 600)); //check binning and range

    mOutList->Add(new TH2D("h_QA_BBC_rate_kaons_HFT_hybridTOF", "BBC_rate_kaons_Vs_N_HFT_hybridTOF_kaons", 400, 0, 2000, 600, 0, 600)); //check binning and range
    mOutList->Add(new TH2D("h_QA_BBC_rate_pions_HFT_hybridTOF", "BBC_rate_pions_Vs_N_HFT_hybridTOF_pions", 400, 0, 2000, 600, 0, 600)); //check binning and range

//    mOutList->Add(new TH2D("h_HFT_ratio_ZDC", "h_HFT_ratio_ZDC", 250, 0, 250, 2, 0, 2));
//    mOutList->Add(new TH2D("h_HFT_ratio_BBC", "h_HFT_ratio_BBC", 2000, 0, 2000, 2, 0, 2));
//
//    mOutList->Add(new TH2D("h_HFT_ratio_ZDC_0406", "h_HFT_ratio_ZDC_0406", 250, 0, 250, 2, 0, 2));
//    mOutList->Add(new TH2D("h_HFT_ratio_BBC_0406", "h_HFT_ratio_BBC_0406", 2000, 0, 2000, 2, 0, 2));
//
//    mOutList->Add(new TH2D("h_HFT_ratio_ZDC_1012", "h_HFT_ratio_ZDC_1012", 250, 0, 250, 2, 0, 2));
//    mOutList->Add(new TH2D("h_HFT_ratio_BBC_1012", "h_HFT_ratio_BBC_1012", 2000, 0, 2000, 2, 0, 2));
//
//    mOutList->Add(new TH2D("h_HFT_ratio_ZDC_3040", "h_HFT_ratio_ZDC_3040", 250, 0, 250, 2, 0, 2));
//    mOutList->Add(new TH2D("h_HFT_ratio_BBC_3040", "h_HFT_ratio_BBC_3040", 2000, 0, 2000, 2, 0, 2));

    mOutList->Add(new TH2I("h_QA_reweight_isNaN", "reweight_isNaNvsRunIndex", 2, 0, 2, RunNumberVector.size() + 1, -1, RunNumberVector.size())); //check for bad weights for refmutCorr
    mOutList->Add(new TH1I("h_QA_nEvents", "Number_of_eventsvsRunIndex", RunNumberVector.size() + 1, -1, RunNumberVector.size())); //number of events in run

    mOutList->Add(new TH2D("h_QA_nTracks", "Number_of_tracksvsRunIndex", 200, 0, 200, RunNumberVector.size() + 1, -1, RunNumberVector.size()));    //total nuber of tracks in event (no cuts)

    mOutList->Add(new TH2D("h_QA_nTracks_HFT_TOF_0406", "Number_of_tracksvsRunIndex", 200, 0, 200, RunNumberVector.size() + 1, -1, RunNumberVector.size()));    //total nuber of tracks in event (no cuts)
    mOutList->Add(new TH2D("h_QA_nTracks_HFT_TOF_1012", "Number_of_tracksvsRunIndex", 200, 0, 200, RunNumberVector.size() + 1, -1, RunNumberVector.size()));    //total nuber of tracks in event (no cuts)
    mOutList->Add(new TH2D("h_QA_nTracks_HFT_TOF_3040", "Number_of_tracksvsRunIndex", 200, 0, 200, RunNumberVector.size() + 1, -1, RunNumberVector.size()));    //total nuber of tracks in event (no cuts)

    mOutList->Add(new TH2D("h_QA_nTracks_HFT_0406", "Number_of_tracksvsRunIndex", 200, 0, 200, RunNumberVector.size() + 1, -1, RunNumberVector.size()));    //total nuber of tracks in event (no cuts)
    mOutList->Add(new TH2D("h_QA_nTracks_HFT_1012", "Number_of_tracksvsRunIndex", 200, 0, 200, RunNumberVector.size() + 1, -1, RunNumberVector.size()));    //total nuber of tracks in event (no cuts)
    mOutList->Add(new TH2D("h_QA_nTracks_HFT_3040", "Number_of_tracksvsRunIndex", 200, 0, 200, RunNumberVector.size() + 1, -1, RunNumberVector.size()));    //total nuber of tracks in event (no cuts)

    mOutList->Add(new TH2D("h_QA_nTracks_TOF_0406", "Number_of_tracksvsRunIndex", 200, 0, 200, RunNumberVector.size() + 1, -1, RunNumberVector.size()));    //total nuber of tracks in event (no cuts)
    mOutList->Add(new TH2D("h_QA_nTracks_TOF_1012", "Number_of_tracksvsRunIndex", 200, 0, 200, RunNumberVector.size() + 1, -1, RunNumberVector.size()));    //total nuber of tracks in event (no cuts)
    mOutList->Add(new TH2D("h_QA_nTracks_TOF_3040", "Number_of_tracksvsRunIndex", 200, 0, 200, RunNumberVector.size() + 1, -1, RunNumberVector.size()));    //total nuber of tracks in event (no cuts)

    mOutList->Add(new TH2D("h_QA_nTracks_0406", "Number_of_tracksvsRunIndex", 200, 0, 200, RunNumberVector.size() + 1, -1, RunNumberVector.size()));    //total nuber of tracks in event (no cuts)
    mOutList->Add(new TH2D("h_QA_nTracks_1012", "Number_of_tracksvsRunIndex", 200, 0, 200, RunNumberVector.size() + 1, -1, RunNumberVector.size()));    //total nuber of tracks in event (no cuts)
    mOutList->Add(new TH2D("h_QA_nTracks_3040", "Number_of_tracksvsRunIndex", 200, 0, 200, RunNumberVector.size() + 1, -1, RunNumberVector.size()));    //total nuber of tracks in event (no cuts)

    mOutList->Add(new TH2D("h_QA_nTracks_TPC", "Number_of_TPC_tracksvsRunIndex", 200, 0, 200, RunNumberVector.size() + 1, -1, RunNumberVector.size())); //number of TPC tracks in run (heve to pass TPC cuts)
    mOutList->Add(new TH2D("h_QA_nTracks_HFT", "Number_of_HFT_tracksvsRunIndex", 200, 0, 200, RunNumberVector.size() + 1, -1, RunNumberVector.size())); //number of HFT tracks (have to pass track->isHFTTrack())
    mOutList->Add(new TH2D("h_QA_nTracks_HFT_TOF", "Number_of_HFT_TOF_tracksvsRunIndex", 200, 0, 200, RunNumberVector.size() + 1, -1, RunNumberVector.size())); //number of HFT tracks (have to pass track->isHFTTrack())
    mOutList->Add(new TH2D("h_QA_nTracks_HFT_PXL1", "Number_of_HFT_PXL1_tracksvsRunIndex", 200, 0, 20, RunNumberVector.size() + 1, -1, RunNumberVector.size())); //number of PXL1 tracks (have to pass track->hasPxl1Hit())
    mOutList->Add(new TH2D("h_QA_nTracks_HFT_PXL2", "Number_of_HFT_PXL2_tracksvsRunIndex", 200, 0, 20, RunNumberVector.size() + 1, -1, RunNumberVector.size())); //number of PXL2 tracks (have to pass track->hasPxl2Hit())
    mOutList->Add(new TH2D("h_QA_nTracks_HFT_IST", "Number_of_HFT_IST_tracksvsRunIndex", 200, 0, 20, RunNumberVector.size() + 1, -1, RunNumberVector.size()));    //number of IST tracks (have to pass track->hasIstHit())
    mOutList->Add(new TH2D("h_QA_nTracks_HFT_SSD", "Number_of_HFT_SSD_tracksvsRunIndex", 200, 0, 20, RunNumberVector.size() + 1, -1, RunNumberVector.size()));    //number of SSD tracks (have to pass track->hasSstHit())
    mOutList->Add(new TH2D("h_QA_nTracks_HFT_IST_or_SSD", "Number_of_HFT_IST_or_SSD_tracksvsRunIndex", 200, 0, 20, RunNumberVector.size() + 1, -1, RunNumberVector.size())); //number of IST or SST tracks (have to pass track->hasIstHit() || track->hasSstHit())
    mOutList->Add(new TH2D("h_QA_nTracks_TOF", "Number_of_TOF_tracksvsRunIndex", 200, 0, 200, RunNumberVector.size() + 1, -1, RunNumberVector.size())); //nuber of TOF tracks (have to have TOF info)
    mOutList->Add(new TH2D("h_QA_nTracks_BEMC", "Number_of_BEMC_tracks_vs_RunIndex", 200, 0, 20, RunNumberVector.size() + 1, -1, RunNumberVector.size())); //number of BEMC tracks

    mOutList->Add(new TH2D("h_QA_nHitsFit", "h_QA_nHitsFit", 60, -0.5, 59.5, RunNumberVector.size() + 1, -1, RunNumberVector.size()));
    mOutList->Add(new TH2D("h_QA_nHFTHits", "h_QA_nHFTHits", 6, -0.5, 5.5, RunNumberVector.size() + 1, -1, RunNumberVector.size()));
    mOutList->Add(new TH2D("h_QA_nHFTHitsTOF", "h_QA_nHFTHitsTOF", 6, -0.5, 5.5, RunNumberVector.size() + 1, -1, RunNumberVector.size()));
    mOutList->Add(new TH2D("h_QA_nHitsFitTOF", "h_QA_nHitsFitTOF", 60, -0.5, 59.5, RunNumberVector.size() + 1, -1, RunNumberVector.size()));
    mOutList->Add(new TH2D("h_QA_nHitsDedx", "h_QA_nHitsDedx",  60, -0.5, 59.5, RunNumberVector.size() + 1, -1, RunNumberVector.size()));
    mOutList->Add(new TH2D("h_QA_nSigmaKaon", "h_QA_nSigmaKaon", 2000, -10, 10, RunNumberVector.size() + 1, -1, RunNumberVector.size()));
    mOutList->Add(new TH2D("h_QA_nSigmaPion", "h_QA_nSigmaPion", 2000, -10, 10, RunNumberVector.size() + 1, -1, RunNumberVector.size()));

    mOutList->Add(new TH2D("h_QA_Beta", "h_QA_Beta", 3000, -1, 2, RunNumberVector.size() + 1, -1, RunNumberVector.size()));
    mOutList->Add(new TH2D("h_QA_OneOverBetaDiffPion", "h_QA_OneOverBetaDiffPion", 1000, -0.5, 0.5, RunNumberVector.size() + 1, -1, RunNumberVector.size()));
    mOutList->Add(new TH2D("h_QA_OneOverBetaDiffKaon", "h_QA_OneOverBetaDiffKaon", 1000, -0.5, 0.5, RunNumberVector.size() + 1, -1, RunNumberVector.size()));

    //_____________________HFT, eta cut_____________________________________________
    mOutList->Add(new TH2D("h_QA_nTracks_TPC_etaCut", "Number_of_TPC_etaCut_tracksvsRunIndex", 200, 0, 20, RunNumberVector.size() + 1, -1, RunNumberVector.size()));
    mOutList->Add(new TH2D("h_QA_nTracks_HFT_etaCut", "Number_of_HFT_etaCut_tracksvsRunIndex", 200, 0, 20, RunNumberVector.size() + 1, -1, RunNumberVector.size())); //number of HFT tracks (have to pass track->isHFTTrack())
    mOutList->Add(new TH2D("h_QA_nTracks_HFT_etaCut_PXL1", "Number_of_HFT_etaCut_PXL1_tracksvsRunIndex", 200, 0, 20, RunNumberVector.size() + 1, -1, RunNumberVector.size())); //number of PXL1 tracks (have to pass track->hasPxl1Hit())
    mOutList->Add(new TH2D("h_QA_nTracks_HFT_etaCut_PXL2", "Number_of_HFT_etaCut_PXL2_tracksvsRunIndex", 200, 0, 20, RunNumberVector.size() + 1, -1, RunNumberVector.size())); //number of PXL2 tracks (have to pass track->hasPxl2Hit())
    mOutList->Add(new TH2D("h_QA_nTracks_HFT_etaCut_IST", "Number_of_HFT_etaCut_IST_tracksvsRunIndex", 200, 0, 20, RunNumberVector.size() + 1, -1, RunNumberVector.size()));    //number of IST tracks (have to pass track->hasIstHit())
    mOutList->Add(new TH2D("h_QA_nTracks_HFT_etaCut_SSD", "Number_of_HFT_etaCut_SSD_tracksvsRunIndex", 200, 0, 20, RunNumberVector.size() + 1, -1, RunNumberVector.size()));    //number of SSD tracks (have to pass track->hasSstHit())
    mOutList->Add(new TH2D("h_QA_nTracks_HFT_etaCut_IST_or_SSD", "Number_of_HFT_etaCut_IST_or_SSD_tracksvsRunIndex", 200, 0, 20, RunNumberVector.size() + 1, -1, RunNumberVector.size())); //number of IST or SST tracks (have to pass track->hasIstHit() || track->hasSstHit())

    //Number of tracks vs. pT - get from integral of pT spectrum
//____general QA histograms (all tracks within (TPC, track quality) cuts, NO TOF and HFT)__________________________________
    mOutList->Add(new TH2F("h_QA_pT", "Transverse_momentum_TPCvsRunIndex", 200, 0, 20, RunNumberVector.size() + 1, -1, RunNumberVector.size())); //check binning and range
    mOutList->Add(new TH2F("h_QA_eta", "Pseudorapidity_TPCvsRunIndex", 400, -2., 2., RunNumberVector.size() + 1, -1, RunNumberVector.size()));
    mOutList->Add(new TH2F("h_QA_phi", "Azimuthal_angle_TPCvsRunIndex", 200, -TMath::Pi(), TMath::Pi(), RunNumberVector.size() + 1, -1, RunNumberVector.size()));

    mOutList->Add(new TH2F("h_QA_vertex_x", "Vertex_x_positionvsRunIndex", 2000, -10, 10, RunNumberVector.size() + 1, -1, RunNumberVector.size()));
    mOutList->Add(new TH2F("h_QA_vertex_z", "Vertex_z_positionvsRunIndex", 2000, -10, 10, RunNumberVector.size() + 1, -1, RunNumberVector.size()));
    mOutList->Add(new TH2F("h_QA_vertex_y", "Vertex_y_positionvsRunIndex", 2000, -10, 10, RunNumberVector.size() + 1, -1, RunNumberVector.size()));

    mOutList->Add(new TH2F("h_QA_DCA_xy_TPC", "DCA_xy_TPCvsRunIndex", 600, -3, 3, RunNumberVector.size() + 1, -1, RunNumberVector.size()));
    mOutList->Add(new TH2F("h_QA_DCA_xy_zoom_TPC", "DCA_xy_zoom_TPCvsRunIndex", 100, -0.1, 0.1, RunNumberVector.size() + 1, -1, RunNumberVector.size()));
    mOutList->Add(new TH2F("h_QA_DCA_z_TPC", "DCA_z_TPCvsRunIndex", 600, -3, 3, RunNumberVector.size() + 1, -1, RunNumberVector.size()));
    mOutList->Add(new TH2F("h_QA_DCA_TPC", "h_QA_DCA_TPC", 600, -3, 3, RunNumberVector.size() + 1, -1, RunNumberVector.size()));
    mOutList->Add(new TH2F("h_QA_DCA_HFT_TOF", "h_QA_DCA_HFT_TOF", 600, -3, 3, RunNumberVector.size() + 1, -1, RunNumberVector.size()));
    mOutList->Add(new TH2F("h_QA_DCA_HFT", "h_QA_DCA_HFT", 600, -3, 3, RunNumberVector.size() + 1, -1, RunNumberVector.size()));
    mOutList->Add(new TH2F("h_QA_DCA_z_zoom_TPC", "DCA_z_zoom_TPCvsRunIndex", 100, -0.1, 0.1, RunNumberVector.size() + 1, -1, RunNumberVector.size()));

//___HFT QA hitograms______________________________________________________________________________________________________
    mOutList->Add(new TH2F("h_QA_pT_HFT", "Transverse_momentum_HFTvsRunIndex", 200, 0, 20, RunNumberVector.size() + 1, -1, RunNumberVector.size())); //check binning and range
    mOutList->Add(new TH2F("h_QA_eta_HFT", "Pseudorapidity_HFTvsRunIndex", 400, -2., 2., RunNumberVector.size() + 1, -1, RunNumberVector.size()));
    mOutList->Add(new TH2F("h_QA_phi_HFT", "Azimuthal_angle_HFTvsRunIndex", 200, -TMath::Pi(), TMath::Pi(), RunNumberVector.size() + 1, -1, RunNumberVector.size()));

    mOutList->Add(new TH2F("h_QA_pT_HFT_TOF", "Transverse_momentum_HFT_TOFvsRunIndex", 200, 0, 20, RunNumberVector.size() + 1, -1, RunNumberVector.size())); //check binning and range
    mOutList->Add(new TH2F("h_QA_eta_HFT_TOF", "Pseudorapidity_HFT_TOFvsRunIndex", 200, -1., 1., RunNumberVector.size() + 1, -1, RunNumberVector.size()));
    mOutList->Add(new TH2F("h_QA_phi_HFT_TOF", "Azimuthal_angle_HFT_TOFvsRunIndex", 200, -TMath::Pi(), TMath::Pi(), RunNumberVector.size() + 1, -1, RunNumberVector.size()));

    mOutList->Add(new TH2F("h_QA_vertex_x_HFT", "Vertex_x_position_HFTvsRunIndex", 100, -0.5, 0.5, RunNumberVector.size() + 1, -1, RunNumberVector.size()));
    mOutList->Add(new TH2F("h_QA_vertex_y_HFT", "Vertex_y_position_HFTvsRunIndex", 100, -0.5, 0.5, RunNumberVector.size() + 1, -1, RunNumberVector.size()));
    mOutList->Add(new TH2F("h_QA_vertex_z_HFT", "Vertex_z_position_HFTvsRunIndex", 100, -0.5, 0.5, RunNumberVector.size() + 1, -1, RunNumberVector.size()));

    mOutList->Add(new TH2F("h_QA_DCA_xy_HFT", "DCA_xy_HFTvsRunIndex", 600, -3, 3, RunNumberVector.size() + 1, -1, RunNumberVector.size()));
    mOutList->Add(new TH2F("h_QA_DCA_xy_zoom_HFT", "DCA_xy_zoom_HFTvsRunIndex", 100, -0.1, 0.1, RunNumberVector.size() + 1, -1, RunNumberVector.size()));
    mOutList->Add(new TH2F("h_QA_DCA_z_HFT", "DCA_z_HFTvsRunIndex", 600, -3, 3, RunNumberVector.size() + 1, -1, RunNumberVector.size()));
    mOutList->Add(new TH2F("h_QA_DCA_z_zoom_HFT", "DCA_z_zoom_HFTvsRunIndex", 100, -0.1, 0.1, RunNumberVector.size() + 1, -1, RunNumberVector.size()));

    mOutList->Add(new TH2F("h_QA_DCA_xy_HFT_TOF", "DCA_xy_HFT_TOFvsRunIndex", 600, -3, 3, RunNumberVector.size() + 1, -1, RunNumberVector.size()));
    mOutList->Add(new TH2F("h_QA_DCA_xy_zoom_HFT_TOF", "DCA_xy_zoom_HFT_TOFvsRunIndex", 100, -0.1, 0.1, RunNumberVector.size() + 1, -1, RunNumberVector.size()));
    mOutList->Add(new TH2F("h_QA_DCA_z_HFT_TOF", "DCA_z_HFT_TOFvsRunIndex", 600, -3, 3, RunNumberVector.size() + 1, -1, RunNumberVector.size()));
    mOutList->Add(new TH2F("h_QA_DCA_z_zoom_HFT_TOF", "DCA_z_zoom_HFT_TOFvsRunIndex", 100, -0.1, 0.1, RunNumberVector.size() + 1, -1, RunNumberVector.size()));

//___HFT QA hitograms, eta cut______________________________________________________________________________________________________
    mOutList->Add(new TH2F("h_QA_pT_HFT_etaCut", "Transverse_momentum_HFT_etaCutvsRunIndex", 200, 0, 20, RunNumberVector.size() + 1, -1, RunNumberVector.size())); //check binning and range
    mOutList->Add(new TH2F("h_QA_eta_HFT_etaCut", "Pseudorapidity_HFT_etaCutvsRunIndex", 200, -1., 1., RunNumberVector.size() + 1, -1, RunNumberVector.size()));
    mOutList->Add(new TH2F("h_QA_phi_HFT_etaCut", "Azimuthal_angle_HFT_etaCutvsRunIndex", 200, -TMath::Pi(), TMath::Pi(), RunNumberVector.size() + 1, -1, RunNumberVector.size()));

    mOutList->Add(new TH2F("h_QA_vertex_x_HFT_etaCut", "Vertex_x_position_HFT_etaCutvsRunIndex", 100, -0.5, 0.5, RunNumberVector.size() + 1, -1, RunNumberVector.size()));
    mOutList->Add(new TH2F("h_QA_vertex_y_HFT_etaCut", "Vertex_y_position_HFT_etaCutvsRunIndex", 100, -0.5, 0.5, RunNumberVector.size() + 1, -1, RunNumberVector.size()));
    mOutList->Add(new TH2F("h_QA_vertex_z_HFT_etaCut", "Vertex_z_position_HFT_etaCutvsRunIndex", 100, -0.5, 0.5, RunNumberVector.size() + 1, -1, RunNumberVector.size()));

    mOutList->Add(new TH2F("h_QA_DCA_xy_HFT_etaCut", "DCA_xy_HFT_etaCutvsRunIndex", 600, -3, 3, RunNumberVector.size() + 1, -1, RunNumberVector.size()));
    mOutList->Add(new TH2F("h_QA_DCA_xy_zoom_HFT_etaCut", "DCA_xy_zoom_HFT_etaCutvsRunIndex", 100, -0.1, 0.1, RunNumberVector.size() + 1, -1, RunNumberVector.size()));
    mOutList->Add(new TH2F("h_QA_DCA_z_HFT_etaCut", "DCA_z_HFT_etaCutvsRunIndex", 600, -3, 3, RunNumberVector.size() + 1, -1, RunNumberVector.size()));
    mOutList->Add(new TH2F("h_QA_DCA_z_zoom_HFT_etaCut", "DCA_z_zoom_HFT_etaCutvsRunIndex", 100, -0.1, 0.1, RunNumberVector.size() + 1, -1, RunNumberVector.size()));

//____TOF QA histograms______________________________________________________________________________________________________
    mOutList->Add(new TH2F("h_QA_pT_TOF", "Transverse_momentum_TOFvsRunIndex", 200, 0, 20, RunNumberVector.size() + 1, -1, RunNumberVector.size())); //check binning and range
    mOutList->Add(new TH2F("h_QA_eta_TOF", "Pseudorapidity_TOFvsRunIndex", 200, -1., 1., RunNumberVector.size() + 1, -1, RunNumberVector.size()));
    mOutList->Add(new TH2F("h_QA_phi_TOF", "Azimuthal_angle_TOFvsRunIndex", 200, -TMath::Pi(), TMath::Pi(), RunNumberVector.size() + 1, -1, RunNumberVector.size()));

    mOutList->Add(new TH2F("h_QA_vertex_x_TOF", "Vertex_x_position_TOFvsRunIndex", 100, -0.5, 0.5, RunNumberVector.size() + 1, -1, RunNumberVector.size()));
    mOutList->Add(new TH2F("h_QA_vertex_y_TOF", "Vertex_y_position_TOFvsRunIndex", 100, -0.5, 0.5, RunNumberVector.size() + 1, -1, RunNumberVector.size()));
    mOutList->Add(new TH2F("h_QA_vertex_z_TOF", "Vertex_z_position_TOFvsRunIndex", 100, -0.5, 0.5, RunNumberVector.size() + 1, -1, RunNumberVector.size()));
    mOutList->Add(new TH2F("h_QA_DCA_TOF", "h_QA_DCA_TOF", 600, -3, 3, RunNumberVector.size() + 1, -1, RunNumberVector.size()));

//___BEMC QA histograms_______________________________________________________________________________________________________
    mOutList->Add(new TH2F("h_QA_BEMC_TOWId", "BEMC_Mached_tower_id_vs_RunIndex", 4801, 0, 4800, RunNumberVector.size() + 1, -1, RunNumberVector.size()));
    mOutList->Add(new TH2F("h_QA_BSMD_nEta", "BSMD_Eta_wires_vs_RunIndex", 11, 0, 10, RunNumberVector.size() + 1, -1, RunNumberVector.size()));
    mOutList->Add(new TH2F("h_QA_BSMD_nPhi", "BSMD_Phi_wires_vs_RunIndex", 11, 0, 10, RunNumberVector.size() + 1, -1, RunNumberVector.size()));
    mOutList->Add(new TH2F("h_pOverE", "E/p_vs_RunIndex", 200, 0, 10, RunNumberVector.size() + 1, -1, RunNumberVector.size()));

    mRunNumber = 0;

    mOutFileBaseName = mOutFileBaseName.ReplaceAll(".root", "");
//    ntp_event = new TNtuple("ntp_event","event tuple","nTrk:nTrk0406:nTrk1012:nTrk3040:nTrkHft:nTrkHft0406:nTrkHft1012:nTrkHft3040:ZDC:BBC:runId:RunIndex");
//    ntp_hft_track = new TNtuple("ntp_hft_track", "hft track tree","pt:zdc:bbc:runId:dca:dca_xy:dca_z");
    return kStOK;
}

// _________________________________________________________
void StPicoQAMaker::ClearHF(Option_t *opt="") {
    return;
}

// _________________________________________________________
int StPicoQAMaker::FinishHF() {
//    ntp_event -> Write(ntp_event->GetName(), TObject::kOverwrite);
    return kStOK;
}

// _________________________________________________________
int StPicoQAMaker::MakeHF() {
//    float dcaxy, dcaz, dca;
//    float isHft;
//    TH1D *phi_hft = static_cast<TH1D*>(mOutList->FindObject("h_phi_hft_track"));
//    TH1D *eta_hft = static_cast<TH1D*>(mOutList->FindObject("h_eta_hft_track"));
//    TH1D *phi = static_cast<TH1D*>(mOutList->FindObject("h_phi_track"));
//    TH1D *eta = static_cast<TH1D*>(mOutList->FindObject("h_eta_track"));

//    float bbc_rate = event->BBCx()/1000.;
//    float zdc_rate = event->ZDCx()/1000.;
//    UInt_t nTracks = mPicoDst->numberOfTracks();
//    for (unsigned short iTrack = 0; iTrack < nTracks; ++iTrack) {
//        isHft = 0;
//        StPicoTrack const* trk = mPicoDst->track(iTrack);
//        if (!trk) continue;
//        if(!mHFCuts->isGoodTrack(trk)) continue;  //TOF match, HFT, NHitFit, dcaMax - hft needs to be set to fault to make all needed comp
//        float pt = trk->gPt();
//        if (pt<0.2) continue;
////        StPicoPhysicalHelix helix = trk->helix(b);
////        dcaxy =helix.geometricSignedDistance(vtx.x(),vtx.y());
////        dcaz =(trk->origin().z() - vtx.z());
////        dca = (vtx-trk->origin()).Mag();
//        phi->Fill(trk->gMom(vtx,b).Phi());
//        eta->Fill(trk->gMom(vtx,b).PseudoRapidity());
//        if(trk->isHFTTrack()) {
//            isHft = 1;
//            phi_hft->Fill(trk->gMom(vtx,b).Phi());
//            eta_hft->Fill(trk->gMom(vtx,b).PseudoRapidity());
////            ntp_hft_track->Fill(pt,zdc_rate,bbc_rate,event->runId(),dca,dcaxy,dcaz);
//        }
////        ntp_track->Fill(isHft,pt,dca,zdc_rate,bbc_rate,event->runId());
//    }

    RunId = mPicoDst->event()->runId();
    int RunIndex = -1; //default value for RunIndex (does not correspond to any RunId)

    for (unsigned int i = 0; i < RunNumberVector.size(); i++) //find corresponding RunIndex to a given RunId
    {
        if (RunNumberVector.at(i) == RunId) {
            RunIndex = i;
            break; //do not continue if coresponding RunId is found
        }
    }

    TH2F *h_gRefmult = static_cast<TH2F *>(mOutList->FindObject("h_gRefmult"));
    TH2F *h_gRefmult_vs_ZDCx = static_cast<TH2F *>(mOutList->FindObject("h_gRefmult_vs_ZDCx"));
    TH2F *h_gRefmult_vs_BBCx = static_cast<TH2F *>(mOutList->FindObject("h_gRefmult_vs_BBCx"));

    TH1F *h_gRefmult_Vz_min6_min4 = static_cast<TH1F *>(mOutList->FindObject("h_gRefmult_Vz_min6_min4"));
    TH1F *h_gRefmult_Vz_min4_min2 = static_cast<TH1F *>(mOutList->FindObject("h_gRefmult_Vz_min4_min2"));
    TH1F *h_gRefmult_Vz_min2_0 = static_cast<TH1F *>(mOutList->FindObject("h_gRefmult_Vz_min2_0"));
    TH1F *h_gRefmult_Vz_0_2 = static_cast<TH1F *>(mOutList->FindObject("h_gRefmult_Vz_0_2"));
    TH1F *h_gRefmult_Vz_2_4 = static_cast<TH1F *>(mOutList->FindObject("h_gRefmult_Vz_2_4"));
    TH1F *h_gRefmult_Vz_4_6 = static_cast<TH1F *>(mOutList->FindObject("h_gRefmult_Vz_4_6"));

    TH2F *h_gRefmult_HFT = static_cast<TH2F *>(mOutList->FindObject("h_gRefmult_HFT"));
    TH2F *h_gRefmult_vs_ZDCx_HFT = static_cast<TH2F *>(mOutList->FindObject("h_gRefmult_vs_ZDCx_HFT"));
    TH2F *h_gRefmult_vs_BBCx_HFT = static_cast<TH2F *>(mOutList->FindObject("h_gRefmult_vs_BBCx_HFT"));

    TH1F *h_gRefmult_Vz_min6_min4_HFT = static_cast<TH1F *>(mOutList->FindObject("h_gRefmult_Vz_min6_min4_HFT"));
    TH1F *h_gRefmult_Vz_min4_min2_HFT = static_cast<TH1F *>(mOutList->FindObject("h_gRefmult_Vz_min4_min2_HFT"));
    TH1F *h_gRefmult_Vz_min2_0_HFT = static_cast<TH1F *>(mOutList->FindObject("h_gRefmult_Vz_min2_0_HFT"));
    TH1F *h_gRefmult_Vz_0_2_HFT = static_cast<TH1F *>(mOutList->FindObject("h_gRefmult_Vz_0_2_HFT"));
    TH1F *h_gRefmult_Vz_2_4_HFT = static_cast<TH1F *>(mOutList->FindObject("h_gRefmult_Vz_2_4_HFT"));
    TH1F *h_gRefmult_Vz_4_6_HFT = static_cast<TH1F *>(mOutList->FindObject("h_gRefmult_Vz_4_6_HFT"));

    TH2F *h_gRefmult_HFT_hybridTOF = static_cast<TH2F *>(mOutList->FindObject("h_gRefmult_HFT_hybridTOF"));
    TH2F *h_gRefmult_vs_ZDCx_HFT_hybridTOF = static_cast<TH2F *>(mOutList->FindObject("h_gRefmult_vs_ZDCx_HFT_hybridTOF"));
    TH2F *h_gRefmult_vs_BBCx_HFT_hybridTOF = static_cast<TH2F *>(mOutList->FindObject("h_gRefmult_vs_BBCx_HFT_hybridTOF"));

    TH1F *h_gRefmult_Vz_min6_min4_HFT_hybridTOF = static_cast<TH1F *>(mOutList->FindObject("h_gRefmult_Vz_min6_min4_HFT_hybridTOF"));
    TH1F *h_gRefmult_Vz_min4_min2_HFT_hybridTOF = static_cast<TH1F *>(mOutList->FindObject("h_gRefmult_Vz_min4_min2_HFT_hybridTOF"));
    TH1F *h_gRefmult_Vz_min2_0_HFT_hybridTOF = static_cast<TH1F *>(mOutList->FindObject("h_gRefmult_Vz_min2_0_HFT_hybridTOF"));
    TH1F *h_gRefmult_Vz_0_2_HFT_hybridTOF = static_cast<TH1F *>(mOutList->FindObject("h_gRefmult_Vz_0_2_HFT_hybridTOF"));
    TH1F *h_gRefmult_Vz_2_4_HFT_hybridTOF = static_cast<TH1F *>(mOutList->FindObject("h_gRefmult_Vz_2_4_HFT_hybridTOF"));
    TH1F *h_gRefmult_Vz_4_6_HFT_hybridTOF = static_cast<TH1F *>(mOutList->FindObject("h_gRefmult_Vz_4_6_HFT_hybridTOF"));

    TH2D *h_QA_Vz = static_cast<TH2D *>(mOutList->FindObject("h_QA_Vz"));
    TH2F *h_QA_Vz_position = static_cast<TH2F *>(mOutList->FindObject("h_QA_Vz_position"));
    TH2F *h_QA_Vx_position = static_cast<TH2F *>(mOutList->FindObject("h_QA_Vx_position"));
    TH2F *h_QA_Vy_position = static_cast<TH2F *>(mOutList->FindObject("h_QA_Vy_position"));

    TH2D *h_QA_VzmVzVPD = static_cast<TH2D *>(mOutList->FindObject("h_QA_VzmVzVPD"));

    TH2D *h_QA_ZDC_rate = static_cast<TH2D *>(mOutList->FindObject("h_QA_ZDC_rate"));
    TH2D *h_QA_BBC_rate = static_cast<TH2D *>(mOutList->FindObject("h_QA_BBC_rate"));
    TH2D *h_QA_ZDC_over_BBC = static_cast<TH2D *>(mOutList->FindObject("h_QA_ZDC_over_BBC"));

    TH2D *h_QA_ZDC_rate_pileUp = static_cast<TH2D *>(mOutList->FindObject("h_QA_ZDC_rate_pileUp"));
    TH2D *h_QA_BBC_rate_pileUp = static_cast<TH2D *>(mOutList->FindObject("h_QA_BBC_rate_pileUp"));

    TH2D *h_QA_ZDC_rate_pileUp_TOF = static_cast<TH2D *>(mOutList->FindObject("h_QA_ZDC_rate_pileUp_TOF"));
    TH2D *h_QA_BBC_rate_pileUp_TOF = static_cast<TH2D *>(mOutList->FindObject("h_QA_BBC_rate_pileUp_TOF"));

    TH2D *h_QA_ZDC_rate_TOF_0406 = static_cast<TH2D *>(mOutList->FindObject("h_QA_ZDC_rate_TOF_0406"));
    TH2D *h_QA_BBC_rate_TOF_0406 = static_cast<TH2D *>(mOutList->FindObject("h_QA_BBC_rate_TOF_0406"));

    TH2D *h_QA_ZDC_rate_TOF_1012 = static_cast<TH2D *>(mOutList->FindObject("h_QA_ZDC_rate_TOF_1012"));
    TH2D *h_QA_BBC_rate_TOF_1012 = static_cast<TH2D *>(mOutList->FindObject("h_QA_BBC_rate_TOF_1012"));

    TH2D *h_QA_ZDC_rate_TOF_3040 = static_cast<TH2D *>(mOutList->FindObject("h_QA_ZDC_rate_TOF_3040"));
    TH2D *h_QA_BBC_rate_TOF_3040 = static_cast<TH2D *>(mOutList->FindObject("h_QA_BBC_rate_TOF_3040"));

    TH2D *h_QA_ZDC_rate_TOF_HFT = static_cast<TH2D *>(mOutList->FindObject("h_QA_ZDC_rate_TOF_HFT"));
    TH2D *h_QA_BBC_rate_TOF_HFT = static_cast<TH2D *>(mOutList->FindObject("h_QA_BBC_rate_TOF_HFT"));

    TH2D *h_QA_ZDC_rate_TOF_HFT_1012 = static_cast<TH2D *>(mOutList->FindObject("h_QA_ZDC_rate_TOF_HFT_1012"));
    TH2D *h_QA_BBC_rate_TOF_HFT_1012 = static_cast<TH2D *>(mOutList->FindObject("h_QA_BBC_rate_TOF_HFT_1012"));

    TH2D *h_QA_ZDC_rate_TOF_HFT_3040 = static_cast<TH2D *>(mOutList->FindObject("h_QA_ZDC_rate_TOF_HFT_3040"));
    TH2D *h_QA_BBC_rate_TOF_HFT_3040 = static_cast<TH2D *>(mOutList->FindObject("h_QA_BBC_rate_TOF_HFT_3040"));

//    TH2D *h_HFT_ratio_ZDC = static_cast<TH2D *>(mOutList->FindObject("h_HFT_ratio_ZDC"));
//    TH2D *h_HFT_ratio_BBC = static_cast<TH2D *>(mOutList->FindObject("h_HFT_ratio_BBC"));
//    TH2D *h_HFT_ratio_ZDC_0406 = static_cast<TH2D *>(mOutList->FindObject("h_HFT_ratio_ZDC_0406"));
//    TH2D *h_HFT_ratio_BBC_0406 = static_cast<TH2D *>(mOutList->FindObject("h_HFT_ratio_BBC_0406"));
//    TH2D *h_HFT_ratio_ZDC_1012 = static_cast<TH2D *>(mOutList->FindObject("h_HFT_ratio_ZDC_1012"));
//    TH2D *h_HFT_ratio_BBC_1012 = static_cast<TH2D *>(mOutList->FindObject("h_HFT_ratio_BBC_1012"));
//    TH2D *h_HFT_ratio_ZDC_3040 = static_cast<TH2D *>(mOutList->FindObject("h_HFT_ratio_ZDC_3040"));
//    TH2D *h_HFT_ratio_BBC_3040 = static_cast<TH2D *>(mOutList->FindObject("h_HFT_ratio_BBC_3040"));

    TH2D *h_QA_BBC_rate_pions = static_cast<TH2D*>(mOutList->FindObject("h_QA_BBC_rate_pions"));
    TH2D *h_QA_BBC_rate_kaons = static_cast<TH2D*>(mOutList->FindObject("h_QA_BBC_rate_kaons"));

    TH2D *h_QA_BBC_rate_kaons_matching = static_cast<TH2D*>(mOutList->FindObject("h_QA_BBC_rate_kaons_matching"));
    TH2D *h_QA_BBC_rate_pions_matching = static_cast<TH2D*>(mOutList->FindObject("h_QA_BBC_rate_pions_matching"));

    TH2D *h_QA_BBC_rate_kaons_HFT = static_cast<TH2D*>(mOutList->FindObject("h_QA_BBC_rate_kaons_HFT"));
    TH2D *h_QA_BBC_rate_pions_HFT = static_cast<TH2D*>(mOutList->FindObject("h_QA_BBC_rate_pions_HFT"));

    TH2D *h_QA_BBC_rate_kaons_HFT_TOF = static_cast<TH2D*>(mOutList->FindObject("h_QA_BBC_rate_kaons_HFT_TOF"));
    TH2D *h_QA_BBC_rate_pions_HFT_TOF = static_cast<TH2D*>(mOutList->FindObject("h_QA_BBC_rate_pions_HFT_TOF"));

    TH2D *h_QA_BBC_rate_kaons_HFT_hybridTOF = static_cast<TH2D*>(mOutList->FindObject("h_QA_BBC_rate_kaons_HFT_hybridTOF"));
    TH2D *h_QA_BBC_rate_pions_HFT_hybridTOF = static_cast<TH2D*>(mOutList->FindObject("h_QA_BBC_rate_pions_HFT_hybridTOF"));

    TH2F *h_QA_reweight_isNaN = static_cast<TH2F *>(mOutList->FindObject("h_QA_reweight_isNaN"));
    TH1I *h_QA_nEvents = static_cast<TH1I *>(mOutList->FindObject("h_QA_nEvents"));

    TH2D *h_QA_nHitsFit = static_cast<TH2D *>(mOutList->FindObject("h_QA_nHitsFit"));
    TH2D *h_QA_nHitsFitTOF = static_cast<TH2D *>(mOutList->FindObject("h_QA_nHitsFitTOF"));
    TH2D *h_QA_nHFTHits = static_cast<TH2D *>(mOutList->FindObject("h_QA_nHFTHits"));
    TH2D *h_QA_nHFTHitsTOF = static_cast<TH2D *>(mOutList->FindObject("h_QA_nHFTHitsTOF"));
    TH2D *h_QA_nSigmaKaon = static_cast<TH2D *>(mOutList->FindObject("h_QA_nSigmaKaon"));
    TH2D *h_QA_nSigmaPion = static_cast<TH2D *>(mOutList->FindObject("h_QA_nSigmaPion"));
    TH2D *h_QA_nHitsDedx = static_cast<TH2D *>(mOutList->FindObject("h_QA_nHitsDedx"));
    TH2D *h_QA_OneOverBetaDiffKaon = static_cast<TH2D *>(mOutList->FindObject("h_QA_OneOverBetaDiffKaon"));
    TH2D *h_QA_OneOverBetaDiffPion = static_cast<TH2D *>(mOutList->FindObject("h_QA_OneOverBetaDiffPion"));
    TH2D *h_QA_Beta = static_cast<TH2D *>(mOutList->FindObject("h_QA_Beta"));

    TH2D *h_QA_nTracks = static_cast<TH2D *>(mOutList->FindObject("h_QA_nTracks"));
    TH2D *h_QA_nTracks_TPC = static_cast<TH2D *>(mOutList->FindObject("h_QA_nTracks_TPC"));
    TH2D *h_QA_nTracks_HFT = static_cast<TH2D *>(mOutList->FindObject("h_QA_nTracks_HFT"));
    TH2D *h_QA_nTracks_HFT_PXL1 = static_cast<TH2D *>(mOutList->FindObject("h_QA_nTracks_HFT_PXL1"));
    TH2D *h_QA_nTracks_HFT_PXL2 = static_cast<TH2D *>(mOutList->FindObject("h_QA_nTracks_HFT_PXL2"));
    TH2D *h_QA_nTracks_HFT_IST = static_cast<TH2D *>(mOutList->FindObject("h_QA_nTracks_HFT_IST"));
    TH2D *h_QA_nTracks_HFT_SSD = static_cast<TH2D *>(mOutList->FindObject("h_QA_nTracks_HFT_SSD"));
    TH2D *h_QA_nTracks_HFT_IST_or_SSD = static_cast<TH2D *>(mOutList->FindObject("h_QA_nTracks_HFT_IST_or_SSD"));
    TH2D *h_QA_nTracks_TOF = static_cast<TH2D *>(mOutList->FindObject("h_QA_nTracks_TOF"));
    TH2D *h_QA_nTracks_BEMC = static_cast<TH2D *>(mOutList->FindObject("h_QA_nTracks_BEMC"));

    TH2D *h_QA_nTracks_HFT_etaCut = static_cast<TH2D *>(mOutList->FindObject("h_QA_nTracks_HFT_etaCut"));
    TH2D *h_QA_nTracks_HFT_etaCut_PXL1 = static_cast<TH2D *>(mOutList->FindObject("h_QA_nTracks_HFT_etaCut_PXL1"));
    TH2D *h_QA_nTracks_HFT_etaCut_PXL2 = static_cast<TH2D *>(mOutList->FindObject("h_QA_nTracks_HFT_etaCut_PXL2"));
    TH2D *h_QA_nTracks_HFT_etaCut_IST = static_cast<TH2D *>(mOutList->FindObject("h_QA_nTracks_HFT_etaCut_IST"));
    TH2D *h_QA_nTracks_HFT_etaCut_SSD = static_cast<TH2D *>(mOutList->FindObject("h_QA_nTracks_HFT_etaCut_SSD"));
    TH2D *h_QA_nTracks_HFT_etaCut_IST_or_SSD = static_cast<TH2D *>(mOutList->FindObject("h_QA_nTracks_HFT_etaCut_IST_or_SSD"));

    TH2D *h_QA_nTracks_HFT_TOF_0406 = static_cast<TH2D *>(mOutList->FindObject("h_QA_nTracks_HFT_TOF_0406"));
    TH2D *h_QA_nTracks_HFT_TOF_1012 = static_cast<TH2D *>(mOutList->FindObject("h_QA_nTracks_HFT_TOF_1012"));
    TH2D *h_QA_nTracks_HFT_TOF_3040 = static_cast<TH2D *>(mOutList->FindObject("h_QA_nTracks_HFT_TOF_3040"));

    TH2D *h_QA_nTracks_HFT_0406 = static_cast<TH2D *>(mOutList->FindObject("h_QA_nTracks_HFT_0406"));
    TH2D *h_QA_nTracks_HFT_1012 = static_cast<TH2D *>(mOutList->FindObject("h_QA_nTracks_HFT_1012"));
    TH2D *h_QA_nTracks_HFT_3040 = static_cast<TH2D *>(mOutList->FindObject("h_QA_nTracks_HFT_3040"));

    TH2D *h_QA_nTracks_TOF_0406 = static_cast<TH2D *>(mOutList->FindObject("h_QA_nTracks_TOF_0406"));
    TH2D *h_QA_nTracks_TOF_1012 = static_cast<TH2D *>(mOutList->FindObject("h_QA_nTracks_TOF_1012"));
    TH2D *h_QA_nTracks_TOF_3040 = static_cast<TH2D *>(mOutList->FindObject("h_QA_nTracks_TOF_3040"));

    TH2D *h_QA_nTracks_0406 = static_cast<TH2D *>(mOutList->FindObject("h_QA_nTracks_0406"));
    TH2D *h_QA_nTracks_1012 = static_cast<TH2D *>(mOutList->FindObject("h_QA_nTracks_1012"));
    TH2D *h_QA_nTracks_3040 = static_cast<TH2D *>(mOutList->FindObject("h_QA_nTracks_3040"));

    TH2D *h_QA_nTracks_HFT_TOF = static_cast<TH2D *>(mOutList->FindObject("h_QA_nTracks_HFT_TOF"));
    TH2D *h_QA_pT_HFT_TOF = static_cast<TH2D *>(mOutList->FindObject("h_QA_pT_HFT_TOF"));
    TH2D *h_QA_eta_HFT_TOF = static_cast<TH2D *>(mOutList->FindObject("h_QA_eta_HFT_TOF"));
    TH2D *h_QA_phi_HFT_TOF = static_cast<TH2D *>(mOutList->FindObject("h_QA_phi_HFT_TOF"));

    TH2D *h_QA_DCA_xy_HFT_TOF = static_cast<TH2D *>(mOutList->FindObject("h_QA_DCA_xy_HFT_TOF"));
    TH2D *h_QA_DCA_xy_zoom_HFT_TOF = static_cast<TH2D *>(mOutList->FindObject("h_QA_DCA_xy_zoom_HFT_TOF"));
    TH2D *h_QA_DCA_z_HFT_TOF = static_cast<TH2D *>(mOutList->FindObject("h_QA_DCA_z_HFT_TOF"));
    TH2D *h_QA_DCA_z_zoom_HFT_TOF = static_cast<TH2D *>(mOutList->FindObject("h_QA_DCA_z_zoom_HFT_TOF"));

    //Number of tracks vs. pT - get from integral of pT spectrum

//____general QA histograms (all tracks within (TPC, track quality) cuts, NO TOF and HFT)__________________________________
    TH2F *h_QA_pT = static_cast<TH2F *>(mOutList->FindObject("h_QA_pT"));
    TH2F *h_QA_eta = static_cast<TH2F *>(mOutList->FindObject("h_QA_eta"));
    TH2F *h_QA_phi = static_cast<TH2F *>(mOutList->FindObject("h_QA_phi"));

    TH2F *h_QA_vertex_x = static_cast<TH2F *>(mOutList->FindObject("h_QA_vertex_x"));
    TH2F *h_QA_vertex_y = static_cast<TH2F *>(mOutList->FindObject("h_QA_vertex_y"));
    TH2F *h_QA_vertex_z = static_cast<TH2F *>(mOutList->FindObject("h_QA_vertex_z"));

    TH2F *h_QA_DCA_xy_TPC = static_cast<TH2F *>(mOutList->FindObject("h_QA_DCA_xy_TPC"));
    TH2F *h_QA_DCA_xy_zoom_TPC = static_cast<TH2F *>(mOutList->FindObject("h_QA_DCA_xy_zoom_TPC"));
    TH2F *h_QA_DCA_z_TPC = static_cast<TH2F *>(mOutList->FindObject("h_QA_DCA_z_TPC"));
    TH2F *h_QA_DCA_TPC = static_cast<TH2F *>(mOutList->FindObject("h_QA_DCA_TPC"));
    TH2F *h_QA_DCA_z_zoom_TPC = static_cast<TH2F *>(mOutList->FindObject("h_QA_DCA_z_zoom_TPC"));

//___HFT QA hitograms______________________________________________________________________________________________________
    TH2F *h_QA_pT_HFT = static_cast<TH2F *>(mOutList->FindObject("h_QA_pT_HFT"));
    TH2F *h_QA_eta_HFT = static_cast<TH2F *>(mOutList->FindObject("h_QA_eta_HFT"));
    TH2F *h_QA_phi_HFT = static_cast<TH2F *>(mOutList->FindObject("h_QA_phi_HFT"));

    TH2F *h_QA_vertex_x_HFT = static_cast<TH2F *>(mOutList->FindObject("h_QA_vertex_x_HFT"));
    TH2F *h_QA_vertex_y_HFT = static_cast<TH2F *>(mOutList->FindObject("h_QA_vertex_y_HFT"));
    TH2F *h_QA_vertex_z_HFT = static_cast<TH2F *>(mOutList->FindObject("h_QA_vertex_z_HFT"));

    TH2F *h_QA_DCA_xy_HFT = static_cast<TH2F *>(mOutList->FindObject("h_QA_DCA_xy_HFT"));
    TH2F *h_QA_DCA_xy_zoom_HFT = static_cast<TH2F *>(mOutList->FindObject("h_QA_DCA_xy_zoom_HFT"));
    TH2F *h_QA_DCA_z_HFT = static_cast<TH2F *>(mOutList->FindObject("h_QA_DCA_z_HFT"));
    TH2F *h_QA_DCA_HFT_TOF = static_cast<TH2F *>(mOutList->FindObject("h_QA_DCA_HFT_TOF"));
    TH2F *h_QA_DCA_HFT = static_cast<TH2F *>(mOutList->FindObject("h_QA_DCA_HFT"));
    TH2F *h_QA_DCA_z_zoom_HFT = static_cast<TH2F *>(mOutList->FindObject("h_QA_DCA_z_zoom_HFT"));

    TH2F *h_QA_DCA_TOF = static_cast<TH2F *>(mOutList->FindObject("h_QA_DCA_TOF"));

//___HFT QA hitograms, eta cut______________________________________________________________________________________________________
    TH2D *h_QA_nTracks_TPC_etaCut = static_cast<TH2D *>(mOutList->FindObject("h_QA_nTracks_TPC_etaCut"));

    TH2F *h_QA_pT_HFT_etaCut = static_cast<TH2F *>(mOutList->FindObject("h_QA_pT_HFT_etaCut"));
    TH2F *h_QA_eta_HFT_etaCut = static_cast<TH2F *>(mOutList->FindObject("h_QA_eta_HFT_etaCut"));
    TH2F *h_QA_phi_HFT_etaCut = static_cast<TH2F *>(mOutList->FindObject("h_QA_phi_HFT_etaCut"));

    TH2F *h_QA_vertex_x_HFT_etaCut = static_cast<TH2F *>(mOutList->FindObject("h_QA_vertex_x_HFT_etaCut"));
    TH2F *h_QA_vertex_y_HFT_etaCut = static_cast<TH2F *>(mOutList->FindObject("h_QA_vertex_y_HFT_etaCut"));
    TH2F *h_QA_vertex_z_HFT_etaCut = static_cast<TH2F *>(mOutList->FindObject("h_QA_vertex_z_HFT_etaCut"));

    TH2F *h_QA_DCA_xy_HFT_etaCut = static_cast<TH2F *>(mOutList->FindObject("h_QA_DCA_xy_HFT_etaCut"));
    TH2F *h_QA_DCA_xy_zoom_HFT_etaCut = static_cast<TH2F *>(mOutList->FindObject("h_QA_DCA_xy_zoom_HFT_etaCut"));
    TH2F *h_QA_DCA_z_HFT_etaCut = static_cast<TH2F *>(mOutList->FindObject("h_QA_DCA_z_HFT_etaCut"));
    TH2F *h_QA_DCA_z_zoom_HFT_etaCut = static_cast<TH2F *>(mOutList->FindObject("h_QA_DCA_z_zoom_HFT_etaCut"));

//____TOF QA histograms______________________________________________________________________________________________________
    TH2F *h_QA_pT_TOF = static_cast<TH2F *>(mOutList->FindObject("h_QA_pT_TOF"));
    TH2F *h_QA_eta_TOF = static_cast<TH2F *>(mOutList->FindObject("h_QA_eta_TOF"));
    TH2F *h_QA_phi_TOF = static_cast<TH2F *>(mOutList->FindObject("h_QA_phi_TOF"));

    TH2F *h_QA_vertex_x_TOF = static_cast<TH2F *>(mOutList->FindObject("h_QA_vertex_x_TOF"));
    TH2F *h_QA_vertex_y_TOF = static_cast<TH2F *>(mOutList->FindObject("h_QA_vertex_y_TOF"));
    TH2F *h_QA_vertex_z_TOF = static_cast<TH2F *>(mOutList->FindObject("h_QA_vertex_z_TOF"));

    TVector3 pVtx = mPicoDst->event()->primaryVertex();
    float b = mPicoDst->event()->bField();
    float grefMult = mPicoDst->event()->refMult();
    h_gRefmult->Fill(grefMult, RunIndex);

    h_QA_nEvents->Fill(RunIndex); //number of 0 filled to this histogram = number of events

    h_QA_ZDC_rate->Fill(mPicoDst->event()->ZDCx() / 1000., RunIndex);
    h_QA_BBC_rate->Fill(mPicoDst->event()->BBCx() / 1000., RunIndex);
    h_QA_ZDC_over_BBC->Fill(mPicoDst->event()->ZDCx()/mPicoDst->event()->BBCx(), RunIndex);

    float ZDC = mPicoDst->event()->ZDCx() / 1000.;
    float BBC = mPicoDst->event()->BBCx() / 1000.;

    float vertex_x_QA = pVtx.x();
    float vertex_y_QA = pVtx.y();
    float vertex_z_QA = pVtx.z();
    float nTrk=0, nTrk0406=0, nTrk1012=0, nTrk3040=0, nTrkHftTof=0, nTrkHftTof0406=0, nTrkHftTof1012=0, nTrkHftTof3040=0, nTrkHft=0, nTrkHft0406=0, nTrkHft1012=0, nTrkHft3040=0, nTrkTof=0, nTrkTof0406=0, nTrkTof1012=0, nTrkTof3040=0;

    if (vertex_z_QA>-6 && vertex_z_QA<=-4) h_gRefmult_Vz_min6_min4->Fill(grefMult);
    if (vertex_z_QA>-4 && vertex_z_QA<=-2) h_gRefmult_Vz_min4_min2->Fill(grefMult);
    if (vertex_z_QA>-2 && vertex_z_QA<=0) h_gRefmult_Vz_min2_0->Fill(grefMult);
    if (vertex_z_QA>0 && vertex_z_QA<=2) h_gRefmult_Vz_0_2->Fill(grefMult);
    if (vertex_z_QA>2 && vertex_z_QA<=4) h_gRefmult_Vz_2_4->Fill(grefMult);
    if (vertex_z_QA>4 && vertex_z_QA<=6) h_gRefmult_Vz_4_6->Fill(grefMult);

    h_gRefmult_vs_ZDCx->Fill(ZDC, grefMult);
    h_gRefmult_vs_BBCx->Fill(BBC, grefMult);

    h_QA_Vz->Fill(vertex_z_QA, RunIndex);
    h_QA_VzmVzVPD->Fill(fabs(vertex_z_QA - mPicoDst->event()->vzVpd()), RunIndex);

//    StPicoKFVertexFitter kfVertexFitter;
//    KFVertex kfVertex = kfVertexFitter.primaryVertexRefit(mPicoDst);
//    h_QA_Vx_position->Fill(kfVertex.GetX(), pVtx.x());
//    h_QA_Vy_position->Fill(kfVertex.GetY(), pVtx.y());
//    h_QA_Vz_position->Fill(kfVertex.GetZ(), pVtx.z());


    Int_t nPions = 0;
    Int_t nPionsTOFMatching = 0;
    Int_t nKaons = 0;
    Int_t nKaonsTOFMatching = 0;
    Int_t nPionsHFT=0;
    Int_t nPionsHFTTOF=0;
    Int_t nPionsHFThybridTOF=0;
    Int_t nKaonsHFT=0;
    Int_t nKaonsHFTTOF=0;
    Int_t nKaonsHFThybridTOF=0;
    Int_t nCommon=0;

    for (unsigned short iTrack = 0; iTrack < mPicoDst->numberOfTracks(); ++iTrack) {
        StPicoTrack const *trk = mPicoDst->track(iTrack);
        if (!trk) continue;
        float pT_QA = trk->gPt();
        if (pT_QA < 0.15) continue;
        float Beta = mHFCuts->getTofBetaBase(trk);

        h_QA_nHitsFit->Fill(trk->nHitsFit(), RunIndex);
        h_QA_nSigmaKaon->Fill(trk->nSigmaKaon(), RunIndex);
        h_QA_nSigmaPion->Fill(trk->nSigmaPion(), RunIndex);
        h_QA_nHitsDedx->Fill(trk->nHitsDedx(),RunIndex);

        if (trk->nHitsFit() < 15) continue;

        if (mHFCuts->isTOFmatched(trk)) {
            if (abs(trk->nSigmaKaon()<3)) h_QA_OneOverBetaDiffKaon->Fill(mHFCuts->getOneOverBeta(trk,Beta,StPicoCutsBase::kKaon), RunIndex);
            if (abs(trk->nSigmaPion()<3)) h_QA_OneOverBetaDiffPion->Fill(mHFCuts->getOneOverBeta(trk,Beta,StPicoCutsBase::kPion), RunIndex);
            h_QA_Beta->Fill(Beta, RunIndex);
            h_QA_nHitsFitTOF->Fill(trk->nHitsFit(), RunIndex);
        }

        TVector3 momentum = trk->gMom();
        float eta_QA = momentum.PseudoRapidity();
        if (eta_QA>1) continue;
        StPicoPhysicalHelix helix = trk->helix(b);
        float dca_xy_QA = float(helix.geometricSignedDistance(pVtx.x(), pVtx.y()));
        float dca_z_QA = trk->origin().z() - vertex_z_QA;
        float dca_QA = (pVtx - trk->origin()).Mag();
        if (abs(dca_QA) > 1.5) continue;
        nTrk+=1;

        if (pT_QA>0.4 && pT_QA<0.6) nTrk0406+=1;
        if (pT_QA>1.0 && pT_QA<1.2) nTrk1012+=1;
        if (pT_QA>3 && pT_QA<4) nTrk3040+=1;

        float phi_QA = momentum.Phi();
        if (mHFCuts->isTOFmatched(trk)) {
            nTrkTof+=1;
            h_QA_ZDC_rate_pileUp_TOF->Fill(mPicoDst->event()->ZDCx() / 1000., RunIndex);
            h_QA_BBC_rate_pileUp_TOF->Fill(mPicoDst->event()->BBCx() / 1000., RunIndex);
            if (pT_QA>0.4 && pT_QA<0.6) {
                nTrkTof0406+=1;
                h_QA_ZDC_rate_TOF_0406->Fill(mPicoDst->event()->ZDCx() / 1000., RunIndex);
                h_QA_BBC_rate_TOF_0406->Fill(mPicoDst->event()->BBCx() / 1000., RunIndex);
            }
            if (pT_QA>1.0 && pT_QA<1.2) {
                nTrkTof1012+=1;
                h_QA_ZDC_rate_TOF_1012->Fill(mPicoDst->event()->ZDCx() / 1000., RunIndex);
                h_QA_BBC_rate_TOF_1012->Fill(mPicoDst->event()->BBCx() / 1000., RunIndex);
            }
            if (pT_QA>3 && pT_QA<4) {
                nTrkTof3040+=1;
                h_QA_ZDC_rate_TOF_3040->Fill(mPicoDst->event()->ZDCx() / 1000., RunIndex);
                h_QA_BBC_rate_TOF_3040->Fill(mPicoDst->event()->BBCx() / 1000., RunIndex);
            }

            if (trk->isHFTTrack()){
                nTrkHftTof+=1;
                h_QA_pT_HFT_TOF->Fill(pT_QA, RunIndex);
                h_QA_eta_HFT_TOF->Fill(eta_QA, RunIndex);
                h_QA_phi_HFT_TOF->Fill(phi_QA, RunIndex);

                h_QA_DCA_xy_HFT_TOF->Fill(dca_xy_QA, RunIndex);
                h_QA_DCA_xy_zoom_HFT_TOF->Fill(dca_xy_QA, RunIndex);
                h_QA_DCA_z_HFT_TOF->Fill(dca_z_QA, RunIndex);
                h_QA_DCA_z_zoom_HFT_TOF->Fill(dca_z_QA, RunIndex);
                h_QA_DCA_HFT_TOF->Fill(dca_QA, RunIndex);

                h_QA_ZDC_rate_TOF_HFT->Fill(mPicoDst->event()->ZDCx() / 1000., RunIndex);
                h_QA_BBC_rate_TOF_HFT->Fill(mPicoDst->event()->BBCx() / 1000., RunIndex);
                if (pT_QA>0.4 && pT_QA<0.6) {
                    nTrkHftTof0406+=1;
                }
                if (pT_QA>1.0 && pT_QA<1.2) {
                    nTrkHftTof1012+=1;
                    h_QA_ZDC_rate_TOF_HFT_1012->Fill(mPicoDst->event()->ZDCx() / 1000., RunIndex);
                    h_QA_BBC_rate_TOF_HFT_1012->Fill(mPicoDst->event()->BBCx() / 1000., RunIndex);
                }
                if (pT_QA>3 && pT_QA<4) {
                    nTrkHftTof3040+=1;
                    h_QA_ZDC_rate_TOF_HFT_3040->Fill(mPicoDst->event()->ZDCx() / 1000., RunIndex);
                    h_QA_BBC_rate_TOF_HFT_3040->Fill(mPicoDst->event()->BBCx() / 1000., RunIndex);
                }
            }
        }

        //add all relevant TPC and track quality cuts here
        h_QA_ZDC_rate_pileUp->Fill(mPicoDst->event()->ZDCx() / 1000., RunIndex); // TPC only
        h_QA_BBC_rate_pileUp->Fill(mPicoDst->event()->BBCx() / 1000., RunIndex);

        //Fill all TPC histograms here
        h_QA_pT->Fill(pT_QA, RunIndex);
        h_QA_eta->Fill(eta_QA, RunIndex);
        h_QA_phi->Fill(phi_QA, RunIndex);

        h_QA_vertex_x->Fill(vertex_x_QA, RunIndex);
        h_QA_vertex_y->Fill(vertex_y_QA, RunIndex);
        h_QA_vertex_z->Fill(vertex_z_QA, RunIndex);

        h_QA_DCA_xy_TPC->Fill(dca_xy_QA, RunIndex);
        h_QA_DCA_xy_zoom_TPC->Fill(dca_xy_QA, RunIndex);
        h_QA_DCA_z_TPC->Fill(dca_z_QA, RunIndex);
        h_QA_DCA_z_zoom_TPC->Fill(dca_z_QA, RunIndex);
        h_QA_DCA_TPC->Fill(dca_QA, RunIndex);

        //---------------FILL HFT QA INFORMATION-------------------------------------------------------
        if (trk->isHFTTrack()) {
            if (pT_QA>0.4 && pT_QA<0.6) nTrkHft0406+=1;
            if (pT_QA>1.0 && pT_QA<1.2) nTrkHft1012+=1;
            if (pT_QA>3 && pT_QA<4) nTrkHft3040+=1;

            //Fill all HFT histograms here
            nTrkHft+=1;
            h_QA_pT_HFT->Fill(pT_QA, RunIndex);
            h_QA_eta_HFT->Fill(eta_QA, RunIndex);
            h_QA_phi_HFT->Fill(phi_QA, RunIndex);

            h_QA_vertex_x_HFT->Fill(vertex_x_QA, RunIndex);
            h_QA_vertex_y_HFT->Fill(vertex_y_QA, RunIndex);
            h_QA_vertex_z_HFT->Fill(vertex_z_QA, RunIndex);

            h_QA_DCA_xy_HFT->Fill(dca_xy_QA, RunIndex);
            h_QA_DCA_xy_zoom_HFT->Fill(dca_xy_QA, RunIndex);
            h_QA_DCA_z_HFT->Fill(dca_z_QA, RunIndex);
            h_QA_DCA_z_zoom_HFT->Fill(dca_z_QA, RunIndex);
            h_QA_DCA_HFT->Fill(dca_QA, RunIndex);
        }
        float nHFTHits = 0;

        if (trk->hasPxl1Hit()) {
            //Fill all PXL1 histograms here
            nHFTHits+=1;
            h_QA_nTracks_HFT_PXL1->Fill(pT_QA, RunIndex);
        }

        if (trk->hasPxl2Hit()) {
            nHFTHits+=1;
            //Fill all PXL2 histograms here
            h_QA_nTracks_HFT_PXL2->Fill(pT_QA, RunIndex);
        }

        if (trk->hasIstHit()) {
            nHFTHits+=1;
            //Fill all IST histograms here
            h_QA_nTracks_HFT_IST->Fill(pT_QA, RunIndex);
        }

        if (trk->hasSstHit()) {
            nHFTHits+=1;
            //Fill all SSD histograms here
            h_QA_nTracks_HFT_SSD->Fill(pT_QA, RunIndex);
        }

        if ((trk->hasIstHit() || trk->hasSstHit())) {
            //Fill all (IST or SSD) histograms here
            h_QA_nTracks_HFT_IST_or_SSD->Fill(pT_QA, RunIndex);
        }

        h_QA_nHFTHits->Fill(nHFTHits,RunIndex);

        if (eta_QA < 0.4) {
            h_QA_nTracks_TPC_etaCut->Fill(pT_QA, RunIndex); //for HFT_etaCut histograms
        }

        if (trk->isHFTTrack() && eta_QA <= 0.4) {
            h_QA_nTracks_HFT_etaCut->Fill(pT_QA, RunIndex);

            h_QA_pT_HFT_etaCut->Fill(pT_QA, RunIndex);
            h_QA_eta_HFT_etaCut->Fill(eta_QA, RunIndex);
            h_QA_phi_HFT_etaCut->Fill(phi_QA, RunIndex);

            h_QA_vertex_x_HFT_etaCut->Fill(vertex_x_QA, RunIndex);
            h_QA_vertex_y_HFT_etaCut->Fill(vertex_y_QA, RunIndex);
            h_QA_vertex_z_HFT_etaCut->Fill(vertex_z_QA, RunIndex);

            h_QA_DCA_xy_HFT_etaCut->Fill(dca_xy_QA, RunIndex);
            h_QA_DCA_xy_zoom_HFT_etaCut->Fill(dca_xy_QA, RunIndex);
            h_QA_DCA_z_HFT_etaCut->Fill(dca_z_QA, RunIndex);
            h_QA_DCA_z_zoom_HFT_etaCut->Fill(dca_z_QA, RunIndex);
        }

        if (trk->hasPxl1Hit() && eta_QA <= 0.4) {
            h_QA_nTracks_HFT_etaCut_PXL1->Fill(pT_QA, RunIndex);
        }

        if (trk->hasPxl2Hit() && eta_QA <= 0.4) {
            h_QA_nTracks_HFT_etaCut_PXL2->Fill(pT_QA, RunIndex);
        }

        if (trk->hasIstHit() && eta_QA <= 0.4) {
            h_QA_nTracks_HFT_etaCut_IST->Fill(pT_QA, RunIndex);
        }

        if (trk->hasSstHit() && eta_QA <= 0.4) {
            h_QA_nTracks_HFT_etaCut_SSD->Fill(pT_QA, RunIndex);
        }

        if ((trk->hasIstHit() || trk->hasSstHit()) && eta_QA <= 0.4) {
            h_QA_nTracks_HFT_etaCut_IST_or_SSD->Fill(pT_QA, RunIndex);
        }

        if (mHFCuts->isTOFmatched(trk)) {
            //Fill all TOF histograms here
            h_QA_pT_TOF->Fill(pT_QA, RunIndex);
            h_QA_eta_TOF->Fill(eta_QA, RunIndex);
            h_QA_phi_TOF->Fill(phi_QA, RunIndex);

            h_QA_vertex_x_TOF->Fill(vertex_x_QA, RunIndex);
            h_QA_vertex_y_TOF->Fill(vertex_y_QA, RunIndex);
            h_QA_vertex_z_TOF->Fill(vertex_z_QA, RunIndex);

            h_QA_DCA_TOF->Fill(dca_QA, RunIndex);

            h_QA_nHFTHitsTOF->Fill(nHFTHits,RunIndex);

        }

        bool isCommon=false;

        if (fabs(trk->nSigmaKaon())<3) {
            nKaons = nKaons + 1;
            if (mHFCuts->isTOFmatched(trk)) {
                nKaonsTOFMatching = nKaonsTOFMatching+1;
            }
            if (trk->isHFTTrack()) nKaonsHFT=nKaonsHFT+1;
            if (trk->isHFTTrack() && mHFCuts->isTOFmatched(trk)) nKaonsHFTTOF=nKaonsHFTTOF+1;
            if (trk->isHFTTrack() && mHFCuts->isHybridTOFKaon(trk)) {
                nKaonsHFThybridTOF=nKaonsHFThybridTOF+1;
                isCommon=true;
            }
        }

        if (fabs(trk->nSigmaPion())<3) {
            nPions = nPions + 1;
            if (mHFCuts->isTOFmatched(trk)) {
                nPionsTOFMatching = nPionsTOFMatching+1;
            }
            if (trk->isHFTTrack()) nPionsHFT=nPionsHFT+1;
            if (trk->isHFTTrack() && mHFCuts->isTOFmatched(trk)) nPionsHFTTOF=nPionsHFTTOF+1;
            if (trk->isHFTTrack() && mHFCuts->isHybridTOFPion(trk)) {
                nPionsHFThybridTOF=nPionsHFThybridTOF+1;
                if (isCommon) isCommon=true;
                else isCommon=false;
            }
        }

        if (isCommon) ++nCommon;

    } // .. end tracks loop

//    ntp_event->Fill(nTrk,nTrk0406,nTrk1012,nTrk3040,nTrkHft,nTrkHft0406,nTrkHft1012,nTrkHft3040,ZDC,BBC,mPicoDst->event()->runId(),RunIndex);

    h_QA_nTracks->Fill(nTrk, RunIndex);
    h_QA_nTracks_TOF->Fill(nTrkTof, RunIndex);
    h_QA_nTracks_HFT->Fill(nTrkHft, RunIndex);
    h_QA_nTracks_HFT_TOF->Fill(nTrkHftTof, RunIndex);

    h_QA_nTracks_TOF_0406->Fill(nTrkTof0406, RunIndex);
    h_QA_nTracks_TOF_1012->Fill(nTrkTof1012, RunIndex);
    h_QA_nTracks_TOF_3040->Fill(nTrkTof3040, RunIndex);

    h_QA_nTracks_HFT_0406->Fill(nTrkHft0406, RunIndex);
    h_QA_nTracks_HFT_1012->Fill(nTrkHft1012, RunIndex);
    h_QA_nTracks_HFT_3040->Fill(nTrkHft3040, RunIndex);

    h_QA_nTracks_HFT_TOF_0406->Fill(nTrkHftTof0406, RunIndex);
    h_QA_nTracks_HFT_TOF_1012->Fill(nTrkHftTof1012, RunIndex);
    h_QA_nTracks_HFT_TOF_3040->Fill(nTrkHftTof3040, RunIndex);

    h_QA_nTracks_0406->Fill(nTrk0406, RunIndex);
    h_QA_nTracks_1012->Fill(nTrk1012, RunIndex);
    h_QA_nTracks_3040->Fill(nTrk3040, RunIndex);

    h_QA_BBC_rate_kaons->Fill(mPicoDst->event()->BBCx()/1000., nKaons);
    h_QA_BBC_rate_pions->Fill(mPicoDst->event()->BBCx()/1000., nPions);
    h_QA_BBC_rate_pions_matching->Fill(mPicoDst->event()->BBCx()/1000., nPionsTOFMatching);
    h_QA_BBC_rate_kaons_matching->Fill(mPicoDst->event()->BBCx()/1000., nKaonsTOFMatching);
    h_QA_BBC_rate_kaons_HFT->Fill(mPicoDst->event()->BBCx()/1000., nKaonsHFT);
    h_QA_BBC_rate_pions_HFT->Fill(mPicoDst->event()->BBCx()/1000., nPionsHFT);
    h_QA_BBC_rate_kaons_HFT_TOF->Fill(mPicoDst->event()->BBCx()/1000., nKaonsHFTTOF);
    h_QA_BBC_rate_pions_HFT_TOF->Fill(mPicoDst->event()->BBCx()/1000., nPionsHFTTOF);
    h_QA_BBC_rate_kaons_HFT_hybridTOF->Fill(mPicoDst->event()->BBCx()/1000., nKaonsHFThybridTOF);
    h_QA_BBC_rate_pions_HFT_hybridTOF->Fill(mPicoDst->event()->BBCx()/1000., nPionsHFThybridTOF);

    if(nTrkHft>1) {
        h_gRefmult_HFT->Fill(mPicoDst->event()->grefMult(), RunIndex);

        if (vertex_z_QA > -6 && vertex_z_QA <= -4) h_gRefmult_Vz_min6_min4_HFT->Fill(grefMult);
        if (vertex_z_QA > -4 && vertex_z_QA <= -2) h_gRefmult_Vz_min4_min2_HFT->Fill(grefMult);
        if (vertex_z_QA > -2 && vertex_z_QA <= 0) h_gRefmult_Vz_min2_0_HFT->Fill(grefMult);
        if (vertex_z_QA > 0 && vertex_z_QA <= 2) h_gRefmult_Vz_0_2_HFT->Fill(grefMult);
        if (vertex_z_QA > 2 && vertex_z_QA <= 4) h_gRefmult_Vz_2_4_HFT->Fill(grefMult);
        if (vertex_z_QA > 4 && vertex_z_QA <= 6) h_gRefmult_Vz_4_6_HFT->Fill(grefMult);

        h_gRefmult_vs_ZDCx_HFT->Fill(ZDC, grefMult);
        h_gRefmult_vs_BBCx_HFT->Fill(BBC, grefMult);
    }

    if((nCommon!=1 && nKaonsHFThybridTOF>0 && nPionsHFThybridTOF>0) || (nCommon==1 && nKaonsHFThybridTOF>1)) {
        h_gRefmult_HFT_hybridTOF->Fill(mPicoDst->event()->grefMult(), RunIndex);

        if (vertex_z_QA > -6 && vertex_z_QA <= -4) h_gRefmult_Vz_min6_min4_HFT_hybridTOF->Fill(grefMult);
        if (vertex_z_QA > -4 && vertex_z_QA <= -2) h_gRefmult_Vz_min4_min2_HFT_hybridTOF->Fill(grefMult);
        if (vertex_z_QA > -2 && vertex_z_QA <= 0) h_gRefmult_Vz_min2_0_HFT_hybridTOF->Fill(grefMult);
        if (vertex_z_QA > 0 && vertex_z_QA <= 2) h_gRefmult_Vz_0_2_HFT_hybridTOF->Fill(grefMult);
        if (vertex_z_QA > 2 && vertex_z_QA <= 4) h_gRefmult_Vz_2_4_HFT_hybridTOF->Fill(grefMult);
        if (vertex_z_QA > 4 && vertex_z_QA <= 6) h_gRefmult_Vz_4_6_HFT_hybridTOF->Fill(grefMult);

        h_gRefmult_vs_ZDCx_HFT_hybridTOF->Fill(ZDC, grefMult);
        h_gRefmult_vs_BBCx_HFT_hybridTOF->Fill(BBC, grefMult);
    }


    return kStOK;
}
