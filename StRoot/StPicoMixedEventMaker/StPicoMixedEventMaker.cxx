#include "TTree.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TList.h"
#include "TString.h" // needed for the Form(...)

#include "StarClassLibrary/StThreeVectorF.hh"
#include "StPicoDstMaker/StPicoDst.h"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoDstMaker/StPicoBTofPidTraits.h"
#include "StPicoMixedEventMaker.h"
#include "StPicoEventMixer.h"
#include "StPicoHFMaker/StHFCuts.h"

#include <vector>
#include <string>

ClassImp(StPicoMixedEventMaker)

// _________________________________________________________
StPicoMixedEventMaker::StPicoMixedEventMaker(char const* name, StPicoDstMaker* picoMaker, StHFCuts* hfCuts,
                                             char const* outputBaseFileName,  char const* inputHFListHFtree = "") :
        StMaker(name),
        mPicoDst(NULL),
        mPicoDstMaker(picoMaker),
        mPicoEvent(NULL),
        mHFCuts(hfCuts),
        mOuputFileBaseName(outputBaseFileName),
        mInputFileName(inputHFListHFtree),
        mEventCounter(0),
        mBufferSize(defaultBufferSize),
        mSETupleSig(NULL),
        mSETupleBack(NULL),
        mMETupleSig(NULL),
        mMETupleBack(NULL),
        mOutputFileTree(NULL),
        mSingePartHists(NULL)
{

    TH1::AddDirectory(false);
    // -- create OutputTree
    mOutputFileTree = new TFile(Form("%s.picoMEtree.root", mOuputFileBaseName.Data()), "RECREATE");
    mOutputFileTree->SetCompressionLevel(1);
    mOutputFileTree->cd();

    const string varList ="pi1_pt:pi1_dca:pi1_nSigma:pi1_nHitFit:pi1_TOFinvbeta:"
                          "k_pt:k_dca:k_nSigma:k_nHitFit:k_TOFinvbeta:"
                          "dcaDaughters:D_rapidity:D_theta:cosTheta:D_decayL:dcaD0ToPv:D_phi:D_eta:D_cosThetaStar:D_pt:D_mass";

    mSETupleSig = new TNtuple("ntp_signal_SE","ntp_signal_SE", varList.data());
    mMETupleSig = new TNtuple("ntp_signal_ME","ntp_signal_ME", varList.data());
    mSETupleBack = new TNtuple("ntp_background_SE","ntp_background_SE", varList.data());
    mMETupleBack = new TNtuple("ntp_background_ME","ntp_background_ME", varList.data());

    mSingePartHists = new TList();
    mSingePartHists->SetOwner(true);
    mSingePartHists->SetName("HFSinglePartHists");

    // create single particle hists
    const std::string evtNames[2] = {"SE", "ME"};
    const std::string partNames[3] = {"pi", "K"};

    for (int i = 0; i < 2; ++i)
    {
        // mSingePartHists->Add(new TH1D(Form("centrality%s", evtNames[i].data()),Form("centrality%s", evtNames[i].data()), 10, -1.5, 8.5));
        // mSingePartHists->Add(new TH1D(Form("centralityCorrection%s", evtNames[i].data()),Form("centrality corrected %s", evtNames[i].data()), 10, -1.5, 8.5));
        // mSingePartHists->Add(new TH1D(Form("refMult%s", evtNames[i].data()), Form("corrected refferernce multiplicity %s", evtNames[i].data()), 100, 0, 800));

        for (int iPart = 0; iPart < 2; ++iPart)
        {
            // eta phi
            mSingePartHists->Add(new TH2D(Form("%sEtaPhi%s",partNames[ iPart ].data(), evtNames[i].data()),
                                          Form("%s Eta phi distribution %s",partNames[ iPart ].data(), evtNames[i].data()),
                                          100, -TMath::Pi(), TMath::Pi(), 100, -1.1, 1.1));

            // phi vs pT
            mSingePartHists->Add(new TH2D(Form("%sPhiPt%s",partNames[ iPart ].data(), evtNames[i].data()),
                                          Form("%s phi vs pT %s",partNames[ iPart ].data(), evtNames[i].data()),
                                          100, 0, 15, 100, -TMath::Pi(), TMath::Pi()));

            // DCA
            mSingePartHists->Add(new TH1D(Form("%sDCA%s",partNames[ iPart ].data(), evtNames[i].data()),
                                          Form("%s DCA %s",partNames[ iPart ].data(), evtNames[i].data()),
                                          200, 0, 0.02));
            // nTracks
            mSingePartHists->Add(new TH1D(Form("%stracks%s",partNames[ iPart ].data(), evtNames[i].data()),
                                          Form("Number of %s tracks %s",partNames[ iPart ].data(), evtNames[i].data()),
                                          100, -0.5, 99.5));
        }
    }

    // loop over all histograms to set Sumw2
    TH1* hist = static_cast<TH1*>(mSingePartHists->First());
    hist->Sumw2();
    while(hist != static_cast<TH1*>(mSingePartHists->Last()))
    {
        hist = static_cast<TH1*>(mSingePartHists->After(hist));
        hist->Sumw2();
    }

}

// _________________________________________________________
StPicoMixedEventMaker::~StPicoMixedEventMaker() {

    for(int iVz =0 ; iVz < 10 ; ++iVz){
        for(int iCentrality = 0 ; iCentrality < 9 ; ++iCentrality){
            delete mPicoEventMixer[iVz][iCentrality];
        }
    }
    mOutputFileTree->Close();
}
// _________________________________________________________
bool StPicoMixedEventMaker::loadEventPlaneCorr(Int_t const run) {
    //needs to implement, will currently break maker
    return false;
}
// _________________________________________________________
Int_t StPicoMixedEventMaker::Init() {
    mOutputFileTree->cd();
    for(int iVz =0 ; iVz < 10 ; ++iVz){
        for(int iCentrality = 0 ; iCentrality < 9 ; ++iCentrality){
            mPicoEventMixer[iVz][iCentrality] = new StPicoEventMixer(Form("Cent_%i_Vz_%i",iCentrality,iVz));
            mPicoEventMixer[iVz][iCentrality]->setEventBuffer(mBufferSize);
            mPicoEventMixer[iVz][iCentrality]->setHFCuts(mHFCuts);
            mPicoEventMixer[iVz][iCentrality]->setSameEvtNtupleSig(mSETupleSig);
            mPicoEventMixer[iVz][iCentrality]->setSameEvtNtupleBack(mSETupleBack);
            mPicoEventMixer[iVz][iCentrality]->setMixedEvtNtupleSig(mMETupleSig);
            mPicoEventMixer[iVz][iCentrality]->setMixedEvtNtupleBack(mMETupleBack);
            mPicoEventMixer[iVz][iCentrality]->setSinglePartHistsList(mSingePartHists);
            mPicoEventMixer[iVz][iCentrality]->setFillSinglePartHists(fillSingleTrackHistos);
        }
    }

    //resetEvent();
    return kStOK;
}

// _________________________________________________________
Int_t StPicoMixedEventMaker::Finish() {
    mOutputFileTree->cd();

    for(int iVz =0 ; iVz < 10 ; ++iVz){
        for(int iCentrality = 0 ; iCentrality < 9 ; ++iCentrality){
            mPicoEventMixer[iVz][iCentrality]->finish();
            //delete mPicoEventMixer[iVz][iCentrality];
        }
    }
    mSingePartHists->Write();

    mSETupleSig -> Write(mSETupleSig->GetName(), TObject::kOverwrite);
    mMETupleSig -> Write(mMETupleSig->GetName(), TObject::kOverwrite);
    mSETupleBack -> Write(mSETupleBack->GetName(), TObject::kOverwrite);
    mMETupleBack -> Write(mMETupleBack->GetName(), TObject::kOverwrite);

    return kStOK;
}
// _________________________________________________________
void StPicoMixedEventMaker::Clear(Option_t* opt) {
}
// _________________________________________________________
Int_t StPicoMixedEventMaker::Make() {

    if(!mPicoDstMaker) {
        LOG_WARN << "No PicoDstMaker! Skipping! "<<endm;
        return kStWarn;
    }

    StPicoDst const* picoDst = mPicoDstMaker->picoDst();
    if (!picoDst) {
        LOG_WARN << "No picoDst ! Skipping! "<<endm;
        return kStWarn;
    }

    if (!mHFCuts->isGoodEvent(picoDst))
        return kStOk;
    StThreeVectorF const pVtx = picoDst->event()->primaryVertex();
    if( fabs(pVtx.z()) >=6.0 )
        return kStOk;

    int const centrality  = 1;

    if(centrality < 0 || centrality >8 ) return kStOk;
    int const vz_bin = (int)((6 +pVtx.z())/1.2) ;

//    if( mPicoEventMixer[vz_bin][centrality] -> addPicoEvent(picoDst, mGRefMultCorrUtil->getWeight()) ==  true )
//        mPicoEventMixer[vz_bin][centrality]->mixEvents();

    if( mPicoEventMixer[vz_bin][centrality] -> addPicoEvent(picoDst, 1 ))
        mPicoEventMixer[vz_bin][centrality]->mixEvents();

    return kStOk;
}
// _________________________________________________________
Int_t StPicoMixedEventMaker::SetCategories() {
    return kStOk;
}
// _________________________________________________________
int StPicoMixedEventMaker::categorize(StPicoDst const * picoDst ) {
    StThreeVectorF pVertex = (picoDst->event())->primaryVertex();
    if( fabs(pVertex.z())>6.0 ) return -99;
    int bin = -6.0 + (pVertex.z()+6.0)/1.2;
    return bin;
}
