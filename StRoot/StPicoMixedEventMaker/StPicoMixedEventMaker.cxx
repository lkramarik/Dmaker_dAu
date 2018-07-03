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
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StPicoMixedEventMaker.h"
#include "StPicoEventMixer.h"
#include "StPicoHFMaker/StHFCuts.h"
#include <vector>
#include <string>

ClassImp(StPicoMixedEventMaker)

static const int m_nmultEdge = 7;
static float constexpr m_multEdge[m_nmultEdge+1] = {0, 4, 8, 12, 16, 20, 24, 200};

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
        mBufferSize(5),
        mSETupleSig(NULL),
        mSETupleBack(NULL),
        mMETupleSig(NULL),
        mMETupleBack(NULL),
        mOutputFileTree(NULL)
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
}

// _________________________________________________________
StPicoMixedEventMaker::~StPicoMixedEventMaker() {

    for(int iVz =0 ; iVz < 10 ; ++iVz){
        for(int iCentrality = 0 ; iCentrality < m_nmultEdge ; ++iCentrality){
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
        for(int iCentrality = 0 ; iCentrality < m_nmultEdge ; ++iCentrality){
            mPicoEventMixer[iVz][iCentrality] = new StPicoEventMixer(Form("Cent_%i_Vz_%i",iCentrality,iVz));
            mPicoEventMixer[iVz][iCentrality]->setEventBuffer(mBufferSize);
            mPicoEventMixer[iVz][iCentrality]->setHFCuts(mHFCuts);
            mPicoEventMixer[iVz][iCentrality]->setSameEvtNtupleSig(mSETupleSig);
            mPicoEventMixer[iVz][iCentrality]->setSameEvtNtupleBack(mSETupleBack);
            mPicoEventMixer[iVz][iCentrality]->setMixedEvtNtupleSig(mMETupleSig);
            mPicoEventMixer[iVz][iCentrality]->setMixedEvtNtupleBack(mMETupleBack);
        }
    }

    //resetEvent();
    return kStOK;
}

// _________________________________________________________
Int_t StPicoMixedEventMaker::Finish() {

    mOutputFileTree->cd();

    for(int iVz =0 ; iVz < 10 ; ++iVz){
        for(int iCentrality = 0 ; iCentrality < m_nmultEdge ; ++iCentrality){
            mPicoEventMixer[iVz][iCentrality]->finish();
        }
    }

    mSETupleSig -> Write(mSETupleSig->GetName(), TObject::kOverwrite);
    mMETupleSig -> Write(mMETupleSig->GetName(), TObject::kOverwrite);
    mSETupleBack -> Write(mSETupleBack->GetName(), TObject::kOverwrite);
    mMETupleBack -> Write(mMETupleBack->GetName(), TObject::kOverwrite);

    mOutputFileTree->Close();

    return kStOK;
}
// _________________________________________________________
void StPicoMixedEventMaker::Clear(Option_t* opt) {
    return;
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

    int* aEventStat = NULL;
    if (!mHFCuts->isGoodEvent(picoDst, aEventStat))
        return kStOk;

    StThreeVectorF const pVtx = picoDst->event()->primaryVertex();


    int multiplicity = mPicoDst->event()->refMult();
    int centrality = getMultIndex(multiplicity);

    if(centrality < 0 || centrality > m_nmultEdge+1 ) return kStOk;
    int const vz_bin = (int)((6 +pVtx.z())/1.2) ;

    if( mPicoEventMixer[vz_bin][centrality] -> addPicoEvent(picoDst, 1)) {
        mPicoEventMixer[vz_bin][centrality] -> mixEvents();
    }
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
// _________________________________________________________
int StPicoMixedEventMaker::getMultIndex(float multiplicity){
    for (int i = 0; i < m_nmultEdge; i++){
        if ((multiplicity >= m_multEdge[i]) && (multiplicity < m_multEdge[i + 1]))
            return i;
    }
    return -1;
}
