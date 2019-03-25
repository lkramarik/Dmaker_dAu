#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoHFMaker/StHFCuts.h"
#include "phys_constants.h"
#include "StPicoKFVertexTools.h"
#include "StiMaker/StKFVerticesCollection.h"
#include "StPicoKFVertexFitter.h"
ClassImp(StPicoKFVertexTools)

// _________________________________________________________
StPicoKFVertexTools::StPicoKFVertexTools(char const* name, StPicoDstMaker* picoMaker, char const* outputBaseFileName) :
        StPicoHFMaker(name, picoMaker, outputBaseFileName),
        mOutFileBaseName(outputBaseFileName){
}

// _________________________________________________________
StPicoKFVertexTools::~StPicoKFVertexTools() {
}

// _________________________________________________________
int StPicoKFVertexTools::InitHF() {
    mOutFileBaseName = mOutFileBaseName.ReplaceAll(".root", "");

    mOutList->Add(new TH1F("hMassUS","hMassUS", 500, 1.6, 2.1));

    ntp_vertex = new TNtuple("ntp_vertex","ntp_vertex","runId:refMult:nGlobTracks:nHftTracks:nD0:StAnnelingChi2Cut:"
                                                       "picoDstVx:picoDstVy:picoDstVz:"
                                                       "picoDstVErrX:picoDstVErrY:picoDstVErrZ:"
                                                       "KFVx:KFVy:KFVz:"
                                                       "KFVErrX:KFVErrY:KFVErrZ");
    return kStOK;
}

// _________________________________________________________
void StPicoKFVertexTools::ClearHF(Option_t *opt="") {
    return;
}

// _________________________________________________________
int StPicoKFVertexTools::FinishHF() {
    ntp_vertex -> Write(ntp_vertex->GetName(), TObject::kOverwrite);
    return kStOK;
}
// _________________________________________________________
int StPicoKFVertexTools::MakeHF() {
    TH1F *hMassUS = static_cast<TH1F*>(mOutList->FindObject("hMassUS"));

    std::vector<unsigned short> mIdxPicoPions;
    std::vector<unsigned short> mIdxPicoKaons;

    int nHftTracks=0;
    int nD0=0;

    // k and pi arrays:
    for (unsigned short iTrack = 0; iTrack < nTracks; ++iTrack) {
        StPicoTrack* trk = mPicoDst->track(iTrack);
        if (trk->isHFTTrack()) nHftTracks++;
        if (abs(trk->gMom().PseudoRapidity())>1) continue;
        if (mHFCuts->isGoodPion(trk)) mIdxPicoPions.push_back(iTrack);
        if (mHFCuts->isGoodKaon(trk)) mIdxPicoKaons.push_back(iTrack);
    }

    std::vector<int> tracksToRemove;

    //pair reconstruction
    for (unsigned short idxPion1 = 0; idxPion1 < mIdxPicoPions.size(); ++idxPion1) {
        StPicoTrack const *pion1 = mPicoDst->track(mIdxPicoPions[idxPion1]);

        for (unsigned short idxKaon = 0; idxKaon < mIdxPicoKaons.size(); ++idxKaon) {
            StPicoTrack const *kaon = mPicoDst->track(mIdxPicoKaons[idxKaon]);

            if((kaon->charge() + pion1->charge() != 0) ) continue;

            StHFPair *pair = new StHFPair(pion1, kaon, mHFCuts->getHypotheticalMass(StPicoCutsBase::kPion), mHFCuts->getHypotheticalMass(StPicoCutsBase::kKaon), mIdxPicoPions[idxPion1], mIdxPicoKaons[idxKaon], mPrimVtx, mBField, kTRUE);

            if ((pair->m()<1.7) || (pair->m()>2)) continue;

            if (mHFCuts->isGoodSecondaryVertexPairPtBin(pair)) {
                nD0+=1;
                tracksToRemove.push_back(mIdxPicoPions[idxPion1]);
                tracksToRemove.push_back(mIdxPicoKaons[idxKaon]);
                hMassUS->Fill(pair->m());
            }
        }
    }



    bool goodEvent=true;

//    if (!(nHftTracks>1)) goodEvent=false;
//    if (!(mPicoEvent->numberOfPxlInnerHits()>0 && mPicoEvent->numberOfPxlOuterHits()>0)) goodEvent=false;
//    if (!(mPicoEvent->BBCx()<950000)) goodEvent=false;

    if (nD0>0) {
        StPicoKFVertexFitter kfVertexFitter;
        KFVertex kfVertex = kfVertexFitter.primaryVertexRefit(mPicoDst, tracksToRemove);
//        KFVertex kfVertex = kfVertexFitter.primaryVertexRefit(mPicoDst);

        const int nNtVars = ntp_vertex->GetNvar();
        Float_t ntVar[nNtVars];
        int ii = 0;

        ntVar[ii++] = mPicoEvent->runId();
        ntVar[ii++] = mPicoEvent->refMult();
        ntVar[ii++] = mPicoEvent->numberOfGlobalTracks();
        ntVar[ii++] = nHftTracks;
        ntVar[ii++] = nD0;
        ntVar[ii++] = StAnneling::Chi2Cut();

        ntVar[ii++] = mPrimVtx.x();
        ntVar[ii++] = mPrimVtx.y();
        ntVar[ii++] = mPrimVtx.z();

        ntVar[ii++] = (mPicoEvent->primaryVertexError()).x();
        ntVar[ii++] = (mPicoEvent->primaryVertexError()).y();
        ntVar[ii++] = (mPicoEvent->primaryVertexError()).z();

        ntVar[ii++] = kfVertex.GetX();
        ntVar[ii++] = kfVertex.GetY();
        ntVar[ii++] = kfVertex.GetZ();

        ntVar[ii++] = kfVertex.GetErrX();
        ntVar[ii++] = kfVertex.GetErrY();
        ntVar[ii++] = kfVertex.GetErrZ();

        ntp_vertex->Fill(ntVar);
    }

    mIdxPicoPions.clear();
    mIdxPicoKaons.clear();
    tracksToRemove.clear();

    return kStOK;
}

