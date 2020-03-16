#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoHFMaker/StHFCuts.h"
#include "phys_constants.h"
#include "StPicoKFVertexTools.h"
#include "StiMaker/StKFVerticesCollection.h"
#include "StPicoKFVertexFitter/StPicoKFVertexFitter.h"
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
    mOutList->Add(new TH1F("hMassUSRefit","hMassUSRefit", 500, 1.6, 2.1));
    mOutList->Add(new TH2F("hPicoPosition","hPicoPosition", 400, -1, 1, 400, 1, 1));

    ntp_vertex = new TNtuple("ntp_vertex","ntp_vertex","runId:refMult:grefMult:nGlobTracks:nHftTracks:nPrimTracks:nD0:StAnnelingChi2Cut:BBC:ZDC:VzVpd:"
                                                       "picoDstVx:picoDstVy:picoDstVz:"
                                                       "picoDstVErrX:picoDstVErrY:picoDstVErrZ:"
                                                       "KFVx:KFVy:KFVz:"
                                                       "KFVErrX:KFVErrY:KFVErrZ:nFlag0");

    ntp_KFReso = new TNtuple("ntp_KFReso","ntp_KFReso","KF1x:KF1y:KF1z:"
                                                       "KF2x:KF2y:KF2z:"
                                                       "KFdiffx:KFdiffy:KFdiffz:"
                                                       "KFdiff:nPrimTracks:nHftTracks:grefMult");
    return kStOK;
}

// _________________________________________________________
void StPicoKFVertexTools::ClearHF(Option_t *opt="") {
    return;
}

// _________________________________________________________
int StPicoKFVertexTools::FinishHF() {
    ntp_vertex -> Write(ntp_vertex->GetName(), TObject::kOverwrite);
    ntp_KFReso -> Write(ntp_KFReso->GetName(), TObject::kOverwrite);
    return kStOK;
}
// _________________________________________________________
int StPicoKFVertexTools::MakeHF() {
    bool goodEvent=true;

//    if (!(mPicoEvent->BBCx()<950000)) return kStOK;
//    if (!(mPrimVtx.x()>-0.28)) return kStOK;
//    if (!(mPrimVtx.x()<-0.13)) return kStOK;
//    if (!(mPrimVtx.y()>-0.28)) return kStOK;
//    if (!(mPrimVtx.y()<-0.13)) return kStOK;

//    if (!(abs(mPrimVtx.x())<0.5)) return kStOK;
//    if (!(abs(mPrimVtx.y())<0.5)) return kStOK;

    TH1F *hMassUS = static_cast<TH1F*>(mOutList->FindObject("hMassUS"));
    TH1F *hMassUSRefit = static_cast<TH1F*>(mOutList->FindObject("hMassUSRefit"));

    TH2F *hPicoPosition = static_cast<TH2F*>(mOutList->FindObject("hPicoPosition"));
    hPicoPosition->Fill(mPrimVtx.x(), mPrimVtx.y());

    std::vector<unsigned short> mIdxPicoPions;
    std::vector<unsigned short> mIdxPicoKaons;
    std::vector<int> primaryTracks;
    std::vector<int> globalTracks;

    int nHftTracks=0;
    int nD0=0;

    // k and pi arrays:
    UInt_t nTracks = mPicoDst->numberOfTracks();
    for (unsigned short iTrack = 0; iTrack < nTracks; ++iTrack) {
        StPicoTrack* trk = mPicoDst->track(iTrack);
        globalTracks.push_back(iTrack);
        if (trk->isPrimary()) primaryTracks.push_back(iTrack);
        if (trk->isHFTTrack()) nHftTracks++;
        if (abs(trk->gMom().PseudoRapidity())>1) continue;
//        if (mHFCuts->isGoodPion(trk)) mIdxPicoPions.push_back(iTrack);
//        if (mHFCuts->isGoodKaon(trk)) mIdxPicoKaons.push_back(iTrack);
    }

//    if (!(nHftTracks>1)) return kStOK;

//    std::vector<int> tracksToRemove;
//
//    //pair reconstruction
//    for (unsigned short idxPion1 = 0; idxPion1 < mIdxPicoPions.size(); ++idxPion1) {
//        StPicoTrack const *pion1 = mPicoDst->track(mIdxPicoPions[idxPion1]);
//
//        for (unsigned short idxKaon = 0; idxKaon < mIdxPicoKaons.size(); ++idxKaon) {
//            StPicoTrack const *kaon = mPicoDst->track(mIdxPicoKaons[idxKaon]);
//
//            if((kaon->charge() + pion1->charge() != 0) ) continue;
//
//            StHFPair *pair = new StHFPair(pion1, kaon, mHFCuts->getHypotheticalMass(StPicoCutsBase::kPion), mHFCuts->getHypotheticalMass(StPicoCutsBase::kKaon), mIdxPicoPions[idxPion1], mIdxPicoKaons[idxKaon], mPrimVtx, mBField, kTRUE);
//
//            if ((pair->m()<1.7) || (pair->m()>2)) continue;
//
//            if (mHFCuts->isGoodSecondaryVertexPairPtBin(pair)) {
//                nD0+=1;
//                tracksToRemove.push_back(mIdxPicoPions[idxPion1]);
//                tracksToRemove.push_back(mIdxPicoKaons[idxKaon]);
//                hMassUS->Fill(pair->m());
//            }
//        }
//    }





    if (goodEvent) {
//    if (nD0>-1) {
        compareFitters(primaryTracks, nD0, nHftTracks);

        //making 2 vertices and comparing:
//        if (primaryTracks.size()>10) {
//            makeKFReso(primaryTracks, nHftTracks);
//            makeKFReso(primaryTracks, nHftTracks);
//        }
//        KFVertex kfVertex = kfVertexFitter.primaryVertexRefit(mPicoDst, tracksToRemove);

        //////////////////////////////////////////////////////////
//        TVector3 newKFVertex(-999., -999., -999.);
//
//        if (kfVertex.GetX()) {
//            newKFVertex.SetXYZ(kfVertex.GetX(), kfVertex.GetY(), kfVertex.GetZ());
//        }
//
//        //pair reconstruction with new KF vertex
//        for (unsigned short idxPion1 = 0; idxPion1 < mIdxPicoPions.size(); ++idxPion1) {
//            StPicoTrack const *pion1 = mPicoDst->track(mIdxPicoPions[idxPion1]);
//
//            for (unsigned short idxKaon = 0; idxKaon < mIdxPicoKaons.size(); ++idxKaon) {
//                StPicoTrack const *kaon = mPicoDst->track(mIdxPicoKaons[idxKaon]);
//                if ((kaon->charge() + pion1->charge() != 0)) continue;
//
//                StHFPair *pair = new StHFPair(pion1, kaon, mHFCuts->getHypotheticalMass(StPicoCutsBase::kPion), mHFCuts->getHypotheticalMass(StPicoCutsBase::kKaon), mIdxPicoPions[idxPion1], mIdxPicoKaons[idxKaon], newKFVertex, mBField, kTRUE);
//                if ((pair->m() < 1.7) || (pair->m() > 2)) continue;
//                if (mHFCuts->isGoodSecondaryVertexPairPtBin(pair)) {
//                    hMassUSRefit->Fill(pair->m());
//                }
//            }
//        }
        ///////////////////////////////////////////////////////////////

    }


    mIdxPicoPions.clear();
    mIdxPicoPions.shrink_to_fit();

    mIdxPicoKaons.clear();
    mIdxPicoKaons.shrink_to_fit();

//    tracksToRemove.clear();
    primaryTracks.clear();
    primaryTracks.shrink_to_fit();

    globalTracks.clear();
    globalTracks.shrink_to_fit();
    return kStOK;
}

// _____________________________________________________________________________
void StPicoKFVertexTools::compareFitters(std::vector<int>&  primaryTracks, int nD0, int nHftTracks) {
    const unsigned int nPrimTracks = primaryTracks.size();

//    StPicoKFVertexFitter kfVertexFitter;
//    KFVertex kfVertex = kfVertexFitter.primaryVertexRefit(mPicoDst);
//    int nFlag0Write =  kfVertexFitter.nFlag0;

    //if you dont want to reconstruct PV
    const Double_t parVertex[6] = {mPrimVtx.x(),mPrimVtx.y(),mPrimVtx.z(), 0, 0, 1000};
    const Double_t covVertex[21] = {1,
                                    0, 1,
                                    0, 0, 100,
                                    0, 0,   0, 100,
                                    0, 0,   0,   0, 100,
                                    0, 0,   0,   0,   0, 100};
    KFVertex kfVertex;
    kfVertex.Create(parVertex,covVertex, 0, 0);
    int nFlag0Write=0;

    const int nNtVars = ntp_vertex->GetNvar();
    Float_t ntVar[nNtVars];
    int ii = 0;

    ntVar[ii++] = mPicoEvent->runId();
    ntVar[ii++] = mPicoEvent->refMult();
    ntVar[ii++] = mPicoEvent->grefMult();
    ntVar[ii++] = mPicoEvent->numberOfGlobalTracks();
    ntVar[ii++] = nHftTracks;
    ntVar[ii++] = nPrimTracks;
    ntVar[ii++] = nD0;
    ntVar[ii++] = StAnneling::Chi2Cut();
    ntVar[ii++] = mPicoEvent->BBCx()/1000;
    ntVar[ii++] = mPicoEvent->ZDCx()/1000;
    ntVar[ii++] = mPicoEvent->vzVpd();

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

    ntVar[ii++] = nFlag0Write;

    ntp_vertex->Fill(ntVar);
}

// _____________________________________________________________________________
void StPicoKFVertexTools::makeKFReso(std::vector<int>&  primaryTracks, int nHftTracks) {
    const int nTestedRefits = 2;
    const unsigned int nPrimTracks = primaryTracks.size();
    std::vector<int> setOfTracks[nTestedRefits];

    for (int m = 0; m < nTestedRefits; ++m) {
        setOfTracks[m].resize(1+primaryTracks.size()/2,-999);
    }

    Float_t testDca[nTestedRefits] = {0, 0};
    int testNumber = 0;
    StPicoTrack *test = new StPicoTrack();

    do {
        testNumber++;
        std::random_shuffle(std::begin(primaryTracks), std::end(primaryTracks));
//        for (int k = 0; k < nTestedRefits; ++k) {
//            setOfTracks[k].clear();
//            setOfTracks[k].shrink_to_fit();
//        }
        for (unsigned int i = 0; i < nPrimTracks/2; ++i) {
            setOfTracks[0][i] = primaryTracks[i];
        }

        for (unsigned int j = nPrimTracks/2; j < nPrimTracks; ++j) {
            setOfTracks[1][j-nPrimTracks/2] = primaryTracks[j];
        }

        for (int l = 0; l < nTestedRefits; ++l) {
            testDca[l]=0;
            for (unsigned int i = 0; i < setOfTracks[l].size(); ++i) {
                if (setOfTracks[l][i]<0) continue;
                test = mPicoDst->track(setOfTracks[l][i]);
                testDca[l] += test->gDCAx(mPrimVtx.x());
                testDca[l] += test->gDCAy(mPrimVtx.y());
                testDca[l] += test->gDCAz(mPrimVtx.z());
            }
            testDca[l] = testDca[l] / (3*setOfTracks[l].size());
        }
    } while (abs(testDca[0]-testDca[1]) > 0.025 && testNumber < 20); //-> one of this is false, do..while ends.

    for (int n = 0; n < nTestedRefits; ++n) {
        setOfTracks[n].erase(std::remove(setOfTracks[n].begin(), setOfTracks[n].end(), -999), setOfTracks[n].end());
    }

    StPicoKFVertexFitter kfVertexFitterSet[nTestedRefits];
    KFVertex kfVertexSet[nTestedRefits];
    for (int k = 0; k < nTestedRefits; ++k) {
        kfVertexSet[k] = kfVertexFitterSet[k].primaryVertexRefitUsingTracks(mPicoDst, setOfTracks[k]);
    }
    const int nNtVars = ntp_KFReso->GetNvar();
    Float_t ntVarKF[nNtVars];
    int ii = 0;
    ntVarKF[ii++] = kfVertexSet[0].GetX();
    ntVarKF[ii++] = kfVertexSet[0].GetY();
    ntVarKF[ii++] = kfVertexSet[0].GetZ();

    ntVarKF[ii++] = kfVertexSet[1].GetX();
    ntVarKF[ii++] = kfVertexSet[1].GetY();
    ntVarKF[ii++] = kfVertexSet[1].GetZ();

    float diffX = kfVertexSet[0].GetX()-kfVertexSet[1].GetX();
    float diffY = kfVertexSet[0].GetY()-kfVertexSet[1].GetY();
    float diffZ = kfVertexSet[0].GetZ()-kfVertexSet[1].GetZ();

    ntVarKF[ii++] = diffX;
    ntVarKF[ii++] = diffY;
    ntVarKF[ii++] = diffZ;

    ntVarKF[ii++] = sqrt(diffX*diffX + diffY*diffY + diffZ*diffZ);
    ntVarKF[ii++] = nPrimTracks;
    ntVarKF[ii++] = nHftTracks;
    ntVarKF[ii++] = mPicoEvent->grefMult();

    ntp_KFReso->Fill(ntVarKF);

    setOfTracks[0].clear();
    setOfTracks[0].shrink_to_fit();
    setOfTracks[1].clear();
    setOfTracks[1].shrink_to_fit();
}
