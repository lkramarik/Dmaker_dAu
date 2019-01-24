#include <algorithm>
#include "StarRoot/MTrack.h"
#include "StiMaker/StKFVerticesCollection.h"
#include "StEvent/StDcaGeometry.h"
#include "StPicoKFVertexFitter.h"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StEvent/StGlobalTrack.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoTrackCovMatrix/StPicoTrackCovMatrix.h"
//#include "KFPTrackVector.h"
//#include "KFParticle.h"
//#include "KFPTrack.h"

using namespace std;

TVector3 StPicoKFVertexFitter::primaryVertexRefit(StPicoDst const* const picoDst, std::vector<int>& tracksToRemove) const {
    // just in case it is not sorted
    std::sort(tracksToRemove.begin(),tracksToRemove.end());

    vector<int> goodTracks;
//    StPicoEvent mPicoEvent = picoDst->event();
    TVector3 Vtx = picoDst->event()->primaryVertex();

    // make a list of good tracks to be used in the KFVertex fit
    for (unsigned int iTrk = 0; iTrk < picoDst->numberOfTracks(); ++iTrk) {
        StPicoTrack* gTrack = (StPicoTrack*)picoDst->track(iTrk);
        if (! gTrack) continue;
        if(std::binary_search(tracksToRemove.begin(), tracksToRemove.end(), iTrk)) continue;
        if(!gTrack->isPrimary()) continue;
//        if(abs(gTrack->gDCAz(Vtx.z()))>3) continue;
//        if(abs(gTrack->gDCAxy(Vtx.x(),Vtx.y()))>1.5) continue;
        goodTracks.push_back(iTrk);
    }

    return primaryVertexRefitUsingTracks(picoDst,goodTracks);
}

TVector3 StPicoKFVertexFitter::primaryVertexRefitUsingTracks(StPicoDst const* const picoDst, std::vector<int>& tracksToUse) const {
    // fill an array of KFParticles
    KFParticle* particles[tracksToUse.size()];

    for (size_t iTrk = 0; iTrk < tracksToUse.size(); ++iTrk) {
        StPicoTrack* gTrack = (StPicoTrack*)picoDst->track(tracksToUse[iTrk]);
        StPicoTrackCovMatrix *cov = picoDst->trackCovMatrix(tracksToUse[iTrk]);
        const StDcaGeometry dcaG = cov->dcaGeometry();

        Double_t xyzp[6], CovXyzp[21]; //ok
        dcaG.GetXYZ(xyzp, CovXyzp); //ok

        Int_t q = 1; if (gTrack->charge() < 0) q = -1; //ok

        MTrack track;
//        KFPTrack track1;
        track.SetParameters(xyzp); //ok
        track.SetCovarianceMatrix(CovXyzp); //ok
        track.SetNDF(1); //ok
        track.SetID(gTrack->id()); //ok
        track.SetCharge(q); //ok

//        Int_t pdg = dcaG.charge() > 0 ? 211 : -211; // assume all tracks are pions.
//
        particles[iTrk] = new KFParticle(track, 211);
//        particles[iTrk] = new KFParticle(track, pdg);
    }

    TArrayC Flag(tracksToUse.size());
    KFVertex aVertex;
//    aVertex.ConstructPrimaryVertex((const KFParticle **) particles, tracksToUse.size(), (Bool_t*) Flag.GetArray(), TMath::Sqrt(StAnneling::Chi2Cut() / 2)); //original
//    aVertex.ConstructPrimaryVertex((const KFParticle **) particles, tracksToUse.size(), (Bool_t*) Flag.GetArray(), 2.); //ok
    aVertex.ConstructPrimaryVertex((const KFParticle **) particles, tracksToUse.size(), (Bool_t*) Flag.GetArray(), StAnneling::Chi2Cut()); //ok

    // clean up
    for(size_t iTrk = 0; iTrk < tracksToUse.size(); ++iTrk) delete particles[iTrk];

    TVector3 kfVertex(-999.,-999.,-999.);

    if (aVertex.GetX()) {
        kfVertex.SetXYZ(aVertex.GetX(), aVertex.GetY(), aVertex.GetZ());
    }

    return kfVertex;

}

//StDcaGeometry &StPicoKFVertexFitter::dcaGeometry(StPicoTrackCovMatrix *cov) const {
//    static StDcaGeometry a;
//    Float_t errMatrix[15];
//    Float16_t* mSigma=cov->sigmas();
//    Int_t ii = 0;
//    for (int i = 0; i < 5; i++) {
//        errMatrix[ii] = mSigma[i]*mSigma[i];
//        for (int j = 0; j < i; j++) {
//            Int_t ij = ii - i - 1 + j + 1;
//            Int_t ij1 = ij - i;
//            errMatrix[ij] = mCorr[ij1]*mSigma[i]*mSigma[j];
//        }
//        ii += i+2;
//    }
//    a.set(params(), errMatrix);
//    return *&a;
//}

//michal
//    Int_t q = 1; if (gTrack->charge() < 0) q = -1;
//    KFPTrack track;
//    if( !GetTrack(dcaG, track, q, index) ) continue;
