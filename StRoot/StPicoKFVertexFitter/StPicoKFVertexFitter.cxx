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
#include "StPicoEvent/StPicoTrackCovMatrix.h"

using namespace std;

TVector3 StPicoKFVertexFitter::primaryVertexRefit(StPicoDst const* const picoDst, std::vector<int>& tracksToRemove) const {
    // just in case it is not sorted
    std::sort(tracksToRemove.begin(),tracksToRemove.end());

    vector<int> goodTracks;

    // make a list of good tracks to be used in the KFVertex fit
    for (unsigned int iTrk = 0; iTrk < picoDst->numberOfTracks(); ++iTrk) {
        StPicoTrack* gTrack = (StPicoTrack*)picoDst->track(iTrk);
        if (! gTrack) continue;
        if(std::binary_search(tracksToRemove.begin(), tracksToRemove.end(), iTrk)) continue;
        goodTracks.push_back(iTrk);
    }

    return primaryVertexRefitUsingTracks(picoDst,goodTracks);
}

TVector3 StPicoKFVertexFitter::primaryVertexRefitUsingTracks(StPicoDst const* const picoDst, std::vector<int>& tracksToUse) const {
    // fill an array of KFParticles
    KFParticle* particles[tracksToUse.size()];

    for (size_t iTrk = 0; iTrk < tracksToUse.size(); ++iTrk) {
        StPicoTrack* gTrack = (StPicoTrack*)picoDst->track(tracksToUse[iTrk]);
//        StDcaGeometry dcaG = gTrack->dcaGeometry();

        StPicoTrackCovMatrix *cov = picoDst->trackCovMatrix(tracksToUse[iTrk]);
//        const StDcaGeometry dcaG = cov->dcaGeometry();
//        Double_t xyzp[6], CovXyzp[21];
//        dcaG.GetXYZ(xyzp, CovXyzp);
        MTrack track;
//        track.SetParameters(xyzp);
        track.SetCovarianceMatrix(cov);
        track.SetNDF(1);
        track.SetID(gTrack->id());
        track.SetCharge(dcaG.charge());

        Int_t pdg = dcaG.charge() > 0 ? 211 : -211; // assume all tracks are pions.

        particles[iTrk] = new KFParticle(track, pdg);
    }

    TArrayC Flag(tracksToUse.size());
    KFVertex aVertex;
    aVertex.ConstructPrimaryVertex((const KFParticle **) particles, tracksToUse.size(), (Bool_t*) Flag.GetArray(), TMath::Sqrt(StAnneling::Chi2Cut() / 2));

    // clean up
    for(size_t iTrk = 0; iTrk < tracksToUse.size(); ++iTrk) delete particles[iTrk];

    TVector3 kfVertex(-999.,-999.,-999.);

    if (aVertex.GetX()) {
        kfVertex.SetXYZ(aVertex.GetX(), aVertex.GetY(), aVertex.GetZ());
    }

    return kfVertex;

}
//michal
//    Int_t q = 1; if (gTrack->charge() < 0) q = -1;
//    KFPTrack track;
//    if( !GetTrack(dcaG, track, q, index) ) continue;
