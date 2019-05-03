#include <algorithm>
#include "StarRoot/MTrack.h"
#include "StarRoot/KFVertex.h"
#include "StiMaker/StKFVerticesCollection.h"
#include "StEvent/StDcaGeometry.h"
#include "StPicoKFVertexFitter.h"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StEvent/StGlobalTrack.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoTrackCovMatrix.h"
//#include "KFPTrackVector.h"
//#include "KFParticle.h"
//#include "KFPTrack.h"

using namespace std;

KFVertex StPicoKFVertexFitter::primaryVertexRefit(StPicoDst const* const picoDst, std::vector<int>& tracksToRemove) const {
    // just in case it is not sorted
    std::sort(tracksToRemove.begin(),tracksToRemove.end());

    vector<int> goodTracks;
    TVector3 Vtx = picoDst->event()->primaryVertex();

    // make a list of good tracks to be used in the KFVertex fit
    for (unsigned int iTrk = 0; iTrk < picoDst->numberOfTracks(); ++iTrk) {
        StPicoTrack* gTrack = (StPicoTrack*)picoDst->track(iTrk);
        if (! gTrack) continue;
        if(std::binary_search(tracksToRemove.begin(), tracksToRemove.end(), iTrk)) continue;
        if(!gTrack->isPrimary()) continue; //no
        goodTracks.push_back(iTrk);
    }

    return primaryVertexRefitUsingTracks(picoDst,goodTracks);
}

KFVertex StPicoKFVertexFitter::primaryVertexRefitUsingTracks(StPicoDst const* const picoDst, std::vector<int>& tracksToUse) const {
    // fill an array of KFParticles
    KFParticle* particles[tracksToUse.size()];
    TVector3 Vtx = picoDst->event()->primaryVertex();

    for (size_t iTrk = 0; iTrk < tracksToUse.size(); ++iTrk) {
        StPicoTrack* gTrack = (StPicoTrack*)picoDst->track(tracksToUse[iTrk]);
        StPicoTrackCovMatrix *cov = picoDst->trackCovMatrix(tracksToUse[iTrk]);

//        const StDcaGeometry dcaG = cov->dcaGeometry();
        Double_t xyzp[6], CovXyzp[21]; //ok
//        dcaG.GetXYZ(xyzp, CovXyzp); //ok
	
	Float_t* p = cov->params();
	for (int ii = 0; ii < 6; ++ii){
		xyzp[ii]=*p;
		cout<<p<<endl;
		p++;
	}        
//			Float_t* pCov = cov->
//        for (int kk = 0; kk < 21; ++kk){
//                CovXyzp[ii]=*p;
//               cout<<p<<endl;
//                p++;
//        }



	Int_t q = 1;
        Int_t pdg = 211;
        if (gTrack->charge() < 0) {
            q=-1;
            pdg=-211;
        }
        MTrack track;
        track.SetParameters(xyzp); //ok
        track.SetCovarianceMatrix(CovXyzp); //ok
        track.SetNDF(1); //ok
        track.SetID(gTrack->id()); //ok
        track.SetCharge(q); //ok

        particles[iTrk] = new KFParticle(track, pdg);
    }

    TArrayC Flag(tracksToUse.size());
    KFVertex aVertex;
    const Double_t parVertex[6] = {Vtx.x(),Vtx.y(),Vtx.z(), 0, 0, 1000};
    const Double_t covVertex[21] = {1,
                              0, 1,
                              0, 0, 100,
                              0, 0,   0, 100,
                              0, 0,   0,   0, 100,
                              0, 0,   0,   0,   0, 100};

    aVertex.Create(parVertex,covVertex, 0, 0);

    aVertex.ConstructPrimaryVertex((const KFParticle **) particles, tracksToUse.size(), (Bool_t*) Flag.GetArray(), 5); //Yuri
//    aVertex.ConstructPrimaryVertex((const KFParticle **) particles, tracksToUse.size(), (Bool_t*) Flag.GetArray(), TMath::Sqrt(StAnneling::Chi2Cut()/2)); //original
//    aVertex.ConstructPrimaryVertex((const KFParticle **) particles, tracksToUse.size(), (Bool_t*) Flag.GetArray(), 3.5); //in makzym stuffs
//    aVertex.ConstructPrimaryVertex((const KFParticle **) particles, tracksToUse.size(), (Bool_t*) Flag.GetArray(), StAnneling::Chi2Cut()); //ok


//Yuri
//    KFVertex     Vertex;
//    const Double_t par[6] = {VGlob.x(),VGlob.y(),VGlob.z(), 0, 0, 1000};
//    const Double_t cov[21] = {1,
//                              0, 1,
//                              0, 0, 100,
//                              0, 0,   0, 100,
//                              0, 0,   0,   0, 100,
//                              0, 0,   0,   0,   0, 100};
//
//    static Double_t Chi2Cut = 10; // Cut on how well particle match to vertex
//
//    Vertex.Create(par,cov, 0, 0);
//    Vertex.SetId(ivx+1);
//
//    UInt_t N = particles.size();
//    if (N < 3) continue;
//    TArrayC Flag(N);
//    KFParticle **parts = new KFParticle*[N];
//
//    for (UInt_t i = 0; i < N; i++) {
//        parts[i] = &particles[i];
//    }
//
//    Vertex.ConstructPrimaryVertex((const KFParticle **) parts, N, (Bool_t*) Flag.GetArray(), Chi2Cut/2);
//    PrPP(Vertex);
//    delete [] parts;
//Yuri end


   // clean up
    for(size_t iTrk = 0; iTrk < tracksToUse.size(); ++iTrk) delete particles[iTrk];


//    TVector3 kfVertex(-999.,-999.,-999.);
//
//    if (aVertex.GetX()) {
//        kfVertex.SetXYZ(aVertex.GetX(), aVertex.GetY(), aVertex.GetZ());
//    }

    return aVertex;
}



