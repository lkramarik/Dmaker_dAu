#include <algorithm>
//#include "StarRoot/MTrack.h"
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

KFVertex StPicoKFVertexFitter::primaryVertexRefit(StPicoDst const* const picoDst, std::vector<int>& tracksToRemove) {
    // just in case it is not sorted
    std::sort(tracksToRemove.begin(),tracksToRemove.end());

    vector<int> goodTracks;
    TVector3 Vtx = picoDst->event()->primaryVertex();

    // make a list of good tracks to be used in the KFVertex fit
    for (unsigned int iTrk = 0; iTrk < picoDst->numberOfTracks(); ++iTrk) {
        StPicoTrack* gTrack = (StPicoTrack*)picoDst->track(iTrk);
        if (!gTrack) continue;
        if (!gTrack->isPrimary()) continue;
        if(std::binary_search(tracksToRemove.begin(), tracksToRemove.end(), iTrk)) continue;
        goodTracks.push_back(iTrk);
    }

    return primaryVertexRefitUsingTracks(picoDst,goodTracks);
}

//________________________________________________________________________________
KFVertex StPicoKFVertexFitter::primaryVertexRefitUsingTracks(StPicoDst const* const picoDst, std::vector<int>& tracksToUse) {
    // fill an array of KFParticles
    KFParticle* particles[tracksToUse.size()];
    TVector3 Vtx = picoDst->event()->primaryVertex();
//    StPicoTrack* gTrack=new StPicoTrack();
//    StPicoTrackCovMatrix *cov=new StPicoTrackCovMatrix();

    for (size_t iTrk = 0; iTrk < tracksToUse.size(); ++iTrk) {
        StPicoTrack* gTrack = (StPicoTrack*)picoDst->track(tracksToUse[iTrk]);
        StPicoTrackCovMatrix *cov = (StPicoTrackCovMatrix*)picoDst->trackCovMatrix(tracksToUse[iTrk]);
        StDcaGeometry dcaG = dcaGeometry(cov);
        Double_t xyzp[6], CovXyzp[21];
        dcaG.GetXYZ(xyzp, CovXyzp);

        /*
	    Int_t q = 1;
        Int_t pdg = 211;
        if (gTrack->charge() < 0) {
            q=-1;
            pdg=-211;
        }
        MTrack track;
        track.SetParameters(xyzp);
        track.SetCovarianceMatrix(CovXyzp);
        track.SetNDF(1);
        track.SetID(gTrack->id());
        track.SetCharge(q);

        particles[iTrk] = new KFParticle(track, pdg);
        */

        Int_t q = 1;
        Int_t pdg = 211;
        if (gTrack->charge() < 0) {
            q=-1;
            pdg=-211;
        }

        particles[iTrk] = new KFParticle();
        particles[iTrk]->Create(xyzp, CovXyzp, q, pdg);
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

    aVertex.ConstructPrimaryVertex((const KFParticle **) particles, tracksToUse.size(), (Bool_t*) Flag.GetArray(), 5); //Yuri - particles, that have worst hi2cut, are not used for final fit
    //if flag is 0, track was not used for PV refit. If all flags are 0, PV is probably not good

    nFlag0=0;
    nFlag1=0;
    for (int i = 0; i < tracksToUse.size(); ++i) {
        if (Flag[i]==0) nFlag0++;
        if (Flag[i]==1) nFlag1++;
    }

//    aVertex.ConstructPrimaryVertex((const KFParticle **) particles, tracksToUse.size(), (Bool_t*) Flag.GetArray(), TMath::Sqrt(StAnneling::Chi2Cut()/2)); //original
//    aVertex.ConstructPrimaryVertex((const KFParticle **) particles, tracksToUse.size(), (Bool_t*) Flag.GetArray(), 3.5); //in makzym stuffs
//    aVertex.ConstructPrimaryVertex((const KFParticle **) particles, tracksToUse.size(), (Bool_t*) Flag.GetArray(), StAnneling::Chi2Cut()); //ok

//Yuri
//    KFVertex     Vertex;
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

    return aVertex;
}

//________________________________________________________________________________
StDcaGeometry StPicoKFVertexFitter::dcaGeometry(StPicoTrackCovMatrix const* cov) const {
    static StDcaGeometry a;
    Float_t errMatrix[15];
    Float_t mSigma[5];
    Float_t mCorr[10];

    const Float_t* sig = cov->sigmas();
    for (int ii = 0; ii < 5; ++ii){
        mSigma[ii]=*sig;
        sig++;
    }

    const Float_t* corr = cov->correlations();
    for (int ii = 0; ii < 10; ++ii){
        mCorr[ii]=*corr;
        corr++;
    }

    Int_t ii = 0;
    for (int i = 0; i < 5; i++) {
        errMatrix[ii] = mSigma[i]*mSigma[i];
        for (int j = 0; j < i; j++) {
            Int_t ij = ii - i - 1 + j + 1;
            Int_t ij1 = ij - i;
            errMatrix[ij] = mCorr[ij1]*mSigma[i]*mSigma[j];
        }
        ii += i+2;
    }

    a.set(cov->params(), errMatrix);
    return *&a;
}