#include <limits>
#include <cmath>

#include "StPicoEvent/StPicoPhysicalHelix.h"
#include "StarClassLibrary/SystemOfUnits.h"

#include "StMixerPair.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StMixerTrack.h"

ClassImp(StMixerPair)

// _________________________________________________________
StMixerPair::StMixerPair(): mLorentzVector(TLorentzVector()), mDecayVertex(TVector3()),
    mPointingAngle(std::numeric_limits<float>::quiet_NaN()), mDecayLength(std::numeric_limits<float>::quiet_NaN()),
    mParticle1Dca(std::numeric_limits<float>::quiet_NaN()), mParticle2Dca(std::numeric_limits<float>::quiet_NaN()),
    mParticle1Mom(TVector3()), mParticle2Mom(TVector3()),
    mDcaDaughters(std::numeric_limits<float>::max()), mCosThetaStar(std::numeric_limits<float>::quiet_NaN()) {
}

// _________________________________________________________
StMixerPair::StMixerPair(StMixerPair const * t) : mLorentzVector(t->mLorentzVector), mDecayVertex(t->mDecayVertex),
    mPointingAngle(t->mPointingAngle), mDecayLength(t->mDecayLength),
    mParticle1Dca(t->mParticle1Dca), mParticle2Dca(t->mParticle2Dca),
    mParticle1Mom(t->mParticle1Mom), mParticle2Mom(t->mParticle2Mom),
    mDcaDaughters(t->mDcaDaughters), mCosThetaStar(t->mCosThetaStar) {
}

// _________________________________________________________
//StMixerPair::StMixerPair(StMixerTrack const& particle1, StMixerTrack const& particle2,
StMixerPair::StMixerPair(StPicoTrack const&  particle1, StPicoTrack const& particle2,
                         float p1MassHypo, float p2MassHypo,
                         TVector3 const& vtx1, TVector3 const& vtx2, float const bField) :  mLorentzVector(TLorentzVector()), mDecayVertex(TVector3()),
    mPointingAngle(std::numeric_limits<float>::quiet_NaN()), mDecayLength(std::numeric_limits<float>::quiet_NaN()),
    mParticle1Dca(std::numeric_limits<float>::quiet_NaN()), mParticle2Dca(std::numeric_limits<float>::quiet_NaN()),
    mParticle1Mom(particle1.gMom()), mParticle2Mom(particle2.gMom()),
    mDcaDaughters(std::numeric_limits<float>::max()), mCosThetaStar(std::numeric_limits<float>::quiet_NaN()) {
    // -- Create pair out of 2 tracks
    //     prefixes code:
    //      p1 means particle 1
    //      p2 means particle 2
    //      pair means particle1-particle2  pair

    TVector3 dVtx = vtx1 -vtx2;

    StPicoPhysicalHelix p1Helix(particle1.gMom(), particle1.origin(),bField*kilogauss, particle1.charge());
    StPicoPhysicalHelix p2Helix(particle2.gMom(), particle2.origin() + dVtx, bField*kilogauss,  particle2.charge());
//    StPhysicalHelixD p1Helix = particle1->helix(bField);
//    StPhysicalHelixD p2Helix = particle2->helix(bField);

    // -- move origins of helices to the primary vertex origin
    p1Helix.moveOrigin(p1Helix.pathLength(vtx1));
    p2Helix.moveOrigin(p2Helix.pathLength(vtx1));

    // -- use straight lines approximation to get point of DCA of particle1-particle2 pair
    TVector3 const p1Mom = p1Helix.momentum(bField * kilogauss);
    TVector3 const p2Mom = p2Helix.momentum(bField * kilogauss);

    StPicoPhysicalHelix const p1StraightLine(p1Mom, p1Helix.origin(), 0, particle1.charge());
    StPicoPhysicalHelix const p2StraightLine(p2Mom, p2Helix.origin(), 0, particle2.charge());

    pair<double, double> const ss = p1StraightLine.pathLengths(p2StraightLine);
    TVector3 const p1AtDcaToP2 = p1StraightLine.at(ss.first);
    TVector3 const p2AtDcaToP1 = p2StraightLine.at(ss.second);

    // -- calculate DCA of particle1 to particle2 at their DCA
    mDcaDaughters = (p1AtDcaToP2 - p2AtDcaToP1).Mag();

    // -- calculate Lorentz vector of particle1-particle2 pair
    TVector3 const p1MomAtDca = p1Helix.momentumAt(ss.first,  bField * kilogauss);
    TVector3 const p2MomAtDca = p2Helix.momentumAt(ss.second, bField * kilogauss);

    TLorentzVector const p1FourMom(p1MomAtDca, sqrt(p1MomAtDca*p1MomAtDca+p1MassHypo*p1MassHypo));
    TLorentzVector const p2FourMom(p2MomAtDca, sqrt(p2MomAtDca*p2MomAtDca+p2MassHypo*p2MassHypo));

    mLorentzVector = p1FourMom + p2FourMom;

    // -- calculate cosThetaStar
//    TVector3 const pairFourMomReverse(-mLorentzVector.Px(), -mLorentzVector.Py(), -mLorentzVector.Pz());
//    TLorentzVector p1FourMomStar = p1FourMom;
//    p1FourMomStar.Boost(pairFourMomReverse);
//    mCosThetaStar = std::cos(p1FourMomStar.Vect().Angle(mLorentzVector.Vect()));

    mCosThetaStar = p2FourMom.Vect().Unit().Dot(mLorentzVector.Vect().Unit()); //same as in FastSim
    if (mCosThetaStar!=mCosThetaStar) mCosThetaStar=-999;

    // -- calculate decay vertex (secondary or tertiary)
    mDecayVertex = (p1AtDcaToP2 + p2AtDcaToP1) * 0.5 ;

    // -- calculate pointing Angle and decay length with respect to primary vertex
    TVector3 const vtxToV0 = mDecayVertex - vtx1;
    mPointingAngle = vtxToV0.Angle(mLorentzVector.Vect());
    mDecayLength = vtxToV0.Mag();

    // -- calculate DCA of tracks to primary vertex
    mParticle1Dca = (p1Helix.origin() - vtx1).Mag();
    mParticle2Dca = (p2Helix.origin() - vtx1).Mag();
}