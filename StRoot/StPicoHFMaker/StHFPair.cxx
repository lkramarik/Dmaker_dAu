#include <limits>
#include <cmath>
#include "StHFPair.h"
#include "StPicoEvent/StPicoPhysicalHelix.h"
#include "StarClassLibrary/SystemOfUnits.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StTypeDefs.h"

ClassImp(StHFPair)

// _________________________________________________________
StHFPair::StHFPair(): mLorentzVector(TLorentzVector()), mDecayVertex(TVector3()),
  mPointingAngle(std::numeric_limits<float>::quiet_NaN()), mDecayLength(std::numeric_limits<float>::quiet_NaN()),
  mParticle1Dca(std::numeric_limits<float>::quiet_NaN()), mParticle2Dca(std::numeric_limits<float>::quiet_NaN()),
  mParticle1Idx(std::numeric_limits<unsigned short>::max()), mParticle2Idx(std::numeric_limits<unsigned short>::max()),
  mDcaDaughters(std::numeric_limits<float>::max()), mCosThetaStar(std::numeric_limits<float>::quiet_NaN()) {
}

// _________________________________________________________
StHFPair::StHFPair(StHFPair const * t) : mLorentzVector(t->mLorentzVector), mDecayVertex(t->mDecayVertex),
   mPointingAngle(t->mPointingAngle), mDecayLength(t->mDecayLength),
   mParticle1Dca(t->mParticle1Dca), mParticle2Dca(t->mParticle2Dca),
   mParticle1Idx(t->mParticle1Idx), mParticle2Idx(t->mParticle2Idx),
   mDcaDaughters(t->mDcaDaughters), mCosThetaStar(t->mCosThetaStar) {
}

// _________________________________________________________
StHFPair::~StHFPair() {
    // destructor

}

// _________________________________________________________
StHFPair::StHFPair(StPicoTrack const * const particle1, StPicoTrack const * const particle2,
		   float p1MassHypo, float p2MassHypo, unsigned short const p1Idx, unsigned short const p2Idx,
		   TVector3 const & vtx, float const bField, bool const useStraightLine) :
  mLorentzVector(TLorentzVector()), mDecayVertex(TVector3()),
  mPointingAngle(std::numeric_limits<float>::quiet_NaN()), mDecayLength(std::numeric_limits<float>::quiet_NaN()),
  mParticle1Dca(std::numeric_limits<float>::quiet_NaN()), mParticle2Dca(std::numeric_limits<float>::quiet_NaN()),
  mParticle1Idx(p1Idx), mParticle2Idx(p2Idx),
  mDcaDaughters(std::numeric_limits<float>::max()), mCosThetaStar(std::numeric_limits<float>::quiet_NaN()) {

  if ((!particle1 || !particle2) || (particle1->id() == particle2->id())) {
    mParticle1Idx = std::numeric_limits<unsigned short>::quiet_NaN();
    mParticle2Idx = std::numeric_limits<unsigned short>::quiet_NaN();
    return;
  }

  StPicoPhysicalHelix p1Helix = particle1->helix(bField);
  StPicoPhysicalHelix p2Helix = particle2->helix(bField);

  // move origins of helices to the primary vertex origin - same as Liang.
  p1Helix.moveOrigin(p1Helix.pathLength(vtx));
  p2Helix.moveOrigin(p2Helix.pathLength(vtx));

  // use straight lines approximation to get point of DCA of kaon-pion pair
  TVector3 const p1Mom = p1Helix.momentum(bField * kilogauss);
  TVector3 const p2Mom = p2Helix.momentum(bField * kilogauss);

  StPicoPhysicalHelix const p1StraightLine(p1Mom, p1Helix.origin(), 0, particle1->charge());
  StPicoPhysicalHelix const p2StraightLine(p2Mom, p2Helix.origin(), 0, particle2->charge());

//  pair<double, double> const ss = (useStraightLine) ? p1StraightLine.pathLengths(p2StraightLine) : p1Helix.pathLengths(p2Helix);
  pair<double, double> const ss = p1StraightLine.pathLengths(p2StraightLine);
  TVector3 const p1AtDcaToP2 = p1StraightLine.at(ss.first);
  TVector3 const p2AtDcaToP1 = p2StraightLine.at(ss.second);

  // -- calculate DCA of particle1 to particle2 at their DCA
  mDcaDaughters = (p1AtDcaToP2 - p2AtDcaToP1).Mag();

  // -- calculate decay vertex (secondary or tertiary)
  mDecayVertex = (p1AtDcaToP2 + p2AtDcaToP1) * 0.5 ;

  // -- calculate Lorentz vector of particle1-particle2 pair, straight line commented 09.01.2018
  TVector3 const p1MomAtDca = p1Helix.momentumAt(ss.first,  bField * kilogauss);
  TVector3 const p2MomAtDca = p2Helix.momentumAt(ss.second, bField * kilogauss);

//  StLorentzVectorF const p1FourMom(p1MomAtDca, p1MomAtDca.massHypothesis(p1MassHypo));
//  StLorentzVectorF const p2FourMom(p2MomAtDca, p2MomAtDca.massHypothesis(p2MassHypo));

  TLorentzVector const p1FourMom(p1MomAtDca, sqrt(p1MomAtDca*p1MomAtDca+p1MassHypo*p1MassHypo));
  TLorentzVector const p2FourMom(p2MomAtDca, sqrt(p2MomAtDca*p2MomAtDca+p2MassHypo*p2MassHypo));

  mLorentzVector = p1FourMom + p2FourMom;

  // -- calculate cosThetaStar
//  TLorentzVector pairFourMomReverse;
//  pairFourMomReverse.SetPxPyPzE(-mLorentzVector.Px(), -mLorentzVector.Py(), -mLorentzVector.Pz(), mLorentzVector.E());
//
//  TLorentzVector p1FourMomStar = p1FourMom;
//  TLorentzVector p2FourMomStar = p2FourMom;
//
//  TVector3 beta = pairFourMomReverse.BoostVector();
//  p1FourMomStar.Boost(beta);
//  p2FourMomStar.Boost(beta);
//  mCosThetaStar = cos(p1FourMomStar.Vect().Angle(mLorentzVector.Vect()));

  mCosThetaStar = p2FourMom.Vect().Unit().Dot(mLorentzVector.Vect().Unit()); //same as in FastSim
  if (mCosThetaStar!=mCosThetaStar) mCosThetaStar=-999;

  TVector3 const vtxToV0 = mDecayVertex - vtx;
  mPointingAngle = vtxToV0.Angle(mLorentzVector.Vect());
  mDecayLength = vtxToV0.Mag();

  mParticle1Dca = (p1Helix.origin() - vtx).Mag();
  mParticle2Dca = (p2Helix.origin() - vtx).Mag();

  if (mParticle1Dca<0 || mParticle2Dca<0) mParticle2Dca=9999;
//  mParticle1Dca = (vtx - particle1->origin()).Mag();
//  mParticle2Dca = (vtx - particle2->origin()).Mag();
}

// _________________________________________________________
float StHFPair::pointingAngle(TVector3 const & vtx2) const{
  // -- Overloaded function recalculates pointing angle given secondary vertex
  TVector3 const vtx2ToTertiary(mDecayVertex - vtx2);
  float const nPointingAngle = vtx2ToTertiary.Angle(mLorentzVector.Vect());
  return nPointingAngle;
}
// _________________________________________________________
float StHFPair::decayLength(TVector3 const & vtx2) const{
  // -- Overloaded function recalculates decayLength given secondary vertex
  TVector3 const vtx2ToTertiary(mDecayVertex - vtx2);
  float const nDecayLength = vtx2ToTertiary.Mag();
  return nDecayLength;
}
// _________________________________________________________
float StHFPair::particle1Dca(StPicoTrack const * p1track, TVector3 const & vtx2, float const bField) const{
  // -- Overloaded function recalculates daughter dca 2 updated vertex
  StPicoPhysicalHelix p1Helix = p1track->helix(bField);
  // -- move origins of helices to the primary vertex origin
  p1Helix.moveOrigin(p1Helix.pathLength(vtx2));
  float const nParticle1Dca = (p1Helix.origin() - vtx2).Mag();
  return nParticle1Dca;
}
// _________________________________________________________
float StHFPair::particle2Dca(StPicoTrack const * p2track, TVector3 const & vtx2, float const bField) const{
  // -- Overloaded function recalculates daughter dca 2 updated vertex
  StPicoPhysicalHelix p2Helix = p2track->helix(bField);
  // -- move origins of helices to the primary vertex origin
  p2Helix.moveOrigin(p2Helix.pathLength(vtx2));
  float const nParticle2Dca = (p2Helix.origin() - vtx2).Mag();
  return nParticle2Dca;
}


