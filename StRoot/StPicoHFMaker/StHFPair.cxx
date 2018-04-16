#include <limits>
#include <cmath>
#include "StHFPair.h"
#include "StarClassLibrary/StPhysicalHelixD.hh"
#include "StarClassLibrary/SystemOfUnits.h"
#include "StPicoEvent/StPicoTrack.h"

ClassImp(StHFPair)

// _________________________________________________________
StHFPair::StHFPair(): mLorentzVector(StLorentzVectorF()), mDecayVertex(StThreeVectorF()),
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
StHFPair::StHFPair(StPicoTrack const * const particle1, StPicoTrack const * const particle2,
		   float p1MassHypo, float p2MassHypo, unsigned short const p1Idx, unsigned short const p2Idx,
		   StThreeVectorF const & vtx, float const bField, bool const useStraightLine) : 
  mLorentzVector(StLorentzVectorF()), mDecayVertex(StThreeVectorF()),
  mPointingAngle(std::numeric_limits<float>::quiet_NaN()), mDecayLength(std::numeric_limits<float>::quiet_NaN()),
  mParticle1Dca(std::numeric_limits<float>::quiet_NaN()), mParticle2Dca(std::numeric_limits<float>::quiet_NaN()),
  mParticle1Idx(p1Idx), mParticle2Idx(p2Idx),
  mDcaDaughters(std::numeric_limits<float>::max()), mCosThetaStar(std::numeric_limits<float>::quiet_NaN()) {

  if ((!particle1 || !particle2) || (particle1->id() == particle2->id())) {
    mParticle1Idx = std::numeric_limits<unsigned short>::quiet_NaN();
    mParticle2Idx = std::numeric_limits<unsigned short>::quiet_NaN();
    return;
  }

  StPhysicalHelixD p1Helix = particle1->helix(bField);
  StPhysicalHelixD p2Helix = particle2->helix(bField);

  // move origins of helices to the primary vertex origin - same as Liang.
  p1Helix.moveOrigin(p1Helix.pathLength(vtx));
  p2Helix.moveOrigin(p2Helix.pathLength(vtx));

  // use straight lines approximation to get point of DCA of kaon-pion pair
  StThreeVectorF const p1Mom = p1Helix.momentum(bField * kilogauss);
  StThreeVectorF const p2Mom = p2Helix.momentum(bField * kilogauss);

  StPhysicalHelixD const p1StraightLine(p1Mom, p1Helix.origin(), 0, particle1->charge());
  StPhysicalHelixD const p2StraightLine(p2Mom, p2Helix.origin(), 0, particle2->charge());

//  pair<double, double> const ss = (useStraightLine) ? p1StraightLine.pathLengths(p2StraightLine) : p1Helix.pathLengths(p2Helix);
  pair<double, double> const ss = p1StraightLine.pathLengths(p2StraightLine);
  StThreeVectorF const p1AtDcaToP2 = p1StraightLine.at(ss.first);
  StThreeVectorF const p2AtDcaToP1 = p2StraightLine.at(ss.second);

  // -- calculate DCA of particle1 to particle2 at their DCA
  mDcaDaughters = (p1AtDcaToP2 - p2AtDcaToP1).mag();

  // -- calculate decay vertex (secondary or tertiary)
  mDecayVertex = (p1AtDcaToP2 + p2AtDcaToP1) * 0.5 ;

  // -- calculate Lorentz vector of particle1-particle2 pair, straight line commented 09.01.2018
  StThreeVectorF const p1MomAtDca = p1Helix.momentumAt(ss.first,  bField * kilogauss);
  StThreeVectorF const p2MomAtDca = p2Helix.momentumAt(ss.second, bField * kilogauss);

  StLorentzVectorF const p1FourMom(p1MomAtDca, p1MomAtDca.massHypothesis(p1MassHypo));
  StLorentzVectorF const p2FourMom(p2MomAtDca, p2MomAtDca.massHypothesis(p2MassHypo));

  mLorentzVector = p1FourMom + p2FourMom;

  // -- calculate cosThetaStar
  StLorentzVectorF const pairFourMomReverse(-mLorentzVector.px(), -mLorentzVector.py(), -mLorentzVector.pz(), mLorentzVector.e());
  StLorentzVectorF const p1FourMomStar = p1FourMom.boost(pairFourMomReverse);
  mCosThetaStar = std::cos(p1FourMomStar.vect().angle(mLorentzVector.vect()));

  StThreeVectorF const vtxToV0 = mDecayVertex - vtx;
  mPointingAngle = vtxToV0.angle(mLorentzVector.vect());
  mDecayLength = vtxToV0.mag();

  mParticle1Dca = (p1Helix.origin() - vtx).mag();
  mParticle2Dca = (p2Helix.origin() - vtx).mag();

//  mParticle1Dca = (vtx - particle1->dcaPoint()).mag();
//  mParticle2Dca = (vtx - particle2->dcaPoint()).mag();
}

// _________________________________________________________
float StHFPair::pointingAngle(StThreeVectorF const & vtx2) const{
  // -- Overloaded function recalculates pointing angle given secondary vertex
  StThreeVectorF const vtx2ToTertiary(mDecayVertex - vtx2);
  float const nPointingAngle = vtx2ToTertiary.angle(mLorentzVector.vect());
  return nPointingAngle;
}
// _________________________________________________________
float StHFPair::decayLength(StThreeVectorF const & vtx2) const{
  // -- Overloaded function recalculates decayLength given secondary vertex
  StThreeVectorF const vtx2ToTertiary(mDecayVertex - vtx2); 
  float const nDecayLength = vtx2ToTertiary.mag();  
  return nDecayLength;
}
// _________________________________________________________
float StHFPair::particle1Dca(StPicoTrack const * p1track, StThreeVectorF const & vtx2, float const bField) const{
  // -- Overloaded function recalculates daughter dca 2 updated vertex
  StPhysicalHelixD p1Helix = p1track->helix(bField);
  // -- move origins of helices to the primary vertex origin
  p1Helix.moveOrigin(p1Helix.pathLength(vtx2));
  float const nParticle1Dca = (p1Helix.origin() - vtx2).mag();
  return nParticle1Dca;
}
// _________________________________________________________
float StHFPair::particle2Dca(StPicoTrack const * p2track, StThreeVectorF const & vtx2, float const bField) const{
  // -- Overloaded function recalculates daughter dca 2 updated vertex
  StPhysicalHelixD p2Helix = p2track->helix(bField);
  // -- move origins of helices to the primary vertex origin
  p2Helix.moveOrigin(p2Helix.pathLength(vtx2));
  float const nParticle2Dca = (p2Helix.origin() - vtx2).mag();
  return nParticle2Dca;
}


