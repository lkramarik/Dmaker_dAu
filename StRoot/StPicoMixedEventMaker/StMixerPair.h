#ifndef StMixerPair_hh
#define StMixerPair_hh

/* **************************************************
 *  Generic class calculating and storing pairs in Event Mixing
 *  Allows to combine:
 *  - two particles, using
 *      StMixerPair(StPicoTrack const * particle1, StPicoTrack const * particle2, ...
 *
 * **************************************************
 *
 *  Initial Authors: 
 *         ** Michael Lomnitz (mrlomnitz@lbl.gov)
 *            Mustafa Mustafa (mmustafa@lbl.gov)
 *
 *  ** Code Maintainer 
 *
 * **************************************************
 */

#include "TObject.h"
#include "TClonesArray.h"
#include "StPicoEvent/StPicoEvent.h"
#include "TVector3.h"
#include "TLorentzVector.h"

class StPicoTrack;
//class StMixerTrack;

class StMixerPair : public TObject
{
public:
    StMixerPair();
    StMixerPair(StMixerPair const*);

//  StMixerPair(StMixerTrack const&  particle1, StMixerTrack const& particle2,
//	   float p1MassHypo, float p2MassHypo,
//	   TVector3 const& vtx1, TVector3 const& vtx2,
//	   float bField);

    StMixerPair(StPicoTrack const&  particle1, StPicoTrack const& particle2,
                float p1MassHypo, float p2MassHypo,
                TVector3 const& vtx1, TVector3 const& vtx2,
                float bField);

    ~StMixerPair() {;}


    TLorentzVector const & lorentzVector() const;
    TVector3 const & decayVertex() const;
    float rapidity()    const;
    float m()    const;
    float pt()   const;
    float eta()  const;
    float phi()  const;
    float pointingAngle() const;
    float decayLength() const;
    float particle1Dca() const;
    float particle2Dca() const;
    TVector3 const & particle1Mom() const;
    TVector3 const & particle2Mom() const;
    float dcaDaughters() const;
    float cosThetaStar() const;
    float v0x() const;
    float v0y() const;
    float v0z() const;
    float px() const;
    float py() const;
    float pz() const;
    float DcaToPrimaryVertex() const;


private:
    StMixerPair(StMixerPair const &);
    StMixerPair& operator=(StMixerPair const &);
    TLorentzVector mLorentzVector;
    TVector3   mDecayVertex;

    float mPointingAngle;
    float mDecayLength;
    float mParticle1Dca;
    float mParticle2Dca;

    TVector3 mParticle1Mom;
    TVector3 mParticle2Mom;

    float mDcaDaughters;
    float mCosThetaStar;

    ClassDef(StMixerPair,1)
};
inline TLorentzVector const& StMixerPair::lorentzVector() const { return mLorentzVector;}
inline float StMixerPair::rapidity()    const { return mLorentzVector.Rapidity();}
inline float StMixerPair::m()    const { return mLorentzVector.M();}
inline float StMixerPair::pt()   const { return mLorentzVector.Perp();}
inline float StMixerPair::eta()  const { return mLorentzVector.PseudoRapidity();}
inline float StMixerPair::phi()  const { return mLorentzVector.Phi();}
inline float StMixerPair::px()   const { return mLorentzVector.Px();}
inline float StMixerPair::py()   const { return mLorentzVector.Py();}
inline float StMixerPair::pz()   const { return mLorentzVector.Pz();}
inline float StMixerPair::pointingAngle() const { return mPointingAngle;}
inline float StMixerPair::decayLength()   const { return mDecayLength;}
inline float StMixerPair::particle1Dca()  const { return mParticle1Dca;}
inline float StMixerPair::particle2Dca()  const { return mParticle2Dca;}
inline TVector3 const & StMixerPair::particle1Mom() const { return mParticle1Mom; }
inline TVector3 const & StMixerPair::particle2Mom() const { return mParticle2Mom; }
inline float StMixerPair::dcaDaughters() const { return mDcaDaughters;}
inline float StMixerPair::cosThetaStar() const { return mCosThetaStar;}
inline TVector3 const & StMixerPair::decayVertex() const { return mDecayVertex;}
inline float StMixerPair::v0x() const { return mDecayVertex.x();}
inline float StMixerPair::v0y() const { return mDecayVertex.y();}
inline float StMixerPair::v0z() const { return mDecayVertex.z();}
inline float StMixerPair::DcaToPrimaryVertex() const { return mDecayLength*sin(mPointingAngle); }

#endif

