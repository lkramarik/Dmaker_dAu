#ifndef StKaonPion_hh
#define StKaonPion_hh
#ifdef __ROOT__

/* **************************************************
 *  A specialized pair class for calculating K-π pair 
 *  lorentz vector and topological decay parameters 
 *  and storing them.
 *
 *  Authors:  Xin Dong (xdong@lbl.gov),
 *          **Mustafa Mustafa (mmustafa@lbl.gov)
 *
 *  ** Code Maintainer
 *
 * **************************************************
 */
#include <cmath>
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TObject.h"
#include "TClonesArray.h"

class StPicoTrack;
class StPicoEvent;

class StKaonPion : public TObject
{
 public:
  StKaonPion();
  StKaonPion(StKaonPion const *);
  StKaonPion(StPicoTrack const * kaon, StPicoTrack const * pion,unsigned short kIdx,unsigned short pIdx,
             TVector3 const & vtx, float bField);
  ~StKaonPion() {}// please keep this non-virtual and NEVER inherit from this class 

  TLorentzVector const & lorentzVector() const;
  float m()    const;
  float pt()   const;
  float eta()  const;
  float phi()  const;
  float pointingAngle() const;
  float decayLength() const;
  float kaonDca() const;
  float pionDca() const;
  unsigned short   kaonIdx() const;
  unsigned short   pionIdx() const;
  float dcaDaughters() const;
  float cosThetaStar() const;
  float perpDcaToVtx() const;
          
 private:
  // disable copy constructor and assignment operator by making them private (once C++11 is available in STAR you can use delete specifier instead)
  StKaonPion(StKaonPion const &);
  StKaonPion& operator=(StKaonPion const &);
  TLorentzVector mLorentzVector; // this owns four float only

  float mPointingAngle;
  float mDecayLength;
  float mKaonDca;
  float mPionDca;

  unsigned short  mKaonIdx; // index of track in StPicoDstEvent
  unsigned short  mPionIdx;

  float mDcaDaughters;
  float mCosThetaStar; 

  ClassDef(StKaonPion,1)
};
inline TLorentzVector const & StKaonPion::lorentzVector() const { return mLorentzVector;}
inline float StKaonPion::m()    const { return mLorentzVector.M();}
inline float StKaonPion::pt()   const { return mLorentzVector.Perp();}
inline float StKaonPion::eta()  const { return mLorentzVector.PseudoRapidity();}
inline float StKaonPion::phi()  const { return mLorentzVector.Phi();}
inline float StKaonPion::pointingAngle() const { return mPointingAngle;}
inline float StKaonPion::decayLength() const { return mDecayLength;}
inline float StKaonPion::kaonDca() const { return mKaonDca;}
inline float StKaonPion::pionDca() const { return mPionDca;}
inline unsigned short   StKaonPion::kaonIdx() const { return mKaonIdx;}
inline unsigned short   StKaonPion::pionIdx() const { return mPionIdx;}
inline float StKaonPion::dcaDaughters() const { return mDcaDaughters;}
inline float StKaonPion::cosThetaStar() const { return mCosThetaStar;}
inline float StKaonPion::perpDcaToVtx() const { return mDecayLength*std::sin(mPointingAngle);}

#endif
#endif

