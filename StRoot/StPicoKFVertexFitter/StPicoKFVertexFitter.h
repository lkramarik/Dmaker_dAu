#ifndef StPicoKFVertexFitter_h
#define StPicoKFVertexFitter_h

/* **************************************************
 *  Class to fit primary vertex using KF vertex maker
 *
 *  Usage:
 *
 * **************************************************
 *  Authors:
 *            **Liang He(he202@purdue.edu)
 *            Mustafa Mustafa (mmustafa@lbl.gov)
 *
 * ** Code Maintainer
 *             Lukas Kramarik
 * **************************************************
 */

#include <vector>
#include "TVector3.h"
#include "KFVertex.h"
#include "StEvent/StDcaGeometry.h"
#include "StPicoEvent//StPicoTrackCovMatrix.h"

//#include "StEvent/StDcaGeometry.h"

class StPicoDst;

class StPicoKFVertexFitter
{
public:
    StPicoKFVertexFitter() {}
    ~StPicoKFVertexFitter() {}

//    StDcaGeometry &dcaGeometry() const;

    KFVertex primaryVertexRefit(StPicoDst const*) const;

    KFVertex primaryVertexRefit(StPicoDst const*, std::vector<int>&) const;

    KFVertex primaryVertexRefitUsingTracks(StPicoDst const*, std::vector<int>&) const;

private:
    StDcaGeometry dcaGeometry(StPicoTrackCovMatrix const*);
};

inline KFVertex StPicoKFVertexFitter::primaryVertexRefit(StPicoDst const*) const
{
    std::vector<int> v;
    return primaryVertexRefit(picoDst,v);
}
#endif
