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

    KFVertex primaryVertexRefit(StPicoDst const*);
    KFVertex primaryVertexRefit(StPicoDst const*, std::vector<int>&);
    KFVertex primaryVertexRefitUsingTracks(StPicoDst const*, std::vector<int>&);

    int nFlag0;
    int nFlag1;


private:
    StDcaGeometry dcaGeometry(StPicoTrackCovMatrix const*) const;


};

inline KFVertex StPicoKFVertexFitter::primaryVertexRefit(StPicoDst const* picoDst) {
    std::vector<int> v;
    return primaryVertexRefit(picoDst,v);
}
#endif
