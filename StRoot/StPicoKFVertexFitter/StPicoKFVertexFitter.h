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

class StPicoDst;

class StPicoKFVertexFitter
{
public:
    StPicoKFVertexFitter() {}
    ~StPicoKFVertexFitter() {}

    TVector3 primaryVertexRefit(StPicoDst const*) const;

    TVector3 primaryVertexRefit(StPicoDst const*,
                                      std::vector<int>& tracksToRemove) const;

    TVector3 primaryVertexRefitUsingTracks(StPicoDst const*,
                                                 std::vector<int>& tracksToUse) const;
};

inline TVector3 StPicoKFVertexFitter::primaryVertexRefit(StPicoDst const* picoDst) const
{
    std::vector<int> v;
    return primaryVertexRefit(picoDst,v);
}
#endif
