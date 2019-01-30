#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoHFMaker/StHFCuts.h"
#include "phys_constants.h"
#include "StPicoKFVertexTools.h"
#include "StPicoKFVertexFitter.h"
ClassImp(StPicoKFVertexTools)

// _________________________________________________________
StPicoKFVertexTools::StPicoKFVertexTools(char const* name, StPicoDstMaker* picoMaker, char const* outputBaseFileName) :
        StPicoHFMaker(name, picoMaker, outputBaseFileName),
        mOutFileBaseName(outputBaseFileName){
}

// _________________________________________________________
StPicoKFVertexTools::~StPicoKFVertexTools() {
}

// _________________________________________________________
int StPicoKFVertexTools::InitHF() {
    mOutFileBaseName = mOutFileBaseName.ReplaceAll(".root", "");
    ntp_vertex = new TNtuple("ntp_vertex","ntp_vertex","runId:"
                                                       "picoDst_vx:picoDst_vy:picoDst_vz:"
                                                       "picoDst_vxErr:picoDst_vyErr:picoDst_vzErr:"
                                                       "KF_vx:KF_vy:KF_vz:"
                                                       "KF_vxErr:KF_vyErr:KF_vzErr");
    return kStOK;
}

// _________________________________________________________
void StPicoKFVertexTools::ClearHF(Option_t *opt="") {
    return;
}

// _________________________________________________________
int StPicoKFVertexTools::FinishHF() {
    ntp_vertex -> Write(ntp_vertex->GetName(), TObject::kOverwrite);
    return kStOK;
}
// _________________________________________________________
int StPicoKFVertexTools::MakeHF() {
    StPicoKFVertexFitter kfVertexFitter;
    KFVertex kfVertex = kfVertexFitter.primaryVertexRefit(mPicoDst);

    const int nNtVars=ntp_vertex->GetNvar();
    Float_t ntVar[nNtVars];
    int ii=0;

    ntVar[ii++]=mPicoEvent->runId();

    ntVar[ii++]=mPrimVtx.x();
    ntVar[ii++]=mPrimVtx.y();
    ntVar[ii++]=mPrimVtx.z();

    ntVar[ii++]=(mPicoEvent->primaryVertexError()).x();
    ntVar[ii++]=(mPicoEvent->primaryVertexError()).y();
    ntVar[ii++]=(mPicoEvent->primaryVertexError()).z();

    ntVar[ii++]=kfVertex.GetX();
    ntVar[ii++]=kfVertex.GetY();
    ntVar[ii++]=kfVertex.GetZ();

    ntVar[ii++]=kfVertex.GetErrX();
    ntVar[ii++]=kfVertex.GetErrY();
    ntVar[ii++]=kfVertex.GetErrZ();

    ntp_vertex->Fill(ntVar);

    return kStOK;
}

