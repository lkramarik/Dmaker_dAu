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
    // constructor
}

// _________________________________________________________
StPicoKFVertexTools::~StPicoKFVertexTools() {
    // destructor
}

// _________________________________________________________
int StPicoKFVertexTools::InitHF() {
    mOutFileBaseName = mOutFileBaseName.ReplaceAll(".root", "");
    ntp_vertex = new TNtuple("ntp_vertex","ntp_vertex","runId:picoDst_vx:picoDst_vy:picoDst_vz:picoDst_vxErr:picoDst_vyErr:picoDst_vzErr:KF_vx:KF_vy:KF_vz:KF_vxErr:KF_vyErr:KF_vzErr");
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
    Float_t runId, picoDst_vx, picoDst_vy, picoDst_vz, picoDst_vxErr, picoDst_vyErr, picoDst_vzErr, KF_vx, KF_vy, KF_vz, KF_vxErr, KF_vyErr, KF_vzErr;

    StPicoKFVertexFitter kfVertexFitter;
    KFVertex kfVertex = kfVertexFitter.primaryVertexRefit(mPicoDst);
    cout<<kfVertex.GetX()<<" "<<kfVertex.GetY()<<" "<<kfVertex.GetZ()<<endl;
    cout<<mPrimVtx.x()<<" "<<mPrimVtx.y()<<" "<<mPrimVtx.z()<<endl;
    cout<<" "<<endl;

    const int nNtVars=ntp_vertex->GetNvar();
    cout<<ntp_vertex<<endl;
    //    Float_t ntVar[];
//
//    runId=mPicoEvent->runId();
//    picoDst_vx=mPrimVtx.x();
//    picoDst_vy=mPrimVtx.y();

    return kStOK;
}

