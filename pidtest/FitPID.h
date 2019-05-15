//
// Created by lukas on 8.11.2018.
//

#ifndef DMAKER_DAU_STFIT_H
#define DMAKER_DAU_STFIT_H

#include "TStyle.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TMath.h"
#include "TNtuple.h"
#include "TList.h"
#include "TLatex.h"
#include "TF1.h"
#include "TCut.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include <iostream>
#include <ctime>
#include <cstdio>
#include <fstream>
#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include "TObject.h"

using namespace std;
using namespace TMath;

class FitPID : public TObject {
public:
    FitPID();
    ~FitPID();

    TString mOutputFileName;
    TH1F* hSig = new TH1F();

    TH1F* projectSubtractBckg(TString, TString, Int_t, Float_t, Float_t, Float_t, Float_t, TString, TCut, TString, TString, bool);
    TH1F* subtractBckg(TH1F*, TH1F*, TString, TFile*, TString, bool);
    void peakFit(TString, TH1F*, Float_t, Float_t, Float_t, Float_t, TString, Float_t, Float_t, TString, Float_t);
    TH1F* peakFitResSub(TH1F*, Float_t, Float_t, Float_t, Float_t, TString, Float_t, Float_t, TString);
    void peakMassFit(TString, TH1F*, Float_t, Float_t, Float_t, Float_t, TString, Float_t, Float_t, TString);
    void setOutputFileName(TString);
    float getMeanError();
    float getMean();
    float getSigmaError();
    float getSigma();
    float getHeight();

    void makeTuple(TString, TCut, bool);
    void setMean(float);
    void setSigma(float);
    void setHeight(float);
private:
    Float_t mMean;
    Float_t mSigma;
    Float_t mSigmaE;
    Float_t mMeanE;
    Float_t mHeight;

ClassDef(FitPID,1)

};


inline void  FitPID::setOutputFileName(TString name) { mOutputFileName = name; }
inline float FitPID::getMean() {return mMean;}
inline float FitPID::getMeanError() {return mMeanE;}
inline float FitPID::getHeight() {return mHeight;}
inline float FitPID::getSigma() {return mSigma;}
inline float FitPID::getSigmaError() {return mSigmaE;}

inline void FitPID::setMean(float a) {mMean = a;}
inline void FitPID::setSigma(float a) {mSigma = a;}
inline void FitPID::setHeight(float a) {mHeight = a;}


#endif //DMAKER_DAU_STFIT_H


