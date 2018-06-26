#include "StMixerHists.h"

StMixerHists::StMixerHists(char* fileBaseName):
  mSE_Vtx(NULL), mME_Vtx(NULL), mSE_LS(NULL), mSE_US_plus(NULL), mSE_US_minus(NULL),
  mME_LS(NULL), mME_US_plus(NULL), mME_US_minus(NULL)
{
  const float massMin = 0.;
  const float massMax = 3.;
  const float massNBin = 300;

  mSE_Vtx = new TH2F(Form("%s_seVtx",fileBaseName),"Vertex pos;vertex x;vertex y",250,-2.5,2.5,250,-2.5,2.5);
  mSE_Vtx->Sumw2();
  mME_Vtx = new TH2F(Form("%s_meVtx",fileBaseName),"Vertex pos;vertex x;vertex y",250,-2.5,2.5,250,-2.5,2.5);
  mME_Vtx->Sumw2();

  mSE_LS = new TH2F(Form("%s_se_ls_mass",fileBaseName),"Same Event wrong sign triplet candidates Invariant mass(K#pi^{}p);p_{T}(K#pi^{}p)(GeV/c),Mass_{K#pi^{p}}(GeV/c^{2})",150,0,15,massNBin,massMin,massMax);
  mSE_LS->Sumw2();
  mSE_US_plus = new TH2F(Form("%s_se_us_plus_mass",fileBaseName),"Same Event good sign #Lambda_{c}^{+} triplet candidates Invariant mass(K^{-}#pi^{+}p^{+});p_{T}(K#pi^{}p)(GeV/c),Mass_{K#pi^{}p}(GeV/c^{2})",150,0,15,massNBin,massMin,massMax);
  mSE_US_minus = new TH2F(Form("%s_se_us_minus_mass",fileBaseName),"Same Event good sign #Lambda_{c}^{-} triplet candidates Invariant mass(K^{+}#pi^{-}p^{-});p_{T}(K#pi^{}p)(GeV/c),Mass_{K#pi^{}p}(GeV/c^{2})",150,0,15,massNBin,massMin,massMax);
  mSE_US_plus->Sumw2();
  mSE_US_minus->Sumw2();
  mME_LS = new TH2F(Form("%s_me_ls_mass",fileBaseName),"Mixed Event wrong sign Invariant mass(K#pi^{}p);p_{T}(K#pi^{}p)(GeV/c),Mass_{K#pi^{}p}(GeV/c^{2})",150,0,15,massNBin,massMin,massMax);
  mME_LS->Sumw2();
  mME_US_plus = new TH2F(Form("%s_me_us_plus_mass",fileBaseName),"Mixed Event good sign #Lambda_{c}^{+} triplet Invariant mass(K^{-}#pi^{+}p^{+});p_{T}(K#pi^{}p)(GeV/c),Mass_{K#pi^{}p}(GeV/c^{2})",150,0,15,massNBin,massMin,massMax);
  mME_US_minus = new TH2F(Form("%s_me_us_minus_mass",fileBaseName),"Mixed Event good sign #Lambda_{c}^{-} triplet Invariant mass(K^{+}#pi^{-}p^{-});p_{T}(K#pi^{}p)(GeV/c),Mass_{K#pi^{}p}(GeV/c^{2})",150,0,15,massNBin,massMin,massMax);
  mME_US_plus->Sumw2();
  mME_US_minus->Sumw2();
}
StMixerHists::~StMixerHists()
{
  delete mSE_Vtx ;
  delete mME_Vtx ;
  delete mSE_LS ;
  delete mSE_US_plus ;
  delete mSE_US_minus ;
  delete mME_LS ;
  delete mME_US_plus ;
  delete mME_US_minus ;
}
void StMixerHists::closeFile()
{
  mSE_Vtx->Write();
  mME_Vtx->Write();
  mSE_LS->Write();
  mSE_US_plus->Write();
  mSE_US_minus->Write();
  mME_LS->Write();
  mME_US_plus->Write();
  mME_US_minus->Write();
}
