#include <iostream>
#include "TProfile.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TMath.h"
#include "TH1.h"

using namespace std;

void flow()
{
 	TFile *file=new TFile("output.root");
 	TProfile *qVec[4];
 	TProfile *refFlow;
    TProfile *dirFlow[5];
    TProfile *corrD[2][5];
    //TProfile *c2;
    TH1D *c2;
    TH1D *d2[5];
    TH1D *v2[5];
    TH1D *allm;
    TProfile *allm2;

    //loading TProfiles
 	TString names[4] = {"cos_B", "cos_F", "sin_B", "sin_F"}; //backward and forward samples
  	float multBin[6] = {0,7,12,16,22,100};
  	for(int m = 0; m < 4; m++)
  	{
	  	TString aa = names[m];
	  	qVec[m] = (TProfile *)file->Get(aa.Data());
  	}
 	float momBins[7] = {0,1,2,3,4,5,10};
  	TString multBinNames[6] = {"0","7","12","16","22","100"};
  	for(int m = 0; m < 5; m++)
  	{
	  	TString aa = "cosD_" + multBinNames[m] + "_" + multBinNames[m+1];
	  	corrD[0][m] = (TProfile *)file->Get(aa.Data());
	  	aa = "sinD_" + multBinNames[m] + "_" + multBinNames[m+1];
	  	corrD[1][m] = (TProfile *)file->Get(aa.Data());
	  	aa = "dirFlow_" + multBinNames[m] + "_" + multBinNames[m+1];
	  	dirFlow[m] = (TProfile *)file->Get(aa.Data());
  	}
  	refFlow = (TProfile *)file->Get("refFlow");

  	allm = new TH1D("allm","All multiplicities v2",6,momBins);
  	allm2 = new TProfile("allm","All multiplicities v2",6,momBins);


	double d,r,c,s,final;
	//c2 = new TProfile("c2","c_2{2}", 5, multBin);
	c2 = new TH1D("c2", "v2 ref = sqrt(c_2{2})", 5, multBin);
	for(int i = 1; i < 6; i++)
		{
			//computing c_2{2} via reference flow
			r = refFlow->GetBinContent(i);
			c = qVec[0]->GetBinContent(i);
			s = qVec[2]->GetBinContent(i);
			final = r - c*c - s*s;
			c2->SetEntries(refFlow->GetEntries());
			c2->SetBinContent(i, (TMath::Sqrt(final)));
			c2->SetBinError(i, refFlow->GetBinError(i));
		}

	for (int i = 0; i < 5; i++)
		{
			//computing v2 via directed flow
			d2[i] = new TH1D(TString::Format("d2_%d", i), "d_2{2}", 6, momBins);
			v2[i] = new TH1D(TString::Format("v2_%d", i), "v_{2};p_{T};v_{2}", 6, momBins);
			d2[i]->SetEntries(dirFlow[i]->GetEntries());
			v2[i]->SetEntries(dirFlow[i]->GetEntries());
			for(int j = 1; j < 7; j++)
			{
				d = dirFlow[i]->GetBinContent(j);
				c = (corrD[0][i]->GetBinContent(j))*(qVec[0]->GetBinContent(i+1));
				s = (corrD[0][i]->GetBinContent(j))*(qVec[2]->GetBinContent(i+1));
				final = d - c - s;
				d2[i]->SetBinContent(j, final);
				d2[i]->SetBinError(j, dirFlow[i]->GetBinError(j));
				v2[i]->SetBinContent(j,(d2[i]->GetBinContent(j))/(c2->GetBinContent(i+1)));
				v2[i]->SetBinError(j,(d2[i]->GetBinError(j))/(c2->GetBinContent(i+1)));
				allm2->Fill(momBins[j-1]+0.5,v2[i]->GetBinContent(j),v2[i]->GetEntries());
			}
			v2[i]->SetBit(TH1::kIsAverage);
			allm->SetBit(TH1::kIsAverage);
			allm->Add(v2[i]);
			if(i==0) allm->Sumw2();
		}	

	

	TFile *output = new TFile("v2_plots.root", "recreate");
	c2->Write();
	for (int i = 0; i < 5; i++)
	{
		d2[i]->Write();
		v2[i]->Write();
	}
	allm->Write();
	allm->Draw();
	allm2->Write();
	new TCanvas;
	allm2->Draw();

 }

 	