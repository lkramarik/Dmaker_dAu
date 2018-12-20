#include <iostream>
#include "TProfile.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TMath.h"
#include "TH1.h"

using namespace std;

void flow()
{
 	TFile *file=new TFile("output_noHybTOF.root");
 	TProfile *qVec[4];
 	TProfile *refFlow;
    TProfile *dirFlow[5];
    TProfile *corrD[2][5];
    TProfile *qVec2[4];
    TProfile *refFlow2;
    TProfile *dirFlow2;
    TProfile *corrD2[2];
    
    //TProfile *c2;
    TH1D *c2;
    TH1D *d2[5];
    TH1D *v2[5];
    
    TH1D *cum;
    TH1D *cum_noC;
    TH1D *d2_all_mult;
    TH1D *v2_all_mult;
    TH1D *d2_all_mult_noC;
    TH1D *v2_all_mult_noC;

    //loading TProfiles
 	TString names[4] = {"cos_B", "cos_F", "sin_B", "sin_F"}; //backward and forward samples
  	float multBin[6] = {0,7,12,16,22,100};
  	for(int m = 0; m < 4; m++)
  	{
	  	TString aa = names[m];
	  	qVec[m] = (TProfile *)file->Get(aa.Data());
	  	aa += "_no_mult";
	  	qVec2[m] = (TProfile *)file->Get(aa.Data());
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

  	corrD2[0] = (TProfile *)file->Get("cosD_no_mult");
  	corrD2[1] = (TProfile *)file->Get("sinD_no_mult");
	dirFlow2 = (TProfile *)file->Get("dirFlow_no_mult");
  	refFlow2 = (TProfile *)file->Get("refFlow_no_mult");

  	
	double d,r,c,s,final,error,cumulant,cumEr;
	//c2 = new TProfile("c2","c_2{2}", 5, multBin);
	c2 = new TH1D("c2", "v2 ref = sqrt(c_2{2})", 5, multBin);
	cum = new TH1D("cum", "v2 ref = sqrt(c_2{2})", 1, 0, 100);
	cum_noC = new TH1D("cum_noC", "v2 ref = sqrt(c_2{2})", 1, 0, 100);
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

	r = refFlow2->GetBinContent(1);
	c = qVec2[0]->GetBinContent(1);
	s = qVec2[2]->GetBinContent(1);
	final = r - c*c - s*s;
	error = TMath::Sqrt( refFlow2->GetBinError(1)*refFlow2->GetBinError(1) + 4*c*c*qVec2[0]->GetBinError(1)*qVec2[0]->GetBinError(1) +  4*s*s*qVec2[1]->GetBinError(1)*qVec2[1]->GetBinError(1));
	cum->SetEntries(refFlow2->GetEntries());
	cum->SetBinContent(1, (TMath::Sqrt(final)));
	//cum->SetBinError(1, refFlow2->GetBinError(1));
	cum->SetBinError(1, error);
	cum_noC->SetEntries(refFlow2->GetEntries());
	cum_noC->SetBinContent(1, (TMath::Sqrt(r)));
	cum_noC->SetBinError(1, refFlow2->GetBinError(1));

	printf("error just from ref flow %f \n", refFlow2->GetBinError(1));
	printf(" fancy error %f \n", error);

	cum->Draw();
	new TCanvas;

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
				s = (corrD[1][i]->GetBinContent(j))*(qVec[2]->GetBinContent(i+1));
				final = d - c - s;
				d2[i]->SetBinContent(j, final);
				d2[i]->SetBinError(j, dirFlow[i]->GetBinError(j));
				v2[i]->SetBinContent(j,(d2[i]->GetBinContent(j))/(c2->GetBinContent(i+1)));
				v2[i]->SetBinError(j,(d2[i]->GetBinError(j))/(c2->GetBinContent(i+1)));
	//			allm2->Fill(momBins[j-1]+0.5,v2[i]->GetBinContent(j),v2[i]->GetEntries());
			}
		}	

	d2_all_mult = new TH1D("d2_all_mult", "d_2{2}", 6, momBins);
	v2_all_mult = new TH1D("v2_all_mult", "v_{2};p_{T};v_{2}", 6, momBins);

	d2_all_mult_noC = new TH1D("d2_all_mult_noC", "d_2{2}", 6, momBins);
	v2_all_mult_noC = new TH1D("v2_all_mult_noC", "v_{2};p_{T};v_{2}", 6, momBins);

	cumulant = cum->GetBinContent(1);
	cumEr = cum->GetBinError(1);

	double dirEr,dirF;
	
	for(int j = 2; j < 6; j++)
			{
				d = dirFlow2->GetBinContent(j);
				c = (corrD2[0]->GetBinContent(j))*(qVec2[0]->GetBinContent(1));
				s = (corrD2[1]->GetBinContent(j))*(qVec2[2]->GetBinContent(1));
				final = d - c - s;
				d2_all_mult->SetBinContent(j, final);
				d2_all_mult->SetBinError(j, dirFlow2->GetBinError(j));
				d2_all_mult_noC->SetBinContent(j, d);
				d2_all_mult_noC->SetBinError(j, dirFlow2->GetBinError(j));
				d2_all_mult->SetBinError(j,  TMath::Sqrt( dirFlow2->GetBinError(j)*dirFlow2->GetBinError(j) + corrD2[0]->GetBinContent(j)*corrD2[0]->GetBinContent(j)*qVec2[0]->GetBinError(1)*qVec2[0]->GetBinError(1) + qVec2[0]->GetBinContent(1)*qVec2[0]->GetBinContent(1)*corrD2[0]->GetBinError(j)*corrD2[0]->GetBinError(j) + corrD2[1]->GetBinContent(j)*corrD2[1]->GetBinContent(j)*qVec2[2]->GetBinError(1)*qVec2[2]->GetBinError(1) + qVec2[2]->GetBinContent(1)*qVec2[2]->GetBinContent(1)*corrD2[1]->GetBinError(j)*corrD2[1]->GetBinError(j) ) );
				printf("fancy error %f \n", d2_all_mult->GetBinError(j));
				printf("normal error from ref flow %f \n", dirFlow2->GetBinError(j));

				dirEr = d2_all_mult->GetBinError(j);
				dirF = d2_all_mult->GetBinContent(j);
				//error = TMath::Sqrt((dirEr*dirEr)/(cumulant*cumulant) +  (dirF*dirF)*(cumEr*cumEr)/(TMath::Power(cumulant,4)));
				v2_all_mult->SetBinContent(j,(d2_all_mult->GetBinContent(j))/(cum->GetBinContent(1)));
				v2_all_mult->SetEntries(dirFlow2->GetEntries());
				error = v2_all_mult->GetBinContent(j)*TMath::Sqrt( (cumEr/cumulant)*(cumEr/cumulant) + (dirEr/dirF)*(dirEr/dirF) );
				//v2_all_mult->SetBinError(j,(d2_all_mult->GetBinError(j))/(cum->GetBinContent(1)));
				v2_all_mult->SetBinError(j,error);
				v2_all_mult_noC->SetEntries(dirFlow2->GetEntries());
				v2_all_mult_noC->SetBinContent(j,(d2_all_mult_noC->GetBinContent(j))/(cum_noC->GetBinContent(1)));
				v2_all_mult_noC->SetBinError(j,(d2_all_mult_noC->GetBinError(j))/(cum_noC->GetBinContent(1)));
			}

	

	TFile *output = new TFile("v2_plots_noHybTOF.root", "recreate");
	c2->Write();
	for (int i = 0; i < 5; i++)
	{
		d2[i]->Write();
		v2[i]->Write();
	}
	cum->Write();
	cum_noC->Write();
	v2_all_mult->Write();
	d2_all_mult->Write();
	v2_all_mult_noC->Write();
	d2_all_mult_noC->Write();

	v2_all_mult->Draw();
	v2_all_mult_noC->SetLineColor(kRed);
	v2_all_mult_noC->Draw("same");
	
 }

 	