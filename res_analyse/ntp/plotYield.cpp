//
// Created by lukas on 10.4.2018.
//

void plotYield(){
    std::fstream inFile("raw_yields.txt", std::ios_base::in);
    std::fstream inFileF("raw_yields_fct.txt", std::ios_base::in);
    std::fstream inFileEff("raw_yields_eff.txt", std::ios_base::in);
    const int nPtBins = 3;
    double x[nPtBins], xMin, xMax, xe[nPtBins], y[nPtBins], ye[nPtBins];
    double xF[nPtBins], yF[nPtBins], xFe[nPtBins], yFe[nPtBins], corrYield[nPtBins], eff[nPtBins], corrYieldErr[nPtBins], effE[nPtBins];

    float NEvents = 124455238;
    float NCharge = 2;
    float deltay = 2;
    float BR = 0.0393 + 0.0001398;
    float BRaerr = sqrt(0.0004*0.0004 + 0.027e-4*0.027e-4);

    for (int i = 0; i < 3; ++i) {
        inFile>>xMin>>xMax>>y[i]>>ye[i];
        x[i]=(xMin+xMax)/2;
        xe[i]=(xMax-xMin)/2;
        cout<<x[i]<<" "<<xe[i]<<" "<<y[i]<<" "<<ye[i]<<endl;

        inFileF>>xMin>>xMax>>yF[i]>>yFe[i];
        xF[i]=(xMin+xMax)/2;
        xFe[i]=(xMax-xMin)/2;
        cout<<xF[i]<<" "<<xFe[i]<<" "<<yF[i]<<" "<<yFe[i]<<endl;
        cout<<" "<<endl;

        cout<<"getting effs. from file..."<<endl;
        inFileEff>>xMin>>xMax>>eff[i]>>effE[i];

//        effE[i]=0.0004;
        cout<<"so...you are working in "<<xMin<<" to "<<xMax<<" and your eff is "<<eff[i]<<endl;
        cout<<"going to calculate inv yield..."<<endl;
        corrYield[i] = y[i]/(2*3.141592*NEvents*NCharge*BR*deltay*x[i]*xe[i]*2*eff[i]);
        corrYieldErr[i] = ye[i]*corrYield[i]/y[i];

    }


    TGraphErrors* gr = new TGraphErrors(nPtBins, x, y, xe, ye);
    gr->SetNameTitle("raw_yield","Raw yield from bin counting");
//    gr->SetMarkerStyle(2);
//    gr->GetXaxis()->SetTitle("p_{T} [GeV/c]");
//    gr->GetYaxis()->SetTitle("D^{0} raw yield");
    gr->SetFillStyle(0);
    gr->SetFillColor(0);

    TGraphErrors* grF = new TGraphErrors(nPtBins, xF, yF, xFe, yFe);
    grF->SetNameTitle("raw_yield_fct","Raw yield from gaus integral");
    grF->SetMarkerStyle(2);
    grF->SetMarkerColor(46);
    grF->SetLineColor(46);
    grF->SetFillStyle(0);
    grF->SetFillColor(0);

    TFile *fOut = new TFile("results.root","update");

    TCanvas *c31 = new TCanvas("c31","c31",1000,1200);
    auto mg = new TMultiGraph();
    mg->SetTitle(" ;p_{T} [GeV/c];D^{0} raw yield");
    mg->Add(gr);
    mg->Add(grF);
    mg->Draw("ap");
    c31->BuildLegend();
    c31->Update();



    TGraphErrors* effGr = new TGraphErrors(nPtBins, x, eff, xe, effE);
    TCanvas *c42 = new TCanvas("c42","c42",1000,1200);
    c42 -> SetLogy();
    auto mgEff = new TMultiGraph();
    mgEff->SetTitle(" ;p_{T} [GeV/c];D^{0} reconstruction eff.");
    mgEff->Add(effGr);
//    mg->Add(grF);
    mgEff->Draw("ap");
//    mgEff->Write("inv_yields");
    c42->SaveAs("eff.png");


    TGraphErrors* grEff = new TGraphErrors(nPtBins, x, corrYield, xe, corrYieldErr);
    TCanvas *c41 = new TCanvas("c41","c41",1000,1200);
    c41->SetLogy();
    auto mgYield = new TMultiGraph();
    mgYield->SetTitle(" ;p_{T} [GeV/c]; d^{2}N/(2#pip_{T}dp_{T}dy) [(GeV/c)^{-2}]");
    mgYield->Add(grEff);
//    mg->Add(grF);
    mgYield->Draw("ap");
    c41->SaveAs("invariant_yield.png");
//    mgYield->Write("inv_yields");

    TFile* fpp = new TFile("published_run10_D0_AuAu_data.root","r"); //2010
    TF1 *ppref = (TF1*)fpp->Get("Levy_pp");

    TGraphErrors* grRAA = new TGraphErrors(nPtBins, x, corrYield, xe, corrYieldErr);

    for (int i=0;i<grRAA->GetN();i++){
        cout<<grRAA->GetY()[i]<<endl;
         grRAA->GetY()[i] /= (0.5*7*ppref->Eval(grRAA->GetX()[i]));
//                grRAA[cent]->GetEY()[i] = yEppf[cent][i]*raappsys[cent]->GetY()[i];
    }


    TCanvas *c51 = new TCanvas("c51","c51",1000,1200);
    auto mgRAA = new TMultiGraph();
//    mgYield->SetTitle(" ;p_{T} [GeV/c];D^{0} inv yield");
    mgRAA->Add(grRAA);
//    mg->Add(grF);
    mgRAA->Draw("ap");
//    mgRAA->Write("raa");

    fOut->Close();

//	CorrYieldErrSys = CorrYield*TMath::Sqrt(pow(simerr, 2) + pow(BRerr , 2)); //incomplete


}
