//
// Created by lukas on 6.12.2018.
//

void graphs(){

    TFile* data = new TFile("results_KK.root" ,"r");
    TString name[3] = {"eff", "mean", "sigma"};
    TString x = "p_T (GeV/c^2)";
    TString title;
    TString y[3] = {"Efficiency", "n#sigma_K sigma", "n#sigma_K sigma"};
    TGraph* gr[3] = {new TGraph(), new TGraph(), new TGraph()};

    for (int i = 0; i < 3; ++i) {
        gr[i] =(TGraph*) data->Get(name[i]);
//        gr[i]->
        title = " ;"+x+";"+y[i];
//        gr[i]->SetTitle(" ;Invariant mass, m_{K#pi} [GeV/c^{2}];Counts");
        gr[i]->SetTitle(title);

    }

    TCanvas *c1 = new TCanvas("c1","c1",1000,900);
    gr[0]->Draw("apze");

    TCanvas *c2 = new TCanvas("c2","c2",1000,900);
    gr[1]->Draw("apze");

    TCanvas *c3 = new TCanvas("c3","c3",1000,900);
    gr[2]->Draw("apze");







}
