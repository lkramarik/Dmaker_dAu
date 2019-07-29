

void refMult() {
    TString input = "qa.2907.gref.root";

    TFile *inFile = new TFile(input, "READ");
    TList *list = (TList*)inFile->Get("picoQAMaker;1");

    TH2F *hgrefZdc = static_cast<TH2F*>(list->FindObject("h_gRefmult_vs_ZDCx"));
    hgrefZdc->Sumw2();

    TProfile *pxy=new TProfile();
    pxy = hgrefZdc->ProfileX("pxy");
    pxy->SetMarkerStyle(2);
    pxy->SetMarkerColor(kBlack);
    pxy->GetXaxis()->SetTitle("ZDC rate [kHz]");
    pxy->GetYaxis()->SetTitle("<grefMult>");
    pxy->SetTitle("");
    pxy->SetStats(0);
    pxy->GetXaxis()->SetLabelFont(42);
    pxy->GetXaxis()->SetTitleFont(42);

    TCanvas *c3 = new TCanvas("c3","c3",1200,900);
    pxy->Draw();

    TString vzInt[]={"min6","min4","min2","0", "2", "4", "6"};
    Int_t vzRange[]={-6, -4, -2, 0, 2, 4, 6};
    Int_t color[]={1,46,8,9,6,7};
    TCanvas *cgRef = new TCanvas("cgRef","cgRef",1200,900);
    cgRef->Divide(3,2,0.01,0.01,0);
    TH1F *hgref[6];
    for (int i = 0; i < 6; ++i) {
        cgRef->cd(i+1);
        hgref[i] = static_cast<TH1F*>(list->FindObject(Form("h_gRefmult_Vz_%s_%s", vzInt[i].Data(), vzInt[i+1].Data())));
        hgref[i]->Sumw2();
        hgref[i]->SetStats(0);
        hgref[i]->SetTitle(Form("%i cm < V_{z} < %i cm", vzRange[i], vzRange[i+1]));
        hgref[i]->GetXaxis()->SetTitle("grefMult");
        hgref[i]->GetYaxis()->SetTitle("Counts");
        hgref[i]->GetYaxis()->SetTitleOffset(1.1);
        hgref[i]->SetMarkerStyle(2);
        hgref[i]->DrawClone();

    }

    TCanvas *cgRefOne = new TCanvas("cgRefOne","cgRefOne",1200,900);
    TLegend *legendPub = new TLegend(0.6, 0.67, 0.77, 0.88, "", "brNDC");

    for (int j = 0; j < 6; ++j) {
        hgref[j]->SetMarkerColor(color[j]);
        hgref[j]->SetLineColor(color[j]);
        hgref[j]->SetTitle("");
        if(j==1) hgref[j]->Draw();
        else hgref[j]->Draw("same");

        legendPub->AddEntry(hgref[j], Form("%i cm < V_{z} < %i cm", vzRange[j], vzRange[j+1]), "pl");
        legendPub->SetFillStyle(0);
        legendPub->SetLineColor(0);
        legendPub->SetTextSize(0.04);
    }
    legendPub->Draw("same");




}