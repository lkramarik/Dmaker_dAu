

void refMult() {
//    TString input = "outputLocal.picoQAAnaMaker.root";
    TString input = "qa.refmult.3001.root";
//    TString input = "qa.2907.gref.root";

    TFile *inFile = new TFile(input, "READ");
    TList *list = (TList*)inFile->Get("picoQAMaker;1");

    TString names[3]={"", "_HFT", "_HFT_hybridTOF"};

    for (int k=0; k<3;++k) {
        TH2F *hgrefZdc = static_cast<TH2F *>(list->FindObject(Form("h_gRefmult_vs_ZDCx%s", names[k].Data())));
        hgrefZdc->Sumw2();

        TProfile *pxy = new TProfile();
        pxy = hgrefZdc->ProfileX("pxy");
        pxy->SetMarkerStyle(2);
        pxy->SetMarkerColor(kBlack);
        pxy->GetXaxis()->SetTitle("ZDC rate [kHz]");
        pxy->GetYaxis()->SetTitle("<grefMult>");
        pxy->SetTitle("");
        pxy->SetStats(0);
        pxy->GetXaxis()->SetLabelFont(42);
        pxy->GetXaxis()->SetRangeUser(50,250);
        pxy->GetYaxis()->SetRangeUser(0,21);
        pxy->GetXaxis()->SetTitleFont(42);

        TCanvas *c3 = new TCanvas("c3", "c3", 1000, 900);
        pxy->Draw();
        c3->SaveAs(Form("gref_vs_ZDC%s.png",names[k].Data()));
        c3->Close();
        delete c3;

        TString vzInt[] = {"min6", "min4", "min2", "0", "2", "4", "6"};
        Int_t vzRange[] = {-6, -4, -2, 0, 2, 4, 6};
        Int_t color[] = {1, 46, 8, 9, 6, 7};
        TCanvas *cgRef = new TCanvas("cgRef", "cgRef", 1200, 900);
        cgRef->Divide(3, 2, 0.01, 0.01, 0);
        TH1F *hgref[6];
        for (int i = 0; i < 6; ++i) {
            cgRef->cd(i + 1);
            gPad->SetLogy();
            hgref[i] = static_cast<TH1F *>(list->FindObject(Form("h_gRefmult_Vz_%s_%s%s", vzInt[i].Data(), vzInt[i + 1].Data(), names[k].Data())));
            hgref[i]->Sumw2();
            hgref[i]->SetStats(0);
            hgref[i]->SetTitle(Form("%i cm < V_{z} < %i cm", vzRange[i], vzRange[i + 1]));
            hgref[i]->GetXaxis()->SetTitle("grefMult");
            hgref[i]->GetYaxis()->SetTitle("Counts");
            hgref[i]->GetYaxis()->SetTitleOffset(1.1);
            hgref[i]->SetMarkerStyle(2);
            hgref[i]->DrawClone();
        }
        cgRef->SetLogy();
        cgRef->SaveAs(Form("grefAll_Vz%s.png",names[k].Data()));
        cgRef->Close();

        TCanvas *cgRefOne = new TCanvas("cgRefOne", "cgRefOne", 1200, 900);
        TLegend *legendPub = new TLegend(0.6, 0.67, 0.77, 0.88, "", "brNDC");

        for (int j = 0; j < 6; ++j) {
            hgref[j]->SetMarkerColor(color[j]);
            hgref[j]->SetLineColor(color[j]);
            hgref[j]->SetTitle("");
            if (j == 1) hgref[j]->Draw();
            else hgref[j]->Draw("same");

            legendPub->AddEntry(hgref[j], Form("%i cm < V_{z} < %i cm", vzRange[j], vzRange[j + 1]), "pl");
            legendPub->SetFillStyle(0);
            legendPub->SetLineColor(0);
            legendPub->SetTextSize(0.04);
        }
        legendPub->Draw("same");
        cgRefOne->SetLogy();
//        cgRefOne->Update();
        cgRefOne->SaveAs(Form("grefOne_Vz%s.png",names[k].Data()));
        cgRefOne->Close();
    }


}