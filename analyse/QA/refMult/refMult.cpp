

void refMult() {
//    TString input = "outputLocal.picoQAAnaMaker.root";
    TString input = "qa.refmult.3001.root";
//    TString input = "qa.2907.gref.root";

    TFile *inFile = new TFile(input, "READ");
    TList *list = (TList*)inFile->Get("picoQAMaker;1");

    TString names[3]={"", "_HFT", "_HFT_hybridTOF"};
    TString legend[3]={"no req. on tracks", "#(HFT)>1", "#(HFT && hybridTOF)>1"};

    TCanvas *c3 = new TCanvas("c3", "c3", 1000, 900);
    TCanvas *c4 = new TCanvas("c4", "c4", 1000, 900);

    TLegend *legendZDC = new TLegend(0.6, 0.67, 0.77, 0.88, "", "brNDC");
    legendZDC->SetFillStyle(0);
    legendZDC->SetLineColor(0);
    legendZDC->SetTextSize(0.04);
    TString textData = "dAu 200 GeV, Run16, trigger VPD-5, |V_{z}|< 6 cm, |V_{z,TPC}-V_{z,VPD}|< 6 cm, SL18f";
    TPaveText *text1 = new TPaveText(0.439, 0.90, 0.91, 0.947, "NDC");
    text1->AddText(textData);
    text1->SetFillColor(0);

    Int_t col[]={1,46,8};
    for (int k=0; k<3;++k) {
        TH2F *hgrefZdc = static_cast<TH2F *>(list->FindObject(Form("h_gRefmult_vs_ZDCx%s", names[k].Data())));
        hgrefZdc->Sumw2();

        TProfile *pxy = new TProfile();
        pxy = hgrefZdc->ProfileX(legend[k]);
        pxy->SetMarkerStyle(2);
        pxy->SetMarkerColor(col[k]);
        pxy->SetLineColor(col[k]);
        pxy->GetXaxis()->SetTitle("ZDC rate [kHz]");
        pxy->GetYaxis()->SetTitle("<grefMult>");
        pxy->SetTitle("");
        pxy->SetStats(0);
        pxy->GetXaxis()->SetLabelFont(42);
        pxy->GetXaxis()->SetRangeUser(50,250);
        pxy->GetYaxis()->SetRangeUser(0,21);
        pxy->GetXaxis()->SetTitleFont(42);
        c3->cd();
        if (k==0) pxy->DrawClone();
        else pxy->DrawClone("same");
        legendZDC->AddEntry(pxy, legend[k], "pl");

        if (k==2) {
//            legendZDC->Draw("same");
            c3->BuildLegend();
            text1->Draw("same");
        }

        TH2F *hgrefOne = static_cast<TH2F *>(list->FindObject(Form("h_gRefmult%s", names[k].Data())));
        hgrefOne->Sumw2();
        TH1D *pxy1 = new TH1D();
        pxy1 = hgrefOne->ProjectionX();
//        pxy1 = hgrefOne->ProjectionY(legend[k]);
        pxy1->SetMarkerStyle(2);
        pxy1->SetMarkerColor(col[k]);
        pxy1->SetLineColor(col[k]);
        pxy1->GetXaxis()->SetTitle("grefMult");
        pxy1->GetYaxis()->SetTitle("1/Entries");
        pxy1->SetTitle("");
        pxy1->SetStats(0);
        pxy1->Scale(1/pxy1->GetEntries());
        pxy1->GetXaxis()->SetLabelFont(42);
        pxy1->GetXaxis()->SetRangeUser(0,100);
//        pxy1->GetYaxis()->SetRangeUser(0,21);
        pxy1->GetXaxis()->SetTitleFont(42);
        pxy1->SetName(legend[k]);

        c4->cd();
        if (k==0) pxy1->DrawClone();
        else pxy1->DrawClone("same");

        if (k==2) {
            c4->BuildLegend();
            text1->Draw("same");
        }


        TString vzInt[] = {"min6", "min4", "min2", "0", "2", "4", "6"};
        Int_t vzRange[] = {-6, -4, -2, 0, 2, 4, 6};
        Int_t color[] = {1, 46, 8, 9, 6, 7};
        TCanvas *cgRef = new TCanvas("cgRef", "cgRef", 1200, 900);
        cgRef->Divide(3, 2, 0.01, 0.01, 0);
        TH1F *hgref[6];
        TF1 *erfF[6];
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
            erfF[i]=new TF1("erf", "[0]*TMath::Erf([1]*(x-[2]))+[0]", 30, 95);
            erfF[i]->SetParName(0, "A");
            erfF[i]->SetParName(1, "sigma");
            erfF[i]->SetParName(2, "h");
            erfF[i]->SetParLimits(0, 3000000, 10000000);
            erfF[i]->SetParLimits(1, -1, 1);
            erfF[i]->SetParLimits(2, -100, 100);

            if (i==0) erfF[i]->SetParameters(30000,0.5,80);
            else erfF[i]->SetParameters(erfF[i-1]->GetParameter(0), erfF[i-1]->GetParameter(1), erfF[i-1]->GetParameter(2));

//            hgref[i]->Fit(erfF[i], "RLMI");
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
            hgref[j]->GetYaxis()->SetTitle("1/Entries");
            hgref[j]->GetXaxis()->SetRangeUser(0,60);
            hgref[j]->SetTitle("");
            hgref[j]->Scale(1/hgref[j]->GetEntries());
            if (j == 1) hgref[j]->Draw();
            else hgref[j]->Draw("same");

            legendPub->AddEntry(hgref[j], Form("%i cm < V_{z} < %i cm", vzRange[j], vzRange[j + 1]), "pl");
            legendPub->SetFillStyle(0);
            legendPub->SetLineColor(0);
            legendPub->SetTextSize(0.04);
        }
        legendPub->Draw("same");
        text1->Draw("same");

//        cgRefOne->SetLogy();
//        cgRefOne->Update();
        cgRefOne->SaveAs(Form("grefOne_Vz%s.png",names[k].Data()));
        cgRefOne->Close();
    }

//    c3->SaveAs(Form("gref_vs_ZDC%s.png",names[k].Data()));
    c3->SaveAs("gref_vs_ZDC.png");
    c3->Close();
    delete c3;

    TFile *inFileSim = new TFile("/home/lukas/work/glauber/ncoll_npart.root", "READ");
    TH1F *hMult = static_cast<TH1F *>(inFileSim->Get("hMult"));
    hMult->Sumw2();
    hMult->Scale(1/hMult->GetEntries());
    hMult->Draw();
    TH2F *hgrefOne = static_cast<TH2F *>(list->FindObject(Form("h_gRefmult%s", names[2].Data())));
    TH1D *pxy1 = new TH1D();
    pxy1 = hgrefOne->ProjectionX();
    pxy1->Scale(1/pxy1->GetEntries());
    pxy1->Draw("same");

}