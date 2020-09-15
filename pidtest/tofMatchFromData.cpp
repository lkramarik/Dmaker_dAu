


void tofMatchFromData(){
    Int_t colors[] = {1,46,8,4,5,6,7,9,41,42,44,45,16};
    TString fileNameData = "/media/lukas/1E4183350724EC9C/tracks.sim.vzvpd3.hotspot.root";
    auto* dataF = new TFile(fileNameData,"r");
    auto* dataOut = new TFile("tracksQA.data.new.root","recreate");
    auto* ntpData = (TNtuple*)dataF -> Get("ntp_tracks");

    double maxPt = 6;

    TString variables[] = {"dca", "nHitsFit"};
    Double_t limMin[] = {0, 0.5};
    Double_t limMax[] = {1.5, 50.5};
    Int_t nBins[] = {100, 51};
    TH1D *his;
    TH1D *hisHft;

    Float_t  pt, dca, isHFT, isTOF, zdc, nHitsFit, refMult, nHftTracks, nTofTracks, particleId, isPrimary, dcaxy, dcaz;
    ntpData->SetBranchAddress("pt", &pt);
    ntpData->SetBranchAddress("dca", &dca);
    ntpData->SetBranchAddress("dcaXy", &dcaxy);
    ntpData->SetBranchAddress("dcaZ", &dcaz);
    ntpData->SetBranchAddress("isHFT", &isHFT);
    ntpData->SetBranchAddress("isTOF", &isTOF);
    ntpData->SetBranchAddress("zdc", &zdc);
    ntpData->SetBranchAddress("nHitsFit", &nHitsFit);
    ntpData->SetBranchAddress("refMult", &refMult);
    ntpData->SetBranchAddress("isPrimary", &isPrimary);
    ntpData->SetBranchAddress("nHftTracks", &nHftTracks);
    ntpData->SetBranchAddress("nTofTracks", &nTofTracks);
    ntpData->SetBranchAddress("particleId", &particleId);

    Double_t maxDca=0.2;

//init histos
    TString histoname;
    TH1D* hnHitsFit = new TH1D("hnHitsFit", "hnHitsFit", 51, 0.5, 50.5);
    TH1D* hnHitsFitHft = new TH1D("hnHitsFitHft", "hnHitsFitHft", 51, 0.5, 50.5);
    TH1D* hDcaAll = new TH1D("hDca_all", "hDca_all", 100, 0, maxDca);
    TH1D* hDcaAllHft = new TH1D("hDca_all_hft", "hDca_all_hft", 100, 0, maxDca);
    TH1D* hDcaxyAll = new TH1D("hDcaxy_all", "hDcaxy_all", 100, -1*maxDca, maxDca);
    TH1D* hDcaxyAllHft = new TH1D("hDcaxy_all_hft", "hDcaxy_all_hft", 100, -1*maxDca, maxDca);
    TH1D* hDcazAll = new TH1D("hDcaz_all", "hDcaz_all", 100, -1*maxDca, maxDca);
    TH1D* hDcazAllHft = new TH1D("hDcaz_all_hft", "hDcaz_all_hft", 100, -1*maxDca, maxDca);

    TH1D* hPtRecoDCA[nDcaRatio];
    TH1D* hPtRecoHftDCA[nDcaRatio];

    TH1D* hDCAPtBins[nPtDca];
    TH1D* hDCAPtBinsHft[nPtDca];
    TH1D* hDCAxyPtBins[nPtDca];
    TH1D* hDCAxyPtBinsHft[nPtDca];
    TH1D* hDCAzPtBins[nPtDca];
    TH1D* hDCAzPtBinsHft[nPtDca];

    for (int j = 0; j < nDcaRatio; ++j) {
        histoname=Form("hPt_dca%.1f", dcaCut[j]);
        hPtRecoDCA[j]=new TH1D(histoname, histoname, 100,0,10);

        histoname=Form("hPt_hft_dca%.1f", dcaCut[j]);
        hPtRecoHftDCA[j]=new TH1D(histoname, histoname, 100,0,10);
    }

    for (int j = 0; j < nPtDca; ++j) {
        histoname=Form("hDca_pt_%.2f_%.2f", ptCut[j], ptCut[j+1]);
        hDCAPtBins[j]=new TH1D(histoname, histoname, 100, 0, maxDca);
        histoname=Form("hDca_hft_pt_%.2f_%.2f", ptCut[j], ptCut[j+1]);
        hDCAPtBinsHft[j]=new TH1D(histoname, histoname, 100, 0, maxDca);

        histoname=Form("hDcaxy_pt_%.2f_%.2f", ptCut[j], ptCut[j+1]);
        hDCAxyPtBins[j]=new TH1D(histoname, histoname, 100, -1*maxDca, maxDca);
        histoname=Form("hDcaxy_hft_pt_%.2f_%.2f", ptCut[j], ptCut[j+1]);
        hDCAxyPtBinsHft[j]=new TH1D(histoname, histoname, 100, -1*maxDca, maxDca);

        histoname=Form("hDcaz_pt_%.2f_%.2f", ptCut[j], ptCut[j+1]);
        hDCAzPtBins[j]=new TH1D(histoname, histoname, 100, -1*maxDca, maxDca);
        histoname=Form("hDcaz_hft_pt_%.2f_%.2f", ptCut[j], ptCut[j+1]);
        hDCAzPtBinsHft[j]=new TH1D(histoname, histoname, 100, -1*maxDca, maxDca);
    }

    const int nVars = sizeof(variables) / sizeof(TString);
    for (int k = 0; k < ntpData->GetEntries()/3; ++k) {
        ntpData->GetEntry(k);
        if (isTOF==0) continue;
        if (isPrimary==0) continue;
        if (particleId!=0) continue;
//        if (nHitsFit<20) continue;
//        if (nTofTracks<1) continue;
        if (dca>1.5) continue;

        hnHitsFit->Fill(nHitsFit);
        hDcaAll->Fill(dca);
        hDcaxyAll->Fill(dcaxy);
        hDcazAll->Fill(dcaz);

        for (int j = 0; j < nPtDca; ++j) {
            if (pt>ptCut[j] && pt<ptCut[j+1]) {
                hDCAPtBins[j]->Fill(dca);
                hDCAxyPtBins[j]->Fill(dcaxy);
                hDCAzPtBins[j]->Fill(dcaz);
            }
        }

        if (isHFT>0) {
            hnHitsFitHft->Fill(nHitsFit);
            hDcaAllHft->Fill(dca);
            hDcaxyAllHft->Fill(dcaxy);
            hDcazAllHft->Fill(dcaz);
            for (int j = 0; j < nPtDca; ++j) {
                if (pt>ptCut[j] && pt<ptCut[j+1]) {
                    hDCAPtBinsHft[j]->Fill(dca);
                    hDCAxyPtBinsHft[j]->Fill(dcaxy);
                    hDCAzPtBinsHft[j]->Fill(dcaz);
                }
            }
            for (int i = 0; i < nDcaRatio; ++i) {
                if (dca<dcaCut[i]) hPtRecoHftDCA[i]->Fill(pt);
            }
        }
    }

    TList *listOutQA = new TList();
    listOutQA->Add(hnHitsFit);
    listOutQA->Add(hnHitsFitHft);
    listOutQA->Add(hDcaAll);
    listOutQA->Add(hDcaAllHft);
    listOutQA->Add(hDcaxyAll);
    listOutQA->Add(hDcaxyAllHft);
    listOutQA->Add(hDcazAll);
    listOutQA->Add(hDcazAllHft);
    for (int l = 0; l < nPtDca; ++l) {
        listOutQA->Add(hDCAPtBins[l]);
        listOutQA->Add(hDCAPtBinsHft[l]);
        listOutQA->Add(hDCAxyPtBins[l]);
        listOutQA->Add(hDCAxyPtBinsHft[l]);
        listOutQA->Add(hDCAzPtBins[l]);
        listOutQA->Add(hDCAzPtBinsHft[l]);
    }

    for (int j = 0; j < nDcaRatio; ++j) {
        listOutQA->Add(hPtRecoDCA[j]);
        listOutQA->Add(hPtRecoHftDCA[j]);
    }

    dataOut->cd();
    listOutQA->Write("hists_QA", 1, 0);

    dataF->Close();
    dataOut->Close();
}