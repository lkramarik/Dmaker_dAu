void fitEffRCF() {
    int bbcMin=0;
    int bbcMax=950;
    int nTof=0;
    float nsigma=3;
    float tofInvBeta=0.03;
    float ptTrackCut=0.;
    TString particle = "pi";

    TString cutComb=Form("bbc%i_%i_nHft%i_nsigma%.1f_tof%.2f_pt%.1f/",bbcMin, bbcMax, nTof, nsigma, tofInvBeta, ptTrackCut);

    TFile *fPID = new TFile("results_total_eff/"+cutComb+"totalEff_"+particle+".root", "READ");
    TGraphErrors *gPID = (TGraphErrors*) fPID->Get("grTotalGraphEffPid_"+particle);

    TF1 *fitF1 = new TF1("eff_fit", "[0]+[1]/x+[2]/x/x", 0.15, 4); //momentum resolution fit
//    TF1 *fitF1 = new TF1("eff_fit", "[0]+[1]/x+[2]*x+[3]/x/x", 0.15, 4); //momentum resolution fit
//    fitF1->SetParameters(1, -0.06, -0.1, 0.02, 0.006);
//    fitF1->SetParLimits(0, 0, 1);
//    fitF1->SetParLimits(1, 0, 1);
//    fitF1->SetParLimits(2, 0, 1);
//    fitF1->SetParLimits(3, -1, 0);

    fitF1->SetParameters(1, -0.06, -0.1, 0.006);

    fitF1->SetParLimits(0, 0, 1);
    fitF1->SetParLimits(1, 0, 1);
    fitF1->SetParLimits(2, -1, 0);

    gPID->Fit(fitF1, "", "", 0.15, 3);
//    gPID->Draw("ap");
    TFile *fOut = new TFile("totalEff_"+particle+".root", "RECREATE");
    fOut->SetCompressionAlgorithm(1);
    fitF1->Write("fTotalGraphEffPid_"+particle);

    fOut->Close();

    fPID->Close();




}

