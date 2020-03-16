
void changeCompression(){
    int bbcMin=0;
    int bbcMax=950;
    int nTof=0;
    float nsigma=3;
    float tofInvBeta=0.03;
    float ptTrackCut=0.;
    TString particle = "K";

    TString cutComb=Form("bbc%i_%i_nHft%i_nsigma%.1f_tof%.2f_pt%.1f/",bbcMin, bbcMax, nTof, nsigma, tofInvBeta, ptTrackCut);
//    TString cutComb1=Form("bbc%i_%i_nHft%i_nsigma%.1f_tof%.2f_pt%.1f/",bbcMin, bbcMax, nTof, nsigma, tofInvBeta, 0.);
    TString pair = particle+particle;

    TFile *fTpcPID = new TFile("tpc_"+cutComb+"rootFiles/results_"+pair+".root", "UPDATE");
    fTpcPID->SetCompressionAlgorithm(1);
    fTpcPID->Close();

    TFile *fTofPID = new TFile("tofPidEff_"+cutComb+"rootFiles/results_"+pair+".root", "READ");
    fTofPID->SetCompressionAlgorithm(1);
    fTofPID->Close();

    TFile *fTofMatch = new TFile(cutComb+"rootFiles/results_"+pair+".root", "READ");
    fTofMatch->SetCompressionAlgorithm(1);
    fTofMatch->Close();
}