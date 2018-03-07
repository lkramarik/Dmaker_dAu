void project_ntp(TString input) {
    TFile* data = new TFile(input ,"r");
    TNtuple* ntpS = (TNtuple*)data -> Get("ntp_signal");
    TNtuple* ntpB = (TNtuple*)data -> Get("ntp_background");

    //    TCut* background = new TCut("flag==4||flag==5");
//    TCut* dcaDaughter = new TCut("dcaMax<0.0100");
//    TCut* k_dca = new TCut("k_dca<0.008");
//    TCut* pi1_dca = new TCut("pi1_dca<0.008");

    TFile *fOut = new TFile("res_"+input,"recreate");

// //     TH1D* h = new TH1D("h", "h", 2000, 0., 0.2);
//     TH1D* h2 = new TH1D("h2", "h2", 2000, 0., 0.2);
    
    TH1D* hS = new TH1D("hS", "hS", 2000, 0.4, 2.4);
    TH1D* hB = new TH1D("hB", "hB", 2000, 0.4, 2.4);

//    string line[9] = {"D_pt","D_mass", "pi1_dca", "k_dca", "D_decayL", "dcaMax", "pi1_nSigma", "k_nSigma"};

//    Double_t minim[] = {0,0.4,0,0,0,0,-3,-3};
//    Double_t maxim[] = {7, 2.4, 2, 2, 2, 2, 3, 3};
//    Double_t nbins[] = {700, 2000, 2000, 2000, 2000, 2000,600, 600};

//    h -> SetTitle(var);
//    h -> GetXaxis() -> Set(nbins[i], minim[i], maxim[i]);

// TMVA all pt
     ntpS -> Project("hS", "D_mass", "");
     ntpB -> Project("hB", "D_mass", "");


//     TMVA liang 1-2
//     ntpS -> Project("hS", "D_mass", "(D_pt>=1)&&(D_pt<2)&&(k_TOFinvbeta>0.)&&(k_TOFinvbeta<0.03)&&(pi1_TOFinvbeta>0.)&&(k_dca>0.0087)&&(pi1_dca>0.0099)&&((D_decayL*sqrt(1-cos(D_theta)*cos(D_theta)))<0.0075)&&(dcaDaughters<0.0093)&&(D_decayL>0.0232)");
//     ntpB -> Project("hB", "D_mass", "(D_pt>=1)&&(D_pt<2)&&(k_TOFinvbeta>0.)&&(k_TOFinvbeta<0.03)&&(pi1_TOFinvbeta>0.)&&(k_dca>0.0087)&&(pi1_dca>0.0099)&&((D_decayL*sqrt(1-cos(D_theta)*cos(D_theta)))<0.0075)&&(dcaDaughters<0.0093)&&(D_decayL>0.0232)");
    
    
// TMVA bit wider 1-2
//     ntpS -> Project("hS", "D_mass", "(D_pt>=1)&&(D_pt<2)&&(k_TOFinvbeta>0.)&&(k_TOFinvbeta<0.03)&&(pi1_TOFinvbeta>0.)&&(k_dca>0.007)&&(pi1_dca>0.007)&&((D_decayL*sqrt(1-cos(D_theta)*cos(D_theta)))<0.0095)&&(dcaDaughters<0.01)&&(D_decayL>0.01)");
//     ntpB -> Project("hB", "D_mass", "(D_pt>=1)&&(D_pt<2)&&(k_TOFinvbeta>0.)&&(k_TOFinvbeta<0.03)&&(pi1_TOFinvbeta>0.)&&(k_dca>0.007)&&(pi1_dca>0.007)&&((D_decayL*sqrt(1-cos(D_theta)*cos(D_theta)))<0.0095)&&(dcaDaughters<0.01)&&(D_decayL>0.01)");
    
// TMVA bit wider 2-3
//    ntpS -> Project("hS", "D_mass", "(D_pt>=2)&&(D_pt<3)&&(k_TOFinvbeta>0.)&&(k_TOFinvbeta<0.03)&&(pi1_TOFinvbeta>0.)&&(k_dca>0.007)&&(pi1_dca>0.007)&&((D_decayL*sqrt(1-cos(D_theta)*cos(D_theta)))<0.0095)&&(dcaDaughters<0.01)&&(D_decayL>0.01)&&(cos(D_theta)>0.5)");
//    ntpB -> Project("hB", "D_mass", "(D_pt>=2)&&(D_pt<3)&&(k_TOFinvbeta>0.)&&(k_TOFinvbeta<0.03)&&(pi1_TOFinvbeta>0.)&&(k_dca>0.007)&&(pi1_dca>0.007)&&((D_decayL*sqrt(1-cos(D_theta)*cos(D_theta)))<0.0095)&&(dcaDaughters<0.01)&&(D_decayL>0.01)");
    

// TMVA bit wider 3-5
//     ntpS -> Project("hS", "D_mass", "(D_pt>=3)&&(D_pt<5)&&(k_TOFinvbeta>0.)&&(k_TOFinvbeta<0.03)&&(pi1_TOFinvbeta>0.)&&(k_dca>0.007)&&(pi1_dca>0.007)&&((D_decayL*sqrt(1-cos(D_theta)*cos(D_theta)))<0.0095)&&(dcaDaughters<0.01)&&(D_decayL>0.01)");
//     ntpB -> Project("hB", "D_mass", "(D_pt>=3)&&(D_pt<5)&&(k_TOFinvbeta>0.)&&(k_TOFinvbeta<0.03)&&(pi1_TOFinvbeta>0.)&&(k_dca>0.007)&&(pi1_dca>0.007)&&((D_decayL*sqrt(1-cos(D_theta)*cos(D_theta)))<0.0095)&&(dcaDaughters<0.01)&&(D_decayL>0.01)");

    
    
    hS -> Rebin(8);
    hB -> Rebin(8);
    
    
    hS -> Write("Dmass_sig");
    hB -> Write("Dmass_back");
    
    hS -> Draw();
    hB -> Draw("same");
//    h -> Draw();

//    data -> Close();
//    fOut -> Close();


}
