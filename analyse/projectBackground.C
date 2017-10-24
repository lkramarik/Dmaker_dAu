void projectBackground(TString input) {
    TFile* data = new TFile(input ,"r");
    TNtuple* ntp = (TNtuple*)data -> Get("ntp");

    TCut* signal = new TCut("flag==1||flag==0");
    TCut* background = new TCut("flag==4||flag==5");
    TCut* dcaDaughter = new TCut("dcaMax<0.0100");
    TCut* k_dca = new TCut("k_dca<0.008");
    TCut* pi1_dca = new TCut("pi1_dca<0.008");

    TList *listOut = new TList();
    TFile *fOut = new TFile("background_properties_"+input,"recreate");

    TH1F* h = new TH1F("h", "h", 10000, 0., 10.);

    string line[8] = {"D_pt","D_mass", "pi1_dca", "k_dca", "D_decayL", "dcaMax", "pi1_nSigma", "k_nSigma"};
    Double_t minim[] = {0,0.4,0,0,0,0,-3,-3};
    Double_t maxim[] = {7, 2.4, 2, 2, 2, 2, 3, 3};
    Double_t nbins[] = {700, 2000, 2000, 2000, 2000, 2000,600, 600};
    for (int i = 0; i < 8; i++){
        std::cout << line[i] << std::endl;
        const char * var = line[i].c_str();
        h -> SetTitle(var);
        h -> GetXaxis() -> Set(nbins[i], minim[i], maxim[i]);
        ntp -> Project("h", var, *background);
        h -> Write(var);
    }
  data -> Close();
  fOut -> Close();

}
