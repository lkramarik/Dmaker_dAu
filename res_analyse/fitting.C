#include "TStyle.h"


void fitting() {
    Float_t fitRMin = 1.65;
    Float_t fitRMax = 2.35;

    TString intLowLow = "1.7";
    TString intLowUp = "1.8";
    TString intUpLow = "1.92";
    TString intUpUp = "2.2";

    int rebin = 5;
//     int rebin = 10;
//    TString folder = "res_ntp/p17ib/";
//    TString input = "res_ntp_dau_p17ib_hpss_all_dca60.root";


    TString folder = "";
//    TString folder = "res_ntp/p17id/compareToTMVA/01/";
//    TString input = "res_ntp_dau_p17id_cos095_dca30.root";
//    TString input = "res_ntp_dau_p17id_cos095_tof003.root";
//    TString input = "res_ntp_dau_p17id_dca30.root";
// TString input = "res_ntp_dau_cos095_dca30_decay60.root";
//    TString input = "res_ntp_dau_p17id_wide.root";
    TString input = "res_MLPoutput.root";

//    TString input = "res_ntp_dau_p17id_costheta095.root";
    cout<<input<<endl;
    TFile* data = new TFile(folder + input ,"r");

//     hInvMassBackMin = (TH1F*)data -> Get("background minus");
//     hInvMassBackPlus = (TH1F*)data -> Get("background plus");
//     hInvMassBackMin -> Sumw2();
//     hInvMassBackPlus -> Sumw2();
    //background from geometric average
    hInvMassBack = (TH1F*)data -> Get("background");

    //background from sum of background minus and plus combos
//    TH1F* hInvMassBack = new TH1F("background", "background", 2000, 0.4, 2.4);
//    hInvMassBack -> Add(hInvMassBackPlus,1);
//    hInvMassBack -> Add(hInvMassBackMin,1);
    hInvMassBack -> Sumw2();

    hInvMassSign = (TH1F*)data -> Get("signal");
    hInvMassSign -> Sumw2();
    hInvMassSign -> Rebin(rebin);
    hInvMassBack -> Rebin(rebin);

    hInvMassSign -> SetMarkerColor(46);
    hInvMassSign -> SetLineColor(46);

    Double_t hBackIntegral = hInvMassBack -> Integral(hInvMassBack -> FindBin(intLowLow.Atof()), hInvMassBack -> FindBin(intLowUp.Atof()), "") + hInvMassBack -> Integral(hInvMassBack -> FindBin(intUpLow.Atof()), hInvMassBack -> FindBin(intUpUp.Atof()), "");
    Double_t hSignIntegral = hInvMassSign -> Integral(hInvMassSign -> FindBin(intLowLow.Atof()), hInvMassSign -> FindBin(intLowUp.Atof()), "") + hInvMassSign -> Integral(hInvMassSign -> FindBin(intUpLow.Atof()), hInvMassSign -> FindBin(intUpUp.Atof()), "");

    cout<<hInvMassBack -> GetNbinsX()<<endl;
    cout<<hInvMassSign -> GetNbinsX() << endl;

    cout<<"Backgroung Integral: "<<hBackIntegral<<endl;
    cout<<"Signal Integral :"<<hSignIntegral<<endl;

    hInvMassBack -> Scale(hSignIntegral/hBackIntegral);

    const int Nbins = 2000/rebin;
    float err[Nbins], errS[Nbins];
    for (int j=0; j<Nbins; j++) {
        err[j] = hInvMassBack -> GetBinError(j);
        errS[j] = hInvMassSign -> GetBinError(j);
    }

//    for (j=0; j<Nbins; j++) {
//        err[j] = err[j]*hSignIntegral/hBackIntegral;
//        hInvMassBack -> SetBinError(j, err[j]);
//    }   //scaling error

//    hInvMassSign->Clone("signal_orig");
    TH1F *hInvMassSignOrig = (TH1F*)hInvMassSign->Clone("signal_orig");
    Double_t nentriesSig = hInvMassSignOrig->Integral(hInvMassSignOrig->FindBin(1.7), hInvMassSignOrig->FindBin(2),"");

    hInvMassSign -> Add(hInvMassBack,-1);
    for (j=0; j<Nbins; j++) {
        hInvMassSign -> SetBinError(j, sqrt(err[j]*err[j] + errS[j]*errS[j]));
    }

    TCanvas *c3 = new TCanvas("c3","c3",1200,900);
    gStyle->SetOptFit(1);
    gStyle->SetStatY(0.899);
    gStyle->SetStatX(0.9);

    TF1 *fun0 = new TF1("fun0","pol1(0)+gaus(2)",fitRMin,fitRMax);
    fun0->SetParameters(1.,1.,1.,1.84,0.01);
    fun0->SetLineColor(2);
    fun0->SetLineStyle(7);
    fun0->SetParName(2,"height");
    fun0->SetParName(3,"mean");
    fun0->SetParName(4,"sigma");

    const float MKSize    = 1.;
    const float rotwthmin = 1.80; // peak mean fitting range
    const float rotwthmax = 1.98; //peak mean fitting range
    double mm[Nbins],ym[Nbins],yme[Nbins],ym1[Nbins];

    for (int ib=0; ib<Nbins; ib++) {
        mm[ib]  = hInvMassSign -> GetBinCenter(ib);
        ym[ib]  = hInvMassSign -> GetBinContent(ib);
        yme[ib] = hInvMassSign -> GetBinError(ib);
    }

    TGraphErrors *gm = new TGraphErrors(Nbins,mm,ym,0,yme);
    gm->SetMarkerStyle(20);
    gm->SetMarkerSize(0.9);
    gm->SetMarkerColor(9);
    gm->SetLineColor(9);
    fun0->SetParLimits(3,rotwthmin,rotwthmax);
    fun0 -> SetLineColor(9);
    gm->Fit(fun0,"OR");

    Double_t mean = fun0->GetParameter(3);
    Double_t sigma = fun0->GetParameter(4);
    cout<<"mean: "<<mean<<endl;
    cout<<"sigma: "<<sigma<<endl;
    Double_t nsigma = 3;


    TF1 *resfunm = new TF1("resfunm","pol1",fitRMin,fitRMax);
    resfunm->SetParameters(fun0->GetParameter(0),fun0->GetParameter(1));
    resfunm->SetLineStyle(7);
    resfunm->SetLineWidth(1);
    //resfunm->Draw("same");
    //fun0->Draw("same");

    TF1 *resfunm1 = new TF1("resfunm1","pol1(0)+gaus(2)",fitRMin,fitRMax);
    resfunm1->SetParameters(0.,0.,fun0->GetParameter(2),fun0->GetParameter(3),fun0->GetParameter(4));
    resfunm1->SetLineColor(4);
    resfunm1->SetLineStyle(7);
    //resfunm1->Draw("same");

    Double_t intLow = hInvMassSignOrig -> FindBin(mean-nsigma*sigma);
    Double_t intUp = hInvMassSignOrig -> FindBin(mean+nsigma*sigma);
    cout<<"Integral: "<<intLow<<" "<<intUp<<endl;
    cout<<"Integral: "<<hInvMassSignOrig->GetBinCenter(intLow)<<" "<<hInvMassSignOrig->GetBinCenter(intUp)<<endl;
    cout<<"Integral: "<<intLow<<" "<<intUp<<endl;
    cout<<"Integral: "<<hInvMassBack->GetBinCenter(intLow)<<" "<<hInvMassBack->GetBinCenter(intUp)<<endl;
    Double_t integral_function_yield = resfunm1->Integral(mean-nsigma*sigma, mean+nsigma*sigma);
//    Double_t integral_function_yield = fun0->Integral(mean-nsigma*sigma, mean+nsigma*sigma);
    cout<<"Integral raw yield: "<<integral_function_yield<<endl;
//    Double_t integral_function_yield_error = resfunm1->IntegralError(mean-nsigma*sigma, mean+nsigma*sigma);
    Double_t integral_function_yield_error = fun0->IntegralError(mean-nsigma*sigma, mean+nsigma*sigma);
    cout<<"Integral raw yield error: "<<integral_function_yield_error<<endl;





    gm->GetYaxis()->SetTitle("Raw Counts");
    gm->GetYaxis()->SetTitleOffset(1.5);
    gm->GetXaxis()->SetTitle("Mass_{K#pi} (GeV/c^{2})");
    gm->SetTitle("");
    gm->PaintStats(fun0);
//    gm->GetXaxis()->SetRangeUser(0.4,2.4);
    gm->GetXaxis()->SetRangeUser(1.72,2.0);
    gm->GetYaxis()->SetRangeUser(-1000.,2000.0);
    gm->Draw("ap");

    for(j=0; j<Nbins; j++) ym1[j] = ym[j] - resfunm->Eval(mm[j]);
    TGraphErrors *gm1 = new TGraphErrors(Nbins,mm,ym1,0,yme);
    gm1->SetMarkerStyle(24);
    gm1->SetMarkerSize(MKSize-0.1);
    gm1->SetMarkerColor(2);
    gm1->SetLineColor(2);
//    gm1->Draw("same,AP");

//    TMultiGraph *mg = new TMultiGraph();
//    mg->Add(gm);
//    mg->Add(gm1);
//    mg->Draw("ap");
//    //mg->GetYaxis()->SetTitle("Raw Counts (#times 10^{3})");
//    mg->GetYaxis()->SetTitle("Raw Counts");
//    //mg->GetYaxis()->SetTitleSize(titsize);
//    mg->GetYaxis()->SetTitleOffset(1.8);
//    //mg->GetYaxis()->SetLabelOffset(0.03);
//    //mg->GetYaxis()->SetLabelSize(0.047);
//    //mg->GetXaxis()->SetNdivisions(208);
//    //  mg->GetXaxis()->CenterTitle();
//    mg->GetXaxis()->SetTitle("Mass_{K#pi} (GeV/c^{2})");


//    c3->Update();

    resfunm->SetLineColor(46);
    resfunm1->SetLineColor(9);

//    resfunm->Draw("same");
//    resfunm1->Draw("same");
    Double_t S = hInvMassSign -> Integral(hInvMassSign->FindBin(mean-nsigma*sigma), hInvMassSign->FindBin(mean+nsigma*sigma) ,"");
    Double_t xtest, ytest;
    cout<<gm1->GetPoint(hInvMassSignOrig->FindBin(mean-nsigma*sigma), xtest, ytest)<<endl;
    cout<<xtest<<endl;
    cout<<ytest<<endl;
    cout<<gm1->GetPoint(hInvMassSignOrig->FindBin(mean+nsigma*sigma), xtest, ytest)<<endl;
    cout<<xtest<<endl;
    cout<<ytest<<endl;
//    Double_t S = gm1 -> Integral(hInvMassSignOrig->FindBin(mean-nsigma*sigma), hInvMassSignOrig->FindBin(mean+nsigma*sigma) );
    Double_t B = hInvMassBack -> Integral(hInvMassBack->FindBin(mean-nsigma*sigma), hInvMassBack->FindBin(mean+nsigma*sigma) ,"");
    Double_t binSize = (2.4-0.4)*rebin/2000;
//    Double_t B  = resfunm->Integral((mean-nsigma*sigma), (mean+nsigma*sigma))/binSize;
    TLine *leftline = new TLine(mean-nsigma*sigma,-1000,mean-nsigma*sigma,2000);
    leftline->SetLineStyle(3);
    leftline->SetLineColor(46);
    TLine *rightline = new TLine(mean+nsigma*sigma,-1000,mean+nsigma*sigma,2000);
    rightline->SetLineStyle(3);
    rightline->SetLineColor(46);

    rightline->Draw("");
    leftline->Draw("same");

    cout<<"S: "<<S<<endl;
    cout<<"B: "<<B<<endl;
    Double_t significance_fit = integral_function_yield/integral_function_yield_error ;
    Double_t significance_bins = S/TMath::Sqrt(S+B);
    cout<<"Significance from fit: "<<abs(significance_fit)<<endl;
    cout<<"Significance from number of entries in bins: "<<significance_bins<<endl;

    TPaveText *text = new TPaveText(0.172,0.735,0.443,0.8300,"brNDC");
    text -> SetTextSize(0.025);
    text->SetTextColor(39);
    text -> SetLineColor(0);
    text -> SetShadowColor(0);
    text -> SetFillColor(0);
//    text -> AddText(input);

    text -> AddText(Form("N_entries in sig. 1.7-2.0 GeV/c^{2}: %g", nentriesSig));
    text -> AddText(Form("Bin size: %g GeV/c^{2}", binSize));
    TString paveSc = "Scaling integral: " + intLowLow + "-" + intLowUp +", " + intUpLow + "-" + intUpUp + " GeV/c^{2}";
    text -> AddText(paveSc);
//    text -> AddText(Form("Significance from fit: %0.3g", significance_fit));
//    text -> AddText(Form("Significance: %0.3g", significance_bins));

    TPaveText *text1 = new TPaveText(0.66,0.91,0.9391,0.967,"brNDC");
    text1 -> SetTextSize(0.04);
    text1 -> SetShadowColor(0);
    text1 -> SetLineColor(0);
    text1 -> SetFillColor(0);
    text1 -> AddText("d+Au 200 GeV");

    TPaveText *text3 = new TPaveText(0.0,0.01,0.27,0.0350,"brNDC");
    text3 -> SetTextSize(0.02);
    text3->SetTextColor(39);
    text3 -> SetLineColor(0);
    text3 -> SetShadowColor(0);
    text3 -> SetFillColor(0);
    text3 -> AddText(input);

    TPaveText *text4 = new TPaveText(0.19,0.848,0.229,0.873,"brNDC");
    text4 -> SetTextSize(0.03);
//    text4->SetTextColor(39);
    text4 -> SetLineColor(0);
    text4 -> SetShadowColor(0);
    text4 -> SetFillColor(0);
    text4 -> AddText(Form("Significance: %0.3g", significance_bins));
    TLatex tx2;
    tx2.SetNDC();
    tx2.SetTextSize(0.04);
    tx2.DrawLatex(0.1,0.93,Form("p_{T}: %3.1f-%3.1f GeV/c", 1, 2));


    text -> Draw("same");
    text1 -> Draw("same");
    text3 -> Draw("same");
    text4 -> Draw("same");
//         gStyle->SetOptFit();

    TCanvas *c4 = new TCanvas("c4","c4",1200,1000);
    TLegend *legend = new TLegend(0.40,0.845, 0.560, 0.897,"","brNDC");
    hInvMassBack -> SetTitle("");
    hInvMassBack -> SetStats(0);
    hInvMassBack -> SetLineColor(46);
    hInvMassBack -> SetMarkerColor(46);
    hInvMassBack -> SetMarkerStyle(2);

    gPad->SetLeftMargin(0.15);

    hInvMassSign -> SetTitle("");
    hInvMassSign -> SetStats(0);
    hInvMassSign -> SetLineColor(9);
    hInvMassSign ->SetMarkerStyle(20);
    hInvMassSign -> SetMarkerColor(9);
    hInvMassSign->GetYaxis()->SetTitle("Raw Counts");
    hInvMassSign->GetYaxis()->SetTitleOffset(1.7);
    hInvMassSign->GetXaxis()->SetTitle("Mass K#pi (GeV/c^{2})");
    hInvMassSign->GetXaxis()->SetLabelFont(42);
    hInvMassSign->GetXaxis()->SetTitleFont(42);
    hInvMassSign->GetYaxis()->SetLabelFont(42);
    hInvMassSign->GetYaxis()->SetTitleFont(42);
    hInvMassSign -> GetXaxis()->SetRangeUser(0.65,2.4);
    hInvMassSign -> Draw("");

//    TLine *leftline1 = new TLine(1.7,hInvMassSign->GetYaxis()->GetXmin(),1.7,2000);
    TLine *leftline1 = new TLine(1.72,hInvMassSign->GetMinimum(),1.72,hInvMassSign->GetMaximum());
    leftline1->SetLineStyle(9);
    leftline1->SetLineColor(28);
    leftline1->Draw("same");

    TLine *leftline2 = new TLine(2,hInvMassSign->GetMinimum(),2,hInvMassSign->GetMaximum());
    leftline2->SetLineStyle(9);
    leftline2->SetLineColor(28);
    leftline2->Draw("same");

//    legend -> AddEntry(hInvMassBack, Form("scaled background (%g-%g & %g-%g GeV/c^{2})", intLowLow.Atof(), intLowUp.Atof(), intUpLow.Atof(), intUpUp.Atof()), "pl");
    legend -> AddEntry(hInvMassSign, "Unlike-sign - like-sign pairs", "pl");
    TLatex tx1;
    tx1.SetNDC();
    tx1.SetTextSize(0.04);
    tx1.DrawLatex(0.1,0.93,Form("p_{T}: %3.1f-%3.1f GeV/c", 0.2, 6.0));
    legend -> SetFillStyle(0);
    legend -> SetLineColor(0);
    legend -> SetTextSize(0.035);
    legend -> Draw("same");
    text1 -> Draw("same");

    //SIGNAL AND BACKGROUND
    TCanvas *c5 = new TCanvas("c5","c5",1200,1000);
    TLegend *legend3 = new TLegend(0.292,0.807, 0.452, 0.877,"","brNDC");
    hInvMassBack -> SetTitle("");
    hInvMassBack -> SetStats(0);
    hInvMassBack -> SetLineColor(46);
    hInvMassBack -> SetMarkerColor(46);
    gPad->SetLeftMargin(0.15);

    hInvMassSignOrig -> SetTitle("");
    hInvMassSignOrig -> SetStats(0);
    hInvMassSignOrig -> SetLineColor(9);
    hInvMassSignOrig ->SetMarkerStyle(20);
    hInvMassSignOrig -> SetMarkerColor(9);
    hInvMassSignOrig->GetYaxis()->SetTitle("Raw Counts");
    hInvMassSignOrig->GetYaxis()->SetTitleOffset(1.7);
    hInvMassSignOrig->GetXaxis()->SetTitle("Mass K#pi (GeV/c^{2})");
    hInvMassSignOrig->GetXaxis()->SetLabelFont(42);
    hInvMassSignOrig->GetXaxis()->SetTitleFont(42);
    hInvMassSignOrig->GetYaxis()->SetLabelFont(42);
    hInvMassSignOrig->GetYaxis()->SetTitleFont(42);
    hInvMassSignOrig -> GetXaxis()->SetRangeUser(1.72,2.);
    hInvMassSignOrig -> Draw("");
    hInvMassBack->Draw("same");

    legend3 -> AddEntry(hInvMassBack, Form("Scaled background (%g-%g & %g-%g GeV/c^{2})", intLowLow.Atof(), intLowUp.Atof(), intUpLow.Atof(), intUpUp.Atof()), "pl");
    legend3 -> AddEntry(hInvMassSign, "Unlike-sign (signal) pairs", "pl");
    TLatex tx2;
    tx2.SetNDC();
    tx2.SetTextSize(0.04);
    tx2.DrawLatex(0.1,0.93,Form("p_{T}: %3.1f-%3.1f GeV/c", 0.2, 6.0));
    legend3 -> SetFillStyle(0);
    legend3 -> SetLineColor(0);
    legend3 -> SetTextSize(0.035);
    legend3 -> Draw("same");
    text1 -> Draw("same");

    c3 -> SaveAs("signal_"+input+".png");
    c4 -> SaveAs("mass_"+input+".png");
    c5 -> SaveAs("mass_zoom_"+input+".png");


//     hInvMassBack -> Draw();
//    hInvMassSignOrig->Draw();
}

void plotDCA() {
    TString folder = "";
//    TString folder = "res_ntp/p17id/";
    TString input = "res_plots_signal.root";
    TFile* data = new TFile(folder + input ,"r");

    TCanvas *c3 = new TCanvas("c3","c3",1000,1200);

//    hEvents = (TH1F*)data -> Get("hcosTheta");
    hEvents = (TH1F*)data -> Get("hk_dca");
    hEvents1 = (TH1F*)data -> Get("hpi1_dca");
    hEvents ->GetYaxis()->SetTitleOffset(1.5);
//    hInvMassBackPlus = (TH1F*)data -> Get("background plus");
    hEvents -> SetTitle("");
    hEvents -> SetStats(0);
//    hEvents -> SetFillColor(18);
//    hEvents1 -> SetFillColor(18);
    hEvents->SetLineWidth(3); //2
    hEvents1->SetLineWidth(3); //2

    hEvents -> SetLineColor(46);
    hEvents1 -> SetLineColor(9);
    hEvents ->GetXaxis()->SetRangeUser(0.,0.04);
//    hEvents ->GetXaxis()-> SetTitle("cos(pointing angle diff.) [-]");
    hEvents ->GetXaxis()-> SetTitle("DCA daughters [cm]");
    hEvents ->GetYaxis()-> SetTitle("Number of pairs");
//    hEvents -> SetFillStyle(3344);

//    gPad->SetLogy();
    gPad->SetLeftMargin(0.15);
    TLegend *legend = new TLegend(0.397,0.82, 0.55, 0.87,"","brNDC");
    legend -> SetFillStyle(0);
    legend -> SetLineColor(0);
    legend -> SetTextSize(0.035);
    legend -> AddEntry(hEvents, "kaons", "pl");
    legend -> AddEntry(hEvents1, "pions", "pl");

    TLine *leftline2 = new TLine(0.003,hEvents->GetMinimum(),0.003,hEvents->GetMaximum()+0.05*hEvents->GetMaximum());
    leftline2->SetLineWidth(3);
    leftline2->SetLineStyle(9);
    leftline2->SetLineColor(28);

//    TPaveText *text1 = new TPaveText(0.66,0.91,0.9391,0.967,"brNDC");
//    text1 -> SetTextSize(0.04);
//    text1 -> SetShadowColor(0);
//    text1 -> SetLineColor(0);
//    text1 -> SetFillColor(0);
//    text1 -> AddText("d+Au 200 GeV");

    TPaveText *text2 = new TPaveText(0.66,0.91,0.9391,0.967,"brNDC");
    text2 -> SetTextSize(0.03);
    text2 -> SetShadowColor(0);
    text2 -> SetLineColor(0);
    text2 -> SetFillColor(0);
    text2->SetTextColor(28);
    text2 -> AddText("cuts I., II.");

//    TLine *leftline2 = new TLine(0.003,hEvents->GetMinimum(),0.003,hEvents->GetMaximum()+0.05*hEvents->GetMaximum());
//    leftline2->SetLineWidth(2);
//    leftline2->SetLineStyle(9);
//    leftline2->SetLineColor(8);

    hEvents -> Draw();
    hEvents1 -> Draw("same");
    legend -> Draw("same");
    leftline2->Draw("same");
//    text1->Draw("same");
    text2->Draw("same");
}

void plot() {
    TString folder = "";
//    TString folder = "res_ntp/p17id/";
    TString input = "res_plots_signal.root";
    TFile* data = new TFile(folder + input ,"r");

    TCanvas *c3 = new TCanvas("c3","c3",1000,1200);

    hEvents = (TH1F*)data -> Get("hcosTheta");
//    hEvents = (TH1F*)data -> Get("hk_dca");
//    hEvents1 = (TH1F*)data -> Get("hpi1_dca");
    hEvents ->GetYaxis()->SetTitleOffset(1.9);
//    hInvMassBackPlus = (TH1F*)data -> Get("background plus");
    hEvents -> SetTitle("");
    hEvents -> SetStats(0);
    hEvents -> SetFillColor(18);
    hEvents->SetLineWidth(3); //2

    hEvents -> SetLineColor(46);
    hEvents ->GetXaxis()->SetRangeUser(0.75,1);
    hEvents ->GetXaxis()-> SetTitle("cos(pointing angle diff.) [-]");
//    hEvents ->GetXaxis()-> SetTitle("DCA daughters [cm]");
    hEvents ->GetYaxis()-> SetTitle("Number of pairs");
//    hEvents -> SetFillStyle(3344);

//    gPad->SetLogy();
    gPad->SetLeftMargin(0.15);
    TLegend *legend = new TLegend(0.397,0.82, 0.55, 0.87,"","brNDC");
    legend -> SetFillStyle(0);
    legend -> SetLineColor(0);
    legend -> SetTextSize(0.035);
//    legend -> AddEntry(hEvents, "kaons", "pl");
//    legend -> AddEntry(hEvents1, "pions", "pl");

    TLine *leftline2 = new TLine(0.8,hEvents->GetMinimum()-0.05*hEvents->GetMinimum(),0.8,hEvents->GetMaximum()+0.04*hEvents->GetMaximum());
    leftline2->SetLineWidth(3);
    leftline2->SetLineStyle(9);
    leftline2->SetLineColor(28);

    TLine *leftline3 = new TLine(0.95,hEvents->GetMinimum()-0.05*hEvents->GetMinimum(),0.95,hEvents->GetMaximum()+0.04*hEvents->GetMaximum());
    leftline3->SetLineWidth(3);
    leftline3->SetLineStyle(9);
    leftline3->SetLineColor(9);

//    TPaveText *text1 = new TPaveText(0.66,0.91,0.9391,0.967,"brNDC");
//    text1 -> SetTextSize(0.04);
//    text1 -> SetShadowColor(0);
//    text1 -> SetLineColor(0);
//    text1 -> SetFillColor(0);
//    text1 -> AddText("d+Au 200 GeV");

    TPaveText *text2 = new TPaveText(0.66,0.91,0.9391,0.967,"brNDC");
    text2 -> SetTextSize(0.03);
    text2 -> SetShadowColor(0);
    text2 -> SetLineColor(0);
    text2 -> SetFillColor(0);
    text2->SetTextColor(28);
    text2 -> AddText("cuts I.");

    TPaveText *text3 = new TPaveText(0.66,0.91,0.9391,0.967,"brNDC");
    text3 -> SetTextSize(0.03);
    text3 -> SetShadowColor(0);
    text3 -> SetLineColor(0);
    text3 -> SetFillColor(0);
    text3->SetTextColor(9);
    text3 -> AddText("cuts II.");
//    TLine *leftline2 = new TLine(0.003,hEvents->GetMinimum(),0.003,hEvents->GetMaximum()+0.05*hEvents->GetMaximum());
//    leftline2->SetLineWidth(2);
//    leftline2->SetLineStyle(9);
//    leftline2->SetLineColor(8);

    hEvents -> Draw();
    legend -> Draw("same");
    leftline2->Draw("same");
    leftline3->Draw("same");
//    text1->Draw("same");
    text2->Draw("same");
    text3->Draw("same");
}



void plotEv() {

    TString folder = "res_ntp/p17id/";
    TString input = "res_ntp_dau_p17id_wide.root";

    TFile* data = new TFile(folder + input ,"r");

    TCanvas *c3 = new TCanvas("c3","c3",1000,1200);

    hEvents = (TH1F*)data -> Get("hEventStat1");
//    hEvents ->GetYaxis()->SetTitleOffset(1.5);
//    hInvMassBackPlus = (TH1F*)data -> Get("background plus");
    hEvents -> SetTitle("");
    hEvents -> SetStats(0);
    hEvents -> SetFillColor(46);
    hEvents->SetLineWidth(2);
    hEvents -> SetLineColor(46);
    hEvents ->GetYaxis()-> SetTitleOffset(1.4);
//    hEvents ->GetXaxis()-> SetTitle("DCA daughters [cm]");
//    hEvents ->GetYaxis()-> SetTitle("Number of pairs");
//    hEvents -> SetFillStyle(3344);

    gPad->SetLeftMargin(0.15);
    hEvents -> Draw();

}

