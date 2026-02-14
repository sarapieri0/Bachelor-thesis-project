#include <TCanvas.h>
#include <TFile.h>
#include <TF1.h>
#include <TLegend.h>

void fit_vbs(){

    TFile *histsTT = TFile::Open("output_diffMassAK4_AK8_loose_ssWW_TT.root", "READ");
    TFile *histsLL = TFile::Open("output_diffMassAK4_AK8_loose_ssWW_LL.root", "READ");
    TFile *hists = TFile::Open("output_diffMassAK4_AK8_loose_ssWW.root", "READ");

    if (!hists || hists->IsZombie()) {
        std::cout<<"Errore: file non trovato\n"<<std::endl;
    }


    TH1F *h_costheta_starLL = nullptr;
    histsLL->GetObject("h_costheta_star", h_costheta_starLL);
    if (!h_costheta_starLL) {
        std::cout<<"Errore: istogramma non trovato\n"<<std::endl;
    }
    h_costheta_starLL->Scale(1.0 / h_costheta_starLL->Integral());


    TH1F *h_costheta_starTT = nullptr;
    histsTT->GetObject("h_costheta_star", h_costheta_starTT);
    if (!h_costheta_starTT) {
        std::cout<<"Errore: istogramma non trovato\n"<<std::endl;
    }
    h_costheta_starTT->Scale(1.0 / h_costheta_starTT->Integral());


    TH1F *h_costheta_star = (TH1F*)h_costheta_starLL->Clone("h_costheta_star");
    h_costheta_star->Scale(0.70);
    h_costheta_star->Add(h_costheta_starTT, 0.30);

    

    TF1 *f1 = new TF1("f1", "[0]*(1-x*x)", -1, 1);
    TF1 *f2 = new TF1("f2", "[0]*x*x + [1]*x + [2]", -1, 1);
    TF1 *f3 = new TF1("f3", "[0]*(-0.0222054*x*x + 0.0255167*x + 0.0564097) + [1]*(-0.0761086*x*x + 0.00190588*x + 0.0744102)", -1, 1);
    f3->SetParameter(0, 0.5);
    f3->SetParameter(1, 0.5);
    f3->SetParLimits(0, 0.0, 1.0);
    f3->SetParLimits(1, 0.0, 1.0);

    h_costheta_starTT->Fit(f2, "R");
    //double fl = f->GetParameter(0);
    //double sigma = f->GetParError(0);
    double chi2  = f2->GetChisquare();


    TCanvas *c_costheta_fit = new TCanvas("c_costheta_fit", "Fit h_costheta_star", 800, 600);
    h_costheta_starTT->Draw("hist");
    f2->Draw("same");
    c_costheta_fit->SaveAs("plots_lp/c_costheta_fit.pdf");


    gROOT->ProcessLine(".q");
}
