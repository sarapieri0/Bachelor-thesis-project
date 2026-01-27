#include "samples_nonPol.h"
#include <TTree.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <iostream>
#include <vector>
#include <TROOT.h>
#include <cmath>
#include <algorithm>





// lista dei branch nel tree : https://cms-xpog.docs.cern.ch/autoDoc/NanoAODv12/2022/2023/doc_DYJetsToLL_M-50_TuneCP5_13p6TeV-madgraphMLM-pythia8_Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2.html

void jet_plots3(
    TString output="jets3.root"
){
    TFile *f=new TFile(output, "recreate");

    Int_t maxEvents = 2000;
    std::string tag = "ssWW";
    

//pt    
    TH1F *h_FJ_pt = new TH1F(("h_FatJet_pt_" + tag).c_str(), "FatJet p_{T};p_{T} [GeV];Events", 50, 170, 1000);
    
    TH1F *h_GFJ_pt = new TH1F(("h_GenJetAK8_pt_" + tag).c_str(), "GenJetAK8 p_{T};p_{T} [GeV];Events", 50, 150, 1000);
    
    
//eta    
    TH1F *h_FJ_eta = new TH1F(("h_FatJet_eta_" + tag).c_str(), "FatJet #eta; #eta Pseudorapidity; Events", 50,-6, 6);
    
    TH1F *h_GFJ_eta = new TH1F(("h_GenJetAK8_eta_" + tag).c_str(), "GenJetAK8 #eta;#eta Pseudorapidity; Events", 50, -6, 6);

//phi
    TH1F *h_FJ_phi = new TH1F(("h_FatJet_phi_" + tag).c_str(), "FatJet #phi; #phi; Events", 40, -4,4);
    
    TH1F *h_GFJ_phi = new TH1F(("h_GenJetAK8_phi_" + tag).c_str(), "GenJetAK8 #phi; #phi; Events", 40, -4, 4);


//mass
    TH1F *h_FJ_mass = new TH1F(("h_FJ_mass_" + tag).c_str(), "FatJet mass; mass [GeV]; Events", 50,0,250);
    
    TH1F *h_GFJ_mass = new TH1F(("h_GFJ_mass_" + tag).c_str(), "GenAK8 FatJet mass;mass [GeV];Events", 50, 0, 250);

//mass softdrop
    TH1F *h_FJ_msd = new TH1F(("h_FJ_msd_" + tag).c_str(), "FatJet SoftDrop mass;mass [GeV]", 50, 0, 250);


//tau
    TH1F *h_tau1 =new TH1F(("h_tau1_" + tag).c_str(), "FatJet #tau_{1}; #tau_{1}; Events", 30, 0, 0.6);
    TH1F *h_tau2 =new TH1F(("h_tau2_" + tag).c_str(), "FatJet #tau_{2}; #tau_{2}; Events", 30, 0, 0.6);
    TH1F *h_tau3= new TH1F(("h_tau3_" + tag).c_str(), "FatJet #tau_{3}; #tau_{3}; Events", 30, 0, 0.6);

    TH1F *h_tau21m= new TH1F(("h_tau21m_" + tag).c_str(), "#tau_{21} (60 < m_{sd} < 100); #tau_{21}; Events", 30, 0, 1);
    TH1F *h_tau21= new TH1F(("h_tau21_" + tag).c_str(), "#tau_{21}; #tau_{21}; Events", 30, 0, 1);
    TH1F *h_tau32m= new TH1F(("h_tau32m_" + tag).c_str(), "#tau_{32} (60 < m_{sd} < 100); #tau_{32}; Events", 30, 0, 1);
    TH1F *h_tau32= new TH1F(("h_tau32_" + tag).c_str(), "#tau_{32}; #tau_{32}; Events", 30, 0, 1);

    TH2F *h_tau21_msd = new TH2F(("h_msd_tau21_" + tag).c_str(), "SoftDrop mass   vs   #tau_{21}; #tau_{21}; SoftDrop", 50, 0, 1,    30, 20, 300);

    TH2F *h_tau32_msd = new TH2F(("h_msd_tau32_" + tag).c_str(), "SoftDrop mass   vs   #tau_{32}; #tau_{32}; SoftDrop", 50, 0, 1,    30, 20, 300);


//subjet
    TH1F *h_SJ_etad_sd = new TH1F(("h_SJ_etad_sd_" + tag).c_str(), "SubJet #eta difference; #eta; Events", 20,-1,1);
    TH1F *h_SJ_etad = new TH1F(("h_SJ_etad_" + tag).c_str(), "SubJet #eta difference; #eta; Events", 30,-1,1);
    
    TH1F *h_SJ_phid_sd = new TH1F(("h_SJ_phid_sd_" + tag).c_str(), "SubJet #phi difference; #phi; Events", 20,-1,1);
    TH1F *h_SJ_phid = new TH1F(("h_SJ_phid_" + tag).c_str(), "SubJet #phi difference; #phi; Events", 30,-1,1);

    //TH1F *h_SJ_ptd_sd = new TH1F(("h_SJ_ptd_sd_" + tag).c_str(), "SubJet p_{t} difference; p_{t}; Events", 30,150,1000);
    //TH1F *h_SJ_ptd = new TH1F(("h_SJ_ptd_" + tag).c_str(), "SubJet p_{t} difference; p_{t}; Events", 30,150,1000);

    TH1F *h_zg_sj = new TH1F(("h_zg_sj_" + tag).c_str(), "Subjets Z_{g}; Z_{g}; Events", 30, 0, 0.5);
    TH1F *h_zg_sjsd = new TH1F(("h_zg_sjsd_" + tag).c_str(), "Subjets Z_{g}; Z_{g}; Events", 30, 0, 0.5);
    TH1F *h_zg_pt200 = new TH1F(("h_zg_pt200_" + tag).c_str(), "Subjets Z_{g}, p_{t} > 200; Z_{g}; Events", 30, 0, 0.5);
    TH1F *h_zg_pt230 = new TH1F(("h_zg_pt230_" + tag).c_str(), "Subjets Z_{g}, p_{t} > 230; Z_{g}; Events", 30, 0, 0.5);
    TH1F *h_zg_pt270 = new TH1F(("h_zg_pt270_" + tag).c_str(), "Subjets Z_{g}, p_{t} > 270; Z_{g}; Events", 30, 0, 0.5);
    TH1F *h_zg_pt300 = new TH1F(("h_zg_pt300_" + tag).c_str(), "Subjets Z_{g}, p_{t} > 300; Z_{g}; Events", 30, 0, 0.5);

    TH1F *h_ptheta_sj = new TH1F(("h_ptheta_sj_" + tag).c_str(), "Subjets p_{#theta}; p_{#theta}; Events", 30, 0, 1);
    TH1F *h_ptheta_sjsd = new TH1F(("h_ptheta_sjsd_" + tag).c_str(), "Subjets p_{#theta}; p_{#theta}; Events", 30, 0, 1);
    TH1F *h_ptheta_pt200 = new TH1F(("h_ptheta_pt200_" + tag).c_str(), "Subjets p_{#theta}, p_{t} > 200; p_{#theta}; Events", 30, 0, 1);
    TH1F *h_ptheta_pt230 = new TH1F(("h_ptheta_pt230_" + tag).c_str(), "Subjets p_{#theta}, p_{t} > 230; p_{#theta}; Events", 30, 0, 1);
    TH1F *h_ptheta_pt270 = new TH1F(("h_ptheta_pt270_" + tag).c_str(), "Subjets p_{#theta}, p_{t} > 270; p_{#theta}; Events", 30, 0, 1);
    TH1F *h_ptheta_pt300 = new TH1F(("h_ptheta_pt300_" + tag).c_str(), "Subjets p_{#theta}, p_{t} > 300; p_{#theta}; Events", 30, 0, 1);



    const Int_t maxFatJets = 100;
    const Int_t maxJets = 100;
    const Int_t maxGenJets = 100;
    const Int_t maxGenFatJets = 100;
    const Int_t maxSubJets = 100;

    Float_t FatJet_pt[maxFatJets], FatJet_eta[maxFatJets], FatJet_phi[maxFatJets], FatJet_mass[maxFatJets], FatJet_msoftdrop[maxFatJets], FatJet_tau1[maxFatJets], FatJet_tau2[maxFatJets], FatJet_tau3[maxFatJets];
    Float_t Jet_pt[maxJets], Jet_eta[maxJets], Jet_phi[maxJets], Jet_mass[maxJets];
    Float_t SubJet_pt[maxSubJets], SubJet_eta[maxSubJets], SubJet_phi[maxSubJets], SubJet_mass[maxSubJets];
    Float_t GenFatJet_partonFlavour[maxGenFatJets];
    Short_t FatJet_subJetIdx1[maxFatJets], FatJet_subJetIdx2[maxFatJets];   
    Float_t GenFatJet_pt[maxGenFatJets], GenFatJet_eta[maxGenFatJets], GenFatJet_phi[maxGenFatJets], GenFatJet_mass[maxGenFatJets];
    Float_t GenJet_pt[maxGenJets], GenJet_eta[maxGenJets], GenJet_phi[maxGenJets], GenJet_mass[maxGenJets];

    Int_t nGenFatJet, nFatJet, nGenJet, nJet, nSubJet;

    TChain *chain = new TChain("Events"); 
    for(auto sample : samples_ssWW){
        chain->Add(sample.c_str()); 
    }
    




    chain->SetBranchAddress("nFatJet", &nFatJet);
    chain->SetBranchAddress("FatJet_mass", FatJet_mass);
    chain->SetBranchAddress("FatJet_msoftdrop", FatJet_msoftdrop);
    chain->SetBranchAddress("FatJet_pt", FatJet_pt);
    chain->SetBranchAddress("FatJet_eta", FatJet_eta);
    chain->SetBranchAddress("FatJet_phi", FatJet_phi);
    chain->SetBranchAddress("FatJet_tau1", FatJet_tau1);
    chain->SetBranchAddress("FatJet_tau2", FatJet_tau2);
    chain->SetBranchAddress("FatJet_tau3", FatJet_tau3);
    chain->SetBranchAddress("FatJet_subJetIdx1", FatJet_subJetIdx1);
    chain->SetBranchAddress("FatJet_subJetIdx2", FatJet_subJetIdx2);


    
    chain->SetBranchAddress("nJet", &nJet);
    chain->SetBranchAddress("Jet_mass", Jet_mass);
    chain->SetBranchAddress("Jet_pt", Jet_pt);
    chain->SetBranchAddress("Jet_eta", Jet_eta);
    chain->SetBranchAddress("Jet_phi", Jet_phi);


    chain->SetBranchAddress("nSubJet", &nSubJet);
    chain->SetBranchAddress("SubJet_pt", SubJet_pt);
    chain->SetBranchAddress("SubJet_eta", SubJet_eta);
    chain->SetBranchAddress("SubJet_phi", SubJet_phi);
    chain->SetBranchAddress("SubJet_mass", SubJet_mass);

    chain->SetBranchAddress("nGenJetAK8", &nGenFatJet);
    chain->SetBranchAddress("GenJetAK8_mass", GenFatJet_mass);
    chain->SetBranchAddress("GenJetAK8_pt", GenFatJet_pt);
    chain->SetBranchAddress("GenJetAK8_eta", GenFatJet_eta);
    chain->SetBranchAddress("GenJetAK8_phi", GenFatJet_phi);

 
    chain->SetBranchAddress("nGenJet", &nGenJet);
    chain->SetBranchAddress("GenJet_mass", GenJet_mass);
    chain->SetBranchAddress("GenJet_pt", GenJet_pt);
    chain->SetBranchAddress("GenJet_eta", GenJet_eta);
    chain->SetBranchAddress("GenJet_phi", GenJet_phi);







    for(int i=0; i< int(chain->GetEntries()); i++){
        chain->GetEntry(i);
        if(i%5000==0) std::cout<<"Processing event "<<i<<std::endl;
        //if(i==maxEvents) break;
        for(Int_t fj=0; fj<nFatJet; fj++){
            h_FJ_mass->Fill(FatJet_mass[fj]);
            h_FJ_msd->Fill(FatJet_msoftdrop[fj]);
            h_FJ_pt->Fill(FatJet_pt[fj]);
            h_FJ_eta->Fill(FatJet_eta[fj]);
            h_FJ_phi->Fill(FatJet_phi[fj]);

            h_tau1->Fill(FatJet_tau1[fj]);
            h_tau2->Fill(FatJet_tau2[fj]);
            h_tau3->Fill(FatJet_tau3[fj]);
            h_tau21->Fill(FatJet_tau2[fj]/FatJet_tau1[fj]);
            h_tau32->Fill(FatJet_tau3[fj]/FatJet_tau2[fj]);
            if(FatJet_msoftdrop[fj]>60 && FatJet_msoftdrop[fj]<100)
            {
                double min_pt = std::min(SubJet_pt[FatJet_subJetIdx1[fj]], SubJet_pt[FatJet_subJetIdx2[fj]]);
                h_zg_sjsd->Fill(min_pt/(SubJet_pt[FatJet_subJetIdx1[fj]] + SubJet_pt[FatJet_subJetIdx2[fj]]));

                double E1, E2;
                E1 = std::sqrt(SubJet_pt[FatJet_subJetIdx1[fj]]*SubJet_pt[FatJet_subJetIdx1[fj]] + SubJet_mass[FatJet_subJetIdx1[fj]]*SubJet_mass[FatJet_subJetIdx1[fj]]);
                E2 = std::sqrt(SubJet_pt[FatJet_subJetIdx2[fj]]*SubJet_pt[FatJet_subJetIdx2[fj]] + SubJet_mass[FatJet_subJetIdx2[fj]]*SubJet_mass[FatJet_subJetIdx2[fj]]);
                h_ptheta_sjsd->Fill(std::abs(E1-E2)/FatJet_pt[fj]);


                h_tau21m->Fill(FatJet_tau2[fj]/FatJet_tau1[fj]);
                h_tau32m->Fill(FatJet_tau3[fj]/FatJet_tau2[fj]);

                h_SJ_etad_sd->Fill(SubJet_eta[FatJet_subJetIdx1[fj]] - SubJet_eta[FatJet_subJetIdx2[fj]]);
                h_SJ_phid_sd->Fill(SubJet_phi[FatJet_subJetIdx1[fj]] - SubJet_phi[FatJet_subJetIdx2[fj]]);

                if(FatJet_pt[fj] > 200){
                    double min_pt = std::min(SubJet_pt[FatJet_subJetIdx1[fj]], SubJet_pt[FatJet_subJetIdx2[fj]]);
                    h_zg_pt200->Fill(min_pt/(SubJet_pt[FatJet_subJetIdx1[fj]] + SubJet_pt[FatJet_subJetIdx2[fj]]));

                    double E1, E2;
                    E1 = std::sqrt(SubJet_pt[FatJet_subJetIdx1[fj]]*SubJet_pt[FatJet_subJetIdx1[fj]] + SubJet_mass[FatJet_subJetIdx1[fj]]*SubJet_mass[FatJet_subJetIdx1[fj]]);
                    E2 = std::sqrt(SubJet_pt[FatJet_subJetIdx2[fj]]*SubJet_pt[FatJet_subJetIdx2[fj]] + SubJet_mass[FatJet_subJetIdx2[fj]]*SubJet_mass[FatJet_subJetIdx2[fj]]);
                    h_ptheta_pt200->Fill(std::abs(E1-E2)/FatJet_pt[fj]);
                }
                if(FatJet_pt[fj] > 230){
                    double min_pt = std::min(SubJet_pt[FatJet_subJetIdx1[fj]], SubJet_pt[FatJet_subJetIdx2[fj]]);
                    h_zg_pt230->Fill(min_pt/(SubJet_pt[FatJet_subJetIdx1[fj]] + SubJet_pt[FatJet_subJetIdx2[fj]]));

                    double E1, E2;
                    E1 = std::sqrt(SubJet_pt[FatJet_subJetIdx1[fj]]*SubJet_pt[FatJet_subJetIdx1[fj]] + SubJet_mass[FatJet_subJetIdx1[fj]]*SubJet_mass[FatJet_subJetIdx1[fj]]);
                    E2 = std::sqrt(SubJet_pt[FatJet_subJetIdx2[fj]]*SubJet_pt[FatJet_subJetIdx2[fj]] + SubJet_mass[FatJet_subJetIdx2[fj]]*SubJet_mass[FatJet_subJetIdx2[fj]]);
                    h_ptheta_pt230->Fill(std::abs(E1-E2)/FatJet_pt[fj]);
                }
                if(FatJet_pt[fj] > 270){
                    double min_pt = std::min(SubJet_pt[FatJet_subJetIdx1[fj]], SubJet_pt[FatJet_subJetIdx2[fj]]);
                    h_zg_pt270->Fill(min_pt/(SubJet_pt[FatJet_subJetIdx1[fj]] + SubJet_pt[FatJet_subJetIdx2[fj]]));

                    double E1, E2;
                    E1 = std::sqrt(SubJet_pt[FatJet_subJetIdx1[fj]]*SubJet_pt[FatJet_subJetIdx1[fj]] + SubJet_mass[FatJet_subJetIdx1[fj]]*SubJet_mass[FatJet_subJetIdx1[fj]]);
                    E2 = std::sqrt(SubJet_pt[FatJet_subJetIdx2[fj]]*SubJet_pt[FatJet_subJetIdx2[fj]] + SubJet_mass[FatJet_subJetIdx2[fj]]*SubJet_mass[FatJet_subJetIdx2[fj]]);
                    h_ptheta_pt270->Fill(std::abs(E1-E2)/FatJet_pt[fj]);
                }
                if(FatJet_pt[fj] > 300){
                    double min_pt = std::min(SubJet_pt[FatJet_subJetIdx1[fj]], SubJet_pt[FatJet_subJetIdx2[fj]]);
                    h_zg_pt300->Fill(min_pt/(SubJet_pt[FatJet_subJetIdx1[fj]] + SubJet_pt[FatJet_subJetIdx2[fj]]));

                    double E1, E2;
                    E1 = std::sqrt(SubJet_pt[FatJet_subJetIdx1[fj]]*SubJet_pt[FatJet_subJetIdx1[fj]] + SubJet_mass[FatJet_subJetIdx1[fj]]*SubJet_mass[FatJet_subJetIdx1[fj]]);
                    E2 = std::sqrt(SubJet_pt[FatJet_subJetIdx2[fj]]*SubJet_pt[FatJet_subJetIdx2[fj]] + SubJet_mass[FatJet_subJetIdx2[fj]]*SubJet_mass[FatJet_subJetIdx2[fj]]);
                    h_ptheta_pt300->Fill(std::abs(E1-E2)/FatJet_pt[fj]);
            }


            }
            


            h_tau21_msd->Fill((FatJet_tau2[fj]/FatJet_tau1[fj]), FatJet_msoftdrop[fj]);
            h_tau32_msd->Fill((FatJet_tau3[fj]/FatJet_tau2[fj]), FatJet_msoftdrop[fj]);

            h_SJ_etad->Fill(SubJet_eta[FatJet_subJetIdx1[fj]] - SubJet_eta[FatJet_subJetIdx2[fj]]);
            h_SJ_phid->Fill(SubJet_phi[FatJet_subJetIdx1[fj]] - SubJet_phi[FatJet_subJetIdx2[fj]]);
            
            double min_pt = std::min(SubJet_pt[FatJet_subJetIdx1[fj]], SubJet_pt[FatJet_subJetIdx2[fj]]);
            h_zg_sj->Fill(min_pt/(SubJet_pt[FatJet_subJetIdx1[fj]] + SubJet_pt[FatJet_subJetIdx2[fj]]));

            double E1, E2;
            E1 = std::sqrt(SubJet_pt[FatJet_subJetIdx1[fj]]*SubJet_pt[FatJet_subJetIdx1[fj]] + SubJet_mass[FatJet_subJetIdx1[fj]]*SubJet_mass[FatJet_subJetIdx1[fj]]);
            E2 = std::sqrt(SubJet_pt[FatJet_subJetIdx2[fj]]*SubJet_pt[FatJet_subJetIdx2[fj]] + SubJet_mass[FatJet_subJetIdx2[fj]]*SubJet_mass[FatJet_subJetIdx2[fj]]);
            h_ptheta_sj->Fill(std::abs(E1-E2)/FatJet_pt[fj]);

        }
        for(Int_t genfj=0; genfj<nGenFatJet; genfj++){
            h_GFJ_mass->Fill(GenFatJet_mass[genfj]);
            h_GFJ_pt->Fill(GenFatJet_pt[genfj]);
            h_GFJ_eta->Fill(GenFatJet_eta[genfj]);
            h_GFJ_phi->Fill(GenFatJet_phi[genfj]);
        }

    }





    /*ISTOGRAMMI PER pt*/
    TCanvas *c_pt_tot_nonp = new TCanvas("c_pt_tot_nonp", "Consituents p_{T}", 2000, 1000);
    c_pt_tot_nonp->Divide(2,1);

    c_pt_tot_nonp->cd(1);
    gPad->SetLeftMargin(0.15);
    h_GFJ_pt->SetTitle("GenJetAK8 p_{t}");
    h_GFJ_pt->Scale(1.0 / h_GFJ_pt->Integral());
    h_GFJ_pt->Draw("hist");
    h_GFJ_pt->SetLineColor(kViolet-3);
    h_GFJ_pt->SetFillColor(kViolet-3);
    h_GFJ_pt->GetXaxis()->SetTitleSize(0.04);
    h_GFJ_pt->GetXaxis()->SetTitleOffset(0.9);


    c_pt_tot_nonp->cd(2);
    gPad->SetLeftMargin(0.15);
    h_FJ_pt->SetTitle("FatJet p_{t}");
    h_FJ_pt->Scale(1.0 / h_FJ_pt->Integral());
    h_FJ_pt->Draw("hist");
    h_FJ_pt->SetLineColor(kViolet-3);
    h_FJ_pt->SetFillColor(kViolet-3);
    h_FJ_pt->GetXaxis()->SetTitleSize(0.04);
    h_FJ_pt->GetXaxis()->SetTitleOffset(0.9);

    c_pt_tot_nonp->SaveAs("plots_pt/c_pt_tot_nonp.pdf");
    


    /*ISTOGRAMMI PER ETA*/
    TCanvas *c_eta_tot_nonp = new TCanvas("c_eta_tot_nonp", "Consituents #eta: pseudorapirity", 2000, 1000);
    c_eta_tot_nonp->Divide(2,1);

    c_eta_tot_nonp->cd(1);
    gPad->SetLeftMargin(0.15);
    h_GFJ_eta->SetTitle("GenJetAK8 #eta");
    h_GFJ_eta->Scale(1.0 / h_GFJ_eta->Integral());
    h_GFJ_eta->Draw("hist");
    h_GFJ_eta->SetLineColor(kViolet-3);
    h_GFJ_eta->SetFillColor(kViolet-3);
    h_GFJ_eta->GetXaxis()->SetTitleSize(0.04);


    c_eta_tot_nonp->cd(2);
    gPad->SetLeftMargin(0.15);
    h_FJ_eta->SetTitle("FatJet #eta");
    h_FJ_eta->Draw("hist");
    h_FJ_eta->Scale(1.0 / h_FJ_eta->Integral());
    h_FJ_eta->SetLineColor(kViolet-3);
    h_FJ_eta->SetFillColor(kViolet-3);
    h_FJ_eta->GetXaxis()->SetTitleSize(0.09);


    c_eta_tot_nonp->SaveAs("plots_eta/c_eta_tot_nonp.pdf");


    /*ISTOGRAMMI PER PHI*/
    TCanvas *c_phi_tot_nonp = new TCanvas("c_phi_tot_nonp", "Consituents #phi", 2000, 1000);
    c_phi_tot_nonp->Divide(2,1);

    c_phi_tot_nonp->cd(1);
    gPad->SetLeftMargin(0.15);
    h_GFJ_phi->Scale(1.0 / h_GFJ_phi->Integral());
    h_GFJ_phi->SetTitle("GenJetAK8 #phi");
    h_GFJ_phi->Draw("hist");
    h_GFJ_phi->SetLineColor(kViolet-3);
    h_GFJ_phi->SetFillColor(kViolet-3);
    h_GFJ_phi->GetXaxis()->SetTitleSize(0.05);
    h_GFJ_phi->GetXaxis()->SetTitleOffset(0.9);


    c_phi_tot_nonp->cd(2);
    gPad->SetLeftMargin(0.15);
    h_FJ_phi->Scale(1.0 / h_GFJ_phi->Integral());
    h_FJ_phi->SetTitle("FatJet #phi");
    h_FJ_phi->Draw("hist");
    h_FJ_phi->SetLineColor(kViolet-3);
    h_FJ_phi->SetFillColor(kViolet-3);
    h_FJ_phi->GetXaxis()->SetTitleSize(0.05);
    h_FJ_phi->GetXaxis()->SetTitleOffset(0.9);

    c_phi_tot_nonp->SaveAs("plots_phi/c_phi_tot_nonp.pdf");
    


    
    /*ISTOGRAMMI PER LE MASSE*/
    TCanvas *c_mass_tot_nonp = new TCanvas("c_mass_tot_nonp", "Consituents Mass", 2000, 1100);
    c_mass_tot_nonp->Divide(2,1);

    c_mass_tot_nonp->cd(1);
    gPad->SetLeftMargin(0.15);
    h_GFJ_mass->SetTitle("GenJetAK8 mass");
    h_GFJ_mass->Scale(1.0 / h_GFJ_mass->Integral());
    h_GFJ_mass->Draw("hist");
    h_GFJ_mass->SetLineColor(kViolet-3);
    h_GFJ_mass->SetFillColor(kViolet-3);
    h_GFJ_mass->GetXaxis()->SetTitleSize(0.05);
    h_GFJ_mass->GetXaxis()->SetTitleOffset(0.9);

    c_mass_tot_nonp->cd(2);
    gPad->SetLeftMargin(0.15);
    h_FJ_mass->SetTitle("FatJet mass");
    h_FJ_mass->Scale(1.0 / h_FJ_mass->Integral());
    h_FJ_mass->Draw("hist");
    h_FJ_mass->SetLineColor(kViolet-3);
    h_FJ_mass->SetFillColor(kViolet-3);
    h_FJ_mass->GetXaxis()->SetTitleSize(0.05);
    h_FJ_mass->GetXaxis()->SetTitleOffset(0.9);
    
    c_mass_tot_nonp->SaveAs("plots_m/c_mass_tot_nonp.pdf");






    /*ISTOGRAMMI TAU 1,2, 3 */
    TCanvas *c_tau_nonp = new TCanvas("c_tau_nonp", "FatJet #tau", 3000, 1100);
    c_tau_nonp->Divide(3,1);

    c_tau_nonp->cd(1);
    gPad->SetLeftMargin(0.15);
    h_tau1->SetTitleSize(0.7);
    h_tau1->Scale(1.0 / h_tau1->Integral());
    h_tau1->Draw("hist");
    h_tau1->SetLineColor(kViolet-3);
    h_tau1->SetFillColor(kViolet-3);
    h_tau1->GetXaxis()->SetTitleSize(0.05);

    c_tau_nonp->cd(2);
    gPad->SetLeftMargin(0.15);
    h_tau2->SetTitleSize(0.7);
    h_tau2->Scale(1.0 / h_tau2->Integral());
    h_tau2->Draw("hist");
    h_tau2->SetLineColor(kViolet-3);
    h_tau2->SetFillColor(kViolet-3);
    h_tau2->GetXaxis()->SetTitleSize(0.05);

    c_tau_nonp->cd(3);
    gPad->SetLeftMargin(0.15);
    h_tau3->SetTitleSize(0.7);
    h_tau3->Scale(1.0 / h_tau3->Integral());
    h_tau3->Draw("hist");
    h_tau3->SetLineColor(kViolet-3);
    h_tau3->SetFillColor(kViolet-3);
    h_tau3->GetXaxis()->SetTitleSize(0.05);

    c_tau_nonp->SaveAs("plots_tau/c_tau_nonp.pdf");



    /*ISTOGRAMMI Tau21, Tau32 confronto con eventi "W" */
    TCanvas *c_tau_ratio_nonp = new TCanvas("c_tau_ratio_nonp", "FatJet #tau ratios", 2000, 1800);
    c_tau_ratio_nonp->Divide(2,2);

    c_tau_ratio_nonp->cd(1);
    h_tau21->SetTitleSize(0.7);
    h_tau21->Scale(1.0 / h_tau21->Integral());
    h_tau21->Draw("hist");
    h_tau21->SetLineColor(kViolet-3);
    h_tau21->SetFillColor(kViolet-3);
    h_tau21->GetXaxis()->SetTitleSize(0.05);


    c_tau_ratio_nonp->cd(2);
    h_tau32->SetTitleSize(0.7);
    h_tau32->Scale(1.0 / h_tau32->Integral());
    h_tau32->Draw("hist");
    h_tau32->SetLineColor(kViolet-3);
    h_tau32->SetFillColor(kViolet-3);
    h_tau32->GetXaxis()->SetTitleSize(0.05);



    c_tau_ratio_nonp->cd(3);
    h_tau21m->Scale(1.0 / h_tau21m->Integral());
    h_tau21m->Draw("hist");
    h_tau21m->SetLineColor(kViolet-3);
    h_tau21m->SetFillColor(kViolet-3);
    h_tau21m->GetXaxis()->SetTitleSize(0.05);



    c_tau_ratio_nonp->cd(4);
    h_tau32m->Scale(1.0 / h_tau32m->Integral());
    h_tau32m->Draw("hist");
    h_tau32m->SetLineColor(kViolet-3);
    h_tau32m->SetFillColor(kViolet-3);
    h_tau32m->GetXaxis()->SetTitleSize(0.05);

    c_tau_ratio_nonp->SaveAs("plots_tau/c_tau_ratio_nonp.pdf");


    TCanvas *c_zg_pt_nonp = new TCanvas("c_zg_pt_nonp", "Z_{g}", 2000, 1800);
    c_zg_pt_nonp->Divide(2,2);

    c_zg_pt_nonp->cd(1);
    h_zg_pt200->Scale(1.0 / h_zg_pt200->Integral());
    h_zg_pt200->Draw("hist");
    gPad->SetLeftMargin(0.15);
    h_zg_pt200->SetLineColor(kViolet-3);
    h_zg_pt200->SetFillColor(kViolet-3);
    h_zg_pt200->GetXaxis()->SetTitleSize(0.04);

    c_zg_pt_nonp->cd(2);
    h_zg_pt230->Scale(1.0 / h_zg_pt230->Integral());
    h_zg_pt230->Draw("hist");
    gPad->SetLeftMargin(0.15);
    h_zg_pt230->SetLineColor(kViolet-3);
    h_zg_pt230->SetFillColor(kViolet-3);
    h_zg_pt230->GetXaxis()->SetTitleSize(0.04);

    c_zg_pt_nonp->cd(3);
    h_zg_pt270->Scale(1.0 / h_zg_pt270->Integral());
    h_zg_pt270->Draw("hist");
    gPad->SetLeftMargin(0.15);
    h_zg_pt270->SetLineColor(kViolet-3);
    h_zg_pt270->SetFillColor(kViolet-3);
    h_zg_pt270->GetXaxis()->SetTitleSize(0.04);

    c_zg_pt_nonp->cd(4);
    h_zg_pt300->Scale(1.0 / h_zg_pt300->Integral());
    h_zg_pt300->Draw("hist");
    gPad->SetLeftMargin(0.15);
    h_zg_pt300->SetLineColor(kViolet-3);
    h_zg_pt300->SetFillColor(kViolet-3);
    h_zg_pt300->GetXaxis()->SetTitleSize(0.04);

    c_zg_pt_nonp->SaveAs("plots_z1/c_zg_pt_nonp.pdf");
    

    
    
    
    
    TCanvas *c_ptheta_pt_nonp = new TCanvas("c_ptheta_pt_nonp", "p_{#theta}", 2000, 1800);
    c_ptheta_pt_nonp->Divide(2,2);

    c_ptheta_pt_nonp->cd(1);
    h_ptheta_pt200->Scale(1.0 / h_ptheta_pt200->Integral());
    h_ptheta_pt200->Draw("hist");
    gPad->SetLeftMargin(0.15);
    h_ptheta_pt200->SetLineColor(kViolet-3);
    h_ptheta_pt200->SetFillColor(kViolet-3);
    h_ptheta_pt200->GetXaxis()->SetTitleSize(0.04);

    c_ptheta_pt_nonp->cd(2);
    h_ptheta_pt230->Scale(1.0 / h_ptheta_pt230->Integral());
    h_ptheta_pt230->Draw("hist");
    gPad->SetLeftMargin(0.15);
    h_ptheta_pt230->SetLineColor(kViolet-3);
    h_ptheta_pt230->SetFillColor(kViolet-3);
    h_ptheta_pt230->GetXaxis()->SetTitleSize(0.04);

    c_ptheta_pt_nonp->cd(3);
    h_ptheta_pt270->Scale(1.0 / h_ptheta_pt270->Integral());
    h_ptheta_pt270->Draw("hist");
    gPad->SetLeftMargin(0.15);
    h_ptheta_pt270->SetLineColor(kViolet-3);
    h_ptheta_pt270->SetFillColor(kViolet-3);
    h_ptheta_pt270->GetXaxis()->SetTitleSize(0.04);

    c_ptheta_pt_nonp->cd(4);
    h_ptheta_pt300->Scale(1.0 / h_ptheta_pt300->Integral());
    h_ptheta_pt300->Draw("hist");
    gPad->SetLeftMargin(0.15);
    h_ptheta_pt300->SetLineColor(kViolet-3);
    h_ptheta_pt300->SetFillColor(kViolet-3);
    h_ptheta_pt300->GetXaxis()->SetTitleSize(0.04);

    c_ptheta_pt_nonp->SaveAs("plots_ptheta/c_ptheta_pt_nonp.pdf");





    //ISTOGRAMMI per SubJet
//eta
    gStyle->SetOptStat(0);
    TCanvas *c_SJ_etad = new TCanvas("c_SJ_etad", "SubJet #eta difference", 800, 600);

    h_SJ_etad->Scale(1.0 / h_SJ_etad->Integral());
    h_SJ_etad->Draw("hist");
    h_SJ_etad->SetLineColor(kViolet-3);
    h_SJ_etad->SetFillColor(kViolet-3);
    h_SJ_etad->SetTitle("SubJets #eta difference");
    h_SJ_etad_sd->Scale(1.0 / h_SJ_etad->Integral());
    h_SJ_etad_sd->Draw("P same");
    h_SJ_etad_sd->SetMarkerStyle(20);
    h_SJ_etad_sd->SetMarkerColor(kBlack);
    h_SJ_etad_sd->SetLineColor(kBlack);
    h_SJ_etad_sd->SetTitleSize(0.05);

    TLegend *leg4 = new TLegend(0.75, 0.75, 0.9, 0.9);
    leg4->AddEntry(h_SJ_etad_sd, "60 < m_{FJ} < 100", "l");
    leg4->AddEntry(h_SJ_etad, "No cuts", "l");
    leg4->Draw();

    c_SJ_etad->SaveAs("plots_sj/c_SJ_etad_ nonp.pdf");


//phi
    /*TCanvas *c_SJ_phid = new TCanvas("c_SJ_phid", "SubJet phi difference", 800, 600);
    c_SJ_phid->Divide(2,1);

    c_SJ_phid->cd(1);
    h_SJ_phid_sdTT->Draw();
    h_SJ_phid_sdTT->SetTitle("SubJet #phi difference: 60 < m_sd < 100");
    h_SJ_phid_sdLL->Draw("same");
    h_SJ_phid_sdLL->SetLineColor(kRed);


    c_SJ_phid->cd(2);
    h_SJ_phidTT->Draw();
    h_SJ_phidTT->SetTitle("SubJet #phi difference: all masses");
    h_SJ_phidLL->Draw("same");
    h_SJ_phidLL->SetLineColor(kRed);

    TLegend *leg5 = new TLegend(0.2, 0.65, 0.35, 0.8);
    leg5->AddEntry(h_SJ_phid_sdTT, "TT", "l");
    leg5->AddEntry(h_SJ_phid_sdLL, "LL", "l");
    leg5->Draw();

    c_SJ_phid->SaveAs("plots_sj/c_SJ_phid_nonp.pdf");*/







    //SCATTER PLOT TAU - MASSA SD
    TCanvas *c_tau_msd_nonp = new TCanvas("c_tau_msd_nonp","Tau vs msd", 2200,1000);
    c_tau_msd_nonp->Divide(2,1);

    c_tau_msd_nonp->cd(1);
    h_tau21_msd->Draw("COLZ"); 
    h_tau21_msd->GetXaxis()->SetTitleSize(0.04);

    c_tau_msd_nonp->cd(2);
    h_tau32_msd->Draw("COLZ"); 
    h_tau32_msd->GetXaxis()->SetTitleSize(0.04);

    c_tau_msd_nonp->SaveAs("plots_tau/c_tau_msd_nonp.pdf");



    /*Soft Drop*/
    TCanvas *c_FJ_msd_nonp = new TCanvas("c_FJ_msd_nonp", "FatJet SoftDrop mass", 800, 600);
    h_FJ_msd->SetTitle("SoftDrop vs FatJet mass");
    h_FJ_msd->Scale(1.0 / h_FJ_msd->Integral());
    h_FJ_msd->Draw("hist");
    h_FJ_msd->SetLineColor(kViolet-3);
    h_FJ_msd->SetFillColor(kViolet-3);
    h_FJ_mass->Scale(1.0 / h_FJ_mass->Integral());
    h_FJ_mass->Draw("same");
    h_FJ_mass->SetLineColor(kBlack);
    h_FJ_mass->SetMarkerColor(kBlack);
    h_FJ_mass->SetMarkerStyle(20);

    TLegend *leg9 = new TLegend(0.5, 0.6, 0.7, 0.8);
    leg9->AddEntry(h_FJ_msd, "SoftDrop", "l");
    leg9->AddEntry(h_FJ_mass, "FatJet", "l");
    leg9->SetTextSize(0.04);
    leg9->SetTextFont(42);
    leg9->Draw();

    c_FJ_msd_nonp->SaveAs("plots_m/c_FJ_msd_nonp.pdf");


    //zg e ptheta
    TCanvas *c_zg_nonp = new TCanvas("c_zg_nonp", "Z_{g}", 800, 600);
    h_zg_sjsd->Scale(1.0 / h_zg_sjsd->Integral());
    h_zg_sjsd->Draw("hist");
    gPad->SetLeftMargin(0.15);
    h_zg_sjsd->SetLineColor(kViolet-3);
    h_zg_sjsd->SetFillColor(kViolet-3);
    h_zg_sjsd->GetXaxis()->SetTitleSize(0.04);
    h_zg_sj->Scale(1.0 / h_zg_sj->Integral());
    h_zg_sj->Draw("P same");
    h_zg_sj->SetMarkerStyle(20);
    h_zg_sj->SetMarkerSize(1.5);
    h_zg_sj->SetMarkerColor(kBlack);
    h_zg_sj->SetLineColor(kBlack);

    TLegend *leg20 = new TLegend(0.7, 0.1, 0.9, 0.25);
    leg20->AddEntry(h_zg_sjsd, "60 < m_{sd} < 100", "l");
    leg20->AddEntry(h_zg_sj, "No cuts", "l");
    leg20->Draw();
    
    c_zg_nonp->SaveAs("plots_z1/c_zg_nonp.pdf");




    TCanvas *c_ptheta_nonp = new TCanvas("c_ptheta_nonp", "p_{#theta}", 800, 600);
    h_ptheta_sjsd->Scale(1.0 / h_ptheta_sjsd->Integral());
    gPad->SetLeftMargin(0.15);
    h_ptheta_sjsd->Draw("hist");
    h_ptheta_sjsd->SetLineColor(kViolet-3);
    h_ptheta_sjsd->SetFillColor(kViolet-3);
    h_ptheta_sj->Scale(1.0 / h_ptheta_sj->Integral());
    h_ptheta_sj->Draw("P same");
    h_ptheta_sj->SetMarkerStyle(20);
    h_ptheta_sj->SetMarkerSize(1.5);
    h_ptheta_sj->SetMarkerColor(kBlack);
    h_ptheta_sj->SetLineColor(kBlack);
    h_ptheta_sj->GetXaxis()->SetTitleOffset(0.09);
    h_ptheta_sj->GetXaxis()->SetTitleSize(0.05);

    TLegend *leg19 = new TLegend(0.15, 0.1, 0.35, 0.25);
    leg19->AddEntry(h_ptheta_sjsd, "60 < m_{sd} < 100", "l");
    leg19->AddEntry(h_ptheta_sj, "No cuts", "l");
    leg19->Draw();
    
    c_ptheta_nonp->SaveAs("plots_ptheta/c_ptheta_nonp.pdf");




    f->Write();
    f->Close();


    gROOT->ProcessLine(".q");



}
