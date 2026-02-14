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
#include "canva.cpp"




// lista dei branch nel tree : https://cms-xpog.docs.cern.ch/autoDoc/NanoAODv12/2022/2023/doc_DYJetsToLL_M-50_TuneCP5_13p6TeV-madgraphMLM-pythia8_Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2.html

void jet_plots3(
    TString output="jets3.root"
){
    TFile *f=new TFile(output, "recreate");

    Int_t maxEvents = 2000;
    std::string tag = "ssWW";
    

//pt    
    TH1F *h_FJ_pt = new TH1F(("h_FatJet_pt_" + tag).c_str(), "FatJet p_{T}  (CMS private work);p_{T} [GeV];Events", 50, 170, 1000);
    
    TH1F *h_GFJ_pt = new TH1F(("h_GenJetAK8_pt_" + tag).c_str(), "GenJetAK8 p_{T}  (CMS private work);p_{T} [GeV];Events", 50, 150, 1000);
    
    
//eta    
    TH1F *h_FJ_eta = new TH1F(("h_FatJet_eta_" + tag).c_str(), "FatJet #eta  (CMS private work); #eta Pseudorapidity; Events", 50,-6, 6);
    
    TH1F *h_GFJ_eta = new TH1F(("h_GenJetAK8_eta_" + tag).c_str(), "GenJetAK8 #eta  (CMS private work);#eta Pseudorapidity; Events", 50, -6, 6);

//phi
    TH1F *h_FJ_phi = new TH1F(("h_FatJet_phi_" + tag).c_str(), "FatJet #phi  (CMS private work); #phi; Events", 40, -4,4);
    
    TH1F *h_GFJ_phi = new TH1F(("h_GenJetAK8_phi_" + tag).c_str(), "GenJetAK8 #phi  (CMS private work); #phi; Events", 40, -4, 4);


//mass
    TH1F *h_FJ_mass = new TH1F(("h_FJ_mass_" + tag).c_str(), "FatJet mass  (CMS private work); mass [GeV]; Events", 50,10,200);
    
    TH1F *h_GFJ_mass = new TH1F(("h_GFJ_mass_" + tag).c_str(), "GenJetAK8 FatJet mass  (CMS private work);mass [GeV];Events", 50, 0, 200);

//mass softdrop
    TH1F *h_FJ_msd = new TH1F(("h_FJ_msd_" + tag).c_str(), "FatJet SoftDrop mass  (CMS private work);mass [GeV]", 50, 10, 200);


//tau
    TH1F *h_tau1 =new TH1F(("h_tau1_" + tag).c_str(), "FatJet #tau_{1}  (CMS private work); #tau_{1}; Events", 30, 0, 0.6);
    TH1F *h_tau2 =new TH1F(("h_tau2_" + tag).c_str(), "FatJet #tau_{2}  (CMS private work); #tau_{2}; Events", 30, 0, 0.6);
    TH1F *h_tau3= new TH1F(("h_tau3_" + tag).c_str(), "FatJet #tau_{3}  (CMS private work); #tau_{3}; Events", 30, 0, 0.6);

    TH1F *h_tau21m= new TH1F(("h_tau21m_" + tag).c_str(), "#tau_{21} (60 < m_{sd} < 100)  (CMS private work); #tau_{21}; Events", 30, 0, 1);
    TH1F *h_tau21= new TH1F(("h_tau21_" + tag).c_str(), "#tau_{21}  (CMS private work); #tau_{21}; Events", 30, 0, 1);
    TH1F *h_tau32m= new TH1F(("h_tau32m_" + tag).c_str(), "#tau_{32} (60 < m_{sd} < 100)  (CMS private work); #tau_{32}; Events", 30, 0, 1);
    TH1F *h_tau32= new TH1F(("h_tau32_" + tag).c_str(), "#tau_{32}  (CMS private work); #tau_{32}; Events", 30, 0, 1);

    TH2F *h_tau21_msd = new TH2F(("h_msd_tau21_" + tag).c_str(), "SoftDrop mass   vs   #tau_{21}    (CMS private work); #tau_{21}; SoftDrop", 50, 0, 1,    30, 20, 300);

    TH2F *h_tau32_msd = new TH2F(("h_msd_tau32_" + tag).c_str(), "SoftDrop mass   vs   #tau_{32}    (CMS private work); #tau_{32}; SoftDrop", 50, 0, 1,    30, 20, 300);


//subjet
    TH1F *h_SJ_etad_sd = new TH1F(("h_SJ_etad_sd_" + tag).c_str(), "SubJet #eta difference; #eta; Events", 20,-1,1);
    TH1F *h_SJ_etad = new TH1F(("h_SJ_etad_" + tag).c_str(), "SubJet #eta difference; #eta; Events", 30,-1,1);
    
    TH1F *h_SJ_phid_sd = new TH1F(("h_SJ_phid_sd_" + tag).c_str(), "SubJet #phi difference; #phi; Events", 20,-1,1);
    TH1F *h_SJ_phid = new TH1F(("h_SJ_phid_" + tag).c_str(), "SubJet #phi difference; #phi; Events", 30,-1,1);

    TH1F *h_SJ_pts1 = new TH1F(("h_SJ_pts1_" + tag).c_str(), "SubJet #1 p_{t} (CMS private work); p_{t} [GeV];Events", 30, 0, 600);
    TH1F *h_SJ_pts2 = new TH1F(("h_SJ_pts2_" + tag).c_str(), "SubJet #2 p_{t} (CMS private work); p_{t} [GeV];Events", 30, 0, 600);

    TH1F *h_SJ_mass1 = new TH1F(("h_SJ_mass1" + tag).c_str(), "SubJet #1 mass (CMS private work); mass [GeV];Events", 30, 0, 200);
    TH1F *h_SJ_mass2 = new TH1F(("h_SJ_mass2" + tag).c_str(), "SubJet #2 mass (CMS private work); mass [GeV];Events", 30, 0, 200);

    TH1F *h_SJ_eta1 = new TH1F(("h_SJ_eta1" + tag).c_str(), "SubJet #1  #eta  (CMS private work); #eta;Events", 30, -4, 4);
    TH1F *h_SJ_eta2 = new TH1F(("h_SJ_eta2" + tag).c_str(), "SubJet #2  #eta  (CMS private work); #eta;Events", 30, -4, 4);

    TH1F *h_SJ_phi1 = new TH1F(("h_SJ_phi1" + tag).c_str(), "SubJet #1  #phi  (CMS private work); #phi;Events", 30, -4, 4);
    TH1F *h_SJ_phi2 = new TH1F(("h_SJ_phi2" + tag).c_str(), "SubJet #2  #phi  (CMS private work); #phi;Events", 30, -4, 4);


    //TH1F *h_SJ_ptd_sd = new TH1F(("h_SJ_ptd_sd_" + tag).c_str(), "SubJet p_{t} difference; p_{t}; Events", 30,150,1000);
    //TH1F *h_SJ_ptd = new TH1F(("h_SJ_ptd_" + tag).c_str(), "SubJet p_{t} difference; p_{t}; Events", 30,150,1000);

    /*TH1F *h_zg_sj = new TH1F(("h_zg_sj_" + tag).c_str(), "Subjets z_{g}; z_{g}; Events", 30, 0.1, 0.5);
    TH1F *h_zg_sjsd = new TH1F(("h_zg_sjsd_" + tag).c_str(), "Subjets z_{g}; z_{g}; Events", 30, 0.1, 0.5);
    TH1F *h_zg_pt200 = new TH1F(("h_zg_pt200_" + tag).c_str(), "Subjets z_{g}, p_{t} > 200; z_{g}; Events", 30, 0.1, 0.5);
    TH1F *h_zg_pt230 = new TH1F(("h_zg_pt230_" + tag).c_str(), "Subjets z_{g}, p_{t} > 230; z_{g}; Events", 30, 0.1, 0.5);
    TH1F *h_zg_pt270 = new TH1F(("h_zg_pt270_" + tag).c_str(), "Subjets z_{g}, p_{t} > 270; z_{g}; Events", 30, 0.1, 0.5);
    TH1F *h_zg_pt300 = new TH1F(("h_zg_pt300_" + tag).c_str(), "Subjets z_{g}, p_{t} > 300; z_{g}; Events", 30, 0.1, 0.5);

    TH1F *h_ptheta_sj = new TH1F(("h_ptheta_sj_" + tag).c_str(), "Subjets p_{#theta}; p_{#theta}; Events", 30, 0, 0.9);
    TH1F *h_ptheta_sjsd = new TH1F(("h_ptheta_sjsd_" + tag).c_str(), "Subjets p_{#theta}; p_{#theta}; Events", 30, 0, 0.9);
    TH1F *h_ptheta_pt200 = new TH1F(("h_ptheta_pt200_" + tag).c_str(), "Subjets p_{#theta}   p_{t} > 200; p_{#theta}; Events", 30, 0, 0.9);
    TH1F *h_ptheta_pt230 = new TH1F(("h_ptheta_pt230_" + tag).c_str(), "Subjets p_{#theta}   p_{t} > 230; p_{#theta}; Events", 30, 0, 0.9);
    TH1F *h_ptheta_pt270 = new TH1F(("h_ptheta_pt270_" + tag).c_str(), "Subjets p_{#theta}   p_{t} > 270; p_{#theta}; Events", 30, 0, 0.9);
    TH1F *h_ptheta_pt300 = new TH1F(("h_ptheta_pt300_" + tag).c_str(), "Subjets p_{#theta}   p_{t} > 300; p_{#theta}; Events", 30, 0, 0.9);*/



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

            h_SJ_eta1->Fill(SubJet_eta[FatJet_subJetIdx1[fj]]);
            h_SJ_eta2->Fill(SubJet_eta[FatJet_subJetIdx2[fj]]);

            h_SJ_phi1->Fill(SubJet_phi[FatJet_subJetIdx1[fj]]);
            h_SJ_phi2->Fill(SubJet_phi[FatJet_subJetIdx2[fj]]);

            h_SJ_pts1->Fill(SubJet_pt[FatJet_subJetIdx1[fj]]);
            h_SJ_pts2->Fill(SubJet_pt[FatJet_subJetIdx2[fj]]);

            h_SJ_mass1->Fill(SubJet_mass[FatJet_subJetIdx1[fj]]);
            h_SJ_mass2->Fill(SubJet_mass[FatJet_subJetIdx2[fj]]);

            if(FatJet_msoftdrop[fj]>60 && FatJet_msoftdrop[fj]<100)
            {
                /*double min_pt = std::min(SubJet_pt[FatJet_subJetIdx1[fj]], SubJet_pt[FatJet_subJetIdx2[fj]]);
                h_zg_sjsd->Fill(min_pt/(SubJet_pt[FatJet_subJetIdx1[fj]] + SubJet_pt[FatJet_subJetIdx2[fj]]));

                double E1, E2;
                E1 = std::sqrt(SubJet_pt[FatJet_subJetIdx1[fj]]*SubJet_pt[FatJet_subJetIdx1[fj]] + SubJet_mass[FatJet_subJetIdx1[fj]]*SubJet_mass[FatJet_subJetIdx1[fj]]);
                E2 = std::sqrt(SubJet_pt[FatJet_subJetIdx2[fj]]*SubJet_pt[FatJet_subJetIdx2[fj]] + SubJet_mass[FatJet_subJetIdx2[fj]]*SubJet_mass[FatJet_subJetIdx2[fj]]);
                h_ptheta_sjsd->Fill(std::abs(E1-E2)/FatJet_pt[fj]);*/


                h_tau21m->Fill(FatJet_tau2[fj]/FatJet_tau1[fj]);
                h_tau32m->Fill(FatJet_tau3[fj]/FatJet_tau2[fj]);

                h_SJ_etad_sd->Fill(SubJet_eta[FatJet_subJetIdx1[fj]] - SubJet_eta[FatJet_subJetIdx2[fj]]);
                h_SJ_phid_sd->Fill(SubJet_phi[FatJet_subJetIdx1[fj]] - SubJet_phi[FatJet_subJetIdx2[fj]]);

                /*if(FatJet_pt[fj] > 200){
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
                    h_ptheta_pt300->Fill(std::abs(E1-E2)/FatJet_pt[fj]);*/
            }


            
            


            h_tau21_msd->Fill((FatJet_tau2[fj]/FatJet_tau1[fj]), FatJet_msoftdrop[fj]);
            h_tau32_msd->Fill((FatJet_tau3[fj]/FatJet_tau2[fj]), FatJet_msoftdrop[fj]);

            h_SJ_etad->Fill(SubJet_eta[FatJet_subJetIdx1[fj]] - SubJet_eta[FatJet_subJetIdx2[fj]]);
            h_SJ_phid->Fill(SubJet_phi[FatJet_subJetIdx1[fj]] - SubJet_phi[FatJet_subJetIdx2[fj]]);
            
            /*double min_pt = std::min(SubJet_pt[FatJet_subJetIdx1[fj]], SubJet_pt[FatJet_subJetIdx2[fj]]);
            h_zg_sj->Fill(min_pt/(SubJet_pt[FatJet_subJetIdx1[fj]] + SubJet_pt[FatJet_subJetIdx2[fj]]));

            double E1, E2;
            E1 = std::sqrt(SubJet_pt[FatJet_subJetIdx1[fj]]*SubJet_pt[FatJet_subJetIdx1[fj]] + SubJet_mass[FatJet_subJetIdx1[fj]]*SubJet_mass[FatJet_subJetIdx1[fj]]);
            E2 = std::sqrt(SubJet_pt[FatJet_subJetIdx2[fj]]*SubJet_pt[FatJet_subJetIdx2[fj]] + SubJet_mass[FatJet_subJetIdx2[fj]]*SubJet_mass[FatJet_subJetIdx2[fj]]);
            h_ptheta_sj->Fill(std::abs(E1-E2)/FatJet_pt[fj]);*/

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
    std::vector<TH1*> pt_nonp;
    pt_nonp.push_back(h_GFJ_pt);
    pt_nonp.push_back(h_FJ_pt);
    Makecanva2x1(c_pt_tot_nonp, pt_nonp, kViolet-4);

    c_pt_tot_nonp->SaveAs("plots_pt/c_pt_tot_nonp.pdf");
    


    /*ISTOGRAMMI PER ETA*/
    TCanvas *c_eta_tot_nonp = new TCanvas("c_eta_tot_nonp", "Consituents #eta: pseudorapirity", 2000, 1000);
    std::vector<TH1*> eta_nonp;
    eta_nonp.push_back(h_GFJ_eta);
    eta_nonp.push_back(h_FJ_eta);
    Makecanva2x1(c_eta_tot_nonp, eta_nonp, kViolet-4);

    c_eta_tot_nonp->SaveAs("plots_eta/c_eta_tot_nonp.pdf");


    /*ISTOGRAMMI PER PHI*/
    TCanvas *c_phi_tot_nonp = new TCanvas("c_phi_tot_nonp", "Consituents #phi", 2000, 1000);
    std::vector<TH1*> phi_nonp;
    phi_nonp.push_back(h_GFJ_phi);
    phi_nonp.push_back(h_FJ_phi);
    Makecanva2x1(c_phi_tot_nonp, phi_nonp, kViolet-4);
    c_phi_tot_nonp->SaveAs("plots_phi/c_phi_tot_nonp.pdf");
    


    
    /*ISTOGRAMMI PER LE MASSE*/
    TCanvas *c_mass_tot_nonp = new TCanvas("c_mass_tot_nonp", "Consituents Mass", 2000, 1100);
    std::vector<TH1*> mass_nonp;
    mass_nonp.push_back(h_GFJ_mass);
    mass_nonp.push_back(h_FJ_mass);
    Makecanva2x1(c_mass_tot_nonp, mass_nonp, kViolet-4);

    c_mass_tot_nonp->SaveAs("plots_m/c_mass_tot_nonp.pdf");






    /*ISTOGRAMMI TAU 1,2, 3 */
    TCanvas *c_tau_nonp = new TCanvas("c_tau_nonp", "FatJet #tau", 3000, 1100);
    std::vector<TH1*> tau_nonp;
    tau_nonp.push_back(h_tau1);
    tau_nonp.push_back(h_tau2);
    tau_nonp.push_back(h_tau3);
    Makecanva3x1(c_tau_nonp, tau_nonp, kViolet-4);
    c_tau_nonp->SaveAs("plots_tau/c_tau_nonp.pdf");



    /*ISTOGRAMMI Tau21, Tau32 confronto con eventi "W" */
    TCanvas *c_tau_ratio_nonp = new TCanvas("c_tau_ratio_nonp", "FatJet #tau ratios", 2000, 1800);
    std::vector<TH1*> tau_ratio_nonp;
    tau_ratio_nonp.push_back(h_tau21);
    tau_ratio_nonp.push_back(h_tau32);
    tau_ratio_nonp.push_back(h_tau21m);
    tau_ratio_nonp.push_back(h_tau32m);
    Makecanva2x2(c_tau_ratio_nonp, tau_ratio_nonp, kViolet-4);
    c_tau_ratio_nonp->SaveAs("plots_tau/c_tau_ratio_nonp.pdf");


    /*TCanvas *c_zg_pt_nonp = new TCanvas("c_zg_pt_nonp", "z_{g}", 2000, 1800);
    std::vector<TH1*> zg_pt_nonp;
    zg_pt_nonp.push_back(h_zg_pt200);
    zg_pt_nonp.push_back(h_zg_pt230);
    zg_pt_nonp.push_back(h_zg_pt270);
    zg_pt_nonp.push_back(h_zg_pt300);
    Makecanva2x2(c_zg_pt_nonp, zg_pt_nonp, kViolet-4);
    c_zg_pt_nonp->SaveAs("plots_z1/c_zg_pt_nonp.pdf");
    

    
    
    
    
    TCanvas *c_ptheta_pt_nonp = new TCanvas("c_ptheta_pt_nonp", "p_{#theta}", 2000, 1800);
    std::vector<TH1*> ptheta_pt_nonp;
    ptheta_pt_nonp.push_back(h_ptheta_pt200);
    ptheta_pt_nonp.push_back(h_ptheta_pt230);
    ptheta_pt_nonp.push_back(h_ptheta_pt270);
    ptheta_pt_nonp.push_back(h_ptheta_pt300);
    Makecanva2x2(c_ptheta_pt_nonp, ptheta_pt_nonp, kViolet-4);
    c_ptheta_pt_nonp->SaveAs("plots_ptheta/c_ptheta_pt_nonp.pdf");*/





    //ISTOGRAMMI per SubJet
//eta
    /*TCanvas *c_SJ_etad = new TCanvas("c_SJ_etad", "SubJet #eta difference", 800, 600);

    h_SJ_etad->Scale(1.0 / h_SJ_etad->Integral());
    h_SJ_etad->Draw("hist");
    h_SJ_etad->SetLineColor(kViolet-4);
    h_SJ_etad->SetFillColor(kViolet-4);
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
    TCanvas *c_SJ_phid = new TCanvas("c_SJ_phid", "SubJet phi difference", 800, 600);
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


    TCanvas *c_sj_mass_nonp = new TCanvas("c_sj_mass_nonp", "Consituents Mass", 1200, 600);
    std::vector<TH1*> sj_mass_nonp;
    sj_mass_nonp.push_back(h_SJ_mass1);
    sj_mass_nonp.push_back(h_SJ_mass2);
    Makecanva2x1(c_mass_tot_nonp, mass_nonp, kViolet+6);
    c_sj_mass_nonp->SaveAs("plots_sj/c_sj_mass_nonp.pdf");

    TCanvas *c_sj_pt_nonp = new TCanvas("c_sj_pt_nonp", "Consituents pt", 1200, 600);
    std::vector<TH1*> sj_pt_nonp;
    sj_pt_nonp.push_back(h_SJ_pts1);
    sj_pt_nonp.push_back(h_SJ_pts2);
    Makecanva2x1(c_sj_pt_nonp, sj_pt_nonp, kViolet+6);
    c_sj_pt_nonp->SaveAs("plots_sj/c_sj_pt_nonp.pdf");

    TCanvas *c_sj_eta_nonp = new TCanvas("c_sj_eta_nonp", "Consituents #eta", 1200, 600);
    std::vector<TH1*> sj_eta_nonp;
    sj_eta_nonp.push_back(h_SJ_eta1);
    sj_eta_nonp.push_back(h_SJ_eta2);
    Makecanva2x1(c_sj_eta_nonp, sj_eta_nonp, kViolet+6);
    c_sj_eta_nonp->SaveAs("plots_sj/c_sj_eta_nonp.pdf");

    TCanvas *c_sj_phi_nonp = new TCanvas("c_sj_phi_nonp", "Consituents #phi", 1200, 600);
    std::vector<TH1*> sj_phi_nonp;
    sj_phi_nonp.push_back(h_SJ_phi1);
    sj_phi_nonp.push_back(h_SJ_phi2);
    Makecanva2x1(c_sj_phi_nonp, sj_eta_nonp, kViolet+6);
    c_sj_phi_nonp->SaveAs("plots_sj/c_sj_phi_nonp.pdf");





    gStyle->SetOptStat(0);


    //SCATTER PLOT TAU - MASSA SD
    TCanvas *c_tau_msd_nonp = new TCanvas("c_tau_msd_nonp","Tau vs msd", 2200,1000);
    c_tau_msd_nonp->Divide(2,1);

    c_tau_msd_nonp->cd(1);
    gPad->SetRightMargin(0.15);
    h_tau21_msd->Scale(1.0/ h_tau21_msd->Integral());
    h_tau21_msd->Draw("COLZ"); 
    h_tau21_msd->GetXaxis()->SetTitleSize(0.05);

    c_tau_msd_nonp->cd(2);
    gPad->SetRightMargin(0.15);
    h_tau32_msd->Scale(1.0/ h_tau32_msd->Integral());
    h_tau32_msd->Draw("COLZ"); 
    h_tau32_msd->GetXaxis()->SetTitleSize(0.05);

    c_tau_msd_nonp->SaveAs("plots_tau/c_tau_msd_nonp.pdf");



    /*Soft Drop*/
    TCanvas *c_FJ_msd_nonp = new TCanvas("c_FJ_msd_nonp", "FatJet SoftDrop mass", 800, 600);
    h_FJ_mass->SetTitle("SoftDrop vs FatJet mass (CMS private work)");
    std::vector<TH1*> FJ_msd_nonp;
    FJ_msd_nonp.push_back(h_FJ_mass);
    FJ_msd_nonp.push_back(h_FJ_msd);
    Makecanva(c_FJ_msd_nonp, FJ_msd_nonp, kViolet-4);
    TLegend *leg9 = new TLegend(0.6, 0.65, 0.8, 0.8);
    leg9->AddEntry(h_FJ_msd, "SoftDrop", "f");
    leg9->AddEntry(h_FJ_mass, "FatJet", "f");
    leg9->Draw();
    c_FJ_msd_nonp->SaveAs("plots_m/c_FJ_msd_nonp.pdf");


    //zg e ptheta
    /*TCanvas *c_zg_nonp = new TCanvas("c_zg_nonp", "z_{g}", 800, 600);
    std::vector<TH1*> zg_nonp;
    zg_nonp.push_back(h_zg_sj);
    zg_nonp.push_back(h_zg_sjsd);
    Makecanva(c_zg_nonp, zg_nonp, kViolet-4);
    TLegend *leg21 = new TLegend(0.6, 0.65, 0.8, 0.85);
    leg21->AddEntry(h_zg_sjsd, "SD mass cut", "l");
    leg21->AddEntry(h_zg_sj, "No cuts", "l");
    leg21->Draw();
    c_zg_nonp->SaveAs("plots_z1/c_zg_nonp.pdf");




    TCanvas *c_ptheta_nonp = new TCanvas("c_ptheta_nonp", "p_{#theta}", 800, 600);
    std::vector<TH1*> ptheta_nonp;
    ptheta_nonp.push_back(h_ptheta_sj);
    ptheta_nonp.push_back(h_ptheta_sjsd);
    Makecanva(c_ptheta_nonp, ptheta_nonp, kViolet-4);
    TLegend *leg1 = new TLegend(0.6, 0.65, 0.8, 0.85);
    leg1->AddEntry(h_ptheta_sjsd, "SD mass cut", "l");
    leg1->AddEntry(h_ptheta_sj, "No cuts", "l");
    leg1->Draw();
    c_ptheta_nonp->SaveAs("plots_ptheta/c_ptheta_nonp.pdf");*/




    f->Write();
    f->Close();


    gROOT->ProcessLine(".q");



}
