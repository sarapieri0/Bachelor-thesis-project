#include "samples.h"
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

void jet_plots2(
    TString output="jets2.root"
){
    TFile *f=new TFile(output, "recreate");

    Int_t maxEvents = 2000;
    std::string tagLL = "ssWWLL";
    std::string tagTT = "ssWWTT";
    

//pt    
    TH1F *h_FJ_ptLL = new TH1F(("h_FatJet_pt_" + tagLL).c_str(), "FatJet p_{T} LL;p_{T} [GeV];Events", 50, 170, 1000);
    TH1F *h_FJ_ptTT = new TH1F(("h_FatJet_pt_" + tagTT).c_str(), "FatJet p_{T} TT;p_{T} [GeV];Events", 50, 170, 1000);
    
    TH1F *h_GFJ_ptLL = new TH1F(("h_GenJetAK8_pt_" + tagLL).c_str(), "GenJetAK8 p_{T} LL;p_{T} [GeV];Events", 50, 150, 1000);
    TH1F *h_GFJ_ptTT = new TH1F(("h_GenJetAK8_pt_" + tagTT).c_str(), "GenJetAK8 p_{T} TT;p_{T} [GeV];Events", 50, 150, 1000);

    
    
//eta    
    TH1F *h_FJ_etaLL = new TH1F(("h_FatJet_eta_" + tagLL).c_str(), "FatJet #eta LL; #eta; Events", 50,-8, 8);
    TH1F *h_FJ_etaTT = new TH1F(("h_FatJet_eta_" + tagTT).c_str(), "FatJet #eta TT; #eta; Events", 50, -8, 8);
    
    TH1F *h_GFJ_etaLL = new TH1F(("h_GenJetAK8_eta_" + tagLL).c_str(), "GenJetAK8 #eta LL;#eta; Events", 50, -8, 8);
    TH1F *h_GFJ_etaTT = new TH1F(("h_GenJetAK8_eta_" + tagTT).c_str(), "GenJetAK8 #eta TT;#eta; Events", 50, -8, 8);

//phi
    TH1F *h_FJ_phiLL = new TH1F(("h_FatJet_phi_" + tagLL).c_str(), "FatJet #phi LL; #phi; Events", 40, -4,4);
    TH1F *h_FJ_phiTT = new TH1F(("h_FatJet_phi_" + tagTT).c_str(), "FatJet #phi TT; #phi; Events", 40, -4, 4);

    TH1F *h_GFJ_phiLL = new TH1F(("h_GenJetAK8_phi_" + tagLL).c_str(), "GenJetAK8 #phi LL; #phi; Events", 40, -4, 4);
    TH1F *h_GFJ_phiTT = new TH1F(("h_GenJetAK8_phi_" + tagTT).c_str(), "GenJetAK8 #phi TT; #phi; Events", 40, -4, 4);


//mass
    TH1F *h_FJ_massLL = new TH1F(("h_FJ_mass_" + tagLL).c_str(), "FatJet mass LL [GeV]", 50,0,300);
    TH1F *h_FJ_massTT = new TH1F(("h_FJ_mass_" + tagTT).c_str(), "FatJet mass TT [GeV]", 50,0,300);

    TH1F *h_GFJ_massLL = new TH1F(("h_GFJ_mass_" + tagLL).c_str(), "GenAK8 FatJet mass LL;mass [GeV];", 50, 0, 300);
    TH1F *h_GFJ_massTT = new TH1F(("h_GFJ_mass_" + tagTT).c_str(), "GenAK8 FatJet mass TT;mass [GeV];", 50, 0, 300);

//mass softdrop
    TH1F *h_FJ_msdLL = new TH1F(("h_FJ_msd_" + tagLL).c_str(), "FatJet SoftDrop mass LL", 50, 0, 250);
    TH1F *h_FJ_msdTT = new TH1F(("h_FJ_msd_" + tagTT).c_str(), "FatJet SoftDrop mass TT", 50, 0, 250);


//tau
    TH1F *h_tau1LL =new TH1F(("h_tau1_" + tagLL).c_str(), "FatJet tau1 LL", 30, 0, 0.6);
    TH1F *h_tau2LL =new TH1F(("h_tau2_" + tagLL).c_str(), "FatJet tau2 LL", 30, 0, 0.6);
    TH1F *h_tau3LL= new TH1F(("h_tau3_" + tagLL).c_str(), "FatJet tau3 LL", 30, 0, 0.6);
    TH1F *h_tau1TT =new TH1F(("h_tau1_" + tagTT).c_str(), "FatJet tau1 TT", 30, 0, 0.6);
    TH1F *h_tau2TT =new TH1F(("h_tau2_" + tagTT).c_str(), "FatJet tau2 TT", 30, 0, 0.6);
    TH1F *h_tau3TT= new TH1F(("h_tau3_" + tagTT).c_str(), "FatJet tau3 TT", 30, 0, 0.6);

    TH1F *h_tau21mLL= new TH1F(("h_tau21m_" + tagLL).c_str(), "#tau_{21} 60< m_sd <100 LL", 30, 0, 1);
    TH1F *h_tau21LL= new TH1F(("h_tau21_" + tagLL).c_str(), "#tau_{21} LL", 30, 0, 1);
    TH1F *h_tau32mLL= new TH1F(("h_tau32m_" + tagLL).c_str(), "#tau_{32} 60< m_sd <100 LL", 30, 0, 1);
    TH1F *h_tau32LL= new TH1F(("h_tau32_" + tagLL).c_str(), "#tau_{32} LL", 30, 0, 1);
    TH1F *h_tau21mTT= new TH1F(("h_tau21m_" + tagTT).c_str(), "#tau_{21} 60< m_sd <100 TT", 30, 0, 1);
    TH1F *h_tau21TT= new TH1F(("h_tau21_" + tagTT).c_str(), "#tau_{21} TT", 30, 0, 1);
    TH1F *h_tau32mTT= new TH1F(("h_tau32m_" + tagTT).c_str(), "#tau_{32} 60< m_sd <100 TT", 30, 0, 1);
    TH1F *h_tau32TT= new TH1F(("h_tau32_" + tagTT).c_str(), "#tau_{32} TT", 30, 0, 1);

    TH2F *h_tau21_msdLL = new TH2F(("h_msd_tau21_" + tagLL).c_str(), "SD mass vs  #tau_{21} LL; #tau_{21}; SoftDrop mass", 50, 0, 1,    30, 50, 300);
    TH2F *h_tau21_msdTT = new TH2F(("h_msd_tau21_" + tagTT).c_str(), "SD mass vs  #tau_{21} TT; #tau_{21}; SoftDrop mass", 50, 0, 1,    30, 50, 300);

    TH2F *h_tau32_msdLL = new TH2F(("h_msd_tau32_" + tagLL).c_str(), "SD mass vs #tau_{32} LL; #tau_{32}; SoftDrop mass", 50, 0, 1,    30, 50, 300);
    TH2F *h_tau32_msdTT = new TH2F(("h_msd_tau32_" + tagTT).c_str(), "SD mass vs #tau_{32} TT; #tau_{32}; SoftDrop mass", 50, 0, 1,    30, 50, 300);


//subjet
    TH1F *h_SJ_etad_sdLL = new TH1F(("h_SJ_etad_sd_" + tagLL).c_str(), "SubJet #eta difference LL; #eta; Events", 20,-1,1);
    TH1F *h_SJ_etad_sdTT = new TH1F(("h_SJ_etad_sd_" + tagTT).c_str(), "SubJet #eta difference TT; #eta; Events", 20,-1,1);
    TH1F *h_SJ_etadLL = new TH1F(("h_SJ_etad_" + tagLL).c_str(), "SubJet #eta difference LL; #eta; Events", 30,-1,1);
    TH1F *h_SJ_etadTT = new TH1F(("h_SJ_etad_" + tagTT).c_str(), "SubJet #eta difference TT; #eta; Events", 30,-1,1);
    
    TH1F *h_SJ_phid_sdLL = new TH1F(("h_SJ_phid_sd_" + tagLL).c_str(), "SubJet #phi difference LL; #phi; Events", 20,-1,1);
    TH1F *h_SJ_phid_sdTT = new TH1F(("h_SJ_phid_sd_" + tagTT).c_str(), "SubJet #phi difference TT; #phi; Events", 20,-1,1);
    TH1F *h_SJ_phidLL = new TH1F(("h_SJ_phid_" + tagLL).c_str(), "SubJet #phi difference LL; #phi; Events", 30,-1,1);
    TH1F *h_SJ_phidTT = new TH1F(("h_SJ_phid_" + tagTT).c_str(), "SubJet #phi difference TT; #phi; Events", 30,-1,1);

    TH1F *h_SJ_ptd_sdLL = new TH1F(("h_SJ_ptd_sd_" + tagLL).c_str(), "SubJet p_{t} difference LL; p_{t}; Events", 30,150,1000);
    TH1F *h_SJ_ptd_sdTT = new TH1F(("h_SJ_ptd_sd_" + tagTT).c_str(), "SubJet p_{t} difference TT; p_{t}; Events", 30,150,1000);
    TH1F *h_SJ_ptdLL = new TH1F(("h_SJ_ptd_" + tagLL).c_str(), "SubJet p_{t} difference LL; p_{t}; Events", 30,150,1000);
    TH1F *h_SJ_ptdTT = new TH1F(("h_SJ_ptd_" + tagTT).c_str(), "SubJet p_{t} difference p_{t}; phi; Events", 30,150,1000);

    TH1F *h_SJ_pts1LL = new TH1F(("h_SJ_pts1_" + tagLL).c_str(), "SubJet #1 p_{t} LL; p_{t};Events", 30, 0, 600);
    TH1F *h_SJ_pts2LL = new TH1F(("h_SJ_pts2_" + tagLL).c_str(), "SubJet #2 p_{t} LL; p_{t};Events", 30, 0, 600);
    TH1F *h_SJ_pts1TT = new TH1F(("h_SJ_pts1_" + tagTT).c_str(), "SubJet #1 p_{t} TT; p_{t};Events", 30, 0, 600);
    TH1F *h_SJ_pts2TT = new TH1F(("h_SJ_pts2_" + tagTT).c_str(), "SubJet #2 p_{t} TT; p_{t};Events", 30, 0, 600);
    
//zg1
    TH1F *h_SJ_z1mtLL= new TH1F(("h_SJ_z1mtLL" + tagLL).c_str(), "Subjet z_{g} LL; z_{g}; Events", 30, 1, 10);
    TH1F *h_SJ_z1mtTT= new TH1F(("h_SJ_z1mtTT" + tagTT).c_str(), "Subjet z_{g} TT; z_{g}; Events", 30, 1, 10);
    TH1F *h_SJ_z1mLL= new TH1F(("h_SJ_z1mLL" + tagLL).c_str(), "Subjet z_{g} LL; z_{g}; Events", 30, 1, 10);
    TH1F *h_SJ_z1mTT= new TH1F(("h_SJ_z1mTT" + tagTT).c_str(), "Subjet z_{g} TT; z_{g}; Events", 30, 1, 10);
    TH1F *h_SJ_z1pt1LL= new TH1F(("h_SJ_z1pt1LL" + tagLL).c_str(), "Subjet z_{g} LL; z_{g}; Events", 30, 1, 10);
    TH1F *h_SJ_z1pt1TT= new TH1F(("h_SJ_z1pt1TT" + tagTT).c_str(), "Subjet z_{g} TT; z_{g}; Events", 30, 1, 10);
    TH1F *h_SJ_z1pt2LL= new TH1F(("h_SJ_z1pt2LL" + tagLL).c_str(), "Subjet z_{g} LL; z_{g}; Events", 30, 1, 10);
    TH1F *h_SJ_z1pt2TT= new TH1F(("h_SJ_z1pt2TT" + tagTT).c_str(), "Subjet z_{g} TT; z_{g}; Events", 30, 1, 10);
    TH1F *h_SJ_z1pt3LL= new TH1F(("h_SJ_z1pt3LL" + tagLL).c_str(), "Subjet z_{g} LL; z_{g}; Events", 30, 1, 10);
    TH1F *h_SJ_z1pt3TT= new TH1F(("h_SJ_z1pt3TT" + tagTT).c_str(), "Subjet z_{g} TT; z_{g}; Events", 30, 1, 10);
    TH1F *h_SJ_z1pt4LL= new TH1F(("h_SJ_z1pt4LL" + tagLL).c_str(), "Subjet z_{g} LL; z_{g}; Events", 30, 1, 10);
    TH1F *h_SJ_z1pt4TT= new TH1F(("h_SJ_z1pt4TT" + tagTT).c_str(), "Subjet z_{g} TT; z_{g}; Events", 30, 1, 10);
//zg2
    TH1F *h_SJ_z2mtLL= new TH1F(("h_SJ_z2mtLL" + tagLL).c_str(), "Subjet z_{g} LL; z_{g}'; Events", 30, 0, 0.5);
    TH1F *h_SJ_z2mtTT= new TH1F(("h_SJ_z2mtTT" + tagTT).c_str(), "Subjet z_{g} TT; z_{g}'; Events", 30, 0, 0.5);
    TH1F *h_SJ_z2mLL= new TH1F(("h_SJ_z2mLL" + tagLL).c_str(), "Subjet z_{g} LL; z_{g}'; Events", 30, 0, 0.5);
    TH1F *h_SJ_z2mTT= new TH1F(("h_SJ_z2mTT" + tagTT).c_str(), "Subjet z_{g} TT; z_{g}'; Events", 30, 0, 0.5);

    TH1F *h_SJ_z2pt1LL= new TH1F(("h_SJ_z2pt1LL" + tagLL).c_str(), "Subjet z_{g} LL; z_{g}'; Events", 30, 0, 0.5);
    TH1F *h_SJ_z2pt1TT= new TH1F(("h_SJ_z2pt1TT" + tagTT).c_str(), "Subjet z_{g} TT; z_{g}'; Events", 30, 0, 0.5);
    TH1F *h_SJ_z2pt2LL= new TH1F(("h_SJ_z2pt2LL" + tagLL).c_str(), "Subjet z_{g} LL; z_{g}'; Events", 30, 0, 0.5);
    TH1F *h_SJ_z2pt2TT= new TH1F(("h_SJ_z2pt2TT" + tagTT).c_str(), "Subjet z_{g} TT; z_{g}'; Events", 30, 0, 0.5);
    TH1F *h_SJ_z2pt3LL= new TH1F(("h_SJ_z2pt3LL" + tagLL).c_str(), "Subjet z_{g} LL; z_{g}'; Events", 30, 0, 0.5);
    TH1F *h_SJ_z2pt3TT= new TH1F(("h_SJ_z2pt3TT" + tagTT).c_str(), "Subjet z_{g} TT; z_{g}'; Events", 30, 0, 0.5);
    TH1F *h_SJ_z2pt4LL= new TH1F(("h_SJ_z2pt4LL" + tagLL).c_str(), "Subjet z_{g} LL; z_{g}'; Events", 30, 0, 0.5);
    TH1F *h_SJ_z2pt4TT= new TH1F(("h_SJ_z2pt4TT" + tagTT).c_str(), "Subjet z_{g} TT; z_{g}'; Events", 30, 0, 0.5);

//ptheta
    TH1F *h_SJ_pemtLL= new TH1F(("h_SJ_pemtLL" + tagLL).c_str(), "Subjet p_{#theta} LL; p_{#theta}; Events", 30, 0, 1);
    TH1F *h_SJ_pemtTT= new TH1F(("h_SJ_pemtTT" + tagTT).c_str(), "Subjet p_{#theta} TT; p_{#theta}; Events", 30, 0, 1);
    TH1F *h_SJ_pemLL= new TH1F(("h_SJ_pemLL" + tagLL).c_str(), "Subjet p_{#theta} LL; p_{#theta}; Events", 30, 0, 1);
    TH1F *h_SJ_pemTT= new TH1F(("h_SJ_pemTT" + tagTT).c_str(), "Subjet p_{#theta} TT; p_{#theta}; Events", 30, 0, 1);
    TH1F *h_SJ_pept1LL= new TH1F(("h_SJ_pept1LL" + tagLL).c_str(), "Subjet p_{#theta} LL; p_{#theta}; Events", 30, 1, 1);
    TH1F *h_SJ_pept1TT= new TH1F(("h_SJ_pept1TT" + tagTT).c_str(), "Subjet p_{#theta} TT; p_{#theta}; Events", 30, 1, 1);
    TH1F *h_SJ_pept2LL= new TH1F(("h_SJ_pept2LL" + tagLL).c_str(), "Subjet p_{#theta} LL; p_{#theta}; Events", 30, 1, 1);
    TH1F *h_SJ_pept2TT= new TH1F(("h_SJ_pept2TT" + tagTT).c_str(), "Subjet p_{#theta} TT; p_{#theta}; Events", 30, 1, 1);
    TH1F *h_SJ_pept3LL= new TH1F(("h_SJ_pept3LL" + tagLL).c_str(), "Subjet p_{#theta} LL; p_{#theta}; Events", 30, 1, 1);
    TH1F *h_SJ_pept3TT= new TH1F(("h_SJ_pept3TT" + tagTT).c_str(), "Subjet p_{#theta} TT; p_{#theta}; Events", 30, 1, 1);
    TH1F *h_SJ_pept4LL= new TH1F(("h_SJ_pept4LL" + tagLL).c_str(), "Subjet p_{#theta} LL; p_{#theta}; Events", 30, 1, 1);
    TH1F *h_SJ_pept4TT= new TH1F(("h_SJ_pept4TT" + tagTT).c_str(), "Subjet p_{#theta} TT; p_{#theta}; Events", 30, 1, 1);



    
 


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

    //Creo due chain, una per LL e una per TT
    TChain *chainLL = new TChain("Events"); 
    for(auto sample : samples_ssWWLL){
        chainLL->Add(sample.c_str()); 
    }

    TChain *chainTT = new TChain("Events");
    for(auto sample : samples_ssWWTT){
        chainTT->Add(sample.c_str());
    }
    




/*BRANCH LL*/
    chainLL->SetBranchAddress("nFatJet", &nFatJet);
    chainLL->SetBranchAddress("FatJet_mass", FatJet_mass);
    chainLL->SetBranchAddress("FatJet_msoftdrop", FatJet_msoftdrop);
    chainLL->SetBranchAddress("FatJet_pt", FatJet_pt);
    chainLL->SetBranchAddress("FatJet_eta", FatJet_eta);
    chainLL->SetBranchAddress("FatJet_phi", FatJet_phi);
    chainLL->SetBranchAddress("FatJet_tau1", FatJet_tau1);
    chainLL->SetBranchAddress("FatJet_tau2", FatJet_tau2);
    chainLL->SetBranchAddress("FatJet_tau3", FatJet_tau3);
    chainLL->SetBranchAddress("FatJet_subJetIdx1", FatJet_subJetIdx1);
    chainLL->SetBranchAddress("FatJet_subJetIdx2", FatJet_subJetIdx2);


    
    chainLL->SetBranchAddress("nJet", &nJet);
    chainLL->SetBranchAddress("Jet_mass", Jet_mass);
    chainLL->SetBranchAddress("Jet_pt", Jet_pt);
    chainLL->SetBranchAddress("Jet_eta", Jet_eta);
    chainLL->SetBranchAddress("Jet_phi", Jet_phi);


    chainLL->SetBranchAddress("nSubJet", &nSubJet);
    chainLL->SetBranchAddress("SubJet_pt", SubJet_pt);
    chainLL->SetBranchAddress("SubJet_eta", SubJet_eta);
    chainLL->SetBranchAddress("SubJet_phi", SubJet_phi);
    chainLL->SetBranchAddress("SubJet_mass", SubJet_mass);    

    chainLL->SetBranchAddress("nGenJetAK8", &nGenFatJet);
    chainLL->SetBranchAddress("GenJetAK8_mass", GenFatJet_mass);
    chainLL->SetBranchAddress("GenJetAK8_pt", GenFatJet_pt);
    chainLL->SetBranchAddress("GenJetAK8_eta", GenFatJet_eta);
    chainLL->SetBranchAddress("GenJetAK8_phi", GenFatJet_phi);

 
    chainLL->SetBranchAddress("nGenJet", &nGenJet);
    chainLL->SetBranchAddress("GenJet_mass", GenJet_mass);
    chainLL->SetBranchAddress("GenJet_pt", GenJet_pt);
    chainLL->SetBranchAddress("GenJet_eta", GenJet_eta);
    chainLL->SetBranchAddress("GenJet_phi", GenJet_phi);


/*BRANCH TT*/
    chainTT->SetBranchAddress("nFatJet", &nFatJet);
    chainTT->SetBranchAddress("FatJet_mass", FatJet_mass);
    chainTT->SetBranchAddress("FatJet_msoftdrop", FatJet_msoftdrop);
    chainTT->SetBranchAddress("FatJet_pt", FatJet_pt);
    chainTT->SetBranchAddress("FatJet_eta", FatJet_eta);
    chainTT->SetBranchAddress("FatJet_phi", FatJet_phi);
    chainTT->SetBranchAddress("FatJet_tau1", FatJet_tau1);
    chainTT->SetBranchAddress("FatJet_tau2", FatJet_tau2);
    chainTT->SetBranchAddress("FatJet_tau3", FatJet_tau3);
    chainTT->SetBranchAddress("FatJet_subJetIdx1", FatJet_subJetIdx1);
    chainTT->SetBranchAddress("FatJet_subJetIdx2", FatJet_subJetIdx2);
    
    chainTT->SetBranchAddress("nJet", &nJet);
    chainTT->SetBranchAddress("Jet_mass", Jet_mass);
    chainTT->SetBranchAddress("Jet_pt", Jet_pt);
    chainTT->SetBranchAddress("Jet_eta", Jet_eta);
    chainTT->SetBranchAddress("Jet_phi", Jet_phi);

    chainTT->SetBranchAddress("nSubJet", &nSubJet);
    chainTT->SetBranchAddress("SubJet_pt", SubJet_pt);
    chainTT->SetBranchAddress("SubJet_eta", SubJet_eta);
    chainTT->SetBranchAddress("SubJet_phi", SubJet_phi);    
    chainTT->SetBranchAddress("SubJet_mass", SubJet_mass);    

    chainTT->SetBranchAddress("nGenJetAK8", &nGenFatJet);
    chainTT->SetBranchAddress("GenJetAK8_mass", GenFatJet_mass);
    chainTT->SetBranchAddress("GenJetAK8_pt", GenFatJet_pt);
    chainTT->SetBranchAddress("GenJetAK8_eta", GenFatJet_eta);
    chainTT->SetBranchAddress("GenJetAK8_phi", GenFatJet_phi);

    chainTT->SetBranchAddress("nGenJet", &nGenJet);
    chainTT->SetBranchAddress("GenJet_mass", GenJet_mass);
    chainTT->SetBranchAddress("GenJet_pt", GenJet_pt);
    chainTT->SetBranchAddress("GenJet_eta", GenJet_eta);
    chainTT->SetBranchAddress("GenJet_phi", GenJet_phi);






    for(int i=0; i< int(chainLL->GetEntries()); i++){
        chainLL->GetEntry(i);
        if(i%5000==0) std::cout<<"Processing event "<<i<<std::endl;
        //if(i==maxEvents) break;
        for(Int_t fj=0; fj<nFatJet; fj++){
            h_FJ_massLL->Fill(FatJet_mass[fj]);
            h_FJ_msdLL->Fill(FatJet_msoftdrop[fj]);
            h_FJ_ptLL->Fill(FatJet_pt[fj]);
            h_FJ_etaLL->Fill(FatJet_eta[fj]);
            h_FJ_phiLL->Fill(FatJet_phi[fj]);

            h_tau1LL->Fill(FatJet_tau1[fj]);
            h_tau2LL->Fill(FatJet_tau2[fj]);
            h_tau3LL->Fill(FatJet_tau3[fj]);
            h_tau21LL->Fill(FatJet_tau2[fj]/FatJet_tau1[fj]);
            h_tau32LL->Fill(FatJet_tau3[fj]/FatJet_tau2[fj]);
            if(FatJet_msoftdrop[fj]>60 && FatJet_msoftdrop[fj]<100)
            {
                h_tau21mLL->Fill(FatJet_tau2[fj]/FatJet_tau1[fj]);
                h_tau32mLL->Fill(FatJet_tau3[fj]/FatJet_tau2[fj]);

                h_SJ_etad_sdLL->Fill(SubJet_eta[FatJet_subJetIdx1[fj]] - SubJet_eta[FatJet_subJetIdx2[fj]]);
                h_SJ_phid_sdLL->Fill(SubJet_phi[FatJet_subJetIdx1[fj]] - SubJet_phi[FatJet_subJetIdx2[fj]]);
                h_SJ_ptd_sdLL->Fill(SubJet_pt[FatJet_subJetIdx1[fj]] - SubJet_pt[FatJet_subJetIdx2[fj]]);
                h_SJ_pts1LL->Fill(SubJet_pt[FatJet_subJetIdx1[fj]]);
                h_SJ_pts2LL->Fill(SubJet_pt[FatJet_subJetIdx2[fj]]);
                
                double max_pt = std::max(SubJet_pt[FatJet_subJetIdx1[fj]], SubJet_pt[FatJet_subJetIdx2[fj]]);
                double min_pt = std::min(SubJet_pt[FatJet_subJetIdx1[fj]], SubJet_pt[FatJet_subJetIdx2[fj]]);
                h_SJ_z1mLL->Fill(max_pt/min_pt);
                h_SJ_z2mLL->Fill(min_pt/(SubJet_pt[FatJet_subJetIdx1[fj]] + SubJet_pt[FatJet_subJetIdx2[fj]]));

                double E1, E2;
                E1 = std::sqrt(SubJet_pt[FatJet_subJetIdx1[fj]]*SubJet_pt[FatJet_subJetIdx1[fj]] + SubJet_mass[FatJet_subJetIdx1[fj]]*SubJet_mass[FatJet_subJetIdx1[fj]]);
                E2 = std::sqrt(SubJet_pt[FatJet_subJetIdx2[fj]]*SubJet_pt[FatJet_subJetIdx2[fj]] + SubJet_mass[FatJet_subJetIdx2[fj]]*SubJet_mass[FatJet_subJetIdx2[fj]]);
                h_SJ_pemLL->Fill(std::abs(E1-E2)/FatJet_pt[fj]);

                if(FatJet_tau2[fj]/FatJet_tau1[fj] < 0.45){
                    double max_pt = std::max(SubJet_pt[FatJet_subJetIdx1[fj]], SubJet_pt[FatJet_subJetIdx2[fj]]);
                    double min_pt = std::min(SubJet_pt[FatJet_subJetIdx1[fj]], SubJet_pt[FatJet_subJetIdx2[fj]]);
                    h_SJ_z1mtLL->Fill(max_pt/min_pt);
                    h_SJ_z2mtLL->Fill(min_pt/(SubJet_pt[FatJet_subJetIdx1[fj]] + SubJet_pt[FatJet_subJetIdx2[fj]]));

                    double E1, E2;
                    E1 = std::sqrt(SubJet_pt[FatJet_subJetIdx1[fj]]*SubJet_pt[FatJet_subJetIdx1[fj]] + SubJet_mass[FatJet_subJetIdx1[fj]]*SubJet_mass[FatJet_subJetIdx1[fj]]);
                    E2 = std::sqrt(SubJet_pt[FatJet_subJetIdx2[fj]]*SubJet_pt[FatJet_subJetIdx2[fj]] + SubJet_mass[FatJet_subJetIdx2[fj]]*SubJet_mass[FatJet_subJetIdx2[fj]]);
                    h_SJ_pemtLL->Fill(std::abs(E1-E2)/FatJet_pt[fj]);
                }

                if(FatJet_pt[fj]>220){
                    double max_pt = std::max(SubJet_pt[FatJet_subJetIdx1[fj]], SubJet_pt[FatJet_subJetIdx2[fj]]);
                    double min_pt = std::min(SubJet_pt[FatJet_subJetIdx1[fj]], SubJet_pt[FatJet_subJetIdx2[fj]]);
                    h_SJ_z1pt1LL->Fill(max_pt/min_pt);
                    h_SJ_z2pt1LL->Fill(min_pt/(SubJet_pt[FatJet_subJetIdx1[fj]] + SubJet_pt[FatJet_subJetIdx2[fj]]));
                    
                    double E1, E2;
                    E1 = std::sqrt(SubJet_pt[FatJet_subJetIdx1[fj]]*SubJet_pt[FatJet_subJetIdx1[fj]] + SubJet_mass[FatJet_subJetIdx1[fj]]*SubJet_mass[FatJet_subJetIdx1[fj]]);
                    E2 = std::sqrt(SubJet_pt[FatJet_subJetIdx2[fj]]*SubJet_pt[FatJet_subJetIdx2[fj]] + SubJet_mass[FatJet_subJetIdx2[fj]]*SubJet_mass[FatJet_subJetIdx2[fj]]);
                    h_SJ_pept1LL->Fill(std::abs(E1-E2)/FatJet_pt[fj]);
                }
                if(FatJet_pt[fj]>250){
                    double max_pt = std::max(SubJet_pt[FatJet_subJetIdx1[fj]], SubJet_pt[FatJet_subJetIdx2[fj]]);
                    double min_pt = std::min(SubJet_pt[FatJet_subJetIdx1[fj]], SubJet_pt[FatJet_subJetIdx2[fj]]);
                    h_SJ_z1pt2LL->Fill(max_pt/min_pt);
                    h_SJ_z2pt2LL->Fill(min_pt/(SubJet_pt[FatJet_subJetIdx1[fj]] + SubJet_pt[FatJet_subJetIdx2[fj]]));

                    double E1, E2;
                    E1 = std::sqrt(SubJet_pt[FatJet_subJetIdx1[fj]]*SubJet_pt[FatJet_subJetIdx1[fj]] + SubJet_mass[FatJet_subJetIdx1[fj]]*SubJet_mass[FatJet_subJetIdx1[fj]]);
                    E2 = std::sqrt(SubJet_pt[FatJet_subJetIdx2[fj]]*SubJet_pt[FatJet_subJetIdx2[fj]] + SubJet_mass[FatJet_subJetIdx2[fj]]*SubJet_mass[FatJet_subJetIdx2[fj]]);
                    h_SJ_pept2LL->Fill(std::abs(E1-E2)/FatJet_pt[fj]);
                }
                if(FatJet_pt[fj]>280){
                    double max_pt = std::max(SubJet_pt[FatJet_subJetIdx1[fj]], SubJet_pt[FatJet_subJetIdx2[fj]]);
                    double min_pt = std::min(SubJet_pt[FatJet_subJetIdx1[fj]], SubJet_pt[FatJet_subJetIdx2[fj]]);
                    h_SJ_z1pt3LL->Fill(max_pt/min_pt);
                    h_SJ_z2pt3LL->Fill(min_pt/(SubJet_pt[FatJet_subJetIdx1[fj]] + SubJet_pt[FatJet_subJetIdx2[fj]]));

                    double E1, E2;
                    E1 = std::sqrt(SubJet_pt[FatJet_subJetIdx1[fj]]*SubJet_pt[FatJet_subJetIdx1[fj]] + SubJet_mass[FatJet_subJetIdx1[fj]]*SubJet_mass[FatJet_subJetIdx1[fj]]);
                    E2 = std::sqrt(SubJet_pt[FatJet_subJetIdx2[fj]]*SubJet_pt[FatJet_subJetIdx2[fj]] + SubJet_mass[FatJet_subJetIdx2[fj]]*SubJet_mass[FatJet_subJetIdx2[fj]]);
                    h_SJ_pept3LL->Fill(std::abs(E1-E2)/FatJet_pt[fj]);
                }
                if(FatJet_pt[fj]>310){
                    double max_pt = std::max(SubJet_pt[FatJet_subJetIdx1[fj]], SubJet_pt[FatJet_subJetIdx2[fj]]);
                    double min_pt = std::min(SubJet_pt[FatJet_subJetIdx1[fj]], SubJet_pt[FatJet_subJetIdx2[fj]]);
                    h_SJ_z1pt4LL->Fill(max_pt/min_pt);
                    h_SJ_z2pt4LL->Fill(min_pt/(SubJet_pt[FatJet_subJetIdx1[fj]] + SubJet_pt[FatJet_subJetIdx2[fj]]));

                    double E1, E2;
                    E1 = std::sqrt(SubJet_pt[FatJet_subJetIdx1[fj]]*SubJet_pt[FatJet_subJetIdx1[fj]] + SubJet_mass[FatJet_subJetIdx1[fj]]*SubJet_mass[FatJet_subJetIdx1[fj]]);
                    E2 = std::sqrt(SubJet_pt[FatJet_subJetIdx2[fj]]*SubJet_pt[FatJet_subJetIdx2[fj]] + SubJet_mass[FatJet_subJetIdx2[fj]]*SubJet_mass[FatJet_subJetIdx2[fj]]);
                    h_SJ_pept4LL->Fill(std::abs(E1-E2)/FatJet_pt[fj]);
                }

            }
            
            h_tau21_msdLL->Fill((FatJet_tau2[fj]/FatJet_tau1[fj]), FatJet_msoftdrop[fj]);
            h_tau32_msdLL->Fill((FatJet_tau3[fj]/FatJet_tau2[fj]), FatJet_msoftdrop[fj]);

            h_SJ_etadLL->Fill(SubJet_eta[FatJet_subJetIdx1[fj]] - SubJet_eta[FatJet_subJetIdx2[fj]]);
            h_SJ_phidLL->Fill(SubJet_phi[FatJet_subJetIdx1[fj]] - SubJet_phi[FatJet_subJetIdx2[fj]]);
            h_SJ_ptdLL->Fill(SubJet_pt[FatJet_subJetIdx1[fj]] - SubJet_pt[FatJet_subJetIdx2[fj]]);
            

        }
        for(Int_t genfj=0; genfj<nGenFatJet; genfj++){
            h_GFJ_massLL->Fill(GenFatJet_mass[genfj]);
            h_GFJ_ptLL->Fill(GenFatJet_pt[genfj]);
            h_GFJ_etaLL->Fill(GenFatJet_eta[genfj]);
            h_GFJ_phiLL->Fill(GenFatJet_phi[genfj]);
        }
    }


    for(int i=0; i< int(chainTT->GetEntries()); i++){
        chainTT->GetEntry(i);
        if(i%5000==0) std::cout<<"Processing event "<<i<<std::endl;
        //if(i==maxEvents) break;
        for(Int_t fj=0; fj<nFatJet; fj++){
            h_FJ_massTT->Fill(FatJet_mass[fj]);
            h_FJ_msdTT->Fill(FatJet_msoftdrop[fj]);
            h_FJ_ptTT->Fill(FatJet_pt[fj]);
            h_FJ_etaTT->Fill(FatJet_eta[fj]);
            h_FJ_phiTT->Fill(FatJet_phi[fj]);

            h_tau1TT->Fill(FatJet_tau1[fj]);
            h_tau2TT->Fill(FatJet_tau2[fj]);
            h_tau3TT->Fill(FatJet_tau3[fj]);
            h_tau21TT->Fill(FatJet_tau2[fj]/FatJet_tau1[fj]);
            h_tau32TT->Fill(FatJet_tau3[fj]/FatJet_tau2[fj]);
            if(FatJet_msoftdrop[fj]>60 && FatJet_msoftdrop[fj]<100)
            {
                h_tau21mTT->Fill(FatJet_tau2[fj]/FatJet_tau1[fj]);
                h_tau32mTT->Fill(FatJet_tau3[fj]/FatJet_tau2[fj]);

                h_SJ_etad_sdTT->Fill(SubJet_eta[FatJet_subJetIdx1[fj]] - SubJet_eta[FatJet_subJetIdx2[fj]]);
                h_SJ_phid_sdTT->Fill(SubJet_phi[FatJet_subJetIdx1[fj]] - SubJet_phi[FatJet_subJetIdx2[fj]]);
                h_SJ_ptd_sdTT->Fill(SubJet_pt[FatJet_subJetIdx1[fj]] - SubJet_pt[FatJet_subJetIdx2[fj]]);

                h_SJ_pts1TT->Fill(SubJet_pt[FatJet_subJetIdx1[fj]]);
                h_SJ_pts2TT->Fill(SubJet_pt[FatJet_subJetIdx2[fj]]);

                double max_pt = std::max(SubJet_pt[FatJet_subJetIdx1[fj]], SubJet_pt[FatJet_subJetIdx2[fj]]);
                double min_pt = std::min(SubJet_pt[FatJet_subJetIdx1[fj]], SubJet_pt[FatJet_subJetIdx2[fj]]);                
                h_SJ_z1mTT->Fill(max_pt/min_pt);
                h_SJ_z2mTT->Fill(min_pt/(SubJet_pt[FatJet_subJetIdx1[fj]] + SubJet_pt[FatJet_subJetIdx2[fj]]));

                double E1, E2;
                E1 = std::sqrt(SubJet_pt[FatJet_subJetIdx1[fj]]*SubJet_pt[FatJet_subJetIdx1[fj]] + SubJet_mass[FatJet_subJetIdx1[fj]]*SubJet_mass[FatJet_subJetIdx1[fj]]);
                E2 = std::sqrt(SubJet_pt[FatJet_subJetIdx2[fj]]*SubJet_pt[FatJet_subJetIdx2[fj]] + SubJet_mass[FatJet_subJetIdx2[fj]]*SubJet_mass[FatJet_subJetIdx2[fj]]);
                h_SJ_pemTT->Fill(std::abs(E1-E2)/FatJet_pt[fj]);


                if(FatJet_tau2[fj]/FatJet_tau1[fj] < 0.45){
                    double max_pt2 = std::max(SubJet_pt[FatJet_subJetIdx1[fj]], SubJet_pt[FatJet_subJetIdx2[fj]]);
                    double min_pt2 = std::min(SubJet_pt[FatJet_subJetIdx1[fj]], SubJet_pt[FatJet_subJetIdx2[fj]]);
                    h_SJ_z1mtTT->Fill(max_pt2/min_pt2);
                    h_SJ_z2mtTT->Fill(min_pt2/(SubJet_pt[FatJet_subJetIdx1[fj]] + SubJet_pt[FatJet_subJetIdx2[fj]]));

                    double E1;
                    E1 = std::sqrt(SubJet_pt[FatJet_subJetIdx1[fj]]*SubJet_pt[FatJet_subJetIdx1[fj]] + SubJet_mass[FatJet_subJetIdx1[fj]]*SubJet_mass[FatJet_subJetIdx1[fj]]);
                    double E2;
                    E2 = std::sqrt(SubJet_pt[FatJet_subJetIdx2[fj]]*SubJet_pt[FatJet_subJetIdx2[fj]] + SubJet_mass[FatJet_subJetIdx2[fj]]*SubJet_mass[FatJet_subJetIdx2[fj]]);
                    h_SJ_pemtTT->Fill(std::abs(E1-E2)/FatJet_pt[fj]);
                }

                if(FatJet_pt[fj]>220){
                    double max_pt3 = std::max(SubJet_pt[FatJet_subJetIdx1[fj]], SubJet_pt[FatJet_subJetIdx2[fj]]);
                    double min_pt3 = std::min(SubJet_pt[FatJet_subJetIdx1[fj]], SubJet_pt[FatJet_subJetIdx2[fj]]);
                    h_SJ_z1pt1TT->Fill(max_pt3/min_pt3);
                    h_SJ_z2pt1TT->Fill(min_pt3/(SubJet_pt[FatJet_subJetIdx1[fj]] + SubJet_pt[FatJet_subJetIdx2[fj]]));
                
                    double E1, E2;
                    E1 = std::sqrt(SubJet_pt[FatJet_subJetIdx1[fj]]*SubJet_pt[FatJet_subJetIdx1[fj]] + SubJet_mass[FatJet_subJetIdx1[fj]]*SubJet_mass[FatJet_subJetIdx1[fj]]);
                    E2 = std::sqrt(SubJet_pt[FatJet_subJetIdx2[fj]]*SubJet_pt[FatJet_subJetIdx2[fj]] + SubJet_mass[FatJet_subJetIdx2[fj]]*SubJet_mass[FatJet_subJetIdx2[fj]]);
                    h_SJ_pept1TT->Fill(std::abs(E1-E2)/FatJet_pt[fj]);
                }
                if(FatJet_pt[fj]>250){
                    double max_pt4 = std::max(SubJet_pt[FatJet_subJetIdx1[fj]], SubJet_pt[FatJet_subJetIdx2[fj]]);
                    double min_pt4 = std::min(SubJet_pt[FatJet_subJetIdx1[fj]], SubJet_pt[FatJet_subJetIdx2[fj]]);
                    h_SJ_z1pt2TT->Fill(max_pt4/min_pt4);
                    h_SJ_z2pt2TT->Fill(min_pt4/(SubJet_pt[FatJet_subJetIdx1[fj]] + SubJet_pt[FatJet_subJetIdx2[fj]]));

                    double E1, E2;
                    E1 = std::sqrt(SubJet_pt[FatJet_subJetIdx1[fj]]*SubJet_pt[FatJet_subJetIdx1[fj]] + SubJet_mass[FatJet_subJetIdx1[fj]]*SubJet_mass[FatJet_subJetIdx1[fj]]);
                    E2 = std::sqrt(SubJet_pt[FatJet_subJetIdx2[fj]]*SubJet_pt[FatJet_subJetIdx2[fj]] + SubJet_mass[FatJet_subJetIdx2[fj]]*SubJet_mass[FatJet_subJetIdx2[fj]]);
                    h_SJ_pept2TT->Fill(std::abs(E1-E2)/FatJet_pt[fj]);

                }
                if(FatJet_pt[fj]>280){
                    double max_pt5 = std::max(SubJet_pt[FatJet_subJetIdx1[fj]], SubJet_pt[FatJet_subJetIdx2[fj]]);
                    double min_pt5 = std::min(SubJet_pt[FatJet_subJetIdx1[fj]], SubJet_pt[FatJet_subJetIdx2[fj]]);
                    h_SJ_z1pt3TT->Fill(max_pt5/min_pt5);
                    h_SJ_z2pt3TT->Fill(min_pt5/(SubJet_pt[FatJet_subJetIdx1[fj]] + SubJet_pt[FatJet_subJetIdx2[fj]]));

                    double E1, E2;
                    E1 = std::sqrt(SubJet_pt[FatJet_subJetIdx1[fj]]*SubJet_pt[FatJet_subJetIdx1[fj]] + SubJet_mass[FatJet_subJetIdx1[fj]]*SubJet_mass[FatJet_subJetIdx1[fj]]);
                    E2 = std::sqrt(SubJet_pt[FatJet_subJetIdx2[fj]]*SubJet_pt[FatJet_subJetIdx2[fj]] + SubJet_mass[FatJet_subJetIdx2[fj]]*SubJet_mass[FatJet_subJetIdx2[fj]]);
                    h_SJ_pept3TT->Fill(std::abs(E1-E2)/FatJet_pt[fj]);

                }
                if(FatJet_pt[fj]>310){
                    double max_pt6 = std::max(SubJet_pt[FatJet_subJetIdx1[fj]], SubJet_pt[FatJet_subJetIdx2[fj]]);
                    double min_pt6 = std::min(SubJet_pt[FatJet_subJetIdx1[fj]], SubJet_pt[FatJet_subJetIdx2[fj]]);
                    h_SJ_z1pt4TT->Fill(max_pt6/min_pt6);
                    h_SJ_z2pt4TT->Fill(min_pt6/(SubJet_pt[FatJet_subJetIdx1[fj]] + SubJet_pt[FatJet_subJetIdx2[fj]]));

                    double E1, E2;
                    E1 = std::sqrt(SubJet_pt[FatJet_subJetIdx1[fj]]*SubJet_pt[FatJet_subJetIdx1[fj]] + SubJet_mass[FatJet_subJetIdx1[fj]]*SubJet_mass[FatJet_subJetIdx1[fj]]);
                    E2 = std::sqrt(SubJet_pt[FatJet_subJetIdx2[fj]]*SubJet_pt[FatJet_subJetIdx2[fj]] + SubJet_mass[FatJet_subJetIdx2[fj]]*SubJet_mass[FatJet_subJetIdx2[fj]]);
                    h_SJ_pept4TT->Fill(std::abs(E1-E2)/FatJet_pt[fj]);
                }
                
            }  
            
            h_tau21_msdTT->Fill((FatJet_tau2[fj]/FatJet_tau1[fj]), FatJet_msoftdrop[fj]);
            h_tau32_msdTT->Fill((FatJet_tau3[fj]/FatJet_tau2[fj]), FatJet_msoftdrop[fj]);

            h_SJ_etadTT->Fill(SubJet_eta[FatJet_subJetIdx1[fj]] - SubJet_eta[FatJet_subJetIdx2[fj]]);
            h_SJ_phidTT->Fill(SubJet_phi[FatJet_subJetIdx1[fj]] - SubJet_phi[FatJet_subJetIdx2[fj]]);
            h_SJ_ptdTT->Fill(SubJet_pt[FatJet_subJetIdx1[fj]] - SubJet_pt[FatJet_subJetIdx2[fj]]);

        }
        for(Int_t genfj=0; genfj<nGenFatJet; genfj++){
            h_GFJ_massTT->Fill(GenFatJet_mass[genfj]);
            h_GFJ_ptTT->Fill(GenFatJet_pt[genfj]);
            h_GFJ_etaTT->Fill(GenFatJet_eta[genfj]);
            h_GFJ_phiTT->Fill(GenFatJet_phi[genfj]);
        }
    }



    

    gStyle->SetOptStat(0);

    /*ISTOGRAMMI PER pt*/
    TCanvas *c_pt_tot = new TCanvas("c_pt_tot", "Consituents p_{T} LL vs TT", 2000, 1100);
    c_pt_tot->Divide(2,1);

    TLegend *leg7 = new TLegend(0.5, 0.65, 0.65, 0.8);
    leg7->AddEntry(h_GJ_ptLL, "L", "l");
    leg7->AddEntry(h_GJ_ptTT, "T", "l");
    leg7->Draw();
//prima gli LL?
    c_pt_tot->cd(1);
    h_GFJ_ptLL->SetTitle("GenJetAK8 p_{t}");
    gPad->SetLeftMargin(0.15);
    h_GFJ_ptLL->Scale(1.0 / h_GFJ_ptLL->Integral());
    h_GFJ_ptLL->Draw("hist");
    h_GFJ_ptLL->SetFillColor(kRed-4);
    h_GFJ_ptTT->Scale(1.0 / h_GFJ_ptTT->Integral());
    h_GFJ_ptTT->SetMarkerColor(kBlack);
    h_GFJ_ptTT->SetMarkerSize(1.5);
    h_GFJ_ptTT->SetMarkerStyle(20);
    h_GFJ_ptTT->Draw("P same");


    c_pt_tot->cd(2);
    gPad->SetLeftMargin(0.15);
    h_FJ_ptLL->SetTitle("FatJet p_{t}");
    h_FJ_ptLL->Scale(1.0 / h_FJ_ptLL->Integral());
    h_FJ_ptLL->Draw("hist");
    h_FJ_ptLL->SetLineColor(kBlue+2);
    h_FJ_ptLL->SetFillColorAlpha(kBlue-7, 0.25);
    h_FJ_ptTT->Scale(1.0 / h_FJ_ptTT->Integral());
    h_FJ_ptTT->SetMarkerColor(kBlack);
    h_FJ_ptTT->SetMarkerSize(1.5);
    h_FJ_ptTT->SetMarkerStyle(20);
    h_FJ_ptTT->Draw("P same");

    c_pt_tot->SaveAs("plots_pt/c_pt_tot.pdf");
    


    /*ISTOGRAMMI PER ETA*/
    TCanvas *c_eta_tot = new TCanvas("c_eta_tot", "Consituents #eta: pseudorapirity LL vs TT", 2000, 1100);
    c_eta_tot->Divide(2,1);

    TLegend *leg6 = new TLegend(0.2, 0.65, 0.35, 0.8);
    leg6->AddEntry(h_GJ_etaTT, "TT", "l");
    leg6->AddEntry(h_GJ_etaLL, "LL", "l");
    leg6->Draw();

    c_eta_tot->cd(1);
    h_GFJ_etaTT->SetTitle("GenJetAK8 #eta");
    h_GFJ_etaTT->Scale(1.0 / h_GFJ_etaTT->Integral());
    h_GFJ_etaTT->Draw("hist");
    h_GFJ_etaTT->SetLineColor(kRed);
    h_GFJ_etaTT->SetFillColorAlpha(kMagenta-9, 0.25);
    h_GFJ_etaLL->Scale(1.0 / h_GFJ_etaLL->Integral());
    h_GFJ_etaLL->SetLineColor(kBlue+2);
    h_GFJ_etaLL->SetFillColorAlpha(kBlue-7, 0.25);
    h_GFJ_etaLL->Draw("hist same");


    c_eta_tot->cd(2);
    h_FJ_etaTT->SetTitle("FatJet #eta");
    h_FJ_etaTT->Draw("hist");
    h_FJ_etaTT->Scale(1.0 / h_FJ_etaTT->Integral());
    h_FJ_etaTT->SetLineColor(kRed);
    h_FJ_etaTT->SetFillColorAlpha(kMagenta-9, 0.25);
    h_FJ_etaLL->Scale(1.0 / h_FJ_etaLL->Integral());
    h_FJ_etaLL->SetLineColor(kBlue+2);
    h_FJ_etaLL->SetFillColorAlpha(kBlue-7, 0.25);
    h_FJ_etaLL->Draw("hist same");

    c_eta_tot->SaveAs("plots_eta/c_eta_tot.pdf");




    /*ISTOGRAMMI PER PHI*/
    TCanvas *c_phi_tot = new TCanvas("c_phi_tot", "Consituents #phi: LL vs TT", 2000, 1100);
    c_phi_tot->Divide(2,1);

    c_phi_tot->cd(1);
    h_GFJ_phiTT->SetTitle("GenJetAK8 #phi");
    h_GFJ_phiTT->Scale(1.0 / h_GFJ_phiTT->Integral());
    h_GFJ_phiTT->Draw("hist");
    h_GFJ_phiTT->SetLineColor(kRed);
    h_GFJ_phiTT->SetFillColorAlpha(kMagenta-9, 0.25);
    h_GFJ_phiLL->Scale(1.0 / h_GFJ_phiLL->Integral());
    h_GFJ_phiLL->SetLineColor(kBlue+2);
    h_GFJ_phiLL->SetFillColorAlpha(kBlue-7, 0.25);
    h_GFJ_phiLL->Draw("hist same");


    c_phi_tot->cd(2);
    h_FJ_phiTT->SetTitle("FatJet #phi");
    h_FJ_phiTT->Scale(1.0 / h_FJ_phiTT->Integral());
    h_FJ_phiTT->Draw("hist");
    h_FJ_phiTT->SetLineColor(kRed);
    h_FJ_phiTT->SetFillColorAlpha(kMagenta-9, 0.25);
    h_FJ_phiLL->Scale(1.0 / h_FJ_phiLL->Integral());
    h_FJ_phiLL->SetLineColor(kBlue+2);
    h_FJ_phiLL->SetFillColorAlpha(kBlue-7, 0.25);
    h_FJ_phiLL->Draw("hist same");

    c_phi_tot->SaveAs("plots_phi/c_phi_tot.pdf");
    


    
    /*ISTOGRAMMI PER LE MASSE*/
    TCanvas *c_mass_tot_pol = new TCanvas("c_mass_tot_pol", "Consituents Mass LL vs TT", 2000, 1100);
    c_mass_tot_pol->Divide(2,1);

    TLegend *leg11 = new TLegend(0.2, 0.65, 0.35, 0.8);
    leg11->AddEntry(h_GJ_massLL, "LL", "l");
    leg11->AddEntry(h_GJ_massTT, "TT", "l");
    leg11->Draw();

    c_mass_tot_pol->cd(1);
    h_GFJ_massLL->SetTitle("GenJetAK8 mass");
    h_GFJ_massLL->Scale(1.0 / h_GFJ_massLL->Integral());
    h_GFJ_massLL->SetLineColor(kBlue+2);
    h_GFJ_massLL->SetFillColorAlpha(kBlue-7, 0.25);
    h_GFJ_massLL->Draw("hist");
    h_GFJ_massTT->Scale(1.0 / h_GFJ_massTT->Integral());
    h_GFJ_massTT->Draw("hist same");
    h_GFJ_massTT->SetLineColor(kRed);
    h_GFJ_massTT->SetFillColorAlpha(kMagenta-9, 0.25);
    

    c_mass_tot_pol->cd(2);
    h_FJ_massLL->SetTitle("FatJet mass");
    h_FJ_massLL->Scale(1.0 / h_FJ_massLL->Integral());
    h_FJ_massLL->SetLineColor(kBlue+2);
    h_FJ_massLL->SetFillColorAlpha(kBlue-7, 0.25);
    h_FJ_massLL->Draw("hist");
    h_FJ_massTT->Scale(1.0 / h_FJ_massTT->Integral());
    h_FJ_massTT->Draw("hist same");
    h_FJ_massTT->SetLineColor(kRed);
    h_FJ_massTT->SetFillColorAlpha(kMagenta-9, 0.25);
    
    c_mass_tot_pol->SaveAs("plots_m/c_mass_tot_pol.pdf");



    /*Soft Drop*/
    TCanvas *c_FJ_msd_pol = new TCanvas("c_FJ_msd_pol", "FatJet SoftDrop mass LL vs TT", 2000, 1100);
    c_FJ_msd_pol->Divide(2,1);

    c_FJ_msd_pol->cd(1);
    h_FJ_msdLL->Scale(1.0 / h_FJ_msdLL->Integral());
    h_FJ_msdLL->Draw("hist");
    h_FJ_msdLL->SetLineColor(kRed);
    h_FJ_msdLL->SetFillColorAlpha(kMagenta-9, 0.25);
    h_FJ_massLL->Scale(1.0 / h_FJ_massLL->Integral());
    h_FJ_massLL->Draw("hist same");
    h_FJ_massLL->SetLineColor(kBlue+2);
    h_FJ_massLL->SetFillColorAlpha(kBlue-0, 0.25);

    TLegend *leg9 = new TLegend(0.5, 0.6, 0.9, 0.8);
    leg9->AddEntry(h_FJ_msdLL, "Mass SoftDrop", "l");
    leg9->AddEntry(h_FJ_massLL, "All masses", "l");
    leg9->SetTextSize(0.04);
    leg9->SetTextFont(42);
    leg9->Draw();
    
    c_FJ_msd_pol->cd(2);
    h_FJ_msdTT->Scale(1.0 / h_FJ_msdTT->Integral());
    h_FJ_msdTT->Draw("hist");
    h_FJ_msdTT->SetLineColor(kRed);
    h_FJ_msdTT->SetFillColorAlpha(kMagenta-9, 0.25);
    h_FJ_massTT->Scale(1.0 / h_FJ_massTT->Integral());
    h_FJ_massTT->Draw("hist same");
    h_FJ_massTT->SetLineColor(kBlue+2);
    h_FJ_massTT->SetFillColorAlpha(kBlue-0, 0.25);

    
    c_FJ_msd_pol->SaveAs("c_FJ_msd_pol.pdf");



    /*ISTOGRAMMI TAU 1,2, 3 */
    gStyle->SetOptStat(0);
    TCanvas *c_tau = new TCanvas("c_tau", "FatJet tau", 3000, 1200);
    c_tau->Divide(3,1);

    c_tau->cd(1);
    h_tau1TT->SetTitle("#tau_{1}");
    h_tau1TT->Scale(1.0 / h_tau1TT->Integral());
    h_tau1TT->SetLineColor(kRed);
    h_tau1TT->SetFillColorAlpha(kMagenta-9, 0.25);
    h_tau1TT->Draw("hist");
    h_tau1LL->Scale(1.0 / h_tau1LL->Integral());
    h_tau1LL->Draw("hist same");
    h_tau1LL->SetLineColor(kBlue+2);
    h_tau1LL->SetFillColorAlpha(kBlue-7, 0.25);

    TLegend *leg1 = new TLegend(0.6, 0.6, 0.8, 0.8);
    leg1->AddEntry(h_tau1TT, "TT", "l");
    leg1->AddEntry(h_tau1LL, "LL", "l");
    leg1->Draw();

    c_tau->cd(2);
    h_tau2TT->SetTitle("#tau_{2}");
    h_tau2TT->Scale(1.0 / h_tau2TT->Integral());
    h_tau2TT->SetLineColor(kRed);
    h_tau2TT->SetFillColorAlpha(kMagenta-9, 0.25);
    h_tau2TT->Draw("hist");
    h_tau2LL->Scale(1.0 / h_tau2LL->Integral());
    h_tau2LL->Draw("hist same");
    h_tau2LL->SetLineColor(kBlue+2);
    h_tau2LL->SetFillColorAlpha(kBlue-7, 0.25);

    c_tau->cd(3);
    h_tau3TT->SetTitle("#tau_{3}");
    h_tau3TT->Scale(1.0 / h_tau3TT->Integral());
    h_tau3TT->SetLineColor(kRed);
    h_tau3TT->SetFillColorAlpha(kMagenta-9, 0.25);
    h_tau3TT->Draw("hist");
    h_tau3LL->Scale(1.0 / h_tau3LL->Integral());
    h_tau3LL->Draw("hist same");
    h_tau3LL->SetLineColor(kBlue+2);
    h_tau3LL->SetFillColorAlpha(kBlue-7, 0.25);

    c_tau->SaveAs("plots_tau/c_tau.pdf");



    /*ISTOGRAMMI Tau21, Tau32 confronto con eventi "W" */
    gStyle->SetOptStat(0);
    TCanvas *c_tau_ratio = new TCanvas("c_tau_ratio", "FatJet tau ratios", 2000, 2000);
    c_tau_ratio->Divide(2,2);

    c_tau_ratio->cd(1);
    h_tau21TT->SetTitle("#tau_{21}");
    h_tau21TT->Scale(1.0 / h_tau21TT->Integral());
    h_tau21TT->SetLineColor(kRed);
    h_tau21TT->SetFillColorAlpha(kMagenta-9, 0.25);
    h_tau21TT->Draw("hist");
    h_tau21LL->Scale(1.0 / h_tau21LL->Integral());
    h_tau21LL->Draw("hist same");
    h_tau21LL->SetLineColor(kBlue+2);
    h_tau21LL->SetFillColorAlpha(kBlue-7, 0.25);

    TLegend *leg2 = new TLegend(0.2, 0.6, 0.4, 0.8);
    leg2->AddEntry(h_tau1TT, "TT", "l");
    leg2->AddEntry(h_tau1LL, "LL", "l");
    leg2->Draw();

    c_tau_ratio->cd(2);
    h_tau32TT->SetTitle("#tau_{32}");
    h_tau32TT->Scale(1.0 / h_tau32TT->Integral());
    h_tau32TT->SetLineColor(kRed);
    h_tau32TT->SetFillColorAlpha(kMagenta-9, 0.25);
    h_tau32TT->Draw("hist");
    h_tau32LL->Scale(1.0 / h_tau32LL->Integral());
    h_tau32LL->Draw("hist same");
    h_tau32LL->SetLineColor(kBlue+2);
    h_tau32LL->SetFillColorAlpha(kBlue-7, 0.25);

    c_tau_ratio->cd(3);
    h_tau21mTT->SetTitle("#tau_{21} 60 < m_{sd} < 100");
    h_tau21mTT->Scale(1.0 / h_tau21mTT->Integral());
    h_tau21mTT->SetLineColor(kRed);
    h_tau21mTT->SetFillColorAlpha(kMagenta-9, 0.25);
    h_tau21mTT->Draw("hist");
    h_tau21mLL->Scale(1.0 / h_tau21mLL->Integral());
    h_tau21mLL->Draw("hist same");
    h_tau21mLL->SetLineColor(kBlue+2);
    h_tau21mLL->SetFillColorAlpha(kBlue-7, 0.25);

    TLegend *leg3 = new TLegend(0.6, 0.6, 0.8, 0.8);
    leg3->AddEntry(h_tau1TT, "TT", "l");
    leg3->AddEntry(h_tau1LL, "LL", "l");
    leg3->Draw();

    c_tau_ratio->cd(4);
    h_tau32mTT->SetTitle("#tau_{32} 60 < m_{sd} < 100");
    h_tau32mTT->Scale(1.0 / h_tau32mTT->Integral());
    h_tau32mTT->SetLineColor(kRed);
    h_tau32mTT->SetFillColorAlpha(kMagenta-9, 0.25);
    h_tau32mTT->Draw("hist");
    h_tau32mLL->Scale(1.0 / h_tau32mLL->Integral());
    h_tau32mLL->Draw("hist same");
    h_tau32mLL->SetLineColor(kBlue+2);
    h_tau32mLL->SetFillColorAlpha(kBlue-7, 0.25);

    c_tau_ratio->SaveAs("plots_tau/c_tau_ratio.pdf");

    
    /*SCATTER PLOT TAU - MASSA SD*/
    gStyle->SetOptStat(0);
    TCanvas *c_tau_msd = new TCanvas("c_tau_msd","Tau vs msd", 2500,1600);
    c_tau_msd->Divide(2,2);

    c_tau_msd->cd(1);
    h_tau21_msdTT->Scale(1.0/ h_tau21_msdTT->Integral());
    h_tau21_msdTT->GetXaxis()->SetTitleSize(0.05);
    h_tau21_msdTT->GetYaxis()->SetTitleSize(0.05);
    h_tau21_msdTT->Draw("COLZ"); 


    c_tau_msd->cd(2);
    h_tau21_msdLL->Scale(1.0/ h_tau21_msdLL->Integral());
    h_tau21_msdLL->GetXaxis()->SetTitleSize(0.05);
    h_tau21_msdLL->GetYaxis()->SetTitleSize(0.05);
    h_tau21_msdLL->Draw("colz");

    c_tau_msd->cd(3);
    h_tau32_msdTT->Scale(1.0/ h_tau32_msdTT->Integral());
    h_tau32_msdTT->GetXaxis()->SetTitleSize(0.05);
    h_tau32_msdTT->GetYaxis()->SetTitleSize(0.05);
    h_tau32_msdTT->Draw("COLZ"); 

    c_tau_msd->cd(4);
    h_tau32_msdLL->Scale(1.0/ h_tau32_msdLL->Integral());
    h_tau21_msdLL->GetXaxis()->SetTitleSize(0.05);
    h_tau21_msdLL->GetYaxis()->SetTitleSize(0.05);
    h_tau32_msdLL->Draw("colz");

    c_tau_msd->SaveAs("plots_tau/c_tau_msd.pdf");
    


    /*ISTOGRAMMI per SubJet*/
//eta
    TCanvas *c_SJ_etad = new TCanvas("c_SJ_etad", "SubJet #eta difference: LL vs TT", 800, 600);
    c_SJ_etad->Divide(2,1);

    c_SJ_etad->cd(1);
    h_SJ_etad_sdTT->Scale(1.0 / h_SJ_etad_sdTT->Integral());
    h_SJ_etad_sdTT->Draw("hist");
    h_SJ_etad_sdTT->SetLineColor(kRed);
    h_SJ_etad_sdTT->SetTitle("SubJet #eta difference: 60 < m_sd < 100");
    h_SJ_etad_sdLL->Scale(1.0 / h_SJ_etad_sdLL->Integral());
    h_SJ_etad_sdLL->Draw("hist same");
    h_SJ_etad_sdLL->SetLineColor(kBlue);

    c_SJ_etad->cd(2);
    h_SJ_etadTT->Scale(1.0 / h_SJ_etadTT->Integral());
    h_SJ_etadTT->Draw("hist");
    h_SJ_etadTT->SetLineColor(kRed);
    h_SJ_etadTT->SetTitle("SubJet #eta difference: all masses");
    h_SJ_etadLL->Scale(1.0 / h_SJ_etadLL->Integral());
    h_SJ_etadLL->Draw("hist same");
    h_SJ_etadLL->SetLineColor(kBlue);

    TLegend *leg4 = new TLegend(0.2, 0.65, 0.35, 0.8);
    leg4->AddEntry(h_SJ_etadTT, "TT", "l");
    leg4->AddEntry(h_SJ_etadLL, "LL", "l");
    leg4->Draw();

    c_SJ_etad->SaveAs("plots_sj/c_SJ_etad.pdf");


//phi
    gStyle->SetOptStat(0);
    TCanvas *c_SJ_phid = new TCanvas("c_SJ_phid", "SubJet phi difference: LL vs TT", 800, 600);
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

    c_SJ_phid->SaveAs("plots_sj/c_SJ_phid.pdf");


//pt
    TCanvas *c_SJ_ptd = new TCanvas("c_SJ_ptd", "SubJet p_{t} difference: LL vs TT", 800, 600);
    c_SJ_ptd->Divide(2,1);

    c_SJ_ptd->cd(1);
    h_SJ_ptd_sdTT->Draw();
    h_SJ_ptd_sdTT->SetTitle("SubJet p_{t} difference: 60 < m_sd < 100");
    h_SJ_ptd_sdLL->Draw("same");
    h_SJ_ptd_sdLL->SetLineColor(kRed);


    c_SJ_ptd->cd(2);
    h_SJ_ptdTT->Draw();
    h_SJ_ptdTT->SetTitle("SubJet p_{t} difference: all masses");
    h_SJ_ptdLL->Draw("same");
    h_SJ_ptdLL->SetLineColor(kRed);

    TLegend *leg8 = new TLegend(0.2, 0.65, 0.35, 0.8);
    leg8->AddEntry(h_SJ_ptd_sdTT, "TT", "l");
    leg8->AddEntry(h_SJ_ptd_sdLL, "LL", "l");
    leg8->Draw();

    c_SJ_ptd->SaveAs("plots_sj/c_SJ_ptd.pdf");




    TCanvas *c_SJ_pts = new TCanvas("c_SJ_pts", "Subjets p_{t} 1 vs 2", 1200, 600);
    c_SJ_pts->Divide(2,1);

    c_SJ_pts->cd(1);
    h_SJ_pts2LL->SetTitle("SubjetLL #1 & #2 p_{t} (60 < m_{sd} < 100)");
    h_SJ_pts2LL->Scale(1.0 / h_SJ_pts2LL->Integral());
    h_SJ_pts2LL->Draw("hist");
    h_SJ_pts2LL->SetLineColor(kMagenta+2);
    h_SJ_pts2LL->SetFillColorAlpha(kMagenta-9, 0.25);

    h_SJ_pts1LL->Scale(1.0 / h_SJ_pts1LL->Integral());
    h_SJ_pts1LL->Draw("hist same");
    h_SJ_pts1LL->SetLineColor(kOrange+9);
    h_SJ_pts1LL->SetFillColorAlpha(kOrange+6, 0.25);

    
    TLegend *leg13 = new TLegend(0.5, 0.65, 0.65, 0.8);
    leg13->AddEntry(h_SJ_pts1LL, "#1", "l");
    leg13->AddEntry(h_SJ_pts2LL, "#2", "l");
    leg13->Draw();


    c_SJ_pts->cd(2);
    h_SJ_pts2TT->SetTitle("SubjetTT #1 & #2 p_{t} (60 < m_{sd} < 100)");
    h_SJ_pts2TT->Scale(1.0 / h_SJ_pts2TT->Integral());
    h_SJ_pts2TT->Draw("hist");
    h_SJ_pts2TT->SetLineColor(kMagenta+2);
    h_SJ_pts2TT->SetFillColorAlpha(kMagenta-9, 0.25);

    h_SJ_pts1TT->Scale(1.0 / h_SJ_pts1TT->Integral());
    h_SJ_pts1TT->Draw("hist same");
    h_SJ_pts1TT->SetLineColor(kOrange+9);
    h_SJ_pts1TT->SetFillColorAlpha(kOrange+6, 0.25);

    

    c_SJ_pts->SaveAs("plots_pt/c_SJ_pts.pdf");
    

//zg1
    TCanvas *c_SJ_z1m = new TCanvas("c_SJ_z1m", "SubJet z_{g}: LL vs TT", 1200, 600);
    c_SJ_z1m->Divide(2,1);

    c_SJ_z1m->cd(1);
    h_SJ_z1mTT->Scale(1.0/h_SJ_z1mTT->Integral());
    h_SJ_z1mTT->SetLineColor(kRed);
    h_SJ_z1mTT->SetFillColorAlpha(kMagenta-9, 0.25);
    h_SJ_z1mTT->Draw("hist");
    h_SJ_z1mTT->SetTitle("SubJet z_{g}: 60 < m_{sd} < 100");
    h_SJ_z1mLL->Scale(1.0/h_SJ_z1mLL->Integral());
    h_SJ_z1mLL->SetLineColor(kBlue+2);
    h_SJ_z1mLL->SetFillColorAlpha(kBlue-7, 0.25);
    h_SJ_z1mLL->Draw("hist same");


    c_SJ_z1m->cd(2);
    h_SJ_z1mtTT->Scale(1.0/h_SJ_z1mtTT->Integral());
    h_SJ_z1mtTT->SetLineColor(kRed);
    h_SJ_z1mtTT->SetFillColorAlpha(kMagenta-9, 0.25);
    h_SJ_z1mtTT->Draw("hist");
    h_SJ_z1mtTT->SetTitle("SubJet z_{g}: 60 < m_{sd} < 100  &  #tau_{21} < 0.45");
    h_SJ_z1mtLL->Scale(1.0/h_SJ_z1mtLL->Integral());
    h_SJ_z1mtLL->SetLineColor(kBlue+2);
    h_SJ_z1mtLL->SetFillColorAlpha(kBlue-7, 0.25);
    h_SJ_z1mtLL->Draw("hist same");

    TLegend *leg12 = new TLegend(0.4, 0.65, 0.55, 0.8);
    leg12->AddEntry(h_SJ_z1mtTT, "TT", "l");
    leg12->AddEntry(h_SJ_z1mtLL, "LL", "l");
    leg12->Draw();

    c_SJ_z1m->SaveAs("plots_z1/c_SJ_z1m.pdf");



    TCanvas *c_SJ_z1pt = new TCanvas("c_SJ_z1pt", "SubJet z_{g}: LL vs TT", 1600, 1200);
    c_SJ_z1pt->Divide(2,2);

    c_SJ_z1pt->cd(1);
    h_SJ_z1pt1TT->Scale(1.0/h_SJ_z1pt1TT->Integral());
    h_SJ_z1pt1TT->SetLineColor(kRed);
    h_SJ_z1pt1TT->SetFillColorAlpha(kMagenta-9, 0.25);
    h_SJ_z1pt1TT->GetXaxis()->SetTitleSize(0.05);
    h_SJ_z1pt1TT->GetYaxis()->SetTitleSize(0.05);
    h_SJ_z1pt1TT->Draw("hist");
    h_SJ_z1pt1TT->SetTitle("SubJet z_{g}: p_{t}>220");
    h_SJ_z1pt1LL->Scale(1.0/h_SJ_z1pt1LL->Integral());
    h_SJ_z1pt1LL->SetLineColor(kBlue+2);
    h_SJ_z1pt1LL->SetFillColorAlpha(kBlue-7, 0.25);
    h_SJ_z1pt1LL->GetXaxis()->SetTitleSize(0.05);
    h_SJ_z1pt1LL->GetYaxis()->SetTitleSize(0.05);
    h_SJ_z1pt1LL->Draw("hist same");

    TLegend *leg14 = new TLegend(0.4, 0.65, 0.55, 0.8);
    leg14->AddEntry(h_SJ_z1pt1TT, "TT", "l");
    leg14->AddEntry(h_SJ_z1pt1LL, "LL", "l");
    leg14->Draw();


    c_SJ_z1pt->cd(2);
    h_SJ_z1pt2TT->Scale(1.0/h_SJ_z1pt2TT->Integral());
    h_SJ_z1pt2TT->SetLineColor(kRed);
    h_SJ_z1pt2TT->SetFillColorAlpha(kMagenta-9, 0.25);
    h_SJ_z1pt2TT->GetXaxis()->SetTitleSize(0.05);
    h_SJ_z1pt2TT->GetYaxis()->SetTitleSize(0.05);
    h_SJ_z1pt2TT->Draw("hist");
    h_SJ_z1pt2TT->SetTitle("SubJet z_{g}: p_{t}>250");
    h_SJ_z1pt2LL->Scale(1.0/h_SJ_z1pt2LL->Integral());
    h_SJ_z1pt2LL->SetLineColor(kBlue+2);
    h_SJ_z1pt2LL->SetFillColorAlpha(kBlue-7, 0.25);
    h_SJ_z1pt2LL->GetXaxis()->SetTitleSize(0.05);
    h_SJ_z1pt2LL->GetYaxis()->SetTitleSize(0.05);
    h_SJ_z1pt2LL->Draw("hist same");


    c_SJ_z1pt->cd(3);
    h_SJ_z1pt3TT->Scale(1.0/h_SJ_z1pt3TT->Integral());
    h_SJ_z1pt3TT->SetLineColor(kRed);
    h_SJ_z1pt3TT->SetFillColorAlpha(kMagenta-9, 0.25);
    h_SJ_z1pt3TT->GetXaxis()->SetTitleSize(0.05);
    h_SJ_z1pt3TT->GetYaxis()->SetTitleSize(0.05);
    h_SJ_z1pt3TT->Draw("hist");
    h_SJ_z1pt3TT->SetTitle("SubJet z_{g}: p_{t}>280");
    h_SJ_z1pt3LL->Scale(1.0/h_SJ_z1pt3LL->Integral());
    h_SJ_z1pt3LL->SetLineColor(kBlue+2);
    h_SJ_z1pt3LL->SetFillColorAlpha(kBlue-7, 0.25);
    h_SJ_z1pt3LL->GetXaxis()->SetTitleSize(0.05);
    h_SJ_z1pt3LL->GetYaxis()->SetTitleSize(0.05);
    h_SJ_z1pt3LL->Draw("hist same");


    c_SJ_z1pt->cd(4);
    h_SJ_z1pt4TT->Scale(1.0/h_SJ_z1pt4TT->Integral());
    h_SJ_z1pt4TT->SetLineColor(kRed);
    h_SJ_z1pt4TT->SetFillColorAlpha(kMagenta-9, 0.25);
    h_SJ_z1pt4TT->GetXaxis()->SetTitleSize(0.05);
    h_SJ_z1pt4TT->GetYaxis()->SetTitleSize(0.05);
    h_SJ_z1pt4TT->Draw("hist");
    h_SJ_z1pt4TT->SetTitle("SubJet z_{g}: p_{t}>310");
    h_SJ_z1pt4LL->Scale(1.0/h_SJ_z1pt4LL->Integral());
    h_SJ_z1pt4LL->SetLineColor(kBlue+2);
    h_SJ_z1pt4LL->SetFillColorAlpha(kBlue-7, 0.25);
    h_SJ_z1pt4LL->GetXaxis()->SetTitleSize(0.05);
    h_SJ_z1pt4LL->GetYaxis()->SetTitleSize(0.05);
    h_SJ_z1pt4LL->Draw("hist same");

    c_SJ_z1pt->SaveAs("plots_z1/c_SJ_z1pt.pdf");

//zg2
    TCanvas *c_SJ_z2m = new TCanvas("c_SJ_z2m", "SubJet z_{g}: LL vs TT", 1200, 600);
    c_SJ_z2m->Divide(2,1);

    c_SJ_z2m->cd(1);
    h_SJ_z2mTT->Scale(1.0/h_SJ_z2mTT->Integral());
    h_SJ_z2mTT->SetLineColor(kRed);
    h_SJ_z2mTT->SetFillColorAlpha(kMagenta-9, 0.25);
    h_SJ_z2mTT->Draw("hist");
    h_SJ_z2mTT->SetTitle("SubJet z_{g}': 60 < m_{sd} < 100");
    h_SJ_z2mLL->Scale(1.0/h_SJ_z2mLL->Integral());
    h_SJ_z2mLL->SetLineColor(kBlue+2);
    h_SJ_z2mLL->SetFillColorAlpha(kBlue-7, 0.25);
    h_SJ_z2mLL->Draw("hist same");

    TLegend *leg15 = new TLegend(0.6, 0.65, 0.75, 0.8);
    leg15->AddEntry(h_SJ_z2mTT, "TT", "l");
    leg15->AddEntry(h_SJ_z2mLL, "LL", "l");
    leg15->Draw();

    c_SJ_z2m->cd(2);
    h_SJ_z2mtTT->Scale(1.0/h_SJ_z2mtTT->Integral());
    h_SJ_z2mtTT->SetLineColor(kRed);
    h_SJ_z2mtTT->SetFillColorAlpha(kMagenta-9, 0.25);
    h_SJ_z2mtTT->Draw("hist");
    h_SJ_z2mtTT->SetTitle("SubJet z_{g}': 60 < m_{sd} < 100  &  #tau_{21} < 0.45");
    h_SJ_z2mtLL->Scale(1.0/h_SJ_z2mtLL->Integral());
    h_SJ_z2mtLL->SetLineColor(kBlue+2);
    h_SJ_z2mtLL->SetFillColorAlpha(kBlue-7, 0.25);
    h_SJ_z2mtLL->Draw("hist same");

    

    c_SJ_z2m->SaveAs("plots_z1/c_SJ_z2m.pdf");




    TCanvas *c_SJ_z2pt = new TCanvas("c_SJ_z2pt", "SubJet z_{g}: LL vs TT", 1600, 1200);
    c_SJ_z2pt->Divide(2,2);

    c_SJ_z2pt->cd(1);
    h_SJ_z2pt1LL->SetTitle("SubJet z_{g}': p_{t}>220");
    h_SJ_z2pt1LL->Scale(1.0 / h_SJ_z2pt1LL->Integral());
    h_SJ_z2pt1LL->SetLineColor(kBlue+2);
    h_SJ_z2pt1LL->SetFillColorAlpha(kBlue-7, 0.25);
    h_SJ_z2pt1LL->GetXaxis()->SetTitleSize(0.05);
    h_SJ_z2pt1LL->GetYaxis()->SetTitleSize(0.05);
    h_SJ_z2pt1LL->Draw("hist");
    h_SJ_z2pt1TT->Scale(1.0 / h_SJ_z2pt1TT->Integral());
    h_SJ_z2pt1TT->SetLineColor(kRed);
    h_SJ_z2pt1TT->SetFillColorAlpha(kMagenta-9, 0.25);
    h_SJ_z2pt1TT->GetXaxis()->SetTitleSize(0.05);
    h_SJ_z2pt1TT->GetYaxis()->SetTitleSize(0.05);
    h_SJ_z2pt1TT->Draw("hist same");
    

    TLegend *leg16 = new TLegend(0.6, 0.65, 0.75, 0.8);
    leg16->AddEntry(h_SJ_z2pt1LL, "LL", "l");
    leg16->AddEntry(h_SJ_z2pt1TT, "TT", "l");
    leg16->Draw();

    c_SJ_z2pt->cd(2);
    h_SJ_z2pt2LL->SetTitle("SubJet z_{g}': p_{t}>250");
    h_SJ_z2pt2LL->Scale(1.0 / h_SJ_z2pt2LL->Integral());
    h_SJ_z2pt2LL->SetLineColor(kBlue+2);
    h_SJ_z2pt2LL->SetFillColorAlpha(kBlue-7, 0.25);
    h_SJ_z2pt2LL->GetXaxis()->SetTitleSize(0.05);
    h_SJ_z2pt2LL->GetYaxis()->SetTitleSize(0.05);
    h_SJ_z2pt2LL->Draw("hist");
    h_SJ_z2pt2TT->Scale(1.0 /h_SJ_z2pt2TT->Integral());
    h_SJ_z2pt2TT->SetLineColor(kRed);
    h_SJ_z2pt2TT->SetFillColorAlpha(kMagenta-9, 0.25);
    h_SJ_z2pt2TT->GetXaxis()->SetTitleSize(0.05);
    h_SJ_z2pt2TT->GetYaxis()->SetTitleSize(0.05);
    h_SJ_z2pt2TT->Draw("hist same");
    

    c_SJ_z2pt->cd(3);
    h_SJ_z2pt3LL->SetTitle("SubJet z_{g}': p_{t}>280");
    h_SJ_z2pt3LL->Scale(1.0 / h_SJ_z2pt3LL->Integral());
    h_SJ_z2pt3LL->SetLineColor(kBlue+2);
    h_SJ_z2pt3LL->SetFillColorAlpha(kBlue-7, 0.25);
    h_SJ_z2pt3LL->GetXaxis()->SetTitleSize(0.05);
    h_SJ_z2pt3LL->GetYaxis()->SetTitleSize(0.05);
    h_SJ_z2pt3LL->Draw("hist");
    h_SJ_z2pt3TT->Scale(1.0 /h_SJ_z2pt3TT->Integral());
    h_SJ_z2pt3TT->SetLineColor(kRed);
    h_SJ_z2pt3TT->SetFillColorAlpha(kMagenta-9, 0.25);
    h_SJ_z2pt3TT->GetXaxis()->SetTitleSize(0.05);
    h_SJ_z2pt3TT->GetYaxis()->SetTitleSize(0.05);
    h_SJ_z2pt3TT->Draw("hist same");


    c_SJ_z2pt->cd(4);
    h_SJ_z2pt4LL->SetTitle("SubJet z_{g}': p_{t}>310");
    h_SJ_z2pt4LL->Scale(1.0 / h_SJ_z2pt4LL->Integral());
    h_SJ_z2pt4LL->SetLineColor(kBlue+2);
    h_SJ_z2pt4LL->SetFillColorAlpha(kBlue-7, 0.25);
    h_SJ_z2pt4LL->GetXaxis()->SetTitleSize(0.05);
    h_SJ_z2pt4LL->GetYaxis()->SetTitleSize(0.05);
    h_SJ_z2pt4LL->Draw("hist");
    h_SJ_z2pt4TT->Scale(1.0 /h_SJ_z2pt4TT->Integral());
    h_SJ_z2pt4TT->SetLineColor(kRed);
    h_SJ_z2pt4TT->SetFillColorAlpha(kMagenta-9, 0.25);
    h_SJ_z2pt4TT->GetXaxis()->SetTitleSize(0.05);
    h_SJ_z2pt4TT->GetYaxis()->SetTitleSize(0.05);
    h_SJ_z2pt4TT->Draw("hist same");
    

    c_SJ_z2pt->SaveAs("plots_z1/c_SJ_z2pt.pdf");

//ptheta
    TCanvas *c_SJ_pem = new TCanvas("c_SJ_pem", "SubJet z_{g}: LL vs TT", 1200, 600);
    c_SJ_pem->Divide(2,1);

    c_SJ_pem->cd(1);
    h_SJ_pemTT->Scale(1.0/h_SJ_pemTT->Integral());
    h_SJ_pemTT->SetLineColor(kRed);
    h_SJ_pemTT->SetFillColorAlpha(kMagenta-9, 0.25);
    h_SJ_pemTT->Draw("hist");
    h_SJ_pemTT->SetTitle("SubJet p_{#theta}: 60 < m_{sd} < 100");
    h_SJ_pemLL->Scale(1.0/h_SJ_pemLL->Integral());
    h_SJ_pemLL->SetLineColor(kBlue+2);
    h_SJ_pemLL->SetFillColorAlpha(kBlue-7, 0.25);
    h_SJ_pemLL->Draw("hist same");

    TLegend *leg17 = new TLegend(0.2, 0.35, 0.4, 0.5);
    leg17->AddEntry(h_SJ_pemTT, "TT", "l");
    leg17->AddEntry(h_SJ_pemLL, "LL", "l");
    leg17->Draw();

    c_SJ_pem->cd(2);
    h_SJ_pemtTT->Scale(1.0/h_SJ_pemtTT->Integral());
    h_SJ_pemtTT->SetLineColor(kRed);
    h_SJ_pemtTT->SetFillColorAlpha(kMagenta-9, 0.25);
    h_SJ_pemtTT->Draw("hist");
    h_SJ_pemtTT->SetTitle("SubJet p_{#theta}: 60 < m_{sd} < 100  &  #tau_{21} < 0.45");
    h_SJ_pemtLL->Scale(1.0/h_SJ_pemtLL->Integral());
    h_SJ_pemtLL->SetLineColor(kBlue+2);
    h_SJ_pemtLL->SetFillColorAlpha(kBlue-7, 0.25);
    h_SJ_pemtLL->Draw("hist same");

    

    c_SJ_pem->SaveAs("plots_ptheta/c_SJ_pem.pdf");    




    TCanvas *c_SJ_pept = new TCanvas("c_SJ_pept", "SubJet z_{g}: LL vs TT", 1600, 1200);
    c_SJ_pept->Divide(2,2);

    c_SJ_pept->cd(1);
    h_SJ_pept1LL->SetTitle("SubJet p_{#theta}: p_{t}>220");
    h_SJ_pept1LL->Scale(1.0 / h_SJ_pept1LL->Integral());
    h_SJ_pept1LL->SetLineColor(kBlue+2);
    h_SJ_pept1LL->SetFillColorAlpha(kBlue-7, 0.25);
    h_SJ_pept1LL->GetXaxis()->SetTitleSize(0.05);
    h_SJ_pept1LL->GetYaxis()->SetTitleSize(0.05);
    h_SJ_pept1LL->Draw("hist");
    h_SJ_pept1TT->Scale(1.0 / h_SJ_pept1TT->Integral());
    h_SJ_pept1TT->SetLineColor(kRed);
    h_SJ_pept1TT->SetFillColorAlpha(kMagenta-9, 0.25);
    h_SJ_pept1TT->GetXaxis()->SetTitleSize(0.05);
    h_SJ_pept1TT->GetYaxis()->SetTitleSize(0.05);
    h_SJ_pept1TT->Draw("hist same");
    

    TLegend *leg18 = new TLegend(0.2, 0.25, 0.3, 0.4);
    leg18->AddEntry(h_SJ_pept1LL, "LL", "l");
    leg18->AddEntry(h_SJ_pept1TT, "TT", "l");
    leg18->Draw();

    c_SJ_pept->cd(2);
    h_SJ_pept2LL->SetTitle("SubJet p_{#theta}: p_{t}>250");
    h_SJ_pept2LL->Scale(1.0 / h_SJ_pept2LL->Integral());
    h_SJ_pept2LL->SetLineColor(kBlue+2);
    h_SJ_pept2LL->SetFillColorAlpha(kBlue-7, 0.25);
    h_SJ_pept2LL->GetXaxis()->SetTitleSize(0.05);
    h_SJ_pept2LL->GetYaxis()->SetTitleSize(0.05);
    h_SJ_pept2LL->Draw("hist");
    h_SJ_pept2TT->Scale(1.0 /h_SJ_pept2TT->Integral());
    h_SJ_pept2TT->SetLineColor(kRed);
    h_SJ_pept2TT->SetFillColorAlpha(kMagenta-9, 0.25);
    h_SJ_pept2TT->GetXaxis()->SetTitleSize(0.05);
    h_SJ_pept2TT->GetYaxis()->SetTitleSize(0.05);
    h_SJ_pept2TT->Draw("hist same");
    

    c_SJ_pept->cd(3);
    h_SJ_pept3LL->SetTitle("SubJet p_{#theta}: p_{t}>280");
    h_SJ_pept3LL->Scale(1.0 / h_SJ_pept3LL->Integral());
    h_SJ_pept3LL->SetLineColor(kBlue+2);
    h_SJ_pept3LL->SetFillColorAlpha(kBlue-7, 0.25);
    h_SJ_pept3LL->GetXaxis()->SetTitleSize(0.05);
    h_SJ_pept3LL->GetYaxis()->SetTitleSize(0.05);
    h_SJ_pept3LL->Draw("hist");
    h_SJ_pept3TT->Scale(1.0 /h_SJ_pept3TT->Integral());
    h_SJ_pept3TT->SetLineColor(kRed);
    h_SJ_pept3TT->SetFillColorAlpha(kMagenta-9, 0.25);
    h_SJ_pept3TT->GetXaxis()->SetTitleSize(0.05);
    h_SJ_pept3TT->GetYaxis()->SetTitleSize(0.05);
    h_SJ_pept3TT->Draw("hist same");


    c_SJ_pept->cd(4);
    h_SJ_pept4LL->SetTitle("SubJet p_{#theta}: p_{t}>310");
    h_SJ_pept4LL->Scale(1.0 / h_SJ_pept4LL->Integral());
    h_SJ_pept4LL->SetLineColor(kBlue+2);
    h_SJ_pept4LL->SetFillColorAlpha(kBlue-7, 0.25);
    h_SJ_pept4LL->GetXaxis()->SetTitleSize(0.05);
    h_SJ_pept4LL->GetYaxis()->SetTitleSize(0.05);
    h_SJ_pept4LL->Draw("hist");
    h_SJ_pept4TT->Scale(1.0 /h_SJ_pept4TT->Integral());
    h_SJ_pept4TT->SetLineColor(kRed);
    h_SJ_pept4TT->SetFillColorAlpha(kMagenta-9, 0.25);
    h_SJ_pept4TT->GetXaxis()->SetTitleSize(0.05);
    h_SJ_pept4TT->GetYaxis()->SetTitleSize(0.05);
    h_SJ_pept4TT->Draw("hist same");
    

    c_SJ_pept->SaveAs("plots_ptheta/c_SJ_pept.pdf");



    f->Write();
    f->Close();


    gROOT->ProcessLine(".q");



}
