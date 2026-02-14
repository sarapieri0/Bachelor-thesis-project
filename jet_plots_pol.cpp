/*Codice per trattare contemporaneamente entrambe le polarizzazioni LL e TT*/


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
#include "canva.cpp"


// lista dei branch nel tree : https://cms-xpog.docs.cern.ch/autoDoc/NanoAODv12/2022/2023/doc_DYJetsToLL_M-50_TuneCP5_13p6TeV-madgraphMLM-pythia8_Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2.html

void jet_plots2(
    TString output="jets2.root"
){
    TFile *f=new TFile(output, "recreate");

    Int_t maxEvents = 2000;
    std::string tagLL = "sWL";
    std::string tagTT = "sWT";
    

//pt    
    TH1F *h_FJ_ptLL = new TH1F(("h_FatJet_pt_" + tagLL).c_str(), "FatJet p_{T}   (CMS private work);p_{T} [GeV];Events", 50, 170, 1000);
    TH1F *h_FJ_ptTT = new TH1F(("h_FatJet_pt_" + tagTT).c_str(), "FatJet p_{T}   (CMS private work);p_{T} [GeV];Events", 50, 170, 1000);
    
    TH1F *h_GFJ_ptLL = new TH1F(("h_GenJetAK8_pt_" + tagLL).c_str(), "GenJetAK8 p_{T}   (CMS private work);p_{T} [GeV];Events", 50, 150, 1000);
    TH1F *h_GFJ_ptTT = new TH1F(("h_GenJetAK8_pt_" + tagTT).c_str(), "GenJetAK8 p_{T}   (CMS private work);p_{T} [GeV];Events", 50, 150, 1000);

    
    
//eta    
    TH1F *h_FJ_etaLL = new TH1F(("h_FatJet_eta_" + tagLL).c_str(), "FatJet #eta   (CMS private work); #eta; Events", 50,-8, 8);
    TH1F *h_FJ_etaTT = new TH1F(("h_FatJet_eta_" + tagTT).c_str(), "FatJet #eta   (CMS private work); #eta; Events", 50, -8, 8);
    
    TH1F *h_GFJ_etaLL = new TH1F(("h_GenJetAK8_eta_" + tagLL).c_str(), "GenJetAK8 #eta   (CMS private work);#eta; Events", 50, -8, 8);
    TH1F *h_GFJ_etaTT = new TH1F(("h_GenJetAK8_eta_" + tagTT).c_str(), "GenJetAK8 #eta   (CMS private work);#eta; Events", 50, -8, 8);

//phi
    TH1F *h_FJ_phiLL = new TH1F(("h_FatJet_phi_" + tagLL).c_str(), "FatJet #phi   (CMS private work); #phi; Events", 40, -4,4);
    TH1F *h_FJ_phiTT = new TH1F(("h_FatJet_phi_" + tagTT).c_str(), "FatJet #phi   (CMS private work); #phi; Events", 40, -4, 4);

    TH1F *h_GFJ_phiLL = new TH1F(("h_GenJetAK8_phi_" + tagLL).c_str(), "GenJetAK8 #phi   (CMS private work); #phi; Events", 40, -4, 4);
    TH1F *h_GFJ_phiTT = new TH1F(("h_GenJetAK8_phi_" + tagTT).c_str(), "GenJetAK8 #phi   (CMS private work); #phi; Events", 40, -4, 4);


//mass
    TH1F *h_FJ_massLL = new TH1F(("h_FJ_mass_" + tagLL).c_str(), "FatJet mass (L)   (CMS private work); Mass [GeV]; Events", 40,10,200);
    TH1F *h_FJ_massTT = new TH1F(("h_FJ_mass_" + tagTT).c_str(), "FatJet mass (T)   (CMS private work); Mass [GeV]; Events", 40,10,200);

    TH1F *h_GFJ_massLL = new TH1F(("h_GFJ_mass_" + tagLL).c_str(), "GenAK8 FatJet mass   (CMS private work);mass [GeV]; Events;", 40, 0, 300);
    TH1F *h_GFJ_massTT = new TH1F(("h_GFJ_mass_" + tagTT).c_str(), "GenAK8 FatJet mass   (CMS private work);mass [GeV]; Events;", 40, 0, 300);

//mass softdrop
    TH1F *h_FJ_msdLL = new TH1F(("h_FJ_msd_" + tagLL).c_str(), "FatJet SoftDrop mass (L)   (CMS private work); mass [GeV]; Events;", 40, 10, 200);
    TH1F *h_FJ_msdTT = new TH1F(("h_FJ_msd_" + tagTT).c_str(), "FatJet SoftDrop mass (T)   (CMS private work); mass [GeV]; Events;", 40, 10, 200);


//tau
    TH1F *h_tau1LL =new TH1F(("h_tau1_" + tagLL).c_str(), "FatJet #tau_{1}   (CMS private work); #tau_{1}; Events", 30, 0, 0.6);
    TH1F *h_tau2LL =new TH1F(("h_tau2_" + tagLL).c_str(), "FatJet #tau_{2}   (CMS private work); #tau_{2}; Events", 30, 0, 0.6);
    TH1F *h_tau3LL= new TH1F(("h_tau3_" + tagLL).c_str(), "FatJet #tau_{3}   (CMS private work); #tau_{3}; Events", 30, 0, 0.6);
    TH1F *h_tau1TT =new TH1F(("h_tau1_" + tagTT).c_str(), "FatJet #tau_{1}   (CMS private work); #tau_{1}; Events", 30, 0, 0.6);
    TH1F *h_tau2TT =new TH1F(("h_tau2_" + tagTT).c_str(), "FatJet #tau_{2}   (CMS private work); #tau_{2}; Events", 30, 0, 0.6);
    TH1F *h_tau3TT= new TH1F(("h_tau3_" + tagTT).c_str(), "FatJet #tau_{3}   (CMS private work); #tau_{3}; Events", 30, 0, 0.6);

    TH1F *h_tau21mLL= new TH1F(("h_tau21m_" + tagLL).c_str(), "#tau_{21} 60< m_sd <100   (CMS private work); #tau_{21}; Events", 30, 0, 1);
    TH1F *h_tau21LL= new TH1F(("h_tau21_" + tagLL).c_str(), "#tau_{21}   (CMS private work); #tau_{21}; Events", 30, 0, 1);
    TH1F *h_tau32mLL= new TH1F(("h_tau32m_" + tagLL).c_str(), "#tau_{32} 60< m_sd <100   (CMS private work); #tau_{32}; Events", 30, 0, 1);
    TH1F *h_tau32LL= new TH1F(("h_tau32_" + tagLL).c_str(), "#tau_{32}   (CMS private work); #tau_{32}; Events", 30, 0, 1);
    TH1F *h_tau21mTT= new TH1F(("h_tau21m_" + tagTT).c_str(), "#tau_{21} 60< m_sd <100   (CMS private work); #tau_{21}; Events", 30, 0, 1);
    TH1F *h_tau21TT= new TH1F(("h_tau21_" + tagTT).c_str(), "#tau_{21}   (CMS private work); #tau_{21}; Events", 30, 0, 1);
    TH1F *h_tau32mTT= new TH1F(("h_tau32m_" + tagTT).c_str(), "#tau_{32} 60< m_sd <100   (CMS private work); #tau_{32}; Events", 30, 0, 1);
    TH1F *h_tau32TT= new TH1F(("h_tau32_" + tagTT).c_str(), "#tau_{32}   (CMS private work); #tau_{32}; Events", 30, 0, 1);

    TH2F *h_tau21_msdLL = new TH2F(("h_msd_tau21_" + tagLL).c_str(), "SD mass vs  #tau_{21} (L)   (CMS private work); #tau_{21}; SoftDrop mass", 50, 0, 1,    30, 50, 300);
    TH2F *h_tau21_msdTT = new TH2F(("h_msd_tau21_" + tagTT).c_str(), "SD mass vs  #tau_{21} (T)   (CMS private work); #tau_{21}; SoftDrop mass", 50, 0, 1,    30, 50, 300);

    TH2F *h_tau32_msdLL = new TH2F(("h_msd_tau32_" + tagLL).c_str(), "SD mass vs #tau_{32} (L)   (CMS private work); #tau_{32}; SoftDrop mass", 50, 0, 1,    30, 50, 300);
    TH2F *h_tau32_msdTT = new TH2F(("h_msd_tau32_" + tagTT).c_str(), "SD mass vs #tau_{32} (T)   (CMS private work); #tau_{32}; SoftDrop mass", 50, 0, 1,    30, 50, 300);


//subjet
    TH1F *h_SJ_etad_sdLL = new TH1F(("h_SJ_etad_sd_" + tagLL).c_str(), "SubJet #eta difference L; #eta; Events", 20,-1,1);
    TH1F *h_SJ_etad_sdTT = new TH1F(("h_SJ_etad_sd_" + tagTT).c_str(), "SubJet #eta difference T; #eta; Events", 20,-1,1);
    TH1F *h_SJ_etadLL = new TH1F(("h_SJ_etad_" + tagLL).c_str(), "SubJet #eta difference L; #eta; Events", 30,-1,1);
    TH1F *h_SJ_etadTT = new TH1F(("h_SJ_etad_" + tagTT).c_str(), "SubJet #eta difference T; #eta; Events", 30,-1,1);
    
    TH1F *h_SJ_phid_sdLL = new TH1F(("h_SJ_phid_sd_" + tagLL).c_str(), "SubJet #phi difference L; #phi; Events", 20,-1,1);
    TH1F *h_SJ_phid_sdTT = new TH1F(("h_SJ_phid_sd_" + tagTT).c_str(), "SubJet #phi difference T; #phi; Events", 20,-1,1);
    TH1F *h_SJ_phidLL = new TH1F(("h_SJ_phid_" + tagLL).c_str(), "SubJet #phi difference L; #phi; Events", 30,-1,1);
    TH1F *h_SJ_phidTT = new TH1F(("h_SJ_phid_" + tagTT).c_str(), "SubJet #phi difference T; #phi; Events", 30,-1,1);

    TH1F *h_SJ_ptd_sdLL = new TH1F(("h_SJ_ptd_sd_" + tagLL).c_str(), "SubJet p_{T} difference L; p_{T}; Events", 30,150,1000);
    TH1F *h_SJ_ptd_sdTT = new TH1F(("h_SJ_ptd_sd_" + tagTT).c_str(), "SubJet p_{T} difference T; p_{T}; Events", 30,150,1000);
    TH1F *h_SJ_ptdLL = new TH1F(("h_SJ_ptd_" + tagLL).c_str(), "SubJet p_{T} difference L; p_{T}; Events", 30,150,1000);
    TH1F *h_SJ_ptdTT = new TH1F(("h_SJ_ptd_" + tagTT).c_str(), "SubJet p_{T} difference T; p_{T}; phi; Events", 30,150,1000);

    TH1F *h_SJ_pts1LL = new TH1F(("h_SJ_pts1_" + tagLL).c_str(), "SubJet #1 p_{T}   (CMS private work); p_{T} [GeV];Events", 30, 0, 600);
    TH1F *h_SJ_pts2LL = new TH1F(("h_SJ_pts2_" + tagLL).c_str(), "SubJet #2 p_{T}   (CMS private work); p_{T} [GeV];Events", 30, 0, 600);
    TH1F *h_SJ_pts1TT = new TH1F(("h_SJ_pts1_" + tagTT).c_str(), "SubJet #1 p_{T}   (CMS private work); p_{T} [GeV];Events", 30, 0, 600);
    TH1F *h_SJ_pts2TT = new TH1F(("h_SJ_pts2_" + tagTT).c_str(), "SubJet #2 p_{T}   (CMS private work); p_{T} [GeV];Events", 30, 0, 600);

    TH1F *h_SJ_mass1LL = new TH1F(("h_SJ_mass1LL" + tagLL).c_str(), "SubJet #1 mass   (CMS private work); mass [GeV];Events", 30, 0, 100);
    TH1F *h_SJ_mass2LL = new TH1F(("h_SJ_mass2TT" + tagLL).c_str(), "SubJet #2 mass   (CMS private work); mass [GeV];Events", 30, 0, 100);
    TH1F *h_SJ_mass1TT = new TH1F(("h_SJ_mass1LL" + tagTT).c_str(), "SubJet #1 mass   (CMS private work); mass [GeV];Events", 30, 0, 100);
    TH1F *h_SJ_mass2TT = new TH1F(("h_SJ_mass2TT" + tagTT).c_str(), "SubJet #2 mass   (CMS private work); mass [GeV];Events", 30, 0, 100);

    TH1F *h_SJ_eta1LL = new TH1F(("h_SJ_eta1LL" + tagLL).c_str(), "SubJet #1 #eta   (CMS private work); #eta;Events", 30, -4, 4);
    TH1F *h_SJ_eta2LL = new TH1F(("h_SJ_eta2LL" + tagLL).c_str(), "SubJet #2 #eta   (CMS private work); #eta;Events", 30, -4, 4);
    TH1F *h_SJ_eta1TT = new TH1F(("h_SJ_eta1TT" + tagTT).c_str(), "SubJet #1 #eta   (CMS private work); #eta;Events", 30, -4, 4);
    TH1F *h_SJ_eta2TT = new TH1F(("h_SJ_eta2TT" + tagTT).c_str(), "SubJet #2 #eta   (CMS private work); #eta;Events", 30, -4, 4);

    TH1F *h_SJ_phi1LL = new TH1F(("h_SJ_phi1LL" + tagLL).c_str(), "SubJet #1 #phi   (CMS private work); #phi;Events", 30, -4, 4);
    TH1F *h_SJ_phi2LL = new TH1F(("h_SJ_phi2LL" + tagLL).c_str(), "SubJet #2 #phi   (CMS private work); #phi;Events", 30, -4, 4);
    TH1F *h_SJ_phi1TT = new TH1F(("h_SJ_phi1TT" + tagTT).c_str(), "SubJet #1 #phi   (CMS private work); #phi;Events", 30, -4, 4);
    TH1F *h_SJ_phi2TT = new TH1F(("h_SJ_phi2TT" + tagTT).c_str(), "SubJet #2 #phi   (CMS private work); #phi;Events", 30, -4, 4);
    
//zg1
    /*TH1F *h_SJ_z1mLL= new TH1F(("h_SJ_z1mtLL" + tagLL).c_str(), "Subjets z_{g}' 60 < m_{sd} < 100; z_{g}'; Events", 30, 1, 10);
    TH1F *h_SJ_z1mTT= new TH1F(("h_SJ_z1mtTT" + tagTT).c_str(), "Subjets z_{g}' 60 < m_{sd} < 100; z_{g}'; Events", 30, 1, 10);
    TH1F *h_SJ_z1mtLL= new TH1F(("h_SJ_z1mLL" + tagLL).c_str(), "Subjets z_{g}'   60 < m_{sd} < 100, #tau_{21}<0.45; z_{g}'; Events", 30, 1, 10);
    TH1F *h_SJ_z1mtTT= new TH1F(("h_SJ_z1mTT" + tagTT).c_str(), "Subjets z_{g}'   60 < m_{sd} < 100, #tau_{21}<0.45; z_{g}'; Events", 30, 1, 10);
    TH1F *h_SJ_z1pt1LL= new TH1F(("h_SJ_z1pt1LL" + tagLL).c_str(), "Subjets z_{g}'   p_{T}>200; z_{g}'; Events", 30, 1, 10);
    TH1F *h_SJ_z1pt1TT= new TH1F(("h_SJ_z1pt1TT" + tagTT).c_str(), "Subjets z_{g}'   p_{T}>200; z_{g}'; Events", 30, 1, 10);
    TH1F *h_SJ_z1pt2LL= new TH1F(("h_SJ_z1pt2LL" + tagLL).c_str(), "Subjets z_{g}'   p_{T}>230; z_{g}'; Events", 30, 1, 10);
    TH1F *h_SJ_z1pt2TT= new TH1F(("h_SJ_z1pt2TT" + tagTT).c_str(), "Subjets z_{g}'   p_{T}>230; z_{g}'; Events", 30, 1, 10);
    TH1F *h_SJ_z1pt3LL= new TH1F(("h_SJ_z1pt3LL" + tagLL).c_str(), "Subjets z_{g}'   p_{T}>270; z_{g}'; Events", 30, 1, 10);
    TH1F *h_SJ_z1pt3TT= new TH1F(("h_SJ_z1pt3TT" + tagTT).c_str(), "Subjets z_{g}'   p_{T}>270; z_{g}'; Events", 30, 1, 10);
    TH1F *h_SJ_z1pt4LL= new TH1F(("h_SJ_z1pt4LL" + tagLL).c_str(), "Subjets z_{g}'   p_{T}>300; z_{g}'; Events", 30, 1, 10);
    TH1F *h_SJ_z1pt4TT= new TH1F(("h_SJ_z1pt4TT" + tagTT).c_str(), "Subjets z_{g}'   p_{T}>300; z_{g}'; Events", 30, 1, 10);*/
//zg2
   /* TH1F *h_SJ_z2mtLL= new TH1F(("h_SJ_z2mtLL" + tagLL).c_str(), "Subjets z_{g}   60 < m_{sd} < 100, #tau_{21}<0.45; z_{g}'; Events", 20, 0.1, 0.5);
    TH1F *h_SJ_z2mtTT= new TH1F(("h_SJ_z2mtTT" + tagTT).c_str(), "Subjets z_{g}   60 < m_{sd} < 100, #tau_{21}<0.45; z_{g}'; Events", 20, 0.1, 0.5);
    TH1F *h_SJ_z2mLL= new TH1F(("h_SJ_z2mLL" + tagLL).c_str(), "Subjets z_{g}   60 < m_{sd} < 100; z_{g}; Events", 20, 0.1, 0.5);
    TH1F *h_SJ_z2mTT= new TH1F(("h_SJ_z2mTT" + tagTT).c_str(), "Subjets z_{g}   60 < m_{sd} < 100; z_{g}; Events", 20, 0.1, 0.5);

    TH1F *h_SJ_z2pt1LL= new TH1F(("h_SJ_z2pt1LL" + tagLL).c_str(), "Subjets z_{g}   p_{T}>200; z_{g}; Events", 30, 0.1, 0.5);
    TH1F *h_SJ_z2pt1TT= new TH1F(("h_SJ_z2pt1TT" + tagTT).c_str(), "Subjets z_{g}   p_{T}>200; z_{g}; Events", 30, 0.1, 0.5);
    TH1F *h_SJ_z2pt2LL= new TH1F(("h_SJ_z2pt2LL" + tagLL).c_str(), "Subjets z_{g}   p_{T}>230; z_{g}; Events", 30, 0.1, 0.5);
    TH1F *h_SJ_z2pt2TT= new TH1F(("h_SJ_z2pt2TT" + tagTT).c_str(), "Subjets z_{g}   p_{T}>230; z_{g}; Events", 30, 0.1, 0.5);
    TH1F *h_SJ_z2pt3LL= new TH1F(("h_SJ_z2pt3LL" + tagLL).c_str(), "Subjets z_{g}   p_{T}>270; z_{g}; Events", 30, 0.1, 0.5);
    TH1F *h_SJ_z2pt3TT= new TH1F(("h_SJ_z2pt3TT" + tagTT).c_str(), "Subjets z_{g}   p_{T}>270; z_{g}; Events", 30, 0.1, 0.5);
    TH1F *h_SJ_z2pt4LL= new TH1F(("h_SJ_z2pt4LL" + tagLL).c_str(), "Subjets z_{g}   p_{T}>200; z_{g}; Events", 30, 0.1, 0.5);
    TH1F *h_SJ_z2pt4TT= new TH1F(("h_SJ_z2pt4TT" + tagTT).c_str(), "Subjets z_{g}   p_{T}>300; z_{g}; Events", 30, 0.1, 0.5);*/

//ptheta
  /*TH1F *h_SJ_pemtLL= new TH1F(("h_SJ_pemtLL" + tagLL).c_str(), "Subjets p_{#theta} ; p_{#theta}; Events", 30, 0, 0.9);
    TH1F *h_SJ_pemtTT= new TH1F(("h_SJ_pemtTT" + tagTT).c_str(), "Subjets p_{#theta}; p_{#theta}; Events", 30, 0, 0.9);
    TH1F *h_SJ_pemLL= new TH1F(("h_SJ_pemLL" + tagLL).c_str(), "Subjets p_{#theta}; p_{#theta}; Events", 30, 0, 0.9);
    TH1F *h_SJ_pemTT= new TH1F(("h_SJ_pemTT" + tagTT).c_str(), "Subjets p_{#theta}; p_{#theta}; Events", 30, 0, 0.9);
    TH1F *h_SJ_pept1LL= new TH1F(("h_SJ_pept1LL" + tagLL).c_str(), "Subjets p_{#theta}   p_{T}>200; p_{#theta}; Events", 30, 0, 0.9);
    TH1F *h_SJ_pept1TT= new TH1F(("h_SJ_pept1TT" + tagTT).c_str(), "Subjets p_{#theta}   p_{T}>200; p_{#theta}; Events", 30, 0, 0.9);
    TH1F *h_SJ_pept2LL= new TH1F(("h_SJ_pept2LL" + tagLL).c_str(), "Subjets p_{#theta}   p_{T}>230; p_{#theta}; Events", 30, 0, 0.9);
    TH1F *h_SJ_pept2TT= new TH1F(("h_SJ_pept2TT" + tagTT).c_str(), "Subjets p_{#theta}   p_{T}>230; p_{#theta}; Events", 30, 0, 0.9);
    TH1F *h_SJ_pept3LL= new TH1F(("h_SJ_pept3LL" + tagLL).c_str(), "Subjets p_{#theta}   p_{T}>270; p_{#theta}; Events", 30, 0, 0.9);
    TH1F *h_SJ_pept3TT= new TH1F(("h_SJ_pept3TT" + tagTT).c_str(), "Subjets p_{#theta}   p_{T}>270; p_{#theta}; Events", 30, 0, 0.9);
    TH1F *h_SJ_pept4LL= new TH1F(("h_SJ_pept4LL" + tagLL).c_str(), "Subjets p_{#theta}   p_{T}>300; p_{#theta}; Events", 30, 0, 0.9);
    TH1F *h_SJ_pept4TT= new TH1F(("h_SJ_pept4TT" + tagTT).c_str(), "Subjets p_{#theta}   p_{T}>300; p_{#theta}; Events", 30, 0, 0.9);*/



    
 


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
                h_SJ_mass1LL->Fill(SubJet_mass[FatJet_subJetIdx1[fj]]);
                h_SJ_mass2LL->Fill(SubJet_mass[FatJet_subJetIdx2[fj]]);
                h_SJ_eta1LL->Fill(SubJet_eta[FatJet_subJetIdx1[fj]]);
                h_SJ_eta2LL->Fill(SubJet_eta[FatJet_subJetIdx2[fj]]);
                h_SJ_phi1LL->Fill(SubJet_phi[FatJet_subJetIdx1[fj]]);
                h_SJ_phi2LL->Fill(SubJet_phi[FatJet_subJetIdx2[fj]]);
                
                /*double max_pt = std::max(SubJet_pt[FatJet_subJetIdx1[fj]], SubJet_pt[FatJet_subJetIdx2[fj]]);
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

                if(FatJet_pt[fj]>200){
                    double max_pt = std::max(SubJet_pt[FatJet_subJetIdx1[fj]], SubJet_pt[FatJet_subJetIdx2[fj]]);
                    double min_pt = std::min(SubJet_pt[FatJet_subJetIdx1[fj]], SubJet_pt[FatJet_subJetIdx2[fj]]);
                    h_SJ_z1pt1LL->Fill(max_pt/min_pt);
                    h_SJ_z2pt1LL->Fill(min_pt/(SubJet_pt[FatJet_subJetIdx1[fj]] + SubJet_pt[FatJet_subJetIdx2[fj]]));
                    
                    double E1, E2;
                    E1 = std::sqrt(SubJet_pt[FatJet_subJetIdx1[fj]]*SubJet_pt[FatJet_subJetIdx1[fj]] + SubJet_mass[FatJet_subJetIdx1[fj]]*SubJet_mass[FatJet_subJetIdx1[fj]]);
                    E2 = std::sqrt(SubJet_pt[FatJet_subJetIdx2[fj]]*SubJet_pt[FatJet_subJetIdx2[fj]] + SubJet_mass[FatJet_subJetIdx2[fj]]*SubJet_mass[FatJet_subJetIdx2[fj]]);
                    h_SJ_pept1LL->Fill(std::abs(E1-E2)/FatJet_pt[fj]);
                }
                if(FatJet_pt[fj]>230){
                    double max_pt = std::max(SubJet_pt[FatJet_subJetIdx1[fj]], SubJet_pt[FatJet_subJetIdx2[fj]]);
                    double min_pt = std::min(SubJet_pt[FatJet_subJetIdx1[fj]], SubJet_pt[FatJet_subJetIdx2[fj]]);
                    h_SJ_z1pt2LL->Fill(max_pt/min_pt);
                    h_SJ_z2pt2LL->Fill(min_pt/(SubJet_pt[FatJet_subJetIdx1[fj]] + SubJet_pt[FatJet_subJetIdx2[fj]]));

                    double E1, E2;
                    E1 = std::sqrt(SubJet_pt[FatJet_subJetIdx1[fj]]*SubJet_pt[FatJet_subJetIdx1[fj]] + SubJet_mass[FatJet_subJetIdx1[fj]]*SubJet_mass[FatJet_subJetIdx1[fj]]);
                    E2 = std::sqrt(SubJet_pt[FatJet_subJetIdx2[fj]]*SubJet_pt[FatJet_subJetIdx2[fj]] + SubJet_mass[FatJet_subJetIdx2[fj]]*SubJet_mass[FatJet_subJetIdx2[fj]]);
                    h_SJ_pept2LL->Fill(std::abs(E1-E2)/FatJet_pt[fj]);
                }
                if(FatJet_pt[fj]>270){
                    double max_pt = std::max(SubJet_pt[FatJet_subJetIdx1[fj]], SubJet_pt[FatJet_subJetIdx2[fj]]);
                    double min_pt = std::min(SubJet_pt[FatJet_subJetIdx1[fj]], SubJet_pt[FatJet_subJetIdx2[fj]]);
                    h_SJ_z1pt3LL->Fill(max_pt/min_pt);
                    h_SJ_z2pt3LL->Fill(min_pt/(SubJet_pt[FatJet_subJetIdx1[fj]] + SubJet_pt[FatJet_subJetIdx2[fj]]));

                    double E1, E2;
                    E1 = std::sqrt(SubJet_pt[FatJet_subJetIdx1[fj]]*SubJet_pt[FatJet_subJetIdx1[fj]] + SubJet_mass[FatJet_subJetIdx1[fj]]*SubJet_mass[FatJet_subJetIdx1[fj]]);
                    E2 = std::sqrt(SubJet_pt[FatJet_subJetIdx2[fj]]*SubJet_pt[FatJet_subJetIdx2[fj]] + SubJet_mass[FatJet_subJetIdx2[fj]]*SubJet_mass[FatJet_subJetIdx2[fj]]);
                    h_SJ_pept3LL->Fill(std::abs(E1-E2)/FatJet_pt[fj]);
                }
                if(FatJet_pt[fj]>300){
                    double max_pt = std::max(SubJet_pt[FatJet_subJetIdx1[fj]], SubJet_pt[FatJet_subJetIdx2[fj]]);
                    double min_pt = std::min(SubJet_pt[FatJet_subJetIdx1[fj]], SubJet_pt[FatJet_subJetIdx2[fj]]);
                    h_SJ_z1pt4LL->Fill(max_pt/min_pt);
                    h_SJ_z2pt4LL->Fill(min_pt/(SubJet_pt[FatJet_subJetIdx1[fj]] + SubJet_pt[FatJet_subJetIdx2[fj]]));

                    double E1, E2;
                    E1 = std::sqrt(SubJet_pt[FatJet_subJetIdx1[fj]]*SubJet_pt[FatJet_subJetIdx1[fj]] + SubJet_mass[FatJet_subJetIdx1[fj]]*SubJet_mass[FatJet_subJetIdx1[fj]]);
                    E2 = std::sqrt(SubJet_pt[FatJet_subJetIdx2[fj]]*SubJet_pt[FatJet_subJetIdx2[fj]] + SubJet_mass[FatJet_subJetIdx2[fj]]*SubJet_mass[FatJet_subJetIdx2[fj]]);
                    h_SJ_pept4LL->Fill(std::abs(E1-E2)/FatJet_pt[fj]);
                }*/

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
                h_SJ_mass1TT->Fill(SubJet_mass[FatJet_subJetIdx1[fj]]);
                h_SJ_mass2TT->Fill(SubJet_mass[FatJet_subJetIdx2[fj]]);
                h_SJ_eta1TT->Fill(SubJet_eta[FatJet_subJetIdx1[fj]]);
                h_SJ_eta2TT->Fill(SubJet_eta[FatJet_subJetIdx2[fj]]);
                h_SJ_phi1TT->Fill(SubJet_phi[FatJet_subJetIdx1[fj]]);
                h_SJ_phi2TT->Fill(SubJet_phi[FatJet_subJetIdx2[fj]]);
                

                /*double max_pt = std::max(SubJet_pt[FatJet_subJetIdx1[fj]], SubJet_pt[FatJet_subJetIdx2[fj]]);
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
                }*/
                
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

    //ISTOGRAMMI PER pt
    TCanvas *c_pt_tot = new TCanvas("c_pt_tot", "Consituents p_{T} LL vs TT", 2000, 1100);
    std::vector<TH1*> pt;
    pt.push_back(h_GFJ_ptLL);
    pt.push_back(h_GFJ_ptTT);
    pt.push_back(h_FJ_ptLL);
    pt.push_back(h_FJ_ptTT);
    Makecanva2x1(c_pt_tot, pt, kAzure-4);
    TLegend *leg7 = new TLegend(0.5, 0.65, 0.7, 0.8);
    leg7->AddEntry(h_GFJ_ptLL, "L pol", "f");
    leg7->AddEntry(h_GFJ_ptTT, "T pol", "f");
    leg7->SetBorderSize(0);
    leg7->SetFillStyle(0);
    leg7->SetMargin(0.25);
    leg7->Draw();
    c_pt_tot->SaveAs("plots_pt/c_pt_tot.pdf");
    


    //ISTOGRAMMI PER ETA
    TCanvas *c_eta_tot = new TCanvas("c_eta_tot", "Consituents #eta: pseudorapirity LL vs TT", 2000, 1100);
    std::vector<TH1*> eta;
    eta.push_back(h_GFJ_etaTT);
    eta.push_back(h_GFJ_etaLL);
    eta.push_back(h_FJ_etaTT);
    eta.push_back(h_FJ_etaLL);
    Makecanva2x1(c_eta_tot, eta, kAzure-4);
    TLegend *leg6 = new TLegend(0.15, 0.65, 0.35, 0.8);
    leg6->AddEntry(h_GFJ_etaTT, "T pol", "f");
    leg6->AddEntry(h_GFJ_etaLL, "L pol", "f");
    leg6->SetBorderSize(0);
    leg6->SetFillStyle(0);
    leg6->SetMargin(0.25);
    leg6->Draw();
    c_eta_tot->SaveAs("plots_eta/c_eta_tot.pdf");




    //ISTOGRAMMI PER PHI
    TCanvas *c_phi_tot = new TCanvas("c_phi_tot", "Consituents #phi: LL vs TT", 2000, 1100);
    std::vector<TH1*> phi;
    phi.push_back(h_GFJ_phiTT);
    phi.push_back(h_GFJ_phiLL);
    phi.push_back(h_FJ_phiTT);
    phi.push_back(h_FJ_phiLL);
    Makecanva2x1(c_phi_tot, phi, kAzure-4);
    TLegend *leg10= new TLegend(0.45, 0.55, 0.6, 0.7);
    leg10->AddEntry(h_GFJ_phiLL, "L pol", "f");
    leg10->AddEntry(h_GFJ_phiTT, "T pol", "f");
    leg10->SetBorderSize(0);
    leg10->Draw();
    c_phi_tot->SaveAs("plots_phi/c_phi_tot.pdf");
    


    
    //ISTOGRAMMI PER LE MASSE
    TCanvas *c_mass_tot_pol = new TCanvas("c_mass_tot_pol", "Consituents Mass LL vs TT", 2000, 1100);
    std::vector<TH1*> mass;
    mass.push_back(h_GFJ_massLL);
    mass.push_back(h_GFJ_massTT);
    mass.push_back(h_FJ_massLL);
    mass.push_back(h_FJ_massTT);
    Makecanva2x1(c_mass_tot_pol, mass, kAzure-4);
    TLegend *leg11 = new TLegend(0.6, 0.65, 0.8, 0.8);
    leg11->AddEntry(h_GFJ_massLL, "L pol", "f");
    leg11->AddEntry(h_GFJ_massTT, "T pol", "f");
    leg11->SetBorderSize(0);
    leg11->SetFillStyle(0);
    leg11->SetMargin(0.25);
    leg11->Draw();
    c_mass_tot_pol->SaveAs("plots_m/c_mass_tot_pol.pdf");



    //Soft Drop
    TCanvas *c_FJ_msd_pol = new TCanvas("c_FJ_msd_pol", "FatJet SoftDrop mass LL vs TT", 2000, 1100);
    std::vector<TH1*> softdrop;
    h_FJ_massTT->SetTitleSize(0.09); 
    h_FJ_massLL->SetTitleSize(0.09);
    h_FJ_massTT->SetTitle("Fatjet mass / soft drop (T)   (CMS private work)");
    h_FJ_massLL->SetTitle("Fatjet mass / soft drop (L)   (CMS private work)");
    softdrop.push_back(h_FJ_massTT);
    softdrop.push_back(h_FJ_msdTT);
    softdrop.push_back(h_FJ_massLL);
    softdrop.push_back(h_FJ_msdLL);
    Makecanva2x1(c_FJ_msd_pol, softdrop, kBlue-4);
    TLegend *leg9 = new TLegend(0.5, 0.6, 0.8, 0.75);
    leg9->AddEntry(h_FJ_msdLL, "SoftDrop mass", "f");
    leg9->AddEntry(h_FJ_massLL, "Mass", "f");
    leg9->SetBorderSize(0);
    leg9->SetFillStyle(0);
    leg9->SetMargin(0.25);
    leg9->Draw();
    c_FJ_msd_pol->SaveAs("c_FJ_msd_pol.pdf");



    //ISTOGRAMMI TAU 1,2, 3 
    TCanvas *c_tau_pol = new TCanvas("c_tau", "FatJet tau", 3000, 1200);
    std::vector<TH1*> tau_pol;
    tau_pol.push_back(h_tau1TT);
    tau_pol.push_back(h_tau1LL);
    tau_pol.push_back(h_tau2TT);
    tau_pol.push_back(h_tau2LL);
    tau_pol.push_back(h_tau3TT);
    tau_pol.push_back(h_tau3LL);
    Makecanva3x1(c_tau_pol, tau_pol, kAzure-4);
    TLegend *leg1 = new TLegend(0.55, 0.7, 0.75, 0.8);
    leg1->AddEntry(h_tau3TT, "T pol", "f");
    leg1->AddEntry(h_tau3LL, "L pol", "f");
    leg1->SetBorderSize(0);
    leg1->SetFillStyle(0);
    leg1->SetMargin(0.25);
    leg1->Draw();
    c_tau_pol->SaveAs("plots_tau/c_tau_pol.pdf");



    //ISTOGRAMMI Tau21, Tau32 confronto con eventi "W" 
    TCanvas *c_tau_ratio_pol = new TCanvas("c_tau_ratio_pol", "FatJet tau ratios", 2000, 2000);

    std::vector<TH1*> tau_ratio_pol;
    tau_ratio_pol.push_back(h_tau21TT);
    tau_ratio_pol.push_back(h_tau21LL);
    tau_ratio_pol.push_back(h_tau32TT);
    tau_ratio_pol.push_back(h_tau32LL);
    tau_ratio_pol.push_back(h_tau21mTT);
    tau_ratio_pol.push_back(h_tau21mLL);
    tau_ratio_pol.push_back(h_tau32mTT);
    tau_ratio_pol.push_back(h_tau32mLL);
    Makecanva2x2(c_tau_ratio_pol, tau_ratio_pol, kAzure-4);
    TLegend *leg2 = new TLegend(0.25, 0.7, 0.45, 0.85);
    leg2->AddEntry(h_tau32mTT, "T pol", "f");
    leg2->AddEntry(h_tau32mLL, "L pol", "f");
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->SetMargin(0.25);
    leg2->Draw();
    c_tau_ratio_pol->SaveAs("plots_tau/c_tau_ratio_pol.pdf");

    
    //SCATTER PLOT TAU - MASSA SD
    TCanvas *c_tau_msd = new TCanvas("c_tau_msd","Tau vs msd", 2000,1600);
    c_tau_msd->Divide(2,2);

    c_tau_msd->cd(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.15);
    h_tau21_msdTT->Scale(1.0/ h_tau21_msdTT->Integral());
    h_tau21_msdTT->GetXaxis()->SetTitleSize(0.05);
    h_tau21_msdTT->GetYaxis()->SetTitleSize(0.05);
    h_tau21_msdTT->Draw("COLZ"); 


    c_tau_msd->cd(2);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.15);
    h_tau21_msdLL->Scale(1.0/ h_tau21_msdLL->Integral());
    h_tau21_msdLL->GetXaxis()->SetTitleSize(0.05);
    h_tau21_msdLL->GetYaxis()->SetTitleSize(0.05);
    h_tau21_msdLL->Draw("colz");

    c_tau_msd->cd(3);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.15);
    h_tau32_msdTT->Scale(1.0/ h_tau32_msdTT->Integral());
    h_tau32_msdTT->GetXaxis()->SetTitleSize(0.05);
    h_tau32_msdTT->GetYaxis()->SetTitleSize(0.05);
    h_tau32_msdTT->Draw("COLZ"); 

    c_tau_msd->cd(4);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.15);
    h_tau32_msdLL->Scale(1.0/ h_tau32_msdLL->Integral());
    h_tau21_msdLL->GetXaxis()->SetTitleSize(0.05);
    h_tau21_msdLL->GetYaxis()->SetTitleSize(0.05);
    h_tau32_msdLL->Draw("colz");

    c_tau_msd->SaveAs("plots_tau/c_tau_msd.pdf");
    


    /*ISTOGRAMMI per SubJet*/
//eta
/*
    TCanvas *c_SJ_etad = new TCanvas("c_SJ_etad", "SubJet #eta difference: LL vs TT", 800, 600);
    std::vector<TH1F* SJ_eta;
    SJ_eta.push_back(h_SJ_eta)

    TLegend *leg4 = new TLegend(0.2, 0.65, 0.35, 0.8);
    leg4->AddEntry(h_SJ_etadTT, "TT", "l");
    leg4->AddEntry(h_SJ_etadLL, "LL", "l");
    leg4->Draw();

    c_SJ_etad->SaveAs("plots_sj/c_SJ_etad.pdf");


//phi
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
    TCanvas *c_SJ_ptd = new TCanvas("c_SJ_ptd", "SubJet p_{T} difference: LL vs TT", 1200, 600);
    std::vector<TH1*> tau_pol;
    tau_pol.push_back(h_tau1TT);
    tau_pol.push_back(h_tau1LL);
    tau_pol.push_back(h_tau2TT);
    tau_pol.push_back(h_tau2LL);
    tau_pol.push_back(h_tau3TT);
    tau_pol.push_back(h_tau3LL);
    Makecanva3x1(c_SJ_ptd, SJ_ptd, kAzure-4);
    c_SJ_ptd->Divide(2,1);

    c_SJ_ptd->cd(1);
    h_SJ_ptd_sdTT->Draw();
    h_SJ_ptd_sdTT->SetTitle("SubJet p_{T} difference: 60 < m_sd < 100");
    h_SJ_ptd_sdLL->Draw("same");
    h_SJ_ptd_sdLL->SetLineColor(kRed);


    c_SJ_ptd->cd(2);
    h_SJ_ptdTT->Draw();
    h_SJ_ptdTT->SetTitle("SubJet p_{T} difference: all masses");
    h_SJ_ptdLL->Draw("same");
    h_SJ_ptdLL->SetLineColor(kRed);

    TLegend *leg8 = new TLegend(0.2, 0.65, 0.35, 0.8);
    leg8->AddEntry(h_SJ_ptd_sdTT, "TT", "l");
    leg8->AddEntry(h_SJ_ptd_sdLL, "LL", "l");
    leg8->Draw();

    c_SJ_ptd->SaveAs("plots_sj/c_SJ_ptd.pdf");*/




    /*TCanvas *c_SJ_pt_pol = new TCanvas("c_SJ_pts", "Subjets p_{T} 1 vs 2", 1200, 600);
    std::vector<TH1*> SJ_pt;
    SJ_pt.push_back(h_SJ_pts1LL);
    SJ_pt.push_back(h_SJ_pts1TT);
    SJ_pt.push_back(h_SJ_pts2LL);
    SJ_pt.push_back(h_SJ_pts2TT);
    Makecanva2x1(c_SJ_pt_pol, SJ_pt, kAzure+6);
    TLegend *leg19 = new TLegend(0.6, 0.7, 0.75, 0.8);
    leg19->AddEntry(h_SJ_pts2TT, "T pol", "f");
    leg19->AddEntry(h_SJ_pts2LL, "L pol", "f");
    leg19->SetBorderSize(0);
    leg19->SetFillStyle(0);
    leg19->SetMargin(0.25);
    leg19->Draw();
    c_SJ_pt_pol->SaveAs("plots_sj/c_SJ_pt_pol.pdf");


    TCanvas *c_SJ_mass_pol = new TCanvas("c_SJ_mass_pol", "Subjets mass 1 vs 2", 1200, 600);
    std::vector<TH1*> SJ_mass;
    SJ_mass.push_back(h_SJ_mass1TT);
    SJ_mass.push_back(h_SJ_mass1LL);
    SJ_mass.push_back(h_SJ_mass2TT);
    SJ_mass.push_back(h_SJ_mass2LL);
    Makecanva2x1(c_SJ_mass_pol, SJ_mass, kAzure+6);
    TLegend *leg20 = new TLegend(0.6, 0.7, 0.75, 0.8);
    leg20->AddEntry(h_SJ_mass2TT, "T pol", "f");
    leg20->AddEntry(h_SJ_mass2LL, "L pol", "f");
    leg20->SetBorderSize(0);
    leg20->SetFillStyle(0);
    leg20->SetMargin(0.25);
    leg20->Draw();
    c_SJ_mass_pol->SaveAs("plots_sj/c_SJ_mass_pol.pdf");


    TCanvas *c_SJ_eta_pol = new TCanvas("c_SJ_eta_pol", "Subjets #eta 1 vs 2", 1200, 600);
    std::vector<TH1*> SJ_eta;
    SJ_eta.push_back(h_SJ_eta1TT);
    SJ_eta.push_back(h_SJ_eta1LL);
    SJ_eta.push_back(h_SJ_eta2TT);
    SJ_eta.push_back(h_SJ_eta2LL);
    Makecanva2x1(c_SJ_eta_pol, SJ_eta, kAzure+6);
    TLegend *leg21 = new TLegend(0.45, 0.4, 0.6, 0.5);
    leg21->AddEntry(h_SJ_eta2TT, "T pol", "f");
    leg21->AddEntry(h_SJ_eta2LL, "L pol", "f");
    leg21->SetBorderSize(0);
    leg21->Draw();
    c_SJ_eta_pol->SaveAs("plots_sj/c_SJ_eta_pol.pdf");


    TCanvas *c_SJ_phi_pol = new TCanvas("c_SJ_phi_pol", "Subjets #phi 1 vs 2", 1200, 600);
    std::vector<TH1*> SJ_phi;
    SJ_phi.push_back(h_SJ_phi1TT);
    SJ_phi.push_back(h_SJ_phi1LL);
    SJ_phi.push_back(h_SJ_phi2TT);
    SJ_phi.push_back(h_SJ_phi2LL);
    Makecanva2x1(c_SJ_phi_pol, SJ_phi, kAzure+6);
    TLegend *leg22 = new TLegend(0.45, 0.5, 0.6, 0.6);
    leg22->AddEntry(h_SJ_phi2TT, "T pol", "f");
    leg22->AddEntry(h_SJ_phi2LL, "L pol", "f");
    leg22->SetBorderSize(0);
    leg22->Draw();
    c_SJ_phi_pol->SaveAs("plots_sj/c_SJ_phi_pol.pdf");*/
    

/*zg1
    TCanvas *c_z1m = new TCanvas("c_z1m", "SubJet z_{g}: L vs T", 1200, 600);
    std::vector<TH1*> z1m;
    z1m.push_back(h_SJ_z1mTT);
    z1m.push_back(h_SJ_z1mLL);
    z1m.push_back(h_SJ_z1mtTT);
    z1m.push_back(h_SJ_z1mtLL);
    Makecanva2x1(c_z1m, z1m, kAzure-4);
    TLegend *leg12 = new TLegend(0.5, 0.7, 0.65, 0.8);
    leg12->AddEntry(h_SJ_z1mtTT, "T pol", "l");
    leg12->AddEntry(h_SJ_z1mtLL, "L pol", "l");
    leg12->Draw();
    c_z1m->SaveAs("plots_z1/c_z1m.pdf");



    TCanvas *c_z1pt = new TCanvas("c_SJ_z1pt", "SubJet z_{g}: LL vs TT", 1600, 1200);
    std::vector<TH1*> z1pt;
    z1pt.push_back(h_SJ_z1pt1LL);
    z1pt.push_back(h_SJ_z1pt1TT);
    z1pt.push_back(h_SJ_z1pt2LL);
    z1pt.push_back(h_SJ_z1pt2TT);
    z1pt.push_back(h_SJ_z1pt3LL);
    z1pt.push_back(h_SJ_z1pt3TT);
    z1pt.push_back(h_SJ_z1pt4LL);
    z1pt.push_back(h_SJ_z1pt4TT);
    Makecanva2x2(c_z1pt, z1pt, kAzure-4);
    TLegend *leg14 = new TLegend(0.4, 0.7, 0.55, 0.8);
    leg14->AddEntry(h_SJ_z1pt1TT, "T pol", "l");
    leg14->AddEntry(h_SJ_z1pt1LL, "L pol", "l");
    leg14->Draw();
    c_z1pt->SaveAs("plots_z1/c_z1pt.pdf");

//zg2
    TCanvas *c_z2m = new TCanvas("c_SJ_z2m", "SubJet z_{g}: LL vs TT", 1200, 600);
    std::vector<TH1*> z2m;
    z2m.push_back(h_SJ_z2mLL);
    z2m.push_back(h_SJ_z2mTT);
    z2m.push_back(h_SJ_z2mtLL);
    z2m.push_back(h_SJ_z2mtTT);
    Makecanva2x1(c_z2m, z2m, kAzure-4);
    TLegend *leg15 = new TLegend(0.2, 0.75, 0.4, 0.85);
    leg15->AddEntry(h_SJ_z2mTT, "T pol", "l");
    leg15->AddEntry(h_SJ_z2mLL, "L pol", "l");
    leg15->Draw();
    c_z2m->SaveAs("plots_z1/c_z2m.pdf");




    TCanvas *c_z2pt = new TCanvas("c_z2pt", "SubJet z_{g}: L vs T", 1600, 1200);
    std::vector<TH1*> z2pt;
    z2pt.push_back(h_SJ_z2pt1LL);
    z2pt.push_back(h_SJ_z2pt1TT);
    z2pt.push_back(h_SJ_z2pt2LL);
    z2pt.push_back(h_SJ_z2pt2TT);
    z2pt.push_back(h_SJ_z2pt3LL);
    z2pt.push_back(h_SJ_z2pt3TT);
    z2pt.push_back(h_SJ_z2pt4LL);
    z2pt.push_back(h_SJ_z2pt4TT);
    Makecanva2x2(c_z2pt, z2pt, kAzure-4);
    TLegend *leg16 = new TLegend(0.6, 0.7, 0.75, 0.8);
    leg16->AddEntry(h_SJ_z2pt1LL, "L pol", "l");
    leg16->AddEntry(h_SJ_z2pt1TT, "T pol", "l");
    leg16->Draw();
    c_z2pt->SaveAs("plots_z1/c_z2pt.pdf");

//ptheta
    TCanvas *c_pthetam_pol = new TCanvas("c_pthetam_pol", "SubJet p_{#theta}: L vs T", 1200, 600);
    std::vector<TH1*> pthetam_pol;
    pthetam_pol.push_back(h_SJ_pemLL);
    pthetam_pol.push_back(h_SJ_pemTT);
    pthetam_pol.push_back(h_SJ_pemtLL);
    pthetam_pol.push_back(h_SJ_pemtTT);
    Makecanva2x1(c_pthetam_pol, pthetam_pol, kAzure-4);
    TLegend *leg17 = new TLegend(0.7, 0.8, 0.9, 0.9);
    leg17->AddEntry(h_SJ_pemTT, "T pol", "l");
    leg17->AddEntry(h_SJ_pemLL, "L pol", "l");
    leg17->Draw();
    c_pthetam_pol->SaveAs("plots_ptheta/c_pthetam_pol.pdf");    




    TCanvas *c_pthetapt_pol = new TCanvas("c_pthetapt_pol", "SubJet z_{g}: L vs T", 1600, 1200);
    std::vector<TH1*> pthetapt_pol;
    pthetapt_pol.push_back(h_SJ_pept1LL);
    pthetapt_pol.push_back(h_SJ_pept1TT);
    pthetapt_pol.push_back(h_SJ_pept2LL);
    pthetapt_pol.push_back(h_SJ_pept2TT);
    pthetapt_pol.push_back(h_SJ_pept3LL);
    pthetapt_pol.push_back(h_SJ_pept3TT);
    pthetapt_pol.push_back(h_SJ_pept4LL);
    pthetapt_pol.push_back(h_SJ_pept4TT);
    Makecanva2x2(c_pthetapt_pol, pthetapt_pol, kAzure-4);
    TLegend *leg18 = new TLegend(0.2, 0.75, 0.35, 0.85);
    leg18->AddEntry(h_SJ_pept1LL, "L pol", "l");
    leg18->AddEntry(h_SJ_pept1TT, "T pol", "l");
    leg18->Draw();
    c_pthetapt_pol->SaveAs("plots_ptheta/c_pthetapt_pol.pdf");*/



    f->Write();
    f->Close();


    gROOT->ProcessLine(".q");



}
