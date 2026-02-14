#include <TTree.h>
#include <TLatex.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <iostream>
#include <cmath>
#include <limits>
#include <string>
#include <vector>
#include "samples.h"
#include "canva.cpp"
#include <TVector3.h>
#include <TVector2.h>
#include "TLorentzVector.h"



std::string tag_cutflow = "diffMassAK4_AK8_loose_";

bool ssWWTT = true;
bool ssWWLL = true;


std::map<std::string, double> good_ele, good_muon, good_jet, good_fatjet;
std::map<std::string, double> V_to_lep, V_to_had, VBS_presel;


good_ele["pt"]  = 15.0;
good_ele["eta"] = 5.0;
good_ele["iso"] = 0.06;

good_muon["pt"]  = 15.0;
good_muon["eta"] = 5.0;
good_muon["iso"] = 0.25;

good_jet["pt"]  = 30.0;
good_jet["eta"] = 5.0;
good_jet["wp_btag"] = 0.0479; // loose working point
//good_jet["wp_btag"] = 0.2431;

good_fatjet["pt"]  = 30.0;
good_fatjet["eta"] = 5.0;
good_fatjet["msd"] = 25.0;

V_to_lep["ele_pt"]  = 30.0;
V_to_lep["muon_pt"] = 30.0;
V_to_lep["MET_pt"]  = 30.0;

V_to_had["m_sd_min"] = 70.0;
V_to_had["m_sd_max"] = 115.0;
V_to_had["mass_min"] = 60;
V_to_had["mass_max"] = 110;
V_to_had["tau21"]    = 0.45;

VBS_presel["m_jj"]       = 500.0;
VBS_presel["pt_leading"] = 50.0;
VBS_presel["deltaEta"]   = 2.5;


// --- Map bin -> description ---
std::map<int, std::string> cutDescriptions = {
    {1, "At least 1 e or 1 #mu"},
    {2, Form("good ele: p_{T} > %.0f , |#eta|< %.1f, isolation WP80 , iso < %.3f, no good #mu", good_ele["pt"], good_ele["eta"], good_ele["iso"])},
    {3, Form("good muon: p_{T} > %.0f , |#eta|< %.1f , tight ID , iso < %.3f, no good e", good_muon["pt"], good_muon["eta"], good_muon["iso"])},
    //{4, Form("good AK4: p_{T}>%.0f , #Delta R(lep,jet)>0.4, |#eta|<%.1f", good_jet["pt"], good_jet["eta"])},
    {4, Form("good AK8: p_{T}>%.0f , m_{SD}>%.0f , #Delta R(l,FatJet)>0.8", good_fatjet["pt"], good_fatjet["msd"])},
    //{6, ""},
    {5, ""},
    {6, Form("VBS kin : |#Delta #eta| > %.1f , m_{jj} > %.2f , p^{lead.}_{T} > %.0f", VBS_presel["deltaEta"], VBS_presel["m_jj"], VBS_presel["pt_leading"])},
    {7, Form("V to lep : p^e_{T} > %.0f or p^{#mu}_{T} > %.0f , MET > %.0f GeV", V_to_lep["ele_pt"], V_to_lep["muon_pt"], V_to_lep["MET_pt"])},
    //{8, ""},
    {8, Form("V to had: %.0f < m_{sd} < %.0f , #tau_{21} < %.2f (if 1 AK8)", V_to_had["m_sd_min"], V_to_had["m_sd_max"], V_to_had["tau21"])},
    {9, ""},
    //{10, Form("b-tagging: WP: %.4f, RobustParT", good_jet["wp_btag"])},

};


const double mW = 80.4; // W boson mass

double deltaPhi(double phi1, double phi2) {
    double dphi = phi1 - phi2;
    while (dphi > M_PI) dphi -= 2*M_PI;
    while (dphi < -M_PI) dphi += 2*M_PI;
    return dphi;
}

double deltaR(double eta1, double phi1, double eta2, double phi2) {
    double dEta = eta1 - eta2;
    double dPhi = deltaPhi(phi1, phi2);
    return std::sqrt(dEta*dEta + dPhi*dPhi);
}

double deltaEta(double eta1, double eta2) { return eta1 - eta2; }

double invariantMass(double pt1, double eta1, double phi1,
                     double pt2, double eta2, double phi2,
                     double m1=0., double m2=0.) 
{
    double theta1 = 2.0 * atan(exp(-eta1));
    double theta2 = 2.0 * atan(exp(-eta2));
    double px1 = pt1*cos(phi1), py1 = pt1*sin(phi1), pz1 = pt1/tan(theta1), E1 = sqrt(px1*px1 + py1*py1 + pz1*pz1 + m1*m1);
    double px2 = pt2*cos(phi2), py2 = pt2*sin(phi2), pz2 = pt2/tan(theta2), E2 = sqrt(px2*px2 + py2*py2 + pz2*pz2 + m2*m2);

    double px = px1 + px2, py = py1 + py2, pz = pz1 + pz2, E = E1 + E2;
    return sqrt(E*E - px*px - py*py - pz*pz);
}
    const int maxJets = 30;
    const int maxFatJets = 30;
    const int maxElectrons = 30;
    const int maxMuons = 30;

    // Branch arrays
    ULong64_t event;
    Short_t FatJet_subJetIdx1[maxFatJets], FatJet_subJetIdx2[maxFatJets], FatJet_electronIdx3SJ[maxFatJets], FatJet_muonIdx3SJ[maxFatJets];
    Float_t SubJet_pt[maxFatJets], SubJet_eta[maxFatJets], SubJet_phi[maxFatJets], SubJet_mass[maxFatJets];
    Float_t Rho_fixedGridRhoFastjetCentral;
    Float_t Rho_fixedGridRhoFastjetAll;
    Float_t Jet_pt[maxJets], Jet_eta[maxJets], Jet_phi[maxJets], Jet_mass[maxJets], Jet_btagRobustParTAK4B[maxJets];
    Float_t FatJet_pt[maxFatJets], FatJet_eta[maxFatJets], FatJet_phi[maxFatJets], FatJet_mass[maxFatJets], FatJet_msoftdrop[maxFatJets], FatJet_tau1[maxFatJets], FatJet_tau2[maxFatJets];
    Float_t Electron_pt[maxElectrons], Electron_eta[maxElectrons], Electron_phi[maxElectrons], Electron_pfRelIso03_all[maxElectrons], Electron_deltaEtaSC[maxElectrons], Electron_mass[maxElectrons];
    Bool_t Electron_mvaIso_WP90[maxElectrons], Electron_mvaIso_WP80[maxElectrons];
    Float_t Muon_pt[maxMuons], Muon_eta[maxMuons], Muon_phi[maxMuons], Muon_pfRelIso04_all[maxMuons], Muon_mass[maxMuons];
    Bool_t Muon_tightId[maxMuons];
    Float_t MET_pt, MET_phi;

    const Int_t maxGenParticle = 10000;
    Float_t GenPart_pt[maxGenParticle], GenPart_eta[maxGenParticle], GenPart_phi[maxGenParticle], GenPart_mass[maxGenParticle];
    Int_t GenPart_status[maxGenParticle];
    UShort_t GenPart_statusFlags[maxGenParticle];
    Short_t GenPart_genPartIdxMother[maxGenParticle];
    Int_t GenPart_pdgId[maxGenParticle];
    Int_t nGenParticle;

    const Int_t maxSubGenJets = 200;
    Int_t nSubGenJetAK8;
    Float_t SubGenJetAK8_pt[maxSubGenJets];
    Float_t SubGenJetAK8_eta[maxSubGenJets];
    Float_t SubGenJetAK8_phi[maxSubGenJets];
    Float_t SubGenJetAK8_mass[maxSubGenJets];

    const Int_t maxGenFatJets = 100;
    Int_t nGenJetAK8;
    Float_t GenJetAK8_pt[maxGenFatJets], GenJetAK8_eta[maxGenFatJets], GenJetAK8_phi[maxGenFatJets], GenJetAK8_mass[maxGenFatJets];
    Short_t GenJetAK8_partonFlavour[maxGenFatJets]; 

    Int_t nJet, nFatJet, nElectron, nMuon;


//GenJet histograms

    TH1F *h_matchedJet_pt = new TH1F(
    "h_matchedJet_pt",
    "Matched GenJetAK8 p_{T} ;p_{T} [GeV]",
    30, 170, 1000
    );

    TH1F *h_matchedJet_mass = new TH1F(
        "h_matchedJet_mass",
        "Matched GenJetAK8 mass ;Mass [GeV]",
        30, 0, 300
    );


    TH1F *h_matchedJet_eta = new TH1F(
        "h_matchedJet_eta",
        "Matched GenJetAK8 #eta ;#eta",
        30, -4, 4
    );


    TH1F *h_matchedJet_phi = new TH1F(
        "h_matchedJet_phi",
        "Matched GenJetAK8 #phi ;#phi",
        30, -4, 4
    );
int pt_lim = 150;

    TH1F *h_zg_subjets = new TH1F("h_zg_subjets", Form("z_{g} subjets p_{t} > %d  ; z_{g} ; Events", pt_lim), 20, 0.1, 0.5);
    TH1F *h_zg_ratio = new TH1F("h_zg_ratio", Form("z_{g}' subjets p_{t} > %d  ; z_{g}' ; Events", pt_lim), 20, 1, 10);    
    TH1F *h_ptheta_subjets = new TH1F("h_ptheta_subjets", Form("p_{#theta} subjets p_{t} > %d; p_{#theta} ; Events", pt_lim), 20, 0, 0.9);

    TH1F *h_Lp = new TH1F("h_Lp1", "Lepton projection (CMS private work); L_{P}; Events", 20, 0, 1.2);
    TH1F *h_cos2D = new TH1F("h_cos2D", "cos2D; cos2D; Events", 20, -1, 1);
    TH1F *h_costheta_star = new TH1F("h_costheta_star", "h_costheta_star; cos#theta*; Events", 40, -1, 1);

    #include "func.cpp"


//_mg5+ms
void process_VBS_nanoAOD() {
    std::string tag[2] = {"ssWW_LL", "ssWW_TT"};
    std::string tag2[2] = {"LL", "TT"};
    float xsec[2] = {0.02143, 0.1789};

    TFile* fGen = TFile::Open("Gen_match.root");
        if (!fGen || fGen->IsZombie()) {
        Error("MyCode", "Impossibile aprire GenJet_match.root");
        return;
    }
    //fGen->ls();

    TH1F* h_cos2D_gen = nullptr;  
    TH1F* h_costheta_star_gen; 
    TH1F* h_Lp_gen; 
    TH1F* h_zg_gen; 
    TH1F* h_ptheta_gen; 
    TH1F* h_zg_ratio_gen;


    for(int isample=0; isample<2; isample++){

        std::vector<std::string> files;

    if(isample==0){ 

        for (auto &f : samples_ssWWLL) files.push_back(f);
        h_cos2D_gen = dynamic_cast<TH1F*>(fGen->Get("h_cos2D_L"));
        if (!h_cos2D_gen) {
            Error("MyCode", "Istogramma h_cos2D_L non trovato");
            return;
        }
        h_costheta_star_gen = dynamic_cast<TH1F*>(fGen->Get("h_costheta_star_L"));
        if (!h_costheta_star_gen) {
            Error("MyCode", "Istogramma h_costheta_star_L non trovato");
            return;
        }
        h_Lp_gen = dynamic_cast<TH1F*>(fGen->Get("h_Lp_L"));
        h_zg_gen = dynamic_cast<TH1F*>(fGen->Get("h_zg_subjetsLL"));
        h_ptheta_gen = dynamic_cast<TH1F*>(fGen->Get("h_ptheta_subjetsLL"));
        h_zg_ratio_gen = dynamic_cast<TH1F*>(fGen->Get("h_zg_ratioLL"));

    } 
    else{ 

        for (auto &f : samples_ssWWTT) files.push_back(f);
        h_cos2D_gen = dynamic_cast<TH1F*>(fGen->Get("h_cos2D_T"));
        if (!h_cos2D_gen) {
            Error("MyCode", "Istogramma h_cos2D_T non trovato");
            return;
        }
        h_costheta_star_gen = dynamic_cast<TH1F*>(fGen->Get("h_costheta_star_T"));
        if (!h_costheta_star_gen) {
            Error("MyCode", "Istogramma h_costheta_star_T non trovato");
            return;
        }
        h_Lp_gen = dynamic_cast<TH1F*>(fGen->Get("h_Lp_T"));
        h_zg_gen = dynamic_cast<TH1F*>(fGen->Get("h_zg_subjetsTT"));
        h_ptheta_gen = dynamic_cast<TH1F*>(fGen->Get("h_ptheta_subjetsTT"));
        h_zg_ratio_gen = dynamic_cast<TH1F*>(fGen->Get("h_zg_ratioTT"));


    }

    TH1F *h_nJet = new TH1F(("h_nJet_" + tag[isample]).c_str(), "Number of clean jets;N_{jet};Events", 15, 0, 15);
    TH1F *h_nFatJet = new TH1F(("h_nFatJet_" + tag[isample]).c_str(), "Number of clean fatjets;N_{fatjet};Events", 10, 0, 10);
    TH1F *h_max_mjj = new TH1F(("h_max_mjj_" + tag[isample]).c_str(), "Max m_{jj} of VBS jets; m_{jj} [GeV]; Events", 50, 0, 2000);
    TH1F *h_deltaEta_vbs = new TH1F(("h_deltaEta_vbs_" + tag[isample]).c_str(), "#Delta#eta of VBS jets;#Delta#eta;Events", 50, 0, 10);
    TH1F *h_MT_lepton = new TH1F(("h_MT_lepton_" + tag[isample]).c_str(), "Lepton-MET MT; MT [GeV]; Events", 50, 0, 500);
    TH2F *h_nJet_vs_nFatJet = new TH2F(("h_nJet_vs_nFatJet_" + tag[isample]).c_str(), "NJet vs NFatJet;N_{jet};N_{fatjet}", 15, 0, 15, 10, 0, 10);


    TH1I *h_singleLepton = new TH1I(("h_singleLepton_" + tag[isample]).c_str(), "h_singleLepton", 2, 0, 2);
    TH1I *h_singleCleanElectron_isowp90 = new TH1I(("h_singleElectron_iswp90_" + tag[isample]).c_str(), "h_singleElectron_iswp90", 2, 0, 2);
    TH1I *h_singleCleanElectron_isowp80 = new TH1I(("h_singleElectron_iswp80_" + tag[isample]).c_str(), "h_singleElectron_iswp80", 2, 0, 2);
    TH1I *h_singleCleanMuon = new TH1I(("h_singleCleanMuon_" + tag[isample]).c_str(), "h_singleCleanMuon", 2, 0, 2);
    TH1I *h_singleCleanLepton = new TH1I(("h_singleCleanLepton_" + tag[isample]).c_str(), "h_singleCleanLepton", 2, 0, 2);
    TH1I *h_CleanFatJet = new TH1I(("h_CleanFatJet_" + tag[isample]).c_str(), "h_CleanFatJet", 2, 0, 2);
    TH1I *h_vbs_topology = new TH1I(("h_vbs_topology_" + tag[isample]).c_str(), "h_vbs_topology", 2, 0, 2);
    TH1I *h_resolvedWjet = new TH1I(("h_resolvedWjet_" + tag[isample]).c_str(), "h_resolvedWjet", 2, 0, 2);
    TH1I *h_boostedWjet = new TH1I(("h_boostedWjet_" + tag[isample]).c_str(), "h_boostedWjet", 2, 0, 2);
    TH1F *h_vbs_mass = new TH1F(("h_vbs_mass_" + tag[isample]).c_str(), "h_vbs_mass", 50, 200, 5000);
    TH1F *h_subjet_zg_ak4 = new TH1F(("h_subjet_zg_ak4_" + tag[isample]).c_str(), "h_subjet_zg_ak4", 30, 0, 0.5);
    TH1F *h_subjet_ratio_ak4 = new TH1F(("h_subjet_ratio_ak4_" + tag[isample]).c_str(), "h_subjet_ratio_ak4", 30, 0, 0.5);
    TH1F *h_deltaR_ak4 = new TH1F(("h_deltaR_ak4_" + tag[isample]).c_str(), "h_deltaR_ak4",50, 0. , 1.5);
    TH1F *h_deltaR_ak8 = new TH1F(("h_deltaR_ak8_" + tag[isample]).c_str(), "h_deltaR_ak8", 50, 0. , 1.5);
    TH1F *h_mass_ak4 = new TH1F(("h_mass_ak4" + tag[isample]).c_str(), "h_mass_ak4", 30, 60,110);
    TH1F *h_subjet_zg_ak8_custom1 = new TH1F(("h_subjet_zg_ak8_custom1_" + tag[isample]).c_str(), "h_subjet_custom_zg_1", 30, 0, 0.5);
    TH1F *h_subjet_zg_ak8_custom2 = new TH1F(("h_subjet_zg_ak8_custom2_" + tag[isample]).c_str(), "h_subjet_custom_zg_2", 30, 0, 0.5);
    TH2F *h_zg_tau21 = new TH2F("h_zg_tau21", "h_zg_tau21", 30, 0.1, 0.5,30, 0, 1);
    TH1F *h_CleanFatJet_pt = new TH1F("h_CleanFatJet_pt", "h_CleanFatJet_pt", 50, 100, 500);

    TH1F *h_rho_central = new TH1F(("rho_central" + tag[isample]).c_str(), "rho_central", 60, 0, 80); 
    TH1F *h_rho_all = new TH1F(("rho_all" + tag[isample]).c_str(), "rho_all", 60, 0, 80); 
    TH1F *h_eta_subjets = new TH1F("h_subjets_eta" , "h_subjets_eta", 40, -5, 5);
    TH1F *h_rho_diff = new TH1F(("rho_diff" + tag[isample]).c_str(), "rho_diff", 60, 0, 80); 



    TH1F *h_subjet_zg_ak8 = new TH1F(("h_subjet_zg_ak8_" + tag[isample]).c_str(), "h_subjet_zg_ak8", 30, 0.1, 0.5);
    TH1F *h_subjet_ratio_ak8 = new TH1F(("h_subjet_ratio_ak8_" + tag[isample]).c_str(), "h_subjet_ratio_ak8", 30, 0, 0.5);
    TH1F *h_ptheta_ak8 = new TH1F(("h_ptheta_ak8_" + tag[isample]).c_str(), "h_ptheta_ak8", 30, 0, 0.9);
    TH1F *h_mass_ak8 = new TH1F(("h_mass_ak8" + tag[isample]).c_str(), "h_mass_ak8", 30, 60, 110);
    TH1F *h_pt_ak8 = new TH1F(("h_pt_ak8" + tag[isample]).c_str(), "h_pt_ak8", 30, 170, 1000);
    TH1F *h_Lp_vbs = new TH1F(("h_Lp_vbs" + tag[isample]).c_str(), "h_Lp_vbs", 30, 0, 1.2);
    TH1F *h_cos2D_vbs = new TH1F(("h_cos2D_vbs" + tag[isample]).c_str(), "h_cos2D_vbs", 30, -1, 1);
    TH1F *h_costheta_star_vbs = new TH1F(("h_costheta_star_vbs" + tag[isample]).c_str(), "h_costheta_star_vbs", 30, 0-1, 1);





    // --- Cutflow histograms ---
TH1F *h_cutflow = new TH1F(("h_cutflow" + tag[isample]).c_str(), "Cutflow; ;Events", 13, 0.5, 10.5);
    double Wmass = 80;
    int total_cleanfat_jets = 0;
    int total_good_electron = 0;
    int total_good_muon = 0;
    int total_clean_jet=0;
int nentries=0;
  //  t->Add("/eos/user/l/ldellape/VBS/VBS_cards/privateMCproduction/jobs_output/MadSpin_WWLL/job_1/MadSpin_WWLL_Run3Summer23wmLHEGS/root/Nanov12/NanoAODv12_Run3Summer23wmLHEGS.root");
  //  t->Add("/eos/user/l/ldellape/VBS/VBS_cards/privateMCproduction/nanov12_samples/MadSpin/WW_TT/ssWW_TT_madspin_nanoAODv12_Run3Summer23wmLHEGS_merged.root");
    for(auto &fileName : files){
    TFile *f = TFile::Open(fileName.c_str());
    //std::cout<<"processing: "<<fileName<<std::endl;
    TTree *t = (TTree*) f->Get("Events");

    size_t pos = fileName.rfind(".root");
    fileName = fileName.substr(0, pos);
    std::cout<<fileName<<std::endl;
    // Maximum sizes


    // Set branch addresses
    t->SetBranchAddress("event", &event);
    t->SetBranchAddress("Jet_pt", Jet_pt); t->SetBranchAddress("nJet", &nJet);
    t->SetBranchAddress("Jet_eta", Jet_eta);
    t->SetBranchAddress("Jet_phi", Jet_phi);
    t->SetBranchAddress("Jet_mass", Jet_mass);
    t->SetBranchAddress("Rho_fixedGridRhoFastjetCentral", &Rho_fixedGridRhoFastjetCentral);
    t->SetBranchAddress("Rho_fixedGridRhoFastjetAll", &Rho_fixedGridRhoFastjetAll);

    t->SetBranchAddress("FatJet_pt", FatJet_pt); t->SetBranchAddress("nFatJet", &nFatJet);
    t->SetBranchAddress("FatJet_eta", FatJet_eta);
    t->SetBranchAddress("FatJet_phi", FatJet_phi);
    t->SetBranchAddress("FatJet_mass", FatJet_mass);
    t->SetBranchAddress("FatJet_msoftdrop", FatJet_msoftdrop);
    t->SetBranchAddress("FatJet_tau1", FatJet_tau1);
    t->SetBranchAddress("FatJet_tau2", FatJet_tau2);
    t->SetBranchAddress("FatJet_subJetIdx1", FatJet_subJetIdx1);
    t->SetBranchAddress("FatJet_subJetIdx2", FatJet_subJetIdx2);
  //  t->SetBranchAddress("FatJet_electronIdx3SJ", FatJet_electronIdx3SJ);
  //  t->SetBranchAddress("FatJet_muonIdx3SJ",FatJet_muonIdx3SJ);
    t->SetBranchAddress("SubJet_pt", SubJet_pt);
    t->SetBranchAddress("SubJet_eta", SubJet_eta);
    t->SetBranchAddress("SubJet_phi", SubJet_phi);
    t->SetBranchAddress("SubJet_mass", SubJet_mass);

    t->SetBranchAddress("Electron_pt", Electron_pt); t->SetBranchAddress("nElectron", &nElectron);
    t->SetBranchAddress("Electron_pfRelIso03_all", &Electron_pfRelIso03_all);
    t->SetBranchAddress("Electron_eta", Electron_eta);
    t->SetBranchAddress("Electron_phi", Electron_phi);
    t->SetBranchAddress("Electron_mvaIso_WP90", Electron_mvaIso_WP90);
    t->SetBranchAddress("Electron_mvaIso_WP80", Electron_mvaIso_WP80);
    t->SetBranchAddress("Electron_deltaEtaSC", Electron_deltaEtaSC);
    t->SetBranchAddress("Electron_mass", Electron_mass);

    t->SetBranchAddress("Muon_pt", Muon_pt); t->SetBranchAddress("nMuon", &nMuon);
    t->SetBranchAddress("Muon_eta", Muon_eta);
    t->SetBranchAddress("Muon_phi", Muon_phi);
    t->SetBranchAddress("Muon_tightId", Muon_tightId);
    t->SetBranchAddress("Muon_pfRelIso04_all", Muon_pfRelIso04_all);
    t->SetBranchAddress("Muon_mass", Muon_mass);


    t->SetBranchAddress("MET_pt", &MET_pt);
    t->SetBranchAddress("MET_phi", &MET_phi);
    t->SetBranchAddress("Jet_btagRobustParTAK4B", Jet_btagRobustParTAK4B);

    t->SetBranchAddress("nGenPart", &nGenParticle);
    t->SetBranchAddress("GenPart_pt", GenPart_pt);
    t->SetBranchAddress("GenPart_eta", GenPart_eta);
    t->SetBranchAddress("GenPart_phi", GenPart_phi);
    t->SetBranchAddress("GenPart_mass", GenPart_mass);
    t->SetBranchAddress("GenPart_status", GenPart_status);
    t->SetBranchAddress("GenPart_statusFlags", GenPart_statusFlags);
    t->SetBranchAddress("GenPart_genPartIdxMother", GenPart_genPartIdxMother);
    t->SetBranchAddress("GenPart_pdgId", GenPart_pdgId);
    t->SetBranchAddress("nSubGenJetAK8", &nSubGenJetAK8);
    t->SetBranchAddress("SubGenJetAK8_pt", SubGenJetAK8_pt);
    t->SetBranchAddress("SubGenJetAK8_eta", SubGenJetAK8_eta);
    t->SetBranchAddress("SubGenJetAK8_phi", SubGenJetAK8_phi);
    t->SetBranchAddress("SubGenJetAK8_mass", SubGenJetAK8_mass);
    t->SetBranchAddress("nGenJetAK8", &nGenJetAK8);
    t->SetBranchAddress("GenJetAK8_pt", GenJetAK8_pt);
    t->SetBranchAddress("GenJetAK8_eta", GenJetAK8_eta);
    t->SetBranchAddress("GenJetAK8_phi", GenJetAK8_phi);
    t->SetBranchAddress("GenJetAK8_mass", GenJetAK8_mass);
    t->SetBranchAddress("GenJetAK8_partonFlavour", GenJetAK8_partonFlavour);




  //  #include "func.cpp"


    // === Histograms ===
// === Histograms ===

    h_cutflow->GetXaxis()->SetBinLabel(1, "Single lepton");
    h_cutflow->GetXaxis()->SetBinLabel(2, "Good Electron");
    h_cutflow->GetXaxis()->SetBinLabel(3, "Good Muon");
    //h_cutflow->GetXaxis()->SetBinLabel(4, "Good Jet");
    h_cutflow->GetXaxis()->SetBinLabel(4, "Good FatJet");
    //h_cutflow->GetXaxis()->SetBinLabel(6, "#geq 4 AK4");
    h_cutflow->GetXaxis()->SetBinLabel(5, "1AK8 + #geq 2AK4");
    h_cutflow->GetXaxis()->SetBinLabel(6, "VBS kinematics");
    h_cutflow->GetXaxis()->SetBinLabel(7, "Lep.Decay");
    //h_cutflow->GetXaxis()->SetBinLabel(10, "Had.Decay (2AK4)");
    h_cutflow->GetXaxis()->SetBinLabel(8, "Had.Decay (1AK8)");
    //h_cutflow->GetXaxis()->SetBinLabel(9, "pre-selected no b-tag ");
    h_cutflow->GetXaxis()->SetBinLabel(9, "pre-selected");

    Long64_t events_tree = t->GetEntries();
    nentries+= events_tree;


    
    for(Long64_t i=0; i<events_tree; ++i){
        t->GetEntry(i);

    
        bool singleLepton = false;
        if((nElectron >= 1  ||  nMuon>=1)){
            singleLepton = true;
            h_singleLepton->Fill(singleLepton);
        }
        else{
            continue;
        }
        h_cutflow->AddBinContent(1); // Passed single lepton

        // good muons and electrons
        bool found_ele_wp80 = false, found_ele_wp90 = false, found_muon = false;

    
    



        // exactly 1 lepton
        int count_electron = 0;
        int count_muon = 0;
        int idx_ele = 0;
        int idx_muon = 0;
        int low_pt_leptons = 0;

        for(int ee=0; ee<nElectron; ee++){
         float etaSC = std::abs(Electron_eta[ee] + Electron_deltaEtaSC[ee]);
         if(etaSC >= 1.4442 && etaSC <= 1.5660) continue;
         if(Electron_mvaIso_WP80[ee] && abs(Electron_eta[ee]) < good_ele["eta"] && Electron_pt[ee] > good_ele["pt"] && Electron_pfRelIso03_all[ee] < good_ele["iso"]){  h_cutflow->AddBinContent(2); count_electron++; idx_ele=ee; found_ele_wp80 = true;}
         if(Electron_mvaIso_WP90[ee] && abs(Electron_eta[ee]) < good_ele["eta"] && Electron_pt[ee] > good_ele["pt"] && Electron_pfRelIso03_all[ee] < good_ele["iso"]) found_ele_wp90 = true; 
            h_singleCleanElectron_isowp80->Fill(found_ele_wp80);
            h_singleCleanElectron_isowp90->Fill(found_ele_wp90);

            if(Electron_pt[ee] < 10) low_pt_leptons++;
        } 
        for(int mu=0; mu<nMuon; mu++){
            if(Muon_pt[mu] > good_muon["pt"] && Muon_tightId[mu] && Muon_pfRelIso04_all[mu] < good_muon["iso"] && abs(Muon_eta[mu]) < good_muon["eta"]) {h_cutflow->AddBinContent(3); count_muon++; idx_muon=mu; found_muon = true; }
            h_singleCleanMuon->Fill(found_muon);

            if(Muon_pt[mu] < 10) low_pt_leptons++;
        }

        bool singleGoodleptons = false;

        if(((count_electron==1 && count_muon==0) || (count_electron==0 && count_muon==1)) ) singleGoodleptons=true;
        if(count_muon==1 && count_electron == 0) ++total_good_muon;
        if(count_electron==1 && count_muon == 0) ++total_good_electron; 
        if(!singleGoodleptons) continue;



        int cleanfatjet = 0;
        std::vector<int> cleanFatJet_idx;
        for(Int_t fj=0; fj<nFatJet; fj++){
            if(FatJet_pt[fj] < good_fatjet["pt"]  || FatJet_msoftdrop[fj] < good_fatjet["msd"] || abs(FatJet_eta[fj]) > good_fatjet["eta"]) continue; 
            if(count_electron==1 && deltaR(Electron_eta[idx_ele], Electron_phi[idx_ele], FatJet_eta[fj], FatJet_phi[fj]) < 0.8) continue;
            if(count_muon==1 && deltaR(Muon_eta[idx_muon], Muon_phi[idx_muon], FatJet_eta[fj], FatJet_phi[fj]) < 0.8) continue;
            ++cleanfatjet;
            ++total_cleanfat_jets;
            cleanFatJet_idx.emplace_back(fj);
        }

        // --- Clean jets
        int cleanjet = 0;
        std::vector<int> cleanjet_idx;
        for(Int_t j=0; j<nJet; j++){
            if(Jet_pt[j] < good_jet["pt"] || Jet_eta[j] > good_jet["eta"]) continue;
            if(count_electron==1 && deltaR(Electron_eta[idx_ele], Electron_phi[idx_ele], Jet_eta[j], Jet_phi[j]) < 0.4) continue;
            if(count_muon==1 && deltaR(Muon_eta[idx_muon], Muon_phi[idx_muon], Jet_eta[j], Jet_phi[j]) < 0.4) continue;           
            ++cleanjet;
            ++total_clean_jet;
            cleanjet_idx.push_back(j);
        }

        h_nJet->Fill(cleanjet);
        h_nFatJet->Fill(cleanfatjet);
        h_nJet_vs_nFatJet->Fill(cleanjet, cleanfatjet);
        //if(cleanjet>0) h_cutflow->AddBinContent(4);
        if(cleanfatjet>0) h_cutflow->AddBinContent(4);

        if(cleanfatjet == 1 && cleanjet >= 2){
            h_cutflow->AddBinContent(5);
        }
        if(cleanfatjet == 0 && cleanjet >= 4){
         //h_cutflow->AddBinContent(6);
        }


        std::pair<double,double> vbs_jets;
        std::pair<int,int> vbs_jets_idx;
        double max_mjj = -1.0;
        int idx_leading_pt = 0;
        for(Int_t j=0;j<cleanjet_idx.size();j++){
            for(Int_t jj=j+1;jj<cleanjet_idx.size();jj++){
                double m_jj = invariantMass(Jet_pt[cleanjet_idx[j]],Jet_eta[cleanjet_idx[j]],Jet_phi[cleanjet_idx[j]],
                                            Jet_pt[cleanjet_idx[jj]],Jet_eta[cleanjet_idx[jj]], Jet_phi[cleanjet_idx[jj]]);
                if(m_jj>max_mjj){
                    max_mjj=m_jj;
                    vbs_jets={deltaEta(Jet_eta[cleanjet_idx[j]],Jet_eta[cleanjet_idx[jj]]), m_jj};
                    vbs_jets_idx={cleanjet_idx[j] , cleanjet_idx[jj]};
                    if(Jet_pt[cleanjet_idx[j]] > Jet_pt[cleanjet_idx[jj]]) idx_leading_pt = cleanjet_idx[j];
                    else{ idx_leading_pt = cleanjet_idx[jj];}
                }
            }
        }
        h_max_mjj->Fill(max_mjj);
        h_deltaEta_vbs->Fill(std::abs(vbs_jets.first));

        bool vbs_topology = ((cleanjet>=2 && cleanfatjet==1) || (cleanjet>=4 && cleanfatjet==0));

        if(!vbs_topology) continue;

        bool vbs_kinematics = (std::abs(vbs_jets.first)> VBS_presel["eta"] && vbs_jets.second > VBS_presel["m_jj"] && Jet_pt[idx_leading_pt] > VBS_presel["pt_leading"]);

        h_vbs_mass->Fill(vbs_jets.second);

        if(!vbs_kinematics) continue;
            h_cutflow->AddBinContent(6); // Passed VBS kinematics


        // --- Leptonic MT
        double MT=0;
        double METx = MET_pt * cos(MET_phi);
        double METy = MET_pt * sin(MET_phi);
        if(found_ele_wp80 || found_ele_wp90){
            MT = sqrt(2*Electron_pt[0]*MET_pt*(1-cos(Electron_phi[0]-MET_phi)));
            if(Electron_pt[idx_ele] > V_to_lep["ele_pt"] && MET_pt > V_to_lep["MET_pt"]) {
                h_cutflow->AddBinContent(7);
            }
        } else if(found_muon){
            MT = sqrt(2*Muon_pt[0]*MET_pt*(1-cos(Muon_phi[0]-MET_phi)));
            if(Muon_pt[idx_muon] > V_to_lep["muon_pt"] && MET_pt > V_to_lep["MET_pt"]) h_cutflow->AddBinContent(7);
        }
        if(MT > 0) h_MT_lepton->Fill(MT);



    
    
        // --- W reconstruction
 
        bool onshell_W = false;
        /*if(cleanfatjet == 0 && cleanjet>=4){
            std::pair<int,int> wjets_idx;
            double closest_mw_diff = 1e6;
            double w_mjj = -1.0;

            for (int j = 0; j < cleanjet_idx.size(); j++) {
                if (j == vbs_jets_idx.first || j == vbs_jets_idx.second) continue;
                for (int jj = j + 1; jj < cleanjet_idx.size(); jj++) {
                    if (jj == vbs_jets_idx.first || jj == vbs_jets_idx.second) continue;
                    double m_jj = invariantMass(Jet_pt[ cleanjet_idx[j]], Jet_eta[ cleanjet_idx[j]], Jet_phi[cleanjet_idx[j]],
                                                Jet_pt[cleanjet_idx[jj]], Jet_eta[cleanjet_idx[jj]], Jet_phi[cleanjet_idx[jj]]);
                    double diff = std::abs(m_jj - Wmass);
                    if (diff < closest_mw_diff) {
                        closest_mw_diff = diff;
                        w_mjj = m_jj;
                        wjets_idx = {j, jj};
                    }
                }
            }
            if(w_mjj > V_to_had["mass_min"] && w_mjj < V_to_had["mass_max"]){ 
                h_cutflow->AddBinContent(10);
               float zg = std::min( Jet_pt[cleanjet_idx[wjets_idx.first]] , Jet_pt[cleanjet_idx[wjets_idx.second]] )/ (Jet_pt[cleanjet_idx[wjets_idx.first]] + Jet_pt[cleanjet_idx[wjets_idx.second]]);
                h_subjet_zg_ak4->Fill(zg);
                h_subjet_ratio_ak4->Fill(  std::min( Jet_pt[cleanjet_idx[wjets_idx.first]] , Jet_pt[cleanjet_idx[wjets_idx.second]] ) / std::max( Jet_pt[cleanjet_idx[wjets_idx.first]] , Jet_pt[cleanjet_idx[wjets_idx.second]] ));
                h_rho_all->Fill(Rho_fixedGridRhoFastjetAll);
                h_rho_central->Fill(Rho_fixedGridRhoFastjetCentral);
                float deltaR_ak4 = deltaR(
                    Jet_eta[cleanjet_idx[wjets_idx.first]] , Jet_phi[cleanjet_idx[wjets_idx.first]], 
                    Jet_eta[cleanjet_idx[wjets_idx.second]] , Jet_phi[cleanjet_idx[wjets_idx.second]] 
                );
                h_deltaR_ak4->Fill(deltaR_ak4);
                h_mass_ak4->Fill(w_mjj);
                onshell_W = true; }
        } */
        if(cleanfatjet==1 && cleanjet>=2) {
            if(FatJet_msoftdrop[cleanFatJet_idx[0]] > V_to_had["m_sd_min"] && FatJet_msoftdrop[cleanFatJet_idx[0]] < V_to_had["m_sd_max"] ){
                h_cutflow->AddBinContent(8);
                h_rho_all->Fill(Rho_fixedGridRhoFastjetAll);
                h_rho_central->Fill(Rho_fixedGridRhoFastjetCentral);
                float zg = std::min( SubJet_pt[FatJet_subJetIdx1[cleanFatJet_idx[0]]] , SubJet_pt[FatJet_subJetIdx2[cleanFatJet_idx[0]]] ) / (SubJet_pt[FatJet_subJetIdx1[cleanFatJet_idx[0]]] +  SubJet_pt[FatJet_subJetIdx2[cleanFatJet_idx[0]]]);
                h_subjet_zg_ak8->Fill(zg);
                h_zg_tau21->Fill(zg,  (FatJet_tau2[cleanFatJet_idx[0]]/FatJet_tau1[cleanFatJet_idx[0]]));
                h_CleanFatJet_pt->Fill(FatJet_pt[cleanFatJet_idx[0]]);
                float deltaR_ak8 = deltaR(
                    SubJet_eta[FatJet_subJetIdx1[cleanFatJet_idx[0]]] , SubJet_phi[FatJet_subJetIdx1[cleanFatJet_idx[0]]],
                    SubJet_eta[FatJet_subJetIdx2[cleanFatJet_idx[0]]] , SubJet_phi[FatJet_subJetIdx2[cleanFatJet_idx[0]]]              
                );
                h_deltaR_ak8->Fill(deltaR_ak8);
                h_eta_subjets->Fill(SubJet_eta[FatJet_subJetIdx1[cleanFatJet_idx[0]]] );
                h_eta_subjets->Fill(SubJet_eta[FatJet_subJetIdx2[cleanFatJet_idx[0]]] );
                h_mass_ak8->Fill( FatJet_msoftdrop[cleanFatJet_idx[0]]);
                h_pt_ak8->Fill( FatJet_pt[cleanFatJet_idx[0]]);
                onshell_W = true;

                if(FatJet_pt[0] >= pt_lim){
                    funzione(t, i);
                }


                h_subjet_ratio_ak8->Fill(std::min( SubJet_pt[FatJet_subJetIdx1[cleanFatJet_idx[0]]] , SubJet_pt[FatJet_subJetIdx2[cleanFatJet_idx[0]]] )/ std::max( SubJet_pt[FatJet_subJetIdx1[cleanFatJet_idx[0]]] , SubJet_pt[FatJet_subJetIdx2[cleanFatJet_idx[0]]] ));

                double E1, E2;
                E1 = std::sqrt(SubJet_pt[FatJet_subJetIdx1[cleanFatJet_idx[0]]]*SubJet_pt[FatJet_subJetIdx1[cleanFatJet_idx[0]]] + SubJet_mass[FatJet_subJetIdx1[cleanFatJet_idx[0]]]*SubJet_mass[FatJet_subJetIdx1[cleanFatJet_idx[0]]]);
                E2 = std::sqrt(SubJet_pt[FatJet_subJetIdx2[cleanFatJet_idx[0]]]*SubJet_pt[FatJet_subJetIdx2[cleanFatJet_idx[0]]] + SubJet_mass[FatJet_subJetIdx2[cleanFatJet_idx[0]]]*SubJet_mass[FatJet_subJetIdx2[cleanFatJet_idx[0]]]);
                h_ptheta_ak8->Fill(std::abs(E1-E2)/FatJet_pt[cleanFatJet_idx[0]]);


                

                int maxSubjets = 40;
            }
        }
        if(!onshell_W) continue;

        bool b_tagged = false;
        for(int jj= 0; jj<cleanjet_idx.size(); jj++){
            if(Jet_btagRobustParTAK4B[cleanjet_idx[jj]] > good_jet["wp_btag"] && abs(Jet_eta[cleanjet_idx[jj]]) < 2.5){
                b_tagged=true;
                break;
            }
        }
        //h_cutflow->AddBinContent(9);
        if(!b_tagged) h_cutflow->AddBinContent(9); 



  
    }
}

    // --- Efficiency histogram ---
    TH1F *h_efficiency = (TH1F*)h_cutflow->Clone("h_efficiency");
    h_efficiency->SetTitle("Cut efficiencies");
    h_efficiency->SetYTitle("Efficiency");

    // --- Save output ---
    TFile* outputFile = new TFile(("output_" + tag_cutflow  + tag[isample] + ".root").c_str(),"RECREATE");
    h_nJet->Write();
    h_nFatJet->Write();
    h_max_mjj->Write();
    h_deltaEta_vbs->Write();
    h_MT_lepton->Write();
    h_nJet_vs_nFatJet->Write();
    h_singleLepton->Write();
    h_singleCleanElectron_isowp90->Write();
    h_singleCleanElectron_isowp80->Write();
    h_singleCleanMuon->Write();
    h_singleCleanLepton->Write();
    h_CleanFatJet->Write();
    h_vbs_topology->Write();
    h_resolvedWjet->Write();
    h_subjet_ratio_ak8->Write();
    h_subjet_ratio_ak4->Write();

    h_subjet_zg_ak4->Write();
    h_eta_subjets->Write();
    h_subjet_zg_ak8->Write();
    h_boostedWjet->Write();
    h_cutflow->Write();
    h_deltaR_ak4->Write();
    h_zg_tau21->Write();
    h_deltaR_ak8->Write();
    h_vbs_mass->Write();
    h_rho_central->Write();
    h_mass_ak4->Write();
    h_mass_ak8->Write();
    h_subjet_zg_ak8_custom1->Write();
    h_subjet_zg_ak8_custom2->Write();
    h_CleanFatJet_pt->Write();
    h_rho_all->Write();
    h_rho_diff = (TH1F*) h_rho_all->Clone("h_rho_diff");
    h_rho_diff->Add(h_rho_central, -1);
    h_rho_diff->Write();
    h_efficiency->Write();

    h_Lp->Write();
    h_costheta_star->Write();
    h_cos2D->Write();
    h_zg_subjets->Write();
    h_ptheta_subjets->Write();
    h_costheta_star_gen->Write();

    outputFile->Close();



    TCanvas *c_cutflow = new TCanvas(("c_cutflow_" + tag[isample]).c_str(), "Cutflow", 1200, 600);
    h_cutflow->SetStats(0);
    h_cutflow->SetFillColor(kAzure-9);
    h_cutflow->SetBarWidth(0.9);
    h_cutflow->SetBarOffset(0.05);
    double ymax = h_cutflow->GetMaximum();
    h_cutflow->SetMaximum(ymax*1.2);
    h_cutflow->Draw("hist bar");
    int nCuts = 9;  // numero reale di bin usati
    h_cutflow->GetXaxis()->SetRange(1, nCuts);

std::map<int, std::string> eff = {
    {1,  Form("#epsilon = %.2f", h_cutflow->GetBinContent(1) / nentries)}, //single lepton
    {2,  Form("#epsilon = %.2f", h_cutflow->GetBinContent(2) / h_cutflow->GetBinContent(1))}, //electrons
    {3,  Form("#epsilon = %.2f", h_cutflow->GetBinContent(3) / h_cutflow->GetBinContent(1))}, //muons
    //{4,  Form("#epsilon = %.2f", h_cutflow->GetBinContent(4) / (h_cutflow->GetBinContent(2) + h_cutflow->GetBinContent(3)))},
    {4,  Form("#epsilon = %.2f", h_cutflow->GetBinContent(4) / (h_cutflow->GetBinContent(2) + h_cutflow->GetBinContent(3)))}, //good fatjet
    //{6,  Form("#epsilon = %.2f", h_cutflow->GetBinContent(6) / h_cutflow->GetBinContent(1))},
    {5,  Form("#epsilon = %.2f", h_cutflow->GetBinContent(5) / h_cutflow->GetBinContent(1))}, //1AK8 + 2AK4
    {6,  Form("#epsilon = %.2f", h_cutflow->GetBinContent(6) / h_cutflow->GetBinContent(1))}, //VBS kinematics
    {7,  Form("#epsilon = %.2f", h_cutflow->GetBinContent(7) / h_cutflow->GetBinContent(6))}, //Lep decay
    //{8, Form("#epsilon = %.2f", h_cutflow->GetBinContent(8) / h_cutflow->GetBinContent(6))},
    {8, Form("#epsilon = %.2f", h_cutflow->GetBinContent(8) / h_cutflow->GetBinContent(6))}, //Had decay
    //{12, Form("")},
    {9, Form("#epsilon = %.2f", h_cutflow->GetBinContent(9) / nentries)} //b jets
};


TLatex latex;
latex.SetNDC();
latex.SetTextSize(0.03);
latex.SetTextAlign(13);  // allineamento a sinistra

int count_text = 0;

for (int bin = 1; bin <= 10; ++bin) {
    if (!cutDescriptions[bin].empty()) {
        ++count_text;
        latex.DrawLatex(0.4, 0.9 - 0.05*count_text,
                        cutDescriptions[bin].c_str());
    }
}

float lumi = nentries/(xsec[isample]*1e3);
if( isample==0){
h_cutflow->SetTitle(Form("cutflow, VBS ssWW LL, entries=%i    CMS simulation (private work)", nentries)); 
}
else if(isample==1){
h_cutflow->SetTitle(Form("cutflow, VBS ssWW TT, entries=%i  CMS simulation (private work)", nentries)); 
}
else if(isample==2){
h_cutflow->SetTitle(Form("cutflow, VBS ssWW TL, entries=%i  CMS simulation (private work)", nentries)); 
}

// Loop 1: draw bin numbers and efficiency
for (int i = 1; i <= h_cutflow->GetNbinsX(); ++i) {
    double x = h_cutflow->GetBinCenter(i);
    double y = h_cutflow->GetBinContent(i);
    double y2 = h_cutflow->GetBinContent(8);


    if (y > 0) {
        TLatex text;
        text.SetTextAlign(22);  // center
        text.SetTextSize(0.04);
        text.DrawLatex(x, y + 0.1*y2, Form("%.0f", y));  // bin content

        TLatex text2;
        text2.SetTextAlign(22);
        text2.SetTextSize(0.04);
        text2.DrawLatex(x, y + 0.1*y2 + 5000, eff[i].c_str()); // vector indexed at i-1
    }
}


// Finalize canvas
c_cutflow->SetTicks(1, 1);
c_cutflow->Update();
c_cutflow->SaveAs(("cuts_" + tag_cutflow +"_" + tag[isample] + ".pdf").c_str());

std::cout<<"total fj: "<<total_cleanfat_jets<<std::endl;
std::cout<<"total jet: "<<total_clean_jet<<std::endl;
std::cout<<"good ele: "<<total_good_electron<<std::endl;
std::cout<<"good muon: "<<total_good_muon<<std::endl;


gStyle->SetOptStat(0);

std::vector<Color_t> colors;
colors.push_back(kRed-3);
colors.push_back(kSpring-6);



TCanvas* c_zg_ratio_vbs = new TCanvas("c_zg_ratio_vbs", "c_zg_ratio_vbs", 1000, 800);
std::vector<TH1*> zg_ratio_vbs;
h_zg_ratio->SetTitle(Form("z_{g}'     p_{T} > %d  %s (CMS private work)", pt_lim, tag2[isample].c_str()));
zg_ratio_vbs.push_back(h_zg_ratio);
zg_ratio_vbs.push_back(h_zg_ratio_gen);
Makecanva(c_zg_ratio_vbs, zg_ratio_vbs, kRed-3);
TLegend *leg7 = new TLegend(0.55, 0.7, 0.88, 0.85);
leg7->AddEntry(h_zg_ratio_gen, "Gen level", "f");
leg7->AddEntry(h_zg_ratio, "Gen level (VBS selection)", "f");
leg7->SetTextSize(0.03);
leg7->SetBorderSize(0);
leg7->SetFillStyle(0);
leg7->SetMargin(0.25);
leg7->Draw();
c_zg_ratio_vbs->SaveAs(("plots_z1/c_zg_ratio_vbs" + tag2[isample] +".pdf").c_str());


TCanvas* c_zg_vbs = new TCanvas("c_zg_vbs", "c_zg_vbs", 1000, 800);
std::vector<TH1*> zg_vbs;
h_zg_subjets->SetTitle(Form("z_{g}     p_{T} > %d  %s (CMS private work)", pt_lim, tag2[isample].c_str()));
zg_vbs.push_back(h_zg_subjets);
zg_vbs.push_back(h_zg_gen);
Makecanva(c_zg_vbs, zg_vbs, kRed-3);
TLegend *leg2 = new TLegend(0.47, 0.74, 0.75, 0.89);
leg2->AddEntry(h_zg_subjets, "Gen level (VBS selection)", "f");
leg2->AddEntry(h_zg_gen, "Gen level", "f");
leg2->SetTextSize(0.025);
leg2->SetBorderSize(0);
leg2->SetFillStyle(0);
leg2->SetMargin(0.25);
leg2->Draw();
c_zg_vbs->SaveAs(("plots_z1/c_zg_vbs" + tag2[isample] +".pdf").c_str());



TCanvas *c_ptheta_vbs = new TCanvas("c_ptheta_vbs_", "c_ptheta_vbs", 1000, 800);
std::vector<TH1*> ptheta_vbs;
h_ptheta_subjets->SetTitle(Form("p_{#theta}     p_{T} > %d  %s (CMS private work)", pt_lim, tag2[isample].c_str()));
ptheta_vbs.push_back(h_ptheta_subjets);
ptheta_vbs.push_back(h_ptheta_gen);
Makecanva(c_ptheta_vbs, ptheta_vbs, kRed-3);
TLegend *leg3 = new TLegend(0.17, 0.18, 0.45, 0.3);
leg3->AddEntry(h_ptheta_subjets, "Gen level (VBS selection)", "f");
leg3->AddEntry(h_ptheta_gen, "Gen level", "f");
leg3->SetTextSize(0.025);
leg3->SetBorderSize(0);
leg3->SetMargin(0.25);
leg3->Draw();
c_ptheta_vbs->SaveAs(("plots_ptheta/c_ptheta_vbs" + tag2[isample] +".pdf").c_str());

TCanvas *c_Lp_vbs = new TCanvas("c_Lp_vbs", "c_Lp_vbs", 1000, 800);
std::vector<TH1*> Lp_vbs;
h_Lp->SetTitle(("L_{p}  " + tag2[isample]+ " (CMS private work)").c_str());
Lp_vbs.push_back(h_Lp);
Lp_vbs.push_back(h_Lp_gen);
Makecanva(c_Lp_vbs, Lp_vbs, kGreen+1);
TLegend *leg4 = new TLegend(0.2, 0.73, 0.45, 0.88);
leg4->AddEntry(h_Lp, "Gen level (VBS selection)", "f");
leg4->AddEntry(h_Lp_gen, "Gen level", "f");
leg4->SetTextSize(0.03);
leg4->SetBorderSize(0);
leg4->SetFillStyle(0);
leg4->SetMargin(0.25);
leg4->Draw();
c_Lp_vbs->SaveAs(("plots_lp/c_Lp_vbs" + tag2[isample] +".pdf").c_str());


TCanvas *c_cos2D_vbs = new TCanvas("c_cos2D_vbs", "c_cos2D_vbs", 1000, 800);
std::vector<TH1*> cos2D_vbs;
h_cos2D->SetTitle(("cos2D  " + tag2[isample]+ " (CMS private work)").c_str());
cos2D_vbs.push_back(h_cos2D);
cos2D_vbs.push_back(h_cos2D_gen);
Makecanva(c_cos2D_vbs, cos2D_vbs, kGreen+1);
TLegend *leg5 = new TLegend(0.35, 0.6, 0.6, 0.8);
leg5->AddEntry(h_cos2D, "Gen level (VBS selection)", "f");
leg5->AddEntry(h_cos2D_gen, "Gen level", "f");
leg5->SetTextSize(0.03);
leg5->SetBorderSize(0);
leg5->SetFillStyle(0);
leg5->SetMargin(0.25);
leg5->Draw();
c_cos2D_vbs->SaveAs(("plots_lp/c_cos2D_vbs_" + tag2[isample] +".pdf").c_str());


TCanvas *c_costheta_star_vbs = new TCanvas("c_costheta_star_vbs", "c_costheta_star_vbs", 1000, 800);
std::vector<TH1*> costheta_star_vbs;
h_costheta_star->SetTitle(("cos#theta*  " + tag2[isample]+ " (CMS private work)").c_str());
costheta_star_vbs.push_back(h_costheta_star);
costheta_star_vbs.push_back(h_costheta_star_gen);
Makecanva(c_costheta_star_vbs, costheta_star_vbs, kGreen+1);
TLegend *leg6 = new TLegend(0.15, 0.68, 0.4, 0.88);
leg6->AddEntry(h_costheta_star, "Gen level (VBS selection)", "f");
leg6->AddEntry(h_costheta_star_gen, "Gen level", "f");
leg6->SetTextFont(42);
leg6->SetTextSize(0.03);
leg6->SetBorderSize(0);
leg6->SetFillStyle(0);
leg6->SetMargin(0.25);
leg6->Draw();
c_costheta_star_vbs->SaveAs(("plots_lp/c_costheta_star_vbs_" + tag2[isample] +".pdf").c_str());





}








gROOT->ProcessLine(".q");


}
