#include <iostream>
#include <vector>
#include <cmath>
#include "samples.h"
#include <utility>
#include <TChain.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TFile.h>
#include <cmath>
#include <TLegend.h>
#include "canva.cpp"
#include "TColor.h"



void GenJetAK8_match(
    TString output="Gen_match.root"
){
    TFile *f=new TFile(output, "recreate");
int pt_lim = 150;

TH1F *h_zg_quarksLL = new TH1F("h_zg_quarksLL", Form("z_{g} quarks p_{T} > %d  (L)   (CMS private work); z_{g} ; Events", pt_lim), 20,  0.1, 0.5);
TH1F *h_zg_subjetsLL = new TH1F("h_zg_subjetsLL", Form("z_{g} subjets p_{T} > %d   (L)   (CMS private work); z_{g} ; Events", pt_lim), 20, 0.1, 0.5);
TH1F *h_zg_ratioLL = new TH1F("h_zg_ratioLL", Form("z_{g}' subjets p_{T} > %d   (L)   (CMS private work); z_{g}' ; Events", pt_lim), 20, 1, 10);
TH1F *h_zg_ratio_quarksLL = new TH1F("h_zg_ratio_quarksLL", "z_{g}' quarks  (L)   (CMS private work); z_{g}' ; Events", 20,  1, 10);
TH1F *h_ptheta_quarksLL = new TH1F("h_ptheta_quarksLL", Form("p_{#theta} quarks p_{T} > %d  (L)   (CMS private work); p_{#theta} ; Events", pt_lim),20, 0, 0.9);
TH1F *h_ptheta_subjetsLL = new TH1F("h_ptheta_subjetsLL", Form("p_{#theta} subjets p_{T} > %d   (L)   (CMS private work); p_{#theta} ; Events", pt_lim), 20, 0, 0.9);

TH1F *h_zg_quarksTT = new TH1F("h_zg_quarksTT", Form("z_{g} quarks p_{T} > %d  (T)   (CMS private work); z_{g} ; Events", pt_lim), 20,  0.1, 0.5);
TH1F *h_zg_subjetsTT = new TH1F("h_zg_subjetsTT", Form("z_{g} subjets p_{T} > %d   (T)   (CMS private work); z_{g} ; Events", pt_lim), 20, 0.1, 0.5);
TH1F *h_zg_ratioTT = new TH1F("h_zg_ratioTT", Form("z_{g}' subjets p_{T} > %d   (T)   (CMS private work); z_{g}' ; Events", pt_lim), 20, 1, 10);
TH1F *h_zg_ratio_quarksTT = new TH1F("h_zg_ratio_quarksTT", "z_{g}' quarks  (T)   (CMS private work); z_{g}' ; Events", 20,  1, 10);
TH1F *h_ptheta_quarksTT = new TH1F("h_ptheta_quarksTT", Form("p_{#theta} quarks p_{T} > %d  (T)   (CMS private work); p_{#theta} ; Events", pt_lim),20, 0, 0.9);
TH1F *h_ptheta_subjetsTT = new TH1F("h_ptheta_subjetsTT", Form("p_{#theta} subjets p_{T} > %d   (T)   (CMS private work); p_{#theta} ; Events", pt_lim), 20, 0, 0.9);




//leptons
TH1F *h_Lp_LL = new TH1F("h_Lp_L", "Lepton projection   (CMS private work); L_{P}; Events", 20, 0, 1.2);
TH1F *h_Lp_TT = new TH1F("h_Lp_T", "Lepton projection   (CMS private work); L_{P}; Events", 20, 0, 1.2);
TH1F *h_cos2D_LL = new TH1F("h_cos2D_L", "Lepton cos2D   (CMS private work); cos2D; Events", 20, -1,1);
TH1F *h_costheta_star_LL = new TH1F("h_costheta_star_L", "Lepton cos #theta*   (CMS private work); cos #theta*; Events", 20, -1,1);
TH1F *h_cos2D_TT = new TH1F("h_cos2D_T", "Lepton cos2D   (CMS private work); cos2D; Events", 20, -1,1);
TH1F *h_costheta_star_TT = new TH1F("h_costheta_star_T", "Lepton cos #theta*   (CMS private work); cos #theta*; Events", 20, -1,1);

TH1F *h_lep_ptLL = new TH1F("h_lep1_ptL", "Lepton p_{T}    (CMS private work); p_{T} [GeV]; Events", 50, 0, 400);
TH1F *h_lep_ptTT = new TH1F("h_lep1_ptT", "Lepton p_{T}    (CMS private work); p_{T} [GeV]; Events", 50, 0, 400);

TH1F *h_lep_etaLL = new TH1F("h_lep_etaLL", "Lepton #eta   (CMS private work); #eta; Events", 50, -6, 6);
TH1F *h_lep_etaTT = new TH1F("h_lep_etaTT", "Lepton #eta   (CMS private work); #eta; Events", 50, -6, 6);

TH1F *h_lep_phiLL = new TH1F("h_lep_phiLL", "Lepton #phi   (CMS private work); #phi; Events", 50, -4, 4);
TH1F *h_lep_phiTT = new TH1F("h_lep_phiTT", "Lepton #phi   (CMS private work); #phi; Events", 50, 4, 4);



TH1F *h_nMatchedJetsLL = new TH1F(
"h_nMatchedJetsLL", 
"Number of matched GenJetAK8 per event (L);N matched jets", 
5, -0.5, 4.5
);
TH1F *h_nMatchedJetsTT = new TH1F(
"h_nMatchedJetsTT", 
"Number of matched GenJetAK8 per event (T);N matched jets", 
5, -0.5, 4.5
);

TH1F *h_dR_q1_jetLL = new TH1F(
    "h_dR_q1_jetLL",
    Form("#DeltaR(q_{1}, GenJetAK8) p_{T} > %d   (CMS private work);#DeltaR(q_{1}, Jet);Events", pt_lim),
    35, 0, 0.8
);
TH1F *h_dR_q1_jetTT = new TH1F(
    "h_dR_q1_jetTT",
    Form("#DeltaR(q_{1}, GenJetAK8) p_{T} > %d   (CMS private work);#DeltaR(q_{1}, Jet); Events", pt_lim),
    35, 0, 0.8
);

TH1F *h_dR_qqLL = new TH1F("h_dR_qqLL", Form("#DeltaR quarks p_{T} > %d  (L)   (CMS private work); #DeltaR; Events", pt_lim), 40, 0, 1);
TH1F *h_dR_qqTT = new TH1F("h_dR_qqTT", Form("#DeltaR quarks p_{T} > %d  (T)   (CMS private work); #DeltaR; Events", pt_lim), 40, 0, 1);

TH1F *h_dR_q2_jetLL = new TH1F(
    "h_dR_q2_jetLL",
    Form("#DeltaR(q_{2}, GenJetAK8) p_{T} > %d   (CMS private work);#DeltaR(q_{2}, Jet);Events", pt_lim),
    35, 0, 0.8
);
TH1F *h_dR_q2_jetTT = new TH1F(
    "h_dR_q2_jetTT",
    Form("#DeltaR(q_{2}, GenJetAK8) p_{T} > %d   (CMS private work);#DeltaR(q_{2}, Jet);Events", pt_lim),
    35, 0, 0.8
);

TH1F *h_matchedJet_ptLL = new TH1F(
    "h_matchedJet_ptLL",
    "Matched GenJetAK8 p_{T} (L)   (CMS private work);p_{T} [GeV]",
    40, 0, 1000
);
TH1F *h_matchedJet_ptTT = new TH1F(
    "h_matchedJet_ptTT",
    "Matched GenJetAK8 p_{T} (T)   (CMS private work);p_{T} [GeV]",
    40, 0, 1000
);

TH1F *h_matchedJet_massLL = new TH1F(
    "h_matchedJet_massLL",
    "Matched GenJetAK8 mass (L)   (CMS private work);Mass [GeV]",
    40, 10, 170
);
TH1F *h_matchedJet_massTT = new TH1F(
    "h_matchedJet_massTT",
    "Matched GenJetAK8 mass (T)   (CMS private work);Mass [GeV]",
    40, 10, 170
);

TH1F *h_matchedJet_etaLL = new TH1F(
    "h_matchedJet_etaLL",
    "Matched GenJetAK8 #eta (L)   (CMS private work);#eta",
    30, -4, 4
);
TH1F *h_matchedJet_etaTT = new TH1F(
    "h_matchedJet_etaTT",
    "Matched GenJetAK8 #eta (T)   (CMS private work);#eta",
    30, -4, 4
);

TH1F *h_matchedJet_phiLL = new TH1F(
    "h_matchedJet_phiLL",
    "Matched GenJetAK8 #phi (L)   (CMS private work);#phi",
    30, -4, 4
);
TH1F *h_matchedJet_phiTT = new TH1F(
    "h_matchedJet_phiTT",
    "Matched GenJetAK8 #phi (T)   (CMS private work);#phi",
    30, -4, 4
);

TH2F *h_Wpt_vs_jetptLL = new TH2F(
    "h_Wpt_vs_jetptLL",
    "Hadronic W p_{T} vs Matched Jet p_{T} (L)  (CMS private work);W p_{T} [GeV];Jet p_{T} [GeV]",
    50, 120, 1000,
    50, 120, 1000
);
TH2F *h_Wpt_vs_jetptTT = new TH2F(
    "h_Wpt_vs_jetptTT",
    "Hadronic W p_{T} vs Matched Jet p_{T} (T)  (CMS private work);W p_{T} [GeV];Jet p_{T} [GeV]",
    50, 120, 1000,
    50, 120, 1000
);

TH2F *h_dRqq_vs_ptLL = new TH2F("h_dRqq_vs_ptLL", "#DeltaR(q_{1},q_{2}) vs Matched Jet p_{T} (L)  (CMS private work); #DeltaR(q_{1},q_{2}); Jet p_{T} [GeV]", 
    50, 0, 1,
    50, 130, 1000
);
TH2F *h_dRqq_vs_ptTT = new TH2F("h_dRqq_vs_ptTT", "#DeltaR(q_{1},q_{2}) vs Matched Jet p_{T} (T)  (CMS private work); #DeltaR(q_{1},q_{2}); Jet p_{T} [GeV]", 
    50, 0, 1,
    50, 130, 1000
);





//subjet
TH1F *h_matchedSubJet1_ptLL = new TH1F(
    "h_matchedSubJet1_ptLL",
    "Matched #1 SubGenJetAK8 p_{T} (L)   (CMS private work);p_{T} [GeV]",
    30, 0, 1000
);
TH1F *h_matchedSubJet1_ptTT = new TH1F(
    "h_matchedSubJet1_ptTT",
    "Matched #1 SubGenJetAK8 p_{T} (T)   (CMS private work);p_{T} [GeV]",
    30, 0, 1000
);

TH1F *h_matchedSubJet1_massLL = new TH1F(
    "h_matchedSubJet1_massLL",
    "Matched #1 SubGenJetAK8 mass (L)   (CMS private work);Mass [GeV]",
    30, 0, 100
);
TH1F *h_matchedSubJet1_massTT = new TH1F(
    "h_matchedSubJet1_massTT",
    "Matched #1 SubGenJetAK8 mass (T)   (CMS private work);Mass [GeV]",
    30, 0, 100
);

TH1F *h_matchedSubJet2_ptLL = new TH1F(
    "h_matchedSubJet2_ptLL",
    "Matched #2 SubGenJetAK8 p_{T} (L)   (CMS private work);p_{T} [GeV]",
    30, 0, 1000
);
TH1F *h_matchedSubJet2_ptTT = new TH1F(
    "h_matchedSubJet2_ptTT",
    "Matched #2 SubGenJetAK8 p_{T} (T)   (CMS private work);p_{T} [GeV]",
    30, 0, 1000
);

TH1F *h_matchedSubJet2_massLL = new TH1F(
    "h_matchedSubJet2_massLL",
    "Matched #2 SubGenJetAK8 mass (L)   (CMS private work);Mass [GeV]",
    30, 0, 100
);
TH1F *h_matchedSubJet2_massTT = new TH1F(
    "h_matchedSubJet2_massTT",
    "Matched #2 SubGenJetAK8 mass (T)   (CMS private work);Mass [GeV]",
    30, 0, 100
);



    Int_t maxEvents = 2000000;
    std::string tagLL = "ssWWLL";
    std::string tagTT = "ssWWTT";

    TChain *chainLL = new TChain("Events");
    for(auto sample : samples_ssWWLL){
        chainLL->Add(sample.c_str());
    }
    TChain *chainTT = new TChain("Events");
    for(auto sample : samples_ssWWTT){
        chainTT->Add(sample.c_str());
    }

    


    const Int_t maxGenParticle = 10000;
    const Int_t maxGenFatJets = 100;

    Int_t nGenParticle;
    Int_t nGenJetAK8;
    const Int_t maxSubGenJets = 200;

    Int_t nSubGenJetAK8;
    Float_t SubGenJetAK8_pt[maxSubGenJets];
    Float_t SubGenJetAK8_eta[maxSubGenJets];
    Float_t SubGenJetAK8_phi[maxSubGenJets];
    Float_t SubGenJetAK8_mass[maxSubGenJets];

    Float_t GenPart_pt[maxGenParticle], GenPart_eta[maxGenParticle], GenPart_phi[maxGenParticle], GenPart_mass[maxGenParticle];
    Int_t GenPart_status[maxGenParticle];
    UShort_t GenPart_statusFlags[maxGenParticle];
    Short_t GenPart_genPartIdxMother[maxGenParticle];
    Int_t GenPart_pdgId[maxGenParticle];

    Float_t GenJetAK8_pt[maxGenFatJets], GenJetAK8_eta[maxGenFatJets], GenJetAK8_phi[maxGenFatJets], GenJetAK8_mass[maxGenFatJets];
    Short_t GenJetAK8_partonFlavour[maxGenFatJets]; 
    chainLL->SetBranchAddress("nGenPart", &nGenParticle);
    chainLL->SetBranchAddress("GenPart_pt", GenPart_pt);
    chainLL->SetBranchAddress("GenPart_eta", GenPart_eta);
    chainLL->SetBranchAddress("GenPart_phi", GenPart_phi);
    chainLL->SetBranchAddress("GenPart_mass", GenPart_mass);
    chainLL->SetBranchAddress("GenPart_status", GenPart_status);
    chainLL->SetBranchAddress("GenPart_statusFlags", GenPart_statusFlags);
    chainLL->SetBranchAddress("GenPart_genPartIdxMother", GenPart_genPartIdxMother);
    chainLL->SetBranchAddress("GenPart_pdgId", GenPart_pdgId);
    chainLL->SetBranchAddress("nSubGenJetAK8", &nSubGenJetAK8);
    chainLL->SetBranchAddress("SubGenJetAK8_pt", SubGenJetAK8_pt);
    chainLL->SetBranchAddress("SubGenJetAK8_eta", SubGenJetAK8_eta);
    chainLL->SetBranchAddress("SubGenJetAK8_phi", SubGenJetAK8_phi);
    chainLL->SetBranchAddress("SubGenJetAK8_mass", SubGenJetAK8_mass);
    chainLL->SetBranchAddress("nGenJetAK8", &nGenJetAK8);
    chainLL->SetBranchAddress("GenJetAK8_pt", GenJetAK8_pt);
    chainLL->SetBranchAddress("GenJetAK8_eta", GenJetAK8_eta);
    chainLL->SetBranchAddress("GenJetAK8_phi", GenJetAK8_phi);
    chainLL->SetBranchAddress("GenJetAK8_mass", GenJetAK8_mass);
    chainLL->SetBranchAddress("GenJetAK8_partonFlavour", GenJetAK8_partonFlavour);


    chainTT->SetBranchAddress("nGenPart", &nGenParticle);
    chainTT->SetBranchAddress("GenPart_pt", GenPart_pt);
    chainTT->SetBranchAddress("GenPart_eta", GenPart_eta);
    chainTT->SetBranchAddress("GenPart_phi", GenPart_phi);
    chainTT->SetBranchAddress("GenPart_mass", GenPart_mass);
    chainTT->SetBranchAddress("GenPart_status", GenPart_status);
    chainTT->SetBranchAddress("GenPart_statusFlags", GenPart_statusFlags);
    chainTT->SetBranchAddress("GenPart_genPartIdxMother", GenPart_genPartIdxMother);
    chainTT->SetBranchAddress("GenPart_pdgId", GenPart_pdgId);
    chainTT->SetBranchAddress("nSubGenJetAK8", &nSubGenJetAK8);
    chainTT->SetBranchAddress("SubGenJetAK8_pt", SubGenJetAK8_pt);
    chainTT->SetBranchAddress("SubGenJetAK8_eta", SubGenJetAK8_eta);
    chainTT->SetBranchAddress("SubGenJetAK8_phi", SubGenJetAK8_phi);
    chainTT->SetBranchAddress("SubGenJetAK8_mass", SubGenJetAK8_mass);
    chainTT->SetBranchAddress("nGenJetAK8", &nGenJetAK8);
    chainTT->SetBranchAddress("GenJetAK8_pt", GenJetAK8_pt);
    chainTT->SetBranchAddress("GenJetAK8_eta", GenJetAK8_eta);
    chainTT->SetBranchAddress("GenJetAK8_phi", GenJetAK8_phi);
    chainTT->SetBranchAddress("GenJetAK8_mass", GenJetAK8_mass);
    chainTT->SetBranchAddress("GenJetAK8_partonFlavour", GenJetAK8_partonFlavour);


//void funzione(TTree *t, Int_t evento){

    for(int iEvent = 0; iEvent < chainLL->GetEntries(); ++iEvent) {
        chainLL->GetEntry(iEvent);
        if(iEvent%10000==0){
            std::cout<<"Processing event L: "<<iEvent<<std::endl;
        }
        if(iEvent == maxEvents) break;
        std::vector<int> Wboson_indices;
        std::vector<int> quark_indices;
        std::pair<int,int> quarks_fromW;
        std::vector<int> leps_indices;
        std::pair<int, int> leps_fromW;
        int hadWboson = -1;
        int lepWboson = -1;
        int first_boson = 0, second_boson = 0;
    

        // tra le particelle generate cercate il bosone W^{+-} e salvarne l'indice dentro l'array 
        for(int gp = 0; gp < nGenParticle; ++gp) {
            //std::cout<<GenPart_status[gp]<<std::endl;
            if(std::abs(GenPart_pdgId[gp]) == 24 && GenPart_statusFlags[gp] == 10497 && GenPart_status[gp] == 62) {
                Wboson_indices.push_back(gp);
            }
        }
        
        if (Wboson_indices.size() < 2) continue;

        // tra i quark cercare quelli che hanno come indice della particella da cui provengono il W 
        for(int gpp = 0; gpp < nGenParticle; ++gpp) {
            if(std::abs(GenPart_pdgId[gpp]) <= 6) {
                if(GenPart_genPartIdxMother[gpp] == Wboson_indices[0]) { quark_indices.push_back(gpp); ++first_boson; }
                if(GenPart_genPartIdxMother[gpp] == Wboson_indices[1]) { quark_indices.push_back(gpp); ++second_boson; }
            }
        }

    
        // primo W è quello che decade adronico
        if(first_boson == 2 && second_boson == 0) {
            quarks_fromW = {quark_indices[0], quark_indices[1]};
            hadWboson = Wboson_indices[0];
            lepWboson = Wboson_indices[1];
            
            
        } 
        // secondo W è quello che decade adronico
        else if(first_boson == 0 && second_boson == 2) {
            quarks_fromW = {quark_indices[0], quark_indices[1]};
            hadWboson = Wboson_indices[1];
            lepWboson = Wboson_indices[0];
            
        }
        // W leptonico ->lv
        for(int gp=0; gp < nGenParticle; gp++){
            //std::cout<<"Indice gp: "<<gp<<std::endl;
            if((std::abs(GenPart_pdgId[gp]) >= 11) && (std::abs(GenPart_pdgId[gp]) <= 16)){
                //std::cout<<"gp ID: "<<std::abs(GenPart_pdgId[gp])<<std::endl;
                if(GenPart_genPartIdxMother[gp] == lepWboson){ leps_indices.push_back(gp);}
            }
        }
        if(leps_indices.size() >= 2){
            leps_fromW = {leps_indices[0], leps_indices[1]};
            //std::cout<<"lep_indices: "<<leps_fromW.first<<", "<<leps_fromW.second<<std::endl;
            //std::cout<<"leps ID: "<<std::abs(GenPart_pdgId[leps_fromW.first])<<", "<<std::abs(GenPart_pdgId[leps_fromW.second])<<std::endl;
            
            if(std::abs(GenPart_pdgId[leps_fromW.first]) == 12 || std::abs(GenPart_pdgId[leps_fromW.first]) == 14 || std::abs(GenPart_pdgId[leps_fromW.first]) == 16){ //cerco neutrino
                double dPhi1 = std::abs(GenPart_phi[leps_fromW.second] - GenPart_phi[lepWboson]);
                if(dPhi1 > M_PI) dPhi1 -= 2*M_PI;
                if(dPhi1 < -M_PI) dPhi1 += 2*M_PI;

                double Lp = (GenPart_pt[lepWboson]*GenPart_pt[leps_fromW.second]*cos(dPhi1))/(std::abs(GenPart_pt[lepWboson]*GenPart_pt[lepWboson]));
                h_Lp_LL->Fill(Lp);
                h_lep_ptLL->Fill(GenPart_pt[leps_fromW.second]);
                h_lep_etaLL->Fill(GenPart_eta[leps_fromW.second]);
                h_lep_phiLL->Fill(GenPart_phi[leps_fromW.second]);

                //boost di Lorentz
                TLorentzVector lep, W;
                lep.SetPtEtaPhiM(GenPart_pt[leps_fromW.second], GenPart_eta[leps_fromW.second], GenPart_phi[leps_fromW.second], GenPart_mass[leps_fromW.second]);
                W.SetPtEtaPhiM(GenPart_pt[lepWboson], GenPart_eta[lepWboson], GenPart_phi[lepWboson], GenPart_mass[lepWboson]);

                TVector3 boostW = -W.BoostVector();
                TLorentzVector lep_star = lep;
                lep_star.Boost(boostW); //leptone nel sdr del W
                TVector2 pt_lep_star(lep_star.Px(), lep_star.Py());
                TVector2 pt_W_lab(W.Px(), W.Py());
                double cos2D = (pt_lep_star*pt_W_lab)/(pt_lep_star.Mod()*pt_W_lab.Mod());
                h_cos2D_LL->Fill(cos2D);

                TVector3 p_lep_star(lep_star.Px(), lep_star.Py(), lep_star.Pz());
                TVector3 p_W_lab(W.Px(), W.Py(), W.Pz());
                double costheta_star = (p_lep_star.Dot(p_W_lab))/(p_lep_star.Mag()*p_W_lab.Mag());
                h_costheta_star_LL->Fill(costheta_star);

            }
            else{
                double dPhi1 = std::abs(GenPart_phi[leps_fromW.first] - GenPart_phi[lepWboson]);
                if(dPhi1 > M_PI) dPhi1 -= 2*M_PI;
                if(dPhi1 < -M_PI) dPhi1 += 2*M_PI;

                double Lp = (GenPart_pt[lepWboson]*GenPart_pt[leps_fromW.first]*cos(dPhi1))/(std::abs(GenPart_pt[lepWboson]*GenPart_pt[lepWboson]));
                h_Lp_LL->Fill(Lp);
                h_lep_ptLL->Fill(GenPart_pt[leps_fromW.first]);
                h_lep_etaLL->Fill(GenPart_eta[leps_fromW.first]);
                h_lep_phiLL->Fill(GenPart_phi[leps_fromW.first]);

                //boost di Lorentz
                TLorentzVector lep, W;
                lep.SetPtEtaPhiM(GenPart_pt[leps_fromW.first], GenPart_eta[leps_fromW.first], GenPart_phi[leps_fromW.first], GenPart_mass[leps_fromW.first]);
                W.SetPtEtaPhiM(GenPart_pt[lepWboson], GenPart_eta[lepWboson], GenPart_phi[lepWboson], GenPart_mass[lepWboson]);

                TVector3 boostW = -W.BoostVector();
                TLorentzVector lep_star = lep;
                lep_star.Boost(boostW); //leptone nel sdr del W
                TVector2 pt_lep_star(lep_star.Px(), lep_star.Py());
                TVector2 pt_W_lab(W.Px(), W.Py());
                double cos2D = (pt_lep_star*pt_W_lab)/(pt_lep_star.Mod()*pt_W_lab.Mod());
                h_cos2D_LL->Fill(cos2D);

                TVector3 p_lep_star(lep_star.Px(), lep_star.Py(), lep_star.Pz());
                TVector3 p_W_lab(W.Px(), W.Py(), W.Pz());
                double costheta_star = (p_lep_star.Dot(p_W_lab))/(p_lep_star.Mag()*p_W_lab.Mag());
                h_costheta_star_LL->Fill(costheta_star);

                }
            //double coss1 = 2*Lp1 -1;
            /*if(coss1 >1 || coss1 <-1){
                h_Lp1_LL->Fill(Lp1);
                h_lep1_massLL->Fill(GenPart_mass[leps_fromW.first]);
                h_lep1_ptLL->Fill(GenPart_pt[leps_fromW.first]);
            }
            double coss2 = 2*Lp2 -1;
            if(coss2 >1 || coss2 <-1){
                h_Lp2_LL->Fill(Lp2);
                h_lep2_massLL->Fill(GenPart_mass[leps_fromW.second]);
                h_lep2_ptLL->Fill(GenPart_pt[leps_fromW.second]);
            }*/

            
            

        }
        
    
   
    int matchedSubJetsThisEventLL = 0;
    int matchedJetsThisEventLL = 0;
    int matchedJetIndexLL=-1;
    int matchedSJIndex_1LL=0;
    int matchedSJIndex_2LL=0;
    std::vector<int> sgjIndex;



    // associazione GenJet - quarks, qui l'associazione non è univoca (non ci sono indici), si 
    // cerca il genjet che contiene i due quark trovati prima dentro 0.8  
    for(int ijet = 0; ijet < nGenJetAK8; ++ijet) {


    // variare questo per vedere come cambiano le distribuzioni 
        if(GenJetAK8_pt[ijet] < pt_lim) continue;
            
        // primo quark
        double dEta1, dPhi1;
        dEta1 = GenPart_eta[quarks_fromW.first] - GenJetAK8_eta[ijet];
        dPhi1 = GenPart_phi[quarks_fromW.first] - GenJetAK8_phi[ijet];
        if(dPhi1 > M_PI) dPhi1 -= 2*M_PI;     //deve stare nel range -pi/2, pi/2
        if(dPhi1 < -M_PI) dPhi1 += 2*M_PI;
        double dR1 = std::sqrt(dEta1*dEta1 + dPhi1*dPhi1);

        //secondo quark
        double dEta2, dPhi2;
        dEta2 = GenPart_eta[quarks_fromW.second] - GenJetAK8_eta[ijet];
        dPhi2 = GenPart_phi[quarks_fromW.second] - GenJetAK8_phi[ijet];
        if(dPhi2 > M_PI) dPhi2 -= 2*M_PI;
        if(dPhi2 < -M_PI) dPhi2 += 2*M_PI;
        double dR2 = std::sqrt(dEta2*dEta2 + dPhi2*dPhi2);

        // distanza tra i due quark
        double deta_qq = std::abs(GenPart_eta[quarks_fromW.first] - GenPart_eta[quarks_fromW.second]);
        double dphi_qq = GenPart_phi[quarks_fromW.first] - GenPart_phi[quarks_fromW.second];
        if(dphi_qq > M_PI) dphi_qq -= 2*M_PI;
        if(dphi_qq < - M_PI) dphi_qq += 2*M_PI;

        // richiesta che i due quark siano dentro il cono e che il pt del jet sia simile a quello del W 
        if(dR1 < 0.8 && dR2 < 0.8 && (std::abs(GenJetAK8_pt[ijet]/GenPart_pt[hadWboson]) < 1.10 && std::abs(GenJetAK8_pt[ijet]/GenPart_pt[hadWboson]) > 0.90)) {
            matchedJetIndexLL = ijet;
            matchedJetsThisEventLL++;
            h_dR_qqLL->Fill(std::sqrt(dphi_qq*dphi_qq + deta_qq*deta_qq));
            h_dR_q1_jetLL->Fill(dR1);
            h_dR_q2_jetLL->Fill(dR2);
            h_matchedJet_ptLL->Fill(GenJetAK8_pt[ijet]);
            h_matchedJet_massLL->Fill(GenJetAK8_mass[ijet]);
            h_matchedJet_etaLL->Fill(GenJetAK8_eta[ijet]);
            h_matchedJet_phiLL->Fill(GenJetAK8_phi[ijet]);
            h_Wpt_vs_jetptLL->Fill(
                GenPart_pt[hadWboson],
                GenJetAK8_pt[ijet]
            );
            h_dRqq_vs_ptLL->Fill(std::sqrt(dphi_qq*dphi_qq + deta_qq*deta_qq), GenJetAK8_pt[ijet]);

            break; 
        }
        
    }
    //h_nMatchedJetsLL->Fill(matchedJetsThisEventLL);

    if(matchedJetIndexLL >= 0){
        for(int isj=0; isj<nSubGenJetAK8; isj++){
            double dEta3 = std::abs(GenJetAK8_eta[matchedJetIndexLL] - SubGenJetAK8_eta[isj]);
            double dPhi3 = std::abs(GenJetAK8_phi[matchedJetIndexLL] - SubGenJetAK8_phi[isj]);
            if(dPhi3 > M_PI) dPhi3 -= 2*M_PI;
            if(dPhi3 < -M_PI) dPhi3 += 2*M_PI;
            double dR3 = std::sqrt(dEta3*dEta3 + dPhi3*dPhi3);

            if(dR3 < 0.8){
                sgjIndex.push_back(isj);
            }

        }
        if(sgjIndex.size() == 2){
            double pt_minv, pt_max;
            pt_max = std::max(SubGenJetAK8_pt[sgjIndex[0]], SubGenJetAK8_pt[sgjIndex[1]]);
            pt_minv = std::min(SubGenJetAK8_pt[sgjIndex[0]], SubGenJetAK8_pt[sgjIndex[1]]);
            h_zg_ratioLL->Fill(pt_max/pt_minv);
            h_zg_subjetsLL->Fill(pt_minv/(SubGenJetAK8_pt[sgjIndex[0]] + SubGenJetAK8_pt[sgjIndex[1]]));
            double E1, E2;
            E1=std::sqrt(SubGenJetAK8_mass[sgjIndex[0]]*SubGenJetAK8_mass[sgjIndex[0]] + SubGenJetAK8_pt[sgjIndex[0]]*SubGenJetAK8_pt[sgjIndex[0]]);
            E2=std::sqrt(SubGenJetAK8_mass[sgjIndex[1]]*SubGenJetAK8_mass[sgjIndex[1]] + SubGenJetAK8_pt[sgjIndex[1]]*SubGenJetAK8_pt[sgjIndex[1]]);
            h_ptheta_subjetsLL->Fill(std::abs(E1 - E2)/GenJetAK8_pt[matchedJetIndexLL]);
        }

        else if(sgjIndex.size() > 2){ //ci sono più di due sottojet: controllo per selezionare i due principali
            double bestSum = 0;
            std::pair<int, int> sgjDefIndex;
            bool foundPair = false;
            for(int i=0; i<sgjIndex.size(); i++){
                for(int j=i+1; j<sgjIndex.size(); j++){
                    double Sum = (SubGenJetAK8_pt[sgjIndex[i]] + SubGenJetAK8_pt[sgjIndex[j]]);
                    double rapp = (SubGenJetAK8_pt[sgjIndex[i]] + SubGenJetAK8_pt[sgjIndex[j]])/GenJetAK8_pt[matchedJetIndexLL];
                    if(Sum > bestSum){ //scelgo la migliore coppia
                        bestSum = Sum;
                        if(rapp > 0.8 && rapp < 1.20){
                            sgjDefIndex = {sgjIndex[i], sgjIndex[j]};
                            foundPair = true;
                        }
                    }                        
                }
            }
            if(foundPair != true) continue;
            double pt_minv, pt_max;
            pt_max = std::max(SubGenJetAK8_pt[sgjIndex[0]], SubGenJetAK8_pt[sgjIndex[1]]);
            pt_minv = std::min(SubGenJetAK8_pt[sgjIndex[0]], SubGenJetAK8_pt[sgjIndex[1]]);
            h_zg_ratioLL->Fill(pt_max/pt_minv);
            h_zg_subjetsLL->Fill(pt_minv/(SubGenJetAK8_pt[sgjDefIndex.first] + SubGenJetAK8_pt[sgjDefIndex.second]));
            double E1, E2;
            E1=std::sqrt(SubGenJetAK8_mass[sgjDefIndex.first]*SubGenJetAK8_mass[sgjDefIndex.first] + SubGenJetAK8_pt[sgjDefIndex.first]*SubGenJetAK8_pt[sgjDefIndex.first]);
            E2=std::sqrt(SubGenJetAK8_mass[sgjDefIndex.second]*SubGenJetAK8_mass[sgjDefIndex.second] + SubGenJetAK8_pt[sgjDefIndex.second]*SubGenJetAK8_pt[sgjDefIndex.second]);
            h_ptheta_subjetsLL->Fill(std::abs(E1 - E2)/GenJetAK8_pt[matchedJetIndexLL]);
                
        }
        else{
            continue;
        }
        
        double pt_min_q = std::min(GenPart_pt[quarks_fromW.first], GenPart_pt[quarks_fromW.second]);
        h_zg_quarksLL->Fill(pt_min_q/(GenPart_pt[quarks_fromW.first] + GenPart_pt[quarks_fromW.second]));
        double pt_max;
        pt_max = std::max(SubGenJetAK8_pt[sgjIndex[0]], SubGenJetAK8_pt[sgjIndex[1]]);
        h_zg_ratio_quarksLL->Fill(pt_max/pt_min_q);

        double E1_q, E2_q;
        E1_q = std::sqrt(GenPart_mass[quarks_fromW.first]*GenPart_mass[quarks_fromW.first] + GenPart_pt[quarks_fromW.first]*GenPart_pt[quarks_fromW.first]);
        E2_q = std::sqrt(GenPart_mass[quarks_fromW.second]*GenPart_mass[quarks_fromW.second] + GenPart_pt[quarks_fromW.second]*GenPart_pt[quarks_fromW.second]);
        h_ptheta_quarksLL->Fill(std::abs(E1_q - E2_q)/ GenPart_pt[hadWboson]);
    }
        


    for(int isjet = 0; isjet < nSubGenJetAK8; ++isjet) {

        // variare questo per vedere come cambiano le distribuzioni 
        if(SubGenJetAK8_pt[isjet] < pt_lim) continue;
            
            // primo quark
            double dEta1, dPhi1;
            dEta1 = GenPart_eta[quarks_fromW.first] - SubGenJetAK8_eta[isjet];
            dPhi1 = GenPart_phi[quarks_fromW.first] - SubGenJetAK8_phi[isjet];
            if(dPhi1 > M_PI) dPhi1 -= 2*M_PI;     //deve stare nel range -pi/2, pi/2
            if(dPhi1 < -M_PI) dPhi1 += 2*M_PI;
            double dR1 = std::sqrt(dEta1*dEta1 + dPhi1*dPhi1);

            //secondo quark
            double dEta2, dPhi2;
            dEta2 = GenPart_eta[quarks_fromW.second] - SubGenJetAK8_eta[isjet];
            dPhi2 = GenPart_phi[quarks_fromW.second] - SubGenJetAK8_phi[isjet];
            if(dPhi2 > M_PI) dPhi2 -= 2*M_PI;
            if(dPhi2 < -M_PI) dPhi2 += 2*M_PI;
            double dR2 = std::sqrt(dEta2*dEta2 + dPhi2*dPhi2);

            // richiesta che il primo quark sia dentro il cono e che il pt del primo subjet sia simile a quello del quark 
            if(dR1 < 0.8 && (std::abs(SubGenJetAK8_pt[isjet]/GenPart_pt[quarks_fromW.first]) < 1.10 && std::abs(SubGenJetAK8_pt[isjet]/GenPart_pt[quarks_fromW.first]) > 0.90)) {
                h_matchedSubJet1_ptLL->Fill(SubGenJetAK8_pt[isjet]);
                h_matchedSubJet1_massLL->Fill(SubGenJetAK8_mass[isjet]);
                matchedSJIndex_1LL=isjet;
                matchedSubJetsThisEventLL++;
                break; 
            }

            // richiesta che il secondo quark sia dentro il cono e che il pt del secondo subjet sia simile a quello del quark 
            if(dR2 < 0.8 && (std::abs(SubGenJetAK8_pt[isjet]/GenPart_pt[quarks_fromW.second]) < 1.10 && std::abs(SubGenJetAK8_pt[isjet]/GenPart_pt[quarks_fromW.second]) > 0.90)) {
                h_matchedSubJet2_ptLL->Fill(SubGenJetAK8_pt[isjet]);
                h_matchedSubJet2_massLL->Fill(SubGenJetAK8_mass[isjet]);
                matchedSJIndex_2LL=isjet;
                matchedSubJetsThisEventLL++;
                break; 
            }

        }
    }
    



    

    for(int iEvent = 0; iEvent < chainTT->GetEntries(); ++iEvent) {
        chainTT->GetEntry(iEvent);
        if(iEvent%10000==0){
            std::cout<<"Processing event T: "<<iEvent<<std::endl;
        }
        if(iEvent == maxEvents) break;
        std::vector<int> Wboson_indices;
        std::vector<int> quark_indices;
        std::pair<int,int> quarks_fromW;
        std::vector<int> leps_indices;
        std::pair<int, int> leps_fromW;
        int hadWboson = -1;
        int lepWboson = -1;
        int first_boson = 0, second_boson = 0;
    

        // tra le particelle generate cercate il bosone W^{+-} e salvarne l'indice dentro l'array 
        for(int gp = 0; gp < nGenParticle; ++gp) {
            //std::cout<<GenPart_status[gp]<<std::endl;
            if(std::abs(GenPart_pdgId[gp]) == 24 && GenPart_statusFlags[gp] == 10497 && GenPart_status[gp] == 62) {
                Wboson_indices.push_back(gp);
            }
        }
        
        // tra i quark cercare quelli che hanno come indice della particella da cui provengono il W 
        for(int gpp = 0; gpp < nGenParticle; ++gpp) {
            if(std::abs(GenPart_pdgId[gpp]) <= 6) {
                if(GenPart_genPartIdxMother[gpp] == Wboson_indices[0]) { quark_indices.push_back(gpp); ++first_boson; }
                if(GenPart_genPartIdxMother[gpp] == Wboson_indices[1]) { quark_indices.push_back(gpp); ++second_boson; }
            }
        }
        if (Wboson_indices.size() < 2) continue;
    
        // primo W è quello che decade adronico
        if(first_boson == 2 && second_boson == 0) {
            quarks_fromW = {quark_indices[0], quark_indices[1]};
            hadWboson = Wboson_indices[0];
            lepWboson = Wboson_indices[1];
        } 
        // secondo W è quello che decade adronico
        else if(first_boson == 0 && second_boson == 2) {
            quarks_fromW = {quark_indices[0], quark_indices[1]};
            hadWboson = Wboson_indices[1];
            lepWboson = Wboson_indices[0];
        }
        // W leptonico ->lv
        for(int gp=0; gp < nGenParticle; gp++){
            if((std::abs(GenPart_pdgId[gp]) >= 11) && (std::abs(GenPart_pdgId[gp]) <= 16)){
                if(GenPart_genPartIdxMother[gp] == lepWboson){ leps_indices.push_back(gp);}
            }
        }
        if(leps_indices.size() >= 2){
            leps_fromW = {leps_indices[0], leps_indices[1]};
            if(std::abs(GenPart_pdgId[leps_fromW.first]) == 12 || std::abs(GenPart_pdgId[leps_fromW.first]) == 14 || std::abs(GenPart_pdgId[leps_fromW.first]) == 16){ //cerco neutrino
                double dPhi1 = std::abs(GenPart_phi[leps_fromW.second] - GenPart_phi[lepWboson]);
                if(dPhi1 > M_PI) dPhi1 -= 2*M_PI;
                if(dPhi1 < -M_PI) dPhi1 += 2*M_PI;

                double Lp = (GenPart_pt[lepWboson]*GenPart_pt[leps_fromW.second]*cos(dPhi1))/(std::abs(GenPart_pt[lepWboson]*GenPart_pt[lepWboson]));
                h_Lp_TT->Fill(Lp);
                
                h_lep_ptTT->Fill(GenPart_pt[leps_fromW.second]);
                h_lep_etaTT->Fill(GenPart_eta[leps_fromW.second]);
                h_lep_phiTT->Fill(GenPart_phi[leps_fromW.second]);

                //boost di Lorentz
                TLorentzVector lep, W;
                lep.SetPtEtaPhiM(GenPart_pt[leps_fromW.second], GenPart_eta[leps_fromW.second], GenPart_phi[leps_fromW.second], GenPart_mass[leps_fromW.second]);
                W.SetPtEtaPhiM(GenPart_pt[lepWboson], GenPart_eta[lepWboson], GenPart_phi[lepWboson], GenPart_mass[lepWboson]);

                TVector3 boostW = -W.BoostVector();
                TLorentzVector lep_star = lep;
                lep_star.Boost(boostW); //leptone nel sdr del W
                TVector2 pt_lep_star(lep_star.Px(), lep_star.Py());
                TVector2 pt_W_lab(W.Px(), W.Py());
                double cos2D = (pt_lep_star*pt_W_lab)/(pt_lep_star.Mod()*pt_W_lab.Mod());
                h_cos2D_TT->Fill(cos2D);

                TVector3 p_lep_star(lep_star.Px(), lep_star.Py(), lep_star.Pz());
                TVector3 p_W_lab(W.Px(), W.Py(), W.Pz());
                double costheta_star = (p_lep_star.Dot(p_W_lab))/(p_lep_star.Mag()*p_W_lab.Mag());
                h_costheta_star_TT->Fill(costheta_star);

            }
            else{
                double dPhi1 = std::abs(GenPart_phi[leps_fromW.first] - GenPart_phi[lepWboson]);
                if(dPhi1 > M_PI) dPhi1 -= 2*M_PI;
                if(dPhi1 < -M_PI) dPhi1 += 2*M_PI;

                double Lp = (GenPart_pt[lepWboson]*GenPart_pt[leps_fromW.first]*cos(dPhi1))/(std::abs(GenPart_pt[lepWboson]*GenPart_pt[lepWboson]));
                h_Lp_TT->Fill(Lp);
                
                h_lep_ptTT->Fill(GenPart_pt[leps_fromW.first]);
                h_lep_etaTT->Fill(GenPart_eta[leps_fromW.first]);
                h_lep_phiTT->Fill(GenPart_phi[leps_fromW.first]);

                //boost di Lorentz
                TLorentzVector lep, W;
                lep.SetPtEtaPhiM(GenPart_pt[leps_fromW.first], GenPart_eta[leps_fromW.first], GenPart_phi[leps_fromW.first], GenPart_mass[leps_fromW.first]);
                W.SetPtEtaPhiM(GenPart_pt[lepWboson], GenPart_eta[lepWboson], GenPart_phi[lepWboson], GenPart_mass[lepWboson]);

                TVector3 boostW = -W.BoostVector();
                TLorentzVector lep_star = lep;
                lep_star.Boost(boostW); //leptone nel sdr del W
                TVector2 pt_lep_star(lep_star.Px(), lep_star.Py());
                TVector2 pt_W_lab(W.Px(), W.Py());
                double cos2D = (pt_lep_star*pt_W_lab)/(pt_lep_star.Mod()*pt_W_lab.Mod());
                h_cos2D_TT->Fill(cos2D);

                TVector3 p_lep_star(lep_star.Px(), lep_star.Py(), lep_star.Pz());
                TVector3 p_W_lab(W.Px(), W.Py(), W.Pz());
                double costheta_star = (p_lep_star.Dot(p_W_lab))/(p_lep_star.Mag()*p_W_lab.Mag());
                h_costheta_star_TT->Fill(costheta_star);
        
            }

        }
    
    int matchedSubJetsThisEventTT = 0;
    int matchedJetsThisEventTT = 0;
    int matchedJetIndexTT=-1;
    int matchedSJIndex_1TT=0;
    int matchedSJIndex_2TT=0;
    std::vector<int> sgjIndex;

    // associazione GenJet - quarks, qui l'associazione non è univoca (non ci sono indici), si 
    // cerca il genjet che contiene i due quark trovati prima dentro 0.8  
    for(int ijet = 0; ijet < nGenJetAK8; ++ijet) {


    // variare questo per vedere come cambiano le distribuzioni 
        if(GenJetAK8_pt[ijet] < pt_lim) continue;
            
        // primo quark
        double dEta1, dPhi1;
        dEta1 = GenPart_eta[quarks_fromW.first] - GenJetAK8_eta[ijet];
        dPhi1 = GenPart_phi[quarks_fromW.first] - GenJetAK8_phi[ijet];
        if(dPhi1 > M_PI) dPhi1 -= 2*M_PI;   
        if(dPhi1 < -M_PI) dPhi1 += 2*M_PI;
        double dR1 = std::sqrt(dEta1*dEta1 + dPhi1*dPhi1);

        //secondo quark
        double dEta2, dPhi2;
        dEta2 = GenPart_eta[quarks_fromW.second] - GenJetAK8_eta[ijet];
        dPhi2 = GenPart_phi[quarks_fromW.second] - GenJetAK8_phi[ijet];
        if(dPhi2 > M_PI) dPhi2 -= 2*M_PI;
        if(dPhi2 < -M_PI) dPhi2 += 2*M_PI;
        double dR2 = std::sqrt(dEta2*dEta2 + dPhi2*dPhi2);

        // distanza tra i due quark
        double deta_qq = std::abs(GenPart_eta[quarks_fromW.first] - GenPart_eta[quarks_fromW.second]);
        double dphi_qq = GenPart_phi[quarks_fromW.first] - GenPart_phi[quarks_fromW.second];
        if(dphi_qq > M_PI) dphi_qq -= 2*M_PI;
        if(dphi_qq < - M_PI) dphi_qq += 2*M_PI;

        // richiesta che i due quark siano dentro il cono e che il pt del jet sia simile a quello del W 
        if(dR1 < 0.8 && dR2 < 0.8 && (std::abs(GenJetAK8_pt[ijet]/GenPart_pt[hadWboson]) < 1.10 && std::abs(GenJetAK8_pt[ijet]/GenPart_pt[hadWboson]) > 0.90)) {
            matchedJetIndexTT = ijet;
            matchedJetsThisEventTT++;
            h_dR_qqTT->Fill(std::sqrt(dphi_qq*dphi_qq + deta_qq*deta_qq));
            h_dR_q1_jetTT->Fill(dR1);
            h_dR_q2_jetTT->Fill(dR2);
            h_matchedJet_ptTT->Fill(GenJetAK8_pt[ijet]);
            h_matchedJet_massTT->Fill(GenJetAK8_mass[ijet]);
            h_matchedJet_etaTT->Fill(GenJetAK8_eta[ijet]);
            h_matchedJet_phiTT->Fill(GenJetAK8_phi[ijet]);
            h_Wpt_vs_jetptTT->Fill(
                GenPart_pt[hadWboson],
                GenJetAK8_pt[ijet]
            );
            h_dRqq_vs_ptTT->Fill(std::sqrt(dphi_qq*dphi_qq + deta_qq*deta_qq), GenJetAK8_pt[ijet]);
            break; 
        }
    }
    h_nMatchedJetsTT->Fill(matchedJetsThisEventTT);

    if(matchedJetIndexTT >= 0){
           
        for(int isj=0; isj<nSubGenJetAK8; isj++){
            double dEta3 = std::abs(GenJetAK8_eta[matchedJetIndexTT] - SubGenJetAK8_eta[isj]);
            double dPhi3 = std::abs(GenJetAK8_phi[matchedJetIndexTT] - SubGenJetAK8_phi[isj]);
            if(dPhi3 > M_PI) dPhi3 -= 2*M_PI;
            if(dPhi3 < -M_PI) dPhi3 += 2*M_PI;
            double dR3 = std::sqrt(dEta3*dEta3 + dPhi3*dPhi3);

            if(dR3 < 0.8){
                sgjIndex.push_back(isj);
            }
        }
        if(sgjIndex.size() == 2){

            double pt_minv, pt_max;
            pt_max = std::max(SubGenJetAK8_pt[sgjIndex[0]], SubGenJetAK8_pt[sgjIndex[1]]);
            pt_minv = std::min(SubGenJetAK8_pt[sgjIndex[0]], SubGenJetAK8_pt[sgjIndex[1]]);
            h_zg_ratioTT->Fill(pt_max/pt_minv);
            h_zg_subjetsTT->Fill(pt_minv/(SubGenJetAK8_pt[sgjIndex[0]] + SubGenJetAK8_pt[sgjIndex[1]]));
            double E1, E2;
            E1=std::sqrt(SubGenJetAK8_mass[sgjIndex[0]]*SubGenJetAK8_mass[sgjIndex[0]] + SubGenJetAK8_pt[sgjIndex[0]]*SubGenJetAK8_pt[sgjIndex[0]]);
            E2=std::sqrt(SubGenJetAK8_mass[sgjIndex[1]]*SubGenJetAK8_mass[sgjIndex[1]] + SubGenJetAK8_pt[sgjIndex[1]]*SubGenJetAK8_pt[sgjIndex[1]]);
            h_ptheta_subjetsTT->Fill(std::abs(E1 - E2)/GenJetAK8_pt[matchedJetIndexTT]);
        }

        else if(sgjIndex.size() > 2){ //ci sono più di due sottojet: controllo per selezionare i due principali
            double bestSum = 0;
            std::pair<int, int> sgjDefIndex;
            bool foundPair = false;
            for(int i=0; i<sgjIndex.size(); i++){
                for(int j=i+1; j<sgjIndex.size(); j++){
                    double Sum = (SubGenJetAK8_pt[sgjIndex[i]] + SubGenJetAK8_pt[sgjIndex[j]]);
                    double rapp = (SubGenJetAK8_pt[sgjIndex[i]] + SubGenJetAK8_pt[sgjIndex[j]])/GenJetAK8_pt[matchedJetIndexTT];
                    if(Sum > bestSum){ //scelgo la migliore coppia
                        bestSum = Sum;
                        if(rapp > 0.8 && rapp < 1.20){
                            sgjDefIndex = {sgjIndex[i], sgjIndex[j]};
                            foundPair = true;
                        }
                    }                        
                }
            }
            if(foundPair != true) continue;
            double pt_minp, pt_max;
            pt_max = std::max(SubGenJetAK8_pt[sgjIndex[0]], SubGenJetAK8_pt[sgjIndex[1]]);
            pt_minp = std::min(SubGenJetAK8_pt[sgjIndex[0]], SubGenJetAK8_pt[sgjIndex[1]]);
            h_zg_ratioTT->Fill(pt_max/pt_minp);
            h_zg_subjetsTT->Fill(pt_minp/(SubGenJetAK8_pt[sgjDefIndex.first] + SubGenJetAK8_pt[sgjDefIndex.second]));
            double E1, E2;
            E1=std::sqrt(SubGenJetAK8_mass[sgjDefIndex.first]*SubGenJetAK8_mass[sgjDefIndex.first] + SubGenJetAK8_pt[sgjDefIndex.first]*SubGenJetAK8_pt[sgjDefIndex.first]);
            E2=std::sqrt(SubGenJetAK8_mass[sgjDefIndex.second]*SubGenJetAK8_mass[sgjDefIndex.second] + SubGenJetAK8_pt[sgjDefIndex.second]*SubGenJetAK8_pt[sgjDefIndex.second]);
            h_ptheta_subjetsTT->Fill(std::abs(E1 - E2)/GenJetAK8_pt[matchedJetIndexTT]);
        }
        else{ //sgjIndex.size() < 2
            continue;
        }

        double pt_min_q;
        pt_min_q = std::min(GenPart_pt[quarks_fromW.first], GenPart_pt[quarks_fromW.second]);
        h_zg_quarksTT->Fill(pt_min_q/(GenPart_pt[quarks_fromW.first] + GenPart_pt[quarks_fromW.second]));
        double pt_max;
        pt_max = std::max(SubGenJetAK8_pt[sgjIndex[0]], SubGenJetAK8_pt[sgjIndex[1]]);
        h_zg_ratio_quarksTT->Fill(pt_max/pt_min_q);

        double E1_q, E2_q;
        E1_q = std::sqrt(GenPart_mass[quarks_fromW.first]*GenPart_mass[quarks_fromW.first] + GenPart_pt[quarks_fromW.first]*GenPart_pt[quarks_fromW.first]);
        E2_q = std::sqrt(GenPart_mass[quarks_fromW.second]*GenPart_mass[quarks_fromW.second] + GenPart_pt[quarks_fromW.second]*GenPart_pt[quarks_fromW.second]);
        h_ptheta_quarksTT->Fill(std::abs(E1_q - E2_q)/GenPart_pt[hadWboson]);

    }

    



    for(int isjet = 0; isjet < nSubGenJetAK8; ++isjet) {

        // variare questo per vedere come cambiano le distribuzioni 
        if(SubGenJetAK8_pt[isjet] < pt_lim) continue;
            
            // primo quark
            double dEta1, dPhi1;
            dEta1 = GenPart_eta[quarks_fromW.first] - SubGenJetAK8_eta[isjet];
            dPhi1 = GenPart_phi[quarks_fromW.first] - SubGenJetAK8_phi[isjet];
            if(dPhi1 > M_PI) dPhi1 -= 2*M_PI;     //deve stare nel range -pi/2, pi/2
            if(dPhi1 < -M_PI) dPhi1 += 2*M_PI;
            double dR1 = std::sqrt(dEta1*dEta1 + dPhi1*dPhi1);

            //secondo quark
            double dEta2, dPhi2;
            dEta2 = GenPart_eta[quarks_fromW.second] - SubGenJetAK8_eta[isjet];
            dPhi2 = GenPart_phi[quarks_fromW.second] - SubGenJetAK8_phi[isjet];
            if(dPhi2 > M_PI) dPhi2 -= 2*M_PI;
            if(dPhi2 < -M_PI) dPhi2 += 2*M_PI;
            double dR2 = std::sqrt(dEta2*dEta2 + dPhi2*dPhi2);

            // richiesta che il primo quark sia dentro il cono e che il pt del primo subjet sia simile a quello del quark 
            if(dR1 < 0.8 && (std::abs(SubGenJetAK8_pt[isjet]/GenPart_pt[quarks_fromW.first]) < 1.10 && std::abs(SubGenJetAK8_pt[isjet]/GenPart_pt[quarks_fromW.first]) > 0.90)) {
                h_matchedSubJet1_ptTT->Fill(SubGenJetAK8_pt[isjet]);
                h_matchedSubJet1_massTT->Fill(SubGenJetAK8_mass[isjet]);
                matchedSJIndex_1TT=isjet;
                matchedSubJetsThisEventTT++;
                break; 
            }

            // richiesta che il secondo quark sia dentro il cono e che il pt del secondo subjet sia simile a quello del quark 
            if(dR2 < 0.8 && (std::abs(SubGenJetAK8_pt[isjet]/GenPart_pt[quarks_fromW.second]) < 1.10 && std::abs(SubGenJetAK8_pt[isjet]/GenPart_pt[quarks_fromW.second]) > 0.90)) {
                h_matchedSubJet2_ptTT->Fill(SubGenJetAK8_pt[isjet]);
                h_matchedSubJet2_massTT->Fill(SubGenJetAK8_mass[isjet]);
                matchedSJIndex_2TT=isjet;
                matchedSubJetsThisEventTT++;
                break; 
            }

        }

    }    

    std::vector<Color_t> colors;
    colors.push_back(kAzure-2);
    colors.push_back(kRed);

    std::vector<Color_t> colors2;
    colors2.push_back(kAzure-1);
    colors2.push_back(kRed+1);

    /*TCanvas *c_matchedJ_pt = new TCanvas("c_matchedJ_pt", "matchedJ pt", 1200, 600);
    std::vector<TH1*> matchedJ_pt;
    matchedJ_pt.push_back(h_matchedJet_ptLL);
    matchedJ_pt.push_back(h_matchedJet_ptTT);
    Makecanva2x1Bicolor(c_matchedJ_pt, matchedJ_pt, colors);
    c_matchedJ_pt->SaveAs("plots_pt/c_matchedJ_pt.pdf");



    TCanvas *c_matchedJ_mass = new TCanvas("c_matchedJ_mass", "matchedJ mass", 1200, 600);
    std::vector<TH1*> matchedJ_mass;
    matchedJ_mass.push_back(h_matchedJet_massLL);
    matchedJ_mass.push_back(h_matchedJet_massTT);
    Makecanva2x1Bicolor(c_matchedJ_mass, matchedJ_mass, colors);
    c_matchedJ_mass->SaveAs("plots_m/c_matchedJ_mass.pdf");



    TCanvas *c_matchedJ_eta = new TCanvas("c_matchedJ_eta", "matchedJ eta", 1200, 600);
    std::vector<TH1*> matchedJ_eta;
    matchedJ_eta.push_back(h_matchedJet_etaLL);
    matchedJ_eta.push_back(h_matchedJet_etaTT);
    Makecanva2x1Bicolor(c_matchedJ_eta, matchedJ_eta, colors);
    c_matchedJ_eta->SaveAs("plots_eta/c_matchedJ_eta.pdf");



    TCanvas *c_matchedJ_phi = new TCanvas("c_matchedJ_phi", "matchedJ phi", 1200, 600);
    std::vector<TH1*> matchedJ_phi;
    matchedJ_phi.push_back(h_matchedJet_phiLL);
    matchedJ_phi.push_back(h_matchedJet_phiTT);
    Makecanva2x1Bicolor(c_matchedJ_phi, matchedJ_phi, colors);
    c_matchedJ_phi->SaveAs("plots_phi/c_matchedJ_phi.pdf");*/



    /*TCanvas *c_W_gjet_pt = new TCanvas("c_W_gjet_pt","W vs jet p_{T}", 1500, 800);
    c_W_gjet_pt->Divide(2,1);
    c_W_gjet_pt->cd(1);
    gPad->SetBottomMargin(0.15);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.15);
    h_Wpt_vs_jetptLL->SetTitleSize(0.07);
    h_Wpt_vs_jetptLL->GetXaxis()->SetTitleSize(0.05);
    h_Wpt_vs_jetptLL->GetYaxis()->SetTitleSize(0.05);
    h_Wpt_vs_jetptLL->Scale(1.0 / h_Wpt_vs_jetptLL->Integral());
    h_Wpt_vs_jetptLL->SetStats(0);
    h_Wpt_vs_jetptLL->Draw("colz");

    c_W_gjet_pt->cd(2);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.15);
    gPad->SetBottomMargin(0.15);
    h_Wpt_vs_jetptTT->SetTitleSize(0.07);
    h_Wpt_vs_jetptTT->GetXaxis()->SetTitleSize(0.05);
    h_Wpt_vs_jetptTT->GetYaxis()->SetTitleSize(0.05);
    h_Wpt_vs_jetptTT->Scale(1.0 / h_Wpt_vs_jetptTT->Integral());
    h_Wpt_vs_jetptTT->SetStats(0);
    h_Wpt_vs_jetptTT->Draw("colz");
    c_W_gjet_pt->SaveAs("plots_pt/c_W_gjet_pt.pdf");*/




    /*TCanvas *c_dRqq = new TCanvas("c_dRqq","#DeltaR qq", 1200, 600);
    std::vector<TH1*> dRqq;
    dRqq.push_back(h_dR_qqLL);
    dRqq.push_back(h_dR_qqTT);
    Makecanva2x1Bicolor(c_dRqq, dRqq, colors);
    c_dRqq->SaveAs("plots_pt/c_dRqq.pdf");




    TCanvas *c_dRqq_vs_pt = new TCanvas("c_dRqq_vs_pt","W vs jet p_{T}", 1500, 800);
    c_dRqq_vs_pt->Divide(2,1);
    c_dRqq_vs_pt->cd(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.15);
    gPad->SetBottomMargin(0.15);
    h_dRqq_vs_ptLL->SetTitleSize(0.05);
    h_dRqq_vs_ptLL->Scale(1.0 / h_dRqq_vs_ptLL->Integral());
    h_dRqq_vs_ptLL->SetStats(0);
    h_dRqq_vs_ptLL->GetYaxis()->SetTitleSize(0.05);
    h_dRqq_vs_ptLL->Draw("colz");

    c_dRqq_vs_pt->cd(2);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.15);
    gPad->SetBottomMargin(0.15);
    h_dRqq_vs_ptTT->SetTitleSize(0.05);
    h_dRqq_vs_ptTT->Scale(1.0 / h_dRqq_vs_ptTT->Integral());
    h_dRqq_vs_ptTT->SetStats(0);
    h_dRqq_vs_ptTT->GetYaxis()->SetTitleSize(0.05);
    h_dRqq_vs_ptTT->Draw("colz");
    c_dRqq_vs_pt->SaveAs("plots_pt/c_dRqq_vs_pt.pdf");*/



    /*TCanvas *c_matchedSJ_pt = new TCanvas("c_matchedSJ_pt", "Matched Subjets pt", 1400, 1000);
    std::vector<TH1*> matchedSJ_pt;
    matchedSJ_pt.push_back(h_matchedSubJet1_ptLL);
    matchedSJ_pt.push_back(h_matchedSubJet1_ptTT);
    matchedSJ_pt.push_back(h_matchedSubJet2_ptLL);
    matchedSJ_pt.push_back(h_matchedSubJet2_ptTT);
    Makecanva2x2Bicolor(c_matchedSJ_pt, matchedSJ_pt, colors);
    c_matchedSJ_pt->SaveAs("plots_sj/c_matchedSJ_pt.pdf");



    TCanvas *c_matchedSJ_mass = new TCanvas("c_matchedSJ_mass", "Matched Subjets mass", 1400, 1000);
    std::vector<TH1*> matchedSJ_mass;
    matchedSJ_mass.push_back(h_matchedSubJet1_massLL);
    matchedSJ_mass.push_back(h_matchedSubJet1_massTT);
    matchedSJ_mass.push_back(h_matchedSubJet2_massLL);
    matchedSJ_mass.push_back(h_matchedSubJet2_massTT);
    Makecanva2x2Bicolor(c_matchedSJ_mass, matchedSJ_mass, colors);
    c_matchedSJ_mass->SaveAs("plots_sj/c_matchedSJ_mass.pdf");*/


    
    TCanvas *c_zg_sgj = new TCanvas("c_zg_sgj","SubGenJet z_{g}", 1200, 600);
    h_zg_subjetsLL->SetTitleSize(0.07);
    h_zg_subjetsTT->SetTitleSize(0.07);
    std::vector<TH1*> zg_sgj;
    zg_sgj.push_back(h_zg_subjetsLL);
    zg_sgj.push_back(h_zg_subjetsTT);
    Makecanva2x1Bicolor(c_zg_sgj, zg_sgj, colors);
    c_zg_sgj->SaveAs("plots_z1/c_zg_sgj.pdf");

    TCanvas *c_zg_sgj_pol = new TCanvas("c_zg_sgj_pol", "c_zg_sgj_pol", 800, 600);
    h_zg_subjetsLL->SetTitleSize(0.07);
    h_zg_subjetsLL->SetTitle(Form("z_{g} subjets  p_{T} > %d   (CMS private work)", pt_lim));
    h_zg_subjetsLL->SetStats(0);
    h_zg_subjetsTT->SetTitleSize(0.07);
    std::vector<TH1*> zg_sgj_pol;
    zg_sgj_pol.push_back(h_zg_subjetsLL);
    zg_sgj_pol.push_back(h_zg_subjetsTT);
    MakecanvaColor(c_zg_sgj_pol, zg_sgj_pol, colors2);
    TLegend *leg2 = new TLegend(0.4, 0.74, 0.68, 0.88);
    leg2->AddEntry(h_zg_subjetsLL, "L pol", "f");
    leg2->AddEntry(h_zg_subjetsTT, "T pol", "f");
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->SetMargin(0.25);
    leg2->Draw();
    c_zg_sgj_pol->SaveAs("plots_z1/c_zg_sgj_pol.pdf");

    /*TCanvas *c_zg_ratio = new TCanvas("c_zg_ratio","SubGenJet z_{g}'", 1200, 600);
    h_zg_ratioLL->SetTitleSize(0.07);
    h_zg_ratioTT->SetTitleSize(0.07);
    std::vector<TH1*> zg_ratio_pol;
    zg_ratio_pol.push_back(h_zg_ratioLL);
    zg_ratio_pol.push_back(h_zg_ratioTT);
    Makecanva2x1Bicolor(c_zg_ratio, zg_ratio_pol, colors);
    c_zg_ratio->SaveAs("plots_z1/c_zg_ratio.pdf");*/


   
   
    TCanvas *c_zg_q = new TCanvas("c_zg_q","GenParticle Quarks z_{g}", 1200, 600);
    std::vector<TH1*> zg_q;
    zg_q.push_back(h_zg_quarksLL);
    zg_q.push_back(h_zg_quarksTT);
    Makecanva2x1Bicolor(c_zg_q, zg_q, colors2);
    c_zg_q->SaveAs("plots_z1/c_zg_q.pdf");

    /*TCanvas *c_zg_ratio_q = new TCanvas("c_zg_ratio_q","GenParticle Quarks z_{g}'", 1200, 600);
    std::vector<TH1*> zg_ratio_q;
    zg_ratio_q.push_back(h_zg_ratio_quarksLL);
    zg_ratio_q.push_back(h_zg_ratio_quarksTT);
    Makecanva2x1Bicolor(c_zg_ratio_q, zg_ratio_q, colors);
    c_zg_ratio_q->SaveAs("plots_z1/c_zg_ratio_q.pdf");*/




    TCanvas *c_ptheta_sj = new TCanvas("c_ptheta_sj","SubGenJet p_{#theta}", 1200, 600);
    h_ptheta_subjetsLL->SetTitleSize(0.07);
    h_ptheta_subjetsTT->SetTitleSize(0.07);
    std::vector<TH1*> ptheta_sj;
    ptheta_sj.push_back(h_ptheta_subjetsLL);
    ptheta_sj.push_back(h_ptheta_subjetsTT);
    Makecanva2x1Bicolor(c_ptheta_sj, ptheta_sj, colors);
    c_ptheta_sj->SaveAs("plots_ptheta/c_ptheta_sj.pdf");

    TCanvas *c_ptheta_sj_pol = new TCanvas("c_ptheta_sj_pol", "c_ptheta_sj_pol", 800, 600);
    h_ptheta_subjetsLL->SetTitleSize(0.07);
    h_ptheta_subjetsLL->SetStats(0);
    h_ptheta_subjetsLL->SetTitle(Form("p_{#theta} subjets  p_{T} > %d   (CMS private work)", pt_lim));
    h_ptheta_subjetsTT->SetTitleSize(0.07);
    std::vector<TH1*> ptheta_sj_pol;
    ptheta_sj_pol.push_back(h_ptheta_subjetsLL);
    ptheta_sj_pol.push_back(h_ptheta_subjetsTT);
    MakecanvaColor(c_ptheta_sj_pol, ptheta_sj_pol, colors2);
    TLegend *leg3 = new TLegend(0.4, 0.74, 0.68, 0.88);
    leg3->AddEntry(h_ptheta_subjetsLL, "L pol", "f");
    leg3->AddEntry(h_ptheta_subjetsTT, "T pol", "f");
    leg3->SetBorderSize(0);
    leg3->SetFillStyle(0);
    leg3->SetMargin(0.25);
    leg3->Draw();
    c_ptheta_sj_pol->SaveAs("plots_ptheta/c_ptheta_sj_pol.pdf");



    TCanvas *c_ptheta_q = new TCanvas("c_ptheta_q","SubGenJet p_{#theta}", 1200, 600);
    std::vector<TH1*> ptheta_q;
    ptheta_q.push_back(h_ptheta_quarksLL);
    ptheta_q.push_back(h_ptheta_quarksTT);
    Makecanva2x1Bicolor(c_ptheta_q, ptheta_q, colors2);
    c_ptheta_q->SaveAs("plots_ptheta/c_ptheta_q.pdf");



    gStyle->SetOptStat(0);

    /*TCanvas *c_dR_q_jet = new TCanvas("c_dR_q_jet", "#DeltaR q_1 - jet", 1200, 600);
    std::vector<TH1*> dR_q_jet;
    h_dR_q1_jetLL->SetTitleSize(0.08);
    h_dR_q2_jetLL->SetTitleSize(0.08);
    dR_q_jet.push_back(h_dR_q1_jetLL);
    dR_q_jet.push_back(h_dR_q1_jetTT);
    dR_q_jet.push_back(h_dR_q2_jetLL);
    dR_q_jet.push_back(h_dR_q2_jetTT);
    Makecanva2x1Bicolor(c_dR_q_jet, dR_q_jet, colors2);
    TLegend *leg9 = new TLegend(0.6, 0.74, 0.88, 0.88);
    leg9->AddEntry(h_dR_q2_jetLL, "L pol", "f");
    leg9->AddEntry(h_dR_q2_jetTT, "T pol", "f");
    leg9->SetBorderSize(0);
    leg9->SetFillStyle(0);
    leg9->SetMargin(0.25);
    leg9->Draw();

    c_dR_q_jet->SaveAs("plots_pt/c_dR_q_jet.pdf");*/




    gStyle->SetOptStat(0);

    /*TCanvas *c_LepP = new TCanvas("c_LepP", "Lepton projetion", 900, 700);
    std::vector<TH1*> LepP;
    LepP.push_back(h_Lp_TT);
    LepP.push_back(h_Lp_LL);
    Makecanva(c_LepP, LepP, kGreen-4);
    TLegend *leg1 = new TLegend(0.2, 0.7, 0.35, 0.8);
    leg1->AddEntry(h_Lp_LL, "L pol", "f");
    leg1->AddEntry(h_Lp_TT, "T pol", "f");
    leg1->SetBorderSize(0);
    leg1->SetFillStyle(0);
    leg1->SetMargin(0.25);
    leg1->Draw();
    c_LepP->SaveAs("plots_lp/c_LepP_gen.pdf");

    

    TCanvas *c_lep_pt = new TCanvas("c_lep_pt", "Lepton p_{T}", 900, 700);
    std::vector<TH1*> lep_pt;
    lep_pt.push_back(h_lep_ptLL);
    lep_pt.push_back(h_lep_ptTT);
    Makecanva(c_lep_pt, lep_pt, kGreen-4);
    TLegend *leg4 = new TLegend(0.5, 0.7, 0.65, 0.8);
    leg4->AddEntry(h_lep_ptLL, "L pol", "f");
    leg4->AddEntry(h_lep_ptTT, "T pol", "f");
    leg4->SetBorderSize(0);
    leg4->SetFillStyle(0);
    leg4->SetMargin(0.25);
    leg4->Draw();
    c_lep_pt->SaveAs("plots_lp/c_lep_pt_gen.pdf");

    TCanvas *c_lep_eta = new TCanvas("c_lep_eta", "Lepton #eta", 900, 700);
    std::vector<TH1*> lep_eta;
    lep_eta.push_back(h_lep_etaTT);
    lep_eta.push_back(h_lep_etaLL);
    Makecanva(c_lep_eta, lep_eta, kGreen-4);
    TLegend *leg5 = new TLegend(0.7, 0.75, 0.85, 0.85);
    leg5->AddEntry(h_lep_etaLL, "L pol", "f");
    leg5->AddEntry(h_lep_etaTT, "T pol", "f");
    leg5->SetBorderSize(0);
    leg5->SetFillStyle(0);
    leg5->SetMargin(0.25);
    leg5->Draw();
    c_lep_eta->SaveAs("plots_lp/c_lep_eta_gen.pdf");


    TCanvas *c_lep_phi = new TCanvas("c_lep_phi", "Lepton #eta", 900, 700);
    std::vector<TH1*> lep_phi;
    lep_phi.push_back(h_lep_phiLL);
    lep_phi.push_back(h_lep_phiTT);
    Makecanva(c_lep_phi, lep_phi, kGreen-4);
    TLegend *leg6 = new TLegend(0.45, 0.4, 0.6, 0.5);
    leg6->AddEntry(h_lep_phiLL, "L pol", "f");
    leg6->AddEntry(h_lep_phiTT, "T pol", "f");
    leg6->SetBorderSize(0);
    leg6->SetMargin(0.25);
    leg6->Draw();
    c_lep_phi->SaveAs("plots_lp/c_lep_phi_gen.pdf");*/


    TCanvas *c_cos2D = new TCanvas("c_cos2D", "Lepton projetions", 900, 700);
    std::vector<TH1*> cos2D;
    cos2D.push_back(h_cos2D_TT);
    cos2D.push_back(h_cos2D_LL);
    Makecanva(c_cos2D, cos2D, kGreen-4);
    TLegend *leg7 = new TLegend(0.45, 0.7, 0.6, 0.8);
    leg7->AddEntry(h_cos2D_LL, "L pol", "f");
    leg7->AddEntry(h_cos2D_TT, "T pol", "f");
    leg7->SetBorderSize(0);
    leg7->SetFillStyle(0);
    leg7->SetMargin(0.25);
    leg7->Draw();
    c_cos2D->SaveAs("plots_lp/c_cos2D_gen.pdf");


    TCanvas *c_costheta_star = new TCanvas("c_costheta_gen", "Lepton projetions", 900, 700);
    std::vector<TH1*> costheta_star;
    costheta_star.push_back(h_costheta_star_TT);
    costheta_star.push_back(h_costheta_star_LL);
    Makecanva(c_costheta_star, costheta_star, kGreen-4);
    TLegend *leg8 = new TLegend(0.2, 0.75, 0.35, 0.85);
    leg8->AddEntry(h_costheta_star_LL, "L pol", "f");
    leg8->AddEntry(h_costheta_star_TT, "T pol", "f");
    leg8->SetBorderSize(0);
    leg8->SetFillStyle(0);
    leg8->SetMargin(0.25);
    leg8->Draw();
    c_costheta_star->SaveAs("plots_lp/c_costheta_star_gen.pdf");




    f->Write();
    f->Close();




    gROOT->ProcessLine(".q");

}
