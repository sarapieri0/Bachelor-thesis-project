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


void GenJetAK8_match(
    TString output="Gen_match.root"
){
    TFile *f=new TFile(output, "recreate");
int pt_lim = 300;

TH1F *h_zg_quarksLL = new TH1F("h_zg_quarksLL", "z_{g} quarks  (L); z_{g} ; Events", 30,  0, 0.5);
TH1F *h_zg_subjetsLL = new TH1F("h_zg_subjetsLL", Form("z_{g} subjets p_{t} > %d   (L); z_{g} ; Events", pt_lim), 30, 0, 0.5);
TH1F *h_ptheta_quarksLL = new TH1F("h_ptheta_quarksLL", "p_{#theta} quarks  (L); p_{#theta} ; Events",30, 0, 1);
TH1F *h_ptheta_subjetsLL = new TH1F("h_ptheta_subjetsLL", Form("p_{#theta} subjets p_{t} > %d   (L); p_{#theta} ; Events", pt_lim), 30, 0, 1);

TH1F *h_zg_quarksTT = new TH1F("h_zg_quarksTT", "z_{g} quarks  (T); z_{g} ; Events", 30,  0, 0.5);
TH1F *h_zg_subjetsTT = new TH1F("h_zg_subjetsTT", Form("z_{g} subjets p_{t} > %d   (T); z_{g} ; Events", pt_lim), 30, 0, 0.5);
TH1F *h_ptheta_quarksTT = new TH1F("h_ptheta_quarksTT", "p_{#theta} quarks  (T); p_{#theta} ; Events",30, 0, 1);
TH1F *h_ptheta_subjetsTT = new TH1F("h_ptheta_subjetsTT", Form("p_{#theta} subjets p_{t} > %d   (T); p_{#theta} ; Events", pt_lim), 30, 0, 1);

TH1F *h_Lp_LL = new TH1F("h_Lp1_L", "Lepton projection (L); L_{P}; Events", 30, 0, 1.2);
TH1F *h_Lp_TT = new TH1F("h_Lp1_T", "Lepton projection (T); L_{P}; Events", 30, 0, 1.2);


//leptons
TH1F *h_lep_ptLL = new TH1F("h_lep1_ptL", "Lepton p_{T} (L); p_{T} [GeV]; Events", 50, 0, 400);
TH1F *h_lep_ptTT = new TH1F("h_lep1_ptT", "Lepton p_{T} (T); p_{T} [GeV]; Events", 50, 0, 400);

TH1F *h_lep_etaLL = new TH1F("h_lep_etaLL", "Lepton #eta (L); #eta; Events", 50, -6, 6);
TH1F *h_lep_etaTT = new TH1F("h_lep_etaTT", "Lepton #eta (T); #eta; Events", 50, -6, 6);

TH1F *h_lep_phiLL = new TH1F("h_lep_phiLL", "Lepton #phi (L); #phi; Events", 50, -4, 4);
TH1F *h_lep_phiTT = new TH1F("h_lep_phiTT", "Lepton #phi (T); #phi; Events", 50, 4, 4);



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
    Form("#DeltaR(quark_{1}, GenJetAK8) (L) p_{t} > %d;#DeltaR", pt_lim),
    50, 0, 1
);
TH1F *h_dR_q1_jetTT = new TH1F(
    "h_dR_q1_jetTT",
    Form("#DeltaR(quark_{1}, GenJetAK8) (T) p_{t} > %d;#DeltaR", pt_lim),
    50, 0, 1
);

TH1F *h_dR_qqLL = new TH1F("h_dR_qqLL", Form("#DeltaR qq p_{t} > %d     (L); #DeltaR", pt_lim), 50, 0, 1);
TH1F *h_dR_qqTT = new TH1F("h_dR_qqTT", Form("#DeltaR qq p_{t} > %d     (T); #DeltaR", pt_lim), 50, 0, 1);

TH1F *h_dR_q2_jetLL = new TH1F(
    "h_dR_q2_jetLL",
    Form("#DeltaR(quark_{2}, GenJetAK8) (L) p_{t} > %d;#DeltaR", pt_lim),
    50, 0, 1
);
TH1F *h_dR_q2_jetTT = new TH1F(
    "h_dR_q2_jetTT",
    Form("#DeltaR(quark_{2}, GenJetAK8) (T) p_{t} > %d;#DeltaR", pt_lim),
    50, 0, 1
);

TH1F *h_matchedJet_ptLL = new TH1F(
    "h_matchedJet_ptLL",
    "Matched GenJetAK8 p_{T} (L);p_{T} [GeV]",
    50, 0, 1000
);
TH1F *h_matchedJet_ptTT = new TH1F(
    "h_matchedJet_ptTT",
    "Matched GenJetAK8 p_{T} (T);p_{T} [GeV]",
    50, 0, 1000
);

TH1F *h_matchedJet_massLL = new TH1F(
    "h_matchedJet_massLL",
    "Matched GenJetAK8 mass (L);Mass [GeV]",
    50, 0, 300
);
TH1F *h_matchedJet_massTT = new TH1F(
    "h_matchedJet_massTT",
    "Matched GenJetAK8 mass (T);Mass [GeV]",
    50, 0, 300
);

TH1F *h_matchedJet_etaLL = new TH1F(
    "h_matchedJet_etaLL",
    "Matched GenJetAK8 #eta (L);#eta",
    40, -4, 4
);
TH1F *h_matchedJet_etaTT = new TH1F(
    "h_matchedJet_etaTT",
    "Matched GenJetAK8 #eta (T);#eta",
    40, -4, 4
);

TH1F *h_matchedJet_phiLL = new TH1F(
    "h_matchedJet_phiLL",
    "Matched GenJetAK8 #phi (L);#phi",
    40, -4, 4
);
TH1F *h_matchedJet_phiTT = new TH1F(
    "h_matchedJet_phiTT",
    "Matched GenJetAK8 #phi (T);#phi",
    40, -4, 4
);

TH2F *h_Wpt_vs_jetptLL = new TH2F(
    "h_Wpt_vs_jetptLL",
    "Hadronic W p_{T} vs Matched Jet p_{T} (L);W p_{T} [GeV];Jet p_{T} [GeV]",
    50, 0, 1000,
    50, 0, 1000
);
TH2F *h_Wpt_vs_jetptTT = new TH2F(
    "h_Wpt_vs_jetptTT",
    "Hadronic W p_{T} vs Matched Jet p_{T} (T);W p_{T} [GeV];Jet p_{T} [GeV]",
    50, 0, 1000,
    50, 0, 1000
);

TH2F *h_dRqq_vs_ptLL = new TH2F("h_dRqq_vs_ptLL", "qq #DeltaR vs Matched Jet p_{t} (L); #DeltaR; Jet p_{t} [GeV]", 
    50, 0, 1,
    50, 0, 1000
);
TH2F *h_dRqq_vs_ptTT = new TH2F("h_dRqq_vs_ptTT", "qq #DeltaR vs Matched Jet p_{t} (T); #DeltaR; Jet p_{t} [GeV]", 
    50, 0, 1,
    50, 0, 1000
);





//subjet
TH1F *h_matchedSubJet1_ptLL = new TH1F(
    "h_matchedSubJet1_ptLL",
    "Matched #1 SubGenJetAK8 p_{T} (L);p_{T} [GeV]",
    30, 0, 1000
);
TH1F *h_matchedSubJet1_ptTT = new TH1F(
    "h_matchedSubJet1_ptTT",
    "Matched #1 SubGenJetAK8 p_{T} (T);p_{T} [GeV]",
    30, 0, 1000
);

TH1F *h_matchedSubJet1_massLL = new TH1F(
    "h_matchedSubJet1_massLL",
    "Matched #1 SubGenJetAK8 mass (L);Mass [GeV]",
    30, 0, 100
);
TH1F *h_matchedSubJet1_massTT = new TH1F(
    "h_matchedSubJet1_massTT",
    "Matched #1 SubGenJetAK8 mass (T);Mass [GeV]",
    30, 0, 100
);

TH1F *h_matchedSubJet2_ptLL = new TH1F(
    "h_matchedSubJet2_ptLL",
    "Matched #2 SubGenJetAK8 p_{T} (L);p_{T} [GeV]",
    30, 0, 1000
);
TH1F *h_matchedSubJet2_ptTT = new TH1F(
    "h_matchedSubJet2_ptTT",
    "Matched #2 SubGenJetAK8 p_{T} (T);p_{T} [GeV]",
    30, 0, 1000
);

TH1F *h_matchedSubJet2_massLL = new TH1F(
    "h_matchedSubJet2_massLL",
    "Matched #2 SubGenJetAK8 mass (L);Mass [GeV]",
    30, 0, 100
);
TH1F *h_matchedSubJet2_massTT = new TH1F(
    "h_matchedSubJet2_massTT",
    "Matched #2 SubGenJetAK8 mass (T);Mass [GeV]",
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


void funzione(TTree *t, Int_t evento){

    for(int iEvent = 0; iEvent < chainLL->GetEntries(); ++iEvent) {
        t->GetEntry(evento);
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
        
    
 /*  
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
    h_nMatchedJetsLL->Fill(matchedJetsThisEventLL);

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
            double pt_minv;
            pt_minv = std::min(SubGenJetAK8_pt[sgjIndex[0]], SubGenJetAK8_pt[sgjIndex[1]]);
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
            double pt_minp;
            pt_minp = std::min(SubGenJetAK8_pt[sgjDefIndex.first], SubGenJetAK8_pt[sgjDefIndex.second]);
            h_zg_subjetsLL->Fill(pt_minp/(SubGenJetAK8_pt[sgjDefIndex.first] + SubGenJetAK8_pt[sgjDefIndex.second]));
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

        }*/
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
        
            }

        }
    /*
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

            double pt_minv;
            pt_minv = std::min(SubGenJetAK8_pt[sgjIndex[0]], SubGenJetAK8_pt[sgjIndex[1]]);
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
            double pt_minp;
            pt_minp = std::min(SubGenJetAK8_pt[sgjDefIndex.first], SubGenJetAK8_pt[sgjDefIndex.second]);
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

        }*/

    }    

    /*
    TCanvas *c_matchedJ_pt = new TCanvas("c_matchedJ_pt", "matchedJ pt", 1400, 800);
    c_matchedJ_pt->Divide(2,1);
    c_matchedJ_pt->cd(1);
    h_matchedJet_ptLL->Scale(1.0 / h_matchedJet_ptLL->Integral());
    h_matchedJet_ptLL->Draw("hist");
    h_matchedJet_ptLL->SetLineColor(kBlue+1);
    h_matchedJet_ptLL->SetFillColor(kBlue+1);

    c_matchedJ_pt->cd(2);
    h_matchedJet_ptTT->Scale(1.0 / h_matchedJet_ptTT->Integral());
    h_matchedJet_ptTT->Draw("hist");
    h_matchedJet_ptTT->SetLineColor(kRed+1);
    h_matchedJet_ptTT->SetFillColor(kRed+1);
    c_matchedJ_pt->SaveAs("plots_pt/c_matchedJ_pt.pdf");



    TCanvas *c_matchedJ_mass = new TCanvas("c_matchedJ_mass", "matchedJ mass", 1400, 800);
    c_matchedJ_mass->Divide(2,1);
    c_matchedJ_mass->cd(1);
    h_matchedJet_massLL->Scale(1.0 / h_matchedJet_massLL->Integral());
    h_matchedJet_massLL->Draw("hist");
    h_matchedJet_massLL->SetLineColor(kBlue+1);
    h_matchedJet_massLL->SetFillColor(kBlue+1);

    c_matchedJ_mass->cd(2);
    h_matchedJet_massTT->Scale(1.0 / h_matchedJet_massTT->Integral());
    h_matchedJet_massTT->Draw("hist");
    h_matchedJet_massTT->SetLineColor(kRed+1);
    h_matchedJet_massTT->SetFillColor(kRed+1);
    c_matchedJ_mass->SaveAs("plots_m/c_matchedJ_mass.pdf");



    TCanvas *c_matchedJ_eta = new TCanvas("c_matchedJ_eta", "matchedJ eta", 1400, 800);
    c_matchedJ_eta->Divide(2,1);
    c_matchedJ_eta->cd(1);
    gPad->SetLeftMargin(0.15);
    h_matchedJet_etaLL->Scale(1.0 / h_matchedJet_etaLL->Integral());
    h_matchedJet_etaLL->Draw("hist");
    h_matchedJet_etaLL->SetLineColor(kBlue+1);
    h_matchedJet_etaLL->SetFillColor(kBlue+1);
    h_matchedJet_etaLL->GetXaxis()->SetTitleSize(0.05);
    gPad->SetBottomMargin(0.15);

    c_matchedJ_eta->cd(2);
    gPad->SetLeftMargin(0.15);
    h_matchedJet_etaTT->Scale(1.0 / h_matchedJet_etaTT->Integral());
    h_matchedJet_etaTT->Draw("hist");
    h_matchedJet_etaTT->SetLineColor(kRed+1);
    h_matchedJet_etaTT->SetFillColor(kRed+1);
    h_matchedJet_etaTT->GetXaxis()->SetTitleSize(0.05);
    gPad->SetBottomMargin(0.15);
    c_matchedJ_eta->SaveAs("plots_eta/c_matchedJ_eta.pdf");



    TCanvas *c_matchedJ_phi = new TCanvas("c_matchedJ_phi", "matchedJ phi", 1400, 800);
    c_matchedJ_phi->Divide(2,1);
    c_matchedJ_phi->cd(1);
    h_matchedJet_phiLL->Scale(1.0 / h_matchedJet_phiLL->Integral());
    h_matchedJet_phiLL->Draw("hist");
    h_matchedJet_phiLL->SetLineColor(kBlue+1);
    h_matchedJet_phiLL->SetFillColor(kBlue+1);
    h_matchedJet_phiLL->GetXaxis()->SetTitleSize(0.05);

    c_matchedJ_phi->cd(2);
    h_matchedJet_phiTT->Scale(1.0 / h_matchedJet_phiTT->Integral());
    h_matchedJet_phiTT->Draw("hist");
    h_matchedJet_phiTT->SetLineColor(kRed+1);
    h_matchedJet_phiTT->SetFillColor(kRed+1);
    h_matchedJet_phiTT->GetXaxis()->SetTitleSize(0.05);
    c_matchedJ_phi->SaveAs("plots_phi/c_matchedJ_phi.pdf");*/



    /*TCanvas *c_W_gjet_pt = new TCanvas("c_W_gjet_pt","W vs jet p_{t}", 1500, 800);
    c_W_gjet_pt->Divide(2,1);
    c_W_gjet_pt->cd(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.15);
    h_Wpt_vs_jetptLL->Scale(1.0 / h_Wpt_vs_jetptLL->Integral());
    h_Wpt_vs_jetptLL->SetStats(0);
    h_Wpt_vs_jetptLL->Draw("colz");

    c_W_gjet_pt->cd(2);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.15);
    h_Wpt_vs_jetptTT->Scale(1.0 / h_Wpt_vs_jetptTT->Integral());
    h_Wpt_vs_jetptTT->SetStats(0);
    h_Wpt_vs_jetptTT->Draw("colz");
    c_W_gjet_pt->SaveAs("plots_pt/c_W_gjet_pt.pdf");*/




    TCanvas *c_dRqq = new TCanvas("c_dRqq","#DeltaR qq", 1300, 800);
    c_dRqq->Divide(2,1);
    c_dRqq->cd(1);
    h_dR_qqLL->Scale(1.0 / h_dR_qqLL->Integral());
    h_dR_qqLL->Draw("hist");
    h_dR_qqLL->SetLineColor(kBlue+1);
    h_dR_qqLL->SetFillColor(kBlue+1);
    h_dR_qqLL->GetXaxis()->SetTitleSize(0.05);

    c_dRqq->cd(2);
    h_dR_qqTT->Scale(1.0 / h_dR_qqTT->Integral());
    h_dR_qqTT->Draw("hist");
    h_dR_qqTT->SetLineColor(kRed+1);
    h_dR_qqTT->SetFillColor(kRed+1);
    h_dR_qqTT->GetXaxis()->SetTitleSize(0.05);

    c_dRqq->SaveAs("plots_pt/c_dRqq.pdf");




    /*TCanvas *c_dRqq_vs_pt = new TCanvas("c_dRqq_vs_pt","W vs jet p_{t}", 1500, 800);
    c_dRqq_vs_pt->Divide(2,1);
    c_dRqq_vs_pt->cd(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.15);
    h_dRqq_vs_ptLL->Scale(1.0 / h_dRqq_vs_ptLL->Integral());
    h_dRqq_vs_ptLL->SetStats(0);
    h_dRqq_vs_ptLL->GetXaxis()->SetTitleSize(0.05);
    h_dRqq_vs_ptLL->GetYaxis()->SetTitleSize(0.05);
    h_dRqq_vs_ptLL->Draw("colz");

    c_dRqq_vs_pt->cd(2);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.15);
    h_dRqq_vs_ptTT->Scale(1.0 / h_dRqq_vs_ptTT->Integral());
    h_dRqq_vs_ptTT->SetStats(0);
    h_dRqq_vs_ptTT->GetXaxis()->SetTitleSize(0.05);
    h_dRqq_vs_ptTT->GetYaxis()->SetTitleSize(0.05);
    h_dRqq_vs_ptTT->Draw("colz");
    c_dRqq_vs_pt->SaveAs("plots_pt/c_dRqq_vs_pt.pdf");*/



    /*TCanvas *c_matchedSJ_pt = new TCanvas("c_matchedSJ_pt", "Matched Subjets pt", 1400, 1000);
    c_matchedSJ_pt->Divide(2,2);

    c_matchedSJ_pt->cd(1);
    gPad->SetBottomMargin(0.15);
    h_matchedSubJet1_ptLL->Scale(1.0/h_matchedSubJet1_ptLL->Integral());
    h_matchedSubJet1_ptLL->Draw("hist");
    h_matchedSubJet1_ptLL->SetLineColor(kBlue+1);
    h_matchedSubJet1_ptLL->SetFillColor(kBlue+1);
    h_matchedSubJet1_ptLL->GetXaxis()->SetTitleSize(0.05);
    //h_matchedSubJet1_ptLL->GetXaxis()->SetTitleOffset(0.09);

    c_matchedSJ_pt->cd(2);
    gPad->SetBottomMargin(0.15);
    h_matchedSubJet1_ptTT->Scale(1.0/h_matchedSubJet1_ptTT->Integral());
    h_matchedSubJet1_ptTT->Draw("hist");
    h_matchedSubJet1_ptTT->SetLineColor(kRed+1);
    h_matchedSubJet1_ptTT->SetFillColor(kRed+1);
    h_matchedSubJet1_ptTT->GetXaxis()->SetTitleSize(0.05);
    //h_matchedSubJet1_ptTT->GetXaxis()->SetTitleOffset(0.09);

    c_matchedSJ_pt->cd(3);
    gPad->SetBottomMargin(0.15);
    h_matchedSubJet2_ptLL->Scale(1.0/h_matchedSubJet2_ptLL->Integral());
    h_matchedSubJet2_ptLL->Draw("hist");
    h_matchedSubJet2_ptLL->SetLineColor(kBlue+1);
    h_matchedSubJet2_ptLL->SetFillColor(kBlue+1);
    h_matchedSubJet2_ptLL->GetXaxis()->SetTitleSize(0.05);
    //h_matchedSubJet2_ptLL->GetXaxis()->SetTitleOffset(0.09);

    c_matchedSJ_pt->cd(4);
    gPad->SetBottomMargin(0.15);
    h_matchedSubJet2_ptTT->Scale(1.0/h_matchedSubJet2_ptTT->Integral());
    h_matchedSubJet2_ptTT->Draw("hist");
    h_matchedSubJet2_ptTT->SetLineColor(kRed+1);
    h_matchedSubJet2_ptTT->SetFillColor(kRed+1);
    h_matchedSubJet2_ptTT->GetXaxis()->SetTitleSize(0.05);
    //h_matchedSubJet2_ptTT->GetXaxis()->SetTitleOffset(0.09);
    c_matchedSJ_pt->SaveAs("plots_sj/c_matchedSJ_pt.pdf");



    TCanvas *c_matchedSJ_mass = new TCanvas("c_matchedSJ_mass", "Matched Subjets mass", 1400, 1000);
    c_matchedSJ_mass->Divide(2,2);

    c_matchedSJ_mass->cd(1);
    h_matchedSubJet1_massLL->Scale(1.0/h_matchedSubJet1_massLL->Integral());
    h_matchedSubJet1_massLL->Draw("hist");
    h_matchedSubJet1_massLL->SetLineColor(kBlue+1);
    h_matchedSubJet1_massLL->SetFillColor(kBlue+1);
    h_matchedSubJet1_massLL->GetXaxis()->SetTitleSize(0.05);
    //h_matchedSubJet1_massLL->GetXaxis()->SetTitleOffset(0.09);

    c_matchedSJ_mass->cd(2);
    h_matchedSubJet1_massTT->Scale(1.0/h_matchedSubJet1_massTT->Integral());
    h_matchedSubJet1_massTT->Draw("hist");
    h_matchedSubJet1_massTT->SetLineColor(kRed+1);
    h_matchedSubJet1_massTT->SetFillColor(kRed+1);
    h_matchedSubJet1_massTT->GetXaxis()->SetTitleSize(0.05);
    //h_matchedSubJet1_massTT->GetXaxis()->SetTitleOffset(0.09);
    c_matchedSJ_pt->SaveAs("plots_sj/c_matchedSJ_pt.pdf");

    c_matchedSJ_mass->cd(3);
    h_matchedSubJet2_massLL->Scale(1.0/h_matchedSubJet2_massLL->Integral());
    h_matchedSubJet2_massLL->Draw("hist");
    h_matchedSubJet2_massLL->SetLineColor(kBlue+1);
    h_matchedSubJet2_massLL->SetFillColor(kBlue+1);
    h_matchedSubJet2_massLL->GetXaxis()->SetTitleSize(0.05);
    //h_matchedSubJet2_massLL->GetXaxis()->SetTitleOffset(0.09);

    c_matchedSJ_mass->cd(4);
    h_matchedSubJet2_massTT->Scale(1.0/h_matchedSubJet2_massTT->Integral());
    h_matchedSubJet2_massTT->Draw("hist");
    h_matchedSubJet2_massTT->SetLineColor(kRed+1);
    h_matchedSubJet2_massTT->SetFillColor(kRed+1);
    h_matchedSubJet2_massTT->GetXaxis()->SetTitleSize(0.05);
    //h_matchedSubJet2_massTT->GetXaxis()->SetTitleOffset(0.09);
    c_matchedSJ_mass->SaveAs("plots_sj/c_matchedSJ_mass.pdf");*/


    
    TCanvas *c_zg_sgj = new TCanvas("c_zg_sgj","SubGenJet z_{g}", 1500, 900);
    c_zg_sgj->Divide(2,1);

    c_zg_sgj->cd(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    h_zg_subjetsLL->Scale(1.0 / h_zg_subjetsLL->Integral());
    h_zg_subjetsLL->Draw("hist");
    h_zg_subjetsLL->SetLineColor(kBlue+1);
    h_zg_subjetsLL->SetFillColor(kBlue+1);
    h_zg_subjetsLL->GetXaxis()->SetTitleSize(0.05);


    c_zg_sgj->cd(2);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    h_zg_subjetsTT->Scale(1.0 / h_zg_subjetsTT->Integral());
    h_zg_subjetsTT->Draw("hist");
    h_zg_subjetsTT->SetLineColor(kRed+1);
    h_zg_subjetsTT->SetFillColor(kRed+1);
    h_zg_subjetsTT->GetXaxis()->SetTitleSize(0.05);
    c_zg_sgj->SaveAs("plots_z1/c_zg_sgj.pdf");


   
   
    TCanvas *c_zg_q = new TCanvas("c_zg_q","GenParticle Quarks z_{g}", 1500, 900);
    c_zg_q->Divide(2,1);

    c_zg_q->cd(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    h_zg_quarksLL->Scale(1.0 / h_zg_quarksLL->Integral());
    h_zg_quarksLL->Draw("hist");
    h_zg_quarksLL->SetLineColor(kBlue+1);
    h_zg_quarksLL->SetFillColor(kBlue+1);
    h_zg_quarksLL->GetXaxis()->SetTitleSize(0.05);

    c_zg_q->cd(2);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    h_zg_quarksTT->Scale(1.0 / h_zg_quarksTT->Integral());
    h_zg_quarksTT->Draw("hist");
    h_zg_quarksTT->SetLineColor(kRed+1);
    h_zg_quarksTT->SetFillColor(kRed+1);
    h_zg_quarksTT->GetXaxis()->SetTitleSize(0.05);

    c_zg_q->SaveAs("plots_z1/c_zg_q.pdf");




    TCanvas *c_ptheta_sj = new TCanvas("c_ptheta_sj","SubGenJet p_{#theta}", 1500, 900);
    c_ptheta_sj->Divide(2,1);

    c_ptheta_sj->cd(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    h_ptheta_subjetsLL->Scale(1.0 / h_ptheta_subjetsLL->Integral());
    h_ptheta_subjetsLL->Draw("hist");
    h_ptheta_subjetsLL->SetLineColor(kBlue+1);
    h_ptheta_subjetsLL->SetFillColor(kBlue+1);
    h_ptheta_subjetsLL->GetXaxis()->SetTitleSize(0.05);


    c_ptheta_sj->cd(2);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    h_ptheta_subjetsTT->Scale(1.0 / h_ptheta_subjetsTT->Integral());
    h_ptheta_subjetsTT->Draw("hist");
    h_ptheta_subjetsTT->SetLineColor(kRed+1);
    h_ptheta_subjetsTT->SetFillColor(kRed+1);
    h_ptheta_subjetsTT->GetXaxis()->SetTitleSize(0.05);

    c_ptheta_sj->SaveAs("plots_ptheta/c_ptheta_sj.pdf");



    TCanvas *c_ptheta_q = new TCanvas("c_ptheta_q","SubGenJet p_{#theta}", 1500, 900);
    c_ptheta_q->Divide(2,1);

    c_ptheta_q->cd(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    h_ptheta_quarksLL->Scale(1.0 / h_ptheta_quarksLL->Integral());
    h_ptheta_quarksLL->Draw("hist");
    h_ptheta_quarksLL->SetLineColor(kBlue+1);
    h_ptheta_quarksLL->SetFillColor(kBlue+1);
    h_ptheta_quarksLL->GetXaxis()->SetTitleSize(0.05);


    c_ptheta_q->cd(2);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    h_ptheta_quarksTT->Scale(1.0 / h_ptheta_quarksTT->Integral());
    h_ptheta_quarksTT->Draw("hist");
    h_ptheta_quarksTT->SetLineColor(kRed+1);
    h_ptheta_quarksTT->SetFillColor(kRed+1);
    h_ptheta_quarksTT->GetXaxis()->SetTitleSize(0.05);

    c_ptheta_q->SaveAs("plots_ptheta/c_ptheta_q.pdf");




    TCanvas *c_dR_q_jet = new TCanvas("c_dR_q_jet", "#DeltaR q_1 - jet", 1300, 1000);
    c_dR_q_jet->Divide(2,2);

    c_dR_q_jet->cd(1);
    h_dR_q1_jetLL->Scale(1.0 / h_dR_q1_jetLL->Integral());
    h_dR_q1_jetLL->Draw("hist");
    h_dR_q1_jetLL->SetLineColor(kBlue+1);
    h_dR_q1_jetLL->SetFillColor(kBlue+1);
    h_dR_q1_jetLL->GetXaxis()->SetTitleSize(0.05);


    c_dR_q_jet->cd(2);
    h_dR_q1_jetTT->Scale(1.0 / h_dR_q1_jetTT->Integral());
    h_dR_q1_jetTT->Draw("hist");
    h_dR_q1_jetTT->SetLineColor(kRed+1);
    h_dR_q1_jetTT->SetFillColor(kRed+1);
    h_dR_q1_jetTT->GetXaxis()->SetTitleSize(0.05);


    c_dR_q_jet->cd(3);
    h_dR_q2_jetLL->Scale(1.0 / h_dR_q2_jetLL->Integral());
    h_dR_q2_jetLL->Draw("hist");
    h_dR_q2_jetLL->SetLineColor(kBlue+1);
    h_dR_q2_jetLL->SetFillColor(kBlue+1);
    h_dR_q2_jetLL->GetXaxis()->SetTitleSize(0.05);


    c_dR_q_jet->cd(4);
    h_dR_q2_jetTT->Scale(1.0 / h_dR_q2_jetTT->Integral());
    h_dR_q2_jetTT->Draw("hist");
    h_dR_q2_jetTT->SetLineColor(kRed+1);
    h_dR_q2_jetTT->SetFillColor(kRed+1);
    h_dR_q2_jetTT->GetXaxis()->SetTitleSize(0.05);

    c_dR_q_jet->SaveAs("plots_pt/c_dR_q_jet.pdf");



    TCanvas *c_NmatchedJets = new TCanvas("c_NmatchedJets", "N matched Jets per event", 1500, 900);
    c_NmatchedJets->Divide(2,1);

    c_NmatchedJets->cd(1);
    h_nMatchedJetsLL->Scale(1.0 / h_nMatchedJetsLL->Integral());
    h_nMatchedJetsLL->Draw("hist");
    h_nMatchedJetsLL->SetLineColor(kBlue+1);
    h_nMatchedJetsLL->SetFillColor(kBlue+1);
    h_nMatchedJetsLL->GetXaxis()->SetTitleSize(0.05);


    c_NmatchedJets->cd(2);
    h_nMatchedJetsTT->Scale(1.0 / h_nMatchedJetsTT->Integral());
    h_nMatchedJetsTT->Draw("hist");
    h_nMatchedJetsTT->SetLineColor(kRed+1);
    h_nMatchedJetsTT->SetFillColor(kRed+1);
    h_nMatchedJetsTT->GetXaxis()->SetTitleSize(0.05);
    c_NmatchedJets->SaveAs("plots_nmatch/c_NmatchedJets.pdf");


    gStyle->SetOptStat(0);
/*
    TCanvas *c_LepP = new TCanvas("c_LepP", "Lepton projetions", 900, 700);
    h_Lp_LL->SetTitle("Lepton projetions");
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    h_Lp_LL->Scale(1.0 / h_Lp_LL->Integral());
    h_Lp_LL->Draw("hist");
    h_Lp_LL->SetLineColor(kGreen-4);
    h_Lp_LL->SetFillColor(kGreen-4);
    h_Lp_TT->Scale(1.0 / h_Lp_TT->Integral());
    h_Lp_TT->SetMarkerColor(kBlack);
    h_Lp_TT->SetLineColor(kBlack);
    h_Lp_TT->SetMarkerStyle(20);
    h_Lp_TT->SetMarkerSize(1.5);
    h_Lp_TT->Draw("P same");
    h_Lp_LL->GetXaxis()->SetTitleSize(0.05);


    TLegend *leg1 = new TLegend(0.75, 0.75, 0.9, 0.9);
    leg1->AddEntry(h_Lp_LL, "L pol", "l");
    leg1->AddEntry(h_Lp_TT, "T pol", "l");
    leg1->Draw();

    c_LepP->SaveAs("plots_lp/c_LepP.pdf");

    

    TCanvas *c_lep_pt = new TCanvas("c_lep_pt", "Leptons p_{t}", 900, 700);

    h_lep_ptLL->SetTitle("Lepton p_{t}");
    gPad->SetBottomMargin(0.15);
    h_lep_ptLL->Scale(1.0 / h_lep_ptLL->Integral());
    h_lep_ptLL->Draw("hist");
    h_lep_ptLL->SetLineColor(kGreen-4);
    h_lep_ptLL->SetFillColor(kGreen-4);
    h_lep_ptLL->GetXaxis()->SetTitleSize(0.05);
    h_lep_ptTT->Scale(1.0 / h_lep_ptTT->Integral());
    h_lep_ptTT->SetMarkerColor(kBlack);
    h_lep_ptTT->SetLineColor(kBlack);
    h_lep_ptTT->SetMarkerStyle(20);
    h_lep_ptTT->SetMarkerSize(1.5);
    h_lep_ptTT->Draw("P same");

    TLegend *leg4 = new TLegend(0.75, 0.75, 0.9, 0.9);
    leg4->AddEntry(h_lep_ptLL, "L pol", "l");
    leg4->AddEntry(h_lep_ptTT, "T pol", "l");
    leg4->Draw();

    c_lep_pt->SaveAs("plots_lp/c_lep_pt.pdf");*/

    TCanvas *c_lep_eta = new TCanvas("c_lep_eta", "Leptons #eta", 900, 700);
    h_lep_etaLL->SetTitle("Lepton pseudorapidity #eta");
    gPad->SetBottomMargin(0.15);
    h_lep_etaLL->Scale(1.0 / h_lep_etaLL->Integral());
    h_lep_etaLL->Draw("hist");
    h_lep_etaLL->SetLineColor(kGreen-4);
    h_lep_etaLL->SetFillColor(kGreen-4);
    h_lep_etaLL->GetXaxis()->SetTitleSize(0.05);
    h_lep_etaTT->Scale(1.0 / h_lep_etaTT->Integral());
    h_lep_etaTT->SetMarkerColor(kBlack);
    h_lep_etaTT->SetLineColor(kBlack);
    h_lep_etaTT->SetMarkerStyle(20);
    h_lep_etaTT->SetMarkerSize(1.5);
    h_lep_etaTT->Draw("P same");

    TLegend *leg5 = new TLegend(0.75, 0.75, 0.9, 0.9);
    leg5->AddEntry(h_lep_etaLL, "L pol", "l");
    leg5->AddEntry(h_lep_etaTT, "T pol", "l");
    leg5->Draw();

    c_lep_eta->SaveAs("plots_lp/c_lep_eta.pdf");


    TCanvas *c_lep_phi = new TCanvas("c_lep_phi", "Leptons #eta", 900, 700);
    h_lep_phiLL->SetTitle("Lepton #phi");
    gPad->SetBottomMargin(0.15);
    h_lep_phiLL->Scale(1.0 / h_lep_phiLL->Integral());
    h_lep_phiLL->Draw("hist");
    h_lep_phiLL->SetLineColor(kGreen-4);
    h_lep_phiLL->SetFillColor(kGreen-4);
    h_lep_phiLL->GetXaxis()->SetTitleSize(0.05);
    h_lep_phiTT->Scale(1.0 / h_lep_phiTT->Integral());
    h_lep_phiTT->SetMarkerColor(kBlack);
    h_lep_phiTT->SetLineColor(kBlack);
    h_lep_phiTT->SetMarkerStyle(20);
    h_lep_phiTT->SetMarkerSize(1.5);
    h_lep_phiTT->Draw("P same");

    TLegend *leg6 = new TLegend(0.75, 0.75, 0.9, 0.9);
    leg6->AddEntry(h_lep_phiLL, "L pol", "l");
    leg6->AddEntry(h_lep_phiTT, "T pol", "l");
    leg6->Draw();

    c_lep_phi->SaveAs("plots_lp/c_lep_phi.pdf");






    f->Write();
    f->Close();




    gROOT->ProcessLine(".q");

}
