void funzione(TTree *tree, Int_t evento){

    int pt_lim=150;

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
            h_Lp->Fill(Lp);

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
            h_cos2D->Fill(cos2D);

            TVector3 p_lep_star(lep_star.Px(), lep_star.Py(), lep_star.Pz());
            TVector3 p_W_lab(W.Px(), W.Py(), W.Pz());
            double costheta_star = (p_lep_star.Dot(p_W_lab))/(p_lep_star.Mag()*p_W_lab.Mag());
            h_costheta_star->Fill(costheta_star);

            //h_lep_pt->Fill(GenPart_pt[leps_fromW.second]);
            //h_lep_eta->Fill(GenPart_eta[leps_fromW.second]);
            //h_lep_phi->Fill(GenPart_phi[leps_fromW.second]);

        }
        else{
            double dPhi1 = std::abs(GenPart_phi[leps_fromW.first] - GenPart_phi[lepWboson]);
            if(dPhi1 > M_PI) dPhi1 -= 2*M_PI;
            if(dPhi1 < -M_PI) dPhi1 += 2*M_PI;

            double Lp = (GenPart_pt[lepWboson]*GenPart_pt[leps_fromW.first]*cos(dPhi1))/(std::abs(GenPart_pt[lepWboson]*GenPart_pt[lepWboson]));
            h_Lp->Fill(Lp);

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
            h_cos2D->Fill(cos2D);

            TVector3 p_lep_star(lep_star.Px(), lep_star.Py(), lep_star.Pz());
            TVector3 p_W_lab(W.Px(), W.Py(), W.Pz());
            double costheta_star = (p_lep_star.Dot(p_W_lab))/(p_lep_star.Mag()*p_W_lab.Mag());
            h_costheta_star->Fill(costheta_star);
           // h_lep_pt->Fill(GenPart_pt[leps_fromW.first]);
           // h_lep_eta->Fill(GenPart_eta[leps_fromW.first]);
           // h_lep_phi->Fill(GenPart_phi[leps_fromW.first]);

            }          
        

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
        if(GenJetAK8_pt[ijet] > pt_lim){
            
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
            //h_dR_qq->Fill(std::sqrt(dphi_qq*dphi_qq + deta_qq*deta_qq));
            //h_dR_q1_jetLL->Fill(dR1);
            //h_dR_q2_jetLL->Fill(dR2);
            h_matchedJet_pt->Fill(GenJetAK8_pt[ijet]);
            h_matchedJet_mass->Fill(GenJetAK8_mass[ijet]);
            h_matchedJet_eta->Fill(GenJetAK8_eta[ijet]);
            h_matchedJet_phi->Fill(GenJetAK8_phi[ijet]);
            /*h_Wpt_vs_jetpt->Fill(
                GenPart_pt[hadWboson],
                GenJetAK8_pt[ijet]
            );*/
            //h_dRqq_vs_ptLL->Fill(std::sqrt(dphi_qq*dphi_qq + deta_qq*deta_qq), GenJetAK8_pt[ijet]);

            break; 
        }
    }
        
    }

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
            pt_minv = std::min(SubGenJetAK8_pt[sgjIndex[0]], SubGenJetAK8_pt[sgjIndex[1]]);
            pt_max = std::max(SubGenJetAK8_pt[sgjIndex[0]], SubGenJetAK8_pt[sgjIndex[1]]);
            h_zg_ratio->Fill(pt_max/pt_minv);
            h_zg_subjets->Fill(pt_minv/(SubGenJetAK8_pt[sgjIndex[0]] + SubGenJetAK8_pt[sgjIndex[1]]));
            double E1, E2;
            E1=std::sqrt(SubGenJetAK8_mass[sgjIndex[0]]*SubGenJetAK8_mass[sgjIndex[0]] + SubGenJetAK8_pt[sgjIndex[0]]*SubGenJetAK8_pt[sgjIndex[0]]);
            E2=std::sqrt(SubGenJetAK8_mass[sgjIndex[1]]*SubGenJetAK8_mass[sgjIndex[1]] + SubGenJetAK8_pt[sgjIndex[1]]*SubGenJetAK8_pt[sgjIndex[1]]);
            h_ptheta_subjets->Fill(std::abs(E1 - E2)/GenJetAK8_pt[matchedJetIndexLL]);
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
            double pt_minp, pt_max;

            pt_minp = std::min(SubGenJetAK8_pt[sgjDefIndex.first], SubGenJetAK8_pt[sgjDefIndex.second]);
            pt_max = std::max(SubGenJetAK8_pt[sgjDefIndex.first], SubGenJetAK8_pt[sgjDefIndex.second]);
            h_zg_ratio->Fill(pt_max/pt_minp);
            h_zg_subjets->Fill(pt_minp/(SubGenJetAK8_pt[sgjDefIndex.first] + SubGenJetAK8_pt[sgjDefIndex.second]));
            double E1, E2;
            E1=std::sqrt(SubGenJetAK8_mass[sgjDefIndex.first]*SubGenJetAK8_mass[sgjDefIndex.first] + SubGenJetAK8_pt[sgjDefIndex.first]*SubGenJetAK8_pt[sgjDefIndex.first]);
            E2=std::sqrt(SubGenJetAK8_mass[sgjDefIndex.second]*SubGenJetAK8_mass[sgjDefIndex.second] + SubGenJetAK8_pt[sgjDefIndex.second]*SubGenJetAK8_pt[sgjDefIndex.second]);
            h_ptheta_subjets->Fill(std::abs(E1 - E2)/GenJetAK8_pt[matchedJetIndexLL]);
                
        }
        
        /*
        double pt_min_q = std::min(GenPart_pt[quarks_fromW.first], GenPart_pt[quarks_fromW.second]);
        h_zg_quarks->Fill(pt_min_q/(GenPart_pt[quarks_fromW.first] + GenPart_pt[quarks_fromW.second]));

        double E1_q, E2_q;
        E1_q = std::sqrt(GenPart_mass[quarks_fromW.first]*GenPart_mass[quarks_fromW.first] + GenPart_pt[quarks_fromW.first]*GenPart_pt[quarks_fromW.first]);
        E2_q = std::sqrt(GenPart_mass[quarks_fromW.second]*GenPart_mass[quarks_fromW.second] + GenPart_pt[quarks_fromW.second]*GenPart_pt[quarks_fromW.second]);
        h_ptheta_quarks->Fill(std::abs(E1_q - E2_q)/ GenPart_pt[hadWboson]);*/
    }

}
