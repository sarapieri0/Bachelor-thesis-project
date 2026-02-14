
void Makecanva(TCanvas* c, const std::vector<TH1*>& hists, Color_t color){
    
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    Int_t nHist = hists.size();
    for(TH1 *h : hists){
        if (!h) {
            Error("MyCode", "Istogramma nullo");
            return;
        }
    }

    if(nHist == 1){
        hists.front()->Scale(1.0 / hists.front()->Integral());
        hists.front()->Draw("hist");
        hists.front()->SetLineColor(color);
        hists.front()->SetFillColor(color);
        hists.front()->GetXaxis()->SetTitleSize(0.05);
    }
    else if(nHist == 2){
        hists.front()->Scale(1.0 / hists.front()->Integral());
        hists.back()->Scale(1.0 / hists.back()->Integral());

        double max1  = hists.front()->GetMaximum();
        double max2 = hists.back()->GetMaximum();

        double ymax = std::max(max1, max2) * 1.2;
        hists.front()->SetMaximum(ymax);

        gStyle->SetOptStat(0);
        hists.front()->Draw("hist");
        hists.front()->SetLineColor(color);
        hists.front()->SetFillColor(color);
        hists.front()->GetXaxis()->SetTitleSize(0.05);

        hists.back()->SetFillStyle(1001);
        hists.back()->SetFillColorAlpha(kBlack, 0.15);
        hists.back()->SetLineColor(kBlack);
        hists.back()->SetMarkerStyle(20);
        hists.back()->SetMarkerSize(1);
        hists.back()->SetMarkerColor(kBlack);
        hists.back()->Draw("hist E same");
        hists.back()->GetXaxis()->SetTitleSize(0.05);
    }
    else{
        return;
    }
    
}

void MakecanvaColor(TCanvas* c, const std::vector<TH1*>& hists, const std::vector<Color_t>& colors){
    for(TH1 *h : hists){
        if (!h) {
            Error("MyCode", "Istogramma nullo");
            return;
        }
    }
    hists.front()->SetFillStyle(0);      
    hists.front()->SetFillColor(0);
    hists.back()->SetFillStyle(0);      
    hists.back()->SetFillColor(0); 

    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);    
    hists.front()->Scale(1.0 / hists.front()->Integral());
    hists.back()->Scale(1.0 / hists.back()->Integral());
    double max1_1 = hists.front()->GetMaximum();
    double max2_1 = hists.back()->GetMaximum();
    double ymax_1 = std::max(max1_1, max2_1) * 1.2;
    hists.front()->SetMaximum(ymax_1);

    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);   
    hists.front()->Draw("hist");
    hists.front()->SetLineColor(colors[0]);
    hists.front()->SetFillStyle(1001);
    hists.front()->SetFillColorAlpha(colors[0], 0.1);
    //hists.front()->SetMarkerSize(1.2);
    //hists.front()->SetMarkerStyle(22);
    //hists.front()->SetMarkerColor(colors[0]);
    hists.front()->GetXaxis()->SetTitleSize(0.05);
    //hists.back()->SetFillStyle(1001);
    hists.back()->Draw("hist E same");
    hists.back()->SetLineColor(colors[1]);
    //hists.back()->SetFillColorAlpha(colors[1], 0.1);
    hists.back()->SetMarkerSize(1.2);
    hists.back()->SetMarkerStyle(22);
    hists.back()->SetMarkerColor(colors[1]);
    hists.back()->GetXaxis()->SetTitleSize(0.05);

    
    
    
}






void Makecanva2x1(TCanvas* c, const std::vector<TH1*>& hists, Color_t color){
    for(TH1 *h : hists){
        if (!h) {
            Error("MyCode", "Istogramma nullo");
            return;
        }
    }
    
    c->Divide(2,1);
    if(hists.size() == 2){
        int i=1;
        for(TH1* h : hists){
            if(!h) continue;
            c->cd(i);
            gPad->SetLeftMargin(0.15);
            gPad->SetBottomMargin(0.15);    
            h->Scale(1.0 / h->Integral());
            h->Draw("hist");
            h->SetLineColor(color);
            h->SetFillColor(color);
            h->GetXaxis()->SetTitleSize(0.05);
            i++;
        }
    }
    else if(hists.size() == 4){
        c->cd(1);
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.15);   
        hists[0]->Scale(1.0 / hists[0]->Integral());
        hists[0]->Draw("hist E");
        hists[0]->SetLineColor(kBlack);
        hists[0]->SetMarkerSize(1.5);
        hists[0]->SetMarkerStyle(22);
        hists[0]->SetMarkerColor(kBlack);
        hists[0]->GetXaxis()->SetTitleSize(0.05);
        hists[1]->Scale(1.0 / hists[1]->Integral());
        hists[1]->SetFillStyle(1001);
        hists[1]->Draw("hist E same");
        hists[1]->SetLineColor(color);
        hists[1]->SetFillColorAlpha(color, 0.1);
        hists[1]->SetMarkerSize(1.2);
        hists[1]->SetMarkerStyle(20);
        hists[1]->SetMarkerColor(color);
        hists[1]->GetXaxis()->SetTitleSize(0.05);


        c->cd(2);
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.15);   
        hists[2]->Scale(1.0 / hists[2]->Integral());
        hists[2]->Draw("hist E");
        hists[2]->SetLineColor(kBlack);
        hists[2]->SetMarkerSize(1.5);
        hists[2]->SetMarkerStyle(22);
        hists[2]->SetMarkerColor(kBlack);
        //hists[0]->SetFillColor(color);
        hists[2]->GetXaxis()->SetTitleSize(0.05);
        hists[3]->Scale(1.0 / hists[3]->Integral());
        hists[3]->SetFillStyle(1001);
        hists[3]->Draw("hist E same");
        hists[3]->SetLineColor(color);
        hists[3]->SetFillColorAlpha(color, 0.1);
        hists[3]->SetMarkerSize(1.2);
        hists[3]->SetMarkerStyle(20);
        hists[3]->SetMarkerColor(color);
        hists[3]->GetXaxis()->SetTitleSize(0.05);
    }
    
    
    else{  return;  }
    
}

void Makecanva2x1Bicolor(TCanvas* c, const std::vector<TH1*>& hists, const std::vector<Color_t>& colors){
    for(TH1 *h : hists){
        if (!h) {
            Error("MyCode", "Istogramma nullo");
            return;
        }
    }

    c->Divide(2,1);
    if(hists.size()==2){
        int i=1;
        for(TH1* h : hists){
            if(!h) continue;
            c->cd(i);
            gPad->SetLeftMargin(0.15);
            gPad->SetBottomMargin(0.15);    
            h->Scale(1.0 / h->Integral());
            h->SetFillStyle(1001);
            h->Draw("hist E same");
            h->SetLineColor(colors[i-1]);
            h->SetFillColorAlpha(colors[i-1], 0.1);
            h->SetMarkerSize(1);
            h->SetMarkerStyle(22);
            h->SetMarkerColor(colors[i-1]);
            h->GetXaxis()->SetTitleSize(0.05);
            i++;
        }
    }
    else if(hists.size()==4){
        c->cd(1);
        hists[0]->Scale(1.0 / hists[0]->Integral());
        hists[1]->Scale(1.0 / hists[1]->Integral());
        double max1_1 = hists[0]->GetMaximum();
        double max2_1 = hists[1]->GetMaximum();
        double ymax_1 = std::max(max1_1, max2_1) * 1.2;
        hists[0]->SetMaximum(ymax_1);

        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.15);   
        hists[0]->Draw("hist");
        hists[0]->SetLineColor(colors[0]);
        hists[0]->SetFillColorAlpha(colors[0], 0.1);
        //hists[0]->SetMarkerSize(1.2);
        //hists[0]->SetMarkerStyle(22);
        //hists[0]->SetMarkerColor(colors[0]);
        hists[0]->GetXaxis()->SetTitleSize(0.05);
        hists[1]->SetFillStyle(1001);
        hists[1]->Draw("hist E same");
        hists[1]->SetLineColor(colors[1]);
        //hists[1]->SetFillColorAlpha(kBlack, 0.1);
        hists[1]->SetMarkerSize(1);
        hists[1]->SetMarkerStyle(20);
        hists[1]->SetMarkerColor(colors[1]);
        hists[1]->GetXaxis()->SetTitleSize(0.05);


        c->cd(2);
        hists[2]->Scale(1.0 / hists[2]->Integral());
        hists[3]->Scale(1.0 / hists[3]->Integral());
        double max1_2 = hists[2]->GetMaximum();
        double max2_2 = hists[3]->GetMaximum();
        double ymax_2 = std::max(max1_2, max2_2) * 1.2;
        hists[2]->SetMaximum(ymax_2);

        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.15);   
        hists[2]->Draw("hist");
        hists[2]->SetLineColor(colors[0]);
        hists[2]->SetFillColorAlpha(colors[0], 0.1);
        //hists[2]->SetMarkerSize(1.2);
        //hists[2]->SetMarkerStyle(22);
        //hists[2]->SetMarkerColor(colors[0]);
        hists[2]->GetXaxis()->SetTitleSize(0.05);
        hists[3]->SetFillStyle(1001);
        hists[3]->Draw("hist E same");
        hists[3]->SetLineColor(colors[1]);
        //hists[3]->SetFillColorAlpha(kBlack, 0.1);
        hists[3]->SetMarkerSize(1);
        hists[3]->SetMarkerStyle(20);
        hists[3]->SetMarkerColor(colors[1]);
        hists[3]->GetXaxis()->SetTitleSize(0.05);

    }
    
    
    
    
}


void Makecanva2x2(TCanvas* c, const std::vector<TH1*>& hists, Color_t color){
    for(TH1 *h : hists){
        if (!h) {
            Error("MyCode", "Istogramma nullo");
            return;
        }
    }

    c->Divide(2,2);
    if(hists.size()==4){
        int i=1;
        for(TH1* h : hists){
            if(!h) continue;
            c->cd(i);
            gPad->SetLeftMargin(0.15);
            gPad->SetBottomMargin(0.15);  
            h->Scale(1.0 / h->Integral());
            h->Draw("hist");
            h->SetLineColor(color);
            h->SetFillColor(color);
            h->GetXaxis()->SetTitleSize(0.05);
            i++;
        }
    }

    else if(hists.size()==8){
        for(TH1 *h : hists){
        if (!h) {
            Error("MyCode", "Istogramma nullo");
            return;
        }
        }
        hists[0]->Scale(1.0 / hists[0]->Integral());
        hists[1]->Scale(1.0 / hists[1]->Integral());
        double max1_1 = hists[0]->GetMaximum();
        double max2_1 = hists[1]->GetMaximum();
        double ymax_1 = std::max(max1_1, max2_1) * 1.2;
        hists[0]->SetMaximum(ymax_1);

        c->cd(1);
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.15);   
        hists[0]->Draw("hist");
        hists[0]->SetLineColor(color);
        hists[0]->SetFillColor(color);
        hists[0]->GetXaxis()->SetTitleSize(0.05);
        hists[1]->SetFillStyle(1001);
        hists[1]->Draw("hist E same");
        hists[1]->SetLineColor(kBlack);
        hists[1]->SetFillColorAlpha(kBlack, 0.15);
        hists[1]->SetMarkerSize(1.2);
        hists[1]->SetMarkerStyle(20);
        hists[1]->GetXaxis()->SetTitleSize(0.05);


        hists[2]->Scale(1.0 / hists[2]->Integral());
        hists[3]->Scale(1.0 / hists[3]->Integral());
        double max1_2 = hists[2]->GetMaximum();
        double max2_2 = hists[3]->GetMaximum();
        double ymax_2 = std::max(max1_2, max2_2) * 1.2;
        hists[2]->SetMaximum(ymax_2);

        c->cd(2);
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.15);   
        hists[2]->Draw("hist");
        hists[2]->SetLineColor(color);
        hists[2]->SetFillColor(color);
        hists[2]->GetXaxis()->SetTitleSize(0.05);
        hists[3]->SetFillStyle(1001);
        hists[3]->Draw("hist E same");
        hists[3]->SetLineColor(kBlack);
        hists[3]->SetFillColorAlpha(kBlack, 0.15);
        hists[3]->SetMarkerSize(1.2);
        hists[3]->SetMarkerStyle(20);
        hists[3]->GetXaxis()->SetTitleSize(0.05);


        hists[4]->Scale(1.0 / hists[4]->Integral());
        hists[5]->Scale(1.0 / hists[5]->Integral());
        double max1_3 = hists[4]->GetMaximum();
        double max2_3 = hists[5]->GetMaximum();
        double ymax_3 = std::max(max1_3, max2_3) * 1.2;
        hists[4]->SetMaximum(ymax_3);

        c->cd(3);
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.15);   
        hists[4]->Draw("hist");
        hists[4]->SetLineColor(color);
        hists[4]->SetFillColor(color);
        hists[4]->GetXaxis()->SetTitleSize(0.05);
        hists[5]->SetFillStyle(1001);
        hists[5]->Draw("hist E same");
        hists[5]->SetLineColor(kBlack);
        hists[5]->SetFillColorAlpha(kBlack, 0.15);
        hists[5]->SetMarkerSize(1.2);
        hists[5]->SetMarkerStyle(20);
        hists[5]->GetXaxis()->SetTitleSize(0.05);


        hists[6]->Scale(1.0 / hists[6]->Integral());
        hists[7]->Scale(1.0 / hists[7]->Integral());
        double max1_4 = hists[6]->GetMaximum();
        double max2_4 = hists[7]->GetMaximum();
        double ymax_4 = std::max(max1_4, max2_4) * 1.2;
        hists[6]->SetMaximum(ymax_4);

        c->cd(4);
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.15);   
        hists[6]->Draw("hist");
        hists[6]->SetLineColor(color);
        hists[6]->SetFillColor(color);
        hists[6]->GetXaxis()->SetTitleSize(0.05);
        hists[7]->SetFillStyle(1001);
        hists[7]->Draw("hist E same");
        hists[7]->SetLineColor(kBlack);
        hists[7]->SetFillColorAlpha(kBlack, 0.15);
        hists[7]->SetMarkerSize(1.2);
        hists[7]->SetMarkerStyle(20);
        //hists[7]->Draw("P0 same");
        hists[7]->GetXaxis()->SetTitleSize(0.05);
    }


    else return;
}

void Makecanva2x2Bicolor(TCanvas* c, const std::vector<TH1*>& hists, const std::vector<Color_t> colors){

    c->Divide(2,2);
    c->cd(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);   
    hists[0]->Scale(1.0 / hists[0]->Integral());
    hists[0]->Draw("hist");
    hists[0]->SetLineColor(colors[0]);
    hists[0]->SetFillColor(colors[0]);
    hists[0]->GetXaxis()->SetTitleSize(0.05);

    c->cd(2);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);   
    hists[1]->Scale(1.0 / hists[1]->Integral());
    hists[1]->Draw("hist");
    hists[1]->SetLineColor(colors[1]);
    hists[1]->SetFillColor(colors[1]);
    hists[1]->GetXaxis()->SetTitleSize(0.05);

    c->cd(3);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);   
    hists[2]->Scale(1.0 / hists[2]->Integral());
    hists[2]->Draw("hist");
    hists[2]->SetLineColor(colors[0]);
    hists[2]->SetFillColor(colors[0]);
    hists[2]->GetXaxis()->SetTitleSize(0.05);

    c->cd(4);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);   
    hists[3]->Scale(1.0 / hists[3]->Integral());
    hists[3]->Draw("hist");
    hists[3]->SetLineColor(colors[1]);
    hists[3]->SetFillColor(colors[1]);
    hists[3]->GetXaxis()->SetTitleSize(0.05);

}







void Makecanva3x1(TCanvas* c, const std::vector<TH1*>& hists, Color_t color){
    c->Divide(3,1);


    if(hists.size() == 3){
        int i=1;
        for(TH1* h : hists){
            if(!h) continue;
            c->cd(i);
            gPad->SetLeftMargin(0.15);
            gPad->SetBottomMargin(0.15);    
            h->Scale(1.0 / h->Integral());
            h->Draw("hist");
            h->SetLineColor(color);
            h->SetFillColor(color);
            h->GetXaxis()->SetTitleSize(0.05);
            i++;
        }
    }

    else if(hists.size() == 6){
        c->cd(1);

        hists[0]->Scale(1.0 / hists[0]->Integral());
        hists[1]->Scale(1.0 / hists[1]->Integral());
        double max1_1 = hists[0]->GetMaximum();
        double max2_1 = hists[1]->GetMaximum();
        double ymax_1 = std::max(max1_1, max2_1) * 1.2;
        hists[0]->SetMaximum(ymax_1);

        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.15);   
        hists[0]->Draw("hist");
        hists[0]->SetLineColor(color);
        hists[0]->SetFillColor(color);
        hists[0]->GetXaxis()->SetTitleSize(0.05);
        hists[1]->SetFillStyle(1001);
        hists[1]->Draw("hist E same");
        hists[1]->SetLineColor(kBlack);
        hists[1]->SetFillColorAlpha(kBlack, 0.15);
        hists[1]->SetMarkerSize(1.2);
        hists[1]->SetMarkerStyle(20);
        hists[1]->GetXaxis()->SetTitleSize(0.05);

        c->cd(2);
        hists[2]->Scale(1.0 / hists[2]->Integral());
        hists[3]->Scale(1.0 / hists[3]->Integral());
        double max1_2 = hists[2]->GetMaximum();
        double max2_2 = hists[3]->GetMaximum();
        double ymax_2 = std::max(max1_2, max2_2) * 1.2;
        hists[2]->SetMaximum(ymax_2);

        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.15);   
        hists[2]->Draw("hist");
        hists[2]->SetLineColor(color);
        hists[2]->SetFillColor(color);
        hists[2]->GetXaxis()->SetTitleSize(0.05);
        hists[3]->SetFillStyle(1001);
        hists[3]->Draw("hist E same");
        hists[3]->SetLineColor(kBlack);
        hists[3]->SetFillColorAlpha(kBlack, 0.15);
        hists[3]->SetMarkerSize(1.2);
        hists[3]->SetMarkerStyle(20);
        hists[3]->GetXaxis()->SetTitleSize(0.05);

        c->cd(3);
        hists[4]->Scale(1.0 / hists[4]->Integral());
        hists[5]->Scale(1.0 / hists[5]->Integral());
        double max1_3 = hists[4]->GetMaximum();
        double max2_3 = hists[5]->GetMaximum();
        double ymax_3 = std::max(max1_3, max2_3) * 1.2;
        hists[4]->SetMaximum(ymax_3);

        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.15);   
        hists[4]->Draw("hist");
        hists[4]->SetLineColor(color);
        hists[4]->SetFillColor(color);
        hists[4]->GetXaxis()->SetTitleSize(0.05);
        hists[5]->SetFillStyle(1001);
        hists[5]->Draw("hist E same");
        hists[5]->SetLineColor(kBlack);
        hists[5]->SetFillColorAlpha(kBlack, 0.15);
        hists[5]->SetMarkerSize(1.2);
        hists[5]->SetMarkerStyle(20);
        hists[5]->GetXaxis()->SetTitleSize(0.05);
    }
}
