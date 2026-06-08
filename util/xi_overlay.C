#include "MultiPlotter.h"

// MultiPlotter *mp = new MultiPlotter();
// 
// void xi_overlay(TFile *f1, TFile *f2){
    // 
    // mp->Add(f1,"Bkg Subtracted");
    // mp->Add(f2,"Bkg Subtracted");
// 
    // mp->SetLegendEntry("Bkg Subtracted","#tilde{#xi} < 30#circ");
    // mp->SetLegendEntry("Bkg Subtracted_perp","#tilde{#xi} > 60#circ");
// 
    // mp->SetLineColor("Bkg Subtracted_perp",kBlue);
    // mp->SetLineColor("Bkg Subtracted",kRed);
    // mp->SetFill(0.4);
    // mp->SetLegendCoordinates(0.78,0.84,0.98,0.98);
    // mp->Draw();
// }

std::string nuclideToLatex(const std::string& nuclide) {
    // Find where the letters begin
    auto it = std::find_if(nuclide.begin(), nuclide.end(), ::isalpha);

    std::string mass(nuclide.begin(), it);   // e.g. "48"
    std::string symbol(it, nuclide.end());   // e.g. "Ca"

    return "^{" + mass + "}" + symbol;
}

void do_format(TH1* hh){
    hh->GetXaxis()->SetTitle("Energy (keV)");
    if (hh->GetBinWidth(1) > 1)
        hh->GetYaxis()->SetTitle(Form("Counts / %d keV",int(hh->GetBinWidth(1))));
    else 
        hh->GetYaxis()->SetTitle("Counts / keV");

    hh->GetXaxis()->CenterTitle();
    hh->GetYaxis()->CenterTitle();

    hh->GetXaxis()->SetLabelFont(43);
    hh->GetYaxis()->SetLabelFont(43);
    
    hh->GetXaxis()->SetLabelSize(23);
    hh->GetYaxis()->SetLabelSize(23);

    hh->GetXaxis()->SetTitleFont(43);
    hh->GetYaxis()->SetTitleFont(43);

    hh->GetXaxis()->SetTitleSize(25);
    hh->GetYaxis()->SetTitleSize(25);

    hh->GetMaximum() > 999 ? hh->GetYaxis()->SetTitleOffset(1.05) : hh->GetYaxis()->SetTitleOffset(0.85);
    hh->GetXaxis()->SetTitleOffset(0.9);

    hh->GetXaxis()->SetNdivisions(705);
}

void xi_overlay(TFile *f1, TFile *f2, std::string species="", std::string beta=""){
    TH1* para = (TH1*) f1->Get("Bkg Subtracted");
    TH1* perp = (TH1*) f2->Get("Bkg Subtracted");
    TLegend *leg = new TLegend(0.785,0.77,0.97,0.98);
    leg->SetTextAlign(22);
    leg->AddEntry(perp,"#tilde{#xi} > 60#circ ","l");
    leg->AddEntry(para,"#tilde{#xi} < 30#circ ","l");

    para->SetLineColor(kRed);
    para->SetFillColorAlpha(kRed,0.4);

    perp->SetLineColor(kBlue);
    perp->SetFillColorAlpha(kBlue,0.4);

    do_format(perp);
    do_format(para);

    TCanvas *canv = new TCanvas("polarization","polarization",600,400);
    // TCanvas *canv = new TCanvas("polarization","polarization",500,450);

    perp->Draw("hist");
    para->Draw("hist same");

    gPad->SetBottomMargin(0.13);
    std::max(para->GetMaximum(),perp->GetMaximum()) > 999 ? gPad->SetLeftMargin(0.13) : gPad->SetLeftMargin(0.11);
    gPad->SetRightMargin(0.03);
    gPad->SetTopMargin(0.02);

    if (species != ""){
        TLatex nuclide;
        nuclide.SetTextAlign(32);
        nuclide.SetTextFont(43);
        nuclide.SetTextSize(25);
        nuclide.DrawLatexNDC(0.96,0.71,nuclideToLatex(species).c_str());
        nuclide.DrawLatexNDC(0.96,0.64,Form("v/c=%s",beta.c_str()));
    }

    canv->Update();

    leg->Draw("same");
}

void drawSubFig(std::string str){
    TText text;
    text.SetTextAlign(12);
    text.SetTextFont(43);
    text.SetTextSize(25);
    double xpos = gPad->GetLeftMargin() > 0.11 ? 0.18 : 0.16;
    text.DrawTextNDC(xpos,0.93,str.c_str());
}

void readFile(std::string filename, std::map<int,std::pair<double,double>> &gmap){
    std::ifstream file(filename);
    std::string line;

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        int first;
        double glo, ghi;
        std::string temp;

        // Read first number
        ss >> first;

        // Skip 7 fields
        for (int i = 0; i < 7; ++i) {
            ss >> temp;
        }

        // Read second-to-last and last number
        ss >> glo >> ghi;

        gmap[first] = std::make_pair(glo,ghi);
    }
}

std::vector<int> get_map_keys(const std::map<int,std::pair<double,double>> & my_map) {
    std::vector<int> keys;
    for (const auto& pair : my_map) {
        keys.push_back(pair.first);
    }
    return keys;
}

void draw_integ_reg(std::string fpara, std::string fperp, int key){
    std::vector<int> keys; keys.push_back(key);
    //get the scale of the histogram
    TH1D *h;
    double ymax = 0;
    TIter iter(gPad->GetListOfPrimitives());
    while(TObject *obj=iter.Next()) {
        if(obj->InheritsFrom(TH1::Class())) {
            h = (TH1D*)obj;
            double tempy = h->GetMaximum();
            if (tempy > ymax) ymax = tempy;
        }
    }

    //read the file
    std::map<int,std::pair<double,double>> para, perp;
    readFile(fpara,para);
    readFile(fperp,perp);

    //draw all peaks if none specified
    if (keys.size() == 0) keys = get_map_keys(para);

    double xwidth = 0;

    //draw the gates
    int Np = keys.size();
    for (int i = 0; i < Np; i++){
        int peak = keys[i];
        std::vector<TLine *> lpara = {new TLine(para[peak].first,0,para[peak].first,ymax), new TLine(para[peak].second,0,para[peak].second,ymax)};
        std::vector<TLine *> lperp = {new TLine(perp[peak].first,0,perp[peak].first,ymax), new TLine(perp[peak].second,0,perp[peak].second,ymax)};
        xwidth = para[peak].second - para[peak].first;
        
        for (int j=0; j < 2; j++){
            lpara[j]->SetLineColor(kRed);
            lperp[j]->SetLineColor(kBlue);
            lpara[j]->SetLineStyle(kDashed);
            lperp[j]->SetLineStyle(kDashed);
            lperp[j]->SetLineWidth(3);
            lpara[j]->SetLineWidth(3);
            lpara[j]->Draw("same");
            lperp[j]->Draw("same");
        }
    }

    TLatex eng;
    eng.SetTextAlign(21);
    eng.SetTextFont(43);
    eng.SetTextSize(25);
    eng.SetTextAngle(90);
    eng.DrawLatex(key-0.1*xwidth,0.9*ymax,Form("%d",key));
    printf("%f %f\n",xwidth,key-0.05*xwidth);

    return;
}

void SetAxisrange(double xlo, double xhi, bool doX){
    if (xlo > xhi) std::swap(xlo,xhi);
    //get the histograms on the canvas
    TH1 *hh;
    TIter iter(gPad->GetListOfPrimitives());
    while(TObject *obj=iter.Next()) {
        if(obj->InheritsFrom(TH1::Class())) {
            hh = (TH1*)obj;
            if (doX) hh->GetXaxis()->SetRangeUser(xlo,xhi);
            else hh->GetYaxis()->SetRangeUser(xlo,xhi);
        }
    }
    gPad->Modified();
    gPad->Update();
    return;
}

void SetXrange(double xlo, double xhi){
    SetAxisrange(xlo,xhi,true);
}

void SetYrange(double xlo, double xhi){
    SetAxisrange(xlo,xhi,false);
}