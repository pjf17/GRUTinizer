#include <vector>
#include <string>

#include "TFile.h"
#include "TH1D.h"
// #include "../include/GH1D.h"

void plot(std::vector<TH1D*> &vec){
    TCanvas *canv = new TCanvas();
    bool first = true;
    for (auto &elem: vec){
        if (first) {
            elem->Draw();
            first = false;
        } else {
            elem->Draw("same");
        }
    }
    return;
}

void pair_check(TFile *f, std::string folder, int bw){
    std::vector<GH1D*> hBlue, hRed, hGold;
    int nbins = 100;
    double xlo = 0;
    double xhi = 0.2;
    TH1D *rRed = new TH1D("redres","Red Pair Resolutions",nbins,xlo,xhi);
    TH1D *rBlue = new TH1D("blueres","Blue Pair Resolutions",nbins,xlo,xhi);
    TH1D *rGold = new TH1D("goldres","Gold Pair Resolutions",nbins,xlo,xhi);
    rRed->SetLineColor(kRed);
    rBlue->SetLineColor(kBlue);
    rGold->SetLineColor(kOrange);

    f->cd(folder.c_str());
    TIter next(gDirectory->GetListOfKeys());
    TKey *key;
    while ((key = (TKey*)next())) {
        std::string hname = std::string(key->GetName());
        GH1D *h = (GH1D*) key->ReadObj();
        h->Rebin(bw);
        // h->GetXaxis()->SetRangeUser(0,600);
        // new TCanvas();
        GGaus *fit = GausFit(h,270,420,"0");
        double res = fit->GetFWHM()/fit->GetCentroid();

        if (hname.find("red_pair_") != std::string::npos){
            hRed.push_back(h);
            rRed->Fill(res);
        }
        if (hname.find("blue_pair_") != std::string::npos){
            hBlue.push_back(h);
            rBlue->Fill(res);
        }
        if (hname.find("gold_pair_") != std::string::npos){
            hGold.push_back(h);
            rGold->Fill(res);
        }
    }

    rRed->Draw();
    rBlue->Draw("same");
    rGold->Draw("same");
    return;
}