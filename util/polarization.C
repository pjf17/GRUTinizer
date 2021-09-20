#include <vector>
#include <string>

#include "TFile.h"
#include "TH1D.h"

void polarization(TFile *f, std::string folder, int bw){
    std::string hpath = folder + "/gamma_corrected_addback_prompt";
    std::vector<TH1D*> hists;
    hists.push_back( (TH1D*) f->Get(Form("%s_%s",hpath.c_str(),"red_pair")));
    hists.back()->SetNameTitle("red","red");
    hists.push_back( (TH1D*) f->Get(Form("%s_%s",hpath.c_str(),"blue_pair")));
    hists.back()->SetNameTitle("blue","blue");
    hists.push_back( (TH1D*) f->Get(Form("%s_%s",hpath.c_str(),"gold_pair")));
    hists.back()->SetNameTitle("gold","gold");
    
    //rebin to 16 keV bins
    for (int i=0; i < 3; i++) {
        hists[i]->Rebin(bw);
    }

    //normalize
    int nRed = 18;
    int nGB = 13;
    for (int i=1; i<3; i++) hists[i]->Scale(nRed*1.0/nGB);

    //add and subtract
    std::vector<TH1D*> sum_hists, diff_hists;
    for (int i=0; i < 3; i++){
        for (int j=i+1; j < 3; j++){
            TH1D *hSum = (TH1D*) hists[i]->Clone(Form("%s+%s",hists[i]->GetName(),hists[j]->GetName()));
            TH1D *hDif = (TH1D*) hists[i]->Clone(Form("%s-%s",hists[i]->GetName(),hists[j]->GetName()));
            hSum->Add(hists[j]);
            hDif->Add(hists[j],-1.0);
            sum_hists.push_back(hSum);
            diff_hists.push_back(hDif);
        }
    }

    TFile *fOut = new TFile(Form("pol_%s.root",folder.c_str()),"RECREATE");
    fOut->cd();
    for (int i=0; i < 3; i++){
        hists[i]->Write();
        sum_hists[i]->Write();
        diff_hists[i]->Write();
    }
    return;
}