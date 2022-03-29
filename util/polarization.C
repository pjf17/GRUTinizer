#include <vector>
#include <string>

#include "TFile.h"
#include "TH1D.h"

void polarization(TFile *f, std::string folder, int bw, bool swapped=false){
    std::string hpath = folder;
    if (swapped) hpath += "/gamma_corrected_swapped_addback_prompt";
    else hpath += "/gamma_corrected_addback_prompt";
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

    //add and subtract
    std::vector<TH1D*> sum_hists, diff_hists, stat_hists;
    for (int i=0; i < 3; i++){
        for (int j=i+1; j < 3; j++){
            double norm = 1.0;
            if (i == 0){
                norm = nRed*1.0/nGB;
            }
            TH1D *hSum = (TH1D*) hists[i]->Clone(Form("%s+%s",hists[i]->GetName(),hists[j]->GetName()));
            TH1D *hDif = (TH1D*) hists[i]->Clone(Form("%s-%s",hists[i]->GetName(),hists[j]->GetName()));
            TH1D *hStat = (TH1D*) hists[i]->Clone(Form("%s-%s/(sigma)",hists[i]->GetName(),hists[j]->GetName()));
            hSum->Add(hists[j],norm);
            hDif->Add(hists[j],-1.0*norm);

            //calc what the expected statistical fluctuation of subtracting the
            //hist would be if it was all random
            int nBins = hists[i]->GetNbinsX();
            for (int bin=1; bin <= nBins; bin++){
                double difOverErr = 0.0;
                if (hists[i]->GetBinContent(bin) > 0 || hists[j]->GetBinContent(bin) > 0) 
                    difOverErr = hDif->GetBinContent(bin)/(TMath::Sqrt(hists[i]->GetBinContent(bin)) + TMath::Sqrt(hists[j]->GetBinContent(bin)*norm));

                hStat->SetBinContent(bin,difOverErr);
            }

            sum_hists.push_back(hSum);
            diff_hists.push_back(hDif);
            stat_hists.push_back(hStat);
        }
    }
    std::string fname = Form("pol_%s.root",folder.c_str());
    if (swapped) fname = Form("pol_swapped_%s.root",folder.c_str());
    TFile *fOut = new TFile(fname.c_str(),"RECREATE");
    fOut->cd();
    for (int i=0; i < 3; i++){
        hists[i]->Write();
        sum_hists[i]->Write();
        diff_hists[i]->Write();
        stat_hists[i]->Write();
    }
    return;
}