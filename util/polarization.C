#include <vector>
#include <string>

#include "TFile.h"
#include "TH1D.h"

int nconfigs = 1;

void polarization(TFile *f, std::string folder, int bw){
    std::vector< std::vector<TH1D*>> hists;
    for (int i=0; i < nconfigs; i++){
        std::string hpath = folder;
        if (i == 1) hpath += "/gamma_corrected_swapped_addback";
        else hpath += "/gamma_corrected_addback";
        std::vector<TH1D*> temp_hists;
        temp_hists.push_back( (TH1D*) f->Get(Form("%s_%s",hpath.c_str(),"red_pair")));
        temp_hists.back()->SetNameTitle("red","red");
        temp_hists.push_back( (TH1D*) f->Get(Form("%s_%s",hpath.c_str(),"blue_pair")));
        temp_hists.back()->SetNameTitle("blue","blue");
        temp_hists.push_back( (TH1D*) f->Get(Form("%s_%s",hpath.c_str(),"gold_pair")));
        temp_hists.back()->SetNameTitle("gold","gold");
        hists.push_back(temp_hists);
    }
    
    //rebin to specified bins
    for (int i=0; i < nconfigs; i++) {
        for (int j=0; j < 3; j++) {
            hists[i][j]->Rebin(bw);
        }
    }

    //normalize
    // double slope[2] = {-8.8686E-6,-3.2959E-06};
    // double incpt[2] = {1.1699,1.1843};
    double norm[2] = {1.16267,1.18164};
    // double norm[2] = {18.0/14,18.0/13};
 
    //add and subtract
    std::vector<std::vector<TH1D*>> sum_hists, diff_hists, stat_hists;
    for (int i=0; i < nconfigs; i++){
        std::vector<TH1D*> temp_sum_hists, temp_diff_hists, temp_stat_hists;
        for (int j=0; j < 3; j++){
            for (int k=j+1; k < 3; k++){
                // double norm = 1.0;
                
                std::string flag = "dflt";
                if (i == 1) flag = "swap"; 
                TH1D *hSum = (TH1D*) hists[i][j]->Clone(Form("%s+%s_%s",hists[i][j]->GetName(),hists[i][k]->GetName(),flag.c_str()));
                TH1D *hDif = (TH1D*) hists[i][j]->Clone(Form("%s-%s_%s",hists[i][j]->GetName(),hists[i][k]->GetName(),flag.c_str()));
                TH1D *hStat = (TH1D*) hists[i][j]->Clone(Form("%s-%s/(sigma)_%s",hists[i][j]->GetName(),hists[i][k]->GetName(),flag.c_str()));
                if (j!=0){
                    hSum->Scale(norm[0]);
                    hDif->Scale(norm[0]);
                }
                hSum->Add(hists[i][k],norm[k-1]);
                hDif->Add(hists[i][k],-1.0*norm[k-1]);

                //calc what the expected statistical fluctuation of subtracting the
                //hist would be if it was all random
                // int nBins = hists[i][j]->GetNbinsX();
                // for (int bin=1; bin <= nBins; bin++){
                //     double binContent = 0.0;
                //     double norm = hists[i][k]->GetBinCenter(bin)*slope[k-1] + incpt[k-1];

                //     binContent = hSum->GetBinContent(bin) + norm*hists[i][k]->GetBinContent(bin);
                //     hSum->SetBinContent(bin,binContent);

                //     binContent = hDif->GetBinContent(bin) - norm*hists[i][k]->GetBinContent(bin);
                //     hDif->SetBinContent(bin,binContent);

                //     binContent = 0.0;
                //     if (hists[i][j]->GetBinContent(bin) > 0 || hists[i][k]->GetBinContent(bin) > 0) 
                //         binContent = hDif->GetBinContent(bin)/TMath::Sqrt(hists[i][j]->GetBinContent(bin) + hists[i][k]->GetBinContent(bin)*norm);

                //     hStat->SetBinContent(bin,binContent);
                // }

                temp_sum_hists.push_back(hSum);
                temp_diff_hists.push_back(hDif);
                temp_stat_hists.push_back(hStat);
            }
        }
        sum_hists.push_back(temp_sum_hists);
        diff_hists.push_back(temp_diff_hists);
        stat_hists.push_back(temp_stat_hists);
    }
    //draw stuff
    for (int j=0; j < 3; j++){
        TCanvas *canv = new TCanvas(diff_hists[0][j]->GetName(),diff_hists[0][j]->GetName());
        canv->Divide(2,2,0.001,0.001);
        // sum_hists[1][j]->GetXaxis()->SetRangeUser(0,800);
        // diff_hists[1][j]->GetXaxis()->SetRangeUser(0,800);
        // stat_hists[1][j]->GetXaxis()->SetRangeUser(0,800);
        sum_hists[0][j]->GetXaxis()->SetRangeUser(0,1600);
        diff_hists[0][j]->GetXaxis()->SetRangeUser(0,1600);
        stat_hists[0][j]->GetXaxis()->SetRangeUser(0,1600);
        for (int i=0; i < nconfigs; i++){
            // canv->cd(2-i);
            // stat_hists[i][j]->Draw();
            canv->cd(2-i);
            diff_hists[i][j]->Draw();
            canv->cd(4-i);
            sum_hists[i][j]->Draw();
        }
        // canv->Write();
    }

    //write out the files
    std::string fname = Form("pol_%s.root",folder.c_str());
    TFile *fOut = new TFile(fname.c_str(),"RECREATE");
    fOut->cd();
    for (int i=0; i < nconfigs; i++){
        for (int j=0; j < 3; j++){
            hists[i][j]->Write();
            sum_hists[i][j]->Write();
            diff_hists[i][j]->Write();
            stat_hists[i][j]->Write();
        }
    }
    return;
}