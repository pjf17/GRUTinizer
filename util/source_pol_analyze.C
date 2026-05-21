#include "TF1.h"
#include "TFitResult.h"
#include "TMath.h"
#include "TH1D.h"
#include "TCanvas.h"
// #include "GH1D.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <string>

class ExclusionPeakFit {
    public:
        ExclusionPeakFit(int target, int deg, double exlo, double exhi, double fl, double fh){
            if (exlo > exhi) std::swap(exlo,exhi);
            if (fl > fh) std::swap(fl,fh);

            target_peak = target;
            degree = deg;
            fitLo = fl;
            fitHi = fh;

            regions[target] = std::make_pair(exlo, exhi); //add target peak to the list of regions
        } 

        ExclusionPeakFit() : target_peak(0), degree(0), fitLo(0), fitHi(0) {}

        void AddRegion(int pk, double exlo, double exhi){
            if (exlo > exhi) std::swap(exlo,exhi);
            if (!regions.count(pk)) regions[pk] = std::make_pair(exlo,exhi);
        }

        double operator() (double *x, double *par){
            bool reject = false;
            for (auto r : regions){
                if (x[0] >= r.second.first && x[0] <= r.second.second) reject = true;
            }
            if (reject) TF1::RejectPoint();
            
            return par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0];
        }

        std::pair<double,double> GetPeakCounts(TH1D *h, int optdeg=-1) {
            // TF1 *fitfunc = new TF1(Form("exfit_%d_%f_%f",degree,fitLo,fitHi),this,fitLo,fitHi,degree);
            int fitdeg = degree;
            if (optdeg != -1) fitdeg = optdeg;
            TF1 *fitfunc = new TF1("tempfit",this,fitLo,fitHi,fitdeg);

            TFitResultPtr fitres;
            int nloops = 0;
            do {
                fitres = h->Fit(fitfunc,"SMLRQN");
                nloops++;
            } while (!fitres->IsValid() && nloops < 4);

            fitParameters = fitfunc->GetParameters();

            //do background subtraction
            double integLo = h->GetBinLowEdge(h->FindBin(regions[target_peak].first)); //make sure the limits are the same limits...
            double integHi = h->GetBinLowEdge(h->FindBin(regions[target_peak].second)+1); //...as your histogram
            double bkg = fitfunc->Integral(integLo,integHi) / h->GetBinWidth(1);
            double dBkg = fitfunc->IntegralError(integLo,integHi,fitres->GetParams(),fitres->GetCovarianceMatrix().GetMatrixArray()) / h->GetBinWidth(1);
                
            //integrate the histogram over the excluded region
            double sigbkg = h->Integral(h->FindBin(regions[target_peak].first),h->FindBin(regions[target_peak].second));
            double dSigbkg = TMath::Sqrt(sigbkg);
            double sig = sigbkg-bkg;
            double dSig = TMath::Sqrt(dBkg*dBkg + dSigbkg*dSigbkg);

            return std::make_pair(sig,dSig);
        }

        int GetPeak() {return target_peak;} 
        double GetFitLo() {return fitLo;}
        double GetFitHi() {return fitHi;}
        double *GetFitParameters() {return fitParameters;}

        std::pair<double, double> bkgCounts(TH1 *h) {
            TSpectrum *s = new TSpectrum();
            TH1D *hbkg = (TH1D*) s->Background(h,12);
            TH1D *hc = (TH1D*) h->Clone("hdiff");
            hc->Add(hbkg);

            //integrate the histogram over the excluded region
            double sig = h->Integral(h->FindBin(regions[target_peak].first),h->FindBin(regions[target_peak].second));

            return std::make_pair(sig,TMath::Sqrt(sig));
        }

    private:
        int degree, target_peak;
        double fitLo, fitHi;
        double *fitParameters;
        std::map<int, std::pair<double,double>> regions;
};

std::pair<double,double> calcAsym(std::pair<double,double> para, std::pair<double,double> perp, std::pair<double,double> correction = std::make_pair(-12345,-12345)){
    if (correction.first != -12345) {
        double tmp_para = para.first * correction.first;
        double tmp_para_err = tmp_para * TMath::Sqrt(std::pow(correction.second/correction.first,2) + std::pow(para.second/para.first,2));
        para.first = tmp_para;
        para.second = tmp_para_err;
    }

    double diff = perp.first - para.first;
    double sum = perp.first + para.first;
    double sumerr = TMath::Sqrt(perp.second*perp.second + para.second*para.second);

    double asym = diff/sum;
    double asymerr = std::abs(asym)*TMath::Sqrt(sumerr*sumerr/sum/sum + sumerr*sumerr/diff/diff);

    return std::make_pair(asym,asymerr);
}

std::pair<double,double> PerpOverPara(std::pair<double,double> para, std::pair<double,double> perp){
    double ratio = perp.first/para.first;
    double ratioErr = ratio * TMath::Sqrt(std::pow(perp.second/perp.first,2) + std::pow(para.second/para.first,2));

    return std::make_pair(ratio,ratioErr);
}

void source_pol_analyze(TFile *f, std::string source) {
    std::vector<std::string> sortMethods = {"MainInt", "ComptonSort", "Tracking"};

    // ================================================================================
    // SOURCE SETTINGS  ===============================================================
    // ================================================================================
    std::map<int, ExclusionPeakFit> correction_regions_56Co = {
        {511, ExclusionPeakFit(511,3,501,522,492,527)},
        {846, ExclusionPeakFit(846,4,836,854,819,875)},
        {1037, ExclusionPeakFit(1037,3,1028,1045,1022,1055)},
        {1238, ExclusionPeakFit(1238,4,1221,1248,1200,1272)},
        {2034, ExclusionPeakFit(2034,3,2022,2043,1974,2060)},
        {2598, ExclusionPeakFit(2598,4,2577,2605,2550,2641)},
        {3253, ExclusionPeakFit(3253,3,3230,3260,3213,3300)}
    };

    std::map<int, ExclusionPeakFit> correction_regions_152Eu = {
        {344, ExclusionPeakFit(344,3,333,352,308,361)},
        {443, ExclusionPeakFit(443,2,435,453,424,460)},
        {779, ExclusionPeakFit(779,4,772,784,750,802)},
        {867, ExclusionPeakFit(867,4,860,874,846,888)},
        {964, ExclusionPeakFit(964,3,953,970,939,980)},
        {1112, ExclusionPeakFit(1112,2,1101,1118,1026,1138)},
        {1408, ExclusionPeakFit(1408,4,1392,1416,1353,1422)},
    };

    //add more regions for 56Co
    correction_regions_56Co[846].AddRegion(861,858,863);
    correction_regions_56Co[846].AddRegion(834,833,835);
    correction_regions_56Co[1238].AddRegion(1260,1255,1264);
    correction_regions_56Co[2034].AddRegion(2014,1999,2021);
    correction_regions_56Co[2598].AddRegion(2615,2606,2623);
    correction_regions_56Co[3253].AddRegion(3272,3260,3281);
    
    //add more regions for 152Eu
    correction_regions_152Eu[344].AddRegion(315,312,318);
    correction_regions_152Eu[344].AddRegion(327,322,332);
    correction_regions_152Eu[779].AddRegion(767,762,771);
    correction_regions_152Eu[779].AddRegion(788,786,791);
    correction_regions_152Eu[779].AddRegion(795,793,797);
    correction_regions_152Eu[1112].AddRegion(1086,1045,1094);
    correction_regions_152Eu[1112].AddRegion(1123,1118,1128);
    correction_regions_152Eu[1408].AddRegion(1378,1376,1380);
    correction_regions_152Eu[1408].AddRegion(1367,1363,1370);
    
    std::map<int, ExclusionPeakFit> correction_regions; 
    if (source == "152Eu") correction_regions = correction_regions_152Eu;
    else if (source == "56Co") correction_regions = correction_regions_56Co;
    else {
        printf("Error: unknown source '%s'. Valid options are '56Co' and '152Eu'.\n", source.c_str());
        return;
    }
    
    std::map<std::string, std::vector<int>> all_gatedpeaks;
    all_gatedpeaks["56Co"]  = {1238, 2598};
    all_gatedpeaks["152Eu"] = {344};

    std::vector<int>   &gatedpeaks         = all_gatedpeaks[source];
    printf("Using source: %s with %zu correction regions and %zu gated peaks.\n",
        source.c_str(), correction_regions.size(), gatedpeaks.size());

    std::map<std::string, ExclusionPeakFit> peakOfInterest;
    peakOfInterest["152Eu"] = ExclusionPeakFit(779,2,773,783,720,805);
    peakOfInterest["152Eu"].AddRegion(742,732,752);
    peakOfInterest["152Eu"].AddRegion(765,759,772);
    peakOfInterest["56Co"] = ExclusionPeakFit(846,2,836,854,819,875);
    // ================================================================================
    // ================================================================================
    // ================================================================================
    
    std::string input_name = std::string(f->GetName());
    input_name = input_name.substr(0, input_name.length() - 5);
    TFile *fitHists = new TFile(Form("fit_%s.root",input_name.c_str()),"RECREATE");

    //final data storage
    std::map<std::string, std::map<int, std::pair<double,double>>> asymmetry_data;

    for (auto method : sortMethods) {
        printf("\n===== Method: %s =====\n", method.c_str());

        fitHists->mkdir(Form("%s", method.c_str()))->cd();

        std::vector<TH1D*> totHists = {
            (TH1D*) f->Get(Form("gretina/all_core_energy_%s_para", method.c_str())),
            (TH1D*) f->Get(Form("gretina/all_core_energy_%s_perp", method.c_str()))
        };

        // std::vector<TH1D*> totHists = {
        //     (TH1D*) f->Get("gretina/all_core_energy"),
        //     (TH1D*) f->Get(Form("gretina/all_core_energy_%s_perp", method.c_str()))
        // };

        // get the counts for the correction histograms
        std::map<int, std::vector<std::pair<double,double>>> correction_counts;
        std::vector<TF1 *> fitfuncs;
        for (auto& r : correction_regions) {
            std::pair<double,double> para_counts = r.second.GetPeakCounts(totHists[0]);
            fitfuncs.push_back(new TF1(Form("fit_para_%d",r.first),"[0] + [1]*x + [2]*x*x + [3]*x*x*x", r.second.GetFitLo(), r.second.GetFitHi()));
            fitfuncs.back()->SetParameters(r.second.GetFitParameters());
            totHists[0]->GetListOfFunctions()->Add(fitfuncs.back());
            
            std::pair<double,double> perp_counts = r.second.GetPeakCounts(totHists[1]);
            fitfuncs.push_back(new TF1(Form("fit_perp_%d",r.first),"[0] + [1]*x + [2]*x*x + [3]*x*x*x", r.second.GetFitLo(), r.second.GetFitHi()));
            fitfuncs.back()->SetParameters(r.second.GetFitParameters());
            totHists[1]->GetListOfFunctions()->Add(fitfuncs.back());

            correction_counts[r.first] = {para_counts, perp_counts};
            // printf("%d %.0f %.0f\n",r.first,para_counts.first,para_counts.second);
        }
        totHists[0]->Write();
        totHists[1]->Write();

        //build the vectors for plotting the correction histograms
        std::vector<double> pkEnergy, perpParaRatio, perpParaRatioErr, asym, asymErr;
        for (auto it = correction_counts.begin(); it != correction_counts.end(); ++it) {
            std::pair<double,double> tempRatio = PerpOverPara(it->second[0], it->second[1]);
            std::pair<double,double> tempAsym = calcAsym(it->second[0], it->second[1]);
            pkEnergy.push_back(it->first);
            perpParaRatio.push_back(tempRatio.first);
            perpParaRatioErr.push_back(tempRatio.second);
            asym.push_back(tempAsym.first);
            asymErr.push_back(tempAsym.second);
        }

        //make the error scatter plots
        TGraphErrors *grRatio  = new TGraphErrors((int) pkEnergy.size(), &pkEnergy[0],  &perpParaRatio[0],  nullptr, &perpParaRatioErr[0]);
        TGraphErrors *grAsym = new TGraphErrors((int) pkEnergy.size(), &pkEnergy[0], &asym[0], nullptr, &asymErr[0]);
        TCanvas *cc = new TCanvas(Form("%s", method.c_str()),Form("%s", method.c_str()));

        //fit the ratio plot to get the correction factor
        TF1 *fconst = new TF1(Form("corr_ratio_%s", method.c_str()), "[0]");
        grRatio->Fit(fconst, "Q");
        // grRatio->Draw();
        std::pair<double,double> correctionFactor = std::make_pair(fconst->GetParameter(0), fconst->GetParError(0));
        printf("Correction factor: %f +/- %f\n", correctionFactor.first, correctionFactor.second);

        grAsym->Fit(fconst, "Q");
        grAsym->Draw("A*");
        printf("Apparative asymmetry: %f +/- %f\n", fconst->GetParameter(0), fconst->GetParError(0));

        grAsym->Write();
        grRatio->Write();

        std::map<int, std::vector<TH1D*>> histos;
        std::map<int, std::vector<std::pair<double,double>>> counts;
        for (auto ss : gatedpeaks) {
            histos[ss] = {
                (TH1D*) f->Get(Form("gretina/all_crrl_%d_core_energy_%s_para", ss, method.c_str())),
                (TH1D*) f->Get(Form("gretina/all_crrl_%d_core_energy_%s_perp", ss, method.c_str()))
            };

            std::vector<std::pair<double,double>> tempcounts;
            for (int i = 0; i < 2; i++) {
                if (ss == 2598) 
                    tempcounts.push_back(peakOfInterest[source].GetPeakCounts(histos[ss][i],3));
                else 
                    tempcounts.push_back(peakOfInterest[source].GetPeakCounts(histos[ss][i]));
                fitfuncs.push_back(new TF1(Form("fit_perp_%d",peakOfInterest[source].GetPeak()),"[0] + [1]*x + [2]*x*x + [3]*x*x*x", peakOfInterest[source].GetFitLo(), peakOfInterest[source].GetFitHi()));
                fitfuncs.back()->SetParameters(peakOfInterest[source].GetFitParameters());
                histos[ss][i]->GetListOfFunctions()->Add(fitfuncs.back());
            }
                
            printf("COUNTS %f +/- %f | %f +/- %f\n", tempcounts[0].first, tempcounts[0].second, tempcounts[1].first, tempcounts[1].second);
            
            counts[ss] = tempcounts;
            histos[ss][0]->Write();
            histos[ss][1]->Write();
        }

        pkEnergy.clear(); asym.clear(); asymErr.clear();
        int nloops = 1;
        for (auto ss : gatedpeaks) {
            std::pair<double,double> tempAsym = calcAsym(counts[ss][0], counts[ss][1], correctionFactor);
            asymmetry_data[method].insert({ss, tempAsym});
            printf("%d %f +/- %f\n", ss, tempAsym.first, tempAsym.second);
            asym.push_back(std::abs(tempAsym.first));
            asymErr.push_back(tempAsym.second);
            pkEnergy.push_back(nloops);
            nloops++;
        }

        grAsym = new TGraphErrors((int) pkEnergy.size(), &pkEnergy[0], &asym[0], nullptr, &asymErr[0]);
        new TCanvas();
        grAsym->SetTitle(Form("Asymmetry - %s", method.c_str()));
        grAsym->Draw("A*");
    }

    printf("\nAsymmetries of the %d keV peak...\nCorrelated with: ",peakOfInterest[source].GetPeak());
    for (auto pk : gatedpeaks) printf("     %4d keV     \t",pk);
    std::cout<<std::endl;
    for (auto method : sortMethods){
        printf("%-11s\t",method.c_str());
        for (auto pk : gatedpeaks) {
            printf("%7.4f +/- %-6.4f\t",asymmetry_data[method][pk].first,asymmetry_data[method][pk].second);
        }
        printf("\n");
    }

    fitHists->Close();
}