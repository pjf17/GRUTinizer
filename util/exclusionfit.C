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

class ExclusionFit {
    public:
        ExclusionFit() {}
        
        ExclusionFit(double pXlo,double pXhi) {regions.push_back(std::make_pair(pXlo,pXhi));}
        
        void AddRegion(double pXlo, double pXhi) {
            if (pXlo > pXhi) std::swap(pXlo,pXhi);

            //check if region is already added
            bool nodupe = true;
            for (auto r : regions){
                if (pXlo*10000 + pXhi == r.first*10000 + r.second) nodupe = false;
            }   
            if (nodupe) regions.push_back(std::make_pair(pXlo,pXhi));
            else printf("%4.0f <> %-4.0f region already added!\n",pXlo,pXhi);
        }

        void LoadRegions(std::string filename){
            std::ifstream input(filename.c_str());
            std::string line;
            double xlo, xhi;
            while (getline(input,line)){
                if (line.find('#') != std::string::npos) continue;
                std::stringstream ss(line);
                ss >> xlo >> xhi;
                AddRegion(xlo,xhi);
                printf("Added %4.0f <> %-4.0f\n",xlo,xhi);
            }
        }

        void EditRegion(double oldXlo, double oldXhi, double newXlo, double newXhi){
            if (newXlo > newXhi) std::swap(newXlo,newXhi);
            if (oldXlo > oldXhi) std::swap(oldXlo,oldXhi);

            //check if region is already added
            bool nodupe = true;
            int nregions = regions.size();
            for (int i=0; i < nregions; i++){
                if (oldXlo*10000 + oldXhi == regions[i].first*10000 + regions[i].second) {
                    nodupe = false;
                    regions[i].first = newXlo;
                    regions[i].second = newXhi;
                    printf("Updated %4.0f <> %-4.0f to %4.0f <> %-4.0f\n",oldXlo,oldXhi,newXlo,newXhi);
                }
            }   
            if (nodupe) printf("Region does not exist\n");
        }

        void OnlyFitInRegions(bool flag = true){
            rejectSwitch = !flag;
        }

        void ListRegions(){
            for (auto r : regions){
                printf("%4.0f <> %-4.0f\n",r.first,r.second);
            }
        }

        void ClearRegions(){regions.clear();}
        
        double operator() (double *x, double *par){
            bool reject = !rejectSwitch; //reject switch is default true
            for (auto r : regions){
                if (x[0] >= r.first && x[0] <= r.second) reject = rejectSwitch;
            }
            if (reject) TF1::RejectPoint();
            // return par[0]*TMath::Exp(par[1]*x[0]) + par[2]*TMath::Exp(par[3]*x[0]);
            return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
        }
        
        void Fit(TH1D *h, int order, double elo, double ehi, double xlo, double xhi){
            AddRegion(elo,ehi);
            Fit(h,order,xlo,xhi);
            regions.pop_back();
        }
        
        void Fit(TH1D *h, int order, double xlo, double xhi, bool useprev=false){
            //copy og hist to sig and bkg hists before fit
            hsub = (TH1D*) h->Clone("bkg_subtracted");
            hbkg = (TH1D*) h->Clone("bkg");
            hbkg->SetLineColor(kRed);
            hsub->SetLineColor(kBlack);

            TF1 *f = new TF1("myfitfunc",this,xlo,xhi,order);
            if (useprev) f->SetParameters(prevpars);
            // else {
            //     f->SetParameter(0,1.90826e+03);
            //     f->SetParameter(1,-0.000828433);
            //     if (order > 2){
            //         f->SetParameter(2,10000);
            //         f->SetParameter(3,-0.01);
            //     }
            // }
            
            TFitResultPtr fitres = h->Fit(f,"SMLR");
            int npars = fitres->NPar();
            for (int p=0; p < npars; p++){
                prevpars[p] = fitres->Parameter(p);
            }

            printf("Chi2  /  ndf: %6.2f / %-3d = %-5.3f\n",fitres->Chi2(),fitres->Ndf(),fitres->Chi2()/fitres->Ndf());

            //make entire fit region bkg and subtract from data hist
            int start = h->FindBin(xlo);
            int end = h->FindBin(xhi);
            for (int b=start; b <= end; b++){
                hbkg->SetBinContent(b,f->Eval(hbkg->GetBinCenter(b)));
            }
            hsub->Add(hbkg,-1);
            double subXmin = hsub->GetXaxis()->GetXmin();
            double subXmax = hsub->GetXaxis()->GetXmax();

            //calculate the counts in each region
            printf("\nExclusion regions within [%4.0f,%-4.0f]\n",xlo,xhi);
            for (auto r : regions){
                //if both of the exclude boundaries are within the fitting range, keep it
                if ( (xlo < r.first && r.first < xhi)  && (xlo < r.second && r.second < xhi) ){
                    printf("=== region %4.0f <> %-4.0f ==============\n",r.first,r.second);
                    //get bkg counts over the excluded region
                    double integLo = h->GetBinLowEdge(h->FindBin(r.first)); //make sure the limits are the same limits...
                    double integHi = h->GetBinLowEdge(h->FindBin(r.second)+1); //...as your histogram
                    double bkg = f->Integral(integLo,integHi) / h->GetBinWidth(1);
                    double dBkg = f->IntegralError(integLo,integHi,fitres->GetParams(),fitres->GetCovarianceMatrix().GetMatrixArray()) / h->GetBinWidth(1);
                    
                    //integrate the histogram over the excluded region
                    double sigbkg = h->Integral(h->FindBin(r.first),h->FindBin(r.second));
                    double dSigbkg = TMath::Sqrt(sigbkg);
                    double sig = sigbkg-bkg;
                    double dSig = TMath::Sqrt(dBkg*dBkg + dSigbkg*dSigbkg);

                    //get the FWHM of the peak
                    double fwhm = getFWHM(hsub,r.first,r.second);

                    //output the results
                    printf("Total Counts: %10.2f +/- %-7.2f\n",sigbkg,dSigbkg);
                    printf("Bkg   Counts: %10.2f +/- %-7.2f\n",bkg,dBkg);
                    printf("Peak  Counts: %10.2f +/- %-7.2f\n",sig,dSig);
                    printf("Peak   FWHM : %10.2f\n\n",fwhm);
                } 
            }
 
            //draw the bkg and signal hists on a new convas
            TCanvas *resultDisplay = new TCanvas("Bkg Subtracted","Bkg Subtracted");
            hsub->Draw();
            hbkg->Draw("same");
        }
    private:
        std::vector<std::pair<double,double>> regions;
        TH1D *hsub, *hbkg = 0;
        double prevpars[4];
        bool rejectSwitch = true;

        double getFWHM(TH1D *h, double lo, double hi){
            TF1 *fhist = ((GH1D*) h)->ConstructTF1();
            double hLo = h->GetXaxis()->GetXmin();
            double hHi = h->GetXaxis()->GetXmax();

            int lobin = h->FindBin(lo);
            int hibin = h->FindBin(hi);
            h->GetXaxis()->SetRange(lobin,hibin);
            double xmax = h->GetBinCenter(h->GetMaximumBin());
            double halfmax = h->GetMaximum()/2.0;
            double xlo = xmax;
            double xhi = xmax;
            while (fhist->Eval(xlo) > halfmax){
                xlo -= 0.001;
            }
            while (fhist->Eval(xhi) > halfmax){
                xhi += 0.001;
            }
            h->GetXaxis()->SetRange(hLo,hHi);
            return xhi-xlo;
        }
};