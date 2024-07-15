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
        
        ExclusionFit(int pRebinFactor) {rebinFactor = pRebinFactor;}

        ExclusionFit(int pRebinFactor, bool pIsExpo) {rebinFactor = pRebinFactor; isExpo = pIsExpo;}
        
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

        void SetExpo() {isExpo = true;}
        void SetPoly() {isExpo = false;}
        void SetRebinFactor (int rbf) {rebinFactor = rbf; currentHistName = "";}
        
        double operator() (double *x, double *par){
            bool reject = !rejectSwitch; //reject switch is default true
            for (auto r : regions){
                if (x[0] >= r.first && x[0] <= r.second) reject = rejectSwitch;
            }
            if (reject) TF1::RejectPoint();
            if (isExpo) return TMath::Exp(par[0] + par[1]*x[0]) + TMath::Exp(par[2] + par[3]*x[0]);
            else        return par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0];
        }
        
        void Fit(TH1D *h, int order, double elo, double ehi, double xlo, double xhi){
            AddRegion(elo,ehi);
            Fit(h,order,xlo,xhi);
            regions.pop_back();
        }
        
        void Fit(TH1D *h, int order, double xlo, double xhi){
            //get current hist name
            bool newName = false;
            if (currentHistName != std::string(h->GetName())) {
                currentHistName = std::string(h->GetName());
                newName = true;
            }

            //rebin hist if needed
            if (rebinFactor > 1 && newName) h->Rebin(rebinFactor);

            //copy og hist to sig and bkg hists before fit
            hsub = (TH1D*) h->Clone("bkg_subtracted");
            hbkg = (TH1D*) h->Clone("bkg");
            hbkg->SetLineColor(kRed);
            hsub->SetLineColor(kBlack);

            TF1 *f = new TF1("myfitfunc",this,xlo,xhi,order);

            if (isExpo) {
                f->SetParameter(0,7);
                f->SetParameter(1,-0.000828433);
                if (order > 2) {
                    f->SetParameter(2,9);
                    f->SetParameter(3,-0.000828433);
                }
            }
            
            TFitResultPtr fitres;
            int nloops = 0;
            do {
                fitres = h->Fit(f,"SMLR");
                nloops++;
                // int npars = fitres->NPar();
                // for (int p=0; p < npars; p++){
                //     prevpars[p] = fitres->Parameter(p);
                // }
                // f->SetParameters(prevpars);
            } while (!fitres->IsValid() && nloops < 4);

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
                    double fwhmErr = 2.355/6*TMath::Sqrt(2)*h->GetBinWidth(1);

                    //output the results
                    printf("Total Counts: %10.2f +/- %-7.2f\n",sigbkg,dSigbkg);
                    printf("Bkg   Counts: %10.2f +/- %-7.2f\n",bkg,dBkg);
                    printf("Peak  Counts: %10.2f +/- %-7.2f\n",sig,dSig);
                    printf("Peak   FWHM : %10.3f +/- %-7.3f\n\n",fwhm,fwhmErr);
                } 
            }
 
            //draw the bkg and signal hists on a new convas
            TCanvas *resultDisplay = new TCanvas("Bkg Subtracted","Bkg Subtracted");
            hsub->Draw();
            hbkg->Draw("same");
            newName = false;
        }
    private:
        std::vector<std::pair<double,double>> regions;
        TH1D *hsub, *hbkg = 0;
        double prevpars[4];
        int rebinFactor = 1;
        bool rejectSwitch = true;
        bool isExpo = false;
        std::string currentHistName = "";

        double interpolate(double x1, double y1, double x2, double y2, double y){
            double m = (y2-y1)/(x2-x1);
            double b = (y1*x2-y2*x1)/(x2-x1);
            return (y-b)/m;
        }

        double getFWHM(TH1D *h, double lo, double hi){
            // TF1 *fhist = ((GH1D*) h)->ConstructTF1();
            double hLo = h->GetXaxis()->GetXmin();
            double hHi = h->GetXaxis()->GetXmax();

            int lobin = h->FindBin(lo);
            int hibin = h->FindBin(hi);
            h->GetXaxis()->SetRange(lobin,hibin);
            double maxBin = h->GetMaximumBin();
            double halfmax = h->GetMaximum()/2.0;
            double loBin = maxBin;
            double hiBin = maxBin;
            while (h->GetBinContent(loBin) > halfmax) loBin--;
            while (h->GetBinContent(hiBin) > halfmax) hiBin++;
            h->GetXaxis()->SetRange(hLo,hHi);
            std::cout<<loBin<<" "<<hiBin<<std::endl;
             
            return interpolate(h->GetBinCenter(hiBin),h->GetBinContent(hiBin),h->GetBinCenter(hiBin-1),h->GetBinContent(hiBin-1),halfmax)- 
                   interpolate(h->GetBinCenter(loBin),h->GetBinContent(loBin),h->GetBinCenter(loBin+1),h->GetBinContent(loBin+1),halfmax);
        }
};