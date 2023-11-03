#include <map>
#include <vector>
#include <string>
#include <iostream>

#include "TFile.h"
#include "TH1.h"
#include "TKey.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TF1.h"
#include "TCanvas.h"
#include "GGaus.h"
#include "GRootCommands.h"
#include "TMath.h"
#include "TFitResult.h"
#include "THStack.h"

#include "MultiPlotter.h"

//PUBLIC FUNCTIONS

void MultiPlotter::Add(TH1* pHist){
    mYMax = 0.0;
    mNHistos++;
    if (mHistos.count(std::string(pHist->GetName()))){
        std::cout<<"Duplicate hist name, add modifier to make it unique."<<std::endl;
        std::string modifier;
        std::cout<<"Modifier: ";
        std::cin>>modifier;
        std::string hname = std::string(pHist->GetName()) + modifier;
        pHist->SetName(hname.c_str());
    }
    std::cout<<"added "<<pHist->GetName()<<std::endl;
    mHistos.insert( std::pair<std::string,TH1*>(pHist->GetName(),pHist) );
}

void MultiPlotter::Add(TFile *f, const char *hname){
    TH1* h = (TH1*) f->Get(hname);
    this->Add(h);
}

void MultiPlotter::Add(TDirectoryFile *f){
    // print out hists you can choose
    std::cout<<"*** Reading "<<f->GetName()<<" ***"<<std::endl;
    TList *hlist = f->GetListOfKeys();
    std::vector<TKey *> sortedList;
    TIter next(hlist);
    while (TObject *obj = next()){
        TKey *key = (TKey *)obj;
        std::string classType(key->GetClassName());
        if (key->IsFolder() || classType.find("H1") != std::string::npos 
            || classType.find("H2") != std::string::npos){
            sortedList.push_back(key);
        }
    }
    //if only one object and it's not a folder, just add it
    if (sortedList.size() == 1 && !(sortedList[0]->IsFolder())){
        TH1 *hout = (TH1*) sortedList[0]->ReadObj();
        this->Add(hout);
    }
    else {
        //sort alphabetically
        std::sort(sortedList.begin(),sortedList.end(), 
            [] (TKey *a, TKey *b) {
                int result = strcmp(a->GetName(),b->GetName());
                if (result >= 0) return false;
                else return true;
            });
        printf("Available objects\n");
        
        int nKeys = (int) sortedList.size();
        for (int i=0; i<nKeys; i++){
            TKey *key = (TKey *) sortedList[i];
            if (key->IsFolder()) 
                printf("--> d%d %s%s%s\n",i,"\033[1;34m",key->GetName(),"\033[m");
            else 
                printf("--> %d %s%s%s\n",i,"",key->GetName(),"");
        }
        
        char input[50];
        while (true){
            std::cout<<"add hist ('q' to quit): ";
            std::cin>>input;
            if (strcmp(input,"q") == 0) break;
            else {
                std::vector<int> nums;
                bool isfolder = ParseInput(input,nums);
                if (isfolder){
                    TKey *key = sortedList[nums.back()];
                    std::cout<<key->GetName()<<std::endl;
                    TDirectoryFile *dir = (TDirectoryFile*) key->ReadObj();
                    this->Add(dir);
                    break;
                } else {
                    int N = nums.size();
                    for (int i=0; i < N; i++){
                        if (nums[i] < nKeys && nums[i] >= 0){
                            TKey *key = sortedList[nums[i]];
                            TH1 *hout = (TH1*) key->ReadObj();
                            this->Add(hout);
                        }
                    }
                }
            }
        }
    }
}

void MultiPlotter::Erase(std::string key){
    if (Exists(key)) mHistos.erase(key);
}

void MultiPlotter::Clear(){
    mHistos.clear();
    mNHistos = 0;
    mYMax = 0.0;
    mCustomColors = false;
}

void MultiPlotter::List(){
    std::map<std::string, TH1*>::iterator it = mHistos.begin();
    std::map<std::string, TH1*>::iterator end = mHistos.end();
    while (it != end){
        std::cout<<it->first<<std::endl;
        it++;
    }
}

TH1* MultiPlotter::GetClone(std::string key){
    return (TH1*) mHistos[key]->Clone();
}

TH1* MultiPlotter::Get(std::string key){
    return mHistos[key];
}

void MultiPlotter::SetLineWidth(int w){ 
    mLineWidth = w;
}

void MultiPlotter::SetLineColor(std::string key, int c){
    if (Exists(key)){
        mCustomColors = true;
        mHistos[key]->SetLineColor(c);   
    }
}

void MultiPlotter::SetYrange(double ylo, double yhi){
    mYlo = ylo;
    mYhi = yhi;
}

void MultiPlotter::SetRange(double xlo, double xhi){
    mXlo = xlo;
    mXhi = xhi;
}

void MultiPlotter::ResetRange(){
    mXlo = -123;
    mXhi = -123;
}

void MultiPlotter::IterateLineStyle(){
    std::map<std::string, TH1*>::iterator it = mHistos.begin();
    std::map<std::string, TH1*>::iterator end = mHistos.end();
    int style = 1;
    while (it != end){
        it->second->SetLineStyle(style);
        style++;
        it++;
    }
}

void MultiPlotter::Norm(std::string mode){
    std::map<std::string, TH1*>::iterator it = mHistos.begin();
    std::map<std::string, TH1*>::iterator end = mHistos.end();

    double norm = 0.0;
    if (mode == "area") norm = it->second->Integral();
    else if (mode == "height") norm = it->second->GetMaximum();
    it++;
    while (it != end){ 
        if (mode == "area") it->second->Scale(norm/it->second->Integral());
        else if (mode == "height") it->second->Scale(norm/it->second->GetMaximum());
        it++;
    }

    //reset the max hist parameter
    mYMax = 0.0;
    mMaxKey = "";
}

void MultiPlotter::Norm(double lo, double hi){
    std::map<std::string, TH1*>::iterator it = mHistos.begin();
    std::map<std::string, TH1*>::iterator end = mHistos.end();
    int binLo = it->second->FindBin(lo);
    int binHi = it->second->FindBin(hi);

    double norm = it->second->Integral(binLo, binHi);
    it++;
    while (it != end){ 
        it->second->Scale(norm/it->second->Integral(binLo, binHi));
        it++;
    }

    //reset the max hist parameter
    mYMax = 0.0;
    mMaxKey = "";
}

void MultiPlotter::Fit(std::string key, TF1 *f, double xlo, double xhi){
    f->SetRange(xlo,xhi);
    if (Exists(key)) mHistos[key]->Fit(f,"R");
}

void MultiPlotter::Fit(TF1 *f, double xlo, double xhi){
    std::map<std::string, TH1*>::iterator it = mHistos.begin();
    std::map<std::string, TH1*>::iterator end = mHistos.end();
    f->SetRange(xlo,xhi);
    while (it != end){
        new TCanvas();
        std::cout<<it->first<<std::endl;
        it->second->Fit(f,"R");
        it++;
    }
}

void MultiPlotter::FitPearson(double xlo, double xhi, double h, double c, double w){
    std::map<std::string, TH1*>::iterator it = mHistos.begin();
    std::map<std::string, TH1*>::iterator end = mHistos.end();
    TF1 *pearsonbkg = new TF1("pearsonbkg","[4] + [5]*x + [0]/TMath::Power( 1 + (TMath::Power(2,1/[3]) - 1)*TMath::Power((2*x-2*[1])/[2],2) ,[3])",xlo,xhi);
    TF1 *bkg = new TF1("bkg","[0] + [1]*x",xlo,xhi);
    pearsonbkg->SetParameter(0,h);
    pearsonbkg->SetParameter(1,c);
    pearsonbkg->SetParameter(2,w);
    pearsonbkg->SetParameter(3,2);
    pearsonbkg->SetParameter(4,500);
    pearsonbkg->SetParameter(5,-2);
    pearsonbkg->SetParLimits(0,0,10000000);

    while (it != end){
        new TCanvas();
        std::cout<<it->first<<std::endl;
        TFitResultPtr fitres = it->second->Fit(pearsonbkg,"RS");
        double area = pearsonbkg->Integral(xlo,xhi) / it->second->GetBinWidth(1);
        double dArea = pearsonbkg->IntegralError(xlo,xhi,fitres->GetParams(),fitres->GetCovarianceMatrix().GetMatrixArray()) / it->second->GetBinWidth(1);
        
        bkg->SetParameters(&pearsonbkg->GetParameters()[4]);
        double bgArea = bkg->Integral(xlo,xhi) / it->second->GetBinWidth(1);
        double bgDArea = bkg->IntegralError(xlo,xhi,&(fitres->GetParams())[4],fitres->GetCovarianceMatrix().GetSub(4,5,4,5).GetMatrixArray()) / it->second->GetBinWidth(1);

        double counts = area - bgArea;
        double dCounts = TMath::Sqrt(dArea*dArea + bgDArea*bgDArea);
        printf("Chi2: %f NDF: %d Chi2/NDF: %f\n",fitres->Chi2(),fitres->Ndf(),fitres->Chi2()/(1.0*fitres->Ndf()));
        printf("Counts\tCountsErr\tArea\tAreaErr\tBkg\tBkgErr\n");
        printf("%f\t%f\t%f\t%f\t%f\t%f\n",counts,dCounts,area,dArea,bgArea,bgDArea);
        it++;
    }
}

void MultiPlotter::FitGaus(double xlo, double xhi, Option_t *opt){
    std::map<std::string, TH1*>::iterator it = mHistos.begin();
    std::map<std::string, TH1*>::iterator end = mHistos.end();
    std::string sOpt = opt;
    sOpt.append("no-print");
    printf("%-30s %6s %4s %3s %7s %5s\n","name","cntrd","fwhm","err","area","err");
    while (it != end){
        GGaus *fitR = GausFit(it->second,xlo,xhi,sOpt.c_str());
        printf("%-30s %6.2f %4.2f %3.2f %7.1f %5.1f\n",it->second->GetName(),fitR->GetCentroid(),fitR->GetFWHM(),fitR->GetFWHMErr(),fitR->GetArea(),fitR->GetAreaErr());
        it++;
    }
}

void MultiPlotter::Draw(std::string key){
    if (Exists(key)) mHistos[key]->Draw("hist");
}

void MultiPlotter::Draw(int ndraw, int noffset){
    doSetLineWidth();
    if (mYMax == 0.0) SortMax();
    std::map<std::string, TH1*>::iterator max = mHistos.find(mMaxKey);
    
    //initialize legend
    double ylo = 0.99 - mNHistos*0.06;
    double xlo = 0.72;
    if (ylo < 0.3) {
        ylo = 0.01;
        xlo = 0.9;
    }
    TLegend *leg = new TLegend(xlo,ylo,0.99,0.99);
    gStyle->SetOptStat(0);

    //draw all histos
    std::map<std::string, TH1*>::iterator it = max;
    std::map<std::string, TH1*>::iterator end = mHistos.end();
    int nloops = 0;
    while ((it != end || nloops == 0) && nloops < ndraw){
        if (it == max && nloops > 0) {
            it++;
            continue;
        }
        
        if (!mCustomColors) it->second->SetLineColor(mColors[nloops%12]);
        leg->AddEntry(it->second,it->second->GetName(),"l");

        if (mXhi != mXlo) it->second->GetXaxis()->SetRangeUser(mXlo,mXhi);
        if (mYhi != mYlo) it->second->GetYaxis()->SetRangeUser(mYlo,mYhi);
        
        if (nloops == 0){
            it->second->Draw("hist");
            it = mHistos.begin();
            for (int j=0; j < noffset; j++) it++;
            nloops++;
        }
        else{
            it->second->Draw("hist same");
            it++;
            nloops++;
        }
    }

    leg->Draw("same");
}

void MultiPlotter::Rebin(int bg){
    std::map<std::string, TH1*>::iterator it = mHistos.begin();
    std::map<std::string, TH1*>::iterator end = mHistos.end();
    while (it != end){
        it->second->Rebin(bg);
        it++;
    }
    mYMax = 0.0;
}

void MultiPlotter::Scroll(std::string control){
    std::map<std::string, TH1*>::iterator begin = mHistos.begin();
    std::map<std::string, TH1*>::iterator end = mHistos.end();
    static std::map<std::string, TH1*>::iterator it = mHistos.begin();
    int xlo = it->second->GetXaxis()->GetFirst();
    int xhi = it->second->GetXaxis()->GetLast();
    int ylo = it->second->GetYaxis()->GetFirst();
    int yhi = it->second->GetYaxis()->GetLast();
    if (control.compare(".")==0 && it!=end) it++;
    if (control.compare(",")==0 && it!=begin) it--;
    it->second->GetXaxis()->SetRange(xlo,xhi);
    it->second->GetYaxis()->SetRange(ylo,yhi);
    it->second->Draw();
}

//pick one histogram and divide all others by it
void MultiPlotter::RatioToHist(std::string key){
    mYMax = 0;
    std::map<std::string, TH1*>::iterator it = mHistos.begin();
    std::map<std::string, TH1*>::iterator div = mHistos.find(key);
    std::map<std::string, TH1*>::iterator end = mHistos.end();
    while (it != end){
        if (it != div) it->second->Divide(div->second);
        it++;
    }
    mHistos.erase(key);
    return;
}

//PRIVATE FUNCTIONS

void MultiPlotter::SortMax(){
    std::map<std::string, TH1*>::iterator it = mHistos.begin();
    std::map<std::string, TH1*>::iterator end = mHistos.end();
    while (it != end){
        double tmpMax = it->second->GetMaximum();
        if (tmpMax > mYMax){
            mYMax = tmpMax;
            mMaxKey = it->first; 
        }
        it++;
    }
    return;
}

bool MultiPlotter::ParseInput(char *input, std::vector<int> &nums){
    if (input[0] == 'd'){
        nums.push_back(atoi(&input[1]));
        return true;
    }   

    //go through comma separated fields
    std::vector<char *> tokens;
    char *tok;
    tok = strtok(input,",");
    while (tok != NULL){
        tokens.push_back(tok);
        tok = strtok(NULL,",");
    }

    //find dashed included fields and add numbers
    int ntokens = tokens.size();
    for (int i=0; i < ntokens; i++){
        int start, end = -1;
        
        tok = strtok(tokens[i],"-");
        start = atoi(tok);
        tok = strtok(NULL,"-");
        if (tok != NULL){
            end = atoi(tok);
            for (int j=start; j < end+1; j++){
                nums.push_back(j);
            }
        } else {
            nums.push_back(start);
        }
    }
    return false;
}

void MultiPlotter::doSetLineWidth(){
    std::map<std::string, TH1*>::iterator it = mHistos.begin();
    std::map<std::string, TH1*>::iterator end = mHistos.end();
    while (it != end){
        it->second->SetLineWidth(mLineWidth);
        it++;
    }
}

bool MultiPlotter::Exists(std::string key){
    bool does_exist = mHistos.count(key);
    if (!does_exist) {
        std::cout<<"Error, histogram not found in list"<<std::endl;
    }
    return does_exist;
}
