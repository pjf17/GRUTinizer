#include <map>
#include <vector>
#include <string>
#include <iostream>

#include "TFile.h"
#include "TH1F.h"
#include "TKey.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TF1.h"
#include "TCanvas.h"

#include "MultiPlotter.h"

//PUBLIC FUNCTIONS

void MultiPlotter::Add(TH1F* pHist){
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
    mHistos.insert( std::pair<std::string,TH1F*>(pHist->GetName(),pHist) );
}

void MultiPlotter::Add(TFile *f, const char *hname){
    TH1F* h = (TH1F*) f->Get(hname);
    this->Add(h);
}

void MultiPlotter::Add(TDirectoryFile *f){
    // print out hists you can choose
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
                        TH1F *hout = (TH1F*) key->ReadObj();
                        this->Add(hout);
                    }
                }
            }
        }
    }
}

void MultiPlotter::Erase(std::string key){
    mHistos.erase(key);
}

void MultiPlotter::Clear(){
    mHistos.clear();
    mNHistos = 0;
    mYMax = 0.0;
    mCustomColors = false;
}

void MultiPlotter::List(){
    std::map<std::string, TH1F*>::iterator it = mHistos.begin();
    std::map<std::string, TH1F*>::iterator end = mHistos.end();
    while (it != end){
        std::cout<<it->first<<std::endl;
        it++;
    }
}

TH1F* MultiPlotter::GetClone(std::string key){
    return (TH1F*) mHistos[key]->Clone();
}

TH1F* MultiPlotter::Get(std::string key){
    return mHistos[key];
}

void MultiPlotter::SetLineWidth(int w){ 
    mLineWidth = w;
}

void MultiPlotter::SetLineColor(std::string key, int c){
    mCustomColors = true;
    mHistos[key]->SetLineColor(c);
}

void MultiPlotter::SetRange(double xlo, double xhi){
    std::map<std::string, TH1F*>::iterator it = mHistos.begin();
    std::map<std::string, TH1F*>::iterator end = mHistos.end();
    while (it != end){
        it->second->GetXaxis()->SetRangeUser(xlo,xhi);
        it++;
    }
}

void MultiPlotter::IterateLineStyle(){
    std::map<std::string, TH1F*>::iterator it = mHistos.begin();
    std::map<std::string, TH1F*>::iterator end = mHistos.end();
    int style = 1;
    while (it != end){
        it->second->SetLineStyle(style);
        style++;
        it++;
    }
}

void MultiPlotter::Norm(std::string mode){
    std::map<std::string, TH1F*>::iterator it = mHistos.begin();
    std::map<std::string, TH1F*>::iterator end = mHistos.end();

    double norm = 0.0;
    if (mode == "area") norm = it->second->Integral();
    else if (mode == "height") norm = it->second->GetMaximum();
    while (it != end){ 
        if (mode == "area") it->second->Scale(norm/it->second->Integral());
        else if (mode == "height") it->second->Scale(norm/it->second->GetMaximum());
        it++;
    }

    //reset the max hist parameter
    mYMax = 0.0;
    mMaxKey = "";
}

void MultiPlotter::Fit(std::string key, TF1 *f){
    mHistos[key]->Fit(f);
}

void MultiPlotter::Fit(TF1 *f){
    std::map<std::string, TH1F*>::iterator it = mHistos.begin();
    std::map<std::string, TH1F*>::iterator end = mHistos.end();
    while (it != end){
        new TCanvas();
        it->second->Fit(f);
        it++;
    }
}

void MultiPlotter::Draw(std::string key){
    mHistos[key]->Draw("hist");
}

void MultiPlotter::Draw(){
    doSetLineWidth();
    if (mYMax == 0.0) SortMax();
    std::map<std::string, TH1F*>::iterator max = mHistos.find(mMaxKey);
    
    //initialize legend
    TLegend *leg = new TLegend(0.72,0.99 - mNHistos*0.06,0.99,0.99);
    gStyle->SetOptStat(0);
    
    //draw all histos
    std::map<std::string, TH1F*>::iterator it = max;
    std::map<std::string, TH1F*>::iterator end = mHistos.end();
    int nloops = 0;
    while (it != end || nloops == 0){
        if (it == max && nloops > 0) {
            it++;
            continue;
        }
        
        if (!mCustomColors) it->second->SetLineColor(mColors[nloops%9]);
        leg->AddEntry(it->second,it->second->GetName(),"l");
        
        if (nloops == 0){
            it->second->Draw("hist");
            it = mHistos.begin();
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

//PRIVATE FUNCTIONS

void MultiPlotter::SortMax(){
    std::map<std::string, TH1F*>::iterator it = mHistos.begin();
    std::map<std::string, TH1F*>::iterator end = mHistos.end();
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
    std::map<std::string, TH1F*>::iterator it = mHistos.begin();
    std::map<std::string, TH1F*>::iterator end = mHistos.end();
    while (it != end){
        it->second->SetLineWidth(mLineWidth);
        it++;
    }
}
