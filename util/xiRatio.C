#include "TFile.h"
#include "TH1D.h"
#include "TList.h"
#include "TKey.h"
#include <vector>
#include <string>

TH1D *loadHist(TFile *f){
    TList *hlist = f->GetListOfKeys();
    TIter next(hlist);
    TKey *key;
    while (TObject *obj = next()){
        key = (TKey *)obj;
        std::string classType(key->GetClassName());
        if (classType.find("H1") != std::string::npos) break;
    }
    return (TH1D*) key->ReadObj();
}

std::vector<std::string> findPolDirectory(TFile *f){
    TList *l = f->GetListOfKeys();
    int N = l->GetEntries();
    std::vector<std::string> output;
    for (int n=0; n < N; n++){
        std::string keyname = std::string(l->At(n)->GetName());
        if (keyname.find("polarization") != std::string::npos) output.push_back(keyname);
    }
    return output;
}

void processDirectory(TFile *f, std::string dirname,std::vector<std::pair<int,std::string>> &names){
    TList *li = f->GetDirectory(dirname.c_str())->GetListOfKeys();
    int N = li->GetEntries();
    for (int i=0 ; i < N; i++) {
        std::string histname = std::string(li->At(i)->GetName());
        int cryNum = std::stoi(histname.substr(histname.size()-2,2));
        if (!(cryNum > 3 && cryNum < 124)) continue;
        names.push_back(std::make_pair(cryNum,histname));
    }

    std::sort(names.begin(),names.end());
    return;
}

GH1D *makeProjection(TFile *f, std::string dir, std::string name, bool axis, double *egate = nullptr){
    // axis == 0 is X
    // axis == 1 is Y
    GH2D *hh = (GH2D*) f->Get(Form("%s/%s",dir.c_str(),name.c_str()));
    if (egate){
        int lobin, hibin;
        if (axis) {
            lobin = hh->GetXaxis()->FindBin(egate[0]);
            hibin = hh->GetXaxis()->FindBin(egate[1]);
            return hh->ProjectionY("_py",lobin,hibin);
        }
        else {
            lobin = hh->GetYaxis()->FindBin(egate[0]);
            hibin = hh->GetYaxis()->FindBin(egate[1]);
            return hh->ProjectionX("_px",lobin,hibin);
        }
    }
    else {
        if (axis) return hh->ProjectionY();
        else return hh->ProjectionX();
    }
}

void xiRatio(TFile *fsource, TFile *fdata, double EsrcLo, double EsrcHi, double EdatLo, double EdatHi) {
    std::vector<std::string> dirsource = findPolDirectory(fsource);
    std::vector<std::string> dirdata = findPolDirectory(fdata);

    if (!(dirsource.size() == 1 && dirdata.size() > 0)) return;
    int dataidx = 0;
    if (dirdata.size() > 1) {
        std::cout<<"Multiple polarization channels, pick one\n";
        for (int i=0 ; i < dirdata.size(); i++) printf("%d %s\n",i,dirdata[i].c_str());
        std::cin >> dataidx;
    }
    std::vector<std::pair<int,std::string>> srcCrystals;
    std::vector<std::pair<int,std::string>> datCrystals;
    processDirectory(fsource,dirsource[0],srcCrystals);
    processDirectory(fdata,dirdata[dataidx],datCrystals);

    if (srcCrystals.size() != datCrystals.size()) {std::cout<<"Error: source and data have different numbers of crystals"; return;}

    double Esrc[2] = {EsrcLo,EsrcHi};
    double Edat[2] = {EdatLo,EdatHi};

    // GH1D *drawSrc = makeProjection(fsource,dirsource[0],srcCrystals[0].second,0,Esrc);
    // GH1D *drawDat = makeProjection(fdata,dirdata[dataidx],datCrystals[0].second,0,Edat);
    // drawSrc->Sumw2(); drawDat->Sumw2();
    GH1D *stemp;
    GH1D *dtemp;
    GH1D *hnorm;
    bool firstReset = true;
    int N = srcCrystals.size();  
    int groupID = -1; 
    int nGroups = -1;
    TGretina *gret = new TGretina();
    for (int i=0; i < N; i++){
        if (srcCrystals[i].first != datCrystals[i].first) break;
        int currentID = srcCrystals[i].first/4-1;
        if (currentID != groupID){
            nGroups++;
            if (i != 0) { 
                dtemp->Divide(stemp);
                if (firstReset) {hnorm = (GH1D*) dtemp->Clone("norm"); firstReset = false;}
                else hnorm->Add(dtemp);
            }
            stemp = makeProjection(fsource,dirsource[0],srcCrystals[i].second,0,Esrc); stemp->Sumw2();
            dtemp = makeProjection(fdata,dirdata[dataidx],datCrystals[i].second,0,Edat); dtemp->Sumw2();
            groupID = currentID;
        }
        else {
            stemp->Add(makeProjection(fsource,dirsource[0],srcCrystals[i].second,0,Esrc));
            dtemp->Add(makeProjection(fdata,dirdata[dataidx],datCrystals[i].second,0,Edat));
        }  
    }
    hnorm->Scale(1.0/nGroups);
    hnorm->Draw();
}

void xiRatio(TFile *fsource, TFile *fdata, int binning=-1){

    std::vector<std::string> dirsource = findPolDirectory(fsource);
    std::vector<std::string> dirdata = findPolDirectory(fdata);

    if (dirsource.size() == 1 && dirdata.size() > 0){
        int dataidx = 0;
        if (dirdata.size() > 1) {
            std::cout<<"Multiple polarization channels, pick one\n";
            for (int i=0 ; i < dirdata.size(); i++) printf("%d %s\n",i,dirdata[i].c_str());
            std::cin >> dataidx;
        }
        std::vector<std::pair<int,std::string>> srcCrystals;
        std::vector<std::pair<int,std::string>> datCrystals;
        processDirectory(fsource,dirsource[0],srcCrystals);
        processDirectory(fdata,dirdata[dataidx],datCrystals);

        if (srcCrystals.size() != datCrystals.size()) {std::cout<<"Error: source and data have different numbers of crystals"; return;}

        //add all projected hists together and draw for energy gate selection
        GH1D *drawSrc = makeProjection(fsource,dirsource[0],srcCrystals[0].second,1);
        GH1D *drawDat = makeProjection(fdata,dirdata[dataidx],datCrystals[0].second,1);
        int N = srcCrystals.size();   
        for (int i=1; i < N; i++){
            if (srcCrystals[i].first != datCrystals[i].first) break;
            drawSrc->Add(makeProjection(fsource,dirsource[0],srcCrystals[i].second,1));
            drawDat->Add(makeProjection(fdata,dirdata[dataidx],datCrystals[i].second,1));
        }
        drawDat->Draw();
        new GCanvas();
        drawSrc->Draw();

        // std::cout<<"Select Energy Gates for the data\n";
        // double Edat[2];
        // std::cin>>Edat[0];std::cin>>Edat[1];
        // if (Edat[0] > Edat[1]) std::swap(Edat[0],Edat[1]);
        
        // std::cout<<"Select Energy Gates for the source\n";
        // double Esrc[2];
        // std::cin>>Esrc[0];std::cin>>Esrc[1];
        // if (Esrc[0] > Esrc[1]) std::swap(Esrc[0],Esrc[1]);
        
        // drawSrc = makeProjection(fsource,dirsource[0],srcCrystals[0].second,0,Esrc);
        // drawDat = makeProjection(fdata,dirdata[dataidx],datCrystals[0].second,0,Edat);
        // for (int i=1; i < N; i++){
        //     if (srcCrystals[i].first != datCrystals[i].first) break;
        //     drawSrc->Add(makeProjection(fsource,dirsource[0],srcCrystals[i].second,0,Esrc));
        //     drawDat->Add(makeProjection(fdata,dirdata[dataidx],datCrystals[i].second,0,Edat));
        // }
        // drawSrc->Draw();
    }

    else {
        std::cout<<"Not doing crystal by crystal\n";
        TH1D *hSource = loadHist(fsource);
        hSource->Sumw2();

        TH1D *hData = loadHist(fdata);
        hData->Sumw2();
        TH1D *div = (TH1D*) hData->Clone();

        if (binning != -1){
            div->Rebin(binning);
            hSource->Rebin(binning);
            hData->Rebin(binning);
        }

        div->Divide(hSource);
        div->Draw();
    }
}