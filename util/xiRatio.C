#include "TFile.h"
#include "TH1D.h"
#include "TList.h"
#include "TKey.h"
#include "TGretina.h"
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

void groupCrystals(std::vector<std::pair<int,std::string>> &crstl, std::string mode){
    int N = crstl.size();
    if (mode.find("quad") != std::string::npos) {
        for (int i=0; i < N; i++){
            crstl[i].first = crstl[i].first/4 - 1; 
        }
        std::sort(crstl.begin(),crstl.end());
    }
    if (mode.find("ring") != std::string::npos) {
        TGretina *gret = new TGretina();
        for (int i=0; i < N; i++){
            crstl[i].first = gret->GetRingNumber(crstl[i].first);
        }
        std::sort(crstl.begin(),crstl.end());
    }
    return;
}

void xiRatio(TFile *fsource, TFile *fdata, double EsrcLo, double EsrcHi, double EdatLo, double EdatHi, std::string grouping="crys", int binning = 1) {
    std::vector<std::string> dirsource = findPolDirectory(fsource);
    std::vector<std::string> dirdata = findPolDirectory(fdata);

    // if (!(dirsource.size() == 1 && dirdata.size() > 0)) return;
    int dataidx = 0;
    if (dirdata.size() > 1) {
        std::cout<<"Multiple polarization channels, pick one\n";
        for (int i=0 ; i < dirdata.size(); i++) printf("%d -> %s\n",i,dirdata[i].c_str());
        std::cout<<"index: ";
        std::cin >> dataidx;
    }
    std::vector<std::pair<int,std::string>> srcCrystals;
    std::vector<std::pair<int,std::string>> datCrystals;
    processDirectory(fsource,dirsource[0],srcCrystals);
    processDirectory(fdata,dirdata[dataidx],datCrystals);

    if (srcCrystals.size() != datCrystals.size()) {std::cout<<"Error: source and data have different numbers of crystals"; return;}

    double Esrc[2] = {EsrcLo,EsrcHi};
    double Edat[2] = {EdatLo,EdatHi};

    GH1D *hnorm;
    GH1D *hsrc;
    GH1D *hdat;
    bool firstReset = true;
    int N = srcCrystals.size();  
    int groupID = -1; 
    int nGroups = 0;

    //split them into the chosen groups
    groupCrystals(srcCrystals,grouping);
    groupCrystals(datCrystals,grouping);

    GH1D *stemp;
    GH1D *dtemp;
    for (int i=0; i < N; i++){
        if (srcCrystals[i].first != datCrystals[i].first) break;
        if (srcCrystals[i].first != groupID){
            nGroups++;
            if (i != 0) { 
                if (firstReset) {
                    hsrc = (GH1D*) stemp->Clone("source");
                    hdat = (GH1D*) dtemp->Clone("data");
                    dtemp->Divide(stemp);
                    hnorm = (GH1D*) dtemp->Clone("norm"); 
                    firstReset = false;
                }
                else {
                    hsrc->Add(stemp);
                    hdat->Add(dtemp);
                    dtemp->Divide(stemp);
                    hnorm->Add(dtemp);
                } 
            }
            stemp = makeProjection(fsource,dirsource[0],srcCrystals[i].second,0,Esrc); stemp->Sumw2();
            dtemp = makeProjection(fdata,dirdata[dataidx],datCrystals[i].second,0,Edat); dtemp->Sumw2();
            groupID = srcCrystals[i].first;
        }
        else {
            stemp->Add(makeProjection(fsource,dirsource[0],srcCrystals[i].second,0,Esrc));
            dtemp->Add(makeProjection(fdata,dirdata[dataidx],datCrystals[i].second,0,Edat));
        }  
    }
    dtemp->Divide(stemp);
    hnorm->Add(dtemp);
    hnorm->Scale(1.0/nGroups/binning);

    TF1 *fitfunc = new TF1("pol","[0]*(1-[1]*TMath::Cos(2*x))",0,TMath::TwoPi());
    hnorm->Fit(fitfunc,"Q");
    double A0 = fitfunc->GetParameter(1);
    double A0_err = fitfunc->GetParError(1);
    printf("%5.3f +/- %5.3f -> %4.2f\n",A0,A0_err,A0/A0_err);

    if (binning != 1){
        hnorm->Rebin(binning);
        hsrc->Rebin(binning);
        hdat->Rebin(binning);
    }
    hnorm->Fit(fitfunc);
    hnorm->Draw();

    new GCanvas();
    double scaling = hsrc->GetEntries()/hdat->GetEntries() * 4.0/5;
    hdat->SetLineColor(kRed);
    hdat->Scale(scaling);
    double ymax = std::max(hdat->GetMaximum(),hsrc->GetMaximum());
    double ymin = std::min(hdat->GetMaximum(),hsrc->GetMaximum());

    hsrc->GetYaxis()->SetRangeUser(ymin*0.8,ymax*1.1);
    if (binning != 1) {hsrc->Draw(); hdat->Draw("same");}
    else {hsrc->Draw("hist"); hdat->Draw("samehist");}
}

void xiRatio(TFile *fsource, TFile *fdata, std::string hname="", std::string fileout="", int binning=-1){

    std::vector<std::string> dirsource = findPolDirectory(fsource);
    std::vector<std::string> dirdata = findPolDirectory(fdata);

    if (dirsource.size() == 1 && dirdata.size() > 0){
        int dataidx = 0;
        if (dirdata.size() > 1) {
            std::cout<<"Multiple polarization channels, pick one\n";
            for (int i=0 ; i < dirdata.size(); i++) printf("%d -> %s\n",i,dirdata[i].c_str());
            std::cout<<"index: ";
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
        GCanvas *canv = new GCanvas();
        canv->Divide(1,2);
        canv->cd(1);
        drawDat->Draw();
        canv->cd(2);
        drawSrc->Draw();
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
        // TF1 *ffit = new TF1("mypolfit","([0] + [1]*x)*(1-[2]*TMath::Cos(2*x))",0,TMath::Pi());
        // div->Fit(ffit);
        // TFile *fout = new TFile(fileout.c_str(),"RECREATE");
        // fout->cd();
        // div->SetNameTitle(hname.c_str(),hname.c_str());
        // div->Write();
        div->Draw();
    }
}