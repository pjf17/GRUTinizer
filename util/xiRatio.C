#include "TFile.h"
#include "TH1D.h"
#include "TList.h"
#include "TKey.h"

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

void xiRatio(TFile *fsource, TFile *fdata, int binning=-1){
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
    // new TCanvas();
    // hData->SetLineColor(kRed);
    // hSource->Draw();
    // hData->Draw("same");
}