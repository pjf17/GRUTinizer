#include "TFile.h"
#include "TString.h"
#include "THStack.h"

const int NCOLORS = 4;
const int NFILLS = 2;
// int colors[NCOLORS] = {kOrange,kViolet,kPink,kAzure};
// int fills[NFILLS] = {1001,69420,3002};
double fills[NFILLS] = {1.00,0.50};
int colors[NCOLORS] = {kViolet-3,kSpring-3,kCyan-3,kOrange-3};
// int colors[NCOLORS] = {kYellow,kRed,kMagenta,kBlue};
// int colors[12] = {kYellow,kMagenta,kRed,kGreen,kBlue,kCyan,kOrange,kSpring,kTeal,kPink,kAzure,kViolet};
// int colors[12] = {kPink-2,kPink,kMagenta,kViolet,kBlue+2,kBlue,kAzure,kCyan+3,kCyan,kTeal,kGreen,kSpring};

THStack *hs;
TH1D *hdata;

void drawStack(TFile *f) {
    //read the file
    TList *hlist = f->GetListOfKeys();
    hs= new THStack("hs","");
    TIter colnext(hlist);
    TH1D *hstopped = nullptr;
    int col = 0;
    while (TObject *obj = colnext()){
        TKey *key = (TKey *)obj;
        std::string classType(key->GetClassName());
        if (classType.find("H1") != std::string::npos){
            TString keyname = TString(key->GetName());
            if (keyname.Contains("Stopped")) {
                TH1D *hh = (TH1D*) key->ReadObj();
                hstopped = (TH1D*) hh->Clone("stopped");
            }
        }
    }

    // Double_t Red[4]    = { 1.00, 0.00, 0.00, 0.00};
    // Double_t Green[4]  = { 0.00, 0.00, 0.00, 0.00};
    // Double_t Blue[4]   = { 0.00, 0.00, 0.00, 1.00};
    // Double_t Length[4] = { 0.00, 0.33, 0.66, 1.00 };
    // Int_t FI = TColor::CreateGradientColorTable(4, Length, Red, Green, Blue, 12);
    bool noStopped = true;
    TIter next(hlist);
    while (TObject *obj = next()){
        TKey *key = (TKey *)obj;
        std::string classType(key->GetClassName());
        if (classType.find("H1") != std::string::npos){
            TString keyname = TString(key->GetName());
            if (keyname.Contains("Peak")) {
                TH1D *hh = (TH1D*) key->ReadObj();
                //cycle through colors & change slightly every repeat
                int cc = colors[col%NCOLORS];// - col/(NCOLORS*2);

                // hh->SetFillColor(cc);
                hh->SetFillColorAlpha(cc,fills[(col/NCOLORS)%NFILLS]);
                // else hh->SetFillStyle(fs); 
                // hh->SetLineColor(kBlack);
                hh->SetLineColor(cc + 5);
                hs->Add(hh);
                col++;
                // if (noStopped && hh->GetBinCenter(hh->GetMaximumBin()) < 600) {
                    // hstopped->SetFillColor(kBlack);
                    // hstopped->SetLineColor(kBlack);
                    // hs->Add(hstopped);
                    // noStopped = false;
                // }
            }
            // if (keyname.Contains("Stopped")) {
                // hstopped = (TH1D*) key->ReadObj();
                // hstopped = (TH1D*) hh->Clone("stopped");
                // hstopped->SetFillColor(kBlack);
                // hstopped->SetLineColor(kBlack);
                // hs->Add(hh);
            // }
            if (keyname.Contains("Exp_Bkg")){
                TH1D *hbk = (TH1D*) key->ReadObj();
                int nBins = hbk->GetNbinsX();
                TH1D *hh = new TH1D("bkg","bkg",nBins,hbk->GetXaxis()->GetXmin(),hbk->GetXaxis()->GetXmax());
                // TH1D *hh = (TH1D*) hbk->Clone("bkg");
                for (int i=1; i <= nBins; i++){
                    double nCounts = hbk->GetBinContent(i);
                    hh->SetBinContent(i,nCounts);
                }
                hh->SetFillColor(kGray+2);
                hh->SetLineColor(kBlack);
                hs->Add(hh);
                
                if (hstopped){
                    hstopped->SetFillColor(kBlack);
                    hstopped->SetLineColor(kBlack);
                    hs->Add(hstopped);
                }
            }
            if (keyname.Contains("gam_dop") || keyname.Contains("gamma_dop"))
                hdata = (TH1D*) key->ReadObj();
        }
    }
    
    // hstopped->SetFillColor(kBlack);
    // hstopped->SetLineColor(kBlack);
    // hs->Add(hstopped);

    // TCanvas *canv = new TCanvas("canv","canv",700,700);
    // TCanvas *canv = new TCanvas("canv","canv",800,400);
    TCanvas *canv = new TCanvas("canv","canv",800,270);
    
    // canv->SetLeftMargin(0.065);
    // canv->SetBottomMargin(0.14);
    ((TF1 *) hdata->GetListOfFunctions()->At(1))->SetLineWidth(2);
    ((TF1 *) hdata->GetListOfFunctions()->At(1))->SetLineColor(kAzure);

    hs->Draw("hist");

    hs->GetXaxis()->SetTitle("Energy [keV]");
    hs->GetXaxis()->SetTitleFont(43);
    hs->GetXaxis()->SetTitleSize(22);
    hs->GetXaxis()->SetTitleOffset(0.84);
    hs->GetXaxis()->SetLabelFont(43);
    hs->GetXaxis()->SetLabelSize(18);
    hs->GetXaxis()->SetLabelOffset(0.01);

    hs->GetYaxis()->SetTitle(Form("Counts / %d keV",(int) hdata->GetBinWidth(1)));
    hs->GetYaxis()->SetTitleFont(43);
    hs->GetYaxis()->SetTitleSize(22);
    hs->GetYaxis()->SetTitleOffset(0.5);
    hs->GetYaxis()->SetLabelFont(43);
    hs->GetYaxis()->SetLabelSize(18);
    hs->GetYaxis()->SetLabelOffset(0.01);

    hdata->SetLineColor(kBlack);
    hdata->Draw("same");

    gPad->SetTopMargin(0.01);
    gPad->SetBottomMargin(0.165);
    gPad->SetLeftMargin(0.085);
    gPad->SetRightMargin(0.02);

    // gPad->SetTopMargin(0.01);
    // gPad->SetBottomMargin(0.18);
    // gPad->SetLeftMargin(0.09);
    // gPad->SetRightMargin(0.02);
    
    canv->Modified();
    canv->Update();
}

void setAxisrange(double xlo, double xhi, bool doX){
    if (xlo > xhi) std::swap(xlo,xhi);
    //get the histograms on the canvas
    if (doX) 
        hs->GetXaxis()->SetRangeUser(xlo,xhi);
    else {
        hs->SetMinimum(xlo);
        hs->SetMaximum(xhi);
    }

    gPad->Modified();
    gPad->Update();
    return;
}

void setXrange(double xlo, double xhi){
    setAxisrange(xlo,xhi,true);
}

void setYrange(double xlo, double xhi){
    setAxisrange(xlo,xhi,false);
}

void drawSubFig(std::string letter, double xx=0.15, double yy=0.88){
    TText text;
    text.SetTextAlign(31);
    text.SetTextFont(43);
    text.SetTextSize(25);
    text.DrawTextNDC(xx,yy,letter.c_str());
}