#include <map>
#include <vector>
#include <string>
#include <iostream>

#include "TFile.h"
#include "TText.h"
#include "TH1.h"
#include "TKey.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TF1.h"
#include "TCanvas.h"
#include "GGaus.h"
#include "GPeak.h"
#include "GRootCommands.h"
#include "TMath.h"
#include "TFitResult.h"
#include "THStack.h"
#include "TAxis.h"

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

THStack *MultiPlotter::CreateStack(){
    THStack *out = new THStack("hs","");
    std::map<std::string, TH1*>::iterator it = mHistos.begin();
    std::map<std::string, TH1*>::iterator end = mHistos.end();
    int nloops = 0;
    while (it != end){
        if (!mCustomColors) it->second->SetFillColor(mColors[nloops%12]);
        out->Add(it->second);
        it++;
        nloops++;
    }
    return out;
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

void MultiPlotter::SetXrange(double xlo, double xhi){
    mXlo = xlo;
    mXhi = xhi;
}

void MultiPlotter::SetLegendEntry(std::string key, std::string label, std::string opt){
    if (mUseDefaultLegend) {
        double ylo = 0.99 - mNHistos*0.06;
        double xlo = 0.72;
        if (ylo < 0.3) {
            ylo = 0.01;
            xlo = 0.9;
        }
        mLeg = new TLegend(xlo,ylo,0.99,0.99); 
    }
    mUseDefaultLegend = false;
    mLeg->AddEntry(mHistos[key],label.c_str(),opt.c_str());
}

void MultiPlotter::ResetRange(){
    mXlo = -123;
    mXhi = -123;
}

void MultiPlotter::SetFill(double alpha){
    if (alpha < 0.0) mDrawFill = false;
    else {
        mFillAlpha = alpha;
        mDrawFill = true;
    }
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

    if (!it->second->GetSumw2()) it->second->Sumw2();

    double norm = 0.0;
    if (mode == "area") norm = it->second->Integral();
    else if (mode == "height") norm = it->second->GetMaximum();
    else if (mode == "unit") {norm = it->second->GetNbinsX(); it->second->Scale(norm/it->second->Integral());}
    it++;
    while (it != end){ 
        double scaleFactor = 0.0;
        if (mode == "area" || mode == "unit") scaleFactor = it->second->Integral();
        else if (mode == "height") scaleFactor = it->second->GetMaximum();
        it->second->Scale(norm/scaleFactor);
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

void MultiPlotter::FitPeak(double xlo, double xhi, Option_t *opt){
    std::map<std::string, TH1*>::iterator it = mHistos.begin();
    std::map<std::string, TH1*>::iterator end = mHistos.end();
    std::string sOpt = opt;
    sOpt.append("no-print");
    printf("%-30s %6s %4s %3s %7s %5s\n","name","cntrd","fwhm","err","area","err");
    while (it != end){
        GPeak *fitR = PhotoPeakFit(it->second,xlo,xhi,sOpt.c_str());
        printf("%-30s %6.2f %4.2f %3.2f %7.1f %5.1f\n",it->second->GetName(),fitR->GetCentroid(),fitR->GetFWHM(),fitR->GetFWHMErr(),fitR->GetArea(),fitR->GetAreaErr());
        it++;
    }
}

void MultiPlotter::Integral(double lo, double hi) {
    std::map<std::string, TH1*>::iterator it = mHistos.begin();
    std::map<std::string, TH1*>::iterator end = mHistos.end();
    while (it != end){
        int xlo = it->second->FindBin(lo);
        int xhi = it->second->FindBin(hi);
        printf("%-30s %6.2f\n",it->second->GetName(),it->second->Integral(xlo,xhi));
        it++;
    }
}

void MultiPlotter::Add(std::string key, double scale){
    if (!Exists(key)) return;
    mYMax = 0;
    std::map<std::string, TH1*>::iterator it = mHistos.begin();
    std::map<std::string, TH1*>::iterator end = mHistos.end();
    TH1 *hadd = (TH1*) mHistos[key]->Clone("htemp");
    mHistos.erase(key);
    while (it != end){
        it->second->Add(hadd,scale);
        it++;
    }
    gFile->Delete("htemp");
    return;
}

// void MultiPlotter::Draw(std::string opt, std::string key){
//     if (Exists(key)) mHistos[key]->Draw("hist");
// }

void MultiPlotter::Draw(int nx, int ny, int wx, int wy, bool mergeX, bool mergeY){
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    // gStyle->SetPadTopMargin(0.05);
    // gStyle->SetPadRightMargin(0.05);
    // gStyle->SetPadBottomMargin(0.16);
    // gStyle->SetPadLeftMargin(0.12);
    // gStyle->SetPadBorderMode(0);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    
    if (nx*ny != mNHistos) {
        printf("ERROR: NUMBER OF HISTOGRAMS DOES NOT MATCH DIMENSIONS\n");
        return;
    }
    
    //list the histograms
    std::map<std::string, TH1*>::iterator it = mHistos.begin();
    std::map<std::string, TH1*>::iterator end = mHistos.end();
    while (it != end){
        std::cout<<it->first<<std::endl;
        it++;
    }

    //show the canvas layout
    std::string outformat = "";
    int nplots = 1;
    for (int yy=0; yy < ny; yy++){
        outformat += "|";
        for (int xx=0; xx < nx; xx++){
            outformat += std::string(Form("%3d |",nplots));
            nplots++;
        }
        outformat += "\n";
    }
    nplots--;
    std::cout<<outformat;
    
    //assign the plots to the layout
    printf("\nAssign the histos to the layout\n");
    it = mHistos.begin();
    std::vector<std::pair<int,TH1*>> drawingstack;
    while (it != end){
        std::cout<<it->first;
        int digit;
        std::cout<<" -> ";
        std::cin>>digit;
        drawingstack.push_back(std::make_pair(digit,it->second));
        it++;
    }
    std::sort(drawingstack.begin(),drawingstack.end());

    //validate
    for (int i=0; i < nplots; i++){
        if (i+1 != drawingstack[i].first) {
            printf("ERROR: USER INPUT INDICES ARE NOT UNIQUE\n");
            return;
        }
    }
    
    TCanvas *c1 = new TCanvas("canv","canv",wx,wy);
    CanvasPartition(c1,nx,ny,0.15,0.15,0.15,0.05);
    std::vector<TPad *> pads;
    int iplt = 0;
    for (Int_t i = 0; i < nx; i++) {
        for (Int_t j = 0; j < ny; j++) {
            TH1 *h = drawingstack[iplt].second;
            c1->cd(0);

            // Get the pads previously created.
            TPad *ipad = (TPad *)c1->FindObject(TString::Format("pad_%d_%d", i, j).Data());
            pads.push_back(ipad);
            pads.back()->Draw();
            pads.back()->SetFillStyle(4000);
            pads.back()->SetFrameFillStyle(4000);
            pads.back()->cd();

            // Size factors
            Float_t xFactor = pads.front()->GetAbsWNDC() / pads.back()->GetAbsWNDC();
            Float_t yFactor = pads.front()->GetAbsHNDC() / pads.back()->GetAbsHNDC();

            TH1F *hFrame = (TH1F *)h->Clone(TString::Format("h_%d_%d", i, j).Data());

            // y axis range
            hFrame->SetMinimum(0.0001); // do not show 0
            hFrame->SetMaximum(1.2 * h->GetMaximum());

            // Format for y axis
            hFrame->GetYaxis()->SetLabelFont(43);
            hFrame->GetYaxis()->SetLabelSize(16);
            hFrame->GetYaxis()->SetLabelOffset(0.02);
            hFrame->GetYaxis()->SetTitleFont(43);
            hFrame->GetYaxis()->SetTitleSize(16);
            hFrame->GetYaxis()->SetTitleOffset(2);

            hFrame->GetYaxis()->CenterTitle();
            hFrame->GetYaxis()->SetNdivisions(505);

            // TICKS Y Axis
            hFrame->GetYaxis()->SetTickLength(xFactor * 0.04 / yFactor);

            // Format for x axis
            hFrame->GetXaxis()->SetLabelFont(43);
            hFrame->GetXaxis()->SetLabelSize(16);
            hFrame->GetXaxis()->SetLabelOffset(0.02);
            hFrame->GetXaxis()->SetTitleFont(43);
            hFrame->GetXaxis()->SetTitleSize(16);
            hFrame->GetXaxis()->SetTitleOffset(1);
            hFrame->GetXaxis()->CenterTitle();
            hFrame->GetXaxis()->SetNdivisions(505);

            // TICKS X Axis
            hFrame->GetXaxis()->SetTickLength(yFactor * 0.06 / xFactor);

            // Draw cloned histogram with individual settings
            hFrame->Draw();

            iplt++;
        //   TText text;
        //   text.SetTextAlign(31);
        //   text.SetTextFont(43);
        //   text.SetTextSize(10);
        //   text.DrawTextNDC(XtoPad(0.9), YtoPad(0.8), gPad->GetName());
        }
    }
    c1->cd();

    /*
    TCanvas *c1 = new TCanvas("canv","canv",wx,wy);
    if (mergeX || mergeY) 
        c1->Divide(nx,ny,0.0,0.0);
    else 
        c1->Divide(nx,ny);

    // c1->cd(0);
    // c1->SetLeftMargin(0.2);
    // c1->SetRightMargin(0.2);
    int iplt = 0;
    for (int yy=0; yy < ny; yy++){
        for (int xx=0; xx < nx; xx++){
            c1->cd(iplt+1);
            if (yy == 0){
                gPad->SetTopMargin(0.1);
            }
            if (yy == ny-1) {
                gPad->SetBottomMargin(0.15);
            }
            // gPad->SetLeftMargin(0.1);
            // gPad->SetRightMargin(0.1);
            if (mXhi != mXlo) drawingstack[iplt].second->GetXaxis()->SetRangeUser(mXlo,mXhi);
            if (mYhi != mYlo) drawingstack[iplt].second->GetYaxis()->SetRangeUser(mYlo,mYhi);
            if (mXLabelSize != 0.04) drawingstack[iplt].second->GetXaxis()->SetLabelSize(mXLabelSize);
            if (mYLabelSize != 0.04) drawingstack[iplt].second->GetYaxis()->SetLabelSize(mYLabelSize);
            if (mXtitle != "") drawingstack[iplt].second->GetXaxis()->SetTitle(mXtitle.c_str());
            if (mYtitle != "") drawingstack[iplt].second->GetYaxis()->SetTitle(mYtitle.c_str());
            drawingstack[iplt].second->Draw();
            iplt++;
        }
    }
    c1->Update();
    */
}

void MultiPlotter::Draw(std::string opt, std::string key){
    if (key.compare("") != 0 && Exists(key)) {
        mHistos[key]->Draw(opt.c_str());
        return;
    }
    int ndraw=100000; 
    int noffset=0;
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
    if (mUseDefaultLegend) mLeg = new TLegend(xlo,ylo,0.99,0.99);
    gStyle->SetOptStat(0);

    //check for an X axis divisible by pi, if so make the axis units of pi
    IsDivisibleByPi();

    //apply typical gamma spec labels y = counts / kev, x = energy [kev]
    if (mGammaSpecLabels) {
        mXtitle = "Energy [keV]";
        int bw = max->second->GetBinWidth(1);
        if (bw > 1) 
            mYtitle = Form("Counts / %d keV",bw);
        else
            mYtitle = "Counts / keV";
    }

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
        if (mDrawFill) it->second->SetFillColorAlpha(it->second->GetLineColor(),mFillAlpha);
        if (mUseDefaultLegend) mLeg->AddEntry(it->second,it->second->GetName(),"l");

        if (mXhi != mXlo) it->second->GetXaxis()->SetRangeUser(mXlo,mXhi);
        if (mYhi != mYlo) it->second->GetYaxis()->SetRangeUser(mYlo,mYhi);

        if (mXtitle != "") it->second->GetXaxis()->SetTitle(mXtitle.c_str());
        if (mYtitle != "") it->second->GetYaxis()->SetTitle(mYtitle.c_str());
        it->second->GetXaxis()->SetLabelSize(mXLabelSize);
        it->second->GetYaxis()->SetLabelSize(mYLabelSize);
        it->second->GetXaxis()->SetTitleOffset(mXtitleOffset);
        it->second->GetYaxis()->SetTitleOffset(mYtitleOffset);
        it->second->GetXaxis()->SetTitleSize(mXtitleSize);
        it->second->GetYaxis()->SetTitleSize(mYtitleSize);
        
        if (nloops == 0){
            it->second->Draw(opt.c_str());
            it = mHistos.begin();
            for (int j=0; j < noffset; j++) it++;
            nloops++;
        }
        else{
            std::string tempOpt = opt + "same";
            it->second->Draw(tempOpt.c_str());
            it++;
            nloops++;
        }
    }

    //draw peak labels if available
    if (mPeaksToLabel.size() > 0){
        int PTL = mPeaksToLabel.size();
        for (int i=0; i < PTL; i++){
            double X = std::stoi(mPeaksToLabel[i]);
            double Y = max->second->GetBinContent(max->second->FindBin(X)) + 0.2*max->second->GetYaxis()->GetXmax();
            TText *tt = new TText(X,Y,mPeaksToLabel[i].c_str());
            tt->SetTextAngle(90);
            tt->SetTextFont(42);
            tt->SetTextSize(0.035);
            tt->Draw("same");
        }
    }

    mLeg->Draw("same");
}

void MultiPlotter::Rebin(int bg, bool bwScale){
    std::map<std::string, TH1*>::iterator it = mHistos.begin();
    std::map<std::string, TH1*>::iterator end = mHistos.end();
    while (it != end){
        it->second->Rebin(bg);
        if (bwScale) it->second->Scale(1.0/bg);
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

//============================================================================================
//PRIVATE FUNCTIONS
//============================================================================================

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

bool MultiPlotter::IsDivisibleByPi(){
    std::map<std::string, TH1*>::iterator it = mHistos.begin();
    double fraction = (it->second->GetXaxis()->GetXmax()-it->second->GetXaxis()->GetXmin())/TMath::Pi();
    if ( fraction - std::floor(fraction) > 0.0001) return false;

    //pick the correct number of pi units 
    int nPi = (int) std::floor(fraction);
    if (fraction - std::floor(fraction) > std::ceil(fraction) - fraction) 
        nPi = (int) std::floor(fraction);

    //do we want divisions of pi/2
    int ndivisions = nPi*2 + nPi*4*100;
    
    // loop over and set all labels to be in units of pi
    std::map<std::string, TH1*>::iterator end = mHistos.end();
    while (it != end){
        TAxis *a = it->second->GetXaxis();
        a->SetNdivisions(-1*ndivisions);
        a->SetLabelOffset(0.02);
        int nlabel = 2*nPi+1;
        for (int l=2; l <= nlabel; l++){
            int num = l-1;
            int denom = 2;
            if (num%2 == 0) {
                num = num/2;
                denom = 0;
            }
            if (num == 1) {
                if (denom != 0) 
                    a->ChangeLabel(l,-1,-1,-1,-1,-1,Form("#frac{#pi}{%d}",denom));
                else 
                    a->ChangeLabel(l,-1,-1,-1,-1,-1,"#pi");
            }
            else if (denom == 0) a->ChangeLabel(l,-1,-1,-1,-1,-1,Form("%d#pi",num));
            else a->ChangeLabel(l,-1,-1,-1,-1,-1,Form("#frac{%d#pi}{%d}",num,denom));
        }
        it++;
    }
    return true;
}

void MultiPlotter::CanvasPartition(TCanvas *C, const Int_t Nx, const Int_t Ny, Float_t lMargin, Float_t rMargin, Float_t bMargin, Float_t tMargin){
   if (!C)
      return;
 
   // Setup Pad layout:
   Float_t vSpacing = 0.0;
   Float_t vStep = (1. - bMargin - tMargin - (Ny - 1) * vSpacing) / Ny;
 
   Float_t hSpacing = 0.0;
   Float_t hStep = (1. - lMargin - rMargin - (Nx - 1) * hSpacing) / Nx;
 
   Float_t vposd, vposu, vmard, vmaru, vfactor;
   Float_t hposl, hposr, hmarl, hmarr, hfactor;
 
   for (Int_t i = 0; i < Nx; i++) {
 
      if (i == 0) {
         hposl = 0.0;
         hposr = lMargin + hStep;
         hfactor = hposr - hposl;
         hmarl = lMargin / hfactor;
         hmarr = 0.0;
      } else if (i == Nx - 1) {
         hposl = hposr + hSpacing;
         hposr = hposl + hStep + rMargin;
         hfactor = hposr - hposl;
         hmarl = 0.0;
         hmarr = rMargin / (hposr - hposl);
      } else {
         hposl = hposr + hSpacing;
         hposr = hposl + hStep;
         hfactor = hposr - hposl;
         hmarl = 0.0;
         hmarr = 0.0;
      }
 
      for (Int_t j = 0; j < Ny; j++) {
 
         if (j == 0) {
            vposd = 0.0;
            vposu = bMargin + vStep;
            vfactor = vposu - vposd;
            vmard = bMargin / vfactor;
            vmaru = 0.0;
         } else if (j == Ny - 1) {
            vposd = vposu + vSpacing;
            vposu = vposd + vStep + tMargin;
            vfactor = vposu - vposd;
            vmard = 0.0;
            vmaru = tMargin / (vposu - vposd);
         } else {
            vposd = vposu + vSpacing;
            vposu = vposd + vStep;
            vfactor = vposu - vposd;
            vmard = 0.0;
            vmaru = 0.0;
         }
 
         C->cd(0);
 
         auto name = TString::Format("pad_%d_%d", i, j);
         auto pad = (TPad *)C->FindObject(name.Data());
         if (pad)
            delete pad;
         pad = new TPad(name.Data(), "", hposl, vposd, hposr, vposu);
         pad->SetLeftMargin(hmarl);
         pad->SetRightMargin(hmarr);
         pad->SetBottomMargin(vmard);
         pad->SetTopMargin(vmaru);
 
         pad->SetFrameBorderMode(0);
         pad->SetBorderMode(0);
         pad->SetBorderSize(0);
 
         pad->Draw();
      }
   }
}