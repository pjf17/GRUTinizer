#ifndef MULTIPLOTTER__H
#define MULTIPLOTTER__H

#include <map>
#include <vector>
#include <string>

#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "THStack.h"
#include "TLegend.h"

class MultiPlotter{
    private:
        //variables
        int mNHistos = 0;
        std::map<std::string, TH1*> mHistos;
        std::map<std::string, TH1*> mHistos_bak;
        TLegend *mLeg = 0;
        bool mUseDefaultLegend = true;

        double mYMax = 0.0;
        std::string mMaxKey = "";
        
        //drawing settings
        int mColors[12] = {kBlack,kRed,kBlue,kGreen,kCyan,kOrange,kViolet,kGray,kYellow+3,kCyan+3,kBlue-3,kRed+3};
        bool mCustomColors = false;
        bool mDrawFill = false;
        bool mGammaSpecLabels = false;
        double mFillAlpha;
        int mLineWidth = 1;
        double mXlo = -123;
        double mXhi = -123;
        double mYlo = -123;
        double mYhi = -123;
        double mXLabelSize = 0.04;
        double mYLabelSize = 0.04;
        std::string mXtitle = "";
        std::string mYtitle = "";
        double mXtitleSize = 0.045;
        double mYtitleSize = 0.045;
        double mXtitleOffset = 1.0;
        double mYtitleOffset = 1.1;

        std::vector<std::string> mPeaksToLabel;

        //functions
        void SortMax();
        bool ParseInput(char *input, std::vector<int> &nums);
        void doSetLineWidth();
        bool Exists(std::string key);
        bool IsDivisibleByPi();
        void CanvasPartition(TCanvas *C, const Int_t Nx = 2, const Int_t Ny = 2, Float_t lMargin = 0.15, Float_t rMargin = 0.05,
            Float_t bMargin = 0.15, Float_t tMargin = 0.05);

    public:
        void Add(TH1* pHist);
        void Add(TFile *f, const char *hname);
        void Add(TDirectoryFile *f);
        void Add(TDirectoryFile *f1,TDirectoryFile *f2,TDirectoryFile *f3=nullptr,
                 TDirectoryFile *f4=nullptr,TDirectoryFile *f5=nullptr) {
                    Add(f1); Add(f2); if(f3) Add(f3); if(f4) Add(f4); if(f5) Add(f5);
                 }
                 
        void Clear();
        THStack *CreateStack();
        void Erase(std::string key);
        void List();
        void ListBackup();

        TH1 *GetClone(std::string key);
        TH1 *Get(std::string key);
        TH1 *GetBackup(std::string key);

        void LabelPeaks(std::vector<std::string> peaks) {mPeaksToLabel = peaks;}
        
        void Rebin(int bw=2, bool bwScale=false);
        void Scroll(std::string control=" ");
        void SetLineWidth(int w);
        void SetLineColor(std::string key, int c);
        void SetFill(double alpha=0.5);
        void SetXrange(double xlo, double xhi);
        void SetYrange(double ylo, double yhi);

        void SetXtitle(std::string title) {mXtitle = title;}
        void SetXtitleSize(double tsize) {mXtitleSize = tsize;}
        void SetXtitleOffset(double toffset) {mXtitleOffset = toffset;}
        void SetXLabelSize(double lsize) {mXLabelSize = lsize;} 
        
        void SetYtitle(std::string title) {mYtitle = title;}
        void SetYtitleSize(double tsize) {mYtitleSize = tsize;}
        void SetYtitleOffset(double toffset) {mYtitleOffset = toffset;}
        void SetYLabelSize(double lsize) {mYLabelSize = lsize;} 

        void SetGammaSpecLabels(bool setting=true) {mGammaSpecLabels = setting;}

        void SetLegendEntry(std::string key, std::string label, std::string opt="l");
        void SetLegendCoordinates(double x1, double y1, double x2, double y2);
        void ResetRange();
        void IterateLineStyle();
        void Integral(double lo=-1, double hi=-1);

        void Norm(std::string mode="area");
        void Norm(double lo, double hi);

        void RatioToHist(std::string key);

        void Fit(std::string key, TF1 *f, double xlo, double xhi);
        void Fit(TF1 *f, double xlo, double xhi);
        void FitPearson(double xlo, double xhi, double h, double c, double w);
        void FitGaus(double xlo, double xhi, Option_t *opt="");
        void FitPeak(double xlo, double xhi, Option_t *opt="");
        void FitExclusion(double exlo, double exhi, double rlo, double rhi);
        
        void Draw(std::string opt="hist",std::string key="");
        void Draw(int nx, int ny, int wx, int wy, bool mergeX=false, bool mergeY=false);
        
        void Add(std::string key, double scale=1.0);
        // void Draw(int ndraw=100000, int noffset=0);

    protected:
        ClassDef(MultiPlotter,4);
};

#endif