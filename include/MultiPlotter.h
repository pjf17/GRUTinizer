#include <map>
#include <vector>
#include <string>

#include "TFile.h"
#include "TH1.h"
#include "TF1.h"

class MultiPlotter{
    private:
        //variables
        int mNHistos = 0;
        std::map<std::string, TH1*> mHistos;

        double mYMax = 0.0;
        std::string mMaxKey = "";
        
        //drawing settings
        int mColors[12] = {kBlack,kRed,kBlue,kGreen,kCyan,kOrange,kViolet,kGray,kYellow+3,kCyan+3,kBlue-3,kRed+3};
        bool mCustomColors = false;
        int mLineWidth = 1;
        double mXlo = -123;
        double mXhi = -123;

        //functions
        void SortMax();
        bool ParseInput(char *input, std::vector<int> &nums);
        void doSetLineWidth();
        bool Exists(std::string key);

    public:
        void Add(TH1* pHist);
        void Add(TFile *f, const char *hname);
        void Add(TDirectoryFile *f);

        void Clear();
        void Erase(std::string key);
        void List();

        TH1 *GetClone(std::string key);
        TH1 *Get(std::string key);
        
        void Rebin(int bg=2);
        void Scroll(std::string control=" ");
        void SetLineWidth(int w);
        void SetLineColor(std::string key, int c);
        void SetRange(double xlo, double xhi);
        void ResetRange();
        void IterateLineStyle();

        void Norm(std::string mode="area");

        void Fit(std::string key, TF1 *f, double xlo, double xhi);
        void Fit(TF1 *f, double xlo, double xhi);
        void FitPearson(double xlo, double xhi, double h, double c, double w);
        void FitGaus(double xlo, double xhi, Option_t *opt="");
        void FitExclusion(double exlo, double exhi, double rlo, double rhi);
        
        void Draw(std::string key);
        void Draw(int ndraw=100000, int noffset=0);
};