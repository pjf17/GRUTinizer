#ifndef GPolAnalyzer_H
#define GPolAnalyzer_H

#include "TObject.h"
#include "TH1.h"
#include "GH1D.h"
#include "TDirectory.h"
#include "TString.h"

#include <vector>

class GPolAnalyzer {
  private:
    TDirectory *fSrce;
    TDirectory *fBeam;
    TString fMassNum;
    bool fDoRand;
    int fEsrceGate[2];
    int fEbeamGate[2];
    TString getHistDirectories(TDirectory *f);
    void getHistNames(TDirectory *f, TString dir, std::vector<std::pair<int,std::string>> &names, bool randMode);
    void groupCrystals(std::vector<std::pair<int,std::string>> &crstl, TString mode);

  public:
    GPolAnalyzer() { fSrce=0; fBeam=0; fMassNum = "";}
    GPolAnalyzer(TH1 *hsource, TH1 *hbeam, int EsrceLo, int EsrceHi, int EbeamLo, int EbeamHi);
    GPolAnalyzer(TH1 *h, int ELo, int EHi);

    void GetXiRatio(GH1D *&hnorm, GH1D *&hsrce, GH1D *&hbeam, int binning = 1, TString opt = "crys"); 
};

#endif