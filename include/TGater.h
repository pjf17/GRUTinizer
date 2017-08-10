#ifndef __TGATER_H__
#define __TGATER_H__


#include <cstdio>
#include <utility>
#include <cmath>

#include <TH2.h>
#include <TH1.h>

#include <TGFrame.h>
#include <TGButton.h>
#include <TGNumberEntry.h>

#include <GRootCommands.h>

class GCanvas;
class TGDoubleHSlider;
class TGCheckButton;


class TGater { 
  public: 
    TGater();
    virtual ~TGater() { }
 
    virtual void Print(Option_t *opt="") const;
    virtual void Clear(Option_t *opt="");
    virtual void Draw(Option_t *opt="");
    virtual void DrawProj(Option_t *opt="");
    

    void SetMatrix(TH2 *mat)            { Clear(); mat->Copy(fMatrix);  }
    void SetPrompt(TH2 *prom)           { prom->Copy(fTimeRandom); }
    void SetTimeRand(TH2 *rand)         { rand->Copy(fTimeRandom); }
    void SetTimeRandScale(double x=sqrt(-1));
    void CheckAddGate();
    void CheckSubGate();
    void SetAddGate(double low,double high);
    void SetSubGate(double low,double high);

    void SetSubtractOption(int);
    void SetSubtractScaling(int);
    
    TH1D *DoProjection(double low,double high);

  //private:
    TH2D fMatrix;
    TH2D fPrompt;
    TH2D fTimeRandom;

    TH1D *fMatrix_total;
    TH1D *fTimeRandom_total;
    
    TH1D *fAddGate;
    TH1D *fSubGate;
    TH1D *fCurrent;

    double fTimeRandScale;
    std::pair<double,double> fAddLimits;
    std::pair<double,double> fSubLimits;
    int fSubGateOptions; //1000,1001,1002
    int fSubGateScaleOption; // 3000,3001;

    GCanvas *fTotalCanvas;
    GCanvas *fProjectCanvas;

  //private:
    TGMainFrame   *fMain;
    
    TGButton      *fPrintButton;
    TGButton      *fDrawButton;
    TGButton      *fProjButton;
    TGNumberEntry *fTRScale;
    TGCheckButton *fLiveUpdate;

    TGDoubleHSlider *fAddSlider;
    TGNumberEntry   *fAddLow;
    TGNumberEntry   *fAddHigh;
    
    //TGHButtonGroup  *fSubOpts; 
    TGDoubleHSlider *fSubSlider;
    TGNumberEntry   *fSubLow;
    TGNumberEntry   *fSubHigh;
    TGNumberEntry   *fSubScale;

  //private:
    void Construct();
    void ConstructCanvas();

    void ResetGates();

    void MarkerAdded();
    void GrabMarkers();
    void UpdateMarkers();
    void UpdateBGMarkers();
  
  //private:
    



  ClassDef(TGater,0) 
};


#endif

