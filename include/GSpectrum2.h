#ifndef GSPECTRUM2_H
#define GSPECTRUM2_H

#include <cstdio>
#include <Rtypes.h>

class TH1;
class TH2;

#include <TNamed.h>


class GSpectrum2:public TNamed {
  public:
    GSpectrum2(TH2 *mat=0);
    virtual ~GSpectrum2() {} 

    virtual void Print(Option_t *opt="") const { printf("blah\n"); }

    void SetTH2(TH2* mat) { fMat=mat;}
    //void DrawBox(int py=100);


    //repurposed root functions.
    void Draw(Option_t *opt="");            // adds obj to gPad.
    void Paint(Option_t *opt="");           // actually draws the obj. 
    int DistancetoPrimitive(int px,int py); // allows mouse to intact with the obj.
    void ExecuteEvent(Int_t event,Int_t px, Int_t py);


  private:
    TH2 *fMat;
    TH1 *fProj;
    TH1 *fBg;

  ClassDef(GSpectrum2,1)
};

#endif


