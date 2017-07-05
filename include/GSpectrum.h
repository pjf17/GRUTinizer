#ifndef GSPECTRUM_H
#define GSPECTRUM_H

#include <cstdio>
#include <Rtypes.h>

class TH1;
class TH2;

#include <TNamed.h>
#include <GObject.h>
#include <Buttons.h> // enums for event types.


class GSpectrum:public TNamed,public GObject {
  public:
    GSpectrum(TH1 *mat=0);
    virtual ~GSpectrum() {} 

    virtual void Print(Option_t *opt="") const { printf("blah\n"); }

    //void DrawBox(int py=100);
    
    //repurposed root functions.
    void Draw(Option_t *opt="");            // adds obj to gPad.
    void Paint(Option_t *opt="");           // actually draws the obj. 
    int DistancetoPrimitive(int px,int py); // allows mouse to intact with the obj.
    void ExecuteEvent(Int_t event,Int_t px, Int_t py);

    
    bool HandleMouseEvent(Int_t event,Int_t px, Int_t py) { return false;}
    bool HandleKeyEvent(Event_t *event,UInt_t *keysym)      { return false;}


  private:
    TH1 *fSelf;

  ClassDef(GSpectrum,1)
};

#endif


