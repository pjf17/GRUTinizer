#ifndef _GOBJECT_H_
#define _GOBJECT_H_

#include <GuiTypes.h>
#include <Rtypes.h>

class GObject {
  public:
    virtual ~GObject() { }

    virtual bool HandleMouseEvent(Int_t event,Int_t px, Int_t py) = 0; 
    virtual bool HandleKeyEvent(Event_t *event,UInt_t *keysym) = 0; 

  protected: 
    GObject() { }
  ClassDef(GObject,0)
};


#endif
