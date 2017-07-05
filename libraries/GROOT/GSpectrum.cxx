
#include<GSpectrum.h>

#include<TROOT.h>
#include<TPad.h>
#include<TRandom.h>
#include<TH2.h>

#include<TVirtualHistPainter.h>

#ifndef kButton1Ctrl
#define kButton1Ctrl 9
#define kButton1CtrlMotion 10
#endif

#ifndef kButton1Alt
#define kButton1Alt 15
#define kButton1AltMotion 16
#endif

GSpectrum::GSpectrum(TH1 *mat):TNamed("GSpectrum",""),fSelf(mat) { 
  if(fSelf) {
    this->SetNameTitle(fSelf->GetName(),fSelf->GetTitle());
  } else {
    fSelf = new TH2D("mat","mat",200,0,200,200,0,200);
    for(int i=0;i<50000;i++){
      fSelf->Fill(gRandom->Gaus(100,5),gRandom->Gaus(100,5));
    }
    Draw();
  }
}

void GSpectrum::Draw(Option_t *opt) {
  printf("%s called.\n",__PRETTY_FUNCTION__);
  TObject *obj=0;
  if(fSelf)
    obj = (TObject*)fSelf;
  else
    return;
  //TString sopt =opt; sopt.ToLower();  // check for same arg
  if(gPad) {
    if(!gPad->IsEditable())
      gROOT->MakeDefCanvas();
    if(obj->TestBit(kCanDelete)) gPad->GetListOfPrimitives()->Remove(this);
    gPad->Clear();
  } else {
     
  }
  AppendPad(opt);
}

void GSpectrum::Paint(Option_t *opt) {
  //printf("%s called.\n",__PRETTY_FUNCTION__);
  if(fSelf) {
    fSelf->Paint("colz");
  }
}


int GSpectrum::DistancetoPrimitive(int px,int py) {
  if(fSelf && fSelf->GetPainter()) {
    return fSelf->GetPainter()->DistancetoPrimitive(px,py);
  }
  return 9999;
}
  

void GSpectrum::ExecuteEvent(Int_t event,Int_t px,Int_t py) {
  bool used = false;
  switch(event) {  
    case kButton1Down: //single left click
      printf("left click\n");
      break;
    case kButton1Shift: //single shift-left click
      printf("left shift click\n");
      break;
    case kButton1Ctrl: //single ctrl-left click
      printf("left ctrl click\n");
      //HandleMousePress(event,x,y);
      break;
    case kButton1Alt: //single alt-left click
      printf("left alt click\n");
      //HandleMousePress(event,x,y);
      break;
    default:
      break;
  }
  printf("event = 0x%08x\t px = 0x%08x\t py = 0x%08x \n\n",event,px,py);
  if(!used && fSelf && fSelf->GetPainter()) {
    fSelf->GetPainter()->ExecuteEvent(event,px,py);
  }

}

//bool GSpectrum::Handle


















