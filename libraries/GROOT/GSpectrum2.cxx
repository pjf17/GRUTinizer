
#include<GSpectrum2.h>

#include<TROOT.h>
#include<TPad.h>
#include<TRandom.h>
#include<TH2.h>

#include<TVirtualHistPainter.h>

GSpectrum2::GSpectrum2(TH2 *mat):TNamed("GSpectrum2",""),fMat(mat),fProj(0),fBg(0) { 
  if(fMat) {
    this->SetNameTitle(fMat->GetName(),fMat->GetTitle());
  } else {
    fMat = new TH2D("mat","mat",200,0,200,200,0,200);
    for(int i=0;i<50000;i++){
      fMat->Fill(gRandom->Gaus(100,5),gRandom->Gaus(100,5));
    }
    Draw();
  }
}



//void GSpectrum2::DrawBox(int py) { 
// 
//  Int_t nbins = 10;
//  gPad->SetDoubleBuffer(0); // turn off double buffer mode
//  gVirtualX->SetDrawMode(TVirtualX::kInvert); // set the drawing mode to XOR mode
// 
//  static int pyold1 = 0;
//  static int pyold2 = 0;
// 
//  float uxmin = gPad->GetUxmin();
//  float uxmax = gPad->GetUxmax();
//  int pxmin   = gPad->XtoAbsPixel(uxmin);
//  int pxmax   = gPad->XtoAbsPixel(uxmax);
//  Float_t upy = gPad->AbsPixeltoY(py);
//  Float_t y   = gPad->PadtoY(upy);
//  Int_t biny1 = fMat->GetYaxis()->FindBin(y);
//  Int_t biny2 = TMath::Min(biny1+nbins-1, fMat->GetYaxis()->GetNbins());
//  Int_t py1   = gPad->YtoAbsPixel(fMat->GetYaxis()->GetBinLowEdge(biny1));
//  Int_t py2   = gPad->YtoAbsPixel(fMat->GetYaxis()->GetBinUpEdge(biny2));
// 
//  if(pyold1 || pyold2) 
//    gVirtualX->DrawBox(pxmin,pyold1,pxmax,pyold2,TVirtualX::kFilled);
//  gVirtualX->DrawBox(pxmin,py1,pxmax,py2,TVirtualX::kFilled);
//  pyold1 = py1;
//  pyold2 = py2;
//}


void GSpectrum2::Draw(Option_t *opt) {
  printf("%s called.\n",__PRETTY_FUNCTION__);
  TObject *obj=0;
  if(fMat)
    obj = (TObject*)fMat;
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

void GSpectrum2::Paint(Option_t *opt) {
  printf("%s called.\n",__PRETTY_FUNCTION__);
  if(fMat) {
    fMat->Paint("colz");
  }
}

void GSpectrum2::ExecuteEvent(Int_t event,Int_t px,Int_t py) {
  printf("%s called.\n",__PRETTY_FUNCTION__);
  printf("event = 0x%08x\t px = 0x%08x\t py = 0x%08x \n",event,px,py);
  if(fMat && fMat->GetPainter()) {
    fMat->GetPainter()->ExecuteEvent(event,px,py);
  }
}

int GSpectrum2::DistancetoPrimitive(int px,int py) {
  if(fMat && fMat->GetPainter()) {
    return fMat->GetPainter()->DistancetoPrimitive(px,py);
  }
  return 9999;
}
  




