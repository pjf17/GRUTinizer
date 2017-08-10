#include <TGater.h>
#include <Globals.h>

#include<GCanvas.h>
#include<TROOT.h>

#include <TGClient.h>
#include <TGButtonGroup.h>
#include <TGLabel.h>
#include <TGSlider.h>
#include <TGDoubleSlider.h>

TGater::TGater():fMatrix_total(0),fTimeRandom_total(0),fAddGate(0),fSubGate(0),fTimeRandScale(1.0),
                 fAddLimits(sqrt(-1),sqrt(-1)),fSubLimits(sqrt(-1),sqrt(-1)),fTotalCanvas(0) {
                   
                   
  Construct();                   
}                   
                   
                   
                   
                   
void TGater::Construct() {                   
                   
  fMain = new TGMainFrame(gClient->GetRoot(),500,500);
  
  TGVerticalFrame *vframe = new TGVerticalFrame(fMain,250,600);
 

//////////////////////
  TGHorizontalFrame *scaleframe = new TGHorizontalFrame(vframe,200,200);
  
  TGLabel *scale_label = new TGLabel(scaleframe,"Time Rand Scale");
  scaleframe->AddFrame(scale_label,new TGLayoutHints(kLHintsExpandX|kLHintsExpandY,1,1,1,1));

  fTRScale = new TGNumberEntry(scaleframe,1.0,5,999,
                                      TGNumberFormat::kNESRealTwo,
                                      TGNumberFormat::kNEAPositive,
                                      TGNumberFormat::kNELNoLimits,
                                      0,99999);
  fTRScale->Connect("ValueSet(Long_t)","TGater",this,"SetTimeRandScale()");
  fTRScale->SetText("TimeRandomScale");
  scaleframe->AddFrame(fTRScale,new TGLayoutHints(kLHintsExpandX|kLHintsExpandY,1,1,1,1));

  vframe->AddFrame(scaleframe,new TGLayoutHints(kLHintsExpandX|kLHintsExpandY,1,1,1,1));
//////////////////////


//////////////////////
  
  TGVerticalFrame *addframe = new TGVerticalFrame(vframe,200,400);
  TGHorizontalFrame *addrow1 = new TGHorizontalFrame(addframe,200,200);
  TGHorizontalFrame *addrow2 = new TGHorizontalFrame(addframe,200,200);
  
  TGLabel *add_label1 = new TGLabel(addrow1,"Add");
  addrow1->AddFrame(add_label1,new TGLayoutHints(kLHintsExpandX|kLHintsExpandY,1,1,1,1));

  fAddLow = new TGNumberEntry(addrow1,1.0,5,999,
                                      TGNumberFormat::kNESRealTwo,
                                      TGNumberFormat::kNEAPositive,
                                      TGNumberFormat::kNELNoLimits,
                                      0,99999);
  fAddLow->Connect("ValueSet(Long_t)","TGater",this,"CheckAddGate()");
  addrow1->AddFrame(fAddLow,new TGLayoutHints(kLHintsExpandX|kLHintsExpandY,1,1,1,1));
  fAddHigh = new TGNumberEntry(addrow1,1.0,5,999,
                                      TGNumberFormat::kNESRealTwo,
                                      TGNumberFormat::kNEAPositive,
                                      TGNumberFormat::kNELNoLimits,
                                      0,99999);
  fAddHigh->Connect("ValueSet(Long_t)","TGater",this,"CheckAddGate()");
  addrow1->AddFrame(fAddHigh,new TGLayoutHints(kLHintsExpandX|kLHintsExpandY,1,1,1,1));
   
  //TGLabel *add_label2 = new TGLabel(addrow2,"   ");
  //addrow2->AddFrame(add_label2,new TGLayoutHints(kLHintsExpandX|kLHintsExpandY,1,1,1,1));
  fLiveUpdate = new TGCheckButton(addrow2,"Live Projections",2000);
  addrow2->AddFrame(fLiveUpdate,new TGLayoutHints(kLHintsExpandX|kLHintsExpandY,1,1,1,1));


  fAddSlider = new TGDoubleHSlider(addrow2,100);
  fAddSlider->Connect("PositionChanged()","TGater",this,"CheckAddGate()");
  addrow2->AddFrame(fAddSlider,new TGLayoutHints(kLHintsExpandX|kLHintsExpandY,1,1,1,1));
  
  addframe->AddFrame(addrow1,new TGLayoutHints(kLHintsExpandX|kLHintsExpandY,1,1,1,1));
  addframe->AddFrame(addrow2,new TGLayoutHints(kLHintsExpandX|kLHintsExpandY,1,1,1,1));

  vframe->AddFrame(addframe,new TGLayoutHints(kLHintsExpandX|kLHintsExpandY,1,1,1,1));
//////////////////////
  

//////////////////////
  TGVerticalFrame *subframe = new TGVerticalFrame(vframe,200,600);
  TGHorizontalFrame *subrow1 = new TGHorizontalFrame(subframe,200,200);
  TGHorizontalFrame *subrow2 = new TGHorizontalFrame(subframe,200,200);
  TGHorizontalFrame *subrow3 = new TGHorizontalFrame(subframe,200,200);
  
  TGLabel *sub_label1 = new TGLabel(subrow1,"Sub");
  subrow1->AddFrame(sub_label1,new TGLayoutHints(kLHintsExpandX|kLHintsExpandY,1,1,1,1));

  fSubLow = new TGNumberEntry(subrow1,1.0,5,999,
                                      TGNumberFormat::kNESRealTwo,
                                      TGNumberFormat::kNEAPositive,
                                      TGNumberFormat::kNELNoLimits,
                                      0,99999);
  fSubLow->Connect("ValueSet(Long_t)","TGater",this,"CheckSubGate()");
  subrow1->AddFrame(fSubLow,new TGLayoutHints(kLHintsExpandX|kLHintsExpandY,1,1,1,1));
  fSubHigh = new TGNumberEntry(subrow1,1.0,5,999,
                                      TGNumberFormat::kNESRealTwo,
                                      TGNumberFormat::kNEAPositive,
                                      TGNumberFormat::kNELNoLimits,
                                      0,99999);
  fSubHigh->Connect("ValueSet(Long_t)","TGater",this,"CheckSubGate()");
  subrow1->AddFrame(fSubHigh,new TGLayoutHints(kLHintsExpandX|kLHintsExpandY,1,1,1,1));
  
  //TGLabel *sub_label2 = new TGLabel(subrow2,"   ");
  //subrow2->AddFrame(sub_label2,new TGLayoutHints(kLHintsExpandX|kLHintsExpandY,1,1,1,1));
  TGHButtonGroup *fSubOpts = new TGHButtonGroup(subrow2,"sub opts");
  TGRadioButton *sub_bgb[3];
  sub_bgb[0] = new TGRadioButton(fSubOpts,new TGHotString("Ignore"),1000);
  sub_bgb[1] = new TGRadioButton(fSubOpts,new TGHotString("Subtract"),1001);
  sub_bgb[1]->SetState(kButtonDown);
  fSubGateOptions = 1001;
  sub_bgb[2] = new TGRadioButton(fSubOpts,new TGHotString("Overlay"),1002);
  fSubOpts->Connect("Clicked(Int_t)","TGater",this,"SetSubtractOption(int)");
  fSubOpts->Show();
  subrow2->AddFrame(fSubOpts,new TGLayoutHints(kLHintsExpandX|kLHintsExpandY,1,1,1,1));

  fSubSlider = new TGDoubleHSlider(subrow2,100);
  fSubSlider->Connect("PositionChanged()","TGater",this,"CheckSubGate()");
  subrow2->AddFrame(fSubSlider,new TGLayoutHints(kLHintsExpandX|kLHintsExpandY,1,1,1,1));
  
  TGHButtonGroup *fSubScaling = new TGHButtonGroup(subrow3,"sub scaling");
  TGRadioButton *sub_scale_bgb[2];
  sub_scale_bgb[0] = new TGRadioButton(fSubScaling,new TGHotString("Auto"),3000);
  sub_scale_bgb[1] = new TGRadioButton(fSubScaling,new TGHotString("Scale"),3001);
  sub_scale_bgb[1]->SetState(kButtonDown);
  fSubGateScaleOption = 3001;
  fSubScaling->Connect("Clicked(Int_t)","TGater",this,"SetSubtractScaling(int)");
  fSubScaling->Show();
  subrow3->AddFrame(fSubScaling,new TGLayoutHints(kLHintsExpandX|kLHintsExpandY,1,1,1,1));
  
  fSubScale = new TGNumberEntry(subrow3,1.0,5,999,
                                      TGNumberFormat::kNESRealTwo,
                                      TGNumberFormat::kNEAPositive,
                                      TGNumberFormat::kNELNoLimits,
                                      0,99999);
  fSubScale->Connect("ValueSet(Long_t)","TGater",this,"DrawProj()");
  subrow3->AddFrame(fSubScale,new TGLayoutHints(kLHintsExpandX|kLHintsExpandY,1,1,1,1));
  
  TGLabel *sub_label3 = new TGLabel(subrow3," ");
  subrow3->AddFrame(sub_label3,new TGLayoutHints(kLHintsExpandX|kLHintsExpandY,1,1,1,1));

  
  
  subframe->AddFrame(subrow1,new TGLayoutHints(kLHintsExpandX|kLHintsExpandY,1,1,1,1));
  subframe->AddFrame(subrow2,new TGLayoutHints(kLHintsExpandX|kLHintsExpandY,1,1,1,1));
  subframe->AddFrame(subrow3,new TGLayoutHints(kLHintsExpandX|kLHintsExpandY,1,1,1,1));

  vframe->AddFrame(subframe,new TGLayoutHints(kLHintsExpandX|kLHintsExpandY,1,1,1,1));
//////////////////////
  

/////////////////////
  TGHorizontalFrame *buttonframe = new TGHorizontalFrame(vframe,200,200);

  fPrintButton = new TGTextButton(buttonframe,"Print");
  fPrintButton->Connect("Clicked()","TGater",this,"Print()");
  buttonframe->AddFrame(fPrintButton,new TGLayoutHints(kLHintsExpandX|kLHintsExpandY,1,1,1,1));

  fDrawButton  = new TGTextButton(buttonframe,"Total");
  fDrawButton->Connect("Clicked()","TGater",this,"Draw()");
  buttonframe->AddFrame(fDrawButton,new TGLayoutHints(kLHintsExpandX|kLHintsExpandY,1,1,1,1));

  fProjButton  = new TGTextButton(buttonframe,"Project");
  fProjButton->Connect("Clicked()","TGater",this,"DrawProj()");
  buttonframe->AddFrame(fProjButton,new TGLayoutHints(kLHintsExpandX|kLHintsExpandY,1,1,1,1));

  vframe->AddFrame(buttonframe,new TGLayoutHints(kLHintsExpandX|kLHintsExpandY,1,1,1,1));
////////////////////


  fMain->AddFrame(vframe,new TGLayoutHints(kLHintsExpandX|kLHintsExpandY,1,1,1,1));
  
  fMain->SetWindowName("gate helper");
  fMain->MapSubwindows();
  fMain->Resize(fMain->GetDefaultSize());
  fMain->MapWindow();
}



void TGater::Print(Option_t *opt) const {
  printf(DCYAN "\t\tThe Gate Helper\n" RESET_COLOR);
  printf("--------------------------\n");
  printf("Prompt Matrix:     %s\n",fMatrix.GetName());
  printf("Time Rand Matrix:  %s\n",fTimeRandom.GetName());
  printf("\n");
  printf("Current Add Gate[0x%08x]: %s\n",fAddGate,fAddGate==0?"":fAddGate->GetName());
  printf("Current Sub Gate[0x%08x]: %s\n",fSubGate,fSubGate==0?"":fSubGate->GetName());
  printf("\n");
  printf("TimeRandom Scale:    %.04f\n", fTimeRandScale);
  printf("Add Gate Markers:    %.02f\t - \t%.02f\n",fAddLimits.first,fAddLimits.second);
  printf("Sub Gate Markers:    %.02f\t - \t%.02f\n",fSubLimits.first,fSubLimits.second);
  printf("--------------------------\n");

}

void TGater::Clear(Option_t *opt) {
  fMatrix_total     = 0;
  fTimeRandom_total = 0;
  fAddGate          = 0;
  fSubGate          = 0;

  fTotalCanvas           = 0;
  fProjectCanvas         = 0;
}


void TGater::Draw(Option_t *opt) { 
  printf("scribble, scribble.\n");
  ConstructCanvas();
  //fTotalCanvas->cd(1);
  if(!fMatrix_total)
    fMatrix_total = fMatrix.ProjectionX();  
  printf("name: %s\n",fMatrix_total->GetName());
  fMatrix_total->Draw();
  //fAddGate = (TH1D*)fMatrix_total->Clone("AddGate");
  //fAddGate->Draw();
  //fTotalCanvas->cd(1);

  ResetGates();
}

void TGater::DrawProj(Option_t *opt) {
  printf("quick, get the crayons!\n");
  if(!fProjectCanvas || !gROOT->GetListOfCanvases()->FindObject(fProjectCanvas)) {
    gROOT->MakeDefCanvas();
    fProjectCanvas   = (GCanvas*)gPad->GetCanvas();
    fProjectCanvas->SetTitle("Projections!");
  }
  fProjectCanvas->cd();
  if(GrabHist(0)) {
    fCurrent = (TH1D*)(GrabHist(0)->Clone("clone"));
  }
  if(fAddGate) {fAddGate->Delete();  }
  fAddGate = DoProjection(fAddLimits.first,fAddLimits.second);
  if(fSubGate) {fSubGate->Delete();  }
  fSubGate = DoProjection(fSubLimits.first,fSubLimits.second);
  //printf("fSubGateOptions = %i\n",fSubGateOptions);

  double scale=1.0;
  switch(fSubGateScaleOption) {
    case 3000:
      scale = (fAddLimits.second-fAddLimits.first)/(fSubLimits.second-fSubLimits.first);
      break;
    case 3001:
      scale = fSubScale->GetNumberEntry()->GetNumber();
      break;
  }
 
  fSubGate->Scale(scale);

  switch(fSubGateOptions) {
    case 1000:
      fAddGate->SetLineColor(kBlue);
      fAddGate->Draw();
      break;
    case 1001:
      fAddGate->Add(fSubGate,-1);
      fAddGate->SetLineColor(kBlue);
      fAddGate->Draw();
      break;
    case 1002:
      fAddGate->Draw();
      fAddGate->SetLineColor(kRed);
      fSubGate->Draw("same");
      break;
  }
  if(fCurrent) { 
    fAddGate->GetXaxis()->SetRangeUser(fCurrent->GetXaxis()->GetBinLowEdge(fCurrent->GetXaxis()->GetFirst()),
                                       fCurrent->GetXaxis()->GetBinUpEdge(fCurrent->GetXaxis()->GetLast()));
  }
  fProjectCanvas->Modified();
  fProjectCanvas->Update();

  if(!fTotalCanvas || !gROOT->GetListOfCanvases()->FindObject(fTotalCanvas)) { Draw(); }
  fTotalCanvas->cd();
}


TH1D *TGater::DoProjection(double low,double high) {
  //printf("%s called\t%.1f to %.1f\n",__PRETTY_FUNCTION__,low,high);
  if(!fMatrix_total) {
    fMatrix.GetXaxis()->UnZoom();
    fMatrix.GetYaxis()->UnZoom();
    fMatrix_total = fMatrix.ProjectionX();  
  }

  int bin_low  = fMatrix_total->GetXaxis()->FindBin(low);
  int bin_high = fMatrix_total->GetXaxis()->FindBin(high);
  static int fNp=0;
  TH1D *proj = fMatrix.ProjectionY(Form("project_%i",fNp++),bin_low,bin_high);
  proj->SetTitle(Form("%.1f to %.1f",low,high));
  return proj;
}

void TGater::ResetGates() {
  if(fMatrix_total) {
    double low  = fMatrix_total->GetXaxis()->GetBinLowEdge(fMatrix_total->GetXaxis()->GetFirst());
    double high = fMatrix_total->GetXaxis()->GetBinUpEdge(fMatrix_total->GetXaxis()->GetLast());
    fAddSlider->SetRange(low,high);
    fSubSlider->SetRange(low,high);
    SetAddGate(low,high);
    SetSubGate(low,high);
  }

}

void TGater::ConstructCanvas() {
  if(!fTotalCanvas || !gROOT->GetListOfCanvases()->FindObject(fTotalCanvas)) {
    gROOT->MakeDefCanvas();
    fTotalCanvas   = (GCanvas*)gPad->GetCanvas();
    //gPad->Divide(1,2);
    //gPad->cd(1);
    gPad->Connect("RangeChanged()","TGater",this,"ResetGates()");
    fTotalCanvas->SetTitle("gate helper"); 
    //fTotalCanvas->Connect("Closed()","TGater",this,"Clear()");
    fTotalCanvas->Connect("AddMarker(int,int,int)","TGater",this,"MarkerAdded()");
    fTotalCanvas->Connect("SetBackgroundMarkers()","TGater",this,"GrabMarkers()");
    fTotalCanvas->Connect("RemoveMarker()","TGater",this,"ResetGates()");
  } 
}

void TGater::MarkerAdded() {
  printf("%s called.\n",__PRETTY_FUNCTION__);
}

void TGater::GrabMarkers() {
  if(!fTotalCanvas || !gROOT->GetListOfCanvases()->FindObject(fTotalCanvas)) 
    return;
  if(fTotalCanvas->GetMarkers().size()<2 || fTotalCanvas->GetBackgroundMarkers().size()<2)
    return;
  SetAddGate(fTotalCanvas->GetMarkers().at(0)->localx,fTotalCanvas->GetMarkers().at(1)->localx);
  SetSubGate(fTotalCanvas->GetBackgroundMarkers().at(0)->localx,fTotalCanvas->GetBackgroundMarkers().at(1)->localx);
  //fTotalCanvas->GetBackgroundMarkers().at(0)->localx
  //fTotalCanvas->GetBackgroundMarkers().at(1)->localx
}

void TGater::UpdateMarkers() {
  if(!fTotalCanvas || !gROOT->GetListOfCanvases()->FindObject(fTotalCanvas)) 
    return;
  if(fTotalCanvas->GetMarkers().size()<2 || fTotalCanvas->GetBackgroundMarkers().size()<2)
    return;
  fTotalCanvas->cd();
  double bin_edge;
  fTotalCanvas->GetMarkers().at(0)->localx = fAddLimits.first;
  fTotalCanvas->GetMarkers().at(0)->binx   = fMatrix_total->GetXaxis()->FindBin(fAddLimits.first);
  bin_edge = fMatrix_total->GetXaxis()->GetBinLowEdge(fTotalCanvas->GetMarkers().at(0)->binx);
  fTotalCanvas->GetMarkers().at(0)->linex->SetX1(bin_edge);
  fTotalCanvas->GetMarkers().at(0)->linex->SetX2(bin_edge);
  
  fTotalCanvas->GetMarkers().at(1)->localx = fAddLimits.second;
  fTotalCanvas->GetMarkers().at(1)->binx   = fMatrix_total->GetXaxis()->FindBin(fAddLimits.second);
  bin_edge = fMatrix_total->GetXaxis()->GetBinLowEdge(fTotalCanvas->GetMarkers().at(1)->binx);
  fTotalCanvas->GetMarkers().at(1)->linex->SetX1(bin_edge);
  fTotalCanvas->GetMarkers().at(1)->linex->SetX2(bin_edge);
  
  fTotalCanvas->RedrawMarkers(); 
}

void TGater::UpdateBGMarkers() {
  if(!fTotalCanvas || !gROOT->GetListOfCanvases()->FindObject(fTotalCanvas)) 
    return;
  if(fTotalCanvas->GetMarkers().size()<2 || fTotalCanvas->GetBackgroundMarkers().size()<2)
    return;
  fTotalCanvas->cd();
  double bin_edge;

  fTotalCanvas->GetBackgroundMarkers().at(0)->localx = fSubLimits.first;
  fTotalCanvas->GetBackgroundMarkers().at(0)->binx   = fMatrix_total->GetXaxis()->FindBin(fSubLimits.first);
  bin_edge = fMatrix_total->GetXaxis()->GetBinLowEdge(fTotalCanvas->GetBackgroundMarkers().at(0)->binx);
  fTotalCanvas->GetBackgroundMarkers().at(0)->linex->SetX1(bin_edge);
  fTotalCanvas->GetBackgroundMarkers().at(0)->linex->SetX2(bin_edge);
  
  fTotalCanvas->GetBackgroundMarkers().at(1)->localx = fSubLimits.second;
  fTotalCanvas->GetBackgroundMarkers().at(1)->binx   = fMatrix_total->GetXaxis()->FindBin(fSubLimits.second);
  bin_edge = fMatrix_total->GetXaxis()->GetBinLowEdge(fTotalCanvas->GetBackgroundMarkers().at(1)->binx);
  fTotalCanvas->GetBackgroundMarkers().at(1)->linex->SetX1(bin_edge);
  fTotalCanvas->GetBackgroundMarkers().at(1)->linex->SetX2(bin_edge);
  
  fTotalCanvas->RedrawMarkers(); 
}




void TGater::SetAddGate(double low,double high) { 
  if(low>high) 
    std::swap(low,high);
  fAddLimits = std::make_pair(low,high); 
  
  fAddLow->SetNumber(low);
  fAddHigh->SetNumber(high);
  //fAddSlider->SetRange(low,high);
  fAddSlider->SetPosition(low,high);

  if(fLiveUpdate->IsOn())
    DrawProj();

  UpdateMarkers();
}

void TGater::CheckAddGate() {
  double low  = fAddLimits.first;
  double high = fAddLimits.second;
  double tol = 0.01;

  //printf("in %s\n",__PRETTY_FUNCTION__);
  //printf("     %.02f \t %.02f\n",low,high);
  //printf("NE:  %.02f \t %.02f\n",fAddLow->GetNumberEntry()->GetNumber(),fAddHigh->GetNumberEntry()->GetNumber());

  if((fabs(fAddLow->GetNumberEntry()->GetNumber()-low)>tol) ||
     (fabs(fAddHigh->GetNumberEntry()->GetNumber()-high)>tol)) {
    SetAddGate(fAddLow->GetNumberEntry()->GetNumber(),fAddHigh->GetNumberEntry()->GetNumber());
    //printf("reset from buttons called!\n");
  }
  float slow,shigh;
  fAddSlider->GetPosition(slow,shigh);
  //printf("SB:  %.02f \t %.02f\n",slow,shigh);
  if((fabs(low-slow)>tol) || (fabs(high-shigh)>tol)) {
    SetAddGate(slow,shigh);
    //printf("reset from slider called!\n");
  }
}

void TGater::SetSubGate(double low,double high) { 
  if(low>high) 
    std::swap(low,high);
  fSubLimits = std::make_pair(low,high);
  
  fSubLow->SetNumber(low);
  fSubHigh->SetNumber(high);
  //fAddSlider->SetRange(low,high);
  fSubSlider->SetPosition(low,high);
 
  if(fLiveUpdate->IsOn())
    DrawProj();

  UpdateBGMarkers();
}

void TGater::CheckSubGate() {
  double low  = fSubLimits.first;
  double high = fSubLimits.second;
  double tol = 0.01;

  //printf("in %s\n",__PRETTY_FUNCTION__);

  if((fabs(fSubLow->GetNumberEntry()->GetNumber()-low)>tol) ||
     (fabs(fSubHigh->GetNumberEntry()->GetNumber()-high)>tol)) {
    SetSubGate(fSubLow->GetNumberEntry()->GetNumber(),fSubHigh->GetNumberEntry()->GetNumber());
    //printf("reset from buttons called!\n");
  }
  float slow,shigh;
  fSubSlider->GetPosition(slow,shigh);
  if((fabs(low-slow)>tol) || (fabs(high-shigh)>tol)) {
    SetSubGate(slow,shigh);
    //printf("reset from slider called!\n");
  }
}


void TGater::SetTimeRandScale(double x) {
  if(std::isnan(x))
    fTimeRandScale = x;
  else 
    fTimeRandScale = fTRScale->GetNumberEntry()->GetNumber();
}


void TGater::SetSubtractOption(int opt) {
  //printf("%s  recieved %i\n",__PRETTY_FUNCTION__,opt);
  fSubGateOptions = opt;
  DrawProj();
}

void TGater::SetSubtractScaling(int opt) {
  fSubGateScaleOption = opt;
  DrawProj();
}


