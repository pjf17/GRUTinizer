/// notes.
//
//  - stange floating point sizes are suppose to % of the thing they are relative to, 
//    howeverIt is not obvious what the relative thing is most of the time.
//  - font 42 =  helvetica-medium-r-normal, scalable. 
//  - font 52 =  helvetica-medium-o-normal, italic, scalable. 
//  - font 62 =  helvetica-bold-r-normal, bold, scalable. 
//  - font 72 =  helvetica-bold-o-normal, bold, italic, scalable. 
//
//  for best results run this before doing anything else or
//  copy every thing between (and including) the curlly brackets to 
//  a file named rootlogon.C in your home directory.

//
// Further magic:
//   hist->SetContour(200), needs to be called directly if this style 
//   isn't loaded when the hist is created.  Calling this makes the 
//   histogram unmeasurably nicer.
//
//
//   I also cannot seem to globally move the color scale.  Setting the 
//   margins below allows for room for the scale to be moved to the
//   right.  To actually move the scale:
//   
//     TPaletteAxis *pal = (TPaletteAxis*)hist->GetListOfFunctions()->FindObject("palette")
//     pal->SetX1NDC(0.86);
//     pal->SetX2NDC(0.91);
//     then resize the canvas a bit for it to take effect.

void NiceMatStyle() {

  //******************//
  //  color gradient  //
  TStyle *style = new TStyle("NiceMat","");
  
  //style->SetPalette(55);
  //  https://root.cern.ch/doc/master/classTColor.html
  //  see TColor::SetPalette
  //
  const Int_t NRGBs = 5;
  const Int_t NCont = 500;
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };

  TColor::CreateGradientColorTable(NRGBs,stops,red,green,blue,NCont);
  style->SetNumberContours(NCont); 
  //
  //******************//

  //misc.
  style->SetOptStat(0);
  style->SetStatTextColor(1);  // 1 = Black
  style->SetDrawBorder(0);
  style->SetFuncColor(kRed);
  style->SetMarkerStyle(20); 
  style->SetOptLogz(1);


  //canvas
  style->SetCanvasBorderSize(0);
  style->SetCanvasBorderMode(0);
  style->SetCanvasBorderMode(0);
  style->SetCanvasColor(0);
  style->SetCanvasDefH(700); //Height of canvas
  style->SetCanvasDefW(900); //Width of canvas
  style->SetCanvasDefX(10);   //POsition on screen
  style->SetCanvasDefY(10);


  //pad
  style->SetPadBorderSize(0);
  style->SetPadBorderMode(0);
  style->SetPadColor(0);
  style->SetPadGridX(false);
  style->SetPadGridY(false);
  style->SetGridColor(0);
  style->SetGridStyle(3);
  style->SetGridWidth(1);

  //frame
  style->SetFrameBorderMode(0);
  style->SetFrameBorderSize(1);
  style->SetFrameFillColor(0);
  style->SetFrameFillStyle(0);
  style->SetFrameLineColor(1);
  style->SetFrameLineStyle(1);
  style->SetFrameLineWidth(1);


  //margins
  //style->SetPadTopMargin(0.05);
  style->SetPadBottomMargin(0.13);
  style->SetPadLeftMargin(0.15);
  style->SetPadRightMargin(0.12);



  //graph title
  style->SetTitleFont(62);
  style->SetTitleColor(1);     // 1 = Black
  style->SetTitleTextColor(1);     
  style->SetTitleFillColor(0); // 0 = White
  style->SetTitleFontSize(0.06);
  style->SetTitleX(0.5);
  //style->SetTitleAlign(13); //align left
  style->SetTitleAlign(23); //align center
  //style->SetTitleAlign(33); //align right
  //style->SetTitleH(0); // Set the height of the title box
  //style->SetTitleW(0); // Set the width of the title box
  //style->SetTitleX(0); // Set the position of the title box
  //style->SetTitleY(0.985); // Set the position of the title box
  //style->SetTitleStyle(Style_t style = 1001);
  style->SetTitleBorderSize(0);


  //axis titles
  style->SetTitleColor(1,"xyz");
  style->SetTitleFont(62,"xyz");
  style->SetTitleSize(0.06,"xyz");
  //style->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  //style->SetTitleYSize(Float_t size = 0.02);
  //style->SetTitleOffset(1.1, "xyz"); // Another way to set the Offset
  //style->SetTitleXOffset(0.9);
  style->SetTitleYOffset(1.25);
  
  //axis labels:
  style->SetLabelColor(1,"xyz");
  style->SetLabelFont(62,"xyz");
  style->SetLabelOffset(0.01,"xyz");
  style->SetLabelSize(0.04,"xyz");

  //the axis:
  style->SetAxisColor(1, "xyz");
  style->SetStripDecimals(true);
  style->SetTickLength(0.03, "xyz");
  //style->SetNdivisions(20, "xyz");
  style->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  style->SetPadTickY(1);


  gROOT->SetStyle("NiceMat");
  gROOT->ForceStyle();


}

