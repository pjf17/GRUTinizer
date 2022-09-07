
#include "TRuntimeObjects.h"

#include <iostream>
#include <map>
#include <cstdio>

#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TRandom.h>
#include <TObject.h>
#include <TLine.h>

#include "TGretina.h"
#include "TS800.h"
#include "TBank29.h"
#include "TS800.h"
#include "GCutG.h"

#include "TChannel.h"
#include "GValue.h"

// extern "C" is needed to prevent name mangling.
// The function signature must be exactly as shown here,
//   or else bad things will happen.
extern "C"
void MakeHistograms(TRuntimeObjects& obj) {
  TS800    *s800    = obj.GetDetector<TS800>();
  TList    *list    = &(obj.GetObjects());
  int numobj = list->GetSize();

  if (s800){
    double crdc_1_x = s800->GetCrdc(0).GetDispersiveX();
    double crdc_2_x = s800->GetCrdc(1).GetDispersiveX();
    
    double crdc_1_y = s800->GetCrdc(0).GetNonDispersiveY();
    double crdc_2_y = s800->GetCrdc(1).GetNonDispersiveY();
    
    double crdc1_slope = GValue::Value("CRDC1_Y_SLOPE");
    double crdc2_slope = GValue::Value("CRDC2_Y_SLOPE");
    double yhigh = 1200;
    double ylow = 0;
    if (! (std::isnan(crdc1_slope) || std::isnan(crdc2_slope))){
      yhigh = 600;
      ylow = -600;
    }
    obj.FillHistogram("ungated", "crdc1 X_Y", 600, -300, 300, crdc_1_x, 1200, ylow, yhigh, crdc_1_y);  
    obj.FillHistogram("ungated", "crdc2 X_Y", 600, -300, 300, crdc_2_x, 1200, ylow, yhigh, crdc_2_y);

    if(numobj!=list->GetSize()){
      list->Sort();
    }
  }
}
