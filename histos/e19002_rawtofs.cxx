
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

  if (!s800){
    return;
  }

  //get the raw TOFS
  // double raw_obj = s800->GetRawOBJ_MESY();
  // double raw_e1 = s800->GetRawE1_MESY();
  // double raw_xf = s800->GetRawXF_MESY();
  
  obj.FillHistogram("ungated", "MTOF_OBJE1", 1000, -10000, 0, s800->GetOBJ_E1Raw());
  obj.FillHistogram("ungated", "MTOF_XFE1", 1000, -6000, 6000, s800->GetXF_E1Raw());
  obj.FillHistogram("ungated", "MTOF_XFE1_vs_MTOF_OBJE1", 2000, 1000, 3000, s800->GetXF_E1Raw(), 2000, -3000,-1000, s800->GetOBJ_E1Raw());
  
  if(numobj!=list->GetSize()){
    list->Sort();
  }
}
