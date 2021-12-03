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
//  InitMap();
  TGretina *gretina = obj.GetDetector<TGretina>();
  TList    *list    = &(obj.GetObjects());
  int numobj = list->GetSize();
  
  if (gretina){
    int nHits = gretina->Size();
    for (int i=0; i < nHits; i++){
      TGretinaHit &hit = gretina->GetGretinaHit(i);
      double phi = hit.GetPhiDeg();
      double theta = 180 - hit.GetThetaDeg();
      int cryID = hit.GetCrystalId();
      obj.FillHistogram("crystals",Form("crystal%d",cryID),180,0,180,theta,360,0,360,phi);
      obj.FillHistogram("crystals",Form("ring%02d",hit.GetRingNumber()),180,0,180,theta,360,0,360,phi);
      obj.FillHistogram("crystals","total",180,0,180,theta,360,0,360,phi);
    }
  }

  if(numobj!=list->GetSize()){
    list->Sort();
  }
}
