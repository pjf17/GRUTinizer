#include "TRuntimeObjects.h"

#include <iostream>
#include <map>
#include <cstdio>
#include <chrono>

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

std::map<int,int> detMap = {
  {24, 0}, {25, 1}, {26, 2}, {27, 3}, {28, 4}, {29, 5}, {30, 6}, {31, 7},
  {32, 8}, {33, 9}, {34,10}, {35,11}, {36,12}, {37,13}, {38,14}, {39,15},
  {44,16}, {45,17}, {46,18}, {47,19}, {48,20}, {49,21}, {50,22}, {51,23},
  {56,24}, {57,25}, {58,26}, {59,27}, {60,28}, {61,29}, {62,30}, {63,31},
  {64,32}, {65,33}, {66,34}, {67,35}, {68,36}, {69,37}, {70,38}, {71,39},
  {76,40}, {77,41}, {78,42}, {79,43}, {80,44}, {81,45}, {82,46}, {83,47}
};

const double TIMESCALE = 1E8; // 10ns
double timeZero = -1;

// extern "C" is needed to prevent name mangling.
// The function signature must be exactly as shown here,
//   or else bad things will happen.
extern "C"
void MakeHistograms(TRuntimeObjects& obj) {
  TGretina *gretina = obj.GetDetector<TGretina>();
  TS800 *s800 = obj.GetDetector<TS800>();
  // TList    *list    = &(obj.GetObjects());
  // int numobj = list->GetSize();
  
  std::string dirname  = "gretina";
  double timeBank29 = 0;

  std::vector<std::chrono::time_point<std::chrono::high_resolution_clock>> points(6);
  
  if (s800){
    unsigned short bits = s800->GetTrigger().GetRegistr();
    for(int j=0;j<16;j++) {
      if(((bits>>j)&0x0001)) {
        obj.FillHistogram("ungated","trig_bit",20,0,20,j);
      }
    }
  }
  
  if (gretina){
    
    int nHits = gretina->Size();
    // obj.FillHistogram(dirname, "gretina_size",50,0,50,nHits);
    for (int i=0; i < nHits; i++){
      TGretinaHit &hit = gretina->GetGretinaHit(i);
      if (hit.GetPad() > 0) continue;
      points[0] = std::chrono::high_resolution_clock::now();
      
      double phi = hit.GetPhi();
      points[1] = std::chrono::high_resolution_clock::now();
      double theta = hit.GetTheta();
      points[2] = std::chrono::high_resolution_clock::now();
      int cryID = hit.GetCrystalId();
      points[3] = std::chrono::high_resolution_clock::now();
      double timestamp = hit.GetTime();
      points[4] = std::chrono::high_resolution_clock::now();
      double core_energy = hit.GetCoreEnergy();
      points[5] = std::chrono::high_resolution_clock::now();
      
      // if (timeZero == -1 && !std::isnan(timestamp)) timeZero = timestamp;
      // obj.FillHistogram(dirname, "gretina_timestamps_t0",4320,0,43200,(timestamp-timeZero)/TIMESCALE); //12 hours, bins 10 seconds wide
      // std::string timeflag = "all";
      
      // if ( (timestamp-timeZero)/TIMESCALE < 30 ) continue;
      // points[2] = std::chrono::high_resolution_clock::now(); 
      
      // if (timeflag != "") {
      //   obj.FillHistogram(dirname, Form("%s_core_energy",timeflag.c_str()), 4096, 0, 4096, core_energy);
        
      //   //gate on 1408
      //   if (1404 < core_energy && core_energy < 1414){
      //     obj.FillHistogram(dirname, Form("%s_theta_vs_phi_1408",timeflag.c_str()),720,0,360, phi*TMath::RadToDeg(), 360,0,180, theta*TMath::RadToDeg());
      //     obj.FillHistogram(dirname, Form("%s_segment_ids_1408",timeflag.c_str()),1692,0,1692,hit.GetSegmentId() + 36*detMap[cryID],4096,0,4096, hit.GetSegmentEng(0));
      //   }
      //   points[3] = std::chrono::high_resolution_clock::now(); 
      //   //make hit map
      //   obj.FillHistogram(dirname, Form("%s_theta_vs_phi",timeflag.c_str()),720,0,360, phi*TMath::RadToDeg(), 360,0,180, theta*TMath::RadToDeg());
        
      //   //energy vs segment number
      //   obj.FillHistogram(dirname, Form("%s_segment_ids",timeflag.c_str()),1692,0,1692,hit.GetSegmentId() + 36*detMap[cryID],4096,0,4096, hit.GetSegmentEng(0));
        
      //   points[4] = std::chrono::high_resolution_clock::now(); 
      //   //summary spectra
      //   obj.FillHistogram(dirname, Form("%s_summary",timeflag.c_str()),48,0,48,detMap[cryID],4096,0,4096, core_energy);
        
      //   //energy vs interaction point
      //   obj.FillHistogram(dirname, Form("%s_energy_vs_nints",timeflag.c_str()),15,0,15,hit.NumberOfInteractions(),1024,0,2048, core_energy);
      // }

      for (int i=1; i < points.size(); ++i){
        auto durr = std::chrono::duration_cast<std::chrono::nanoseconds>(points[i] - points[0]);
        obj.FillHistogram("timing", Form("point_%d",i),5000,0,50000,(int) durr.count());

        if (i == points.size()-1) {
          obj.FillHistogram("timing",Form("point_%d_vs_cryID",i),48,0,48,detMap[cryID],5000,0,50000,(int) durr.count());
        }
      }
    }
  }

  // if(numobj!=list->GetSize()){
    // list->Sort();
  // }
}
