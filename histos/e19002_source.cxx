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

std::vector<std::pair<double,double>> gamma_gates = {
  std::make_pair(1400,1412),
  std::make_pair(1103,1117),
  std::make_pair(954,970),
  std::make_pair(771,782),
  std::make_pair(336,348),
  std::make_pair(118,125)
};

// extern "C" is needed to prevent name mangling.
// The function signature must be exactly as shown here,
//   or else bad things will happen.
extern "C"
void MakeHistograms(TRuntimeObjects& obj) {
//  InitMap();
  TGretina *gretina = obj.GetDetector<TGretina>();
  // TBank29  *bank29  = obj.GetDetector<TBank29>();
  TList    *list    = &(obj.GetObjects());
  int numobj = list->GetSize();

  std::string dirname  = "gretina";
  
  if (gretina){
    //ADDBACK STUFF
    int nABHits = gretina->AddbackSize();
    for (int i=0; i < nABHits; i++){
      TGretinaHit abhit = gretina->GetAddbackHit(i);
      double abEnergy_corrected = abhit.GetDoppler(0.0);
      obj.FillHistogram(dirname, "gamma_corrected_addback_prompt", 8192,0,8192, abEnergy_corrected);
    }

    //OTHER STUFF
    for (unsigned int i = 0; i < gretina->Size(); i++){
      TGretinaHit &hit = gretina->GetGretinaHit(i);
      double energy = hit.GetDoppler(0.0);
      obj.FillHistogram(dirname, "gamma_singles", 8192,0,8192, energy);
      
      TVector3 position = hit.GetPosition();
      double theta = position.Theta();
        
      obj.FillHistogram(dirname, "core_energy_vs_theta", 8192,0,8192, hit.GetCoreEnergy(), 100, 0, 4, theta);
      obj.FillHistogram(dirname, "gamma_singles_prompt", 8192,0,8192, energy);
    }

    //CRYID Hits
    dirname = "crystal-hits";
    int nHits = gretina->Size();
    for (int i=0; i < nHits; i++){
      TGretinaHit hit = gretina->GetGretinaHit(i);
      int cryID = hit.GetCrystalId();
      double energy = hit.GetDoppler(0.0);
      obj.FillHistogram(dirname,"crystal-hits-total",100,0,100,cryID);
      
      int nGammaGates = (int) gamma_gates.size();
      for (int j=0; j < nGammaGates; j++){
        if (energy > gamma_gates[j].first && energy < gamma_gates[j].second){
          obj.FillHistogram(dirname,Form("cryHit(%4.0f,%4.0f)",gamma_gates[j].first,gamma_gates[j].second),100,0,100,cryID);
        }
      }
    }
  }
  
  if(numobj!=list->GetSize()){
    list->Sort();
  }
}
