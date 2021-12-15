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

std::vector<std::pair<int,int>> redPairs = {
  std::make_pair(46,44),
  std::make_pair(46,48),
  std::make_pair(48,49),
  std::make_pair(50,51),
  std::make_pair(51,47),
  std::make_pair(56,57),
  std::make_pair(57,61),
  std::make_pair(59,58),
  std::make_pair(60,58),
  std::make_pair(60,62),
  std::make_pair(64,62),
  std::make_pair(64,65),
  std::make_pair(69,65),
  std::make_pair(63,67),
  std::make_pair(66,67),
  std::make_pair(66,68),
  std::make_pair(70,68),
  std::make_pair(78,76)
};

std::vector<std::pair<int,int>> goldPairs = {
  std::make_pair(44,45),
  std::make_pair(46,47),
  std::make_pair(48,51),
  std::make_pair(49,50),
  std::make_pair(56,59),
  std::make_pair(57,58),
  std::make_pair(61,60),
  std::make_pair(62,63),
  std::make_pair(64,67),
  std::make_pair(65,66),
  std::make_pair(69,68),
  std::make_pair(70,71),
  std::make_pair(78,79),
};

std::vector<std::pair<int,int>> bluePairs = {
  std::make_pair(44,47),
  std::make_pair(45,46),
  std::make_pair(46,51),
  std::make_pair(48,50),
  std::make_pair(56,58),
  std::make_pair(57,60),
  std::make_pair(60,63),
  std::make_pair(61,62),
  std::make_pair(62,67),
  std::make_pair(64,66),
  std::make_pair(65,68),
  std::make_pair(68,71),
  std::make_pair(76,79),
};

bool PairHit(const TGretinaHit& abhit, std::vector<std::pair<int, int>> &pairs) {
  int cryId1 = abhit.GetCrystalId();
  int cryId2 = abhit.GetNeighbor().GetCrystalId();
  bool hit = false;
  
  for (auto &p : pairs){
    if ( (cryId1 == p.first && cryId2 == p.second) 
        || (cryId2 == p.first && cryId1 == p.second) ) {
        hit = true;
        break;
    }
  }
  return hit;
}

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
    //SINGLES
    int gSize = gretina->Size();
    for (int i=0; i < gSize; i++){
      TGretinaHit &hit = gretina->GetGretinaHit(i);
      double core_energy = hit.GetCoreEnergy();
      double theta = hit.GetTheta();
      int cryID = hit.GetCrystalId();
      //deal with the crystal ID gaps
      // if (cryID >= 40 && cryID < 52 ){
      //   cryID -= 4;
      // } else if (cryID >= 52 && cryID < 72 ) {
      //   cryID -=8;
      // } else if (cryID == 76) {
      //   cryID -=12;
      // } else if (cryID > 77){
      //   cryID -=13;
      // }
      // cryID -= 24;

      //Make all the singles spectra
      obj.FillHistogram(dirname, "gamma_singles", 8192,0,8192, core_energy);
      obj.FillHistogram(dirname, "gamma_vs_theta", 8192,0,8192, core_energy, 100, 0, 2.5, theta);
      obj.FillHistogram(dirname, "gamma_vs_crystalID", 56, 24, 80, cryID, 8192,0,8192, core_energy);
    }

    //NNADDBACK
    //loop over multiplicity
    for (int n=0; n<2; n++){
      //loop over hits for each multiplicity spectrum
      int nnSize = gretina->NNAddbackSize(n);
      for (int i=0; i < nnSize; i++){

        //get hit and hit data 
        TGretinaHit nnhit = gretina->GetNNAddbackHit(n,i);
        int cryID = nnhit.GetCrystalId();
        int ringNum = nnhit.GetRingNumber();
        double nncore_energy = nnhit.GetCoreEnergy();
        
        //exclude the ng spectrum (n==3)
        if (n < 3){
          obj.FillHistogram(dirname, "gamma_addback_prompt", 8192,0,8192, nncore_energy);
          // //GAMMA GAMMA CORRELATION
          // for (int j=i+1; j < nnSize; j++){
          //   if (i==j) continue;
          //   TGretinaHit nnhit2 = gretina->GetNNAddbackHit(n,j);
          //   double nnEnergy_corrected2 = nnhit.GetDopplerYta(s800->AdjustedBeta(GValue::Value("BETA")), s800->GetYta(), &track);
          //   if (prompt_timing_gate->IsInside(timeBank29-nnhit.GetTime(), nnEnergy_corrected)){
          //     obj.FillHistogram(dirname, "gamma_gamma", 8192,0,8192, nnEnergy_corrected2, 8192,0,8192, nnEnergy_corrected);
          //   }
          // }
        }

        if (n == 1){
          //POLARIZATION
          if ( PairHit(nnhit,redPairs) ){
            obj.FillHistogram(dirname,"gamma_addback_red_pair", 8192,0,8192, nncore_energy);
          }

          if ( PairHit(nnhit,goldPairs) ){
            obj.FillHistogram(dirname,"gamma_addback_gold_pair", 8192,0,8192, nncore_energy);
          }

          if ( PairHit(nnhit,bluePairs) ){
            obj.FillHistogram(dirname,"gamma_addback_blue_pair", 8192,0,8192, nncore_energy);
          }

          //N1 Neighbor correlations
          TGretinaHit nnhit1 = nnhit.GetInitialHit();
          TGretinaHit nnhit2 = nnhit.GetNeighbor();
          double singleCrystalEnergy = nnhit2.GetCoreEnergy();
          if ( nnhit1.GetCrystalPosition().Theta() < nnhit2.GetCrystalPosition().Theta() ){
            singleCrystalEnergy = nnhit1.GetCoreEnergy();
          }
          obj.FillHistogram(dirname, Form("total_energy_vs_single_hit_ring%02d_cryID%d",ringNum,cryID),4096,0,8192,singleCrystalEnergy,4096,0,8192,nncore_energy);
        }

        char *multiplicity = Form("%d",n);
        if (n == 3) multiplicity = Form("g");
        obj.FillHistogram(dirname, Form("gamma_n%s_prompt",multiplicity), 8192,0,8192, nncore_energy);
        obj.FillHistogram(dirname, Form("gamma_n%s_ring%02d_crystal%d_prompt",multiplicity,ringNum,cryID),8192,0,8192, nncore_energy);
      }
    }
  }
  
  if(numobj!=list->GetSize()){
    list->Sort();
  }
}
