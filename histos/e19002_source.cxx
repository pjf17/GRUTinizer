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

std::vector<std::map<std::string,std::pair<int,int>>> scatterGroups = {
  {{"A-A",std::make_pair(65,69)}, {"A-B",std::make_pair(65,68)}, {"B-B",std::make_pair(66,68)}},
  {{"A-A",std::make_pair(63,67)}, {"A-B",std::make_pair(62,67)}, {"B-B",std::make_pair(62,64)}},
  {{"A-A",std::make_pair(57,61)}, {"A-B",std::make_pair(57,60)}, {"B-B",std::make_pair(58,60)}},
  {{"A-A",std::make_pair(47,51)}, {"A-B",std::make_pair(46,51)}, {"B-B",std::make_pair(46,48)}}
};

bool checkScatterType(const TGretinaHit& abhit, std::pair<int,int> crystals){
  int cryId1 = abhit.GetCrystalId();
  int cryId2 = abhit.GetNeighbor().GetCrystalId();
  return (cryId1 == crystals.first && cryId2 == crystals.second) || (cryId1 == crystals.second && cryId2 == crystals.first);
}

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
  std::make_pair(69,70)
};

std::vector<std::pair<int,int>> TwoQuadPairs = {
  std::make_pair(46,48),
  std::make_pair(47,51),
  std::make_pair(46,51),
  std::make_pair(61,57),
  std::make_pair(60,58),
  std::make_pair(57,60),
  std::make_pair(62,64),
  std::make_pair(63,67),
  std::make_pair(62,67),
  std::make_pair(65,69),
  std::make_pair(66,68),
  std::make_pair(65,68)
};

std::vector<std::pair<int,int>> OneQuadPairs = {
  std::make_pair(44,45),
  std::make_pair(46,47),
  std::make_pair(45,46),
  std::make_pair(44,47),
  std::make_pair(44,46),
  std::make_pair(48,51),
  std::make_pair(49,50),
  std::make_pair(49,48),
  std::make_pair(50,51),
  std::make_pair(48,50),
  std::make_pair(56,59),
  std::make_pair(57,58),
  std::make_pair(56,57),
  std::make_pair(58,59),
  std::make_pair(56,58),
  std::make_pair(61,60),
  std::make_pair(62,63),
  std::make_pair(62,61),
  std::make_pair(60,63),
  std::make_pair(60,62),
  std::make_pair(64,67),
  std::make_pair(65,66),
  std::make_pair(65,64),
  std::make_pair(66,67),
  std::make_pair(64,66),
  std::make_pair(69,68),
  std::make_pair(70,71),
  std::make_pair(69,70),
  std::make_pair(68,71),
  std::make_pair(70,68),
  std::make_pair(78,76),
  std::make_pair(78,79),
  std::make_pair(76,79)
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
      double phi = hit.GetPhiDeg();
      int cryID = hit.GetCrystalId();

      obj.FillHistogram(dirname, "core_energy", 8192,0,8192, core_energy);
      obj.FillHistogram(dirname, "core_energy_vs_theta", 8192,0,8192, core_energy, 100, 0, 2.5, theta);
      obj.FillHistogram(dirname, "core_energy_vs_crystalID", 56, 24, 80, cryID, 8192,0,8192, core_energy);
      obj.FillHistogram(dirname, "gretina_theta_vs_phi",720,0,360,phi,360,0,180,theta*TMath::RadToDeg());
    }

    //NNADDBACK
    //loop over multiplicity
    for (int n=0; n<4; n++){
      //loop over hits for each multiplicity spectrum
      int nnSize = gretina->NNAddbackSize(n);
      for (int i=0; i < nnSize; i++){

        //get hit and hit data 
        TGretinaHit nnhit = gretina->GetNNAddbackHit(n,i);
        int cryID = nnhit.GetCrystalId();
        int ringNum = nnhit.GetRingNumber();
        double core_energy = nnhit.GetCoreEnergy();
        
        //make sure hits are prompt
        //exclude the ng spectrum (n==3)
        if (n < 3){
          obj.FillHistogram(dirname, "core_energy_addback", 8192,0,8192, core_energy);
        }

        if (n == 1) {
          //GAMMA GAMMA CORRELATION
          // for (int j=0; j < nnSize; j++){
          //   if (i==j) continue;
          //   TGretinaHit nnhit2 = gretina->GetNNAddbackHit(n,j);
          //   double nnEnergy_corrected2 = nnhit2.GetDopplerYta(s800->AdjustedBeta(GValue::Value("BETA")), s800->GetYta(), &track);
          //   obj.FillHistogram(dirname, "gamma_gamma", 8192,0,8192, nnEnergy_corrected2, 8192,0,8192, nnEnergy_corrected);
          // }

          //POLARIZATION
          if ( PairHit(nnhit,redPairs) ){
            obj.FillHistogram(dirname,"gamma_corrected_addback_red_pair", 8192,0,8192, core_energy);
          }

          if ( PairHit(nnhit,goldPairs) ){
            obj.FillHistogram(dirname,"gamma_corrected_addback_gold_pair", 8192,0,8192, core_energy);
          }

          if ( PairHit(nnhit,bluePairs) ){
            obj.FillHistogram(dirname,"gamma_corrected_addback_blue_pair", 8192,0,8192, core_energy);
          }

          //SCATTER TYPE
          int ngroups = (int) scatterGroups.size();
          for (int i=0; i < ngroups; i++){
            std::map<std::string, std::pair<int,int>>::iterator it = scatterGroups[i].begin();
            std::map<std::string, std::pair<int,int>>::iterator end = scatterGroups[i].end();
            while (it != end){
                if (checkScatterType(nnhit,it->second)){
                  obj.FillHistogram(dirname,Form("gamma_corrected_n1_grp%d_%s",i+1,it->first.c_str()),8192,0,8192, core_energy);
                }
                it++;
            }
          }

          //totals
          if (cryID > 40){ //only use 90 degree quads
            int id1 = nnhit.GetCrystalId();
            int id2 = nnhit.GetNeighbor().GetCrystalId();
            if (id1%2 == 1 && id2%2 == 1){
              obj.FillHistogram(dirname,Form("gamma_n1_A-A"),1600,0,1600, core_energy);
            }
            else if (id1%2 == 0 && id2%2 == 0){
              obj.FillHistogram(dirname,Form("gamma_n1_B-B"),1600,0,1600, core_energy);
            }
            else{
              obj.FillHistogram(dirname,Form("gamma_n1_A-B"),1600,0,1600, core_energy);
            }
          }

          if ( PairHit(nnhit,TwoQuadPairs) ){
            obj.FillHistogram(dirname,"gamma_n1_qd2_pair", 8192,0,8192, core_energy);
          }
          if ( PairHit(nnhit,OneQuadPairs) ){
            obj.FillHistogram(dirname,"gamma_n1_qd1_pair", 8192,0,8192, core_energy);
          }
          // int id1 = nnhit.GetCrystalId();
          // int id2 = nnhit.GetNeighbor().GetCrystalId();
          // if ((id1 == 65 && id2 == 69) || (id2 == 65 && id1 == 69)) 
          //   obj.FillHistogram(dirname,"gamma_corrected_n1_A-A_scatter",8192,0,8192, core_energy); //both are type A
          // if ((id1 == 66 && id2 == 68) || (id2 == 66 && id1 == 68)) 
          //   obj.FillHistogram(dirname,"gamma_corrected_n1_B-B_scatter",8192,0,8192, core_energy); //both are type B
          // if ((id1 == 65 && id2 == 68) || (id1 == 65 && id2 == 68))
          //   obj.FillHistogram(dirname,"gamma_corrected_n1_A-B_scatter",8192,0,8192, core_energy); //type A and B
        
          // double singleCrystalEnergy = nnhit2.GetCoreEnergy();
          // // if ( nnhit1.GetCrystalPosition().Theta() > nnhit2.GetCrystalPosition().Theta() ){
          // //   singleCrystalEnergy = nnhit1.GetCoreEnergy();
          // // }
          // obj.FillHistogram(dirname, Form("total_energy_vs_single_hit_ring%02d_cryID%d",ringNum,cryID),4096,0,8192,singleCrystalEnergy,4096,0,8192,nnEnergy_corrected);
          
          // obj.FillHistogram(dirname, Form("total_energy_swapped_vs_single_hit_ring%02d_cryID%d",ringNum,cryID),4096,0,8192,nnhit1.GetCoreEnergy(),4096,0,8192,swappedEnergy);                  
        }

        char *multiplicity = Form("%d",n);
        if (n == 3) multiplicity = Form("g");
        obj.FillHistogram(dirname, Form("gamma_corrected_n%s",multiplicity), 8192,0,8192, core_energy);
        if (n == 0){
          obj.FillHistogram(dirname, Form("gamma_corrected_n%s_vs_cryID",multiplicity),56, 24, 80, cryID, 8192,0,8192, core_energy);
          obj.FillHistogram(dirname, Form("gamma_corrected_n%s_cr%d",multiplicity,cryID), 8192,0,8192, core_energy);
        }
      }
    }
  }
  
  if(numobj!=list->GetSize()){
    list->Sort();
  }
}
