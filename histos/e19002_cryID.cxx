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

//24
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
  std::make_pair(78,76),
  std::make_pair(78,80),
  std::make_pair(80,81),
  std::make_pair(45,81),
  std::make_pair(79,83),
  std::make_pair(82,83),
  std::make_pair(44,82)
};

//16
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
  std::make_pair(76,77),
  std::make_pair(80,83),
  std::make_pair(81,82)
};

//18
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
  std::make_pair(69,70),
  std::make_pair(77,78),
  std::make_pair(78,83),
  std::make_pair(80,82),
  std::make_pair(44,81)
};

//16
std::vector<std::pair<int,int>> OneQuadPlus = {
  std::make_pair(44,47),
  std::make_pair(45,46),
  std::make_pair(60,63),
  std::make_pair(61,62),
  std::make_pair(68,71),
  std::make_pair(76,79),
  std::make_pair(69,70),
  std::make_pair(49,50),
  std::make_pair(56,59),
  std::make_pair(57,58),
  std::make_pair(64,67),
  std::make_pair(65,66),
  std::make_pair(48,51),
  std::make_pair(77,78),
  std::make_pair(80,83),
  std::make_pair(81,82),
};

//24
std::vector<std::pair<int,int>> OneQuadDefault = {
  std::make_pair(44,45),
  std::make_pair(46,47),
  std::make_pair(44,46),
  std::make_pair(48,49),
  std::make_pair(50,51),
  std::make_pair(48,50),
  std::make_pair(56,57),
  std::make_pair(58,59),
  std::make_pair(56,58),
  std::make_pair(60,61),
  std::make_pair(62,63),
  std::make_pair(60,62),
  std::make_pair(64,65),
  std::make_pair(66,67),
  std::make_pair(64,66),
  std::make_pair(68,69),
  std::make_pair(70,71),
  std::make_pair(68,70),
  std::make_pair(76,78),
  std::make_pair(78,79),
  std::make_pair(77,76),
  std::make_pair(80,82),
  std::make_pair(81,80),
  std::make_pair(82,83),
};


//18
std::vector<std::pair<int,int>> TwoQuadPairs = {
  std::make_pair(46,48),
  std::make_pair(46,51),
  std::make_pair(47,51),
  std::make_pair(61,57),
  std::make_pair(58,60),
  std::make_pair(57,60),
  std::make_pair(62,64),
  std::make_pair(63,67),
  std::make_pair(62,67),
  std::make_pair(65,69),
  std::make_pair(66,68),
  std::make_pair(65,68),
  std::make_pair(78,80),
  std::make_pair(78,83),
  std::make_pair(79,83),
  std::make_pair(45,81),
  std::make_pair(44,81),
  std::make_pair(44,82),
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
  TList    *list    = &(obj.GetObjects());
  int numobj = list->GetSize();
  
  if (gretina){
    int nHits = gretina->NNAddbackSize();
    for (int i=0; i < nHits; i++){
      TGretinaHit hit = gretina->GetNNAddbackHit(i);
      double phi = hit.GetPhiDeg();
      double theta = 180 - hit.GetThetaDeg();
      int cryID = hit.GetCrystalId();
      
      std::string polColor = "blank";
      if (PairHit(hit,redPairs)) polColor = "red";
      else if (PairHit(hit,goldPairs)) polColor = "gold";
      else if (PairHit(hit,bluePairs)) polColor = "blue";

      std::string quadType = "blank";
      if (PairHit(hit,OneQuadPlus)) quadType = "qd1+";
      else if (PairHit(hit,OneQuadDefault)) quadType = "qd1";
      else if (PairHit(hit,TwoQuadPairs)) quadType = "qd2";

      obj.FillHistogram("crystals",Form("crystal%d",cryID),180,0,180,theta,360,0,360,phi);
      // obj.FillHistogram("crystals",Form("ring%02d",hit.GetRingNumber()),180,0,180,theta,360,0,360,phi);
      obj.FillHistogram("crystals",Form("%s",polColor.c_str()),180,0,180,theta,360,0,360,phi);
      obj.FillHistogram("crystals",Form("%s",quadType.c_str()),180,0,180,theta,360,0,360,phi);
      obj.FillHistogram("crystals",Form("%s_%s",polColor.c_str(),quadType.c_str()),180,0,180,theta,360,0,360,phi);
      obj.FillHistogram("crystals","total",180,0,180,theta,360,0,360,phi);
    }
  }

  if(numobj!=list->GetSize()){
    list->Sort();
  }
}
