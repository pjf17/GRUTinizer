#include "TRuntimeObjects.h"

#include <iostream>
#include <map>
#include <cstdio>
#include <vector>

#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TRandom.h>
#include <TObject.h>
#include <TLine.h>
#include <TVector3.h>

#include "TGretina.h"
#include "TS800.h"
#include "TBank29.h"
#include "TS800.h"
#include "GCutG.h"

#include "TChannel.h"
#include "GValue.h"

std::map<int,int> detMap = {
  {26, 0}, {30, 1}, {34, 2}, {38, 3}, {25, 4}, {29, 5}, {33, 6}, {37, 7},
  {27, 8}, {31, 9}, {35,10}, {39,11}, {24,12}, {28,13}, {32,14}, {36,15},
  {47,16}, {63,17}, {71,18}, {79,19}, {51,20}, {59,21}, {67,22}, {83,23},
  {50,24}, {58,25}, {66,26}, {82,27}, {44,28}, {60,29}, {68,30}, {76,31},
  {46,32}, {62,33}, {70,34}, {78,35}, {48,36}, {56,37}, {64,38}, {80,39},
  {49,40}, {57,41}, {65,42}, {81,43}, {45,44}, {61,45}, {69,46}, {77,47}
};

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

std::map<int,int> SCATTERPAIRS;

void makemap(){
  int index = 0;
  for (auto &p : redPairs){
    int label1 = p.first*100 + p.second;
    int label2 = p.first + p.second*100;
    SCATTERPAIRS.insert(std::pair<int,int>(label1,index));
    SCATTERPAIRS.insert(std::pair<int,int>(label2,index));
    index++;
  }
  for (auto &p : bluePairs){
    int label1 = p.first*100 + p.second;
    int label2 = p.first + p.second*100;
    SCATTERPAIRS.insert(std::pair<int,int>(label1,index));
    SCATTERPAIRS.insert(std::pair<int,int>(label2,index));
    index++;
  }
  for (auto &p : goldPairs){
    int label1 = p.first*100 + p.second;
    int label2 = p.first + p.second*100;
    SCATTERPAIRS.insert(std::pair<int,int>(label1,index));
    SCATTERPAIRS.insert(std::pair<int,int>(label2,index));
    index++;
  }
}

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

void LoadGates(TList *gates_list, std::map<std::string,std::vector<GCutG*>> &gates){
  TIter iter(gates_list);
  std::cout << "loading gates:" <<std::endl;
  while(TObject *obj = iter.Next()) {
    GCutG *gate = (GCutG*)obj;
    std::string tag = gate->GetTag();
    gates[tag].push_back(gate);
  }
  for (std::map<std::string,std::vector<GCutG*>>::iterator it=gates.begin(); it!=gates.end(); ++it){
    int ngate = it->second.size();
    for (int i=0; i < ngate; i++) std::cout<<it->first<<" << "<<it->second[i]->GetName()<<std::endl;
  }
  return;
}

void CheckGates(std::vector<GCutG*> gates, std::vector<unsigned short> &passed, double x, double y){
  unsigned short ngates = gates.size();
  for (unsigned short i=0; i < ngates; i++){
    if (gates.at(i)->IsInside(x,y)) passed.push_back(i);
  }
  return;
}

bool gates_loaded = false;
std::map<std::string,std::vector<GCutG*>> gates;

const double TIMESCALE = 1E8; // 10ns

// extern "C" is needed to prevent name mangling.
// The function signature must be exactly as shown here,
//   or else bad things will happen.
double timeZero = -1;

extern "C"
void MakeHistograms(TRuntimeObjects& obj) {
//  InitMap();
  TGretina *gretina = obj.GetDetector<TGretina>();
  TBank29  *bank29  = obj.GetDetector<TBank29>();
  TList    *list    = &(obj.GetObjects());
  int numobj = list->GetSize();

  //load in the gates
  if (!gates_loaded) {
    LoadGates(&(obj.GetGates()),gates);
    gates_loaded = true;
  }
  

  //Use this spectrum for the time-energy cut for GRETINA
  if(bank29 && gretina) {
    for(unsigned int i=0;i<gretina->Size();i++) {
      //Time-energy cut
      TGretinaHit &hit = gretina->GetGretinaHit(i);
      obj.FillHistogram("Bank29","Gretina_dop_t0_Bank29_time",
          600,-600,600,bank29->Timestamp()-hit.GetTime(),
          2500,0,10000, hit.GetCoreEnergy());
    }//loop over gretina hits
  }//bank29 and gretina exist

  std::string dirname  = "gretina";
  makemap();
  double timeBank29 = 0;

  if (gretina){
    if (bank29) timeBank29 = bank29->Timestamp();
    //SINGLES
    int gSize = gretina->Size();
    for (int i=0; i < gSize; i++){
      TGretinaHit &hit = gretina->GetGretinaHit(i);
      double core_energy = hit.GetCoreEnergy();
      double theta = hit.GetTheta();
      double phi = hit.GetPhiDeg();
      int cryID = hit.GetCrystalId();
      double timestamp = hit.GetTime();

      bool prompt = false; 
      if (gates["prompt"][0]) prompt = gates["prompt"][0]->IsInside(timeBank29-timestamp, core_energy); 
      if (timeZero == -1 && !std::isnan(timestamp)) timeZero = timestamp;
      
      std::string timeflag = "";
      double timethresh = 15*60; //seconds
      if (bank29 && prompt) timeflag = "prompt";
      else if (bank29) timeflag = "not-prompt";
      
      if ((timestamp-timeZero)/TIMESCALE < timethresh) timeflag = "gtTime"; //timestamp is in 10ns convert to seconds

      obj.FillHistogram(dirname, "prompt_gretina_timestamps_t0",3600,0,3600,(timestamp-timeZero)/TIMESCALE);

      if (timeflag != "") {
        obj.FillHistogram(dirname, Form("%s_core_energy",timeflag.c_str()), 8192,0,8192, core_energy);
        obj.FillHistogram(dirname, Form("%s_core_energy_vs_theta",timeflag.c_str()), 180, 0, 180, theta*TMath::RadToDeg(), 4000,0,4000, core_energy);
        obj.FillHistogram(dirname, Form("%s_core_energy_vs_crystalID",timeflag.c_str()), 48, 0, 48, detMap[cryID], 8192,0,8192, core_energy);
        // obj.FillHistogram(dirname, Form("%s_core_energy_vs_crystal%02d",timeflag.c_str(),detMap[cryID]), 8192,0,8192, core_energy);
        obj.FillHistogram(dirname, Form("%s_gretina_theta_vs_phi",timeflag.c_str()),720,0,360,phi,360,0,180,theta*TMath::RadToDeg());
        if (hit.NumberOfInteractions() > 1){
          double xi = hit.GetXi();
          
          // if (hit.GetScatterAngle() > 1.2*TMath::ACos(1-511/core_energy)) {
          //   obj.FillHistogram(dirname, Form("%s_energy_vs_xi_nugated",timeflag.c_str()),360,0,TMath::TwoPi(),xi,1024,0,2048,core_energy);
          //   double thAngles[6] = {40,55,73,95,120,149}; 
          //   for (int th=0; th < 5; th++){
          //     if (theta*TMath::RadToDeg() >= thAngles[th] && theta*TMath::RadToDeg() < thAngles[th+1]) 
          //       obj.FillHistogram(dirname, Form("%s_energy_vs_xi_nugated_%3.0f-%3.0f",timeflag.c_str(),thAngles[th],thAngles[th+1]),360,0,TMath::TwoPi(),xi,1024,0,2048,core_energy);
          //   }
          // }
          obj.FillHistogram(dirname, Form("%s_energy_vs_xi",timeflag.c_str()),360,0,TMath::TwoPi(),xi,1024,0,2048,core_energy);
          if (theta*TMath::RadToDeg() >= 55 && theta*TMath::RadToDeg() <= 100)
            obj.FillHistogram(dirname, Form("%s_energy_vs_xi_theta_gate",timeflag.c_str()),360,0,TMath::TwoPi(),xi,1024,0,2048,core_energy);

          if (cryID < 40)
            obj.FillHistogram(dirname, Form("%s_IP_fwd_pass_Edop_vs_xi",timeflag.c_str()),360,0,TMath::TwoPi(),xi,2048,0,2048,core_energy);
          else 
            obj.FillHistogram(dirname, Form("%s_IP_90deg_pass_Edop_vs_xi",timeflag.c_str()),360,0,TMath::TwoPi(),xi,2048,0,2048,core_energy);
        }
      }
    }

    //NNADDBACK
    //loop over multiplicity
    // for (int n=0; n<4; n++){
    //   //loop over hits for each multiplicity spectrum
    //   int nnSize = gretina->NNAddbackSize(n);
    //   for (int i=0; i < nnSize; i++){

    //     //get hit and hit data 
    //     TGretinaHit nnhit = gretina->GetNNAddbackHit(n,i);
    //     int cryID = nnhit.GetCrystalId();
    //     int ringNum = nnhit.GetRingNumber();
    //     double core_energy = nnhit.GetCoreEnergy();
    //     double timestamp = nnhit.GetTime();
    //     int nInteractions = nnhit.NumberOfInteractions();
    //     double theta = nnhit.GetTheta()*TMath::RadToDeg();
        
    //     bool prompt = false; 
    //     if (prompt_timing_gate) prompt = prompt_timing_gate->IsInside(timeBank29-timestamp, core_energy);
    //     if (timeZero == -1 && !std::isnan(timestamp)) timeZero = timestamp;

    //     std::string timeflag = "";
    //     if ((bank29 && prompt)) timeflag = "prompt";
    //     else if ((timestamp-timeZero)/1000000000 < 90) timeflag = "gtTime";

    //     if (timeflag != ""){
    //       //exclude the ng spectrum (n==3)
    //       if (n < 3){
    //         obj.FillHistogram(dirname, Form("%s_core_energy_addback",timeflag.c_str()), 8192,0,8192, core_energy);
    //       }

    //       char *multiplicity = Form("%d",n);
    //       if (n == 3) multiplicity = Form("g");
    //       obj.FillHistogram(dirname, Form("%s_addback_n%s",timeflag.c_str(),multiplicity), 8192,0,8192, core_energy);
    //       obj.FillHistogram(dirname, Form("%s_addback_n%s_vs_crystalID",timeflag.c_str(),multiplicity), 48, 0, 48, detMap[cryID], 8192,0,8192, core_energy);
    //     }
    //   }
    // }
  }
  
  if(numobj!=list->GetSize()){
    list->Sort();
  }
}
