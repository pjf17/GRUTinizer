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

std::vector<GCutG*> incoming_gates = {};
std::vector<GCutG*> outgoing_gates = {};
std::vector<GCutG*> isoline_gates = {};
GCutG *prompt_timing_gate=0;
GCutG *afp_gate=0;
int gates_loaded=0;

void LoadGates(TRuntimeObjects &obj){
  TList *gates = &(obj.GetGates());
  TIter iter(gates);
  std::cout << "loading gates:" <<std::endl;
  while(TObject *obj = iter.Next()) {
    GCutG *gate = (GCutG*)obj;
    std::string tag = gate->GetTag();
    if(!tag.compare("incoming")) {
      incoming_gates.push_back(gate);
      std::cout << "\t incoming: << " << gate->GetName() << std::endl;
    } else if(!tag.compare("outgoing")) {
      outgoing_gates.push_back(gate);
      std::cout << "\t outgoing: << " << gate->GetName() << std::endl;
    } else if(!tag.compare("isoline")){
      isoline_gates.push_back(gate);
      std::cout << "\t isoline: << " << gate->GetName() << std::endl;
    } else if(!tag.compare("prompt")){
      prompt_timing_gate = new GCutG(*gate);
      std::cout << "\t prompt_timing_gate: << " << gate->GetName() << std::endl;
    } else if(!tag.compare("afp")){
      afp_gate = new GCutG(*gate);
      std::cout << "\t afp_gate: << " << gate->GetName() << std::endl;
    } else {
      std::cout << "\t unknown: << " << gate->GetName() << std::endl;
    }
    gates_loaded++;
  }
  std::cout << "outgoing size: " << outgoing_gates.size() << std::endl;
}

double azimuthalCompton(const TGretinaHit &hit, const TVector3 *beam){
  TVector3 interaction1 = hit.GetIntPosition(0);
  TVector3 interaction2 = hit.GetIntPosition(1);
  TVector3 comptonNorm = interaction1.Cross(interaction2);
  TVector3 reactionNorm = beam->Cross(interaction1);
  TVector3 basisNorm = interaction1.Cross(reactionNorm);
  double angle = reactionNorm.Angle(comptonNorm);
  if (basisNorm.Angle(comptonNorm) > TMath::PiOver2()) angle = TMath::TwoPi() - angle;
  return angle;
}

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

  TList *gates = &(obj.GetGates());

  if(gates_loaded!=gates->GetSize()) {
    LoadGates(obj);
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
      if (prompt_timing_gate) prompt = prompt_timing_gate->IsInside(timeBank29-timestamp, core_energy); 
      if (timeZero == -1 && !std::isnan(timestamp)) timeZero = timestamp;
      
      std::string timeflag = "";
      if (bank29 && prompt) timeflag = "prompt";
      else if (bank29) timeflag = "not-prompt";
      else if ((timestamp-timeZero)/1000000000 < 90) timeflag = "gtTime";

      obj.FillHistogram(dirname, "prompt_gretina_timestamps_t0",500000,0,5000,(timestamp-timeZero)/1000000000);

      if (timeflag != "") {
        obj.FillHistogram(dirname, Form("%s_core_energy",timeflag.c_str()), 8192,0,8192, core_energy);
        obj.FillHistogram(dirname, Form("%s_core_energy_vs_theta",timeflag.c_str()), 8192,0,8192, core_energy, 100, 0, 2.5, theta);
        obj.FillHistogram(dirname, Form("%s_core_energy_vs_crystalID",timeflag.c_str()), 48, 0, 48, detMap[cryID], 8192,0,8192, core_energy);
        obj.FillHistogram(dirname, Form("%s_gretina_theta_vs_phi",timeflag.c_str()),720,0,360,phi,360,0,180,theta*TMath::RadToDeg());
      }
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
        double timestamp = nnhit.GetTime();
        int nInteractions = nnhit.NumberOfInteractions();
        double theta = nnhit.GetTheta()*TMath::RadToDeg();
        
        bool prompt = false; 
        if (prompt_timing_gate) prompt = prompt_timing_gate->IsInside(timeBank29-timestamp, core_energy);
        if (timeZero == -1 && !std::isnan(timestamp)) timeZero = timestamp;

        std::string timeflag = "";
        if ((bank29 && prompt)) timeflag = "prompt";
        else if ((timestamp-timeZero)/1000000000 < 90) timeflag = "gtTime";

        if (timeflag != ""){
          //exclude the ng spectrum (n==3)
          if (n < 3){
            obj.FillHistogram(dirname, Form("%s_core_energy_addback",timeflag.c_str()), 8192,0,8192, core_energy);
            if (nInteractions > 1){
              TVector3 *track = new TVector3(0,0,1);
              double aziCompt = azimuthalCompton(nnhit,track);
              if (50 < theta && theta < 75)
                obj.FillHistogram(dirname, Form("%s_azmthl_compton_theta_cut",timeflag.c_str()),360,0,TMath::TwoPi(),aziCompt,1024,0,2048,core_energy);
              else 
                obj.FillHistogram(dirname, Form("%s_azmthl_compton_anti_theta",timeflag.c_str()),360,0,TMath::TwoPi(),aziCompt,1024,0,2048,core_energy);
              
              //everything
              obj.FillHistogram(dirname, Form("%s_azmthl_compton",timeflag.c_str()),360,0,TMath::TwoPi(),aziCompt,1024,0,2048,core_energy);
            }
          }

          char *multiplicity = Form("%d",n);
          if (n == 3) multiplicity = Form("g");
          obj.FillHistogram(dirname, Form("%s_addback_n%s",timeflag.c_str(),multiplicity), 8192,0,8192, core_energy);
          obj.FillHistogram(dirname, Form("%s_addback_n%s_vs_crystalID",timeflag.c_str(),multiplicity), 48, 0, 48, detMap[cryID], 8192,0,8192, core_energy);
        }
      }
    }
  }
  
  if(numobj!=list->GetSize()){
    list->Sort();
  }
}
