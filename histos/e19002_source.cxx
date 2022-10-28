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

      //bool prompt = prompt_timing_gate->IsInside(timeBank29-timestamp, core_energy); 

      if ((timestamp-timeZero)/1000000000 < 90){
        if (timeZero == -1 && !std::isnan(timestamp)) timeZero = timestamp;
        obj.FillHistogram(dirname, "core_energy", 8192,0,8192, core_energy);
        obj.FillHistogram(dirname, "core_energy_vs_theta", 8192,0,8192, core_energy, 100, 0, 2.5, theta);
        obj.FillHistogram(dirname, "core_energy_vs_crystalID", 48, 0, 48, detMap[cryID], 8192,0,8192, core_energy);
        obj.FillHistogram(dirname, "gretina_theta_vs_phi",720,0,360,phi,360,0,180,theta*TMath::RadToDeg());
        obj.FillHistogram(dirname, "gretina_timestamps_t0",500000,0,5000,(timestamp-timeZero)/1000000000);
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

        //bool prompt = prompt_timing_gate->IsInside(timeBank29-nnhit.GetTime(), core_energy);
        if ((nnhit.GetTime()-timeZero)/1000000000 < 90){
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

            //totals
            int id1 = nnhit.GetCrystalId();
            int id2 = nnhit.GetNeighbor().GetCrystalId();

            if ( PairHit(nnhit,TwoQuadPairs) ){
              obj.FillHistogram(dirname,"gamma_n1_qd2_pair", 8192,0,8192, core_energy);
            }
            if ( PairHit(nnhit,OneQuadPlus) ){
              obj.FillHistogram(dirname,"gamma_n1_qd1_preferred_pair", 8192,0,8192, core_energy);
            }
            if ( PairHit(nnhit,OneQuadDefault) ){
              obj.FillHistogram(dirname,"gamma_n1_qd1_pair", 8192,0,8192, core_energy);
            }

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
        }
      }
    }
  }
  
  if(numobj!=list->GetSize()){
    list->Sort();
  }
}
