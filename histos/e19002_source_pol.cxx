#include "TRuntimeObjects.h"

#include <iostream>
#include <map>
#include <cstdio>
#include <vector>
#include <algorithm>
#include <queue>

#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TRandom3.h>
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

std::map<int,int> detMapRing = {
  {26, 0}, {30, 1}, {34, 2}, {38, 3}, {25, 4}, {29, 5}, {33, 6}, {37, 7},
  {27, 8}, {31, 9}, {35,10}, {39,11}, {24,12}, {28,13}, {32,14}, {36,15},
  {47,16}, {63,17}, {71,18}, {79,19}, {51,20}, {59,21}, {67,22}, {83,23},
  {50,24}, {58,25}, {66,26}, {82,27}, {44,28}, {60,29}, {68,30}, {76,31},
  {46,32}, {62,33}, {70,34}, {78,35}, {48,36}, {56,37}, {64,38}, {80,39},
  {49,40}, {57,41}, {65,42}, {81,43}, {45,44}, {61,45}, {69,46}, {77,47}
};

int thetaGate(double angle){
  if (angle < 0.933) return 0;
  else if (0.933 <= angle && angle < 1.094) return 1;
  else if (1.094 <= angle && angle < 1.275) return 2;
  else if (1.275 <= angle && angle <	1.485)  return 3;
  else if (1.485 <= angle && angle <	1.73)  return 4;
  else if (1.73	<= angle)  return 5;
}

double lastPointPenalty(double x){
  double val = TMath::TanH((x/100-1.1)/0.5)*exp(-3*(x/100-1.1)+1);
  if (val < 0) val = 0;
  return val + 1;
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

const char * polCode(double angle){
  angle = angle*TMath::RadToDeg();
  if (angle < 30) return "para";
  else if (angle > 60) return "perp";
  else return "mid";
}

bool gates_loaded = false;
std::map<std::string,std::vector<GCutG*>> gates;

const double TIMESCALE = 1E8; // 10ns

// extern "C" is needed to prevent name mangling.
// The function signature must be exactly as shown here,
//   or else bad things will happen.
double timeZero = -1;
const int QUEUESIZE = 20;
std::map<std::string,std::vector<std::array<TVector3,2>>> mixing_queue;

extern "C"
void MakeHistograms(TRuntimeObjects& obj) {
//  InitMap();
  TGretina *gretina = obj.GetDetector<TGretina>();
  TBank29  *bank29  = obj.GetDetector<TBank29>();
  TList    *list    = &(obj.GetObjects());
  int numobj = list->GetSize();

  // static TRandom3 *rand_gen = new TRandom3(59953614);

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
  double timeBank29 = 0;

  if (gretina){
    if (bank29) timeBank29 = bank29->Timestamp();
    //SINGLES
    int gSize = gretina->Size();
    std::map<std::string,int> correlations;
    correlations.insert(std::make_pair("1238",-1)); //56Co
    correlations.insert(std::make_pair("2598",-1)); //56Co
    correlations.insert(std::make_pair("344",-1)); //152Eu

    std::vector<std::map<std::string,TGretinaHit>> ghits(gSize);
    // std::vector<int> roll_indices;
    // roll_indices.reserve(gSize);
    for (int i=0; i < gSize; i++) {
      TGretinaHit &hit = gretina->GetGretinaHit(i);
      if (hit.GetPad() > 0) continue;
      
      if ( 1235 < hit.GetCoreEnergy() && hit.GetCoreEnergy() < 1242 ) correlations["1238"] = i;
      if ( 2592 < hit.GetCoreEnergy() && hit.GetCoreEnergy() < 2602 ) correlations["2598"] = i;
      if ( 339 < hit.GetCoreEnergy() && hit.GetCoreEnergy() < 348 ) correlations["344"] = i;

      std::map<std::string,TGretinaHit> sortTypes;
      sortTypes["MainInt"] = TGretinaHit(hit);
      sortTypes["ComptonSort"] = TGretinaHit(hit);
      sortTypes["Tracking"] = TGretinaHit(hit);
      
      //run sorting algorithms
      sortTypes["ComptonSort"].PICCSort();
      sortTypes["Tracking"].TrackingSort();

      //add to list
      ghits[i] = sortTypes;
      // roll_indices.push_back(i);
    }
    // int nAvailIdx = roll_indices.size();

    for (int i=0; i < gSize; i++){
      if (ghits[i].empty()) continue;
      
      // //find a good random index for the correlation normalization
      // int rand_idx = -1;
      // if (nAvailIdx > 1){
      //   do {
      //     rand_idx = roll_indices[rand_gen->Integer(nAvailIdx)];
      //   } while (rand_idx == i);
      // }
      
      //grab generic info
      double core_energy = ghits[i]["MainInt"].GetCoreEnergy();
      double theta = ghits[i]["MainInt"].GetTheta();
      double phi = ghits[i]["MainInt"].GetPhi();
      int cryID = ghits[i]["MainInt"].GetCrystalId();
      int Ninteractions = ghits[i]["MainInt"].NumberOfInteractions();
      double timestamp = ghits[i]["MainInt"].GetTime();
      
      // if (cryID < 40) continue;

      bool prompt = false; 
      if (gates.count("prompt") > 0) prompt = gates["prompt"][0]->IsInside(timeBank29-timestamp, core_energy);
      if (timeZero == -1 && !std::isnan(timestamp)) timeZero = timestamp;
      obj.FillHistogram(dirname, "gretina_timestamps_t0",1440,0,43200,(timestamp-timeZero)/TIMESCALE); //12 hours, bins in 30 seconds
      std::string timeflag = "all";
      
      if ( (timestamp-timeZero)/TIMESCALE < 30 ) continue;
      
      if (timeflag != "") {
        //make hit map
        obj.FillHistogram(dirname, Form("%s_theta_vs_phi",timeflag.c_str()),360,0,360, phi*TMath::RadToDeg(), 180,0,180, theta*TMath::RadToDeg());
        //make total spectrum
        obj.FillHistogram(dirname, Form("%s_core_energy",timeflag.c_str()), 8192,0,8192, core_energy);
        obj.FillHistogram(dirname, Form("%s_core_energy_summary",timeflag.c_str()),48,0,48,detMapRing[cryID],4096,0,4096, core_energy);
    
        //loop over sort types
        for (auto stHit = ghits[i].begin(); stHit != ghits[i].end(); ++stHit) {
          bool overOneInt = Ninteractions > 1;
          if (overOneInt){
            //fill total, uncorrelated hists
            obj.FillHistogram(dirname, 
              Form("%s_core_energy_%s_%s",timeflag.c_str(),stHit->first.c_str(),polCode(stHit->second.GetXiSimple(nullptr))),
              4096,0,4096,core_energy);
            obj.FillHistogram(dirname, 
              Form("%s_energy_vs_xi_%s",timeflag.c_str(),stHit->first.c_str()),
              180,0,180,stHit->second.GetXi(nullptr)*TMath::RadToDeg(),1024,0,2048,core_energy);
          }

          //loop over correlation
          for (auto it = correlations.begin(); it != correlations.end(); ++it) {
            if (it->second != -1 && it->second != i) { //only do if correlated hit exists and is not self
              TGretinaHit hitTemp = gretina->GetGretinaHit(it->second);
              if (stHit->first == "ComptonSort") hitTemp.PICCSort();
              else if (stHit->first == "Tracking") hitTemp.TrackingSort();
              TVector3 g1dir = hitTemp.GetPosition();

              //for 152Eu GRETINA response, get 744 gamma in correlation with 344
              // if (it->first == "1238" || it->first == "2598"){
              if (it->first == "344"){
                // if ( 842 < core_energy && core_energy < 850 ) {
                if ( 776 < core_energy && core_energy < 782 ) {
                  mixing_queue[stHit->first].push_back({stHit->second.GetPosition(),g1dir});
                }
                //back of queue is current hit, calculate the 4x fold mixing with the whole queue
                int current_size = mixing_queue[stHit->first].size();
                if (current_size > QUEUESIZE) { 
                  for (int mq=0; mq < current_size-1; mq++){
                    for (int cidx=0; cidx < 2; ++cidx){
                      for (int oidx=0; oidx < 2; ++oidx){
                        double mixing_theta = mixing_queue[stHit->first].back()[cidx].Angle(mixing_queue[stHit->first][mq][oidx])*TMath::RadToDeg();
                        obj.FillHistogram(dirname, Form("%s_polar_%s_response_%s",timeflag.c_str(),stHit->first.c_str(),it->first.c_str()), 180,0,180, mixing_theta);
                      }
                    }
                  }
                  // remove first element
                  mixing_queue[stHit->first].erase(mixing_queue[stHit->first].begin());
                }
              }

              //calculate polar angles
              obj.FillHistogram(dirname, Form("%s_crrl_%s_polar_%s",timeflag.c_str(),it->first.c_str(),stHit->first.c_str()), 180,0,180, g1dir.Angle(stHit->second.GetPosition())*TMath::RadToDeg(), 2048,0,2048,core_energy);

              //calculate xi angles
              if (overOneInt){
                double xi = stHit->second.GetXi(&g1dir)*TMath::RadToDeg();
                double xiSimple = stHit->second.GetXiSimple(&g1dir);
                obj.FillHistogram(dirname, Form("%s_crrl_%s_core_energy_%s_%s",timeflag.c_str(),it->first.c_str(),stHit->first.c_str(),polCode(xiSimple)), 4096, 0, 4096, core_energy);
                obj.FillHistogram(dirname, Form("%s_crrl_%s_core_energy_vs_xi_%s",timeflag.c_str(),it->first.c_str(),stHit->first.c_str()),180,0,180,xi,1024,0,2048,core_energy);
              }
            } 
            
          // } else if (it->second != -1) {
          //   obj.FillHistogram(dirname, Form("%s_crrl_%s_core_energy_remain",timeflag.c_str(),it->first.c_str()), 4096, 0, 4096, core_energy);
          // }
          }
        }
      }
    }
  }

  if(numobj!=list->GetSize()){
    list->Sort();
  }
}
