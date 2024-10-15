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

std::map<int,int> detMap12 = {
  {26, 0}, {30, 1}, {34, 2}, {38, 3}, {25, 4}, {29, 5}, {33, 6}, {37, 7},
  {27, 8}, {31, 9}, {35,10}, {39,11}, {24,12}, {28,13}, {32,14}, {36,15},
  {47,16}, {63,17}, {71,18}, {79,19}, {51,20}, {59,21}, {67,22}, {83,23},
  {50,24}, {58,25}, {66,26}, {82,27}, {44,28}, {60,29}, {68,30}, {76,31},
  {46,32}, {62,33}, {70,34}, {78,35}, {48,36}, {56,37}, {64,38}, {80,39},
  {49,40}, {57,41}, {65,42}, {81,43}, {45,44}, {61,45}, {69,46}, {77,47}
};

std::map<int,int> detMap = {
  {26, 0}, {30, 1}, {34, 2}, {38, 3}, {25, 4}, {29, 5}, {33, 6}, {37, 7},
  {27, 8}, {31, 9}, {35,10}, {39,11}, {24,12}, {28,13}, {32,14}, {36,15},
  {47,16}, {63,17}, {71,18}, {79,19}, {51,20}, {59,21}, {67,22},
  {50,23}, {58,24}, {66,25}, {44,26}, {60,27}, {68,28}, {76,29},
  {46,30}, {62,31}, {70,32}, {78,33}, {48,34}, {56,35}, {64,36},
  {49,37}, {57,38}, {65,39}, {45,40}, {61,41}, {69,42}
};

std::map<int,int> detMapAB = {
  {24, 0}, {28, 1}, {32, 2}, {36, 3},
  {26, 4}, {30, 5}, {34, 6}, {38, 7}, 
  {46, 8}, {62, 9}, {70,10}, 
  {50,11}, {58,12}, {66,13}, {44,14}, 
  {60,15}, {68,16}, {76,17},
  {78,18}, {48,19}, {56,20}, {64,21},
  {25,22}, {29,23}, {33,24}, {37,25},
  {27,26}, {31,27}, {35,28}, {39,29}, 
  {45,30}, {61,31}, {69,32},
  {47,33}, {63,34}, {71,35}, {79,36}, 
  {49,37}, {57,38}, {65,39}, 
  {51,40}, {59,41}, {67,42},
};

std::map<int,int> holeMap = {
  {5,0}, {6,1}, {7,2}, {8,3}, 
  {10,4}, {11,5}, {13,6}, {14,7}, 
  {15,8}, {16,9}, {18,10}, {19,11}
};

std::map<int,int> quadType = {
  {5,1}, {6,1}, {7,1}, {8,1}, 
  {10,3}, {11,2}, {13,2}, {14,3}, 
  {15,2}, {16,3}, {18,3}, {19,2}
};

double lastPointPenalty(double x){
  double val = TMath::TanH((x/100-1.1)/0.5)*exp(-3*(x/100-1.1)+1);
  if (val < 0) val = 0;
  return val + 1;
}

void comptonSort(const TGretinaHit &ghit, int &FP, int &SP) {
  double FOM = 1e10;
  int N = ghit.NumberOfInteractions();
  double ipFactor = 308.76*std::pow(N,-0.542);
  for (int fp=0; fp < N; fp++){
    double E = ghit.GetCoreEnergy();
    double E1 = ghit.GetSegmentEng(fp);
    double er = 511.0/E * E1/(E - E1);
    for (int sp=0; sp < N; sp++){
      if (fp == sp) continue;
      double cosp = TMath::Cos(ghit.GetScatterAngle(fp,sp));
      double E2 = ghit.GetSegmentEng(sp);
      double x = er + cosp;
      double kn = pow((E - E1)/E,2)*((E-E1)/E + E/(E-E1) - pow(TMath::Sin(ghit.GetScatterAngle(fp,sp)),2) );
      double ffom = std::pow(std::abs(1-x),2.0/3)*ghit.GetAlpha(fp,sp)*kn*std::pow(E/E1,3)*ghit.GetLocalPosition(fp).Z(); //*std::pow(E/E2,1.0/3);
      ffom *= std::pow(E/E2,2)*TMath::Sqrt(ghit.GetLocalPosition(sp).Z());
      // ffom /= TMath::Sqrt(abs(E1-125)*abs(E2-125))/E;
      ffom *= lastPointPenalty(E1)*lastPointPenalty(E2);

      if (ffom < FOM) {
        FOM = ffom;
        FP = fp;
        SP = sp;
      }
      // if (fp == sp) continue;
      // double E2 = ghit.GetSegmentEng(sp);
      // double cosp = TMath::Cos(ghit.GetScatterAngle(fp,sp));
      // double x = er + cosp;
      // double kn = pow((E - E1)/E,2)*((E-E1)/E + E/(E-E1) - pow(TMath::Sin(ghit.GetScatterAngle(fp,sp)),2) );
      // // double kn = pow((E - E1)/E,2)*((E-E1)/E + E/(E-E1) - pow(TMath::Sin(ghit.GetScatterAngle(fp,sp)),2) );
      // // double ffom = std::pow(std::abs(1-x),2.0/3)*ghit.GetAlpha(fp,sp)*kn*E/E1;
      // double ffom = std::pow(std::abs(1-x),2.0/3)*ghit.GetAlpha(fp,sp)*std::pow(E/E1,2)*std::pow(E/E2,1.0/3)*kn;
      // // double ffom = std::pow(std::abs(1-x),2.0/3)*ghit.GetAlpha(fp,sp)*std::pow(E/E1,2)*kn;

      // if (ffom < FOM) {
      //   FOM = ffom;
      //   FP = fp;
      //   SP = sp;
      // }
    }
  }
  return;
}

void comptonSortTest(const TGretinaHit &ghit, int &FP, int &SP) {
  double FOM = 1e10;
  FP = 0; SP = 1;
  int N = ghit.NumberOfInteractions();
  double E = ghit.GetCoreEnergy();
  //scale the interaction points so they match the core energy
  double scaleFactor = 0;
  for (int i=0; i < N; i++) scaleFactor += ghit.GetSegmentEng(i);
  scaleFactor = E/scaleFactor;

  for (int fp=0; fp < N; fp++){
    double E1 = ghit.GetSegmentEng(fp)*scaleFactor;
    double er = 511.0/E * E1/(E - E1);

    for (int sp=0; sp < N; sp++){
      if (fp == sp) continue;
      
      double cosp = TMath::Cos(ghit.GetScatterAngle(fp,sp));
      double E2 = ghit.GetSegmentEng(sp)*scaleFactor;
      double x = er + cosp;
      double kn = pow((E - E1)/E,2)*((E-E1)/E + E/(E-E1) - pow(TMath::Sin(ghit.GetScatterAngle(fp,sp)),2) );
      double ffom = std::pow(std::abs(1-x),2.0/3) + std::pow(E1/E - 1/(1+511/E/(1-cosp)),2);
      ffom *= lastPointPenalty(E1)*lastPointPenalty(E2)*std::pow(E/E1,2)*ghit.GetLocalPosition(fp).Z()*ghit.GetAlpha(fp,sp)*kn;
      ffom *= std::pow(E/E2,2) * ghit.GetLocalPosition(sp).Z();

      if (ffom < FOM) {
        FOM = ffom;
        FP = fp;
        SP = sp;
      }
    }
  }
  return;
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
  double timeBank29 = 0;

  if (gretina){
    if (bank29) timeBank29 = bank29->Timestamp();
    //SINGLES
    int gSize = gretina->Size();
    for (int i=0; i < gSize; i++){
      TGretinaHit &hit = gretina->GetGretinaHit(i);
      if (hit.GetPad() != 0) continue;
      double core_energy = hit.GetCoreEnergy();
      double theta = hit.GetTheta();
      double phi = hit.GetPhi();
      int cryID = hit.GetCrystalId();
      double timestamp = hit.GetTime();
      double decompNormChi2 = hit.GetDecompNormChi2();

      bool prompt = false; 
      if (gates.count("prompt") > 0) prompt = gates["prompt"][0]->IsInside(timeBank29-timestamp, core_energy);
      if (timeZero == -1 && !std::isnan(timestamp)) timeZero = timestamp;
      obj.FillHistogram(dirname, "gretina_timestamps_t0",3600,0,7200,(timestamp-timeZero)/TIMESCALE);

      std::string timeflag = "all";
      double timethresh = 15*60; //seconds
      // if (bank29 && prompt) timeflag = "tgated";
      // else if (bank29) timeflag = "not-tgated";
      
      // if ((timestamp-timeZero)/TIMESCALE < timethresh) timeflag = "gtTime"; //timestamp is in 10ns convert to seconds
      if ( (timestamp-timeZero)/TIMESCALE < 30 ) continue;

      if (timeflag != "") {
        obj.FillHistogram(dirname, Form("%s_core_energy",timeflag.c_str()), 8192,0,8192, core_energy);
        obj.FillHistogram(dirname, Form("%s_core_energy_vs_theta",timeflag.c_str()), 180, 0, 180, theta*TMath::RadToDeg(), 4000,0,4000, core_energy);
        obj.FillHistogram(dirname, Form("%s_core_energy_vs_crystalID",timeflag.c_str()), 48, 0, 48, detMap12[cryID], 8192,0,8192, core_energy);
        // obj.FillHistogram(dirname, Form("%s_core_energy_vs_chi2",timeflag.c_str()), 8192,0,8192, core_energy);
        obj.FillHistogram(dirname, Form("%s_gretina_map",timeflag.c_str()),720,0,360,phi*TMath::RadToDeg(),360,0,180,theta*TMath::RadToDeg());
        obj.FillHistogram(dirname, "pad_vs_gretina_timestamps_t0",360,0,3600,(timestamp-timeZero)/TIMESCALE,10,0,10,hit.GetPad());

        if (hit.NumberOfInteractions() > 1){
          double nu = hit.GetScatterAngle();

          int myFP = 0;
          int mySP = 1;
          // comptonSortTest(hit,myFP,mySP);
          double xi = hit.GetXi(nullptr,myFP,mySP);
          obj.FillHistogram(dirname, Form("%s_energy_vs_xi",timeflag.c_str()),360,0,TMath::TwoPi(),xi,2048,0,2048,core_energy);
          double egateLo = 774; //842;       
          double egateHi = 784; //851;       
          if (egateLo < core_energy && core_energy < egateHi) {
            obj.FillHistogram(dirname, Form("%s_gretina_map_E%3.0f",timeflag.c_str(),(egateHi+egateLo)/2),720,0,360,hit.GetPhi(myFP)*TMath::RadToDeg(),360,0,180,hit.GetTheta(myFP)*TMath::RadToDeg());
            // obj.FillHistogram(dirname, Form("%s_gretina_map_CRYS%d_E%3.0f",timeflag.c_str(),detMap12[cryID],(egateHi+egateLo)/2),720,0,360,hit.GetPhi(myFP)*TMath::RadToDeg(),360,0,180,hit.GetTheta(myFP)*TMath::RadToDeg());
            obj.FillHistogram(dirname, Form("%s_chi2norm_E%3.0f",timeflag.c_str(),(egateHi+egateLo)/2),500,0,50,decompNormChi2);
            obj.FillHistogram(dirname, Form("%s_xi_vs_chi2norm_E%3.0f",timeflag.c_str(),(egateHi+egateLo)/2),200,0,50,decompNormChi2,90,0,TMath::TwoPi(),xi);
            obj.FillHistogram(dirname, Form("%s_xi_vs_holenumber_E%3.0f",timeflag.c_str(),(egateHi+egateLo)/2),12,0,12,holeMap[hit.GetHoleNumber()],180,0,TMath::TwoPi(),xi);
            obj.FillHistogram(dirname, Form("%s_chi2norm_vs_holenumber_E%3.0f",timeflag.c_str(),(egateHi+egateLo)/2),12,0,12,holeMap[hit.GetHoleNumber()],200,0,50,decompNormChi2);
            obj.FillHistogram(dirname, Form("%s_phi_vs_holenumber_E%3.0f",timeflag.c_str(),(egateHi+egateLo)/2),12,0,12,holeMap[hit.GetHoleNumber()],360,0,TMath::TwoPi(),hit.GetPhi(myFP));
            // obj.FillHistogram(dirname, Form("%s_chi2norm_HN%d_E%3.0f",timeflag.c_str(),holeMap[hit.GetHoleNumber()],(egateHi+egateLo)/2),500,0,50,decompNormChi2);
            // obj.FillHistogram(dirname, Form("%s_phi_HN%d_E%3.0f",timeflag.c_str(),holeMap[hit.GetHoleNumber()],(egateHi+egateLo)/2),360,0,360,hit.GetPhi(myFP)*TMath::RadToDeg());
            // obj.FillHistogram(dirname, Form("%s_phi_HN%d_ID%d_E%3.0f",timeflag.c_str(),holeMap[hit.GetHoleNumber()],cryID,(egateHi+egateLo)/2),360,0,360,hit.GetPhi(myFP)*TMath::RadToDeg());

            // obj.FillHistogram(dirname, Form("%s_Local_YvsX_tp%d_E%3.0f",timeflag.c_str(),cryID%2,(egateHi+egateLo)/2),50,-50,50,hit.GetLocalPosition(myFP).X(),50,-50,50,hit.GetLocalPosition(myFP).Y());
            // obj.FillHistogram(dirname, Form("%s_Local_YvsX_cryid%d_E%3.0f",timeflag.c_str(),cryID,(egateHi+egateLo)/2),50,-50,50,hit.GetLocalPosition(myFP).X(),50,-50,50,hit.GetLocalPosition(myFP).Y());
            // if (detMap12[cryID] < 16) obj.FillHistogram(dirname, Form("%s_xi_type%d_fwd_E%3.0f",timeflag.c_str(),cryID%2,(egateHi+egateLo)/2),180,0,TMath::TwoPi(),xi);
            // obj.FillHistogram(dirname, Form("%s_Local_YvsX_CRYID%d_E%3.0f",timeflag.c_str(),cryID,(egateHi+egateLo)/2),100,-50,50,hit.GetLocalPosition(myFP).X(),100,-50,50,hit.GetLocalPosition(myFP).Y());
            obj.FillHistogram(dirname, Form("%s_xi_vs_crystalIDmap_E%3.0f",timeflag.c_str(),(egateHi+egateLo)/2),48,0,48,detMap12[cryID],36,0,360,xi*TMath::RadToDeg());
            obj.FillHistogram(dirname, Form("%s_segmentID_vs_crystalID_E%3.0f",timeflag.c_str(),(egateHi+egateLo)/2),48,0,48,detMap[cryID],36,0,36,hit.GetSegmentId());
            // obj.FillHistogram(dirname, Form("%s_xi_vs_crystalID_E%3.0f",timeflag.c_str(),(egateHi+egateLo)/2),120,0,120,cryID,72,0,TMath::TwoPi(),xi);
            // obj.FillHistogram(dirname, Form("%s_xi_qdtype%d_E%3.0f",timeflag.c_str(),quadType[hit.GetHoleNumber()],(egateHi+egateLo)/2),360,0,TMath::TwoPi(),xi);
            // int ihn = holeMap[hit.GetHoleNumber()];
            // if (ihn == 1 || ihn == 2 || (5 < ihn && ihn < 10) )
            //   obj.FillHistogram(dirname, Form("%s_xi_balanced_E847",timeflag.c_str()),360,0,TMath::TwoPi(),xi);
          }
          // if (hit.NumberOfInteractions() < 4) obj.FillHistogram(dirname, Form("%s_energy_vs_xi<4intp",timeflag.c_str()),360,0,TMath::TwoPi(),xi,1024,0,2048,core_energy);
          
          if (theta*TMath::RadToDeg() >= 55 && theta*TMath::RadToDeg() <= 100)
            obj.FillHistogram(dirname, Form("%s_energy_vs_xi_theta_gate",timeflag.c_str()),360,0,TMath::TwoPi(),xi,1024,0,2048,core_energy);

          // if (cryID < 40)
          //   obj.FillHistogram(dirname, Form("%s_fwd_Ecore_vs_xi",timeflag.c_str()),180,0,TMath::Pi(),xi,2048,0,2048,core_energy);
          // else 
          //   obj.FillHistogram(dirname, Form("%s_90deg_Ecore_vs_xi",timeflag.c_str()),180,0,TMath::Pi(),xi,2048,0,2048,core_energy);
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
    //       obj.FillHistogram(dirname, Form("%s_addback_n%s_vs_crystalID",timeflag.c_str(),multiplicity), 48, 0, 48, detMap12[cryID], 8192,0,8192, core_energy);
    //     }
    //   }
    // }
  }
  
  if(numobj!=list->GetSize()){
    list->Sort();
  }
}
