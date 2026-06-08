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

int thetaGate(double angle){
  //attempt 1
  // if (0.701 <= angle && angle < 0.933) return 0;
  // else if (0.933 <= angle && angle < 1.094) return 1;
  // else if (1.094 <= angle && angle < 1.275) return 2;
  // else if (1.275 <= angle && angle <	1.485)  return 3;
  // else if (1.485 <= angle && angle <	1.73)  return 4;
  // else if (1.73	<= angle && angle < 2.02)  return 5;
  //attempt 2
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

std::string polCode(double angle){
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
    std::map<std::string,int> correlations;
    correlations.insert(std::make_pair("1238",-1));
    correlations.insert(std::make_pair("2598",-1));
    correlations.insert(std::make_pair("344",-1));
    for (int i=0; i < gSize; i++) {
      TGretinaHit &hit = gretina->GetGretinaHit(i);
      if (hit.GetPad() > 0) continue;
      if ( 1235 < hit.GetCoreEnergy() && hit.GetCoreEnergy() < 1242 ) correlations["1238"] = i;
      if ( 2592 < hit.GetCoreEnergy() && hit.GetCoreEnergy() < 2602 ) correlations["2598"] = i;
      if ( 339 < hit.GetCoreEnergy() && hit.GetCoreEnergy() < 348 ) correlations["344"] = i;
    }

    for (int i=0; i < gSize; i++){
      TGretinaHit &hit = gretina->GetGretinaHit(i);
      if (hit.GetPad() > 0) continue;
      TGretinaHit hitMain;
      TGretinaHit hitTrack;
      hit.Copy(hitMain);
      hit.Copy(hitTrack);
      hitTrack.TrackingSort();
      hit.PICCSort();
      double core_energy = hit.GetCoreEnergy();
      double theta = hit.GetTheta();
      double phi = hit.GetPhi();
      int cryID = hit.GetCrystalId();
      double timestamp = hit.GetTime();
      double decompNormChi2 = hit.GetDecompNormChi2();
      double dop_energy = hit.GetDoppler(0.3346);

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
        // bool bin75 = hit.GetTheta()*TMath::RadToDeg() > 65 && hit.GetTheta()*TMath::RadToDeg() < 85;
        // bool bin100 = hit.GetTheta()*TMath::RadToDeg() > 90 && hit.GetTheta()*TMath::RadToDeg() < 115;
        // if (bin75) {
          // obj.FillHistogram(dirname, Form("%s_core_energy_bin75",timeflag.c_str()),10000,0,10000, core_energy);
          // obj.FillHistogram(dirname, Form("%s_theta_bin75",timeflag.c_str()),180,0,180, hit.GetTheta()*TMath::RadToDeg());
        // }
        // if (bin100) {
          // obj.FillHistogram(dirname, Form("%s_core_energy_bin100",timeflag.c_str()),10000,0,10000, core_energy);
          // obj.FillHistogram(dirname, Form("%s_theta_bin100",timeflag.c_str()),180,0,180, hit.GetTheta()*TMath::RadToDeg());
        // }
    
        // if (hit.NumberOfInteractions() > 1){
        //   for (auto it = correlations.begin(); it != correlations.end(); ++it) {
        //     if (it->second != -1 && it->second != i) {
        //       TVector3 g1dir = (gretina->GetGretinaHit(it->second)).GetPosition();
        //       double xi = hit.GetXi(&g1dir)*TMath::RadToDeg();
        //       double xiMain = hitMain.GetXi(&g1dir)*TMath::RadToDeg();
        //       obj.FillHistogram(dirname, Form("%s_crrl_%s_core_energy_vs_xi",timeflag.c_str(),it->first.c_str()),180,0,180,xi,1024,0,1024,core_energy);
        //       obj.FillHistogram(dirname, Form("%s_crrl_%s_core_energy_vs_ximain",timeflag.c_str(),it->first.c_str()),180,0,180,xiMain,1024,0,1024,core_energy);

        //       obj.FillHistogram(dirname, Form("%s_crrl_%s_core_energy",timeflag.c_str(),it->first.c_str()), 4096, 0, 4096, core_energy);
        //       obj.FillHistogram(dirname, Form("%s_crrl_%s_core_energy_%s",timeflag.c_str(),it->first.c_str(),polCode(hit.GetXiSimple(&g1dir)).c_str()),4096,0,4096,core_energy);
        //       obj.FillHistogram(dirname, Form("%s_crrl_%s_core_energy_main_%s",timeflag.c_str(),it->first.c_str(),polCode(hitMain.GetXiSimple(&g1dir)).c_str()),4096,0,4096,core_energy);

        //     } else if (it->second != -1) {
        //       obj.FillHistogram(dirname, Form("%s_crrl_%s_core_energy_remain",timeflag.c_str(),it->first.c_str()), 4096, 0, 4096, core_energy);
        //     }
        //   }
        // }

        obj.FillHistogram(dirname, Form("%s_core_energy_%d",timeflag.c_str(),thetaGate(hit.GetTheta())), 8192,0,8192, core_energy);
        // obj.FillHistogram(dirname, Form("%s_core_energy_%d_vs_theta",timeflag.c_str(),thetaGate(hit.GetTheta())), 180, 0, 180, theta*TMath::RadToDeg(), 8192,0,8192, core_energy);
        // obj.FillHistogram(dirname, Form("%s_core_energy_vs_theta",timeflag.c_str()), 180, 0, 180, theta*TMath::RadToDeg(), 4000,0,4000, core_energy);
        // obj.FillHistogram(dirname, Form("%s_core_energy_vs_crystalID",timeflag.c_str()), 48, 0, 48, detMap[cryID], 8192,0,8192, core_energy);
        // obj.FillHistogram(dirname, Form("%s_core_energy_vs_nint",timeflag.c_str()), 12,0,12,hit.NumberOfInteractions(), 8192,0,8192, core_energy);
        // obj.FillHistogram(dirname, Form("%s_core_energy_vs_chi2",timeflag.c_str()), 8192,0,8192, core_energy);
        // obj.FillHistogram(dirname, Form("%s_gretina_map",timeflag.c_str()),720,0,360,phi*TMath::RadToDeg(),360,0,180,theta*TMath::RadToDeg());
        // obj.FillHistogram(dirname, "pad_vs_gretina_timestamps_t0",360,0,3600,(timestamp-timeZero)/TIMESCALE,10,0,10,hit.GetPad());
        
        if (hit.NumberOfInteractions() > 1){
          int nXiPlotEbins = 4096;
          double XiPlotLo = 0;
          double XiPlotHi = 4096;
          // obj.FillHistogram(dirname, Form("%s_core_energy_nInt>1",timeflag.c_str()), nXiPlotEbins,XiPlotLo,XiPlotHi,core_energy);
          obj.FillHistogram(dirname, Form("%s_core_energy_%s",timeflag.c_str(),polCode(hit.GetXiSimple(nullptr)).c_str()),8192,0,8192,core_energy);
          obj.FillHistogram(dirname, Form("%s_core_energy_main_%s",timeflag.c_str(),polCode(hitMain.GetXiSimple(nullptr)).c_str()),8192,0,8192,core_energy);
          // int myFP = 0;
          // int mySP = 1;
          // comptonSortTest(hit,myFP,mySP);
          // double xi = hit.GetXi(nullptr,myFP,mySP);
          // bool goodXi;
          // double xi = hit.GetXi(nullptr);
          bool perp = hit.GetXiSimple(nullptr)*TMath::RadToDeg() > 60;
          bool para = hit.GetXiSimple(nullptr)*TMath::RadToDeg() < 30;
          if (para) obj.FillHistogram(dirname, Form("%s_core_energy_para",timeflag.c_str()),8192,0,8192,core_energy);
          if (perp) obj.FillHistogram(dirname, Form("%s_core_energy_perp",timeflag.c_str()),8192,0,8192,core_energy);
          // if (perp) {
          //   obj.FillHistogram(dirname, Form("%s_core_energy_perp",timeflag.c_str()),8192,0,8192,core_energy);
          //   // if (hit.NumberOfInteractions() < 5) 
          //     // obj.FillHistogram(dirname, Form("%s_core_energy_perp_%d",timeflag.c_str(),hit.NumberOfInteractions()),8192,0,8192,core_energy);
          // }
          // if (para) {
          //   obj.FillHistogram(dirname, Form("%s_core_energy_para",timeflag.c_str()),8192,0,8192,core_energy);
          //   // if (hit.NumberOfInteractions() < 5) 
          //     // obj.FillHistogram(dirname, Form("%s_core_energy_para_%d",timeflag.c_str(),hit.NumberOfInteractions()),8192,0,8192,core_energy);
          // }
          // obj.FillHistogram(dirname, Form("%s_energy_vs_xi_90",timeflag.c_str()),90,0,90,hit.GetXiSimple(nullptr)*TMath::RadToDeg(),2048,0,2048,core_energy);
          obj.FillHistogram(dirname, Form("%s_energy_vs_xi",timeflag.c_str()),180,0,180,hit.GetXi(nullptr)*TMath::RadToDeg(),2048,0,2048,core_energy);
          obj.FillHistogram(dirname, Form("%s_energy_vs_xi_main",timeflag.c_str()),180,0,180,hitMain.GetXi(nullptr)*TMath::RadToDeg(),2048,0,2048,core_energy);
          // obj.FillHistogram(dirname, Form("%s_energy_vs_alpha",timeflag.c_str()),3000,0,30,hitMain.GetAlpha()*TMath::RadToDeg(),2048,0,2048,core_energy);
          // obj.FillHistogram(dirname, Form("%s_energy_xi_bin_%d",timeflag.c_str(),(int) (xi/TMath::Pi()*5)),2048,0,2048,core_energy);
          // obj.FillHistogram("polarization", Form("%s_energy_vs_xi_%d",timeflag.c_str(),cryID),180,0,TMath::Pi(),xi,nXiPlotEbins,XiPlotLo,XiPlotHi,core_energy);
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
