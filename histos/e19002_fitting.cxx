#include "TRuntimeObjects.h"

#include <iostream>
#include <map>
#include <cstdio>

#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TObject.h>
#include <TLine.h>

#include "TGretina.h"
#include "TS800.h"
#include "TBank29.h"
#include "TS800.h"
#include "GCutG.h"
#include "TS800Sim.h"
#include "TGretSim.h"

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

bool efficiencyCorrection(TRandom3 *rand, const TGretinaHit &hit1, int nNeighborHits=-1){
  double thresh_param1 = GValue::Value(Form("DET%i_THRESH1",detMap[hit1.GetCrystalId()]));
  double thresh_param2 = GValue::Value(Form("DET%i_THRESH2",detMap[hit1.GetCrystalId()]));
  double EfficiencyCorrection = (1+TMath::TanH((hit1.GetCoreEnergy()-thresh_param1)/thresh_param2))/2;
  double_t rDraw = rand->Uniform();
  bool value = rDraw <= EfficiencyCorrection;
  
  for (int n = 0; n < nNeighborHits; n++){
    value = value && efficiencyCorrection(rand,hit1.GetNeighbor(n));
  }

  return value;
}

double GetAfp(double crdc_1_x,double  crdc_2_x){
  return TMath::ATan( (crdc_2_x - crdc_1_x)/1073.0 );
}

int gates_loaded=0;
GCutG *crystal_xy, *crystal_zy = nullptr;
// extern "C" is needed to prevent name mangling.
// The function signature must be exactly as shown here,
//   or else bad things will happen.
extern "C"
void MakeHistograms(TRuntimeObjects& obj) {
//  InitMap();
  TS800Sim *s800sim = obj.GetDetector<TS800Sim>();
  TGretSim *gretsim = obj.GetDetector<TGretSim>();
  TGretina *gretina = obj.GetDetector<TGretina>();
  
  TList *gates = &(obj.GetGates());
  if (gates_loaded!=gates->GetSize()){
    TIter iter(gates);
    while(TObject *obj = iter.Next()) {
      GCutG *gate = (GCutG*)obj;
      std::string tag = gate->GetTag();
      if(!tag.compare("crystal-xy")) {
        crystal_xy = gate;
      } 
      if(!tag.compare("crystal-border-zy")) {
        crystal_zy = gate;
      }
      gates_loaded++;
    }
  }
  
  TList    *list    = &(obj.GetObjects());
  int numobj = list->GetSize();

  static TRandom3 *rand_gen = new TRandom3(59953614);

  bool stopped = false;

  if (!gretina){
    return;
  }
  if (!s800sim || !s800sim->Size()){
    stopped = true;
  }

  std::string dirname = "gretsim";
  TVector3 track;
  double yta=0;
  if (!stopped){
    track = s800sim->Track(0,0);
    yta = s800sim->GetS800SimHit(0).GetYTA();
  }

  TGretSimHit simHit = gretsim->GetGretinaSimHit(0);
  // obj.FillHistogram("simcheck","X",1000,-5,5,simHit.GetX());
  // obj.FillHistogram("simcheck","Y",1000,-5,5,simHit.GetY());
  // obj.FillHistogram("simcheck","Z",1000,-5,5,simHit.GetZ());
  // obj.FillHistogram("simcheck","phi",1000,-180,720,simHit.GetPhi()*TMath::RadToDeg());
  // obj.FillHistogram("simcheck","theta",1000,-360,360,simHit.GetTheta()*TMath::RadToDeg());
  // obj.FillHistogram("simcheck","beta",300,0.2,0.6,simHit.GetBeta());
  // obj.FillHistogram("simcheck","energy",2000,0,2000,simHit.GetEn());

  double gammaEn = GValue::Value("FEP_EN");
  if (std::isnan(gammaEn)) gammaEn = simHit.GetEn();
  double beta = GValue::Value("BETA");
  double simBeta = simHit.GetBeta();
  bool isFEP = simHit.IsFEP();
  obj.FillHistogram("ucgretina","beta_sim",200,0.3,0.6,simBeta);
  //S800 coordinates
  if (!stopped){
    obj.FillHistogram("s800sim","ata", 600,-0.1,0.1, s800sim->GetS800SimHit(0).GetATA());
    obj.FillHistogram("s800sim","bta", 600,-0.1,0.1, s800sim->GetS800SimHit(0).GetBTA());
    obj.FillHistogram("s800sim","yta", 1000,-0.003,0.003, s800sim->GetS800SimHit(0).GetYTA());
    obj.FillHistogram("s800sim","dta", 1000,-0.5,0.5, s800sim->GetS800SimHit(0).GetDTA());
  } else {
    obj.FillHistogram("gretsim","Sim Energies",10000,0,10000,simHit.GetEn());
  }
  double SIGMA = (2.1*TMath::Exp(-0.1*gammaEn/1000.0) + 60.0*TMath::Exp(-10.2*gammaEn/1000.0));
  if(SIGMA > 3.8) {
    SIGMA = 3.8;
  }
  if(gammaEn < 150.) {
    SIGMA = 6.0;
  }

  //SINGLES
  int gSize = gretina->Size();
  for (int i=0; i < gSize; i++){
    TGretinaHit &hit = gretina->GetGretinaHit(i);

    double core_energy = hit.GetCoreEnergy();
    int cryID = hit.GetCrystalId();
    int ringNum = hit.GetRingNumber();
    // if (cryID == 77) continue;

    double theta = hit.GetTheta();
    double phi = hit.GetPhi();

    TVector3 local_pos(hit.GetLocalPosition(0));
    double smear_x = local_pos.X() + rand_gen->Gaus(0, SIGMA); 
    double smear_y = local_pos.Y() + rand_gen->Gaus(0, SIGMA);
    double smear_z = local_pos.Z() + rand_gen->Gaus(0, SIGMA);
    if (gates && crystal_xy && crystal_zy){
      bool inXY = crystal_xy->IsInside(smear_x,smear_y);
      bool inZY = crystal_zy->IsInside(smear_z,smear_y);
      obj.FillHistogram("positionsmear","inout",2,0,2,int(inXY && inZY));
      if (inXY && inZY){
        obj.FillHistogram("positionsmear",Form("X_vs_Y_smear_inside"),200,-100,100,smear_x,200,-100,100,smear_y); 
        obj.FillHistogram("positionsmear",Form("Z_vs_Y_smear_inside"),200,-100,100,smear_z,200,-100,100,smear_y);
      } else {
        if (!inXY) obj.FillHistogram("positionsmear",Form("X_vs_Y_smear_outside"),200,-100,100,smear_x,200,-100,100,smear_y); 
        if (!inZY) obj.FillHistogram("positionsmear",Form("Z_vs_Y_smear_outside"),200,-100,100,smear_z,200,-100,100,smear_y);
      }
    }
    hit.SetPosition(0,smear_x,smear_y,smear_z); //this resets the positions also in NNaddback

    double energy_track_yta_dta;
    double energy_track;

    if (!stopped){
      energy_track = hit.GetDoppler(beta, &track);
      energy_track_yta_dta = hit.GetDopplerYta(s800sim->AdjustedBeta(beta), yta, &track);
    } 
    else{
      energy_track = energy_track_yta_dta = hit.GetCoreEnergy();
    }

    //efficiency correction
    if (efficiencyCorrection(rand_gen,hit)){
      obj.FillHistogram(dirname,"HitPhi_v_HitTheta",180,0,180,180-theta*TMath::RadToDeg(),360,0,360,phi*TMath::RadToDeg());
      // obj.FillHistogram(dirname,Form("HitPhi_v_HitTheta_Smeared_c%d",detMap[cryID]),180,0,180,180-theta_smear*TMath::RadToDeg(),360,0,360,phi_smear*TMath::RadToDeg());
      obj.FillHistogram(dirname,"CoreEnergy",10000,0,10000,core_energy);

      obj.FillHistogram(dirname,"Theta_vs_Energy_btyd",4000,0,4000,energy_track_yta_dta,100,0,3,hit.GetTheta());

      // obj.FillHistogram(dirname,"EmissionAngle_vs_DetectedAngle",180,0,180,hit.GetTheta()*TMath::RadToDeg(),180,0,180,simHit.GetTheta()*TMath::RadToDeg());
      
      //fitting hists
      obj.FillHistogram(dirname,"gretina_B&T&Y&D",10000,0,10000,energy_track_yta_dta);
      if (gates && crystal_xy && crystal_zy && crystal_xy->IsInside(smear_x,smear_y) && crystal_zy->IsInside(smear_z,smear_y)) 
        obj.FillHistogram(dirname,"gretina_B&T&Y&D_inside",10000,0,10000,energy_track_yta_dta);
      else obj.FillHistogram(dirname,"gretina_B&T&Y&D_outside",10000,0,10000,energy_track_yta_dta);

      //SUMMARY SPECTRUM
      obj.FillHistogram(dirname,"dop_btyd_summary",48,0,48,detMap[cryID],3000,0,3000,energy_track_yta_dta);

      //organization
      // if(detMap[cryID] < 16) {
      //   obj.FillHistogram(dirname,"gretina_B&T_Fwd",10000,0,10000,energy_track);
      // } else {
      //   obj.FillHistogram(dirname,"gretina_B&T_90Deg",10000,0,10000,energy_track);
      // }

      if (isFEP){ //full energy peak event
        obj.FillHistogram(dirname,"gretina_B&T&Y&D_fep",10000,0,10000,energy_track_yta_dta);
            
      //   if(detMap[cryID] < 17) {
      //     obj.FillHistogram(dirname,"gretina_B&T_fep_Fwd",10000,0,10000,energy_track);
      //   } else {
      //     obj.FillHistogram(dirname,"gretina_B&T_fep_90Deg",10000,0,10000,energy_track);
      //   }
      } else {
        obj.FillHistogram(dirname,"gretina_B&T&Y&D_bg",10000,0,10000,energy_track_yta_dta);

      //   if(detMap[cryID] < 17) {
      //     obj.FillHistogram(dirname,"gretina_B&T_bg_Fwd",10000,0,10000,energy_track);
      //   } else {
      //     obj.FillHistogram(dirname,"gretina_B&T_bg_90Deg",10000,0,10000,energy_track);
      //   }
      }
    }
  }

  //NNADDBACK
  //loop over multiplicity
  dirname = "addback";
  for (int n=0; n<4; n++){
    //loop over hits for each multiplicity spectrum
    int nnSize = gretina->NNAddbackSize(n);
    for (int i=0; i < nnSize; i++){

      //get hit and hit data 
      TGretinaHit nnhit = gretina->GetNNAddbackHit(n,i);
      int ringNum = nnhit.GetRingNumber();
      int cryID = nnhit.GetCrystalId();
      double core_energy = nnhit.GetCoreEnergy();
      if (cryID == 77) continue;

      if (efficiencyCorrection(rand_gen,nnhit,nnhit.GetNNeighborHits())){
        double energy_track_yta_dta;
        double energy_track_yta;
        double energy_track;
        
        // TVector3 local_pos(nnhit.GetLocalPosition(0));
        // double smear_x = local_pos.X() + rand_gen->Gaus(0, SIGMA*0.9); 
        // double smear_y = local_pos.Y() + rand_gen->Gaus(0, SIGMA*0.9);
        // double smear_z = local_pos.Z() + rand_gen->Gaus(0, SIGMA*0.9);
        // nnhit.SetPosition(0,smear_x,smear_y,smear_z);

        if (!stopped){
          energy_track = nnhit.GetDoppler(beta, &track);
          energy_track_yta = nnhit.GetDopplerYta(beta, yta, &track);
          energy_track_yta_dta = nnhit.GetDopplerYta(s800sim->AdjustedBeta(beta), yta, &track);
        } else {
          energy_track = energy_track_yta = energy_track_yta_dta = nnhit.GetCoreEnergy();
        }

        //compare uncorrected energy reported from geant with core energy. if FEP they should be the same
        bool isNNFEP = fabs(simHit.GetEn() - nnhit.GetCoreEnergy()) < 1.5;

        char *multiplicity = Form("%d",n);
        if (n == 3) multiplicity = Form("g");
        obj.FillHistogram(dirname, Form("gretina_n%s_CoreEnergy",multiplicity), 10000,0,10000, core_energy);
        obj.FillHistogram(dirname, Form("gretina_n%s_B&T",multiplicity), 10000,0,10000, energy_track);
        obj.FillHistogram(dirname, Form("gretina_n%s_B&T&Y",multiplicity), 10000,0,10000, energy_track_yta);
        obj.FillHistogram(dirname, Form("gretina_n%s_B&T&Y&D",multiplicity), 10000,0,10000, energy_track_yta_dta);
        if (isNNFEP){
          obj.FillHistogram(dirname, Form("gretina_n%s_B&T&Y&D_fep",multiplicity), 10000,0,10000, energy_track_yta_dta);
        } else {
          obj.FillHistogram(dirname, Form("gretina_n%s_B&T&Y&D_bg",multiplicity), 10000,0,10000, energy_track_yta_dta);
        }
        
        //exclude the ng spectrum (n==3)
        if (n < 3){
          obj.FillHistogram(dirname, "gretina_ab_CoreEnergy", 10000,0,10000, core_energy);
          obj.FillHistogram(dirname, "gretina_ab_B&T", 10000,0,10000, energy_track);
          obj.FillHistogram(dirname, "gretina_ab_B&T&Y", 10000,0,10000, energy_track_yta);
          obj.FillHistogram(dirname, "gretina_ab_B&T&Y&D", 10000,0,10000, energy_track_yta_dta);
          if (isNNFEP){
            obj.FillHistogram(dirname, "gretina_ab_B&T&Y&D_fep", 10000,0,10000, energy_track_yta_dta);
          } else {
            obj.FillHistogram(dirname, "gretina_ab_B&T&Y&D_bg", 10000,0,10000, energy_track_yta_dta);
          }
        }

        if (n == 1) {
          //POLARIZATION
          std::string swaptype = "pol";
          for (int t=0; t < 1; t++){
            if (t==1) {
              swaptype = "pol_swapped";
              TGretinaHit swap = nnhit.GetNeighbor();
              swap.NNAdd(nnhit.GetInitialHit());
              energy_track = swap.GetDoppler(beta, &track);
              energy_track_yta = swap.GetDopplerYta(beta, yta, &track);
              energy_track_yta_dta = swap.GetDopplerYta(s800sim->AdjustedBeta(beta), yta, &track);
            }
            std::string polColor = "blank";
            if (PairHit(nnhit,redPairs)) polColor = "red";
            else if (PairHit(nnhit,goldPairs)) polColor = "gold";
            else if (PairHit(nnhit,bluePairs)) polColor = "blue";

            std::string quadType = "blank";
            if (PairHit(nnhit,OneQuadPlus)) quadType = "qd1+";
            else if (PairHit(nnhit,OneQuadDefault)) quadType = "qd1";
            else if (PairHit(nnhit,TwoQuadPairs)) quadType = "qd2";

            if ( polColor.compare("blank") != 0 ){
              obj.FillHistogram(dirname,Form("gretina_%s_%s_B&T&Y&D",swaptype.c_str(),polColor.c_str()), 8192,0,8192, energy_track_yta_dta);
              if (isNNFEP) obj.FillHistogram(dirname,Form("gretina_%s_%s_B&T&Y&D_fep",swaptype.c_str(),polColor.c_str()), 8192,0,8192, energy_track_yta_dta);
              else obj.FillHistogram(dirname,Form("gretina_%s_%s_B&T&Y&D_bg",swaptype.c_str(),polColor.c_str()), 8192,0,8192, energy_track_yta_dta);
              obj.FillHistogram(dirname,Form("gretina_%s_%s",swaptype.c_str(),quadType.c_str()), 8192,0,8192, energy_track_yta_dta);
            }
          }
        }
      }
    }
  }

  if(numobj!=list->GetSize()){
    list->Sort();
  }
}
