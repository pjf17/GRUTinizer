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

std::map<int,int> detMap = {{26,0}, {30,1}, {34,2}, {38,3}, {25,4}, {29,5}, {33,6}, {37,7}, {27,8}, {31,9}, {35,10}, {39,11},
			    {24,12},{28,13},{32,14},{36,15},{47,16},{51,17},{59,18},{63,19},{67,20},{71,21}, {79,22}, {44,23},
			    {50,24}, {58,25}, {60,26}, {66,27}, {68,28}, {76,29}, {46,30}, {48,31}, {56,32}, {62,33},
			    {64,34}, {70,35}, {78,36}, {45,37}, {49,38}, {57,39}, {61,40}, {65,41}, {69,42}, {77,43}};

std::vector<std::pair<int,int>> redPairs = {
  //1qAB
  std::make_pair(48,49),
  std::make_pair(50,51),
  std::make_pair(56,57),
  std::make_pair(59,58),
  std::make_pair(64,65),
  std::make_pair(66,67),
  //1qBB
  std::make_pair(46,44),
  std::make_pair(60,62),
  std::make_pair(70,68),
  std::make_pair(78,76),
  //2qAA
  std::make_pair(51,47),
  std::make_pair(63,67),
  std::make_pair(57,61),
  std::make_pair(69,65),
  //2qBB
  std::make_pair(60,58),
  std::make_pair(66,68),
  std::make_pair(64,62),
  std::make_pair(46,48)
};

std::vector<std::pair<int,int>> goldPairs = {
  std::make_pair(49,50),
  std::make_pair(56,59),
  std::make_pair(57,58),
  std::make_pair(64,67),
  std::make_pair(65,66),
  std::make_pair(44,45),
  std::make_pair(46,47),
  std::make_pair(48,51),
  std::make_pair(61,60),
  std::make_pair(62,63),
  std::make_pair(69,68),
  std::make_pair(70,71),
  std::make_pair(78,79)
};

std::vector<std::pair<int,int>> bluePairs = {
  //1qAB
  std::make_pair(44,47),
  std::make_pair(45,46),
  std::make_pair(60,63),
  std::make_pair(61,62),
  std::make_pair(68,71),
  std::make_pair(76,79),
  std::make_pair(69,70),
  //1qBB
  std::make_pair(48,50),
  std::make_pair(56,58),
  std::make_pair(64,66),
  //2qAB
  std::make_pair(57,60),
  std::make_pair(62,67),
  std::make_pair(65,68),
  std::make_pair(46,51)
};

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
  std::make_pair(65,68)
};

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
  std::make_pair(65,66)
};

std::vector<std::pair<int,int>> OneQuadDefault = {
  std::make_pair(44,45),
  std::make_pair(46,47),
  std::make_pair(44,46),
  std::make_pair(48,51),
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
  std::make_pair(78,79)
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

// extern "C" is needed to prevent name mangling.
// The function signature must be exactly as shown here,
//   or else bad things will happen.
extern "C"
void MakeHistograms(TRuntimeObjects& obj) {
//  InitMap();
  TS800Sim *s800sim = obj.GetDetector<TS800Sim>();
  TGretSim *gretsim = obj.GetDetector<TGretSim>();
  TGretina *gretina = obj.GetDetector<TGretina>();
  TList    *list    = &(obj.GetObjects());
  int numobj = list->GetSize();

  TRandom3 *rand_gen = new TRandom3(59953614);

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
  double gammaEn = GValue::Value("FEP_EN");
  double beta = GValue::Value("BETA");
  double simBeta = simHit.GetBeta();
  bool isFEP = simHit.IsFEP();
  obj.FillHistogram("ucgretina","beta_sim",200,0.2,0.4,simBeta);
  //S800 coordinates
  if (!stopped){
    obj.FillHistogram("s800sim","ata", 600,-0.1,-0.1, s800sim->GetS800SimHit(0).GetATA());
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
    double theta = hit.GetTheta();
    double phi = hit.GetPhi();
    int cryID = hit.GetCrystalId();
    int number = detMap[cryID];
    if (cryID == 77) continue;

    // TVector3 local_pos(hit.GetLocalPosition(0));
    // double smear_x = local_pos.X() + rand_gen->Gaus(0, SIGMA);
    // double smear_y = local_pos.Y() + rand_gen->Gaus(0, SIGMA);
    // double smear_z = local_pos.Z() + rand_gen->Gaus(0, SIGMA);
    // hit.SetPosition(0,smear_x,smear_y,smear_z);

    double energy_track_yta_dta;
    double energy_track_yta;
    double energy_track;
    double energy_beta;

    if (!stopped){
      energy_beta = hit.GetDoppler(beta);
      energy_track = hit.GetDoppler(beta, &track);
      energy_track_yta = hit.GetDopplerYta(beta, yta, &track);
      energy_track_yta_dta = hit.GetDopplerYta(s800sim->AdjustedBeta(beta), yta, &track);
    } 
    else{
      energy_track = energy_track_yta = energy_track_yta_dta = hit.GetDoppler(beta);
    }

    double energy_track_yta_dta_sim;
    double energy_track_yta_sim;
    double energy_track_sim;

    if (!stopped){
      energy_track_sim = hit.GetDoppler(simBeta, &track);
      energy_track_yta_sim = hit.GetDopplerYta(simBeta, yta, &track);
      energy_track_yta_dta_sim = hit.GetDopplerYta(s800sim->AdjustedBeta(simBeta), yta, &track);
    } 
    else{
      energy_track_sim = energy_track_yta_sim = energy_track_yta_dta_sim = hit.GetDoppler(simBeta);
    }

    //efficiency correction
    if (efficiencyCorrection(rand_gen,hit)){
      obj.FillHistogram(dirname,"HitTheta_v_HitPhi",360,0,360,theta*TMath::RadToDeg(),360,0,360,phi*TMath::RadToDeg());
      obj.FillHistogram(dirname,"CoreEnergy",10000,0,10000,core_energy);

      obj.FillHistogram(dirname,"gretina_summary_B&T",36,1,37,number,8000,0,4000,energy_track);
      obj.FillHistogram(dirname,"gretina_summary_B&T&Y&D",36,1,37,number,8000,0,4000,energy_track_yta_dta);
      
      //fitting hists
      obj.FillHistogram(dirname,"gretina_B",10000,0,10000,energy_beta);
      obj.FillHistogram(dirname,"gretina_B&T",10000,0,10000,energy_track);
      obj.FillHistogram(dirname,"gretina_B&T&Y",10000,0,10000,energy_track_yta);
      obj.FillHistogram(dirname,"gretina_B&T&Y&D",10000,0,10000,energy_track_yta_dta);

      obj.FillHistogram(dirname,"gretina_B&T_sim",10000,0,10000,energy_track_sim);
      obj.FillHistogram(dirname,"gretina_B&T&Y_sim",10000,0,10000,energy_track_yta_sim);
      obj.FillHistogram(dirname,"gretina_B&T&Y&D_sim",10000,0,10000,energy_track_yta_dta_sim);

      //organization
      // if(detMap[cryID] < 17) {
      //   obj.FillHistogram(dirname,"gretina_B&T_Fwd",10000,0,10000,energy_track);
      // } else {
      //   obj.FillHistogram(dirname,"gretina_B&T_90Deg",10000,0,10000,energy_track);
      // }

      if (isFEP){ //full energy peak event
        obj.FillHistogram(dirname,"gretina_B&T_fep",10000,0,10000,energy_track);
            
      //   if(detMap[cryID] < 17) {
      //     obj.FillHistogram(dirname,"gretina_B&T_fep_Fwd",10000,0,10000,energy_track);
      //   } else {
      //     obj.FillHistogram(dirname,"gretina_B&T_fep_90Deg",10000,0,10000,energy_track);
      //   }
      } else {
        obj.FillHistogram(dirname,"gretina_B&T_bg",10000,0,10000,energy_track);

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
    int nnSize = gretina->NNAddbackSize(n,false);
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
        if (!stopped){
          energy_track = nnhit.GetDoppler(beta, &track);
          energy_track_yta = nnhit.GetDopplerYta(beta, yta, &track);
          energy_track_yta_dta = nnhit.GetDopplerYta(s800sim->AdjustedBeta(beta), yta, &track);
        } else {
          energy_track = energy_track_yta = energy_track_yta_dta = nnhit.GetDoppler(beta);
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
          obj.FillHistogram(dirname, Form("gretina_n%s_B&T_fep",multiplicity), 10000,0,10000, energy_track);
        } else {
          obj.FillHistogram(dirname, Form("gretina_n%s_B&T_bg",multiplicity), 10000,0,10000, energy_track);
        }
        
        //exclude the ng spectrum (n==3)
        if (n < 3){
          obj.FillHistogram(dirname, "gretina_ab_CoreEnergy", 10000,0,10000, core_energy);
          obj.FillHistogram(dirname, "gretina_ab_B&T", 10000,0,10000, energy_track);
          obj.FillHistogram(dirname, "gretina_ab_B&T&Y", 10000,0,10000, energy_track_yta);
          obj.FillHistogram(dirname, "gretina_ab_B&T&Y&D", 10000,0,10000, energy_track_yta_dta);
          if (isNNFEP){
            obj.FillHistogram(dirname, "gretina_ab_B&T_fep", 10000,0,10000, energy_track);
          } else {
            obj.FillHistogram(dirname, "gretina_ab_B&T_bg", 10000,0,10000, energy_track);
          }
        }

        if (n == 1) {
          //POLARIZATION
          if ( PairHit(nnhit,redPairs) ){
            obj.FillHistogram(dirname,"gretina_pol_red_B&T&Y&D", 10000,0,10000, energy_track_yta_dta);
            if(isNNFEP){
              obj.FillHistogram(dirname,"gretina_pol_red_B&T_fep", 10000,0,10000, energy_track);
            } else {
              obj.FillHistogram(dirname,"gretina_pol_red_B&T_bg", 10000,0,10000, energy_track);
            }
          }
          if ( PairHit(nnhit,goldPairs) ){
            obj.FillHistogram(dirname,"gamma_pol_gold_B&T&Y&D", 10000,0,10000, energy_track_yta_dta);
            if(isNNFEP){
              obj.FillHistogram(dirname,"gamma_pol_gold_B&T_fep", 10000,0,10000, energy_track);
            } else {
              obj.FillHistogram(dirname,"gamma_pol_gold_B&T_bg", 10000,0,10000, energy_track);
            }
          }
          if ( PairHit(nnhit,bluePairs) ){
            obj.FillHistogram(dirname,"gamma_pol_blue_B&T&Y&D", 10000,0,10000, energy_track_yta_dta);
            if(isNNFEP){
              obj.FillHistogram(dirname,"gamma_pol_blue_B&T_fep", 10000,0,10000, energy_track);
            } else {
              obj.FillHistogram(dirname,"gamma_pol_blue_B&T_bg", 10000,0,10000, energy_track);
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
