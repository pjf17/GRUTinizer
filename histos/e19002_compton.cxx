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

int scatterType(const TGretinaHit& abhit){
  int cryId1 = abhit.GetCrystalId();
  int cryId2 = abhit.GetNeighbor().GetCrystalId();
  int out = 0; //one A and one B
  if (cryId1%2 == 0 && cryId2%2 == 0) out = -1; //both are B type
  else if (cryId1%2 != 0 && cryId2%2 != 0) out = 1; //both are A type
  return out;
}

double calcScatteringAngle(const TGretinaHit &hit1, const TGretinaHit &hit2){
  TVector3 pos1 = hit1.GetPosition();
  TVector3 diff = hit2.GetPosition() - pos1;
  return TMath::Cos(diff.Angle(pos1));
}

double calcComptonAngle(double E1, double E2){
  double argument = 1 - 511/E2 + 511/(E1 + E2);
  // if (argument < -1){
  //   int phase = std::floor(std::fabs(argument)); 
  //   return TMath::ACos(argument + phase)*180/TMath::Pi() + phase*180;
  // } else {
  //   return TMath::ACos(argument)*180/TMath::Pi();
  // }
  return argument;
}

double GetAfp(double crdc_1_x,double  crdc_2_x){
  return TMath::ATan( (crdc_2_x - crdc_1_x)/1073.0 );
}

std::vector<GCutG*> gates_2D = {};
int gates_loaded=0;

void LoadGates(TRuntimeObjects &obj){
  TList *gates = &(obj.GetGates());
  TIter iter(gates);
  std::cout << "loading gates:" <<std::endl;
  while(TObject *obj = iter.Next()) {
    GCutG *gate = (GCutG*)obj;
    // std::string tag = gate->GetTag();
    gates_2D.push_back(gate);
    gates_loaded++;
  }
  std::cout << "gates size: " << gates_2D.size() << std::endl;
}

int THETAMAX = 360;

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

  TList *gates = &(obj.GetGates());
  if(gates_loaded!=gates->GetSize()) {
    LoadGates(obj);
  }

  bool stopped = false;

  if (!gretina){
    return;
  }
  if (!s800sim || !s800sim->Size()){
    stopped = true;
  }

  std::string dirname("basicsim");
  TVector3 track;
  double yta;
  if (!stopped){
    track = s800sim->Track(0,0);
    yta = s800sim->GetS800SimHit(0).GetYTA();
  }
    
  if (gretina){
    double FEP = GValue::Value("FEP_EN");
    TGretSimHit simHit = gretsim->GetGretinaSimHit(0);
    double gammaEn = simHit.GetEn();
    double gammaEnDop = simHit.GetDoppler();
    double gammaBeta = simHit.GetBeta();

    //SINGLES
    int gSize = gretina->Size();
    for (int i=0; i < gSize; i++){
      TGretinaHit &hit = gretina->GetGretinaHit(i);

      double energy_corrected = hit.GetDopplerYta(s800sim->AdjustedBeta(GValue::Value("BETA")), yta, &track);
      double energy = hit.GetDoppler(GValue::Value("BETA"));
      double core_energy = hit.GetCoreEnergy();
      double theta = hit.GetTheta();
      double phi = hit.GetPhi();
      int cryID = hit.GetCrystalId();
      
      obj.FillHistogram(dirname, "core_energy_prompt", 8192,0,8192, core_energy);
      obj.FillHistogram(dirname, "gamma_singles_prompt", 8192,0,8192, energy);
      obj.FillHistogram(dirname, "gamma_corrected_singles_prompt", 8192,0,8192, energy_corrected);
      obj.FillHistogram(dirname, "gamma_corrected_vs_theta_prompt", 8192,0,8192, energy_corrected, 100, 0, 2.5, theta);
      obj.FillHistogram(dirname, "gamma_corrected_vs_crystalID_prompt", 56, 24, 80, cryID, 8192,0,8192, energy_corrected);
      obj.FillHistogram(dirname, "core_energy_vs_theta_prompt", 8192,0,8192, hit.GetCoreEnergy(), 100, 0, 2.5, theta);
      
      //crystal angles
      obj.FillHistogram(dirname, Form("theta_vs_crystalID_ring%d",gretina->GetRingNumber(cryID)),56,24,80,cryID,360,0,180,hit.GetThetaDeg());
      obj.FillHistogram(dirname, "gretina_theta_vs_phi",720,0,360,phi*TMath::RadToDeg(),360,0,180, theta*TMath::RadToDeg());

    }

    //NNADDBACK
    int nnSize = gretina->NNAddbackSize();
    for (int i=0; i < nnSize; i++){
      dirname = "basicsim";

      //get hit and hit data 
      TGretinaHit nnhit = gretina->GetNNAddbackHit(i);
      int cryID = nnhit.GetCrystalId();
      int ringNum = nnhit.GetRingNumber();
      double nnEnergy_corrected = nnhit.GetDopplerYta(s800sim->AdjustedBeta(GValue::Value("BETA")), yta, &track);
      int n = nnhit.GetABDepth();
      
      //exclude the ng spectrum (n==3)
      if (n < 3){
        obj.FillHistogram(dirname, "gamma_corrected_addback_prompt", 8192,0,8192, nnEnergy_corrected);
      }

      if (n == 1) {
        //N1 Neighbor correlations
        TGretinaHit nnhit1 = nnhit.GetInitialHit();
        TGretinaHit nnhit2 = nnhit.GetNeighbor();
        int layer1 = nnhit1.GetSegmentLayer();
        int layer2 = nnhit2.GetSegmentLayer();

        int hit1Ints = nnhit1.NumberOfInteractions();
        int hit2Ints = nnhit2.NumberOfInteractions();

        //Number of interactions
        obj.FillHistogram("interactions","nnhit1_interactions",12,0,12,hit1Ints);
        obj.FillHistogram("interactions","nnhit2_interactions",12,0,12,hit2Ints);
        
        //compton analysis
        double E1 = nnhit1.GetCoreEnergy();
        double E2 = nnhit2.GetCoreEnergy();
        double Vtheta = calcScatteringAngle(nnhit1,nnhit2); 
        double Vphi = calcScatteringAngle(nnhit2,nnhit1);
        double Etheta = calcComptonAngle(E1,E2);
        double Ephi = calcComptonAngle(E2,E1);


        if (abs(gammaEn - (E1+E2)) < 1.5) {
          dirname = "FEP_gated";
          obj.FillHistogram(dirname,"gammaEn_vs_rawSum",500,0,500,E1+E2,500,0,500,gammaEn);

          double thetaDiff = nnhit1.GetThetaDeg() - nnhit2.GetThetaDeg();
          double dE = gammaEn * gammaBeta * TMath::Sin(nnhit1.GetTheta()) / TMath::Sqrt(1 - gammaBeta*gammaBeta) * thetaDiff*TMath::DegToRad();
          if (abs(thetaDiff) < 1.5) dE = 0;

          //POLARIZATION
          if ( PairHit(nnhit,redPairs) ){
            obj.FillHistogram(dirname,"gamma_corrected_addback_prompt_red_pair", 1200,0,1200, nnEnergy_corrected);
            obj.FillHistogram(dirname, Form("predicted_vs_actual_red_ring%02d",ringNum),1200,0,1200,nnEnergy_corrected,1200,0,1200, gammaEnDop + dE);
            obj.FillHistogram(dirname, Form("gamma_red_ring%02d",ringNum), 1200,0,1200, nnEnergy_corrected);
            obj.FillHistogram(dirname, Form("gamma_calc_red_ring%02d",ringNum), 1200,0,1200, gammaEnDop + dE);
            obj.FillHistogram(dirname, Form("theta_vs_crystalID_red_pair_ring%02d",ringNum),56,24,80,nnhit1.GetCrystalId(),360,0,180,nnhit1.GetThetaDeg());
            obj.FillHistogram(dirname, Form("theta_vs_crystalID_red_pair_ring%02d",ringNum),56,24,80,nnhit2.GetCrystalId(),360,0,180,nnhit2.GetThetaDeg());
            obj.FillHistogram(dirname, Form("thetadiff_red_ring%02d",ringNum),240,-60,60,thetaDiff);
          }

          if ( PairHit(nnhit,goldPairs) ){
            obj.FillHistogram(dirname,"gamma_corrected_addback_prompt_gold_pair", 1200,0,1200, nnEnergy_corrected);
            obj.FillHistogram(dirname, Form("predicted_vs_actual_gold_ring%02d",ringNum),1200,0,1200,nnEnergy_corrected,1200,0,1200, gammaEnDop + dE);
            obj.FillHistogram(dirname,Form("gamma_gold_ring%02d",ringNum), 1200,0,1200, nnEnergy_corrected);
            obj.FillHistogram(dirname,Form("gamma_calc_gold_ring%02d",ringNum), 1200,0,1200, gammaEnDop + dE);
            obj.FillHistogram(dirname, Form("theta_vs_crystalID_gold_pair_ring%02d",ringNum),56,24,80,nnhit1.GetCrystalId(),360,0,180,nnhit1.GetThetaDeg());
            obj.FillHistogram(dirname, Form("theta_vs_crystalID_gold_pair_ring%02d",ringNum),56,24,80,nnhit2.GetCrystalId(),360,0,180,nnhit2.GetThetaDeg());
            obj.FillHistogram(dirname, Form("thetadiff_gold_ring%02d",ringNum),240,-60,60,thetaDiff);
          }

          if ( PairHit(nnhit,bluePairs) ){
            obj.FillHistogram(dirname,"gamma_corrected_addback_prompt_blue_pair", 1200,0,1200, nnEnergy_corrected);
            obj.FillHistogram(dirname, Form("predicted_vs_actual_blue_ring%02d",ringNum),1200,0,1200,nnEnergy_corrected,1200,0,1200, gammaEnDop + dE);
            obj.FillHistogram(dirname,Form("gamma_blue_ring%02d",ringNum), 1200,0,1200, nnEnergy_corrected);
            obj.FillHistogram(dirname,Form("gamma_calc_blue_ring%02d",ringNum), 1200,0,1200, gammaEnDop + dE);
            obj.FillHistogram(dirname, Form("theta_vs_crystalID_blue_pair_ring%02d",ringNum),56,24,80,nnhit1.GetCrystalId(),360,0,180,nnhit1.GetThetaDeg());
            obj.FillHistogram(dirname, Form("theta_vs_crystalID_blue_pair_ring%02d",ringNum),56,24,80,nnhit2.GetCrystalId(),360,0,180,nnhit2.GetThetaDeg());
            obj.FillHistogram(dirname, Form("thetadiff_blue_ring%02d",ringNum),240,-60,60,thetaDiff);
          }
          
          // obj.FillHistogram("scattering-angles", Form("pos_angles_n%d_ring%02d_prompt",n,ringNum),THETAMAX,0,THETAMAX,Vtheta);
          // obj.FillHistogram("scattering-angles", Form("eng_angles_n%d_ring%02d_prompt",n,ringNum),THETAMAX,0,THETAMAX,Etheta);

          if (layer1 - layer2 > 1){
            std::swap(nnhit1,nnhit2);
          }
          nnhit1.NNAdd(nnhit2);
          double improvedEnergy = nnhit1.GetDopplerYta(s800sim->AdjustedBeta(GValue::Value("BETA")), yta, &track);
          
          //make swapped energy
          nnhit1 = nnhit.GetInitialHit();
          nnhit2 = nnhit.GetNeighbor();
          nnhit2.NNAdd(nnhit1);
          double swappedEnergy = nnhit2.GetDopplerYta(s800sim->AdjustedBeta(GValue::Value("BETA")), yta, &track);

          if (ringNum == 1 || ringNum == 4 || ringNum == 12){
            obj.FillHistogram("layer-differences", Form("gamma_corrected_n%d_vs_L1-L2_ring%02d_prompt",n,ringNum),12,-6,6,layer1-layer2, 1600,0,1600,nnEnergy_corrected);
            obj.FillHistogram("layer-differences", Form("gamma_corrected_n%d_swapped_vs_L1-L2_ring%02d_prompt",n,ringNum),12,-6,6,layer1-layer2, 1600,0,1600,swappedEnergy);
            obj.FillHistogram("layer-differences", Form("gamma_corrected_n%d_improve_vs_L1-L2_ring%02d_prompt",n,ringNum),12,-6,6,layer1-layer2, 1600,0,1600,improvedEnergy);
            obj.FillHistogram("layer-differences", Form("E1_vs_L1-L2_ring%02d_prompt",ringNum),12,-6,6,layer1-layer2, 1400,0,1400,E1);
            obj.FillHistogram("layer-differences", Form("E2_vs_L1-L2_ring%02d_prompt",ringNum),12,-6,6,layer1-layer2, 1200,0,1200,E2);
          }

          //GATED STUFF
          // std::vector<std::pair<std::string,bool>> EnergyGates;
          // // EnergyGates.push_back(std::make_pair("defaultPeakCut",fabs(nnEnergy_corrected - FEP) < 1.5 * TMath::Sqrt(FEP)));
          // EnergyGates.push_back(std::make_pair("swappedPeakCut",fabs(swappedEnergy - FEP) < 1.2 * TMath::Sqrt(FEP)));
          // int nEnergyGates = (int) EnergyGates.size();

          std::vector<std::pair<std::string,bool>> InteractionGates;
          // InteractionGates.push_back(std::make_pair("E1_E2_int>1",hit1Ints>1 && hit2Ints>1));
          InteractionGates.push_back(std::make_pair("E1_E2_int=1",hit1Ints==1 && hit2Ints==1));
          int nInteractionGates = (int) InteractionGates.size();

          // std::vector<std::pair<std::string,bool>> LineGates;
          // LineGates.push_back(std::make_pair("thetaCut",Etheta - Vtheta < 0.15));
          // LineGates.push_back(std::make_pair("phiCut",Ephi - Vphi < 0.05));

          // for (int eg=0; eg < nEnergyGates; eg++){
          //   if (EnergyGates[eg].second){
          //     if (InteractionGates[0].second){
          //       if (LineGates[1].second){
          //         dirname = EnergyGates[eg].first + "_" + LineGates[1].first + "_" + InteractionGates[0].first;
          //         obj.FillHistogram(dirname, Form("%s_gated_Ephi_vs_Vphi",dirname.c_str()),200,-1,1,Vphi,1000,-1,25,Ephi);
          //         obj.FillHistogram(dirname, Form("%s_gated_Etheta_vs_Vtheta",dirname.c_str()),200,-1,1,Vtheta,1000,-25,1,Etheta);
          //       }

          //       if (LineGates[0].second){
          //         dirname = EnergyGates[eg].first + "_" + LineGates[0].first + "_" + InteractionGates[0].first;
          //         obj.FillHistogram(dirname, Form("%s_gated_Ephi_vs_Vphi",dirname.c_str()),200,-1,1,Vphi,1000,-1,25,Ephi);
          //         obj.FillHistogram(dirname, Form("%s_gated_Etheta_vs_Vtheta",dirname.c_str()),200,-1,1,Vtheta,1000,-25,1,Etheta);
          //       }

          //       if (LineGates[0].second){
          //         if (LineGates[1].second){
          //           dirname = EnergyGates[eg].first + "_" + LineGates[0].first + "_" + LineGates[1].first + "_" + InteractionGates[0].first;
          //           //Hit vs Angle
          //           obj.FillHistogram(dirname, Form("%s_gated_E2_vs_Vtheta",dirname.c_str()),200,-1,1,Vtheta,1200,0,1200,E2);
          //           obj.FillHistogram(dirname, Form("%s_gated_E1_vs_Vtheta",dirname.c_str()),200,-1,1,Vtheta,1600,0,1600,E1);

          //           //Hit vs Hit
          //           obj.FillHistogram(dirname, Form("%s_gated_E1_vs_E2",dirname.c_str()),1200,0,1200,E2,1600,0,1600,E1);

          //           //Angle calculation correlation
          //           obj.FillHistogram(dirname, Form("%s_gated_Etheta_vs_Vtheta",dirname.c_str()),200,-1,1,Vtheta,1000,-25,1,Etheta);
          //           obj.FillHistogram(dirname, Form("%s_gated_Ephi_vs_Vphi",dirname.c_str()),200,-1,1,Vphi,1000,-1,25,Ephi);
          //           obj.FillHistogram(dirname, Form("%s_gated_Ephi-Vphi_vs_Etheta-Vtheta",dirname.c_str()),1100,-50,5,Etheta-Vtheta,700,-3,4,Ephi-Vphi);

          //           //Hit energy difference vs angle calculation difference
          //           obj.FillHistogram(dirname, Form("%s_gated_E1-E2/sum_vs_Etheta-Vtheta",dirname.c_str()),1500,-30,5,Etheta-Vtheta,1000,0,1,(E1-E2)/(E1+E2));

          //           if (ringNum == 1 || ringNum == 4 || ringNum == 12){
          //             obj.FillHistogram(dirname, Form("%s_gated_defaultE_vs_Vtheta_ring%02d",dirname.c_str(),ringNum),200,-1,1,Vtheta, 1600,0,1600,nnEnergy_corrected);
          //             obj.FillHistogram(dirname, Form("%s_gated_swappedE_vs_Vtheta_ring%02d",dirname.c_str(),ringNum),200,-1,1,Vtheta, 1600,0,1600,swappedEnergy);
                      
          //             //Hit vs Angle
          //             obj.FillHistogram(dirname, Form("%s_gated_E2_vs_Vtheta_ring%02d",dirname.c_str(),ringNum),200,-1,1,Vtheta,1200,0,1200,E2);
          //             obj.FillHistogram(dirname, Form("%s_gated_E1_vs_Vtheta_ring%02d",dirname.c_str(),ringNum),200,-1,1,Vtheta,1600,0,1600,E1);

          //             //Hit vs Hit
          //             obj.FillHistogram(dirname, Form("%s_gated_E1_vs_E2_ring%02d",dirname.c_str(),ringNum),1200,0,1200,E2,1600,0,1600,E1);
          //           }
          //         } 
          //         else {
          //           dirname = EnergyGates[eg].first + "_" + LineGates[0].first + "_NOT_" + LineGates[1].first + "_" + InteractionGates[0].first;
          //           //Hit vs Angle
          //           obj.FillHistogram(dirname, Form("%s_gated_E2_vs_Vtheta",dirname.c_str()),200,-1,1,Vtheta,1200,0,1200,E2);
          //           obj.FillHistogram(dirname, Form("%s_gated_E1_vs_Vtheta",dirname.c_str()),200,-1,1,Vtheta,1600,0,1600,E1);

          //           //Hit vs Hit
          //           obj.FillHistogram(dirname, Form("%s_gated_E1_vs_E2",dirname.c_str()),1200,0,1200,E2,1600,0,1600,E1);

          //           //Angle calculation correlation
          //           obj.FillHistogram(dirname, Form("%s_gated_Etheta_vs_Vtheta",dirname.c_str()),200,-1,1,Vtheta,1000,-25,1,Etheta);
          //           obj.FillHistogram(dirname, Form("%s_gated_Ephi_vs_Vphi",dirname.c_str()),200,-1,1,Vphi,1000,-1,25,Ephi);
          //           obj.FillHistogram(dirname, Form("%s_gated_Ephi-Vphi_vs_Etheta-Vtheta",dirname.c_str()),1100,-50,5,Etheta-Vtheta,700,-3,4,Ephi-Vphi);

          //           //Hit energy difference vs angle calculation difference
          //           obj.FillHistogram(dirname, Form("%s_gated_E1-E2/sum_vs_Etheta-Vtheta",dirname.c_str()),1500,-30,5,Etheta-Vtheta,1000,0,1,(E1-E2)/(E1+E2));

          //           if (ringNum == 1 || ringNum == 4 || ringNum == 12){
          //             obj.FillHistogram(dirname, Form("%s_gated_defaultE_vs_Vtheta_ring%02d",dirname.c_str(),ringNum),200,-1,1,Vtheta, 1600,0,1600,nnEnergy_corrected);
          //             obj.FillHistogram(dirname, Form("%s_gated_swappedE_vs_Vtheta_ring%02d",dirname.c_str(),ringNum),200,-1,1,Vtheta, 1600,0,1600,swappedEnergy);

          //             //Hit vs Angle
          //             obj.FillHistogram(dirname, Form("%s_gated_E2_vs_Vtheta_ring%02d",dirname.c_str(),ringNum),200,-1,1,Vtheta,1200,0,1200,E2);
          //             obj.FillHistogram(dirname, Form("%s_gated_E1_vs_Vtheta_ring%02d",dirname.c_str(),ringNum),200,-1,1,Vtheta,1600,0,1600,E1);

          //             //Hit vs Hit
          //             obj.FillHistogram(dirname, Form("%s_gated_E1_vs_E2_ring%02d",dirname.c_str(),ringNum),1200,0,1200,E2,1600,0,1600,E1);
          //           }
          //         }
          //       } 
          //       else {
          //         if (LineGates[1].second){
          //           dirname = EnergyGates[eg].first + "_NOT_" + LineGates[0].first + "_" + LineGates[1].first + "_" + InteractionGates[0].first;
          //           //Hit vs Angle
          //           obj.FillHistogram(dirname, Form("%s_gated_E2_vs_Vtheta",dirname.c_str()),200,-1,1,Vtheta,1200,0,1200,E2);
          //           obj.FillHistogram(dirname, Form("%s_gated_E1_vs_Vtheta",dirname.c_str()),200,-1,1,Vtheta,1600,0,1600,E1);

          //           //Hit vs Hit
          //           obj.FillHistogram(dirname, Form("%s_gated_E1_vs_E2",dirname.c_str()),1200,0,1200,E2,1600,0,1600,E1);

          //           //Angle calculation correlation
          //           obj.FillHistogram(dirname, Form("%s_gated_Etheta_vs_Vtheta",dirname.c_str()),200,-1,1,Vtheta,1000,-25,1,Etheta);
          //           obj.FillHistogram(dirname, Form("%s_gated_Ephi_vs_Vphi",dirname.c_str()),200,-1,1,Vphi,1000,-1,25,Ephi);
          //           obj.FillHistogram(dirname, Form("%s_gated_Ephi-Vphi_vs_Etheta-Vtheta",dirname.c_str()),1100,-50,5,Etheta-Vtheta,700,-3,4,Ephi-Vphi);

          //           //Hit energy difference vs angle calculation difference
          //           obj.FillHistogram(dirname, Form("%s_gated_E1-E2/sum_vs_Etheta-Vtheta",dirname.c_str()),1500,-30,5,Etheta-Vtheta,1000,0,1,(E1-E2)/(E1+E2));

          //           if (ringNum == 1 || ringNum == 4 || ringNum == 12){
          //             obj.FillHistogram(dirname, Form("%s_gated_defaultE_vs_Vtheta_ring%02d",dirname.c_str(),ringNum),200,-1,1,Vtheta, 1600,0,1600,nnEnergy_corrected);
          //             obj.FillHistogram(dirname, Form("%s_gated_swappedE_vs_Vtheta_ring%02d",dirname.c_str(),ringNum),200,-1,1,Vtheta, 1600,0,1600,swappedEnergy);

          //             //Hit vs Angle
          //             obj.FillHistogram(dirname, Form("%s_gated_E2_vs_Vtheta_ring%02d",dirname.c_str(),ringNum),200,-1,1,Vtheta,1200,0,1200,E2);
          //             obj.FillHistogram(dirname, Form("%s_gated_E1_vs_Vtheta_ring%02d",dirname.c_str(),ringNum),200,-1,1,Vtheta,1600,0,1600,E1);

          //             //Hit vs Hit
          //             obj.FillHistogram(dirname, Form("%s_gated_E1_vs_E2_ring%02d",dirname.c_str(),ringNum),1200,0,1200,E2,1600,0,1600,E1);
          //           }
          //         } 
          //         else {
          //           dirname = EnergyGates[eg].first + "_NOT_" + LineGates[0].first + "_NOT_" + LineGates[1].first + "_" + InteractionGates[0].first;
          //           //Hit vs Angle
          //           obj.FillHistogram(dirname, Form("%s_gated_E2_vs_Vtheta",dirname.c_str()),200,-1,1,Vtheta,1200,0,1200,E2);
          //           obj.FillHistogram(dirname, Form("%s_gated_E1_vs_Vtheta",dirname.c_str()),200,-1,1,Vtheta,1600,0,1600,E1);

          //           //Hit vs Hit
          //           obj.FillHistogram(dirname, Form("%s_gated_E1_vs_E2",dirname.c_str()),1200,0,1200,E2,1600,0,1600,E1);

          //           //Angle calculation correlation
          //           obj.FillHistogram(dirname, Form("%s_gated_Etheta_vs_Vtheta",dirname.c_str()),200,-1,1,Vtheta,1000,-25,1,Etheta);
          //           obj.FillHistogram(dirname, Form("%s_gated_Ephi_vs_Vphi",dirname.c_str()),200,-1,1,Vphi,1000,-1,25,Ephi);
          //           obj.FillHistogram(dirname, Form("%s_gated_Ephi-Vphi_vs_Etheta-Vtheta",dirname.c_str()),1100,-50,5,Etheta-Vtheta,700,-3,4,Ephi-Vphi);

          //           //Hit energy difference vs angle calculation difference
          //           obj.FillHistogram(dirname, Form("%s_gated_E1-E2/sum_vs_Etheta-Vtheta",dirname.c_str()),1500,-30,5,Etheta-Vtheta,1000,0,1,(E1-E2)/(E1+E2));

          //           if (ringNum == 1 || ringNum == 4 || ringNum == 12){
          //             obj.FillHistogram(dirname, Form("%s_gated_defaultE_vs_Vtheta_ring%02d",dirname.c_str(),ringNum),200,-1,1,Vtheta, 1600,0,1600,nnEnergy_corrected);
          //             obj.FillHistogram(dirname, Form("%s_gated_swappedE_vs_Vtheta_ring%02d",dirname.c_str(),ringNum),200,-1,1,Vtheta, 1600,0,1600,swappedEnergy);

          //             //Hit vs Angle
          //             obj.FillHistogram(dirname, Form("%s_gated_E2_vs_Vtheta_ring%02d",dirname.c_str(),ringNum),200,-1,1,Vtheta,1200,0,1200,E2);
          //             obj.FillHistogram(dirname, Form("%s_gated_E1_vs_Vtheta_ring%02d",dirname.c_str(),ringNum),200,-1,1,Vtheta,1600,0,1600,E1);

          //             //Hit vs Hit
          //             obj.FillHistogram(dirname, Form("%s_gated_E1_vs_E2_ring%02d",dirname.c_str(),ringNum),1200,0,1200,E2,1600,0,1600,E1);
          //           }
          //         }
          //       }
          //     } 
          //   }
          // }

          // // 2D gates
          // for (int g=0; g < gates_loaded; g++){
          //   std::string tagname(gates_2D[g]->GetTag());
          //   std::string gatename(gates_2D[g]->GetName());
          //   gatename = "NOPEAKCUT_" + gatename;
          //   double varX, varY = -12345;
          //   if (tagname.compare("Ephi-Vphi_Etheta-Vtheta") == 0){
          //     varX = Etheta-Vtheta;
          //     varY = Ephi-Vphi;
          //   }
          //   // if (tagname.compare("E1-E2_Etheta-Vtheta") == 0){
          //   //   varX = Etheta-Vtheta;
          //   //   varY = (E1-E2)/(E1+E2);
          //   // }
          //   if (varX != -12345 && varY != -12345 && gates_2D[g]->IsInside(varX,varY)){
          //     dirname = gatename;

          //     //E angle calculation
          //     obj.FillHistogram(dirname, Form("%s_gated_Etheta",dirname.c_str()),1000,-25,1,Etheta);
          //     //Angle calculation correlation
          //     obj.FillHistogram(dirname, Form("%s_gated_Etheta_vs_Vtheta",dirname.c_str()),200,-1,1,Vtheta,1000,-25,1,Etheta);
          //     obj.FillHistogram(dirname, Form("%s_gated_Ephi_vs_Vphi",dirname.c_str()),200,-1,1,Vphi,1000,-1,25,Ephi);
          //     obj.FillHistogram(dirname, Form("%s_gated_Ephi-Vphi_vs_Etheta-Vtheta",dirname.c_str()),1100,-50,5,Etheta-Vtheta,700,-3,4,Ephi-Vphi);
          //     obj.FillHistogram(dirname, Form("%s_gated_E1-E2/sum_vs_Vtheta",dirname.c_str()),200,-1,1,Vtheta,1000,0,1,(E1-E2)/(E1+E2));
          //     //Hit energy difference vs angle calculation difference
          //     obj.FillHistogram(dirname, Form("%s_gated_E1-E2/sum_vs_Etheta-Vtheta",dirname.c_str()),1500,-30,5,Etheta-Vtheta,1000,0,1,(E1-E2)/(E1+E2));

          //     //Ring by Ring Hists
          //     if (ringNum == 1 || ringNum == 4 || ringNum == 12){
          //       // Doppler corrected vs Angle
          //       obj.FillHistogram(dirname, Form("%s_gated_defaultE_vs_Vtheta_ring%02d_prompt",dirname.c_str(),ringNum),200,-1,1,Vtheta, 1600,0,1600,nnEnergy_corrected);
          //       obj.FillHistogram(dirname, Form("%s_gated_swappedE_vs_Vtheta_ring%02d_prompt",dirname.c_str(),ringNum),200,-1,1,Vtheta, 1600,0,1600,swappedEnergy);
          //       // obj.FillHistogram(dirname, Form("%s_gated_defaultE_vs_Etheta_ring%02d_prompt",dirname.c_str(),ringNum),1000,-25,1,Etheta, 1600,0,1600,nnEnergy_corrected);
          //       // obj.FillHistogram(dirname, Form("%s_gated_swappedE_vs_Etheta_ring%02d_prompt",dirname.c_str(),ringNum),1000,-25,1,Etheta, 1600,0,1600,swappedEnergy);
          //       // E1 vs E2 hit
          //       obj.FillHistogram(dirname, Form("%s_gated_E1_vs_E2_ring%02d",dirname.c_str(),ringNum),1600,0,1600,E2,1600,0,1600,E1);
          //     }

          //     //Interaction Gates
          //     for (int i=0; i < nInteractionGates; i++){
          //       if (InteractionGates[i].second){
          //         dirname = gatename + "_" + InteractionGates[i].first;
          //         //E angle calculation
          //         obj.FillHistogram(dirname, Form("%s_gated_Etheta",dirname.c_str()),1000,-25,1,Etheta);
          //         //Angle calculation correlation
          //         obj.FillHistogram(dirname, Form("%s_gated_Etheta_vs_Vtheta",dirname.c_str()),200,-1,1,Vtheta,1000,-25,1,Etheta);
          //         obj.FillHistogram(dirname, Form("%s_gated_Ephi_vs_Vphi",dirname.c_str()),200,-1,1,Vphi,1000,-1,25,Ephi);
          //         obj.FillHistogram(dirname, Form("%s_gated_Ephi-Vphi_vs_Etheta-Vtheta",dirname.c_str()),1100,-50,5,Etheta-Vtheta,700,-3,4,Ephi-Vphi);
          //         obj.FillHistogram(dirname, Form("%s_gated_E1-E2/sum_vs_Vtheta",dirname.c_str()),200,-1,1,Vtheta,1000,0,1,(E1-E2)/(E1+E2));
          //         //Hit energy difference vs angle calculation difference
          //         obj.FillHistogram(dirname, Form("%s_gated_E1-E2/sum_vs_Etheta-Vtheta",dirname.c_str()),1500,-30,5,Etheta-Vtheta,1000,0,1,(E1-E2)/(E1+E2));

          //         //Ring by Ring Hists
          //         if (ringNum == 1 || ringNum == 4 || ringNum == 12){
          //           // Doppler corrected vs Angle
          //           obj.FillHistogram(dirname, Form("%s_gated_defaultE_vs_Vtheta_ring%02d_prompt",dirname.c_str(),ringNum),200,-1,1,Vtheta, 1600,0,1600,nnEnergy_corrected);
          //           obj.FillHistogram(dirname, Form("%s_gated_swappedE_vs_Vtheta_ring%02d_prompt",dirname.c_str(),ringNum),200,-1,1,Vtheta, 1600,0,1600,swappedEnergy);
          //           // obj.FillHistogram(dirname, Form("%s_gated_defaultE_vs_Etheta_ring%02d_prompt",dirname.c_str(),ringNum),1000,-25,1,Etheta, 1600,0,1600,nnEnergy_corrected);
          //           // obj.FillHistogram(dirname, Form("%s_gated_swappedE_vs_Etheta_ring%02d_prompt",dirname.c_str(),ringNum),1000,-25,1,Etheta, 1600,0,1600,swappedEnergy);
          //           // E1 vs E2 hit
          //           obj.FillHistogram(dirname, Form("%s_gated_E1_vs_E2_ring%02d",dirname.c_str(),ringNum),1600,0,1600,E2,1600,0,1600,E1);
          //         }
          //       }
          //     }
          //   }
          // }

          // //Hit vs Angle
          // obj.FillHistogram(dirname, Form("%s_gated_E2_vs_Vtheta",dirname.c_str()),200,-1,1,Vtheta,1200,0,1200,E2);
          // obj.FillHistogram(dirname, Form("%s_gated_E1_vs_Vtheta",dirname.c_str()),200,-1,1,Vtheta,1600,0,1600,E1);

          // //Angle calculation correlation
          // obj.FillHistogram(dirname, Form("%s_gated_Etheta_vs_Vtheta",dirname.c_str()),200,-1,1,Vtheta,1000,-25,1,Etheta);
          // obj.FillHistogram(dirname, Form("%s_gated_Ephi_vs_Vphi",dirname.c_str()),200,-1,1,Vphi,1000,-1,25,Ephi);
          // obj.FillHistogram(dirname, Form("%s_gated_Ephi-Vphi_vs_Etheta-Vtheta",dirname.c_str()),1100,-50,5,Etheta-Vtheta,700,-3,4,Ephi-Vphi);

          // //Hit energy difference vs angle calculation difference
          // obj.FillHistogram(dirname, Form("%s_gated_E1-E2/sum_vs_Etheta-Vtheta",dirname.c_str()),1500,-30,5,Etheta-Vtheta,1000,0,1,(E1-E2)/(E1+E2));

          if (ringNum == 1 || ringNum == 4 || ringNum == 12){
            // Doppler corrected vs Angle
            obj.FillHistogram(dirname, Form("%s_gated_defaultE_vs_Vtheta_ring%02d_prompt",dirname.c_str(),ringNum),200,-1,1,Vtheta, 1600,0,1600,nnEnergy_corrected);
            obj.FillHistogram(dirname, Form("%s_gated_swappedE_vs_Vtheta_ring%02d_prompt",dirname.c_str(),ringNum),200,-1,1,Vtheta, 1600,0,1600,swappedEnergy);
            // obj.FillHistogram(dirname, Form("%s_gated_defaultE_vs_Etheta_ring%02d_prompt",dirname.c_str(),ringNum),1000,-25,1,Etheta, 1600,0,1600,nnEnergy_corrected);
            // obj.FillHistogram(dirname, Form("%s_gated_swappedE_vs_Etheta_ring%02d_prompt",dirname.c_str(),ringNum),1000,-25,1,Etheta, 1600,0,1600,swappedEnergy);
            
            // Hit vs Angle
            obj.FillHistogram(dirname, Form("%s_gated_E2_vs_Vtheta_ring%02d",dirname.c_str(),ringNum),200,-1,1,Vtheta,1200,0,1200,E2);
            obj.FillHistogram(dirname, Form("%s_gated_E1_vs_Vtheta_ring%02d",dirname.c_str(),ringNum),200,-1,1,Vtheta,1600,0,1600,E1);
          }

          //gate on angles above E/2
          double angleThresh = TMath::ACos(1 - 511/((E1+E2)/2) + 511/(E1+E2));
          if (TMath::Cos(Vtheta) > angleThresh) {
            dirname = "FEP_gated_>thresh_angle_gated";

            // //Hit vs Angle
            // obj.FillHistogram(dirname, Form("%s_gated_E2_vs_Vtheta",dirname.c_str()),200,-1,1,Vtheta,1200,0,1200,E2);
            // obj.FillHistogram(dirname, Form("%s_gated_E1_vs_Vtheta",dirname.c_str()),200,-1,1,Vtheta,1600,0,1600,E1);

            // //Angle calculation correlation
            // obj.FillHistogram(dirname, Form("%s_gated_Etheta_vs_Vtheta",dirname.c_str()),200,-1,1,Vtheta,1000,-25,1,Etheta);
            // obj.FillHistogram(dirname, Form("%s_gated_Ephi_vs_Vphi",dirname.c_str()),200,-1,1,Vphi,1000,-1,25,Ephi);
            // obj.FillHistogram(dirname, Form("%s_gated_Ephi-Vphi_vs_Etheta-Vtheta",dirname.c_str()),1100,-50,5,Etheta-Vtheta,700,-3,4,Ephi-Vphi);

            // //Hit energy difference vs angle calculation difference
            // obj.FillHistogram(dirname, Form("%s_gated_E1-E2/sum_vs_Etheta-Vtheta",dirname.c_str()),1500,-30,5,Etheta-Vtheta,1000,0,1,(E1-E2)/(E1+E2));

            if (ringNum == 1 || ringNum == 4 || ringNum == 12){
              // Doppler corrected vs Angle
              obj.FillHistogram(dirname, Form("%s_gated_defaultE_vs_Vtheta_ring%02d_prompt",dirname.c_str(),ringNum),200,-1,1,Vtheta, 1600,0,1600,nnEnergy_corrected);
              obj.FillHistogram(dirname, Form("%s_gated_swappedE_vs_Vtheta_ring%02d_prompt",dirname.c_str(),ringNum),200,-1,1,Vtheta, 1600,0,1600,swappedEnergy);
              // obj.FillHistogram(dirname, Form("%s_gated_defaultE_vs_Etheta_ring%02d_prompt",dirname.c_str(),ringNum),1000,-25,1,Etheta, 1600,0,1600,nnEnergy_corrected);
              // obj.FillHistogram(dirname, Form("%s_gated_swappedE_vs_Etheta_ring%02d_prompt",dirname.c_str(),ringNum),1000,-25,1,Etheta, 1600,0,1600,swappedEnergy);
              
              // Hit vs Angle
              obj.FillHistogram(dirname, Form("%s_gated_E2_vs_Vtheta_ring%02d",dirname.c_str(),ringNum),200,-1,1,Vtheta,1200,0,1200,E2);
              obj.FillHistogram(dirname, Form("%s_gated_E1_vs_Vtheta_ring%02d",dirname.c_str(),ringNum),200,-1,1,Vtheta,1600,0,1600,E1);
            }
          } else {
            dirname = "FEP_gated_<thresh_angle_gated";

            // //Hit vs Angle
            // obj.FillHistogram(dirname, Form("%s_gated_E2_vs_Vtheta",dirname.c_str()),200,-1,1,Vtheta,1200,0,1200,E2);
            // obj.FillHistogram(dirname, Form("%s_gated_E1_vs_Vtheta",dirname.c_str()),200,-1,1,Vtheta,1600,0,1600,E1);

            // //Angle calculation correlation
            // obj.FillHistogram(dirname, Form("%s_gated_Etheta_vs_Vtheta",dirname.c_str()),200,-1,1,Vtheta,1000,-25,1,Etheta);
            // obj.FillHistogram(dirname, Form("%s_gated_Ephi_vs_Vphi",dirname.c_str()),200,-1,1,Vphi,1000,-1,25,Ephi);
            // obj.FillHistogram(dirname, Form("%s_gated_Ephi-Vphi_vs_Etheta-Vtheta",dirname.c_str()),1100,-50,5,Etheta-Vtheta,700,-3,4,Ephi-Vphi);

            // //Hit energy difference vs angle calculation difference
            // obj.FillHistogram(dirname, Form("%s_gated_E1-E2/sum_vs_Etheta-Vtheta",dirname.c_str()),1500,-30,5,Etheta-Vtheta,1000,0,1,(E1-E2)/(E1+E2));

            if (ringNum == 1 || ringNum == 4 || ringNum == 12){
              // Doppler corrected vs Angle
              obj.FillHistogram(dirname, Form("%s_gated_defaultE_vs_Vtheta_ring%02d_prompt",dirname.c_str(),ringNum),200,-1,1,Vtheta, 1600,0,1600,nnEnergy_corrected);
              obj.FillHistogram(dirname, Form("%s_gated_swappedE_vs_Vtheta_ring%02d_prompt",dirname.c_str(),ringNum),200,-1,1,Vtheta, 1600,0,1600,swappedEnergy);
              // obj.FillHistogram(dirname, Form("%s_gated_defaultE_vs_Etheta_ring%02d_prompt",dirname.c_str(),ringNum),1000,-25,1,Etheta, 1600,0,1600,nnEnergy_corrected);
              // obj.FillHistogram(dirname, Form("%s_gated_swappedE_vs_Etheta_ring%02d_prompt",dirname.c_str(),ringNum),1000,-25,1,Etheta, 1600,0,1600,swappedEnergy);
              
              // Hit vs Angle
              obj.FillHistogram(dirname, Form("%s_gated_E2_vs_Vtheta_ring%02d",dirname.c_str(),ringNum),200,-1,1,Vtheta,1200,0,1200,E2);
              obj.FillHistogram(dirname, Form("%s_gated_E1_vs_Vtheta_ring%02d",dirname.c_str(),ringNum),200,-1,1,Vtheta,1600,0,1600,E1);
            }
          }

          // //gate on interactions
          // for (int i=0; i < nInteractionGates; i++){
          //   if (InteractionGates[i].second){
          //     dirname = InteractionGates[i].first;
              
          //     //Hit vs Angle
          //     obj.FillHistogram(dirname, Form("%s_gated_E2_vs_Vtheta",dirname.c_str()),200,-1,1,Vtheta,1200,0,1200,E2);
          //     obj.FillHistogram(dirname, Form("%s_gated_E1_vs_Vtheta",dirname.c_str()),200,-1,1,Vtheta,1600,0,1600,E1);
              
          //     //Angle calculation correlation
          //     obj.FillHistogram(dirname, Form("%s_gated_Etheta_vs_Vtheta",dirname.c_str()),200,-1,1,Vtheta,1000,-25,1,Etheta);
          //     obj.FillHistogram(dirname, Form("%s_gated_Ephi_vs_Vphi",dirname.c_str()),200,-1,1,Vphi,1000,-1,25,Ephi);
          //     obj.FillHistogram(dirname, Form("%s_gated_Ephi-Vphi_vs_Etheta-Vtheta",dirname.c_str()),1100,-50,5,Etheta-Vtheta,700,-3,4,Ephi-Vphi);

          //     //Hit energy difference vs angle calculation difference
          //     obj.FillHistogram(dirname, Form("%s_gated_E1-E2/sum_vs_Etheta-Vtheta",dirname.c_str()),1500,-30,5,Etheta-Vtheta,1000,0,1,(E1-E2)/(E1+E2));

          //     //Ring By Ring Hits
          //     if (ringNum == 1 || ringNum == 4 || ringNum == 12){
          //       // Doppler corrected vs Angle
          //       obj.FillHistogram(dirname, Form("%s_gated_defaultE_vs_Vtheta_ring%02d_prompt",dirname.c_str(),ringNum),200,-1,1,Vtheta, 1600,0,1600,nnEnergy_corrected);
          //       obj.FillHistogram(dirname, Form("%s_gated_swappedE_vs_Vtheta_ring%02d_prompt",dirname.c_str(),ringNum),200,-1,1,Vtheta, 1600,0,1600,swappedEnergy);
          //       // obj.FillHistogram(dirname, Form("%s_gated_defaultE_vs_Etheta_ring%02d_prompt",dirname.c_str(),ringNum),1000,-25,1,Etheta, 1600,0,1600,nnEnergy_corrected);
          //       // obj.FillHistogram(dirname, Form("%s_gated_swappedE_vs_Etheta_ring%02d_prompt",dirname.c_str(),ringNum),1000,-25,1,Etheta, 1600,0,1600,swappedEnergy);

          //       // Hit vs Angle
          //       obj.FillHistogram(dirname, Form("%s_gated_E2_vs_Vtheta_ring%02d",dirname.c_str(),ringNum),200,-1,1,Vtheta,1200,0,1200,E2);
          //       obj.FillHistogram(dirname, Form("%s_gated_E1_vs_Vtheta_ring%02d",dirname.c_str(),ringNum),200,-1,1,Vtheta,1600,0,1600,E1);
          //       // E1 vs E2 hit
          //       obj.FillHistogram(dirname, Form("%s_gated_E1_vs_E2_ring%02d",dirname.c_str(),ringNum),1600,0,1600,E2,1600,0,1600,E1);
          //     }
          //   }
          // }
        }
      }

      if (n==0){
        if (gretsim->GetGretinaSimHit(0).IsFEP()){
          dirname = "basicsim";
          obj.FillHistogram(dirname,"gamma_corrected_n0_prompt_fep", 1600,0,1600, nnEnergy_corrected);
          obj.FillHistogram("crystal-specific", Form("gamma_corrected_n0_ring%02d_crystal%d_prompt_fep",ringNum,cryID),1600,0,1600, nnEnergy_corrected);
        }
      }
      
    }
  }

  if(numobj!=list->GetSize()){
    list->Sort();
  }
}
