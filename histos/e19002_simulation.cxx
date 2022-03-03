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

double calcScatteringAngle(const TGretinaHit &hit1, const TGretinaHit &hit2){
  TVector3 pos1 = hit1.GetPosition();
  TVector3 diff = hit2.GetPosition() - pos1;
  return 180.0/TMath::Pi() * diff.Angle(pos1);
}

double calcComptonAngle(double E1, double E2){
  double argument = 1 - 511/E2 + 511/(E1 + E2);
  // if (fabs(argument) > 1) argument = 1 - 511/E1 + 511/(E1 + E2);
  return TMath::ACos(argument)*180/TMath::Pi();
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
    
    //SINGLES
    int gSize = gretina->Size();
    for (int i=0; i < gSize; i++){
      TGretinaHit &hit = gretina->GetGretinaHit(i);

      double energy_corrected = hit.GetDopplerYta(s800sim->AdjustedBeta(GValue::Value("BETA")), yta, &track);
      double energy = hit.GetDoppler(GValue::Value("BETA"));
      double core_energy = hit.GetCoreEnergy();
      double theta = hit.GetTheta();
      int cryID = hit.GetCrystalId();
      
      obj.FillHistogram(dirname, "core_energy_prompt", 8192,0,8192, core_energy);
      obj.FillHistogram(dirname, "gamma_singles_prompt", 8192,0,8192, energy);
      obj.FillHistogram(dirname, "gamma_corrected_singles_prompt", 8192,0,8192, energy_corrected);
      obj.FillHistogram(dirname, "gamma_corrected_vs_theta_prompt", 8192,0,8192, energy_corrected, 100, 0, 2.5, theta);
      obj.FillHistogram(dirname, "gamma_corrected_vs_crystalID_prompt", 56, 24, 80, cryID, 8192,0,8192, energy_corrected);
      obj.FillHistogram(dirname, "core_energy_vs_theta_prompt", 8192,0,8192, hit.GetCoreEnergy(), 100, 0, 2.5, theta);
      
      //crystal angles
      obj.FillHistogram(dirname, Form("theta_vs_crystalID_ring%d",gretina->GetRingNumber(cryID)),56,24,80,cryID,180,0,180,hit.GetThetaDeg());
      
    }

    //NNADDBACK
    //loop over multiplicity
    for (int n=0; n<2; n++){
      //loop over hits for each multiplicity spectrum
      int nnSize = gretina->NNAddbackSize(n);
      for (int i=0; i < nnSize; i++){
      dirname = "basicsim";

        //get hit and hit data 
        TGretinaHit nnhit = gretina->GetNNAddbackHit(n,i);
        int cryID = nnhit.GetCrystalId();
        int ringNum = nnhit.GetRingNumber();
        double nnEnergy_corrected = nnhit.GetDopplerYta(s800sim->AdjustedBeta(GValue::Value("BETA")), yta, &track);
        
        //exclude the ng spectrum (n==3)
        if (n < 3){
          obj.FillHistogram(dirname, "gamma_corrected_addback_prompt", 8192,0,8192, nnEnergy_corrected);
        }

        if (n == 1) {
          // // GAMMA GAMMA CORRELATION
          // for (int j=0; j < nnSize; j++){
          //   if (i==j) continue;
          //   TGretinaHit nnhit2 = gretina->GetNNAddbackHit(n,j);
          //   double nnEnergy_corrected2 = nnhit2.GetDopplerYta(s800sim->AdjustedBeta(GValue::Value("BETA")), yta, &track);
          //   obj.FillHistogram(dirname, "gamma_gamma", 8192,0,8192, nnEnergy_corrected2, 8192,0,8192, nnEnergy_corrected);
          // }

          //POLARIZATION
          // if ( PairHit(nnhit,redPairs) ){
          //   obj.FillHistogram(dirname,"gamma_corrected_addback_prompt_red_pair", 8192,0,8192, nnEnergy_corrected);
          // }

          // if ( PairHit(nnhit,goldPairs) ){
          //   obj.FillHistogram(dirname,"gamma_corrected_addback_prompt_gold_pair", 8192,0,8192, nnEnergy_corrected);
          // }

          // if ( PairHit(nnhit,bluePairs) ){
          //   obj.FillHistogram(dirname,"gamma_corrected_addback_prompt_blue_pair", 8192,0,8192, nnEnergy_corrected);
          // }

          //N1 Neighbor correlations
          TGretinaHit nnhit1 = nnhit.GetInitialHit();
          TGretinaHit nnhit2 = nnhit.GetNeighbor();
          int hit1Ints = nnhit1.NumberOfInteractions();
          int hit2Ints = nnhit2.NumberOfInteractions();

          //Number of interactions
          obj.FillHistogram("interactions","nnhit1_interactions",12,0,12,hit1Ints);
          obj.FillHistogram("interactions","nnhit2_interactions",12,0,12,hit2Ints);
          
          //compton analysis
          double hiHit = nnhit1.GetCoreEnergy();
          double loHit = nnhit2.GetCoreEnergy();
          double posAngle = calcScatteringAngle(nnhit1,nnhit2); 
          double compPosAngle = calcScatteringAngle(nnhit2,nnhit1);
          double engAngle = calcComptonAngle(hiHit,loHit);
          double compEngAngle = calcComptonAngle(loHit,hiHit);
          
          // if (fabs(posAngle - engAngle) > 10 && hit1Ints==1 && hit2Ints==1){
          //   std::swap(nnhit1,nnhit2);
          // }
          // nnhit1.NNAdd(nnhit2);
          // double newEnergy = nnhit1.GetDopplerYta(s800sim->AdjustedBeta(GValue::Value("BETA")), yta, &track);
          // obj.FillHistogram("crystal-specific", Form("gamma_corrected_newEnergy_n%d_ring%02d",n,ringNum),8192,0,8192, newEnergy);
          // obj.FillHistogram("scattering-angles", Form("gamma_corrected_newEnergy_sorted_n%d_vs_scangle_ring%02d",n,ringNum),180,0,180,posAngle, 2000,0,2000,newEnergy);

          obj.FillHistogram("scattering-angles", Form("pos_angles_n%d_ring%02d_prompt",n,ringNum),180,0,180,posAngle);
          obj.FillHistogram("scattering-angles", Form("eng_angles_n%d_ring%02d_prompt",n,ringNum),180,0,180,engAngle);
          // obj.FillHistogram("scattering-angles", Form("high_energy_hit_vs_scangle_n%d_ring%02d_crystal%d_prompt",n,ringNum,cryID),180,0,180,posAngle,2000,0,2000,hiHit);
          // obj.FillHistogram("scattering-angles", Form("low_energy_hit_vs_scangle_n%d_ring%02d_crystal%d_prompt",n,ringNum,cryID),180,0,180,posAngle,2000,0,2000,loHit);
          // obj.FillHistogram("scattering-angles", Form("tot_energy_hit_vs_scangle_n%d_ring%02d_crystal%d_prompt",n,ringNum,cryID),180,0,180,posAngle,2000,0,2000,loHit+hiHit);
          
          //make swapped energy
          nnhit1 = nnhit.GetInitialHit();
          nnhit2 = nnhit.GetNeighbor();
          nnhit2.NNAdd(nnhit1);
          double swappedEnergy = nnhit2.GetDopplerYta(s800sim->AdjustedBeta(GValue::Value("BETA")), yta, &track);

          // obj.FillHistogram(dirname, Form("gamma_corrected_n%d_prompt",n), 8192,0,8192, nnEnergy_corrected);
          // obj.FillHistogram(dirname, Form("gamma_corrected_n%d_swapped_prompt",n), 8192,0,8192, swappedEnergy);
          obj.FillHistogram("crystal-specific", Form("gamma_corrected_n%d_ring%02d_prompt",n,ringNum),1600,0,1600, nnEnergy_corrected);
          obj.FillHistogram("crystal-specific", Form("gamma_corrected_n%d_swapped_ring%02d_prompt",n,ringNum),1600,0,1600, swappedEnergy);
          obj.FillHistogram("scattering-angles", Form("gamma_corrected_n%d_vs_posAngle_ring%02d_prompt",n,ringNum),180,0,180,posAngle, 1600,0,1600,nnEnergy_corrected);
          obj.FillHistogram("scattering-angles", Form("gamma_corrected_n%d_swapped_vs_posAngle_ring%02d_prompt",n,ringNum),180,0,180,posAngle, 1600,0,1600,swappedEnergy);
          // obj.FillHistogram("scattering-angles", Form("gamma_corrected_n%d_vs_engAngle_ring%02d_prompt",n,ringNum),180,0,180,engAngle, 1600,0,1600,nnEnergy_corrected);
          // obj.FillHistogram("scattering-angles", Form("gamma_corrected_n%d_swapped_vs_engAngle_ring%02d_prompt",n,ringNum),180,0,180,engAngle, 1600,0,1600,swappedEnergy);

          //GATED STUFF
          std::vector<std::pair<std::string,bool>> EnergyGates;
          EnergyGates.push_back(std::make_pair("defaultE",fabs(nnEnergy_corrected - FEP) < 1.5 * TMath::Sqrt(FEP)));
          EnergyGates.push_back(std::make_pair("swappedE",fabs(swappedEnergy - FEP) < 1.5 * TMath::Sqrt(FEP)));
          int nEnergyGates = (int) EnergyGates.size();

          std::vector<std::pair<std::string,bool>> InteractionGates;
          // InteractionGates.push_back(std::make_pair("hiHit_int>1",hit1Ints>1));
          // InteractionGates.push_back(std::make_pair("hiHit_int=1",hit1Ints==1));
          // InteractionGates.push_back(std::make_pair("loHit_int>1",hit2Ints>1));
          // InteractionGates.push_back(std::make_pair("loHit_int=1",hit2Ints==1));
          InteractionGates.push_back(std::make_pair("bothHit_int>1",hit1Ints>1 && hit2Ints>1));
          InteractionGates.push_back(std::make_pair("bothHit_int=1",hit1Ints==1 && hit2Ints==1));
          int nInteractionGates = (int) InteractionGates.size();

          //2D gates
          // for (int g=0; g < gates_loaded; g++){
          //   std::string gatename(gates_2D[g]->GetName());
          //   double varX, varY = -12345;
          //   if (gatename.compare("swappedE_posAngle") == 0){
          //     varX = posAngle;
          //     varY = swappedEnergy;
          //   }
          //   if (gatename.compare("engAngle_posAngle") == 0){
          //     varX = posAngle;
          //     varY = engAngle;
          //   }
          //   if (gates_2D[g]->IsInside(varX,varY) && varX != -12345 && varY != -12345){
          //     obj.FillHistogram(gatename, Form("%s_gated_engAngle_vs_posAngle_ring%02d",gatename.c_str(),ringNum),180,0,180,posAngle,180,0,180,engAngle);
          //     obj.FillHistogram(gatename, Form("%s_gated_loHit_vs_posAngle_ring%02d",gatename.c_str(),ringNum),180,0,180,posAngle,1200,0,1200,loHit);
          //     obj.FillHistogram(gatename, Form("%s_gated_hiHit_vs_posAngle_ring%02d",gatename.c_str(),ringNum),180,0,180,posAngle,1600,0,1600,hiHit);
          //     obj.FillHistogram(gatename, Form("%s_gated_engAngle_vs_posAngle",gatename.c_str()),180,0,180,posAngle,180,0,180,engAngle);
          //     obj.FillHistogram(gatename, Form("%s_gated_loHit_vs_posAngle",gatename.c_str()),180,0,180,posAngle,1200,0,1200,loHit);
          //     obj.FillHistogram(gatename, Form("%s_gated_hiHit_vs_posAngle",gatename.c_str()),180,0,180,posAngle,1600,0,1600,hiHit);
          //     obj.FillHistogram(gatename, Form("%s_gated_hiHit",gatename.c_str()),1600,0,1600,hiHit);
          //     obj.FillHistogram(gatename, Form("%s_gated_loHit",gatename.c_str()),1600,0,1600,loHit);
          //   }
          // }
          
          //gate on energy
          for (int eg=0; eg < nEnergyGates; eg++){
            if (EnergyGates[eg].second){
              dirname = EnergyGates[eg].first;
              // obj.FillHistogram(dirname, Form("%s_gated_engAngle_vs_posAngle_ring%02d",dirname.c_str(),ringNum),180,0,180,posAngle,180,0,180,engAngle);
              // obj.FillHistogram(dirname, Form("%s_gated_loHit_vs_posAngle_ring%02d",dirname.c_str(),ringNum),180,0,180,posAngle,1200,0,1200,loHit);
              // obj.FillHistogram(dirname, Form("%s_gated_hiHit_vs_posAngle_ring%02d",dirname.c_str(),ringNum),180,0,180,posAngle,1600,0,1600,hiHit);
              // obj.FillHistogram(dirname, Form("%s_gated_loHit_vs_posAngle",dirname.c_str()),180,0,180,posAngle,1200,0,1200,loHit);
              // obj.FillHistogram(dirname, Form("%s_gated_hiHit_vs_posAngle",dirname.c_str()),180,0,180,posAngle,1600,0,1600,hiHit);
              obj.FillHistogram(dirname, Form("%s_gated_engAngle_vs_posAngle",dirname.c_str()),180,0,180,posAngle,180,0,180,engAngle);
              if (ringNum == 1 || ringNum == 4 || ringNum == 8 || ringNum == 12){
                obj.FillHistogram(dirname, Form("%s_gated_defaultE_vs_posAngle_ring%02d_prompt",dirname.c_str(),ringNum),180,0,180,posAngle, 1600,0,1600,nnEnergy_corrected);
                obj.FillHistogram(dirname, Form("%s_gated_swappedE_vs_posAngle_ring%02d_prompt",dirname.c_str(),ringNum),180,0,180,posAngle, 1600,0,1600,swappedEnergy);
              }

              //Hit energy difference vs angle calculation difference
              obj.FillHistogram(dirname, Form("%s_gated_hiHit-loHit_vs_engAngle-posAngle",dirname.c_str()),360,-180,180,engAngle-posAngle,1200,0,1200,hiHit-loHit);
              obj.FillHistogram(dirname, Form("%s_gated_hiHit-loHit/sum_vs_engAngle-posAngle",dirname.c_str()),360,-180,180,engAngle-posAngle,1000,0,1,(hiHit-loHit)/(hiHit+loHit));
              obj.FillHistogram(dirname, Form("%s_gated_hiHit-loHit/sum_vs_(engAngle-posAngle)/sum",dirname.c_str()),360,-1,1,(engAngle-posAngle)/(posAngle+engAngle),1000,0,1,(hiHit-loHit)/(hiHit+loHit));

              //gate on interactions
              for (int i=0; i < nInteractionGates; i++){
                if (InteractionGates[i].second){
                  dirname = EnergyGates[eg].first + "_" + InteractionGates[i].first;
                  // obj.FillHistogram(dirname, Form("%s_gated_engAngle_vs_posAngle_ring%02d",dirname.c_str(),ringNum),180,0,180,posAngle,180,0,180,engAngle);
                  
                  //Hit vs Angle
                  obj.FillHistogram(dirname, Form("%s_gated_loHit_vs_posAngle",dirname.c_str()),180,0,180,posAngle,1200,0,1200,loHit);
                  obj.FillHistogram(dirname, Form("%s_gated_hiHit_vs_posAngle",dirname.c_str()),180,0,180,posAngle,1600,0,1600,hiHit);
                  if (ringNum == 1 || ringNum == 4 || ringNum == 8 || ringNum == 12){
                    obj.FillHistogram(dirname, Form("%s_gated_loHit_vs_posAngle_ring%02d",dirname.c_str(),ringNum),180,0,180,posAngle,1200,0,1200,loHit);
                    obj.FillHistogram(dirname, Form("%s_gated_hiHit_vs_posAngle_ring%02d",dirname.c_str(),ringNum),180,0,180,posAngle,1600,0,1600,hiHit);
                  }
                  
                  //Angle calculation correlation
                  obj.FillHistogram(dirname, Form("%s_gated_engAngle_vs_posAngle",dirname.c_str()),180,0,180,posAngle,180,0,180,engAngle);
                  obj.FillHistogram(dirname, Form("%s_gated_abs(engAngle-posAngle)",dirname.c_str()),180,0,180,posAngle -engAngle);

                  //Hit energy difference vs angle calculation difference
                  obj.FillHistogram(dirname, Form("%s_gated_hiHit-loHit_vs_engAngle-posAngle",dirname.c_str()),360,-180,180,engAngle-posAngle,1200,0,1200,hiHit-loHit);
                  obj.FillHistogram(dirname, Form("%s_gated_hiHit-loHit/sum_vs_engAngle-posAngle",dirname.c_str()),360,-180,180,engAngle-posAngle,1000,0,1,(hiHit-loHit)/(hiHit+loHit));
                  obj.FillHistogram(dirname, Form("%s_gated_hiHit-loHit/sum_vs_(engAngle-posAngle)/sum",dirname.c_str()),360,-1,1,(engAngle-posAngle)/(posAngle+engAngle),1000,0,1,(hiHit-loHit)/(hiHit+loHit));

                  //Doppler corrected vs Angle
                  if (ringNum == 1 || ringNum == 4 || ringNum == 8 || ringNum == 12){
                    obj.FillHistogram(dirname, Form("%s_gated_defaultE_vs_posAngle_ring%02d_prompt",dirname.c_str(),ringNum),180,0,180,posAngle, 1600,0,1600,nnEnergy_corrected);
                    obj.FillHistogram(dirname, Form("%s_gated_swappedE_vs_posAngle_ring%02d_prompt",dirname.c_str(),ringNum),180,0,180,posAngle, 1600,0,1600,swappedEnergy);
                  }
                  
                  //2D gates
                  for (int g=0; g < gates_loaded; g++){
                    std::string tagname(gates_2D[g]->GetTag());
                    std::string gatename(gates_2D[g]->GetName());
                    double varX, varY = -12345;
                    if (tagname.compare("engAngle_posAngle") == 0){
                      varX = posAngle;
                      varY = engAngle;
                    }
                    if (varX != -12345 && varY != -12345 && gates_2D[g]->IsInside(varX,varY)){
                      dirname = EnergyGates[eg].first + "_" + InteractionGates[i].first + "_" + gatename;

                      //Hit vs Angle
                      obj.FillHistogram(dirname, Form("%s_gated_loHit_vs_posAngle",dirname.c_str()),180,0,180,posAngle,1200,0,1200,loHit);
                      obj.FillHistogram(dirname, Form("%s_gated_hiHit_vs_posAngle",dirname.c_str()),180,0,180,posAngle,1600,0,1600,hiHit);
                      if (ringNum == 1 || ringNum == 4 || ringNum == 8 || ringNum == 12){
                        obj.FillHistogram(dirname, Form("%s_gated_loHit_vs_posAngle_ring%02d",dirname.c_str(),ringNum),180,0,180,posAngle,1200,0,1200,loHit);
                        obj.FillHistogram(dirname, Form("%s_gated_hiHit_vs_posAngle_ring%02d",dirname.c_str(),ringNum),180,0,180,posAngle,1600,0,1600,hiHit);
                      }

                      //Angle Calculation correlations
                      obj.FillHistogram(dirname, Form("%s_gated_engAngle_vs_complPosAngle",dirname.c_str()),180,0,180,compPosAngle,180,0,180,engAngle);
                      obj.FillHistogram(dirname, Form("%s_gated_complEngAngle_vs_complPosAngle",dirname.c_str()),180,0,180,compPosAngle,180,0,180,compEngAngle);
                      obj.FillHistogram(dirname, Form("%s_gated_abs(engAngle-posAngle)",dirname.c_str()),180,0,180,fabs(posAngle -engAngle));
                      
                      //Hit energy difference vs angle calculation difference
                      obj.FillHistogram(dirname, Form("%s_gated_hiHit-loHit_vs_engAngle-posAngle",dirname.c_str()),360,-180,180,engAngle-posAngle,1200,0,1200,hiHit-loHit);
                      obj.FillHistogram(dirname, Form("%s_gated_hiHit-loHit/sum_vs_engAngle-posAngle",dirname.c_str()),360,-180,180,engAngle-posAngle,1000,0,1,(hiHit-loHit)/(hiHit+loHit));
                      obj.FillHistogram(dirname, Form("%s_gated_hiHit-loHit/sum_vs_(engAngle-posAngle)/sum",dirname.c_str()),360,-1,1,(engAngle-posAngle)/(posAngle+engAngle),1000,0,1,(hiHit-loHit)/(hiHit+loHit));

                      // Doppler corrected vs Angle
                      if (ringNum == 1 || ringNum == 4 || ringNum == 8 || ringNum == 12){
                        obj.FillHistogram(dirname, Form("%s_gated_defaultE_vs_posAngle_ring%02d_prompt",dirname.c_str(),ringNum),180,0,180,posAngle, 1600,0,1600,nnEnergy_corrected);
                        obj.FillHistogram(dirname, Form("%s_gated_swappedE_vs_posAngle_ring%02d_prompt",dirname.c_str(),ringNum),180,0,180,posAngle, 1600,0,1600,swappedEnergy);
                      }
                    }
                  }
                }
              }
            }
          }
        }

        // char *multiplicity = Form("%d",n);
        // if (n == 3) multiplicity = Form("g");
        // obj.FillHistogram(dirname, Form("gamma_corrected_n%s_prompt",multiplicity), 8192,0,8192, nnEnergy_corrected);
        // obj.FillHistogram(dirname, Form("gamma_corrected_n%s_ring%02d_crystal%d_prompt",multiplicity,ringNum,cryID),8192,0,8192, nnEnergy_corrected);
        if (n==0){
          if (gretsim->GetGretinaSimHit(0).IsFEP()){
            dirname = "basicsim";
            obj.FillHistogram(dirname,"gamma_corrected_n0_prompt_fep", 1600,0,1600, nnEnergy_corrected);
            obj.FillHistogram("crystal-specific", Form("gamma_corrected_n0_ring%02d_crystal%d_prompt_fep",ringNum,cryID),1600,0,1600, nnEnergy_corrected);
          }
        }
        
      }
    }
  }

  if(numobj!=list->GetSize()){
    list->Sort();
  }
}
