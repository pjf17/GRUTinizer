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

std::vector<std::pair<int,int>> OneQuadPairs = {
  std::make_pair(44,45),
  std::make_pair(46,47),
  std::make_pair(45,46),
  std::make_pair(44,47),
  std::make_pair(44,46),
  std::make_pair(48,51),
  std::make_pair(49,50),
  std::make_pair(48,49),
  std::make_pair(50,51),
  std::make_pair(48,50),
  std::make_pair(56,59),
  std::make_pair(57,58),
  std::make_pair(56,57),
  std::make_pair(58,59),
  std::make_pair(56,58),
  std::make_pair(60,61),
  std::make_pair(62,63),
  std::make_pair(61,62),
  std::make_pair(60,63),
  std::make_pair(60,62),
  std::make_pair(64,67),
  std::make_pair(65,66),
  std::make_pair(64,65),
  std::make_pair(66,67),
  std::make_pair(64,66),
  std::make_pair(68,69),
  std::make_pair(70,71),
  std::make_pair(69,70),
  std::make_pair(68,71),
  std::make_pair(68,70),
  std::make_pair(76,78),
  std::make_pair(76,79),
  std::make_pair(78,79)
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
  index++;
  for (auto &p : bluePairs){
    int label1 = p.first*100 + p.second;
    int label2 = p.first + p.second*100;
    SCATTERPAIRS.insert(std::pair<int,int>(label1,index));
    SCATTERPAIRS.insert(std::pair<int,int>(label2,index));
    index++;
  }
  index++;
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
  double yta=0;
  if (!stopped){
    track = s800sim->Track(0,0);
    yta = s800sim->GetS800SimHit(0).GetYTA();
  }
  
  makemap();
  if (gretina){
    TGretSimHit simHit = gretsim->GetGretinaSimHit(0);
    double gammaEn = simHit.GetEn();
    bool isFEP = simHit.IsFEP();
    // double gammaEnDop = simHit.GetDoppler();
    // double gammaBeta = simHit.GetBeta();

    //SINGLES
    int gSize = gretina->Size();
    for (int i=0; i < gSize; i++){
      TGretinaHit &hit = gretina->GetGretinaHit(i);

      // double energy_corrected = hit.GetDopplerYta(s800sim->AdjustedBeta(GValue::Value("BETA")), yta, &track);
      // double energy = hit.GetDoppler(GValue::Value("BETA"));
      double core_energy = hit.GetCoreEnergy();
      double theta = hit.GetTheta();
      double phi = hit.GetPhi();
      int cryID = hit.GetCrystalId();
      if (cryID == 77) continue;
      
      obj.FillHistogram(dirname, "core_energy_prompt", 8192,0,8192, core_energy);
      // obj.FillHistogram(dirname, "gamma_singles_prompt", 8192,0,8192, energy);
      // obj.FillHistogram(dirname, "gamma_corrected_singles_prompt", 8192,0,8192, energy_corrected);
      // obj.FillHistogram(dirname, "gamma_corrected_vs_theta_prompt", 8192,0,8192, energy_corrected, 100, 0, 2.5, theta);
      // obj.FillHistogram(dirname, "gamma_corrected_vs_crystalID_prompt", 56, 24, 80, cryID, 8192,0,8192, energy_corrected);
      obj.FillHistogram(dirname, "core_energy_vs_theta_prompt", 8192,0,8192, hit.GetCoreEnergy(), 100, 0, 2.5, theta);
      
      //crystal angles
      // obj.FillHistogram(dirname, Form("theta_vs_crystalID_ring%d",gretina->GetRingNumber(cryID)),56,24,80,cryID,360,0,180,hit.GetThetaDeg());
      obj.FillHistogram(dirname, "gretina_theta_vs_phi",720,0,360,phi*TMath::RadToDeg(),360,0,180, theta*TMath::RadToDeg());

    }

    //NNADDBACK
    //loop over multiplicity
    for (int n=0; n<2; n++){
      //loop over hits for each multiplicity spectrum
      int nnSize = gretina->NNAddbackSize(n,false);
      for (int i=0; i < nnSize; i++){
        dirname = "basicsim";

        //get hit and hit data 
        TGretinaHit nnhit = gretina->GetNNAddbackHit(n,i);
        int cryID = nnhit.GetCrystalId();
        int ringNum = nnhit.GetRingNumber();
        double gEnergy = 0;
        if (stopped) gEnergy = nnhit.GetCoreEnergy();
        else gEnergy = nnhit.GetDopplerYta(s800sim->AdjustedBeta(GValue::Value("BETA")), yta, &track);

        double thresh = 100;
        if (nnhit.GetCoreEnergy() < thresh) continue;
        
        //exclude the ng spectrum (n==3)
        if (n < 3){
          obj.FillHistogram(dirname, "gamma_corrected_addback_prompt", 1600,0,1600, gEnergy);
        }

        if (n == 1) {
          if (nnhit.GetNeighbor().GetCoreEnergy() < thresh) continue;
          bool isN1FEP = fabs(gammaEn - gEnergy) < 1.5;

          //SCATTER TYPE
          // int ngroups = (int) scatterGroups.size();
          // for (int i=0; i < ngroups; i++){
          //   std::map<std::string, std::pair<int,int>>::iterator it = scatterGroups[i].begin();
          //   std::map<std::string, std::pair<int,int>>::iterator end = scatterGroups[i].end();
          //   while (it != end){
          //       if (checkScatterType(nnhit,it->second)){
          //         if (isN1FEP)
          //           obj.FillHistogram(dirname,Form("gamma_n1_FEP_grp%d_%s",i+1,it->first.c_str()),1600,0,1600, gEnergy);
          //         else 
          //           obj.FillHistogram(dirname,Form("gamma_n1_COMPT_grp%d_%s",i+1,it->first.c_str()),1600,0,1600, gEnergy);
          //         obj.FillHistogram(dirname,Form("gamma_n1_FEP+COMPT_grp%d_%s",i+1,it->first.c_str()),1600,0,1600, gEnergy);
          //       }
          //       it++;
          //   }
          // }

          //totals
          int id1 = nnhit.GetCrystalId();
          int id2 = nnhit.GetNeighbor().GetCrystalId();
          int quadType = 0;
          std::string scatType = "00";
          if (id1 > 43 && id2 > 43){ //only use 90 degree quads
            if ( PairHit(nnhit,TwoQuadPairs) ){
              quadType = 2;
              obj.FillHistogram(dirname,"gamma_n1_qd2_pair", 1600,0,1600, gEnergy);
              if(isN1FEP){
                obj.FillHistogram(dirname,"gamma_n1_FEP_qd2_pair", 1600,0,1600, gEnergy);
              }
            }
            if ( PairHit(nnhit,OneQuadPairs) ){
              quadType = 1;
              obj.FillHistogram(dirname,"gamma_n1_qd1_pair", 1600,0,1600, gEnergy);
              if(isN1FEP){
                obj.FillHistogram(dirname,"gamma_n1_FEP_qd1_pair", 1600,0,1600, gEnergy);
              }
            }

            if (id1%2 == 1 && id2%2 == 1){
              scatType = "AA";
              obj.FillHistogram(dirname,Form("gamma_n1_A-A"),1600,0,1600, gEnergy);
              if (isN1FEP){ 
                obj.FillHistogram(dirname,Form("gamma_n1_q%dAA_FEP",quadType),1600,0,1600, gEnergy);
                obj.FillHistogram(dirname,Form("gamma_n1_AA_FEP"),1600,0,1600, gEnergy);
                obj.FillHistogram(dirname,"crystal-map_A-A_FEP",360,0,360,nnhit.GetPhiDeg(),180,0,180,nnhit.GetThetaDeg());
                obj.FillHistogram(dirname,"crystal-map_A-A_FEP",360,0,360,nnhit.GetNeighbor().GetPhiDeg(),180,0,180,nnhit.GetNeighbor().GetThetaDeg());
              }
            }
            else if (id1%2 == 0 && id2%2 == 0){
              scatType = "BB";
              obj.FillHistogram(dirname,Form("gamma_n1_B-B"),1600,0,1600, gEnergy);
              if (isN1FEP){
                obj.FillHistogram(dirname,Form("gamma_n1_q%dBB_FEP",quadType),1600,0,1600, gEnergy);
                obj.FillHistogram(dirname,Form("gamma_n1_BB_FEP"),1600,0,1600, gEnergy);
                obj.FillHistogram(dirname,"crystal-map_B-B_FEP",360,0,360,nnhit.GetPhiDeg(),180,0,180,nnhit.GetThetaDeg());
                obj.FillHistogram(dirname,"crystal-map_B-B_FEP",360,0,360,nnhit.GetNeighbor().GetPhiDeg(),180,0,180,nnhit.GetNeighbor().GetThetaDeg());
              }
            }
            else{
              scatType = "AB";
              obj.FillHistogram(dirname,Form("gamma_n1_A-B"),1600,0,1600, gEnergy);
              if (isN1FEP){
                obj.FillHistogram(dirname,Form("gamma_n1_q%dAB_FEP",quadType),1600,0,1600, gEnergy);
                obj.FillHistogram(dirname,Form("gamma_n1_AB_FEP"),1600,0,1600, gEnergy);
                obj.FillHistogram(dirname,"crystal-map_A-B_FEP",360,0,360,nnhit.GetPhiDeg(),180,0,180,nnhit.GetThetaDeg());
                obj.FillHistogram(dirname,"crystal-map_A-B_FEP",360,0,360,nnhit.GetNeighbor().GetPhiDeg(),180,0,180,nnhit.GetNeighbor().GetThetaDeg());
                obj.FillHistogram(dirname,"gamma_n1_A-B_FEP_first_hit",2,0,2,id1%2);
              }
            }
          }

          std::string color = "nothing";
           //POLARIZATION
          if ( PairHit(nnhit,redPairs) ){
            color = "red";
            obj.FillHistogram(dirname,"gamma_n1_red_pair", 8192,0,8192, gEnergy);
            if(isN1FEP){
              obj.FillHistogram(dirname,"gamma_n1_red_pair_FEP", 8192,0,8192, gEnergy);
              obj.FillHistogram(dirname,Form("gamma_n1_red_pair_%s_FEP",scatType.c_str()), 2000,0,2000, gEnergy);
              obj.FillHistogram(dirname,Form("gamma_n1_red_pair_q%d_FEP",quadType), 2000,0,2000, gEnergy);
              obj.FillHistogram(dirname,Form("gamma_n1_red_pair_q%d%s_FEP",quadType,scatType.c_str()), 2000,0,2000, gEnergy);
              if (scatType == "AB"){
                obj.FillHistogram(dirname,"crystal-map_red_FEP",360,0,360,nnhit.GetPhiDeg(),180,0,180,nnhit.GetThetaDeg());
                obj.FillHistogram(dirname,"crystal-map_red_FEP",360,0,360,nnhit.GetNeighbor().GetPhiDeg(),180,0,180,nnhit.GetNeighbor().GetThetaDeg());
              }
            }
          }
          if ( PairHit(nnhit,goldPairs) ){
            color = "gold";
            obj.FillHistogram(dirname,"gamma_n1_gold_pair", 8192,0,8192, gEnergy);
            if(isN1FEP){
              obj.FillHistogram(dirname,"gamma_n1_gold_pair_FEP", 8192,0,8192, gEnergy);
              obj.FillHistogram(dirname,Form("gamma_n1_gold_pair_%s_FEP",scatType.c_str()), 2000,0,2000, gEnergy);
              obj.FillHistogram(dirname,Form("gamma_n1_gold_pair_q%d_FEP",quadType), 2000,0,2000, gEnergy);
              obj.FillHistogram(dirname,Form("gamma_n1_gold_pair_q%d%s_FEP",quadType,scatType.c_str()), 2000,0,2000, gEnergy);
              if (scatType == "AB"){
                obj.FillHistogram(dirname,"crystal-map_gold_FEP",360,0,360,nnhit.GetPhiDeg(),180,0,180,nnhit.GetThetaDeg());
                obj.FillHistogram(dirname,"crystal-map_gold_FEP",360,0,360,nnhit.GetNeighbor().GetPhiDeg(),180,0,180,nnhit.GetNeighbor().GetThetaDeg());
              }
            }
          }
          if ( PairHit(nnhit,bluePairs) ){
            color = "blue";
            obj.FillHistogram(dirname,"gamma_n1_blue_pair", 8192,0,8192, gEnergy);
            if(isN1FEP){
              obj.FillHistogram(dirname,"gamma_n1_blue_pair_FEP", 8192,0,8192, gEnergy);
              obj.FillHistogram(dirname,Form("gamma_n1_blue_pair_%s_FEP",scatType.c_str()), 2000,0,2000, gEnergy);
              obj.FillHistogram(dirname,Form("gamma_n1_blue_pair_q%d_FEP",quadType), 2000,0,2000, gEnergy);
              obj.FillHistogram(dirname,Form("gamma_n1_blue_pair_q%d%s_FEP",quadType,scatType.c_str()), 2000,0,2000, gEnergy);
              if (scatType == "AB"){
                obj.FillHistogram(dirname,"crystal-map_blue_FEP",360,0,360,nnhit.GetPhiDeg(),180,0,180,nnhit.GetThetaDeg());
                obj.FillHistogram(dirname,"crystal-map_blue_FEP",360,0,360,nnhit.GetNeighbor().GetPhiDeg(),180,0,180,nnhit.GetNeighbor().GetThetaDeg());
              }
            }
          }
          if (isN1FEP && id1 > 43 && id2 > 43){
            int pairCombo = id1*100 + id2;
            obj.FillHistogram(dirname,"n1_90degRing_pair_counts",100,0,50,SCATTERPAIRS[pairCombo]);
            obj.FillHistogram(dirname,Form("n1_90degRing_pair_counts_%s",color.c_str()),100,0,50,SCATTERPAIRS[pairCombo]);
            obj.FillHistogram(dirname,Form("n1_90degRing_pair_counts_%s_q%d",color.c_str(),quadType),100,0,50,SCATTERPAIRS[pairCombo]);
            obj.FillHistogram(dirname,Form("n1_90degRing_pair_counts_%s",scatType.c_str()),100,0,50,SCATTERPAIRS[pairCombo]);
            obj.FillHistogram(dirname,Form("n1_90degRing_pair_counts_%s_q%d",scatType.c_str(),quadType),100,0,50,SCATTERPAIRS[pairCombo]);
          }
        }

        char *multiplicity = Form("%d",n);
        if (n == 3) multiplicity = Form("g");
        obj.FillHistogram(dirname, Form("gamma_n%s",multiplicity), 1600,0,1600, gEnergy);

        if (n==0){
          // if (isFEP){
          //   dirname = "basicsim";
          //   obj.FillHistogram(dirname,"gamma_n0_prompt_fep", 1600,0,1600, gEnergy);
          //   // obj.FillHistogram("crystal-specific", Form("gamma_corrected_n0_ring%02d_crystal%d_prompt_fep",ringNum,cryID),1600,0,1600, gEnergy);
          // }
          if (cryID > 40){
            if (cryID%2 == 1) {
              obj.FillHistogram(dirname,"gamma_n0_A", 1600,0,1600, gEnergy);
              obj.FillHistogram(dirname, "crystal-map_A",360,0,360,nnhit.GetPhiDeg(),180,0,180,nnhit.GetThetaDeg());
              if (isFEP) {
                obj.FillHistogram(dirname,"gamma_n0_A_FEP", 1600,0,1600, gEnergy);
                obj.FillHistogram(dirname, "crystal-map_A_FEP",360,0,360,nnhit.GetPhiDeg(),180,0,180,nnhit.GetThetaDeg());
              }
            }
            else{
              obj.FillHistogram(dirname,"gamma_n0_B", 1600,0,1600, gEnergy);
              obj.FillHistogram(dirname, "crystal-map_B",360,0,360,nnhit.GetPhiDeg(),180,0,180,nnhit.GetThetaDeg());
              if (isFEP) {
                obj.FillHistogram(dirname,"gamma_n0_B_FEP", 1600,0,1600, gEnergy);
                obj.FillHistogram(dirname, "crystal-map_B_FEP",360,0,360,nnhit.GetPhiDeg(),180,0,180,nnhit.GetThetaDeg());
              }
            }
          }
          obj.FillHistogram(dirname, Form("gamma_corrected_n%s_vs_cryID",multiplicity),56, 24, 80, cryID, 1600,0,1600, gEnergy);
          // obj.FillHistogram(dirname, Form("gamma_corrected_n%s_cr%d",multiplicity,cryID), 1600,0,1600, gEnergy);
        }
      }
    }
  }

  if(numobj!=list->GetSize()){
    list->Sort();
  }
}
