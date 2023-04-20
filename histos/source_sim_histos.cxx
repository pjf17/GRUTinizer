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

double azimuthalCompton(const TGretinaHit &hit, const TVector3 *beam, double &dp, bool doSwap=false){
  TVector3 interaction1 = hit.GetIntPosition(0);
  TVector3 interaction2 = hit.GetIntPosition(1);
  if (doSwap) std::swap(interaction1,interaction2);
  TVector3 comptonPlaneNorm = interaction1.Cross(interaction2);
  TVector3 reactionPlaneNorm = beam->Cross(interaction1);
  TVector3 basisNorm = interaction1.Cross(reactionPlaneNorm);
  dp = reactionPlaneNorm.Dot(comptonPlaneNorm)/reactionPlaneNorm.Mag()/comptonPlaneNorm.Mag();
  double angle = reactionPlaneNorm.Angle(comptonPlaneNorm);
  if (basisNorm.Angle(comptonPlaneNorm) > TMath::PiOver2()) angle = TMath::TwoPi() - angle;
  return angle;
}

void EnergySmear(TGretinaHit &hit, TRandom3 *rand){
  double E = hit.GetCoreEnergy();
  double sigma = (2.1*TMath::Exp(-0.1*E/1000.0) + 60.0*TMath::Exp(-10.2*E/1000.0));
  if(sigma > 3.8) {
    sigma = 3.8;
  }
  if(E < 150.) {
    sigma = 6.0;
  }
  TVector3 local_pos(hit.GetLocalPosition(0));
  double smear_x = local_pos.X() + rand->Gaus(0, sigma); 
  double smear_y = local_pos.Y() + rand->Gaus(0, sigma);
  double smear_z = local_pos.Z() + rand->Gaus(0, sigma);
  hit.SetPosition(0,smear_x,smear_y,smear_z);
}

double altDoppler(double beta, const TVector3 &beam, double yta, const TGretinaHit &hit, int pos){
  double gamma = 1./(sqrt(1.-pow(beta,2.)));
  TVector3 gret_pos = hit.GetIntPosition(pos);
  gret_pos.SetY(gret_pos.Y() + yta);
  return hit.GetCoreEnergy()*gamma*(1 - beta*TMath::Cos(gret_pos.Angle(beam)));
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
  static TRandom3 *rand_gen = new TRandom3(59953614);
  bool stopped = false;

  if (gretsim){
    int nhits = gretsim->Size();
    if (nhits == 2){
      TGretSimHit simHit = gretsim->GetGretinaSimHit(0);
      TGretSimHit simOtherHit = gretsim->GetGretinaSimHit(1);
      double simTheta = simHit.GetTheta();
      double simAngle = simHit.fPosit.Angle(simOtherHit.fPosit);
      obj.FillHistogram("GEANT", "simtheta_dist", 360,0,180,simTheta*TMath::RadToDeg());
      obj.FillHistogram("GEANT", "sim_angular_correlation", 360,0,180,simAngle*TMath::RadToDeg());
    }
  }

  if (!gretina || gretsim->Size() < 2){
    return;
  }
  
  if (!s800sim || !s800sim->Size() || GValue::Value("BETA") == 0.0){
    stopped = true;
  }
  
  std::string dirname("basicsim");
  TVector3 track;
  double yta, ata, bta, dta=0;
  if (!stopped){
    track = s800sim->Track(0,0);
    yta = s800sim->GetS800SimHit(0).GetYTA();
    ata = s800sim->GetS800SimHit(0).GetATA();
    bta = s800sim->GetS800SimHit(0).GetBTA();
    dta = s800sim->GetS800SimHit(0).GetDTA();
    //S800 coordinates
    obj.FillHistogram("s800sim","bta_vs_ata", 600,-0.1,0.1, ata, 600,-0.1,0.1, bta);
    obj.FillHistogram("s800sim","yta", 1000,-0.003,0.003, yta);
    obj.FillHistogram("s800sim","dta", 1000,-0.5,0.5, dta);
  } 
  else track = TVector3(0,0,1);

  bool atabtaGate = std::abs(ata) < 0.002 && std::abs(bta) < 0.002;
  bool ytaGate = std::abs(yta) < 0.0002;
  if (gretina){
    TGretSimHit simHit = gretsim->GetGretinaSimHit(1);
    // TGretSimHit simHit = gretsim->GetGretinaSimHit(1);
    // track = gretsim->GetGretinaSimHit(0).fPosit;

    //SINGLES
    int gSize = gretina->Size();
    for (int i=0; i < gSize; i++){
      TGretinaHit &hit = gretina->GetGretinaHit(i);
      // hit.SortSegments();
      // EnergySmear(hit,rand_gen);
      double energy_corrected = rand_gen->Gaus(hit.GetCoreEnergy(),hit.GetCoreEnergy()*0.0027/2.35);
      if (!stopped) energy_corrected = hit.GetDopplerYta(s800sim->AdjustedBeta(GValue::Value("BETA")), yta, &track);
      bool is1332FEP = abs(simHit.GetEn() - hit.GetCoreEnergy()) < 1; 

      double core_energy = hit.GetCoreEnergy();
      double theta = hit.GetTheta();
      double phi = hit.GetPhi();
      int nInteractions = hit.NumberOfInteractions();
      int cryID = hit.GetCrystalId();
      int hole = hit.GetHoleNumber();
      obj.FillHistogram(dirname,"gretina-map",180,0,180,theta*TMath::RadToDeg(),360,0,360,phi*TMath::RadToDeg());
      // obj.FillHistogram(dirname,Form("gretina-map_cr%d",cryID),180,0,180,theta*TMath::RadToDeg(),360,0,360,phi*TMath::RadToDeg());
      
      obj.FillHistogram(dirname, "core_energy", 4096,0,4096, core_energy);
      obj.FillHistogram(dirname, "core_energy_dop", 4096,0,4096, hit.GetDoppler(0.35));
      obj.FillHistogram(dirname, "core_energy_vs_theta",180,0,180,theta*TMath::RadToDeg(),4096,0,4096, core_energy);
      
      // if (is1332FEP){
      //   obj.FillHistogram(dirname, "gamma_corrected_singles_FEP_theta_dist", 360,0,180,theta*TMath::RadToDeg());
      //   obj.FillHistogram(dirname, Form("gamma_corrected_singles_FEP_type%d",cryID%2), 2048,0,2048, energy_corrected);
      // }
      obj.FillHistogram(dirname, "gamma_corrected_singles_vs_nInteraction",12,0,12,nInteractions, 2048,0,2048, energy_corrected);

      if (nInteractions > 1){
        double dummy;
        double xi = azimuthalCompton(hit,&track,dummy);
        if (!simHit.IsFEP()) obj.FillHistogram(dirname, "gamma_corrected_singles_FEP_vs_xi_1000counts",360,0,360,xi*TMath::RadToDeg(), 2048,0,2048, hit.GetDoppler(0.4474));
        if (!simHit.IsFEP()) obj.FillHistogram(dirname, "gamma_corrected_singles_FEP_vs_xi_2000counts",360,0,360,xi*TMath::RadToDeg(), 2048,0,2048, hit.GetDoppler(0.4474));
        if (!simHit.IsFEP()) obj.FillHistogram(dirname, "gamma_corrected_singles_FEP_vs_xi_4000counts",360,0,360,xi*TMath::RadToDeg(), 2048,0,2048, hit.GetDoppler(0.4474));
        if (!simHit.IsFEP()) obj.FillHistogram(dirname, "gamma_corrected_singles_FEP_vs_xi_8000counts",360,0,360,xi*TMath::RadToDeg(), 2048,0,2048, hit.GetDoppler(0.4474));
        // obj.FillHistogram(dirname,"gamma_corrected_singles_vs_xi",360,0,360,xi*TMath::RadToDeg(),2048,0,2048,energy_corrected);
        // obj.FillHistogram(dirname,"gamma_corrected_singles_xi",360,0,360,xi*TMath::RadToDeg());
        // if (is1332FEP){
        //   // obj.FillHistogram(dirname,Form("gamma_corrected_singles_FEP_xi_qd%d",hole),360,0,360,xi*TMath::RadToDeg());
        //   if (theta*TMath::RadToDeg() >= 80 && theta*TMath::RadToDeg() <= 100){
        //     obj.FillHistogram(dirname,"gamma_corrected_singles_FEP_vs_xi_theta_gate",360,0,360,xi*TMath::RadToDeg(),2048,0,2048,energy_corrected);
        //     obj.FillHistogram(dirname,"gamma_corrected_singles_FEP_xi_theta_gate",360,0,360,xi*TMath::RadToDeg()); 
        //   }
        //   for (int i=60; i<=100; i+=10){
        //     int hw = 5; //half width
        //     if (theta*TMath::RadToDeg() >= i-hw && theta*TMath::RadToDeg() <= i+hw){
        //       obj.FillHistogram(dirname,Form("xi_theta_%d_%d",i-hw,i+hw),360,0,360,xi*TMath::RadToDeg());
        //     }
        //   } 
        // }
      }
    }
  }

  if(numobj!=list->GetSize()){
    list->Sort();
  }
}
