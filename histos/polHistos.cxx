#include "TRuntimeObjects.h"

#include <iostream>
#include <map>
#include <cstdio>
#include <unordered_set>

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


double GetAfp(double crdc_1_x,double  crdc_2_x){
  return TMath::ATan( (crdc_2_x - crdc_1_x)/1073.0 );
}

void EnergySmear(TGretinaHit &hit, TRandom3 *rand){
  int nPoints = hit.NumberOfInteractions();
  // if (nPoints > 1) nPoints = 2; 
  for (int i=0; i<nPoints; i++){
    double E = hit.GetSegmentEng(i);
    double sigma = (2.1*TMath::Exp(-0.1*E/1000.0) + 60.0*TMath::Exp(-10.2*E/1000.0));
    if(sigma > 3.8) {
      sigma = 3.8;
    }
    if(sigma > 6) {
      sigma = 6.0;
    }
    TVector3 local_pos(hit.GetLocalPosition(i));
    double smear_x = local_pos.X() + rand->Gaus(0, sigma); 
    double smear_y = local_pos.Y() + rand->Gaus(0, sigma);
    double smear_z = local_pos.Z() + rand->Gaus(0, sigma);
    hit.SetPosition(i,smear_x,smear_y,smear_z);
  }
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

bool gates_loaded = false;
std::map<std::string,std::vector<GCutG*>> gates;

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

  //load in the gates
  if (!gates_loaded) {
    LoadGates(&(obj.GetGates()),gates);
    gates_loaded = true;
  }

  static TRandom3 *rand_gen = new TRandom3(59953614);
  // static int nEvents = 0;
  bool stopped = false;

  if (gretsim){
    int nhits = gretsim->Size();
    obj.FillHistogram("GEANT", "GEANT_gamma_hits",20,0,20,nhits);
    for (int i=0; i < nhits; i++){
      TGretSimHit simHit = gretsim->GetGretinaSimHit(i);
      double simTheta = simHit.GetTheta();
      obj.FillHistogram("GEANT", "simtheta_dist", 360,0,180,simTheta*TMath::RadToDeg());
      obj.FillHistogram("GEANT", "simBeta", 100,0.3,0.5,simHit.GetBeta());
    }
  }

  if (!gretina || gretsim->Size() == 0){
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
    TGretSimHit simHit = gretsim->GetGretinaSimHit(0);
    double simTheta = simHit.GetTheta();
    obj.FillHistogram(dirname,"GEANT_energies",1000,500,1500,simHit.GetDoppler(0,&track));

    //SINGLES
    int gSize = gretina->Size();
    for (int i=0; i < gSize; i++){
      TGretinaHit &hit = gretina->GetGretinaHit(i);
      // hit.SortSegments();
      // hit.ReverseSegments();
      // EnergySmear(hit,rand_gen);
      
      int nInteractions = hit.NumberOfInteractions();
      double energy_corrected = rand_gen->Gaus(hit.GetCoreEnergy(),hit.GetCoreEnergy()*0.0027/2.35);
      if (!stopped) energy_corrected = hit.GetDopplerYta(s800sim->AdjustedBeta(simHit.GetBeta()), yta, &track);
      
      double core_energy = hit.GetCoreEnergy();
      double theta = hit.GetTheta();
      double phi = hit.GetPhi();
      int cryID = hit.GetCrystalId();
      
      obj.FillHistogram(dirname,"gretina-map",180,0,180,theta*TMath::RadToDeg(),360,0,360,phi*TMath::RadToDeg());
      obj.FillHistogram(dirname, "core_energy", 8192,0,8192, core_energy);
      obj.FillHistogram(dirname, "theta_vs_core_energy",8192,0,8192, core_energy,100,0.5,2.5,theta);
      obj.FillHistogram(dirname, "gam_dop_sgl_vs_theta", 180, 0, 180,theta*TMath::RadToDeg(), 4096,0,4096, energy_corrected);
      
      if (nInteractions > 1) obj.FillHistogram(dirname,"gamma_corrected_singles_nInt>2",4096,0,4096, energy_corrected);
      if (simHit.IsFEP()) {
        obj.FillHistogram(dirname, "FEP_theta_dist",180,0,180,theta*TMath::RadToDeg());
        obj.FillHistogram(dirname, "core_energy_FEP", 4096,0,4096, core_energy);
        obj.FillHistogram(dirname, "gamma_corrected_singles_FEP_theta_dist", 360,0,180,theta*TMath::RadToDeg());
        obj.FillHistogram(dirname, "gamma_corrected_singles_vs_nInteraction",12,0,12,nInteractions, 4096,0,4096, energy_corrected);
        obj.FillHistogram(dirname, "gamma_corrected_singles_vs_Efrac",200,0,1,hit.GetSegmentEng(0)/core_energy, 4096,0,4096, energy_corrected);
        
        if (nInteractions > 1){
          double xi = 0.0; 
          if (stopped) xi = hit.GetXi();
          else xi = hit.GetXi(&track);

          double scatterAngle = hit.GetScatterAngle();
          double diffEratio = 511.0/core_energy * hit.GetSegmentEng(0)/(core_energy - hit.GetSegmentEng(0));

          obj.FillHistogram(dirname, "nopolgate_FEP_scatterAngle",180,0,TMath::Pi(),scatterAngle);
          obj.FillHistogram(dirname, "core_energy_FEP_vs_xi",360,0,TMath::TwoPi(),xi, 4096,0,4096, core_energy);
          obj.FillHistogram(dirname, "nopolgate_FEP_Edop_vs_xi",360,0,TMath::TwoPi(),xi,4096,0,4096,energy_corrected);
          obj.FillHistogram(dirname, "nopolgate_FEP_theta_vs_xi",360,0,TMath::TwoPi(),xi,360,0,TMath::Pi(),theta);
          obj.FillHistogram(dirname, "nopolgate_FEP_phi_vs_xi",360,0,TMath::TwoPi(),xi,360,0,TMath::TwoPi(),phi);

          // check if all interaction points are > 1.5 degrees away
          bool allBad = true;
          for (int ip=0; ip < nInteractions; ip++){
            if (std::abs(simTheta-hit.GetTheta(ip))*TMath::RadToDeg() < 1.5) {
              allBad = false;
              break;
            }
          }
          if (allBad) obj.FillHistogram(dirname,"evt_noGoodPoints_E",300,1200,1500,energy_corrected);
          else obj.FillHistogram(dirname,"evt_goodPoints_E",300,1200,1500,energy_corrected);
        }
      } 
    }
  }

  if(numobj!=list->GetSize()){
    list->Sort();
  }
}
