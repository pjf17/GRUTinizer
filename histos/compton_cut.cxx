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

double altDoppler(double beta, double yta, const TVector3 &beam, const TGretinaHit &hit, int pos){
  double gamma = 1./(sqrt(1.-pow(beta,2.)));
  TVector3 gret_pos = hit.GetIntPosition(pos);
  gret_pos.SetY(gret_pos.Y() + yta);
  return hit.GetCoreEnergy()*gamma*(1 - beta*TMath::Cos(gret_pos.Angle(beam)));
}

double thetaCM(double theta,double beta){
  double cosT = TMath::Cos(theta);
  return TMath::ACos((cosT - beta)/(1 - beta*cosT));
}

int testXi(double xi, double gatewidth){
    //0 is parallel 1 is perpendicular
    int flag = -1;
    double gw = gatewidth/2;
    if ((xi >= 90-gw && xi <= 90+gw) || (xi >= 270-gw && xi <= 270+gw)){
        flag = 1;
    }
    else if ((xi >= 360-gw || xi <= gw) || (xi >= 180-gw && xi <= 180+gw)){
        flag = 0;
    }
    return flag;
}

int getClosestPoint(const TGretinaHit &hit){
  int nInteractions = hit.NumberOfInteractions();
  double length = hit.GetIntPosition(0).Mag();
  int shortestPoint = 0;
  for (int i=1; i < nInteractions; i++){
    double comp = hit.GetIntPosition(i).Mag();
    if (length > comp){
      length = comp;
      shortestPoint = i;
    }
  }
  return shortestPoint;
}

int getTrueFirst(const TGretinaHit &hitSm, const TGretinaHit &hitTr){
  int nInteractions = hitSm.NumberOfInteractions();
  double diff = 200;
  double target = hitTr.GetSegmentEng(0);
  int p =0;
  for (int i=0; i < nInteractions;i++){
    double comp = std::abs(hitSm.GetSegmentEng(i)-target);
    if (comp < diff) {
      diff = comp;
      p = i;
    }
  }
  return p;
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
  static int nEvents = 0;
  bool stopped = false;

  if (gretsim){
    int nhits = gretsim->Size();
    obj.FillHistogram("GEANT", "GEANT_gamma_hits",20,0,20,nhits);
    for (int i=0; i < nhits; i++){
      TGretSimHit simHit = gretsim->GetGretinaSimHit(i);
      double simTheta = simHit.GetTheta();
      obj.FillHistogram("GEANT", "simtheta_dist", 360,0,180,simTheta*TMath::RadToDeg());
      obj.FillHistogram("GEANT", "simtheta_cm_dist", 360,0,180,thetaCM(simTheta,simHit.GetBeta())*TMath::RadToDeg());
      if (simTheta*TMath::RadToDeg() > 40 && simTheta*TMath::RadToDeg() <120)
        obj.FillHistogram("GEANT", "simtheta_cm_limited_dist", 360,0,180,thetaCM(simTheta,simHit.GetBeta())*TMath::RadToDeg());
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
    // double gammaEnDop = simHit.GetDoppler();
    // double gammaBeta = simHit.GetBeta();

    //SINGLES
    int gSize = gretina->Size();
    for (int i=0; i < gSize; i++){
      TGretinaHit &hit = gretina->GetGretinaHit(i);
      hit.SortSegments();
      EnergySmear(hit,rand_gen);

      for (int i=0; i < 6; i++){
        double CUT = 2 - i*0.5;
        
        TGretinaHit hitCopy;
        hit.Copy(hitCopy);
        // hitCopy.ComptonSort(CUT);
        
        double energy_corrected = energy_corrected = hitCopy.GetDopplerYta(s800sim->AdjustedBeta(GValue::Value("BETA")), yta, &track);
        
        if (simHit.IsFEP() && hit.NumberOfInteractions() > 1) {
          obj.FillHistogram(dirname,Form("gamma_corrected_singles_FEP_cut%2.2f",CUT),2048,0,2048, energy_corrected);
        } 
        if (hit.NumberOfInteractions() > 1){
          obj.FillHistogram(dirname,Form("gamma_corrected_singles_cut%2.2f",CUT),2048,0,2048, energy_corrected);
        }
      }
    }
  }

  if(numobj!=list->GetSize()){
    list->Sort();
  }
}
