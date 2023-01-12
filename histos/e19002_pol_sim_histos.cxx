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

double azimuthalCompton(const TGretinaHit &hit, const TVector3 *beam){
  TVector3 interaction1 = hit.GetIntPosition(0);
  TVector3 interaction2 = hit.GetIntPosition(1);
  TVector3 comptonPlaneNorm = interaction1.Cross(interaction2);
  TVector3 reactionPlaneNorm = beam->Cross(interaction1);
  TVector3 basisNorm = interaction1.Cross(reactionPlaneNorm);
  double angle = reactionPlaneNorm.Angle(comptonPlaneNorm);
  if (basisNorm.Angle(comptonPlaneNorm) > TMath::PiOver2()) angle = TMath::TwoPi() - angle;
  return angle;
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
  
  if (gretina){
    TGretSimHit simHit = gretsim->GetGretinaSimHit(0);
    double gammaEn = simHit.GetEn();
    bool isFEP = simHit.IsFEP();
    double gammaEnDop = simHit.GetDoppler();
    double SIGMA = (2.1*TMath::Exp(-0.1*gammaEn/1000.0) + 60.0*TMath::Exp(-10.2*gammaEn/1000.0));
    if(SIGMA > 3.8) {
      SIGMA = 3.8;
    }
    if(gammaEn < 150.) {
      SIGMA = 6.0;
    }
    // double gammaBeta = simHit.GetBeta();

    //SINGLES
    int gSize = gretina->Size();
    for (int i=0; i < gSize; i++){
      TGretinaHit &hit = gretina->GetGretinaHit(i);

      TVector3 local_pos(hit.GetLocalPosition(0));
      double smear_x = local_pos.X() + rand_gen->Gaus(0, SIGMA); 
      double smear_y = local_pos.Y() + rand_gen->Gaus(0, SIGMA);
      double smear_z = local_pos.Z() + rand_gen->Gaus(0, SIGMA);
      hit.SetPosition(0,smear_x,smear_y,smear_z);

      double energy_corrected = hit.GetDopplerYta(s800sim->AdjustedBeta(GValue::Value("BETA")), yta, &track);
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
    for (int n=0; n<3; n++){
      //loop over hits for each multiplicity spectrum
      int nnSize = gretina->NNAddbackSize(n,true);
      for (int i=0; i < nnSize; i++){
        dirname = "basicsim";

        //get hit and hit data 
        TGretinaHit nnhit = gretina->GetNNAddbackHit(n,i);
        int cryID = nnhit.GetCrystalId();
        int ringNum = nnhit.GetRingNumber();
        double theta = nnhit.GetThetaDeg();
        double phi = nnhit.GetPhiDeg();
        int nInteractions = nnhit.NumberOfInteractions();
        // get the energy depending on whether it's a sim stopped beam
        double gEnergy = 0;
        if (stopped) gEnergy = nnhit.GetCoreEnergy();
        else gEnergy = nnhit.GetDopplerYta(s800sim->AdjustedBeta(GValue::Value("BETA")), yta, &track);

        //only do FEP events
        if (!(fabs(gammaEn - nnhit.GetCoreEnergy()) < 1.5)) continue;

        double thresh = 90;
        if (nnhit.GetCoreEnergy() < thresh) continue;

        if (nInteractions > 1){
          double aziCompt = azimuthalCompton(nnhit,&track)*TMath::RadToDeg();
          obj.FillHistogram(dirname, "azmthl_compton",360,0,360,aziCompt,1024,0,2048,gEnergy);
          if (aziCompt > 90 && aziCompt < 270){
            obj.FillHistogram(dirname,"Xi_center_gamma_vs_theta",180,0,180,theta,8192,0,8192,gEnergy);
            obj.FillHistogram(dirname,"Xi_center_gamma_vs_phi",360,0,360,phi,8192,0,8192,gEnergy);
          } else {
            obj.FillHistogram(dirname,"Xi_extreme_gamma_vs_theta",180,0,180,theta,8192,0,8192,gEnergy);
            obj.FillHistogram(dirname,"Xi_extreme_gamma_vs_phi",360,0,360,phi,8192,0,8192,gEnergy);
          }
          if (theta > 100) 
            obj.FillHistogram(dirname, "azmthl_compton_theta>100",360,0,360,aziCompt,1024,0,2048,gEnergy);
          else if (theta < 50) 
            obj.FillHistogram(dirname, "azmthl_compton_theta<50",360,0,360,aziCompt,1024,0,2048,gEnergy);
          else if (theta > 75 && theta < 95) 
            obj.FillHistogram(dirname, "azmthl_compton_theta_flat",360,0,360,aziCompt,1024,0,2048,gEnergy);
          else 
            obj.FillHistogram(dirname, "azmthl_compton_theta_varied",360,0,360,aziCompt,1024,0,2048,gEnergy);
        }
        
        //exclude the ng spectrum (n==3)
        if (n < 3){
          obj.FillHistogram(dirname, "gamma_corrected_addback_prompt", 1600,0,1600, gEnergy);
        }

        char *multiplicity = Form("%d",n);
        if (n == 3) multiplicity = Form("g");
        obj.FillHistogram(dirname, Form("gamma_n%s",multiplicity), 1600,0,1600, gEnergy);
      }
    }
  }

  if(numobj!=list->GetSize()){
    list->Sort();
  }
}
