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

// int LoadGates(TRuntimeObjects &obj, std::map<std::string, std::vector<GCutG*>> loaded_gates){
//   TList *gates = &(obj.GetGates());
//   TIter iter(gates);
//   std::cout << "loading gates:" <<std::endl;
//   while(TObject *obj = iter.Next()) {
//     GCutG *gate = (GCutG*)obj;
//     std::string tag = gate->GetTag();
//     loaded_gates[tag].push_back(gate);
//   }
//   std::cout << "outgoing size: " << outgoing_gates.size() << std::endl;
// }

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

  // TList *gates = &(obj.GetGates());
  // if(gates_loaded!=gates->GetSize()) {
    // LoadGates(obj);
  // }
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
    // double gammaEn = simHit.GetEn();
    // obj.FillHistogram(dirname,"GEANT_gretsim_size",20,0,20,gretsim->Size());
    // obj.FillHistogram(dirname,"GEANT_s800sim_size",20,0,20,s800sim->Size());
    // obj.FillHistogram(dirname,"GEANT_phi_vs_theta",360,0,360,simHit.GetPhi()*TMath::RadToDeg(),180,0,180,simHit.GetTheta()*TMath::RadToDeg());
    // bool isFEP = simHit.IsFEP();
    // double gammaEnDop = simHit.GetDoppler();
    // double gammaBeta = simHit.GetBeta();

    //SINGLES
    int gSize = gretina->Size();
    for (int i=0; i < gSize; i++){
      TGretinaHit &hit = gretina->GetGretinaHit(i);
      TGretinaHit hitCopy;
      hit.Copy(hitCopy);
      hit.SortSegments();
      EnergySmear(hit,rand_gen);
      double CUT = 0.0; 
      hit.ComptonSort(CUT);
      int nInteractions = hit.NumberOfInteractions();
      int nSegments = hit.GetNSegments();
      double energy_corrected = rand_gen->Gaus(hit.GetCoreEnergy(),hit.GetCoreEnergy()*0.0027/2.35);
      if (!stopped) energy_corrected = hit.GetDopplerYta(s800sim->AdjustedBeta(GValue::Value("BETA")), yta, &track);
      
      double core_energy = hit.GetCoreEnergy();
      double theta = hit.GetTheta();
      double phi = hit.GetPhi();
      int cryID = hit.GetCrystalId();
      int ringNum = hit.GetRingNumber();
      
      obj.FillHistogram(dirname,"gretina-map",180,0,180,theta*TMath::RadToDeg(),360,0,360,phi*TMath::RadToDeg());
      
      // double dopplerFactor = std::abs((theta-simTheta)*(gammaBeta*TMath::Sin(simTheta)/(1-gammaBeta*TMath::Cos(simTheta))));
      obj.FillHistogram(dirname, "core_energy", 8192,0,8192, core_energy);
      // obj.FillHistogram(dirname, "gamma_singles_prompt", 8192,0,8192, energy);
      obj.FillHistogram(dirname, "gam_dop_sgl_vs_theta", 180, 0, 180,theta*TMath::RadToDeg(), 2048,0,2048, energy_corrected);
      obj.FillHistogram(dirname, "gam_dop_sgl_vs_theta_cm", 180, 0, 180,thetaCM(theta,GValue::Value("BETA"))*TMath::RadToDeg(), 2048,0,2048, energy_corrected);
      
      if (simHit.IsFEP()) {
        obj.FillHistogram(dirname, "core_energy_FEP", 2048,0,2048, energy_corrected);
        obj.FillHistogram(dirname, "gamma_corrected_singles_FEP_theta_dist", 360,0,180,theta*TMath::RadToDeg());
        obj.FillHistogram(dirname, "gamma_corrected_singles_vs_nInteraction",12,0,12,nInteractions, 2048,0,2048, energy_corrected);
        obj.FillHistogram(dirname, "gamma_corrected_singles_vs_Efrac",200,0,1,hit.GetSegmentEng(0)/core_energy, 2048,0,2048, energy_corrected);
        obj.FillHistogram(dirname, "gamma_corrected_singles_vs_nSegments",12,0,12,nSegments, 2048,0,2048, energy_corrected);
        obj.FillHistogram(dirname, "gamma_corrected_singles_nSegments_vs_nInteraction",12,0,12,nInteractions,12,0,12,nSegments);
        
        if (nInteractions > 1){
          nEvents++;
          double xi = hit.GetXi(&track);
          double alpha = hit.GetAlpha();
          double scatterAngle = hit.GetScatterAngle();
          double nu = hit.GetScatterAngle(1,0); //swap the first two points
          double E2 = altDoppler(s800sim->AdjustedBeta(GValue::Value("BETA")),yta,track,hit,getTrueFirst(hit,hitCopy));

          double diffEratio = 511.0/core_energy * hit.GetSegmentEng(0)/(core_energy - hit.GetSegmentEng(0));
          obj.FillHistogram(dirname, "gamma_corrected_singles_FEP_nInt>2",2048,0,2048, energy_corrected);
          obj.FillHistogram(dirname,"gamma_corrected_singles_FEP_vs_xi",360,0,360,xi*TMath::RadToDeg(),2048,0,2048,energy_corrected);
          obj.FillHistogram(dirname,"gamma_corrected_singles_FEP_vs_scatterAngle",180,0,180,scatterAngle*TMath::RadToDeg(),2048,0,2048,energy_corrected);
          obj.FillHistogram(dirname,"gamma_corrected_singles_FEP_vs_nuprime",180,0,180,nu*TMath::RadToDeg(),2048,0,2048,energy_corrected);
          obj.FillHistogram(dirname,"gamma_corrected_singles_FEP_vs_E-E1/E1",1000,0,10,diffEratio,2048,0,2048,energy_corrected);
          obj.FillHistogram(dirname,"gamma_corrected_singles_FEP_cosScatterAngle_vs_diffEratio",1000,0,10,diffEratio,200,-1,1,TMath::Cos(scatterAngle));

          //print relative coords to first point
          // bool passed = TMath::Cos(scatterAngle) > - diffEratio;
          // TVector3 v1 = hitCopy.GetIntPosition(0);
          // TVector3 v2 = TVector3(0,-v1.Z()/v1.Y(),1);
          // TVector3 v3 = v1.Cross(v2);
          // double mtxData[9] = {v1.Unit().X(),v1.Unit().Y(),v1.Unit().Z(),v2.Unit().X(),v2.Unit().Y(),v2.Unit().Z(),v3.Unit().X(),v3.Unit().Y(),v3.Unit().Z()};
          // TMatrixT <double> mtx = TMatrixT<double> (3,3,mtxData,"F");
          // mtx.Invert();
          // TVector3 pos;
          // if (nEvents < 10001){
          //   static FILE *pFile;
          //   pFile = fopen("scatter_coords.dat","a");
          //   if (nEvents == 1) fprintf(pFile,"evt,nint,cut,eng,x,y,z\n");
          //   for (int np=0; np < nInteractions; np++){
          //     pos = mtx*(hitCopy.GetIntPosition(np));
          //     fprintf(pFile,"%d,%d,%d,%f,%f,%f,%f\n",nEvents-1,nInteractions,passed,hitCopy.GetSegmentEng(np),pos.X(),pos.Y(),pos.Z());
          //   }
          //   fclose(pFile);
          // } 

          //NU GATED
          // if (nu > 1.5*TMath::ACos(1-1.0/(nInteractions-1) * 511.0/core_energy)){
          if (diffEratio > 2){
            obj.FillHistogram(dirname,"gamma_corrected_singles_COMPT>2GATE_FEP_vs_xi",360,0,360,xi*TMath::RadToDeg(),2048,0,2048,energy_corrected);
            obj.FillHistogram(dirname,"gamma_corrected_singles_COMPT>2GATE_FEP_vs_E-E1/E1",1000,0,10,diffEratio,2048,0,2048,energy_corrected);
            obj.FillHistogram(dirname,"gamma_corrected_singles_COMPT>2GATE_FEP_cosScatterAngle_vs_E-E1/E1",1000,0,10,diffEratio,200,-1,1,TMath::Cos(scatterAngle));
          } else if (TMath::Cos(scatterAngle) > 0.5 - diffEratio) {
            obj.FillHistogram(dirname,"gamma_corrected_singles_COMPTGATE_FEP_vs_xi",360,0,360,xi*TMath::RadToDeg(),2048,0,2048,energy_corrected);
            obj.FillHistogram(dirname,"gamma_corrected_singles_COMPTGATE_FEP_vs_E-E1/E1",1000,0,10,diffEratio,2048,0,2048,energy_corrected);
            obj.FillHistogram(dirname,"gamma_corrected_singles_COMPTGATE_FEP_vs_AltDop",2048,0,2048,E2,2048,0,2048,energy_corrected);
          } else {
            obj.FillHistogram(dirname,"gamma_corrected_singles_COMPTNOTGATE_FEP_vs_xi",360,0,360,xi*TMath::RadToDeg(),2048,0,2048,energy_corrected);
            obj.FillHistogram(dirname,"gamma_corrected_singles_COMPTNOTGATE_FEP_vs_E-E1/E1",1000,0,10,diffEratio,2048,0,2048,energy_corrected);
            obj.FillHistogram(dirname,"gamma_corrected_singles_COMPTNOTGATE_FEP_vs_AltDop",2048,0,2048,E2,2048,0,2048,energy_corrected);
          }
        }
      }
    }

    //NNADDBACK
    //loop over multiplicity
    // for (int n=0; n<3; n++){
    //   //loop over hits for each multiplicity spectrum
    //   int nnSize = gretina->NNAddbackSize(n);
    //   for (int i=0; i < nnSize; i++){
    //     dirname = "basicsim";

    //     //get hit and hit data 
    //     TGretinaHit nnhit = gretina->GetNNAddbackHit(n,i);
    //     bool isFEP = (fabs(gammaEn - nnhit.GetCoreEnergy()) < 1.5);
    //     double trueTheta = nnhit.GetThetaDeg();
    //     double truePhi = nnhit.GetPhiDeg();
    //     double dummy;
    //     double trueXi = azimuthalCompton(nnhit,&track,dummy)*TMath::RadToDeg();
    //     double trueEnergy = nnhit.GetDopplerYta(s800sim->AdjustedBeta(GValue::Value("BETA")), yta, &track);
    //     int nInteractions = nnhit.NumberOfInteractions();
    //     // if (std::abs(trueEnergy - 1018) > 4){
    //     //   double rel_beta = s800sim->AdjustedBeta(GValue::Value("BETA"));
    //     //   double rel_gamma = 1./(sqrt(1.-pow(rel_beta,2.)));
    //     //   TVector3 gret_pos = nnhit.GetIntPosition(2);
    //     //   gret_pos.SetY(gret_pos.Y() + yta);
    //     //   trueEnergy = nnhit.GetCoreEnergy()*rel_gamma *(1 - rel_beta*TMath::Cos(gret_pos.Angle(track)));
    //     // }
    //     if (std::abs(trueEnergy - 1018) > 4 && isFEP && nInteractions > 1 && n < 3){
    //       obj.FillHistogram(dirname,"bad_events_theta_vs_phi",720,0,360,truePhi,360,0,180,trueTheta);
    //       obj.FillHistogram(dirname,"bad_events_TrueEnergy_vs_TrueXi",360,0,360,trueXi,2048,0,2048,trueEnergy);
    //       obj.FillHistogram(dirname,"bad_events_energy_vs_ninteractions>1",10,0,10,nInteractions,2048,0,2048,trueEnergy);
    //     }
        
    //     EnergySmear(nnhit,rand_gen);
        
    //     int cryID = nnhit.GetCrystalId();
    //     int ringNum = nnhit.GetRingNumber();
    //     double theta = nnhit.GetThetaDeg();
    //     double phi = nnhit.GetPhiDeg();
    //     // get the energy depending on whether it's a sim stopped beam
    //     double gEnergy = 0;
    //     if (stopped) gEnergy = nnhit.GetCoreEnergy();
    //     else gEnergy = nnhit.GetDopplerYta(s800sim->AdjustedBeta(GValue::Value("BETA")), yta, &track);

    //     //only do FEP events
    //     if (!isFEP) continue;

    //     double thresh = 90;
    //     if (nnhit.GetCoreEnergy() < thresh) continue;
        
    //     if (nInteractions == 2) {
    //       obj.FillHistogram(dirname, "gamma_corrected_2interactions_prompt", 1600,0,1600, gEnergy);
    //       obj.FillHistogram(dirname, Form("gamma_corrected_2int_n%d",n), 1600,0,1600, gEnergy);
    //       if (n == 0){
    //         double int1E = nnhit.GetSegmentEng(0);
    //         double int2E = nnhit.GetSegmentEng(1);
    //         obj.FillHistogram(dirname,"ip1_vs_ip2_energy", 1600,0,1600, int2E, 1600,0,1600, int1E);
    //         if (int1E > int2E) obj.FillHistogram(dirname, "gamma_corrected_ip1>ip2", 1600,0,1600, gEnergy);
    //         else obj.FillHistogram(dirname, "gamma_corrected_ip1<ip2", 1600,0,1600, gEnergy);
    //       }
    //     }
    //     if (nInteractions > 1 && n < 3){
    //       double dp;
    //       double xi = azimuthalCompton(nnhit,&track,dp)*TMath::RadToDeg();
    //       // obj.FillHistogram(dirname,"energy_vs_ninteractions>1",16,0,16,nInteractions,2048,0,2048,gEnergy);
    //       obj.FillHistogram(dirname,"energy_vs_dot_product",440,-1.1,1.1,dp,2048,0,2048,gEnergy);
    //       obj.FillHistogram(dirname,"dp_vs_xi",360,0,360,xi,400,-2,2,dp); 
    //       obj.FillHistogram(dirname,"energy_vs_xi",360,0,360,xi,1024,0,2048,gEnergy);
    //       obj.FillHistogram(dirname,"beta_vs_xi",360,0,360,xi,200,0.3,0.5,gammaBeta);
    //       // obj.FillHistogram(dirname,"trueTheta_vs_trueXi",360,0,360,trueXi,180,0,180,trueTheta);
    //       // obj.FillHistogram(dirname,"theta_diff_vs_xi_diff",400,-400,400,trueXi - xi,500,-50,50,trueTheta - theta);
    //       // obj.FillHistogram(dirname,"energy_diff_vs_theta_diff",500,-50,50,theta - trueTheta,600,-300,300,gEnergy - trueEnergy);
    //       obj.FillHistogram(dirname,"energy_vs_TrueXi",360,0,360,trueXi,2048,0,2048,gEnergy);
    //       obj.FillHistogram(dirname,"TrueEnergy_vs_TrueXi",360,0,360,trueXi,2048,0,2048,trueEnergy);
    //       obj.FillHistogram(dirname,"(xi-TrueXi)_vs_(theta-TrueTheta)",200,-50,50,theta-trueTheta,760,-380,380,xi-trueXi);
    //       // obj.FillHistogram(dirname,"phi_vs_xi",360,0,360,xi,360,0,360,phi);
    //       // obj.FillHistogram(dirname,"xi_vs_trueXi",360,0,360,trueXi,360,0,360,xi);

    //     } else {
    //       obj.FillHistogram(dirname, "gamma_corrected_oneIntPoint_prompt", 1600,0,1600, gEnergy);
    //     }
        
    //     //exclude the ng spectrum (n==3)
    //     if (n < 3){
    //       obj.FillHistogram(dirname, "gamma_corrected_addback_prompt", 1600,0,1600, gEnergy);
    //     }

    //     char *multiplicity = Form("%d",n);
    //     if (n == 3) multiplicity = Form("g");
    //     obj.FillHistogram(dirname, Form("gamma_n%s",multiplicity), 1600,0,1600, gEnergy);
    //   }
    // }
  }

  if(numobj!=list->GetSize()){
    list->Sort();
  }
}
