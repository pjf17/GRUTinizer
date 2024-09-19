#include "TRuntimeObjects.h"

#include <iostream>
#include <map>
#include <cstdio>
#include <unordered_set>

#include <TF1.h>
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

std::map<int,int> detMapRing = {
  {26, 0}, {30, 1}, {34, 2}, {38, 3}, {25, 4}, {29, 5}, {33, 6}, {37, 7},
  {27, 8}, {31, 9}, {35,10}, {39,11}, {24,12}, {28,13}, {32,14}, {36,15},
  {47,16}, {63,17}, {71,18}, {79,19}, {51,20}, {59,21}, {67,22}, {83,23},
  {50,24}, {58,25}, {66,26}, {82,27}, {44,28}, {60,29}, {68,30}, {76,31},
  {46,32}, {62,33}, {70,34}, {78,35}, {48,36}, {56,37}, {64,38}, {80,39},
  {49,40}, {57,41}, {65,42}, {81,43}, {45,44}, {61,45}, {69,46}, {77,47}
};

std::map<int,int> holeMap = {
  {5,0}, {6,1}, {7,2}, {8,3}, 
  {10,4}, {11,5}, {13,6}, {14,7}, 
  {15,8}, {16,9}, {18,10}, {19,11}
};

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

double thetaCM(double theta,double beta){
  double cosT = TMath::Cos(theta);
  return TMath::ACos((cosT - beta)/(1 - beta*cosT));
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

bool checkEnergyTheta(TF1* fdop, TF1* fres, double energy, double theta){
  return fdop->Eval(theta) + fres->Eval(theta) > energy && fdop->Eval(theta) - fres->Eval(theta) < energy;
}

double dopplerbroad(double *x, double *p){
  double sum = pow(p[1]*sin(x[0])/(1-p[1]*cos(x[0]))*p[2],2) + 
                pow( (-p[1]+cos(x[0]))/((1-p[1]*p[1])*(1-p[1]*cos(x[0])))*p[3] ,2) +
                pow(p[4]*cos(x[0]),2);

  return p[0]*p[5]*TMath::Sqrt(sum);
}

void comptonSort(const TGretinaHit &ghit, int &FP, int &SP) {
  double FOM = 1e10;
  FP = 0; SP = 1;
  int N = ghit.NumberOfInteractions();
  double E = ghit.GetCoreEnergy();
  for (int fp=0; fp < N; fp++){
    double E1 = ghit.GetSegmentEng(fp);
    double er = 511.0/E * E1/(E - E1);
    for (int sp=0; sp < N; sp++){
      if (fp == sp) continue;
      double cosp = TMath::Cos(ghit.GetScatterAngle(fp,sp));
      // double E2 = ghit.GetSegmentEng(sp);
      // double x = er + cosp;
      // double kn = pow((E - E1)/E,2)*((E-E1)/E + E/(E-E1) - pow(TMath::Sin(ghit.GetScatterAngle(fp,sp)),2) );
      // double ffom = std::pow(std::abs(1-x),2.0/3)*ghit.GetAlpha(fp,sp)*kn*std::pow(E/E1,2); // std::pow(E/E2,1.0/3)
      // if (N==3) ffom *= std::pow(E2/E,1.0/3);
      // if (N > 3) ffom *= std::pow(E/E2,1.0/3)*ghit.GetLocalPosition(fp).Z()/(TMath::Sqrt(abs(E1-125)*abs(E2-125))/E);
      double ffom = std::pow(1-er - cosp,2);
      if (ffom < FOM) {
        FOM = ffom;
        FP = fp;
        SP = sp;
      }
    }
  }
  return;
}

double lastPointPenalty(double x){
  double val = TMath::TanH((x/100-1.1)/0.5)*exp(-3*(x/100-1.1)+1);
  if (val < 0) val = 0;
  return val + 1;
}

bool comptonAllowed(const TGretinaHit &ghit){
  int N = ghit.NumberOfInteractions();
  double E = ghit.GetCoreEnergy();
  bool allowed = true;
  for (int i=0; i < N-1; i++){
    if (ghit.GetSegmentEng(i) > 2.0*E/(2+511.0/E)){
      allowed = false;
      break;
    }
    E -= ghit.GetSegmentEng(i);
  }
  return allowed;
}

double comptonChi2(const TGretinaHit &ghit, int fp, int sp){
  static int ncalls = 0;
  static double E = 0;
  if (ncalls == 0) {
    E = ghit.GetCoreEnergy();
  }
  double Edep = ghit.GetSegmentEng(fp);
  double cosE = 1 - 511.0/E * Edep/(E - Edep);
  double cosA = TMath::Cos(ghit.GetScatterAngle(fp,sp));
  double chi2 = pow(cosE-cosA,2);

  double Edep2 = ghit.GetSegmentEng(sp);
  double cosE2 = 1 - 511.0/E * Edep/(E - Edep);
  double cosA2 = TMath::Cos(ghit.GetScatterAngle(sp,fp));

  if (chi2 > pow(cosE2-cosA2,2)) {
    chi2 = pow(cosE2-cosA2,2);
    E -= ghit.GetSegmentEng(sp);
  }
  else 
    E -= ghit.GetSegmentEng(fp);
  
  ncalls++;

  if (ncalls == ghit.NumberOfInteractions()-1) ncalls = 0;
  return chi2;
}

void comptonSortReal(const TGretinaHit &ghit, int &FP, int &SP) {
  double FOM = 1e10;
  int N = ghit.NumberOfInteractions();

  for (int fp=0; fp < N; fp++){
    double E = ghit.GetCoreEnergy();
    double E1 = ghit.GetSegmentEng(fp);
    double er = 511.0/E * E1/(E - E1);
    for (int sp=0; sp < N; sp++){
      if (fp == sp) continue;
      double cosp = TMath::Cos(ghit.GetScatterAngle(fp,sp));
      double E2 = ghit.GetSegmentEng(sp);
      double x = er + cosp;
      double kn = pow((E - E1)/E,2)*((E-E1)/E + E/(E-E1) - pow(TMath::Sin(ghit.GetScatterAngle(fp,sp)),2) );
      double ffom = std::pow(std::abs(1-x),2.0/3)*ghit.GetAlpha(fp,sp)*kn*std::pow(E/E1,3); //*std::pow(E/E2,1.0/3);
      if (N > 2) ffom *= ghit.GetLocalPosition(fp).Z();
      ffom *= std::pow(E/E2,0.5);
      // ffom /= TMath::Sqrt(abs(E1-125)*abs(E2-125))/E;
      ffom *= lastPointPenalty(E1)*lastPointPenalty(E2);
      // double ffom = std::pow(cosp - (1-er),2);

      if (ffom < FOM) {
        FOM = ffom;
        FP = fp;
        SP = sp;
      }
    }
  }
  return;
}

void comptonSortTest2(const TGretinaHit &ghit, int &FP, int &SP) {
  double FOM = 1e10;
  FP = 0; SP = 1;
  int N = ghit.NumberOfInteractions();
  double E = ghit.GetCoreEnergy();

  for (int fp=0; fp < N; fp++){
    double E1 = ghit.GetSegmentEng(fp);
    double er = 511.0/E * E1/(E - E1);

    for (int sp=0; sp < N; sp++){
      if (fp == sp) continue;
      
      double cosp = TMath::Cos(ghit.GetScatterAngle(fp,sp));
      double E2 = ghit.GetSegmentEng(sp);
      double x = er + cosp;
      double kn = pow((E - E1)/E,2)*((E-E1)/E + E/(E-E1) - pow(TMath::Sin(ghit.GetScatterAngle(fp,sp)),2) );
      // double ffom = std::pow(std::abs(1-x),2.0/3)*ghit.GetAlpha(fp,sp)*kn*std::pow(E/E1,2) + std::pow(E1/E - 1/(1+511/E/(1-cosp)),2)*ghit.GetLocalPosition(fp).Z();

      // double ffom = std::pow(std::abs(1-x),2.0/3);
      // double ffom = std::pow(E1/E - 1/(1+511/E/(1-cosp)),2);
      // double ffom = std::pow(std::abs(1-x),2.0/3) + std::pow(E1/E - 1/(1+511/E/(1-cosp)),2);
      double ffom = std::pow(1-x,2);
      // double ffom = std::pow(std::abs(E1/E - 1/(1+511/E/(1-cosp))),2.0/3);
      ffom *= lastPointPenalty(E1)*lastPointPenalty(E2)*std::pow(E/E1,3)*ghit.GetLocalPosition(fp).Z()*ghit.GetAlpha(fp,sp)*kn;
      ffom *= std::pow(E/E2,2) * ghit.GetLocalPosition(sp).Z();
      // if (N > 2) ffom *= ghit.GetLocalPosition(fp).Z();

      // if (N > 3) ffom += E/E2;
      if (ffom < FOM) {
        FOM = ffom;
        FP = fp;
        SP = sp;
      }
    }
  }
  return;
}

void comptonSortTest(const TGretinaHit &ghit, int &FP, int &SP) {
  double FOM = 1e10;
  FP = 0; SP = 1;
  int N = ghit.NumberOfInteractions();
  double E = ghit.GetCoreEnergy();
  
  //find the last point
  int lastpoint = -1;
  if (N > 2){
    int lp1 = -1;
    int lp2 = -1;
    double hypExp = 0;
    for (int i=0; i < N; i++){
      double x = ghit.GetSegmentEng(i)/100;
      double temphypExp = TMath::TanH((x-3.9809*std::pow(N,-1.08))/0.5)*exp(-4*(x-3.9809*std::pow(N,-1.08)) + 2);
      if (temphypExp > hypExp){
        lp2 = lp1;
        lp1 = i;
        hypExp = temphypExp;
      }
    }
    if (lp2 == -1) lastpoint = lp1;
    else if (std::abs(ghit.GetSegmentEng(lp1) - ghit.GetSegmentEng(lp2)) < 80) {
      if (ghit.GetLocalPosition(lp1).Z() > ghit.GetLocalPosition(lp2).Z()) lastpoint = lp1;
      else lastpoint = lp2;
    }
    else lastpoint = lp1;
  }

  //loop through data and exclude the last point
  for (int fp=0; fp < N; fp++){
    if (fp == lastpoint) continue;

    double E1 = ghit.GetSegmentEng(fp);
    double er = 511.0/E * E1/(E - E1);

    for (int sp=0; sp < N; sp++){
      if (fp == sp || sp == lastpoint) continue;
      
      double cosp = TMath::Cos(ghit.GetScatterAngle(fp,sp));
      double E2 = ghit.GetSegmentEng(sp);
      double x = er + cosp;
      double kn = pow((E - E1)/E,2)*((E-E1)/E + E/(E-E1) - pow(TMath::Sin(ghit.GetScatterAngle(fp,sp)),2) );
      double ffom = std::pow(std::abs(1-x),2.0/3)*ghit.GetAlpha(fp,sp)*kn*std::pow(E/E1*ghit.GetLocalPosition(fp).Z()/100,2);
      // double ipFactor = 308.76*std::pow(N,-0.542);
      // if (N > 2) ffom /= TMath::Sqrt(abs(E1-ipFactor)*abs(E2-ipFactor))/E;
      // if (N-1 == 3) ffom *= 1.0/(E/E2*ghit.GetLocalPosition(sp).Z()/100);
      if (N-1 > 2) ffom *= std::pow(E/E2*ghit.GetLocalPosition(sp).Z()/100,1.0/3);

      if (ffom < FOM) {
        FOM = ffom;
        FP = fp;
        SP = sp;
      }
    }
  }
  return;
}

void comptonSortGate(const TGretinaHit &ghit, std::vector<int> mypassed, int &FP, int &SP) {
  double FOM = 1e10;
  int N = ghit.NumberOfInteractions();
  for (int fp=0; fp < (int) mypassed.size(); fp++){
    double E = ghit.GetCoreEnergy();
    double E1 = ghit.GetSegmentEng(mypassed[fp]);
    double er = 511.0/E * E1/(E - E1);
    for (int sp=0; sp < N; sp++){
      if (mypassed[fp] == sp) continue;
      double cosp = TMath::Cos(ghit.GetScatterAngle(mypassed[fp],sp));
      double x = er + cosp;
      double y = ghit.GetAlpha(mypassed[fp],sp);
      // double kn = pow((E - E1)/E,2)*((E-E1)/E + E/(E-E1) - pow(TMath::Sin(ghit.GetScatterAngle(mypassed[fp],sp)),2) );
      double kn = pow((E - E1)/E1,2)*((E-E1)/E + E/(E-E1) - pow(TMath::Sin(ghit.GetScatterAngle(mypassed[fp],sp)),2) );
      //NOT THIS ONE
      double ffom = std::pow(std::abs(1-x),2.0/3)*y*kn /*std::pow(E/E1,2)*kn*/;
      if (ffom < FOM) {
        FOM = ffom;
        FP = mypassed[fp];
        SP = sp;
      }
    }
  }
  return;
}

bool segmentationCheck(const TGretinaHit &ghit) {
  bool segmentation = false;
  // int nInt = ghit.NumberOfInteractions();
  double zSeg[6] = {8,22,38,56,76,90};

  for (int i=0 ; i < 2; i++){
    double x = ghit.GetLocalPosition(i).X();
    double y = ghit.GetLocalPosition(i).Y();

    //phi segmentation check
    bool iseg = false;
    double width = 1.0;
    int crystalType = (1-ghit.GetCrystalId()%2); //a or b type
    for (int j=-1; j < 2; j++) {
      double angle = j*TMath::Pi()/3 + crystalType*TMath::Pi()/6;
      iseg = iseg || std::abs(y - std::tan(angle)*x) < width/2/std::cos(angle);
    }

    //Z segmentation check
    double z = ghit.GetLocalPosition(i).Z();
    for (int j=0; j < 6; j++) iseg = iseg || abs(zSeg[j]-z) < width;

    if (iseg) {
      segmentation = true;
      break;
    }
  }

  return segmentation;
}

double calcE1Compt(const TGretinaHit &ghit, int fp, int sp){
  double E = ghit.GetCoreEnergy();
  double E1 = ghit.GetSegmentEng(fp);
  return E1/E - 1./(1 + 511/E/(1-TMath::Cos(ghit.GetScatterAngle(fp,sp))));
}

double calcEnergyRatio(double E, double Edep){
  return 511.0/E * Edep/(E - Edep);
}

double calcEnergyCos(double E, double Edep){
  return 1 - 511.0/E * Edep/(E - Edep);
}

double calcKN(const TGretinaHit &ghit, int fp, int sp){
  double E = ghit.GetCoreEnergy();
  double E1 = ghit.GetSegmentEng(fp);
  return pow((E - E1)/E,2)*((E-E1)/E + E/(E-E1) - pow(TMath::Sin(ghit.GetScatterAngle(fp,sp)),2) );  // *TMath::Sin(ghit.GetScatterAngle(fp,sp));
}

std::map<int,int> buildCorrectMap(const TGretinaHit &htime,const TGretinaHit &hmain) {
  int N = htime.Size();
  std::map<int,int> outmap;
  for (int i=0; i < N; i++){
    for (int j=0; j < N; j++){
      if (htime.GetSegmentEng(i) == hmain.GetSegmentEng(j) && htime.GetTheta(i) == hmain.GetTheta(j) && htime.GetPhi(i) == hmain.GetPhi(j)) {
        outmap.insert(std::pair<int,int> (j,i));
        break;
      }
    }
  }

  return outmap;
}

bool gates_loaded = false;
std::map<std::string,std::vector<GCutG*>> gates;

int EVTnumber = -1;

// extern "C" is needed to prevent name mangling.
// The function signature must be exactly as shown here,
//   or else bad things will happen.
extern "C"
void MakeHistograms(TRuntimeObjects& obj) {
  EVTnumber++;
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

  if (!gretina || gretsim->Size() == 0){
    return;
  }
  
  if (!s800sim || !s800sim->Size() || GValue::Value("BETA") == 0.0){
    stopped = true;
  }

  // int nhits = gretsim->Size();
  // obj.FillHistogram("GEANT", "GEANT_gamma_hits",20,0,20,nhits);
  // for (int i=0; i < nhits; i++){
  //   TGretSimHit simHit = gretsim->GetGretinaSimHit(i);
  //   double simTheta = simHit.GetTheta();
  //   obj.FillHistogram("GEANT", "simtheta_dist", 360,0,180,simTheta*TMath::RadToDeg());
  //   obj.FillHistogram("GEANT", "simtheta_cm_dist", 360,0,180,thetaCM(simTheta,simHit.GetBeta())*TMath::RadToDeg());
  //   obj.FillHistogram("GEANT", "simBeta", 100,0.3,0.5,simHit.GetBeta());
  //   if (simTheta*TMath::RadToDeg() > 40 && simTheta*TMath::RadToDeg() <120)
  //     obj.FillHistogram("GEANT", "simtheta_cm_limited_dist", 360,0,180,thetaCM(simTheta,simHit.GetBeta())*TMath::RadToDeg());
  // }
  
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

  TGretSimHit simHit = gretsim->GetGretinaSimHit(0);
  double simTheta = simHit.GetTheta();
  obj.FillHistogram(dirname,"GEANT_energies",1000,500,1500,simHit.GetDoppler(0,&track));
  // double gammaEn = simHit.GetEn();
  // obj.FillHistogram(dirname,"GEANT_gretsim_size",20,0,20,gretsim->Size());
  // obj.FillHistogram(dirname,"GEANT_s800sim_size",20,0,20,s800sim->Size());
  // bool isFEP = simHit.IsFEP();
  // double gammaEnDop = simHit.GetDoppler();
  // double gammaBeta = simHit.GetBeta();
  /*
  //ADDBACK
  int abSize = gretina->NNAddbackSize(0);
  for (int i=0; i < abSize; i++){
    TGretinaHit nnhit = gretina->GetNNAddbackHit(0,i);
    if (nnhit.NumberOfInteractions() > 1){
      double r0 = 0;
      for (int ip=0; ip < nnhit.NumberOfInteractions(); ip++){
        double xx = nnhit.GetLocalPosition(ip).X();
        double yy = nnhit.GetLocalPosition(ip).Y();
        double ipPhi = nnhit.GetPhi(ip);
        double ipTheta = nnhit.GetTheta(ip);

        if (simHit.IsFEP()) {
          obj.FillHistogram(dirname, "n0_local_pos_xy_FEP",200,-50,50,xx,200,-50,50,yy);
          obj.FillHistogram(dirname, "n0_theta_vs_phi_FEP",720,0,TMath::TwoPi(),ipPhi,360,0,TMath::Pi(),ipTheta);
        }
        else {
          obj.FillHistogram(dirname, "n0_local_pos_xy_!FEP",200,-50,50,xx,200,-50,50,yy);
          obj.FillHistogram(dirname, "n0_theta_vs_phi_!FEP",720,0,TMath::TwoPi(),ipPhi,360,0,TMath::Pi(),ipTheta);
        }
        if (ip == 0) r0 = TMath::Sqrt(xx*xx + yy*yy);
        else {
          if (ip == 1 && 19.5 < r0 && r0 < 20.5){
            if (simHit.IsFEP()) {
              obj.FillHistogram(dirname, "n0_local_pos_xy_FEP_r0",200,-50,50,xx,200,-50,50,yy);
              obj.FillHistogram(dirname, "n0_theta_vs_phi_FEP_r0",720,0,TMath::TwoPi(),ipPhi,360,0,TMath::Pi(),ipTheta);
            }
            else { 
              obj.FillHistogram(dirname, "n0_local_pos_xy_!FEP_r0",200,-50,50,xx,200,-50,50,yy);
              obj.FillHistogram(dirname, "n0_theta_vs_phi_!FEP_r0",720,0,TMath::TwoPi(),ipPhi,360,0,TMath::Pi(),ipTheta);
            }
          }
        }
      }
    }
  }
  */
  
  //SINGLES
  int gSize = gretina->Size();
  for (int i=0; i < gSize; i++){
    TGretinaHit &hit = gretina->GetGretinaHit(i);
    TGretinaHit hitCopy;
    hit.Copy(hitCopy);
    hitCopy.SortSegments();
    std::map<int,int> truePoints = buildCorrectMap(hit,hitCopy);
    // EnergySmear(hit,rand_gen);
    EnergySmear(hitCopy,rand_gen);
    
    int nInteractions = hit.NumberOfInteractions();
    // double energy_corrected = rand_gen->Gaus(hit.GetCoreEnergy(),hit.GetCoreEnergy()*0.0027/2.35);
    // if (!stopped) energy_corrected = hit.GetDoppler(simHit.GetBeta(),&track);
    double core_energy = hit.GetCoreEnergy();

    double true_energy = core_energy;
    double time_energy = core_energy;
    double main_energy = core_energy;
    if (!stopped){
      true_energy = simHit.GetDoppler(0.0);
      time_energy = hit.GetDopplerYta(s800sim->AdjustedBeta(GValue::Value("BETA")),yta,&track);
      main_energy = hitCopy.GetDopplerYta(simHit.GetBeta(), yta, &track);
    }
    
    double theta = hit.GetTheta();
    double phi = hit.GetPhi();
    int cryID = hit.GetCrystalId();
    int ringNum = hit.GetRingNumber();

    double xi = 0.0; 
    if (stopped) xi = hit.GetXi();
    else xi = hit.GetXi(&track);
    
    obj.FillHistogram(dirname, "core_energy", 8192,0,8192, core_energy);
    obj.FillHistogram(dirname, "gam_dop_sgl",4096,0,4096, time_energy);
    
    if (simHit.IsFEP()) {
      // obj.FillHistogram(dirname, "FEP_spec_true_energy_vs_theta",1000,0,TMath::Pi(),theta,4000,1100,1500,true_energy);
      // obj.FillHistogram(dirname, "FEP_spec_time_sorted_vs_theta",1000,0,TMath::Pi(),theta,4000,1100,1500,time_energy);
      // obj.FillHistogram(dirname, "FEP_spec_main_sorted_vs_theta",1000,0,TMath::Pi(),theta,4000,1100,1500,main_energy);
      obj.FillHistogram(dirname,"FEP_nInteractions",16,0,16,nInteractions);
      if (nInteractions > 1){
        int myFP = -1;
        int mySP = -1;
        comptonSortTest2(hitCopy,myFP,mySP);
        bool segmentation = segmentationCheck(hit);

        // obj.FillHistogram(dirname,Form("gretina_map_FEP_hn%d_cry%d",holeMap[hit.GetHoleNumber()],cryID),360,0,360,phi*TMath::RadToDeg());
        // obj.FillHistogram(dirname,Form("gretina_map_FEP_hn%d",holeMap[hit.GetHoleNumber()]),360,0,360,phi*TMath::RadToDeg());

        obj.FillHistogram(dirname,"gretina_map_FEP_nInt>1",900,0,360,phi*TMath::RadToDeg(),450,0,180,theta*TMath::RadToDeg());
        if (segmentation) 
          obj.FillHistogram(dirname,"gretina_map_FEP_nInt>1_seg",900,0,360,phi*TMath::RadToDeg(),450,0,180,theta*TMath::RadToDeg());
        else 
          obj.FillHistogram(dirname,"gretina_map_FEP_nInt>1_!seg",900,0,360,phi*TMath::RadToDeg(),450,0,180,theta*TMath::RadToDeg());
        // double eSum = 0;
        // for (int nip=0; nip < nInteractions; nip++) eSum += hit.GetSegmentEng(nip);
        // obj.FillHistogram(dirname,"Esum_vs_CoreEnergy",1000,1000,2000,eSum,1000,1000,2000,core_energy);

        // compton and klein nishina quantities
        // > 1 means that the scattering angle cos is > the energy cos
        double trueCompt = calcEnergyRatio(core_energy,hit.GetSegmentEng(0)) + TMath::Cos(hit.GetScatterAngle());
        double mainCompt = calcEnergyRatio(core_energy,hitCopy.GetSegmentEng(0)) + TMath::Cos(hitCopy.GetScatterAngle());
        double algoCompt = calcEnergyRatio(core_energy,hitCopy.GetSegmentEng(myFP)) + TMath::Cos(hitCopy.GetScatterAngle(myFP,mySP));
        double trueE1Compt = calcE1Compt(hit,0,1);
        double trueKn = calcKN(hit,0,1);
        double mainKn = calcKN(hitCopy,0,1);
        double algoKn = calcKN(hitCopy,myFP,mySP);
        double truePosCos = TMath::Cos(hit.GetScatterAngle(0,1));
        double trueEngCos = calcEnergyCos(core_energy,hit.GetSegmentEng(0));

        bool allowedByCompton = hit.GetSegmentEng(0)/core_energy > 1.0/(1+511.0/998/(1-truePosCos)) && hit.GetSegmentEng(0)/core_energy < 1.0/(1+511.0/1834/(1-truePosCos));

        //Guess Matrices
        int algoGoodFP = truePoints[myFP] == 0;
        int algoGoodSP = truePoints[mySP] == 1;
        int mainGoodFP = truePoints[0] == 0;
        int mainGoodSP = truePoints[1] == 1;
        obj.FillHistogram(dirname, "guess_matrix_algo",2,0,2,algoGoodFP,2,0,2,algoGoodSP);
        obj.FillHistogram(dirname, "guess_matrix_main",2,0,2,mainGoodFP,2,0,2,mainGoodSP);

        // bool allowedByCompton = trueEngCos > truePosCos - 0.02 && trueEngCos < truePosCos + 0.02;
        // bool allowedByCompton = trueEngCos < truePosCos - 0.1;
        // bool testRegion = trueEngCos < -0.4 && trueEngCos > -0.5 && truePosCos > -0.9 && truePosCos < -0.85;

        // if (hit.GetSegmentEng(0)/core_energy > 0.2 && hit.GetSegmentEng(0)/core_energy < 0.3 && truePosCos > -0.6 && truePosCos < -0.4){
        // // if (allowedByCompton && testRegion) {
        //   std::cout<<"Evt # "<<EVTnumber<<" E: "<<simHit.GetEn()<<std::endl;
        // } 
        
        if (nInteractions == 2 && allowedByCompton) {
          obj.FillHistogram(dirname, "Algo_Chi2",1000,0,5,std::pow(calcEnergyCos(core_energy,hit.GetSegmentEng(myFP)) - TMath::Cos(hitCopy.GetScatterAngle(myFP,mySP)),2));
          double theChi2 = std::pow(calcEnergyCos(core_energy,hit.GetSegmentEng(0)) - TMath::Cos(hit.GetScatterAngle(0,1)),2);
          double swapChi2 = std::pow(calcEnergyCos(core_energy,hit.GetSegmentEng(1)) - TMath::Cos(hit.GetScatterAngle(1,0)),2);
          bool doswap = theChi2 > swapChi2;
          // bool bothRight = truePoints[myFP] == 0 && truePoints[mySP] == 1;

          // if (bothRight && doswap) {
          //   obj.FillHistogram(dirname,"AlgoChi2_alpha_vs_compt_FEP",400,-1,3,trueCompt,200,0,20,hit.GetAlpha()*TMath::RadToDeg());
          //   obj.FillHistogram(dirname,"AlgoChi2_kn_vs_compt_FEP",400,-1,3,trueCompt,200,0,2,trueKn);
          //   obj.FillHistogram(dirname,"AlgoChi2_engCos_vs_scatCos_FEP",200,-1,1,truePosCos,300,-5,1,trueEngCos);
          // }
          
          int theFP = 1;
          int theSP = 1;
          if (doswap) {theFP = 0; theSP = 0;}
          obj.FillHistogram(dirname, "guess_matrix_2int_min_chi2",2,0,2,theFP,2,0,2,theSP);
          obj.FillHistogram(dirname, "guess_matrix_2int_algo",2,0,2,algoGoodFP,2,0,2,algoGoodSP);
        }

        // alpha vs compt
        obj.FillHistogram(dirname,"surf_TRUE_alpha_vs_compt_FEP",400,-1,3,trueCompt,200,0,20,hit.GetAlpha()*TMath::RadToDeg());
        obj.FillHistogram(dirname,"surf_MAIN_alpha_vs_compt_FEP",400,-1,3,mainCompt,200,0,20,hitCopy.GetAlpha()*TMath::RadToDeg());
        obj.FillHistogram(dirname,"surf_ALGO_alpha_vs_compt_FEP",400,-1,3,algoCompt,200,0,20,hitCopy.GetAlpha(myFP,mySP)*TMath::RadToDeg());
        // kn vs compt
        obj.FillHistogram(dirname,"surf_TRUE_kn_vs_compt_FEP",400,-1,3,trueCompt,200,0,2,trueKn);
        obj.FillHistogram(dirname,"surf_MAIN_kn_vs_compt_FEP",400,-1,3,mainCompt,200,0,2,mainKn);
        obj.FillHistogram(dirname,"surf_ALGO_kn_vs_compt_FEP",400,-1,3,algoCompt,200,0,2,algoKn);
        // e1Compt
        obj.FillHistogram(dirname,"surf_TRUE_E1compt_vs_compt_FEP",400,-1,3,trueCompt,400,-2,2,trueE1Compt);
        obj.FillHistogram(dirname,"surf_TRUE_kn_vs_E1compt_FEP",400,-1,1,trueE1Compt,200,0,2,trueKn);
        obj.FillHistogram(dirname,"surf_TRUE_alpha_vs_E1compt_FEP",400,-1,1,trueE1Compt,200,0,20,hit.GetAlpha()*TMath::RadToDeg());
        //eratio vs scatCos
        obj.FillHistogram(dirname,"surf_TRUE_engCos_vs_scatCos_FEP",400,-1,1,truePosCos,400,-1,1,trueEngCos);
        obj.FillHistogram(dirname,"surf_TRUE_E1_vs_engCos_FEP",300,-5,1,trueEngCos,200,0,1,hit.GetSegmentEng(0)/core_energy);
        obj.FillHistogram(dirname,"surf_TRUE_E1_vs_scatterCos_FEP",400,-1,1,truePosCos,400,0,1,hit.GetSegmentEng(0)/core_energy);
        //theta & phi vs compt
        obj.FillHistogram(dirname,"surf_TRUE_theta_vs_compt_FEP",400,-1,3,trueCompt,360,0,180,theta*TMath::RadToDeg());
        obj.FillHistogram(dirname,"surf_TRUE_phi_vs_compt_FEP",400,-1,3,trueCompt,720,0,360,phi*TMath::RadToDeg());

        //e1 & e2 vs compt
        obj.FillHistogram(dirname,"surf_TRUE_E1_vs_compt_FEP",400,-1,3,trueCompt,500,0,1,hit.GetSegmentEng(0)/core_energy);
        obj.FillHistogram(dirname,"surf_TRUE_E2_vs_compt_FEP",400,-1,3,trueCompt,500,0,1,hit.GetSegmentEng(1)/core_energy);

        //INTERACTION POINT ANALYSIS
        if (nInteractions < 5){
          obj.FillHistogram(dirname, Form("guess_matrix_algo_%dint",nInteractions),2,0,2,algoGoodFP,2,0,2,algoGoodSP);
          obj.FillHistogram(dirname, Form("guess_matrix_main_%dint",nInteractions),2,0,2,mainGoodFP,2,0,2,mainGoodSP);

          //matrix breakdown
          obj.FillHistogram(dirname, Form("brkdn_%dint_F%dS%d_TRUE_kn_vs_compt",nInteractions,algoGoodFP,algoGoodSP),400,-1,3,trueCompt,200,0,2,trueKn);
          obj.FillHistogram(dirname, Form("brkdn_%dint_F%dS%d_TRUE_alpha_vs_compt",nInteractions,algoGoodFP,algoGoodSP),400,-1,3,trueCompt,200,0,20,hit.GetAlpha()*TMath::RadToDeg());
          if (!(algoGoodFP == 1 && algoGoodSP == 1)){
            obj.FillHistogram(dirname, Form("brkdn_%dint_F%dS%d_ALGO_kn_vs_compt",nInteractions,algoGoodFP,algoGoodSP),400,-1,3,algoCompt,200,0,2,algoKn);
            obj.FillHistogram(dirname, Form("brkdn_%dint_F%dS%d_ALGO_alpha_vs_compt",nInteractions,algoGoodFP,algoGoodSP),400,-1,3,algoCompt,200,0,20,hitCopy.GetAlpha(myFP,mySP)*TMath::RadToDeg());
          }

          for (int pp=0; pp < nInteractions; pp++) {
            obj.FillHistogram(dirname, Form("IP_energy_vs_idx_%dint",nInteractions), 7,0,7, pp, 1800,0,1800,hit.GetSegmentEng(pp)); 
            obj.FillHistogram(dirname, Form("IP_esort_idx_vs_idx_%dint",nInteractions), 7,0,7, pp, 7,0,7, truePoints[pp]); 
            obj.FillHistogram(dirname, Form("IP_localz_vs_idx_%dint",nInteractions), 7,0,7, pp, 100,0,100,hit.GetLocalPosition(pp).Z());
            if (gates.count("e1_cos") && gates["e1_cos"][0]->IsInside(truePosCos,hit.GetSegmentEng(0)/core_energy)) {
              obj.FillHistogram(dirname, Form("GATED_IP_energy_vs_idx_%dint",nInteractions), 7,0,7, pp, 1800,0,1800,hit.GetSegmentEng(pp));
              obj.FillHistogram(dirname, Form("GATED_IP_localz_vs_idx_%dint",nInteractions), 7,0,7, pp, 100,0,100,hit.GetLocalPosition(pp).Z());
            }
          }
          if (gates.count("e1_cos") && gates["e1_cos"][0]->IsInside(truePosCos,hit.GetSegmentEng(0)/core_energy)) {
            obj.FillHistogram(dirname, "GATED_gretina_map",900,0,360,phi*TMath::RadToDeg(),450,0,180,theta*TMath::RadToDeg());
          }
          //region gates
          // if (gates["algo_kn_compt"][0] && gates["algo_kn_compt"][0]->IsInside(trueCompt,trueKn)) {
          //   obj.FillHistogram(dirname,"GATED_surf_alpha_vs_compt_true_FEP",400,-1,3,trueCompt,200,0,20,hit.GetAlpha()*TMath::RadToDeg());
          //   for (int pp=0; pp < nInteractions; pp++){
          //     obj.FillHistogram(dirname, Form("GATED_energy_vs_idx_%dint",nInteractions), 7,0,7, pp, 1800,0,1800,hit.GetSegmentEng(pp));
          //     obj.FillHistogram(dirname, Form("GATED_localz_vs_idx_%dint",nInteractions), 7,0,7, pp, 100,0,100,hit.GetLocalPosition(pp).Z());
          //     obj.FillHistogram(dirname, Form("GATED_esort_idx_vs_idx_%dint",nInteractions), 7,0,7, pp, 7,0,7, truePoints[pp]);
          //   }
          // }
        }

        //energy spectra + vs xi
        double algo_energy = core_energy;
        if (!stopped) algo_energy = hitCopy.GetDopplerYta(s800sim->AdjustedBeta(GValue::Value("BETA")),yta,&track,myFP);
        double algo_xi = hitCopy.GetXi(&track,myFP,mySP);
        // obj.FillHistogram(dirname, "gam_dop_algo_fep",1500,0,1500, algo_energy);
        // obj.FillHistogram(dirname, "gam_dop_time_fep",1500,0,1500, time_energy);
        // obj.FillHistogram(dirname, "gam_dop_main_fep",1500,0,1500, main_energy);
        obj.FillHistogram(dirname, "gam_dop_algo_fep_vs_xi",360,0,TMath::TwoPi(),algo_xi,1500,0,1500,algo_energy);
        obj.FillHistogram(dirname, "True_xi_vs_Algo_xi",360,0,TMath::TwoPi(),algo_xi,360,0,TMath::TwoPi(),xi);
        obj.FillHistogram(dirname, "algo_xi_vs_hole_number",12,0,12,holeMap[hit.GetHoleNumber()],180,0,TMath::TwoPi(),algo_xi);
        obj.FillHistogram(dirname, "algo_xi_vs_crystalID",48,0,48,detMapRing[cryID],90,0,180,algo_xi*TMath::RadToDeg());
        // if (segmentation)
        //   obj.FillHistogram(dirname, "True_xi_vs_hole_number_seg",12,0,12,holeMap[hit.GetHoleNumber()],180,0,TMath::TwoPi(),algo_xi);
        // else 
        //   obj.FillHistogram(dirname, "True_xi_vs_hole_number_!seg",12,0,12,holeMap[hit.GetHoleNumber()],180,0,TMath::TwoPi(),algo_xi);
        // if (hit.GetHoleNumber() < 10 || (hit.GetHoleNumber() > 13 && hit.GetHoleNumber() < 17)){
        //   obj.FillHistogram(dirname, "True_xi_vs_Algo_xi_e12501_quads",360,0,TMath::TwoPi(),algo_xi,360,0,TMath::TwoPi(),xi);
        // }
      }
    }

    if (nInteractions > 1){
      // obj.FillHistogram(dirname, "core_energy_vs_xi",360,0,TMath::TwoPi(),xi, 4096,0,4096, core_energy);
      obj.FillHistogram(dirname, "nopolgate_Edop_vs_xi",360,0,TMath::TwoPi(),xi,1500,0,1500,time_energy);
      if (simHit.IsFEP()) {
        obj.FillHistogram(dirname, "nopolgate_Edop_vs_xi_FEP",360,0,TMath::TwoPi(),xi,4000,1100,1500,time_energy);
        // if (nInteractions < 10) 
        //   obj.FillHistogram(dirname, Form("nopolgate_Edop_vs_xi_FEP_nint%d",nInteractions),360,0,TMath::TwoPi(),xi,4000,1100,1500,time_energy);
      }
    }

    // //IP core energy theta
    // for (int ip=0; ip < nInteractions; ip++) {
    //   double ip_theta = hit.GetTheta(ip);
    //   obj.FillHistogram(dirname, "IP_Ecore_vs_theta",120,0.6,2.1,ip_theta,1850,0,1850,core_energy);
    //   // obj.FillHistogram(dirname, "IP_Ecore_vs_theta_fine",600,0.6,2.1,ip_theta,4000,1000,1800,core_energy);
    //   // if (ip == 0 && simHit.IsFEP()) obj.FillHistogram(dirname, "IP_Ecore_vs_theta_fine_fp_fep",600,0.6,2.1,ip_theta,4000,1000,1800,core_energy);
    //   // if (simHit.IsFEP()) 
    //   //   obj.FillHistogram(dirname, "IP_Ecore_FEP_vs_theta",120,0.6,2.1,ip_theta,1850,0,1850,core_energy);
    //   // else 
    //   //   obj.FillHistogram(dirname, "IP_Ecore_!FEP_vs_theta",120,0.6,2.1,ip_theta,1850,0,1850,core_energy);
    // }
    
    /*
    if (nInteractions > 1) {
      TF1 *fDOP = new TF1("fdop","[0]*TMath::Sqrt(1-[1]*[1])/(1 - [1]*TMath::Cos(x))",0,TMath::Pi());
      TF1 *fRES = new TF1("fresolution",dopplerbroad,0,TMath::Pi(),6);
      fDOP->SetParameter(1,0.4266);
      fRES->SetParameter(1,0.4266);
      fRES->SetParameter(2,GValue::Value("DOPP_BROAD_DTHETA"));
      fRES->SetParameter(3,GValue::Value("DOPP_BROAD_DBETA"));
      fRES->SetParameter(4,GValue::Value("DOPP_BROAD_BEAMSPOT"));
      fRES->SetParameter(5,GValue::Value("DOPP_BROAD_SCALE"));

      double H_elo = 1309;
      double H_ehi = 1329;
      int step = 20;
      int H_nbins = int(H_ehi-H_elo)/step;
      double CentroidRestEnergy = H_elo;
      
      std::vector<int> passedIP;
      
      while (CentroidRestEnergy < H_ehi){
        CentroidRestEnergy += step;
        fDOP->SetParameter(0,CentroidRestEnergy);
        fRES->SetParameter(0,CentroidRestEnergy);
        
        for (int ip=0; ip < nInteractions; ip++){
          if (checkEnergyTheta(fDOP,fRES,core_energy,hitCopy.GetTheta(ip))){
            passedIP.push_back(ip);
          }
        }

        // if (CentroidRestEnergy == 1329.0 && simHit.IsFEP()){
        //   double energyThresh = 0;
        //   while (energyThresh < 500) {
        //     bool passedGate = false;
        //     for (int nps =0 ; nps < passedIP.size(); nps++){
        //       if (hit.GetSegmentEng(passedIP[nps]) > energyThresh) {
        //         passedGate = true;
        //         break;
        //       }
        //     }
        //     if (passedGate) {
        //       obj.FillHistogram(dirname, Form("RE1329_pass_xi_thresh%f",energyThresh),90,0,TMath::TwoPi(),xi);
        //       obj.FillHistogram(dirname, Form("RE1329_pass_coreE_vs_theta_thresh%f",energyThresh),600,0.6,2.1,hit.GetTheta(0),4000,1000,1800,core_energy);
        //     }
        //     else {
        //       obj.FillHistogram(dirname, Form("RE1329_fail_xi_thresh%f",energyThresh),90,0,TMath::TwoPi(),xi);
        //       obj.FillHistogram(dirname, Form("RE1329_fail_coreE_vs_theta_thresh%f",energyThresh),600,0.6,2.1,hit.GetTheta(0),4000,1000,1800,core_energy);
        //     }
        //     energyThresh += 50;
        //   }
        // }

        if (passedIP.empty()) 
          obj.FillHistogram(dirname, "CentroidEnergy_vs_DopReconNotPassed",H_nbins,H_elo,H_ehi,main_energy,H_nbins,H_elo-0.5,H_ehi-0.5,CentroidRestEnergy);
        
        else {
          double E_pass_dop = hitCopy.GetDopplerYta(simHit.GetBeta(), 0.0, &track, passedIP[0]);
          obj.FillHistogram(dirname, "CentroidEnergy_vs_DopRecon",H_nbins,H_elo,H_ehi,E_pass_dop,H_nbins,H_elo-0.5,H_ehi-0.5,CentroidRestEnergy);

          //NEW SECOND INTERACTION POINT
          double xiOrder, xiMain, xiCompt = hitCopy.GetXi(&track);
          int spOrder, spMain, spCompt = 1;
          if (passedIP[0] != 0) {
            xiMain = hitCopy.GetXi(&track,passedIP[0],0);
            spMain = 0;

            xiOrder = hitCopy.GetXi(&track,passedIP[0],passedIP[0]+1);
            spOrder = passedIP[0]+1;
          }

          int myFP = -1;
          int mySP = -1;
          comptonSortGate(hitCopy,passedIP,myFP,mySP);
          // double FOM = 0;
          // for (int fp=0; fp < passedIP.size(); fp++){
          //   double er = 511.0/hitCopy.GetCoreEnergy() * hitCopy.GetSegmentEng(passedIP[fp])/(hitCopy.GetCoreEnergy() - hitCopy.GetSegmentEng(passedIP[fp]));
          //   for (int sp=0; sp < nInteractions; sp++){
          //     if (passedIP[fp] == sp) continue;
          //     double cosp = TMath::Cos(hitCopy.GetScatterAngle(passedIP[fp],sp));
          //     double x = er + cosp;
          //     double y = hitCopy.GetAlpha(passedIP[fp],sp);
          //     double ffom = std::abs(1-x)*y;
          //     if (FOM < ffom) {
          //       FOM = ffom;
          //       // best_points = std::make_pair(fp,sp);
          //       myFP = passedIP[fp];
          //       mySP = sp;
          //     }
          //   }
          // }
          xiCompt = hitCopy.GetXi(&track,myFP,mySP);
          
          // int fpCompt = passedIP[0];
          // for (int fpidx = 0; fpidx < (int) passedIP.size(); fpidx++){
          //   double FoM = 0;
          //   double cosE = 1 - 511.0/hitCopy.GetCoreEnergy() * hitCopy.GetSegmentEng(passedIP[fpidx])/(hitCopy.GetCoreEnergy() - hitCopy.GetSegmentEng(passedIP[fpidx]));
          //   if (cosE < -1) break;
          //   for (int sp=0; sp < nInteractions; sp++){
          //     if (sp == passedIP[fpidx]) continue;
          //     double cosP = TMath::Cos(hitCopy.GetScatterAngle(passedIP[fpidx],sp));
          //     if (FoM < std::pow(cosP - cosE,2)){
          //       FoM = std::pow(cosP - cosE,2);
          //       fpCompt = passedIP[fpidx];
          //       spCompt = sp;
          //     }
          //   }
          //   if ((1-cosE) + TMath::Cos(hitCopy.GetScatterAngle(fpCompt,spCompt)) > -0.1) break;
          // }
          // xiCompt = hitCopy.GetXi(&track,fpCompt,spCompt);

          // obj.FillHistogram(dirname, "CentroidEnergy_vs_TrueXi",90,0,TMath::TwoPi(),xi,H_nbins,H_elo-0.5,H_ehi-0.5,CentroidRestEnergy);
          // obj.FillHistogram(dirname, "CentroidEnergy_vs_OrderXi",90,0,TMath::TwoPi(),xiOrder,H_nbins,H_elo-0.5,H_ehi-0.5,CentroidRestEnergy);
          // obj.FillHistogram(dirname, "CentroidEnergy_vs_MainXi",90,0,TMath::TwoPi(),xiMain,H_nbins,H_elo-0.5,H_ehi-0.5,CentroidRestEnergy);
          // obj.FillHistogram(dirname, "CentroidEnergy_vs_ComptXi",90,0,TMath::TwoPi(),xiCompt,H_nbins,H_elo-0.5,H_ehi-0.5,CentroidRestEnergy);

          if (simHit.IsFEP()){
            // obj.FillHistogram(dirname, "TrueXi_vs_OrderXi",90,0,TMath::TwoPi(),xiOrder,90,0,TMath::TwoPi(),xi);
            // obj.FillHistogram(dirname, "TrueXi_vs_MainXi",90,0,TMath::TwoPi(),xiMain,90,0,TMath::TwoPi(),xi);
            obj.FillHistogram(dirname, "TrueXi_vs_ComptXi",360,0,TMath::TwoPi(),xiCompt,360,0,TMath::TwoPi(),xi);

            // int firstPointRight = !(abs(hit.GetTheta() - hitCopy.GetTheta(fpCompt)) > 0.0 || abs(hit.GetPhi() - hitCopy.GetPhi(fpCompt)) > 0.0);
            // int secndPointRight = !(abs(hit.GetTheta(1) - hitCopy.GetTheta(spCompt)) > 0.0 || abs(hit.GetPhi(1) - hitCopy.GetPhi(spCompt)) > 0.0);

            // obj.FillHistogram(dirname,"Guess_mtx_compt_sp_vs_fp",2,0,2,firstPointRight,2,0,2,secndPointRight);
            // obj.FillHistogram(dirname, Form("TrueXi_vs_ComptXi_fp%d_sp%d",firstPointRight,secndPointRight),90,0,TMath::TwoPi(),xiCompt,90,0,TMath::TwoPi(),xi);
          }
          
          // if (simHit.IsFEP()) {
          //   obj.FillHistogram(dirname, "CentroidEnergy_vs_TrueXi_FEP",90,0,TMath::TwoPi(),xi,H_nbins,H_elo-0.5,H_ehi-0.5,CentroidRestEnergy);
          //   obj.FillHistogram(dirname, "CentroidEnergy_vs_CalcXi_FEP",90,0,TMath::TwoPi(),xi,H_nbins,H_elo-0.5,H_ehi-0.5,CentroidRestEnergy);
          // }
          // else {
          //   obj.FillHistogram(dirname, "CentroidEnergy_vs_TrueXi_!FEP",90,0,TMath::TwoPi(),xi,H_nbins,H_elo-0.5,H_ehi-0.5,CentroidRestEnergy);
          // }
        }
        passedIP.clear();
      }
    }*/
    
    /*
    if (nInteractions > 1){
      // ECORE THETA INTERACTION POINT GATES
      obj.FillHistogram(dirname, "IP_Ecore_FEP_vs_theta_TRUEFP",120,0.6,2.1,hitCopy.GetTheta(),4096,0,4096,core_energy);
      for (int ip=0; ip < nInteractions; ip++)
        obj.FillHistogram(dirname, "IP_Ecore_FEP_vs_theta",120,0.6,2.1,hit.GetTheta(ip),4096,0,4096,core_energy);
      
      for (auto ipgate : gates["intpnt"]){
        std::string ipgname = std::string(ipgate->GetName());
        std::vector<int> IPCorepass;
        for (int ip=0; ip < nInteractions; ip++){
          if (ipgate->IsInside(hit.GetTheta(ip),core_energy) && hit.GetSegmentEng(ip) > 100){
            IPCorepass.push_back(ip);
          }
        }

        if (IPCorepass.size() > 0){
          obj.FillHistogram(dirname, Form("IP_%s_npass_vs_totalPoints",ipgname.c_str()),11,1,12,nInteractions,9,1,10,(int) IPCorepass.size());
          double E_pass_dop = hit.GetDopplerYta(s800sim->AdjustedBeta(simHit.GetBeta()), yta, &track, IPCorepass[0]);
          
          double IPxi = xi;
          if (IPCorepass[0]+1 < nInteractions) IPxi = hit.GetXi(&track,IPCorepass[0],IPCorepass[0]+1);
          else if (IPCorepass[0] != 0) IPxi = hit.GetXi(&track,IPCorepass[0],0);
          
          obj.FillHistogram(dirname, Form("IP_%s_pass_Edop_vs_totPoint",ipgname.c_str()),9,2,11,nInteractions,18,0,360,IPxi*TMath::RadToDeg());
          obj.FillHistogram(dirname, Form("IP_%s_pass_Edop_vs_totPass",ipgname.c_str()),9,1,10,(int) IPCorepass.size(),18,0,360,IPxi*TMath::RadToDeg());
          obj.FillHistogram(dirname, Form("IP_%s_pass_theta_vs_xi",ipgname.c_str()),360,0,TMath::TwoPi(),IPxi,180,0,TMath::Pi(),theta);
          obj.FillHistogram(dirname, Form("IP_%s_pass_Edop_vs_xi",ipgname.c_str()),360,0,TMath::TwoPi(),IPxi,4096,0,4096,E_pass_dop);
          
          if (theta*TMath::RadToDeg() > 55 && theta*TMath::RadToDeg() < 100)
            obj.FillHistogram(dirname, Form("IP_%s_pass_Edop_vs_theta55-100",ipgname.c_str()),360,0,TMath::TwoPi(),IPxi,4096,0,4096,E_pass_dop);
          else 
            obj.FillHistogram(dirname, Form("IP_%s_pass_Edop_vs_thetaNOT55-100",ipgname.c_str()),360,0,TMath::TwoPi(),IPxi,4096,0,4096,E_pass_dop);

          obj.FillHistogram(dirname, Form("IP_%s_Edop_pass",ipgname.c_str()),4096,0,4096,E_pass_dop);
        } 
        else {
          obj.FillHistogram(dirname, Form("IP_%s_Edop_NOpass",ipgname.c_str()),4096,0,4096,energy_corrected);
          obj.FillHistogram(dirname, Form("IP_%s_Edop_NOpass_vs_xi",ipgname.c_str()),360,0,TMath::TwoPi(),xi,4096,0,4096,energy_corrected);
        }
      }

      // obj.FillHistogram(dirname,"gamma_corrected_singles_FEP_nInt>2",4096,0,4096, energy_corrected);
      // obj.FillHistogram(dirname,"gamma_corrected_singles_FEP_vs_xi",360,0,360,xi*TMath::RadToDeg(),4096,0,4096,energy_corrected);
      // obj.FillHistogram(dirname,"gamma_corrected_singles_FEP_vs_scatterAngle",180,0,180,scatterAngle*TMath::RadToDeg(),4096,0,4096,energy_corrected);
      // obj.FillHistogram(dirname,"gamma_corrected_singles_FEP_vs_E-E1/E1",1000,0,10,diffEratio,4096,0,4096,energy_corrected);
      // obj.FillHistogram(dirname,"gamma_corrected_singles_FEP_cosScatterAngle_vs_diffEratio",1000,0,10,diffEratio,200,-1,1,TMath::Cos(scatterAngle));

      // print relative coords to first point
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
      
    }*/
  }

  if(numobj!=list->GetSize()){
    list->Sort();
  }
}
