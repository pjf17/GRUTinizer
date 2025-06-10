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

std::map<int,int> detMap = {
  {26, 0}, {30, 1}, {34, 2}, {38, 3}, {25, 4}, {29, 5}, {33, 6}, {37, 7},
  {27, 8}, {31, 9}, {35,10}, {39,11}, {24,12}, {28,13}, {32,14}, {36,15},
  {47,16}, {63,17}, {71,18}, {79,19}, {51,20}, {59,21}, {67,22}, {83,23},
  {50,24}, {58,25}, {66,26}, {82,27}, {44,28}, {60,29}, {68,30}, {76,31},
  {46,32}, {62,33}, {70,34}, {78,35}, {48,36}, {56,37}, {64,38}, {80,39},
  {49,40}, {57,41}, {65,42}, {81,43}, {45,44}, {61,45}, {69,46}, {77,47}
};

int thetaGate(double angle){
  //attmpet 1
  // if (angle <	0.919) return 0;
  // else if (0.919 < angle && angle < 1.055) return 1;
  // else if (1.055 < angle && angle < 1.189) return 2;
  // else if (1.189 < angle && angle <	1.31)  return 3;
  // else if (1.31	< angle && angle < 1.446)  return 4;
  // else if (1.446 < angle && angle < 1.576) return 5;
  // else if (1.576 < angle && angle < 1.726) return 6;
  // else if (1.726 < angle && angle < 3) return 7;
  //attempt 2
  // if (0.701 <= angle && angle < 0.933) return 0;
  // else if (0.933 <= angle && angle < 1.094) return 1;
  // else if (1.094 <= angle && angle < 1.275) return 2;
  // else if (1.275 <= angle && angle <	1.597)  return 3;
  // else if (1.597	<= angle && angle < 2.02)  return 4;
  if (0.701 <= angle && angle < 0.933) return 0;
  else if (0.933 <= angle && angle < 1.094) return 1;
  else if (1.094 <= angle && angle < 1.275) return 2;
  else if (1.275 <= angle && angle <	1.485)  return 3;
  else if (1.485 <= angle && angle <	1.73)  return 4;
  else if (1.73	<= angle && angle < 2.02)  return 5;
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

int gates_loaded=0;
GCutG *crystal_xy, *crystal_zy = nullptr;
// extern "C" is needed to prevent name mangling.
// The function signature must be exactly as shown here,
//   or else bad things will happen.
extern "C"
void MakeHistograms(TRuntimeObjects& obj) {
//  InitMap();
  TS800Sim *s800sim = obj.GetDetector<TS800Sim>();
  TGretSim *gretsim = obj.GetDetector<TGretSim>();
  TGretina *gretina = obj.GetDetector<TGretina>();
  
  TList *gates = &(obj.GetGates());
  if (gates_loaded!=gates->GetSize()){
    TIter iter(gates);
    while(TObject *obj = iter.Next()) {
      GCutG *gate = (GCutG*)obj;
      std::string tag = gate->GetTag();
      if(!tag.compare("crystal-xy")) {
        crystal_xy = gate;
      } 
      if(!tag.compare("crystal-border-zy")) {
        crystal_zy = gate;
      }
      gates_loaded++;
    }
  }
  
  TList    *list    = &(obj.GetObjects());
  int numobj = list->GetSize();

  static TRandom3 *rand_gen = new TRandom3(59953614);

  bool stopped = false;

  if (!gretina){
    return;
  }
  if ((!s800sim || !s800sim->Size())){
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
  // obj.FillHistogram("simcheck","X",1000,-5,5,simHit.GetX());
  // obj.FillHistogram("simcheck","Y",1000,-5,5,simHit.GetY());
  // obj.FillHistogram("simcheck","Z",1000,-5,5,simHit.GetZ());
  // obj.FillHistogram("simcheck","phi",1000,-180,720,simHit.GetPhi()*TMath::RadToDeg());
  // obj.FillHistogram("simcheck","theta",1000,-360,360,simHit.GetTheta()*TMath::RadToDeg());
  // obj.FillHistogram("simcheck","beta",300,0.2,0.6,simHit.GetBeta());
  // obj.FillHistogram("simcheck","energy",2000,0,2000,simHit.GetEn());

  double gammaEn = GValue::Value("FEP_EN");
  if (std::isnan(gammaEn)) gammaEn = simHit.GetEn();
  double beta = GValue::Value("BETA");
  double simBeta = simHit.GetBeta();
  bool isFEP = simHit.IsFEP();
  obj.FillHistogram("ucgretina","beta_sim",200,0.3,0.6,simBeta);
  //S800 coordinates
  if (!stopped){
    obj.FillHistogram("s800sim","ata", 600,-0.1,0.1, s800sim->GetS800SimHit(0).GetATA());
    obj.FillHistogram("s800sim","bta", 600,-0.1,0.1, s800sim->GetS800SimHit(0).GetBTA());
    obj.FillHistogram("s800sim","yta", 1000,-0.003,0.003, s800sim->GetS800SimHit(0).GetYTA());
    obj.FillHistogram("s800sim","dta", 1000,-0.5,0.5, s800sim->GetS800SimHit(0).GetDTA());
  } else {
    obj.FillHistogram("gretsim","Sim Energies",10000,0,10000,simHit.GetEn());
  }
  // double SIGMA = (2.1*TMath::Exp(-0.1*gammaEn/1000.0) + 60.0*TMath::Exp(-10.2*gammaEn/1000.0));
  double SIGMA = 1.2*(2.21*TMath::Exp(-0.135*gammaEn/1000.0) + 60.0*TMath::Exp(-10.2*gammaEn/1000.0)); //0.92
  // double SIGMA = 0.9*(2.21*TMath::Exp(-0.135*gammaEn/1000.0) + 60.0*TMath::Exp(-10.2*gammaEn/1000.0));
  // double SIGMA = 2.1*TMath::Exp(-0.15*gammaEn/1000.0) + 60.0*TMath::Exp(-10.2*gammaEn/1000.0);
  if(SIGMA > 6.0) {
    SIGMA = 6.0;
  }
  // if(gammaEn < 150.) {
  //   SIGMA = 6.0;
  // }
  // std::vector<double> sigmas = {1.66,1.68,1.7,1.72};
  //SINGLES
  int gSize = gretina->Size();
  for (int i=0; i < gSize; i++){
    TGretinaHit &hit = gretina->GetGretinaHit(i);

    double core_energy = hit.GetCoreEnergy();
    int cryID = hit.GetCrystalId();
    int ringNum = hit.GetRingNumber();
    if (cryID == 77) continue;

    double theta = hit.GetTheta();
    double phi = hit.GetPhi();
    hit.SortSegments();
    // hit.ComptonSort();
    // if (rand_gen->Uniform(0,100) < 24) { //4
    //   //pick a random interaction point
    //   int ip = std::floor(rand_gen->Uniform(1,hit.NumberOfInteractions()-1));
    //   bool result = hit.SwapSegments(0,ip);
    //   if (!result) std::cout<<"AAAAAAAAAAAAAAHHHHHHG!!!"<<std::endl;
    // }

    // double sigMax = 3.0;
    // double sigMin = 2.0;
    // double sigStep = 0.1;
    // double sig = sigMin;
    // int sigBins = (sigMax-sigMin)/sigStep + 1;
    // TVector3 local_pos_loop(hit.GetLocalPosition(0));
    // while (sig < sigMax) {
    //   double smear_x = local_pos_loop.X() + rand_gen->Gaus(0, sig); 
    //   double smear_y = local_pos_loop.Y() + rand_gen->Gaus(0, sig);
    //   double smear_z = local_pos_loop.Z() + rand_gen->Gaus(0, sig);
    //   hit.SetPosition(0,smear_x,smear_y,smear_z);
    //   obj.FillHistogram(dirname,"dopEn_Sig",sigBins*2,sigMin,sigMax,sig,6000,0,6000,hit.GetDoppler(beta));
    //   sig += sigStep;
    // }

    TVector3 local_pos(hit.GetLocalPosition(0));
    double smear_x = local_pos.X() + rand_gen->Gaus(0, SIGMA); 
    double smear_y = local_pos.Y() + rand_gen->Gaus(0, SIGMA);
    double smear_z = local_pos.Z() + rand_gen->Gaus(0, SIGMA);
    hit.SetPosition(0,smear_x,smear_y,smear_z); //this resets the positions also in NNaddback

    // if (gates && crystal_xy && crystal_zy){
    //   bool inXY = crystal_xy->IsInside(smear_x,smear_y);
    //   bool inZY = crystal_zy->IsInside(smear_z,smear_y);
    //   obj.FillHistogram("positionsmear","inout",2,0,2,int(inXY && inZY));
    //   if (inXY && inZY){
    //     obj.FillHistogram("positionsmear",Form("X_vs_Y_smear_inside"),200,-100,100,smear_x,200,-100,100,smear_y); 
    //     obj.FillHistogram("positionsmear",Form("Z_vs_Y_smear_inside"),200,-100,100,smear_z,200,-100,100,smear_y);
    //   } else {
    //     if (!inXY) obj.FillHistogram("positionsmear",Form("X_vs_Y_smear_outside"),200,-100,100,smear_x,200,-100,100,smear_y); 
    //     if (!inZY) obj.FillHistogram("positionsmear",Form("Z_vs_Y_smear_outside"),200,-100,100,smear_z,200,-100,100,smear_y);
    //   }
    // }
    double dop_corrected;

    if (!stopped || !std::isnan(beta)){
      // dop_corrected = hit.GetDoppler(s800sim->AdjustedBeta(beta), &track); 
      dop_corrected = hit.GetDoppler(beta, &track); 
      // dop_corrected = hit.GetDopplerYta(beta); 
    } 
    else{
      dop_corrected = hit.GetCoreEnergy();
    }

    obj.FillHistogram(dirname,"HitPhi_v_HitTheta",180,0,180,180-theta*TMath::RadToDeg(),360,0,360,phi*TMath::RadToDeg());
    // obj.FillHistogram(dirname,Form("HitPhi_v_HitTheta_Smeared_c%d",detMap[cryID]),180,0,180,180-theta_smear*TMath::RadToDeg(),360,0,360,phi_smear*TMath::RadToDeg());
    obj.FillHistogram(dirname,"CoreEnergy",10000,0,10000,core_energy);

    obj.FillHistogram(dirname,"Theta_vs_Energy_btyd",4000,0,4000,dop_corrected,100,0,3,hit.GetTheta());
    if (hit.NumberOfInteractions() > 1 && isFEP) obj.FillHistogram(dirname,"Fep_dopEn_vs_xi",360,0,TMath::TwoPi(),hit.GetXi(),4000,0,4000,dop_corrected);

    // obj.FillHistogram(dirname,"EmissionAngle_vs_DetectedAngle",180,0,180,hit.GetTheta()*TMath::RadToDeg(),180,0,180,simHit.GetTheta()*TMath::RadToDeg());
    
    bool bin75 = hit.GetTheta()*TMath::RadToDeg() > 65 && hit.GetTheta()*TMath::RadToDeg() < 85;
    bool bin100 = hit.GetTheta()*TMath::RadToDeg() > 90 && hit.GetTheta()*TMath::RadToDeg() < 115;
    //fitting hists
    obj.FillHistogram(dirname,"dopEn",10000,0,10000,dop_corrected);

    if (bin75) obj.FillHistogram(dirname,"dopEn_bin75",10000,0,10000,dop_corrected);
    if (bin100) obj.FillHistogram(dirname,"dopEn_bin100",10000,0,10000,dop_corrected);
    obj.FillHistogram(dirname,Form("dopEn_%d",thetaGate(hit.GetTheta())),10000,0,10000,dop_corrected);
    double xi = hit.GetXi(&track);
    double simpleXi = xi;
    if (simpleXi > TMath::Pi()/2) simpleXi = TMath::Pi() - simpleXi;
    bool perp = simpleXi*TMath::RadToDeg() > 60 && hit.NumberOfInteractions() > 1;
    bool para = simpleXi*TMath::RadToDeg() < 30 && hit.NumberOfInteractions() > 1;
    if (hit.NumberOfInteractions() > 1) {
      if (perp) obj.FillHistogram(dirname, "dopEn_xi_perp",10000,0,10000,dop_corrected);
      if (para) obj.FillHistogram(dirname, "dopEn_xi_para",10000,0,10000,dop_corrected);
      if (para || perp) obj.FillHistogram(dirname, "dopEn_xi_psum",10000,0,10000,dop_corrected);
      // obj.FillHistogram(dirname,Form("dopEn_xi_%d",(int) (xi/TMath::Pi()*5)),10000,0,10000,dop_corrected); 
      obj.FillHistogram(dirname,"dopEn_nint>1",10000,0,10000,dop_corrected);
    }
    // obj.FillHistogram(dirname,"dopEn_vs_theta",180,0,TMath::Pi(),theta,10000,0,10000,dop_corrected);
    if (detMap[cryID] < 4) obj.FillHistogram(dirname,"dopEn_fwd",10000,0,10000,dop_corrected);
    if (detMap[cryID] > 43) obj.FillHistogram(dirname,"dopEn_bkwd",10000,0,10000,dop_corrected);
    // if (cryID > 40) obj.FillHistogram(dirname,"dopEn_90qds_B&T&Y&D",10000,0,10000,dop_corrected);
    // if (gates && crystal_xy && crystal_zy && crystal_xy->IsInside(smear_x,smear_y) && crystal_zy->IsInside(smear_z,smear_y)) 
    //   obj.FillHistogram(dirname,"dopEn_inside",10000,0,10000,energy_track_yta_dta);
    // else obj.FillHistogram(dirname,"dopEn_outside",10000,0,10000,energy_track_yta_dta);

    //angular distributions
    // obj.FillHistogram("polarization",Form("dopEn_vs_theta_%d",cryID),180,0,TMath::Pi(),hit.GetTheta(),10000,0,10000,dop_corrected);
    // if (isFEP) obj.FillHistogram(dirname, "fep_ringnum",6,4,16,hit.GetRingNumber());
    // obj.FillHistogram("angdist", "dopEn_vs_ringnum",6,4,16,hit.GetRingNumber(),10000,0,10000,dop_corrected);
    // if (hit.NumberOfInteractions() > 1) obj.FillHistogram("poldist", "dopEn_vs_ringnum",6,4,16,hit.GetRingNumber(),10000,0,10000,dop_corrected);

    //SUMMARY SPECTRUM
    // obj.FillHistogram(dirname,"dop_btyd_summary",48,0,48,detMap[cryID],3000,0,3000,dop_corrected);    
    
    if (isFEP){ //full energy peak event
      obj.FillHistogram(dirname,"dopEn_fep",10000,0,10000,dop_corrected); 
      if (bin75) obj.FillHistogram(dirname,"dopEn_bin75_fep",10000,0,10000,dop_corrected);
      if (bin100) obj.FillHistogram(dirname,"dopEn_bin100_fep",10000,0,10000,dop_corrected);
      if (perp) obj.FillHistogram(dirname, "dopEn_xi_perp_fep",10000,0,10000,dop_corrected);
      if (para) obj.FillHistogram(dirname, "dopEn_xi_para_fep",10000,0,10000,dop_corrected);
      if (para || perp) obj.FillHistogram(dirname, "dopEn_xi_psum_fep",10000,0,10000,dop_corrected);
      obj.FillHistogram(dirname,Form("theta_%d_fep",thetaGate(hit.GetTheta())),10000,0,10000,dop_corrected);
      // if (cryID > 40) obj.FillHistogram(dirname,"dopEn_90qds_fep",10000,0,10000,dop_corrected);     
    } else {
      obj.FillHistogram(dirname,"dopEn_bg",10000,0,10000,dop_corrected);
      if (bin75 && gSize == 1) obj.FillHistogram(dirname,"dopEn_bin75_bg",10000,0,10000,dop_corrected);
      if (bin100 && gSize == 1) obj.FillHistogram(dirname,"dopEn_bin100_bg",10000,0,10000,dop_corrected);
      if (perp) obj.FillHistogram(dirname, "dopEn_xi_perp_bg",10000,0,10000,dop_corrected);
      if (para) obj.FillHistogram(dirname, "dopEn_xi_para_bg",10000,0,10000,dop_corrected);
      if (para || perp) obj.FillHistogram(dirname, "dopEn_xi_psum_bg",10000,0,10000,dop_corrected);
      obj.FillHistogram(dirname,Form("dopEn_%d_bg",thetaGate(hit.GetTheta())),10000,0,10000,dop_corrected);
      // if (cryID > 40) obj.FillHistogram(dirname,"dopEn_90qds_bg",10000,0,10000,dop_corrected);
    }
  }

  //NNADDBACK
  //loop over multiplicity
  /*
  dirname = "addback";
  int nnSize = gretina->NNAddbackSize();
  for (int i=0; i < nnSize; i++){

    //get hit and hit data 
    TGretinaHit nnhit = gretina->GetNNAddbackHit(i);
    int ringNum = nnhit.GetRingNumber();
    int cryID = nnhit.GetCrystalId();
    double core_energy = nnhit.GetCoreEnergy();
    // if (cryID == 77) continue;

    if (efficiencyCorrection(rand_gen,nnhit,nnhit.GetNNeighborHits())){
      double energy_track_yta_dta;
      double energy_track_yta;
      double energy_track;
      
      // TVector3 local_pos(nnhit.GetLocalPosition(0));
      // double smear_x = local_pos.X() + rand_gen->Gaus(0, SIGMA); 
      // double smear_y = local_pos.Y() + rand_gen->Gaus(0, SIGMA);
      // double smear_z = local_pos.Z() + rand_gen->Gaus(0, SIGMA);
      // nnhit.SetPosition(0,smear_x,smear_y,smear_z);

      if (!stopped){
        energy_track = nnhit.GetDoppler(beta, &track);
        energy_track_yta = nnhit.GetDopplerYta(beta, yta, &track);
        energy_track_yta_dta = nnhit.GetDopplerYta(s800sim->AdjustedBeta(beta), yta, &track);
      } else {
        energy_track = energy_track_yta = energy_track_yta_dta = nnhit.GetCoreEnergy();
      }

      //compare uncorrected energy reported from geant with core energy. if FEP they should be the same
      bool isNNFEP = fabs(simHit.GetEn() - nnhit.GetCoreEnergy()) < 1.5;

      // obj.FillHistogram(dirname, Form("gretina_n%s_summary",multiplicity),48,0,48,detMap[cryID],4000,0,4000,energy_track_yta_dta);
      obj.FillHistogram(dirname, Form("gretina_n%d_CoreEnergy",nnhit.GetABDepth()), 10000,0,10000, core_energy);
      obj.FillHistogram(dirname, Form("gretina_n%d_B&T",nnhit.GetABDepth()), 10000,0,10000, energy_track);
      obj.FillHistogram(dirname, Form("gretina_n%d_B&T&Y",nnhit.GetABDepth()), 10000,0,10000, energy_track_yta);
      obj.FillHistogram(dirname, Form("gretina_n%d_B&T&Y&D",nnhit.GetABDepth()), 10000,0,10000, energy_track_yta_dta);
      if (isNNFEP){
        obj.FillHistogram(dirname, Form("gretina_n%d_B&T&Y&D_fep",nnhit.GetABDepth()), 10000,0,10000, energy_track_yta_dta);
      } else {
        obj.FillHistogram(dirname, Form("gretina_n%d_B&T&Y&D_bg",nnhit.GetABDepth()), 10000,0,10000, energy_track_yta_dta);
      }
      
      //exclude the ng spectrum (n==3)
      if (nnhit.GetABDepth() < 3){
        // obj.FillHistogram(dirname, "gretina_ab_CoreEnergy", 10000,0,10000, core_energy);
        // obj.FillHistogram(dirname, "gretina_ab_B&T", 10000,0,10000, energy_track);
        // obj.FillHistogram(dirname, "gretina_ab_B&T&Y", 10000,0,10000, energy_track_yta);
        // obj.FillHistogram(dirname, "gretina_ab_B&T&Y&D", 10000,0,10000, energy_track_yta_dta);
        // if (isNNFEP){
        //   obj.FillHistogram(dirname, "gretina_ab_B&T&Y&D_fep", 10000,0,10000, energy_track_yta_dta);
        // } else {
        //   obj.FillHistogram(dirname, "gretina_ab_B&T&Y&D_bg", 10000,0,10000, energy_track_yta_dta);
        // }
      }

      if (nnhit.GetABDepth() == 1) {
        if (cryID > 40)
          obj.FillHistogram(dirname,"gretina_pol_B&T&Y&D", 8192,0,8192, energy_track_yta_dta);
        //POLARIZATION
        std::string swaptype = "pol";
        for (int t=0; t < 1; t++){
          if (t==1) {
            swaptype = "pol_swapped";
            TGretinaHit swap = nnhit.GetNeighbor();
            swap.NNAdd(nnhit.GetInitialHit());
            energy_track = swap.GetDoppler(beta, &track);
            energy_track_yta = swap.GetDopplerYta(beta, yta, &track);
            energy_track_yta_dta = swap.GetDopplerYta(s800sim->AdjustedBeta(beta), yta, &track);
          }
          std::string polColor = "blank";
          if (PairHit(nnhit,redPairs)) polColor = "red";
          else if (PairHit(nnhit,goldPairs)) polColor = "gold";
          else if (PairHit(nnhit,bluePairs)) polColor = "blue";

          std::string quadType = "blank";
          if (PairHit(nnhit,OneQuadPlus)) quadType = "qd1+";
          else if (PairHit(nnhit,OneQuadDefault)) quadType = "qd1";
          else if (PairHit(nnhit,TwoQuadPairs)) quadType = "qd2";

          if ( polColor.compare("blank") != 0 ){
            obj.FillHistogram(dirname,Form("gretina_%s_%s_B&T&Y&D",swaptype.c_str(),polColor.c_str()), 8192,0,8192, energy_track_yta_dta);
            // obj.FillHistogram(dirname,Form("gretina_%s_%s",swaptype.c_str(),quadType.c_str()), 8192,0,8192, energy_track_yta_dta);
            // obj.FillHistogram(dirname,Form("gretina_%s_%s_%s_B&T&Y&D",swaptype.c_str(),polColor.c_str(),quadType.c_str()), 8192,0,8192, energy_track_yta_dta);
            if (isNNFEP) {
              obj.FillHistogram(dirname,Form("gretina_%s_%s_B&T&Y&D_fep",swaptype.c_str(),polColor.c_str()), 8192,0,8192, energy_track_yta_dta);
              // obj.FillHistogram(dirname,Form("gretina_%s_%s_B&T&Y&D_fep",swaptype.c_str(),quadType.c_str()), 8192,0,8192, energy_track_yta_dta);
              // obj.FillHistogram(dirname,Form("gretina_%s_%s_%s_B&T&Y&D_fep",swaptype.c_str(),polColor.c_str(),quadType.c_str()), 8192,0,8192, energy_track_yta_dta);
            }
            else {
              obj.FillHistogram(dirname,Form("gretina_%s_%s_B&T&Y&D_bg",swaptype.c_str(),polColor.c_str()), 8192,0,8192, energy_track_yta_dta);
              // obj.FillHistogram(dirname,Form("gretina_%s_%s_B&T&Y&D_bg",swaptype.c_str(),quadType.c_str()), 8192,0,8192, energy_track_yta_dta);
              // obj.FillHistogram(dirname,Form("gretina_%s_%s_%s_B&T&Y&D_bg",swaptype.c_str(),polColor.c_str(),quadType.c_str()), 8192,0,8192, energy_track_yta_dta);
            }
          }
        }
      }
    }
  }*/

  if(numobj!=list->GetSize()){
    list->Sort();
  }
}
