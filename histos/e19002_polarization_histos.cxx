#include "TRuntimeObjects.h"

#include <iostream>
#include <map>
#include <cstdio>
#include <vector>
#include <string>

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

double GetAfp(double crdc_1_x,double  crdc_2_x){
  return TMath::ATan( (crdc_2_x - crdc_1_x)/1073.0 );
}

int thetaGate(double angle){
  //attempt 1
  // if (0.701 <= angle && angle < 0.933) return 0;
  // else if (0.933 <= angle && angle < 1.094) return 1;
  // else if (1.094 <= angle && angle < 1.275) return 2;
  // else if (1.275 <= angle && angle <	1.485)  return 3;
  // else if (1.485 <= angle && angle <	1.73)  return 4;
  // else if (1.73	<= angle && angle < 2.02)  return 5;
  //attempt 2
  if (angle < 0.933) return 0;
  else if (0.933 <= angle && angle < 1.094) return 1;
  else if (1.094 <= angle && angle < 1.275) return 2;
  else if (1.275 <= angle && angle <	1.485)  return 3;
  else if (1.485 <= angle && angle <	1.73)  return 4;
  else if (1.73	<= angle)  return 5;
}

double lastPointPenalty(double x){
  double val = TMath::TanH((x/100-1.1)/0.5)*exp(-3*(x/100-1.1)+1);
  if (val < 0) val = 0;
  return val + 1;
}

void comptonSortParts(const TGretinaHit &ghit, int &FP, int &SP, int PART=1) {
  int N = ghit.NumberOfInteractions();
  if (N < 2) return;
  
  double FOM = 1e10;
  FP = 0; 
  SP = 1;
  double E = ghit.GetCoreEnergy();
  
  //scale the interaction points so they match the core energy
  double scaleFactor = 0;
  for (int i=0; i < N; i++) scaleFactor += ghit.GetSegmentEng(i);
  scaleFactor = E/scaleFactor;

  //find the best first two interaction points that satisfy the minimization function
  for (int fp=0; fp < N; fp++){
    double E1 = ghit.GetSegmentEng(fp)*scaleFactor;
    double er = 511.0/E * E1/(E - E1);

    for (int sp=0; sp < N; sp++){
      if (fp == sp) continue;
      
      double cosp = TMath::Cos(ghit.GetScatterAngle(fp,sp));
      double E2 = ghit.GetSegmentEng(fp)*scaleFactor;
      double x = er + cosp;
      double kn = pow((E - E1)/E,2)*((E-E1)/E + E/(E-E1) - pow(TMath::Sin(ghit.GetScatterAngle(fp,sp)),2));
      
      double ffom = std::pow(std::abs(1-x),2.0/3)*ghit.GetAlpha(fp,sp)*std::pow(E/E1,3); //compton scattering
      if (PART>=2) ffom *= kn;
      if (PART>=3) ffom *= std::pow(E/E2,2)*lastPointPenalty(E1)*lastPointPenalty(E2);
      if (PART>=4) ffom *= ghit.GetLocalPosition(fp).Z()*ghit.GetLocalPosition(sp).Z();

      if (ffom < FOM) {
        FOM = ffom;
        FP = fp;
        SP = sp;
      }
    }
  }
}

//Get the Ion Chamber DE depending on whether IC_DE_XTILT is set
double GetGoodICE(TS800 *s800){
  static int ncalls = 0;
  double value = 0;
  double crdc_1_x = s800->GetCrdc(0).GetDispersiveX();
  double crdc_1_y = s800->GetCrdc(0).GetNonDispersiveY();
  
  double xtilt = GValue::Value("IC_DE_XTILT");
  double x0tilt = GValue::Value("IC_DE_X0TILT");
  double ytilt = GValue::Value("IC_DE_YTILT");
  if (!std::isnan(xtilt) && !std::isnan(x0tilt) && !std::isnan(ytilt)){
    value = s800->GetIonChamber().GetdE(crdc_1_x, crdc_1_y);
  } else {
    value = s800->GetIonChamber().GetAve();
    if (ncalls == 0){
      std::cout<<"XTILT, X0TILT, YTILT NOT SET SWITCHING TO GETAVE()\n";
      ncalls++;
    }
  }
  
  return value;
}

bool checkEnergyTheta(TF1* fdop, TF1* fres, double energy, double theta){
  return fdop->Eval(theta) + fres->Eval(theta) > energy && fdop->Eval(theta) - fres->Eval(theta) < energy;
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

void CheckGates(std::vector<GCutG*> gates, std::vector<unsigned short> &passed, double x, double y){
  unsigned short ngates = gates.size();
  for (unsigned short i=0; i < ngates; i++){
    if (gates.at(i)->IsInside(x,y)) passed.push_back(i);
  }
  return;
}

TVector3 *randomBeam(TRandom3 *rand){
  double theta = TMath::ACos(2*rand->Uniform()-1);
  double phi = rand->Uniform()*TMath::TwoPi();
  return new TVector3(TMath::Sin(theta)*TMath::Cos(phi),TMath::Sin(theta)*TMath::Sin(phi),TMath::Cos(theta));
}

bool gates_loaded = false;
bool messageGiven = false;
std::map<std::string,std::vector<GCutG*>> gates;

// extern "C" is needed to prevent name mangling.
// The function signature must be exactly as shown here,
//   or else bad things will happen.
extern "C"
void MakeHistograms(TRuntimeObjects& obj) {
  TGretina *gretina = obj.GetDetector<TGretina>();
  TBank29  *bank29  = obj.GetDetector<TBank29>();
  TS800    *s800    = obj.GetDetector<TS800>();
  TList    *list    = &(obj.GetObjects());
  
  if (!s800){
    return;
  }
  
  int numobj = list->GetSize();

  //load in the gates
  if (!gates_loaded) {
    LoadGates(&(obj.GetGates()),gates);
    gates_loaded = true;
  }
  
  //Use this spectrum for the time-energy cut for GRETINA
  if(bank29 && gretina) {
    for(unsigned int i=0;i<gretina->Size();i++) {
      //Time-energy cut
      TGretinaHit &hit = gretina->GetGretinaHit(i);
      obj.FillHistogram("Bank29","Gretina_dop_t0_Bank29_time",
          600,-600,600,bank29->Timestamp()-hit.GetTime(),
          2500,0,10000, hit.GetCoreEnergy());
    }//loop over gretina hits
  }//bank29 and gretina exist

  static TRandom3 *rand_gen = new TRandom3(249783);

  //---------------------------------------------------------------
  //UNGATED

  unsigned short bits = s800->GetTrigger().GetRegistr();
  bool badevent = false;
  for(int j=0;j<16;j++) {
    if(((bits>>j)&0x0001)) {
      if (j > 1) badevent = true;
      if (!badevent) obj.FillHistogram("ungated","trig_bit",20,0,20,j);
    }
  }
  if (badevent) return;

  //MAKE RAW TOF HISTS
  double raw_obj = s800->GetRawOBJ_MESY();
  double raw_e1 = s800->GetRawE1_MESY();
  double raw_xf = s800->GetRawXF_MESY();
  
  obj.FillHistogram("ungated", "MTOF_OBJE1", 5000, -10000, 0, raw_obj - raw_e1);
  obj.FillHistogram("ungated", "MTOF_XFE1", 4000, -4000, 6000, raw_xf - raw_e1);

  //MAKE INCOMING PID
  double tof_obje1 = s800->GetMTof().GetCorrelatedObjE1(); 
  double tof_xfpe1 = s800->GetMTof().GetCorrelatedXfpE1();
  obj.FillHistogram("ungated", "incoming_pid", 600, -2800, -2200, tof_obje1, 600, 1300, 1900, tof_xfpe1);                                              

  //CRDC PLOTS
  double crdc_1_x = s800->GetCrdc(0).GetDispersiveX();
  double crdc_2_x = s800->GetCrdc(1).GetDispersiveX();
  double crdc_1_y = s800->GetCrdc(0).GetNonDispersiveY();
  double crdc_2_y = s800->GetCrdc(1).GetNonDispersiveY();
  double afp = GetAfp(crdc_1_x, crdc_2_x);
  
  double ylow = -200;
  double yhigh = 200;
  double ybins = 400;
  
  double yslope = GValue::Value("CRDC1_Y_SLOPE");
  if (std::isnan(yslope) || yslope == 0){
    ylow = 0;
    yhigh = 1500;
    ybins = 1500;
  }
  obj.FillHistogram("ungated", "crdc1 X_Y", 600, -300, 300, crdc_1_x, ybins, ylow, yhigh, crdc_1_y);  
  obj.FillHistogram("ungated", "crdc2 X_Y", 600, -300, 300, crdc_2_x, ybins, ylow, yhigh, crdc_2_y);

  double s800time = s800->GetTimestamp()/1E8/3600;
  obj.FillHistogram("ungated", "crdc1Y_vs_timestamp",120,0,2,s800time, ybins, ylow, yhigh, crdc_1_y);  
  obj.FillHistogram("ungated", "crdc2Y_vs_timestamp",120,0,2,s800time, ybins, ylow, yhigh, crdc_2_y);

  //---------------------------------------------------------------
  //GATED
  std::vector<unsigned short> incoming_passed;
  CheckGates(gates["incoming"],incoming_passed,tof_obje1,tof_xfpe1);
  
  //GET TOFs AND DE FOR PID AND CORRELATION PLOTS
  double tof_obje1_corr = 0;
  double ic_energy = GetGoodICE(s800);
  
  //obje1 xfp correlation hist binning parameters
  double xocor_lowbin = -2028;
  double xocor_highbin = 8192;
  int xocor_nbins = 4096;
  double xfp_obj = tof_xfpe1 - tof_obje1;
  double xfp_obj_shift = GValue::Value("TOFXFP_OBJ_SHIFT");
  if (!std::isnan(xfp_obj_shift)){
    xfp_obj -= xfp_obj_shift;
    xocor_lowbin = -128;
    xocor_highbin = 128;
    xocor_nbins = 256;
  }

  std::vector<unsigned short> outgoing_passed;
  
  //check if proper TOF GValues are set to work with outoing PID
  
  if (std::isnan(GValue::Value("OBJ_MTOF_CORR_AFP")) || std::isnan(GValue::Value("OBJ_MTOF_CORR_XFP"))){
    if (!messageGiven) {
      std::cout<<"OBJ_MTOF_CORR_AFP OR OBJ_MTOF_CORR_XFP NOT SET, SKIPPING OUTGOING GATES"<<std::endl;
      messageGiven = true;
    }
  } 
  else {
    tof_obje1_corr = s800->GetMTofObjE1();
    CheckGates(gates["outgoing"],outgoing_passed,tof_obje1_corr,ic_energy);
  }

  std::string dirname  = "";
  for (auto ind_out : outgoing_passed){
    dirname = Form("%s_gated", gates["outgoing"].at(ind_out)->GetName());
    obj.FillHistogram(dirname, "incoming_pid", 600, -2800, -2200, tof_obje1, 600, 1300, 1900, tof_xfpe1);   
  }

  //---------------------------------------------------------------
  //INCOMING
  for (auto ind_in: incoming_passed){
    dirname = Form("%s_gated", gates["incoming"].at(ind_in)->GetName());
    obj.FillHistogram(dirname, "outgoing_pid", 600, -2800, -2200, tof_obje1_corr, 2048, 0, 4096, ic_energy);                                                
    obj.FillHistogram(dirname, "yta", 500,-100,100,s800->GetYta());                                                

    //time energy spectrum
    if (bank29 && gretina){
      int gSize = gretina->Size();
      double timeBank29 = bank29->Timestamp();
      for (int i=0; i < gSize; i++){
        TGretinaHit &hit = gretina->GetGretinaHit(i);
        obj.FillHistogram(dirname,"core_energy_t0_Bank29_time",
            1600,-400,400,timeBank29-hit.GetTime(),
            2500,0,10000,hit.GetCoreEnergy());
      }
    }
    
    //---------------------------------------------------------------
    //OUTGOING
    for (auto ind_out : outgoing_passed){
      dirname = Form("%s_%s_gated", gates["incoming"].at(ind_in)->GetName(), gates["outgoing"].at(ind_out)->GetName());

      //CRDC PLOTS
      obj.FillHistogram(dirname, "crdc1_XvsY", 600, -300, 300, crdc_1_x, 1000, -500, 500, crdc_1_y);  
      obj.FillHistogram(dirname, "crdc2_XvsY", 600, -300, 300, crdc_2_x, 1000, -500, 500, crdc_2_y);
      
      //CORRELATION PLOTS
      obj.FillHistogram(dirname, "corrobje1_crdc1x", 3000, -6000, 0, tof_obje1_corr,600, -300, 300, crdc_1_x);
      obj.FillHistogram(dirname, "corrobje1_afp", 3000, -6000, 0, tof_obje1_corr,1000, -0.1, 0.1, afp);

      //GAMMA RAY ANALYSIS
      if (gretina){
        TVector3 track = s800->Track();
        if (bank29){
          double timeBank29 = bank29->Timestamp();

          //get Beta
          double outgoingBeta = GValue::Value(Form("BETA_%s",gates["outgoing"].at(ind_out)->GetName()));
          if (std::isnan(outgoingBeta)) outgoingBeta = GValue::Value("BETA");
          // if (crdc_1_x > 6) continue; 
          
          //SINGLES
          int gSize = gretina->Size();
          std::vector<int> prompt_gammas;
          double total_corrected_energy = 0;
          double total_core_energy = 0;
          for (int i=0; i < gSize; i++){
            TGretinaHit &hit = gretina->GetGretinaHit(i);
            if (hit.GetPad() > 0) continue;
            TGretinaHit hitMain;
            hit.Copy(hitMain);
            hit.PICCSort();
            double energy_corrected = hit.GetDopplerYta(outgoingBeta, s800->GetYta(), &track);
            double energy_corrected_main = hitMain.GetDopplerYta(outgoingBeta, s800->GetYta(), &track);
            double core_energy = hit.GetCoreEnergy();
            // double energy_corrected_main = hitMain.GetDopplerYta(outgoingBeta, s800->GetYta(), &track);
            // double energy = hit.GetDoppler(outgoingBeta);
            // double theta = hit.GetTheta();
            // double phi = hit.GetPhi();
            int cryID = hit.GetCrystalId();
            int nInteractions = hit.NumberOfInteractions();
            // std::vector<bool> momentum_gates;
            // for (int mg=0; mg < 5; mg++){
            //   momentum_gates.push_back(crdc_1_x > 80 - mg*40);
            // }
            
            //PROMPT GATE
            bool tgate = false;
            if (gates["prompt"].size() > 0) tgate = gates["prompt"][0]->IsInside(timeBank29-hit.GetTime(), core_energy);
            else std::cout<<"NO PROMPT GATE LOADED\n";

            bool bin75 = hit.GetTheta()*TMath::RadToDeg() > 65 && hit.GetTheta()*TMath::RadToDeg() < 85;
            bool bin100 = hit.GetTheta()*TMath::RadToDeg() > 90 && hit.GetTheta()*TMath::RadToDeg() < 115;
            
            if (tgate) {
              prompt_gammas.push_back(i);
              obj.FillHistogram(dirname, "crdcx1",1000,-500,500, crdc_1_x);
              obj.FillHistogram(dirname, "gam_dop_sgl_prompt",10000,0,10000, energy_corrected);
              obj.FillHistogram(dirname, "gam_core_sgl_prompt",10000,0,10000, core_energy);
              if (nInteractions < 5 && nInteractions > 1) obj.FillHistogram(dirname, Form("gam_dop_sgl_prompt_nint%d",nInteractions),4096,0,4096, energy_corrected);
              obj.FillHistogram(dirname, "gam_dop_sgl_prompt_vs_theta",180,0,TMath::Pi(),hit.GetTheta(),4096,0,4096, energy_corrected);
              
              if (bin75) {
                obj.FillHistogram(dirname, "gam_dop_sgl_bin75",10000,0,10000, energy_corrected);
                obj.FillHistogram(dirname, "theta_bin75",180,0,180, hit.GetTheta()*TMath::RadToDeg());
              }
              if (bin100) {
                obj.FillHistogram(dirname, "gam_dop_sgl_bin100",10000,0,10000, energy_corrected);
                obj.FillHistogram(dirname, "theta_bin100",180,0,180, hit.GetTheta()*TMath::RadToDeg());
              }
              obj.FillHistogram(dirname, Form("gam_dop_sgl_%d",thetaGate(hit.GetTheta())),10000,0,10000, energy_corrected);
              obj.FillHistogram(dirname, "gam_sngl_main",4096,0,4096, energy_corrected_main);
              // for (int mg=0; mg < 5; mg++) {
                // if (momentum_gates[mg]) obj.FillHistogram(dirname, Form("gam_dop_sgl_%d_%d",mg,thetaGate(hit.GetTheta())),10000,0,10000, energy_corrected);
              // }
              // obj.FillHistogram(dirname, Form("gam_dop_vs_theta_%d",thetaGate(hit.GetTheta())),80,0,TMath::Pi(),hit.GetTheta(),4096,0,4096, energy_corrected);
            }

            if (tgate && nInteractions > 1){
              int nXiPlotEbins = 10000;
              double XiPlotLo = 0;
              double XiPlotHi = 10000;
              obj.FillHistogram(dirname, "gam_sngl_nint>1",nXiPlotEbins,XiPlotLo,XiPlotHi, energy_corrected);

              obj.FillHistogram(dirname, "gam_sngl_cosxiprime_vs_cosxi",200,-1,1,TMath::Cos(hit.GetXi(&track,false,0,1)),200,-1,1,TMath::Cos(hit.GetXi(&track,false,1,0)));
              double xi = hitMain.GetXiSimple(&track); //EDIT THIS BACK!!!!!
              bool perp = xi*TMath::RadToDeg() > 60;
              bool para = xi*TMath::RadToDeg() < 30;
              if (perp) obj.FillHistogram(dirname, "prompt_gamma_dop_xi_perp",10000,0,10000,energy_corrected);
              if (para) obj.FillHistogram(dirname, "prompt_gamma_dop_xi_para",10000,0,10000,energy_corrected); //THESE ARE THE GOOD ONES, SWITCH BACK EVENTUALLY
              // if (perp) obj.FillHistogram(dirname, "prompt_gamma_dop_xi_perp",10000,0,10000,energy_corrected_main);
              // if (para) obj.FillHistogram(dirname, "prompt_gamma_dop_xi_para",10000,0,10000,energy_corrected_main);

              double ximain = hitMain.GetXiSimple(&track);
              if (ximain*TMath::RadToDeg() > 60) 
                obj.FillHistogram(dirname, "prompt_gamma_dop_xi_main_perp",10000,0,10000,energy_corrected_main);
              if (ximain*TMath::RadToDeg() < 30) 
                obj.FillHistogram(dirname, "prompt_gamma_dop_xi_main_para",10000,0,10000,energy_corrected_main);

              obj.FillHistogram(dirname, "COMP_gam_dop_sgl_prompt_NINT>1",4096,0,4096, energy_corrected);
              obj.FillHistogram(dirname, "COMP_gam_dop_sgl_prompt_main_NINT>1",4096,0,4096, energy_corrected_main);
              int partfp = 0, partsp = 1;
              comptonSortParts(hitMain,partfp,partsp,4);
              if (partfp != 0) {
                obj.FillHistogram(dirname, "COMP_gam_dop_sgl_prompt_NINT>1_diff",4096,0,4096, energy_corrected);
                obj.FillHistogram(dirname, "COMP_gam_dop_sgl_prompt_main_NINT>1_diff",4096,0,4096, energy_corrected_main);
              }

              // for (int mg=0; mg < 5; mg++) {
                // if (momentum_gates[mg]) {
                  // if (perp) obj.FillHistogram(dirname, Form("prompt_gamma_dop_xi_%d_perp",mg),10000,0,10000,energy_corrected);
                  // if (para) obj.FillHistogram(dirname, Form("prompt_gamma_dop_xi_%d_para",mg),10000,0,10000,energy_corrected);
                // }
              // }
              // if (para || perp) obj.FillHistogram(dirname, "prompt_gamma_dop_xi_psum",10000,0,10000,energy_corrected);

              obj.FillHistogram(dirname, "gam_sngl_vs_xi_main",180,0,180,hitMain.GetXi(&track)*TMath::RadToDeg(),4096,0,4096, energy_corrected_main);
              obj.FillHistogram(dirname, "gam_sngl_main_nint>1",4096,0,4096, energy_corrected_main);
              obj.FillHistogram(dirname, "gam_sngl_vs_xi",180,0,180,hit.GetXi(&track)*TMath::RadToDeg(),4096,0,4096, energy_corrected);
              // if (hitMain.GetXi(&track)*TMath::RadToDeg() > 90) 
                // obj.FillHistogram(dirname, "gam_sngl_xi>90",4096,0,4096, energy_corrected_main);
              // else
                // obj.FillHistogram(dirname, "gam_sngl_xi<90",4096,0,4096, energy_corrected_main);
              // 
              // if (hit.GetXiSimple(&track)*TMath::RadToDeg() > 75 && hit.GetAlpha()*TMath::RadToDeg() < 2 && nInteractions < 3)
                // obj.FillHistogram(dirname, "gam_dop_sngl_super_filtered",4096,0,4096, energy_corrected);
              // for (int inint=3; inint <= 5; inint++){
              //   if (!(nInteractions < inint)) continue;
              //   if (perp) obj.FillHistogram(dirname, Form("prompt_gamma_dop_xi_perp_<%d",inint),10000,0,10000,energy_corrected);
              //   if (para) obj.FillHistogram(dirname, Form("prompt_gamma_dop_xi_para_<%d",inint),10000,0,10000,energy_corrected);
              // }
              // if (nInteractions < 5) {
              //   if (perp) obj.FillHistogram(dirname, Form("prompt_gamma_dop_xi_perp_%d",nInteractions),10000,0,10000,energy_corrected);
              //   if (para) obj.FillHistogram(dirname, Form("prompt_gamma_dop_xi_para_%d",nInteractions),10000,0,10000,energy_corrected);
              //   obj.FillHistogram(dirname, Form("gam_sngl_vs_xi_nint%d",nInteractions),180,0,TMath::Pi(),xi,4096,0,4096, energy_corrected);
              // }
              // obj.FillHistogram(dirname, Form("gam_dop_sgl_xi_%d",(int) (xi/TMath::Pi()*5)),10000,0,10000, energy_corrected);
              // obj.FillHistogram(Form("polarization_%s",gates["outgoing"].at(ind_out)->GetName()), Form("gam_sngl_vs_xi_%d",cryID),180,0,TMath::Pi(),xi,nXiPlotEbins,XiPlotLo,XiPlotHi, energy_corrected);
              
              // for (int ri=0; ri<400; ri++) {
              //   double randXi = hit.GetXi(randomBeam(rand_gen));
              //   obj.FillHistogram(dirname, "gam_sngl_vs_xirand",180,0,TMath::Pi(),randXi,4096,0,4096, energy_corrected);
              //   obj.FillHistogram(Form("polarization_%s",gates["outgoing"].at(ind_out)->GetName()), Form("gam_sngl_vs_xirand_%d",cryID),180,0,TMath::Pi(),randXi,4096,0,4096, energy_corrected);
              // }
            } 
          } 
          
          if (prompt_gammas.size() == 1 || prompt_gammas.size() == 2) {
            TGretinaHit &hit = gretina->GetGretinaHit(prompt_gammas[0]);
            double energy_corrected = hit.GetDopplerYta(outgoingBeta, s800->GetYta(), &track);
            bool bin75 = hit.GetTheta()*TMath::RadToDeg() > 65 && hit.GetTheta()*TMath::RadToDeg() < 85;
            bool bin100 = hit.GetTheta()*TMath::RadToDeg() > 90 && hit.GetTheta()*TMath::RadToDeg() < 115;
            if (bin75) {
              obj.FillHistogram(dirname, Form("gam_dop_sgl_bin75_mult%d",prompt_gammas.size()),10000,0,10000, energy_corrected);
            }

            if (bin100) {
              obj.FillHistogram(dirname, Form("gam_dop_sgl_bin100_mult%d",prompt_gammas.size()),10000,0,10000, energy_corrected);
            }
          }
          
          //NNADDBACK
          /*
          int nnSize = gretina->NNAddbackSize();
          std::vector<int> goodNN;
          for (int i=0; i < nnSize; i++) if (gretina->GetNNAddbackHit(i).GetABDepth() < 1) goodNN.push_back(i);
          for (auto i : goodNN){
            TGretinaHit nnhit = gretina->GetNNAddbackHit(i);
            if (nnhit.GetPad() > 0) continue;
            nnhit.ComptonSort();
            int cryID = nnhit.GetCrystalId();
            double nnEnergy_corrected = nnhit.GetDopplerYta(outgoingBeta, s800->GetYta(), &track);
            double nnCore_energy = nnhit.GetCoreEnergy();
            int nInteractions = nnhit.NumberOfInteractions();

            //make sure hits are prompt
            bool tgate = false;
            if (gates["prompt"].size() > 0) tgate = gates["prompt"][0]->IsInside(timeBank29-nnhit.GetTime(), nnCore_energy);
            if (!tgate) continue;

            //only gate on events with certain interaction point numbers
            double xi = nnhit.GetXi(&track);
            if (nInteractions > 1 && nInteractions < 4 && xi != -100) {
              obj.FillHistogram(dirname, "gam_n0_nint>1&nint<4",8192,0,8192, nnEnergy_corrected);
              
              obj.FillHistogram(dirname, "gam_n0_vs_xi",180,0,TMath::Pi(),xi,4096,0,4096, nnEnergy_corrected);
              obj.FillHistogram(Form("polarization_%s",gates["outgoing"].at(ind_out)->GetName()), Form("gam_n0_vs_xi_%d",cryID),180,0,TMath::Pi(),xi,4096,0,4096, nnEnergy_corrected);
              
              for (int ri=0; ri<200; ri++) {
                double randXi = nnhit.GetXi(randomBeam(rand_gen));
                if (randXi == -100) continue;

                obj.FillHistogram(dirname, "gam_n0_vs_xirand",180,0,TMath::Pi(),randXi,4096,0,4096, nnEnergy_corrected);
                obj.FillHistogram(Form("polarization_%s",gates["outgoing"].at(ind_out)->GetName()), Form("gam_n0_vs_xirand_%d",cryID),180,0,TMath::Pi(),randXi,4096,0,4096, nnEnergy_corrected);
              }
            }
          }*/
        }
      }
    }
  }
  if(numobj!=list->GetSize()){
    list->Sort();
  }
}
