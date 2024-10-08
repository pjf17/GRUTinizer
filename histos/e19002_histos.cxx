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
#include <TRandom.h>
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

double dopplerbroad(double *x, double *p){
  double sum = pow(p[1]*sin(x[0])/(1-p[1]*cos(x[0]))*p[2],2) + 
                pow( (-p[1]+cos(x[0]))/((1-p[1]*p[1])*(1-p[1]*cos(x[0])))*p[3] ,2) +
                pow(p[4]*cos(x[0]),2);

  return p[0]*p[5]*TMath::Sqrt(sum);
}

double GetAfp(double crdc_1_x,double  crdc_2_x){
  return TMath::ATan( (crdc_2_x - crdc_1_x)/1073.0 );
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

double thetaCorrelation(double xi, double alpha, double theta, double beta, double E){
  alpha *= TMath::DegToRad();
  theta *= TMath::DegToRad();
  double cosThetaPrime = TMath::Cos(alpha)*TMath::Cos(theta) + TMath::Sin(alpha)*TMath::Sin(theta)*TMath::Cos(xi);
  return E*(1 - beta*cosThetaPrime)/(1 - beta*TMath::Cos(theta));
}

double lastPointPenalty(double x){
  double val = TMath::TanH((x/100-1.1)/0.5)*exp(-3*(x/100-1.1)+1);
  if (val < 0) val = 0;
  return val + 1;
}

void comptonSort(const TGretinaHit &ghit, int &FP, int &SP) {
  double FOM = 1e10;
  int N = ghit.NumberOfInteractions();
  double ipFactor = 308.76*std::pow(N,-0.542);
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
      double ffom = std::pow(std::abs(1-x),2.0/3)*ghit.GetAlpha(fp,sp)*kn*std::pow(E/E1,3)*ghit.GetLocalPosition(fp).Z(); //*std::pow(E/E2,1.0/3);
      ffom *= std::pow(E/E2,2)*TMath::Sqrt(ghit.GetLocalPosition(sp).Z());
      // ffom /= TMath::Sqrt(abs(E1-125)*abs(E2-125))/E;
      ffom *= lastPointPenalty(E1)*lastPointPenalty(E2);

      if (ffom < FOM) {
        FOM = ffom;
        FP = fp;
        SP = sp;
      }
      // if (fp == sp) continue;
      // double E2 = ghit.GetSegmentEng(sp);
      // double cosp = TMath::Cos(ghit.GetScatterAngle(fp,sp));
      // double x = er + cosp;
      // double kn = pow((E - E1)/E,2)*((E-E1)/E + E/(E-E1) - pow(TMath::Sin(ghit.GetScatterAngle(fp,sp)),2) );
      // // double kn = pow((E - E1)/E,2)*((E-E1)/E + E/(E-E1) - pow(TMath::Sin(ghit.GetScatterAngle(fp,sp)),2) );
      // // double ffom = std::pow(std::abs(1-x),2.0/3)*ghit.GetAlpha(fp,sp)*kn*E/E1;
      // double ffom = std::pow(std::abs(1-x),2.0/3)*ghit.GetAlpha(fp,sp)*std::pow(E/E1,2)*std::pow(E/E2,1.0/3)*kn;
      // // double ffom = std::pow(std::abs(1-x),2.0/3)*ghit.GetAlpha(fp,sp)*std::pow(E/E1,2)*kn;

      // if (ffom < FOM) {
      //   FOM = ffom;
      //   FP = fp;
      //   SP = sp;
      // }
    }
  }
  return;
}

void comptonSortTest(const TGretinaHit &ghit, int &FP, int &SP) {
  double FOM = 1e10;
  FP = 0; SP = 1;
  int N = ghit.NumberOfInteractions();
  double E = ghit.GetCoreEnergy();
  //scale the interaction points so they match the core energy
  double scaleFactor = 0;
  for (int i=0; i < N; i++) scaleFactor += ghit.GetSegmentEng(i);
  scaleFactor = E/scaleFactor;

  for (int fp=0; fp < N; fp++){
    double E1 = ghit.GetSegmentEng(fp)*scaleFactor;
    double er = 511.0/E * E1/(E - E1);

    for (int sp=0; sp < N; sp++){
      if (fp == sp) continue;
      
      double cosp = TMath::Cos(ghit.GetScatterAngle(fp,sp));
      double E2 = ghit.GetSegmentEng(sp)*scaleFactor;
      double x = er + cosp;
      double kn = pow((E - E1)/E,2)*((E-E1)/E + E/(E-E1) - pow(TMath::Sin(ghit.GetScatterAngle(fp,sp)),2) );
      double ffom = std::pow(std::abs(1-x),2.0/3) + std::pow(E1/E - 1/(1+511/E/(1-cosp)),2);
      ffom *= lastPointPenalty(E1)*lastPointPenalty(E2)*std::pow(E/E1,2)*ghit.GetLocalPosition(fp).Z()*ghit.GetAlpha(fp,sp)*kn;
      ffom *= std::pow(E/E2,2) * ghit.GetLocalPosition(sp).Z();

      if (ffom < FOM) {
        FOM = ffom;
        FP = fp;
        SP = sp;
      }
    }
  }
  return;
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
    obj.FillHistogram(dirname, "outgoing_pid_uncorrected", 600, -2800, -2200, tof_obje1, 2048, 0, 4096, ic_energy);
    obj.FillHistogram(dirname, "outgoing_pid", 600, -2800, -2200, tof_obje1_corr, 2048, 0, 4096, ic_energy);                                                
    
    //CRDC DE CORRECTION
    // double ic_ave = s800->GetIonChamber().GetAve();
    // for (auto ind_iso : isoline_passed){
    //   dirname = Form("%s_%s_gated", gates["incoming"].at(ind_in)->GetName(), isoline_gates.at(ind_iso)->GetName());
    //   obj.FillHistogram(dirname, "crdc1x_dE", 2048, 0, 4096, ic_energy, 600, -300, 300, crdc_1_x);
    //   obj.FillHistogram(dirname, "crdc1x_Ave", 2048, 0, 4096, ic_ave, 600, -300, 300, crdc_1_x);
    //   obj.FillHistogram(dirname, "crdc1y_dE", 2048, 0, 4096, ic_energy, 600, -300, 300, crdc_1_y);
    //   obj.FillHistogram(dirname, "crdc1y_Ave", 2048, 0, 4096, ic_ave, 600, -300, 300, crdc_1_y);
    // }

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
      obj.FillHistogram(dirname, "obje1_crdc1x", 3000, -6000, 0, tof_obje1,600, -300, 300, crdc_1_x);                                     
      obj.FillHistogram(dirname, "corrobje1_afp", 3000, -6000, 0, tof_obje1_corr,1000, -0.1, 0.1, afp);
      obj.FillHistogram(dirname, "obje1_afp", 3000, -6000, 0, tof_obje1,1000, -0.1, 0.1, afp);
      obj.FillHistogram(dirname, "corrobje1_tofxfpobj", 3000, -6000, 0, tof_obje1_corr,xocor_nbins, xocor_lowbin, xocor_highbin, xfp_obj);
      obj.FillHistogram(dirname, "obje1_tofxfpobj", 3000, -6000, 0, tof_obje1,xocor_nbins, xocor_lowbin, xocor_highbin, xfp_obj);

      //GAMMA RAY ANALYSIS
      if (gretina){
        TVector3 track = s800->Track();
        if (bank29){
          double timeBank29 = bank29->Timestamp();

          //get Beta
          double outgoingBeta = GValue::Value(Form("BETA_%s",gates["outgoing"].at(ind_out)->GetName()));
          if (std::isnan(outgoingBeta)) outgoingBeta = GValue::Value("BETA");
          
          //SINGLES
          int gSize = gretina->Size();
          int nPromptGamma = 0;
          double total_corrected_energy = 0;
          double total_core_energy = 0;
          for (int i=0; i < gSize; i++){
            TGretinaHit &hit = gretina->GetGretinaHit(i);
            hit.ComptonSort();
            double energy_corrected = hit.GetDopplerYta(outgoingBeta, s800->GetYta(), &track);
            double energy = hit.GetDoppler(outgoingBeta);
            double core_energy = hit.GetCoreEnergy();
            double theta = hit.GetTheta();
            double phi = hit.GetPhi();
            int cryID = hit.GetCrystalId();
            int nInteractions = hit.NumberOfInteractions();
            double xi = hit.GetXi(&track);
            
            //PROMPT GATE
            bool tgate = false;
            if (gates["prompt"].size() > 0) tgate = gates["prompt"][0]->IsInside(timeBank29-hit.GetTime(), core_energy);
            else std::cout<<"NO PROMPT GATE LOADED\n";

            if (tgate){
              nPromptGamma++;
              total_corrected_energy += energy_corrected;
              total_core_energy += core_energy;

              obj.FillHistogram(dirname, "gretina_theta_vs_phi",360,0,360,phi*TMath::RadToDeg(),180,0,180,theta*TMath::RadToDeg());
              // obj.FillHistogram(dirname, Form("gretina_theta_vs_phi_rn%02d",hit.GetRingNumber()),360,0,360,phi*TMath::RadToDeg(),180,0,180,theta*TMath::RadToDeg());
              
              //summary spectra
              obj.FillHistogram(dirname, "summary_gam_dop_sgl_prompt",48,0,48,detMapRing[cryID],4096,0,4096, energy_corrected);
              obj.FillHistogram(dirname, "summary_core_energy_prompt",48,0,48,detMapRing[cryID],4096,0,4096, core_energy);
              // if ( !((1002 < energy_corrected && energy_corrected < 1034) || (1064 < energy_corrected && energy_corrected < 1094) ||
              //      (712 < energy_corrected && energy_corrected < 764)) )
              
              obj.FillHistogram(dirname, "gam_dop_sgl_prompt",8192,0,8192, energy_corrected);
              obj.FillHistogram(dirname, "gam_dop_sgl_prompt_vs_nInteraction",10,0,10,nInteractions,1024,0,4096, energy_corrected);
              obj.FillHistogram(dirname, "gam_dop_sgl_prompt_vs_nInteraction",10,0,10,nInteractions,1024,0,4096, energy_corrected);
              // obj.FillHistogram(dirname, Form("gam_dop_sgl_prompt_rn%02d",hit.GetRingNumber()),4096,0,4096, energy_corrected);

              // for (int j=0; j < gSize; j++) {
              //   if (i==j) continue;
              //   TGretinaHit &hit2 = gretina->GetGretinaHit(j);
              //   if ( !gates["prompt"][0]->IsInside(timeBank29-hit2.GetTime(), hit2.GetCoreEnergy()) ) continue;
              //   double energy_corrected2 = hit2.GetDopplerYta(s800->AdjustedBeta(outgoingBeta), s800->GetYta(), &track);
              //   obj.FillHistogram(dirname, "gamma_gamma",2048,0,4096,energy_corrected2,2048,0,4096,energy_corrected);
              // }

              if (nInteractions > 1){
                // int myFP, mySP;
                // comptonSortTest(hit,myFP,mySP);
                double myxi = hit.GetXi(&track);
                // double new_energy = hit.GetDopplerYta(s800->AdjustedBeta(outgoingBeta), s800->GetYta(),&track,myFP);
                obj.FillHistogram(dirname, "gam_sngl_vs_xi",360,0,TMath::TwoPi(),myxi,4096,0,4096, energy_corrected);
                // if (nInteractions < 4) obj.FillHistogram(dirname, "new_gam_sngl_vs_new_xi<4intp",360,0,TMath::TwoPi(),myxi,4096,0,4096, new_energy);
                // obj.FillHistogram(dirname, "new_gam_sngl",4096,0,4096, new_energy);

                // if (myFP != 0) {
                //   obj.FillHistogram(dirname, "new!=main_new",4096,0,4096, new_energy);
                //   obj.FillHistogram(dirname, "new!=main_main",4096,0,4096, energy_corrected);
                // }

                // double scaleFactor = 0;
                // for (int ipx=0; ipx < nInteractions; ipx++) scaleFactor += hit.GetSegmentEng(ipx);
                // scaleFactor = core_energy/scaleFactor;
                // obj.FillHistogram(dirname, "new_gam_sngl_vs_firstIPEratio",4096,0,4096, energy_corrected,100,0,1,hit.GetSegmentEng(myFP)/core_energy*scaleFactor);
                // for (int ipx=0; ipx < nInteractions; ipx++){ 
                //   obj.FillHistogram(dirname, "new_gam_sngl_vs_IPEratio",4096,0,4096, energy_corrected,100,0,1,hit.GetSegmentEng(ipx)/core_energy*scaleFactor);
                // }

                // if (theta*TMath::RadToDeg() >= 55 && theta*TMath::RadToDeg() <= 100)
                //   obj.FillHistogram(dirname, "new_gam_sngl_vs_new_xi_theta_gate",360,0,TMath::TwoPi(),myxi,4096,0,4096, new_energy);
              } 
              // else {
              //   obj.FillHistogram(dirname, "new_gam_sngl",4096,0,4096, energy_corrected);
              // }

              /*
              //Ecore theta plot for all IPs
              for (int ip=0; ip < nInteractions; ip++)
                obj.FillHistogram(dirname, "IP_Ecore_vs_theta",60,0.6,2.1,hit.GetTheta(ip),2048,0,4096,core_energy);
              
              TF1 *fDOP = new TF1("fdop","[0]*TMath::Sqrt(1-[1]*[1])/(1 - [1]*TMath::Cos(x))",0,TMath::Pi());
              TF1 *fRES = new TF1("fresolution",dopplerbroad,0,TMath::Pi(),6);
              fDOP->SetParameter(1,outgoingBeta);
              fRES->SetParameter(1,outgoingBeta);
              fRES->SetParameter(2,GValue::Value("DOPP_BROAD_DTHETA"));
              fRES->SetParameter(3,GValue::Value("DOPP_BROAD_DBETA"));
              fRES->SetParameter(4,GValue::Value("DOPP_BROAD_BEAMSPOT"));
              fRES->SetParameter(5,GValue::Value("DOPP_BROAD_SCALE"));

              // double CentroidRestEnergy = std::floor(core_energy*(1-BETA*TMath::Cos(0.65))/TMath::Sqrt(1-BETA*BETA));
              // double ECentMax = std::ceil(core_energy*(1-BETA*TMath::Cos(2.05))/TMath::Sqrt(1-BETA*BETA));
              // double EMAXCONSIDERED = 4000;
              double H_elo = GValue::Value("DOPP_BROAD_SCAN_LO");;
              double H_ehi = GValue::Value("DOPP_BROAD_SCAN_HI");;
              int H_nbins = int(H_ehi-H_elo);
              double CentroidRestEnergy = H_elo;
              double CRE_step = GValue::Value("DOPP_BROAD_SCAN_STEP");

              // fDOP->SetParameter(0,1017);
              // fRES->SetParameter(0,1017);
              std::vector<int> passedIP;   
              // //stats test
              // for (int ip=0; ip < nInteractions; ip++){
              //   if (checkEnergyTheta(fDOP,fRES,core_energy,hit.GetTheta(ip))){
              //     passedIP.push_back(ip);
              //   }
              // }
              // if (!passedIP.empty()) {
              //   obj.FillHistogram(dirname, "CentroidEnergy1017_DopReconPassed",300,800,1100,hit.GetDopplerYta(s800->AdjustedBeta(outgoingBeta), s800->GetYta(), &track, passedIP[0]));
              // }
              // passedIP.clear();
              // while (CentroidRestEnergy < ECentMax && CentroidRestEnergy < EMAXCONSIDERED){
              while (CentroidRestEnergy < H_ehi){
                CentroidRestEnergy += CRE_step;
                fDOP->SetParameter(0,CentroidRestEnergy);
                fRES->SetParameter(0,CentroidRestEnergy);
                
                for (int ip=0; ip < nInteractions; ip++){
                  if (checkEnergyTheta(fDOP,fRES,core_energy,hit.GetTheta(ip))){
                    passedIP.push_back(ip);
                  }
                }

                obj.FillHistogram(dirname, "CentroidEnergy_vs_DopReconTot",H_nbins,H_elo,H_ehi,energy_corrected,H_nbins,H_elo,H_ehi,CentroidRestEnergy);

                if (passedIP.empty()) {
                  obj.FillHistogram(dirname, "CentroidEnergy_vs_DopReconNotPassed",H_nbins,H_elo,H_ehi,energy_corrected,H_nbins,H_elo,H_ehi,CentroidRestEnergy);
                  // for (int j=0; j < gSize; j++) {
                  //   if (i==j) continue;
                  //   TGretinaHit &hit2 = gretina->GetGretinaHit(j);
                  //   if ( !gates["prompt"][0]->IsInside(timeBank29-hit2.GetTime(), hit2.GetCoreEnergy()) ) continue;
                  //   double energy_corrected2 = hit2.GetDopplerYta(s800->AdjustedBeta(outgoingBeta), s800->GetYta(), &track);
                  //   obj.FillHistogram(dirname, "CentroidEnergy_vs_CoincidentGammaNotPassed",2048,0,4096,energy_corrected2,H_nbins,H_elo,H_ehi,CentroidRestEnergy);
                  // }
                } 
                
                else {
                  double E_pass_dop = hit.GetDopplerYta(s800->AdjustedBeta(outgoingBeta), s800->GetYta(), &track, passedIP[0]);
                  obj.FillHistogram(dirname, "CentroidEnergy_vs_DopReconPassed",H_nbins,H_elo,H_ehi,E_pass_dop,H_nbins,H_elo,H_ehi,CentroidRestEnergy);
                  //Coincident gamma rays
                  // for (int j=0; j < gSize; j++) {
                  //   if (i==j) continue;
                  //   TGretinaHit &hit2 = gretina->GetGretinaHit(j);
                  //   if ( !gates["prompt"][0]->IsInside(timeBank29-hit2.GetTime(), hit2.GetCoreEnergy()) ) continue;
                  //   double energy_corrected2 = hit2.GetDopplerYta(s800->AdjustedBeta(outgoingBeta), s800->GetYta(), &track);
                  //   obj.FillHistogram(dirname, "CentroidEnergy_vs_CoincidentGamma",2048,0,4096,energy_corrected2,H_nbins,H_elo,H_ehi,CentroidRestEnergy);
                  // }

                  //NEW SECOND INTERACTION POINT
                  if (nInteractions > 1){
                    // double IPxi_main = xi;
                    // double IPxi_order = xi;
                    // int secondPoint_compt = 1;
                    int myFP, mySP;
                    comptonSortGate(hit,passedIP,myFP,mySP);
                    double newxi = hit.GetXi(&track,myFP,mySP);
                    double E_new_pass_dop = hit.GetDopplerYta(s800->AdjustedBeta(outgoingBeta), s800->GetYta(), &track, myFP);
                    obj.FillHistogram(dirname, "CentroidEnergy_vs_test_fp_energy",H_nbins,H_elo,H_ehi,E_new_pass_dop,H_nbins,H_elo,H_ehi,CentroidRestEnergy);
                    obj.FillHistogram(dirname, "CentroidEnergy_vs_newxi",360,0,TMath::TwoPi(),newxi,H_nbins,H_elo,H_ehi,CentroidRestEnergy);

                    // double IPxi_compt = hit.GetXi(&track,myFP,mySP);
                    // if (passedIP[0] != 0) IPxi_main = hit.GetXi(&track,passedIP[0],0);
                    // if (passedIP[0]+1 < nInteractions) IPxi_order = hit.GetXi(&track,passedIP[0],passedIP[0]+1);
                    // else if (passedIP[0] != 0) IPxi_order = IPxi_main;

                    // obj.FillHistogram(dirname, "CentroidEnergy_vs_Xi_main",180,0,TMath::TwoPi(),IPxi_main,H_nbins,H_elo,H_ehi,CentroidRestEnergy);
                    // obj.FillHistogram(dirname, "CentroidEnergy_vs_Xi_order",360,0,TMath::TwoPi(),IPxi_order,H_nbins,H_elo,H_ehi,CentroidRestEnergy);
                    // if (hit.GetTheta(passedIP[0])*TMath::RadToDeg() > 50 && hit.GetTheta(passedIP[0])*TMath::RadToDeg() < 100) 
                      // obj.FillHistogram(dirname, "CentroidEnergy_vs_Xi_order_theta_gate",360,0,TMath::TwoPi(),IPxi_order,H_nbins,H_elo,H_ehi,CentroidRestEnergy);
                    // obj.FillHistogram(dirname, "CentroidEnergy_vs_Xi_compt",180,0,TMath::TwoPi(),IPxi_compt,H_nbins,H_elo,H_ehi,CentroidRestEnergy);
                    // obj.FillHistogram(dirname, "CentroidEnergy_vs_Xi_Coarse",360,0,TMath::TwoPi(),IPxi, int(EMAXCONSIDERED),0,EMAXCONSIDERED,CentroidRestEnergy);
                  }
                }
                passedIP.clear();
              }
              &&&&&&*/

              if (nInteractions > 1){
                double plXi = xi;
                double xiMax = TMath::TwoPi();
                obj.FillHistogram(dirname, "gam_dop_sgl_IP>1_prompt_vs_xi",360,0,xiMax,plXi, 1500,0,3000, energy_corrected);
                obj.FillHistogram(dirname, "gam_dop_sgl_IP>1_prompt_vs_theta",180,0,180,theta*TMath::RadToDeg(), 2048,0,4096, energy_corrected);
              }
            } 
          }

          //count multiplicity
          obj.FillHistogram(dirname, "gamma_multi_vs_prompt_gamma_multi", 20,0,20,nPromptGamma,20,0,20,gSize);
          if (total_corrected_energy > 0) obj.FillHistogram(dirname,"total_corrected_energy_vs_prompt_multi",20,0,20,nPromptGamma,2048,0,8192,total_corrected_energy);
          if (total_core_energy > 0) obj.FillHistogram(dirname,"total_core_energy_vs_prompt_multi",20,0,20,nPromptGamma,2048,0,8192,total_core_energy);

          //NNADDBACK
          //loop over multiplicity
          for (int n=0; n<4; n++){
            //loop over hits for each multiplicity spectrum
            int nnSize = gretina->NNAddbackSize(n);
            for (int i=0; i < nnSize; i++){
              //get hit and hit data 
              TGretinaHit nnhit = gretina->GetNNAddbackHit(n,i);
              // nnhit.ComptonSort();
              // int cryID = nnhit.GetCrystalId();
              // int ringNum = nnhit.GetRingNumber();
              double nnEnergy_corrected = nnhit.GetDopplerYta(outgoingBeta, s800->GetYta(), &track);
              double nnCore_energy = nnhit.GetCoreEnergy();
              // double theta = nnhit.GetThetaDeg();
              // double phi = nnhit.GetPhiDeg();
              int nInteractions = nnhit.NumberOfInteractions();
              
              //make sure hits are prompt
              bool tgate = false;
              if (gates["prompt"].size() > 0) tgate = gates["prompt"][0]->IsInside(timeBank29-nnhit.GetTime(), nnCore_energy);
              if (!tgate) continue;
              
              char *multiplicity = Form("%d",n);
              if (n == 3) multiplicity = Form("g");
              obj.FillHistogram(dirname, Form("gamma_corrected_n%s_prompt",multiplicity), 8192,0,8192, nnEnergy_corrected);
              
              if (n < 3) {
                obj.FillHistogram(dirname, "gamma_corrected_addback_prompt", 8192,0,8192, nnEnergy_corrected);
                // GAMMA GAMMA CORRELATION
                for (int j=0; j < nnSize; j++){
                  if (i==j) continue;
                  TGretinaHit nnhit2 = gretina->GetNNAddbackHit(n,j);
                  double nnEnergy_corrected2 = nnhit2.GetDopplerYta(outgoingBeta, s800->GetYta(), &track);

                  bool tgate2 = false;
                  if (gates["prompt"].size() > 0) tgate2 = gates["prompt"][0]->IsInside(timeBank29-nnhit2.GetTime(), nnhit2.GetCoreEnergy());
                  if (!tgate2) continue;
                  
                  obj.FillHistogram(dirname, "gamma_gamma", 2048,0,8192, nnEnergy_corrected2, 2048,0,8192, nnEnergy_corrected);
                }
              }
              
            }
          }
        } 
      }
      
    }
  }
  if(numobj!=list->GetSize()){
    list->Sort();
  }
}
