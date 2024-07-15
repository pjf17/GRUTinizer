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

bool gates_loaded = false;
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
  for(int j=0;j<16;j++) {
    if(((bits>>j)&0x0001))
      obj.FillHistogram("ungated","trig_bit",20,0,20,j);
  }
  
  //MAKE RAW TOF HISTS
  double raw_obj = s800->GetRawOBJ_MESY();
  double raw_e1 = s800->GetRawE1_MESY();
  double raw_xf = s800->GetRawXF_MESY();
  
  obj.FillHistogram("ungated", "MTOF_OBJE1", 5000, -10000, 0, raw_obj - raw_e1);
  obj.FillHistogram("ungated", "MTOF_XFE1", 4000, -4000, 6000, raw_xf - raw_e1);

  //MAKE INCOMING PID
  double tof_obje1 = s800->GetMTof().GetCorrelatedObjE1(); 
  double tof_xfpe1 = s800->GetMTof().GetCorrelatedXfpE1();
  obj.FillHistogram("ungated", "incoming_pid", 500, -5000, -3000, tof_obje1,
                                               500,  2000,  5000, tof_xfpe1);                                               

  //CRDC PLOTS
  double crdc_1_x = s800->GetCrdc(0).GetDispersiveX();
  double crdc_2_x = s800->GetCrdc(1).GetDispersiveX();
  double crdc_1_y = s800->GetCrdc(0).GetNonDispersiveY();
  double crdc_2_y = s800->GetCrdc(1).GetNonDispersiveY();
  double afp = GetAfp(crdc_1_x, crdc_2_x);
  
  double ylow = -500;
  double yhigh = 500;
  double ybins = 500;
  
  double yslope = GValue::Value("CRDC1_Y_SLOPE");
  if (std::isnan(yslope) || yslope == 0){
    ylow = 0;
    yhigh = 1500;
    ybins = 1500;
  }
  obj.FillHistogram("ungated", "crdc1 X_Y", 600, -300, 300, crdc_1_x, ybins, ylow, yhigh, crdc_1_y);  
  obj.FillHistogram("ungated", "crdc2 X_Y", 600, -300, 300, crdc_2_x, ybins, ylow, yhigh, crdc_2_y);

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
    std::cout<<"OBJ_MTOF_CORR_AFP OR OBJ_MTOF_CORR_XFP NOT SET, SKIPPING OUTGOING GATES"<<std::endl;
  } 
  else {
    tof_obje1_corr = s800->GetMTofObjE1();
    CheckGates(gates["outgoing"],outgoing_passed,tof_obje1_corr,ic_energy);
  }

  std::string dirname  = "";
  for (auto ind_out : outgoing_passed){
    dirname = Form("%s_gated", gates["outgoing"].at(ind_out)->GetName());
    obj.FillHistogram(dirname, "incoming_pid", 500, -5000, -3000, tof_obje1,
                                               500,  2000,  5000, tof_xfpe1);   
  }

  //---------------------------------------------------------------
  //INCOMING
  for (auto ind_in: incoming_passed){
    dirname = Form("%s_gated", gates["incoming"].at(ind_in)->GetName());
    obj.FillHistogram(dirname, "outgoing_pid_uncorrected", 2000, -5000, -3000, tof_obje1, 2048, 0, 4096, ic_energy);
    obj.FillHistogram(dirname, "outgoing_pid", 2000, -5000, -3000, tof_obje1_corr, 2048, 0, 4096, ic_energy);                                                
    
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
            double energy_corrected = hit.GetDopplerYta(s800->AdjustedBeta(outgoingBeta), s800->GetYta(), &track);
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
              // ECORE THETA INTERACTION POINT GATES
              for (int ip=0; ip < nInteractions; ip++)
                obj.FillHistogram(dirname, "IP_Ecore_FEP_vs_theta",120,0.6,2.1,hit.GetTheta(ip),4096,0,4096,core_energy);
              
              for (auto ipgate : gates["EnergyTheta"]){
                std::string ipgname = std::string(ipgate->GetName());
                std::vector<int> IPCorepass;
                for (int ip=0; ip < nInteractions; ip++){
                  if (ipgate->IsInside(hit.GetTheta(ip),core_energy) && hit.GetSegmentEng(ip) > 100){
                    IPCorepass.push_back(ip);
                  }
                }

                if (IPCorepass.size() > 0){
                  obj.FillHistogram(dirname, Form("IP_%s_npass_vs_totalPoints",ipgname.c_str()),11,1,12,nInteractions,9,1,10,(int) IPCorepass.size());
                  double E_pass_dop = hit.GetDopplerYta(s800->AdjustedBeta(outgoingBeta), s800->GetYta(), &track, IPCorepass[0]);
                  
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
            } 
          }

          //count multiplicity
          obj.FillHistogram(dirname, "gamma_multi_vs_prompt_gamma_multi", 20,0,20,nPromptGamma,20,0,20,gSize);
          if (total_corrected_energy > 0) obj.FillHistogram(dirname,"total_corrected_energy_vs_prompt_multi",20,0,20,nPromptGamma,2048,0,8192,total_corrected_energy);
          if (total_core_energy > 0) obj.FillHistogram(dirname,"total_core_energy_vs_prompt_multi",20,0,20,nPromptGamma,2048,0,8192,total_core_energy);

          // //NNADDBACK
          // int nABpromptGamma = 0;
          // int nTotalABgamma = 0;
          // double tot_energy_corrected_AB = 0;
          // double tot_energy_core_AB = 0;

          // //loop over multiplicity
          // for (int n=0; n<4; n++){
          //   //loop over hits for each multiplicity spectrum
          //   int nnSize = gretina->NNAddbackSize(n);
          //   nTotalABgamma += nnSize;
          //   for (int i=0; i < nnSize; i++){
          //     //get hit and hit data 
          //     TGretinaHit nnhit = gretina->GetNNAddbackHit(n,i);
          //     // int cryID = nnhit.GetCrystalId();
          //     // int ringNum = nnhit.GetRingNumber();
          //     double nnEnergy_corrected = nnhit.GetDopplerYta(s800->AdjustedBeta(outgoingBeta), s800->GetYta(), &track);
          //     double nnCore_energy = nnhit.GetCoreEnergy();
          //     // double theta = nnhit.GetThetaDeg();
          //     // double phi = nnhit.GetPhiDeg();
          //     int nInteractions = nnhit.NumberOfInteractions();
              
          //     bool tgate = false;
          //     if (gates["prompt"].size() > 0) tgate = gates["prompt"][0]->IsInside(timeBank29-nnhit.GetTime(), nnCore_energy);

          //     //make sure hits are prompt
          //     if (tgate){
          //       nABpromptGamma++;
          //       tot_energy_corrected_AB += nnEnergy_corrected;
          //       tot_energy_core_AB += nnCore_energy;
          //       //exclude the ng spectrum (n==3)
          //       // obj.FillHistogram(dirname,"crystal-map",180,0,180,theta,360,0,360,phi);
          //       char *multiplicity = Form("%d",n);
          //       if (n == 3) multiplicity = Form("g");
          //       obj.FillHistogram(dirname, Form("gamma_corrected_n%s_prompt",multiplicity), 8192,0,8192, nnEnergy_corrected);
                
          //       if (n < 3) {
          //         obj.FillHistogram(dirname, "gamma_corrected_addback_prompt", 8192,0,8192, nnEnergy_corrected);
          //         // GAMMA GAMMA CORRELATION
          //         for (int j=0; j < nnSize; j++){
          //           if (i==j) continue;
          //           TGretinaHit nnhit2 = gretina->GetNNAddbackHit(n,j);
          //           double nnEnergy_corrected2 = nnhit2.GetDopplerYta(s800->AdjustedBeta(outgoingBeta), s800->GetYta(), &track);

          //           bool tgate2 = false;
          //           if (gates["prompt"].size() > 0) tgate2 = gates["prompt"][0]->IsInside(timeBank29-nnhit2.GetTime(), nnhit2.GetCoreEnergy());
                    
          //           if (tgate2){
          //             obj.FillHistogram(dirname, "gamma_gamma", 2048,0,4096, nnEnergy_corrected2, 2048,0,4096, nnEnergy_corrected);
          //             // if (n==0) obj.FillHistogram(dirname, "gamma_gamma_n0", 1024,0,4096, nnEnergy_corrected2, 1024,0,4096, nnEnergy_corrected);
          //             // if (n==1) obj.FillHistogram(dirname, "gamma_gamma_n1", 1024,0,4096, nnEnergy_corrected2, 1024,0,4096, nnEnergy_corrected);
          //           }
          //         }
          //       }
          //     }
          //   }
          // }

          // //count ab multiplicity
          // obj.FillHistogram(dirname, "gamma_multi_vs_prompt_gamma_multi_AB", 20,0,20,nABpromptGamma,20,0,20,nTotalABgamma);
          // if (tot_energy_corrected_AB > 0) obj.FillHistogram(dirname, "total_corrected_energy_vs_prompt_multi_AB",20,0,20,nPromptGamma,2048,0,8192,tot_energy_corrected_AB);
          // if (tot_energy_core_AB > 0) obj.FillHistogram(dirname, "total_core_energy_vs_prompt_multi_AB",20,0,20,nPromptGamma,2048,0,8192,tot_energy_core_AB);
        } 
      }
    }
  }
  if(numobj!=list->GetSize()){
    list->Sort();
  }
}
