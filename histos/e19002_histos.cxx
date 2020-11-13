
#include "TRuntimeObjects.h"

#include <iostream>
#include <map>
#include <cstdio>

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

std::vector<GCutG*> incoming_gates = {};
std::vector<GCutG*> outgoing_gates = {};
std::vector<GCutG*> isoline_gates = {};
GCutG *prompt_timing_gate=0;
GCutG *afp_gate=0;
int gates_loaded=0;

double GetAfp(double crdc_1_x,double  crdc_2_x){
  return TMath::ATan( (crdc_2_x - crdc_1_x)/1073.0 );
}

//Get the corresponding OBJE1 TOF based on whether there are AFP and XFP corrections
bool GetGoodMTOFObjE1(TS800 *s800, double &obje1){
  bool flag = true; 
  if(std::isnan(GValue::Value("OBJ_MTOF_CORR_AFP")) || 
    std::isnan(GValue::Value("OBJ_MTOF_CORR_XFP")) ) {
      obje1 = -123;
      flag = false;
  } else {
    obje1 = s800->GetMTofObjE1();
  }
  return flag;
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

void CheckGates(std::vector<unsigned short> &incoming_passed, std::vector<unsigned short> &outgoing_passed, std::vector<unsigned short> &isoline_passed);

void LoadGates(TRuntimeObjects &obj){
  TList *gates = &(obj.GetGates());
  TIter iter(gates);
  std::cout << "loading gates:" <<std::endl;
  while(TObject *obj = iter.Next()) {
    GCutG *gate = (GCutG*)obj;
    std::string tag = gate->GetTag();
    if(!tag.compare("incoming")) {
      incoming_gates.push_back(gate);
      std::cout << "\t incoming: << " << gate->GetName() << std::endl;
    } else if(!tag.compare("outgoing")) {
      outgoing_gates.push_back(gate);
      std::cout << "\t outgoing: << " << gate->GetName() << std::endl;
    } else if(!tag.compare("isoline")){
      isoline_gates.push_back(gate);
      std::cout << "\t isoline: << " << gate->GetName() << std::endl;
    } else if(!tag.compare("prompt")){
      prompt_timing_gate = new GCutG(*gate);
      std::cout << "\t prompt_timing_gate: << " << gate->GetName() << std::endl;
    } else if(!tag.compare("afp")){
      afp_gate = new GCutG(*gate);
      std::cout << "\t afp_gate: << " << gate->GetName() << std::endl;
    } else {
      std::cout << "\t unknown: << " << gate->GetName() << std::endl;
    }
    gates_loaded++;
  }
  std::cout << "outgoing size: " << outgoing_gates.size() << std::endl;
}

void CheckGates(TS800 *s800, std::vector<unsigned short> &incoming_passed, std::vector<unsigned short> &outgoing_passed, std::vector<unsigned short> &isoline_passed){
  for(unsigned short i=0;i<incoming_gates.size();i++) {
    if(incoming_gates.at(i)->IsInside(s800->GetMTof().GetCorrelatedObjE1(), 
                                      s800->GetMTof().GetCorrelatedXfpE1())){
      incoming_passed.push_back(i);
    }
  }

  double corr_obje1 = s800->GetMTofObjE1();
  double ic = GetGoodICE(s800);
  
  for(unsigned int i=0;i<outgoing_gates.size();i++) {
    if(outgoing_gates.at(i)->IsInside(corr_obje1,ic)){
      outgoing_passed.push_back(i);
    }
  }

  for(unsigned int i=0;i<isoline_gates.size();i++) {
    if(isoline_gates.at(i)->IsInside(corr_obje1,ic)){
      isoline_passed.push_back(i);
    }
  }
}

// extern "C" is needed to prevent name mangling.
// The function signature must be exactly as shown here,
//   or else bad things will happen.
extern "C"
void MakeHistograms(TRuntimeObjects& obj) {
//  InitMap();
  TGretina *gretina = obj.GetDetector<TGretina>();
  TBank29  *bank29  = obj.GetDetector<TBank29>();
  TS800    *s800    = obj.GetDetector<TS800>();
  TList    *list    = &(obj.GetObjects());
  int numobj = list->GetSize();


  TList *gates = &(obj.GetGates());
  if (!s800){
    return;
  }

  if(gates_loaded!=gates->GetSize()) {
    LoadGates(obj);
  }

  //Use this spectrum for the time-energy cut for GRETINA
  if(bank29 && gretina) {
    for(unsigned int i=0;i<gretina->Size();i++) {
      //Time-energy cut
      TGretinaHit &hit = gretina->GetGretinaHit(i);
      obj.FillHistogram("Bank29","Gretina_dop_t0_Bank29_time",
          600,-600,600,bank29->Timestamp()-hit.GetTime(),
          2500,0,10000, hit.GetDoppler(GValue::Value("BETA")));
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
  
  obj.FillHistogram("ungated", "MTOF_OBJE1", 3000, -6000, 0, raw_obj - raw_e1);
  obj.FillHistogram("ungated", "MTOF_XFE1", 3000, -2000, 4000, raw_xf - raw_e1);

  //MAKE INCOMING PID
  double tof_obje1 = s800->GetMTof().GetCorrelatedObjE1(); 
  double tof_xfpe1 = s800->GetMTof().GetCorrelatedXfpE1();
  obj.FillHistogram("ungated", "incoming_pid", 2500, -5000, 0, tof_obje1,
                                               2000, 0, 4000, tof_xfpe1);                                              

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
    yhigh = 1200;
    ybins = 1200;
  }
  obj.FillHistogram("ungated", "crdc1 X_Y", 600, -300, 300, crdc_1_x, ybins, ylow, yhigh, crdc_1_y);  
  obj.FillHistogram("ungated", "crdc2 X_Y", 600, -300, 300, crdc_2_x, ybins, ylow, yhigh, crdc_2_y);

  //GET TOFs FOR PID AND CORRELATION PLOTS
  double tof_obje1_corr;
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
  if (GetGoodMTOFObjE1(s800, tof_obje1_corr)){    

    //OUTGOING PID
    obj.FillHistogram("ungated", "outgoing_pid", 2000, -4000, 0, tof_obje1_corr,
                                                2048, 0, 4096, ic_energy);
    obj.FillHistogram("ungated", "outgoing_pid_uncorrected", 2000, -4000, 0, tof_obje1,
                                                2048, 0, 4096, ic_energy);  
  } 
  
  //---------------------------------------------------------------
  //GATED

  std::vector<unsigned short> incoming_passed;
  std::vector<unsigned short> outgoing_passed;
  std::vector<unsigned short> isoline_passed;
  CheckGates(s800, incoming_passed, outgoing_passed, isoline_passed);
  std::string dirname  = "";

  for (auto ind_out : outgoing_passed){
    dirname = Form("%s_gated", outgoing_gates.at(ind_out)->GetName());
    obj.FillHistogram(dirname, "incoming_pid", 2500, -5000, 0, tof_obje1,
                                               2000, 0, 4000, tof_xfpe1);
  }

  for (auto ind_in : incoming_passed){
    dirname = Form("%s_gated", incoming_gates.at(ind_in)->GetName());
    double ic_ave = s800->GetIonChamber().GetAve();
    if (GetGoodMTOFObjE1(s800, tof_obje1_corr)){
      obj.FillHistogram(dirname, "outgoing_pid", 4000, -4000, 0, tof_obje1_corr, 2048, 0, 4096, ic_energy);                                                
    }
    
    for (auto ind_iso : isoline_passed){
      dirname = Form("%s_%s_gated", incoming_gates.at(ind_in)->GetName(), isoline_gates.at(ind_iso)->GetName());
      obj.FillHistogram(dirname, "crdc1x_dE", 2048, 0, 4096, ic_energy, 600, -300, 300, crdc_1_x);
      obj.FillHistogram(dirname, "crdc1x_Ave", 2048, 0, 4096, ic_ave, 600, -300, 300, crdc_1_x);
      obj.FillHistogram(dirname, "crdc1y_dE", 2048, 0, 4096, ic_energy, 600, -300, 300, crdc_1_y);
      obj.FillHistogram(dirname, "crdc1y_Ave", 2048, 0, 4096, ic_ave, 600, -300, 300, crdc_1_y);
    }
    for (auto ind_out : outgoing_passed){
      for (int a=0; a < 2; a++){
        if (a==1 || afp_gate->IsInside(tof_obje1_corr,afp)){
          dirname = Form("%s_%s", incoming_gates.at(ind_in)->GetName(), outgoing_gates.at(ind_out)->GetName());
          if (a==0) dirname += "_afp";
          dirname += "_gated";

          //CRDC PLOTS
          obj.FillHistogram(dirname, "crdc1_XvsY", 600, -300, 300, crdc_1_x, 1000, -500, 500, crdc_1_y);  
          obj.FillHistogram(dirname, "crdc2_XvsY", 600, -300, 300, crdc_2_x, 1000, -500, 500, crdc_2_y);
          
          //CORRELATION PLOTS
          if (GetGoodMTOFObjE1(s800, tof_obje1_corr)){
            obj.FillHistogram(dirname, "corrobje1_crdc1x", 3000, -3000, 0, tof_obje1_corr,
                                                  600, -300, 300, crdc_1_x);

            obj.FillHistogram(dirname, "corrobje1_afp", 3000, -3000, 0, tof_obje1_corr,
                                                        1000, -0.1, 0.1, afp);

            obj.FillHistogram(dirname, "corrobje1_tofxfpobj", 3000, -3000, 0, tof_obje1_corr,
                                            xocor_nbins, xocor_lowbin, xocor_highbin, xfp_obj);
          }

          if (gretina){
            //ADDBACK STUFF
            TVector3 track = s800->Track();
            int nABevents = gretina->AddbackSize();
            for (int i=0; i < nABevents; i++){
              if (prompt_timing_gate && bank29){
                TGretinaHit abHit = gretina->GetAddbackHit(i);
                TVector3 track1 = s800->Track();
                double abEnergy_corrected = abHit.GetDopplerYta(s800->AdjustedBeta(GValue::Value("BETA")), s800->GetYta(), &track1);
                double time = bank29->Timestamp()-abHit.GetTime();
                
                if (prompt_timing_gate->IsInside(time, abEnergy_corrected)){
                  obj.FillHistogram(dirname, "gamma_corrected_addback_prompt", 8192,0,8192, abEnergy_corrected);
                  
                  //BETA CALIBRATION
                  /*
                  double minBeta = 0.3;
                  double maxBeta = 0.4;
                  double stepSize = 0.001;
                  int nbinsBeta = (maxBeta - minBeta)/stepSize;
                  double beta = minBeta;

                  for (int i=0; i < nbinsBeta; i++){
                    double energyBeta = abHit.GetDoppler(beta);
                    obj.FillHistogram(dirname, "Energy_vs_Beta", nbinsBeta, minBeta, maxBeta, beta, 8192, 0, 8192, energyBeta);
                    beta+=stepSize;
                  }
                  */
                  //GAMMA GAMMA CORRELATION PLOT
                  for (int j=0; j < nABevents; j++){
                    if (i==j) continue;
                    if (prompt_timing_gate && bank29){
                      TGretinaHit abHit2 = gretina->GetAddbackHit(j);
                      TVector3 track2 = s800->Track();
                      double abEnergy_corrected2 = abHit2.GetDopplerYta(s800->AdjustedBeta(GValue::Value("BETA")), s800->GetYta(), &track2);
                      double time2 = bank29->Timestamp()-abHit2.GetTime();
                      if (prompt_timing_gate->IsInside(time2, abEnergy_corrected2)){
                        obj.FillHistogram(dirname, "gamma_gamma", 8192,0,8192, abEnergy_corrected2, 8192,0,8192, abEnergy_corrected);
                      }
                    }
                  }
                }
              }
            }
            
            //OTHER STUFF
            for (unsigned int i = 0; i < gretina->Size(); i++){
              TGretinaHit &hit = gretina->GetGretinaHit(i);
              double energy_corrected = hit.GetDopplerYta(s800->AdjustedBeta(GValue::Value("BETA")), s800->GetYta(), &track);
              double energy = hit.GetDoppler(GValue::Value("BETA"));
              obj.FillHistogram(dirname, "gamma_singles", 8192,0,8192, energy);
              
              if (bank29){
                obj.FillHistogram(dirname,"gamma_corrected_t0_Bank29_time",
                    600,-600,600,bank29->Timestamp()-hit.GetTime(),
                    2500,0,10000,energy_corrected);
              }
              if (prompt_timing_gate && bank29){
                double time = bank29->Timestamp()-hit.GetTime();
                if (prompt_timing_gate->IsInside(time, energy_corrected)){
                  TVector3 position = hit.GetPosition();
                  double theta = position.Theta();
                  
                  obj.FillHistogram(dirname, "gamma_corrected_vs_theta", 8192,0,8192, energy_corrected, 100, 0, 4, theta);
                  obj.FillHistogram(dirname, "core_energy_vs_theta", 8192,0,8192, hit.GetCoreEnergy(), 100, 0, 4, theta);
                  obj.FillHistogram(dirname, "gamma_corrected_singles_prompt", 8192,0,8192, energy_corrected);
                  obj.FillHistogram(dirname, "gamma_singles_prompt", 8192,0,8192, energy);
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
