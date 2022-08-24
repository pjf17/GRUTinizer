
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
GCutG *prompt_timing_gate=0;
int gates_loaded=0;

double GetAfp(double crdc_1_x,double  crdc_2_x){
  return TMath::ATan( (crdc_2_x - crdc_1_x)/1073.0 );
}

//Get the corresponding OBJE1 TOF based on whether there are AFP and XFP corrections
double GetGoodMTOFObjE1(TS800 *s800){
  double value = 0;
  if(std::isnan(GValue::Value("OBJ_MTOF_CORR_AFP")) || 
    std::isnan(GValue::Value("OBJ_MTOF_CORR_XFP")) ) {
    value = s800->GetMTof().GetCorrelatedObjE1();

  } else {
    value = s800->GetMTofObjE1();
  }
  return value;
}

//Get the Ion Chamber DE depending on whether IC_DE_XTILT is set
double GetGoodIC_DE(TS800 *s800){
  double value = 0;
  double crdc_1_x = s800->GetCrdc(0).GetDispersiveX();
  
  if (!std::isnan(GValue::Value("IC_DE_XTILT"))){
    double crdc_1_y = s800->GetCrdc(0).GetNonDispersiveY();
    value = s800->GetIonChamber().GetdE(crdc_1_x, crdc_1_y);
  } else{
    value = s800->GetIonChamber().GetAve();
  }
  
  return value;
}

void CheckGates(std::vector<unsigned short> &incoming_passed, std::vector<unsigned short> &outgoing_passed);

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
    } else if(!tag.compare("prompt")){
      prompt_timing_gate = new GCutG(*gate);
      std::cout << "\t prompt_timing_gate: << " << gate->GetName() << std::endl;
    } else {
      std::cout << "\t unknown: << " << gate->GetName() << std::endl;
    }
    gates_loaded++;
  }
  std::cout << "outgoing size: " << outgoing_gates.size() << std::endl;
}

void CheckGates(TS800 *s800, std::vector<unsigned short> &incoming_passed, std::vector<unsigned short> &outgoing_passed){
  for(unsigned short i=0;i<incoming_gates.size();i++) {
    if(incoming_gates.at(i)->IsInside(s800->GetMTof().GetCorrelatedObjE1(), 
                                      s800->GetMTof().GetCorrelatedXfpE1())){
      incoming_passed.push_back(i);
    }
  }

  double corr_obj = GetGoodMTOFObjE1(s800);
  double ic_de = GetGoodIC_DE(s800);
  
  for(unsigned int i=0;i<outgoing_gates.size();i++) {
    if(outgoing_gates.at(i)->IsInside(corr_obj,ic_de)){
      outgoing_passed.push_back(i);
    }
  }
}

// extern "C" is needed to prevent name mangling.
// The function signature must be exactly as shown here,
//   or else bad things will happen.
extern "C"
void MakeHistograms(TRuntimeObjects& obj) {
//  InitMap();
  //TGretina *gretina = obj.GetDetector<TGretina>();
  //TBank29  *bank29  = obj.GetDetector<TBank29>();
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

  std::string dirname  = "";

  //UPDATE OR REMOVE THESE FOR YOUR ANALYSIS!!!!
  //GValue::AddValue(new GValue("TARGET_MTOF_OBJE1CHN15",-2000));
  //GValue::AddValue(new GValue("TARGET_MTOF_XFPE1CHN15", 2000));
  //GValue::AddValue(new GValue("CRDC1_X_OFFSET",-281.940));
  //GValue::AddValue(new GValue("CRDC1_X_SLOPE",2.540));
  //GValue::AddValue(new GValue("CRDC2_X_OFFSET",-281.940));
  //GValue::AddValue(new GValue("CRDC2_X_SLOPE",2.540));
  //GValue::AddValue(new GValue("CRDC1_Y_OFFSET", 88.9471));
  //GValue::AddValue(new GValue("CRDC1_Y_SLOPE", -0.184187));
  //GValue::AddValue(new GValue("OBJ_MTOF_CORR_AFP", 1500));
  //GValue::AddValue(new GValue("OBJ_MTOF_CORR_XFP", 0.2));
  //GValue::AddValue(new GValue("BETA", 0.4));

  /*
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
  */

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
  
  obj.FillHistogram("ungated", "MTOF_OBJE1", 2000, -5000,-3000, raw_obj - raw_e1);
  obj.FillHistogram("ungated", "MTOF_XFE1", 2000, 1000, 6000, raw_xf - raw_e1);

  //MAKE INCOMING PID
  double tof_obje1 = s800->GetMTof().GetCorrelatedObjE1(); 
  double tof_xfpe1 = s800->GetMTof().GetCorrelatedXfpE1();
  obj.FillHistogram("ungated", "incoming_pid", 2000, -5000, -3000, tof_obje1,
                                               2000, 1000, 6000, tof_xfpe1);                                              

  //MAKE OUTGOING PID
  double ic_ave = s800->GetIonChamber().GetAve();
  double tof_obje1_corr = GetGoodMTOFObjE1(s800);
  obj.FillHistogram("ungated", "outgoing_pid", 1000, -5000, -3000, tof_obje1_corr,
                                               1024, 0, 4096, ic_ave);

  //CORRELATION PLOTS
  double crdc_1_x = s800->GetCrdc(0).GetDispersiveX();
  double crdc_2_x = s800->GetCrdc(1).GetDispersiveX();
  double afp = GetAfp(crdc_1_x, crdc_2_x);
  obj.FillHistogram("ungated", "corrobje1_crdc1x", 1000, -5000, -3000, tof_obje1_corr,
                                               600, -300, 300, crdc_1_x);
  obj.FillHistogram("ungated", "corrobje1_afp", 1000, -5000, -3000, tof_obje1_corr,
                                               1000, -0.1, 0.1, afp);
  double xfp_obj = tof_xfpe1-tof_obje1;
  obj.FillHistogram("ungated", "corrobje1_tofxfpobj", 1000, -5000, -3000, tof_obje1_corr,
                                               2048, -1024, 1024, xfp_obj-GValue::Value("TOFXFP_OBJ_SHIFT"));
  
  //CRDC Coordinates
  double crdc_1_y = s800->GetCrdc(0).GetNonDispersiveY();
  double crdc_2_y = s800->GetCrdc(1).GetNonDispersiveY();

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
  
  //---------------------------------------------------------------
  //GATED

  std::vector<unsigned short> incoming_passed;
  std::vector<unsigned short> outgoing_passed;
  CheckGates(s800, incoming_passed, outgoing_passed);
  
  for (auto ind_out : outgoing_passed){
    dirname = Form("%s_gated", outgoing_gates.at(ind_out)->GetName());
    obj.FillHistogram(dirname, "incoming_pid", 1000, -6000, -2000, tof_obje1,
                                               1000, 0, 4000, tof_xfpe1);
  }
  
  for (auto ind_in : incoming_passed){
    dirname = Form("%s_gated", incoming_gates.at(ind_in)->GetName());
    obj.FillHistogram(dirname, "outgoing_pid", 1000, -6000, -2000, tof_obje1_corr,
                                               2048, 0, 4096, ic_ave);
  }
  
  if(numobj!=list->GetSize()){
    list->Sort();
  }
}
