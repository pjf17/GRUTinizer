
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
GCutG *crdc_gate=0;
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
    } else if(!tag.compare("crdc")){
      crdc_gate = new GCutG(*gate);
      std::cout << "\t crdc_gate: << " << gate->GetName() << std::endl;
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

  double corr_obje1 = s800->GetMTofObjE1();
  double icE = GetGoodICE(s800);
  
  for(unsigned int i=0;i<outgoing_gates.size();i++) {
    if(outgoing_gates.at(i)->IsInside(corr_obje1,icE)){
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

  std::vector<unsigned short> incoming_passed;
  std::vector<unsigned short> outgoing_passed;
  CheckGates(s800, incoming_passed, outgoing_passed);
  std::string dirname  = "";

  for (auto ind_in : incoming_passed){
    dirname = Form("%s_gated", incoming_gates.at(ind_in)->GetName());

    double tof_obje1_corr;
    double ic_energy = GetGoodICE(s800);

    if (GetGoodMTOFObjE1(s800, tof_obje1_corr)){
      //OUTGOING PID
      obj.FillHistogram(dirname, "outgoing_pid", 
        4000, -4000, 0, tof_obje1_corr, 2048, 0, 4096, ic_energy); 
    }

    for (auto ind_out : outgoing_passed){
      dirname = Form("%s_%s_gated", incoming_gates.at(ind_in)->GetName(), outgoing_gates.at(ind_out)->GetName());
      obj.FillHistogram(dirname, "outgoing_pid", 
        4000, -4000, 0, tof_obje1_corr, 2048, 0, 4096, ic_energy);
    }

  }
  
  if(numobj!=list->GetSize()){
    list->Sort();
  }
}
