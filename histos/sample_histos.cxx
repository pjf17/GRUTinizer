
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

void CheckGates(std::vector<unsigned short> &incoming_passed, std::vector<unsigned short> &outgoing_passed);

void LoadGates(TRuntimeObjects &obj){
  TList *gates = &(obj.GetGates());
  TIter iter(gates);
  while(TObject *obj = iter.Next()) {
    GCutG *gate = (GCutG*)obj;
    std::string tag = gate->GetTag();
    if(!tag.compare("incoming")) {
      incoming_gates.push_back(gate);
      std::cout << "incoming: << " << gate->GetName() << std::endl;
    } else if(!tag.compare("outgoing")) {
      outgoing_gates.push_back(gate);
      std::cout << "outgoing: << " << gate->GetName() << std::endl;
    } else if(!tag.compare("prompt")){
      prompt_timing_gate = new GCutG(*gate);
      std::cout << "prompt_timing_gate: << " << gate->GetName() << std::endl;
    } else {
      std::cout << "unknown: << " << gate->GetName() << std::endl;
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

  double corr_obj = s800->GetMTofObjE1();
  double crdc_1_x = s800->GetCrdc(0).GetDispersiveX();
  double ic_de;
  if (!std::isnan(GValue::Value("IC_DE_XTILT"))){
    double crdc_1_y = s800->GetCrdc(0).GetNonDispersiveY();
    ic_de = s800->GetIonChamber().GetdE(crdc_1_x, crdc_1_y);
  }
  else{
    ic_de = s800->GetIonChamber().GetAve();
  }
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

  std::string dirname  = "";

  //UPDATE OR REMOVE THESE FOR YOUR ANALYSIS!!!!
  GValue::AddValue(new GValue("TARGET_MTOF_OBJE1",-2000));
  GValue::AddValue(new GValue("TARGET_MTOF_XFPE1", 2500));
  GValue::AddValue(new GValue("CRDC1_X_OFFSET",-281.940));
  GValue::AddValue(new GValue("CRDC1_X_SLOPE",2.540));
  GValue::AddValue(new GValue("CRDC2_X_OFFSET",-281.940));
  GValue::AddValue(new GValue("CRDC2_X_SLOPE",2.540));
  GValue::AddValue(new GValue("CRDC1_Y_OFFSET", 88.9471));
  GValue::AddValue(new GValue("CRDC1_Y_SLOPE", -0.184187));
  GValue::AddValue(new GValue("OBJ_MTOF_CORR_AFP", 1500));
  GValue::AddValue(new GValue("OBJ_MTOF_CORR_XFP", 0.2));
  GValue::AddValue(new GValue("BETA", 0.4));

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


  unsigned short bits = s800->GetTrigger().GetRegistr();
  for(int j=0;j<16;j++) {
    if(((bits>>j)&0x0001))
      obj.FillHistogram("ungated","trig_bit",20,0,20,j);
  }

  double tof_obje1 = s800->GetMTof().GetCorrelatedObjE1(); 
  double tof_xfpe1 =   s800->GetMTof().GetCorrelatedXfpE1();
  obj.FillHistogram("ungated", "incoming_pid", 1000, -4000, 0, tof_obje1,
                                               1000, 0, 4000, tof_xfpe1);

  double ic_ave = s800->GetIonChamber().GetAve();;
  double tof_obje1_corr = s800->GetMTofObjE1(); 
  obj.FillHistogram("ungated", "outgoing_pid", 1000, -4000, 0, tof_obje1_corr,
                                               1024, 0, 4096, ic_ave);

  double crdc_1_x = s800->GetCrdc(0).GetDispersiveX();
  double crdc_2_x = s800->GetCrdc(1).GetDispersiveX();
  double afp = GetAfp(crdc_1_x, crdc_2_x);
  obj.FillHistogram("ungated", "corrobje1_crdc1x", 1000, -4000, 0, tof_obje1_corr,
                                               1024, 0, 4096, crdc_1_x);
  obj.FillHistogram("ungated", "corrobje1_afp", 1000, -4000, 0, tof_obje1_corr,
                                               1024, 0, 4096, afp);
  double xfp_obj = tof_xfpe1-tof_obje1;
  obj.FillHistogram("ungated", "corrobje1_tofxfpobj", 1000, -4000, 0, tof_obje1_corr,
                                               2048, 0, 8192, xfp_obj);

  std::vector<unsigned short> incoming_passed;
  std::vector<unsigned short> outgoing_passed;
  CheckGates(s800, incoming_passed, outgoing_passed);
  for (auto ind_inc : incoming_passed){
    dirname = Form("%s_gated", incoming_gates.at(ind_inc)->GetName());
    obj.FillHistogram(dirname, "outgoing_pid", 1000, -4000, 0, tof_obje1_corr,
        1024, 0, 4096, ic_ave);
    if (gretina){
      for (auto ind_out : outgoing_passed){
        gretina->CleanHits();
        dirname = Form("%s_%s_gated", incoming_gates.at(ind_inc)->GetName(), 
            outgoing_gates.at(ind_out)->GetName()); 
        for (unsigned int i = 0; i < gretina->Size(); i++){
          TGretinaHit &hit = gretina->GetGretinaHit(i);
          double energy = hit.GetDoppler(GValue::Value("BETA"));
          obj.FillHistogram(dirname, "gamma_singles", 8192,0,8192, energy);
          if (prompt_timing_gate){
            if (bank29){
              double time = bank29->Timestamp()-hit.GetTime();
              if (prompt_timing_gate->IsInside(time, energy)){
                obj.FillHistogram(dirname, "gamma_singles_prompt", 8192,0,8192, energy);
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
