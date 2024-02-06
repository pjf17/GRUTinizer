// To make a new filter, copy this file under a new name in the "filter" directory.
// The "FilterCondition" function should return a boolean value.
// The boolean indicates whether the event should be kept or not.

#include <map>
#include <vector>

#include "TRuntimeObjects.h"
#include "GCutG.h"
#include "TFile.h"
#include "TS800.h"
#include "TGretina.h"
#include "TBank29.h"

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

bool gates_loaded = false;
std::map<std::string,std::vector<GCutG*>> gates;

extern "C"
bool FilterCondition(TRuntimeObjects& obj) {
  TS800 *s800 = obj.GetDetector<TS800>();
  TGretina *gretina = obj.GetDetector<TGretina>();
  TBank29 *bank29  = obj.GetDetector<TBank29>();

  if (!s800 || !gretina || !bank29){
    return false;
  }

  //load in the gates
  if (!gates_loaded) {
    LoadGates(&(obj.GetGates()),gates);
    gates_loaded = true;
  }

  bool gamma_gate = false;
  int nGret = gretina->Size();
  for (int i=0 ; i < nGret; i++){
    if (gretina->GetGretinaHit(i).GetCoreEnergy() > 100) {
      gamma_gate = true;
      break;
    }
  }

  double tof_xfpe1_corr = s800->GetTof().GetTacXFP() + GValue::Value("XFP_CORR_AFP") * s800->GetAFP() +  GValue::Value("XFP_CORR_XFP") * s800->GetCrdc(0).GetDispersiveX();

  if (gamma_gate && !gates["dirbeam"][0]->IsInside(tof_xfpe1_corr,s800->GetCrdc(0).GetDispersiveX())){ 
    return true;
  }
  
  return false;
}
