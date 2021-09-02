
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


constexpr unsigned int HODO_TIME_CUT = 200;

std::vector<GCutG*> incoming_gates = {};
std::vector<GCutG*> outgoing_gates = {};
std::vector<GCutG*> x_tof_gates = {};
std::vector<GCutG*> x_tof_contam_gates = {};
GCutG *prompt_timing_gate=0;
GCutG *x_de_gate=0;
GCutG *gg_timing_gate=0;

int gates_loaded=0;
double GetAfp(double crdc_1_x,double  crdc_2_x){
  return TMath::ATan( (crdc_2_x - crdc_1_x)/1073.0 );
}


void HandleInverseMap(TRuntimeObjects &obj, const std::string &dirname);
void HandleTotalHodoscope(TRuntimeObjects &obj, const std::string &dirname);
void HandleS800Gated(TRuntimeObjects &obj,GCutG *incoming=0, GCutG *outgoing=0, GCutG *x_tof = 0);
void HandleHodoscopeGretina(TRuntimeObjects &obj, const std::string &dirname, 
                            const std::vector<double> &gret_energies, const std::vector<double> &hodo_energies); 
void CheckGates(std::vector<unsigned short> &incoming_passed, std::vector<unsigned short> &outgoing_passed, 
                std::vector<unsigned short> &x_tof_passed, std::vector<unsigned short> &x_tof_contam_passed);
void HandleMultiplicityAndMatrices(TRuntimeObjects &obj, const std::string &dirname, const std::vector<double> &gret_energies);
void HandleGretinaGated(TRuntimeObjects &obj,GCutG *incoming=0, GCutG *outgoing=0, GCutG *timing=0, GCutG *x_tof=0);
void HandleGretinaGatedAddback(TRuntimeObjects &obj, const std::string &dirname, TGretina *gretina_ab, GCutG *timing, 
                               const TVector3 &track_shifted, double yta, const std::vector<double> &hodo_energies);
bool CheckHodoYrastTransition(const std::vector<double> &hodo_energies);
bool CheckHodoTime(TRuntimeObjects &obj);

void HandleInverseMap(TRuntimeObjects &obj, const std::string &dirname){
  TS800    *s800    = obj.GetDetector<TS800>();
  double ata = s800->GetAta()*1000;
  double bta = s800->GetBta()*1000;
  double yta = s800->GetYta();
  double dta = s800->GetDta();
  double crdc_1_x = s800->GetCrdc(0).GetDispersiveX();
  obj.FillHistogram(dirname,"ata",400,-200,200, ata);
  obj.FillHistogram(dirname,"bta",400,-200,200, bta);
  obj.FillHistogram(dirname,"dta",400,-0.1,0.1, dta);
  obj.FillHistogram(dirname,"dta_x",400,-0.1,0.1, dta,
                                    800, -400, 400, crdc_1_x);
  obj.FillHistogram(dirname,"yta",160,-40,40, yta);
  obj.FillHistogram(dirname,"azita",360,0,360, s800->Azita()*TMath::RadToDeg());
  obj.FillHistogram(dirname,"theta",1800,0,90, s800->Track().Theta()*TMath::RadToDeg());
  obj.FillHistogram(dirname,"ata_bta",400,-200,200,ata,
                                      400, -200, 200,bta);
  obj.FillHistogram(dirname,"dta_ata",400,-0.1,0.1,dta,
                                      400, -200, 200,ata);
  obj.FillHistogram(dirname,"yta_bta",160,-40,40,yta,
                                      400, -200, 200,ata);
}

void HandleTotalHodoscope(TRuntimeObjects &obj, const std::string &dirname){
  TS800    *s800    = obj.GetDetector<TS800>();
  THodoscope &hodo = s800->GetHodoscope();
  bool found_good_time = false;//found at least 1 time greater than HODO_TIME_CUT
  std::vector<double> hodo_energies;

  if (s800->GetMHodoscopeSize() && s800->GetME1Size()){
    for (int i = 0; i < s800->GetMHodoscopeSize(); i++){
      for (int j = 0; j < s800->GetME1Size(); j++){
        double hodo_time = (s800->GetMHodoscope(i) - s800->GetME1Up(j))*0.0625;
        if (hodo_time > HODO_TIME_CUT){
          found_good_time = true;
        }
        obj.FillHistogram(dirname, "hodo_time_alltimeloop", 5000, -1000, 4000, hodo_time);
        for (auto &hit: hodo.GetHodoHits()){
          double energy = hit.GetEnergy();
          obj.FillHistogram(dirname,"hodo_time_energy_alltimeloop", 5000, -1000, 4000, hodo_time,
              4096,0,4096, energy);
          hodo_energies.push_back(energy);
        }//for hit in hodo hits
      }//loop over E1
    }//loop over hodo
  }//times exist

  if (found_good_time){
    for (auto energy: hodo_energies){
      obj.FillHistogram(dirname,Form("hodo_energy_tcut%d_noloop",HODO_TIME_CUT), 512 ,0,4096, energy);
    }
  }

  for (auto &hit : hodo.GetHodoHits()){
    int channel = hit.GetChannel();
    obj.FillHistogram(dirname,"hodo_charge",4096,0,4096, hit.GetCharge());
    obj.FillHistogram(dirname,"hodo_channel",100,0,100,channel);
    obj.FillHistogram(dirname,"hodo_charge_channel_summary", 32, 0, 32,channel,
        4096,0,4096, hit.GetCharge());
    obj.FillHistogram(dirname,"hodo_energy_channel_summary", 32, 0, 32,channel,
        4096,0,4096, hit.GetEnergy());
  }
}

void HandleS800Gated(TRuntimeObjects &obj, GCutG *incoming, GCutG *outgoing, GCutG *x_tof) {
  TS800    *s800    = obj.GetDetector<TS800>();
  if(!s800)
    return;

  std::string dirname = "s800_gated";

  if (incoming){
    dirname += Form("_%s", incoming->GetName());
  }
  if (outgoing){
    dirname += Form("_%s", outgoing->GetName());
  }
  if (x_tof){//already passing only passed cuts
    dirname += Form("_%s", x_tof->GetName());
  }

  HandleInverseMap(obj, dirname);
  HandleTotalHodoscope(obj, dirname);

  double crdc_1_x = s800->GetCrdc(0).GetDispersiveX();
  double crdc_1_y = s800->GetCrdc(0).GetNonDispersiveY();
  double crdc_2_x = s800->GetCrdc(1).GetDispersiveX();
  double crdc_2_y = s800->GetCrdc(1).GetNonDispersiveY();
  obj.FillHistogram(dirname,"CRDC1_X",800,-400,400,crdc_1_x);
  obj.FillHistogram(dirname,"CRDC2_X",800,-400,400,crdc_2_x);
  obj.FillHistogram(dirname,"CRDC1_Y",800,-400,400,crdc_1_y);
  obj.FillHistogram(dirname,"CRDC2_Y",800,-400,400,crdc_2_y);

  obj.FillHistogram(dirname, "CRDC1_XY", 800, -400, 400, crdc_1_x,
                                         800,-400,400, crdc_2_x);

  double corr_obj = s800->GetMTofObjE1();
  obj.FillHistogram(dirname,"corr_obj_vs_crdc_1_x",2000,-4000,0,corr_obj,
      600,-300,300,crdc_1_x);

  obj.FillHistogram(dirname,"corr_obj_vs_crdc_2_x",2000,-4000,0,corr_obj,
      600,-300,300,crdc_2_x);
  double ic_de;
  if (!std::isnan(GValue::Value("IC_DE_XTILT"))){
    ic_de = s800->GetIonChamber().GetdE(crdc_1_x, crdc_1_y);
  }
  else{
    ic_de = s800->GetIonChamber().GetAve();
  }
  obj.FillHistogram(dirname,"x_de", 600,-300,300,crdc_1_x,
      4096, 0, 4096, ic_de);

  double afp = GetAfp(crdc_1_x, crdc_2_x);
  obj.FillHistogram(dirname,"time_afp",2000,-4000,0,corr_obj,
                                       1000,-.1,.1, afp);

  double bfp = s800->GetBFP();//mrad
  obj.FillHistogram(dirname,"time_bfp",2000,-4000,0,corr_obj,
                                       1000,-.1,.1, bfp);

  obj.FillHistogram(dirname, "pid", 2000,-4000,0, corr_obj,
                                    4096, 0, 4096, ic_de);

  unsigned short bits = s800->GetTrigger().GetRegistr();
  for(int j=0;j<16;j++) {
    if(((bits>>j)&0x0001))
      obj.FillHistogram(dirname,"trig_bit",20,0,20,j);
    if (crdc_1_x < 160){
      obj.FillHistogram(dirname,"trig_bit_crdc1xlt160",20,0,20,j);
    }
  }

  obj.FillHistogram(dirname, "trig_reg",20,0,20, bits);
  double obj_tof = s800->GetMTof().GetCorrelatedObjE1();
  double xfp_tof = s800->GetMTof().GetCorrelatedXfpE1();
  obj.FillHistogram(dirname, "incoming_pid", 1000, -2000, -1000, obj_tof,
                                            4000, 2000, 6000, xfp_tof);
  obj.FillHistogram(dirname, "xfp_vs_crdc_1_x",4000, 2000, 6000, xfp_tof,
                                               800,-400,400, crdc_1_x) ;
  obj.FillHistogram(dirname, "xfp_vs_crdc_2_x",4000, 2000, 6000, xfp_tof,
                                               800,-400,400, crdc_2_x) ;

  return;
}

bool CheckHodoTime(TRuntimeObjects &obj){
  TS800    *s800    = obj.GetDetector<TS800>();
  if (s800->GetMHodoscopeSize() && s800->GetME1Size()){
    for (int i = 0; i < s800->GetMHodoscopeSize(); i++){
      for (int j = 0; j < s800->GetME1Size(); j++){
        if((s800->GetMHodoscope(i) - s800->GetME1Up(j))*0.0625 > HODO_TIME_CUT){
          return true;
        }
      }//loop over E1 times
    }//loop over Hodo times
  }//both times exist
  return false;
}
bool CheckHodoYrastTransition(const std::vector<double> &hodo_energies){
  const double PERCENT_DIFF = 0.1;
  for (auto h_nrg : hodo_energies){
    if ((h_nrg > 183-PERCENT_DIFF*183 && h_nrg < 183+PERCENT_DIFF*183) ||
        (h_nrg > 448-PERCENT_DIFF*448 && h_nrg < 448+PERCENT_DIFF*448) ||
        (h_nrg > 970-PERCENT_DIFF*970 && h_nrg < 970+PERCENT_DIFF*970) ||
        (h_nrg > 1260-PERCENT_DIFF*1260 && h_nrg < 1260+PERCENT_DIFF*1260)){
      return true;
    }//check for yrast transition
  }//loop over hodo energies
  return false;
}



void HandleHodoscopeGretina(TRuntimeObjects &obj, const std::string &dirname, 
                            const std::vector<double> &gret_energies, const std::vector<double> &hodo_energies){ 

  bool found_good_time = CheckHodoTime(obj);
  bool found_yrast_transition = CheckHodoYrastTransition(hodo_energies);

  for (auto corr_gret_energy : gret_energies){
    for (auto hodo_energy : hodo_energies){
      obj.FillHistogram(dirname,"gretina_hodoscope_full", 512,0,4096, hodo_energy,
          5000,0,10000, corr_gret_energy);
      if (found_good_time){
        obj.FillHistogram(dirname,Form("gretina_hodo_energy_tcut_%dns",HODO_TIME_CUT), 
            512, 0, 4096, hodo_energy, 
            5000, 0, 10000, corr_gret_energy);
      }
    }
    
    if (found_yrast_transition && found_good_time){
      obj.FillHistogram(dirname, Form("gretina_energy_hodoyrast_tcut_%dns", HODO_TIME_CUT),
                          5000, 0, 10000, corr_gret_energy);
    }
    if (found_yrast_transition){
      obj.FillHistogram(dirname, "gretina_energy_hodoyrast", 5000, 0, 10000, corr_gret_energy);
    }   
  }//loop over gretina energies
}

void HandleGretinaGatedAddback(TRuntimeObjects &obj, const std::string &dirname, TGretina *gretina_ab,
                               GCutG *timing, const TVector3 &track_shifted, double yta,
                               const std::vector<double> &hodo_energies){

  TS800    *s800    = obj.GetDetector<TS800>();
  TBank29    *bank29    = obj.GetDetector<TBank29>();
  std::vector<double> ab_hits_energy;
  for(int x=0;x<gretina_ab->AddbackSize();x++) {
    const TGretinaHit &ab_hit = gretina_ab->GetAddbackHit(x);
    double energy = ab_hit.GetDoppler(GValue::Value("BETA"));
    double time = bank29->Timestamp()-ab_hit.GetTime();
    if (timing->IsInside(time, energy)){
      double energy_track_yta_dta_shifts = ab_hit.GetDopplerYta(s800->AdjustedBeta(GValue::Value("BETA")), yta, &track_shifted);
      int number = ab_hit.GetNumber();
      ab_hits_energy.push_back(energy_track_yta_dta_shifts);
      obj.FillHistogram(dirname,"gretina_ab_allcorr_summary_prompt", 50, 0, 50, number,
          10000,0,10000, energy_track_yta_dta_shifts);
      obj.FillHistogram(dirname,"gretina_ab_allcorr_prompt", 10000,0,10000, energy_track_yta_dta_shifts);
    }//inside timing gate
  }//loop over addback hits


  std::string dirname_hodo(dirname);
  dirname_hodo += "_hodoscope";
  HandleHodoscopeGretina(obj, dirname_hodo, ab_hits_energy, hodo_energies);

  std::string dirname_mgg(dirname);
  dirname_mgg += "_multcuts_ggmatrices";
  HandleMultiplicityAndMatrices(obj, dirname_mgg, ab_hits_energy);

}

void HandleMultiplicityAndMatrices(TRuntimeObjects &obj, const std::string &dirname, 
                                             const std::vector<double> &gret_energies){
  if (gret_energies.size() == 1){
    obj.FillHistogram(dirname,"mult_1", 5000,0,10000,gret_energies.at(0));
  }
  if (gret_energies.size() == 2){
    obj.FillHistogram(dirname,"mult_2", 5000,0,10000,gret_energies.at(0));
    obj.FillHistogram(dirname,"mult_2", 5000,0,10000,gret_energies.at(1));
  }
  if (gret_energies.size() > 2){
    for (auto &nrg : gret_energies){
      obj.FillHistogram(dirname,"mult_greater_than_3", 5000,0,10000, nrg);
    }
  }

  TS800    *s800    = obj.GetDetector<TS800>();
  double dta = s800->GetDta();
  for (unsigned int i = 0; i < gret_energies.size(); i++){
    for (unsigned int j = i+1; j < gret_energies.size(); j++){
      obj.FillHistogram(dirname,"gg_matrix", 6000,0,6000, gret_energies.at(i),
                                          6000,0,6000, gret_energies.at(j));
      obj.FillHistogram(dirname,"gg_matrix", 6000,0,6000, gret_energies.at(j),
                                          6000,0,6000, gret_energies.at(i));

      obj.FillHistogram(dirname, "dta_gretina_matrix", 
          1000, -1, 1, dta, 
          6000,0,6000, gret_energies.at(i));


      if(gret_energies.size() == 2){
        obj.FillHistogram(dirname,"mult2_gg_matrix", 6000,0,6000,gret_energies.at(i),
                                            6000,0,6000,gret_energies.at(j));
        obj.FillHistogram(dirname,"mult2_gg_matrix", 6000,0,6000,gret_energies.at(j),
                                            6000,0,6000,gret_energies.at(i));   
      }
    }
  }


}
void HandleGretinaGated(TRuntimeObjects &obj,GCutG *incoming, GCutG *outgoing, GCutG *timing, GCutG *x_tof) {

  TGretina *gretina = obj.GetDetector<TGretina>();
  TGretina *gretina_ab = new TGretina(*gretina);
  gretina_ab->Clear();
  TS800    *s800    = obj.GetDetector<TS800>();
  TBank29  *bank29  = obj.GetDetector<TBank29>();

  std::string dirname("gretina");
  if (incoming){
    dirname += Form("_%s", incoming->GetName());
  }
  if (outgoing){
    dirname += Form("_%s", outgoing->GetName());
  }
  if (timing){
    dirname += Form("_%s", timing->GetName());
  }
  if (x_tof){
    dirname += Form("_%s", x_tof->GetName());
  }
  double ata_shift = GValue::Value("GRETINA_ATA_OFFSET");
  double bta_shift = GValue::Value("GRETINA_BTA_OFFSET");
  TVector3 track_shifted = s800->Track(ata_shift, bta_shift);
  double yta = s800->GetYta();
  double afp = s800->GetAFP()*1000;//mrad
  double crdc_1_x = s800->GetCrdc(0).GetDispersiveX();

  //for gg_matrix and multiplicity spectra
  std::vector<double> nonab_hits_energy;
  for(unsigned int x = 0; x < gretina->Size(); x++){
    double energy = gretina->GetGretinaHit(x).GetDoppler(GValue::Value("BETA"));
    double energy_stopped = gretina->GetGretinaHit(x).GetDoppler(0.0);
    double time = bank29->Timestamp()-gretina->GetGretinaHit(x).GetTime();
    obj.FillHistogram(dirname,"gretina_stopped", 10000,0,10000, energy_stopped);
    if (timing->IsInside(time, energy)){
      gretina_ab->InsertHit(gretina->GetGretinaHit(x));
      int number = gretina->GetGretinaHit(x).GetNumber();
      double energy_track_yta_dta_shifts = gretina->GetGretinaHit(x).GetDopplerYta(s800->AdjustedBeta(GValue::Value("BETA")), yta, &track_shifted);
      nonab_hits_energy.push_back(energy_track_yta_dta_shifts);
      obj.FillHistogram(dirname,"gretina_summary_allcorr_prompt", 50, 0, 50, number,
          10000,0,10000, energy_track_yta_dta_shifts);
      obj.FillHistogram(dirname,"gretina_allcorr_prompt", 10000,0,10000, energy_track_yta_dta_shifts);
      obj.FillHistogram(dirname,"gretina_allcorr_theta_energy", 10000,0,10000, energy_track_yta_dta_shifts,
                                                                   2000, 0,2, gretina->GetGretinaHit(x).GetTheta());
      obj.FillHistogram(dirname,"gretina_stopped_prompt", 10000,0,10000, energy_stopped);

      //handle forward detectors
      if (number < 16){
        obj.FillHistogram(dirname,"energy_vs_afp_forward", 10000,0,10000, energy_track_yta_dta_shifts,
                                                             3600, -180, 180, afp);
        obj.FillHistogram(dirname,"energy_vs_theta_forward", 10000,0,10000, energy_track_yta_dta_shifts,
                                                             2000, 0, 2, gretina->GetGretinaHit(x).GetTheta());
      }
      else if (number >= 16){
        obj.FillHistogram(dirname,"energy_vs_afp_90deg", 10000,0,10000, energy_track_yta_dta_shifts,
                                                             3600, -180, 180, afp);
        obj.FillHistogram(dirname,"energy_vs_theta_90deg", 10000,0,10000, energy_track_yta_dta_shifts,
                                                             2000, 0, 2, gretina->GetGretinaHit(x).GetTheta());
      }

      if (crdc_1_x <= 160){
        obj.FillHistogram(dirname,"gretina_allcorr_prompt_crdcxbelow160", 10000,0,10000, energy_track_yta_dta_shifts);
      }
      else{
        obj.FillHistogram(dirname,"gretina_allcorr_prompt_crdcxabove160", 10000,0,10000, energy_track_yta_dta_shifts);
      }
    }//inside timing gate
  }//loop over gretina hits

  std::vector<double> hodo_energies;
  for (const auto &hh : s800->GetHodoscope().GetHodoHits()){
    hodo_energies.push_back(hh.GetEnergy());
  }
 
  // std::string dirname_hodo(dirname); 
  // dirname_hodo += "_hodoscope";
  //HandleHodoscopeGretina(obj, dirname_hodo, nonab_hits_energy, hodo_energies);
  
  //std::string dirname_ab(dirname);
  //dirname_ab += "_addback";
  //HandleGretinaGatedAddback(obj, dirname_ab, gretina_ab, timing, track_shifted, yta, hodo_energies); 

  //std::string dirname_mgg(dirname);
  //dirname_mgg += "_multcuts_ggmatrices";
  //HandleMultiplicityAndMatrices(obj, dirname_mgg, nonab_hits_energy);
  return;
}


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
    } else if(!tag.compare("gretina")) {
      std::string cutname = gate->GetName();
      if(!cutname.compare("timing")) {
        gg_timing_gate=new GCutG(*gate);
      }
    } else if(!tag.compare("prompt")){
      prompt_timing_gate = new GCutG(*gate);
      std::cout << "timing: << " << gate->GetName() << std::endl;
    } else if(!tag.compare("x_de")){
      x_de_gate = new GCutG(*gate);
      std::cout << "x_de: << " << gate->GetName() << std::endl;
    } else if(!tag.compare("x_tof")){
      x_tof_gates.push_back(gate);
      std::cout << "x_tof: << " << gate->GetName() << std::endl;
    } else if(!tag.compare("x_tof_contam")){
      x_tof_contam_gates.push_back(gate);
      std::cout << "x_tof_contam: << " << gate->GetName() << std::endl;
    } 
    else {
      std::cout << "unknown: << " << gate->GetName() << std::endl;
    }
    gates_loaded++;
  }
  std::cout << "outgoing size: " << outgoing_gates.size() << std::endl;
}

void CheckGates(TS800 *s800, std::vector<unsigned short> &incoming_passed, std::vector<unsigned short> &outgoing_passed,
                std::vector<unsigned short> &x_tof_passed, std::vector<unsigned short> &x_tof_contam_passed){
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
  
  for (unsigned int i = 0; i < x_tof_gates.size(); i++){
    if (x_tof_gates.at(i)->IsInside(corr_obj, crdc_1_x)){
      x_tof_passed.push_back(i);
    }
  }

  for (unsigned int i = 0; i < x_tof_contam_gates.size(); i++){
    if (!x_tof_contam_gates.at(i)->IsInside(corr_obj, crdc_1_x)){
      x_tof_contam_passed.push_back(i);
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

  if(bank29 && gretina) {
    for(unsigned int i=0;i<gretina->Size();i++) {
      TGretinaHit &hit = gretina->GetGretinaHit(i);
      obj.FillHistogram("Bank29","Gretina_dop_t0_Bank29_time",
          600,-600,600,bank29->Timestamp()-hit.GetTime(),
          5000,0,10000, hit.GetDoppler(GValue::Value("BETA")));
    }//loop over gretina hits
  }//bank29 and gretina exist

  std::vector<unsigned short> incoming_passed;
  std::vector<unsigned short> outgoing_passed;
  std::vector<unsigned short> x_tof_passed;
  std::vector<unsigned short> x_tof_contam_passed;

  unsigned short bits = s800->GetTrigger().GetRegistr();
  for(int j=0;j<16;j++) {
    if(((bits>>j)&0x0001))
      obj.FillHistogram("ungated","trig_bit",20,0,20,j);
  }

  obj.FillHistogram("ungated", "trig_reg",20,0,20, bits);
  
  double rawobj_e1 = s800->GetMTof().GetCorrelatedObjE1();
  double rawxfp_e1 = s800->GetMTof().GetCorrelatedXfpE1();
  obj.FillHistogram("ungated", "incoming_pid", 4000, -4000, 0, rawobj_e1,
                                               4000, 0, 4000, rawxfp_e1);

  double ic_ave = s800->GetIonChamber().GetAve();
  obj.FillHistogram("ungated", "outgoing_pid_rawtof", 4000, -4000, 0, rawobj_e1,
                                               4096, 0, 4096, ic_ave);
  double corr_obj = s800->GetMTofObjE1Chn15();
  obj.FillHistogram("ungated", "outgoing_pid_rawtof", 4000, -4000, 0, corr_obj,
                                               4096, 0, 4096, ic_ave);

  double ic_de;
  double crdc_1_x = s800->GetCrdc(0).GetDispersiveX();
  double crdc_1_y = s800->GetCrdc(0).GetNonDispersiveY();
  if (!std::isnan(GValue::Value("IC_DE_XTILT"))){
    ic_de = s800->GetIonChamber().GetdE(crdc_1_x, crdc_1_y);
  }
  else{
    ic_de = s800->GetIonChamber().GetAve();
  }
  obj.FillHistogram("ungated", "pid", 2000,-4000,0, corr_obj,
                                    4096, 0, 4096, ic_de);

  if(incoming_gates.at(0)->IsInside(s800->GetMTof().GetCorrelatedObjE1(), 
        s800->GetMTof().GetCorrelatedXfpE1())){
    obj.FillHistogram("inCu71", "pid", 2000,-4000,0, corr_obj,
        4096, 0, 4096, ic_de);

  }
  CheckGates(s800, incoming_passed, outgoing_passed, x_tof_passed, x_tof_contam_passed);


  if(numobj!=list->GetSize()){
    list->Sort();
  }
}