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

std::vector<std::pair<int,int>> redPairs = {
  std::make_pair(46,44),
  std::make_pair(46,48),
  std::make_pair(48,49),
  std::make_pair(50,51),
  std::make_pair(51,47),
  std::make_pair(56,57),
  std::make_pair(57,61),
  std::make_pair(59,58),
  std::make_pair(60,58),
  std::make_pair(60,62),
  std::make_pair(64,62),
  std::make_pair(64,65),
  std::make_pair(69,65),
  std::make_pair(63,67),
  std::make_pair(66,67),
  std::make_pair(66,68),
  std::make_pair(70,68),
  std::make_pair(78,76)
};

std::vector<std::pair<int,int>> goldPairs = {
  std::make_pair(44,45),
  std::make_pair(46,47),
  std::make_pair(48,51),
  std::make_pair(49,50),
  std::make_pair(56,59),
  std::make_pair(57,58),
  std::make_pair(61,60),
  std::make_pair(62,63),
  std::make_pair(64,67),
  std::make_pair(65,66),
  std::make_pair(69,68),
  std::make_pair(70,71),
  std::make_pair(78,79),
};

std::vector<std::pair<int,int>> bluePairs = {
  std::make_pair(44,47),
  std::make_pair(45,46),
  std::make_pair(46,51),
  std::make_pair(48,50),
  std::make_pair(56,58),
  std::make_pair(57,60),
  std::make_pair(60,63),
  std::make_pair(61,62),
  std::make_pair(62,67),
  std::make_pair(64,66),
  std::make_pair(65,68),
  std::make_pair(68,71),
  std::make_pair(76,79),
};

bool PairHit(const TGretinaHit& abhit, std::vector<std::pair<int, int>> &pairs) {
  int cryId1 = abhit.GetCrystalId();
  int cryId2 = abhit.GetNeighbor().GetCrystalId();
  bool hit = false;
  
  for (auto &p : pairs){
    if ( (cryId1 == p.first && cryId2 == p.second) 
        || (cryId2 == p.first && cryId1 == p.second) ) {
        hit = true;
        break;
    }
  }
  return hit;
}

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

  //INCOMING
  for (auto ind_in : incoming_passed){
    dirname = Form("%s_gated", incoming_gates.at(ind_in)->GetName());
    double ic_ave = s800->GetIonChamber().GetAve();
    if (GetGoodMTOFObjE1(s800, tof_obje1_corr)){
      obj.FillHistogram(dirname, "outgoing_pid", 4000, -4000, 0, tof_obje1_corr, 2048, 0, 4096, ic_energy);                                                
    }

    //CRDC DE CORRECTION
    for (auto ind_iso : isoline_passed){
      dirname = Form("%s_%s_gated", incoming_gates.at(ind_in)->GetName(), isoline_gates.at(ind_iso)->GetName());
      obj.FillHistogram(dirname, "crdc1x_dE", 2048, 0, 4096, ic_energy, 600, -300, 300, crdc_1_x);
      obj.FillHistogram(dirname, "crdc1x_Ave", 2048, 0, 4096, ic_ave, 600, -300, 300, crdc_1_x);
      obj.FillHistogram(dirname, "crdc1y_dE", 2048, 0, 4096, ic_energy, 600, -300, 300, crdc_1_y);
      obj.FillHistogram(dirname, "crdc1y_Ave", 2048, 0, 4096, ic_ave, 600, -300, 300, crdc_1_y);
    }
    
    //OUTGOING
    for (auto ind_out : outgoing_passed){
      dirname = Form("%s_%s_gated", incoming_gates.at(ind_in)->GetName(), outgoing_gates.at(ind_out)->GetName());

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
        TVector3 track = s800->Track();
        if (bank29){
          double timeBank29 = bank29->Timestamp();
          
          //SINGLES
          int gSize = gretina->Size();
          for (int i=0; i < gSize; i++){
            TGretinaHit &hit = gretina->GetGretinaHit(i);
            double energy_corrected = hit.GetDopplerYta(s800->AdjustedBeta(GValue::Value("BETA")), s800->GetYta(), &track);
            double energy = hit.GetDoppler(GValue::Value("BETA"));
            double core_energy = hit.GetCoreEnergy();
            double theta = hit.GetTheta();
            int cryID = hit.GetCrystalId();
            //deal with the crystal ID gaps
            // if (cryID >= 40 && cryID < 52 ){
            //   cryID -= 4;
            // } else if (cryID >= 52 && cryID < 72 ) {
            //   cryID -=8;
            // } else if (cryID == 76) {
            //   cryID -=12;
            // } else if (cryID > 77){
            //   cryID -=13;
            // }
            // cryID -= 24;

            //time energy spectrum
            if (bank29){
              obj.FillHistogram(dirname,"gamma_corrected_t0_Bank29_time",
                  600,-600,600,timeBank29-hit.GetTime(),
                  2500,0,10000,energy_corrected);
            }

            //Make all the singles spectra
            obj.FillHistogram(dirname, "gamma_singles", 8192,0,8192, energy);
            //prompt
            if (prompt_timing_gate && prompt_timing_gate->IsInside(timeBank29-hit.GetTime(), energy_corrected)){
              obj.FillHistogram(dirname, "gamma_singles_prompt", 8192,0,8192, energy);
              obj.FillHistogram(dirname, "gamma_corrected_singles_prompt", 8192,0,8192, energy_corrected);
              obj.FillHistogram(dirname, "gamma_corrected_vs_theta_prompt", 8192,0,8192, energy_corrected, 100, 0, 2.5, theta);
              obj.FillHistogram(dirname, "gamma_corrected_vs_crystalID_prompt", 56, 24, 80, cryID, 8192,0,8192, energy_corrected);
              obj.FillHistogram(dirname, "core_energy_vs_theta_prompt", 8192,0,8192, hit.GetCoreEnergy(), 100, 0, 2.5, theta);
            //off prompt
            } else {
              obj.FillHistogram(dirname, "gamma_singles_off-prompt", 8192,0,8192, energy);
              obj.FillHistogram(dirname, "core_energy_off-prompt", 8192,0,8192, core_energy);
              obj.FillHistogram(dirname, "core_energy_vs_theta_off-prompt", 8192,0,8192, core_energy, 100, 0, 2.5, theta);
              obj.FillHistogram(dirname, "core_energy_vs_crystalID_off-prompt", 56, 24, 80, cryID, 8192,0,8192, core_energy);
              obj.FillHistogram(dirname, "core_energy_vs_theta_off-prompt", 8192,0,8192, core_energy, 100, 0, 2.5, theta);
            }
          }

          //NNADDBACK
          //loop over multiplicity
          for (int n=0; n<4; n++){
            //loop over hits for each multiplicity spectrum
            int nnSize = gretina->NNAddbackSize(n);
            for (int i=0; i < nnSize; i++){

              //get hit and hit data 
              TGretinaHit nnhit = gretina->GetNNAddbackHit(n,i);
              int cryID = nnhit.GetCrystalId();
              int ringNum = nnhit.GetRingNumber();
              double nnEnergy_corrected = nnhit.GetDopplerYta(s800->AdjustedBeta(GValue::Value("BETA")), s800->GetYta(), &track);
              
              //make sure hits are prompt
              if (prompt_timing_gate && prompt_timing_gate->IsInside(timeBank29-nnhit.GetTime(), nnEnergy_corrected)){
                //exclude the ng spectrum (n==3)
                if (n < 3){
                  obj.FillHistogram(dirname, "gamma_corrected_addback_prompt", 8192,0,8192, nnEnergy_corrected);
                  // //GAMMA GAMMA CORRELATION
                  // for (int j=i+1; j < nnSize; j++){
                  //   if (i==j) continue;
                  //   TGretinaHit nnhit2 = gretina->GetNNAddbackHit(n,j);
                  //   double nnEnergy_corrected2 = nnhit.GetDopplerYta(s800->AdjustedBeta(GValue::Value("BETA")), s800->GetYta(), &track);
                  //   if (prompt_timing_gate->IsInside(timeBank29-nnhit.GetTime(), nnEnergy_corrected)){
                  //     obj.FillHistogram(dirname, "gamma_gamma", 8192,0,8192, nnEnergy_corrected2, 8192,0,8192, nnEnergy_corrected);
                  //   }
                  // }
                }

                if (n == 1){
                  //POLARIZATION
                  if ( PairHit(nnhit,redPairs) ){
                    obj.FillHistogram(dirname,"gamma_corrected_addback_prompt_red_pair", 8192,0,8192, nnEnergy_corrected);
                  }

                  if ( PairHit(nnhit,goldPairs) ){
                    obj.FillHistogram(dirname,"gamma_corrected_addback_prompt_gold_pair", 8192,0,8192, nnEnergy_corrected);
                  }

                  if ( PairHit(nnhit,bluePairs) ){
                    obj.FillHistogram(dirname,"gamma_corrected_addback_prompt_blue_pair", 8192,0,8192, nnEnergy_corrected);
                  }

                  //CRDC 2 XY correlation
                  obj.FillHistogram(dirname, "n1_addback_crdc2_y", 8192,0,8192, nnEnergy_corrected, 400, -200, 200, crdc_2_y);
                  obj.FillHistogram(dirname, "n1_addback_crdc2_y", 8192,0,8192, nnEnergy_corrected, 600, -300, 300, crdc_2_x);
                }

                char *multiplicity = Form("%d",n);
                if (n == 3) multiplicity = Form("g");
                obj.FillHistogram(dirname, Form("gamma_corrected_n%s_prompt",multiplicity), 8192,0,8192, nnEnergy_corrected);
                obj.FillHistogram(dirname, Form("gamma_corrected_n%s_ring%02d_crystal%d_prompt",multiplicity,ringNum,cryID),8192,0,8192, nnEnergy_corrected);
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
