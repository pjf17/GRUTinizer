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

std::map<int,int> detMap = {
  {24, 0}, {25, 1}, {26, 2}, {27, 3}, {28, 4}, {29, 5}, {30, 6}, {31, 7},
  {32, 8}, {33, 9}, {34,10}, {35,11}, {36,12}, {37,13}, {38,14}, {39,15},
  {44,16}, {45,17}, {46,18}, {47,19}, {48,20}, {49,21}, {50,22}, {51,23},
  {56,24}, {57,25}, {58,26}, {59,27}, {60,28}, {61,29}, {62,30}, {63,31},
  {64,32}, {65,33}, {66,34}, {67,35}, {68,36}, {69,37}, {70,38}, {71,39},
  {76,40}, {77,41}, {78,42}, {79,43}, {80,44}, {81,45}, {82,46}, {83,47}
};

std::map<int,int> detMapRing = {
  {26, 0}, {30, 1}, {34, 2}, {38, 3}, {25, 4}, {29, 5}, {33, 6}, {37, 7},
  {27, 8}, {31, 9}, {35,10}, {39,11}, {24,12}, {28,13}, {32,14}, {36,15},
  {47,16}, {63,17}, {71,18}, {79,19}, {51,20}, {59,21}, {67,22}, {83,23},
  {50,24}, {58,25}, {66,26}, {82,27}, {44,28}, {60,29}, {68,30}, {76,31},
  {46,32}, {62,33}, {70,34}, {78,35}, {48,36}, {56,37}, {64,38}, {80,39},
  {49,40}, {57,41}, {65,42}, {81,43}, {45,44}, {61,45}, {69,46}, {77,47}
};

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

  //---------------------------------------------------------------
  //GATED

  std::vector<unsigned short> incoming_passed;
  std::vector<unsigned short> outgoing_passed;
  std::vector<unsigned short> isoline_passed;
  CheckGates(s800, incoming_passed, outgoing_passed, isoline_passed);
  std::string dirname  = "";

  //INCOMING
  for (auto ind_in : incoming_passed){
    
    //OUTGOING
    for (auto ind_out : outgoing_passed){
      dirname = Form("%s_%s_gated", incoming_gates.at(ind_in)->GetName(), outgoing_gates.at(ind_out)->GetName());
      
      if (gretina){
        if (bank29 && prompt_timing_gate){
          double yta = s800->GetYta();
          double dta = s800->GetDta();
          double ata = s800->GetAta();
          double bta = s800->GetBta();
          obj.FillHistogram("s800","ata",500,-0.1,0.1,ata);
          obj.FillHistogram("s800","bta",500,-0.1,0.1,bta);
          obj.FillHistogram("s800","yta",1000,-20,20,yta);
          obj.FillHistogram("s800","dta",1000,-1,1,dta);

          TVector3 track = s800->Track();
          double timeBank29 = bank29->Timestamp();
          int nGretina = gretina->Size();
          double BETA = GValue::Value("BETA");

          //BETA CORRECTION
          for (int g=0; g < nGretina; g++){
            TGretinaHit &hit = gretina->GetGretinaHit(g);
            int cryID = hit.GetCrystalId();
            int ringnum = hit.GetRingNumber();
            
            if (!isnan(GValue::Value("BETA_SCAN_STEP"))){
              //loop to find optimal beta
              double betaMin = GValue::Value("BETA_SCAN_MIN");
              double betaMax = GValue::Value("BETA_SCAN_MAX");
              double betaStep = GValue::Value("BETA_SCAN_STEP");
              int nBetaBins = (betaMax - betaMin)/betaStep;
              double beta = betaMin;
              for (int i=0; i < nBetaBins; i++){
                double energy = hit.GetDoppler(beta); 
                if (prompt_timing_gate->IsInside(timeBank29-hit.GetTime(),energy)){
                  obj.FillHistogram(dirname,"Energy_vs_beta",nBetaBins,betaMin,betaMax,beta,4000,0,4000,energy);
                  obj.FillHistogram(dirname,Form("Theta_vs_Energy_beta%f",beta),300,900,1200,energy,100,0,3,hit.GetTheta());
                  obj.FillHistogram(dirname,Form("summary_beta%f",beta),48,0,48,detMapRing[cryID],300,900,1200,energy);
                }
                beta += betaStep;
              }
            } else {
              //finer doppler corrections after finding a good beta
              double energy_b = hit.GetDoppler(BETA);
              double energy_bt = hit.GetDoppler(BETA,&track);
              double energy_bty = hit.GetDopplerYta(BETA,s800->GetYta(),&track); 
              double energy_btyd = hit.GetDopplerYta(s800->AdjustedBeta(BETA),s800->GetYta(),&track); 
              if (prompt_timing_gate->IsInside(timeBank29-hit.GetTime(),energy_btyd)){
                //PHI CORRELATION
                // double phiTa = s800->Azita();//TMath::ATan(TMath::Sin(ata)/(TMath::Sin(bta)*-1.0));
                // if (phiTa < 0) phiTa += TMath::TwoPi();
                double phi = (TMath::TwoPi() - s800->Azita() - hit.GetPhi())*TMath::RadToDeg();
                if (phi < 0) phi += 360;

                obj.FillHistogram(dirname,"Phi_vs_Energy",4000,0,4000,energy_b,360,0,360,phi);
                obj.FillHistogram(dirname,"Phi_vs_Energy_corrected",4000,0,4000,energy_bt,360,0,360,phi);
              
                // //YTA CORRELATION
                // obj.FillHistogram(dirname,Form("Yta_vs_Energy_r%02d_c%d",ringnum,cryID),700,600,1300,energy_b,200,-20,20,yta);
                // obj.FillHistogram(dirname,Form("Yta_vs_Energy_corrected_r%02d_c%d",ringnum,cryID),700,600,1300,energy_bty,200,-20,20,yta);

                // //DTA CORRELATION
                // obj.FillHistogram(dirname,Form("Dta_vs_Energy_r%02d_c%d",ringnum,cryID),700,600,1300,energy_b,50,-0.06,-0.02,dta);
                // obj.FillHistogram(dirname,Form("Dta_vs_Energy_corrected_r%02d_c%d",ringnum,cryID),700,600,1300,energy_btyd,50,-0.06,-0.02,dta);

                //SUMMARY SPECTRUM
                obj.FillHistogram(dirname,"Doppler_summary",48,0,48,detMapRing[cryID],2000,0,2000,energy_btyd);

                //Spectrum
                obj.FillHistogram(dirname,"gamma_singles_corrected",4000,0,4000,energy_btyd);
                obj.FillHistogram(dirname,Form("gamma_singles_corrected_i%02d",detMapRing[cryID]),4000,0,4000,energy_btyd);
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
