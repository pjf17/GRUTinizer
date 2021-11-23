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
  std::make_pair(78,76),
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

bool PairHit(const TGretinaHit& one, const TGretinaHit &two, std::vector<std::pair<int, int>> &pairs) {
  int cryId1 = one.GetCrystalId();
  int cryId2 = two.GetCrystalId();
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

  for (auto ind_in: incoming_passed){
    for (auto ind_out : outgoing_passed){
      dirname = Form("%s_%s_gated", incoming_gates.at(ind_in)->GetName(), outgoing_gates.at(ind_out)->GetName());

      if (gretina){    
        if (prompt_timing_gate && bank29){
          double timeBank29 = bank29->Timestamp(); 

          //polarization
          int nHits = gretina->Size();
          
          std::vector<TGretinaHit> gHits;
          for (int i=0; i < nHits; i++){
            gHits.push_back(gretina->GetGretinaHit(i));
          }

          std::sort(gHits.begin(), gHits.end(),
            [](const TGretinaHit& a, const TGretinaHit& b) {
              return a.GetCoreEnergy() > b.GetCoreEnergy();
            });

          TVector3 track = s800->Track();

          for (int i=0; i < nHits; i++){
            //SINGLES
            double energy_corrected = gHits[i].GetDopplerYta(s800->AdjustedBeta(GValue::Value("BETA")), s800->GetYta(), &track);
            obj.FillHistogram(dirname,"singles_crystal_hits",100,0,100,gHits[i].GetCrystalId());
            if (prompt_timing_gate->IsInside(timeBank29-gHits[i].GetTime(),energy_corrected)){
              obj.FillHistogram(dirname,"gamma_corrected_singles_prompt", 8192,0,8192, energy_corrected);
              // if (gHits[i].GetCrystalId() >= 44 && gHits[i].GetCrystalId() <= 79){
              //   obj.FillHistogram(dirname,Form("corrected_prompt_singles_cryID%d",gHits[i].GetCrystalId()), 8192,0,8192, energy_corrected);
              // }
            }

            for (int j=i+1; j < nHits; j++){
              TGretinaHit sum = TGretinaHit(gHits[i]);
              sum.Add(gHits[j]);

              double tot_energy = sum.GetDopplerYta(s800->AdjustedBeta(GValue::Value("BETA")), s800->GetYta(), &track);
              
              //ensure the hit was prompt and that the hits were close in time
              if (prompt_timing_gate->IsInside(timeBank29-sum.GetTime(),tot_energy) && 
                (std::abs(gHits[i].GetTime() - gHits[j].GetTime()) < 44.0) ){
                
                int cryID1 = gHits[i].GetCrystalId(); 
                int cryID2 = gHits[j].GetCrystalId();
                if (cryID1 < cryID2){
                  std::swap(cryID1,cryID2);
                }
                
                if ( PairHit(gHits[i],gHits[j],redPairs) ){
                  obj.FillHistogram(dirname,"gamma_corrected_addback_prompt_red_pair", 8192,0,8192, tot_energy);
                  obj.FillHistogram(dirname,Form("red_pair_%d_%d",cryID1,cryID2), 8192,0,8192, tot_energy);
                }

                if ( PairHit(gHits[i],gHits[j],goldPairs) ){
                  obj.FillHistogram(dirname,"gamma_corrected_addback_prompt_gold_pair", 8192,0,8192, tot_energy);
                  obj.FillHistogram(dirname,Form("gold_pair_%d_%d",cryID1,cryID2), 8192,0,8192, tot_energy);
                }

                if ( PairHit(gHits[i],gHits[j],bluePairs) ){
                  obj.FillHistogram(dirname,"gamma_corrected_addback_prompt_blue_pair", 8192,0,8192, tot_energy);
                  obj.FillHistogram(dirname,Form("blue_pair_%d_%d",cryID1,cryID2), 8192,0,8192, tot_energy);
                }
              }
            } 
          }

          //ADDBACK STUFF
          int nABHits = gretina->AddbackSize();
          
          for (int i=0; i < nABHits; i++){
            TGretinaHit abhit = gretina->GetAddbackHit(i);
            double abEnergy_corrected = abhit.GetDopplerYta(s800->AdjustedBeta(GValue::Value("BETA")), s800->GetYta(), &track);
            
            if (prompt_timing_gate->IsInside(timeBank29-abhit.GetTime(), abEnergy_corrected)){
              obj.FillHistogram(dirname, "gamma_corrected_addback_prompt", 8192,0,8192, abEnergy_corrected);
              obj.FillHistogram(dirname,"addback_crystal_hits",100,0,100,abhit.GetCrystalId());
              // if (gHits[i].GetCrystalId() >= 44 && gHits[i].GetCrystalId() <= 79){
              //   obj.FillHistogram(dirname, Form("corrected_prompt_addback_cryID%d",gHits[i].GetCrystalId()), 8192,0,8192, abEnergy_corrected);  
              // }  
            }
          }

          //NNADDBACK
          for (int n=0; n<2; n++){
            int nnSize = gretina->NNAddbackSize(n);
            for (int i=0; i < nnSize; i++){
              TGretinaHit nnhit = gretina->GetNNAddbackHit(n,i);
              double nnEnergy_corrected = nnhit.GetDopplerYta(s800->AdjustedBeta(GValue::Value("BETA")), s800->GetYta(), &track);
              gretina->GetRingNumber(nnhit);
              if (prompt_timing_gate->IsInside(timeBank29-nnhit.GetTime(), nnEnergy_corrected)){
                obj.FillHistogram(dirname, Form("gamma_corrected_n%d_prompt",n), 8192,0,8192, nnEnergy_corrected);
                obj.FillHistogram(dirname, Form("gamma_corrected_n%d_crystal%d_prompt",n,nnhit.GetCrystalId()), 8192,0,8192, nnEnergy_corrected);
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
