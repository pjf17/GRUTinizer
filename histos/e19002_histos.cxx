#include "TRuntimeObjects.h"

#include <iostream>
#include <map>
#include <cstdio>
#include <vector>
#include <string>

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
          
          //SINGLES
          int gSize = gretina->Size();
          for (int i=0; i < gSize; i++){
            TGretinaHit &hit = gretina->GetGretinaHit(i);
            hit.ComptonSort();
            double energy_corrected = hit.GetDopplerYta(s800->AdjustedBeta(GValue::Value("BETA")), s800->GetYta(), &track);
            double energy = hit.GetDoppler(GValue::Value("BETA"));
            double core_energy = hit.GetCoreEnergy();
            double theta = hit.GetTheta();
            int cryID = hit.GetCrystalId();
            int nInteractions = hit.NumberOfInteractions();

            //time energy spectrum
            if (bank29){
              obj.FillHistogram(dirname,"gamma_corrected_t0_Bank29_time",
                  600,-600,600,timeBank29-hit.GetTime(),
                  2500,0,10000,energy_corrected);
            }

            //Make all the singles spectra
            obj.FillHistogram(dirname, "gamma_singles", 8192,0,8192, energy);
            
            //PROMPT GATE
            bool tgate = false;
            if (gates["prompt"].size() > 0) tgate = gates["prompt"][0]->IsInside(timeBank29-hit.GetTime(), energy_corrected);
            else std::cout<<"NO PROMPT GATE LOADED\n";

            if (tgate){
              if (nInteractions > 1){
                double EPsum = 0;
                for (int ip=0; ip < nInteractions; ip++) EPsum+= hit.GetSegmentEng(ip);
                obj.FillHistogram(dirname, "Ecore_vs_intPSum",1024,0,2048,EPsum,1024,0,2048,core_energy);
                // for (int in=0; in < nInteractions; in++){
                //   double th = hit.GetTheta(in)*TMath::RadToDeg();
                //   double eng = hit.GetDopplerYta(s800->AdjustedBeta(GValue::Value("BETA")), s800->GetYta(), &track, in);
                //   obj.FillHistogram(dirname,"theta_vs_gam_dop_sgl_prompt_INTP",180,0,180,th,4096,0,4096,eng);
                //   obj.FillHistogram(dirname,"theta_vs_core_eng_INTP",180,0,180,th,4096,0,4096,core_energy);
                // }
                double xi = hit.GetXi();
                double nu = hit.GetScatterAngle();
                double alpha = hit.GetAlpha();
                double Eratio = 511.0/core_energy * hit.GetSegmentEng(0)/(core_energy - hit.GetSegmentEng(0));

                // if (energy_corrected > 1000 && energy_corrected < 1036){
                //   obj.FillHistogram(dirname, "1018-cos_vs_Eratio",1000,0,10,Eratio,200,-1,1,TMath::Cos(nu));
                // }

                double plXi = xi;
                double xiMax = TMath::TwoPi();
                obj.FillHistogram(dirname, "gam_dop_sgl_prompt_vs_xi",360,0,xiMax,plXi, 2048,0,4096, energy_corrected);
                
                if (Eratio > 2) obj.FillHistogram(dirname, "eratio>2",2048,0,4096, energy_corrected);
                else obj.FillHistogram(dirname, "gam_dop_sgl_prompt_vs_phase_reg",15,-1.5,3.5,TMath::Cos(nu)+Eratio,2048,0,4096, energy_corrected);
                // if (nu > 1.2*TMath::ACos(1-511/core_energy)) {
                if (TMath::Cos(nu) > 0.25 - Eratio) {
                  obj.FillHistogram(dirname, "gam_dop_sgl_prompt_vs_xi_scatter_gated",360,0,xiMax,plXi, 2048,0,4096, energy_corrected);
                  if (theta*TMath::RadToDeg() > 55 && theta*TMath::RadToDeg() < 100)
                    obj.FillHistogram(dirname, "gam_dop_sgl_prompt_vs_xi_scatter_gated_theta55-100",360,0,xiMax,plXi, 2048,0,4096, energy_corrected);
                  else
                    obj.FillHistogram(dirname, "gam_dop_sgl_prompt_vs_xi_scatter_gated_theta_NOT_55-100",360,0,xiMax,plXi, 2048,0,4096, energy_corrected);
                  // double thAngles[6] = {40,55,73,95,120,149}; 
                  // for (int th=0; th < 5; th++){
                  //   if (theta*TMath::RadToDeg() >= thAngles[th] && theta*TMath::RadToDeg() < thAngles[th+1]) 
                  //     obj.FillHistogram(dirname, Form("gam_dop_sgl_prompt_vs_xi_scatter_gated_theta%3.0f-%3.0f",thAngles[th],thAngles[th+1]),360,0,TMath::TwoPi(),xi, 1024,0,2048, energy_corrected);
                  // }
                } 
                else {
                  obj.FillHistogram(dirname, "gam_dop_sgl_prompt_vs_xi_scatter_not_gated",360,0,xiMax,plXi, 2048,0,4096, energy_corrected);
                }
                for (auto gate : gates["pol"]){
                  if (gate->IsInside(xi,energy_corrected)){
                    obj.FillHistogram(dirname, Form("xi_%s_gated",gate->GetName()),360,0,TMath::TwoPi(),xi);
                  }
                }
              }
              
              obj.FillHistogram(dirname, "core_energy_prompt", 8192,0,8192, core_energy);
              obj.FillHistogram(dirname, "core_energy_vs_theta_prompt", 360, 0, 180, theta*TMath::RadToDeg(), 4000,0,4000, hit.GetCoreEnergy());
              obj.FillHistogram(dirname, "gam_dop_sgl_prompt", 8192,0,8192, energy_corrected);
              obj.FillHistogram(dirname, "gam_dop_sgl_prompt_summary", 48, 0, 48, detMapRing[cryID], 2048,0,2048, energy_corrected);
              obj.FillHistogram(dirname, "gam_dop_sgl_prompt_vs_theta", 180, 0, 180, theta*TMath::RadToDeg(), 4000,0,4000, energy_corrected);

              // if (cryID > 40) 
              //   obj.FillHistogram(dirname, "gamma_corrected_singles_prompt_90qds", 8192,0,8192, energy_corrected);
              // else 
              //   obj.FillHistogram(dirname, "gamma_corrected_singles_prompt_fwdqds", 8192,0,8192, energy_corrected);
            } 
          }

          //NNADDBACK
          //loop over multiplicity
          // for (int n=0; n<4; n++){
          //   //loop over hits for each multiplicity spectrum
          //   int nnSize = gretina->NNAddbackSize(n);
          //   for (int i=0; i < nnSize; i++){

          //     //get hit and hit data 
          //     TGretinaHit nnhit = gretina->GetNNAddbackHit(n,i);
          //     // int cryID = nnhit.GetCrystalId();
          //     // int ringNum = nnhit.GetRingNumber();
          //     double nnEnergy_corrected = nnhit.GetDopplerYta(s800->AdjustedBeta(GValue::Value("BETA")), s800->GetYta(), &track);
          //     // double nnCore_energy = nnhit.GetCoreEnergy();
          //     // double theta = nnhit.GetThetaDeg();
          //     // double phi = nnhit.GetPhiDeg();
          //     int nInteractions = nnhit.NumberOfInteractions();
              
          //     bool tgate = false;
          //     if (gates["prompt"].size() > 0) tgate = gates["prompt"][0]->IsInside(timeBank29-nnhit.GetTime(), nnEnergy_corrected);

          //     //make sure hits are prompt
          //     if (tgate){
          //       //exclude the ng spectrum (n==3)
          //       // obj.FillHistogram(dirname,"crystal-map",180,0,180,theta,360,0,360,phi);
          //       char *multiplicity = Form("%d",n);
          //       if (n == 3) multiplicity = Form("g");
          //       obj.FillHistogram(dirname, Form("gamma_corrected_n%s_prompt",multiplicity), 8192,0,8192, nnEnergy_corrected);
                
          //       if (n < 3)
          //         obj.FillHistogram(dirname, "gamma_corrected_addback_prompt", 8192,0,8192, nnEnergy_corrected);
          //     }
          //   }
          // }
                  // obj.FillHistogram(dirname, "gamma_corrected_addback_prompt_vs_nInteractions",20,0,20,nInteractions,1024,0,2048,nnEnergy_corrected);
                  // obj.FillHistogram(dirname, "core_energy_addback_prompt", 8192,0,8192,nnCore_energy);
          //         if (nInteractions > 1){
          //           // TVector3 pos1 = nnhit.GetIntPosition(0);
          //           // TVector3 pos2 = nnhit.GetIntPosition(1);
          //           // double angleSpread = pos1.Angle(pos2)*TMath::RadToDeg();
          //           // obj.FillHistogram(dirname, "gamma_corrected_addback_prompt_vs_angleSpread",150,0,30,angleSpread,1024,0,2048,nnEnergy_corrected);
          //           double aziCompt = azimuthalCompton(nnhit,&track)*TMath::RadToDeg();
          //           obj.FillHistogram(dirname, "energy_corrected_vs_xi",360,0,360,aziCompt,2048,0,2048,nnEnergy_corrected);
          //           if (aziCompt > 90 && aziCompt < 270){
          //             obj.FillHistogram(dirname,"Xi_center_gamma_vs_theta",180,0,180,theta,8192,0,8192,nnEnergy_corrected);
          //             obj.FillHistogram(dirname,"Xi_center_gamma_vs_phi",360,0,360,phi,8192,0,8192,nnEnergy_corrected);
          //           } else {
          //             obj.FillHistogram(dirname,"Xi_extreme_gamma_vs_theta",180,0,180,theta,8192,0,8192,nnEnergy_corrected);
          //             obj.FillHistogram(dirname,"Xi_extreme_gamma_vs_phi",360,0,360,phi,8192,0,8192,nnEnergy_corrected);
          //           }
          //           if (theta > 100) 
          //             obj.FillHistogram(dirname, "azmthl_compton_theta>100",360,0,360,aziCompt,1024,0,2048,nnEnergy_corrected);
          //           else if (theta < 50) 
          //             obj.FillHistogram(dirname, "azmthl_compton_theta<50",360,0,360,aziCompt,1024,0,2048,nnEnergy_corrected);
          //           else if (theta > 75 && theta < 95) 
          //             obj.FillHistogram(dirname, "azmthl_compton_theta_flat",360,0,360,aziCompt,1024,0,2048,nnEnergy_corrected);
          //           else 
          //             obj.FillHistogram(dirname, "azmthl_compton_theta_varied",360,0,360,aziCompt,1024,0,2048,nnEnergy_corrected);
          //         }
          //       }

          //       //if (n == 1) {
          //         //GAMMA GAMMA CORRELATION
          //         // for (int j=0; j < nnSize; j++){
          //         //   if (i==j) continue;
          //         //   TGretinaHit nnhit2 = gretina->GetNNAddbackHit(n,j);
          //         //   double nnEnergy_corrected2 = nnhit2.GetDopplerYta(s800->AdjustedBeta(GValue::Value("BETA")), s800->GetYta(), &track);
          //         //   if (prompt_timing_gate->IsInside(timeBank29-nnhit.GetTime(), nnEnergy_corrected)){
          //         //     obj.FillHistogram(dirname, "gamma_gamma", 8192,0,8192, nnEnergy_corrected2, 8192,0,8192, nnEnergy_corrected);
          //         //   }
          //         // }

          //         //make swapped spectra
          //         // TGretinaHit nnhit1 = nnhit.GetInitialHit();
          //         // TGretinaHit nnhit2 = nnhit.GetNeighbor();
          //         // nnhit2.NNAdd(nnhit1);
          //         // double swappedEnergy = nnhit2.GetDopplerYta(s800->AdjustedBeta(GValue::Value("BETA")), s800->GetYta(), &track);

          //         // //POLARIZATION
          //         // std::string polColor = "blank";
          //         // if (PairHit(nnhit,redPairs)) polColor = "red";
          //         // else if (PairHit(nnhit,goldPairs)) polColor = "gold";
          //         // else if (PairHit(nnhit,bluePairs)) polColor = "blue";

          //         // std::string quadType = "blank";
          //         // if (PairHit(nnhit,OneQuadPlus)) quadType = "qd1+";
          //         // else if (PairHit(nnhit,OneQuadDefault)) quadType = "qd1";
          //         // else if (PairHit(nnhit,TwoQuadPairs)) quadType = "qd2";

          //         // if ( polColor.compare("blank") != 0 ){
          //         //   obj.FillHistogram(dirname,Form("ab_prompt_%s_%s_pair",polColor.c_str(),quadType.c_str()), 8192,0,8192, nnEnergy_corrected);
          //         //   obj.FillHistogram(dirname,Form("ab_prompt_%s_pair",polColor.c_str()), 8192,0,8192, nnEnergy_corrected);
          //         //   obj.FillHistogram(dirname,Form("ab_prompt_%s_pair",quadType.c_str()), 8192,0,8192, nnEnergy_corrected);
          //         //   obj.FillHistogram(dirname,Form("ab_swapped_prompt_%s_%s_pair",polColor.c_str(),quadType.c_str()), 8192,0,8192, swappedEnergy);
          //         // }

          //       //}
          //     }
          //   }
          // }
        } 
      }
    }
  }
  if(numobj!=list->GetSize()){
    list->Sort();
  }
}
