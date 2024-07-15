
#include "TRuntimeObjects.h"

#include <iostream>
#include <map>
#include <cstdio>

#include <TF1.h>
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

bool checkEnergyTheta(TF1* fdop, TF1* fres, double energy, double theta){
  return fdop->Eval(theta) + fres->Eval(theta) > energy && fdop->Eval(theta) - fres->Eval(theta) < energy;
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
  
  if (!gretina || !s800 || !bank29){
    return;
  }

  bool gamma_gate = false;
  int nGret = gretina->Size();
  for (int i=0 ; i < nGret; i++){
    if (gretina->GetGretinaHit(i).GetCoreEnergy() > 100) {
      gamma_gate = true;
      break;
    }
  }

  if (!gamma_gate) return;
  
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
          2500,0,10000, hit.GetCoreEnergy());
    }//loop over gretina hits
  }//bank29 and gretina exist

  //---------------------------------------------------------------
  //UNGATED

  unsigned short bits = s800->GetTrigger().GetRegistr();
  int trig_bit = -1;
  for(int j=0;j<16;j++) {
    if(((bits>>j)&0x0001)){
      obj.FillHistogram("ungated","trig_bit",20,0,20,j);
      trig_bit = j;
    }
  }
  
  //MAKE RAW TOF HISTS
  //mesytec
  obj.FillHistogram("ungated", "MTOF_OBJE1", 1000, -10000, 0, s800->GetOBJ_E1Raw());
  obj.FillHistogram("ungated", "MTOF_XFE1", 1000, -6000, 6000, s800->GetXF_E1Raw());
  obj.FillHistogram("ungated", "MTOF_XFE1_vs_MTOF_OBJE1", 3000, -3000, 3000, s800->GetXF_E1Raw(), 1500, -3000, 0, s800->GetOBJ_E1Raw());

  //TAC
  double raw_obje1 = s800->GetTof().GetTacOBJ();
  double raw_xfpe1 = s800->GetTof().GetTacXFP();
  obj.FillHistogram("ungated", "RAW_TOF_OBJE1", 5000, 0, 5000, raw_obje1);
  obj.FillHistogram("ungated", "RAW_TOF_XFE1", 5000, 0, 5000, raw_xfpe1);
  
  //Correlated TOFs 
  double tof_obje1 = s800->GetMTof().GetCorrelatedObjE1(); 
  double tof_xfpe1 = s800->GetMTof().GetCorrelatedXfpE1();
  //MAKE INCOMING PID
  obj.FillHistogram("ungated", "TOF_OBJE1", 5000, -10000, 5000, tof_obje1);
  obj.FillHistogram("ungated", "TOF_XFPE1", 5000, -5000, 10000, tof_xfpe1);

  obj.FillHistogram("ungated", "incoming_pid", 900, 2100, 3000, tof_xfpe1, 900, -2300, -1400, tof_obje1);
  obj.FillHistogram("ungated", Form("incoming_pid_trigbit%d",trig_bit), 900, 2100, 3000, tof_xfpe1, 900, -2300, -1400, tof_obje1);

  //CRDC Values
  double crdc_1_x = s800->GetCrdc(0).GetDispersiveX();
  double crdc_2_x = s800->GetCrdc(1).GetDispersiveX();
  double afp = GetAfp(crdc_1_x, crdc_2_x);
  
  //CRDC Coordinates
  double crdc_1_y = s800->GetCrdc(0).GetNonDispersiveY();
  double crdc_2_y = s800->GetCrdc(1).GetNonDispersiveY();

  double ylow = -500;
  double yhigh = 500;
  double ybins = 1000;
  
  double yslope = GValue::Value("CRDC1_Y_SLOPE");
  if (std::isnan(yslope) || yslope == 0){
    ylow = 0;
    yhigh = 3000;
    ybins = 3000;
  }
  obj.FillHistogram("ungated", "crdc1 X_Y", 600, -300, 300, crdc_1_x, ybins, ylow, yhigh, crdc_1_y);  
  obj.FillHistogram("ungated", "crdc2 X_Y", 600, -300, 300, crdc_2_x, ybins, ylow, yhigh, crdc_2_y);

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
  CheckGates(gates["incoming"],incoming_passed,tof_xfpe1,tof_obje1);

  //GET TOFs AND DE FOR PID AND CORRELATION PLOTS
  double tof_obje1_corr = 0;
  double tof_xfpe1_corr = 0;
  double ic_energy = GetGoodICE(s800);

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
  for (auto ind_in: incoming_passed){
    dirname = Form("%s_gated", gates["incoming"].at(ind_in)->GetName());
    obj.FillHistogram(dirname, "crdc1 X_Y", 600, -300, 300, crdc_1_x, ybins, ylow, yhigh, crdc_1_y);  
    obj.FillHistogram(dirname, "crdc2 X_Y", 600, -300, 300, crdc_2_x, ybins, ylow, yhigh, crdc_2_y);

    //Correlation plots
    obj.FillHistogram(dirname, "corrobje1_crdc1x", 800, -2200, -1400, tof_obje1_corr,600, -300, 300, crdc_1_x);
    obj.FillHistogram(dirname, "obje1_crdc1x", 800, -2200, -1400, tof_obje1,600, -300, 300, crdc_1_x);                                     
    obj.FillHistogram(dirname, "corrobje1_afp", 800, -2200, -1400, tof_obje1_corr,1000, -0.1, 0.1, afp);
    obj.FillHistogram(dirname, "obje1_afp", 800, -2200, -1400, tof_obje1,1000, -0.1, 0.1, afp);
    obj.FillHistogram(dirname, "corrobje1_tofxfpobj", 800, -2200, -1400, tof_obje1_corr,xocor_nbins, xocor_lowbin, xocor_highbin, xfp_obj);
    obj.FillHistogram(dirname, "obje1_tofxfpobj", 800, -2200, -1400, tof_obje1,xocor_nbins, xocor_lowbin, xocor_highbin, xfp_obj);
    
    //MAKE OUTGOING PID
    double ic_ave = s800->GetIonChamber().GetAve();
    obj.FillHistogram(dirname, "outgoing_pid_uncorrected", 800, -2200, -1400, tof_obje1, 2048, 0, 4096, ic_ave);
    obj.FillHistogram(dirname, "outgoing_pid", 800, -2200, -1400, tof_obje1_corr, 2048, 0, 4096, ic_energy);

    //CRDC DE correction
    if (gates.count("isoline") && gates["isoline"][0]->IsInside(tof_obje1_corr,ic_ave)) {
      dirname = Form("%s_%s_gated", gates["incoming"].at(ind_in)->GetName(), gates["isoline"][0]->GetName());
      obj.FillHistogram(dirname, "crdc1x_dE", 2048, 0, 4096, ic_energy, 600, -300, 300, crdc_1_x);
      obj.FillHistogram(dirname, "crdc1x_Ave", 2048, 0, 4096, ic_ave, 600, -300, 300, crdc_1_x);
      obj.FillHistogram(dirname, "crdc1y_dE", 2048, 0, 4096, ic_energy, 600, -300, 300, crdc_1_y);
      obj.FillHistogram(dirname, "crdc1y_Ave", 2048, 0, 4096, ic_ave, 600, -300, 300, crdc_1_y);
    }

    //energy time
    for(unsigned int i=0;i<gretina->Size();i++) {
      //Time-energy cut
      TGretinaHit &hit = gretina->GetGretinaHit(i);
      obj.FillHistogram(dirname,"Gretina_dop_t0_Bank29_time",
          1200,-300,300,bank29->Timestamp()-hit.GetTime(),
          1000,0,10000, hit.GetCoreEnergy());
    }//loop over gretina hits
    
    for (auto ind_out : outgoing_passed){
      //ADDITIONAL TIME OF FLIGHT GATES
      int nTOFgates = 0;
      if (gates.count("crdc_tof")) nTOFgates = gates["crdc_tof"].size();

      int tof_gate_idx = -1;
      for (int tfg=0; tfg < nTOFgates; tfg++){
        std::string tof_gate_name = std::string(gates["crdc_tof"][tfg]->GetName());
        tof_gate_name = tof_gate_name.substr(tof_gate_name.rfind("_")+1,tof_gate_name.length()-1-tof_gate_name.rfind("_"));

        std::string outgoing_name = std::string(gates["outgoing"].at(ind_out)->GetName());
        if (outgoing_name.find(tof_gate_name) != std::string::npos) tof_gate_idx = tfg;
      }

      if (tof_gate_idx == -1)
        dirname = Form("%s_%s_gated",gates["incoming"].at(ind_in)->GetName(),gates["outgoing"].at(ind_out)->GetName());
      else if (gates["crdc_tof"][tof_gate_idx]->IsInside(tof_obje1_corr,crdc_1_x))
        dirname = Form("%s_%s_%s_gated",gates["incoming"].at(ind_in)->GetName(),gates["outgoing"].at(ind_out)->GetName(),gates["crdc_tof"][tof_gate_idx]->GetName());
      else
        continue;

      //Correlation plots
      obj.FillHistogram(dirname, "corrobje1_crdc1x", 800, -2200, -1400, tof_obje1_corr,600, -300, 300, crdc_1_x);
      obj.FillHistogram(dirname, "obje1_crdc1x", 800, -2200, -1400, tof_obje1,600, -300, 300, crdc_1_x);                                     
      obj.FillHistogram(dirname, "corrobje1_afp", 800, -2200, -1400, tof_obje1_corr,1000, -0.1, 0.1, afp);
      obj.FillHistogram(dirname, "obje1_afp", 800, -2200, -1400, tof_obje1,1000, -0.1, 0.1, afp);
      obj.FillHistogram(dirname, "corrobje1_tofxfpobj", 800, -2200, -1400, tof_obje1_corr,xocor_nbins, xocor_lowbin, xocor_highbin, xfp_obj);
      obj.FillHistogram(dirname, "obje1_tofxfpobj", 800, -2200, -1400, tof_obje1,xocor_nbins, xocor_lowbin, xocor_highbin, xfp_obj);

      obj.FillHistogram(dirname, "outgoing_pid", 800, -2200, -1400, tof_obje1_corr, 2048, 0, 4096, ic_energy);

      //energy time
      for(unsigned int i=0;i<gretina->Size();i++) {
        //Time-energy cut
        TGretinaHit &hit = gretina->GetGretinaHit(i);
        obj.FillHistogram(dirname,"Gretina_dop_t0_Bank29_time",
            1200,-300,300,bank29->Timestamp()-hit.GetTime(),
            1000,0,10000, hit.GetCoreEnergy());
      }//loop over gretina hits

      //GAMMA ANALYSIS
      TVector3 track = s800->Track();
      double timeBank29 = bank29->Timestamp();

      //SINGLES
      int gSize = gretina->Size();
      double BETA = GValue::Value("BETA");
      for (int i=0; i < gSize; i++){
        TGretinaHit &hit = gretina->GetGretinaHit(i);
        double theta = hit.GetTheta();
        double phi = hit.GetPhi();
        double core_energy = hit.GetCoreEnergy();
        int cryID = hit.GetCrystalId();
        int ringnum = hit.GetRingNumber();
        int nInteractions = hit.NumberOfInteractions();

        //PROMPT GATE
        bool tgate = false;
        if (gates["prompt"].size() > 0) tgate = gates["prompt"][0]->IsInside(timeBank29-hit.GetTime(), core_energy);
        else std::cout<<"NO PROMPT GATE LOADED\n";

        if (tgate){
          obj.FillHistogram(dirname,"gretina_theta_vs_phi",360,0,360,phi*TMath::RadToDeg(),180,0,180,theta*TMath::RadToDeg());
          obj.FillHistogram(dirname,"core_energy_summary",48,0,48,cryID,4000,0,4000,core_energy);
          /*
          if (!isnan(GValue::Value("BETA_SCAN_STEP"))){
            // beta scan parameters
            double betaMin = GValue::Value("BETA_SCAN_MIN");
            double betaMax = GValue::Value("BETA_SCAN_MAX");
            double betaStep = GValue::Value("BETA_SCAN_STEP");
            int nBetaBins = (betaMax - betaMin)/betaStep;
            double beta = betaMin;

            //histogram parameters
            double scanElo = GValue::Value("BETA_SCAN_ELO");
            double scanEhi = GValue::Value("BETA_SCAN_EHI");
            int nEbins = int (scanEhi - scanElo); 

            //loop to find optimal beta
            for (int i=0; i < nBetaBins; i++){
              double energy = hit.GetDopplerYta(beta,s800->GetYta(),&track); 
              obj.FillHistogram(dirname,"Energy_vs_beta",nBetaBins,betaMin,betaMax,beta,4000,0,4000,energy);
              obj.FillHistogram(dirname,Form("Theta_vs_Energy_beta%f",beta),nEbins,scanElo,scanEhi,energy,100,0,3,hit.GetTheta());
              obj.FillHistogram(dirname,Form("summary_beta%f",beta),48,0,48,cryID,nEbins,scanElo,scanEhi,energy);
              beta += betaStep;
            }
          } else {
            //finer doppler corrections after finding a good beta
            double energy_b = hit.GetDoppler(BETA);
            double energy_bt = hit.GetDoppler(BETA,&track);
            double energy_bty = hit.GetDopplerYta(BETA,s800->GetYta(),&track); 
            double energy_btyd = hit.GetDopplerYta(s800->AdjustedBeta(BETA),s800->GetYta(),&track); 

            obj.FillHistogram(dirname,"prompt_gam_b",4096,0,8192,energy_b);
            obj.FillHistogram(dirname,"prompt_gam_bt",4096,0,8192,energy_bt);
            obj.FillHistogram(dirname,"prompt_gam_bty",4096,0,8192,energy_bty);
            obj.FillHistogram(dirname,"prompt_gam_btyd",4096,0,8192,energy_btyd);
            
            //PHI CORRELATION
            double phi = (TMath::TwoPi() - s800->Azita() - hit.GetPhi())*TMath::RadToDeg();
            if (phi < 0) phi += 360;

            obj.FillHistogram(dirname,"Phi_vs_Energy",4000,0,4000,energy_b,360,0,360,phi);
            obj.FillHistogram(dirname,"Phi_vs_Energy_corrected",4000,0,4000,energy_bt,360,0,360,phi);
          
            // //YTA CORRELATION
            // obj.FillHistogram(dirname,Form("Yta_vs_Energy_r%02d_c%d",ringnum,cryID),200,480,680,energy_b,200,-20,20,s800->GetYta());
            // obj.FillHistogram(dirname,Form("Yta_vs_Energy_corrected_r%02d_c%d",ringnum,cryID),200,480,680,energy_bty,200,-20,20,s800->GetYta());

            // //DTA CORRELATION
            // obj.FillHistogram(dirname,Form("Dta_vs_Energy_r%02d_c%d",ringnum,cryID),200,480,680,energy_b,50,-0.06,0.06,s800->GetDta());
            // obj.FillHistogram(dirname,Form("Dta_vs_Energy_corrected_r%02d_c%d",ringnum,cryID),200,480,680,energy_btyd,50,-0.06,0.06,s800->GetDta());

            //SUMMARY SPECTRUM
            obj.FillHistogram(dirname,"Doppler_summary",48,0,48,detMapRing[cryID],1024,0,4096,energy_btyd);

            //Spectrum
            obj.FillHistogram(dirname,"gamma_singles_corrected",2048,0,8192,energy_btyd);
            obj.FillHistogram(dirname,Form("gamma_singles_corrected_i%02d",detMapRing[cryID]),1024,0,4096,energy_btyd);
          }
          */
          double energy_corrected = hit.GetDopplerYta(s800->AdjustedBeta(BETA),s800->GetYta(),&track);
          obj.FillHistogram(dirname,"gamma_singles_beta_yta_dta_phi",4096,0,8192,energy_corrected);
          obj.FillHistogram(dirname,"gamma_singles_beta_yta_dta_phi_vs_theta",180,0,TMath::Pi(),theta,4096,0,8192,energy_corrected);
          obj.FillHistogram(dirname,"gamma_singles_beta_phi",4096,0,8192,hit.GetDoppler(BETA,&track));

          if (nInteractions > 1){
            double xi = hit.GetXi(&track);
            obj.FillHistogram(dirname,"gamma_singles_corrected",4096,0,8192,energy_corrected);

            TF1 *fDOP = new TF1("fdop","[0]*TMath::Sqrt(1-[1]*[1])/(1 - [1]*TMath::Cos(x))",0,TMath::Pi());
            TF1 *fRES = new TF1("fresolution","[0]*[1]*sin(x)/(1-[1]*cos(x))*0.036",0,TMath::Pi());
            fDOP->SetParameter(1,BETA);
            fRES->SetParameter(1,BETA);

            double CentroidRestEnergy = std::floor(core_energy*(1-BETA*TMath::Cos(0.65))/TMath::Sqrt(1-BETA*BETA));
            double ECentMax = std::ceil(core_energy*(1-BETA*TMath::Cos(2.05))/TMath::Sqrt(1-BETA*BETA));
            double EMAXCONSIDERED = 4000;

            std::vector<int> passedIP;
            while (CentroidRestEnergy < ECentMax && CentroidRestEnergy < EMAXCONSIDERED){
              CentroidRestEnergy ++;
              fDOP->SetParameter(0,CentroidRestEnergy);
              fRES->SetParameter(0,CentroidRestEnergy);
              
              for (int ip=0; ip < nInteractions; ip++){
                if (hit.GetSegmentEng(ip) > 100 && checkEnergyTheta(fDOP,fRES,core_energy,hit.GetTheta(ip))){
                  passedIP.push_back(ip);
                }
              }
              if (passedIP.empty()) continue;
              //NEW SECOND INTERACTION POINT
              double IPxi = xi;
              if (passedIP[0]+1 < nInteractions) IPxi = hit.GetXi(&track,passedIP[0],passedIP[0]+1);
              else if (passedIP[0] != 0) IPxi = hit.GetXi(&track,passedIP[0],0);
              obj.FillHistogram(dirname, "CentroidEnergy_vs_Xi",18,0,TMath::TwoPi(),IPxi, int(EMAXCONSIDERED),0,EMAXCONSIDERED,CentroidRestEnergy);
              obj.FillHistogram(dirname, "CentroidEnergy_vs_Xi_Coarse",360,0,TMath::TwoPi(),IPxi, int(EMAXCONSIDERED),0,EMAXCONSIDERED,CentroidRestEnergy);
              passedIP.clear();
            }

            for (int ip=0; ip < nInteractions; ip++)
              obj.FillHistogram(dirname, "IP_Ecore_vs_theta",40,0.6,2.1,hit.GetTheta(ip),2048,0,2048,core_energy);
            

            
            // ECORE THETA INTERACTION POINT GATES
            std::map<std::string,std::vector<int>> IPpass;
            for (auto ipgate : gates["EnergyTheta"]){
              std::string ipgname = std::string(ipgate->GetName());
              for (int ip=0; ip < nInteractions; ip++){
                if (ipgname.find(gates["outgoing"].at(ind_out)->GetName()) != std::string::npos && 
                    ipgate->IsInside(hit.GetTheta(ip),core_energy) && hit.GetSegmentEng(ip) > 100){
                  IPpass[ipgname].push_back(ip);
                }
              }
            }

            //loop over gates again
            for (auto ipgate : gates["EnergyTheta"]){
              std::string ipgname = std::string(ipgate->GetName());
              if (ipgname.find(gates["outgoing"].at(ind_out)->GetName()) == std::string::npos) continue;
              if (IPpass.count(ipgname)) {              
                //NEW FIRST INTERACTION POINT
                double IP_theta = hit.GetTheta(IPpass[ipgname][0]);
                double IP_phi = hit.GetPhi(IPpass[ipgname][0]);
                bool thetacut = IP_theta*TMath::RadToDeg() > 55 && IP_theta*TMath::RadToDeg() < 100;
                double E_pass_dop = hit.GetDopplerYta(s800->AdjustedBeta(GValue::Value("BETA")), s800->GetYta(), &track, IPpass[ipgname][0]);

                obj.FillHistogram(dirname, Form("IP_%s_npass_vs_totalPoints",ipgname.c_str()),11,1,12,nInteractions,9,1,10,(int) IPpass[ipgname].size());

                //NEW SECOND INTERACTION POINT
                double IPxi = xi;
                if (IPpass[ipgname][0]+1 < nInteractions) IPxi = hit.GetXi(&track,IPpass[ipgname][0],IPpass[ipgname][0]+1);
                else if (IPpass[ipgname][0] != 0) IPxi = hit.GetXi(&track,IPpass[ipgname][0],0);

                obj.FillHistogram(dirname, Form("IP_%s_Edop_vs_xi",ipgname.c_str()),360,0,TMath::TwoPi(),IPxi,4096,0,4096,E_pass_dop);
                if (thetacut) obj.FillHistogram(dirname, Form("IP_%s_Edop_vs_xi_thct",ipgname.c_str()),360,0,TMath::TwoPi(),IPxi,4096,0,4096,E_pass_dop);
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

  return;
}
