
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
  int trig_bit = -1;
  for(int j=0;j<16;j++) {
    if(((bits>>j)&0x0001)){
      obj.FillHistogram("ungated","trig_bit",20,0,20,j);
      trig_bit = j;
    }
  }

  bool gamma_gate = false;
  if (gretina) {
    int nGret = gretina->Size();
    for (int i=0 ; i < nGret; i++){
      if (gretina->GetGretinaHit(i).GetCoreEnergy() > 100) {
        gamma_gate = true;
        break;
      }
    }
  }
  
  //MAKE RAW TOF HISTS
  // double raw_obj = s800->GetRawOBJ_MESY();
  // double raw_e1 = s800->GetRawE1_MESY();
  // double raw_xf = s800->GetRawXF_MESY();
  double raw_obje1 = s800->GetTof().GetTacOBJ();
  double raw_xfpe1 = s800->GetTof().GetTacXFP();
  
  obj.FillHistogram("ungated", "RAW_TOF_OBJE1", 5000, 0, 5000, raw_obje1);
  obj.FillHistogram("ungated", "RAW_TOF_XFE1", 5000, 0, 5000, raw_xfpe1);

  //MAKE INCOMING PID
  // double tof_obje1 = s800->GetMTof().GetCorrelatedObjE1(); 
  // double tof_xfpe1 = s800->GetMTof().GetCorrelatedXfpE1();
  obj.FillHistogram("ungated", "incoming_pid", 3000, 0, 3000, raw_obje1, 2000, 0, 2000, raw_xfpe1);                                              
  if (gamma_gate) obj.FillHistogram("ungated", "incoming_pid_gammagate", 3000, 0, 3000, raw_obje1, 2000, 0, 2000, raw_xfpe1);                                              
  obj.FillHistogram("ungated", Form("incoming_pid_trigbit%d",trig_bit), 3000, 0, 3000, raw_obje1, 2000, 0, 2000, raw_xfpe1);                                              

  //MAKE OUTGOING PID
  double ic_ave = s800->GetIonChamber().GetAve();
  double tof_obje1_corr = s800->GetTofE1_TAC(GValue::Value("OBJ_CORR_AFP"),GValue::Value("OBJ_CORR_XFP"));
  double tof_xfpe1_corr = s800->GetTof().GetTacXFP() + GValue::Value("XFP_CORR_AFP") * s800->GetAFP() +  GValue::Value("XFP_CORR_XFP") * s800->GetCrdc(0).GetDispersiveX();
  obj.FillHistogram("ungated", "outgoing_pid_uncorrected", 3000, 0, 3000, raw_obje1, 1024, 0, 2048, ic_ave);
  obj.FillHistogram("ungated", "outgoing_pid_obj", 3000, 0, 3000, tof_obje1_corr, 1500, 0, 1500, ic_ave);
  obj.FillHistogram("ungated", "outgoing_pid_xfp", 2000, 0, 2000, tof_xfpe1_corr, 1500, 0, 1500, ic_ave);

  //CORRELATION PLOTS
  double crdc_1_x = s800->GetCrdc(0).GetDispersiveX();
  double crdc_2_x = s800->GetCrdc(1).GetDispersiveX();
  double afp = GetAfp(crdc_1_x, crdc_2_x);
  obj.FillHistogram("ungated", "corrobje1_crdc1x", 3000, 0, 3000, tof_obje1_corr, 600, -300, 300, crdc_1_x);
  obj.FillHistogram("ungated", "corrobje1_afp", 3000, 0, 3000, tof_obje1_corr, 1000, -0.1, 0.1, afp);
  
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
  
  //---------------------------------------------------------------
  //GATED
  std::vector<unsigned short> incoming_passed;
  CheckGates(gates["incoming"],incoming_passed,raw_obje1,raw_xfpe1);

  std::vector<unsigned short> outgoing_passed;
  CheckGates(gates["outgoing"],outgoing_passed,tof_xfpe1_corr,ic_ave);
  
  std::string dirname  = "";

  bool nodirbeam = false; 
  // if (gates.count("dirbeam")) nodirbeam = !gates["dirbeam"][0]->IsInside(tof_xfpe1_corr,crdc_1_x);

  // bool xtofgate = false; 
  // if (gates.count("dirbeam")) xtofgate = gates["crdcx1_tof"][0]->IsInside(tof_xfpe1_corr,crdc_1_x);

  if (gamma_gate){
    for (auto ind_in : incoming_passed){
      dirname = Form("%s_gated", gates["incoming"].at(ind_in)->GetName());
      //time energy spectrum
      if (bank29){
        int gSize = gretina->Size();
        double timeBank29 = bank29->Timestamp();
        for (int i=0; i < gSize; i++){
          TGretinaHit &hit = gretina->GetGretinaHit(i);
          obj.FillHistogram(dirname,"core_energy_t0_Bank29_time",
              3200,-400,400,timeBank29-hit.GetTime(),
              2000,0,10000,hit.GetCoreEnergy());
        }
      }
      
      for (int rep=0; rep < 2; rep++){
        dirname = Form("%s_gated", gates["incoming"].at(ind_in)->GetName());
        if (rep == 1 && nodirbeam) dirname = Form("%s_gated_dirBeamOut", gates["incoming"].at(ind_in)->GetName());
        else if (rep == 1) continue;
        // obj.FillHistogram(dirname, "outgoing_pid_obj", 500, 1600, 2100, tof_obje1_corr, 2500, 0, 2500, ic_ave);
        // obj.FillHistogram(dirname, "outgoing_pid_obj_uncorrected", 500, 1600, 2100, raw_obje1, 2500, 0, 2500, ic_ave);
        // obj.FillHistogram(dirname, "corrobje1_crdc1x", 3000, 0, 3000, tof_obje1_corr, 600, -300, 300, crdc_1_x);
        // obj.FillHistogram(dirname, "corrobje1_afp", 3000, 0, 3000, tof_obje1_corr, 1000, -0.1, 0.1, afp);
        // obj.FillHistogram(dirname, "obje1_crdc1x", 3000, 0, 3000, raw_obje1, 600, -300, 300, crdc_1_x);
        // obj.FillHistogram(dirname, "obje1_afp", 3000, 0, 3000, raw_obje1, 1000, -0.1, 0.1, afp);
        // obj.FillHistogram(dirname, "crdc1 X_Y", 600, -300, 300, crdc_1_x, ybins, ylow, yhigh, crdc_1_y);  
        // obj.FillHistogram(dirname, "crdc2 X_Y", 600, -300, 300, crdc_2_x, ybins, ylow, yhigh, crdc_2_y);

        obj.FillHistogram(dirname, "outgoing_pid_xfp", 3000, 0, 3000, tof_xfpe1_corr, 2500, 0, 2500, ic_ave);
        obj.FillHistogram(dirname, "outgoing_pid_xfp_uncorrected", 3000, 0, 3000, raw_xfpe1, 2500, 0, 2500, ic_ave);
        obj.FillHistogram(dirname, "corrxfpe1_crdc1x", 3000, 0, 3000, tof_xfpe1_corr, 600, -300, 300, crdc_1_x);
        obj.FillHistogram(dirname, "corrxfpe1_afp", 3000, 0, 3000, tof_xfpe1_corr, 1000, -0.1, 0.1, afp);
        obj.FillHistogram(dirname, "crdc1x_afp", 1200, -300, 300, crdc_1_x, 5000, -0.1, 0.1, afp);
        obj.FillHistogram(dirname, "crdc2x_afp", 1200, -300, 300, crdc_2_x, 5000, -0.1, 0.1, afp);
        obj.FillHistogram(dirname, "xfpe1_crdc1x", 3000, 0, 3000, raw_xfpe1, 600, -300, 300, crdc_1_x);
        obj.FillHistogram(dirname, "xfpe1_afp", 3000, 0, 3000, raw_xfpe1, 1000, -0.1, 0.1, afp);
        obj.FillHistogram(dirname, "crdc1 X_Y", 600, -300, 300, crdc_1_x, ybins, ylow, yhigh, crdc_1_y);  
        obj.FillHistogram(dirname, "crdc2 X_Y", 600, -300, 300, crdc_2_x, ybins, ylow, yhigh, crdc_2_y);
      }

      // if (!nodirbeam) continue;
      for (auto ind_out : outgoing_passed){
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
        else if (gates["crdc_tof"][tof_gate_idx]->IsInside(tof_xfpe1_corr,crdc_1_x))
          dirname = Form("%s_%s_%s_gated",gates["incoming"].at(ind_in)->GetName(),gates["outgoing"].at(ind_out)->GetName(),gates["crdc_tof"][tof_gate_idx]->GetName());
        else
          continue;
        
        obj.FillHistogram(dirname, "outgoing_pid_xfp", 3000, 0, 3000, tof_xfpe1_corr, 2500, 0, 2500, ic_ave);
        obj.FillHistogram(dirname, "corrxfpe1_crdc1x", 3000, 0, 3000, tof_xfpe1_corr, 600, -300, 300, crdc_1_x);
        obj.FillHistogram(dirname, "corrxfpe1_afp", 3000, 0, 3000, tof_xfpe1_corr, 1000, -0.1, 0.1, afp);
        obj.FillHistogram(dirname, "crdc1x_afp", 1200, -300, 300, crdc_1_x, 5000, -0.1, 0.1, afp);

        // GAMMA RAY ANALYSIS
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
            // obj.FillHistogram(dirname,Form("gretina_theta_vs_phi_r%02d",ringnum),360,0,360,phi*TMath::RadToDeg(),180,0,180,theta*TMath::RadToDeg());
            obj.FillHistogram(dirname,"core_energy_summary",48,0,48,cryID,4000,0,4000,core_energy);

            double yta = s800->GetYta();
            double dta = s800->GetDta();
            double ata = s800->GetAta();
            double bta = s800->GetBta();
            obj.FillHistogram("s800","ata",500,-0.1,0.1,ata);
            obj.FillHistogram("s800","bta",500,-0.1,0.1,bta);
            obj.FillHistogram("s800","yta",1000,-20,20,yta);
            obj.FillHistogram("s800","dta",1000,-1,1,dta);

            double energy_corrected = hit.GetDopplerYta(GValue::Value("BETA"), s800->GetYta(), &track);
            double xi = hit.GetXi();

            obj.FillHistogram(dirname, "prompt_gamma_dop",8192,0,8192,energy_corrected);

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
              //PHI CORRELATION
              // double phiTa = s800->Azita();//TMath::ATan(TMath::Sin(ata)/(TMath::Sin(bta)*-1.0));
              // if (phiTa < 0) phiTa += TMath::TwoPi();
              double phi = (TMath::TwoPi() - s800->Azita() - hit.GetPhi())*TMath::RadToDeg();
              if (phi < 0) phi += 360;

              obj.FillHistogram(dirname,"Phi_vs_Energy",4000,0,4000,energy_b,360,0,360,phi);
              obj.FillHistogram(dirname,"Phi_vs_Energy_corrected",4000,0,4000,energy_bt,360,0,360,phi);
            
              // //YTA CORRELATION
              // obj.FillHistogram(dirname,Form("Yta_vs_Energy_r%02d_c%d",ringnum,cryID),200,480,680,energy_b,200,-20,20,yta);
              // obj.FillHistogram(dirname,Form("Yta_vs_Energy_corrected_r%02d_c%d",ringnum,cryID),200,480,680,energy_bty,200,-20,20,yta);

              // //DTA CORRELATION
              // obj.FillHistogram(dirname,Form("Dta_vs_Energy_r%02d_c%d",ringnum,cryID),200,480,680,energy_b,50,-0.06,0.06,dta);
              // obj.FillHistogram(dirname,Form("Dta_vs_Energy_corrected_r%02d_c%d",ringnum,cryID),200,480,680,energy_btyd,50,-0.06,0.06,dta);

              //SUMMARY SPECTRUM
              obj.FillHistogram(dirname,"Doppler_summary",48,0,48,cryID,2000,0,2000,energy_btyd);

              //Spectrum
              obj.FillHistogram(dirname,"gamma_singles_corrected",4000,0,4000,energy_btyd);
              obj.FillHistogram(dirname,Form("gamma_singles_corrected_i%02d",detMapRing[cryID]),4000,0,4000,energy_btyd);
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

                // if (IPpass[ipgname].size() < 5){
                //   obj.FillHistogram(dirname, Form("IP_%s_xi_%dpass",ipgname.c_str(), (int) IPpass[ipgname].size()),360,0,TMath::TwoPi(),IPxi);
                // }

                // if (nInteractions < 5 && thetacut)
                //   obj.FillHistogram(dirname, Form("IP_%s_xi_%dINTPNT",ipgname.c_str(), nInteractions),360,0,TMath::TwoPi(),IPxi);

                obj.FillHistogram(dirname, Form("IP_%s_Edop_vs_xi",ipgname.c_str()),360,0,TMath::TwoPi(),IPxi,4096,0,4096,E_pass_dop);
                if (thetacut) obj.FillHistogram(dirname, Form("IP_%s_Edop_vs_xi_thct",ipgname.c_str()),360,0,TMath::TwoPi(),IPxi,4096,0,4096,E_pass_dop);
              }
            }
          }
        }

        // if (xtofgate) {
        //   obj.FillHistogram(dirname, "outgoing_pid_xfp_GATED", 3000, 0, 3000, tof_xfpe1_corr, 2500, 0, 2500, ic_ave);
        //   obj.FillHistogram(dirname, "corrxfpe1_crdc1x_GATED", 3000, 0, 3000, tof_xfpe1_corr, 600, -300, 300, crdc_1_x);
        //   obj.FillHistogram(dirname, "corrxfpe1_afp_GATED", 3000, 0, 3000, tof_xfpe1_corr, 1000, -0.1, 0.1, afp);
        //   obj.FillHistogram(dirname, "crdc1x_afp_GATED", 1200, -300, 300, crdc_1_x, 5000, -0.1, 0.1, afp);
        // } else {
        //   obj.FillHistogram(dirname, "outgoing_pid_xfpNOTGATED", 3000, 0, 3000, tof_xfpe1_corr, 2500, 0, 2500, ic_ave);
        //   obj.FillHistogram(dirname, "corrxfpe1_crdc1xNOTGATED", 3000, 0, 3000, tof_xfpe1_corr, 600, -300, 300, crdc_1_x);
        //   obj.FillHistogram(dirname, "corrxfpe1_afpNOTGATED", 3000, 0, 3000, tof_xfpe1_corr, 1000, -0.1, 0.1, afp);
        //   obj.FillHistogram(dirname, "crdc1x_afpNOTGATED", 1200, -300, 300, crdc_1_x, 5000, -0.1, 0.1, afp);
        // }
      }
    }
  }
  
  if(numobj!=list->GetSize()){
    list->Sort();
  }

  return;
}
