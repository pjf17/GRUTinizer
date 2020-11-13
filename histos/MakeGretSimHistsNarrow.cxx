
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
#include "TGretSim.h"
#include "TS800.h"
#include "TS800Sim.h"
#include "GCutG.h"

#include "TChannel.h"
#include "GValue.h"
#include "TRandom.h"

class BinCenters {
    public:
        BinCenters(int bins, double low, double high)
            : bins(bins), low(low), high(high) { }

        class iterator {
            public:
                iterator(BinCenters& axis, int binnum)
                    : axis(axis), binnum(binnum) { }

                double operator*() const;

                bool operator==(const iterator& other) const {
                    return
                        (&axis == &other.axis) &&
                        (binnum == other.binnum);
                }

                bool operator!=(const iterator& other) const {
                    return !(*this == other);
                }

                iterator& operator++() {
                    binnum++;
                    return *this;
                }

            private:
                BinCenters& axis;
                int binnum;
        };
        friend class BinCenters::iterator;

        iterator begin() {
            return iterator(*this, 0);
        }

        iterator end() {
            return iterator(*this, bins);
        }

    private:
        int bins;
        double low;
        double high;
};

double BinCenters::iterator::operator*() const {
    return axis.low + (binnum+0.5) * (axis.high-axis.low)/axis.bins;
}

std::map<int,int> detMap = {{26,0}, {30,1}, {34,2}, {38,3}, {25,4}, {29,5}, {33,6}, {37,7}, {27,8}, {31,9}, {35,10}, {39,11},
			    {24,12},{28,13},{32,14},{36,15},{47,16},{51,17},{59,18},{63,19},{67,20},{71,21}, {79,22}, {44,23},
			    {50,24}, {58,25}, {60,26}, {66,27}, {68,28}, {76,29}, {46,30}, {48,31}, {56,32}, {62,33},
			    {64,34}, {70,35}, {78,36}, {45,37}, {49,38}, {57,39}, {61,40}, {65,41}, {69,42}, {77,43}};

TRandom *rand_gen = 0;
void HandleS800Sim(TRuntimeObjects &obj);
void HandleGretSim(TRuntimeObjects &obj);
void HandleGretSimAddback(TRuntimeObjects &obj, const std::string &dirname, TGretina *gretina_ab,
                               const TVector3 &track_shifted, double yta);
void HandleDopplerCorrectionHists (TRuntimeObjects &obj, TGretina *gretina, const std::string &dirname);
void HandleMultiplicityAndMatrices(TRuntimeObjects &obj, const std::string &dirname, 
                                             const std::vector<double> &gret_energies);
void HandleS800Sim(TRuntimeObjects &obj){
  TS800Sim    *s800sim    = obj.GetDetector<TS800Sim>();
  if(!s800sim || !s800sim->Size()){
    return;
  }

  std::string dirname("s800sim");
  obj.FillHistogram(dirname,"ata", 600,-100,100, s800sim->GetS800SimHit(0).GetATA());
  obj.FillHistogram(dirname,"bta", 600,-100,100, s800sim->GetS800SimHit(0).GetBTA());
  obj.FillHistogram(dirname,"yta", 1000,-3,3, s800sim->GetS800SimHit(0).GetYTA());
  obj.FillHistogram(dirname,"dta", 1000,-0.5,0.5, s800sim->GetS800SimHit(0).GetDTA());

}

void HandleGretSimAddback(TRuntimeObjects &obj, const std::string &dirname, TGretina *gretina_ab){

  TS800Sim    *s800sim    = obj.GetDetector<TS800Sim>();
  TGretSim    *gretsim    = obj.GetDetector<TGretSim>();
  double yta = s800sim->GetS800SimHit(0).GetYTA();
  TVector3 track = s800sim->Track(0, 0);

  std::vector<double> ab_hits_energy;
  for(int x=0;x<gretina_ab->AddbackSize();x++) {
    const TGretinaHit &ab_hit = gretina_ab->GetAddbackHit(x);
    double energy_track_yta_dta = ab_hit.GetDopplerYta(s800sim->AdjustedBeta(GValue::Value("BETA")), yta, &track);
    int number = ab_hit.GetNumber();
    ab_hits_energy.push_back(energy_track_yta_dta);
    obj.FillHistogram(dirname,"gretina_ab_summary_allcorr", 50, 0, 50, number,
        10000,0,10000, energy_track_yta_dta);
    obj.FillHistogram(dirname,"gretina_ab_allcorr", 10000,0,10000, energy_track_yta_dta);

    if (gretsim->GetGretinaSimHit(0).fIsFull){ //full energy peak event
      obj.FillHistogram(dirname,"gretina_ab_allcorr_fep", 10000,0,10000, energy_track_yta_dta);
    }
    else{
      obj.FillHistogram(dirname,"gretina_ab_allcorr_bg", 10000,0,10000, energy_track_yta_dta);
    }

  }//loop over addback hits


  std::string dirname_dc(dirname);
  dirname_dc += "_doppler";
  HandleDopplerCorrectionHists(obj, gretina_ab, dirname_dc);

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

  for (unsigned int i = 0; i < gret_energies.size(); i++){
    for (unsigned int j = i+1; j < gret_energies.size(); j++){
      obj.FillHistogram(dirname,"gg_matrix", 2500,0,10000, gret_energies.at(i),
                                          2500,0,10000, gret_energies.at(j));
      obj.FillHistogram(dirname,"gg_matrix", 2500,0,10000, gret_energies.at(j),
                                             2500,0,10000, gret_energies.at(i));



      if(gret_energies.size() == 2){
        obj.FillHistogram(dirname,"mult2_gg_matrix", 2500,0,10000,gret_energies.at(i),
                                            2500,0,10000,gret_energies.at(j));
        obj.FillHistogram(dirname,"mult2_gg_matrix", 2500,0,10000,gret_energies.at(j),
                                            2500,0,10000,gret_energies.at(i));   
      }
    }
  }

}

void HandleDopplerCorrectionHists (TRuntimeObjects &obj, TGretina *gretina, const std::string &dirname){
  
  TS800Sim    *s800sim    = obj.GetDetector<TS800Sim>();
  double ata = s800sim->GetS800SimHit(0).GetATA();
  //std::cout << "ATA: " << ata << std::endl;
  double bta = s800sim->GetS800SimHit(0).GetBTA();
  //std::cout << "BTA: " << bta << std::endl;
  double yta = s800sim->GetS800SimHit(0).GetYTA();
  double dta = s800sim->GetS800SimHit(0).GetDTA();
  double azita = s800sim->Azita(ata,bta)*TMath::RadToDeg();//GRUTinizer expects YTA in mm
  TVector3 track = s800sim->Track(0, 0);

  std::string dirname_ytadbd(dirname);
  dirname_ytadbd += "_ytadetbydet";
  std::string dirname_dtadbd(dirname);
  dirname_dtadbd += "_dtadetbydet";
  std::string histname;

  double beta = GValue::Value("BETA");
  //double beta = 0.408;
  
  for(unsigned int x=0;x<gretina->Size();x++) {
    
    //int number = gretina->GetGretinaHit(x).GetNumber();
    int number = detMap[gretina->GetGretinaHit(x).GetCrystalId()];
    double energy = gretina->GetGretinaHit(x).GetDoppler(beta);
    double energy_track = gretina->GetGretinaHit(x).GetDoppler(beta, &track);
    double energy_track_yta = gretina->GetGretinaHit(x).GetDopplerYta(beta, yta, &track);
    double energy_track_yta_dta = gretina->GetGretinaHit(x).GetDopplerYta(s800sim->AdjustedBeta(beta),yta,&track);

    obj.FillHistogram("doppler_corrections","nocorr_summary",36,1,37,number,8000,0,4000,energy);
    obj.FillHistogram("doppler_corrections","trackcorr_summary",36,1,37,number,8000,0,4000,energy_track);
    obj.FillHistogram("doppler_corrections","trackytacorr_summary",36,1,37,number,8000,0,4000,energy_track_yta);
    obj.FillHistogram("doppler_corrections","allcorr_summary",36,1,37,number,4000,0,8000,energy_track_yta_dta);

    double reac_phi = (track.Cross(TVector3(0.0,0.0,1.0))).Phi();
    if(reac_phi < 0)
      {reac_phi += TMath::TwoPi();}

    double det_phi = ((gretina->GetGretinaHit(x).GetPosition()).Cross(TVector3(0.0,0.0,1.0))).Phi();
    if(det_phi < 0)
      {det_phi += TMath::TwoPi();}

    double phi1 = reac_phi - det_phi;
    if(phi1 < 0)
    {phi1 += TMath::TwoPi();}
    
    double gret_phi = gretina->GetGretinaHit(x).GetPhiDeg();
    double plane_angle = (360. - azita) - gret_phi;

    if (plane_angle < 0){
      plane_angle += 360;
    }

    //First the Phi Correction Plots
    obj.FillHistogram(dirname,"before_phi_corr",360,0,360,plane_angle,5000,0,2500,energy);
    obj.FillHistogram(dirname,"after_phi_corr",360,0,360,plane_angle,5000,0,2500,energy_track);

    obj.FillHistogram(dirname,"before_phi_corr1",360,0,360,phi1*TMath::RadToDeg(),5000,0,2500,energy);
    obj.FillHistogram(dirname,"after_phi_corr1",360,0,360,phi1*TMath::RadToDeg(),5000,0,2500,energy_track);
    
    //Then the Yta correction plot after Phi
    histname = Form("trackcorr_yta_num%02d",number);
    obj.FillHistogram(dirname_ytadbd,histname,60,1200,1320,energy_track,60,-15,15,yta);
    
    histname = Form("trackytacorr_yta_num%02d",number);
    obj.FillHistogram(dirname_ytadbd,histname,60,1200,1320,energy_track_yta,60,-15,15,yta);

    //Then the Dta correction plots after Yta+Phi
    histname = Form("trackytacorr_dta_num%02d",number);
    obj.FillHistogram(dirname_dtadbd,histname, 60,1200,1320,energy_track_yta,48,-0.06,0.06, dta);
    
    histname = Form("trackytadtacorr_dta_num%02d",number);
    obj.FillHistogram(dirname_dtadbd,histname, 60,1200,1320,energy_track_yta_dta,48,-0.06,0.06,dta);
  }//loop over gretina hits
}
void HandleGretSim(TRuntimeObjects &obj){

  TGretina *gretina = obj.GetDetector<TGretina>();
  TS800Sim    *s800sim    = obj.GetDetector<TS800Sim>();
  TGretSim *gretsim = obj.GetDetector<TGretSim>();
  bool stopped = false;

  if (!gretina){
    return;
  }
  if (!s800sim || !s800sim->Size()){
    stopped = true;
  }
  TGretina *gretina_ab = new TGretina(*gretina);
  gretina_ab->Clear();

  std::string dirname("gretsim");
  TVector3 track;
  double yta;
  if (!stopped){
    track = s800sim->Track(0,0);
    yta = s800sim->GetS800SimHit(0).GetYTA();
  }
  //for gg_matrix and multiplicity spectra
  std::vector<double> nonab_hits_energy;

  double beta = GValue::Value("BETA");
  //double SIGMA = 7.0;
  //double SIGMA = 3.0;
  
  int En = GValue::Value("FEP_EN"); //keV
  //if(En==471 || En == 547 || En ==680 || En==1873 || En==2157) {
  //double SIGMA = 2.086*TMath::Exp(-0.0574*En/1000.0) + 60.0277*TMath::Exp(-10.2297*En/1000.0);
  //double SIGMA = (2.1*TMath::Exp(-0.1*En/1000.0) + 60.0*TMath::Exp(-10.2*En/1000.0))*0.5;
  double SIGMA = (2.1*TMath::Exp(-0.1*En/1000.0) + 60.0*TMath::Exp(-10.2*En/1000.0));
  if(SIGMA > 3.8) {
    SIGMA = 3.8;
  }
  if(En < 150.) {
    SIGMA = 6.0;
  }
  //}
  
  for(unsigned int x = 0; x < gretina->Size(); x++){
    TVector3 local_pos(gretina->GetGretinaHit(x).GetLocalPosition(0));

    double smear_x = local_pos.X() + rand_gen->Gaus(0, SIGMA);
    double smear_y = local_pos.Y() + rand_gen->Gaus(0, SIGMA);
    double smear_z = local_pos.Z() + rand_gen->Gaus(0, SIGMA);
    gretina->GetGretinaHit(x).SetPosition(0,smear_x, smear_y, smear_z);


    gretina_ab->InsertHit(gretina->GetGretinaHit(x));
    //int number = gretina->GetGretinaHit(x).GetNumber();
    int number = detMap[gretina->GetGretinaHit(x).GetCrystalId()];
    double energy_track_yta_dta;
    double energy_track_yta;
    double energy_track;
    if (!stopped){
      energy_track = gretina->GetGretinaHit(x).GetDoppler(beta, &track);
      energy_track_yta = gretina->GetGretinaHit(x).GetDopplerYta(beta, yta, &track);
      energy_track_yta_dta = gretina->GetGretinaHit(x).GetDopplerYta(s800sim->AdjustedBeta(beta), yta, &track);
    }
    else{
      //energy_track = gretina->GetGretinaHit(x).GetDoppler(beta, &track);
      //energy_track_yta = gretina->GetGretinaHit(x).GetDopplerYta(beta, yta, &track);
      energy_track = gretina->GetGretinaHit(x).GetDoppler(beta);
      energy_track_yta = gretina->GetGretinaHit(x).GetDoppler(beta);
      energy_track_yta_dta = gretina->GetGretinaHit(x).GetDoppler(beta);
    }
    nonab_hits_energy.push_back(energy_track_yta_dta);
    //double EfficiencyCorrection = (1+TMath::TanH((gretina->GetGretinaHit(x).GetCoreEnergy()-185.1)/82.5))/2;
    //double EfficiencyCorrection = (1+TMath::TanH((gretina->GetGretinaHit(x).GetCoreEnergy()-47.4)/5.3))/2;
    obj.FillHistogram(dirname,"HitTheta_v_DetMap",100,0,100,detMap[gretina->GetGretinaHit(x).GetCrystalId()],
                                                  226,0,3.2,gretina->GetGretinaHit(x).GetTheta());
    obj.FillHistogram(dirname,"HitPhi_v_DetMap",100,0,100,detMap[gretina->GetGretinaHit(x).GetCrystalId()],
                                                  452,0,6.3,gretina->GetGretinaHit(x).GetPhi());
    double thresh_param1 = GValue::Value(Form("DET%i_THRESH1",detMap[gretina->GetGretinaHit(x).GetCrystalId()]));
    double thresh_param2 = GValue::Value(Form("DET%i_THRESH2",detMap[gretina->GetGretinaHit(x).GetCrystalId()]));
    double EfficiencyCorrection = (1+TMath::TanH((gretina->GetGretinaHit(x).GetCoreEnergy()-thresh_param1)/thresh_param2))/2;
    double_t rDraw = rand_gen->Uniform();
    if(rDraw <= EfficiencyCorrection) {

    obj.FillHistogram(dirname,"CoreEnergy",10000,0,10000,gretina->GetGretinaHit(x).GetCoreEnergy());
    
    //obj.FillHistogram(dirname,"HitTheta_v_DetMap",36,1,37,number,
    //		      360,0,360,gretina->GetGretinaHit(x).GetTheta()*TMath::RadToDeg());
    //obj.FillHistogram(dirname,"HitPhi_v_DetMap",36,1,37,number,
    //		      360,0,360,gretina->GetGretinaHit(x).GetPhiDeg());

    obj.FillHistogram(dirname,"HitTheta_v_HitPhi",
    		      360,0,360,gretina->GetGretinaHit(x).GetTheta()*TMath::RadToDeg(),
    		      360,0,360,gretina->GetGretinaHit(x).GetPhi()*TMath::RadToDeg());

    obj.FillHistogram(dirname,"gretina_summary_B&T",36,1,37,number,8000,0,4000,energy_track);
    obj.FillHistogram(dirname,"gretina_summary_B&T&Y&D",36,1,37,number,8000,0,4000,energy_track_yta_dta);
    //obj.FillHistogram(dirname,"gretina_summary_B&T&Y",36,1,37,number,8000,0,4000,energy_track_yta);
    //obj.FillHistogram(dirname,"gretina_summary_allcorr",36,1,37,number,8000,0,4000,energy_track_yta_dta);
    
    //obj.FillHistogram(dirname,"gretina_allcorr",10000,0,10000,energy_track_yta_dta);
    obj.FillHistogram(dirname,"gretina_B&T",10000,0,10000,energy_track);
    obj.FillHistogram(dirname,"gretina_B&T&Y",10000,0,10000,energy_track_yta);
    obj.FillHistogram(dirname,"gretina_B&T&Y&D",10000,0,10000,energy_track_yta_dta);
    //obj.FillHistogram(dirname,"gretina_B&T&Y",10000,0,10000,energy_track_yta);

    if(detMap[gretina->GetGretinaHit(x).GetCrystalId()] < 17) {
      obj.FillHistogram(dirname,"gretina_B&T_Fwd",10000,0,10000,energy_track);

      //obj.FillHistogram(dirname,"HitTheta_v_HitPhi_Fwd",
      //	      360,0,360,gretina->GetGretinaHit(x).GetTheta()*TMath::RadToDeg(),
      //	      360,0,360,gretina->GetGretinaHit(x).GetPhi()*TMath::RadToDeg());
    }
    else {
      obj.FillHistogram(dirname,"gretina_B&T_90Deg",10000,0,10000,energy_track);

      //obj.FillHistogram(dirname,"HitTheta_v_HitPhi_90Deg",
      //	      360,0,360,gretina->GetGretinaHit(x).GetTheta()*TMath::RadToDeg(),
      //	      360,0,360,gretina->GetGretinaHit(x).GetPhi()*TMath::RadToDeg());
    }
    /*
    for(double b : BinCenters(100,0.30,0.40)) {
    obj.FillHistogram(dirname,"gretina_B&T_BetaScan",100,0.30,0.40,b,
    			4000,4000,8000,gretina->GetGretinaHit(x).GetDoppler(b,&track));
    }
    */
    if (gretsim->GetGretinaSimHit(0).fIsFull){ //full energy peak event
      
      //obj.FillHistogram(dirname,"gretina_allcorr_fep",10000,0,10000,energy_track_yta_dta);
      //obj.FillHistogram(dirname,"gretina_B&T&Y_fep",10000,0,10000,energy_track_yta);
      obj.FillHistogram(dirname,"gretina_B&T_fep",10000,0,10000,energy_track);
      
      if(detMap[gretina->GetGretinaHit(x).GetCrystalId()] < 17) {
	obj.FillHistogram(dirname,"gretina_B&T_fep_Fwd",10000,0,10000,energy_track);
      }
      else {
	obj.FillHistogram(dirname,"gretina_B&T_fep_90Deg",10000,0,10000,energy_track);
      }
      
    }
    else{
      
      //obj.FillHistogram(dirname,"gretina_allcorr_bg",10000,0,10000,energy_track_yta_dta);
      //obj.FillHistogram(dirname,"gretina_B&T&Y_bg",10000,0,10000,energy_track_yta);
      obj.FillHistogram(dirname,"gretina_B&T_bg",10000,0,10000,energy_track);

      if(detMap[gretina->GetGretinaHit(x).GetCrystalId()] < 17) {
	obj.FillHistogram(dirname,"gretina_B&T_bg_Fwd",10000,0,10000,energy_track);
      }
      else {
	obj.FillHistogram(dirname,"gretina_B&T_bg_90Deg",10000,0,10000,energy_track);
      }
    }
  }//Random draw efficiency correction
  }//loop over gretina hits

//std::string dirname_ab(dirname);
//dirname_ab += "_addback";
//HandleGretSimAddback(obj, dirname_ab, gretina_ab); 

//std::string dirname_dc(dirname);
//dirname_dc += "_doppler";
//HandleDopplerCorrectionHists(obj, gretina, dirname_dc);
//
//std::string dirname_mgg(dirname);
//dirname_mgg += "_multcuts_ggmatrices";
//HandleMultiplicityAndMatrices(obj, dirname_mgg, nonab_hits_energy);
  return;
}

// extern "C" is needed to prevent name mangling.
// The function signature must be exactly as shown here,
//   or else bad things will happen.
extern "C"
void MakeHistograms(TRuntimeObjects& obj) {
//  InitMap();
  TGretina *gretina = obj.GetDetector<TGretina>();
  //TS800Sim    *s800sim    = obj.GetDetector<TS800Sim>();
  TList    *list    = &(obj.GetObjects());
  int numobj = list->GetSize();

  if (!gretina){
    return;
  }

  if (!rand_gen){
    rand_gen = new TRandom(50747227);
  }
  std::string dirname  = "";
  
  HandleS800Sim(obj);
  HandleGretSim(obj);

  if(numobj!=list->GetSize()){
    list->Sort();
  }
}
