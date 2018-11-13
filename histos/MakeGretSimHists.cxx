
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
  obj.FillHistogram(dirname,"dta", 1000,-0.07,0.07, s800sim->GetS800SimHit(0).GetDTA());

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
  double bta = s800sim->GetS800SimHit(0).GetBTA();
  double yta = s800sim->GetS800SimHit(0).GetYTA();
  double dta = s800sim->GetS800SimHit(0).GetDTA();
  double azita = s800sim->Azita(ata,bta);//GRUTinizer expects YTA in mm
  TVector3 track = s800sim->Track(0, 0);

  std::string dirname_ytadbd(dirname);
  dirname_ytadbd += "_ytadetbydet";
  std::string dirname_dtadbd(dirname);
  dirname_dtadbd += "_dtadetbydet";
  std::string histname;
  for(unsigned int x=0;x<gretina->Size();x++) {
    int number = gretina->GetGretinaHit(x).GetNumber();
    double energy = gretina->GetGretinaHit(x).GetDoppler(GValue::Value("BETA"));
    double energy_track = gretina->GetGretinaHit(x).GetDoppler(GValue::Value("BETA"), &track);
    double energy_track_yta = gretina->GetGretinaHit(x).GetDopplerYta(GValue::Value("BETA"), yta, &track);
    double energy_track_yta_dta = gretina->GetGretinaHit(x).GetDopplerYta(s800sim->AdjustedBeta(GValue::Value("BETA")), yta, &track);

    obj.FillHistogram("doppler_corrections","allcorr_summary", 50, 0, 50, number,
        4000,0,4000, energy_track_yta_dta);
    double gret_phi = gretina->GetGretinaHit(x).GetPhiDeg();
    double plane_angle = (360. - azita) - gret_phi;

    if (plane_angle < 0){
      plane_angle += 360;
    }

    //First the Phi Correction Plots
    obj.FillHistogram(dirname,"before_phi_corr", 360,0, 360,plane_angle,
        200,1160,1360, energy);
    obj.FillHistogram(dirname,"after_phi_corr", 360, 0, 360,plane_angle,
        200,1160,1360, energy_track);
    //Then the Yta correction plot after Phi
    histname = Form("trackcorr_yta_num%02d",number);
    obj.FillHistogram(dirname_ytadbd,histname,60,1200,1320, energy_track, 
        60, -15, 15, yta);
    histname = Form("trackytacorr_yta_num%02d",number);
    obj.FillHistogram(dirname_ytadbd,histname,60,1200,1320, energy_track_yta, 
        60, -15, 15, yta);

    //Then the Dta correction plots after Yta+Phi
    histname = Form("trackytacorr_dta_num%02d",number);
    obj.FillHistogram(dirname_dtadbd,histname, 60,1200,1320, energy_track_yta,
        48,-0.06,0.06, dta);
    histname = Form("trackytadtacorr_dta_num%02d",number);
    obj.FillHistogram(dirname_dtadbd,histname, 60,1200,1320, energy_track_yta_dta,
        48,-0.06,0.06, dta);
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
  //stopped=true;
  TGretina *gretina_ab = new TGretina(*gretina);
  gretina_ab->Clear();

  std::string dirname("gretsim");
  TVector3 track;
  double yta;
  if (!stopped){
    track = s800sim->Track(0, 0);
    yta = s800sim->GetS800SimHit(0).GetYTA();
  }
  //for gg_matrix and multiplicity spectra
  std::vector<double> nonab_hits_energy;
  double SIGMA = GValue::Value("SIGMA");//mm
  std::vector<double> betas;
  const double START_BETA = 0.300;
  const double END_BETA = 0.360;
  const double STEP_SIZE = 0.0001;
  const int N_BETA_BINS = static_cast<int>((END_BETA-START_BETA)/STEP_SIZE);

  double cur_beta = START_BETA;
  while (cur_beta <= END_BETA){
    betas.push_back(cur_beta);
    cur_beta += STEP_SIZE;
  }
  for(unsigned int x = 0; x < gretina->Size(); x++){
    TVector3 local_pos(gretina->GetGretinaHit(x).GetLocalPosition(0));


    double smear_x = local_pos.X() + rand_gen->Gaus(0, SIGMA);
    double smear_y = local_pos.Y() + rand_gen->Gaus(0, SIGMA);
    double smear_z = local_pos.Z() + rand_gen->Gaus(0, SIGMA);
    gretina->GetGretinaHit(x).SetPosition(0,smear_x, smear_y, smear_z);


    gretina_ab->InsertHit(gretina->GetGretinaHit(x));
    int number = gretina->GetGretinaHit(x).GetNumber();
    double energy_track_yta_dta;
    if (!stopped){
      energy_track_yta_dta = gretina->GetGretinaHit(x).GetDopplerYta(s800sim->AdjustedBeta(GValue::Value("BETA")), yta, &track);
    }
    else{
      energy_track_yta_dta = gretina->GetGretinaHit(x).GetDoppler(GValue::Value("BETA"));
    }
    nonab_hits_energy.push_back(energy_track_yta_dta);
    obj.FillHistogram(dirname,"gretina_summary_allcorr", 50, 0, 50, number,
        10000,0,10000, energy_track_yta_dta);
    obj.FillHistogram(dirname,"gretina_allcorr", 10000,0,10000, energy_track_yta_dta);
    if (gretsim->GetGretinaSimHit(0).fIsFull){ //full energy peak event

      obj.FillHistogram(dirname,"gretina_allcorr_fep", 10000,0,10000, energy_track_yta_dta);
      for (unsigned int  i = 0; i < betas.size(); i++){
        if (!stopped){
          double cur_energy = gretina->GetGretinaHit(x).GetDopplerYta(s800sim->AdjustedBeta(betas.at(i)), yta, &track);
          obj.FillHistogram(dirname, "gretina_beta_scan", 
              N_BETA_BINS, START_BETA, END_BETA, betas.at(i),
              10000, 0, 10000, cur_energy);
        }
      }
    }
    else{
      obj.FillHistogram(dirname,"gretina_allcorr_bg", 10000,0,10000, energy_track_yta_dta);
    }
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
  TS800Sim    *s800sim    = obj.GetDetector<TS800Sim>();
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
