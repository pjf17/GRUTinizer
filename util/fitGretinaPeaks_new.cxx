#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "GH1D.h"
#include "TF1.h"
#include "TF1Sum.h"
#include "TFile.h"
#include "TH1.h"
#include "TFitResult.h"
#include "TNtuple.h"
#include "TMath.h"

// const std::string INPUT_HIST = "inCl45_S40_gated/gam_dop_sgl_prompt";
// const std::string INPUT_HIST = "inCl45_S40_gated/prompt_gamma_dop_xi_para";
// const std::string MODE = "gretsim/dopEn";
// const std::string MODE = "gretsim/dopEn_xi_para";

// const std::string INPUT_HIST = "inCl45_Ar46_gated/gam_dop_sgl_bin100";
// const std::string INPUT_HIST = "inCl45_Ar46_gated/gam_dop_sgl_bin75";
// std::string INPUT_HIST = "inCl45_Ar46_gated/prompt_gamma_dop_xi_psum";
// const std::string MODE = "gretsim/dopEn_bin100";
// const std::string MODE = "gretsim/dopEn_bin75";
// std::string MODE = "gretsim/dopEn_xi_psum";

std::string INPUT_HIST = "";
std::string MODE = "";

double getEff(double energy) {
  return (4.532*pow(energy+100.,-0.621)*10.75/8.)*(1+TMath::TanH((energy-185.1)/82.5))/2; //This is for GRETINA with 11 quads and one disabled detector. If you want accurate efficiency corrections you will need to find your own, but it's not important to the actual fitting
}

int getEnergy(const std::string &name) {
  return std::stoi(name.substr(name.find("t")+1,name.find(".")-name.find("t")-1));
}

void printResults(std::vector<double> fep_counts, std::vector<double> fep_counts_unc, std::vector<int> energies){
  
  std::vector<double> intensities;
  std::vector<double> intensities_unc;
  
  std::cout << "============================== COUNTS ==============================\n";
  for (unsigned int i=0;i<energies.size();i++){
    std::cout << energies.at(i) << "\t" << fep_counts.at(i) << "\t" << fep_counts_unc.at(i) << "\n";
    double eff = getEff(energies.at(i));
    intensities.push_back(fep_counts.at(i)/eff);
    intensities_unc.push_back(intensities.at(i)*sqrt(pow(fep_counts_unc.at(i)/fep_counts.at(i), 2.) + 0.021*0.021));
  }
  std::cout << "\n\n\n";
  
  std::cout << "============================== INTENSITIES ==============================\n";
  for (unsigned int i = 0; i < intensities.size(); i++){
    std::cout << energies.at(i) << "\t" << intensities.at(i) << "\t" << intensities_unc.at(i) << "\n";
  }
  
}

TF1 *constructBackground(std::string param_list, int nrebin) {
  //We generally use exponential or double exponential backgrounds. They are constructed here. Generally you only need the double exponential for very low energy regions of your spectra. The last two parameters (in the tanh component) model the reduced efficiency at low energies due to the thresholding on the GRETINA detectors. They should be fit to your raw, uncorrected spectra and then fixed for all of the actual fitting.
  TNtuple* t = new TNtuple("t","t","parameters");
  t->ReadFile(param_list.c_str());
  t->Draw("parameters","","goff");

  TF1 *bg;
  if(t->GetSelectedRows() == 5) {
    
    bg = new TF1("Single_Exp","([0]*TMath::Exp([1]*x))*(1+TMath::TanH((x-[2])/[3]))/2",0,10000);
    
    bg->FixParameter(2,t->GetV1()[2]);
    bg->FixParameter(3,t->GetV1()[3]);

    for (int i=0 ; i < 2; i++) {
      if (t->GetV1()[4] == i+1 || t->GetV1()[4] > 10) 
        bg->FixParameter(i,t->GetV1()[i]);
      else 
        bg->SetParameter(i,t->GetV1()[i]);
    }

    // bg->SetParameter(0,t->GetV1()[0]);
    // bg->SetParameter(1,t->GetV1()[1]);

    // if(t->GetV1()[4]) {
    //   bg->FixParameter(0,t->GetV1()[0]);
    //   bg->FixParameter(1,t->GetV1()[1]);
    //   bg->FixParameter(2,t->GetV1()[2]);
    //   bg->FixParameter(3,t->GetV1()[3]);
    // }
  }
  
  else if(t->GetSelectedRows() == 7) {
    
    bg = new TF1("Double_Exp","([0]*TMath::Exp([1]*x)+[2]*TMath::Exp([3]*x))*(1+TMath::TanH((x-[4])/[5]))/2",0,10000);
    
    bg->SetParameter(0,t->GetV1()[0]);
    bg->SetParameter(1,t->GetV1()[1]);
    bg->SetParameter(2,t->GetV1()[2]);
    bg->SetParameter(3,t->GetV1()[3]);
    bg->FixParameter(4,t->GetV1()[4]);
    bg->FixParameter(5,t->GetV1()[5]);

    if(t->GetV1()[6]) {
    bg->FixParameter(0,t->GetV1()[0]);
    bg->FixParameter(1,t->GetV1()[1]);
    bg->FixParameter(2,t->GetV1()[2]);
    bg->FixParameter(3,t->GetV1()[3]);
    bg->FixParameter(4,t->GetV1()[4]);
    bg->FixParameter(5,t->GetV1()[5]);
    }
    
  }
  else {
    bg = new TF1("Single_Exp","[0]*TMath::Exp([1]*x)",0,10000);
    bg->SetParameter(0,4348.);
    bg->SetParameter(1,-0.00484);
    std::cout << "Background model TNtuple does not have 5 or 7 rows. "
	      << "Ignoring it and using a single exponential with default values."
	      << std::endl;
  }
   
  bg->SetRange(0,10000);
  bg->SetNpx(50000/nrebin);
  
  return bg;
}


TFitResultPtr fitAllPeaks(GH1D* data_hist, TF1Sum &fullSum, const std::vector<TF1*> &fit_funcs, int fit_low_x, int fit_high_x, int nrebin) {
  // fullSum.AddRegion(3230,3260); 
  // fullSum.AddRegion(3845,3865); 

  fullSum.AddRegion(fit_low_x,fit_high_x);
  fullSum.FitInRegions();
  for (unsigned int i=0;i<fit_funcs.size();i++){
    fullSum.AddTF1(fit_funcs.at(i)); 
  }

  fullSum.GetFunc()->SetRange(0,10000);
  fullSum.GetFunc()->SetNpx(50000/nrebin);
  int count = 0;
  TFitResultPtr r;
  std::cout<<"Fitting "<<INPUT_HIST<<std::endl;
  while (1) {
    if (count == 0) 
      r = data_hist->Fit(fullSum.GetFunc(),"LSQ 0","",60,fit_high_x);
    else
      r = data_hist->Fit(fullSum.GetFunc(),"LESQ 0","",60,fit_high_x);
    std::cout << "Fit with r->Status() = " << r->Status()  << " r->IsValid() = " <<  r->IsValid() << std::endl;
    if (count >= 1 && r->IsValid()) {
      break;
    } 
    count++;
  }
  r->Print();
  std::cout<<"Done!"<<std::endl;
  // for (int i = 0; i < fullSum.GetFunc()->GetNpar(); i++) {
  //   std::cout << i << "\t" << fullSum.GetFunc()->GetParName(i) << "\t" << fullSum.GetFunc()->GetParameter(i) << "\t" << "+/-" << "\t"
	//       << fullSum.GetFunc()->GetParError(i) << "\n";
  // }
  return r;
}

void fitGretinaPeaks(std::string data_file_name, std::string output_fn, std::string peak_input, std::string bg_line_input, std::string exp_bg_params, int fit_low_x, int fit_high_x, int nrebin) {
  
  TFile *data_file = new TFile(data_file_name.c_str(), "read");
  if(data_file->IsZombie()) {
    std::cout << "Data File " << data_file_name << " does not exist!" << std::endl;
    return;
  }

  //The name of the histogram you are fitting needs to be hard-coded here, so you'll have to edit and recompile every time.  
  GH1D* data_hist((GH1D*)data_file->Get(INPUT_HIST.c_str()));

  if(!data_hist) {
    std::cout << "Histogram not properly retrieved from data file." << std::endl;
    return;
  }
  data_hist->Sumw2();
  data_hist->Rebin(nrebin);

  TNtuple* tpl = new TNtuple("input","input","energy:shift:parameter:ip:ic:fix:low:high:intlow:inthigh");
  tpl->ReadFile(peak_input.c_str());

  std::vector<int> energies;
  std::map<int,double> shifts;
  std::map<int,double> params;
  std::map<int,bool> peaks;
  std::map<int,bool> comps;
  std::map<int,bool> fix;
  std::map<int,double> llims;
  std::map<int,double> hlims;
  std::map<int,double> integral_llims;
  std::map<int,double> integral_hlims;

  float en; //energy
  float sh; //shift
  float pa; //parameter
  float ip; //include fep
  float ic; //include compt
  float fx; //fix parameter
  float pl; //parameter limit low
  float ph; //parameter limit high
  float il; //integration limit low
  float ih; //integration limit high

  tpl->SetBranchAddress("energy",&en);
  tpl->SetBranchAddress("shift",&sh);
  tpl->SetBranchAddress("parameter",&pa);
  tpl->SetBranchAddress("ip",&ip);
  tpl->SetBranchAddress("ic",&ic);
  tpl->SetBranchAddress("fix",&fx);
  tpl->SetBranchAddress("low",&pl);
  tpl->SetBranchAddress("high",&ph);
  tpl->SetBranchAddress("intlow",&il);
  tpl->SetBranchAddress("inthigh",&ih);

  std::cout << "\nPeak Info: " << "\nEnergy\tShift\tParam\tPeak\tComp\tFix\tLow\tHigh\tIntLow\tIntHigh" << std::endl;
  
  for(int i=0;i<tpl->GetEntries();i++) {

    tpl->GetEntry(i);

    if((bool)ip || (bool)ic) {
      
      energies.push_back((int)en);
      shifts[energies.back()] = (int)sh;
      params[energies.back()] = (double)pa; 
      peaks[energies.back()] = (bool)ip;
      comps[energies.back()] = (bool)ic;
      fix[energies.back()] = (bool)fx;
      llims[energies.back()] = (double)pl;
      hlims[energies.back()] = (double)ph;
      integral_llims[energies.back()] = (double)il;
      integral_hlims[energies.back()] = (double)ih;
    
      std::cout << energies.back() << "\t" << shifts[energies.back()] << "\t" << params[energies.back()] << "\t"
		            << peaks[energies.back()] << "\t" << comps[energies.back()] << "\t" << fix[energies.back()]
		            << "\t" << llims[energies.back()] << "\t" << hlims[energies.back()]
                << "\t" << integral_llims[energies.back()] << "\t" << integral_hlims[energies.back()];
      
      if(en < fit_low_x || en > fit_high_x) {
	      std::cout << "\tOutside fit range (" << fit_low_x << "," << fit_high_x << ")!";
      }
      std::cout << "\n";

      if(!fix[energies.back()] && (!peaks[energies.back()] || !comps[energies.back()])) {
	      std::cout << "\tThe " << energies.back() << " keV peak has a free scaling parameter, but either the peak("
                  << peaks[energies.back()] << ") or compton(" << comps[energies.back()] << ") is not included!!"
                  << std::endl;
      }
    }
  }

  std::cout << "\n";
  
  std::vector<GH1D*> fit_hists;
  std::vector<GH1D*> fep_hists;
  std::vector<GH1D*> com_hists;
  
  double NARROW_SCALE = 0.7;
  double WIDE_SCALE = 0.3;
  //Here we read in the simulated histograms. We use a narrow and a wide component for each histogram to more accurately simulate the GRETINA response.
  
  for(unsigned int i=0;i<energies.size();i++) {
    TFile* f = new TFile(Form("hist%d.root",energies.at(i)),"read");
    if(!f->IsZombie()) {
      if(peaks[energies.at(i)] && comps[energies.at(i)]) {
        fit_hists.push_back((GH1D*)f->Get(Form("%s",MODE.c_str())));
      }
      else if(peaks[energies.at(i)] && !comps[energies.at(i)]) {
        fit_hists.push_back((GH1D*)f->Get(Form("%s_fep",MODE.c_str())));
      }
      else if(!peaks[energies.at(i)] && comps[energies.at(i)]) {
	      fit_hists.push_back((GH1D*)f->Get(Form("%s_bg",MODE.c_str())));
      }
      fep_hists.push_back((GH1D*)(((TH1*)f->Get(Form("%s_fep",MODE.c_str())))->Clone()));
      com_hists.push_back((GH1D*)(((TH1*)f->Get(Form("%s_bg",MODE.c_str())))->Clone()));
      
      fit_hists.back()->Sumw2();

      // TFile* fw = new TFile(Form("wider/hist%d.root",energies.at(i)),"read");
      // if(!fw->IsZombie()) {
      //   // if (energies.at(i) < 900) {
      //   //   NARROW_SCALE = 0.7;
      //   //   WIDE_SCALE = 0.3;
      //   // }
      //   fit_hists.back()->Scale(NARROW_SCALE);
      //   fep_hists.back()->Scale(NARROW_SCALE);
      //   com_hists.back()->Scale(NARROW_SCALE);
	
      //   TH1* hw;
      //   TH1* hw_fep;
      //   TH1* hw_com;
      //   if(peaks[energies.at(i)] && comps[energies.at(i)]) {
      //     hw = (GH1D*)fw->Get(Form("%s",MODE.c_str()));
      //   }
      //   else if(peaks[energies.at(i)] && !comps[energies.at(i)]) {
      //     hw = (GH1D*)fw->Get(Form("%s_fep",MODE.c_str()));
      //   }
      //   else if(!peaks[energies.at(i)] && comps[energies.at(i)]) {
      //     hw = (GH1D*)fw->Get(Form("%s_bg",MODE.c_str()));
      //   }

      //   hw_fep = (GH1D*)(((TH1*)fw->Get(Form("%s_fep",MODE.c_str())))->Clone());
      //   hw_com = (GH1D*)(((TH1*)fw->Get(Form("%s_bg",MODE.c_str())))->Clone());
      //   hw->Sumw2();

      //   hw->Scale(WIDE_SCALE);
      //   hw_fep->Scale(WIDE_SCALE);
      //   hw_com->Scale(WIDE_SCALE);

      //   fit_hists.back()->Add(hw);
      //   fep_hists.back()->Add(hw_fep);
      //   com_hists.back()->Add(hw_com);
      // }
      // else {
      //   std::cout << "The histogram file for the wide component of the " << energies.at(i)
      //             << " keV peak was not found! It will not be included!" << std::endl;
      // }

      //Do not change this naming convention. It will break the next loop.
      fit_hists.back()->SetName(Form("hist%i",energies.at(i)));
      fep_hists.back()->SetName(Form("fep_hist%i",energies.at(i)));
      com_hists.back()->SetName(Form("com_hist%i",energies.at(i)));

      fit_hists.back()->GetXaxis()->SetLimits(fit_hists.back()->GetXaxis()->GetXmin()+shifts[energies.at(i)],
				              fit_hists.back()->GetXaxis()->GetXmax()+shifts[energies.at(i)]);

      fep_hists.back()->GetXaxis()->SetLimits(fep_hists.back()->GetXaxis()->GetXmin()+shifts[energies.at(i)],
				              fep_hists.back()->GetXaxis()->GetXmax()+shifts[energies.at(i)]);

      com_hists.back()->GetXaxis()->SetLimits(com_hists.back()->GetXaxis()->GetXmin()+shifts[energies.at(i)],
				              com_hists.back()->GetXaxis()->GetXmax()+shifts[energies.at(i)]);
      fit_hists.back()->Rebin(nrebin);
      fep_hists.back()->Rebin(nrebin);
      com_hists.back()->Rebin(nrebin); 
    }
    else {
      std::cout << "The histogram file for the " << energies.at(i) << " keV peak was not found!" << std::endl;
    }  
  }

  std::vector<TF1*> fit_funcs;
  for(unsigned int i=0;i<fit_hists.size();i++) {

    fit_funcs.push_back(fit_hists.at(i)->ConstructTF1());
    fit_funcs.back()->SetParameter(0,params[getEnergy(fit_hists.at(i)->GetName())]);
    
    if(fix[getEnergy(fit_hists.at(i)->GetName())]) {
      fit_funcs.back()->FixParameter(0,params[getEnergy(fit_hists.at(i)->GetName())]);
    }
    
    else {
      fit_funcs.back()->SetParLimits(0,llims[getEnergy(fit_hists.at(i)->GetName())],hlims[getEnergy(fit_hists.at(i)->GetName())]);
    }
    
    
    fit_funcs.back()->SetName(Form("func%d",getEnergy(fit_hists.at(i)->GetName()))); 
  }

  TNtuple* bg_tpl = new TNtuple("bg_input","bg_input","energy:parameter:include");
  bg_tpl->ReadFile(bg_line_input.c_str());
  bg_tpl->Draw("energy:parameter:include","","goff");

  std::vector<int> bg_energies;
  std::map<int,double> bg_params;

  std::cout << "\nStopped Line Info: " << "\nEnergy\tScale" << std::endl;
  
  for(int i=0;i<bg_tpl->GetSelectedRows();i++) {
    if(bg_tpl->GetV3()[i]) {
      
      bg_energies.push_back((int)bg_tpl->GetV1()[i]);
      bg_params[bg_energies.back()] = bg_tpl->GetV2()[i];

      std::cout << bg_energies.back() << "\t" << bg_params[bg_energies.back()];

      if(bg_energies.back() < fit_low_x || bg_energies.back() > fit_high_x) {
        std::cout << "\tOutside fit range (" << fit_low_x << "," << fit_high_x << ")!";
      }
      std::cout << std::endl;
     
    } 
  }
  std::cout << "\n";

  std::vector<GH1D*> bg_hists;
  for(unsigned int i=0;i<bg_energies.size();i++) {

    TFile* f = new TFile(Form("bg_hists/hist%d.root",bg_energies.at(i)),"read");
    if(!f->IsZombie()) {
      bg_hists.push_back((GH1D*)f->Get(Form("%s",MODE.c_str())));
      bg_hists.back()->Sumw2();
      bg_hists.back()->Scale(bg_params[bg_energies.at(i)]);
    }
    else {
      std::cout << "Histogram file for " << bg_energies.at(i) << " stopped line was not found! It will not be included."
                << std::endl;
    }
  }

  for(unsigned int i=1;i<bg_hists.size();i++) {
    bg_hists.at(0)->Add(bg_hists.at(i));
  }
      
  if(bg_hists.size() > 0) {
    bg_hists.at(0)->Rebin(nrebin);
    fit_funcs.push_back(bg_hists.at(0)->ConstructTF1());
    //fit_funcs.back()->SetParameter(0,0.01);
    fit_funcs.back()->SetParameter(0,0.002);
    // fit_funcs.back()->SetParameter(0,bg_tpl->GetV2()[bg_tpl->GetSelectedRows()-1]);
    // if((bool)bg_tpl->GetV1()[bg_tpl->GetSelectedRows()-1]) {
    //   fit_funcs.back()->FixParameter(0,bg_tpl->GetV2()[bg_tpl->GetSelectedRows()-1]);
    // }
    //else {
    //fit_funcs.back()->SetParLimits(0,0,1);
    //}
    fit_funcs.back()->SetName("StoppedLineFit");
  }

  fit_funcs.push_back(constructBackground(exp_bg_params,nrebin));
  int nbkgPar = fit_funcs.back()->GetNpar();
  GH1D* hr = (GH1D*) data_hist->Clone("Residuals");
  TF1Sum fSum, fBg; //, fBgPK;
  GH1D *hResids = (GH1D*) data_hist->Clone("Bkg Subtracted");
  TFitResultPtr res = fitAllPeaks(data_hist,fSum,fit_funcs,fit_low_x,fit_high_x,nrebin);
  data_hist = (GH1D*) data_hist->Clone(Form("%s_fit",data_hist->GetName()));
  GH1D* hTotbg;
  for(unsigned int i=0;i<com_hists.size();i++) {
    if (i == 0){
      hTotbg = (GH1D*) com_hists.at(i)->Clone();
      hTotbg->SetNameTitle("htotbg","Total Background");
      hTotbg->Scale(res->Parameter(i));
    }
    else hTotbg->Add(com_hists.at(i),res->Parameter(i));
    fBg.AddTF1(com_hists.at(i)->ConstructTF1());
  }
  if (bg_hists.size() > 0) {
    fBg.AddTF1(bg_hists.at(0)->ConstructTF1());
    hTotbg->Add(bg_hists.at(0),res->Parameter(com_hists.size()));
  }
  fBg.AddTF1(constructBackground(exp_bg_params,nrebin));
  fBg.GetFunc()->SetLineColor(kBlue);
  fBg.GetFunc()->SetParameters(fSum.GetFunc()->GetParameters());
  fBg.GetFunc()->SetNpx(10000/nrebin);
  // hResids->Add(fBg.GetFunc(),-1);
  data_hist->GetListOfFunctions()->Add(fSum.GetFunc());
  data_hist->GetListOfFunctions()->Add(fBg.GetFunc());

  //alternate counts/error estimation
  //get the exponential background info
  TMatrixDSym mat = res->GetCovarianceMatrix().GetSub(res->NPar()-nbkgPar,res->NPar()-3,res->NPar()-nbkgPar,res->NPar()-3);

  std::vector<double> ppar;
  for (int bp=0; bp < nbkgPar-2; bp++){
    if (bp%2 == 1) {
      ppar.push_back(res->Parameter(res->NPar()-nbkgPar + bp)*1000);
      mat[bp][bp] *= 1000000;
    }
    else 
      ppar.push_back(res->Parameter(res->NPar()-nbkgPar + bp));
  }

  TF1 *fitexp;
  if (nbkgPar == 4) 
    fitexp = new TF1("ffitexp","[0]*TMath::Exp([1]*x)*1000",0,10);
  else 
    fitexp = new TF1("ffitexp","([0]*TMath::Exp([1]*x)+[2]*TMath::Exp([3]*x))*1000",0,10);
  fitexp->SetParameters(&ppar[0]);

  int offset = 0;
  if(bg_hists.size()){
    // bg_hists.at(0)->Scale(fSum.GetFunc()->GetParameter(energies.size()));
    bg_hists.at(0)->SetName("Stopped_Lines");
    bg_hists.at(0)->SetTitle("Stopped_Lines");
    offset=1;
  } 

  TF1* bg; //*bgup, *bgdn;
  if(nbkgPar == 4) {
    bg = new TF1("exp_bg","([0]*TMath::Exp([1]*x))*(1+TMath::TanH((x-[2])/[3]))/2",0,10000);
    // bgup = new TF1("exp_bg","([0]*TMath::Exp([1]*x))*(1+TMath::TanH((x-[2])/[3]))/2",0,10000);
    // bgdn = new TF1("exp_bg","([0]*TMath::Exp([1]*x))*(1+TMath::TanH((x-[2])/[3]))/2",0,10000);
  }
  else if (nbkgPar == 6) {
    bg = new TF1("exp_bg","([0]*TMath::Exp([1]*x)+[2]*TMath::Exp([3]*x))*(1+TMath::TanH((x-[4])/[5]))/2",0,10000);
    // bgup = new TF1("exp_bg","([0]*TMath::Exp([1]*x)+[2]*TMath::Exp([3]*x))*(1+TMath::TanH((x-[4])/[5]))/2",0,10000);
    // bgdn = new TF1("exp_bg","([0]*TMath::Exp([1]*x)+[2]*TMath::Exp([3]*x))*(1+TMath::TanH((x-[4])/[5]))/2",0,10000);
  }

  for (int i=0; i < nbkgPar; i++) {
    bg->SetParameter(i,fSum.GetFunc()->GetParameter(energies.size()+offset+i));
    // bgup->SetParameter(i,fSum.GetFunc()->GetParameter(energies.size()+offset+i) + fSum.GetFunc()->GetParError(energies.size()+offset+i));
    // bgdn->SetParameter(i,fSum.GetFunc()->GetParameter(energies.size()+offset+i) - fSum.GetFunc()->GetParError(energies.size()+offset+i));
  }
  
  bg->SetNpx(10000/nrebin);
  // bgup->SetNpx(10000/nrebin);
  // bgdn->SetNpx(10000/nrebin);
  TH1D* hbg = (TH1D*) bg->GetHistogram();
  // TH1D* hbgup = (TH1D*) bgup->GetHistogram();
  // TH1D* hbgdn = (TH1D*) bgdn->GetHistogram();
  hbg->SetName("Exp_Bkg");
  hbg->SetTitle("Exp_Bkg");
  hTotbg->Add(hbg);
  hResids->Add(hTotbg,-1);

  // TFile *fExpFitData = new TFile("expfitdata.root","recreate");
  // fitexp->Write();
  // mat.Write();
  // fExpFitData->Close();
  FILE *countOut;
  TString outfilename = Form("counts_%s.dat",INPUT_HIST.c_str());
  outfilename.ReplaceAll('/','_');
  countOut = fopen(outfilename.Data(),"w");
  bool using_stopped_lines = bg_hists.size() > 0;
  for (unsigned int i=0;i<energies.size();i++){
    double ilo = integral_llims[energies[i]];
    double ihi = integral_hlims[energies[i]];
    if (ilo == -1 || ihi == -1) continue;
    // double residIntegral = hResids->Integral(hResids->FindBin(ilo),hResids->FindBin(ihi));

    //get the integral of all background elements in the region
    double bkgCounts = 0;
    double bkgCountsErr2 = 0;
    if (using_stopped_lines) {
      double bg_integral = bg_hists[0]->Integral(bg_hists[0]->FindBin(ilo),bg_hists[0]->FindBin(ihi)); 
      bkgCounts += bg_integral*res->Parameter(energies.size());
      bkgCountsErr2 += bg_integral*bg_integral*res->Error(energies.size())*res->Error(energies.size());
      // bkgCounts += bg_integral;
      // bkgCountsErr2 += bg_integral*bg_integral*res->Error(energies.size())*res->Error(energies.size())/res->Parameter(energies.size())/res->Parameter(energies.size());
    }

    std::vector<int> grouped;
    for (unsigned int j=0;j<energies.size();j++) {
      double integral = com_hists[j]->Integral(com_hists[j]->FindBin(ilo),com_hists[j]->FindBin(ihi)); 
      if (ilo == integral_llims[energies[j]] && ihi == integral_hlims[energies[j]]) grouped.push_back(j);
      // if (i != j)
      //   integral = fit_hists[j]->Integral(fit_hists[j]->FindBin(ilo),fit_hists[j]->FindBin(ihi));
      // else 
      //   integral = com_hists[j]->Integral(com_hists[j]->FindBin(ilo),com_hists[j]->FindBin(ihi));
      if (res->Parameter(j) < 1E-6) continue;
      bkgCounts += integral*res->Parameter(j);
      bkgCountsErr2 += integral*integral*res->Error(j)*res->Error(j);
      // printf("%d %d %f %f %f %f\n",energies[i],energies[j],bkgCounts,bkgCountsErr2,res->Parameter(j),res->Error(j));
    }
    //add exp 
    // bkgCounts += fitexp->Integral(ilo/1000,ihi/1000)/nrebin;
    double tmpbkg = hbg->Integral(hbg->FindBin(ilo),hbg->FindBin(ihi));
    // double tmpbkgup = std::abs(hbgup->Integral(hbgup->FindBin(ilo),hbgup->FindBin(ihi)) - tmpbkg);
    // double tmpbkgdn = std::abs(hbgdn->Integral(hbgdn->FindBin(ilo),hbgdn->FindBin(ihi)) - tmpbkg);
    // printf("%f %f %f %f %f\n",bkgCounts,bkgCountsErr2,tmpbkg,tmpbkgup,tmpbkgdn);

    bkgCounts += tmpbkg;
    // bkgCountsErr2 += std::pow((tmpbkgup+tmpbkgdn)/2,2);
    bkgCountsErr2 += std::pow(fitexp->IntegralError(ilo/1000,ihi/1000,&ppar[0],mat.GetMatrixArray())/nrebin,2);
    double dataCounts = data_hist->Integral(data_hist->FindBin(ilo),data_hist->FindBin(ihi));
    double peakCounts = dataCounts-bkgCounts;
    double peakCountsErr = TMath::Sqrt(dataCounts + bkgCountsErr2);
    // peakCounts = residIntegral;

    if (grouped.size() > 1) {
      i += grouped.size()-1;
      double parSum = 0;
      double parSumErr2 = 0;
      for (int k=0; k < grouped.size(); k++) {
        parSum += res->GetParams()[grouped[k]];
        parSumErr2 += res->GetErrors()[grouped[k]]*res->GetErrors()[grouped[k]];
      }

      for (int k=0; k < grouped.size(); k++){
        double fracCounts = peakCounts*res->GetParams()[grouped[k]]/parSum;
        double fracCountsErr = peakCountsErr*res->GetParams()[grouped[k]]/parSum;
        // double fracCountsErr = fracCounts * TMath::Sqrt( std::pow(peakCountsErr/peakCounts,2) + std::pow(res->GetErrors()[grouped[k]]/res->GetParams()[grouped[k]],2) + parSumErr2/parSum/parSum);
        fprintf(countOut,"REGION: [%3.0f,%3.0f] PEAK: %d --> DATA: %3.0f +/- %3.0f\tBKG: %3.0f +/- %3.0f\tCOUNTS: %3.0f +/- %3.0f\n",
          ilo,ihi,energies[grouped[k]],dataCounts,TMath::Sqrt(dataCounts),bkgCounts,TMath::Sqrt(bkgCountsErr2),fracCounts,fracCountsErr);
      }
    }
    else {
      fprintf(countOut,"REGION: [%3.0f,%3.0f] PEAK: %d --> DATA: %3.0f +/- %3.0f\tBKG: %3.0f +/- %3.0f\tCOUNTS: %3.0f +/- %3.0f\n",
        ilo,ihi,energies[i],dataCounts,TMath::Sqrt(dataCounts),bkgCounts,TMath::Sqrt(bkgCountsErr2),peakCounts,peakCountsErr);
    }
  }
  fclose(countOut);

  fSum.GetFunc()->SetNpx(50000/nrebin);
  std::vector<double> fep_counts;
  std::vector<double> fep_counts_unc;
  std::vector<double> fep_subt_counts;
  std::vector<double> fep_subt_counts_unc;
  std::vector<int> ens;
  for (unsigned int i=0;i<fep_hists.size();i++){
    fep_counts.push_back(fep_hists.at(i)->Integral()*fSum.GetFunc()->GetParameter(i));
    fep_counts_unc.push_back(fep_counts.at(i)*(fSum.GetFunc()->GetParError(i)/fSum.GetFunc()->GetParameter(i)));
    ens.push_back(getEnergy(fit_hists.at(i)->GetName()));

    // double subtCounts, subtCountsUnc;
    // bkgSubtractedCounts(ens.back(),fep_hists.at(i),data_hist,fBg.GetFunc(),res,subtCounts,subtCountsUnc);
    // fep_subt_counts.push_back(subtCounts);
    // fep_subt_counts_unc.push_back(subtCountsUnc);
  }
  // printResults(fep_counts,fep_counts_unc,ens);
  
  //print updated peak info
  std::ofstream outFile("new_input_pars.txt");
  outFile<<"\n\n===========================NEW PEAK INPUT===============================\n";
  int npar = energies.size();
  for (int pp=0; pp < npar; pp++){
    outFile << energies[pp] << "\t" << shifts[energies[pp]] << "\t" << res->Parameter(pp) << "\t"
              << peaks[energies[pp]] << "\t" << comps[energies[pp]] << "\t" << fix[energies[pp]]
		            << "\t" << llims[energies[pp]] << "\t" << hlims[energies[pp]]
                << "\t" << integral_llims[energies[pp]] << "\t" << integral_hlims[energies[pp]] << std::endl;
  }
  outFile.close();

  // TH1* hf = fSum.GetFunc()->GetHistogram();
  // hf->Rebin( (int) hf->GetNbinsX()/data_hist->GetNbinsX() );
  
  // TH1D* hbg = new TH1D("Exp_Bkg","Exp_Bkg",10000/nrebin,0,10000);
  // for (int idx=0; idx < 5000000; idx++) hbg->Fill(bg->GetRandom());
  // hbg->Scale(hbg_scale->Integral()/hbg->Integral());
  
  
  // hr->Sumw2();
  hr->Add(fSum.GetFunc(),-1);
  for (int bb=0; bb <= hr->GetNbinsX(); bb++){
    hr->SetBinError(bb,TMath::Sqrt(2*data_hist->GetBinContent(bb)));
  }

  data_hist->GetXaxis()->SetTitle("Energy (keV)");
  data_hist->GetYaxis()->SetTitle("Counts");
  TFile *outfile = new TFile(output_fn.c_str(), "recreate");
  outfile->cd();
  
  data_hist->Write();
  // hf->Write();
  
  if(bg_hists.size()) {
    bg_hists.at(0)->Scale(res->Parameter(energies.size()));
    bg_hists.at(0)->Write();
  }
  
  hbg->Write();
  hResids->Write();
  hTotbg->Write();
  
  hr->Write();

  for(unsigned int i=0;i<fit_hists.size();i++) {
    
    fit_hists.at(i)->Scale(fSum.GetFunc()->GetParameter(i));
    fit_hists.at(i)->SetName(Form("Peak%02i",i));
    fit_hists.at(i)->SetTitle(Form("hist%i",energies.at(i)));
    fit_hists.at(i)->GetXaxis()->SetLimits(0,10000);
    fit_hists.at(i)->Write();

    fep_hists.at(i)->Scale(fSum.GetFunc()->GetParameter(i));
    fep_hists.at(i)->SetName(Form("FEP%02i",i));
    fep_hists.at(i)->SetTitle(Form("FEP_hist%i",energies.at(i)));
    fep_hists.at(i)->GetXaxis()->SetLimits(0,10000);
    fep_hists.at(i)->Write();

    com_hists.at(i)->Scale(fSum.GetFunc()->GetParameter(i));
    com_hists.at(i)->SetName(Form("COM%02i",i));
    com_hists.at(i)->SetTitle(Form("COM_hist%i",energies.at(i)));
    com_hists.at(i)->GetXaxis()->SetLimits(0,10000);
    com_hists.at(i)->Write();
    
  }
  
  outfile->Close();
  data_file->Close();
  return;
}


int main(int argc, char **argv) { 
  
  std::string USAGE("fitGretinaPeaks input_data_file data_hist_path sim_hist_path output_file_name peak_input bg_line_input bg_model_params x_low x_high nrebin\n");
  int nrebin = 4;

  if (argc > 9) {
    std::string input_data_file(argv[1]);
    INPUT_HIST = std::string(argv[2]);
    MODE = std::string(argv[3]);
    std::string output_fn(argv[4]);
    std::string peak_input(argv[5]);
    std::string bg_line_input(argv[6]);
    std::string bg_model_params(argv[7]);
    int fit_low_x = std::stoi(argv[8]);
    int fit_high_x = std::stoi(argv[9]);
    if (argc == 11) nrebin = std::stoi(argv[10]);
    
    fitGretinaPeaks(input_data_file,output_fn,peak_input,bg_line_input,bg_model_params,fit_low_x,fit_high_x,nrebin);
    return 0;
  }

  else {
    std::cout << "USAGE: " << USAGE;
    return -1;
  }
}
