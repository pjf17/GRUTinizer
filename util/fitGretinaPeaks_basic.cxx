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

const std::string INPUT_HIST = "ab_prompt_red_pair";
// const std::string INPUT_HIST = "inBeam_Fe64_gated/gamma_corrected_singles_prompt";
const std::string MODE = "addback/gretina_pol_red";
// const std::string MODE = "gretsim/gretina";

const int REBIN_FACTOR = 4; //What binning do you want to use on your histograms for the fit (8000 and 10000 must be divisible by this number)

double getEff(double energy) {
  return (4.532*pow(energy+100.,-0.621)*10.75/8.)*(1+TMath::TanH((energy-185.1)/82.5))/2; //This is for GRETINA with 11 quads and one disabled detector. If you want accurate efficiency corrections you will need to find your own, but it's not important to the actual fitting
}

int getEnergy(const std::string &name) {
  return std::stoi(name.substr(4,4));
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

TF1 *constructBackground(std::string param_list) {
  //We generally use exponential or double exponential backgrounds. They are constructed here. Generally you only need the double exponential for very low energy regions of your spectra. The last two parameters (in the tanh component) model the reduced efficiency at low energies due to the thresholding on the GRETINA detectors. They should be fit to your raw, uncorrected spectra and then fixed for all of the actual fitting.
  TNtuple* t = new TNtuple("t","t","parameters");
  t->ReadFile(param_list.c_str());
  t->Draw("parameters","","goff");

  TF1 *bg;
  if(t->GetSelectedRows() == 5) {
    
    bg = new TF1("Single_Exp","([0]*TMath::Exp([1]*x))*(1+TMath::TanH((x-[2])/[3]))/2",0,10000);
    
    bg->SetParameter(0,t->GetV1()[0]);
    bg->SetParameter(1,t->GetV1()[1]);
    bg->FixParameter(2,t->GetV1()[2]);
    bg->FixParameter(3,t->GetV1()[3]);

    if(t->GetV1()[4]) {
      bg->SetParameter(0,t->GetV1()[0]);
      bg->FixParameter(1,t->GetV1()[1]);
      bg->FixParameter(2,t->GetV1()[2]);
      bg->FixParameter(3,t->GetV1()[3]);
    }
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
  bg->SetNpx(50000/REBIN_FACTOR);
  
  return bg;
}


TF1Sum fitAllPeaks(GH1D* data_hist, const std::vector<TF1*> &fit_funcs, int fit_low_x, int fit_high_x) {
  TF1Sum fullSum;
  for (unsigned int i=0;i<fit_funcs.size();i++){
    fullSum.AddTF1(fit_funcs.at(i)); 
  }

  fullSum.GetFunc()->SetRange(0,10000);
  fullSum.GetFunc()->SetNpx(50000/REBIN_FACTOR);
  
  int count = 0;
  while (1) {
    TFitResultPtr r(data_hist->Fit(fullSum.GetFunc(),"LMES","",fit_low_x,fit_high_x));
    r->Print();
    std::cout << "Fit with r->Status() = " << r->Status()  << " r->IsValid() = " <<  r->IsValid() << std::endl;
    count++;
    if (count >= 2 && r->IsValid()) {
      break;
    } 
  }
  
  for (int i = 0; i < fullSum.GetFunc()->GetNpar(); i++) {
    std::cout << fullSum.GetFunc()->GetParName(i) << "\t" << fullSum.GetFunc()->GetParameter(i) << "\t" << "+/-" << "\t"
	      << fullSum.GetFunc()->GetParError(i) << "\n";
  }
  return fullSum;
}

void fitGretinaPeaks(std::string data_file_name, std::string output_fn, std::string peak_input, std::string bg_line_input, 
		     std::string exp_bg_params, int fit_low_x, int fit_high_x) {
  
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
  data_hist->Rebin(REBIN_FACTOR);

  TNtuple* tpl = new TNtuple("input","input","energy:shift:parameter:ip:ic:fix:low:high");
  tpl->ReadFile(peak_input.c_str());

  std::vector<int> energies;
  std::map<int,double> shifts;
  std::map<int,double> params;
  std::map<int,bool> peaks;
  std::map<int,bool> comps;
  std::map<int,bool> fix;
  std::map<int,double> llims;
  std::map<int,double> hlims;

  float en; //energy
  float sh; //shift
  float pa; //parameter
  float ip; //include fep
  float ic; //include compt
  float fx; //fix parameter
  float pl; //parameter limit low
  float ph; //parameter limit high

  tpl->SetBranchAddress("energy",&en);
  tpl->SetBranchAddress("shift",&sh);
  tpl->SetBranchAddress("parameter",&pa);
  tpl->SetBranchAddress("ip",&ip);
  tpl->SetBranchAddress("ic",&ic);
  tpl->SetBranchAddress("fix",&fx);
  tpl->SetBranchAddress("low",&pl);
  tpl->SetBranchAddress("high",&ph);

  std::cout << "\nPeak Info: " << "\nEnergy\tShift\tParam\tPeak\tComp\tFix\tLow\tHigh" << std::endl;
  
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
    
      std::cout << energies.back() << "\t" << shifts[energies.back()] << "\t" << params[energies.back()] << "\t"
		            << peaks[energies.back()] << "\t" << comps[energies.back()] << "\t" << fix[energies.back()]
		            << "\t" << llims[energies.back()] << "\t" << hlims[energies.back()];
      
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
  
  const double NARROW_SCALE = 0.7;
  const double WIDE_SCALE = 0.3;
  //Here we read in the simulated histograms. We use a narrow and a wide component for each histogram to more accurately simulate the GRETINA response.

  for(unsigned int i=0;i<energies.size();i++) {

    TFile* f = new TFile(Form("hist%04d.root",energies.at(i)),"read");
    if(!f->IsZombie()) {
      if(peaks[energies.at(i)] && comps[energies.at(i)]) {
        fit_hists.push_back((GH1D*)f->Get(Form("%s_B&T&Y&D",MODE.c_str())));
      }
      else if(peaks[energies.at(i)] && !comps[energies.at(i)]) {
        fit_hists.push_back((GH1D*)f->Get(Form("%s_B&T_fep",MODE.c_str())));
      }
      else if(!peaks[energies.at(i)] && comps[energies.at(i)]) {
	      fit_hists.push_back((GH1D*)f->Get(Form("%s_B&T_bg",MODE.c_str())));
      }
      fep_hists.push_back((GH1D*)(((TH1*)f->Get(Form("%s_B&T_fep",MODE.c_str())))->Clone()));
      com_hists.push_back((GH1D*)(((TH1*)f->Get(Form("%s_B&T_bg",MODE.c_str())))->Clone()));
      
      fit_hists.back()->Sumw2();

      TFile* fw = new TFile(Form("wider/hist%04d.root",energies.at(i)),"read");
      if(!fw->IsZombie()) {

        fit_hists.back()->Scale(NARROW_SCALE);
        fep_hists.back()->Scale(NARROW_SCALE);
        com_hists.back()->Scale(NARROW_SCALE);
	
        TH1* hw;
        TH1* hw_fep;
        TH1* hw_com;
        if(peaks[energies.at(i)] && comps[energies.at(i)]) {
          hw = (GH1D*)fw->Get(Form("%s_B&T&Y&D",MODE.c_str()));
        }
        else if(peaks[energies.at(i)] && !comps[energies.at(i)]) {
          hw = (GH1D*)fw->Get(Form("%s_B&T_fep",MODE.c_str()));
        }
        else if(!peaks[energies.at(i)] && comps[energies.at(i)]) {
          hw = (GH1D*)fw->Get(Form("%s_B&T_bg",MODE.c_str()));
        }

        hw_fep = (GH1D*)(((TH1*)fw->Get(Form("%s_B&T_fep",MODE.c_str())))->Clone());
        hw_com = (GH1D*)(((TH1*)fw->Get(Form("%s_B&T_bg",MODE.c_str())))->Clone());
        hw->Sumw2();

        hw->Scale(WIDE_SCALE);
        hw_fep->Scale(WIDE_SCALE);
        hw_com->Scale(WIDE_SCALE);

        fit_hists.back()->Add(hw);
        fep_hists.back()->Add(hw_fep);
        com_hists.back()->Add(hw_com);
      }
      else {
        std::cout << "The histogram file for the wide component of the " << energies.at(i)
                  << " keV peak was not found! It will not be included!" << std::endl;
      }

      //Do not change this naming convention. It will break the next loop.
      fit_hists.back()->SetName(Form("hist%04i",energies.at(i)));
      fep_hists.back()->SetName(Form("fep_hist%04i",energies.at(i)));
      com_hists.back()->SetName(Form("com_hist%04i",energies.at(i)));

      fit_hists.back()->GetXaxis()->SetLimits(fit_hists.back()->GetXaxis()->GetXmin()+shifts[energies.at(i)],
				              fit_hists.back()->GetXaxis()->GetXmax()+shifts[energies.at(i)]);

      fep_hists.back()->GetXaxis()->SetLimits(fep_hists.back()->GetXaxis()->GetXmin()+shifts[energies.at(i)],
				              fep_hists.back()->GetXaxis()->GetXmax()+shifts[energies.at(i)]);

      com_hists.back()->GetXaxis()->SetLimits(com_hists.back()->GetXaxis()->GetXmin()+shifts[energies.at(i)],
				              com_hists.back()->GetXaxis()->GetXmax()+shifts[energies.at(i)]);
      fit_hists.back()->Rebin(REBIN_FACTOR);
      fep_hists.back()->Rebin(REBIN_FACTOR);
      com_hists.back()->Rebin(REBIN_FACTOR); 
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
    /*
    else {
      fit_funcs.back()->SetParLimits(0,llims[getEnergy(fit_hists.at(i)->GetName())],
				       hlims[getEnergy(fit_hists.at(i)->GetName())]
				    );
    }
    */
    
    fit_funcs.back()->SetName(Form("func%04d",getEnergy(fit_hists.at(i)->GetName()))); 
  }

  TNtuple* bg_tpl = new TNtuple("bg_input","bg_input","energy:parameter:include");
  bg_tpl->ReadFile(bg_line_input.c_str());
  bg_tpl->Draw("energy:parameter:include","","goff");

  std::vector<int> bg_energies;
  std::map<int,double> bg_params;

  std::cout << "\nStopped Line Info: " << "\nEnergy\tScale" << std::endl;
  
  for(int i=0;i<bg_tpl->GetSelectedRows()-1;i++) {
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

    TFile* f = new TFile(Form("bg_hists/hist%04d.root",bg_energies.at(i)),"read");
    if(!f->IsZombie()) {
      bg_hists.push_back((GH1D*)f->Get(Form("%s_B&T&Y&D",MODE.c_str())));
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
      
  if(bg_hists.size()) {
    bg_hists.at(0)->Rebin(REBIN_FACTOR);
    fit_funcs.push_back(bg_hists.at(0)->ConstructTF1());
    //fit_funcs.back()->SetParameter(0,0.01);
    fit_funcs.back()->SetParameter(0,bg_tpl->GetV2()[bg_tpl->GetSelectedRows()-1]);
    if((bool)bg_tpl->GetV1()[bg_tpl->GetSelectedRows()-1]) {
      fit_funcs.back()->FixParameter(0,bg_tpl->GetV2()[bg_tpl->GetSelectedRows()-1]);
    }
    //else {
    //fit_funcs.back()->SetParLimits(0,0,1);
    //}
    fit_funcs.back()->SetName("StoppedLineFit");
  }

  fit_funcs.push_back(constructBackground(exp_bg_params));
  int npars = fit_funcs.back()->GetNpar();
  
  TF1Sum fSum(fitAllPeaks(data_hist,fit_funcs,fit_low_x,fit_high_x));

  fSum.GetFunc()->SetNpx(50000/REBIN_FACTOR);
  std::vector<double> fep_counts;
  std::vector<double> fep_counts_unc;
  std::vector<int> ens;
  for (unsigned int i=0;i<fep_hists.size();i++){
    fep_counts.push_back(fep_hists.at(i)->Integral()*fSum.GetFunc()->GetParameter(i));
    fep_counts_unc.push_back(fep_counts.at(i)*(fSum.GetFunc()->GetParError(i)/fSum.GetFunc()->GetParameter(i)));
    ens.push_back(getEnergy(fit_hists.at(i)->GetName()));
  }
  printResults(fep_counts,fep_counts_unc,ens);

  TH1* hf = fSum.GetFunc()->GetHistogram();
  
  int offset = 0;
  if(bg_hists.size()) {
    bg_hists.at(0)->Scale(fSum.GetFunc()->GetParameter(energies.size()));
    bg_hists.at(0)->SetName("Stopped_Lines");
    bg_hists.at(0)->SetTitle("Stopped_Lines");
    offset=1;
  }

  TF1* bg;
  if(npars == 2) {
    bg = new TF1("exp_bg","([0]*TMath::Exp([1]*x))*(1+TMath::TanH((x-[2])/[3]))/2",0,10000);
    bg->SetParameter(0,fSum.GetFunc()->GetParameter(energies.size()+offset));
    bg->SetParameter(1,fSum.GetFunc()->GetParameter(energies.size()+offset+1));
    bg->SetParameter(2,fSum.GetFunc()->GetParameter(energies.size()+offset+2));
    bg->SetParameter(3,fSum.GetFunc()->GetParameter(energies.size()+offset+3));
  }
  else {
    bg = new TF1("exp_bg","([0]*TMath::Exp([1]*x)+[2]*TMath::Exp([3]*x))*(1+TMath::TanH((x-[4])/[5]))/2",0,10000);
    bg->SetParameter(0,fSum.GetFunc()->GetParameter(energies.size()+offset));
    bg->SetParameter(1,fSum.GetFunc()->GetParameter(energies.size()+offset+1));
    bg->SetParameter(2,fSum.GetFunc()->GetParameter(energies.size()+offset+2));
    bg->SetParameter(3,fSum.GetFunc()->GetParameter(energies.size()+offset+3));
    bg->SetParameter(4,fSum.GetFunc()->GetParameter(energies.size()+offset+4));
    bg->SetParameter(5,fSum.GetFunc()->GetParameter(energies.size()+offset+5));
  }
  
  bg->SetNpx(50000/REBIN_FACTOR);
  TH1* hbg = bg->GetHistogram();
  hbg->SetName("Exp_Bg");
  hbg->SetTitle("Exp_Bg");
  
  //TH1* hr = (TH1*)data_hist.Clone();
  //hr->Sumw2();
  //hr->Add(hf,-1);

  //only necessary if data_hist and hf don't have same number of bins. If they do, use the above three lines
  TH1* hr = new TH1D("Residuum","Residuals;Energy (keV);Counts",
		     data_hist->GetNbinsX(),data_hist->GetXaxis()->GetXmin(),data_hist->GetXaxis()->GetXmax()
		    );

  TH1* he = new TH1D("Error","Error;Energy (keV);Counts",
		     data_hist->GetNbinsX(),data_hist->GetXaxis()->GetXmin(),data_hist->GetXaxis()->GetXmax()
		    );

  TH1* he1 = new TH1D("Error1","Error1;Energy (keV);Counts",
		     data_hist->GetNbinsX(),data_hist->GetXaxis()->GetXmin(),data_hist->GetXaxis()->GetXmax()
		    );

  TH1* he2 = new TH1D("Error2","Error2;Energy (keV);Counts",
		     data_hist->GetNbinsX(),data_hist->GetXaxis()->GetXmin(),data_hist->GetXaxis()->GetXmax()
		    );

  TH1* he21 = new TH1D("Error21","Error21;Energy (keV);Counts",
		     data_hist->GetNbinsX(),data_hist->GetXaxis()->GetXmin(),data_hist->GetXaxis()->GetXmax()
		    );

  TH1* h0 = new TH1D("Zero","Zero",
		     data_hist->GetNbinsX(),data_hist->GetXaxis()->GetXmin(),data_hist->GetXaxis()->GetXmax()
		    );
  
  for(int i=1;i<data_hist->GetNbinsX()+1;i++) {
    
    hr->SetBinContent(i,data_hist->GetBinContent(i) - hf->GetBinContent(i));
    hr->SetBinError(i,data_hist->GetBinError(i));

    he->SetBinContent(i,data_hist->GetBinError(i));
    he->SetBinError(i,0);

    he1->SetBinContent(i,-data_hist->GetBinError(i));
    he1->SetBinError(i,0);

    he2->SetBinContent(i,2*data_hist->GetBinError(i));
    he2->SetBinError(i,0);

    he21->SetBinContent(i,-2*data_hist->GetBinError(i));
    he21->SetBinError(i,0);

    h0->SetBinContent(i,0);
    h0->SetBinError(i,0);
    
  }

  data_hist->GetXaxis()->SetTitle("Energy (keV)");
  data_hist->GetYaxis()->SetTitle("Counts");
  
  TFile *outfile = new TFile(output_fn.c_str(), "recreate");
  
  data_hist->Write();
  hf->Write();
  
  if(bg_hists.size()) {
    bg_hists.at(0)->Write();
  }
  
  hbg->Write();
  
  hr->Write();
  he->Write();
  he1->Write();
  he2->Write();
  he21->Write();
  h0->Write();

  for(unsigned int i=0;i<fit_hists.size();i++) {
    
    fit_hists.at(i)->Scale(fSum.GetFunc()->GetParameter(i));
    fit_hists.at(i)->SetName(Form("Peak%02i",i));
    fit_hists.at(i)->SetTitle(Form("hist%04i",energies.at(i)));
    fit_hists.at(i)->Write();

    fep_hists.at(i)->Scale(fSum.GetFunc()->GetParameter(i));
    fep_hists.at(i)->SetName(Form("FEP%02i",i));
    fep_hists.at(i)->SetTitle(Form("FEP_hist%04i",energies.at(i)));
    fep_hists.at(i)->Write();

    com_hists.at(i)->Scale(fSum.GetFunc()->GetParameter(i));
    com_hists.at(i)->SetName(Form("COM%02i",i));
    com_hists.at(i)->SetTitle(Form("COM_hist%04i",energies.at(i)));
    com_hists.at(i)->Write();
    
  }
  
  outfile->Close();
  return;
}


int main(int argc, char **argv) { 
  
  std::string USAGE("fitGretinaPeaks input_data_file output_file_name peak_input bg_line_input bg_model_params x_low x_high\n");

  if (argc == 8) {
    
    std::string input_data_file(argv[1]);
    std::string output_fn(argv[2]);
    std::string peak_input(argv[3]);
    std::string bg_line_input(argv[4]);
    std::string bg_model_params(argv[5]);
    int fit_low_x = std::stoi(argv[6]);
    int fit_high_x = std::stoi(argv[7]);
    
    fitGretinaPeaks(input_data_file,output_fn,peak_input,bg_line_input,bg_model_params,fit_low_x,fit_high_x);
    
    return 0;
  }

  else {
    std::cout << "USAGE: " << USAGE;
    return -1;
  }
}
