#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <tuple>

#include "GH1D.h"
#include "TF1.h"
#include "TF1Sum.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TH1.h"
#include "TMatrixD.h"

//offset from start of string for peaks to get energy.
const int CHAR_OFFSET = 28;

double getEfficiency(double energy){
  return 4.532*pow(energy+100., -0.621)*8.75/8.;
}
void printResults(std::vector<double> fep_counts, std::vector<double> fep_counts_unc, std::vector<int> energies,
                  std::vector<bool> isStoppedLine){
  std::vector<double> intensities;
  std::vector<double> intensities_unc;
  std::cout << "============================== COUNTS ==============================\n";
  std::array<int,7> source_energies({511, 596, 718, 834, 844, 1014, 1039});
  int prev_stopped_lines = 0;
  for (unsigned int i = 0; i < energies.size(); i++){
    if (isStoppedLine.at(i)){
      prev_stopped_lines++;
      continue;
    }
    double eff = getEfficiency(energies.at(i));
    intensities.push_back(fep_counts.at(i-prev_stopped_lines)/eff);
    intensities_unc.push_back(intensities.back()*sqrt(pow(fep_counts_unc.at(i-prev_stopped_lines)/fep_counts.at(i-prev_stopped_lines), 2.) + 0.021*0.021));
    std::cout << energies.at(i) << "\t" << fep_counts.at(i-prev_stopped_lines) << "\t" << fep_counts_unc.at(i-prev_stopped_lines) << "\n";
  }

  std::cout << "\n\n\n";
  std::cout << "============================== INTENSITIES ==============================\n";
  prev_stopped_lines = 0;
  for (unsigned int i = 0; i < energies.size(); i++){
    if (isStoppedLine.at(i)){
      prev_stopped_lines++;
      continue;
    }
    std::cout << energies.at(i) << "\t" << intensities.at(i-prev_stopped_lines) << "\t" << intensities_unc.at(i-prev_stopped_lines) << "\n";
  }
}

bool isStopped(std::string line){
  const std::string stopped("stopped");

  if (line.find(stopped) != std::string::npos){
    return true;
  }
  return false;
}

std::vector<bool> findStoppedLines(std::vector<std::string> hist_file_names, std::vector<int> energies){
  std::vector<bool> isStoppedLine;
  for (unsigned int i = 0; i < hist_file_names.size(); i++){
    if(isStopped(hist_file_names.at(i))){
      std::cout << "Found stopped line with energy: " << energies.at(i) << "\n";
      isStoppedLine.push_back(true);
    }
    else{
      isStoppedLine.push_back(false);
    }
  }
  return isStoppedLine;
}
std::vector<std::string> parseInputFile(const std::string &input_fn, int fit_low_x, int fit_high_x){
  std::vector<std::string> hist_file_names;
  std::ifstream in_file(input_fn.c_str()); 
  std::string line;
  int energy;

  while (std::getline(in_file,line)){
    //because of this line, the peak energy MUST be 0 padded to length 4
    energy = std::stoi(line.substr(CHAR_OFFSET, 4));
    if (energy > fit_low_x && energy < fit_high_x){ 
      hist_file_names.push_back(line);
    }
  }

  return hist_file_names;

}

std::vector<int> getEnergies(const std::vector<std::string> &hist_file_names){
  std::vector<int> energies;
  for(auto &fn : hist_file_names){
    energies.push_back(std::stoi(fn.substr(CHAR_OFFSET, 4)));
  }

  return energies;
}

/*
 *This function will accept a filename containing a list of files in the form:
 *hist0275.root, where 0275 is the 0-padded peak energy of a Gretina GEANT 
 *simulated peak at 275 keV. 
 *
 */
TF1 *constructBackground(std::string background_type){


  TF1 *bg = 0;
  if (background_type == "se"){
    bg = new TF1("bg","[0]*TMath::Exp([1]*x)",0,10000);
    bg->SetParameter(0, 808.352);
    bg->SetParameter(1, -0.00306607);
  }

  else if (background_type == "de"){
    bg = new TF1("bg","[0]*TMath::Exp([1]*x)+[2]*TMath::Exp([3]*x)",0,10000);
    bg->SetParameter(0, 221.342);
    bg->SetParameter(1, -0.000697268);
    bg->SetParameter(2, 808.352);
    bg->SetParameter(3, -0.00306607);
  }

  else{
    std::cout << "Unknown background type. Reverting to single exponential\n";
    bg = new TF1("bg","[0]*TMath::Exp([1]*x)",0,10000);
    bg->SetParameter(0, 808.352);
    bg->SetParameter(1, -0.00306607);
  }
  bg->SetRange(0,10000);
  bg->SetNpx(10000);
  return bg;
}

TF1Sum fitAllPeaks(GH1D &data_hist, const std::vector<TF1*> &fit_funcs, int fit_low_x, int fit_high_x,
                   std::string background_type){
  TF1Sum fullSum;
  for (unsigned int i = 0; i < fit_funcs.size(); i++){
    fullSum.AddTF1(fit_funcs.at(i)); 
  }

  TF1 *de = constructBackground(background_type);
  de->SetNpx(10000/data_hist.GetXaxis()->GetBinWidth(1));
  fullSum.AddTF1(de);
  fullSum.GetFunc()->SetRange(0,10000);
  fullSum.GetFunc()->SetNpx(10000/data_hist.GetXaxis()->GetBinWidth(1));
  bool successful_fit = false;
  unsigned int count = 0;
  while (!successful_fit){
    TFitResultPtr r(data_hist.Fit(fullSum.GetFunc(), "MELS", "", fit_low_x, fit_high_x));
    if (r->IsValid() && count >= 1){
      successful_fit = true;
      std::cout << "Fit succeeded with r->Status() = " << r->Status()  << " r->IsValid() = " <<  r->IsValid() << std::endl;
      TMatrixD cov = r->GetCovarianceMatrix();
      TMatrixD cor = r->GetCorrelationMatrix();
      cov.Print();
      cor.Print();
    }
    else{ 
      std::cout << "Fit failed with r->Status() = " << r->Status()  << " r->IsValid() = " <<  r->IsValid() << std::endl;
    }
    count++;
  }

  for (int i = 0; i < fullSum.GetFunc()->GetNpar(); i++){
    std::cout << fullSum.GetFunc()->GetParName(i) << "\t" << fullSum.GetFunc()->GetParameter(i) << "\t"<<"+/-"<<"\t"<< fullSum.GetFunc()->GetParError(i) << "\n";
  }

  std::vector<double> residuals;
  GH1D *f = (GH1D*)fullSum.GetFunc()->GetHistogram();
  GH1D *sub = (GH1D*)data_hist.Clone("sub");
  sub->Add(f,-1);
  residuals.push_back(sub->Integral(sub->FindBin(300), sub->FindBin(800)));
  residuals.push_back(sub->Integral(sub->FindBin(800), sub->FindBin(1300)));
  residuals.push_back(sub->Integral(sub->FindBin(1300), sub->FindBin(2000)));
  residuals.push_back(sub->Integral(sub->FindBin(2000), sub->FindBin(3000)));
  std::cout << "\n\n==========RESIDUALS==========\n\n";
  std::cout << "(300,800): " << residuals.at(0) << "\n";
  std::cout << "(800,1300): " << residuals.at(1) << "\n";
  std::cout << "(1300,2000): " << residuals.at(2) << "\n";
  std::cout << "(2000,3000): " << residuals.at(3) << "\n\n";
  return fullSum;
}

std::map<int, std::pair<double,double>> getParMap(std::string element_name){
  using  Pair = std::pair<double,double>;
  using  ParMap = std::pair<int, Pair>;
  std::map<int, Pair> parmap;
  if (element_name == "cu"){
    //cu71-1p
    //stopped lines
    //before contam removal
//  parmap.insert(ParMap(0,Pair(1.0,0)));//overall scaling of summed stopped lines
//  parmap.insert(ParMap(511,Pair(0.04763,0)));
//  parmap.insert(ParMap(596,Pair(0.0101114,0)));
//  parmap.insert(ParMap(718,Pair(0.0069071,0)));
//  parmap.insert(ParMap(834,Pair(0.00670763,0)));
//  parmap.insert(ParMap(844,Pair(0.0045335,0)));
//  parmap.insert(ParMap(1014,Pair(0.00334314,0)));
//  parmap.insert(ParMap(1039,Pair(0.00861912,0)));

    
    //after contam removal
    parmap.insert(ParMap(0,Pair(1.0,0)));//overall scaling of summed stopped lines
    parmap.insert(ParMap(511,Pair(0.0435298,0)));
    parmap.insert(ParMap(596,Pair(0.0123906,0)));
    parmap.insert(ParMap(718,Pair(0.00677142,0)));
    parmap.insert(ParMap(834,Pair(0.0113627,0)));
    parmap.insert(ParMap(844,Pair(0.00971463,0)));
    parmap.insert(ParMap(1014,Pair(0.00713203,0)));
    parmap.insert(ParMap(1039,Pair(0.00805667,0)));

    parmap.insert(ParMap(448,Pair(0.00661639,0)));
    parmap.insert(ParMap(631,Pair(0.0136572,1.77925)));
    parmap.insert(ParMap(640,Pair(0.0323438,0)));
    parmap.insert(ParMap(664,Pair(0.00907501,0)));
    parmap.insert(ParMap(676,Pair(0.001,0.0)));
    parmap.insert(ParMap(683,Pair(0.01,0)));
    parmap.insert(ParMap(915,Pair(0.0179,0)));
    parmap.insert(ParMap(930,Pair(0.00800364,1.79118)));
    parmap.insert(ParMap(958,Pair(0.009703,0)));
    parmap.insert(ParMap(970,Pair(0.221035,0)));
    parmap.insert(ParMap(1075,Pair(0.0204767,3.11125)));
    parmap.insert(ParMap(1198,Pair(0.0145478,-0.5)));
    parmap.insert(ParMap(1247,Pair(0.0341441,-2.1644)));
    parmap.insert(ParMap(1260,Pair(0.587616,0.520159)));
    parmap.insert(ParMap(1321,Pair(0.00953202,0)));
    parmap.insert(ParMap(1428,Pair(0.00340795,1.45695)));
    parmap.insert(ParMap(1440,Pair(0.00845154,0.874094)));
    parmap.insert(ParMap(1467,Pair(0.00498927,4)));
    parmap.insert(ParMap(1583,Pair(0.0529319,0)));
    parmap.insert(ParMap(1666,Pair(0.0028249,0)));
    parmap.insert(ParMap(1682,Pair(0.013533,0)));
    parmap.insert(ParMap(1787,Pair(0.0180573,0)));
    parmap.insert(ParMap(1868,Pair(0.0395677,0)));
    parmap.insert(ParMap(1955,Pair(0.0253539,4)));
    parmap.insert(ParMap(2026,Pair(0.0160012,0)));
    parmap.insert(ParMap(2071,Pair(0.00948323,0.676568)));
    parmap.insert(ParMap(2114,Pair(0.0101164,0.764034)));
    parmap.insert(ParMap(2342,Pair(0.00747437,0)));
    parmap.insert(ParMap(2399,Pair(0.0116595,0.560432)));
    parmap.insert(ParMap(2697,Pair(0.0226159,2.43902)));
    parmap.insert(ParMap(2758,Pair(0.0196746,0)));
    parmap.insert(ParMap(2980,Pair(0.0115789,0)));
    parmap.insert(ParMap(3036,Pair(0.0185723,4)));
    parmap.insert(ParMap(3846,Pair(0.00712318,4)));
    //boosted lines
//  parmap.insert(ParMap(275,Pair(0.00661966,0.75801)));
//  parmap.insert(ParMap(385,Pair(0.0191657,0)));
//  parmap.insert(ParMap(479,Pair(0.00748137,0)));
//  parmap.insert(ParMap(609,Pair(0.0508895,0)));
//  parmap.insert(ParMap(631,Pair(0.0130008,0.939047)));
//  parmap.insert(ParMap(640,Pair(0.0334761,0)));
//  parmap.insert(ParMap(664,Pair(0.00991149,0)));
//  parmap.insert(ParMap(676,Pair(0.0073168,0.808131)));
//  parmap.insert(ParMap(683,Pair(0.00495003,0)));
//  parmap.insert(ParMap(915,Pair(0.011477,0)));
//  parmap.insert(ParMap(958,Pair(0.00754714,0)));
//  parmap.insert(ParMap(970,Pair(0.209511,0)));
//  parmap.insert(ParMap(1075,Pair(0.0125948,2.5)));
//  parmap.insert(ParMap(1247,Pair(0.0340944,1.00)));
//  parmap.insert(ParMap(1260,Pair(0.550913,1.648782)));
//  parmap.insert(ParMap(1321,Pair(0.0108471,0)));
//  parmap.insert(ParMap(1428,Pair(0.00439928,1.99965)));
//  parmap.insert(ParMap(1440,Pair(0.0080098,1.59399)));
//  parmap.insert(ParMap(1467,Pair(0.00606613,2.5)));
//  parmap.insert(ParMap(1583,Pair(0.0548508,0)));
//  parmap.insert(ParMap(1666,Pair(0.00451966,-0.5)));
//  parmap.insert(ParMap(1682,Pair(0.0146917,0)));
//  parmap.insert(ParMap(1787,Pair(0.0195832,0)));
//  parmap.insert(ParMap(1868,Pair(0.0407686,0)));
//  parmap.insert(ParMap(1955,Pair(0.0279031,2.5)));
//  parmap.insert(ParMap(2026,Pair(0.0167558,1.09732e-05)));
//  parmap.insert(ParMap(2071,Pair(0.00986898,1.99994)));
//  parmap.insert(ParMap(2114,Pair(0.0105043,0.793158)));
//  parmap.insert(ParMap(2342,Pair(0.00752715,-0.408376)));
//  parmap.insert(ParMap(2399,Pair(0.0118713,0.589449)));
//  parmap.insert(ParMap(2697,Pair(0.0225466,1.26662)));
//  parmap.insert(ParMap(2758,Pair(0.0196882,-0.5)));
//  parmap.insert(ParMap(2980,Pair(0.0115273,-0.5)));
//  parmap.insert(ParMap(3036,Pair(0.0187725,2.5)));
//  parmap.insert(ParMap(3846,Pair(0.0079845,2.00023))); 
  }

  else if (element_name == "ni"){
    //ni71-1n
    //stopped lines
  //parmap.insert(ParMap(0,Pair(1.0,0)));
  //parmap.insert(ParMap(511,Pair(0.00389334,0)));
  //parmap.insert(ParMap(597,Pair(0.00204799,0)));
  //parmap.insert(ParMap(836,Pair(0.00226323,0)));
  //parmap.insert(ParMap(1039,Pair(0.00158093,0)));
    parmap.insert(ParMap(0,Pair(1.0,0)));
    parmap.insert(ParMap(511,Pair(0.0037307,0)));
    parmap.insert(ParMap(597,Pair(0.00244869,0)));
    parmap.insert(ParMap(836,Pair(0.00243733,0)));
    parmap.insert(ParMap(1039,Pair(0.00182415,0)));
    //boosted lines
    parmap.insert(ParMap(230,Pair(0.00350921,0)));
    parmap.insert(ParMap(385,Pair(0.0079055,0)));
    parmap.insert(ParMap(424,Pair(0.00168021,0.0)));
    parmap.insert(ParMap(444,Pair(0.00136974,0)));
    parmap.insert(ParMap(609,Pair(0.00547365,0)));
    parmap.insert(ParMap(626,Pair(0.00144228,0)));
    parmap.insert(ParMap(640,Pair(0.00335833,0)));
    parmap.insert(ParMap(648,Pair(0.00218808,0.0)));
    parmap.insert(ParMap(660,Pair(0.00384201,0.0)));
    parmap.insert(ParMap(676,Pair(0.00509541,0)));
    parmap.insert(ParMap(683,Pair(0.00209193,0)));
    parmap.insert(ParMap(872,Pair(0.00262068,0.0)));
    parmap.insert(ParMap(912,Pair(0.00450256,2.0)));
    parmap.insert(ParMap(930,Pair(0.00297641,2.0)));
    parmap.insert(ParMap(942,Pair(0.00297641,2.0)));
    parmap.insert(ParMap(946,Pair(0.00297641,2.0)));
    parmap.insert(ParMap(953,Pair(0.0719071,2.0)));
    parmap.insert(ParMap(958,Pair(0.0719071,0.0)));
    parmap.insert(ParMap(970,Pair(0.0131901,2.0)));
    parmap.insert(ParMap(1075,Pair(0.000864805,0)));
    parmap.insert(ParMap(1189,Pair(0.004083,0)));
    parmap.insert(ParMap(1225,Pair(0.00434377,0)));
    parmap.insert(ParMap(1243,Pair(0.00493533,0)));
    parmap.insert(ParMap(1260,Pair(0.0403839,0)));
    parmap.insert(ParMap(1678,Pair(0.0002,0)));
    parmap.insert(ParMap(1868,Pair(0.0037445,0)));
    parmap.insert(ParMap(1916,Pair(0.0002,0)));
    parmap.insert(ParMap(1957,Pair(0.00223584,0)));
    parmap.insert(ParMap(2105,Pair(0.00127498,0)));
  }

  else if (element_name == "zn"){
    //zn72-2p
    //stopped lines
    //before contam removal
    //parmap.insert(ParMap(0,Pair(1.0,0)));
    //parmap.insert(ParMap(511,Pair(0.00265228,0)));
    //parmap.insert(ParMap(596,Pair(0.000667017,0)));
    //parmap.insert(ParMap(718,Pair(0.000407026,0)));
    //parmap.insert(ParMap(834,Pair(0.00050978,0)));
    //parmap.insert(ParMap(1039,Pair(0.000592849,0)));
    //after contam removal
      parmap.insert(ParMap(0,Pair(1.0,0)));
      parmap.insert(ParMap(511,Pair(0.00297622,0)));
      parmap.insert(ParMap(596,Pair(0.000558640,0)));
      parmap.insert(ParMap(834,Pair(0.000174398,0)));

      parmap.insert(ParMap(385,Pair(0.0010504,-0.499998)));
      parmap.insert(ParMap(406,Pair(0.000462085,0.0112393)));
      parmap.insert(ParMap(420,Pair(0.000830714,-0.00966671)));
      parmap.insert(ParMap(609,Pair(0.00282452,-0.499999)));
      parmap.insert(ParMap(626,Pair(0.00054261,3.99999)));
      parmap.insert(ParMap(640,Pair(0.00256204,1.44612)));
      parmap.insert(ParMap(660,Pair(0.00132496,2.84014)));
      parmap.insert(ParMap(676,Pair(0.00164968,2.73191)));
      parmap.insert(ParMap(710,Pair(0.00135391,2.74687)));
      parmap.insert(ParMap(932,Pair(0.000906392,-0.499984)));
      parmap.insert(ParMap(953,Pair(0.000975905,-0.386969)));
      parmap.insert(ParMap(970,Pair(0.0089368,-0.112842)));
      parmap.insert(ParMap(1075,Pair(0.0179312,-0.499848)));
      parmap.insert(ParMap(1212,Pair(0.00149868,1.35044)));
      parmap.insert(ParMap(1243,Pair(0.00175887,-0.499993)));
      parmap.insert(ParMap(1260,Pair(0.0237596,0.0889338)));
      parmap.insert(ParMap(1442,Pair(0.000551466,0.363908)));
      parmap.insert(ParMap(1583,Pair(5.69544e-05,0.0993606)));
      parmap.insert(ParMap(1666,Pair(8.75671e-06,0.190195)));
      parmap.insert(ParMap(1682,Pair(0.000586391,3.99559)));
      parmap.insert(ParMap(1868,Pair(0.00212205,0.000270743)));
      parmap.insert(ParMap(1955,Pair(0.000532261,-0.499875)));
      parmap.insert(ParMap(2035,Pair(0.00209021,0.0250933)));
  //parmap.insert(ParMap(0,Pair(1.0,0)));
  //parmap.insert(ParMap(511,Pair(0.00307745,0)));
  //parmap.insert(ParMap(596,Pair(0.000625677,0)));
  //parmap.insert(ParMap(1039,Pair(0.000623252,0)));
    //boosted lines
//  parmap.insert(ParMap(385,Pair(0.000994279,0)));
//  parmap.insert(ParMap(609,Pair(0.00547365,0)));
//  parmap.insert(ParMap(626,Pair(0.00144228,0)));
//  parmap.insert(ParMap(640,Pair(0.00335833,0)));
//  parmap.insert(ParMap(648,Pair(0.00218808,0.0)));
//  parmap.insert(ParMap(660,Pair(0.00384201,0.0)));
//  parmap.insert(ParMap(676,Pair(0.00509541,0)));
//  parmap.insert(ParMap(683,Pair(0.00209193,0)));
//  parmap.insert(ParMap(912,Pair(0.00450256,2.0)));
//  parmap.insert(ParMap(930,Pair(0.00297641,2.0)));
//  parmap.insert(ParMap(942,Pair(0.00297641,2.0)));
//  parmap.insert(ParMap(946,Pair(0.00297641,2.0)));
//  parmap.insert(ParMap(953,Pair(0.0719071,2.0)));
//  parmap.insert(ParMap(958,Pair(0.0719071,2.0)));
//  parmap.insert(ParMap(970,Pair(0.0131901,2.0)));
//  parmap.insert(ParMap(1075,Pair(0.000864805,0)));
//  parmap.insert(ParMap(1225,Pair(0.00434377,0)));
//  parmap.insert(ParMap(1243,Pair(0.00493533,0)));
//  parmap.insert(ParMap(1260,Pair(0.0403839,0)));
//  parmap.insert(ParMap(1678,Pair(0.0002,0)));
//  parmap.insert(ParMap(1868,Pair(0.0037445,0)));
//  parmap.insert(ParMap(1916,Pair(0.0002,0)));
//  parmap.insert(ParMap(1957,Pair(0.00223584,0)));
//  parmap.insert(ParMap(2105,Pair(0.00127498,0)));
  }
  else{
    std::cout << "Unknown element name!\n";
  }
  return parmap;
}

void fitGretinaPeaks(std::string input_file_list, std::string input_wider_list, std::string data_file_name, 
    std::string data_hist_dir, int fit_low_x, int fit_high_x, std::string background_type,
    std::string element_name){

  const int REBIN_FACTOR = 2;
  std::vector<std::string> hist_file_names = parseInputFile(input_file_list, fit_low_x, fit_high_x);
  std::cout << "Found " << hist_file_names.size() << " input files!" << std::endl;

  TFile *data_file = new TFile(data_file_name.c_str(), "read");
  if (data_file == nullptr){
    std::cout << "Failed to open data file\n";
    return;
  }
  GH1D data_hist(*((TH1*)data_file->Get(Form("%s/gretina_allcorr_prompt", data_hist_dir.c_str()))));
  data_hist.Sumw2();
  data_hist.Rebin(REBIN_FACTOR);
  std::vector<int> energies = getEnergies(hist_file_names);
  std::vector<bool> isStoppedLine = findStoppedLines(hist_file_names, energies);

  std::vector<TFile*> hist_files;
  for (unsigned int i = 0; i <  hist_file_names.size(); i++){
    hist_files.push_back(new TFile(hist_file_names.at(i).c_str(), "read"));
    if (hist_files.back() == nullptr){
      std::cout << "Failed to open hist_files_wider file\n";
      return;
    }
  }
  std::vector<GH1D> fit_hists;
  std::vector<GH1D> fep_hists;
  const double NARROW_SCALE = 0.8;
  for (unsigned int i = 0; i < hist_files.size(); i++){
    fit_hists.push_back(*((TH1*)hist_files.at(i)->Get("gretsim/gretina_allcorr")));
    fep_hists.push_back(*((TH1*)hist_files.at(i)->Get("gretsim/gretina_allcorr_fep")));
    fit_hists.back().Sumw2();
    fit_hists.back().Scale(NARROW_SCALE);//
    fep_hists.back().Scale(NARROW_SCALE);//
  }

  std::vector<std::string> hist_file_names_wider = parseInputFile(input_wider_list, fit_low_x, fit_high_x);
  std::vector<TFile*> hist_files_wider;
  for (unsigned int i = 0; i <  hist_file_names_wider.size(); i++){
    hist_files_wider.push_back(new TFile(hist_file_names_wider.at(i).c_str(), "read"));
    if (hist_files_wider.back() == nullptr){
      std::cout << "Failed to open hist_files_wider file\n";
      return;
    }
  }
  std::vector<GH1D> fit_hists_wider;
  std::vector<GH1D> fep_hists_wider;
  const double WIDE_SCALE = 0.2;
  for (unsigned int i = 0; i < hist_files_wider.size(); i++){
    fit_hists_wider.push_back(*((TH1*)hist_files_wider.at(i)->Get("gretsim/gretina_allcorr")));
    fep_hists_wider.push_back(*((TH1*)hist_files_wider.at(i)->Get("gretsim/gretina_allcorr_fep")));
    fit_hists_wider.back().Sumw2();
    fit_hists_wider.back().Scale(WIDE_SCALE);
    fep_hists_wider.back().Scale(WIDE_SCALE);
  }

  if (fit_hists.size() != fit_hists_wider.size()){
    std::cout << "Fit func vector sizes are different!";
    return;
  }

  //Combine good and bad decomp. components
  for (unsigned int i = 0; i < fit_hists.size(); i++){
    fit_hists.at(i).Add(&fit_hists_wider.at(i)); 
    fep_hists.at(i).Add(&fep_hists_wider.at(i));
    fep_hists.at(i).SetDirectory(0);
  }


  std::vector<TF1*> fit_funcs;
  auto parmap = getParMap(element_name); 

  GH1D *stopped_hist_ptr = 0;
  for (unsigned int i = 0; i < fit_hists.size(); i++){
    fit_hists.at(i).Rebin(data_hist.GetXaxis()->GetBinWidth(1));

    if (isStoppedLine.at(i)){
      auto r = parmap.find(energies.at(i));
      if (r != parmap.end()){
        fit_hists.at(i).Scale(r->second.first);
      }
      if (!stopped_hist_ptr){
        stopped_hist_ptr = new GH1D(fit_hists.at(i)); 

      }
      else{
        stopped_hist_ptr->Add((GH1D*)(&fit_hists.at(i)));
      }
    }
    else{
      fit_funcs.push_back(fit_hists.at(i).ConstructTF1Shift()); 
      auto r = parmap.find(energies.at(i));
      if(r != parmap.end()){
        if(r->second.second < 1e-02){
          fit_funcs.back()->FixParameter(1,0);
        }

        if(energies.at(i) == 448 && element_name == "cu"){
          fit_funcs.back()->FixParameter(0, r->second.first);
        }
        else{
          fit_funcs.back()->SetParameter(1, r->second.second);
          fit_funcs.back()->SetParLimits(1, -4.0, 4.0);//shift
        }
        //FIXING FOR NI-71
        if (energies.at(i) > 6000){
          fit_funcs.back()->FixParameter(0, r->second.first);
        }
        else{
          fit_funcs.back()->SetParameter(0, r->second.first);
          fit_funcs.back()->SetParLimits(0, 0, 1);//scale
        }
      }
      else{
        fit_funcs.back()->SetParLimits(1, -0.5, 4.0);//shift
        fit_funcs.back()->SetParLimits(0, 0, 1);//scale
      }
      fit_funcs.back()->SetName(Form("hist%04d", energies.at(i)));
    }
  }

  GH1D stopped_hist;
  if (stopped_hist_ptr != nullptr){
    stopped_hist = *stopped_hist_ptr;
    fit_funcs.push_back(stopped_hist.ConstructTF1Shift());
    auto r = parmap.find(0);//I let 0 be thei overall scaling for the stopped lines!
    if (r != parmap.end()){
      if (r->second.second < 1e-02){
        fit_funcs.back()->FixParameter(1,0);
      }
      else{
        fit_funcs.back()->SetParameter(1,r->second.second);
        fit_funcs.back()->SetParLimits(1,-0.5,4.0);
      }
      //  fit_funcs.back()->SetParameter(0,r->second.first);
      fit_funcs.back()->FixParameter(0, 1);
      //  fit_funcs.back()->SetParLimits(0, 0, 10);

    }

    else{
      fit_funcs.back()->SetParLimits(1, -0.5, 4.0);//shift
    }
    //fit_funcs.back()->SetParLimits(0, 0, 1);//scale
    fit_funcs.back()->SetName(Form("stopped_lines"));
  }

  TF1Sum fSum(fitAllPeaks(data_hist, fit_funcs, fit_low_x, fit_high_x, background_type));
  fSum.GetFunc()->SetNpx(10000/data_hist.GetXaxis()->GetBinWidth(1));
  std::vector<double> fep_counts;
  std::vector<double> fep_counts_unc;
  int cur_par = 0;
  for (unsigned int i = 0; i < fep_hists.size(); i++){
    if (isStoppedLine.at(i)){
      continue;
    }
    fep_counts.push_back(fep_hists.at(i).Integral()*fSum.GetFunc()->GetParameter(cur_par));
    fep_counts_unc.push_back(fep_counts.back()*(fSum.GetFunc()->GetParError(cur_par)/fSum.GetFunc()->GetParameter(cur_par)));
    fep_hists.at(i).Scale(fSum.GetFunc()->GetParameter(cur_par));//to save all the final hists
    cur_par += 2;
  }
  printResults(fep_counts, fep_counts_unc, energies, isStoppedLine);
  TFile *outfile = new TFile(Form("%s_%04d_%04d_%s.root", element_name.c_str(), fit_low_x, fit_high_x,background_type.c_str()), "recreate");
  std::cout << "saving output fit to " << outfile->GetName() << "\n";

  //Clone fep hists for saving
  std::vector<GH1D *> scaled_fep_hists;
  std::cout << "saving scaled fep hists\n";
  for (unsigned int i = 0; i < fep_hists.size(); i++){
    scaled_fep_hists.push_back((GH1D*)fep_hists.at(i).Clone(Form("%s_%04d_scaled",fep_hists.at(i).GetName(), energies.at(i))));
    scaled_fep_hists.back()->SetDirectory(outfile);
    scaled_fep_hists.back()->Write();
  }
  if (stopped_hist_ptr != nullptr){
    stopped_hist.SetDirectory(outfile);
    stopped_hist.SetName("stopped_hist");
    stopped_hist.Write();
  }

  fSum.GetFunc()->GetHistogram()->Write();
  data_hist.Write();
  outfile->Close();
  return;
}


int main(int argc, char **argv){ 

  std::string USAGE("fitGretinaPeaks input_goodres_list input_badres_list input_data_file data_hist_dir leftHighRange rightLowRange background_type element_name\n");

  if (argc == 9){
    std::cout << "using good res input file: " << argv[1] << "\n";
    std::string input_thinner_list(argv[1]);
    std::cout << "using wide input file: " << argv[2] << "\n";
    std::string input_wider_list(argv[2]);
    std::cout << "using data file: " << argv[3] << "\n";
    std::string input_data_file(argv[3]);
    std::cout << "assuming spectrum in folder: " << argv[4] << "\n";
    std::string data_hist_dir(argv[4]);
    std::cout << "fitting region: (" << argv[5] << ","<<argv[6]<<")\n";
    int fit_low_x = std::stoi(argv[5]);
    int fit_high_x = std::stoi(argv[6]);
    std::cout << "Background Type: " << argv[7] <<"\n";
    std::cout << "Element name: " << argv[8] <<"\n";
    
    fitGretinaPeaks(input_thinner_list,input_wider_list, input_data_file, data_hist_dir, fit_low_x,fit_high_x, argv[7], argv[8]);
    return 0;
  }

  else {
    std::cout << USAGE;
    return -1;
  }
}
