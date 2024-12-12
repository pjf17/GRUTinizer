#include "GPolAnalyzer.h"
#include "TGretina.h"
#include "GH1D.h"
#include "GH2D.h"

#include "TString.h"
#include "TH1.h"
#include "TF1.h"
#include "TFile.h"
#include "TList.h"
#include "TKey.h"
#include "TDirectory.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"

#include <string>
#include <iostream>

GPolAnalyzer::GPolAnalyzer(TH1 *hsource, TH1 *hbeam, int EsrceLo, int EsrceHi, int EbeamLo, int EbeamHi) {
  fDoRand = false;
  fSrce = hsource->GetDirectory()->GetMotherDir();
  fBeam = hbeam->GetDirectory()->GetMotherDir();

  std::string dirname = std::string(hbeam->GetDirectory()->GetName());
  fMassNum = TString(dirname.substr(dirname.find("_gated")-4,4));
  if(EsrceHi < EsrceLo) std::swap(EsrceLo,EsrceHi);
  fEsrceGate[0] = EsrceLo;
  fEsrceGate[1] = EsrceHi;
  if(EbeamHi < EbeamLo) std::swap(EbeamLo,EbeamHi);
  fEbeamGate[0] = EbeamLo;
  fEbeamGate[1] = EbeamHi;
}
    
GPolAnalyzer::GPolAnalyzer(TH1 *h, int ELo, int EHi){
  fDoRand = true;
  fSrce = h->GetDirectory()->GetMotherDir();
  fBeam = h->GetDirectory()->GetMotherDir();

  std::string dirname = std::string(h->GetDirectory()->GetName());
  fMassNum = TString(dirname.substr(dirname.find("_gated")-4,4));
  if(EHi < ELo) std::swap(ELo,EHi);
  fEbeamGate[0] = ELo;
  fEbeamGate[1] = EHi;
  fEsrceGate[0] = ELo;
  fEsrceGate[1] = EHi;
}

void GPolAnalyzer::GetXiRatio(GH1D *&hnorm, GH1D *&hsrce, GH1D *&hbeam, int binning, TString opt) {
  bool linFit = opt.Contains("lin");
  if (linFit) opt.ReplaceAll("lin","");
  bool saveFile = opt.Contains("save");
  if (saveFile) opt.ReplaceAll("save","");

  //find the directories for the polarization histograms, must have "polarization" in the name
  TString BeamDir = getHistDirectories(fBeam);
  TString SrceDir;
  if (fDoRand) 
    SrceDir = BeamDir;
  else 
    SrceDir = getHistDirectories(fSrce);
  
  printf("\n==========LINEAR POLARIZATION RESULTS==========\n");
  // printf("Beam hist directory: %s\nNorm hist directory: %s\n",BeamDir.Data(),SrceDir.Data());
  //read the hist names from the directories, must end in the crystal number
  std::vector<std::pair<int,std::string>> srceCrystals;
  std::vector<std::pair<int,std::string>> beamCrystals;
  getHistNames(fSrce,SrceDir,srceCrystals,fDoRand);
  getHistNames(fBeam,BeamDir,beamCrystals,false);

  //split them into the chosen groups
  groupCrystals(srceCrystals,opt);
  groupCrystals(beamCrystals,opt);

  //loop over the histograms
  int N = srceCrystals.size();  
  int nGroups = 0;
  bool readNew = true;
  double scalingRatio = 0;
  double energyLo = 0;
  double energyHi = 0;

  GH1D *stemp;
  GH1D *dtemp;
  GH2D *h2d;

  // ***NOTE*** if srceCrystals[i].first > 1000 then we're just grouping all crystals together
  for (int i=0; i < N; i++){
    if (srceCrystals[i].first != beamCrystals[i].first) {std::cout<<"BEAM AND SOURCE CRYSTAL MISMATCH\n"; break;}
    
    //read in and overwrite the temp histograms
    if (readNew) {
      //source
      h2d = (GH2D*) fSrce->Get(Form("%s/%s",SrceDir.Data(),srceCrystals[i].second.c_str()));
      stemp = h2d->ProjectionX("_px",fEsrceGate[0],fEsrceGate[1]);
      stemp->Sumw2();

      //beam
      h2d = (GH2D*) fBeam->Get(Form("%s/%s",BeamDir.Data(),beamCrystals[i].second.c_str()));
      dtemp = h2d->ProjectionX("_px",fEbeamGate[0],fEbeamGate[1]);
      dtemp->Sumw2();

      energyHi = h2d->GetYaxis()->GetBinLowEdge(fEbeamGate[1]+1);
      energyLo = h2d->GetYaxis()->GetBinLowEdge(fEbeamGate[0]);
      
      scalingRatio += stemp->Integral()/dtemp->Integral();
      readNew = false;
    } 

    //if current group is the same as the last, add to the temp histograms
    if (i > 0 && (srceCrystals[i].first == srceCrystals[i-1].first || srceCrystals[i].first > 1000)) {
      //source
      h2d = (GH2D*) fSrce->Get(Form("%s/%s",SrceDir.Data(),srceCrystals[i].second.c_str()));
      stemp->Add(h2d->ProjectionX("_px",fEsrceGate[0],fEsrceGate[1]));

      //beam
      h2d = (GH2D*) fBeam->Get(Form("%s/%s",BeamDir.Data(),beamCrystals[i].second.c_str()));
      dtemp->Add(h2d->ProjectionX("_px",fEbeamGate[0],fEbeamGate[1])); 
    }
    
    //if a new group is coming, add the queue to the output
    if (i == N-1 || (srceCrystals[i].first != srceCrystals[i+1].first && srceCrystals[i].first < 1000) ) {
      //initialize output hists 
      if (nGroups == 0) {
        hsrce = (GH1D*) stemp->Clone("source"); 
        hbeam = (GH1D*) dtemp->Clone("data");
        dtemp->Divide(stemp);
        hnorm = (GH1D*) dtemp->Clone("norm");
      }
      else {
        hsrce->Add(stemp);
        hbeam->Add(dtemp);
        dtemp->Divide(stemp);
        hnorm->Add(dtemp);
      }

      readNew = true;
      nGroups++;
    }
  }

  if (nGroups > 1){
    scalingRatio /= nGroups;
    hnorm->Scale(1.0/nGroups*scalingRatio);
  }
  
  //fit the normalized distribution
  TF1 *fitfunc;
  int parA0;
  if (linFit) {
    fitfunc = new TF1("mypol","([0] + [1]*x)*(1-[2]*TMath::Cos(2*x))",0,TMath::TwoPi());
    parA0 = 2;
  }
  else {
    fitfunc = new TF1("mypol","[0]*(1-[1]*TMath::Cos(2*x))",0,TMath::TwoPi());
    parA0 = 1;
  }

  //fit with 1 degree binning (best for determining asym)
  TFitResultPtr r1;
  r1 = hnorm->Fit(fitfunc,"SQN0");

  if (binning != 1){
    hnorm->Rebin(binning);
    hsrce->Rebin(binning);
    hbeam->Rebin(binning);
    hnorm->Scale(1.0/binning);
  }
  
  //fit again to show the function on chosen binning
  TFitResultPtr rbw;
  rbw = hnorm->Fit(fitfunc,"SQ");

  //output results
  printf("Grouping mode: %s\n",opt.Data());
  printf("Gate => |%3.1f < %3.1f > %-3.1f| \n",energyLo,(energyLo+energyHi)/2,energyHi);
  printf("Total Events bewteen gates: %3.1f\n\n",hbeam->Integral());
  printf("Bin Width:      1 deg      |        %d deg   \n",binning);
  printf("Const  = %5.4f +/- %5.4f | %5.4f +/- %5.4f\n",r1->Parameter(0),r1->ParError(0),rbw->Parameter(0),rbw->ParError(0));
  if (linFit) printf("Slope  = %5.4f +/- %5.4f | %5.4f +/- %5.4f\n",r1->Parameter(1),r1->ParError(1),rbw->Parameter(1),rbw->ParError(1));
  printf("A0     = %5.4f +/- %5.4f | %5.4f +/- %5.4f\n",r1->Parameter(parA0),r1->ParError(parA0),rbw->Parameter(parA0),rbw->ParError(parA0));
  printf("A0/err =      %5.4f       |      %5.4f \n\n",r1->Parameter(parA0)/r1->ParError(parA0),rbw->Parameter(parA0)/rbw->ParError(parA0));

  if (saveFile) {
    TString outfilename = Form("lpl_%s_%3.0f_%3.0f_bw%d_%s.root",fMassNum.Data(),energyLo,energyHi,binning,opt.Data());
    TFile *foutput = new TFile(outfilename,"RECREATE");
    foutput->cd();
    hnorm->Write();
    hsrce->Write();
    hbeam->Write();
    foutput->Close();
    printf("Created %s\n",outfilename.Data());
  }

  //scale down the source hist to match data hist
  double scaling = 1.0/(std::floor(hsrce->Integral()/hbeam->Integral())); 
  hsrce->Scale(scaling);
  hbeam->SetLineColor(kRed);
  double ymax = std::max(hbeam->GetMaximum(),hsrce->GetMaximum());
  double ymin = std::min(hbeam->GetMaximum(),hsrce->GetMaximum());
  hsrce->GetYaxis()->SetRangeUser(ymin*0.9,ymax*1.1);
  hbeam->GetYaxis()->SetRangeUser(ymin*0.9,ymax*1.1);
}

void GPolAnalyzer::getHistNames(TDirectory *f, TString dir, std::vector<std::pair<int,std::string>> &names, bool randMode) {
  //get the histograms in the directory
  TDirectory* input = (TDirectory*) f->Get(dir.Data());
  TList *l = input->GetListOfKeys();
  int N = l->GetEntries();
  for (int i=0 ; i < N; i++) {
    std::string histname = std::string(l->At(i)->GetName());
    TString cn = TString(histname.substr(histname.size()-2,2));
    int cryNum = cn.Atoi();
    if (!(cryNum > 3 && cryNum < 124)) continue;
    bool isRand = histname.find("rand") != std::string::npos;
    if (randMode && isRand) names.push_back(std::make_pair(cryNum,histname));
    if (!randMode && !isRand) names.push_back(std::make_pair(cryNum,histname));
  }

  //sort by crystal number
  std::sort(names.begin(),names.end());
}

TString GPolAnalyzer::getHistDirectories(TDirectory *f){
  TList *l = f->GetListOfKeys();
  int N = l->GetEntries();
  std::vector<TString> poldirs;
  for (int n=0; n < N; n++){
    TString keyname = TString(l->At(n)->GetName());
    if (keyname.Contains("polarization")) poldirs.push_back(keyname);
  }

  int npol = poldirs.size();
  if (npol > 1) {
    for (int i=0; i < npol; i++){
      if (poldirs[i].Contains(fMassNum)){
        return poldirs[i];
      }
    }
  }
  else return poldirs[0];
}

void GPolAnalyzer::groupCrystals(std::vector<std::pair<int,std::string>> &crstl, TString mode){
  int N = crstl.size();
  if (mode.Contains("quad")) {
    for (int i=0; i < N; i++){
      crstl[i].first = crstl[i].first/4 - 1; 
    }
    std::sort(crstl.begin(),crstl.end());
  }
  if (mode.Contains("ring")) {
    TGretina *gret = new TGretina();
    for (int i=0; i < N; i++){
      crstl[i].first = gret->GetRingNumber(crstl[i].first);
    }
    std::sort(crstl.begin(),crstl.end());
  }
  if (mode.Contains("all")) {
    for (int i=0; i < N; i++){
      crstl[i].first += 1000;
    }
    std::sort(crstl.begin(),crstl.end());
  }
  return;
}
