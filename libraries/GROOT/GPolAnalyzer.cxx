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

double calcQ(double E){
  double Qp = 1.0/(E/511 + 511/(E+511));
  double p0 = 0.131;
  double p1 = 3.25E-4;
  return (p0 + p1*E)*Qp;
}

GPolAnalyzer::GPolAnalyzer(TH1 *hsource, TH1 *hbeam, int EsrceLo, int EsrceHi, int EbeamLo, int EbeamHi, int EbkgLo, int EbkgHi) {
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
  if(EbkgHi < EbkgLo) std::swap(EbkgLo,EbkgHi);
  fEbkgGate[0] = EbkgLo;
  fEbkgGate[1] = EbkgHi;
}
    
GPolAnalyzer::GPolAnalyzer(TH1 *h, int ELo, int EHi, int EbkgLo, int EbkgHi){
  fDoRand = true;
  fSrce = h->GetDirectory()->GetMotherDir();
  fBeam = h->GetDirectory()->GetMotherDir();

  std::string dirname = std::string(h->GetDirectory()->GetName());
  fMassNum = TString(dirname.substr(dirname.find("_gated")-4,4));
  if(EHi < ELo) std::swap(ELo,EHi);
  fEbeamGate[0] = ELo;
  fEbeamGate[1] = EHi;
  fEbkgGate[0] = EbkgLo;
  fEbkgGate[1] = EbkgHi;
  fEsrceGate[0] = ELo;
  fEsrceGate[1] = EHi;
}

void GPolAnalyzer::GetXiRatio(GH1D *&hnorm, GH1D *&hsrce, GH1D *&hbeam, int binning, TString opt) {
  bool linFit = opt.Contains("lin");
  if (linFit) opt.ReplaceAll("lin","");
  bool saveFile = opt.Contains("save");
  if (saveFile) opt.ReplaceAll("save","");

  //background subtracion
  bool doBkgSubtract = false;
  if (!(fEbkgGate[0] == fEbkgGate[1] && fEbkgGate[0] == -1)) doBkgSubtract = true;

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

  //scaling for background subtraction
  double bkgScale = 0; 
  if (doBkgSubtract) 
    bkgScale = (double(fEbeamGate[1])-double(fEbeamGate[0]))/(double(fEbkgGate[1])-double(fEbkgGate[0]));
  GH1D *bkgtemp;
  GH1D *hbkg;

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

      //background
      if (doBkgSubtract) {
        bkgtemp = h2d->ProjectionX("_px",fEbkgGate[0],fEbkgGate[1]);
        bkgtemp->Sumw2();
      }

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

      //background
      if (doBkgSubtract) {
        bkgtemp->Add(h2d->ProjectionX("_px",fEbkgGate[0],fEbkgGate[1])); 
      }
    }
    
    //if a new group is coming, add the queue to the output
    if (i == N-1 || (srceCrystals[i].first != srceCrystals[i+1].first && srceCrystals[i].first < 1000) ) {
      //rebin
      if (binning != 1){
        stemp->Rebin(binning);
        dtemp->Rebin(binning);
        if (doBkgSubtract) bkgtemp->Rebin(binning);
      }

      //initialize output hists 
      if (nGroups == 0) {
        hsrce = (GH1D*) stemp->Clone("source"); 
        hbeam = (GH1D*) dtemp->Clone("data");
        if (doBkgSubtract) {
          hbkg = (GH1D*) bkgtemp->Clone("bkg");
          dtemp->Add(bkgtemp,-1*bkgScale);
        }
        dtemp->Divide(stemp);
        
        hnorm = (GH1D*) dtemp->Clone("norm");
      }
      else {
        hsrce->Add(stemp);
        hbeam->Add(dtemp);
        if (doBkgSubtract) {
          hbkg->Add(bkgtemp);
          dtemp->Add(bkgtemp,-1*bkgScale);
        }
        dtemp->Divide(stemp);
        
        hnorm->Add(dtemp);
      }

      readNew = true;
      nGroups++;
    }
  }

  if (doBkgSubtract) {
    hbkg->Scale(bkgScale);
    hbkg->SetLineColor(kGreen);
  }

  if (nGroups > 1){
    // scalingRatio /= nGroups;
    double Ndiff = hbeam->Integral();
    if (doBkgSubtract) Ndiff -= hbkg->Integral();
    hnorm->Scale(1.0/nGroups*hsrce->Integral()/Ndiff);
  } else {
    double Ndiff = hbeam->Integral();
    if (doBkgSubtract) Ndiff -= hbkg->Integral();
    hnorm->Scale(hsrce->Integral()/Ndiff);
  }
  
  //fit the normalized distribution
  TF1 *fitfunc;
  int parA0;
  if (linFit) {
    fitfunc = new TF1("mypol","([0] + [1]*x)*(1-[2]*TMath::Cos(2*x))",0,TMath::Pi());
    parA0 = 2;
  }
  else {
    fitfunc = new TF1("mypol","[0]*(1-[1]*TMath::Cos(2*x))",0,TMath::Pi());
    parA0 = 1;
    fitfunc->FixParameter(0,1.0);
  }

  // TFitResultPtr r1;
  // r1 = hnorm->Fit(fitfunc,"SQN0");

  // if (binning != 1){
  //   hnorm->Rebin(binning);
  //   hsrce->Rebin(binning);
  //   hbeam->Rebin(binning);
  //   if (doBkgSubtract) hbkg->Rebin(binning);
  //   hnorm->Scale(1.0/binning);
  // }

  // hsrce = hbkg;

  //fit again to show the function on chosen binning
  TFitResultPtr rbw;
  rbw = hnorm->Fit(fitfunc,"SQN0");
  hnorm->GetListOfFunctions()->Add(fitfunc);

  TFitResultPtr r1 = rbw;

  double midPoint = (energyLo+energyHi)/2;
  
  //output results
  printf("Grouping mode: %s\n",opt.Data());
  printf("|   Gate    |  MidPnt  |   Events  |\n");
  printf("| %4.0f:%4.0f |  %6.1f  |  %7.0f  |\n\n",energyLo,energyHi,midPoint,hbeam->Integral());
  printf("Bin Width  =       1 deg       |      %d deg   \n",binning);
  printf("Chi2/ndf   = %4.1f/%d = %3.2f  | %4.3f/%d = %3.2f\n",r1->Chi2(),r1->Ndf(),r1->Chi2()/r1->Ndf(),rbw->Chi2(),rbw->Ndf(),rbw->Chi2()/rbw->Ndf());
  printf("Const      = %5.4f +/- %5.4f | %5.4f +/- %5.4f\n",r1->Parameter(0),r1->ParError(0),rbw->Parameter(0),rbw->ParError(0));
  if (linFit) printf("Slope      = %5.4f +/- %5.4f | %5.4f +/- %5.4f\n",r1->Parameter(1),r1->ParError(1),rbw->Parameter(1),rbw->ParError(1));
  printf("A0         = %5.4f +/- %5.4f | %5.4f +/- %5.4f\n",r1->Parameter(parA0),r1->ParError(parA0),rbw->Parameter(parA0),rbw->ParError(parA0));
  printf("A0/err     =      %5.4f       |      %5.4f \n",r1->Parameter(parA0)/r1->ParError(parA0),rbw->Parameter(parA0)/rbw->ParError(parA0));
  printf("P          = %5.4f +/- %5.4f | %5.4f +/- %5.4f\n\n",2*r1->Parameter(parA0)/calcQ(midPoint),2*r1->ParError(parA0)/calcQ(midPoint),
                                                          2*rbw->Parameter(parA0)/calcQ(midPoint),2*rbw->ParError(parA0)/calcQ(midPoint));

  // printf("Beam BinError/BinContent = %4.3f/%4.3f = %4.3f%%\n",hbeam->GetBinError(2),hbeam->GetBinContent(2),hbeam->GetBinError(2)/hbeam->GetBinContent(2)*100);
  // printf("Source BinError/BinContent = %4.3f/%4.3f = %4.3f%%\n",hsrce->GetBinError(2),hsrce->GetBinContent(2),hsrce->GetBinError(2)/hsrce->GetBinContent(2)*100);
  // printf("Norm BinError/BinContent = %4.3f/%4.3f = %4.3f%%\n",hnorm->GetBinError(2),hnorm->GetBinContent(2),hnorm->GetBinError(2)/hnorm->GetBinContent(2)*100);

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
  /*
  double scaling = 1.0/(std::floor(hsrce->Integral()/hbeam->Integral())); 
  hsrce->Scale(scaling);
  hbeam->SetLineColor(kRed);
  double ymax = std::max(hbeam->GetMaximum(),hsrce->GetMaximum());
  double ymin = std::min(hbeam->GetMaximum(),hsrce->GetMaximum());
  hsrce->GetYaxis()->SetRangeUser(ymin*0.9,ymax*1.1);
  hbeam->GetYaxis()->SetRangeUser(0.0,hbeam->GetMaximum()*1.1);
  */
  double ymax = std::max(hbeam->GetMaximum(),hsrce->GetMaximum());
  hbeam->GetYaxis()->SetRangeUser(0.0,ymax*1.2);
  hsrce->GetYaxis()->SetRangeUser(0.0,ymax*1.2);
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
