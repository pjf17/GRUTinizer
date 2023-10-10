#include "TMath.h"
#include "TF1.h"

double LegendreDoppler (double *x, double *par){
  double cosThCM = (std::cos(x[0]*TMath::DegToRad()) - par[3])/(1 - par[3]*std::cos(x[0]*TMath::DegToRad()));
  double leg2 = 0.5*(3*pow(cosThCM,2) - 1);
  double leg4 = 0.125*(35*pow(cosThCM,4) - 30*pow(cosThCM,2) +3);
  return par[0]*(1 + par[1]*leg2 + par[2]*leg4)*(1-par[3]*par[3])/pow(par[3]*TMath::Cos(x[0]*TMath::DegToRad())-1,2);
}

double polarization(double *x, double *par){
  double cosx = std::cos(x[0]*TMath::DegToRad());
  double leg2 = 0.5*(3*pow(cosx,2) - 1);
  double leg4 = 0.125*(35*pow(cosx,4) - 30*pow(cosx,2) +3);
  double leg22 = 3*(1 - pow(cosx,2));
  double leg42 = 15./2*(7*pow(cosx,2) - 1)*(1 - pow(cosx,2));
  return (1./2*par[0]*leg22 - 1./12*par[1]*leg42)/(1 + par[0]*leg2 + par[1]*leg4);
}

double thetaCM(double theta,double beta){
  double cosT = TMath::Cos(theta);
  return TMath::ACos((cosT - beta)/(1 - beta*cosT));
}

TF1 *flegendre = new TF1("fitleg",LegendreDoppler,45,95,4);
TF1 *fPol = new TF1("fpol",polarization,0,180,2);