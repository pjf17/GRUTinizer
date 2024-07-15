double dopplerbroadHI(double *x, double *p){
  double sum = pow(p[1]*sin(x[0])/(1-p[1]*cos(x[0]))*p[2],2) + 
                pow( (-p[1]+cos(x[0]))/((1-p[1]*p[1])*(1-p[1]*cos(x[0])))*p[3] ,2) +
                pow(p[4]*cos(x[0]),2);

  return p[0]*TMath::Sqrt(1-p[1]*p[1])/(1 - p[1]*TMath::Cos(x[0])) + p[0]*p[5]*TMath::Sqrt(sum);
}

double dopplerbroadLO(double *x, double *p){
  double sum = pow(p[1]*sin(x[0])/(1-p[1]*cos(x[0]))*p[2],2) + 
                pow( (-p[1]+cos(x[0]))/((1-p[1]*p[1])*(1-p[1]*cos(x[0])))*p[3] ,2) +
                pow(p[4]*cos(x[0]),2);

  return p[0]*TMath::Sqrt(1-p[1]*p[1])/(1 - p[1]*TMath::Cos(x[0])) - p[0]*p[5]*TMath::Sqrt(sum);
}

TF1 *fdopHI = new TF1("fdopHI",dopplerbroadHI,0,TMath::Pi(),6);
TF1 *fdopLO = new TF1("fdopLO",dopplerbroadLO,0,TMath::Pi(),6);

void draw_doppler_gate(double centroidEnergy,double beta,double dtheta,double dbeta,double beamspot, double scale){

    fdopHI->SetParameter(0,centroidEnergy);
    fdopHI->SetParameter(1,beta);
    fdopHI->SetParameter(2,dtheta);
    fdopHI->SetParameter(3,dbeta);
    fdopHI->SetParameter(4,beamspot);
    fdopHI->SetParameter(5,scale);

    fdopLO->SetParameter(0,centroidEnergy);
    fdopLO->SetParameter(1,beta);
    fdopLO->SetParameter(2,dtheta);
    fdopLO->SetParameter(3,dbeta);
    fdopLO->SetParameter(4,beamspot);
    fdopLO->SetParameter(5,scale);

    fdopHI->Draw("same");
    fdopLO->Draw("same");
}