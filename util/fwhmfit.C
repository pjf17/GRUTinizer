double dopplerbroad(double *x, double *p){
  double radX = x[0]*TMath::DegToRad();
  double sum = pow(p[0]*sin(radX)/(1-p[0]*cos(radX))*p[1],2) + 
                pow( (-p[0]+cos(radX))/((1-p[0]*p[0])*(1-p[0]*cos(radX)))*p[2] ,2) +
                pow(p[3]*cos(radX),2);

  return TMath::Sqrt(sum);
}

void readDataFile(std::string filename, std::vector<double> &vx, std::vector<double> &vy, std::vector<double> &vxerr, std::vector<double> &vyerr){
    std::ifstream input(filename.c_str());
    std::string line;
    double x, y, xerr, yerr;
    while (getline(input,line)){
        if (line.find('#') != std::string::npos) continue;
        std::stringstream ss(line);
        ss >> x >> y >> xerr >> yerr;
        vx.push_back(x*TMath::RadToDeg()); vy.push_back(y);
        vxerr.push_back(xerr); vyerr.push_back(yerr);
    }

    return;
}

void fwhmfit(std::string filename, double beta) {
    std::vector<double> X, Y, Xerr, Yerr;
    readDataFile(filename,X,Y,Xerr,Yerr);

    TGraphErrors *gr = new TGraphErrors((int) X.size(),&X[0],&Y[0],&Xerr[0],&Yerr[0]);
    TF1 *fitfunc = new TF1("fitfunc",dopplerbroad,0,180,4);
    fitfunc->FixParameter(0,beta);
    fitfunc->SetParameter(1,0.002);
    fitfunc->SetParameter(2,0.001);
    fitfunc->SetParameter(3,0.001);

    gr->Fit(fitfunc);
    gr->Draw("A*");
}