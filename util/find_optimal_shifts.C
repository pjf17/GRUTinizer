#include <map>
#include <vector>
#include <fstream>
#include <sstream>
#include "TMinuit.h"

TGretina *gret = new TGretina();
std::vector<std::pair<int,double>> energies;

std::map<int,int> detMapRing = {
  {26, 0}, {30, 1}, {34, 2}, {38, 3}, {25, 4}, {29, 5}, {33, 6}, {37, 7},
  {27, 8}, {31, 9}, {35,10}, {39,11}, {24,12}, {28,13}, {32,14}, {36,15},
  {47,16}, {63,17}, {71,18}, {79,19}, {51,20}, {59,21}, {67,22}, {83,23},
  {50,24}, {58,25}, {66,26}, {82,27}, {44,28}, {60,29}, {68,30}, {76,31},
  {46,32}, {62,33}, {70,34}, {78,35}, {48,36}, {56,37}, {64,38}, {80,39},
  {49,40}, {57,41}, {65,42}, {81,43}, {45,44}, {61,45}, {69,46}, {77,47}
};

std::map<int,int> invDetMap = {
  { 0,26}, { 1,30}, { 2,34}, { 3,38}, { 4,25}, { 5,29}, { 6,33}, { 7,37},
  { 8,27}, { 9,31}, {10,35}, {11,39}, {12,24}, {13,28}, {14,32}, {15,36},
  {16,47}, {17,63}, {18,71}, {19,79}, {20,51}, {21,59}, {22,67}, {23,83},
  {24,50}, {25,58}, {26,66}, {27,82}, {28,44}, {29,60}, {30,68}, {31,76},
  {32,46}, {33,62}, {34,70}, {35,78}, {36,48}, {37,56}, {38,64}, {39,80},
  {40,49}, {41,57}, {42,65}, {43,81}, {44,45}, {45,61}, {46,69}, {47,77}
};

double readInitFile(std::string filename){
    std::ifstream inFile(filename);
    std::string line;
    getline(inFile,line);
    double beta = std::stod(line); 
    while (std::getline(inFile,line)){
        std::stringstream ss = std::stringstream(line);
        int idx; double eng;
        ss >> idx >> eng;
        // std::cout<<idx<<" "<<eng<<std::endl;
        energies.push_back(std::make_pair(idx,eng));
    }
    inFile.close();
    return beta;
}

double invDoppler(double E, double theta, double beta){
    double gamma = 1./(sqrt(1.-pow(beta,2.)));
    return E / (gamma *(1 - beta*TMath::Cos(theta)));
}

void labFrame(std::vector<std::pair<int,double>> &eng,double beta){
    TGretina *gret = new TGretina();
    // double x = TMath::Sin(-0.005236);
    // double y = x;
    double x = 0.0;
    double y = 0.0;
    TVector3 beam = TVector3(x,-y,TMath::Sqrt(1 - x*x - y*y));
    int ncrystals = (int) energies.size();
    for (int c=0; c < ncrystals; c++){
        TVector3 crstl_pos = gret->GetCrystalPosition(invDetMap[eng[c].first]);
        // TVector3 crstl_pos = gret->GetCrystalPosition(eng[c].first);
        eng[c].second = invDoppler(eng[c].second,crstl_pos.Angle(beam),beta);
    }
    return;
}

double dopplerCorrect(int idx, double beta, double ataShift, double btaShift, double xshift, double yshift, double zshift){
    TVector3 track = TVector3(ataShift,-btaShift,sqrt(1-ataShift*ataShift-btaShift*btaShift));
    TVector3 gret_pos = gret->GetCrystalPosition(idx);
    gret_pos.SetXYZ(gret_pos.x() - xshift,
                    gret_pos.y() - yshift,
                    gret_pos.z() - zshift);
    
    double gamma = 1./(sqrt(1.-pow(beta,2.)));
    return gamma*(1 - beta*TMath::Cos(gret_pos.Angle(track)));
}

//MINIMIZER FUNCTIONS
void fStdev(int &npar, double *gin, double &f, double *par, int iflag){
    double enAvg = 0;
    int ncrystals = (int) energies.size();
    std::vector<double> dop_en;
    //make list and find avg
    for (int c=0; c < ncrystals; c++){
        dop_en.push_back(energies[c].second*dopplerCorrect(invDetMap[energies[c].first],par[0],par[1],par[2],par[3],par[4],par[5]));
        // dop_en.push_back(energies[c].second*dopplerCorrect(energies[c].first,par[0],par[1],par[2],par[3],par[4],par[5]));
        enAvg += dop_en.back();
    }
    
    enAvg /= ncrystals;

    //calc stdev
    double stdev = 0;
    for (int c=0; c < ncrystals; c++){
        stdev += (enAvg - dop_en[c])*(enAvg - dop_en[c]);
    }
    stdev /= ncrystals;
    f = TMath::Sqrt(stdev);
}

//for when you know the value of the gamma ray
void fKnown(int &npar, double *gin, double &f, double *par, int iflag){
    int ncrystals = (int) energies.size();
    std::vector<double> dop_en;
    double chi2 = 0;
    
    //calc chi2
    for (int c=0; c < ncrystals; c++){
        //par[6] is the known energy
        double diff = par[6] - energies[c].second*dopplerCorrect(invDetMap[energies[c].first],par[0],par[1],par[2],par[3],par[4],par[5]);
        chi2 += diff*diff;
    }
    f = chi2;
}

void printFitenergies(std::string filename, double *par){
    FILE *fp; 
    filename = filename.substr(0,filename.find("."));
    filename = filename + "-eng.txt"; 
    fp = fopen(filename.c_str(),"w");
    int ncrystals = (int) energies.size();
    for (int c=0; c < ncrystals; c++){
        fprintf(fp,"%f\n",energies[c].second*dopplerCorrect(invDetMap[energies[c].first],par[0],par[1],par[2],par[3],par[4],par[5]));
        // printf("%f\n",energies[c].second*dopplerCorrect(energies[c].first,par[0],par[1],par[2],par[3],par[4],par[5]));
    }
    fclose(fp);
    return;
}

void outputValFile(std::string filename, double *par) {
    FILE *fp; 
    filename = filename.substr(0,filename.find("."));
    filename = filename + "-out.val"; 
    fp = fopen(filename.c_str(),"w");
    std::vector<std::string> variable_names = {"BETA","ATA_SHIFT","BTA_SHIFT","TARGET_X_OFFSET","TARGET_Y_OFFSET","TARGET_Z_OFFSET"};
    for (int i=0; i < 6; i++) fprintf(fp,"%s {\n  Value: %f\n}\n\n",variable_names[i].c_str(),par[i]);
    fclose(fp);
    return;
}

void drawBeforeAfter(double oldbeta, double *par) {
   std::vector<double> x; 
   std::vector<double> enBefore; 
   std::vector<double> enAfter; 
   int ncrystals = (int) energies.size();
   for (int c=0; c < ncrystals; c++){
    x.push_back(c);
    enBefore.push_back(energies[c].second*dopplerCorrect(invDetMap[energies[c].first],oldbeta,0,0,0,0,0));
    enAfter.push_back(energies[c].second*dopplerCorrect(invDetMap[energies[c].first],par[0],par[1],par[2],par[3],par[4],par[5]));
   }

   TGraph *grBefore = new TGraph((int) enBefore.size(), &x[0], &enBefore[0]);
   TGraph *grAfter =  new TGraph((int) enAfter.size(), &x[0], &enAfter[0]);

   TCanvas *canv = new TCanvas();
   canv->SetGrid();
   grAfter->SetMarkerColor(kRed);
   grAfter->SetLineColor(kRed);
   grBefore->Draw("AL*");
   grAfter->Draw("L*same");
}

void find_optimal_shifts(std::string filename, bool varyBeta = false, double knownEnergy = -1){
    double oldbeta = readInitFile(filename);
    labFrame(energies,oldbeta);

    //initialize
    int npars = 6;
    TMinuit *min;
    if (knownEnergy == -1) 
        min = new TMinuit(npars);
    else 
        min = new TMinuit(npars+1);

    min->DefineParameter(0,"beta",oldbeta,varyBeta*0.001,0.3,0.5);
    min->DefineParameter(1,"ata_shift",0.001,1E-3,-0.5,0.5);
    min->DefineParameter(2,"bta_shift",0.001,1E-3,-0.5,0.5);
    min->DefineParameter(3,"x_shift",0,0.01,-10.0,10.0);
    min->DefineParameter(4,"y_shift",0,0.01,-10.0,10.0);
    min->DefineParameter(5,"z_shift",0,0.01,-10.0,10.0);
    
    if (knownEnergy == -1)
        min->SetFCN(fStdev);
    else {
        min->DefineParameter(6,"known_energy",knownEnergy,0,0,8000);
        min->SetFCN(fKnown);
    }

    //do minimization
    double pars[6], parerrs[6];
    min->Migrad();
    for (int p=0; p < npars; p++){
        min->GetParameter(p,pars[p],parerrs[p]);
    }

    printFitenergies(filename,pars);
    outputValFile(filename,pars);
    drawBeforeAfter(oldbeta,pars);
    energies.clear();
    return;
}