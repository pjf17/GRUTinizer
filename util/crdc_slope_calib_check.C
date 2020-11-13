#include <cmath>
#include <TF1.h>
#include <TFile.h>
#include <TCanvas.h>

TFile *outfile;
const std::string VAL_FILE_DIR = "/mnt/analysis/pecan-2015/farris/e19002/config/crdcY/"; 

GH1D *GetYHist(std::string filename, int crdc){
  TFile *f = new TFile(filename.c_str(),"READ");
  
  std::string inHistName = Form("ungated/crdc%d X_Y",crdc);
  GH2D *hist2d = (GH2D*) f->Get(inHistName.c_str());

  if (!hist2d){
    std::cout<<"error reading hist"<<std::endl;
    return NULL;
  }

  outfile->cd();
  
  GH1D *hy = new GH1D(*hist2d->ProjectionY());
  std::string histname = Form("run%s_crdc%d_Y",filename.substr(4,4).c_str(),crdc);
  hy->SetNameTitle(histname.c_str(),histname.c_str());
  
  delete hist2d;
  f->Close();
  delete f;
  
  return hy;
}

double GetMean(GH1D *hy){
  if (!hy) {
    return sqrt(-1);
  }
  
  double low = -600;
  double high = 600;

  std::string fitname(hy->GetName());
  fitname += "_fit";

  TF1 *fitfunc = new TF1(fitname.c_str(),"gaus",low,high);
  hy->Fit(fitfunc,"QR0","",low, high);

  const Int_t kNotDraw = 1<<9;
  hy->GetFunction(fitname.c_str())->ResetBit(kNotDraw);

  hy->Write("",TObject::kOverwrite);
  
  return fitfunc->GetParameter(1);
}

void CheckYDists(std::vector<std::string> allfiles){
    outfile = new TFile("crdc_calib_hists_check.root","RECREATE");
    ofstream statusOut; 
    statusOut.open("calib_check.txt");

    double good = 0.01, okay = 0.1, fine = 0.25;
    for (auto filename : allfiles){
        statusOut<<filename<<std::endl; 
        for (int i=0; i < 2; i++){
            statusOut<<"\tcrdc"<<i+1<<": ";
            //get mean and check if its nan
            double mean = GetYHist(filename,i+1)->GetMean();
            mean = std::abs(mean);
            std::cout<<mean<<std::endl;
            if (std::isnan(mean)){
                statusOut<<"Error";
                break;
            }

            //check if mean is close to zero
            if (mean < good){
                statusOut<<"good";
            } else if (mean < okay){
                statusOut<<"okay";
            } else if (mean < fine){
                statusOut<<"fine";
            } else {
                statusOut<<"BAD";
            }
            statusOut<<std::endl;
        } 
    }
    statusOut.close();
    return;
}

bool ReadFiles(std::vector<std::string> &allfiles){
  ifstream input("list_o_hists.txt");

  if (!input.is_open()){
    std::cout<<"Error: cannot open files list\n";
    return false;
  }

  std::string line;
  while ( getline(input,line) ){
    allfiles.push_back(line);
  }
  
  return true;
}

void crdc_slope_calib_check(){
  std::vector<std::string> files;
  if (!ReadFiles(files)) return;
  CheckYDists(files);
  return;
}
