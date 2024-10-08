#include <cmath>
#include <TF1.h>
#include <TFile.h>
#include <TCanvas.h>

//NOTE FOR THIS TO WORK PROPERLY YOU NEED TO CREATE A TXT FILE WITH ALL THE HISTOGRAM NAMES IN IT
//TO DO THIS JUST ENTER 'ls hist*.root > name_of_your_choice.txt' IN THE DIRECTORY WITH YOUR HISTS

TFile *outfile;

GH1D *GetYHist(std::string filename, std::string histdir, int crdc){
  TFile *f = new TFile(filename.c_str(),"READ");
  
  std::string inHistName = Form("%s/crdc%d X_Y",histdir.c_str(),crdc);
  GH2D *hist2d = (GH2D*) f->Get(inHistName.c_str());
  // printf("%p\n",&hist2d);
  if (!hist2d){
    std::cout<<"error reading hist in "<<filename<<std::endl;
    return NULL;
  }

  std::string dirname = "run" + filename.substr(4,4);
  outfile->cd(dirname.c_str());
  
  GH1D *hy = hist2d->ProjectionY("_py",0,-1,"");
  std::string histname = Form("crdc%d_Y",crdc);
  std::string histtitle = dirname + histname;
  hy->SetNameTitle(histname.c_str(),histtitle.c_str());
  
  delete hist2d;
  f->Close();
  delete f;

  return hy;
}

double GetMean(GH1D *hy/* double low, double high*/){
  if (!hy) {
    return sqrt(-1);
  }
  
  std::string fitname(hy->GetName());
  fitname += "_fit";

  double mean = hy->GetMean();

  if (hy->GetEntries() > 200){
    TF1 *fitfunc = new TF1(fitname.c_str(),"gaus",mean*0.5,mean*1.5);
    hy->Fit(fitfunc,"QR0","",mean*0.5,mean*1.5);
    const Int_t kNotDraw = 1<<9;
    hy->GetFunction(fitname.c_str())->ResetBit(kNotDraw);
    hy->Write("",TObject::kOverwrite);
    mean = fitfunc->GetParameter(1);
  }
  
  return mean;
}

void GetNewSlopes(std::vector<std::string> &allfiles, std::string histdir, std::string val_file_dir){
  outfile = new TFile("crdc_calib_hists.root","RECREATE");
  int nloops = 1;
  int nfiles = (int) allfiles.size();
  for (auto filename : allfiles){
    std::cout<<"\rProcessing "<<nloops<<"/"<<nfiles<<std::flush;
    //create run directory for crdcY 1 & 2
    std::string dirname = "run" + filename.substr(4,4);
    if (!outfile->GetDirectory(dirname.c_str())){
      outfile->mkdir(dirname.c_str());
    }

    //create .val file
    std::string valfile = val_file_dir + "/" + dirname + ".val";
    ofstream out_valfile(valfile.c_str());
  
    for (int i=0; i < 2; i++){
      //get mean and check if its nan
      double mean = GetMean(GetYHist(filename,histdir,i+1));
      // double mean = GetYHist(filename,i+1)->GetMean();
      if (std::isnan(mean)){
        std::cout<<"Error in crdc "<<i+1<< "mean\n";
        return;
      }
      // get initial slope
      // std::string gSlope = Form("CRDC%d_Y_SLOPE",i+1);
      // double old_slope = GValue::Value(gSlope.c_str());
      // if (std::isnan(old_slope)){
	    //   std::cout<<"No old slope for crdc"<<i+1<<", setting to 1"<<std::endl;
      //   old_slope = 1.0;
      // }

      //get offset
      std::string gOffset = Form("CRDC%d_Y_OFFSET",i+1);
      double offset = GValue::Value(gOffset.c_str());
      if (std::isnan(offset)){
	      std::cout<<"No offset for crdc"<<i+1<<std::endl;
	      break;
      }

      //undo old calibration
      // double uncal_mean = (mean - offset)/old_slope;
      
      //find new slope such that the new mean is zero
      double new_slope = -(offset/mean);
      out_valfile << "CRDC" << i+1 << "_Y_SLOPE {\n  value: " << new_slope << "\n}\n\n";
    }
  
    out_valfile.close();
    nloops++;
  }
  std::cout<<std::endl;

  outfile->Write();
  return;
}

void CheckYDists(std::vector<std::string> allfiles, std::string histdir){
    outfile = new TFile("crdc_calib_hists_check.root","RECREATE");

    TH1D *hmeans = new TH1D("means","means",500,-2.5,2.5);
    TH2D *hindex = new TH2D("run_summary","run_summary",350,0,350,500,-2.5,2.5);
    int nloops = 1;
    int nfiles = (int) allfiles.size();
    for (auto filename : allfiles){
        //output to user
        std::cout<<"\rProcessing "<<nloops<<"/"<<nfiles<<std::flush;
        nloops++;

        //create run directory for crdcY 1 & 2
        std::string dirname = "run" + filename.substr(4,4);
        if (!outfile->GetDirectory(dirname.c_str())){
          outfile->mkdir(dirname.c_str());
        }

        int runNum = atoi(filename.substr(4,4).c_str());

        for (int i=0; i < 2; i++){
            //get mean and check if its nan
            double mean = GetMean(GetYHist(filename,histdir,i+1)/*,-60,60*/);
            // double mean = GetYHist(filename,i+1)->GetMean();
            if (std::isnan(mean)){
                std::cout<<"Error with hist"<<runNum<<std::endl;
                break;
            }
            hmeans->Fill(mean);
            hindex->Fill(runNum,mean);
        } 
    }
    std::cout<<std::endl;
    outfile->Write();
    return;
}

bool ReadFiles(std::string filelist, std::vector<std::string> &allfiles){
  ifstream input(filelist.c_str());

  if (!input.is_open()){
    std::cout<<"Error: cannot open files list\n";
    return false;
  }

  std::string line;
  while ( getline(input,line) ){
    allfiles.push_back(line);
  }
  
  input.close();
  return true;
}

void crdc_slope_calib(std::string mode, std::string val_file_dir = "", std::string histdir="ungated", std::string filelist="list_o_hists.txt"){
  std::vector<std::string> files;
  if (!ReadFiles(filelist, files)) return;
  
  if (histdir.compare("ungated") == 0)
    std::cout<<"No histogram directory specified, assuming ungated directory"<<std::endl;

  if (mode == "calib"){
    if (val_file_dir == ""){
      std::cout<<"Must enter .val file output path"<<std::endl;
    } else {
      GetNewSlopes(files,histdir,val_file_dir);
    }
  } else if (mode == "check"){
    CheckYDists(files,histdir);
  } else {
    std::cout<<"Enter either calib, or check as your mode\n";
  }
  return;
}
