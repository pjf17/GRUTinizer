#include <vector>
#include <string>

#include "TH1.h"
#include "GCanvas.h"
#include "GH1D.h"
#include "TFile.h"


double calculateCenterOfMass(GH1D *h, double x_low, double x_high){
    //Calculate center of mass as sum over bins of (x*bin_content)/sum over bins

    double numerator = 0;
    double denominator = 0;

    int start_bin = h->FindBin(x_low);
    int end_bin = h->FindBin(x_high);

    for (int bin = start_bin; bin < end_bin; bin++){
        numerator += h->GetBinContent(bin)*h->GetBinCenter(bin);
        denominator += h->GetBinContent(bin);
    }

    if (denominator){
        return numerator/denominator;
    }

    else {
        std::cout << "Failed to calculate center of mass; denominator is 0" << std::endl;
        return -1;
    }
}

void check_y_correction(std::vector<std::string> input_file_names, 
                        std::string hist_name_crdc1, std::string hist_name_crdc2){
    //input_file_names are the files containing some histogram (hist_name)
    //of the CRDC Y Data. 
    std::vector<TFile*> files;
    std::vector<GH1D*> crdc1_hists;
    std::vector<GH1D*> crdc2_hists;

    for (auto &file_name : input_file_names){
        files.push_back(new TFile(file_name.c_str(), "read")); 
    }

    for (auto &file : files) {
        crdc1_hists.push_back((GH1D*)((GH2D*)file->Get(hist_name_crdc1.c_str()))); 
        crdc2_hists.push_back((GH1D*)((GH2D*)file->Get(hist_name_crdc2.c_str()))); 
    }


    GCanvas *crdc1_can = new GCanvas("crdc1_can", "crdc1_can", 1000,800);
    crdc1_can->cd();
    double crdc1_max = 0;
    for (auto &hist : crdc1_hists){
        hist->Rebin(2);
        hist->Draw("same");
        hist->SetLineWidth(3);
        if (crdc1_max < hist->GetMaximum()){
            crdc1_max = hist->GetMaximum();
        }

    }
    TLine l1(0,0,0,crdc1_max);
    l1.Draw();
    GCanvas *crdc2_can = new GCanvas("crdc2_can", "crdc2_can", 1000,800);
    crdc2_can->cd();
    double crdc2_max = 0;
    for (auto &hist : crdc2_hists){
        hist->Draw("same");
        hist->SetLineWidth(3);
        hist->Rebin(2);
        if (crdc2_max < hist->GetMaximum()){
            crdc2_max = hist->GetMaximum();
        }
    }
    TLine l2(0,0,0,crdc2_max);
    l2.Draw();
    
    //Plot all spectra on top of each other to show centered at 0
}

//  void run_first_set(){
//      std::vector<std::string> file_names = {
//          "hist0012.root",
//          "hist0013.root",
//          "hist0014.root",
//          "hist0015.root",
//          "hist0016.root",
//          "hist0017.root",
//          "hist0018.root",
//          "hist0019.root",
//          "hist0020.root"
//      };

//      double fit_low_x = -50;
//      double fit_high_x = 50;
//      std::string gvalue_file_name = "/projects/e15127/cfg_files/val_files/cu71_1p/cu71_1p_runs_012_020_init.val";
//      std::string hist_name_crdc1 = "s800_gated/CRDC1_Y_incut_top_blob_ni70fixed";
//      std::string hist_name_crdc2 = "s800_gated/CRDC2_Y_incut_top_blob_ni70fixed";
//      do_crdc_y_shift(file_names, hist_name_crdc1, hist_name_crdc2, gvalue_file_name, fit_low_x, fit_high_x);
//  }
//void run_all_cu71(){
//    std::vector<std::string> file_names = {
//        "hist0012.root",
//        "hist0013.root",
//        "hist0014.root",
//        "hist0015.root",
//        "hist0016.root",
//        "hist0017.root",
//        "hist0018.root",
//        "hist0019.root",
//        "hist0020.root",
//        "hist0096.root",
//        "hist0097.root",
//        "hist0098.root",
//        "hist0099.root",
//        "hist0100.root",
//        "hist0101.root",
//        "hist0103.root",
//        "hist0104.root",
//        "hist0106.root",
//        "hist0107.root",
//        "hist0108.root",
//        "hist0109.root",
//        "hist0110.root",
//        "hist0111.root",
//        "hist0112.root",
//        "hist0113.root",
//        "hist0114.root",
//        "hist0115.root",
//        "hist0116.root",
//        "hist0117.root",
//        "hist0118.root",
//        "hist0119.root",
//        "hist0121.root",
//        "hist0122.root",
//        "hist0123.root",
//        "hist0124.root",
//        "hist0125.root",
//        "hist0126.root",
//        "hist0128.root",
//        "hist0129.root",
//        "hist0130.root",
//        "hist0131.root",
//        "hist0132.root",
//        "hist0133.root",
//        "hist0134.root",
//        "hist0135.root",
//        "hist0136.root",
//        "hist0137.root",
//        "hist0138.root",
//        "hist0140.root",
//        "hist0141.root",
//        "hist0142.root",
//    };
//
//    plot_obj(file_names);
//}
void run_second_set(){

    std::vector<std::string> file_names = {
        "hist0096.root",
        "hist0097.root",
        "hist0098.root",
        "hist0099.root",
        "hist0100.root",
        "hist0101.root",
        "hist0103.root",
        "hist0104.root",
        "hist0106.root",
        "hist0107.root",
        "hist0108.root",
        "hist0109.root",
        "hist0110.root",
        "hist0111.root",
        "hist0112.root",
        "hist0113.root",
        "hist0114.root",
        "hist0115.root",
        "hist0116.root",
        "hist0117.root",
        "hist0118.root",
        "hist0119.root",
        "hist0121.root",
        "hist0122.root",
        "hist0123.root",
        "hist0124.root",
        "hist0125.root",
        "hist0126.root",
        "hist0128.root",
        "hist0129.root",
        "hist0130.root",
        "hist0131.root",
        "hist0132.root",
        "hist0133.root",
        "hist0134.root",
        "hist0135.root",
        "hist0136.root",
        "hist0137.root",
        "hist0138.root",
        "hist0140.root",
        "hist0141.root",
        "hist0142.root",
    };

    std::string hist_name_crdc1 = "s800_gated/CRDC1_Y_incut_top_blob_ni70_tight";
    std::string hist_name_crdc2 = "s800_gated/CRDC2_Y_incut_top_blob_ni70_tight";
    check_y_correction(file_names, hist_name_crdc1, hist_name_crdc2);
}

