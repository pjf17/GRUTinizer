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

void do_crdc_y_shift(std::vector<std::string> input_file_names, std::string hist_name_crdc1, 
        std::string hist_name_crdc2, std::string gvalue_file_name, double fit_low_x, double fit_high_x){
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

    TF1 *mygaus = new TF1("mygaus","gaus");
    std::vector<double> com_crdc1;
    for (auto &hist : crdc1_hists){
        hist->Rebin(2);
            hist->Fit(mygaus, "ME", "", fit_low_x, fit_high_x);
        com_crdc1.push_back(mygaus->GetParameter(1));
        //com_crdc1.push_back(calculateCenterOfMass(hist, fit_low_x, fit_high_x));
    }
    std::vector<double> com_crdc2;
    for (auto &hist : crdc2_hists){
        hist->Rebin(2);
        hist->Fit(mygaus, "ME", "", fit_low_x, fit_high_x);
        com_crdc2.push_back(mygaus->GetParameter(1));
        //com_crdc2.push_back(calculateCenterOfMass(hist, fit_low_x, fit_high_x));
    }

    std::cout << "Input File Size: " << input_file_names.size() << std::endl;
    std::cout << "COM CRDC1 Size: " << com_crdc1.size() << std::endl;
    std::cout << "COM CRDC2 Size: " << com_crdc2.size() << std::endl;

    //Read in current GValues
    GValue::ReadValFile(gvalue_file_name.c_str());

    double crdc1_y_slope = GValue::Value("CRDC1_Y_SLOPE");
    double crdc2_y_slope = GValue::Value("CRDC2_Y_SLOPE");

    double crdc1_y_offset = GValue::Value("CRDC1_Y_OFFSET");
    double crdc2_y_offset = GValue::Value("CRDC2_Y_OFFSET");

    //Update Y_SLOPE variables to center runs at 0.
    std::vector<double> new_slopes_crdc1;
    double base_com_crdc1 = com_crdc1.at(0);
    double base_tac_crdc1 = (base_com_crdc1-crdc1_y_offset)/crdc1_y_slope;
    for (int i = 1; i < com_crdc1.size(); i++){
        double tac_crdc1 = (com_crdc1.at(i)-crdc1_y_offset)/ crdc1_y_slope;
        double new_slope  = (base_tac_crdc1/tac_crdc1) * crdc1_y_slope;
        new_slopes_crdc1.push_back(new_slope);
    }
    std::cout << "Finished getting new slopes for CRDC1" << std::endl;
    std::vector<double> new_slopes_crdc2;
    double base_com_crdc2 = com_crdc2.at(0);
    double base_tac_crdc2 = (base_com_crdc2-crdc2_y_offset)/crdc2_y_slope;
    for (int i = 1; i < com_crdc2.size(); i++){
        double tac_crdc2 = (com_crdc2.at(i)-crdc2_y_offset)/ crdc2_y_slope;
        double new_slope  = (base_tac_crdc2/tac_crdc2) * crdc2_y_slope;
        new_slopes_crdc2.push_back(new_slope);
    }
    std::cout << "Finished getting new slopes for CRDC2" << std::endl;

    std::cout << "Rechecking file name size: " << input_file_names.size() << std::endl;
    for (int i = 1; i < input_file_names.size(); i++){
        std::string new_fn = "cu71_1p_run"+input_file_names.at(i).substr(4,4)+".val";
        std::cout << "\tCreating new file: " << new_fn << "\n";
        GValue::SetReplaceValue("CRDC1_Y_SLOPE", new_slopes_crdc1.at(i-1));
        GValue::SetReplaceValue("CRDC2_Y_SLOPE", new_slopes_crdc2.at(i-1));
        GValue::WriteValFile(new_fn);
    }
    std::cout << "Finished writing new files" << std::endl;


    //Plot all spectra on top of each other to show centered at 0
    //  GCanvas *c = new GCanvas();
    //  TLine l1(0,0,0,c->GetMaximum());

    //  c->cd();
    //  for (auto &hist : hists){
    //      hist->Draw("same");
    //      hist->SetLineWidth(3);
    //      l1.Draw();
    //  }
}

void run_first_set(){
    std::vector<std::string> file_names = {
        "hist0012.root",
        "hist0013.root",
        "hist0014.root",
        "hist0015.root",
        "hist0016.root",
        "hist0017.root",
        "hist0018.root",
        "hist0019.root",
        "hist0020.root"
    };

    double fit_low_x = -50;
    double fit_high_x = 50;
    std::string gvalue_file_name = "/projects/e15127/cfg_files/val_files/cu71_1p/cu71_1p_runs_012_020_init.val";
    std::string hist_name_crdc1 = "s800_gated/CRDC1_Y_inCu71_outNi70";
    std::string hist_name_crdc2 = "s800_gated/CRDC2_Y_inCu71_outNi70";
    do_crdc_y_shift(file_names, hist_name_crdc1, hist_name_crdc2, gvalue_file_name, fit_low_x, fit_high_x);
}
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

    //    plot_obj(file_names);
    double fit_low_x = -80;
    double fit_high_x = 80;
    std::string gvalue_file_name = "/projects/e15127/cfg_files/val_files/cu71_1p/cu71_1p_runs_096_138_init.val";

    std::string hist_name_crdc1 = "s800_gated/CRDC1_Y_in_spectcl_cu_blob_ni70_loose";
    std::string hist_name_crdc2 = "s800_gated/CRDC2_Y_in_spectcl_cu_blob_ni70_loose";
    do_crdc_y_shift(file_names, hist_name_crdc1, hist_name_crdc2, gvalue_file_name, fit_low_x, fit_high_x);
}

void run_ni71(){

    std::vector<std::string> file_names = {
        "hist0039.root",
        "hist0040.root",
        "hist0041.root",
        "hist0042.root",
        "hist0043.root",
        "hist0044.root",
        "hist0045.root",
        "hist0046.root",
        "hist0047.root",
        "hist0048.root",
        "hist0049.root",
        "hist0050.root",
        "hist0053.root",
    };

    //    plot_obj(file_names);
    double fit_low_x = -80;
    double fit_high_x = 80;
    std::string gvalue_file_name = "/projects/e15127/cfg_files/val_files/ni71_1n/ni71_1n_runs_039_053_init.val";
    std::string hist_name_crdc1 = "s800_gated/CRDC1_Y_incut_bot_blob_ni70";
    std::string hist_name_crdc2 = "s800_gated/CRDC2_Y_incut_bot_blob_ni70";
    do_crdc_y_shift(file_names, hist_name_crdc1, hist_name_crdc2, gvalue_file_name, fit_low_x, fit_high_x);
}

void run_zn72(){

    std::vector<std::string> file_names = {
        "hist0061.root",
        "hist0063.root",
        "hist0064.root",
        "hist0065.root",
        "hist0066.root",
        "hist0067.root",
        "hist0068.root",
        "hist0069.root",
        "hist0070.root",
        "hist0071.root",
        "hist0072.root",
        "hist0073.root",
        "hist0074.root",
        "hist0075.root",
        "hist0076.root",
        "hist0077.root",
        "hist0078.root",
        "hist0079.root",
        "hist0080.root",
        "hist0081.root",
        "hist0082.root",
        "hist0083.root",
        "hist0084.root",
        "hist0085.root",
        "hist0086.root"
    };

    //    plot_obj(file_names);
    double fit_low_x = -80;
    double fit_high_x = 80;
    std::string gvalue_file_name = "/projects/e15127/cfg_files/val_files/zn72_2p/zn72_2p_runs_061_086_init.val";
    std::string hist_name_crdc1 = "s800_gated/CRDC1_Y_incut_blob_ni70";
    std::string hist_name_crdc2 = "s800_gated/CRDC2_Y_incut_blob_ni70";
    do_crdc_y_shift(file_names, hist_name_crdc1, hist_name_crdc2, gvalue_file_name, fit_low_x, fit_high_x);
}

void run_zn72_e15020(){

    std::vector<std::string> file_names = {
        "hist0138.root",
        "hist0139.root",
        "hist0140.root",
        "hist0141.root",
        "hist0142.root",
        "hist0143.root",
        "hist0144.root",
        "hist0145.root"
    };

    //    plot_obj(file_names);
    double fit_low_x = -80;
    double fit_high_x = 80;
    std::string gvalue_file_name = "/projects/e15127/cfg_files/val_files/zn72_2p_e15020/zn72_2p_runs_138_146.val";
    std::string hist_name_crdc1 = "s800_gated/CRDC1_Y_in_zn72_out_Ni70_looser";
    std::string hist_name_crdc2 = "s800_gated/CRDC2_Y_in_zn72_out_Ni70_looser";
    do_crdc_y_shift(file_names, hist_name_crdc1, hist_name_crdc2, gvalue_file_name, fit_low_x, fit_high_x);
}
