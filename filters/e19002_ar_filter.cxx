// To make a new filter, copy this file under a new name in the "filter" directory.
// The "FilterCondition" function should return a boolean value.
// The boolean indicates whether the event should be kept or not.

#include "TRuntimeObjects.h"
#include "GCutG.h"
#include "TS800.h"
#include "TFile.h"

double GetGoodICE(TS800 *s800){
  static int ncalls = 0;
  double value = 0;
  double crdc_1_x = s800->GetCrdc(0).GetDispersiveX();
  double crdc_1_y = s800->GetCrdc(0).GetNonDispersiveY();
  
  double xtilt = GValue::Value("IC_DE_XTILT");
  double x0tilt = GValue::Value("IC_DE_X0TILT");
  double ytilt = GValue::Value("IC_DE_YTILT");
  if (!std::isnan(xtilt) && !std::isnan(x0tilt) && !std::isnan(ytilt)){
    value = s800->GetIonChamber().GetdE(crdc_1_x, crdc_1_y);
  } else {
    value = s800->GetIonChamber().GetAve();
    if (ncalls == 0){
      std::cout<<"XTILT, X0TILT, YTILT NOT SET SWITCHING TO GETAVE()\n";
      ncalls++;
    }
  }
  
  return value;
}

int gates_loaded = 0;
std::vector<GCutG*> outgoing_cuts;

std::string CUT_OF_INTEREST_NAME = "allAr";
std::string CUT_FILE_NAME = "/mnt/analysis/pecan-2015/farris/e19002/config/allAr_filter.cuts";

TFile *cut_file = 0;
GCutG *cut_of_interest = 0;
extern "C"
bool FilterCondition(TRuntimeObjects& obj) {

    if (!cut_file){
        cut_file = new TFile(CUT_FILE_NAME.c_str(), "read");
        if (!cut_file){
            std::cout << "Failed to load cut file: " << CUT_FILE_NAME << std::endl;
            exit(1);
        }

        if (!cut_of_interest){
            cut_of_interest = (GCutG*)cut_file->Get(CUT_OF_INTEREST_NAME.c_str());
        }

        if (!cut_of_interest){
            std::cout << "Failed to load cut of interest!" << std::endl;
            exit(1);
        }
    }

    TS800    *s800    = obj.GetDetector<TS800>();
    if (!s800){
        return false;
    }
    double ic = GetGoodICE(s800);

    if (cut_of_interest->IsInside(s800->GetMTofObjE1(),ic)){ 
        return true;
    }
    return false;
}
