// To make a new filter, copy this file under a new name in the "filter" directory.
// The "FilterCondition" function should return a boolean value.
// The boolean indicates whether the event should be kept or not.

#include "TRuntimeObjects.h"
#include "GCutG.h"
#include "TS800.h"
#include "TFile.h"

int gates_loaded = 0;
std::vector<GCutG*> outgoing_cuts;

std::string CUT_OF_INTEREST_NAME = "outAllAr";
std::string CUT_FILE_NAME = "/mnt/analysis/pecan-2015/farris/e19002/ar46/ar46_cuts.cuts";

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
    
    //if (cut_of_interest->IsInside(s800->GetObjE1Correlated_MTDC(0,0),s800->GetXfpE1Correlated_MTDC())){
    //if (cut_of_interest->IsInside(s800->GetMTofObjE1(),s800->GetIonChamber().GetSum())){
    if (cut_of_interest->IsInside(s800->GetMTofObjE1Chn15(),s800->GetIonChamber().GetAve())){
        return true;
    }
    return false;
}
