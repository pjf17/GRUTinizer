#include "TRuntimeObjects.h"

#include "TObject.h"
#include "TS800.h"

// extern "C" is needed to prevent name mangling.
// The function signature must be exactly as shown here,
//   or else bad things will happen.
extern "C"
void MakeHistograms(TRuntimeObjects& obj) {

  TS800 *s800       = obj.GetDetector<TS800>();
  
  std::string histname = "";
  std::string dirname  = "";

  if(s800){

    double ic_sum = s800->GetIonChamber().GetAve();

    // std::cout << "Dispersive X value: " << s800->GetCrdc(0).GetDispersiveX() << std::endl; 
    // std::cout << "NonDispersive Y value: " << s800->GetCrdc(0).GetNonDispersiveY() << std::endl;

    obj.FillHistogram("s800","CRDC1Y",10000,-5000,5000,s800->GetCrdc(0).GetNonDispersiveY());
    obj.FillHistogram("s800","CRDC2Y",10000,-5000,5000,s800->GetCrdc(1).GetNonDispersiveY());
    obj.FillHistogram("s800","CRDC1X",800,-400,400,s800->GetCrdc(0).GetDispersiveX());
    obj.FillHistogram("s800","CRDC2X",800,-400,400,s800->GetCrdc(1).GetDispersiveX());

    dirname = "MaskCal_gated";
    if(ic_sum > 200){
      histname = "CRDC1_Gated";
      obj.FillHistogram(dirname,histname,
			520,-10,250,s800->GetCrdc(0).GetDispersiveX(),
			500,0,4000,s800->GetCrdc(0).GetNonDispersiveY());
      
      histname = "CRDC2_Gated";
      obj.FillHistogram(dirname,histname,
			520,-10,250,s800->GetCrdc(1).GetDispersiveX(),
			500,0,4000,s800->GetCrdc(1).GetNonDispersiveY());

      histname = "CRDC1_Cal_Gated";
      obj.FillHistogram(dirname,histname,
		        800,-400,400,s800->GetCrdc(0).GetDispersiveX(),
		        3000,-200,200,s800->GetCrdc(0).GetNonDispersiveY());

      histname = "CRDC2_Cal_Gated";
      obj.FillHistogram(dirname,histname,
		        800,-400,400,s800->GetCrdc(1).GetDispersiveX(),
		        3000,-200,200,s800->GetCrdc(1).GetNonDispersiveY());
    }

    histname = "IC_Energy";   
    obj.FillHistogram(histname,1000,-100,4000,ic_sum);               

    dirname = "MaskCal";
    histname = "CRDC1";
    obj.FillHistogram(dirname,histname,
		      520,-10,250,s800->GetCrdc(0).GetDispersiveX(),
		      500,0,4000,s800->GetCrdc(0).GetNonDispersiveY());
      
    histname = "CRDC2";
    obj.FillHistogram(dirname,histname,
		      520,-10,250,s800->GetCrdc(1).GetDispersiveX(),
		      500,0,4000,s800->GetCrdc(1).GetNonDispersiveY());

    histname = "CRDC1_Cal";
    obj.FillHistogram(dirname,histname,
		      520,-400,400,s800->GetCrdc(0).GetDispersiveX(),
		      3000,-200,200,s800->GetCrdc(0).GetNonDispersiveY());

    histname = "CRDC2_Cal";
    obj.FillHistogram(dirname,histname,
		      800,-400,400,s800->GetCrdc(1).GetDispersiveX(),
		      3000,-200,200,s800->GetCrdc(1).GetNonDispersiveY());

    dirname = "GetYOffset";
    histname = "CRDC1_Y_vs_S800Timestamp";
    obj.FillHistogram(dirname,histname,
		      10000,0,6000,s800->GetTimestamp()/1e8,
		      800,-400,400,s800->GetCrdc(0).GetNonDispersiveY());
    
    histname = "CRDC2_Y_vs_S800Timestamp";
    obj.FillHistogram(dirname,histname,
		      10000,0,6000,s800->GetTimestamp()/1e8,
		      800,-400,400,s800->GetCrdc(1).GetNonDispersiveY());
  
    histname = "CRDC1_Y_vs_S800Timestamp_UnCal";
    obj.FillHistogram(dirname,histname,
		      10000,0,6000,s800->GetTimestamp()/1e8,
		      800,-4000,4000,s800->GetCrdc(0).GetNonDispersiveY());
    
    histname = "CRDC2_Y_vs_S800Timestamp_UnCal";
    obj.FillHistogram(dirname,histname,
		      10000,0,6000,s800->GetTimestamp()/1e8,
		      800,-4000,4000,s800->GetCrdc(1).GetNonDispersiveY());
    
      dirname = "InverseMap";
     
      histname = "S800_YTA";
      obj.FillHistogram(dirname,histname,
			1000,-50,50,s800->GetYta());
      
      histname = "S800_DTA";
      obj.FillHistogram(dirname,histname,
			1000,-0.2,0.2,s800->GetDta());
      
      histname = "ATA";
      obj.FillHistogram(dirname,histname,
                        1000,-0.2,0.2,s800->GetAta());
      
      histname = "BTA";
      obj.FillHistogram(dirname,histname,
                        1000,-0.2,0.2,s800->GetBta());
      
      histname = "ATA_vs_BTA";
      obj.FillHistogram(dirname,histname,
			1000,-0.2,0.2,s800->GetAta(),
			1000,-0.2,0.2,s800->GetBta());
	    
    
  }


}
