#include "TRuntimeObjects.h"

#include "TObject.h"
#include "TS800.h"
#include "TOldSega_e05017.h"
#include "TMath.h"

#include "GCutG.h"


std::vector<GCutG*> incoming_cuts   = {0};
std::vector<GCutG*> outgoing_cuts   = {0};
std::vector<GCutG*> timeenergy_cuts; // = {0};
int gates_loaded=0;

bool OutgoingBeam(TRuntimeObjects& obj,GCutG *incoming) {
  TS800    *s800    = obj.GetDetector<TS800>();
  if(!s800)
    return false;

  std::string histname;
  std::string dirname;
  if(incoming)
    dirname = Form("outgoing_%s",incoming->GetName());
  else
    dirname = "outgoing";


  double objtime = s800->GetTof().GetOBJ();
  double xfptime = s800->GetTof().GetXFP();
  if(incoming) {
    if(!incoming->IsInside(objtime,xfptime))
      return false;
  }

  double obj_corr   = s800->GetCorrTOF_OBJ(); 
  double ic_sum     = s800->GetIonChamber().GetAve();
  double afp        = s800->GetAFP();
  double xfp        = s800->GetXFP(0);

  histname = "AFP_vs_OBJTOF";
  obj.FillHistogram(dirname,histname,
  		    1000,-1000,1000,obj_corr,
        	    1000,-0.5,0.5,afp); // check units of AFP
 
  histname = "XFP_vs_OBJTOF";
  obj.FillHistogram(dirname,histname,
        	    1000,-1000,1000,obj_corr,
        	    350,-50,300,xfp);
 
  histname = "IC_vs_OBJTOF_PID";
  obj.FillHistogram(dirname,histname,
        	    1000,-1000,1000,obj_corr,
        	    1000,0,2000,ic_sum);

  obj.FillHistogram(dirname,"CRDC1Y",5000,-5000,5000,s800->GetCrdc(0).GetNonDispersiveY());
  obj.FillHistogram(dirname,"CRDC2Y",5000,-5000,5000,s800->GetCrdc(1).GetNonDispersiveY());
  obj.FillHistogram(dirname,"CRDC1X",400,-400,400,s800->GetCrdc(0).GetDispersiveX());
  obj.FillHistogram(dirname,"CRDC2X",400,-400,400,s800->GetCrdc(1).GetDispersiveX());
  obj.FillHistogram(dirname,"TrigBit",30,0,30,s800->GetTrigger().GetRegistr());
  obj.FillHistogram(dirname,"S800_YTA",1000,-50,50,s800->GetYta());
  obj.FillHistogram(dirname,"S800_DTA",1000,-0.2,0.2,s800->GetDta());
  obj.FillHistogram(dirname,"ATA_vs_BTA",1000,-0.2,0.2,s800->GetAta(),
			                 1000,-0.2,0.2,s800->GetBta());

  return true;
}

bool IncomingBeam(TRuntimeObjects& obj,GCutG *outgoing) {
  TS800    *s800    = obj.GetDetector<TS800>();
  if(!s800)
    return false;

  std::string histname;
  std::string dirname;
  if(outgoing)
    dirname = Form("incoming_%s",outgoing->GetName());
  else
    dirname = "incoming";
  
  double ic_sum      = s800->GetIonChamber().GetAve();
  double obj_corr    = s800->GetCorrTOF_OBJ();
  if(outgoing) {
    if(!outgoing->IsInside(obj_corr,ic_sum))
      return false;
  }

  double objtime = s800->GetTof().GetOBJ();
  double xfptime = s800->GetTof().GetXFP();
  histname = "IncomingPID";
  obj.FillHistogram(dirname,histname,
                     1000,-1000,5000,objtime,
                     1000,-1000,5000,xfptime);
  return true;
}

int HandleUngatedOldSega(TRuntimeObjects& obj) {

  //TS800     *s800 = obj.GetDetector<TS800>();
  TOldSega  *sega = obj.GetDetector<TOldSega>();
  //if(!s800 || !sega)
    //return false;

  std::string dirname;
  std::string histname;

  for(unsigned int i=0;i<sega->Size();i++) {
    TOldSegaHit hit = sega->GetSegaHit(i);

    obj.FillHistogram("Ungated","Det_Charge",
                    4096,0,4096,hit.Charge(),
                    30,0,30,hit.GetDetId());
    obj.FillHistogram("Ungated","Det_Energy",
                    4096,0,4096,hit.GetEnergy(),
                    30,0,30,hit.GetDetId());

    if(hit.Size()>0){
      histname = Form("SegDoppler_Det%02i",hit.GetDetId());
      obj.FillHistogram("Ungated",histname,
                        4096,0,4096,hit.GetDoppler(0.4),
                        40,0,40,hit.GetSegId(0));
    } 
    
    for(int j=0; j<hit.Size(); j++) {
      histname = Form("SegCharge_Det%02i",hit.GetDetId());
      obj.FillHistogram("Ungated",histname,
                        4096,0,4096,hit.GetSegChg(j),
                        40,0,40,hit.GetSegId(j));
    } 
  }

  return 0;
}

int HandleOldSega(TRuntimeObjects& obj,GCutG *incoming,GCutG *outgoing) {

  TS800     *s800 = obj.GetDetector<TS800>();
  TOldSega  *sega = obj.GetDetector<TOldSega>();
  if(!s800 || !sega)
    return false;

  std::string dirname;
  std::string histname;
  std::string histname2;

  dirname = Form("sega_%s",outgoing->GetName());

  obj.FillHistogram(dirname,"CRDC1Y",5000,-5000,5000,s800->GetCrdc(0).GetNonDispersiveY());
  obj.FillHistogram(dirname,"CRDC2Y",5000,-5000,5000,s800->GetCrdc(1).GetNonDispersiveY());

   //this enforces that the outgoing PID blob (AX_from_BY_incoming) is in the incoming gate (BY_incoming)
  const char *outgoing_name = outgoing->GetName();
  const char *toRemove = outgoing->GetTitle(); //going to remove AX from name
  size_t length = strlen(toRemove); //length of AX (can be 23Al or 26P, for example)
  outgoing_name += length+6; //remove AX_from_  
  if(strcmp(outgoing_name,incoming->GetName()) != 0)
    return false;

  double beta = GValue::Value(Form("BETA_%s",outgoing->GetTitle()));
  TVector3 track = s800->Track();
  TVector3 junk(0,0,1);
  //printf("beta = %f\n",beta);
  //track.Print();
  double calorimeter = 0;

  for(unsigned int i=0;i<sega->Size();i++) {
    TOldSegaHit hit = sega->GetSegaHit(i);
    //double z_shift = GValue::Value(Form("Z_SHIFT_DET%i",hit.GetDetId()));
    //hit.GetPosition().Print();
    //printf("Angle no Track: %f\n", (hit.GetPosition()).Angle(junk));
    //printf("Angle with Track: %f\n", (hit.GetPosition()).Angle(track));
    //printf("Cos(Angle) with Track: %f\n", TMath::Cos((hit.GetPosition()).Angle(track)));

    //double corr_time = hit.GetTime() - ((s800->GetTrigger().GetS800Source() + s800->GetTof().GetOBJ())*(0.1/0.25)); 

    //printf("corr_time = %f\n ",corr_time);
    //printf("GetTime = %f\n ",hit.GetTime());
    //printf("S800Source = %f\n ",s800->GetTrigger().GetS800Source());
    //printf("S800TOF = %f\n ",s800->GetTof().GetOBJ());

    if(i == 0){
      calorimeter = 0;
    }
    calorimeter += hit.GetDoppler(beta,&track);
    if(i == (sega->Size() - 1)){
      obj.FillHistogram(dirname,"Doppler_Calorimeter",
                        4096,0,4096,calorimeter);
    }

    if(sega->Size()==1){
      obj.FillHistogram(dirname,"Doppler_Singles",
                        4096,0,4096,hit.GetDoppler(beta,&track));
    }
    if(sega->Size()>1){
      for(unsigned int j=i+1;j<sega->Size();j++) {
        TOldSegaHit hit2 = sega->GetSegaHit(j);
        if(hit2.Charge()<50 || hit.Charge()<50)
          break;
        obj.FillHistogram(dirname,"Doppler_Doppler",
                          1024,0,4096,hit.GetDoppler(beta,&track),
                          1024,0,4096,hit2.GetDoppler(beta,&track));
        obj.FillHistogram(dirname,"Doppler_Doppler",
                          1024,0,4096,hit2.GetDoppler(beta,&track),
                          1024,0,4096,hit.GetDoppler(beta,&track));
        if(sega->Size()==2){
          obj.FillHistogram(dirname,"Doppler_Doppler_Mult2",
                          1024,0,4096,hit.GetDoppler(beta,&track),
                          1024,0,4096,hit2.GetDoppler(beta,&track));
          obj.FillHistogram(dirname,"Doppler_Doppler_Mult2",
                          1024,0,4096,hit2.GetDoppler(beta,&track),
                          1024,0,4096,hit.GetDoppler(beta,&track));
        }
      }
    }

    obj.FillHistogram(dirname,"Energy_Time",
                    //1000,0,30000,corr_time,
                    1000,0,30000,hit.GetTime(),
                    1024,0,4096,hit.GetEnergy());

    obj.FillHistogram(dirname,"Det_Doppler_NoTrack",
                    4096,0,4096,hit.GetDoppler(beta),
                    30,0,30,hit.GetDetId());

    obj.FillHistogram(dirname,"Det_Doppler_Track",
                    4096,0,4096,hit.GetDoppler(beta,&track),
                    30,0,30,hit.GetDetId()); 

    obj.FillHistogram(dirname,Form("Det_Doppler_Track_Theta_det%02i",hit.GetDetId()),
                    180,0,180,hit.GetPosition().Theta()*TMath::RadToDeg(),
                    4096,0,4096,hit.GetDoppler(beta,&track));
    
    obj.FillHistogram(dirname,"Det_Doppler_Track_Phi",
                    360,-180,180,hit.GetPosition().Phi()*TMath::RadToDeg(),
                    4096,0,4096,hit.GetDoppler(beta,&track));

    //TVector3 trackrot = track;
    //obj.FillHistogram(dirname,"Det_Doppler_Track0",
    //                4096,0,4096,hit.GetDoppler(beta,&trackrot),
    //                30,0,30,hit.GetDetId());
    //trackrot.RotateZ(TMath::Pi()/2);
    //obj.FillHistogram(dirname,"Det_Doppler_Track1",
    //                4096,0,4096,hit.GetDoppler(beta,&trackrot),
    //                30,0,30,hit.GetDetId());
    //trackrot.RotateZ(TMath::Pi()/2);
    //obj.FillHistogram(dirname,"Det_Doppler_Track2",
    //                4096,0,4096,hit.GetDoppler(beta,&trackrot),
    //                30,0,30,hit.GetDetId());
    //trackrot.RotateZ(TMath::Pi()/2);
    //obj.FillHistogram(dirname,"Det_Doppler_Track3",
    //                4096,0,4096,hit.GetDoppler(beta,&trackrot),
    //                30,0,30,hit.GetDetId());

    histname = Form("Doppler_ZShift_Track_Det%02i",hit.GetDetId());
    histname2 = Form("Doppler_ZShift_NoTrack_Det%02i",hit.GetDetId());
    for(int z_i=0;z_i<601;z_i++){
      double z_use = -3 + z_i*0.01;
      TVector3 vec(0,0,z_use);
      TVector3 vecx(z_use,0,0); // x shift
      TVector3 vecy(0,z_use,0); // y shift
      obj.FillHistogram(dirname,histname,
                      601,-3,3,z_use,
                      4096,0,4096,hit.GetEnergy()*(1/(sqrt(1-pow(beta,2))))*(1 - beta*TMath::Cos((hit.GetPosition() + vec).Angle(track)))
);
      obj.FillHistogram(dirname,"Doppler_ZShift_Track",
                      601,-3,3,z_use,
                      4096,0,4096,hit.GetEnergy()*(1/(sqrt(1-pow(beta,2))))*(1 - beta*TMath::Cos((hit.GetPosition() + vec).Angle(track)))
);

      obj.FillHistogram(dirname,"Doppler_XShift_Track",
                      601,-3,3,z_use,
                      4096,0,4096,hit.GetEnergy()*(1/(sqrt(1-pow(beta,2))))*(1 - beta*TMath::Cos((hit.GetPosition() + vecx).Angle(track)))
);

      obj.FillHistogram(dirname,"Doppler_YShift_Track",
                      601,-3,3,z_use,
                      4096,0,4096,hit.GetEnergy()*(1/(sqrt(1-pow(beta,2))))*(1 - beta*TMath::Cos((hit.GetPosition() + vecy).Angle(track)))
);

     if(hit.GetDetId()<11){
       obj.FillHistogram(dirname,"Doppler_ZShift_Track_37",
                      601,-3,3,z_use,
                      4096,0,4096,hit.GetEnergy()*(1/(sqrt(1-pow(beta,2))))*(1 - beta*TMath::Cos((hit.GetPosition() + vec).Angle(track)))
);

       obj.FillHistogram(dirname,"Doppler_XShift_Track_37",
                      601,-3,3,z_use,
                      4096,0,4096,hit.GetEnergy()*(1/(sqrt(1-pow(beta,2))))*(1 - beta*TMath::Cos((hit.GetPosition() + vecx).Angle(track)))
);
  
       obj.FillHistogram(dirname,"Doppler_YShift_Track_37",
                      601,-3,3,z_use,
                      4096,0,4096,hit.GetEnergy()*(1/(sqrt(1-pow(beta,2))))*(1 - beta*TMath::Cos((hit.GetPosition() + vecy).Angle(track)))
);
     }
     if(hit.GetDetId()>=11){
       obj.FillHistogram(dirname,"Doppler_ZShift_Track_90",
                      601,-3,3,z_use,
                      4096,0,4096,hit.GetEnergy()*(1/(sqrt(1-pow(beta,2))))*(1 - beta*TMath::Cos((hit.GetPosition() + vec).Angle(track)))
);

       obj.FillHistogram(dirname,"Doppler_XShift_Track_90",
                      601,-3,3,z_use,
                      4096,0,4096,hit.GetEnergy()*(1/(sqrt(1-pow(beta,2))))*(1 - beta*TMath::Cos((hit.GetPosition() + vecx).Angle(track)))
);

       obj.FillHistogram(dirname,"Doppler_YShift_Track_90",
                      601,-3,3,z_use,
                      4096,0,4096,hit.GetEnergy()*(1/(sqrt(1-pow(beta,2))))*(1 - beta*TMath::Cos((hit.GetPosition() + vecy).Angle(track)))
);
     }

      obj.FillHistogram(dirname,histname2,
                      601,-3,3,z_use,
                      4096,0,4096,hit.GetEnergy()*(1/(sqrt(1-pow(beta,2))))*(1 - beta*TMath::Cos((hit.GetPosition() + vec).Angle(junk)))
);
      obj.FillHistogram(dirname,"Doppler_ZShift_NoTrack",
                      601,-3,3,z_use,
                      4096,0,4096,hit.GetEnergy()*(1/(sqrt(1-pow(beta,2))))*(1 - beta*TMath::Cos((hit.GetPosition() + vec).Angle(junk)))
);
    }
    
    //histname = Form("Doppler_ThetaShift_NoTrack_Det%02i",hit.GetDetId());
    //histname2 = Form("Doppler_ThetaShift_Track_Det%02i",hit.GetDetId());
    //for(int t_i=0;t_i<501;t_i++){
    //  double t_use = -0.1 + t_i*0.0004;
    //   obj.FillHistogram(dirname,histname,
    //                  501,-0.1,0.1,t_use,
    //                  4096,0,4096,hit.GetEnergy()*(1/(sqrt(1-pow(beta,2))))*(1 - beta*TMath::Cos((hit.GetPosition()).Angle(junk) + t_use))
//);
    //  obj.FillHistogram(dirname,"Doppler_ThetaShift_NoTrack",
    //                  501,-0.1,0.1,t_use,
    //                  4096,0,4096,hit.GetEnergy()*(1/(sqrt(1-pow(beta,2))))*(1 - beta*TMath::Cos((hit.GetPosition()).Angle(junk) + t_use))
//);
    //  obj.FillHistogram(dirname,histname2,
    //                  501,-0.1,0.1,t_use,
    //                  4096,0,4096,hit.GetEnergy()*(1/(sqrt(1-pow(beta,2))))*(1 - beta*TMath::Cos((hit.GetPosition()).Angle(track) + t_use))
//);
    //  obj.FillHistogram(dirname,"Doppler_ThetaShift_Track",
    //                  501,-0.1,0.1,t_use,
    //                  4096,0,4096,hit.GetEnergy()*(1/(sqrt(1-pow(beta,2))))*(1 - beta*TMath::Cos((hit.GetPosition()).Angle(track) + t_use))
//);
    //  if(hit.GetDetId()<11){  
    //  obj.FillHistogram(dirname,"Doppler_ThetaShift_Track_37",
    //                  501,-0.1,0.1,t_use,
    //                  4096,0,4096,hit.GetEnergy()*(1/(sqrt(1-pow(beta,2))))*(1 - beta*TMath::Cos((hit.GetPosition()).Angle(track) + t_use))
//);
    //  }
    //  if(hit.GetDetId()>=11){  
    //  obj.FillHistogram(dirname,"Doppler_ThetaShift_Track_90",
    //                  501,-0.1,0.1,t_use,
    //                  4096,0,4096,hit.GetEnergy()*(1/(sqrt(1-pow(beta,2))))*(1 - beta*TMath::Cos((hit.GetPosition()).Angle(track) + t_use))
//);
    //  }
    //}

    obj.FillHistogram(dirname,"Det_Energy",
                    4096,0,4096,hit.GetEnergy(),
                    30,0,30,hit.GetDetId());

    for(int beta_i=0;beta_i<500;beta_i++){
      double beta_use = 0.001*double(beta_i);
      histname = "FindBeta_Doppler";
      obj.FillHistogram(dirname,histname,
		        500,0,0.499,beta_use,
			4096,0,4096,hit.GetDoppler(beta_use,&track));
      if(hit.GetDetId()<11){
        histname = "FindBeta_Doppler_37";
        obj.FillHistogram(dirname,histname,
		          500,0,0.499,beta_use,
			  4096,0,4096,hit.GetDoppler(beta_use,&track));
      }
      if(hit.GetDetId()>=11){
        histname = "FindBeta_Doppler_90";
        obj.FillHistogram(dirname,histname,
		          500,0,0.499,beta_use,
			  4096,0,4096,hit.GetDoppler(beta_use,&track));
      } 
    }

    //for(int j=0; j<hit.Size(); j++) {
    if(hit.Size()>0){
      histname = Form("SegDoppler_Det%02i",hit.GetDetId());
      obj.FillHistogram(dirname,histname,
                        4096,0,4096,hit.GetDoppler(beta,&track),
                        40,0,40,hit.GetSegId(0));
      histname = Form("SegEnergy_Det%02i",hit.GetDetId());
      obj.FillHistogram(dirname,histname,
                        4096,0,4096,hit.GetEnergy(),
                        40,0,40,hit.GetSegId(0));
    } 
  }  

  return 0;
}

// extern "C" is needed to prevent name mangling.
// The function signature must be exactly as shown here,
//   or else bad things will happen.
extern "C"
void MakeHistograms(TRuntimeObjects& obj) {
 
  TList *list  = &(obj.GetObjects());
  TList *gates = &(obj.GetGates());
  int numobj = list->GetSize();

  HandleUngatedOldSega(obj); 

  if(gates_loaded!=gates->GetSize()) {
    TIter iter(gates);
    while(TObject *obj = iter.Next()) {
      GCutG *gate = (GCutG*)obj;
      std::string tag = gate->GetTag();
      if(!tag.compare("incoming")) {
        incoming_cuts.push_back(gate);
      } else if(!tag.compare("outgoing")) {
        outgoing_cuts.push_back(gate);
      } else if(!tag.compare("timeenergy")){
        timeenergy_cuts.push_back(gate);
      }
      gates_loaded++;
    }
  }

  //printf("incoming.size() == %i\n",incoming_cuts.size());
  int incoming_passed=-1;
  int outgoing_passed=-1;
  for(unsigned int x=0;x<incoming_cuts.size();x++) {
    bool passed = OutgoingBeam(obj,incoming_cuts.at(x)); 
    if(x!=0 && passed) {
      incoming_passed = x;
      break;
    }
  }
  for(unsigned int x=0;x<outgoing_cuts.size();x++) {
    bool passed = IncomingBeam(obj,outgoing_cuts.at(x)); 
    if(x!=0 && passed) {
      outgoing_passed = x;
      break;
    }
  } 


  if(incoming_passed>0 && outgoing_passed>0) {
    HandleOldSega(obj,incoming_cuts.at(incoming_passed),outgoing_cuts.at(outgoing_passed));
  }



  if(numobj!=list->GetSize())
    list->Sort();

}//end MakeHistograms
