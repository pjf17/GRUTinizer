#include <fstream>
#include <string>
#include <sstream>

#include "TS800Sim.h"

#include "TGEBEvent.h"

TS800Sim::TS800Sim(){
  Clear();
}

TS800Sim::~TS800Sim() {}

void TS800Sim::Copy(TObject& obj) const {
  TDetector::Copy(obj);

  TS800Sim& s800sim = (TS800Sim&)obj;
  s800sim.s800sim_hits = s800sim_hits; 
}

void TS800Sim::InsertHit(const TDetectorHit& hit){
  s800sim_hits.emplace_back((TS800SimHit&)hit);
}

int TS800Sim::BuildHits(std::vector<TRawEvent>& raw_data){
  for(auto& event : raw_data){
    TGEBEvent& geb = (TGEBEvent&)event;
    BuildFrom(geb);
  }
  return 0;
}

void TS800Sim::Print(Option_t *opt) const { }

void TS800Sim::Clear(Option_t *opt) {
  TDetector::Clear(opt);
  s800sim_hits.clear(); 
}

TVector3 TS800Sim::Track(double sata,double sbta) const {
  //Dividing by 1000 because GRUTinizer uses radians for ATA/BTA
  double ata = TMath::Sin(GetS800SimHit(0).GetATA()+sata);
  double bta = TMath::Sin(GetS800SimHit(0).GetBTA()+sbta);
  TVector3 track(ata,-bta,sqrt(1-ata*ata-bta*bta));
  return track;//.Unit();
}

float TS800Sim::Azita(float ata, float bta) const{
  float xsin = TMath::Sin(ata);//ATA and BTA are in mrad for the simulation!!
  float ysin = TMath::Sin(bta);//ATA and BTA are in mrad for the simulation!!
  float azita = 0.0;
  if(xsin>0 && ysin>0){
    azita = TMath::ATan(ysin/xsin);
  } else if(xsin<0 && ysin>0){
    azita = TMath::Pi()-TMath::ATan(ysin/TMath::Abs(xsin));
  } else if(xsin<0 && ysin<0){
    azita = TMath::Pi()+TMath::ATan(TMath::Abs(ysin)/TMath::Abs(xsin));
  } else if(xsin>0 && ysin<0){
    azita = 2.0*TMath::Pi()-TMath::ATan(TMath::Abs(ysin)/xsin);
  } else{
    azita = 0;
  }
  return azita;
}

float TS800Sim::AdjustedBeta(float beta) const {
  const TS800SimHit &s800simhit = GetS800SimHit(0);
  double gamma = 1.0/(sqrt(1.-beta*beta));
  double dp_p = gamma/(1.+gamma) * s800simhit.GetDTA();
  beta *=(1.+dp_p/(gamma*gamma));
  return beta;
}

void TS800Sim::BuildFrom(TGEBEvent &event){
  const char* data = event.GetPayload();
  //  std::cout << " IN HERE !!!! " << std::endl;
  TRawEvent::G4S800* s800pack = (TRawEvent::G4S800*)data;
  //  event.Print("all");
  //  std::cout << *s800pack;

  SetTimestamp(event.GetTimestamp());
  
  //  std::cout << "TS800Sim::BuildFrom: DTA =" << s800pack->GetDTA() << std::endl;

  TS800SimHit hit;
  hit.fATA = s800pack->GetATA(); 
  hit.fBTA = s800pack->GetBTA(); 
  hit.fDTA = s800pack->GetDTA(); 
  hit.fYTA = s800pack->GetYTA(); 
  s800sim_hits.push_back(hit);
  //exit(1);


}
