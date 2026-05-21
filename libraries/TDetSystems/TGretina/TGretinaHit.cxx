#include "TGretinaHit.h"

#include <algorithm>
#include <ctime>
#include <cmath>
#include <set>

#include <TRandom.h>

#include "GValue.h"
#include "TGEBEvent.h"
#include "TGretina.h"
#include "TS800.h"

double lastPointPenalty(double x){
  double val = TMath::TanH((x/100-1.1)/0.5)*exp(-3*(x/100-1.1)+1);
  if (val < 0) val = 0;
  return val + 1;
}

TGretinaHit::TGretinaHit(){ Clear(); }

TGretinaHit::~TGretinaHit(){ }

void TGretinaHit::Copy(TObject &rhs) const {
  TDetectorHit::Copy(rhs);
  ((TGretinaHit&)rhs).fWalkCorrection = fWalkCorrection;
  ((TGretinaHit&)rhs).fPad            = fPad;
  ((TGretinaHit&)rhs).fTOffset           = fTOffset;
  ((TGretinaHit&)rhs).fCrystalId      = fCrystalId;
  ((TGretinaHit&)rhs).fCoreEnergy     = fCoreEnergy;
  ((TGretinaHit&)rhs).fCoreCharge[0]  = fCoreCharge[0];
  ((TGretinaHit&)rhs).fCoreCharge[1]  = fCoreCharge[1];
  ((TGretinaHit&)rhs).fCoreCharge[2]  = fCoreCharge[2];
  ((TGretinaHit&)rhs).fCoreCharge[3]  = fCoreCharge[3];
  ((TGretinaHit&)rhs).fNumberOfInteractions = fNumberOfInteractions;
  ((TGretinaHit&)rhs).fSegments       = fSegments;
  ((TGretinaHit&)rhs).fSingles      = fSingles;
  ((TGretinaHit&)rhs).fDecompChi2  = fDecompChi2;
  ((TGretinaHit&)rhs).fDecompNormChi2  = fDecompNormChi2;
  ((TGretinaHit&)rhs).fAB  = fAB;
}

Float_t  TGretinaHit::GetCoreEnergy() const {
  double cal_offst = GValue::Value(Form("GRETINA_%d_ECALIB_a0",fCrystalId));
  double cal_slope = GValue::Value(Form("GRETINA_%d_ECALIB_a1",fCrystalId));
  if (std::isnan(cal_offst) || std::isnan(cal_slope)){
    cal_offst = 0;
    cal_slope = 1;
  }

  return fCoreEnergy*cal_slope + cal_offst;
}

Float_t TGretinaHit::GetCoreEnergy(int i) const {
  float charge = (float)GetCoreCharge(i) + gRandom->Uniform();
  //board_id=; //card  number : 0x0030  information not available here.
  TChannel *channel = TChannel::GetChannel(Address()+(i<<4));
  //printf("GetAddress() + i = 0x%08x\n",GetAddress()+i); 
  if(!channel)
    return charge;
  return channel->CalEnergy(charge);
}

const char *TGretinaHit::GetName() const {
  TChannel *channel = TChannel::GetChannel(Address());
  if(channel) 
    return channel->GetName();
  std::string name = "";
  return name.c_str();
}

Int_t TGretinaHit::GetRingNumber() const {
  return TGretina::GetRingNumber(fCrystalId);
}

TVector3 TGretinaHit::GetCrystalPosition()  const { 
  return TGretina::GetCrystalPosition(fCrystalId); 
}

TGretinaHit TGretinaHit::GetNeighbor(int i) const {
  i+=1; //ensure we never return the original hit
  if (i < (int) fSingles.size() && i > 0){
    return fSingles[i];
  } else {
    std::cout<<"ERROR RETURNING NEIGHBOR HIT"<<std::endl;
    return fSingles[0];
  }
}

TGretinaHit TGretinaHit::GetInitialHit() const {
  TGretinaHit output = fSingles[0];
  output.fSingles.clear();
  output.fSetFirstSingles = false;
  return output;
}

//TVector3 TGretinaHit::GetSegmentPosition(int i)  const { 
//  if(fSegmentNumber.size()) 
//    return TGretina::GetSegmentPosition(fCrystalId,fSegmentNumber.at(0));
//  else
//    return TGretina::GetSegmentPosition(fCrystalId,0);
//}



void TGretinaHit::BuildFrom(TSmartBuffer& buf){
  const TRawEvent::GEBBankType1& raw = *(const TRawEvent::GEBBankType1*)buf.GetData();

  //std::cout << "GretinaHit: " << raw << std::endl;
  //std::cout << raw << std::endl;
  Clear();

  SetTimestamp(raw.timestamp);
  fWalkCorrection = raw.t0;
  fCrystalId = raw.crystal_id;
  fCoreEnergy = raw.tot_e;

  //fAddress = (1<<24) + (fCrystalId<<16);
  //fAddress = (1<<24) + ( raw.board_id );
  int board_id = ((fCrystalId/4) << 8) ;  //hole  number : 0x1f00
  //    board_id =                       ;  //card  number : 0x0030  information not available here.
  board_id += ((fCrystalId%4) << 6) ;  //x-tal number : 0x00c0
  board_id += 9;                       //chan  number : 0x000f  information not available here(assume core).
  SetAddress((1<<24) + board_id);


  for(int i=0; i<4; i++){
    fCoreCharge[i] = raw.core_e[i];
  }

  fNumberOfInteractions = raw.num;
  fPad = raw.pad;
  fDecompChi2 = raw.chisq;
  fDecompNormChi2 = raw.norm_chisq;

  //std::cout<< "pad\t" << fPad << "\tints\t" << fNumberOfInteractions << std::endl; 
  for(int i=0; i<fNumberOfInteractions; i++) {
    try {
      interaction_point pnt; 
      pnt.fSeg = raw.intpts[i].seg;
      pnt.fX   = raw.intpts[i].x;
      pnt.fY   = raw.intpts[i].y;
      pnt.fZ   = raw.intpts[i].z;
      pnt.fEng = raw.intpts[i].e;
      fSegments.push_back(pnt);
    } catch(...) {
      std::cout << "in try catch block!" << std::endl;
    }

  }
  fTOffset = raw.intpts[MAX_INTPTS-1].z;
  //std::cout << "[13].z :  " << raw.intpts[13].z << std::endl;
  //std::cout << "[14].z :  " << raw.intpts[14].z << std::endl;
  //std::cout << "[15].z :  " << raw.intpts[15].z << std::endl;
  //std::cout << "fTOffset :  " << fTOffset             << std::endl;
  
  if (std::isnan(GValue::Value("NO_E_SORT"))) std::sort(fSegments.begin(),fSegments.end());
  // std::reverse(fSegments.begin(),fSegments.end());
  //  Print("all");
  // }
}

TVector3 TGretinaHit::GetIntPosition(unsigned int i) const {
  if(i<fSegments.size()){
    double xoffset = GValue::Value("GRETINA_X_OFFSET");
    if(std::isnan(xoffset))
      xoffset=0.00;
    double yoffset = GValue::Value("GRETINA_Y_OFFSET");
    if(std::isnan(yoffset))
      yoffset=0.00;
    double zoffset = GValue::Value("GRETINA_Z_OFFSET");
    if(std::isnan(zoffset))
      zoffset=0.00;
    TVector3 v = TGretina::CrystalToGlobal(fCrystalId,fSegments.at(i).fX + xoffset,
        fSegments.at(i).fY + yoffset,
        fSegments.at(i).fZ + zoffset);
    return v; 
  } else {
    return TDetectorHit::BeamUnitVec;
  }
}

void TGretinaHit::SetPosition(unsigned int i, double x, double y, double z) {
  fSegments.at(i).fX = x;
  fSegments.at(i).fY = y;
  fSegments.at(i).fZ = z;
}

TVector3 TGretinaHit::GetLocalPosition(unsigned int i) const {
  if(i<fSegments.size()){
    return TVector3(fSegments.at(i).fX,
        fSegments.at(i).fY,
        fSegments.at(i).fZ);
  } else {
    return TVector3(0,0,1);
  }
}


double TGretinaHit::GetDopplerYta(double beta, double yta, const TVector3 *vec, int point, int EngRange) const {
  if(Size()<1)
    return 0.0;
  if(vec==0) {
    vec = &BeamUnitVec;
  }
  if (point == -1 || point > fNumberOfInteractions) point = 0;
  //Target offsets determine new reference point in lab frame
  double xoffset = GValue::Value("TARGET_X_OFFSET");
  if(std::isnan(xoffset))
    xoffset=0.00;
  double yoffset = GValue::Value("TARGET_Y_OFFSET");
  if(std::isnan(yoffset))
    yoffset=0.00;
  double zoffset = GValue::Value("TARGET_Z_OFFSET");
  if(std::isnan(zoffset))
      zoffset=0.00; 
  return GetDopplerYta(beta, yta, xoffset, yoffset, zoffset, vec, point, EngRange);
}

double TGretinaHit::GetDopplerYta(double beta, double yta, double xoffset, double yoffset, double zoffset, const TVector3 *vec, int point, int EngRange) const {
  if(Size()<1)
    return 0.0;
  if(vec==0) {
    vec = &BeamUnitVec;
  }
  double tmp = 0.0;
  double gamma = 1./(sqrt(1.-pow(beta,2.)));
  TVector3 gret_pos = GetIntPosition(point);

  //Target offsets determine new reference point in lab frame
  gret_pos.SetX(gret_pos.X() - xoffset);
  gret_pos.SetY(gret_pos.Y() - (yoffset - yta));
  gret_pos.SetZ(gret_pos.Z() - zoffset);

  double cal_offst = GValue::Value(Form("GRETINA_%d_ECALIB_a0",fCrystalId));
  double cal_slope = GValue::Value(Form("GRETINA_%d_ECALIB_a1",fCrystalId));
  if (std::isnan(cal_offst) || std::isnan(cal_slope)){
    cal_offst = 0;
    cal_slope = 1;
  }

  if(EngRange>0) 
    tmp = GetCoreEnergy(EngRange)*gamma *(1 - beta*TMath::Cos(gret_pos.Angle(*vec)));
  else
    tmp = (fCoreEnergy*cal_slope + cal_offst)*gamma *(1 - beta*TMath::Cos(gret_pos.Angle(*vec)));
  return tmp;
}
/*
   double TGretinaHit::GetDoppler_dB(double beta, const TVector3 *vec,double Dta){
   if(Size()<1)
   return 0.0;
   if(vec==0) {
   vec = &BeamUnitVec;
   }
   double tmp = 0.0;
   double gamma = 1.0/(sqrt(1.0-pow(beta,2.0)));
// Do beta correction here.
double dp_p = gamma/(1.0+gamma)*Dta;
beta *= (1.0+dp_p/(gamma*gamma));
double TheGamma = 1.0/TMath::Sqrt(1.0-beta*beta);
tmp = fCoreEnergy*TheGamma *(1.0 - beta*TMath::Cos(GetPosition().Angle(*vec)));
return tmp;
}
*/


/*
   double TGretinaHit::GetDoppler(const TS800 *s800,bool doDTAcorr,int EngRange) {
   if(!s800 || Size()<1)
   return 0.0;
   double beta  = GValue::Value("BETA");
   if(std::isnan(beta))
   return 0.0;
   double gata =  GValue::Value("GRETINA_ATA_OFFSET");
   if(std::isnan(gata))
   gata = 0.0;
   else 
   gata = gata*TMath::DegToRad();
   double gbta =  GValue::Value("GRETINA_BTA_OFFSET");
   if(std::isnan(gbta))
   gbta = 0.0;
   else 
   gbta = gata*TMath::DegToRad();
   if(doDTAcorr){
   double gamma = 1.0/(sqrt(1.-beta*beta));
   double dp_p = gamma/(1.+gamma) * s800->GetDta();
   beta *=(1.+dp_p/(gamma*gamma));
   }
   TVector3 track = s800->Track(gata,gbta);  //(TMath::Sin(s800->GetAta()),-TMath::Sin(s800->GetBta()),1);
   if(EngRange>-1)
   return GetDoppler(EngRange,beta,&track);
   return GetDoppler(beta,&track);
   }
   */


//TVector3 TGretinaHit::GetCrystalPosition(int i) const {
//  std::cerr << __PRETTY_FUNCTION__ << " NOT IMPLEMENTED YET" << std::endl;
//}

/*
   void TGretinaHit::SortHits(){
// sets are sorted, so this will sort all properties together.
//
// !! When multiple interactions are assigned to a single
//    segment, the first interaction is currenlty assigned
//    to that segment!
//

std::vector<interaction_point> ips;
for(int i=0; i<fNumberOfInteractions; i++){
ips.push_back(interaction_point(fSegmentNumber[i],
fGlobalInteractionPosition[i],
fLocalInteractionPosition[i],
fInteractionEnergy[i],
fInteractionFraction[i]));
}
//printf("ips.size() == %i\n",ips.size());
// Fill all interaction points
//
fSegmentNumber.clear();
fGlobalInteractionPosition.clear();
fLocalInteractionPosition.clear();
fInteractionEnergy.clear();
fInteractionFraction.clear();
//
fNumberOfInteractions = 0;

std::sort(ips.begin(), ips.end());

for(auto& point : ips){
if(fNumberOfInteractions >= MAXHPGESEGMENTS){
break;
}
fSegmentNumber.push_back(point.segnum);
fGlobalInteractionPosition.push_back(point.pos);
fLocalInteractionPosition.push_back(point.local_pos);
fInteractionEnergy.push_back(point.energy);
fInteractionFraction.push_back(point.energy_fraction);
fNumberOfInteractions++;
}

//Print("all");
// Because they are now sorted
fFirstInteraction = 0;
if(fNumberOfInteractions>1)
fSecondInteraction = 1;
}
*/

// TODO: Handle interactions points better
//       Right now, the "first interaction point" is the one with the highest energy,
//       and the "second" is the one with the second highest energy.
//       First and second may be assigned across crystal boundaries.
void TGretinaHit::Add(const TGretinaHit& rhs) {

  // qStash all interaction points
  std::set<interaction_point> ips;
  // Copy other information to self if needed
  double my_core_energy = fCoreEnergy;
  if(fCoreEnergy < rhs.fCoreEnergy) {
    for(unsigned int i=0; i<rhs.fSegments.size(); i++){
      ips.insert(rhs.fSegments[i]);
    }

    for(unsigned int i=0; i<fSegments.size(); i++){
      ips.insert(fSegments[i]);
    }
    rhs.Copy(*this);
    fCoreEnergy += my_core_energy;
  } else {
    for(unsigned int i=0; i<fSegments.size(); i++){
      ips.insert(fSegments[i]);
    }
    for(unsigned int i=0; i<rhs.fSegments.size(); i++){
      ips.insert(rhs.fSegments[i]);
    }
    fCoreEnergy += rhs.fCoreEnergy;
  }

  // Fill all interaction points
  fNumberOfInteractions = 0;
  fSegments.clear();
  for(auto& point : ips){
    if(fNumberOfInteractions >= MAXHPGESEGMENTS){
      break;
    }
    fSegments.push_back(point);
    fNumberOfInteractions++;
  }
}

void TGretinaHit::NNAdd(const TGretinaHit& rhs) {
  //copy original hit into singles
  if (!fSetFirstSingles){
    TGretinaHit myCopy(*this);
    fSingles.push_back(myCopy);
    fSetFirstSingles = true;
  }

  double cal_slope = GValue::Value(Form("GRETINA_%d_ECALIB_a1",fCrystalId));
  if (std::isnan(cal_slope)) fCoreEnergy += rhs.GetCoreEnergy();
  else fCoreEnergy += rhs.GetCoreEnergy()/cal_slope; 
  
  fSingles.push_back(rhs);

  // Fill all interaction points
  for(unsigned int i=0; i<rhs.fSegments.size(); i++){
    if(fNumberOfInteractions >= MAXHPGESEGMENTS){
      break;
    }
    fSegments.push_back(rhs.fSegments[i]);
    fNumberOfInteractions++;
  }
}

void TGretinaHit::NNSwap(int index){
  TVector3 pos = GetNeighbor(index).GetPosition();
  this->SetPosition(0,pos.X(),pos.Y(),pos.Z());
  return;
}
/*
   TGretinaHit& TGretinaHit::operator+=(const TGretinaHit& rhs) {
   AddToSelf(rhs);
   return *this;
   }

   TGretinaHiti& TGretinaHit::operator+(TGretinaHit lhs,const TGretinaHit& rhs) {

   lhs += rhs;
   return lhs;
   }
   */


TVector3 TGretinaHit::GetLastPosition() const {
  if(fSegments.size()<1)
    return TDetectorHit::BeamUnitVec;
  return GetIntPosition(fSegments.size()-1);
}


void TGretinaHit::Print(Option_t *opt) const {

  std::cout << "TGretinaHit:" <<  std::endl;
  //std::cout << "\tAddress:        \t0x" << std::hex << fAddress << std::dec << std::endl;
  printf("\tAddress:        \t0x%08x\n",Address());
  std::cout << "\tHole:           \t" << GetHoleNumber()               << std::endl;
  std::cout << "\tCrystalNum      \t" << GetCrystalNumber()            << std::endl;
  std::cout << "\tCrystalId:      \t" << GetCrystalId()                << std::endl;

  std::cout << "\tLocal Timestamp:\t" << Timestamp()                   << std::endl;
  std::cout << "\tCorrected time: \t" << GetTime()                     << std::endl;
  std::cout << "\tDecomp Energy:  \t" << fCoreEnergy                   << std::endl;
  if(!strcmp(opt,"all")) {
    for(int i=0;i<4;i++)
      std::cout << "\tCharge[" << i << "]:       \t" << fCoreCharge[i] << std::endl;
  }
  std::cout << "\tErrorCode:       \t" << fPad                         << std::endl;
  std::cout << "\tInteractions:    \t" << fSegments.size()             << std::endl;

  if(!strcmp(opt,"all")) {
    for(int i=0;i<fNumberOfInteractions;i++) {
      printf("\t\t");
      fSegments.at(i).Print();
    }
  }
  std::cout << "------------------------------"  << std::endl;

}

void TGretinaHit::Clear(Option_t *opt) {
  TDetectorHit::Clear(opt);

  fWalkCorrection = sqrt(-1);

  fCrystalId      = -1;
  fCoreEnergy     = sqrt(-1);
  fCoreCharge[0]  = -1;
  fCoreCharge[1]  = -1;
  fCoreCharge[2]  = -1;
  fCoreCharge[3]  = -1;

  fTOffset = sqrt(-1);
  fPad  = 0;

  fNumberOfInteractions = 0;

  fSegments.clear();
}

Int_t TGretinaHit::GetNSegments() const {
  int segs[fNumberOfInteractions];
  for (int i=0; i< fNumberOfInteractions; i++){
    segs[i] = fSegments[i].fSeg;
  }
  size_t len = sizeof(segs) / sizeof(segs[0]);
  std::unordered_set<int> s(segs,segs+len);
  return (Int_t) s.size();
}

void TGretinaHit::TrimSegments(int type) {
  // 0: drop multiple ident int pnts.  
  // 1: make into wedge "data"
  if(type==0) {
    std::set<interaction_point,intpnt_compare> pset;
    for(auto x=fSegments.begin();x!=fSegments.end();x++) {
      pset.insert(*x);
    }
    fSegments.clear();
    for(auto x=pset.begin();x!=pset.end();x++) {
      fSegments.push_back(*x);
    }
    // std::sort(fSegments.begin(),fSegments.end());
    fNumberOfInteractions = fSegments.size();
  } else if (type==1) {
    std::set<interaction_point,intpnt_compare_wedge> pset;
    for(auto x=fSegments.begin();x!=fSegments.end();x++) {
      pset.insert(*x);
    }
    fSegments.clear();
    for(auto x=pset.begin();x!=pset.end();x++) {
      fSegments.push_back(*x);
      fSegments.back().fSeg = (fSegments.back().fSeg%6);
    }
    // std::sort(fSegments.begin(),fSegments.end());
    fNumberOfInteractions = fSegments.size();
  }
}

void TGretinaHit::ScaleIntEng(){
  if (fNumberOfInteractions == 1) return;
  double sum = 0;
  for (int i=0; i < fNumberOfInteractions; i++){
    sum += fSegments[i].fEng;
  }
  for (int i=0; i < fNumberOfInteractions; i++){
    fSegments[i].fEng *= fCoreEnergy/sum;
  }
}

double TGretinaHit::GetAlpha(int p1, int p2) const {
  if (fNumberOfInteractions > 1 && p1 < fNumberOfInteractions && p2 < fNumberOfInteractions) {
    return GetIntPosition(p1).Angle(GetIntPosition(p2));
  }
  else return -100;
}

double TGretinaHit::GetScatterAngle(int p1, int p2) const {
  if (fNumberOfInteractions > 1 && p1 < fNumberOfInteractions && p2 < fNumberOfInteractions) {
    TVector3 diff = GetIntPosition(p2) - GetIntPosition(p1);
    return GetIntPosition(p1).Angle(diff);
  }
  else return -100;
}

double TGretinaHit::GetScatterCosine(int p1, int p2) const {
  if (fNumberOfInteractions > 1 && p1 < fNumberOfInteractions && p2 < fNumberOfInteractions) {
    TVector3 v1 = GetIntPosition(p1);
    TVector3 diff = GetIntPosition(p2) - v1;
    return v1.Dot(diff)/diff.Mag()/v1.Mag();
  }
  else return -2;
}

double TGretinaHit::GetXiSimple(const TVector3 *beam, int p1, int p2) const{
  if (fNumberOfInteractions == 1 || p1 >= fNumberOfInteractions || p2 >= fNumberOfInteractions) 
    return -1;

  if (!beam) beam = new TVector3(0,0,1);

  //get interaction points
  TVector3 interaction1 = GetIntPosition(p1); 
  TVector3 interaction2 = GetIntPosition(p2);

  //get angle between plane norms
  TVector3 comptonPlaneNorm = interaction1.Cross(interaction2);
  TVector3 reactionPlaneNorm = beam->Cross(interaction1);
  double xi = reactionPlaneNorm.Angle(comptonPlaneNorm);

  //domain only matters from 0 to 90
  if (xi > TMath::Pi()/2) xi = TMath::Pi() - xi;
  
  return xi;
}

double TGretinaHit::GetXi(const TVector3 *beam, bool extend, int p1, int p2) const{
  if (fNumberOfInteractions > 1 && p1 < fNumberOfInteractions && p2 < fNumberOfInteractions) {
    if (!beam) beam = new TVector3(0,0,1);

    //get interaction points
    TVector3 interaction1 = GetIntPosition(p1); 
    TVector3 interaction2 = GetIntPosition(p2);

    //calculate xi
    TVector3 comptonPlaneNorm = interaction1.Cross(interaction2);
    TVector3 reactionPlaneNorm = beam->Cross(interaction1);
    double xi = reactionPlaneNorm.Angle(comptonPlaneNorm);
    
    if (extend){
      TVector3 basisNorm = reactionPlaneNorm.Cross(interaction1);
      if (basisNorm.Angle(comptonPlaneNorm) > TMath::PiOver2()) xi = TMath::TwoPi() - xi;
    }

    return xi;
  }
}

double TGretinaHit::GetXi(double rpLo, double rpHi, double cpLo, double cpHi, bool &gateCondition, const TVector3 *beam, int p1, int p2) const{
  if (fNumberOfInteractions > 1 && p1 < fNumberOfInteractions && p2 < fNumberOfInteractions) {
    if (!beam) beam = new TVector3(0,0,1);

    //calculate the phase of the crystal wrt the reaction plane
    /*
    TVector3 pos = TGretina::CrystalToGlobal(fCrystalId,0,0,0);
    TVector3 ref = (TVector3(0,0,1)).Cross(pos);
    TVector3 xax = TGretina::CrystalToGlobal(fCrystalId,1,0,0) - pos;
    TVector3 yax = TGretina::CrystalToGlobal(fCrystalId,0,1,0) - pos;
    double phase = ref.Angle(xax);
    if (ref.Angle(yax) > TMath::PiOver2()) phase = TMath::TwoPi() - phase;
    */

    //get interaction points and rotate
    TVector3 interaction1 = GetIntPosition(p1); 
    TVector3 interaction2 = GetIntPosition(p2);

    gateCondition = true;
    if (GetScatterAngle() < cpLo*TMath::DegToRad() || GetScatterAngle() > cpHi*TMath::DegToRad()) gateCondition = false;
    if (beam->Angle(interaction1) < rpLo*TMath::DegToRad() || beam->Angle(interaction1) > rpHi*TMath::DegToRad()) gateCondition = false;

    // if (beam->Angle(interaction1)*TMath::RadToDeg() < 5 || beam->Angle(interaction1)*TMath::RadToDeg() > 175) return -10; 

    // interaction1 = interaction1 + pos;
    // interaction2 = interaction2 + pos;

    //calculate xi
    TVector3 comptonPlaneNorm = interaction1.Cross(interaction2);
    TVector3 reactionPlaneNorm = beam->Cross(interaction1);
    TVector3 basisNorm = reactionPlaneNorm.Cross(interaction1);

    double xi = reactionPlaneNorm.Angle(comptonPlaneNorm);
    // if (basisNorm.Angle(comptonPlaneNorm) > TMath::PiOver2()) xi = TMath::TwoPi() - xi;
    // xi -= phase;
    // if (fCrystalId%2) xi -= 120*TMath::DegToRad();
    // while (xi < 0) xi += TMath::TwoPi();

    return xi;
  }
  // if (fNumberOfInteractions > 1 && p1 < fNumberOfInteractions && p2 < fNumberOfInteractions) {

  //   TVector3 tr = TGretina::CrystalToGlobal(fCrystalId,0,0,0);
  //   TVector3 vX = TGretina::CrystalToGlobal(fCrystalId,1,0,0)-tr;
  //   TVector3 vY = TGretina::CrystalToGlobal(fCrystalId,0,1,0)-tr;
  //   TVector3 vZ = TGretina::CrystalToGlobal(fCrystalId,0,0,1)-tr;
  //   double mtxData[9] = {vX.Unit().X(),vX.Unit().Y(),vX.Unit().Z(),vY.Unit().X(),vY.Unit().Y(),vY.Unit().Z(),vZ.Unit().X(),vZ.Unit().Y(),vZ.Unit().Z()};
  //   TMatrixT <double> mtx = TMatrixT<double> (3,3,mtxData,"F");
  //   mtx.Invert();

  //   TVector3 b = mtx*TVector3(0,0,1);
  //   TVector3 origin = mtx*tr;

  //   TVector3 diff = GetLocalPosition(p2) - GetLocalPosition(p1);
  //   TVector3 gamma = GetLocalPosition(p1) + origin;
  //   TVector3 comptonPlaneNorm = gamma.Cross(diff);
  //   TVector3 reactionPlaneNorm = b.Cross(gamma);
  //   // TVector3 xax = TVector3(1,0,0);
  //   double xi = reactionPlaneNorm.Angle(comptonPlaneNorm);

  //   return xi;
  // }
  else return -10;
}

double TGretinaHit::GetXiChris(const TVector3 *beam, int p1, int p2) const{
  if (fNumberOfInteractions > 1 && p1 < fNumberOfInteractions && p2 < fNumberOfInteractions) {
    if (!beam) beam = new TVector3(0,0,1);

    TVector3 vZ = GetIntPosition(p1);
    TVector3 vY = vZ.Cross(*beam);
    TVector3 vX = vY.Cross(vZ);
    TVector3 r2 = vZ - GetIntPosition(p2);

    double mtxData[9] = {vX.Unit().X(),vX.Unit().Y(),vX.Unit().Z(),vY.Unit().X(),vY.Unit().Y(),vY.Unit().Z(),vZ.Unit().X(),vZ.Unit().Y(),vZ.Unit().Z()};
    TMatrixT <double> mtx = TMatrixT<double> (3,3,mtxData,"F");
    mtx.Invert();
    r2 = mtx*r2; //x' = A(inv)x
    double xi = (TMath::Pi() - r2.Phi()) + TMath::Pi();
    if (xi > TMath::TwoPi()) xi = xi - TMath::TwoPi();
    return xi;
  }
  else return -10;
}
 
double TGretinaHit::TrackingSort(){
  if (fNumberOfInteractions > 7) return -1; //7 is hardcoded in as the max
  if (fNumberOfInteractions < 2) return 0;
  
  std::array<int,7> indices;
  std::array<int,7> best_order;
  for (int i=0; i < fNumberOfInteractions; i++) {
    indices[i] = i;
    best_order[i] = i;
  }

  double Etot = GetCoreEnergy();
  double Esum = 0; //sum int point energies to scale them
  
  //store the interaction point positions and energies
  std::array<TVector3,7> pos;
  std::array<double,7> eng;
  for (int i=0; i < fNumberOfInteractions; i++) {
    pos[i] = GetIntPosition(i);
    eng[i] = fSegments[i].fEng;
    Esum += eng[i];
  }
  
  // scale energies in case Esum =/= Etot
  double scaleFactor = Etot/Esum;
  for (int i=0; i < fNumberOfInteractions; i++) eng[i] *= scaleFactor;

  // go over all permutations to find the optimal order
  double chi2 = 1e22;
  do {
    double Eg = Etot;
    double temp_chi2 = 0;
    for (int i=0; i < fNumberOfInteractions-1; i++) {
      //get the "energy" cosine
      double Eprime = Eg - eng[indices[i]];
      double cosE = 1 - 511/Eg * (Eg/Eprime - 1);
      Eg = Eprime;
      
      //get the scattering angle cos
      TVector3 posdiff = pos[indices[i+1]] - pos[indices[i]];
      double cosV = posdiff.Dot(pos[indices[i]])/pos[indices[i]].Mag()/posdiff.Mag();
      
      //compare the two via chi2
      double cosDiff = cosE - cosV;
      temp_chi2 += cosDiff*cosDiff;
    }

    // temp_chi2 = TMath::Sqrt(temp_chi2)/(fNumberOfInteractions-1); // Torben's FOM
    temp_chi2 = temp_chi2/fNumberOfInteractions; // Dirk's FOM

    if (temp_chi2 < chi2){
      chi2 = temp_chi2;
      for (int i=0; i < fNumberOfInteractions; i++) best_order[i] = indices[i];
    }

  } while (std::next_permutation(indices.begin(), indices.begin() + fNumberOfInteractions));

  std::array<interaction_point, 7> tmp;
  int point_order = 0;
  for (int i = 0; i < fNumberOfInteractions; ++i) tmp[i] = fSegments[best_order[i]];
  for (int i = 0; i < fNumberOfInteractions; ++i) fSegments[i] = tmp[i];

  return chi2;
}

void TGretinaHit::PICCSort(){
  if (fNumberOfInteractions < 2) return;
  
  double FOM = 1e10;
  int FP = 0; 
  int SP = 1;
  double E = GetCoreEnergy();
  
  //scale the interaction points so they match the core energy
  double scaleFactor = 0;
  for (int i=0; i < fNumberOfInteractions; i++) scaleFactor += fSegments[i].fEng;
  scaleFactor = E/scaleFactor;

  //find the best first two interaction points that satisfy the minimization function
  for (int fp=0; fp < fNumberOfInteractions; fp++){
    double E1 = fSegments[fp].fEng*scaleFactor;
    double cosE = 1 - 511.0/E * E1/(E - E1);
    double z1 = GetLocalPosition(fp).Z();
    // double lp1 = lastPointPenalty(E1);
    double er = (E - E1)/E; 

    for (int sp=0; sp < fNumberOfInteractions; sp++){
      if (fp == sp) continue;
      double scatter_angle = GetScatterAngle(fp,sp);
      double cosP = TMath::Cos(scatter_angle);
      double sinP = TMath::Sin(scatter_angle);
      double E2 = fSegments[sp].fEng*scaleFactor;
      double z2 = GetLocalPosition(sp).Z();
      double diff = std::abs(cosE - cosP);
      double kn = er*er*(er + 1./er - sinP*sinP); 
      // double lp2 = lastPointPenalty(E2);
      double ffom = std::pow(diff,2.0/3)*GetAlpha(fp,sp)*kn //compton scattering
              *std::pow(E/E1,3)*std::pow(E/E2,2) //energies
              *z1*z2; //depth

      if (ffom < FOM) {
        FOM = ffom;
        FP = fp;
        SP = sp;
      }
    }
  }

  //place the interaction points
  interaction_point ipFP = fSegments[FP];
  interaction_point ipSP = fSegments[SP];
  if (SP > FP){
    fSegments.erase(fSegments.begin() + SP);
    fSegments.erase(fSegments.begin() + FP);
  } else {
    fSegments.erase(fSegments.begin() + FP);
    fSegments.erase(fSegments.begin() + SP);
  }
  fSegments.insert(fSegments.begin(),ipSP);
  fSegments.insert(fSegments.begin(),ipFP);

  return;
}

std::map<int,int> TGretinaHit::EquivalentPointMap(const TGretinaHit &comp) {
  std::map<int,int> outmap;

  if (fNumberOfInteractions != comp.NumberOfInteractions()) return outmap;

  for (int i=0; i < fNumberOfInteractions; i++){
    for (int j=0; j < fNumberOfInteractions; j++){
      TVector3 posdiff = GetLocalPosition(i) - comp.GetLocalPosition(j);
      if (GetSegmentEng(i) == comp.GetSegmentEng(j) && posdiff.Mag() == 0) {
        outmap.insert(std::pair<int,int> (j,i)); //keys are the comp indeces
        break;
      }
    }
  }

  if (fNumberOfInteractions != (int) outmap.size()) {
    outmap.clear();
    return outmap;
  }

  return outmap;
}