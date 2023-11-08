#ifndef TGRETINAHIT_H
#define TGRETINAHIT_H

#include <TObject.h>
#include <Rtypes.h>
#include <TVector3.h>
#include <TMath.h>
#include <TChain.h>

#include <cmath>

#include "TDetectorHit.h"

#define MAXHPGESEGMENTS 36

class TSmartBuffer;
class TS800;

class interaction_point { 
  public:
  interaction_point():fSeg(-1),fX(sqrt(-1)),fY(sqrt(-1)),fZ(sqrt(-1)),fEng(sqrt(-1)),fFrac(sqrt(-1)) { }
  interaction_point(int seg,float x,float y,float z,float energy,float fraction=100.0)
    : fSeg(seg),fX(x),fY(y),fZ(z),fEng(energy),fFrac(fraction) { }
  virtual ~interaction_point() { }
  int   fSeg;
  float fX;
  float fY;
  float fZ;
  float fEng;
  float fFrac;

  bool operator<(const interaction_point &other) const {
    if(fEng!=other.fEng) {
      return fEng>other.fEng;
    } 
    if(fSeg==other.fSeg) {
      return fFrac>other.fFrac;
    }
    return fSeg<other.fSeg;
  }

  void Print(Option_t *opt="") const { 
    printf("Seg[%02i]\tWedge[%i]\tEng: % 4.2f / % 4.2f  \t(X,Y,Z) % 3.2f % 3.2f % 3.2f\n",
            fSeg,(fSeg%6),fEng,fFrac,fX,fY,fZ);
  
  };
  ClassDef(interaction_point,1)
};

#ifndef __CINT__ 

struct intpnt_compare {
  bool operator()(const interaction_point &p1,const interaction_point &p2) {
    return p1.fSeg < p2.fSeg;
  }
};

struct intpnt_compare_wedge {
  bool operator()(const interaction_point &p1,const interaction_point &p2) {
    return (p1.fSeg%6) < (p2.fSeg%6);
  }
};
#endif

class TGretinaHit : public TDetectorHit {

public:
  TGretinaHit();
  TGretinaHit(const TGretinaHit& hit){ hit.Copy(*this); }
  ~TGretinaHit();

  void Copy(TObject& obj) const;

  void BuildFrom(TSmartBuffer& raw);

  TGretinaHit GetNeighbor(int i=0) const;
  TGretinaHit GetInitialHit() const;
  Int_t GetNNeighborHits() const {return fSingles.size()-1;} //don't count first hit

  Double_t GetTime()               const { return (double)Timestamp() + (double)fWalkCorrection; } 
  Float_t  GetT0()                 const { return fWalkCorrection; }
  Float_t  GetTFit()               const { return fWalkCorrection - fTOffset; }
  Float_t  GetTOffset()            const { return fTOffset; }

  Int_t    GetRingNumber() const;
  Int_t    GetCrystalId()          const { return fCrystalId;      }
  Int_t    GetHoleNumber()         const { return fCrystalId/4-1;  }
  Int_t    GetCrystalNumber()      const { return fCrystalId%4;    }
  Float_t  GetCoreEnergy()         const; 
  Int_t    GetCoreCharge(int i)    const { return fCoreCharge[i];  }
  Float_t  GetCoreEnergy(int i)    const;
  virtual Int_t Charge()           const { return GetCoreCharge(3); }
  Int_t GetPad() const { return fPad; }

  const char *GetName() const;

  void  Print(Option_t *opt="") const;
  void  Clear(Option_t *opt="");
  
  Int_t Size()  const { return fSegments.size();  }
  
  double GetX() const { return GetPosition().X(); }
  double GetY() const { return GetPosition().Y(); }
  double GetZ() const { return GetPosition().Z(); }
  
  double GetPhi(int point = 0) const {
    double phi = GetIntPosition(point).Phi();
    if(phi<0) {
      return TMath::TwoPi()+phi;
    } else {
      return phi;
    }
  }
  double GetTheta(int point = 0)    const { return GetIntPosition(point).Theta(); }
  double GetPhiDeg(int point = 0)   const { return GetPhi(point)*TMath::RadToDeg(); }
  double GetThetaDeg(int point = 0) const { return GetTheta(point)*TMath::RadToDeg(); }
  double GetAlpha(int p1=0, int p2=1) const; //get the angle between two interaction points
  double GetScatterAngle(int p1=0, int p2=1) const; //get the polar compton scattering angle 
  double GetXi(const TVector3 *beam=nullptr,int p1=0, int p2=1) const; //get the azimuthal compton scattering angle
  double GetXiChris(const TVector3 *beam=nullptr,int p1=0, int p2=1) const; //get the azimuthal compton scattering angle

  Int_t Compare(const TObject *obj) const { 
    TGretinaHit *other = (TGretinaHit*)obj;
    if(this->GetCoreEnergy()>other->GetCoreEnergy())
      return -1;
    else if(this->GetCoreEnergy()<other->GetCoreEnergy())
      return 1;  //sort largest to smallest.
    return 0;
  }
  
  bool HasInteractions() { return !fSegments.empty(); }
  
  bool operator<(const TGretinaHit &rhs) const { return fCoreEnergy > rhs.fCoreEnergy; }

  double GetDoppler(double beta,const TVector3 *vec=0,int EngRange=-1) const {
    if(Size()<1)
      return 0.0;
    if(vec==0) {
      vec = &BeamUnitVec;
    }
    double tmp = 0.0;
    double gamma = 1/(sqrt(1-pow(beta,2)));
    if(EngRange>0) 
      tmp = GetCoreEnergy(EngRange)*gamma *(1 - beta*TMath::Cos(GetPosition().Angle(*vec)));
    else
      tmp = fCoreEnergy*gamma *(1 - beta*TMath::Cos(GetPosition().Angle(*vec)));
    return tmp;
  } 
  
  double GetDopplerYta(double beta , double yta, const TVector3 *vec=0, int point = -1, int EngRange =-1) const;
  double GetDopplerYta(double beta , double yta, double target_x_shift, double target_y_shift, double target_z_shift, const TVector3 *vec=0, int point = -1, int EngRange =-1) const;
  double GetDoppler(const TS800 *s800,bool doDTAcorr=false,int EngRange=-1);
  double GetDoppler_dB(double beta,const TVector3 *vec=0, double Dta=0);



  Int_t    NumberOfInteractions()        const { return fNumberOfInteractions; }
  Int_t    GetNSegments()                const; //{ return (int)fSegments.size(); }
  Int_t    GetSegmentId(int i=-1)        const { if(i>=GetNSegments()||GetNSegments()==0) return -1;
                                                 if(i==-1) return fSegments.at(0).fSeg;
                                                           return fSegments.at(i).fSeg;  }
  Float_t  GetSegmentEng(const int &i)   const { return fSegments.at(i).fEng;  }
  Int_t    GetSegmentLayer(int i=-1)     const { return std::floor(GetSegmentId(i)/6.0);}

  TVector3 GetIntPosition(unsigned int i)   const;  // position of the ith segment, Global coor.
  TVector3 GetLocalPosition(unsigned int i) const;  // position of the ith segment, Local coor.
  TVector3 GetPosition()                    const { return GetIntPosition(0); }
  TVector3 GetLastPosition()                const;

  TVector3 GetCrystalPosition()           const; 

  //Set segment position; useful for smearing simulated Gretina data
  void SetPosition(unsigned int i, double x, double y, double z);
                                                
  void Add(const TGretinaHit& other);
  void NNAdd(const TGretinaHit& other);
  void NNSwap(int index);
  void SetCoreEnergy(float temp) const { fCoreEnergy = temp; }

  void TrimSegments(int type); // 0: drop multiple ident int pnts.  1: make into wedge "data"
  bool IsClean() const { return !fPad; }

  void ScaleIntEng();
  void SortSegments() { std::sort(fSegments.begin(),fSegments.end());}
  void ReverseSegments() { std::reverse(fSegments.begin(),fSegments.end());}
  void ComptonSort(double cut=0.0);
  void SimTracking(const double realTheta);

private:
  void SortHits();
/* All possible decomp information and
 * where is is stored:
 * -------------------
  Int_t     type;       // endiness identifier; droppped.
  Int_t     crystal_id; //                                      -> TGreintaHit.fCrystalId
  Int_t     num;        // number of interactions of error code -> TGreintaHit.fNumberOfINteractions
  Float_t   tot_e;      // energy used for decomp               -> TGretinaHit.fCoreEnergy   
  Int_t     core_e[4];  // charge reported at dig for each gain -> TGretinaHit.fCoreCharge[4]
  Long_t    timestamp;  // timestamp for the hit                -> TDetectorHit.fTimestamp
  Long_t    trig_time;  // currently unsed (?)
  Float_t   t0;         // rise time as reported by decomp      -> TGretinaHit.fWalkCorrectrion
  Float_t   cfd;        // cfd as reported by the dig.          -> TDetectorHit.fTime
  Float_t   chisq;      // chisq value reported by decomp
  Float_t   norm_chisq; //
  Float_t   baseline;   //
  Float_t   prestep;    // avg trace value before step
  Float_t   poststep;   // avg trave value after step
  Int_t     pad;        // decomp error code.
 *------------------
*/
  Int_t           fCrystalId;
  Int_t           fCoreCharge[4];
  Int_t   fPad;
  Int_t   fNumberOfInteractions; 
  
  mutable Float_t fCoreEnergy;
  Float_t         fWalkCorrection;   //also called t0.
  Float_t         fTOffset; //  t0 = toffset + tFit

  std::vector<TGretinaHit> fSingles;
  bool fSetFirstSingles = false;

  std::vector<interaction_point> fSegments;
  ClassDef(TGretinaHit,5)
};


#endif
