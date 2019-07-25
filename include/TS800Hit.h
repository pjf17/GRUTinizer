#ifndef TSEIGHTHUNDRADHIT_H
#define TSEIGHTHUNDRADHIT_H

#include <TObject.h>
#include <TClass.h>
#include <iostream>

#include <TDetectorHit.h>
#include <TMath.h>
#include "TRandom.h"
#include "GValue.h"

class TF1;

#define MAXCRDC 513

class TS800Channel : public TDetectorHit {
  public:
    TS800Channel()                     { Clear();    }
    TS800Channel(short value)          { Set(value); }
    TS800Channel(unsigned short value) { Set(value); }
    ~TS800Channel()                    {}

    void Set(short value)          { fValue = (unsigned short)value; }
    void Set(unsigned short value) { fValue = value; }

    virtual short GetId()      const { return (((fValue)&0xf000)>>12);  }
    short GetValue()           const { return fValue&0x0fff;  }

    virtual void Clear(Option_t *opt="")       { TDetectorHit::Clear(opt); fValue = 0; }
    virtual void Print(Option_t *opt="") const { printf("[%i] = %i\n",GetId(),GetValue());}
    virtual void Copy(TObject &obj)      const { TDetectorHit::Copy(obj); ((TS800Channel&)obj).fValue = fValue; }

    virtual int  Charge() const { return GetValue(); }

  private:
    short fValue;

  ClassDef(TS800Channel,1);
};

class TTrigger : public TDetectorHit {
  public:
    TTrigger();
    ~TTrigger();

    void SetRegistr(unsigned short reg) { fregistr=reg; }
    void SetS800Source(short sou)       { fs800source=sou; }
    void SetExternalSource1(short sou)  { fexternalsource1=sou; }
    void SetExternalSource2(short sou)  { fexternalsource2=sou; }
    void SetSecondarySource(short sou)  { fsecondarysource=sou; }

    unsigned short GetRegistr() const { return fregistr; }
    short GetS800Source()       const { return fs800source; }
    short GetExternalSource1()  const { return fexternalsource1; }
    short GetExternalSource2()  const { return fexternalsource2; }
    short GetSecondarySource()  const { return fsecondarysource; }

    virtual void Copy(TObject &)         const;
    virtual void Print(Option_t *opt="") const;
    virtual void Clear(Option_t *opt="");

  private:
    virtual int Charge() const { return 0; }

    unsigned short fregistr;
    short fs800source;
    short fexternalsource1;
    short fexternalsource2;
    short fsecondarysource;

  ClassDef(TTrigger,1)
};


class TTof : public TDetectorHit { // S800 Time of Flight
  public:
    TTof();
    ~TTof();

    void SetRF(short rf)                { frf=rf; }
    void SetOBJ(short obj)              { fobj=obj; }
    void SetXFP(short xfp)              { fxfp=xfp; }
    void SetSI(short si)                { fsi=si; }
    void SetTacOBJ(short obj)           { ftac_obj=obj; }
    void SetTacXFP(short xfp)           { ftac_xfp=xfp; }

    short GetRF()                         { return frf;     }
    short GetOBJ()                        { return fobj;    }  // tdc!
    short GetXFP()                        { return fxfp;    }  // tdc!
    short GetSI()                         { return fsi;     }
    short GetTacOBJ()                     { return ftac_obj;}  // tac!
    short GetTacXFP()                     { return ftac_xfp;}  // tac!


    virtual void Copy(TObject &)         const;
    virtual void Print(Option_t *opt="") const;
    virtual void Clear(Option_t *opt="");

  private:
    virtual int Charge() const { return 0; }

    short frf;
    short fobj;
    short fxfp;
    short fsi;
    short ftac_obj;
    short ftac_xfp;

  ClassDef(TTof,1);
};

class TCrdc : public TDetectorHit {
  public:
    TCrdc();
    ~TCrdc();

    short GetId()    { return fId;   }
    short GetAnode() { return anode; }
    short GetTime()  { return time;  }
    float GetTimeRand()  { return ((float)(time)+gRandom->Uniform());  }

    int  Size()        const { return channel.size(); }
    int  GetNSamples() const { return sample.size(); }

    void SetId(short id)    { fId = id;  }
    void SetAnode(short an) {anode = an; }
    void SetTime(short ti)  {time = ti;  }


    int  Address(int i) const { return TDetectorHit::Address() + channel.at(i); }

    void AddPoint(int chan,int samp,int dat) { channel.push_back(chan);
                                               sample.push_back(samp);
                                               data.push_back(dat);    }
    int GetChannel(int i) const    { if(i>=Size()) return -1; return channel.at(i);    }
    int GetSample(int i)  const    { if(i>=Size()) return -1; return sample.at(i);     }
    int GetData(int i)    const    { if(i>=Size()) return -1; return data.at(i);       }

    int GetWidth();

    float GetDispersiveX() const;     
    float GetNonDispersiveY();  
    int GetMaxPad() const;
    int GetMaxPadSum() const;

    virtual void Copy(TObject&) const;
    virtual void Print(Option_t *opt="") const;
    virtual void Clear(Option_t *opt="");

    virtual void DrawChannels(Option_t *opt="",bool calibrate=true) const;
    virtual void DrawHit(Option_t *opt="") const;

    int Sum() const { int result = 0; for(unsigned int x=0;x<data.size();x++) result +=data[x]; return result; }

  private:
    virtual int Charge() const { return 0; }

    bool IsGoodSample(int i) const;
    
    short fId;
    std::vector<int> channel;
    std::vector<int> sample;
    std::vector<int> data;

    unsigned short anode;
    unsigned short time;

    mutable bool has_cached_dispersive_x; //!
    mutable double cached_dispersive_x; //!


    static TF1 *fgaus;


  ClassDef(TCrdc,1)
};


class TScintillator : public TDetectorHit {
  public:
    TScintillator();
    ~TScintillator();

    void SetID(int id)         { fID=id; }
    void SetdE_Up(float de)    { fdE_up     = de; }
    void SetdE_Down(float de)  { fdE_down   = de; }
    void SetTime_Up(float t)   { fTime_up   = t; }  // tdc
    void SetTime_Down(float t) { fTime_down = t; }  // tdc

    int GetID()         { return fID;        }
    float GetEUp()      { return fdE_up;     }
    float GetEDown()    { return fdE_down;   }
    float GetTimeUp()   { return fTime_up;   }
    float GetTimeDown() { return fTime_down; }

    virtual void Copy(TObject&) const;
    virtual void Print(Option_t *opt="") const;
    virtual void Clear(Option_t *opt="");

  private:
    virtual int Charge() const { return 0; }
    int fID;
    float fdE_up;
    float fdE_down;
    float fTime_up;
    float fTime_down;

    ClassDef(TScintillator,1)
};

class TIonChamber : public TDetectorHit {
  public:
    TIonChamber();
    ~TIonChamber();

    void Set(int ch, int data);

    int GetChannel(int i) const { if(i>=Size()) return -1; return fChan.at(i); }
    int GetData(int i)    const { if(i>=Size()) return -1; return fData.at(i); }
    int Size() const { return fChan.size(); }
    float GetdE(TCrdc *);
    float GetdE(double crdc_1_x, double crdc_1_y);
    float GetSum() const;
    float GetAve();

    int  Address(int i) const { return TDetectorHit::Address() + GetChannel(i); }


    float GetCalData(int i) const {
      TChannel *c = TChannel::GetChannel(Address(GetChannel(i)));
      if(c){
        return c->CalEnergy(GetData(i));
      }else{
       return (float)GetData(i);
      }
    }

    virtual void Copy(TObject&) const;
    virtual void Print(Option_t *opt="") const;
    virtual void Clear(Option_t *opt="");
    int Charge() const { int sum=0;for(int i=0;i<Size();i++)sum+=GetData(i);return sum;}
  private:

    std::vector<int> fChan;
    std::vector<int> fData;
    ClassDef(TIonChamber,1)
};

class THodoHit : public TDetectorHit{
  public:
    THodoHit() { Clear(); }
    THodoHit(const THodoHit&);

    virtual void Copy(TObject &) const;
    virtual void Print(Option_t *opt="") const;
    virtual void Clear(Option_t *opt="");


    void SetChannel(int channel) { fChannel = channel;} 
    
    int GetCharge() const { return Charge(); }
    int GetChannel() const { return fChannel; }
  private:
    int fChannel;

    ClassDef(THodoHit, 1); 
};

class THodoscope : public TDetectorHit {
  public:
    THodoscope();
    virtual ~THodoscope();
    THodoscope(const THodoscope&);

    virtual void Copy(TObject&) const;
    virtual void Print(Option_t *opt="") const;
    virtual void Clear(Option_t *opt="");
    void InsertHit(const TDetectorHit&);

    int  Address(int i) const { return TDetectorHit::Address() + GetChannel(i); }

    int GetChannel(int i) const { return GetHodoHit(i).GetChannel();} 
    TDetectorHit& GetHit(int i);
    THodoHit& GetHodoHit(int i);
    const THodoHit& GetHodoHit(int i) const;
    std::vector<THodoHit> GetHodoHits() const { return hodo_hits; } ;
    
    size_t Size() const { return hodo_hits.size(); } 

    std::vector<THodoHit> hodo_hits;

    ClassDef(THodoscope, 1); 

};


class TMTof : public TDetectorHit {
  public:
    TMTof();
    ~TMTof();
    TMTof(const TMTof&);

    virtual void Copy(TObject&) const;
    virtual void Print(Option_t *opt="") const;
    virtual void Clear(Option_t *opt="");

    int E1UpSize()       const { return fE1Up.size();   }
    int E1DownSize()     const { return fE1Down.size();  }
    int XfpSize()        const { return fXfp.size();      }
    int ObjSize()        const { return fObj.size();       }
    int RfSize()         const { return fRf.size();        }
    int Crdc1AnodeSize() const { return fCrdc1Anode.size();}
    int Crdc2AnodeSize() const { return fCrdc2Anode.size();}
    int HodoSize()       const { return fHodoscope.size(); }
    int RefSize()        const { return fRef.size(); }

    //Determines correlated time-of-flights for multi-hit TDC's based on
    //TARGET_MTOF_# GValues. Values are set in TMTof as fCorrelatedXfp,
    //fCorrelatedXfp, etc.  Note that if the GValues are not set, the first
    //value in each time-of-flight is taken.
    double  GetCorrelatedXfpE1()  const;   //!
    double  GetCorrelatedObjE1()  const;   //!
    double  GetCorrelatedXfpE1Chn15()  const;   //!
    double  GetCorrelatedObjE1Chn15()  const;   //!



  private:
    mutable double fCorrelatedXFPE1;   //!
    mutable double fCorrelatedOBJE1;   //!

  public:
    std::vector<unsigned short> fE1Up;         // Channel 0
    std::vector<unsigned short> fE1Down;       // Channel 1
    std::vector<unsigned short> fXfp;          // Channel 2
    std::vector<unsigned short> fObj;          // Channel 3
    std::vector<unsigned short> fRf;           // Channel 5
    std::vector<unsigned short> fCrdc1Anode;   // Channel 6
    std::vector<unsigned short> fCrdc2Anode;   // Channel 7
    std::vector<unsigned short> fHodoscope;    // Channel 12
    std::vector<unsigned short> fRef;          // Channel 15, same as E1Up (different cable.)
    virtual Int_t Charge() const  {return 0;}
  ClassDef(TMTof,1)
};


class TS800Hit : public TDetectorHit {
  public:
    TS800Hit()   {  }
    ~TS800Hit()  {  }
    virtual Int_t Charge() const { return 0; }

    virtual void Clear(Option_t *opt ="")       { TDetectorHit::Clear(opt); }
    virtual void Print(Option_t *opt ="") const {  }
    virtual void Copy(TObject &obj)       const { TDetectorHit::Copy(obj);  }

  private:

  ClassDef(TS800Hit,1);
};
#endif
