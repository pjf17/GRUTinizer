#ifndef TBANK29_H
#define TBANK29_H

#include <TObject.h>
#include <TMath.h>

#include <TClonesArray.h>

#include "TDetector.h"
#include "TMode3.h"

class TBank29 : public TDetector {

public:
  TBank29();
  ~TBank29();

  virtual void Copy(TObject& obj) const;
  virtual void Clear(Option_t *opt = "");
  virtual void Print(Option_t *opt = "") const;

  virtual void          InsertHit(const TDetectorHit& hit);
  virtual TDetectorHit& GetHit(int i)        { return channels.at(i); }
  virtual TMode3Hit&    GetMode3Hit(int i)   { return channels.at(i); }


  Long_t GetTimestamp() { return Timestamp(); }

private:
  virtual int BuildHits(std::vector<TRawEvent>& raw_data);

  std::vector<TMode3Hit> channels;

  ClassDef(TBank29,1);
};





#endif
