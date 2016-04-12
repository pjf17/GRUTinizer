#ifndef _TCAESAR_H_
#define _TCAESAR_H_

#define N_RINGS 10
#define MAX_DETS 24
#define MAX_VSN 12
#define MAX_CHN 16
#define MAX_NEIGHBORS 10

#include "TDetector.h"
#include "TCaesarHit.h"
#include "GValue.h"
#include "TS800.h"

//For easy parsing of vsn map
#include "TEnv.h" 
//For easy parsing of detector positions
#include <fstream>

class TCaesar : public TDetector {



public:
  TCaesar();
  virtual ~TCaesar();

  void Copy(TObject& obj) const;
  virtual void InsertHit(const TDetectorHit&);

  virtual void Print(Option_t* opt="") const; 
  virtual void Clear(Option_t* opt = "");

  virtual TDetectorHit& GetHit(int i);
  //Is there any reason the below are virtual? There is no GetCaesarHit in other
  //classes.
  virtual TCaesarHit&   GetCaesarHit(int i);
  virtual const TCaesarHit& GetCaesarHit(int i) const;
  virtual std::vector<TCaesarHit> GetCaesarHits() const  { return caesar_hits; };
  

  int GetULM() const { return fUlm; }
  void SetULM(int ulm) { fUlm = ulm; }


  
  //static TVector3 GetPosition(int detnum, int ring_num, int sector_num);
  int Size() const { return caesar_hits.size(); }

  double GetEnergyDC(TCaesarHit hit);
  double GetEnergyDC(int ring, int det, double energy);
  double GetCorrTime(TCaesarHit hit, TS800 *s800);


  static TVector3 GetPotition(const TCaesarHit*) const { return GetPosition(hit->GetRingNumber(),hit->GetDetectorNumber()); }
  static TVector3 GetPosition(int ring,int det) const;


  static char const ring_names[]; 
  static int  const det_per_ring[]; 

  static double detector_positions[N_RINGS][MAX_DETS][3];
  static int num_neighbors[N_RINGS][MAX_DETS];

  //neighbors is an array where neighbors[ring][det][i][#]
  //gives you the i-th neighbor (maximum num_neighbors)
  //of the given (ring,det) combo. # refers to either ring (0)
  //or det (1).
  static int neighbors[N_RINGS][MAX_DETS][MAX_NEIGHBORS][2];

  static int vsnchn_ring_map_energy[MAX_VSN][MAX_CHN];
  static int vsnchn_ring_map_time[MAX_VSN][MAX_CHN];
  static int vsnchn_det_map_energy[MAX_VSN][MAX_CHN];
  static int vsnchn_det_map_time[MAX_VSN][MAX_CHN];

  static bool filled_map;
  static bool filled_neighbor_map;
  static bool filled_det_pos;

private:
  void SetCharge(int vsn, int channel, int data);
  void SetTime(int vsn, int channel, int data);
  void ReadDetectorPositions(std::string in_file_name);
  void ReadVSNMap(std::string in_file_name);
  void ReadNeighborMap(std::string in_file_name);
  TCaesarHit& GetHit_VSNChannel(int vsn, int channel);

  virtual int  BuildHits(std::vector<TRawEvent>& raw_data);

  void Build_Single_Read(TSmartBuffer buf);

  int fUlm;
  std::string vsn_file_name      = std::string(getenv("GRUTSYS")) + "/config/caesar/VSNMap.dat";
  std::string det_pos_file_name  = std::string(getenv("GRUTSYS")) + "/config/caesar/CaesarPos.dat";
  std::string neighbor_file_name = std::string(getenv("GRUTSYS")) + "/config/caesar/CaesarNeighbors.dat";
  std::vector<TCaesarHit> caesar_hits;

  ClassDef(TCaesar,1);
};

#endif /* _TCAESAR_H_ */
