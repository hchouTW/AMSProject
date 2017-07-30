#ifndef __EcalHadron_H__
#define __EcalHadron_H__

#include <root.h>
#include <root_setup.h> 

#include <iostream>
#include <fstream>

class EcalHadron : public TObject {

 private:

  static double ECALPAR[3][8]; //[18]; 
  static double EcalEdepCell[18][72];
  static double EcalEdepPla[18];
  static double EcalS13Pla[18];
  static double EcalS35Pla[18];

 public:

  static int EcalApex;
  static double EcalEdepRaw;
  static double EcalEdep;
  static double EcalRigidity;

  static void InitParameters();
  static void InitEvent();
  static void SetShowerProperties(EcalShowerR* shower);
  static void SetApex();
  static void SetRigidity(int Z);
  static void Build(EcalShowerR* shower, int Z);
};

#endif // __EcalHadron_H__

