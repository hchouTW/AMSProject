// ---- Before the event loop ----

#include “EcalHadron.h”
...
EcalHadron::InitParameters();


// ---- In the event loop ----

// 1. Get a pointer to EcalShowerR* Shower
//    I usually take the one that best matches the TrTrack, OR the highest edep one.

// 2. Determine/guess an integer charge Z

// 3. Build the ecal energy: EcalHadron::Build(Shower, Z); 

// 4. Get Ecal energy (actually rigidity, GV units) and shower Apex: 
   double R = EcalHadron::EcalRigidity; // GV
   int Apex = EcalHadron::EcalApex;     // [1 to 19]

// I usually reject high-apex events (say, Apex>12) i.e. CRs showering in the bottom or MIP-like.


Nicola
