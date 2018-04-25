#if defined(_PGTRACK_) || defined(__ROOTSHAREDLIBRARY__)
#ifndef __TRACKLibs_InterfaceAms_H__
#define __TRACKLibs_InterfaceAms_H__

#include "root.h"
#include "TrTrack.h"
#include "TrCharge.h"
#include "Tofrec02_ihep.h"


namespace TrackSys {
namespace InterfaceAms {

// TODO :
static std::vector<HitStTRK> GetHitStTRK(TrTrackR& track, Int_t patt, Bool_t hasQ = false);
static std::vector<HitStTOF> GetHitStTOF(BetaHR& betaH, Bool_t hasT = true, Bool_t hasQ = false);


} // namespace InterfaceAms
} // namespace TrackSys


#endif // __TRACKLibs_InterfaceAms_H__
#endif // _PGTRACK_ __ROOTSHAREDLIBRARY__ 
