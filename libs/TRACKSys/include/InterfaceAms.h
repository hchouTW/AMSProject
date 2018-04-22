#if defined(_PGTRACK_) || defined(__ROOTSHAREDLIBRARY__)
#ifndef __TRACKLibs_InterfaceAms_H__
#define __TRACKLibs_InterfaceAms_H__


namespace TrackSys {
namespace InterfaceAms {

// TODO :
static std::vector<HitStTRK> GetHitStTRK(const TrTrackR& track, Int_t patt, Bool_t hasQ = false);
static std::vector<HitStTOF> GetHitStTOF(const BetaHR& betaH, Bool_t hasQ = false, Bool_t hasT = false);


} // namespace InterfaceAms
} // namespace TrackSys


#endif // __TRACKLibs_InterfaceAms_H__
#endif // _PGTRACK_ __ROOTSHAREDLIBRARY__ 
