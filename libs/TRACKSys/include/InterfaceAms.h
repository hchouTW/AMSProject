#if defined(_PGTRACK_) || defined(__ROOTSHAREDLIBRARY__)
#ifndef __TRACKLibs_InterfaceAms_H__
#define __TRACKLibs_InterfaceAms_H__

#include "root.h"
#include "TrTrack.h"
#include "TrCharge.h"
#include "Tofrec02_ihep.h"


namespace TrackSys {
namespace InterfaceAms {

/*
class AMSEvent {
    public :
        AMSEvent(AMSEventR* event) {}
        ~AMSEvent() {}
    
    protected :
        HitStTRK              get_HitStTRK(TrRecHitR& hit, Bool_t hasQ = false);
        std::vector<HitStTOF> get_HitStTOF(BetaHR& betaH, Bool_t hasQ = false);
        HitStRICH             get_HitStRICH(RichRingR& rich);

    private :
        // TRK
        std::vector<HitStTRK> trHitIn;
        HitStTRK              trHitL1;
        HitStTRK              trHitL9;

        // TOF
        std::vector<HitStTOF> tfHit;

        // RICH
        HitStRICH rhHit;
};
*/




enum class TrackerPatt { MaxSpan, Inner, InnerL1, InnerL9, FullSpan };
constexpr Int_t QOptTracker = TrClusterR::kTotSign2017 | TrClusterR::kSimAsym | TrClusterR::kSimSignal | TrClusterR::kLoss | TrClusterR::kAngle;
constexpr Int_t QOptTOF     = TofRecH::kThetaCor | TofRecH::kBirkCor | TofRecH::kReAttCor | TofRecH::kDAWeight | TofRecH::kQ2Q;

// TODO :
static std::vector<HitStTRK> GetHitStTRK(TrTrackR& track, const TrackerPatt& patt = TrackerPatt::MaxSpan, Bool_t hasQ = false);
static std::vector<HitStTOF> GetHitStTOF(BetaHR& betaH, Bool_t hasQ = false);
static HitStRICH             GetHitStRICH(RichRingR& rich);


} // namespace InterfaceAms
} // namespace TrackSys


#endif // __TRACKLibs_InterfaceAms_H__
#endif // _PGTRACK_ __ROOTSHAREDLIBRARY__ 
