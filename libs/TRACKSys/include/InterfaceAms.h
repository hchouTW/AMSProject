#if defined(_PGTRACK_) || defined(__ROOTSHAREDLIBRARY__)
#ifndef __TRACKLibs_InterfaceAms_H__
#define __TRACKLibs_InterfaceAms_H__

#include "root.h"
#include "TrTrack.h"
#include "TrCharge.h"
#include "Tofrec02_ihep.h"


namespace TrackSys {
namespace InterfaceAms {


class Event {
    private :
        static constexpr Int_t QOptTracker = TrClusterR::kTotSign2017 | TrClusterR::kSimAsym | TrClusterR::kSimSignal | TrClusterR::kLoss | TrClusterR::kAngle;
        static constexpr Int_t QOptTOF     = TofRecH::kThetaCor | TofRecH::kBirkCor | TofRecH::kReAttCor | TofRecH::kDAWeight | TofRecH::kQ2Q;

    public :
        enum class TrackerPatt { MaxSpan, Inner, InnerL1, InnerL9, FullSpan };

    public :
        Event(AMSEventR* event, Bool_t withQ = false);
        ~Event() { init(); }

        TrFitPar get(const PartInfo& info = PartInfo(PartType::Proton), const TrackerPatt& trPatt = TrackerPatt::MaxSpan, Bool_t withTOF = false, Bool_t withRICH = false);
        
        const Bool_t& status() const { return status_; }
        const Short_t& going() const { return going_; }

    protected :
        void init();

        Bool_t bulid_HitStTRK();
        Bool_t bulid_HitStTOF();
        Bool_t bulid_HitStRICH();

    private :
        Bool_t  status_;
        Bool_t  withQ_;
        Short_t going_; // 0, no  1, up-going  -1, down-going

        // AMS Event
        AMSEventR* event_;
        TrTrackR*  trtk_;
        BetaHR*    btah_;
        RichRingR* rich_;

        // TRK
        std::vector<HitStTRK> trHitIn_;
        HitStTRK              trHitL1_;
        HitStTRK              trHitL9_;

        // TOF
        std::vector<HitStTOF> tfHit_;

        // RICH
        HitStRICH  rhHit_;
};












/*
enum class TrackerPatt { MaxSpan, Inner, InnerL1, InnerL9, FullSpan };
constexpr Int_t QOptTracker = TrClusterR::kTotSign2017 | TrClusterR::kSimAsym | TrClusterR::kSimSignal | TrClusterR::kLoss | TrClusterR::kAngle;
constexpr Int_t QOptTOF     = TofRecH::kThetaCor | TofRecH::kBirkCor | TofRecH::kReAttCor | TofRecH::kDAWeight | TofRecH::kQ2Q;

// TODO :
static std::vector<HitStTRK> GetHitStTRK(TrTrackR& track, const TrackerPatt& patt = TrackerPatt::MaxSpan, Bool_t hasQ = false);
static std::vector<HitStTOF> GetHitStTOF(BetaHR& betaH, Bool_t hasQ = false);
static HitStRICH             GetHitStRICH(RichRingR& rich);
*/

} // namespace InterfaceAms
} // namespace TrackSys


#endif // __TRACKLibs_InterfaceAms_H__
#endif // _PGTRACK_ __ROOTSHAREDLIBRARY__ 
