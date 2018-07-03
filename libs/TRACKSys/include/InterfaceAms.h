#if defined(_PGTRACK_) || defined(__ROOTSHAREDLIBRARY__)
#ifndef __TRACKLibs_InterfaceAms_H__
#define __TRACKLibs_InterfaceAms_H__

#include "root.h"
#include "TrTrack.h"
#include "TrCharge.h"
#include "Tofrec02_ihep.h"


namespace TrackSys {
namespace InterfaceAms {


class TkOpt {
    public :
        static constexpr Int_t QOptISS = TrClusterR::kAsym | TrClusterR::kGain | TrClusterR::kLoss | TrClusterR::kMIP | TrClusterR::kAngle;
        static constexpr Int_t QOptMC  = TrClusterR::kTotSign2017 | TrClusterR::kSimAsym | TrClusterR::kSimSignal | TrClusterR::kLoss | TrClusterR::kAngle;
    
    public :
        enum Patt { MaxSpan, Inner, InnerL1, InnerL9, FullSpan };

    public :
        TkOpt(Patt patt = Patt::Inner, Bool_t dedx = false) : used_(true), patt_(patt), dedx_(dedx) {}
        ~TkOpt() {}

        inline const Bool_t& used() const { return used_; }
        inline const Patt&   patt() const { return patt_; }
        inline const Bool_t& dedx() const { return dedx_; }

        inline Bool_t reqL1() const { return (Patt::InnerL1 == patt_ || Patt::FullSpan == patt_); }
        inline Bool_t reqL9() const { return (Patt::InnerL9 == patt_ || Patt::FullSpan == patt_); }
        
        inline Bool_t useL1() const { return (Patt::InnerL1 == patt_ || Patt::FullSpan == patt_ || Patt::MaxSpan == patt_); }
        inline Bool_t useL9() const { return (Patt::InnerL9 == patt_ || Patt::FullSpan == patt_ || Patt::MaxSpan == patt_); }

    private :
        Bool_t used_;
        Patt   patt_;
        Bool_t dedx_;
};


class TfOpt {
    public :
        static constexpr Int_t QOpt = TofRecH::kThetaCor | TofRecH::kBirkCor | TofRecH::kReAttCor | TofRecH::kDAWeight | TofRecH::kQ2Q;

    public :
        TfOpt() : time_(false), dedx_(false), used_(false) {}
        TfOpt(Bool_t dedx) : time_(true), dedx_(dedx) { used_ = (time_ || dedx_); }
        ~TfOpt() {}

        inline const Bool_t& used() const { return used_; }
        inline const Bool_t& time() const { return time_; }
        inline const Bool_t& dedx() const { return dedx_; }

    private :
        Bool_t used_;
        Bool_t time_;
        Bool_t dedx_;
};


class RhOpt {
    public :
        enum Rad { NO, AGL, NAF };

    public :
        RhOpt(Rad rad = Rad::NO) : rad_(rad) { used_ = (rad_ != Rad::NO); }
        ~RhOpt() {}

        inline const Bool_t& used() const { return used_; }
        inline const Rad&    rad()  const { return rad_; }

    private :
        Bool_t used_;
        Rad    rad_;
};


class Event {
    public :
        Event() {}
        ~Event() {}

    public : 
        inline static void SetArg(
            const TrFitPar::Orientation& ortt = TrFitPar::Orientation::kDownward,
            const VirtualHitSt::NoiseController& noise_ctler = VirtualHitSt::NoiseController::ON,
            const Bool_t& sw_mscat = PhyArg::OptMscat(), const Bool_t& sw_eloss = PhyArg::OptEloss())
        { ArgOrtt = ortt; ArgNoiseCtler = noise_ctler; ArgSwMscat = sw_mscat; ArgSwEloss = sw_eloss; }

    protected :
        static TrFitPar::Orientation         ArgOrtt;
        static VirtualHitSt::NoiseController ArgNoiseCtler;
        static Bool_t                        ArgSwMscat;
        static Bool_t                        ArgSwEloss;

    public :
        static Bool_t   Status() { return (Ev != nullptr && StatusTk); }

        static void     Clear();
        static Bool_t   Load(AMSEventR* event = nullptr, UInt_t ipart = 0, Bool_t rebuild = true);
        static TrFitPar Get(const PartInfo& info = PartInfo(PartType::Proton), const TkOpt& tkOpt = TkOpt(), const TfOpt& tfOpt = TfOpt(), const RhOpt& rhOpt = RhOpt());
   
    protected :
        // AMS Event
        static UInt_t     RunID;
        static UInt_t     EvID;
        static UInt_t     PtID;
        static AMSEventR* Ev;
        static TrTrackR*  Trtk;
        static BetaHR*    Btah;
        static RichRingR* Rich;
        static Bool_t     StatusTk;
        static Bool_t     StatusTf;
        static Bool_t     StatusRh;
    
    protected :
        static void   Init();
        static Bool_t BulidHitStTRK();
        static Bool_t BulidHitStTOF();
        static Bool_t BulidHitStRICH();
        
    protected :
        // TRK
        static std::vector<HitStTRK> TkHitIn;
        static HitStTRK              TkHitL1;
        static HitStTRK              TkHitL9;
        
        static std::vector<HitStTRK> TkHitInQ;
        static HitStTRK              TkHitL1Q;
        static HitStTRK              TkHitL9Q;

        // TOF
        static std::vector<HitStTOF> TfHitT;
        static std::vector<HitStTOF> TfHitQ;
        static std::vector<HitStTOF> TfHitTQ;

        // RICH
        static HitStRICH RhHit;
};


} // namespace InterfaceAms
} // namespace TrackSys


namespace TrackSys {

using AmsTkOpt = InterfaceAms::TkOpt;
using AmsTfOpt = InterfaceAms::TfOpt;
using AmsRhOpt = InterfaceAms::RhOpt;
using AmsEvent = InterfaceAms::Event;

} // namespace TrackSys


#endif // __TRACKLibs_InterfaceAms_H__
#endif // _PGTRACK_ __ROOTSHAREDLIBRARY__ 
