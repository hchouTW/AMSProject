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
        enum Patt { MaxSpan, Inner, InnerL1, InnerL9, FullSpan };

    public :
        TkOpt(Patt patt = Patt::Inner, Bool_t mesc = false, Bool_t dedx = false) : patt_(patt), mesc_(mesc), dedx_(dedx) { used_ = (mesc_ || dedx_); }
        ~TkOpt() {}

        inline const Bool_t& used() const { return used_; }
        inline const Patt&   patt() const { return patt_; }
        inline const Bool_t& mesc() const { return mesc_; }
        inline const Bool_t& dedx() const { return dedx_; }

        inline Bool_t reqL1() const { return (Patt::InnerL1 == patt_ || Patt::FullSpan == patt_); }
        inline Bool_t reqL9() const { return (Patt::InnerL9 == patt_ || Patt::FullSpan == patt_); }
        
        inline Bool_t useL1() const { return (Patt::InnerL1 == patt_ || Patt::FullSpan == patt_ || Patt::MaxSpan == patt_); }
        inline Bool_t useL9() const { return (Patt::InnerL9 == patt_ || Patt::FullSpan == patt_ || Patt::MaxSpan == patt_); }

    private :
        Bool_t used_;
        Patt   patt_;
        Bool_t mesc_;
        Bool_t dedx_;
};


class TfOpt {
    public :
        static constexpr Int_t QOpt = TofRecH::kThetaCor | TofRecH::kBirkCor | TofRecH::kReAttCor | TofRecH::kDAWeight | TofRecH::kQ2Q;

    public :
        TfOpt(Bool_t used = false, Bool_t dedx = false) : used_(used) { dedx_ = (used_ && dedx); }
        ~TfOpt() {}

        inline const Bool_t& used() const { return used_; }
        inline const Bool_t& dedx() const { return dedx_; }

    private :
        Bool_t used_;
        Bool_t dedx_;
};


class RhOpt {
    public :
        RhOpt(Bool_t used = false) : used_(used) {}
        ~RhOpt() {}

        inline const Bool_t& used() const { return used_; }

    private :
        Bool_t used_;
};


class Event {
    public :
        Event() {}
        ~Event() {}

    public :
        static void   Clear();
        
        static Bool_t Status() { return (Ev != nullptr && StatusTk); }
        static Bool_t Load(AMSEventR* event = nullptr, UInt_t ipart = 0);
        static Bool_t HasTRK()  { return StatusTk; }
        static Bool_t HasTOF()  { return StatusTf; }
        static Bool_t HasRICH() { return StatusRh; }
        static Bool_t HasRICH_AGL() { return (StatusRh && (HitStRICH::Radiator::AGL == RhHit.rad())); }
        static Bool_t HasRICH_NAF() { return (StatusRh && (HitStRICH::Radiator::NAF == RhHit.rad())); }

        static const RichObjAms& GetRichObj() { return RichObj; }

    public : 
        static TrFitPar GetTrFitPar(const PartInfo& info = PartInfo(PartType::Proton), const TrFitPar::Orientation& ortt = TrFitPar::Orientation::kDownward, const Bool_t& sw_mscat = PhyArg::OptMscat(), const Bool_t& sw_eloss = PhyArg::OptEloss(), const TkOpt& tkOpt = TkOpt(), const TfOpt& tfOpt = TfOpt(), const RhOpt& rhOpt = RhOpt());

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
    
        static RichObjAms RichObj;

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
        
        static std::vector<HitStTRK> TkHitInQ_NOxy;
        static HitStTRK              TkHitL1Q_NOxy;
        static HitStTRK              TkHitL9Q_NOxy;

        // TOF
        static std::vector<HitStTOF> TfHitT;
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
