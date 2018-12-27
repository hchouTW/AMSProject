#ifndef __TRACKLibs_TrFitPar_C__
#define __TRACKLibs_TrFitPar_C__


#include "Sys.h"
#include "Math.h"
#include "CooMeas.h"
#include "TmeMeas.h"
#include "IonEloss.h"
#include "IonTrEloss.h"
#include "PartInfo.h"
#include "PhySt.h"
#include "MagEnv.h"
#include "MatEnv.h"
#include "Prop.h"
#include "HitSt.h"
#include "TrFitPar.h"


namespace TrackSys {


TrFitPar& TrFitPar::operator=(const TrFitPar& rhs) {
    if (this != &rhs) {
        sw_mscat_    = rhs.sw_mscat_;
        sw_eloss_    = rhs.sw_eloss_;
        info_        = rhs.info_;
        ortt_        = rhs.ortt_;
        hits_TRK_    = rhs.hits_TRK_;
        hits_TOF_    = rhs.hits_TOF_;
        hits_RICH_   = rhs.hits_RICH_;
        hits_TRD_    = rhs.hits_TRD_;
        onlycx_nseq_ = rhs.onlycx_nseq_;
        onlycy_nseq_ = rhs.onlycy_nseq_;
        onlyc_nseq_  = rhs.onlyc_nseq_;
        nseq_        = rhs.nseq_;
        nseg_        = rhs.nseg_;
        nmes_        = rhs.nmes_;
        nmes_cx_     = rhs.nmes_cx_;
        nmes_cy_     = rhs.nmes_cy_;
        nmes_ib_     = rhs.nmes_ib_;
        nmes_TRKq_   = rhs.nmes_TRKq_;
        nmes_TOFt_   = rhs.nmes_TOFt_;
        nmes_TOFq_   = rhs.nmes_TOFq_;
        nmes_RICHib_ = rhs.nmes_RICHib_;
        nmes_TRDel_  = rhs.nmes_TRDel_;
        is_check_    = rhs.is_check_;
        
        hits_.clear();
        if (is_check_) {
            for (auto&& hit : hits_TRK_ ) hits_.push_back(&hit); 
            for (auto&& hit : hits_TOF_ ) hits_.push_back(&hit); 
            for (auto&& hit : hits_RICH_) hits_.push_back(&hit); 
            for (auto&& hit : hits_TRD_)  hits_.push_back(&hit); 
            
            if (ortt_ == Orientation::kDownward) VirtualHitSt::Sort(hits_, VirtualHitSt::Orientation::kDownward);
            else                                 VirtualHitSt::Sort(hits_, VirtualHitSt::Orientation::kUpward);
        }
        else zero();

        timer_ = rhs.timer_;
    }
   
    return *this;
}
   

void TrFitPar::zero() {
    hits_.clear();

    onlycx_nseq_ = 0;
    onlycy_nseq_ = 0;
    onlyc_nseq_  = 0;
    nseq_        = 0;
    nseg_        = 0;
    nmes_        = 0;
    nmes_cx_     = 0;
    nmes_cy_     = 0;
    nmes_ib_     = 0;
    nmes_TRKq_   = 0;
    nmes_TOFt_   = 0;
    nmes_TOFq_   = 0;
    nmes_RICHib_ = 0;
    nmes_TRDel_  = 0;
    
    is_check_ = false;

    timer_.clear();
}


void TrFitPar::clear() {
    sw_mscat_ = false;
    sw_eloss_ = false;
    
    info_ = PartInfo(PartType::Proton);
    ortt_ = Orientation::kDownward;
   
    hits_TRK_.clear();
    hits_TOF_.clear();
    hits_RICH_.clear();
    hits_TRD_.clear();

    zero();
}
        
void TrFitPar::set_info(const PartInfo& info, const Orientation& ortt, const Bool_t& sw_mscat, const Bool_t& sw_eloss) {
    sw_mscat_ = sw_mscat;
    sw_eloss_ = sw_eloss;
    info_ = info;
    ortt_ = ortt;

    zero();
}


Bool_t TrFitPar::sort_hits() {
    zero();

    if (ortt_ == Orientation::kDownward) Hit<HitStTRK>::Sort(hits_TRK_, VirtualHitSt::Orientation::kDownward);
    else                                 Hit<HitStTRK>::Sort(hits_TRK_, VirtualHitSt::Orientation::kUpward);
    if (ortt_ == Orientation::kDownward) Hit<HitStTOF>::Sort(hits_TOF_, VirtualHitSt::Orientation::kDownward);
    else                                 Hit<HitStTOF>::Sort(hits_TOF_, VirtualHitSt::Orientation::kUpward);
    if (ortt_ == Orientation::kDownward) Hit<HitStRICH>::Sort(hits_RICH_, VirtualHitSt::Orientation::kDownward);
    else                                 Hit<HitStRICH>::Sort(hits_RICH_, VirtualHitSt::Orientation::kUpward);
    if (ortt_ == Orientation::kDownward) Hit<HitStTRD>::Sort(hits_TRD_, VirtualHitSt::Orientation::kDownward);
    else                                 Hit<HitStTRD>::Sort(hits_TRD_, VirtualHitSt::Orientation::kUpward);
 
    for (auto&& hit : hits_TRK_ ) hits_.push_back(&hit); 
    for (auto&& hit : hits_TOF_ ) hits_.push_back(&hit); 
    for (auto&& hit : hits_RICH_) hits_.push_back(&hit); 
    for (auto&& hit : hits_TRD_ ) hits_.push_back(&hit); 
    
    if (ortt_ == Orientation::kDownward) VirtualHitSt::Sort(hits_, VirtualHitSt::Orientation::kDownward);
    else                                 VirtualHitSt::Sort(hits_, VirtualHitSt::Orientation::kUpward);
    
    Short_t onlycx_nseq =  0;
    Short_t onlycy_nseq =  0;
    Short_t onlyc_nseq  =  0;
    Short_t nseq        =  0;
    Short_t nseg        = -1;
    for (auto&& hit : hits_) {
        hit->set_type(info_);
        onlycx_nseq += hit->set_onlycx_seqID(onlycx_nseq);
        onlycy_nseq += hit->set_onlycy_seqID(onlycy_nseq);
        onlyc_nseq  += hit->set_onlyc_seqID(onlyc_nseq);
        nseq        += hit->set_seqID(nseq);
        
        if (hit->scx()) nmes_cx_++;
        if (hit->scy()) nmes_cy_++;
        if (hit->scx() || hit->scy()) nseg++;
    }

    for (auto&& hit : hits_TRK_) {
        if (hit.sq()) nmes_TRKq_++;
    }
        
    for (auto&& hit : hits_TOF_) {
        if (hit.st()) nmes_TOFt_++;
        if (hit.sq()) nmes_TOFq_++;
    }
    
    for (auto&& hit : hits_RICH_) {
        if (hit.sib()) nmes_RICHib_++;
    }
    
    for (auto&& hit : hits_TRD_) {
        if (hit.sel()) nmes_TRDel_++;
    }
    
    nmes_ib_ = nmes_TRKq_ + nmes_TOFt_ + nmes_TOFq_ + nmes_RICHib_ + nmes_TRDel_;
    
    onlycx_nseq_ = onlycx_nseq;
    onlycy_nseq_ = onlycy_nseq;
    onlyc_nseq_  = onlyc_nseq;
    nseq_        = nseq;
    nseg_        = nseg;
    nmes_        = (nmes_cx_ + nmes_cy_ + nmes_ib_);

    return true;
}

Bool_t TrFitPar::check_hits() {
    if (is_check_) return is_check_;
    is_check_ = sort_hits();
    return is_check_;
}


Bool_t TrFitPar::survival_test_and_modify(PhySt& part, Bool_t with_eloss) {
    if (Numc::EqualToZero(part.mom())) return false;
    Double_t finalZ = hits_.back()->cz();
   
    Bool_t   succ = false;
    Short_t  survival_iter = 0;
    Double_t survival_fact = 1.0;
    while (!succ && survival_iter <= SURVIVAL_LMTN) {
        PhySt ppst(part);
        ppst.arg().reset(false, with_eloss);
        ppst.set_eta(part.eta() * survival_fact);
        succ = (PropMgnt::PropToZ(finalZ, ppst) && ppst.bta() > SURVIVAL_BETA);
        if (!succ) survival_fact *= SURVIVAL_FACT;
        survival_iter++;
    }
    
    if (succ && Numc::Compare(survival_fact, Numc::ONE<>) < 0)
        part.set_eta(part.eta() * survival_fact);
    return succ;
}


} // namespace TrackSys


#endif // __TRACKLibs_TrFitPar_C__
