#ifndef __TRACKLibs_PhyFit_C__
#define __TRACKLibs_PhyFit_C__


#include "Sys.h"
#include "Math.h"
#include "IonEloss.h"
#include "GmIonEloss.h"
#include "PartInfo.h"
#include "PhySt.h"
#include "MagEnv.h"
#include "MatEnv.h"
#include "Prop.h"
#include "HitSt.h"
#include "PhyFit.h"


namespace TrackSys {
        
    
TrFitPar& TrFitPar::operator=(const TrFitPar& rhs) {
    if (this != &rhs) {
        sw_mscat_     = rhs.sw_mscat_;
        sw_eloss_     = rhs.sw_eloss_;
        info_         = rhs.info_;
        ortt_         = rhs.ortt_;
        hits_TRK_     = rhs.hits_TRK_;
        hits_TOF_     = rhs.hits_TOF_;
        hits_RICH_    = rhs.hits_RICH_;
        hits_TRD_     = rhs.hits_TRD_;
        onlyc_nseq_   = rhs.onlyc_nseq_;
        nseq_         = rhs.nseq_;
        nseg_         = rhs.nseg_;
        nmes_         = rhs.nmes_;
        nmes_cx_      = rhs.nmes_cx_;
        nmes_cy_      = rhs.nmes_cy_;
        nmes_ib_      = rhs.nmes_ib_;
        nmes_TRKqx_   = rhs.nmes_TRKqx_;
        nmes_TRKqy_   = rhs.nmes_TRKqy_;
        nmes_TOFt_    = rhs.nmes_TOFt_;
        nmes_TOFq_    = rhs.nmes_TOFq_;
        nmes_RICHib_  = rhs.nmes_RICHib_;
        nmes_TRDel_   = rhs.nmes_TRDel_;
        is_check_     = rhs.is_check_;
        
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
    }
   
    return *this;
}
    
TrFitPar::TrFitPar(const PartInfo& info, const Orientation& ortt, const Bool_t& sw_mscat, const Bool_t& sw_eloss) {
    clear();

    sw_mscat_ = sw_mscat;
    sw_eloss_ = sw_eloss;
    info_ = info;
    ortt_ = ortt;
}

void TrFitPar::zero() {
    hits_.clear();

    onlyc_nseq_  = 0;
    nseq_        = 0;
    nseg_        = 0;
    nmes_        = 0;
    nmes_cx_     = 0;
    nmes_cy_     = 0;
    nmes_ib_     = 0;
    nmes_TRKqx_  = 0;
    nmes_TRKqy_  = 0;
    nmes_TOFt_   = 0;
    nmes_TOFq_   = 0;
    nmes_RICHib_ = 0;
    nmes_TRDel_  = 0;
    
    is_check_ = false;
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
    
    Short_t onlyc_nseq =  0;
    Short_t nseq       =  0;
    Short_t nseg       = -1;
    for (auto&& hit : hits_) {
        hit->set_type(info_);
        onlyc_nseq += hit->set_onlyc_seqID(onlyc_nseq);
        nseq       += hit->set_seqID(nseq);
        
        if (hit->scx()) nmes_cx_++;
        if (hit->scy()) nmes_cy_++;
        if (hit->scx() || hit->scy()) nseg++;
    }

    for (auto&& hit : hits_TRK_) {
        if (hit.sqx()) nmes_TRKqx_++;
        if (hit.sqy()) nmes_TRKqy_++;
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
    
    nmes_ib_ = nmes_TRKqx_ + nmes_TRKqy_ + nmes_TOFt_ + nmes_TOFq_ + nmes_RICHib_ + nmes_TRDel_;
    
    onlyc_nseq_ = onlyc_nseq;
    nseq_       = nseq;
    nseg_       = nseg;
    nmes_       = (nmes_cx_ + nmes_cy_ + nmes_ib_);

    return true;
}

Bool_t TrFitPar::check_hits() {
    if (is_check_) return is_check_;
    sort_hits();
   
    Bool_t passed = (nmes_cx_ > LMTN_CX && nmes_cy_ > LMTN_CY);
    if (passed) is_check_ = true;
   
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
        ppst.arg().reset(sw_mscat_, with_eloss);
        ppst.set_eta(part.eta() * survival_fact);
        succ = (PropMgnt::PropToZ(finalZ, ppst) && ppst.bta() > SURVIVAL_BETA);
        if (!succ) survival_fact *= SURVIVAL_FACT;
        survival_iter++;
    }
    
    if (succ && Numc::Compare(survival_fact, Numc::ONE<>) < 0)
        part.set_eta(part.eta() * survival_fact);
    return succ;
}


SimpleTrFit& SimpleTrFit::operator=(const SimpleTrFit& rhs) {
    if (this != &rhs) {
        dynamic_cast<TrFitPar&>(*this) = dynamic_cast<const TrFitPar&>(rhs);
        succ_     = rhs.succ_;
        part_     = rhs.part_;
        TOFt_sft_ = rhs.TOFt_sft_;
        ndof_     = rhs.ndof_;
        ndof_cx_  = rhs.ndof_cx_;
        ndof_cy_  = rhs.ndof_cy_;
        ndof_ib_  = rhs.ndof_ib_;
        nchi_     = rhs.nchi_;
        nchi_cx_  = rhs.nchi_cx_;
        nchi_cy_  = rhs.nchi_cy_;
        nchi_ib_  = rhs.nchi_ib_;
    }
    return *this;
}


void SimpleTrFit::clear() {
    succ_ = false;
    part_.reset(info_);
    part_.arg().reset(sw_mscat_, sw_eloss_);

    TOFt_sft_ = 0;

    ndof_    = 0;
    ndof_cx_ = 0;
    ndof_cy_ = 0;
    ndof_ib_ = 0;
        
    nchi_    = 0;
    nchi_cx_ = 0;
    nchi_cy_ = 0;
    nchi_ib_ = 0;
}


SimpleTrFit::SimpleTrFit(const TrFitPar& fitPar, Bool_t advanced) : TrFitPar(fitPar) {
    SimpleTrFit::clear();
    if (!check_hits()) return;
    ndof_cx_ = (nmes_cx_ > LMTN_CX) ? (nmes_cx_ - LMTN_CX) : 0;
    ndof_cy_ = (nmes_cy_ > LMTN_CY) ? (nmes_cy_ - LMTN_CY) : 0;
    ndof_ib_ = nmes_ib_ - (nmes_TOFt_ >= LMTN_TOF_T);
    ndof_    = (ndof_cx_ + ndof_cy_ + ndof_ib_);
    if (ndof_cx_            <= Numc::ONE<Short_t>) { SimpleTrFit::clear(); return; }
    if ((ndof_cy_+ndof_ib_) <= Numc::ONE<Short_t>) { SimpleTrFit::clear(); return; }

    succ_ = (analyticalFit() ? simpleFit() : false);
    if (advanced && succ_) succ_ = advancedSimpleFit();
    if (!succ_) { SimpleTrFit::clear(); TrFitPar::clear(); }

    //if (!succ_) CERR("FAILURE === SimpleTrFit\n"); // testcode
}


Bool_t SimpleTrFit::analyticalFit() {
    // Global
    Double_t prefit_ortt = ((ortt_ == Orientation::kDownward) ? Numc::NEG<> : Numc::ONE<>);
    
    // Linear Fit on X
    // Equation of Motion
    // X  = PX + TX * (Z - RefZ)
    // dZ = Z - RefZ
    // UX = TX * UZ
    // Chisq = (X - Xm)^2
    // | PX |   | Sum(1)     Sum(dZ)   |^-1   | Sum(Xm)    |
    // |    | = |                      |    * |            |
    // | TX |   | Sum(dZ)    Sum(dZ^2) |      | Sum(dZ*Xm) |
    Double_t prefit_pz = hits_.at(0)->cz();
    Double_t prefit_px = Numc::ZERO<>;
    Double_t prefit_tx = Numc::ZERO<>;
    {
        SMtxSymD<2> mtx;
        SVecD<2>    res;
        for (auto&& hit : hits_) {
            if (!hit->scx()) continue;
            Double_t ex  = hit->ecx();
            Double_t err = (Numc::ONE<> / ex / ex);
            Double_t dz1 = hit->cz() - prefit_pz;
            Double_t dz2 = dz1 * dz1;
            mtx(0, 0) += err * Numc::ONE<>;
            mtx(0, 1) += err * dz1;
            mtx(1, 1) += err * dz2;
            res(0)    += err * hit->cx();
            res(1)    += err * dz1 * hit->cx();
        }
        if (!mtx.Invert()) return false;
        SVecD<2>&& rsl = mtx * res;
        prefit_px = rsl(0);
        prefit_tx = rsl(1);
    }
    

    // Curve Fit on Y
    // Equation of Motion
    // Y   = PY + UY * S + ETA * CRS * (0.5 * S^2)
    // U   = UY + ETA * CRS * S;
    // CRS = Lambda * (U x M) at mid-point
    // Yi  = PX + UY * Sum(Si) + ETA * (0.5 * SUM(CRSi*Si^2) + SUM(CRSj*Sj*Si)j<i)
    // Yi  = PX + UY * Au + ETA * Ae
    // Chisq = (Y - Ym)^2
    // | PY  |   | 1          Sum(Au)       Sum(Ae)    |^-1   | Sum(Ym)    |
    // |     |   |                                     |      |            |
    // | UY  | = | Sum(Au)    Sum(Au^2)     Sum(Au*Ae) |    * | Sum(Au*Ym) |
    // |     |   |                                     |      |            |
    // | ETA |   | Sum(Ae)    Sum(Au*Ae)    Sum(Ae^2)  |      | Sum(Ae*Ym) |
    Double_t prefit_py = Numc::ZERO<>;
    Double_t prefit_uy = Numc::ZERO<>;
    Double_t prefit_ea = 0.001;
    {
        const Double_t PROP_FACT = 2.99792458e-04;
        Double_t Lambda = PROP_FACT * part_.info().chrg_to_atomic_mass(); 
        
        std::vector<UInt_t> mapID;
        for (UInt_t ih = 0; ih < hits_.size(); ++ih)
            if (hits_.at(ih)->scy()) mapID.push_back(ih);

        std::vector<Double_t> stp(mapID.size(), Numc::ZERO<>);
        std::vector<Double_t> crs(mapID.size(), Numc::ZERO<>);

        const UInt_t nstp = 5;
        for (UInt_t id = 1; id < mapID.size(); ++id) {
            UInt_t ih = mapID.at(id-1);
            UInt_t jh = mapID.at(id);
            SVecD<3>&& ref_l = (hits_.at(jh)->c() - hits_.at(ih)->c());
            SVecD<3>&& ref_u = LA::Unit(ref_l);
            Double_t   ref_s = LA::Mag(ref_l);
            
            SVecD<3> mfldv;
            for (UInt_t it = 0; it < nstp; ++it) {
                Double_t stp = ((static_cast<Double_t>(it) + Numc::HALF) / static_cast<Double_t>(nstp));
                SVecD<3>&& refm = ((Numc::ONE<> - stp) * hits_.at(ih)->c() + stp * hits_.at(jh)->c());
                MagFld&&   mfld = MagMgnt::Get(refm);
                mfldv += mfld();
            }
            mfldv /= static_cast<Double_t>(nstp);
            Double_t mucrs = Lambda * (ref_u(2) * mfldv(0) - ref_u(0) * mfldv(2));

            stp.at(id) = ref_s;
            crs.at(id) = mucrs;
        }
        
        SMtxSymD<3> mtx;
        SVecD<3>    res;
        Double_t    cur_Au = Numc::ZERO<>;
        Double_t    cur_Ae = Numc::ZERO<>;
        for (UInt_t id = 0; id < mapID.size(); ++id) {
            VirtualHitSt* hit = hits_.at(mapID.at(id));
            Double_t ey  = hit->ecy();
            Double_t err = (Numc::ONE<> / ey / ey);
            
            cur_Au += stp.at(id);
            cur_Ae += Numc::HALF * crs.at(id) * stp.at(id) * stp.at(id);
            for (UInt_t jd = 0; jd < id; ++jd)
                cur_Ae += crs.at(jd) * stp.at(jd) * stp.at(id);

            mtx(0, 0) += err;
            mtx(0, 1) += err * cur_Au;
            mtx(0, 2) += err * cur_Ae;
            mtx(1, 1) += err * cur_Au * cur_Au;
            mtx(1, 2) += err * cur_Au * cur_Ae;
            mtx(2, 2) += err * cur_Ae * cur_Ae;
            res(0)    += err * hit->cy();
            res(1)    += err * cur_Au * hit->cy();
            res(2)    += err * cur_Ae * hit->cy();
        }

        // Scale
        Double_t scl = LA::Mag(res);
        if (Numc::EqualToZero(scl)) scl = Numc::ONE<>;
        mtx /= scl;
        res /= scl;

        Bool_t no_mag = (
                Numc::EqualToZero(res(2)) && 
                Numc::EqualToZero(mtx(0, 2)) && 
                Numc::EqualToZero(mtx(1, 2)) && 
                Numc::EqualToZero(mtx(2, 2)));

        if (no_mag) {
            SVecD<2> res_nomag(res(0), res(1));
            SMtxSymD<2>&& mtx_nomag = mtx.Sub<SMtxSymD<2>>(0, 0);
            if (!mtx_nomag.Invert()) return false;
            SVecD<2>&& rsl_nomag = mtx_nomag * res_nomag;
            prefit_py = rsl_nomag(0);
            prefit_uy = rsl_nomag(1);
        }
        else {
            if (!mtx.Invert()) return false;
            SVecD<3>&& rsl = mtx * res;
            prefit_py = rsl(0);
            prefit_uy = rsl(1);
            prefit_ea = rsl(2);
        }
        
        // transport to first hit
        UInt_t hitID = mapID.at(0);
        if (hitID != 0) {
            VirtualHitSt* hit = hits_.at(hitID);
            Double_t dz = (prefit_pz - hit->cz());
            Double_t ty = prefit_ortt * prefit_uy / std::sqrt(Numc::ONE<> - prefit_uy * prefit_uy);
            if (!Numc::Valid(ty)) ty = Numc::ZERO<>;
            prefit_py = prefit_py + (ty * dz); 
        }
    }
   
    // Merge Fitting Result
    Double_t prefit_uz = prefit_ortt * std::fabs((Numc::ONE<> - prefit_uy * prefit_uy) / (Numc::ONE<> + prefit_tx * prefit_tx));
    Double_t prefit_ux = prefit_tx * prefit_uz;
    part_.set_state_with_cos(
        prefit_px, prefit_py, prefit_pz,
        prefit_ux, prefit_uy, prefit_uz
    );
    part_.set_eta(prefit_ea);

    // Eta Limit Setting
    if (part_.eta_abs() > LMTU_ETA)
        part_.set_eta(part_.eta_sign() * LMTU_ETA);
    survival_test_and_modify(part_, false);
    
    return true;
}


Bool_t SimpleTrFit::simpleFit() {
    // Turn Off (mscat, eloss)
    part_.arg().reset(false, false);
    
    Bool_t succ    = false;
    Bool_t preSucc = false;
    Bool_t curSucc = false;

    Double_t       curLmRhoDen = Numc::ONE<>;
    Double_t       lambda = LAMBDA0;
    PhySt          rltSt(part_);
    SVecD<DIMG>    curGrdG;
    SMtxSymD<DIMG> curCvGG;

    Short_t updIter = 0;
    Short_t curIter = 0;
    while (curIter <= LMTU_ITER && !succ) {
        Bool_t resetTOF = true;
        HitStTOF::SetOffsetTime(Numc::ZERO<>);
        HitStTOF::SetOffsetPath(Numc::ZERO<>);
        HitStTOF::SetTimeShiftCorr(false);
        
        Double_t chi_cx = 0;
        Double_t chi_cy = 0;
        Double_t chi_ib = 0;
        
        SVecD<DIMG>    grdG;
        SMtxSymD<DIMG> cvGG;

        Int_t cnt_nhit = 0;
        PhySt ppst(rltSt);
        SMtxD<DIMG>&& ppjb = SMtxId();
        for (auto&& hit : hits_) {
            PhyJb curjb;
            if (!PropMgnt::PropToZ(hit->cz(), ppst, nullptr, &curjb)) break;
            ppjb = curjb.gg() * ppjb;
        
            // Hit Status: Setting TOF reference time and path
            if (resetTOF && Hit<HitStTOF>::IsSame(hit)) { // set reference
                HitStTOF::SetOffsetPath(ppst.path());
                HitStTOF::SetOffsetTime(ppst.time()-Hit<HitStTOF>::Cast(hit)->orgt());
                resetTOF = false;
            }
            hit->cal(ppst);

            SVecD<2>       rsC;
            SMtxD<2, DIMG> jbC;
            if (hit->scx()) rsC(0) = (hit->cx() - ppst.cx()) / hit->ecx();
            if (hit->scy()) rsC(1) = (hit->cy() - ppst.cy()) / hit->ecy();
            for (Short_t it = 0; it < DIMG; ++it) {
                if (hit->scx()) jbC(0, it) += (Numc::NEG<> / hit->ecx()) * ppjb(0, it);
                if (hit->scy()) jbC(1, it) += (Numc::NEG<> / hit->ecy()) * ppjb(1, it);
            }
           
            grdG += LA::Transpose(jbC) * rsC;
            cvGG += LA::SimilarityT(jbC, SMtxSymD<2>(SMtxId()));
            if (hit->scx()) chi_cx += rsC(0) * rsC(0);
            if (hit->scy()) chi_cy += rsC(1) * rsC(1);
            
            // TRK
            HitStTRK* hitTRK = Hit<HitStTRK>::Cast(hit);
            if (hitTRK != nullptr) {
                SVecD<2> rsTRK;
                SVecD<2> jbTRK;

                if (hitTRK->sqx()) rsTRK(0) = hitTRK->nrmqx();
                if (hitTRK->sqy()) rsTRK(1) = hitTRK->nrmqy();
                    
                if (hitTRK->sqx()) jbTRK(0) = hitTRK->divqx_eta() * ppjb(4, 4);
                if (hitTRK->sqy()) jbTRK(1) = hitTRK->divqy_eta() * ppjb(4, 4);
                
                grdG(4)    += (jbTRK(0) * rsTRK(0) + jbTRK(1) * rsTRK(1));
                cvGG(4, 4) += (jbTRK(0) * jbTRK(0) + jbTRK(1) * jbTRK(1));
                if (hitTRK->sqx()) chi_ib += hitTRK->chiqx() * hitTRK->chiqx();
                if (hitTRK->sqy()) chi_ib += hitTRK->chiqy() * hitTRK->chiqy();
            }

            // TOF
            HitStTOF* hitTOF = Hit<HitStTOF>::Cast(hit);
            if (hitTOF != nullptr) {
                SVecD<2> rsTOF;
                SVecD<2> jbTOF;

                if (hitTOF->st()) rsTOF(0) = hitTOF->nrmt();
                if (hitTOF->sq()) rsTOF(1) = hitTOF->nrmq();
                
                if (hitTOF->st()) jbTOF(0) = hitTOF->divt_eta() * ppjb(4, 4);
                if (hitTOF->sq()) jbTOF(1) = hitTOF->divq_eta() * ppjb(4, 4);
                
                grdG(4)    += (jbTOF(0) * rsTOF(0) + jbTOF(1) * rsTOF(1));
                cvGG(4, 4) += (jbTOF(0) * jbTOF(0) + jbTOF(1) * jbTOF(1));
                if (hitTOF->st()) chi_ib += hitTOF->chit() * hitTOF->chit();
                if (hitTOF->sq()) chi_ib += hitTOF->chiq() * hitTOF->chiq();
            }
            
            // RICH
            HitStRICH* hitRICH = Hit<HitStRICH>::Cast(hit);
            if (hitRICH != nullptr) {
                Double_t rsRICH = Numc::ZERO<>;
                Double_t jbRICH = Numc::ZERO<>;

                if (hitRICH->sib()) rsRICH = hitRICH->nrmib();
                if (hitRICH->sib()) jbRICH = hitRICH->divib_eta() * ppjb(4, 4);

                grdG(4)    += (jbRICH * rsRICH);
                cvGG(4, 4) += (jbRICH * jbRICH);

                if (hitRICH->sib()) chi_ib += hitRICH->chiib() * hitRICH->chiib();
            }

            cnt_nhit++;
        }
        if (cnt_nhit != hits_.size()) break;
        Double_t chi  = (chi_cx + chi_cy + chi_ib);
        Double_t nchi = (chi / static_cast<Double_t>(ndof_));

        Bool_t isSucc   = false;
        Bool_t isUpdate = false;
        if (curIter != 0) {
            Double_t lmRho = (nchi_ - nchi) / curLmRhoDen;
            Bool_t   isLmt = (Numc::Compare(lambda, LMTU_LAMBDA) >= 0);
            Double_t convg = std::sqrt(Numc::ONE<> + lambda);
            isSucc = (Numc::Compare(std::fabs((nchi_ - nchi) / (nchi_ + nchi + CONVG_TOLERANCE)) * convg, CONVG_TOLERANCE) <= 0);

            if (Numc::Compare(lmRho, CONVG_EPSILON) < 0) {
                lambda = std::min(lambda*LAMBDA_UP_FAC, LMTU_LAMBDA); 
                grdG   = curGrdG;
                cvGG   = curCvGG;
                rltSt  = part_;
                if      (isSucc) updIter++;
                else if (isLmt)  break;
            }
            else {
                lambda   = std::max(lambda/LAMBDA_DN_FAC, LMTL_LAMBDA);
                nchi_cx_ = ((ndof_cx_ > 0) ? (chi_cx / static_cast<Double_t>(ndof_cx_)) : 0);
                nchi_cy_ = ((ndof_cy_ > 0) ? (chi_cy / static_cast<Double_t>(ndof_cy_)) : 0);
                nchi_ib_ = ((ndof_ib_ > 0) ? (chi_ib / static_cast<Double_t>(ndof_ib_)) : 0);
                nchi_      = nchi;
                part_      = rltSt;
                isUpdate   = true;
                updIter++;
            }
        }
        else { nchi_ = nchi; }

        SMtxSymD<DIMG> lmCvGG(cvGG);
        SVecD<DIMG>&&  diagCvGG = (lambda * cvGG.Diagonal());
        lmCvGG.SetDiagonal(SVecD<DIMG>(lmCvGG.Diagonal() + diagCvGG));

        Bool_t isNomag = (Numc::EqualToZero(lmCvGG(4, 4)) && Numc::EqualToZero(grdG(4))); // Fast Check
        if (isNomag) {
            grdG(4) = Numc::ZERO<>;
            SMtxSymD<4>&& lmCvGG_nomag = lmCvGG.Sub<SMtxSymD<4>>(0, 0);
            if (!lmCvGG_nomag.Invert()) break;
            lmCvGG = std::move(SMtxSymD<DIMG>());
            for (Short_t ielem = 0; ielem < 4; ++ielem)
                for (Short_t jelem = ielem; jelem < 4; ++jelem)
                    lmCvGG(ielem, jelem) = lmCvGG_nomag(ielem, jelem);
        }
        else {
            if (!lmCvGG.Invert()) break;
        }

        SVecD<DIMG>&& rslG = (lmCvGG * grdG);
       
        curLmRhoDen = Numc::ZERO<>;
        for (Short_t it = 0; it < DIMG; ++it)
            curLmRhoDen += (rslG(it) * (diagCvGG(it)*rslG(it) + grdG(it)));
        
        if (curIter == 0 || isUpdate) {
            curGrdG = grdG;
            curCvGG = cvGG;
        }

        rltSt.set_state_with_uxy(
            rltSt.cx() - rslG(0),
            rltSt.cy() - rslG(1),
            rltSt.cz(),
            rltSt.ux() - rslG(2),
            rltSt.uy() - rslG(3),
            ((ortt_ == Orientation::kDownward) ? Numc::NEG<Short_t> : Numc::ONE<Short_t>)
        );
        rltSt.set_eta(rltSt.eta() - rslG(4));
        
        preSucc = curSucc;
        curSucc = (isSucc && updIter >= LMTL_ITER);
        succ    = (preSucc && curSucc);
        
        if (!succ) curIter++;
    }
    
    // Turn Back (mscat, eloss)
    part_.arg().reset(sw_mscat_, sw_eloss_);
   
    //if (!succ) CERR("FAIL. IT %2d %2d (RIG %14.8f MASS %14.8f CHI %14.8f) LAMBDA %14.8f\n", curIter, updIter, part_.rig(), part_.mass(), nchi_, lambda);
    //else       CERR("SUCC. IT %2d %2d (RIG %14.8f MASS %14.8f CHI %14.8f) LAMBDA %14.8f\n", curIter, updIter, part_.rig(), part_.mass(), nchi_, lambda);
    
    return succ;
}


Bool_t SimpleTrFit::advancedSimpleFit() {
    if (!survival_test_and_modify(part_, sw_eloss_)) return false;
    Bool_t  opt_tsft = (nmes_TOFt_ >= LMTN_TOF_T);
    
    std::vector<double> params({ part_.cx(), part_.cy(), part_.ux(), part_.uy(), part_.eta() });
    if (opt_tsft) { TOFt_sft_ = Numc::ZERO<>; params.push_back(TOFt_sft_); } // TOF Shift Time
    Short_t parIDtsft = (opt_tsft ? (params.size() - 1) : -1);

    // CeresSolver: Problem
    ceres::GradientProblem problem(new VirtualSimpleTrFit(dynamic_cast<TrFitPar&>(*this), part_));
    
    // CeresSolver: Options
    ceres::GradientProblemSolver::Options options;
    options.line_search_direction_type = ceres::LBFGS;
    options.use_approximate_eigenvalue_bfgs_scaling = true;
    options.nonlinear_conjugate_gradient_type = ceres::HESTENES_STIEFEL;

    // CeresSolver: Summary
    ceres::GradientProblemSolver::Summary summary;
    ceres::Solve(options, problem, params.data(), &summary);
    if (!summary.IsSolutionUsable()) return false;

    Double_t partcz = part_.cz();
    Short_t  signuz = Numc::Compare(part_.uz());
    part_.set_state_with_uxy(params.at(0), params.at(1), partcz, params.at(2), params.at(3), signuz);
    part_.set_eta(params.at(4));
    if (opt_tsft) TOFt_sft_ = params.at(parIDtsft);
    
    return true;
}


bool VirtualSimpleTrFit::Evaluate(const double* parameters, double* cost, double* gradient) const {
    if (numOfRes_ <= 0 || numOfPar_ <= 0) return false;
    cost[0] = Numc::ZERO<>;
    
    Bool_t hasGrad = (gradient != nullptr);
    if (hasGrad) std::fill_n(gradient, numOfPar_, Numc::ZERO<>);
    
    // Reset TOF Time and Path
    Bool_t resetTOF = true;
    HitStTOF::SetOffsetTime(Numc::ZERO<>);
    HitStTOF::SetOffsetPath(Numc::ZERO<>);
    HitStTOF::SetTimeShiftCorr(opt_tsft_);

    // TOF time shift
    Double_t TOFt_sft = (opt_tsft_ ? parameters[parIDtsft_] : Numc::ZERO<>);

    // Particle Status
    PhySt ppst(part_);
    ppst.arg().clear();
    ppst.set_state_with_uxy(parameters[0], parameters[1], part_.cz(), parameters[2], parameters[3], Numc::Compare(part_.uz()));
    ppst.set_eta(parameters[4]);
    
    // Matrix (Rs, Jb)
    PhyJb::SMtxDGG&& jbGG = SMtxId();
    std::vector<PhyJb::SMtxDGL> jbGL(nseg_);
    
    ceres::Vector&& rsC  = ceres::Vector::Zero(onlyc_nseq_);
    ceres::Matrix&& jbCG = ceres::Matrix::Zero(onlyc_nseq_, numOfPar_);     // Jacb for coord
    ceres::Matrix&& jbCW = ceres::Matrix::Zero(onlyc_nseq_, nseg_ * DIML_); // Jacb for interaction 
    ceres::Matrix&& cvCC = ceres::Matrix::Zero(onlyc_nseq_, onlyc_nseq_); 
    
    Double_t        costG = Numc::ZERO<>;
    ceres::Vector&& grdG  = ceres::Vector::Zero(numOfPar_);
    
    Short_t cnt_nhit =  0;
    Short_t cnt_nseg = -1;
    PhySt                       nearPpst = ppst;
    PhyJb::SMtxDGG              nearJbGG = jbGG;
    std::vector<PhyJb::SMtxDGL> nearJbGL = jbGL;
    for (auto&& hit : hits_) {
        // Interaction Local Parameters
        Bool_t  hasLoc  = (cnt_nseg >= 0);
        Bool_t  isInner = (cnt_nseg >= 0 && cnt_nseg < nseg_);
        Short_t itnseg  = (cnt_nseg == nseg_) ? (nseg_ - Numc::ONE<Short_t>) : cnt_nseg;
        Bool_t  hasCxy  = (hit->scx() || hit->scy());

        if (isInner && !hasCxy) {
            ppst = nearPpst;
            jbGG = nearJbGG;
            jbGL = nearJbGL;
        }
       
        // Propagate
        PhyJb curjb;
        if (!PropMgnt::PropToZ(hit->cz(), ppst, nullptr, &curjb)) break;
        ppst.symbk();
            
        // Hit Status: Setting TOF reference time and path
        if (resetTOF && Hit<HitStTOF>::IsSame(hit)) { // set reference
            HitStTOF::SetOffsetPath(ppst.path());
            HitStTOF::SetOffsetTime(ppst.time()-Hit<HitStTOF>::Cast(hit)->orgt());
            resetTOF = false;
        }
        hit->cal(ppst);
        
        // Update Jacb
        jbGG = curjb.gg() * jbGG;
        if (opt_int_) {
            if (hasLoc && hasCxy) jbGL.at(cnt_nseg) = curjb.gl();
            for (Short_t is = 0; is < cnt_nseg; ++is)
                jbGL.at(is) = curjb.gg() * jbGL.at(is);
        }
        
        // Coord
        if (hit->scx()) rsC(hit->onlyc_seqIDcx()) += (hit->cx() - ppst.cx());
        if (hit->scy()) rsC(hit->onlyc_seqIDcy()) += (hit->cy() - ppst.cy());
        if (hit->scx()) cvCC(hit->onlyc_seqIDcx(), hit->onlyc_seqIDcx()) += (hit->ecx() * hit->ecx());
        if (hit->scy()) cvCC(hit->onlyc_seqIDcy(), hit->onlyc_seqIDcy()) += (hit->ecy() * hit->ecy());
        if (hasCxy) {
            if (hasGrad) {
                for (Short_t it = 0; it < PhyJb::DIMG; ++it) {
                    if (hit->scx()) jbCG(hit->onlyc_seqIDcx(), it) += Numc::NEG<> * jbGG(0, it);
                    if (hit->scy()) jbCG(hit->onlyc_seqIDcy(), it) += Numc::NEG<> * jbGG(1, it);
                }
            } // hasGrad
            if (opt_int_) {
                for (Short_t is = 0; is <= itnseg; ++is) {
                    for (Short_t it = 0; it < PhyJb::DIML; ++it) {
                        if (hit->scx()) jbCW(hit->onlyc_seqIDcx(), is*DIML_+it) += Numc::NEG<> * jbGL.at(is)(0, it);
                        if (hit->scy()) jbCW(hit->onlyc_seqIDcy(), is*DIML_+it) += Numc::NEG<> * jbGL.at(is)(1, it);
                    }
                }
            } // hasGrad
        } // hasCxy
            
        // TRK
        HitStTRK* hitTRK = Hit<HitStTRK>::Cast(hit);
        if (hitTRK != nullptr) {
            if (hitTRK->sqx()) costG += hitTRK->chiqx() * hitTRK->chiqx();
            if (hitTRK->sqy()) costG += hitTRK->chiqy() * hitTRK->chiqy();
                
            if (hitTRK->sqx()) grdG(4) += (hitTRK->divqx_eta() * jbGG(4, 4)) * hitTRK->nrmqx();
            if (hitTRK->sqy()) grdG(4) += (hitTRK->divqy_eta() * jbGG(4, 4)) * hitTRK->nrmqy();
        }
        
        // TOF
        HitStTOF* hitTOF = Hit<HitStTOF>::Cast(hit);
        if (hitTOF != nullptr) {
            if (hitTOF->st()) costG += hitTOF->chit() * hitTOF->chit();
            if (hitTOF->sq()) costG += hitTOF->chiq() * hitTOF->chiq();
            
            if (hitTOF->st()) grdG(4) += (hitTOF->divt_eta() * jbGG(4, 4)) * hitTOF->nrmt();
            if (hitTOF->sq()) grdG(4) += (hitTOF->divq_eta() * jbGG(4, 4)) * hitTOF->nrmq();
            
            if (opt_tsft_ && hitTOF->st()) grdG(parIDtsft_) += hitTOF->divt_sft(); // TOF time shift
        }
        
        // RICH
        HitStRICH* hitRICH = Hit<HitStRICH>::Cast(hit);
        if (hitRICH != nullptr) {
            if (hitRICH->sib()) costG += hitRICH->chiib() * hitRICH->chiib();
            if (hitRICH->sib()) grdG(4) += (hitRICH->divib_eta() * jbGG(4, 4)) * hitRICH->nrmib();
        }
        
        if (hasCxy) {
            nearPpst = ppst;
            nearJbGG = jbGG;
            nearJbGL = jbGL;
            cnt_nseg++;
        }
        cnt_nhit++;
    }
    if (cnt_nhit != hits_.size()) return false;
    if (cnt_nseg != nseg_) return false;
    
    cvCC += (jbCW * jbCW.transpose());
    ceres::Matrix&& wtCC = cvCC.inverse();

    Double_t        costC = Numc::ZERO<>;
    ceres::Vector&& grdC  = ceres::Vector::Zero(numOfPar_);
    costC += (rsC.transpose() * wtCC * rsC);
    if (hasGrad) grdC += (jbCG.transpose() * wtCC * rsC);

    Double_t chiC = std::sqrt(costC / static_cast<Double_t>(nmes_cx_ + nmes_cy_ - Numc::THREE<Short_t>));
    Robust robust(Robust::Opt::ON, Numc::FOUR<long double>);
    std::array<long double, 3>&& mini = robust.minimizer(chiC);
    costC *= (mini.at(0) * mini.at(0));
    if (hasGrad) grdC *= (mini.at(1) * mini.at(2));

    costG += costC;
    if (hasGrad) grdG += grdC;

    Short_t ndof = (numOfRes_ - numOfPar_);
    cost[0] = (costG / static_cast<Double_t>(ndof));
    if (hasGrad) for (int it = 0; it < numOfPar_; ++it) gradient[it] = (grdG(it) / static_cast<Double_t>(ndof));
    
    return true;
}


PhyTrFit& PhyTrFit::operator=(const PhyTrFit& rhs) {
    if (this != &rhs) {
        dynamic_cast<TrFitPar&>(*this) = dynamic_cast<const TrFitPar&>(rhs);
    
        mu_opt_ = rhs.mu_opt_;
        succ_   = rhs.succ_;
        part_   = rhs.part_;
        args_   = rhs.args_;
        stts_   = rhs.stts_;

        TOFt_sft_ = rhs.TOFt_sft_;
      
        ndof_    = rhs.ndof_;
        nchi_    = rhs.nchi_;
        quality_ = rhs.quality_;

        ndof_tt_ = rhs.ndof_tt_;
        ndof_cx_ = rhs.ndof_cx_;
        ndof_cy_ = rhs.ndof_cy_;
        ndof_ib_ = rhs.ndof_ib_;
        nchi_tt_ = rhs.nchi_tt_;
        nchi_cx_ = rhs.nchi_cx_;
        nchi_cy_ = rhs.nchi_cy_;
        nchi_ib_ = rhs.nchi_ib_;

        err_ = rhs.err_;
    }
    return *this;
}


void PhyTrFit::clear() {
    mu_opt_ = MuOpt::kFixed;
    succ_   = false;
    part_.reset(info_);
    part_.arg().reset(sw_mscat_, sw_eloss_);
    args_.clear();
    stts_.clear();

    TOFt_sft_ = 0;

    ndof_.at(0) = 0;
    ndof_.at(1) = 0;
    nchi_.at(0) = 0;
    nchi_.at(1) = 0;
    quality_.at(0) = 0;
    quality_.at(1) = 0;

    ndof_tt_ = 0;
    ndof_cx_ = 0;
    ndof_cy_ = 0;
    ndof_ib_ = 0;
        
    nchi_tt_ = 0;
    nchi_cx_ = 0;
    nchi_cy_ = 0;
    nchi_ib_ = 0;

    err_.fill(Numc::ZERO<>);
}


PhyTrFit::PhyTrFit(const TrFitPar& fitPar, const MuOpt& mu_opt) : TrFitPar(fitPar) {
    PhyTrFit::clear();
    mu_opt_ = mu_opt;
    if (!check_hits()) return;
    ndof_cx_ = (nmes_cx_ > LMTN_CX) ? (nmes_cx_ - LMTN_CX) : 0;
    ndof_cy_ = (nmes_cy_ > LMTN_CY) ? (nmes_cy_ - LMTN_CY) : 0;
    ndof_ib_ = nmes_ib_ - (nmes_TOFt_ >= LMTN_TOF_T);
    ndof_tt_ = (ndof_cx_ + ndof_cy_ + ndof_ib_);
    if (MuOpt::kFree == mu_opt_) { ndof_ib_ -= 1; ndof_tt_ -= 1; }
    if (ndof_cx_            <= Numc::ONE<Short_t>) { PhyTrFit::clear(); return; }
    if ((ndof_cy_+ndof_ib_) <= Numc::ONE<Short_t>) { PhyTrFit::clear(); return; }
    if (MuOpt::kFree == mu_opt_) {
        if (ndof_cy_ <= Numc::ONE<Short_t>) { PhyTrFit::clear(); return; }
        if (ndof_ib_ <= Numc::ONE<Short_t>) { PhyTrFit::clear(); return; }
    }
    ndof_.at(0) = ndof_cx_;
    ndof_.at(1) = ndof_cy_ + ndof_ib_;

    // Fitting
    if (MuOpt::kFixed == mu_opt_) succ_ = physicalTrFit();
    else                          succ_ = physicalMuFit();
    
    if (!succ_) { PhyTrFit::clear(); TrFitPar::clear(); }
    
    //if (!succ_) CERR("FAILURE === PhyTrFit\n"); // testcode
}


Bool_t PhyTrFit::simpleFit() {
    SimpleTrFit simple(dynamic_cast<TrFitPar&>(*this), true); // testcode
    if (simple.status()) part_ = simple.part();
    part_.arg().reset(sw_mscat_, sw_eloss_);
    Bool_t succ = (simple.status() && survivalTestAndModify());
    return succ;
}


Bool_t PhyTrFit::physicalFit(const MuOpt& mu_opt, Double_t fluc_eta, Double_t fluc_igb, Bool_t with_mu_est) {
    Bool_t opt_mu   = (MuOpt::kFree == mu_opt);
    Bool_t opt_int  = sw_mscat_;
    Bool_t opt_tsft = (nmes_TOFt_ >= LMTN_TOF_T);
    
    Short_t DIMG = PhyJb::DIMG + opt_mu + opt_tsft;
    Short_t DIML = PhyJb::DIML;

    Double_t eta = part_.eta();
    const Bool_t is_fluc_eta = (MuOpt::kFree == mu_opt && Numc::Compare(fluc_eta) > 0);
    if (is_fluc_eta) {
        const Int_t niter = 5; Int_t iter = 0;
        do {
            Double_t rndm = Rndm::NormalGaussian() + MU_FLUC_BASE * Numc::TWO<> * (Rndm::DecimalUniform() - Numc::HALF);
            if (std::fabs(rndm) > Numc::TWO<>) { iter++; continue; }
            eta = part_.eta() * (Numc::ONE<> + fluc_eta * rndm);
            iter++;
        } while ((iter < niter) && Numc::EqualToZero(eta));
        if (iter >= niter) eta = part_.eta();
    }
    
    Double_t igb = part_.igmbta();
    const Bool_t is_fluc_igb = (MuOpt::kFree == mu_opt && Numc::Compare(fluc_igb) > 0);
    if (is_fluc_igb) {
        const Int_t niter = 5; Int_t iter = 0;
        do {
            Double_t rndm = Rndm::NormalGaussian() + MU_FLUC_BASE * Numc::TWO<> * (Rndm::DecimalUniform() - Numc::HALF);
            if (std::fabs(rndm) > Numc::TWO<>) { iter++; continue; }
            igb = part_.igmbta() * (Numc::ONE<> + fluc_igb * rndm);
            iter++;
        } while ((iter < niter) && (igb < LMTL_INV_GB || igb > LMTU_INV_GB));
        if (iter >= niter) igb = part_.igmbta();
    }

    // Gobal Parameters
    TOFt_sft_ = Numc::ZERO<>;
    std::vector<double> params_glb({ part_.cx(), part_.cy(), part_.ux(), part_.uy(), eta });

    Short_t parIDeta  =  4;
    Short_t parIDigb  = -1;
    Short_t parIDtsft = -1;
    if (opt_mu)   { parIDigb  = 5;        params_glb.push_back(igb); }
    if (opt_tsft) { parIDtsft = DIMG - 1; params_glb.push_back(TOFt_sft_); } // TOF Shift Time  

    // Interaction Parameters
    args_ = std::move(std::vector<PhyArg>(nseg_, PhyArg(sw_mscat_, sw_eloss_)));
    std::vector<double> params_int(nseg_*DIML, Numc::ZERO<>);

    // CeresSolver: Cost Function
    ceres::CostFunction* cost_function = new VirtualPhyTrFit(dynamic_cast<TrFitPar&>(*this), part_, opt_mu);

    // CeresSolver: Problem
    ceres::Problem problem;
    if (opt_int) problem.AddResidualBlock(cost_function, nullptr, params_glb.data(), params_int.data());
    else         problem.AddResidualBlock(cost_function, nullptr, params_glb.data());
    problem.SetParameterLowerBound(params_glb.data(), 2, -1.0);
    problem.SetParameterUpperBound(params_glb.data(), 2,  1.0);
    problem.SetParameterLowerBound(params_glb.data(), 3, -1.0);
    problem.SetParameterUpperBound(params_glb.data(), 3,  1.0);
    if (opt_mu) problem.SetParameterLowerBound(params_glb.data(), parIDigb, LMTL_INV_GB);
    if (opt_mu) problem.SetParameterUpperBound(params_glb.data(), parIDigb, LMTU_INV_GB);
    
    // CeresSolver: Options
    ceres::Solver::Options options;
    //options.line_search_direction_type = ceres::LBFGS;
    //options.check_gradients = true;

    // CeresSolver: Summary
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    if (!summary.IsSolutionUsable()) return false;
   
    // Result (Global)
    Double_t partcz = part_.cz();
    Short_t  signuz = Numc::Compare(part_.uz());
    Double_t partmu = (opt_mu ? std::fabs(params_glb.at(parIDigb) / params_glb.at(parIDeta)) : part_.mu());
    if (opt_mu) part_.reset(partmu);
    else        part_.arg().clear();
    part_.set_state_with_uxy(params_glb.at(0), params_glb.at(1), partcz, params_glb.at(2), params_glb.at(3), signuz);
    part_.set_eta(params_glb.at(parIDeta));
    
    if (opt_tsft) { TOFt_sft_ = params_glb.at(parIDtsft); } // TOF Shift Time

    // Result (Interaction)
    if (opt_int) {
    for (Short_t is = 0; is < nseg_; ++is) {
        args_.at(is).set_mscat(
            params_int.at(is*DIML+0), 
            params_int.at(is*DIML+1), 
            params_int.at(is*DIML+2), 
            params_int.at(is*DIML+3));
    }}

    // Result (Mu Fitting)
    //std::cerr << summary.FullReport() << std::endl;
    Bool_t   rw_err_mu = (MuOpt::kFixed == mu_opt && with_mu_est && evolve(MuOpt::kFree) && Numc::Valid(err_.at(6)) && Numc::Compare(err_.at(6)) > 0);
    Double_t err_mu    = (rw_err_mu ? err_.at(6) : Numc::ZERO<>);
    
    Bool_t succ = evolve(mu_opt);
    if (rw_err_mu) err_.at(6) = err_mu;
    return succ;
}


Bool_t PhyTrFit::physicalTrFit() {
    if (!simpleFit()) return false;
    return true; // testcode

    if (sw_mscat_) {
        Bool_t sw_mscat = sw_mscat_;
        Bool_t sw_eloss = sw_eloss_;
        resetPhyArg(false, sw_eloss_);
        Bool_t presucc = physicalFit(MuOpt::kFixed);
        resetPhyArg(sw_mscat, sw_eloss);
        if (!presucc) return false;
    }

    Bool_t succ = physicalFit(MuOpt::kFixed);
    return succ;
}


Bool_t PhyTrFit::physicalMuFit() {
    Short_t chrg = std::abs(info_.chrg());
    if (chrg <= Numc::ZERO<Short_t> || chrg >= PartListMassQ.size()) return false;
    if (PartListMassQ.at(chrg).size() == 0) return false;

    class PartElem {
        public :
            PartElem() : succ(false), fluc_eta(0), fluc_igb(0), qlt(0) {}
            PartElem(const PhySt& _part, Double_t _fluc_eta, Double_t _fluc_igb, Double_t _qlt) : succ(true), part(_part), fluc_eta(_fluc_eta), fluc_igb(_fluc_igb), qlt(_qlt) {}
            Bool_t              succ;
            PhySt               part;
            Double_t            fluc_eta;
            Double_t            fluc_igb;
            Double_t            qlt;
    };
   
    // List of Particle Mass (Init)
    PartElem condElem;
    for (auto&& mass : PartListMassQ.at(chrg)) {
        info_.reset(chrg, mass);
        if (!(simpleFit() ? physicalFit(MuOpt::kFixed, Numc::ZERO<>, Numc::ZERO<>, false) : false)) continue;
        
        Bool_t firstTime = (!condElem.succ);
        if (!firstTime && Numc::Compare(quality_.at(1), condElem.qlt) > 0) continue;
        if (!evolve(MuOpt::kFree)) continue;
        
        Double_t fluc_eta = std::fabs(err_.at(4) / part_.eta());
        Double_t fluc_igb = (err_.at(5) / part_.igmbta());
        condElem = std::move(PartElem(part_, fluc_eta, fluc_igb, quality_.at(1)));
    }
    if (!condElem.succ) return false;

    info_ = condElem.part.info();
    part_ = condElem.part;
    Double_t fluc_eta = condElem.fluc_eta;
    Double_t fluc_igb = condElem.fluc_igb;
    for (Short_t iter = 1; iter <= LMT_MU_ITER; ++iter) {
        Double_t wgt_eta = fluc_eta * fluc_eta;
        Double_t wgt_igb = fluc_igb * fluc_igb;
        Double_t rat_eta = wgt_eta / (wgt_eta + wgt_igb);
        Double_t rat_igb = wgt_igb / (wgt_eta + wgt_igb);

        if (!physicalFit(MuOpt::kFree,
                         MU_FLUC * rat_eta * fluc_eta,
                         MU_FLUC * rat_igb * fluc_igb)) return false;
        fluc_eta = std::fabs(err_.at(4) / part_.eta());
        fluc_igb = (err_.at(5) / part_.igmbta());
    }

    return true;
}


Bool_t PhyTrFit::evolve(const MuOpt& mu_opt) {
    Bool_t opt_mu   = (MuOpt::kFree == mu_opt);
    Bool_t opt_int  = (sw_mscat_);
    Bool_t opt_tsft = (nmes_TOFt_ >= LMTN_TOF_T);
    
    Short_t DIMG = PhyJb::DIMG + opt_mu + opt_tsft;
    Short_t DIML = PhyJb::DIML;
   
    // Jacb
    Bool_t hasJacbInt = opt_int;
    
    // Number of Res and Par
    Short_t numOfRes = (nseq_ + opt_int*nseg_*DIML);
    Short_t numOfPar = (DIMG  + opt_int*nseg_*DIML);

    Short_t parIDeta  =  4;
    Short_t parIDigb  = -1;
    Short_t parIDtsft = -1;
    if (opt_mu)   parIDigb  = 5;
    if (opt_tsft) parIDtsft = DIMG - 1;
    
    // Final State
    std::vector<PhySt> stts;

    // Reset TOF Time and Path
    Bool_t resetTOF = true;
    HitStTOF::SetOffsetTime(Numc::ZERO<>);
    HitStTOF::SetOffsetPath(Numc::ZERO<>);
    HitStTOF::SetTimeShiftCorr(opt_tsft);
    
    // TOF time shift
    Double_t TOFt_sft = (opt_tsft ? TOFt_sft_ : Numc::ZERO<>);
   
    Double_t chi_cx = 0;
    Double_t chi_cy = 0;
    Double_t chi_ib = 0;

    // Particle Status
    PhySt ppst(part_);
    
    // Matrix (Jb)
    PhyJb::SMtxDGG&& jbGG = SMtxId();
    std::vector<PhyJb::SMtxDGL> jbGL(nseg_);
    ceres::Matrix jb = ceres::Matrix::Zero(numOfRes, numOfPar);
    
    // Interaction Local Parameters
    if (opt_int) {
    for (Short_t is = 0; is < nseg_; ++is) { // Interaction
        SVecD<5> ichi, idiv;
        args_.at(is).cal_chi_and_div(ichi, idiv);
        
        chi_cx += ichi(0) * ichi(0); // tauu
        chi_cx += ichi(2) * ichi(2); // taul
        chi_cy += ichi(1) * ichi(1); // rhou
        chi_cy += ichi(3) * ichi(3); // rhol

        for (Short_t it = 0; it < DIML; ++it)
            jb(nseq_ + is*DIML+it, DIMG + is*DIML+it) += idiv(it);
    }} // Interaction
    
    Short_t cnt_nhit =  0;
    Short_t cnt_nseg = -1;
    PhySt nearPpst                       = ppst;
    PhyJb::SMtxDGG nearJbGG              = jbGG;
    std::vector<PhyJb::SMtxDGL> nearJbGL = jbGL;
    for (auto&& hit : hits_) {
        // Interaction Local Parameters
        Bool_t  hasLoc  = (cnt_nseg >= 0);
        Bool_t  isInner = (cnt_nseg >= 0 && cnt_nseg < nseg_);
        Short_t itnseg  = (cnt_nseg == nseg_) ? (nseg_ - Numc::ONE<Short_t>) : cnt_nseg;
        Bool_t  hasCxy  = (hit->scx() || hit->scy());
        
        if (isInner && !hasCxy) {
            ppst = nearPpst;
            jbGG = nearJbGG;
            jbGL = nearJbGL;
        }
        
        // Propagate
        PhyJb curjb;
        if (isInner) ppst.arg() = args_.at(cnt_nseg);
        else         ppst.arg().clear();
        if (!PropMgnt::PropToZ(hit->cz(), ppst, nullptr, &curjb)) break;
        ppst.symbk();
       
        // Hit Status: Setting TOF reference time and path
        if (resetTOF && Hit<HitStTOF>::IsSame(hit)) { // set reference
            const HitStTOF* firstHitTOF = Hit<HitStTOF>::Cast(hit);
            if (firstHitTOF->st()) {
                HitStTOF::SetOffsetPath(ppst.path());
                HitStTOF::SetOffsetTime((ppst.time() - firstHitTOF->orgt()) + TOFt_sft);
                resetTOF = false;
            }
        }
        hit->cal(ppst);
        
        // Update Jacb
        jbGG = curjb.gg() * jbGG;
        if (opt_int) {
            if (hasLoc && hasCxy) jbGL.at(cnt_nseg) = curjb.gl();
            for (Short_t is = 0; is < cnt_nseg; ++is)
                jbGL.at(is) = curjb.gg() * jbGL.at(is);
        }

        // State
        if (hasCxy) stts.push_back(ppst);

        // Coord
        if (hit->scx()) chi_cx += hit->chicx() * hit->chicx(); 
        if (hit->scy()) chi_cy += hit->chicy() * hit->chicy();
        if (hasCxy) {
            for (Short_t it = 0; it < PhyJb::DIMG; ++it) {
                if (hit->scx()) jb(hit->seqIDcx(), it) += hit->divcx() * jbGG(0, it);
                if (hit->scy()) jb(hit->seqIDcy(), it) += hit->divcy() * jbGG(1, it);
            }
            if (hasJacbInt) {
                for (Short_t is = 0; is <= itnseg; ++is) {
                    for (Short_t it = 0; it < PhyJb::DIML; ++it) {
                        if (hit->scx()) jb(hit->seqIDcx(), DIMG + is*DIML+it) += hit->divcx() * jbGL.at(is)(0, it);
                        if (hit->scy()) jb(hit->seqIDcy(), DIMG + is*DIML+it) += hit->divcy() * jbGL.at(is)(1, it);
                    }
                }
            } // Local
        } // hasCxy
         
        // TRK
        HitStTRK* hitTRK = Hit<HitStTRK>::Cast(hit);
        if (hitTRK != nullptr) {
            if (hitTRK->sqx()) chi_ib += hitTRK->chiqx() * hitTRK->chiqx();
            if (hitTRK->sqy()) chi_ib += hitTRK->chiqy() * hitTRK->chiqy();
            if (hitTRK->sqx()) {
                if (opt_mu) jb(hitTRK->seqIDqx(), parIDigb) += hitTRK->divqx_igb() * jbGG(4, 4);
                else        jb(hitTRK->seqIDqx(), parIDeta) += hitTRK->divqx_eta() * jbGG(4, 4);
            }
            if (hitTRK->sqy()) {
                if (opt_mu) jb(hitTRK->seqIDqy(), parIDigb) += hitTRK->divqy_igb() * jbGG(4, 4);
                else        jb(hitTRK->seqIDqy(), parIDeta) += hitTRK->divqy_eta() * jbGG(4, 4);
            }    
        }

        // TOF
        HitStTOF* hitTOF = Hit<HitStTOF>::Cast(hit);
        if (hitTOF != nullptr) {
            if (hitTOF->st()) chi_ib += hitTOF->chit() * hitTOF->chit();
            if (hitTOF->sq()) chi_ib += hitTOF->chiq() * hitTOF->chiq();
            if (hitTOF->st()) {
                if (opt_mu) jb(hitTOF->seqIDt(), parIDigb) += hitTOF->divt_igb() * jbGG(4, 4);
                else        jb(hitTOF->seqIDt(), parIDeta) += hitTOF->divt_eta() * jbGG(4, 4);
            }
            if (hitTOF->sq()) {
                if (opt_mu) jb(hitTOF->seqIDq(), parIDigb) += hitTOF->divq_igb() * jbGG(4, 4);
                else        jb(hitTOF->seqIDq(), parIDeta) += hitTOF->divq_eta() * jbGG(4, 4);
            }
            if (opt_tsft && hitTOF->st()) jb(hitTOF->seqIDt(), parIDtsft) += hitTOF->divt_sft(); // TOF time shift
        }
        
        // RICH
        HitStRICH* hitRICH = Hit<HitStRICH>::Cast(hit);
        if (hitRICH != nullptr) {
            if (hitRICH->sib()) chi_ib += hitRICH->chiib() * hitRICH->chiib();
            if (hitRICH->sib()) {
                if (opt_mu) jb(hitRICH->seqIDib(), parIDigb) += hitRICH->divib_igb() * jbGG(4, 4);
                else        jb(hitRICH->seqIDib(), parIDeta) += hitRICH->divib_eta() * jbGG(4, 4);
            }
        }
        
        // TRD
        HitStTRD* hitTRD = Hit<HitStTRD>::Cast(hit);
        if (hitTRD != nullptr) {
            if (hitTRD->sel()) chi_ib += hitTRD->nrmel() * hitTRD->nrmel();
            if (hitTRD->sel()) {
                if (opt_mu) jb(hitTRD->seqIDel(), parIDigb) += hitTRD->divel_igb() * jbGG(4, 4);
                else        jb(hitTRD->seqIDel(), parIDeta) += hitTRD->divel_eta() * jbGG(4, 4);
            }
        }
        
        if (hasCxy) {
            nearPpst = ppst;
            nearJbGG = jbGG;
            nearJbGL = jbGL;
            cnt_nseg++;
        }
        cnt_nhit++;
    }
    if (cnt_nhit != hits_.size()) return false;
    if (cnt_nseg != nseg_) return false;

    std::vector<Double_t> errs(numOfPar, Numc::ZERO<>);
    ceres::Matrix cov  = (jb.transpose() * jb);
    ceres::Vector diag = cov.inverse().diagonal();
    for (UInt_t it = 0; it < diag.rows(); ++it) {
        Double_t err = std::sqrt(diag(it));
        if (!Numc::Valid(diag(it)) || !Numc::Valid(err)) continue;
        errs.at(it) = err;
    }
    
    Double_t chi       = (chi_cx + chi_cy + chi_ib);
    Double_t nchi_tt   = (chi / static_cast<Double_t>(ndof_tt_));
    Double_t nchi_cyib = ((ndof_.at(1) > 0) ? ((chi_cy + chi_ib) / static_cast<Double_t>(ndof_.at(1))) : 0);

    for (auto&& stt : stts) stt.arg().clear();
    for (UInt_t it = 0; it < args_.size(); ++it) {
        PhyArg& arg = args_.at(it);
        stts.at(it).arg().set_mscat(arg.tauu(), arg.rhou(), arg.taul(), arg.rhol());
    }
    if (!hits_.at(0)->scx() && !hits_.at(0)->scy()) stts.insert(stts.begin(), part_);

    nchi_cx_ = ((ndof_cx_ > 0) ? (chi_cx / static_cast<Double_t>(ndof_cx_)) : 0);
    nchi_cy_ = ((ndof_cy_ > 0) ? (chi_cy / static_cast<Double_t>(ndof_cy_)) : 0);
    nchi_ib_ = ((ndof_ib_ > 0) ? (chi_ib / static_cast<Double_t>(ndof_ib_)) : 0);
    nchi_tt_ = nchi_tt;
    
    nchi_.at(0) = nchi_cx_;
    nchi_.at(1) = nchi_cyib;
    quality_.at(0) = NormQuality(nchi_.at(0), ndof_.at(0));
    quality_.at(1) = NormQuality(nchi_.at(1), ndof_.at(1));
    
    // Mu Error
    Double_t errMu = Numc::ZERO<>;
    if (opt_mu) {
        Double_t reEta = std::fabs(errs.at(4) / part_.eta());
        Double_t reIgb = (errs.at(parIDigb) / part_.igmbta());
        if (Numc::Compare(reIgb, Numc::HUNDRED<>) > 0) {
            errs.at(parIDigb) = Numc::HUNDRED<> * part_.igmbta();
            reIgb = Numc::HUNDRED<>;
        }
        Double_t reMu = std::hypot(reEta, reIgb);
        errMu = reMu * part_.mu();
        if (!Numc::Valid(errMu) || Numc::Compare(errMu) < 0) errMu = Numc::ZERO<>;
    }

    err_.fill(Numc::ZERO<>);
    err_.at(0) = errs.at(0);
    err_.at(1) = errs.at(1);
    err_.at(2) = errs.at(2);
    err_.at(3) = errs.at(3);
    err_.at(4) = errs.at(4);
    if (opt_mu) err_.at(5) = errs.at(parIDigb);
    if (opt_mu) err_.at(6) = errMu;

    stts_ = stts;
    info_ = part_.info();
 
    //CERR("MASS %14.8f RIG %14.8f QLT %14.8f %14.8f\n", part_.mass(), part_.rig(), quality_.at(0), quality_.at(1));
    return true;
}


bool VirtualPhyTrFit::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const {
    if (numOfRes_ <= 0 || numOfParGlb_ <= 0) return false;
    std::fill_n(residuals, numOfRes_, Numc::ZERO<>);
    Bool_t hasJacbGlb = (            jacobians != nullptr && jacobians[0] != nullptr);
    Bool_t hasJacbInt = (opt_int_ && jacobians != nullptr && jacobians[1] != nullptr);
    if (hasJacbGlb) std::fill_n(jacobians[0], numOfRes_ * numOfParGlb_, Numc::ZERO<>);
    if (hasJacbInt) std::fill_n(jacobians[1], numOfRes_ * numOfParInt_, Numc::ZERO<>);
    
    // Reset TOF Time and Path
    Bool_t resetTOF = true;
    HitStTOF::SetOffsetTime(Numc::ZERO<>);
    HitStTOF::SetOffsetPath(Numc::ZERO<>);
    HitStTOF::SetTimeShiftCorr(opt_tsft_);

    // TOF time shift
    Double_t TOFt_sft = (opt_tsft_ ? parameters[0][parIDtsft_] : Numc::ZERO<>);

    // Particle Status
    PhySt ppst(part_);
    Double_t partmu = (opt_mu_ ? std::fabs(parameters[0][parIDigb_] / parameters[0][parIDeta_]) : part_.mu());
    if (opt_mu_) ppst.reset(partmu);
    else         ppst.arg().clear();
    ppst.set_state_with_uxy(parameters[0][0], parameters[0][1], part_.cz(), parameters[0][2], parameters[0][3], Numc::Compare(part_.uz()));
    ppst.set_eta(parameters[0][parIDeta_]);

    // Interaction
    std::vector<PhyArg> args(nseg_, PhyArg(sw_mscat_, sw_eloss_));
    if (opt_int_) {
    for (Short_t is = 0; is < nseg_; ++is) {
        args.at(is).set_mscat(
            parameters[1][is*DIML_+0], 
            parameters[1][is*DIML_+1], 
            parameters[1][is*DIML_+2], 
            parameters[1][is*DIML_+3]);
    }}

    // Matrix (Rs, Jb)
    PhyJb::SMtxDGG&& jbGG = SMtxId();
    std::vector<PhyJb::SMtxDGL> jbGL(nseg_);
    ceres::Vector rs = ceres::Vector::Zero(numOfRes_);
    ceres::Matrix jb = ceres::Matrix::Zero(numOfRes_, numOfParGlb_ + numOfParInt_);

    if (opt_int_) {
    for (Short_t is = 0; is < nseg_; ++is) { // Interaction
        SVecD<5> inrm, idiv;
        args.at(is).cal_nrm_and_div(inrm, idiv);

        for (Short_t it = 0; it < DIML_; ++it)
            rs(nseq_ + is*DIML_+it) += inrm(it);

        if (hasJacbInt) {
            for (Short_t it = 0; it < DIML_; ++it)
                jb(nseq_ + is*DIML_+it, DIMG_ + is*DIML_+it) += idiv(it);
        } // hasJacbInt
    }} // Interaction
    
    Short_t cnt_nhit =  0;
    Short_t cnt_nseg = -1;
    PhySt                       nearPpst = ppst;
    PhyJb::SMtxDGG              nearJbGG = jbGG;
    std::vector<PhyJb::SMtxDGL> nearJbGL = jbGL;
    for (auto&& hit : hits_) {
        // Interaction Local Parameters
        Bool_t  hasLoc  = (cnt_nseg >= 0);
        Bool_t  isInner = (cnt_nseg >= 0 && cnt_nseg < nseg_);
        Short_t itnseg  = (cnt_nseg == nseg_) ? (nseg_ - Numc::ONE<Short_t>) : cnt_nseg;
        Bool_t  hasCxy  = (hit->scx() || hit->scy());

        if (isInner && !hasCxy) {
            ppst = nearPpst;
            jbGG = nearJbGG;
            jbGL = nearJbGL;
        }
       
        // Propagate
        PhyJb curjb;
        if (isInner) { // Internal Region
            ppst.arg().set_mscat(
                args.at(cnt_nseg).tauu(),
                args.at(cnt_nseg).rhou(),
                args.at(cnt_nseg).taul(),
                args.at(cnt_nseg).rhol());
        }
        else ppst.arg().clear(); // External Region
        if (!PropMgnt::PropToZ(hit->cz(), ppst, nullptr, ((hasJacbGlb)?&curjb:nullptr))) break;
        ppst.symbk();

        // Hit Status: Setting TOF reference time and path
        if (resetTOF && Hit<HitStTOF>::IsSame(hit)) { // set reference
            const HitStTOF* firstHitTOF = Hit<HitStTOF>::Cast(hit);
            if (firstHitTOF->st()) {
                HitStTOF::SetOffsetPath(ppst.path());
                HitStTOF::SetOffsetTime((ppst.time() - firstHitTOF->orgt()) + TOFt_sft);
                resetTOF = false;
            }
        }
        hit->cal(ppst);

        // Update Jacb
        if (hasJacbGlb) jbGG = curjb.gg() * jbGG;
        if (hasJacbInt) {
            if (hasLoc && hasCxy) jbGL.at(cnt_nseg) = curjb.gl();
            for (Short_t is = 0; is < cnt_nseg; ++is)
                jbGL.at(is) = curjb.gg() * jbGL.at(is);
        }
      
        // Coord
        if (hit->scx()) rs(hit->seqIDcx()) += hit->nrmcx();
        if (hit->scy()) rs(hit->seqIDcy()) += hit->nrmcy();
        if (hasCxy) {
            if (hasJacbGlb) {
                for (Short_t it = 0; it < PhyJb::DIMG; ++it) {
                    if (hit->scx()) jb(hit->seqIDcx(), it) += hit->divcx() * jbGG(0, it);
                    if (hit->scy()) jb(hit->seqIDcy(), it) += hit->divcy() * jbGG(1, it);
                }
            } // hasJacbGlb
            if (hasJacbInt && hasLoc) {
                for (Short_t is = 0; is <= itnseg; ++is) {
                    for (Short_t it = 0; it < PhyJb::DIML; ++it) {
                        if (hit->scx()) jb(hit->seqIDcx(), DIMG_ + is*DIML_+it) += hit->divcx() * jbGL.at(is)(0, it);
                        if (hit->scy()) jb(hit->seqIDcy(), DIMG_ + is*DIML_+it) += hit->divcy() * jbGL.at(is)(1, it);
                    }
                }
            } // hasJacbInt
        } // hasCxy

        // TRK
        HitStTRK* hitTRK = Hit<HitStTRK>::Cast(hit);
        if (hitTRK != nullptr) {
            if (hitTRK->sqx()) rs(hitTRK->seqIDqx()) += hitTRK->nrmqx();
            if (hitTRK->sqy()) rs(hitTRK->seqIDqy()) += hitTRK->nrmqy();
            if (hasJacbGlb && hitTRK->sqx()) {
                if (opt_mu_) jb(hitTRK->seqIDqx(), parIDigb_) += hitTRK->divqx_igb() * jbGG(4, 4);
                else         jb(hitTRK->seqIDqx(), parIDeta_) += hitTRK->divqx_eta() * jbGG(4, 4);
            } // hasJacbGlb
            if (hasJacbGlb && hitTRK->sqy()) {
                if (opt_mu_) jb(hitTRK->seqIDqy(), parIDigb_) += hitTRK->divqy_igb() * jbGG(4, 4);
                else         jb(hitTRK->seqIDqy(), parIDeta_) += hitTRK->divqy_eta() * jbGG(4, 4);
            } // hasJacbGlb
        }
        
        // TOF
        HitStTOF* hitTOF = Hit<HitStTOF>::Cast(hit);
        if (hitTOF != nullptr) {
            if (hitTOF->sq()) rs(hitTOF->seqIDq()) += hitTOF->nrmq();
            if (hitTOF->st()) rs(hitTOF->seqIDt()) += hitTOF->nrmt();
            if (hasJacbGlb && hitTOF->st()) {
                if (opt_mu_) jb(hitTOF->seqIDt(), parIDigb_) += hitTOF->divt_igb() * jbGG(4, 4);
                else         jb(hitTOF->seqIDt(), parIDeta_) += hitTOF->divt_eta() * jbGG(4, 4);
            } // hasJacbGlb
            if (hasJacbGlb && hitTOF->sq()) {
                if (opt_mu_) jb(hitTOF->seqIDq(), parIDigb_) += hitTOF->divq_igb() * jbGG(4, 4);
                else         jb(hitTOF->seqIDq(), parIDeta_) += hitTOF->divq_eta() * jbGG(4, 4);
            } // hasJacbGlb
            if (hasJacbGlb && opt_tsft_ && hitTOF->st()) jb(hitTOF->seqIDt(), parIDtsft_) += hitTOF->divt_sft(); // TOF time shift
        }
        
        // RICH
        HitStRICH* hitRICH = Hit<HitStRICH>::Cast(hit);
        if (hitRICH != nullptr) {
            if (hitRICH->sib()) rs(hitRICH->seqIDib()) += hitRICH->nrmib();
            if (hasJacbGlb && hitRICH->sib()) {
                if (opt_mu_) jb(hitRICH->seqIDib(), parIDigb_) += hitRICH->divib_igb() * jbGG(4, 4);
                else         jb(hitRICH->seqIDib(), parIDeta_) += hitRICH->divib_eta() * jbGG(4, 4);
            } // hasJacbGlb
        }
        
        // TRD
        HitStTRD* hitTRD = Hit<HitStTRD>::Cast(hit);
        if (hitTRD != nullptr) {
            if (hitTRD->sel()) rs(hitTRD->seqIDel()) += hitTRD->nrmel();
            if (hasJacbGlb && hitTRD->sel()) {
                if (opt_mu_) jb(hitTRD->seqIDel(), parIDigb_) += hitTRD->divel_igb() * jbGG(4, 4);
                else         jb(hitTRD->seqIDel(), parIDeta_) += hitTRD->divel_eta() * jbGG(4, 4);
            } // hasJacbGlb
        }

        if (hasCxy) {
            nearPpst = ppst;
            nearJbGG = jbGG;
            nearJbGL = jbGL;
            cnt_nseg++;
        }
        cnt_nhit++;
    }
    if (cnt_nhit != hits_.size()) return false;
    if (cnt_nseg != nseg_) return false;
    
    for (Short_t it = 0; it < numOfRes_; ++it) {
        if (!Numc::Valid(rs(it))) rs(it) = Numc::ZERO<>;
        residuals[it] = rs(it);
    }
    
    if (hasJacbGlb) {
        for (Short_t it = 0; it < numOfRes_; ++it) {
        for (Short_t jt = 0; jt < numOfParGlb_; ++jt) {
            if (!Numc::Valid(jb(it, jt))) jb(it, jt) = Numc::ZERO<>;
            jacobians[0][it * numOfParGlb_ + jt] = jb(it, jt); 
        }}
    }
    if (hasJacbInt) {
        for (Short_t it = 0; it < numOfRes_; ++it) {
        for (Short_t jt = 0; jt < numOfParInt_; ++jt) {
            if (!Numc::Valid(jb(it, DIMG_ + jt))) jb(it, DIMG_ + jt) = Numc::ZERO<>;
            jacobians[1][it * numOfParInt_ + jt] = jb(it, DIMG_ + jt); 
        }}
    }

    return true;
}
    

PhySt PhyTrFit::interpolate_to_z(Double_t zcoo) const {
    PhySt nullst = part_; nullst.reset(part_.info());
    if (!succ_ || stts_.size() == 0) return nullst;

    // Find Index
    Int_t idx = -1;
    if (ortt_ == Orientation::kDownward)
        for (UInt_t it = 0; it < stts_.size(); ++it)
            if (idx == -1 && Numc::Compare(zcoo, stts_.at(it).cz()) <= 0) idx = it;
    else // kUpward
        for (UInt_t it = 0; it < stts_.size(); ++it)
            if (idx == -1 && Numc::Compare(zcoo, stts_.at(it).cz()) >= 0) idx = it;

    // Set Particle
    PhySt ppst = ((idx == -1) ? stts_.at(0) : stts_.at(idx));
    if (idx == -1) ppst.arg().clear();
    PhyArg arg = ppst.arg();
        
    if (!PropMgnt::PropToZ(zcoo, ppst)) return nullst;
    ppst.symbk();

    ppst.arg().clear();
    ppst.arg().set_mscat(arg.tauu(), arg.rhou(), arg.taul(), arg.rhol());

    return ppst;
}


MatFld PhyTrFit::get_mat(Double_t zbd1, Double_t zbd2) const {
    if (!succ_ || stts_.size() == 0) return MatFld();
    
    // Set Boundary
    Short_t zbd_sign = Numc::Compare(zbd1, zbd2);
    if (zbd_sign == 0) {
        Double_t zcoo = zbd1;
        PhySt&& st = interpolate_to_z(zcoo);
        if (Numc::EqualToZero(st.mom())) return MatFld();
        else return MatMgnt::Get(Numc::ZERO<>, st);
    }
    Double_t zsat = ((ortt_ == Orientation::kDownward && zbd_sign > 0) || (ortt_ == Orientation::kUpward && zbd_sign < 0)) ? zbd1 : zbd2;
    Double_t zend = ((ortt_ == Orientation::kDownward && zbd_sign > 0) || (ortt_ == Orientation::kUpward && zbd_sign < 0)) ? zbd2 : zbd1;

    // Find Index
    Int_t isat = -1, iend = -1;
    if (ortt_ == Orientation::kDownward) {
        for (UInt_t it = 0; it < stts_.size(); ++it) {
            if (isat == -1 && Numc::Compare(zsat, stts_.at(it).cz()) <= 0) isat = it;
            if (iend == -1 && Numc::Compare(zend, stts_.at(it).cz()) <= 0) iend = it;
        }
    }
    else { // kUpward
        for (UInt_t it = 0; it < stts_.size(); ++it) {
            if (isat == -1 && Numc::Compare(zsat, stts_.at(it).cz()) >= 0) isat = it;
            if (iend == -1 && Numc::Compare(zend, stts_.at(it).cz()) >= 0) iend = it;
        }
    }

    // Material
    std::list<MatFld> mflds;
    for (Int_t it = isat; it <= iend; ++it) {
        PhySt ppst = (it == isat) ? interpolate_to_z(zsat) : stts_.at(it);
        if (Numc::EqualToZero(ppst.mom())) return MatFld();

        MatFld fld;
        Double_t zcoo = (it == iend) ? zend : stts_.at(it+1).cz();
        if (!PropMgnt::PropToZ(zcoo, ppst, &fld)) return MatFld();
        mflds.push_back(fld);
    }
    
    MatFld&& mfld = MatFld::Merge(mflds);
    return mfld;
}


Double_t PhyTrFit::NormQuality(Double_t nchi, Short_t ndof) {
    if (Numc::Compare(nchi) < 0 || ndof <= Numc::ZERO<Short_t>) return Numc::ZERO<>;
    Double_t chi = nchi * static_cast<Double_t>(ndof);
    if (Numc::EqualToZero(chi)) return Numc::ZERO<>;
    if (ndof <= Numc::TWO<Short_t>) return std::sqrt(nchi);
    Double_t qmin  = static_cast<Double_t>(ndof - Numc::TWO<Short_t>);
    Double_t sign  = static_cast<Double_t>(Numc::Compare(chi - qmin));
    Double_t qfunc = (chi - qmin) - qmin * std::log(chi / qmin);
    if (!Numc::Valid(qfunc)) return Numc::ZERO<>;
    Double_t xfunc = sign * std::sqrt(qfunc / static_cast<Double_t>(ndof));
    if (Numc::Valid(xfunc)) return xfunc;
    return Numc::ZERO<>;
}


} // namespace TrackSys


#endif // __TRACKLibs_PhyFit_C__
