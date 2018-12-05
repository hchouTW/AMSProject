#ifndef __TRACKLibs_PhyBtaFit_C__
#define __TRACKLibs_PhyBtaFit_C__


#include "Sys.h"
#include "Math.h"
#include "TmeMeas.h"
#include "IonEloss.h"
#include "GmIonEloss.h"
#include "PartInfo.h"
#include "PhySt.h"
#include "MagEnv.h"
#include "MatEnv.h"
#include "Prop.h"
#include "HitSt.h"
#include "TrFitPar.h"
#include "SimpleBtaFit.h"
#include "PhyBtaFit.h"


namespace TrackSys {
        
        
TrFitPar PhyBtaFit::bulidFitPar(const TrFitPar& fitPar) {
    TrFitPar btapar(fitPar.info(), fitPar.ortt(), fitPar.sw_mscat(), fitPar.sw_eloss());
    if (fitPar.hitsTOF().size()  != 0) btapar.add_hit( fitPar.hitsTOF()  );
    if (fitPar.hitsRICH().size() != 0) btapar.add_hit( fitPar.hitsRICH() );
    if (fitPar.hitsTRD().size()  != 0) btapar.add_hit( fitPar.hitsTRD()  );
    for (auto&& hitTRK : fitPar.hitsTRK()) {
        if (!hitTRK.sq()) continue;
        HitStTRK hit(false, false, hitTRK.lay(), hitTRK.isInnTr());
        hit.set_coo(hitTRK.cx(), hitTRK.cy(), hitTRK.cz());
        hit.set_q(hitTRK.q(), hitTRK.qx(), hitTRK.qy(), fitPar.info().chrg());
        btapar.add_hit(hit);
    }
    btapar.check();
    return btapar;
}

    
PhyBtaFit& PhyBtaFit::operator=(const PhyBtaFit& rhs) {
    if (this != &rhs) {
        dynamic_cast<TrFitPar&>(*this) = dynamic_cast<const TrFitPar&>(rhs);
        succ_ = rhs.succ_;
        part_ = rhs.part_;
        tsft_ = rhs.tsft_;
        rerr_ = rhs.rerr_;
      
        ndof_    = rhs.ndof_;
        nchi_    = rhs.nchi_;
        quality_ = rhs.quality_;
    }
    return *this;
}


void PhyBtaFit::clear() {
    succ_ = false;
    part_.reset(info_);
    part_.arg().reset(sw_mscat_, sw_eloss_);
    tsft_ = 0;
    rerr_ = 0;

    ndof_    = 0;
    nchi_    = 0;
    quality_ = 0;
}


PhyBtaFit::PhyBtaFit(const TrFitPar& fitPar, const PhySt& refSt) : TrFitPar(bulidFitPar(fitPar)) {
    PhyBtaFit::clear();
    if (!check_hits()) { PhyBtaFit::clear(); TrFitPar::clear(); return; }
    if (nmes_TOFt_ <= LMTN_TOF_T) { PhyBtaFit::clear(); TrFitPar::clear(); return; }
    
    ndof_ = nmes_ib_ - (Numc::ONE<Short_t> + (nmes_TOFt_ >= LMTN_TOF_T));
    if (ndof_ <= Numc::ONE<Short_t>) { PhyBtaFit::clear(); TrFitPar::clear(); return; }
    
    // check particle type
    if (refSt.info().type() != PartType::Fixed && refSt.info().type() != info_.type()) { PhyBtaFit::clear(); TrFitPar::clear(); return; }
    else if (!Numc::EqualToZero(refSt.mass()-info_.mass()) || !Numc::EqualToZero(refSt.chrg()-info_.chrg())) { PhyBtaFit::clear(); TrFitPar::clear(); return; }
    
    succ_ = (simpleFit(refSt) ? physicalFit() : false);
    if (!succ_) { PhyBtaFit::clear(); TrFitPar::clear(); }
   
    //if (!succ_) CERR("FAILURE === PhyBtaFit\n");
}


Bool_t PhyBtaFit::simpleFit(const PhySt& refSt) {
    SimpleBtaFit btafit(dynamic_cast<TrFitPar&>(*this), refSt);
    if (!btafit.status()) return false;

    // init state
    part_ = btafit.interpolate_to_z(hits_.at(0)->cz());
    if (Numc::EqualToZero(part_.mom())) return false;
    if (Numc::Compare(part_.uz() * refSt.uz()) < 0) return false;
    if (!survival_test_and_modify(part_, sw_eloss_)) return false;
    tsft_ = btafit.tsft();

    return true;
}


Bool_t PhyBtaFit::physicalFit() {
    if (Numc::EqualToZero(part_.mom())) return false;
    Bool_t opt_tsft = (nmes_TOFt_ > LMTN_TOF_T);
 
    // Gobal Parameters
    std::vector<double> params_glb({ part_.ibta() });
    if (opt_tsft) { params_glb.push_back(tsft_); } // time shift

    // CeresSolver: Cost Function
    ceres::CostFunction* cost_function = new VirtualPhyBtaFit(dynamic_cast<TrFitPar&>(*this), part_);

    // CeresSolver: Problem
    ceres::Problem problem;
    problem.AddResidualBlock(cost_function, nullptr, params_glb.data());
    problem.SetParameterLowerBound(params_glb.data(), parIDibta, LMTL_IBTA);
    problem.SetParameterUpperBound(params_glb.data(), parIDibta, LMTU_IBTA);

    // CeresSolver: Options
    ceres::Solver::Options options;
    options.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT;
    //options.trust_region_strategy_type = ceres::DOGLEG;
    options.max_num_iterations = 50;
    //options.max_solver_time_in_seconds = 5.0;

    // CeresSolver: Summary
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    if (!summary.IsSolutionUsable()) return false;
    //if (ceres::NO_CONVERGENCE == summary.termination_type) return false;
 
    // Result (Global)
    Double_t igmbta  = std::sqrt((params_glb.at(parIDibta) + Numc::ONE<>) * (params_glb.at(parIDibta) - Numc::ONE<>));
    Double_t parteta = part_.eta_sign() * (igmbta / part_.mu());
    part_.set_eta(parteta);
    
    if (opt_tsft) { tsft_ = params_glb.at(parIDtsft); } // TOF Shift Time

    Bool_t succ = evolve();
    return succ;
}


Bool_t PhyBtaFit::evolve() {
    Bool_t opt_tsft = (nmes_TOFt_ > LMTN_TOF_T);
   
    // Number of Res and Par
    Short_t numOfRes = (nseq_);
    Short_t numOfPar = (Numc::ONE<Short_t> + opt_tsft);

    // Final State
    std::vector<PhySt> stts;

    // Reset TOF Time and Path
    Bool_t resetTOF = true;
    HitStTOF::SetOffsetPathTime();
    
    // time shift
    Double_t tsft = (opt_tsft ? tsft_ : Numc::ZERO<>);
   
    Double_t chi = 0;

    // Particle Status
    PhySt ppst(part_);
    
    // Matrix (Jb)
    Double_t jbBB = Numc::ONE<>;
    ceres::Vector rs = ceres::Vector::Zero(numOfRes);
    ceres::Matrix jb = ceres::Matrix::Zero(numOfRes, numOfPar);
    
    Short_t cnt_nhit = 0;
    for (auto&& hit : hits_) {
        // Propagate
        Double_t preIbta = ppst.ibta();
        if (!PropMgnt::PropToZ(hit->cz(), ppst)) break;
        ppst.symbk();
       
        // Hit Status: Setting TOF reference time and path
        if (resetTOF && Hit<HitStTOF>::IsSame(hit)) { // set reference
            const HitStTOF* firstHitTOF = Hit<HitStTOF>::Cast(hit);
            if (firstHitTOF->st()) {
                HitStTOF::SetOffsetPathTime(
                    ppst.path(), 
                    (ppst.time() - firstHitTOF->orgt()) + tsft, 
                    true
                );
                resetTOF = false;
            }
        }
        hit->cal(ppst);
        
        // Update Jacb
        Double_t curJbBB = (ppst.ibta() / preIbta);
        jbBB = curJbBB * jbBB;
         
        // TRK
        HitStTRK* hitTRK = Hit<HitStTRK>::Cast(hit);
        if (hitTRK != nullptr) {
            if (hitTRK->sq()) chi += hitTRK->chiq() * hitTRK->chiq();
            if (hitTRK->sq()) rs(hitTRK->seqIDq()) += hitTRK->nrmq();
            if (hitTRK->sq()) jb(hitTRK->seqIDq(), parIDibta) += hitTRK->divq_ibta() * jbBB;
        }

        // TOF
        HitStTOF* hitTOF = Hit<HitStTOF>::Cast(hit);
        if (hitTOF != nullptr) {
            if (hitTOF->st()) chi += hitTOF->chit() * hitTOF->chit();
            if (hitTOF->sq()) chi += hitTOF->chiq() * hitTOF->chiq();
            if (hitTOF->st()) rs(hitTOF->seqIDt()) += hitTOF->nrmt();
            if (hitTOF->sq()) rs(hitTOF->seqIDq()) += hitTOF->nrmq();
            if (hitTOF->st()) jb(hitTOF->seqIDt(), parIDibta) += hitTOF->divt_ibta() * jbBB;
            if (hitTOF->sq()) jb(hitTOF->seqIDq(), parIDibta) += hitTOF->divq_ibta() * jbBB;
            if (hitTOF->st() && opt_tsft) jb(hitTOF->seqIDt(), parIDtsft) += hitTOF->divtsft(); // TOF time shift
        }
        
        // RICH
        HitStRICH* hitRICH = Hit<HitStRICH>::Cast(hit);
        if (hitRICH != nullptr) {
            if (hitRICH->sib()) chi += hitRICH->chiib() * hitRICH->chiib();
            if (hitRICH->sib()) rs(hitRICH->seqIDib()) += hitRICH->nrmib();
            if (hitRICH->sib()) jb(hitRICH->seqIDib(), parIDibta) += hitRICH->divib_ibta() * jbBB;
        }
        
        // TRD
        HitStTRD* hitTRD = Hit<HitStTRD>::Cast(hit);
        if (hitTRD != nullptr) {
            if (hitTRD->sel()) chi+= hitTRD->nrmel() * hitTRD->nrmel();
            if (hitTRD->sel()) rs(hitTRD->seqIDel()) += hitTRD->nrmel();
            if (hitTRD->sel()) jb(hitTRD->seqIDel(), parIDibta) += hitTRD->divel_igb() * jbBB;
        }
        
        cnt_nhit++;
    }
    if (cnt_nhit != hits_.size()) return false;

    Double_t ndof = static_cast<Double_t>(ndof_);
    nchi_    = (chi / ndof);
    quality_ = Numc::NormQuality(nchi_, ndof);

    ceres::Vector grd = (jb.transpose() * rs);
    ceres::Matrix hes = (jb.transpose() * jb);
    ceres::Matrix cov = hes.inverse();
    ceres::Vector diag = cov.diagonal();
    ceres::Vector step = cov * grd;
    Double_t errIbta = std::sqrt(diag(parIDibta));
    Double_t errTsft = (opt_tsft ? std::sqrt(diag(parIDtsft)) : Numc::ZERO<>);
    if (!Numc::Valid(errIbta)) errIbta = Numc::ZERO<>;
    if (!Numc::Valid(errTsft)) errTsft = Numc::ZERO<>;

    rerr_ = (errIbta / part_.ibta());

    return true;
}


bool VirtualPhyBtaFit::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const {
    if (numOfRes_ <= 0 || numOfPar_ <= 0) return false;
    std::fill_n(residuals, numOfRes_, Numc::ZERO<>);
    Bool_t hasJacb = (jacobians != nullptr && jacobians[0] != nullptr);
    if (hasJacb) std::fill_n(jacobians[0], numOfRes_ * numOfPar_, Numc::ZERO<>);
    
    // Reset TOF Time and Path
    Bool_t resetTOF = true;
    HitStTOF::SetOffsetPathTime();

    // time shift
    Double_t tsft = (opt_tsft_ ? parameters[0][parIDtsft] : Numc::ZERO<>);

    // Particle Status
    Double_t igmbta  = std::sqrt((parameters[0][parIDibta] + Numc::ONE<>) * (parameters[0][parIDibta] - Numc::ONE<>));
    Double_t parteta = part_.eta_sign() * (igmbta / part_.mu());
    PhySt ppst(part_);
    ppst.arg().clear();
    ppst.set_eta(parteta);

    Double_t jbBB = Numc::ONE<>;
    ceres::Vector rs = ceres::Vector::Zero(numOfRes_);
    ceres::Matrix jb = ceres::Matrix::Zero(numOfRes_, numOfPar_);

    Short_t cnt_nhit =  0;
    for (auto&& hit : hits_) {
        // Propagate
        Double_t preIbta = ppst.ibta();
        if (!PropMgnt::PropToZ(hit->cz(), ppst)) break;
        ppst.symbk();

        // Hit Status: Setting TOF reference time and path
        if (resetTOF && Hit<HitStTOF>::IsSame(hit)) { // set reference
            const HitStTOF* firstHitTOF = Hit<HitStTOF>::Cast(hit);
            if (firstHitTOF->st()) {
                HitStTOF::SetOffsetPathTime(
                    ppst.path(),
                    (ppst.time() - firstHitTOF->orgt()) + tsft,
                    true
                );
                resetTOF = false;
            }
        }
        hit->cal(ppst);

        // Update Jacb
        Double_t curJbBB = (ppst.ibta() / preIbta);
        if (hasJacb) jbBB = curJbBB * jbBB;

        // TRK
        HitStTRK* hitTRK = Hit<HitStTRK>::Cast(hit);
        if (hitTRK != nullptr) {
            if (hitTRK->sq()) rs(hitTRK->seqIDq()) += hitTRK->nrmq();
            if (hasJacb && hitTRK->sq()) jb(hitTRK->seqIDq(), parIDibta) += hitTRK->divq_ibta() * jbBB;
        }
        
        // TOF
        HitStTOF* hitTOF = Hit<HitStTOF>::Cast(hit);
        if (hitTOF != nullptr) {
            if (hitTOF->sq()) rs(hitTOF->seqIDq()) += hitTOF->nrmq();
            if (hitTOF->st()) rs(hitTOF->seqIDt()) += hitTOF->nrmt();
            if (hasJacb && hitTOF->st()) jb(hitTOF->seqIDt(), parIDibta) += hitTOF->divt_ibta() * jbBB;
            if (hasJacb && hitTOF->sq()) jb(hitTOF->seqIDq(), parIDibta) += hitTOF->divq_ibta() * jbBB;
            if (hasJacb && hitTOF->st() && opt_tsft_) jb(hitTOF->seqIDt(), parIDtsft) += hitTOF->divtsft(); // TOF time shift
        }
        
        // RICH
        HitStRICH* hitRICH = Hit<HitStRICH>::Cast(hit);
        if (hitRICH != nullptr) {
            if (hitRICH->sib()) rs(hitRICH->seqIDib()) += hitRICH->nrmib();
            if (hasJacb && hitRICH->sib()) jb(hitRICH->seqIDib(), parIDibta) += hitRICH->divib_ibta() * jbBB;
        }
        
        // TRD
        HitStTRD* hitTRD = Hit<HitStTRD>::Cast(hit);
        if (hitTRD != nullptr) {
            if (hitTRD->sel()) rs(hitTRD->seqIDel()) += hitTRD->nrmel();
            if (hasJacb && hitTRD->sel()) jb(hitTRD->seqIDel(), parIDibta) += hitTRD->divel_igb() * jbBB;
        }

        cnt_nhit++;
    }
    if (cnt_nhit != hits_.size()) return false;

    for (Short_t it = 0; it < numOfRes_; ++it) {
        if (!Numc::Valid(rs(it))) rs(it) = Numc::ZERO<>;
        residuals[it] = rs(it);
    }
    
    if (hasJacb) {
        for (Short_t it = 0; it < numOfRes_; ++it) {
        for (Short_t jt = 0; jt < numOfPar_; ++jt) {
            if (!Numc::Valid(jb(it, jt))) jb(it, jt) = Numc::ZERO<>;
            jacobians[0][it * numOfPar_ + jt] = jb(it, jt); 
        }}
    }

    return true;
}


PhySt PhyBtaFit::interpolate_to_z(Double_t zcoo) const {
    PhySt nullst = part_; nullst.reset(part_.info());
    if (!succ_) return nullst;
   
    PhySt ppst(part_);
    if (!PropMgnt::PropToZ(zcoo, ppst)) return nullst;
    ppst.symbk();
    return ppst;
};


MatFld PhyBtaFit::get_mat(Double_t zbd1, Double_t zbd2) const {
    if (!succ_) return MatFld();
    
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
    
    MatFld mfld;
    PhySt ppst(part_);
    if (!PropMgnt::PropToZ(zsat, ppst)) return MatFld();
    if (!PropMgnt::PropToZ(zend, ppst, &mfld)) return MatFld();
    return mfld;
};


} // namespace TrackSys


#endif // __TRACKLibs_PhyBtaFit_C__
