#ifndef __TRACKLibs_PhyTrFit_C__
#define __TRACKLibs_PhyTrFit_C__


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
#include "TrFitPar.h"
#include "SimpleTrFit.h"
#include "PhyTrFit.h"


namespace TrackSys {
        

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
    ndof_.at(0) = ndof_cx_;
    ndof_.at(1) = ndof_cy_ + ndof_ib_;
    if (MuOpt::kFree == mu_opt_) { ndof_ib_ -= 1; ndof_tt_ -= 1; }
    if (ndof_.at(0) <= Numc::ONE<Short_t>) { PhyTrFit::clear(); return; }
    if (ndof_.at(1) <= Numc::ONE<Short_t>) { PhyTrFit::clear(); return; }
    if (MuOpt::kFree == mu_opt_) {
        if (ndof_cy_ <= Numc::ONE<Short_t>) { PhyTrFit::clear(); return; }
        if (ndof_ib_ <= Numc::ONE<Short_t>) { PhyTrFit::clear(); return; }
    }


    // Fitting
    if (MuOpt::kFixed == mu_opt_) succ_ = physicalTrFit();
    else                          succ_ = physicalMuFit();
    
    if (!succ_) { PhyTrFit::clear(); TrFitPar::clear(); }
    
    //if (!succ_) CERR("FAILURE === PhyTrFit\n"); // testcode
}


Bool_t PhyTrFit::simpleFit() {
    SimpleTrFit simple(dynamic_cast<TrFitPar&>(*this), true);
    if (simple.status()) {
        part_     = simple.part();
        TOFt_sft_ = simple.tsft();
        args_     = simple.args();
    }
    else {
        TOFt_sft_ = Numc::ZERO<>;
        args_     = std::move(std::vector<PhyArg>(nseg_, PhyArg(sw_mscat_, sw_eloss_)));
    }
    return simple.status();
}


Bool_t PhyTrFit::physicalFit(Bool_t fixedGlb, const MuOpt& mu_opt, Double_t fluc_eta, Double_t fluc_igb, Bool_t with_mu_est) {
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
    //if (!fixedGlb) TOFt_sft_ = Numc::ZERO<>;
    std::vector<double> params_glb({ part_.cx(), part_.cy(), part_.ux(), part_.uy(), eta });

    Short_t parIDeta  =  4;
    Short_t parIDigb  = -1;
    Short_t parIDtsft = -1;
    if (opt_mu)   { parIDigb  = 5;        params_glb.push_back(igb); }
    if (opt_tsft) { parIDtsft = DIMG - 1; params_glb.push_back(TOFt_sft_); } // TOF Shift Time  

    // Interaction Parameters
    //if (fixedGlb) args_ = std::move(std::vector<PhyArg>(nseg_, PhyArg(sw_mscat_, sw_eloss_))); // testocde
    std::vector<double> params_int(nseg_*DIML, Numc::ZERO<>);
    for (Short_t is = 0; is < nseg_; ++is) {
        params_int.at(is*PhyJb::DIML+0) = args_.at(is).tauu();
        params_int.at(is*PhyJb::DIML+1) = args_.at(is).rhou();
        params_int.at(is*PhyJb::DIML+2) = args_.at(is).taul();
        params_int.at(is*PhyJb::DIML+3) = args_.at(is).rhol();
    }

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

    // testcode
    //if (fixedGlb) problem.SetParameterBlockConstant(params_glb.data());

    // CeresSolver: Options
    ceres::Solver::Options options;
    //options.check_gradients = true;
    options.max_solver_time_in_seconds = 10.0;
    options.max_num_iterations = 500;
    options.max_num_line_search_direction_restarts = 50;
    options.max_num_line_search_step_size_iterations = 100;

    // CeresSolver: Summary
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    if (!summary.IsSolutionUsable()) return false;
    //std::cerr << summary.FullReport() << std::endl;
   
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
    Bool_t succ = physicalFit(false, MuOpt::kFixed);
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
        if (!(simpleFit() ? physicalFit(false, MuOpt::kFixed, Numc::ZERO<>, Numc::ZERO<>, false) : false)) continue;
        
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

        if (!physicalFit(false, MuOpt::kFree,
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
    if (!(hits_.at(0)->scx() || hits_.at(0)->scy())) stts.insert(stts.begin(), part_);

    nchi_cx_ = ((ndof_cx_ > 0) ? (chi_cx / static_cast<Double_t>(ndof_cx_)) : 0);
    nchi_cy_ = ((ndof_cy_ > 0) ? (chi_cy / static_cast<Double_t>(ndof_cy_)) : 0);
    nchi_ib_ = ((ndof_ib_ > 0) ? (chi_ib / static_cast<Double_t>(ndof_ib_)) : 0);
    nchi_tt_ = nchi_tt;
    
    nchi_.at(0) = nchi_cx_;
    nchi_.at(1) = nchi_cyib;
    quality_.at(0) = TrFitPar::NormQuality(nchi_.at(0), ndof_.at(0));
    quality_.at(1) = TrFitPar::NormQuality(nchi_.at(1), ndof_.at(1));
    
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


} // namespace TrackSys


#endif // __TRACKLibs_PhyTrFit_C__
