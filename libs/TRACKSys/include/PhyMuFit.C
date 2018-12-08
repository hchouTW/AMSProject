#ifndef __TRACKLibs_PhyMuFit_C__
#define __TRACKLibs_PhyMuFit_C__


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
#include "SimpleTrFit.h"
#include "PhyTrFit.h"
#include "SimpleMuScan.h"
#include "PhyMuFit.h"


namespace TrackSys {
        

PhyMuFit& PhyMuFit::operator=(const PhyMuFit& rhs) {
    if (this != &rhs) {
        dynamic_cast<TrFitPar&>(*this) = dynamic_cast<const TrFitPar&>(rhs);
        succ_ = rhs.succ_;
        part_ = rhs.part_;
        tsft_ = rhs.tsft_;
        
        args_ = rhs.args_;
        stts_ = rhs.stts_;
      
        ndof_    = rhs.ndof_;
        nchi_    = rhs.nchi_;
        quality_ = rhs.quality_;

        ndof_tt_ = rhs.ndof_tt_;
        nchi_tt_ = rhs.nchi_tt_;
    }
    return *this;
}


void PhyMuFit::clear() {
    succ_   = false;
    part_.reset(info_);
    part_.arg().reset(sw_mscat_, sw_eloss_);
    tsft_ = 0;
    
    args_.clear();
    stts_.clear();

    ndof_.fill(0);
    nchi_.fill(0);
    quality_.fill(0);

    ndof_tt_ = 0;
    nchi_tt_ = 0;
}


PhyMuFit::PhyMuFit(const TrFitPar& fitPar) : TrFitPar(fitPar) {
    PhyMuFit::clear();
    if (!check_hits()) { PhyMuFit::clear(); TrFitPar::clear(); return; }
    if (nmes_cx_ <= LMTN_CX) { PhyMuFit::clear(); TrFitPar::clear(); return; }
    if (nmes_cy_ <= LMTN_CY) { PhyMuFit::clear(); TrFitPar::clear(); return; }
    if (nmes_TOFt_ <= LMTN_TOF_T) { PhyMuFit::clear(); TrFitPar::clear(); return; }
    
    if (!(simpleScan(fitPar) && check_hits())) { PhyMuFit::clear(); TrFitPar::clear(); return; }
   
    Short_t ndof_cx = (nmes_cx_ > LMTN_CX) ? (nmes_cx_ - LMTN_CX) : 0;
    Short_t ndof_cy = (nmes_cy_ > LMTN_CY) ? (nmes_cy_ - LMTN_CY) : 0;
    Short_t ndof_ib = nmes_ib_ - (Numc::ONE<Short_t> + (nmes_TOFt_ >= LMTN_TOF_T));
    ndof_tt_ = (ndof_cx + ndof_cy + ndof_ib);
    
    ndof_.at(0) = ndof_cx;
    ndof_.at(1) = ndof_cy;
    ndof_.at(2) = ndof_ib;
    if (ndof_.at(0) <= Numc::ONE<Short_t>) { PhyMuFit::clear(); TrFitPar::clear(); return; }
    if (ndof_.at(1) <= Numc::ONE<Short_t>) { PhyMuFit::clear(); TrFitPar::clear(); return; }
    if (ndof_.at(2) <= Numc::ONE<Short_t>) { PhyMuFit::clear(); TrFitPar::clear(); return; }

    succ_ = physicalFit();
    if (!succ_) { PhyMuFit::clear(); TrFitPar::clear(); return; }
    
    // testcode
    //COUT("PHY MASS %14.8f RIG %14.8f QLT %14.8f %14.8f\n", part_.mass(), part_.rig(), quality_.at(0), quality_.at(1));
    
    //if (!succ_) CERR("FAILURE === PhyMuFit\n");
}


Bool_t PhyMuFit::simpleScan(const TrFitPar& fitPar) {
    SimpleMuScan muscan(fitPar);
    if (!muscan.status()) {
        args_ = std::move(std::vector<PhyArg>(nseg_, PhyArg(sw_mscat_, sw_eloss_)));
        return false;
    }
    else {
        info_ = muscan.part().info();
        part_ = muscan.part();
        tsft_ = muscan.tsft();
        args_ = muscan.args();
    }
    return true;
}


Bool_t PhyMuFit::physicalFit() {
    if (Numc::EqualToZero(part_.mom())) return false;
    Bool_t  opt_loc  = sw_mscat_;
    Bool_t  opt_tsft = (nmes_TOFt_ > LMTN_TOF_T);
    Short_t DIMG     = (PhyJb::DIMG+1) + opt_tsft;
    
    // Gobal Parameters
    std::vector<double> params_glb({ part_.cx(), part_.cy(), part_.ux(), part_.uy(), part_.eta(), part_.ibta() });
    if (opt_tsft) { params_glb.push_back(tsft_); } // time shift

    // Local Parameters
    std::vector<double> params_loc(nseg_*PhyJb::DIML, Numc::ZERO<>);
    if (opt_loc && args_.size() == nseg_) {
    for (Short_t is = 0; is < nseg_; ++is) {
        params_loc.at(is*PhyJb::DIML+0) = args_.at(is).tauu();
        params_loc.at(is*PhyJb::DIML+1) = args_.at(is).rhou();
        params_loc.at(is*PhyJb::DIML+2) = args_.at(is).taul();
        params_loc.at(is*PhyJb::DIML+3) = args_.at(is).rhol();
    }}
    
    // CeresSolver: Cost Function
    ceres::CostFunction* cost_function = new VirtualPhyMuFit(dynamic_cast<TrFitPar&>(*this), part_);

    // CeresSolver: Problem
    ceres::Problem problem;
    if (opt_loc) problem.AddResidualBlock(cost_function, nullptr, params_glb.data(), params_loc.data());
    else         problem.AddResidualBlock(cost_function, nullptr, params_glb.data());
    problem.SetParameterLowerBound(params_glb.data(), 2, -1.0);
    problem.SetParameterUpperBound(params_glb.data(), 2,  1.0);
    problem.SetParameterLowerBound(params_glb.data(), 3, -1.0);
    problem.SetParameterUpperBound(params_glb.data(), 3,  1.0);
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
    //std::cerr << summary.FullReport() << std::endl;
   
    // Result (Global)
    Double_t partcz = part_.cz();
    Short_t  signuz = Numc::Compare(part_.uz());
    Double_t igmbta = std::sqrt((params_glb.at(parIDibta) + Numc::ONE<>) * (params_glb.at(parIDibta) - Numc::ONE<>));
    Double_t partmu = std::fabs(igmbta / params_glb.at(parIDeta));
    if (Numc::Compare(partmu) <= 0) return false;

    info_.reset(partmu);
    part_.reset(partmu);
    part_.set_state_with_uxy(params_glb.at(0), params_glb.at(1), partcz, params_glb.at(2), params_glb.at(3), signuz);
    part_.set_eta(params_glb.at(parIDeta));
    
    if (opt_tsft) { tsft_ = params_glb.at(parIDtsft); } // TOF Shift Time

    // Result (Local)
    if (opt_loc) {
    for (Short_t is = 0; is < nseg_; ++is) {
        args_.at(is).set_mscat(
            params_loc.at(is*PhyJb::DIML+0), 
            params_loc.at(is*PhyJb::DIML+1), 
            params_loc.at(is*PhyJb::DIML+2), 
            params_loc.at(is*PhyJb::DIML+3));
    }}

    Bool_t succ = evolve();
    return succ;
}


Bool_t PhyMuFit::evolve() {
    Bool_t  opt_loc  = (sw_mscat_);
    Bool_t  opt_tsft = (nmes_TOFt_ > LMTN_TOF_T);
    Short_t DIMG     = (PhyJb::DIMG+1) + opt_tsft;
   
    // Number of Res and Par
    Short_t numOfRes = (nseq_ + opt_loc*nseg_*PhyJb::DIML);
    Short_t numOfPar = (DIMG  + opt_loc*nseg_*PhyJb::DIML);

    // Final State
    std::vector<PhySt> stts;

    // Reset TOF Time and Path
    Bool_t resetTOF = true;
    HitStTOF::SetOffsetPathTime();
    
    // time shift
    Double_t tsft = (opt_tsft ? tsft_ : Numc::ZERO<>);
   
    Double_t chi_cx = 0;
    Double_t chi_cy = 0;
    Double_t chi_ib = 0;

    // Particle Status
    PhySt ppst(part_);
    
    // Matrix (Jb)
    Double_t         jbBB = Numc::ONE<>;
    PhyJb::SMtxDGG&& jbGG = SMtxId();
    std::vector<PhyJb::SMtxDGL> jbGL(nseg_);
    std::vector<PhyJb::SMtxDGE> jbGE(nseg_);
    ceres::Matrix jb = ceres::Matrix::Zero(numOfRes, numOfPar);
    
    // Interaction Local Parameters
    if (opt_loc) {
    for (Short_t is = 0; is < nseg_; ++is) { // Interaction
        SVecD<PhyJb::DIMI> ichi, idiv;
        args_.at(is).cal_chi_and_div(ichi, idiv);
        
        chi_cx += ichi(0) * ichi(0); // tauu
        chi_cx += ichi(2) * ichi(2); // taul
        chi_cy += ichi(1) * ichi(1); // rhou
        chi_cy += ichi(3) * ichi(3); // rhol

        for (Short_t it = 0; it < PhyJb::DIML; ++it)
            jb(nseq_ + is*PhyJb::DIML+it, DIMG + is*PhyJb::DIML+it) += idiv(it);
    }} // Interaction
    
    Short_t cnt_nhit =  0;
    Short_t cnt_nseg = -1;
    PhySt nearPpst                       = ppst;
    Double_t       nearJbBB              = jbBB;
    PhyJb::SMtxDGG nearJbGG              = jbGG;
    std::vector<PhyJb::SMtxDGL> nearJbGL = jbGL;
    std::vector<PhyJb::SMtxDGE> nearJbGE = jbGE;
    for (auto&& hit : hits_) {
        // Interaction Local Parameters
        Bool_t  hasLoc  = (cnt_nseg >= 0);
        Bool_t  isInner = (cnt_nseg >= 0 && cnt_nseg < nseg_);
        Short_t itnseg  = (cnt_nseg == nseg_) ? (nseg_ - Numc::ONE<Short_t>) : cnt_nseg;
        Bool_t  hasCxy  = (hit->scx() || hit->scy());
        
        if (isInner && !hasCxy) {
            ppst = nearPpst;
            jbBB = nearJbBB;
            jbGG = nearJbGG;
            jbGL = nearJbGL;
            jbGE = nearJbGE;
        }
        
        // Propagate
        PhyJb curjb;
        if (isInner) ppst.arg() = args_.at(cnt_nseg);
        else         ppst.arg().clear();
        Double_t preIbta = ppst.ibta();
        if (!PropMgnt::PropToZ(hit->cz(), ppst, nullptr, &curjb)) break;
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
        jbGG = curjb.gg() * jbGG;
        if (opt_loc) {
            if (hasLoc && hasCxy) jbGL.at(cnt_nseg) = curjb.gl();
            for (Short_t is = 0; is < cnt_nseg; ++is)
                jbGL.at(is) = curjb.gg() * jbGL.at(is);
            
            if (hasLoc && hasCxy) jbGE.at(cnt_nseg) = curjb.ge();
            for (Short_t is = 0; is < cnt_nseg; ++is)
                jbGE.at(is) = curjb.gg() * jbGE.at(is);
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
            if (opt_loc) {
                for (Short_t is = 0; is <= itnseg; ++is) {
                    for (Short_t it = 0; it < PhyJb::DIML; ++it) {
                        if (hit->scx()) jb(hit->seqIDcx(), DIMG + is*PhyJb::DIML+it) += hit->divcx() * jbGL.at(is)(0, it);
                        if (hit->scy()) jb(hit->seqIDcy(), DIMG + is*PhyJb::DIML+it) += hit->divcy() * jbGL.at(is)(1, it);
                    }
                }
            } // Local
        } // hasCxy
         
        // TRK
        HitStTRK* hitTRK = Hit<HitStTRK>::Cast(hit);
        if (hitTRK != nullptr) {
            if (hitTRK->sq()) chi_ib += hitTRK->chiq() * hitTRK->chiq();
            if (hitTRK->sq()) jb(hitTRK->seqIDq(), parIDibta) += hitTRK->divq_ibta() * jbBB;
        }

        // TOF
        HitStTOF* hitTOF = Hit<HitStTOF>::Cast(hit);
        if (hitTOF != nullptr) {
            if (hitTOF->st()) chi_ib += hitTOF->chit() * hitTOF->chit();
            if (hitTOF->sq()) chi_ib += hitTOF->chiq() * hitTOF->chiq();
            if (hitTOF->st()) jb(hitTOF->seqIDt(), parIDibta) += hitTOF->divt_ibta() * jbBB;
            if (hitTOF->sq()) jb(hitTOF->seqIDq(), parIDibta) += hitTOF->divq_ibta() * jbBB;
            if (hitTOF->st() && opt_tsft) jb(hitTOF->seqIDt(), parIDtsft) += hitTOF->divtsft(); // TOF time shift
        }
        
        // RICH
        HitStRICH* hitRICH = Hit<HitStRICH>::Cast(hit);
        if (hitRICH != nullptr) {
            if (hitRICH->sib()) chi_ib += hitRICH->chiib() * hitRICH->chiib();
            if (hitRICH->sib()) jb(hitRICH->seqIDib(), parIDibta) += hitRICH->divib_ibta() * jbBB;
        }
        
        // TRD
        HitStTRD* hitTRD = Hit<HitStTRD>::Cast(hit);
        if (hitTRD != nullptr) {
            if (hitTRD->sel()) chi_ib += hitTRD->nrmel() * hitTRD->nrmel();
            if (hitTRD->sel()) jb(hitTRD->seqIDel(), parIDeta) += hitTRD->divel_eta() * jbGG(4, 4);
        }
        
        if (hasCxy) {
            nearPpst = ppst;
            nearJbBB = jbBB;
            nearJbGG = jbGG;
            nearJbGL = jbGL;
            nearJbGE = jbGE;
            cnt_nseg++;
        }
        cnt_nhit++;
    }
    if (cnt_nhit != hits_.size()) return false;
    if (cnt_nseg != nseg_) return false;

    Double_t chi     = (chi_cx + chi_cy + chi_ib);
    Double_t ndof_tt = static_cast<Double_t>(ndof_tt_);
    Double_t nchi_tt = (chi / ndof_tt);

    nchi_.at(0) = ((ndof_.at(0) > 0) ? (chi_cx / static_cast<Double_t>(ndof_.at(0))) : 0);
    nchi_.at(1) = ((ndof_.at(1) > 0) ? (chi_cy / static_cast<Double_t>(ndof_.at(1))) : 0);
    nchi_.at(2) = ((ndof_.at(2) > 0) ? (chi_ib / static_cast<Double_t>(ndof_.at(2))) : 0);
    nchi_tt_    = nchi_tt;
    
    quality_.at(0) = Numc::NormQuality(nchi_.at(0), ndof_.at(0));
    quality_.at(1) = Numc::NormQuality(nchi_.at(1), ndof_.at(1));
    quality_.at(2) = Numc::NormQuality(nchi_.at(2), ndof_.at(2));
   
    // Local States
    for (auto&& stt : stts) stt.arg().clear();
    for (UInt_t it = 0; it < args_.size(); ++it) {
        PhyArg& arg = args_.at(it);
        stts.at(it).arg().set_mscat(arg.tauu(), arg.rhou(), arg.taul(), arg.rhol());
        stts.at(it).arg().set_eloss(arg.elion(), arg.elbrm());
    }
    if (!(hits_.at(0)->scx() || hits_.at(0)->scy())) stts.insert(stts.begin(), part_);
    stts_ = stts;
    
    return true;
}


bool VirtualPhyMuFit::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const {
    if (numOfRes_ <= 0 || numOfParGlb_ <= 0) return false;
    std::fill_n(residuals, numOfRes_, Numc::ZERO<>);
    Bool_t hasJacbGlb = (            jacobians != nullptr && jacobians[0] != nullptr);
    Bool_t hasJacbLoc = (opt_loc_ && jacobians != nullptr && jacobians[1] != nullptr);
    if (hasJacbGlb) std::fill_n(jacobians[0], numOfRes_ * numOfParGlb_, Numc::ZERO<>);
    if (hasJacbLoc) std::fill_n(jacobians[1], numOfRes_ * numOfParLoc_, Numc::ZERO<>);
    
    // Reset TOF Time and Path
    Bool_t resetTOF = true;
    HitStTOF::SetOffsetPathTime();

    // time shift
    Double_t tsft = (opt_tsft_ ? parameters[0][parIDtsft] : Numc::ZERO<>);

    // Particle Status
    Double_t partcz = part_.cz();
    Short_t  signuz = Numc::Compare(part_.uz());
    Double_t igmbta = std::sqrt((parameters[0][parIDibta] + Numc::ONE<>) * (parameters[0][parIDibta] - Numc::ONE<>));
    Double_t partmu = std::fabs(igmbta / parameters[0][parIDeta]);
    if (Numc::Compare(partmu) <= 0) return false;
    
    PhySt ppst(part_);
    ppst.arg().clear();
    ppst.reset(partmu);
    ppst.set_state_with_uxy(parameters[0][0], parameters[0][1], partcz, parameters[0][2], parameters[0][3], signuz);
    ppst.set_eta(parameters[0][parIDeta]);

    // Interaction Local Parameters
    std::vector<PhyArg> args(nseg_, PhyArg(sw_mscat_, sw_eloss_));
    if (opt_loc_) {
    for (Short_t is = 0; is < nseg_; ++is) {
        args.at(is).set_mscat(
            parameters[1][is*PhyJb::DIML+0], 
            parameters[1][is*PhyJb::DIML+1], 
            parameters[1][is*PhyJb::DIML+2], 
            parameters[1][is*PhyJb::DIML+3]);
    }}

    // Matrix (Rs, Jb)
    Double_t         jbBB = Numc::ONE<>;
    PhyJb::SMtxDGG&& jbGG = SMtxId();
    std::vector<PhyJb::SMtxDGL> jbGL(nseg_);
    std::vector<PhyJb::SMtxDGE> jbGE(nseg_);
    ceres::Vector rs = ceres::Vector::Zero(numOfRes_);
    ceres::Matrix jb = ceres::Matrix::Zero(numOfRes_, numOfParGlb_ + numOfParLoc_);

    // Interaction Local Parameters
    if (opt_loc_) {
    for (Short_t is = 0; is < nseg_; ++is) {
        SVecD<PhyJb::DIMI> inrm, idiv;
        args.at(is).cal_nrm_and_div(inrm, idiv);

        for (Short_t it = 0; it < PhyJb::DIML; ++it)
            rs(nseq_ + is*PhyJb::DIML+it) += inrm(it);

        if (hasJacbLoc) {
        for (Short_t it = 0; it < PhyJb::DIML; ++it) {
            jb(nseq_ + is*PhyJb::DIML+it, DIMG_ + is*PhyJb::DIML+it) += idiv(it);
        }}
    }} // Interaction
    
    Short_t cnt_nhit =  0;
    Short_t cnt_nseg = -1;
    PhySt                       nearPpst = ppst;
    Double_t                    nearJbBB = jbBB;
    PhyJb::SMtxDGG              nearJbGG = jbGG;
    std::vector<PhyJb::SMtxDGL> nearJbGL = jbGL;
    std::vector<PhyJb::SMtxDGE> nearJbGE = jbGE;
    for (auto&& hit : hits_) {
        // Interaction Local Parameters
        Bool_t  hasLoc  = (cnt_nseg >= 0);
        Bool_t  isInner = (cnt_nseg >= 0 && cnt_nseg < nseg_);
        Short_t itnseg  = (cnt_nseg == nseg_) ? (nseg_ - Numc::ONE<Short_t>) : cnt_nseg;
        Bool_t  hasCxy  = (hit->scx() || hit->scy());

        if (isInner && !hasCxy) {
            ppst = nearPpst;
            jbBB = nearJbBB;
            jbGG = nearJbGG;
            jbGL = nearJbGL;
            jbGE = nearJbGE;
        }
       
        // Propagate
        PhyJb curjb;
        if (isInner) { // Internal Region
            ppst.arg().set_mscat(
                args.at(cnt_nseg).tauu(),
                args.at(cnt_nseg).rhou(),
                args.at(cnt_nseg).taul(),
                args.at(cnt_nseg).rhol());
            ppst.arg().set_eloss(
                args.at(cnt_nseg).elion(),
                args.at(cnt_nseg).elbrm());
        }
        else ppst.arg().clear(); // External Region
        Double_t preIbta = ppst.ibta();
        if (!PropMgnt::PropToZ(hit->cz(), ppst, nullptr, ((hasJacbGlb)?&curjb:nullptr))) break;
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
        if (hasJacbGlb) jbBB = curJbBB * jbBB;
        if (hasJacbGlb) jbGG = curjb.gg() * jbGG;
        if (hasJacbLoc) {
            if (hasLoc && hasCxy) jbGL.at(cnt_nseg) = curjb.gl();
            for (Short_t is = 0; is < cnt_nseg; ++is)
                jbGL.at(is) = curjb.gg() * jbGL.at(is);
            
            if (hasLoc && hasCxy) jbGE.at(cnt_nseg) = curjb.ge();
            for (Short_t is = 0; is < cnt_nseg; ++is)
                jbGE.at(is) = curjb.gg() * jbGE.at(is);
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
            if (hasJacbLoc && hasLoc) {
                for (Short_t is = 0; is <= itnseg; ++is) {
                    for (Short_t it = 0; it < PhyJb::DIML; ++it) {
                        if (hit->scx()) jb(hit->seqIDcx(), DIMG_ + is*PhyJb::DIML+it) += hit->divcx() * jbGL.at(is)(0, it);
                        if (hit->scy()) jb(hit->seqIDcy(), DIMG_ + is*PhyJb::DIML+it) += hit->divcy() * jbGL.at(is)(1, it);
                    }
                }
            } // hasJacbInt
        } // hasCxy

        // TRK
        HitStTRK* hitTRK = Hit<HitStTRK>::Cast(hit);
        if (hitTRK != nullptr) {
            if (hitTRK->sq()) rs(hitTRK->seqIDq()) += hitTRK->nrmq();
            if (hasJacbGlb && hitTRK->sq()) jb(hitTRK->seqIDq(), parIDibta) += hitTRK->divq_ibta() * jbBB;
        }
        
        // TOF
        HitStTOF* hitTOF = Hit<HitStTOF>::Cast(hit);
        if (hitTOF != nullptr) {
            if (hitTOF->sq()) rs(hitTOF->seqIDq()) += hitTOF->nrmq();
            if (hitTOF->st()) rs(hitTOF->seqIDt()) += hitTOF->nrmt();
            if (hasJacbGlb && hitTOF->st()) jb(hitTOF->seqIDt(), parIDibta) += hitTOF->divt_ibta() * jbBB;
            if (hasJacbGlb && hitTOF->sq()) jb(hitTOF->seqIDq(), parIDibta) += hitTOF->divq_ibta() * jbBB;
            if (hasJacbGlb && hitTOF->st() && opt_tsft_) jb(hitTOF->seqIDt(), parIDtsft) += hitTOF->divtsft(); // TOF time shift
        }
        
        // RICH
        HitStRICH* hitRICH = Hit<HitStRICH>::Cast(hit);
        if (hitRICH != nullptr) {
            if (hitRICH->sib()) rs(hitRICH->seqIDib()) += hitRICH->nrmib();
            if (hasJacbGlb && hitRICH->sib()) jb(hitRICH->seqIDib(), parIDibta) += hitRICH->divib_ibta() * jbBB;
        }
        
        // TRD
        HitStTRD* hitTRD = Hit<HitStTRD>::Cast(hit);
        if (hitTRD != nullptr) {
            if (hitTRD->sel()) rs(hitTRD->seqIDel()) += hitTRD->nrmel();
            if (hasJacbGlb && hitTRD->sel()) jb(hitTRD->seqIDel(), parIDeta) += hitTRD->divel_eta() * jbGG(4, 4);
        }

        if (hasCxy) {
            nearPpst = ppst;
            nearJbBB = jbBB;
            nearJbGG = jbGG;
            nearJbGL = jbGL;
            nearJbGE = jbGE;
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
    if (hasJacbLoc) {
        for (Short_t it = 0; it < numOfRes_; ++it) {
        for (Short_t jt = 0; jt < numOfParLoc_; ++jt) {
            if (!Numc::Valid(jb(it, DIMG_ + jt))) jb(it, DIMG_ + jt) = Numc::ZERO<>;
            jacobians[1][it * numOfParLoc_ + jt] = jb(it, DIMG_ + jt); 
        }}
    }

    return true;
}
    

PhySt PhyMuFit::interpolate_to_z(Double_t zcoo) const {
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
    ppst.arg().set_eloss(arg.elion(), arg.elbrm());

    return ppst;
}


MatFld PhyMuFit::get_mat(Double_t zbd1, Double_t zbd2) const {
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


#endif // __TRACKLibs_PhyMuFit_C__
