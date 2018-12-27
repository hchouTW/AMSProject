#ifndef __TRACKLibs_SimpleTrFit_C__
#define __TRACKLibs_SimpleTrFit_C__


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
#include "SimpleTrFit.h"


namespace TrackSys {
        

SimpleTrFit& SimpleTrFit::operator=(const SimpleTrFit& rhs) {
    if (this != &rhs) {
        dynamic_cast<TrFitPar&>(*this) = dynamic_cast<const TrFitPar&>(rhs);
        succ_ = rhs.succ_;
        part_ = rhs.part_;
        
        args_ = rhs.args_;
        stts_ = rhs.stts_;

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
    }
    return *this;
}


void SimpleTrFit::clear() {
    succ_ = false;
    part_.reset(info_);
    part_.arg().reset(sw_mscat_, sw_eloss_);
 
    args_.clear();
    stts_.clear();

    ndof_.fill(0);
    nchi_.fill(0);
    quality_.fill(0);
    
    ndof_tt_ = 0;
    ndof_cx_ = 0;
    ndof_cy_ = 0;
    ndof_ib_ = 0;
        
    nchi_tt_ = 0;
    nchi_cx_ = 0;
    nchi_cy_ = 0;
    nchi_ib_ = 0;
}


SimpleTrFit::SimpleTrFit(const TrFitPar& fitPar, Bool_t withLocal) : TrFitPar(fitPar) {
    SimpleTrFit::clear();
    if (!check_hits()) { SimpleTrFit::clear(); TrFitPar::clear(); return; }
    if (nmes_cx_ <= LMTN_CX) { SimpleTrFit::clear(); TrFitPar::clear(); return; }
    if (nmes_cy_ <= LMTN_CY) { SimpleTrFit::clear(); TrFitPar::clear(); return; }

    ndof_cx_ = (nmes_cx_ > LMTN_CX) ? (nmes_cx_ - LMTN_CX) : 0;
    ndof_cy_ = (nmes_cy_ > LMTN_CY) ? (nmes_cy_ - LMTN_CY) : 0;
    ndof_ib_ = nmes_ib_ - (nmes_TOFt_ >= LMTN_TOF_T);
    ndof_tt_ = (ndof_cx_ + ndof_cy_ + ndof_ib_);
    
    ndof_.at(0) = ndof_cx_;
    ndof_.at(1) = ndof_cy_ + ndof_ib_;
    if (ndof_.at(0) <= Numc::ONE<Short_t>) { SimpleTrFit::clear(); TrFitPar::clear(); return; }
    if (ndof_.at(1) <= Numc::ONE<Short_t>) { SimpleTrFit::clear(); TrFitPar::clear(); return; }

    timer_.start();
    
    // Initial Status
    succ_ = (analyticalFit() ? simpleFit() : false);
    if (succ_) succ_ = advancedSimpleCooFit();
   
    // Simple Fit
    if (succ_) succ_ = advancedSimpleFit();
    if (withLocal && (succ_ && sw_mscat_)) succ_ = localSimpleFit();
    
    if (succ_) succ_ = evolve();
    if (!succ_) { SimpleTrFit::clear(); TrFitPar::clear(); }

    timer_.stop();

    //if (!succ_) CERR("FAILURE === SimpleTrFit\n");
}


Bool_t SimpleTrFit::analyticalFit() {
    // Global
    Double_t prefit_ortt = ((ortt_ == Orientation::kDownward) ? Numc::NEG<> : Numc::ONE<>);
        
    Int_t firstHitID = -1;
    for (UInt_t ih = 0; ih < hits_.size(); ++ih) {
    if (hits_.at(ih)->scy()) { 
        firstHitID = ih; break; 
    }}
    if (firstHitID < 0) return false;

    // Linear Fit on X
    // Equation of Motion
    // X  = PX + TX * (Z - RefZ)
    // dZ = Z - RefZ
    // UX = TX * UZ
    // Chisq = (X - Xm)^2
    // | PX |   | Sum(1)     Sum(dZ)   |^-1   | Sum(Xm)    |
    // |    | = |                      |    * |            |
    // | TX |   | Sum(dZ)    Sum(dZ^2) |      | Sum(dZ*Xm) |
    Double_t prefit_pz = hits_.at(firstHitID)->cz();
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
        for (UInt_t ih = firstHitID; ih < hits_.size(); ++ih)
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
                Double_t stp = ((static_cast<Double_t>(it) + Numc::ONE_TO_TWO) / static_cast<Double_t>(nstp));
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
            cur_Ae += Numc::ONE_TO_TWO * crs.at(id) * stp.at(id) * stp.at(id);
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
        
    // Set to first Hit
    part_.arg().reset(false, false);
    Bool_t succ = PropMgnt::PropToZ(hits_.at(0)->cz(), part_);
    part_.arg().reset(sw_mscat_, sw_eloss_);
    if (succ) {
        part_.set_state_with_cos(
            part_.cx(), part_.cy(), hits_.at(0)->cz(),
            part_.ux(), part_.uy(), part_.uz()
        );
    } 
    
    return succ;
}


Bool_t SimpleTrFit::simpleFit() {
    // Turn Off (mscat, eloss)
    part_.arg().reset(false, false);
    
    Bool_t succ    = false;
    Bool_t preSucc = false;
    Bool_t curSucc = false;

    Double_t              curLmRhoDen = Numc::ONE<>;
    Double_t              lambda = LAMBDA0;
    PhySt                 rltSt(part_);
    SVecD<PhyJb::DIMG>    curGrdG;
    SMtxSymD<PhyJb::DIMG> curCvGG;

    Short_t updIter = 0;
    Short_t curIter = 0;
    while (curIter <= LMTU_ITER && !succ) {
        Double_t chi_cx = 0;
        Double_t chi_cy = 0;
        
        SVecD<PhyJb::DIMG>    grdG;
        SMtxSymD<PhyJb::DIMG> cvGG;

        Int_t cnt_nhit = 0;
        PhySt ppst(rltSt);
        SMtxD<PhyJb::DIMG>&& ppjb = SMtxId();
        for (auto&& hit : hits_) {
            if (!(hit->scx() || hit->scy())) { cnt_nhit++; continue; }

            // Propagate
            PhyJb curjb;
            if (!PropMgnt::PropToZ(hit->cz(), ppst, nullptr, &curjb)) break;
            ppst.symbk();
       
            // Hit
            hit->cal(ppst);
        
            // Update
            ppjb = curjb.gg() * ppjb;

            // Coord
            SVecD<2>              rsC;
            SMtxD<2, PhyJb::DIMG> jbC;
            if (hit->scx()) {
                rsC(0) = (hit->cx() - ppst.cx()) / hit->ecx();
                for (Short_t it = 0; it < PhyJb::DIMG; ++it)
                    jbC(0, it) += (Numc::NEG<> / hit->ecx()) * ppjb(0, it);
                chi_cx += rsC(0) * rsC(0);
            }
            if (hit->scy()) {
                rsC(1) = (hit->cy() - ppst.cy()) / hit->ecy();
                for (Short_t it = 0; it < PhyJb::DIMG; ++it)
                    jbC(1, it) += (Numc::NEG<> / hit->ecy()) * ppjb(1, it);
                chi_cy += rsC(1) * rsC(1);
            }
            grdG += LA::Transpose(jbC) * rsC;
            cvGG += LA::SimilarityT(jbC, SMtxSymD<2>(SMtxId()));
            
            cnt_nhit++;
        }
        if (cnt_nhit != hits_.size()) break;
        Double_t chi_tt  = (chi_cx + chi_cy);
        Double_t nchi_tt = (chi_tt / static_cast<Double_t>(ndof_cx_ + ndof_cy_));
        
        Bool_t isSucc   = false;
        Bool_t isUpdate = false;
        if (curIter != 0) {
            Double_t lmRho = (nchi_tt_ - nchi_tt) / curLmRhoDen;
            Bool_t   isLmt = (Numc::Compare(lambda, LMTU_LAMBDA) >= 0);
            Double_t convg = std::sqrt(Numc::ONE<> + lambda);
            isSucc = (Numc::Compare(std::fabs((nchi_tt_ - nchi_tt) / (nchi_tt_ + nchi_tt + CONVG_TOLERANCE)) * convg, CONVG_TOLERANCE) <= 0);

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
                nchi_tt_ = nchi_tt;
                part_    = rltSt;
                isUpdate = true;
                updIter++;
            }
        }
        else { nchi_tt_ = nchi_tt; }

        SMtxSymD<PhyJb::DIMG> lmCvGG(cvGG);
        SVecD<PhyJb::DIMG>&&  diagCvGG = (lambda * cvGG.Diagonal());
        lmCvGG.SetDiagonal(SVecD<PhyJb::DIMG>(lmCvGG.Diagonal() + diagCvGG));

        Bool_t isNomag = (Numc::EqualToZero(lmCvGG(4, 4)) && Numc::EqualToZero(grdG(4))); // Fast Check
        if (isNomag) {
            grdG(4) = Numc::ZERO<>;
            SMtxSymD<4>&& lmCvGG_nomag = lmCvGG.Sub<SMtxSymD<4>>(0, 0);
            if (!lmCvGG_nomag.Invert()) break;
            lmCvGG = std::move(SMtxSymD<PhyJb::DIMG>());
            for (Short_t ielem = 0; ielem < 4; ++ielem)
                for (Short_t jelem = ielem; jelem < 4; ++jelem)
                    lmCvGG(ielem, jelem) = lmCvGG_nomag(ielem, jelem);
        }
        else {
            if (!lmCvGG.Invert()) break;
        }

        SVecD<PhyJb::DIMG>&& rslG = (lmCvGG * grdG);
       
        curLmRhoDen = Numc::ZERO<>;
        for (Short_t it = 0; it < PhyJb::DIMG; ++it)
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
    nchi_tt_ = Numc::ZERO<>;

    //if (!succ) CERR("FAIL. IT %2d %2d (RIG %14.8f MASS %14.8f CHI %14.8f) LAMBDA %14.8f\n", curIter, updIter, part_.rig(), part_.mass(), nchi_tt_, lambda);
    //else       CERR("SUCC. IT %2d %2d (RIG %14.8f MASS %14.8f CHI %14.8f) LAMBDA %14.8f\n", curIter, updIter, part_.rig(), part_.mass(), nchi_tt_, lambda);
    
    return succ;
}


Bool_t SimpleTrFit::advancedSimpleCooFit() {
    if (!survival_test_and_modify(part_, sw_eloss_)) return false;

    std::vector<double> params({ part_.cx(), part_.cy(), part_.ux(), part_.uy(), part_.eta() });

    // CeresSolver: Problem
    ceres::GradientProblem problem(new VirtualSimpleTrCooFit(dynamic_cast<TrFitPar&>(*this), part_));
    
    // CeresSolver: Options
    ceres::GradientProblemSolver::Options options;
    options.nonlinear_conjugate_gradient_type = ceres::FLETCHER_REEVES;
    //options.nonlinear_conjugate_gradient_type = ceres::POLAK_RIBIERE;
    //options.nonlinear_conjugate_gradient_type = ceres::HESTENES_STIEFEL;
    options.max_num_iterations = 30;
    //options.max_solver_time_in_seconds = 5.0;

    // CeresSolver: Summary
    ceres::GradientProblemSolver::Summary summary;
    ceres::Solve(options, problem, params.data(), &summary);
    if (!summary.IsSolutionUsable()) return false;
    //if (ceres::NO_CONVERGENCE == summary.termination_type) return false;
    //std::cout << summary.FullReport() << std::endl;

    part_.set_state_with_uxy(params.at(0), params.at(1), part_.cz(), params.at(2), params.at(3), Numc::Compare(part_.uz()));
    part_.set_eta(params.at(4));

    return true;
}


Bool_t SimpleTrFit::advancedSimpleFit() {
    std::vector<double> params({ part_.cx(), part_.cy(), part_.ux(), part_.uy(), part_.eta() });

    // CeresSolver: Problem
    ceres::GradientProblem problem(new VirtualSimpleTrFit(dynamic_cast<TrFitPar&>(*this), part_));
    
    // CeresSolver: Options
    ceres::GradientProblemSolver::Options options;
    options.nonlinear_conjugate_gradient_type = ceres::FLETCHER_REEVES;
    //options.nonlinear_conjugate_gradient_type = ceres::POLAK_RIBIERE;
    //options.nonlinear_conjugate_gradient_type = ceres::HESTENES_STIEFEL;
    options.max_num_iterations = 50;
    //options.max_solver_time_in_seconds = 5.0;

    // CeresSolver: Summary
    ceres::GradientProblemSolver::Summary summary;
    ceres::Solve(options, problem, params.data(), &summary);
    if (!summary.IsSolutionUsable()) return false;
    //if (ceres::NO_CONVERGENCE == summary.termination_type) return false;
    //std::cout << summary.FullReport() << std::endl;
    
    part_.set_state_with_uxy(params.at(0), params.at(1), part_.cz(), params.at(2), params.at(3), Numc::Compare(part_.uz()));
    part_.set_eta(params.at(4));

    return true;
}


Bool_t SimpleTrFit::localSimpleFit() {
    if (Numc::EqualToZero(part_.mom())) return false;
   
    // Interaction Local Parameters
    args_.clear();
    std::vector<double> params(nseg_ * PhyJb::DIML, Numc::ZERO<>);
    
    // CeresSolver: Cost Function
    ceres::CostFunction* cost_function = new VirtualSimpleTrLocFit(dynamic_cast<TrFitPar&>(*this), part_);
    
    // CeresSolver: Problem
    ceres::Problem problem;
    problem.AddResidualBlock(cost_function, nullptr, params.data());
    
    // CeresSolver: Options
    ceres::Solver::Options options;
    options.max_num_iterations = 25;
    //options.max_solver_time_in_seconds = 5.0;

    // CeresSolver: Summary
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    if (!summary.IsSolutionUsable()) return false;
    //if (ceres::NO_CONVERGENCE == summary.termination_type) return false;
    //std::cerr << summary.FullReport() << std::endl;

    args_ = std::move(std::vector<PhyArg>(nseg_, PhyArg(sw_mscat_, sw_eloss_)));
    for (Short_t is = 0; is < nseg_; ++is) {
        args_.at(is).set_mscat(
            params.at(is*PhyJb::DIML+0), 
            params.at(is*PhyJb::DIML+1), 
            params.at(is*PhyJb::DIML+2), 
            params.at(is*PhyJb::DIML+3));
    }
    
    return true;
}


Bool_t SimpleTrFit::evolve() {
    if (Numc::EqualToZero(part_.mom())) return false;
    Bool_t opt_loc   = sw_mscat_; 
    Short_t numOfPar = PhyJb::DIMG;

    // Reset TOF Time and Path
    Bool_t resetTOF = true;
    HitStTOF::SetOffsetPathTime();
    
    // Particle Status
    PhySt ppst(part_);
    ppst.arg().clear();

    // Matrix (Rs, Jb)
    PhyJb::SMtxDGG&& jbGG = SMtxId();
    std::vector<PhyJbMS> jbMS(nseg_);
    
    ceres::Vector&& rsCx  = ceres::Vector::Zero(onlycx_nseq_);                       // Residual
    ceres::Matrix&& cvCCx = ceres::Matrix::Zero(onlycx_nseq_, onlycx_nseq_);         // Convariance
    ceres::Matrix&& jbCWx = ceres::Matrix::Zero(onlycx_nseq_, nseg_ * PhyJbMS::DIM); // Jacb for interaction 
    
    ceres::Vector&& rsCy  = ceres::Vector::Zero(onlycy_nseq_);                       // Residual
    ceres::Matrix&& cvCCy = ceres::Matrix::Zero(onlycy_nseq_, onlycy_nseq_);         // Convariance
    ceres::Matrix&& jbCWy = ceres::Matrix::Zero(onlycy_nseq_, nseg_ * PhyJbMS::DIM); // Jacb for interaction 
    
    Double_t chi_cx = Numc::ZERO<>;
    Double_t chi_cy = Numc::ZERO<>;
    Double_t chi_ib = Numc::ZERO<>;
   
    Short_t cnt_nhit =  0;
    Short_t cnt_nseg = -1;
    PhySt                nearPpst = ppst;
    PhyJb::SMtxDGG       nearJbGG = jbGG;
    std::vector<PhyJbMS> nearJbMS = jbMS;
    for (auto&& hit : hits_) {
        // Interaction Local Parameters
        Bool_t  hasLoc  = (cnt_nseg >= 0);
        Bool_t  isInner = (cnt_nseg >= 0 && cnt_nseg < nseg_);
        Short_t itnseg  = (cnt_nseg == nseg_) ? (nseg_ - Numc::ONE<Short_t>) : cnt_nseg;
        Bool_t  hasCxy  = (hit->scx() || hit->scy());

        if (isInner && !hasCxy) {
            ppst = nearPpst;
            jbGG = nearJbGG;
            jbMS = nearJbMS;
        }
       
        // Propagate
        PhyJb curjb;
        if (!PropMgnt::PropToZ(hit->cz(), ppst, nullptr, &curjb)) break;
        ppst.symbk();
            
        // Hit Status
        if (resetTOF && Hit<HitStTOF>::IsSame(hit)) { // set reference
            HitStTOF::SetOffsetPathTime(
                ppst.path(),
                (ppst.time() - Hit<HitStTOF>::Cast(hit)->orgt())
            );
            resetTOF = false;
        }
        hit->cal(ppst);
        
        // Update Jacb
        jbGG = curjb.gg() * jbGG;
        if (opt_loc) {
            if (hasLoc && hasCxy) jbMS.at(cnt_nseg) = curjb.ms();
            for (Short_t is = 0; is < cnt_nseg; ++is)
                jbMS.at(is).multiplied(curjb.len());
        }
        
        // Coord
        if (hit->scx()) {
            rsCx(hit->onlycx_seqID()) += (hit->cx() - ppst.cx());
            cvCCx(hit->onlycx_seqID(), hit->onlycx_seqID()) += (hit->ecx() * hit->ecx());
            if (opt_loc) {
            for (Short_t is = 0; is <= itnseg; ++is) {
            for (Short_t it = 0; it < PhyJbMS::DIM; ++it) {
                jbCWx(hit->onlycx_seqID(), is*PhyJbMS::DIM+it) += Numc::NEG<> * jbMS.at(is)(0, it);
            }}}
        }
        if (hit->scy()) {
            rsCy(hit->onlycy_seqID()) += (hit->cy() - ppst.cy());
            cvCCy(hit->onlycy_seqID(), hit->onlycy_seqID()) += (hit->ecy() * hit->ecy());
            if (opt_loc) {
            for (Short_t is = 0; is <= itnseg; ++is) {
            for (Short_t it = 0; it < PhyJbMS::DIM; ++it) {
                jbCWy(hit->onlycy_seqID(), is*PhyJbMS::DIM+it) += Numc::NEG<> * jbMS.at(is)(0, it);
            }}}
        }

        // TRK
        HitStTRK* hitTRK = Hit<HitStTRK>::Cast(hit);
        if (hitTRK != nullptr) {
            if (hitTRK->sq()) chi_ib += hitTRK->chiq() * hitTRK->chiq();
        }
        
        // TOF
        HitStTOF* hitTOF = Hit<HitStTOF>::Cast(hit);
        if (hitTOF != nullptr) {
            if (hitTOF->st()) chi_ib += hitTOF->chit() * hitTOF->chit();
            if (hitTOF->sq()) chi_ib += hitTOF->chiq() * hitTOF->chiq();
        }
        
        // RICH
        HitStRICH* hitRICH = Hit<HitStRICH>::Cast(hit);
        if (hitRICH != nullptr) {
            if (hitRICH->sib()) chi_ib += hitRICH->chiib() * hitRICH->chiib();
        }
        
        // TRD
        HitStTRD* hitTRD = Hit<HitStTRD>::Cast(hit);
        if (hitTRD != nullptr) {
            if (hitTRD->sel()) chi_ib += hitTRD->chiel() * hitTRD->chiel();
        }
        
        if (hasCxy) {
            nearPpst = ppst;
            nearJbGG = jbGG;
            nearJbMS = jbMS;
            cnt_nseg++;
        }
        cnt_nhit++;
    }
    if (cnt_nhit != hits_.size()) return false;
    if (cnt_nseg != nseg_) return false;
    
    cvCCx += (jbCWx * jbCWx.transpose());
    ceres::Matrix&& wtCCx = cvCCx.inverse();
    chi_cx += (rsCx.transpose() * wtCCx * rsCx);
    
    cvCCy += (jbCWy * jbCWy.transpose());
    ceres::Matrix&& wtCCy = cvCCy.inverse();
    chi_cy += (rsCy.transpose() * wtCCy * rsCy);

    Double_t nchi_cx = ((ndof_cx_ > 0) ? (chi_cx / static_cast<Double_t>(ndof_cx_)) : 0);
    Double_t nchi_cy = ((ndof_cy_ > 0) ? (chi_cy / static_cast<Double_t>(ndof_cy_)) : 0);
    Double_t nchi_ib = ((ndof_ib_ > 0) ? (chi_ib / static_cast<Double_t>(ndof_ib_)) : 0);
   
    Double_t nchi_cyib = ((ndof_.at(1) > 0) ? ((chi_cy + chi_ib) / static_cast<Double_t>(ndof_.at(1))) : 0);
    if (!Numc::Valid(nchi_cyib) || Numc::Compare(nchi_cyib) <= 0) nchi_cyib = Numc::ZERO<>;
   
    Double_t ndof_tt = static_cast<Double_t>(ndof_tt_);
    Double_t nchi_tt = ((chi_cx + chi_cy + chi_ib) / ndof_tt);
    if (!Numc::Valid(nchi_tt) || Numc::Compare(nchi_tt) <= 0) nchi_tt = Numc::ZERO<>;

    nchi_tt_ = nchi_tt;
    nchi_cx_ = nchi_cx;
    nchi_cy_ = nchi_cy;
    nchi_ib_ = nchi_ib;
    
    nchi_.at(0)    = nchi_cx;
    nchi_.at(1)    = nchi_cyib;
    quality_.at(0) = Numc::NormQuality(nchi_.at(0), ndof_.at(0));
    quality_.at(1) = Numc::NormQuality(nchi_.at(1), ndof_.at(1));
    
    // Local States
    stts_.clear();
    Bool_t succLocal = false;
    if (args_.size() == nseg_) {
        PhySt ppstts(part_);
        std::vector<PhySt> stts;

        Short_t cnt_nseg = -1;
        for (auto&& hit : hits_) {
            if (!(hit->scx() || hit->scy())) continue;
            if (cnt_nseg != -1) ppstts.arg() = args_.at(cnt_nseg);
            else                ppstts.arg().clear();
            if (!PropMgnt::PropToZ(hit->cz(), ppstts)) break;
            ppstts.symbk();
            stts.push_back(ppstts);
            cnt_nseg++;
        }
        succLocal = (cnt_nseg == nseg_ && stts.size() == (nseg_ + 1));
        
        if (!succLocal) { args_.clear(); stts_.clear(); return false; }
        else {
            stts_ = stts;
            for (auto&& stt : stts_) stt.arg().clear();
            for (UInt_t it = 0; it < args_.size(); ++it) {
                PhyArg& arg = args_.at(it);
                stts_.at(it).arg().set_mscat(arg.tauu(), arg.rhou(), arg.taul(), arg.rhol());
            }
            if (!(hits_.at(0)->scx() || hits_.at(0)->scy())) stts_.insert(stts_.begin(), part_);
        }
    }

    return true;
}


bool VirtualSimpleTrCooFit::Evaluate(const double* parameters, double* cost, double* gradient) const {
    if (numOfRes_ <= 0 || numOfPar_ <= 0) return false;
    cost[0] = Numc::ZERO<>;
    
    Bool_t hasGrd = (gradient != nullptr);
    if (hasGrd) std::fill_n(gradient, numOfPar_, Numc::ZERO<>);

    // Particle Status
    PhySt ppst(part_);
    ppst.arg().clear();
    ppst.set_state_with_uxy(parameters[0], parameters[1], part_.cz(), parameters[2], parameters[3], Numc::Compare(part_.uz()));
    ppst.set_eta(parameters[4]);
    
    // Matrix (Rs, Jb)
    PhyJb::SMtxDGG&& jbGG = SMtxId();
    std::vector<PhyJb::SMtxDGL> jbGL(nseg_);
    
    ceres::Vector&& rsC  = ceres::Vector::Zero(onlyc_nseq_);                    // Residual
    ceres::Matrix&& cvCC = ceres::Matrix::Zero(onlyc_nseq_, onlyc_nseq_);       // Convariance
    ceres::Matrix&& jbCG = ceres::Matrix::Zero(onlyc_nseq_, numOfPar_);         // Jacb for coord
    ceres::Matrix&& jbCL = ceres::Matrix::Zero(onlyc_nseq_, nseg_*PhyJb::DIML); // Jacb for interaction 
    
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
        if (!(hit->scx() || hit->scy())) { cnt_nhit++; continue; }

        if (isInner) {
            ppst = nearPpst;
            jbGG = nearJbGG;
            jbGL = nearJbGL;
        }
       
        // Propagate
        PhyJb curjb;
        if (!PropMgnt::PropToZ(hit->cz(), ppst, nullptr, &curjb)) break;
        ppst.symbk();
            
        // Hit Status
        hit->cal(ppst);
        
        // Update Jacb
        jbGG = curjb.gg() * jbGG;
        if (opt_loc_) {
            if (hasLoc) jbGL.at(cnt_nseg) = curjb.gl();
            for (Short_t is = 0; is < cnt_nseg; ++is)
                jbGL.at(is) = curjb.gg() * jbGL.at(is);
        }
        
        // Coord
        if (hit->scx()) {
            rsC(hit->onlyc_seqIDcx()) += (hit->cx() - ppst.cx());
            cvCC(hit->onlyc_seqIDcx(), hit->onlyc_seqIDcx()) += (hit->ecx() * hit->ecx());
            if (hasGrd) {
            for (Short_t it = 0; it < PhyJb::DIMG; ++it) {
                jbCG(hit->onlyc_seqIDcx(), it) += Numc::NEG<> * jbGG(0, it);
            }}
            if (opt_loc_) {
            for (Short_t is = 0; is <= itnseg; ++is) {
            for (Short_t it = 0; it < PhyJb::DIML; ++it) {
                jbCL(hit->onlyc_seqIDcx(), is*PhyJb::DIML+it) += Numc::NEG<> * jbGL.at(is)(0, it);
            }}}
        }
        if (hit->scy()) {
            rsC(hit->onlyc_seqIDcy()) += (hit->cy() - ppst.cy());
            cvCC(hit->onlyc_seqIDcy(), hit->onlyc_seqIDcy()) += (hit->ecy() * hit->ecy());
            if (hasGrd) {
            for (Short_t it = 0; it < PhyJb::DIMG; ++it) {
                jbCG(hit->onlyc_seqIDcy(), it) += Numc::NEG<> * jbGG(1, it);
            }}
            if (opt_loc_) {
            for (Short_t is = 0; is <= itnseg; ++is) {
            for (Short_t it = 0; it < PhyJb::DIML; ++it) {
                jbCL(hit->onlyc_seqIDcy(), is*PhyJb::DIML+it) += Numc::NEG<> * jbGL.at(is)(1, it);
            }}}
        }
        
        nearPpst = ppst;
        nearJbGG = jbGG;
        nearJbGL = jbGL;
        cnt_nseg++;
        cnt_nhit++;
    }
    if (cnt_nhit != hits_.size()) return false;
    if (cnt_nseg != nseg_) return false;
    
    cvCC += (jbCL * jbCL.transpose());
    ceres::Matrix&& wtCC = cvCC.inverse();
    Double_t costC = (rsC.transpose() * wtCC * rsC);
    
    Double_t costG = costC / static_cast<Double_t>(numOfDof_);
    ceres::Vector&& grdG = ceres::Vector::Zero(numOfPar_);
    ceres::Matrix&& hesG = ceres::Matrix::Zero(numOfPar_, numOfPar_);
    if (hasGrd) grdG += (jbCG.transpose() * wtCC * rsC);
    if (hasGrd) hesG += (jbCG.transpose() * wtCC * jbCG);

    if (!Numc::Valid(costG) || Numc::Compare(costG) <= 0) costG = Numc::ZERO<>;
    if (hasGrd) {
    for (int it = 0; it < numOfPar_; ++it) {
        if (!Numc::Valid(grdG(it))) grdG(it) = Numc::ZERO<>;
    }}
    if (hasGrd) {
    for (int it = 0; it < numOfPar_; ++it) {
    for (int jt = 0; jt < numOfPar_; ++jt) {
        if (!Numc::Valid(hesG(it, jt))) hesG(it, jt) = Numc::ZERO<>;
    }}}
    
    // Newton Method
    ceres::Vector&& grdGnwt = ceres::Vector::Zero(numOfPar_);
    if (hasGrd) {
        grdGnwt = hesG.inverse() * grdG;
        for (int it = 0; it < numOfPar_; ++it) {
            if (!Numc::Valid(grdGnwt(it))) grdGnwt(it) = Numc::ZERO<>;
        }
    }

    cost[0] = costG;
    if (hasGrd) {
    for (int it = 0; it < numOfPar_; ++it) {
        gradient[it] = grdGnwt(it);
    }}
    
    return true;
}


bool VirtualSimpleTrFit::Evaluate(const double* parameters, double* cost, double* gradient) const {
    if (numOfRes_ <= 0 || numOfPar_ <= 0) return false;
    cost[0] = Numc::ZERO<>;
    
    Bool_t hasGrd = (gradient != nullptr);
    if (hasGrd) std::fill_n(gradient, numOfPar_, Numc::ZERO<>);
    
    // Reset TOF Time and Path
    Bool_t resetTOF = true;
    HitStTOF::SetOffsetPathTime();

    // Particle Status
    PhySt ppst(part_);
    ppst.arg().clear();
    ppst.set_state_with_uxy(parameters[0], parameters[1], part_.cz(), parameters[2], parameters[3], Numc::Compare(part_.uz()));
    ppst.set_eta(parameters[4]);
    
    // Matrix (Rs, Jb)
    PhyJb::SMtxDGG&& jbGG = SMtxId();
    std::vector<PhyJb::SMtxDGL> jbGL(nseg_);
    
    ceres::Vector&& rsC  = ceres::Vector::Zero(onlyc_nseq_);                    // Residual
    ceres::Matrix&& cvCC = ceres::Matrix::Zero(onlyc_nseq_, onlyc_nseq_);       // Convariance
    ceres::Matrix&& jbCG = ceres::Matrix::Zero(onlyc_nseq_, numOfPar_);         // Jacb for coord
    ceres::Matrix&& jbCL = ceres::Matrix::Zero(onlyc_nseq_, nseg_*PhyJb::DIML); // Jacb for interaction 
    
    Double_t        costIb = Numc::ZERO<>;
    ceres::Vector&& grdIb  = ceres::Vector::Zero(numOfPar_);
    ceres::Matrix&& hesIb  = ceres::Matrix::Zero(numOfPar_, numOfPar_);

    Short_t  cnt_nhit =  0;
    Short_t  cnt_nseg = -1;
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
            
        // Hit Status
        if (resetTOF && Hit<HitStTOF>::IsSame(hit)) { // set reference
            HitStTOF::SetOffsetPathTime(
                ppst.path(),
                (ppst.time() - Hit<HitStTOF>::Cast(hit)->orgt())
            );
            resetTOF = false;
        }
        hit->cal(ppst);
        
        // Update Jacb
        jbGG = curjb.gg() * jbGG;
        if (opt_loc_) {
            if (hasLoc && hasCxy) jbGL.at(cnt_nseg) = curjb.gl();
            for (Short_t is = 0; is < cnt_nseg; ++is)
                jbGL.at(is) = curjb.gg() * jbGL.at(is);
        }
        
        // Coord
        if (hit->scx()) {
            rsC(hit->onlyc_seqIDcx()) += (hit->cx() - ppst.cx());
            cvCC(hit->onlyc_seqIDcx(), hit->onlyc_seqIDcx()) += (hit->ecx() * hit->ecx());
            if (hasGrd) {
            for (Short_t it = 0; it < PhyJb::DIMG; ++it) {
                jbCG(hit->onlyc_seqIDcx(), it) += Numc::NEG<> * jbGG(0, it);
            }}
            if (opt_loc_) {
            for (Short_t is = 0; is <= itnseg; ++is) {
            for (Short_t it = 0; it < PhyJb::DIML; ++it) {
                jbCL(hit->onlyc_seqIDcx(), is*PhyJb::DIML+it) += Numc::NEG<> * jbGL.at(is)(0, it);
            }}}
        }
        if (hit->scy()) {
            rsC(hit->onlyc_seqIDcy()) += (hit->cy() - ppst.cy());
            cvCC(hit->onlyc_seqIDcy(), hit->onlyc_seqIDcy()) += (hit->ecy() * hit->ecy());
            if (hasGrd) {
            for (Short_t it = 0; it < PhyJb::DIMG; ++it) {
                jbCG(hit->onlyc_seqIDcy(), it) += Numc::NEG<> * jbGG(1, it);
            }}
            if (opt_loc_) {
            for (Short_t is = 0; is <= itnseg; ++is) {
            for (Short_t it = 0; it < PhyJb::DIML; ++it) {
                jbCL(hit->onlyc_seqIDcy(), is*PhyJb::DIML+it) += Numc::NEG<> * jbGL.at(is)(1, it);
            }}}
        }
        
        // TRK
        HitStTRK* hitTRK = Hit<HitStTRK>::Cast(hit);
        if (hitTRK != nullptr) {
            if (hitTRK->sq()) costIb += hitTRK->chiq() * hitTRK->chiq();
            if (hasGrd && hitTRK->sq()) grdIb(4) += (hitTRK->divq_eta() * jbGG(4, 4)) * hitTRK->nrmq();
            if (hasGrd && hitTRK->sq()) hesIb(4, 4) += (hitTRK->divq_eta() * jbGG(4, 4)) * (hitTRK->divq_eta() * jbGG(4, 4));
        }
        
        // TOF
        HitStTOF* hitTOF = Hit<HitStTOF>::Cast(hit);
        if (hitTOF != nullptr) {
            if (hitTOF->st()) costIb += hitTOF->chit() * hitTOF->chit();
            if (hitTOF->sq()) costIb += hitTOF->chiq() * hitTOF->chiq();
            
            if (hasGrd && hitTOF->st()) grdIb(4) += (hitTOF->divt_eta() * jbGG(4, 4)) * hitTOF->nrmt();
            if (hasGrd && hitTOF->sq()) grdIb(4) += (hitTOF->divq_eta() * jbGG(4, 4)) * hitTOF->nrmq();

            if (hasGrd && hitTOF->st()) hesIb(4, 4) += (hitTOF->divt_eta() * jbGG(4, 4)) * (hitTOF->divt_eta() * jbGG(4, 4));
            if (hasGrd && hitTOF->sq()) hesIb(4, 4) += (hitTOF->divq_eta() * jbGG(4, 4)) * (hitTOF->divq_eta() * jbGG(4, 4));
        }
        
        // RICH
        HitStRICH* hitRICH = Hit<HitStRICH>::Cast(hit);
        if (hitRICH != nullptr) {
            if (hitRICH->sib()) costIb += hitRICH->chiib() * hitRICH->chiib();
            if (hasGrd && hitRICH->sib()) grdIb(4) += (hitRICH->divib_eta() * jbGG(4, 4)) * hitRICH->nrmib();
            if (hasGrd && hitRICH->sib()) hesIb(4, 4) += (hitRICH->divib_eta() * jbGG(4, 4)) * (hitRICH->divib_eta() * jbGG(4, 4));
        }
        
        // TRD
        HitStTRD* hitTRD = Hit<HitStTRD>::Cast(hit);
        if (hitTRD != nullptr) {
            if (hitTRD->sel()) costIb += hitTRD->chiel() * hitTRD->chiel();
            if (hasGrd && hitTRD->sel()) grdIb(4) += (hitTRD->divel_eta() * jbGG(4, 4)) * hitTRD->nrmel();
            if (hasGrd && hitTRD->sel()) hesIb(4, 4) += (hitTRD->divel_eta() * jbGG(4, 4)) * (hitTRD->divel_eta() * jbGG(4, 4));
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
    
    Double_t costC = Numc::ZERO<>;
    ceres::Vector&& grdC = ceres::Vector::Zero(numOfPar_);
    ceres::Matrix&& hesC = ceres::Matrix::Zero(numOfPar_, numOfPar_);
    
    cvCC += (jbCL * jbCL.transpose());
    ceres::Matrix&& wtCC = cvCC.inverse();
    costC += (rsC.transpose() * wtCC * rsC);
    if (hasGrd) grdC += (jbCG.transpose() * wtCC * rsC);
    if (hasGrd) hesC += (jbCG.transpose() * wtCC * jbCG);

    Double_t        costG = Numc::ZERO<>;
    ceres::Vector&& grdG  = ceres::Vector::Zero(numOfPar_);
    ceres::Matrix&& hesG  = ceres::Matrix::Zero(numOfPar_, numOfPar_);
    costG = (costC + costIb) / static_cast<Double_t>(numOfDof_);
    if (hasGrd) grdG = (grdC + grdIb);
    if (hasGrd) hesG = (hesC + hesIb);
    if (!Numc::Valid(costG) || Numc::Compare(costG) <= 0) costG = Numc::ZERO<>;
    if (hasGrd) {
    for (int it = 0; it < numOfPar_; ++it) {
        if (!Numc::Valid(grdG(it))) grdG(it) = Numc::ZERO<>;
    }}
    if (hasGrd) {
    for (int it = 0; it < numOfPar_; ++it) {
    for (int jt = 0; jt < numOfPar_; ++jt) {
        if (!Numc::Valid(hesG(it, jt))) hesG(it, jt) = Numc::ZERO<>;
    }}}
    
    // Newton Method
    ceres::Vector&& grdGnwt = ceres::Vector::Zero(numOfPar_);
    if (hasGrd) {
        grdGnwt = hesG.inverse() * grdG;
        for (int it = 0; it < numOfPar_; ++it) {
            if (!Numc::Valid(grdGnwt(it))) grdGnwt(it) = Numc::ZERO<>;
        }
    }

    cost[0] = costG;
    if (hasGrd) {
    for (int it = 0; it < numOfPar_; ++it) {
        gradient[it] = grdGnwt(it);
    }}
    
    return true;
}


bool VirtualSimpleTrLocFit::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const {
    if (numOfRes_ <= 0 || numOfPar_ <= 0) return false;
    std::fill_n(residuals, numOfRes_, Numc::ZERO<>);
    
    Bool_t hasJacb = (jacobians != nullptr && jacobians[0] != nullptr);
    if (hasJacb) std::fill_n(jacobians[0], numOfRes_ * numOfPar_, Numc::ZERO<>);
    
    // Particle Status
    PhySt ppst(part_);
    ppst.arg().clear();
    
    // Interaction Local Parameters
    std::vector<PhyArg> args(nseg_, PhyArg(sw_mscat_, sw_eloss_));
    for (Short_t is = 0; is < nseg_; ++is) {
        args.at(is).set_mscat(
            parameters[0][is*PhyJb::DIML+0], 
            parameters[0][is*PhyJb::DIML+1], 
            parameters[0][is*PhyJb::DIML+2], 
            parameters[0][is*PhyJb::DIML+3]);
    }
    
    // Matrix (Rs, Jb)
    PhyJb::SMtxDGG&& jbGG = SMtxId();
    std::vector<PhyJb::SMtxDGL> jbGL(nseg_);

    ceres::Vector rs = ceres::Vector::Zero(numOfRes_);
    ceres::Matrix jb = ceres::Matrix::Zero(numOfRes_, numOfPar_);
    
    // Interaction Local Parameters
    for (Short_t is = 0; is < nseg_; ++is) {
        SVecD<6> inrm, idiv;
        if (hasJacb) args.at(is).cal_nrm_and_div(inrm, idiv);
        else         args.at(is).cal_nrm(inrm);
            
        for (Short_t it = 0; it < PhyJb::DIML; ++it)
            rs(onlyc_nseq_ + is*PhyJb::DIML+it) += inrm(it);

        if (hasJacb) {
        for (Short_t it = 0; it < PhyJb::DIML; ++it) {
            jb(onlyc_nseq_ + is*PhyJb::DIML+it, is*PhyJb::DIML+it) += idiv(it);
        }}
    }
    
    Short_t dof1st   =  0;
    Short_t cnt_nhit =  0;
    Short_t cnt_nseg = -1;
    for (auto&& hit : hits_) {
        // Interaction Local Parameters
        Bool_t  hasLoc  = (cnt_nseg >= 0);
        Bool_t  isInner = (cnt_nseg >= 0 && cnt_nseg < nseg_);
        Short_t itnseg  = (cnt_nseg == nseg_) ? (nseg_ - Numc::ONE<Short_t>) : cnt_nseg;
        if (!(hit->scx() || hit->scy())) { cnt_nhit++; continue; }
        if (dof1st != 0) dof1st = (hit->scx() + hit->scy());

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
        if (!PropMgnt::PropToZ(hit->cz(), ppst, nullptr, ((hasJacb)?&curjb:nullptr))) break;
        ppst.symbk();
       
        // Hit Status
        hit->cal(ppst);
       
        // Update
        jbGG = curjb.gg() * jbGG;
        if (hasLoc) jbGL.at(cnt_nseg) = curjb.gl();
        for (Short_t is = 0; is < cnt_nseg; ++is)
            jbGL.at(is) = curjb.gg() * jbGL.at(is);
        
        // Coord
        if (hit->scx()) {
            Double_t resx = ((hit->cx() - ppst.cx()) / hit->ecx());
            Double_t divx = (Numc::NEG<> / hit->ecx());
            rs(hit->onlyc_seqIDcx()) += resx;
            if (hasJacb) {
                for (Short_t is = 0; is <= itnseg; ++is) {
                for (Short_t it = 0; it < PhyJb::DIML; ++it) {
                    jb(hit->onlyc_seqIDcx(), is*PhyJb::DIML+it) += divx * jbGL.at(is)(0, it);
                }}
            }
        }
        if (hit->scy()) {
            Double_t resy = ((hit->cy() - ppst.cy()) / hit->ecy());
            Double_t divy = (Numc::NEG<> / hit->ecy());
            rs(hit->onlyc_seqIDcy()) += resy;
            if (hasJacb) {
                for (Short_t is = 0; is <= itnseg; ++is) {
                for (Short_t it = 0; it < PhyJb::DIML; ++it) {
                    jb(hit->onlyc_seqIDcy(), is*PhyJb::DIML+it) += divy * jbGL.at(is)(1, it);
                }}
            }
        }
       
        cnt_nseg++;
        cnt_nhit++;
    }
    if (cnt_nhit != hits_.size()) return false;
    if (cnt_nseg != nseg_) return false;
        
    for (int it = 0; it < numOfRes_; ++it) {
        if (!Numc::Valid(rs(it))) rs(it) = Numc::ZERO<>;
        residuals[it] = rs[it];
    }
    
    if (hasJacb) {
    for (int it = 0; it < numOfRes_; ++it) {
    for (int jt = 0; jt < numOfPar_; ++jt) {
        if (!Numc::Valid(jb(it, jt))) jb(it, jt) = Numc::ZERO<>;
        jacobians[0][it*numOfPar_+jt] = jb(it, jt);
    }}}

    return true;
}


PhySt SimpleTrFit::interpolate_to_z(Double_t zcoo) const {
    PhySt nullst = part_; nullst.reset(part_.info());
    if (!succ_) return nullst;
    
    Bool_t isLocal = (args_.size() == nseg_ && stts_.size() != 0);
    if (!isLocal) {
        PhySt  ppst(part_);
        Bool_t succ = PropMgnt::PropToZ(zcoo, ppst);
        if (!succ) return nullst;
        return ppst;
    }
    else {
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

    return nullst;
}


} // namespace TrackSys


#endif // __TRACKLibs_SimpleTrFit_C__
