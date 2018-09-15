#ifndef __TRACKLibs_SimpleTrFit_C__
#define __TRACKLibs_SimpleTrFit_C__


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


namespace TrackSys {
        

SimpleTrFit& SimpleTrFit::operator=(const SimpleTrFit& rhs) {
    if (this != &rhs) {
        dynamic_cast<TrFitPar&>(*this) = dynamic_cast<const TrFitPar&>(rhs);
        succ_     = rhs.succ_;
        part_     = rhs.part_;
        TOFt_sft_ = rhs.TOFt_sft_;
        ndof_     = rhs.ndof_;
        nchi_     = rhs.nchi_;
        quality_  = rhs.quality_;
        ndof_tt_  = rhs.ndof_tt_;
        ndof_cx_  = rhs.ndof_cx_;
        ndof_cy_  = rhs.ndof_cy_;
        ndof_ib_  = rhs.ndof_ib_;
        nchi_tt_  = rhs.nchi_tt_;
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
}


SimpleTrFit::SimpleTrFit(const TrFitPar& fitPar) : TrFitPar(fitPar) {
    SimpleTrFit::clear();
    if (!check_hits()) return;
    ndof_cx_ = (nmes_cx_ > LMTN_CX) ? (nmes_cx_ - LMTN_CX) : 0;
    ndof_cy_ = (nmes_cy_ > LMTN_CY) ? (nmes_cy_ - LMTN_CY) : 0;
    ndof_ib_ = nmes_ib_ - (nmes_TOFt_ >= LMTN_TOF_T);
    ndof_tt_ = (ndof_cx_ + ndof_cy_ + ndof_ib_);
    ndof_.at(0) = ndof_cx_;
    ndof_.at(1) = ndof_cy_ + ndof_ib_;
    if (ndof_.at(0) <= Numc::ONE<Short_t>) { SimpleTrFit::clear(); return; }
    if (ndof_.at(1) <= Numc::ONE<Short_t>) { SimpleTrFit::clear(); return; }

    succ_ = (analyticalFit() ? simpleFit() : false);
    if (succ_) succ_ = advancedSimpleFit();
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
        Double_t chi_tt  = (chi_cx + chi_cy + chi_ib);
        Double_t nchi_tt = (chi_tt / static_cast<Double_t>(ndof_tt_));
        
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
                nchi_cx_ = ((ndof_cx_ > 0) ? (chi_cx / static_cast<Double_t>(ndof_cx_)) : 0);
                nchi_cy_ = ((ndof_cy_ > 0) ? (chi_cy / static_cast<Double_t>(ndof_cy_)) : 0);
                nchi_ib_ = ((ndof_ib_ > 0) ? (chi_ib / static_cast<Double_t>(ndof_ib_)) : 0);
                Double_t nchi_cyib = ((ndof_.at(1) > 0) ? ((chi_cy + chi_ib) / static_cast<Double_t>(ndof_.at(1))) : 0);
                nchi_.at(0)    = nchi_cx_;
                nchi_.at(1)    = nchi_cyib;
                quality_.at(0) = TrFitPar::NormQuality(nchi_.at(0), ndof_.at(0));
                quality_.at(1) = TrFitPar::NormQuality(nchi_.at(1), ndof_.at(1));
                part_    = rltSt;
                isUpdate = true;
                updIter++;
            }
        }
        else { nchi_tt_ = nchi_tt; }

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
   
    //if (!succ) CERR("FAIL. IT %2d %2d (RIG %14.8f MASS %14.8f CHI %14.8f) LAMBDA %14.8f\n", curIter, updIter, part_.rig(), part_.mass(), nchi_tt_, lambda);
    //else       CERR("SUCC. IT %2d %2d (RIG %14.8f MASS %14.8f CHI %14.8f) LAMBDA %14.8f\n", curIter, updIter, part_.rig(), part_.mass(), nchi_tt_, lambda);
    
    return succ;
}


Bool_t SimpleTrFit::advancedSimpleFit() {
    if (!survival_test_and_modify(part_, sw_eloss_)) return false;
    Bool_t opt_tsft = (nmes_TOFt_ >= LMTN_TOF_T);

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
    
    return evolve();
}


Bool_t SimpleTrFit::evolve() {
    if (Numc::EqualToZero(part_.mom())) return false;
    Bool_t opt_int   = sw_mscat_; 
    Bool_t opt_tsft  = (nmes_TOFt_ >= LMTN_TOF_T);
    Short_t numOfPar = PhyJb::DIMG + opt_tsft;

    // Reset TOF Time and Path
    Bool_t resetTOF = true;
    HitStTOF::SetOffsetTime(Numc::ZERO<>);
    HitStTOF::SetOffsetPath(Numc::ZERO<>);
    HitStTOF::SetTimeShiftCorr(opt_tsft);
    
    // TOF time shift
    Double_t TOFt_sft = (opt_tsft ? TOFt_sft_ : Numc::ZERO<>);
    
    // Particle Status
    PhySt ppst(part_);

    // Matrix (Rs, Jb)
    PhyJb::SMtxDGG&& jbGG = SMtxId();
    std::vector<PhyJb::SMtxDGL> jbGL(nseg_);
    
    ceres::Vector&& rsCx  = ceres::Vector::Zero(onlycx_nseq_);
    ceres::Matrix&& jbCWx = ceres::Matrix::Zero(onlycx_nseq_, nseg_ * DIML); // Jacb for interaction 
    ceres::Matrix&& cvCCx = ceres::Matrix::Zero(onlycx_nseq_, onlycx_nseq_); 
    
    ceres::Vector&& rsCy  = ceres::Vector::Zero(onlycy_nseq_);
    ceres::Matrix&& jbCWy = ceres::Matrix::Zero(onlycy_nseq_, nseg_ * DIML); // Jacb for interaction 
    ceres::Matrix&& cvCCy = ceres::Matrix::Zero(onlycy_nseq_, onlycy_nseq_); 
  
    Double_t chi_cx = Numc::ZERO<>;
    Double_t chi_cy = Numc::ZERO<>;
    Double_t chi_ib = Numc::ZERO<>;
   
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
            HitStTOF::SetOffsetTime((ppst.time() - Hit<HitStTOF>::Cast(hit)->orgt()) + TOFt_sft);
            resetTOF = false;
        }
        hit->cal(ppst);
        
        // Update Jacb
        jbGG = curjb.gg() * jbGG;
        if (opt_int) {
            if (hasLoc && hasCxy) jbGL.at(cnt_nseg) = curjb.gl();
            for (Short_t is = 0; is < cnt_nseg; ++is)
                jbGL.at(is) = curjb.gg() * jbGL.at(is);
        }
        
        // Coord
        if (hit->scx()) {
            rsCx(hit->onlycx_seqID()) += (hit->cx() - ppst.cx());
            cvCCx(hit->onlycx_seqID(), hit->onlycx_seqID()) += (hit->ecx() * hit->ecx());
            if (opt_int) {
                for (Short_t is = 0; is <= itnseg; ++is) {
                for (Short_t it = 0; it < PhyJb::DIML; ++it) {
                    if (hit->scx()) jbCWx(hit->onlycx_seqID(), is*DIML+it) += Numc::NEG<> * jbGL.at(is)(0, it);
                }}
            }
        }
        if (hit->scy()) {
            rsCy(hit->onlycy_seqID()) += (hit->cy() - ppst.cy());
            cvCCy(hit->onlycy_seqID(), hit->onlycy_seqID()) += (hit->ecy() * hit->ecy());
            if (opt_int) {
                for (Short_t is = 0; is <= itnseg; ++is) {
                for (Short_t it = 0; it < PhyJb::DIML; ++it) {
                    if (hit->scy()) jbCWy(hit->onlycy_seqID(), is*DIML+it) += Numc::NEG<> * jbGL.at(is)(0, it);
                }}
            }
        }
        
        // TRK
        HitStTRK* hitTRK = Hit<HitStTRK>::Cast(hit);
        if (hitTRK != nullptr) {
            if (hitTRK->sqx()) chi_ib += hitTRK->chiqx() * hitTRK->chiqx();
            if (hitTRK->sqy()) chi_ib += hitTRK->chiqy() * hitTRK->chiqy();
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
    
    cvCCx += (jbCWx * jbCWx.transpose());
    ceres::Matrix&& wtCCx = cvCCx.inverse();
    chi_cx += (rsCx.transpose() * wtCCx * rsCx);
    
    cvCCy += (jbCWy * jbCWy.transpose());
    ceres::Matrix&& wtCCy = cvCCy.inverse();
    chi_cy += (rsCy.transpose() * wtCCy * rsCy);

    Double_t nchi_cx = ((ndof_cx_ > 0) ? (chi_cx / static_cast<Double_t>(ndof_cx_)) : 0);
    Double_t nchi_cy = ((ndof_cy_ > 0) ? (chi_cy / static_cast<Double_t>(ndof_cy_)) : 0);
    Double_t nchi_ib = ((ndof_ib_ > 0) ? (chi_ib / static_cast<Double_t>(ndof_ib_)) : 0);
    
    //Robust robust(Robust::Opt::ON, Numc::FIVE<long double>);
    //long double crrx = robust.minimizer(std::sqrt(nchi_cx)).at(0); crrx = (crrx * crrx);
    //long double crry = robust.minimizer(std::sqrt(nchi_cy)).at(0); crry = (crry * crry);
    //chi_cx *= crrx;
    //chi_cy *= crry;
    //nchi_cx *= crrx;
    //nchi_cy *= crry;
    
    Double_t nchi_cyib = ((ndof_.at(1) > 0) ? ((chi_cy + chi_ib) / static_cast<Double_t>(ndof_.at(1))) : 0);
    if (!Numc::Valid(nchi_cyib) || Numc::Compare(nchi_cyib) <= 0) nchi_cyib = Numc::ZERO<>;
    
    Double_t nchi_tt = ((chi_cx + chi_cy + chi_ib) / static_cast<Double_t>(ndof_tt_));
    if (!Numc::Valid(nchi_tt) || Numc::Compare(nchi_tt) <= 0) nchi_tt = Numc::ZERO<>;

    nchi_tt_ = nchi_tt;
    nchi_cx_ = nchi_cx;
    nchi_cy_ = nchi_cy;
    nchi_ib_ = nchi_ib;
    
    nchi_.at(0)    = nchi_cx_;
    nchi_.at(1)    = nchi_cyib;
    quality_.at(0) = TrFitPar::NormQuality(nchi_.at(0), ndof_.at(0));
    quality_.at(1) = TrFitPar::NormQuality(nchi_.at(1), ndof_.at(1));

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
            HitStTOF::SetOffsetTime((ppst.time() - Hit<HitStTOF>::Cast(hit)->orgt()) + TOFt_sft);
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

    //Double_t chiCxy = std::sqrt(costC / static_cast<Double_t>(nmes_cx_ + nmes_cy_ - Numc::FIVE<Short_t>));
    //Robust robust(Robust::Opt::ON, Numc::FIVE<long double>);
    //std::array<long double, 3>&& mini = robust.minimizer(chiCxy);
    //costC *= (mini.at(0) * mini.at(0));
    //if (hasGrad) grdC *= (mini.at(1) * mini.at(2));

    costG += costC;
    if (hasGrad) grdG += grdC;

    Short_t ndof = (numOfRes_ - numOfPar_);
    cost[0] = (costG / static_cast<Double_t>(ndof));
    if (hasGrad) for (int it = 0; it < numOfPar_; ++it) gradient[it] = (grdG(it) / static_cast<Double_t>(ndof));
    
    return true;
}


PhySt SimpleTrFit::interpolate_to_z(Double_t zcoo) const {
    PhySt nullst = part_; nullst.reset(part_.info());
    if (!succ_) return nullst;
    
    PhySt  ppst = part_;
    Bool_t succ = PropMgnt::PropToZ(zcoo, ppst);
    if (!succ) return nullst;
    return ppst;
}


} // namespace TrackSys


#endif // __TRACKLibs_SimpleTrFit_C__
