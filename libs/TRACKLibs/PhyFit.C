#ifndef __TRACKLibs_PhyFit_C__
#define __TRACKLibs_PhyFit_C__


namespace TrackSys {

    
TrFitPar::TrFitPar(const PartType& type, const Orientation& ortt, Bool_t sw_mscat, Bool_t sw_eloss) {
    clear();

    sw_mscat_ = sw_mscat;
    sw_eloss_ = sw_eloss;
    type_ = type;
    ortt_ = ortt;
}


void TrFitPar::clear() {
    sw_mscat_ = false;
    sw_eloss_ = false;
    
    type_ = PartType::Proton;
    ortt_ = Orientation::kDownward;
    hits_.clear();
    nhtx_ = 0;
    nhty_ = 0;

    is_check_ = false;
}


Bool_t TrFitPar::checkHit() {
    if (is_check_) return true;
    if (hits_.size() == 0) return false;
    if (ortt_ == Orientation::kDownward) HitSt::Sort(hits_, HitSt::Orientation::kDownward);
    else                                 HitSt::Sort(hits_, HitSt::Orientation::kUpward);
    
    Short_t nx = 0, ny = 0;
    for (auto&& hit : hits_) {
        if (hit.sx()) nx++;
        if (hit.sy()) ny++;
    }
    if (nx < LMTL_NHIT_X && ny < LMTL_NHIT_Y) return false;
   
    Short_t seq = 0;
    for (auto&& hit : hits_) {
        hit.set_err(type_);
        hit.set_seqID(seq);
        seq += (hit.sx() + hit.sy());
    }

    nhtx_ = nx;
    nhty_ = ny;

    is_check_ = true;
    return true;
}
        

SimpleTrFit::SimpleTrFit(TrFitPar& fitPar) : TrFitPar(fitPar) {
    SimpleTrFit::clear();
    if (!checkHit()) return;
   
    ndfx_ = nhtx_ - 2;
    ndfy_ = nhty_ - 3;
    ndof_ = ndfx_ + ndfy_;
    succ_ = (analyticalFit() ? simpleFit() : false);
    if (!succ_) SimpleTrFit::clear();
}


void SimpleTrFit::clear() {
    succ_ = false;
    part_.reset(type_);
    part_.arg().reset(false, false);

    ndfx_ = 0;
    ndfy_ = 0;
    chix_ = 0;
    chiy_ = 0;

    ndof_ = 0;
    nchi_ = 0;
}


Bool_t SimpleTrFit::analyticalFit() {
    // Linear Fit on X
    // Equation of Motion
    // X  = PX + TX * (Z - RefZ)
    // dZ = Z - RefZ
    // UX = TX * UZ
    // Chisq = (X - Xm)^2
    // | PX |   | Sum(1)     Sum(dZ)   |^-1   | Sum(Xm)    |
    // |    | = |                      |    * |            |
    // | TX |   | Sum(dZ)    Sum(dZ^2) |      | Sum(dZ*Xm) |
    Double_t prefit_pz = hits_.at(0).cz();
    Double_t prefit_px = MGMath::ZERO;
    Double_t prefit_tx = MGMath::ZERO;
    {
        SMtxSymD<2> mtx;
        SVecD<2>    res;
        for (auto&& hit : hits_) {
            if (!hit.sx()) continue;
            Double_t ex  = hit.ex();
            Double_t err = (MGMath::ONE / ex / ex);
            Double_t dz1 = hit.cz() - prefit_pz;
            Double_t dz2 = dz1 * dz1;
            mtx(0, 0) += err * MGMath::ONE;
            mtx(0, 1) += err * dz1;
            mtx(1, 1) += err * dz2;
            res(0)    += err * hit.cx();
            res(1)    += err * dz1 * hit.cx();
        }
        if (!mtx.Invert()) return false;
        SVecD<2>&& rsl = mtx * res;
        prefit_px = rsl(0);
        prefit_tx = rsl(1);
        
        for (auto&& hit : hits_) {
            if (hit.sx()) continue;
            Double_t dz = hit.cz() - prefit_pz;
            Double_t px = prefit_px + prefit_tx * dz;
            hit.set_dummy_x(px);
        }
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
    Double_t prefit_py = MGMath::ZERO;
    Double_t prefit_uy = MGMath::ZERO;
    Double_t prefit_ea = 0.001;
    {
        const Double_t PROP_FACT = 2.99792458e-04;
        Double_t Lambda = PROP_FACT * part_.info().chrg_to_mass(); 

        std::vector<Double_t> stp(hits_.size(), MGMath::ZERO);
        std::vector<Double_t> crs(hits_.size(), MGMath::ZERO);

        const Int_t nstp = 3;
        for (Int_t ih = 1; ih < hits_.size(); ++ih) {
            SVecD<3>&& ref_l = (hits_.at(ih).c() - hits_.at(ih-1).c());
            SVecD<3>&& ref_u = LA::Unit(ref_l);
            Double_t   ref_s = LA::Mag(ref_l);
            
            SVecD<3> mfldv;
            for (Int_t it = 0; it < nstp; ++it) {
                Double_t stp = ((static_cast<Double_t>(it) + MGMath::HALF) / static_cast<Double_t>(nstp));
                SVecD<3>&& ref_m = ((MGMath::ONE - stp) * hits_.at(ih-1).c() + stp * hits_.at(ih).c());
                MagFld&&   mfld  = MagMgnt::Get(ref_m);
                mfldv += mfld();
            }
            
            mfldv /= static_cast<Double_t>(nstp);
            Double_t mucrs = Lambda * (ref_u(2) * mfldv(0) - ref_u(0) * mfldv(2));

            stp.at(ih) = ref_s;
            crs.at(ih) = mucrs;
        }
        
        SMtxSymD<3> mtx;
        SVecD<3>    res;
        Double_t    cur_Au = MGMath::ZERO;
        Double_t    cur_Ae = MGMath::ZERO;
        for (Int_t ih = 0; ih < hits_.size(); ++ih) {
            HitSt& hit = hits_.at(ih);
            Double_t ey  = hit.ey();
            Double_t err = (MGMath::ONE / ey / ey);
            
            cur_Au += stp.at(ih);
            cur_Ae += MGMath::HALF * crs.at(ih) * stp.at(ih) * stp.at(ih);
            for (Int_t jh = 0; jh < ih; ++jh)
                cur_Ae += crs.at(jh) * stp.at(jh) * stp.at(ih);

            mtx(0, 0) += err * MGMath::ONE;
            mtx(0, 1) += err * cur_Au;
            mtx(0, 2) += err * cur_Ae;
            mtx(1, 1) += err * cur_Au * cur_Au;
            mtx(1, 2) += err * cur_Au * cur_Ae;
            mtx(2, 2) += err * cur_Ae * cur_Ae;
            res(0)    += err * hit.cy();
            res(1)    += err * cur_Au * hit.cy();
            res(2)    += err * cur_Ae * hit.cy();
        }

        if (!mtx.Invert()) return false;
        SVecD<3>&& rsl = mtx * res;
        prefit_py = rsl(0);
        prefit_uy = rsl(1);
        prefit_ea = rsl(2);
    }
    
   
    // Merge Fitting Result
    Double_t prefit_ortt = ((ortt_ == Orientation::kDownward) ? MGMath::NEG : MGMath::ONE);
    Double_t prefit_uz = prefit_ortt * std::fabs((MGMath::ONE - prefit_uy * prefit_uy) / (MGMath::ONE + prefit_tx * prefit_tx));
    Double_t prefit_ux = prefit_tx * prefit_uz;
    part_.set_state_with_cos(
        prefit_px, prefit_py, prefit_pz,
        prefit_ux, prefit_uy, prefit_uz
    );
    part_.set_eta(prefit_ea);
    
    return true;
}


Bool_t SimpleTrFit::simpleFit() {
    Bool_t succ    = false;
    Bool_t preSucc = false;
    Bool_t curSucc = false;

    Double_t    curLmRhoDen = MGMath::ONE;
    Double_t    lambda = LAMBDA0;
    PhySt       rltSt(part_);
    SVecD<5>    curGrdG;
    SMtxSymD<5> curCvGG;

    Int_t  updIter = 0;
    Int_t  curIter = 0;
    while (curIter <= LMTU_ITER && !succ) {
        Double_t chix = MGMath::ZERO;
        Double_t chiy = MGMath::ZERO;
        SVecD<5>    grdG;
        SMtxSymD<5> cvGG;

        Int_t cnt_nhit = 0;
        PhySt ppst(rltSt);
        PhyJb::SMtxDGG&& ppjb = SMtxId();
        for (auto&& hit : hits_) {
            PhyJb curjb;
            if (!PropMgnt::PropToZ(hit.cz(), ppst, nullptr, &curjb)) break;
            ppjb = curjb.gg() * ppjb;

            Double_t mex = (hit.sx() ? hit.ex(hit.cx() - ppst.cx()) : MGMath::ZERO);
            Double_t mey = (hit.sy() ? hit.ey(hit.cy() - ppst.cy()) : MGMath::ZERO);

            SMtxSymD<2> cvM;
            cvM(0, 0) = (hit.sx() ? (MGMath::ONE / mex / mex) : MGMath::ZERO);
            cvM(1, 1) = (hit.sy() ? (MGMath::ONE / mey / mey) : MGMath::ZERO);
            
            SVecD<2> rsM;
            rsM(0) = (hit.sx() ? cvM(0, 0) * (hit.cx() - ppst.cx()) : MGMath::ZERO);
            rsM(1) = (hit.sy() ? cvM(1, 1) * (hit.cy() - ppst.cy()) : MGMath::ZERO);
            
            PhyJb::SMtxDXYG&& subJbF = PhyJb::SubXYG(ppjb);
            grdG += LA::Transpose(subJbF) * rsM;
            cvGG += LA::SimilarityT(subJbF, cvM);

            if (hit.sx()) { chix += rsM(0) * (hit.cx() - ppst.cx()); }
            if (hit.sy()) { chiy += rsM(1) * (hit.cy() - ppst.cy()); }

            cnt_nhit++;
            if (!hit.sx()) hit.set_dummy_x(ppst.cx());
        }
        if (cnt_nhit != hits_.size()) break;
        Double_t chi  = (chix + chiy);
        Double_t nchi = ((chi) / static_cast<Double_t>(ndof_));

        Bool_t isSucc   = false;
        Bool_t isUpdate = false;
        if (curIter != 0) {
            Double_t lmRho = (nchi_ - nchi) / curLmRhoDen;
            Bool_t   isLmt = (MGNumc::Compare(lambda, LMTU_LAMBDA) >= 0);
            Double_t convg = std::sqrt(MGMath::ONE + lambda);
            isSucc = (MGNumc::Compare(std::fabs((nchi_ - nchi) / (nchi_ + nchi + CONVG_TOLERANCE)) * convg, CONVG_TOLERANCE) <= 0);

            if (MGNumc::Compare(lmRho, CONVG_EPSILON) < 0) {
                lambda = std::min(lambda*LAMBDA_UP_FAC, LMTU_LAMBDA); 
                grdG   = curGrdG;
                cvGG   = curCvGG;
                rltSt  = part_;
                if      (isSucc) updIter++;
                else if (isLmt)  break;
            }
            else {
                lambda   = std::max(lambda/LAMBDA_DN_FAC, LMTL_LAMBDA);
                chix_    = chix;
                chiy_    = chiy;
                nchi_    = nchi;
                part_    = rltSt;
                isUpdate = true;
                updIter++;
            }
        }
        else { nchi_ = nchi; }

        SMtxSymD<5> lmCvGG(cvGG);
        SVecD<5>&&  diagCvGG = (lambda * cvGG.Diagonal());
        lmCvGG.SetDiagonal(SVecD<5>(lmCvGG.Diagonal() + diagCvGG));

        if (!lmCvGG.Invert()) break;
        SVecD<5>&& rslG = (lmCvGG * grdG);
       
        curLmRhoDen = MGMath::ZERO;
        for (Int_t p = 0; p < 5; ++p)
            curLmRhoDen += (rslG(p) * (diagCvGG(p)*rslG(p) + grdG(p)));
        
        if (curIter == 0 || isUpdate) {
            curGrdG = grdG;
            curCvGG = cvGG;
        }

        rltSt.set_state_with_uxy(
            rltSt.cx() + rslG(0),
            rltSt.cy() + rslG(1),
            rltSt.cz(),
            rltSt.ux() + rslG(2),
            rltSt.uy() + rslG(3),
            ((ortt_ == Orientation::kDownward) ? -1 : 1)
        );
        rltSt.set_eta(rltSt.eta() + rslG(4));
        
        preSucc = curSucc;
        curSucc = (isSucc && updIter >= LMTL_ITER);
        succ    = (preSucc && curSucc);
        
        if (!succ) curIter++;
    }
    //if (!succ) std::cout << Form("FAIL. IT %2d %2d (RIG %14.8f CHI %14.8f) LAMBDA %14.8f\n", curIter, updIter, part_.rig(), nchi_, lambda);
    //else       std::cout << Form("SUCC. IT %2d %2d (RIG %14.8f CHI %14.8f) LAMBDA %14.8f\n", curIter, updIter, part_.rig(), nchi_, lambda);
    
    return succ;
}


#ifdef __CeresSolver__
PhyTrFit::PhyTrFit(TrFitPar& fitPar) : TrFitPar(fitPar) {
    if (!checkHit()) { clear(); return; }
    clear();
    
    ndfx_ = nhtx_ - 2;
    ndfy_ = nhty_ - 3;
    ndof_ = ndfx_ + ndfy_;
    succ_ = physicalFit();
    if (!succ_) clear();
}


void PhyTrFit::clear() {
    succ_ = false;
    part_.reset(type_);
    part_.arg().reset(sw_mscat_, sw_eloss_);
    args_.clear();
    stts_.clear();

    map_hits_.clear();
    map_stts_.clear();

    ndfx_ = 0;
    ndfy_ = 0;
    chix_ = 0;
    chiy_ = 0;
    chit_ = 0;
    chir_ = 0;

    ndof_ = 0;
    nchi_ = 0;
}


Bool_t PhyTrFit::simpleFit() {
    SimpleTrFit simple(dynamic_cast<TrFitPar&>(*this));
    if (!simple.status()) return false;
    part_ = simple.part();
    return true;
}


Bool_t PhyTrFit::physicalFit() {
    if (!simpleFit()) return false;
    part_.arg().reset(sw_mscat_, sw_eloss_);
    std::vector<double> interaction_parameters((numOfHit()-1)*PhyJb::DIM_L, 0.);
    std::vector<double> parameters({ part_.cx(), part_.cy(), part_.ux(), part_.uy(), part_.eta() });
    parameters.insert(parameters.end(), interaction_parameters.begin(), interaction_parameters.end());

    ceres::CostFunction* cost_function = new VirtualPhyTrFit(dynamic_cast<TrFitPar&>(*this), part_);

    ceres::Problem problem;
    problem.AddResidualBlock(cost_function, NULL, parameters.data());
    problem.SetParameterLowerBound(parameters.data(), 2, -1.0);
    problem.SetParameterUpperBound(parameters.data(), 2,  1.0);
    problem.SetParameterLowerBound(parameters.data(), 3, -1.0);
    problem.SetParameterUpperBound(parameters.data(), 3,  1.0);
    
    ceres::Solver::Options options;
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    if (!summary.IsSolutionUsable()) { std::cerr << summary.FullReport() << std::endl; return false; }

    interaction_parameters = std::vector<double>(parameters.begin()+PhyJb::DIM_G, parameters.end());
    part_.set_state_with_uxy(parameters.at(0), parameters.at(1), part_.cz(), parameters.at(2), parameters.at(3), MGNumc::Compare(part_.uz()));
    part_.set_eta(parameters.at(4));

    stts_ = std::vector<PhySt>(numOfHit());
    args_ = std::vector<PhyArg>(numOfHit()-1, PhyArg(sw_mscat_, sw_eloss_));
    for (Int_t it = 0; it < numOfHit()-1; ++it) {
        args_.at(it).set_mscat(
                interaction_parameters.at(it*PhyJb::DIM_L+0), 
                interaction_parameters.at(it*PhyJb::DIM_L+1), 
                interaction_parameters.at(it*PhyJb::DIM_L+2), 
                interaction_parameters.at(it*PhyJb::DIM_L+3));
    }
    
    if (!evolve()) return false;
    return true;
}


Bool_t PhyTrFit::evolve() {
    Double_t chix = 0;
    Double_t chiy = 0;
    Double_t chit = 0;
    Double_t chir = 0;

    PhySt ppst(part_);
    Int_t cnt_nhit = 0;
    for (auto&& hit : hits_) {
        if (cnt_nhit != 0) ppst.arg() = args_.at(cnt_nhit-1);
        if (!PropMgnt::PropToZ(hit.cz(), ppst)) break;
        
        if (cnt_nhit != 0) args_.at(cnt_nhit-1).setvar_mat(ppst.arg().mat(), ppst.arg().nrl(), ppst.arg().ela());
        PhyArg curArg = ppst.arg();
        
        ppst.symbk();
        stts_.at(cnt_nhit) = ppst;

        SVecD<2> rsm((hit.sx() ? (hit.cx() - ppst.cx()) : MGMath::ZERO), (hit.sy() ? (hit.cy() - ppst.cy()) : MGMath::ZERO));
        SVecD<2> rse((hit.sx() ?         hit.ex(rsm(0)) : MGMath::ZERO), (hit.sy() ?         hit.ey(rsm(1)) : MGMath::ZERO));
        SVecD<2> rs2((hit.sx() ?      (rsm(0) / rse(0)) : MGMath::ZERO), (hit.sy() ?      (rsm(1) / rse(1)) : MGMath::ZERO));

        if (hit.sx()) chix += (rs2(0) * rs2(0));
        if (hit.sy()) chiy += (rs2(1) * rs2(1));

        if (cnt_nhit != 0) {
            SVecD<PhyJb::DIM_L> intm(-curArg.tauu(), -curArg.rhou(), -curArg.taul(), -curArg.rhol());
            SVecD<PhyJb::DIM_L> inte(curArg.etauu(), curArg.erhou(), curArg.etaul(), curArg.erhol());
            SVecD<PhyJb::DIM_L> int2(intm(0)/inte(0), intm(1)/inte(1), intm(2)/inte(2), intm(3)/inte(3));
            
            chit += (int2(0) * int2(0) + int2(2) * int2(2));
            chir += (int2(1) * int2(1) + int2(3) * int2(3));
        }

        cnt_nhit++;
    }
    if (cnt_nhit != hits_.size()) return false;

    stts_.back().arg().reset();
    for (Int_t it = 0; it < stts_.size()-1; ++it) {
        stts_.at(it).arg() = args_.at(it);
        stts_.at(it).arg().zero();
    }
    
    for (Int_t it = 0; it < hits_.size(); ++it) {
        map_hits_[hits_.at(it).lay()] = &hits_.at(it);
        map_stts_[hits_.at(it).lay()] = &stts_.at(it);
    }

    Double_t chi  = (chix + chiy + chit + chir);
    Double_t nchi = ((chi) / static_cast<Double_t>(ndof_));

    chix_ = chix;
    chiy_ = chiy;
    chit_ = chit;
    chir_ = chir;
    nchi_ = nchi;
    succ_ = true;

    return true;
}


bool VirtualPhyTrFit::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const {
    std::fill_n(residuals, numOfRes_, 0.);
    Bool_t hasJacb = (jacobians != nullptr && jacobians[0] != nullptr);
    if (hasJacb) std::fill_n(jacobians[0], numOfRes_*numOfPar_, 0.);

    PhySt ppst(part_);
    ppst.set_state_with_uxy(parameters[0][0], parameters[0][1], part_.cz(), parameters[0][2], parameters[0][3], MGNumc::Compare(part_.uz()));
    ppst.set_eta(parameters[0][4]);
    ppst.arg().reset();

    Int_t cnt_nhit = 0;
    PhyJb::SMtxDGG&& jbGG = SMtxId();
    std::vector<PhyJb::SMtxDGL> jbGL(hits_.size()-1);
    Eigen::VectorXd rs = Eigen::VectorXd::Zero(numOfRes_);
    Eigen::MatrixXd jb = Eigen::MatrixXd::Zero(numOfRes_, numOfPar_);
    for (auto&& hit : hits_) {
        PhyJb curjb;
        if (!PropMgnt::PropToZ(hit.cz(), ppst, nullptr, ((hasJacb)?&curjb:nullptr)))  break;
        if (cnt_nhit != 0)
            ppst.arg().set_mscat(
                parameters[0][PhyJb::DIM_G+(cnt_nhit-1)*PhyJb::DIM_L+0], 
                parameters[0][PhyJb::DIM_G+(cnt_nhit-1)*PhyJb::DIM_L+1], 
                parameters[0][PhyJb::DIM_G+(cnt_nhit-1)*PhyJb::DIM_L+2], 
                parameters[0][PhyJb::DIM_G+(cnt_nhit-1)*PhyJb::DIM_L+3]);
        PhyArg curArg = ppst.arg();
        ppst.symbk();

        SVecD<2> rsm((hit.sx() ? (hit.cx() - ppst.cx()) : MGMath::ZERO), (hit.sy() ? (hit.cy() - ppst.cy()) : MGMath::ZERO));
        SVecD<2> rse((hit.sx() ?         hit.ex(rsm(0)) : MGMath::ZERO), (hit.sy() ?         hit.ey(rsm(1)) : MGMath::ZERO));
       
        if (hit.sx()) rs(hit.seqIDx()) += rsm(0) / rse(0);
        if (hit.sy()) rs(hit.seqIDy()) += rsm(1) / rse(1);
        
        if (hasJacb) {
            jbGG = curjb.gg() * jbGG;
            for (Int_t it = 0; it < PhyJb::DIM_G; ++it) {
                if (hit.sx()) jb(hit.seqIDx(), it) += -jbGG(0, it) / rse(0);
                if (hit.sy()) jb(hit.seqIDy(), it) += -jbGG(1, it) / rse(1);
            }
        }    

        SVecD<PhyJb::DIM_L> intm(-curArg.tauu(), -curArg.rhou(), -curArg.taul(), -curArg.rhol());
        SVecD<PhyJb::DIM_L> inte(curArg.etauu(), curArg.erhou(), curArg.etaul(), curArg.erhol());
        
        if (cnt_nhit != 0) {
            rs(numOfSeq()+(cnt_nhit-1)*PhyJb::DIM_L+0) += intm(0) / inte(0);
            rs(numOfSeq()+(cnt_nhit-1)*PhyJb::DIM_L+1) += intm(1) / inte(1);
            rs(numOfSeq()+(cnt_nhit-1)*PhyJb::DIM_L+2) += intm(2) / inte(2);
            rs(numOfSeq()+(cnt_nhit-1)*PhyJb::DIM_L+3) += intm(3) / inte(3);
        }
        
        if (cnt_nhit != 0 && hasJacb) {
            jbGL.at(cnt_nhit-1) = curjb.gl();
            for (Int_t it = 0; it < cnt_nhit-1; ++it)
                jbGL.at(it) = PhyJb::Multiply(curjb.gg(), jbGL.at(it));

            for (Int_t it = 0; it < cnt_nhit; ++it) {
                for (Int_t jl = 0; jl < PhyJb::DIM_L; ++jl) {
                    if (hit.sx()) jb(hit.seqIDx(), PhyJb::DIM_G+it*PhyJb::DIM_L+jl) += -jbGL.at(it)(0, jl) / rse(0);
                    if (hit.sy()) jb(hit.seqIDy(), PhyJb::DIM_G+it*PhyJb::DIM_L+jl) += -jbGL.at(it)(1, jl) / rse(1);
                }
            }
            
            jb(numOfSeq()+(cnt_nhit-1)*PhyJb::DIM_L+0, PhyJb::DIM_G+(cnt_nhit-1)*PhyJb::DIM_L+0) += -1.0 / inte(0);
            jb(numOfSeq()+(cnt_nhit-1)*PhyJb::DIM_L+1, PhyJb::DIM_G+(cnt_nhit-1)*PhyJb::DIM_L+1) += -1.0 / inte(1);
            jb(numOfSeq()+(cnt_nhit-1)*PhyJb::DIM_L+2, PhyJb::DIM_G+(cnt_nhit-1)*PhyJb::DIM_L+2) += -1.0 / inte(2);
            jb(numOfSeq()+(cnt_nhit-1)*PhyJb::DIM_L+3, PhyJb::DIM_G+(cnt_nhit-1)*PhyJb::DIM_L+3) += -1.0 / inte(3);
        }

        cnt_nhit++;
    }
    if (cnt_nhit != hits_.size()) return false;
    
    for (Int_t it = 0; it < numOfRes_; ++it)
        residuals[it] = rs(it);
    
    if (hasJacb)
        for (Int_t it = 0; it < numOfRes_; ++it)
            for (Int_t jt = 0; jt < numOfPar_; ++jt)
                jacobians[0][it * numOfPar_ + jt] = jb(it, jt); 

    return true;
}
#endif


} // namespace TrackSys


#endif // __TRACKLibs_PhyFit_C__
