#ifndef __TRACKLibs_PhyFit_C__
#define __TRACKLibs_PhyFit_C__


namespace TrackSys {
        
    
void TrFitPar::print() const {
    std::string outstr;
    outstr += "TrFitPar::Print()\n";
    for (auto&& hit : hits_TRK_) outstr += STR("TRK L%d COO %14.8f %14.8f %14.8f CERR %14.8f %14.8f\n", hit.lay(), hit.cx(), hit.cy(), hit.cz(), hit.cex(), hit.cey());
    for (auto&& hit : hits_TOF_) outstr += STR("TOF L%d COO %14.8f %14.8f %14.8f CERR %14.8f %14.8f\n", hit.lay(), hit.cx(), hit.cy(), hit.cz(), hit.cex(), hit.cey());
    outstr += "\n";
    COUT("%s", outstr.c_str());
}
    
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
    nseq_ = 0;
    nhtx_ = 0;
    nhty_ = 0;

    is_check_ = false;
    rlt_check_ = -1;
}

Bool_t TrFitPar::rebuildHit() {
    hits_.clear();
    if (ortt_ == Orientation::kDownward) Hit<HitStTRK>::Sort(hits_TRK_, VirtualHitSt::Orientation::kDownward);
    else                                 Hit<HitStTRK>::Sort(hits_TRK_, VirtualHitSt::Orientation::kUpward);
    if (ortt_ == Orientation::kDownward) Hit<HitStTOF>::Sort(hits_TOF_, VirtualHitSt::Orientation::kDownward);
    else                                 Hit<HitStTOF>::Sort(hits_TOF_, VirtualHitSt::Orientation::kUpward);
    
    if (hits_TRK_.size() == 0 && hits_TOF_.size() == 0) return false;

    for (auto&& hit : hits_TRK_) hits_.push_back(&hit); 
    for (auto&& hit : hits_TOF_) hits_.push_back(&hit); 
    if (ortt_ == Orientation::kDownward) VirtualHitSt::Sort(hits_, VirtualHitSt::Orientation::kDownward);
    else                                 VirtualHitSt::Sort(hits_, VirtualHitSt::Orientation::kUpward);
    if (hits_.size() == 0) return false;

    return true;
}

Short_t TrFitPar::checkHit() {
    if (is_check_) return rlt_check_;
    nseq_ = 0; nhtx_ = 0; nhty_ = 0;
    
    rlt_check_ = -1;
    if (!rebuildHit()) return rlt_check_;
    
    Short_t nx = 0, ny = 0;
    for (auto&& hit : hits_) {
        if (hit->csx()) nx++;
        if (hit->csy()) ny++;
    }
    rlt_check_ = 0;
    if (nx < LMTL_NHIT_X && ny < LMTL_NHIT_Y) return rlt_check_;
   
    Short_t nseq = 0;
    for (auto&& hit : hits_) {
        hit->set_type(type_);
        nseq += hit->set_seqID(nseq);
    }

    nseq_ = nseq;
    nhtx_ = nx;
    nhty_ = ny;

    is_check_ = true;
    rlt_check_ = 1;
    return rlt_check_;
}
        

SimpleTrFit::SimpleTrFit(TrFitPar& fitPar) : TrFitPar(fitPar) {
    SimpleTrFit::clear();
    if (recheckHit() <= 0) return;
   
    ndfx_ = nhtx_ - 2;
    ndfy_ = nhty_ - 3;
    ndof_ = ndfx_ + ndfy_;
    
    succ_ = (analyticalFit() ? simpleFit() : false);
    if (!succ_) { SimpleTrFit::clear(); COUT("SimpleTrFit::FAIL.\n"); }
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
            if (!hit->csx()) continue;
            Double_t ex  = hit->cex();
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
        Double_t Lambda = PROP_FACT * part_.info().chrg_to_mass(); 
        
        std::vector<UInt_t> mapID;
        for (UInt_t ih = 0; ih < hits_.size(); ++ih)
            if (hits_.at(ih)->csy()) mapID.push_back(ih);

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
            Double_t ey  = hit->cey();
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
    
    return true;
}


Bool_t SimpleTrFit::simpleFit() {
    Bool_t succ    = false;
    Bool_t preSucc = false;
    Bool_t curSucc = false;

    Double_t    curLmRhoDen = Numc::ONE<>;
    Double_t    lambda = LAMBDA0;
    PhySt       rltSt(part_);
    SVecD<5>    curGrdG;
    SMtxSymD<5> curCvGG;

    UInt_t updIter = 0;
    UInt_t curIter = 0;
    while (curIter <= LMTU_ITER && !succ) {
        Double_t chix = Numc::ZERO<>;
        Double_t chiy = Numc::ZERO<>;
        SVecD<5>    grdG;
        SMtxSymD<5> cvGG;

        UInt_t cnt_nhit = 0;
        PhySt ppst(rltSt);
        PhyJb::SMtxDGG&& ppjb = SMtxId();
        for (auto&& hit : hits_) {
            PhyJb curjb;
            if (!PropMgnt::PropToZ(hit->cz(), ppst, nullptr, &curjb)) break;
            ppjb = curjb.gg() * ppjb;
            hit->cal(ppst);
            
            SVecD<2> rsM;
            rsM(0) = (hit->seqIDcx()>=0 ? hit->cnrmx() * hit->cdivx() : Numc::ZERO<>);
            rsM(1) = (hit->seqIDcy()>=0 ? hit->cnrmy() * hit->cdivy() : Numc::ZERO<>);
            
            SMtxSymD<2> cvM;
            cvM(0, 0) = (hit->seqIDcx()>=0 ? (hit->cdivx() * hit->cdivx()) : Numc::ZERO<>);
            cvM(1, 1) = (hit->seqIDcy()>=0 ? (hit->cdivy() * hit->cdivy()) : Numc::ZERO<>);
            
            PhyJb::SMtxDXYG&& subJbF = PhyJb::SubXYG(ppjb);
            grdG += LA::Transpose(subJbF) * rsM;
            cvGG += LA::SimilarityT(subJbF, cvM);
            
            if (hit->seqIDcx()>=0) chix += hit->cnrmx() * hit->cnrmx(); 
            if (hit->seqIDcy()>=0) chiy += hit->cnrmy() * hit->cnrmy();
            
            cnt_nhit++;
        }
        if (cnt_nhit != hits_.size()) break;
        Double_t chi  = (chix + chiy);
        Double_t nchi = ((chi) / static_cast<Double_t>(ndof_));

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

        Bool_t isNomag = (Numc::EqualToZero(lmCvGG(4, 4)) && Numc::EqualToZero(grdG(4))); // Fast Check
        if (isNomag) {
            SMtxSymD<4>&& lmCvGG_nomag = lmCvGG.Sub<SMtxSymD<4>>(0, 0);
            if (!lmCvGG_nomag.Invert()) break;
            lmCvGG = std::move(SMtxSymD<5>());
            for (UInt_t ielem = 0; ielem < 4; ++ielem)
                for (UInt_t jelem = ielem; jelem < 4; ++jelem)
                    lmCvGG(ielem, jelem) = lmCvGG_nomag(ielem, jelem);
        }
        else {
            if (!lmCvGG.Invert()) break;
        }
        SVecD<5>&& rslG = (lmCvGG * grdG);
       
        curLmRhoDen = Numc::ZERO<>;
        for (UInt_t p = 0; p < 5; ++p)
            curLmRhoDen += (rslG(p) * (diagCvGG(p)*rslG(p) + grdG(p)));
        
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
            ((ortt_ == Orientation::kDownward) ? -1 : 1)
        );
        rltSt.set_eta(rltSt.eta() - rslG(4));
        
        preSucc = curSucc;
        curSucc = (isSucc && updIter >= LMTL_ITER);
        succ    = (preSucc && curSucc);
        
        if (!succ) curIter++;
    }
    //if (!succ) COUT("FAIL. IT %2d %2d (RIG %14.8f CHI %14.8f) LAMBDA %14.8f\n", curIter, updIter, part_.rig(), nchi_, lambda);
    //else       COUT("SUCC. IT %2d %2d (RIG %14.8f CHI %14.8f) LAMBDA %14.8f\n", curIter, updIter, part_.rig(), nchi_, lambda);
    
    return succ;
}


PhyTrFit::PhyTrFit(TrFitPar& fitPar) : TrFitPar(fitPar) {
    if (recheckHit() <= 0) { clear(); return; }
    clear();
    
    ndfx_ = nhtx_ - 2;
    ndfy_ = nhty_ - 3;
    ndof_ = ndfx_ + ndfy_;
    succ_ = physicalFit();
    if (!succ_) { PhyTrFit::clear(); COUT("PhyTrFit::FAIL.\n"); }
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
    if (!summary.IsSolutionUsable()) return false;

    interaction_parameters = std::vector<double>(parameters.begin()+PhyJb::DIM_G, parameters.end());
    part_.set_state_with_uxy(parameters.at(0), parameters.at(1), part_.cz(), parameters.at(2), parameters.at(3), Numc::Compare(part_.uz()));
    part_.set_eta(parameters.at(4));

    stts_ = std::vector<PhySt>(numOfHit());
    args_ = std::vector<PhyArg>(numOfHit()-1, PhyArg(sw_mscat_, sw_eloss_));
    for (Int_t it = 0; it < numOfHit()-1; ++it) {
        args_.at(it).set_mscat(
                interaction_parameters.at(it*PhyJb::DIM_L+0), 
                interaction_parameters.at(it*PhyJb::DIM_L+1), 
                interaction_parameters.at(it*PhyJb::DIM_L+2), 
                interaction_parameters.at(it*PhyJb::DIM_L+3));
        //args_.at(it).set_eloss(
        //        interaction_parameters.at(it*PhyJb::DIM_L+4));
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
    UInt_t cnt_nhit = 0;
    for (auto&& hit : hits_) {
        if (cnt_nhit != 0) ppst.arg() = args_.at(cnt_nhit-1);
        if (!PropMgnt::PropToZ(hit->cz(), ppst)) break;
        
        if (cnt_nhit != 0) args_.at(cnt_nhit-1).setvar_mat(ppst.arg().mat(), ppst.arg().nrl(), ppst.arg().ela());
        PhyArg curArg = ppst.arg();
        ppst.symbk();

        stts_.at(cnt_nhit) = ppst;
        hit->cal(ppst);

        if (hit->seqIDcx()>=0) chix += hit->cnrmx() * hit->cnrmx(); 
        if (hit->seqIDcy()>=0) chiy += hit->cnrmy() * hit->cnrmy();
        
        if (cnt_nhit != 0) {
            SVecD<PhyJb::DIM_L> intm(-curArg.tauu(), -curArg.rhou(), -curArg.taul(), -curArg.rhol());
            SVecD<PhyJb::DIM_L> inte(curArg.etauu(), curArg.erhou(), curArg.etaul(), curArg.erhol());
            SVecD<PhyJb::DIM_L> int2(intm(0)/inte(0), intm(1)/inte(1), intm(2)/inte(2), intm(3)/inte(3));
            //SVecD<PhyJb::DIM_L> intm(-curArg.tauu(), -curArg.rhou(), -curArg.taul(), -curArg.rhol(), -curArg.elion());
            //SVecD<PhyJb::DIM_L> inte(curArg.etauu(), curArg.erhou(), curArg.etaul(), curArg.erhol(), curArg.eelion());
            //SVecD<PhyJb::DIM_L> int2(intm(0)/inte(0), intm(1)/inte(1), intm(2)/inte(2), intm(3)/inte(3), intm(4)/inte(4));
            
            chit += (int2(0) * int2(0) + int2(2) * int2(2));
            chir += (int2(1) * int2(1) + int2(3) * int2(3));
        }

        cnt_nhit++;
    }
    if (cnt_nhit != hits_.size()) return false;

    stts_.back().arg().reset();
    for (UInt_t it = 0; it < stts_.size()-1; ++it) {
        stts_.at(it).arg() = args_.at(it);
        stts_.at(it).arg().zero();
    }
    
    for (UInt_t it = 0; it < hits_.size(); ++it) {
        if (!Hit<HitStTRK>::IsSame(hits_.at(it))) continue;
        map_hits_[hits_.at(it)->lay()] =  hits_.at(it);
        map_stts_[hits_.at(it)->lay()] = &stts_.at(it);
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

    Bool_t resetTOF = true;
    HitStTOF::SetOffsetTime(Numc::ZERO<>);
    HitStTOF::SetOffsetPath(Numc::ZERO<>);

    PhySt ppst(part_);
    ppst.set_state_with_uxy(parameters[0][0], parameters[0][1], part_.cz(), parameters[0][2], parameters[0][3], Numc::Compare(part_.uz()));
    ppst.set_eta(parameters[0][4]);
    ppst.arg().reset();

    UInt_t cnt_nhit = 0;
    PhyJb::SMtxDGG&& jbGG = SMtxId();
    std::vector<PhyJb::SMtxDGL> jbGL(hits_.size()-1);
    Eigen::VectorXd rs = Eigen::VectorXd::Zero(numOfRes_);
    Eigen::MatrixXd jb = Eigen::MatrixXd::Zero(numOfRes_, numOfPar_);
    for (auto&& hit : hits_) {
        PhyJb curjb;
        if (!PropMgnt::PropToZ(hit->cz(), ppst, nullptr, ((hasJacb)?&curjb:nullptr)))  break;
        if (cnt_nhit != 0) {
            ppst.arg().set_mscat(
                parameters[0][PhyJb::DIM_G+(cnt_nhit-1)*PhyJb::DIM_L+0], 
                parameters[0][PhyJb::DIM_G+(cnt_nhit-1)*PhyJb::DIM_L+1], 
                parameters[0][PhyJb::DIM_G+(cnt_nhit-1)*PhyJb::DIM_L+2], 
                parameters[0][PhyJb::DIM_G+(cnt_nhit-1)*PhyJb::DIM_L+3]);
            //ppst.arg().set_eloss(
            //    parameters[0][PhyJb::DIM_G+(cnt_nhit-1)*PhyJb::DIM_L+4]);
        }
        PhyArg curArg = ppst.arg();
        ppst.symbk();

        if (resetTOF && Hit<HitStTOF>::IsSame(hit)) { // set reference
            HitStTOF::SetOffsetPath(ppst.path());
            HitStTOF::SetOffsetTime(ppst.time()-Hit<HitStTOF>::Cast(hit)->t());
            resetTOF = false;
        }
        hit->cal(ppst);
       
        if (hit->seqIDcx()>=0) rs(hit->seqIDcx()) += hit->cnrmx();
        if (hit->seqIDcy()>=0) rs(hit->seqIDcy()) += hit->cnrmy();
        if (hasJacb) {
            jbGG = curjb.gg() * jbGG;
            for (UInt_t it = 0; it < PhyJb::DIM_G; ++it) {
                if (hit->seqIDcx()>=0) jb(hit->seqIDcx(), it) += hit->cdivx() * jbGG(0, it);
                if (hit->seqIDcy()>=0) jb(hit->seqIDcy(), it) += hit->cdivy() * jbGG(1, it);
            }
        }    

        // TRK
        HitStTRK* hitTRK = Hit<HitStTRK>::Cast(hit);
        if (hitTRK != nullptr) {
            if (hitTRK->seqIDax()>=0) rs(hitTRK->seqIDax()) += hitTRK->anrmx();
            if (hitTRK->seqIDay()>=0) rs(hitTRK->seqIDay()) += hitTRK->anrmy();
            if (hasJacb) {
            for (UInt_t it = 0; it < PhyJb::DIM_G; ++it) {
                if (hitTRK->seqIDax()>=0) jb(hitTRK->seqIDax(), it) += hitTRK->adivx() * jbGG(4, it);
                if (hitTRK->seqIDay()>=0) jb(hitTRK->seqIDay(), it) += hitTRK->adivy() * jbGG(4, it);
            }}
        }
        
        // TOF
        HitStTOF* hitTOF = Hit<HitStTOF>::Cast(hit);
        if (hitTOF != nullptr) {
            if (hitTOF->seqIDq()>=0) rs(hitTOF->seqIDq()) += hitTOF->qnrm();
            if (hasJacb) {
            for (UInt_t it = 0; it < PhyJb::DIM_G; ++it) {
                if (hitTOF->seqIDq()>=0) jb(hitTOF->seqIDq(), it) += hitTOF->qdiv() * jbGG(4, it);
                if (hitTOF->seqIDt()>=0) jb(hitTOF->seqIDt(), it) += hitTOF->tdiv() * jbGG(4, it);
            }}
        }

        // Interaction
        SVecD<PhyJb::DIM_L> intm(-curArg.tauu(), -curArg.rhou(), -curArg.taul(), -curArg.rhol());
        SVecD<PhyJb::DIM_L> inte(curArg.etauu(), curArg.erhou(), curArg.etaul(), curArg.erhol());
        //SVecD<PhyJb::DIM_L> intm(-curArg.tauu(), -curArg.rhou(), -curArg.taul(), -curArg.rhol(), -curArg.elion());
        //SVecD<PhyJb::DIM_L> inte(curArg.etauu(), curArg.erhou(), curArg.etaul(), curArg.erhol(), curArg.eelion());
        
        if (cnt_nhit != 0) {
            for (UInt_t it = 0; it < PhyJb::DIM_L; ++it)
                rs(numOfSeq()+(cnt_nhit-1)*PhyJb::DIM_L+it) += intm(it) / inte(it);
        }
        
        if (cnt_nhit != 0 && hasJacb) {
            jbGL.at(cnt_nhit-1) = curjb.gl();
            for (UInt_t it = 0; it < cnt_nhit-1; ++it)
                jbGL.at(it) = curjb.gg() * jbGL.at(it);

            for (UInt_t it = 0; it < cnt_nhit; ++it) {
                for (UInt_t jl = 0; jl < PhyJb::DIM_L; ++jl) {
                    if (hit->seqIDcx()>=0) jb(hit->seqIDcx(), PhyJb::DIM_G+it*PhyJb::DIM_L+jl) += hit->cdivx() * jbGL.at(it)(0, jl);
                    if (hit->seqIDcy()>=0) jb(hit->seqIDcy(), PhyJb::DIM_G+it*PhyJb::DIM_L+jl) += hit->cdivy() * jbGL.at(it)(1, jl);
                }
            }
            
            for (UInt_t it = 0; it < PhyJb::DIM_L; ++it)
                jb(numOfSeq()+(cnt_nhit-1)*PhyJb::DIM_L+it, PhyJb::DIM_G+(cnt_nhit-1)*PhyJb::DIM_L+it) += -1.0 / inte(it);
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


} // namespace TrackSys


#endif // __TRACKLibs_PhyFit_C__
