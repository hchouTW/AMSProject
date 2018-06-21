#ifndef __TRACKLibs_PhyFit_C__
#define __TRACKLibs_PhyFit_C__


namespace TrackSys {
        
    
TrFitPar& TrFitPar::operator=(const TrFitPar& rhs) {
    if (this != &rhs) {
        noise_ctler_  = rhs.noise_ctler_;
        sw_mscat_     = rhs.sw_mscat_;
        sw_eloss_     = rhs.sw_eloss_;
        info_         = rhs.info_;
        ortt_         = rhs.ortt_;
        hits_TRK_     = rhs.hits_TRK_;
        hits_TOF_     = rhs.hits_TOF_;
        hits_RICH_    = rhs.hits_RICH_;
        hits_TRD_     = rhs.hits_TRD_;
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
    
TrFitPar::TrFitPar(const PartInfo& info, const Orientation& ortt, const Bool_t& sw_mscat, const Bool_t& sw_eloss, const VirtualHitSt::NoiseController& noise_ctler) {
    clear();

    noise_ctler_ = noise_ctler;
    sw_mscat_ = sw_mscat;
    sw_eloss_ = sw_eloss;
    info_ = info;
    ortt_ = ortt;
}

void TrFitPar::zero() {
    hits_.clear();

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
    noise_ctler_ = VirtualHitSt::NoiseController::OFF;

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
    
    Short_t nseq =  0;
    Short_t nseg = -1;
    for (auto&& hit : hits_) {
        hit->set_type(info_);
        nseq += hit->set_seqID(nseq);
        
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
    
    nseq_ = nseq;
    nseg_ = ((sw_mscat_ && nseg > 0) ? nseg : 0);
    nmes_ = (nmes_cx_ + nmes_cy_ + nmes_ib_);

    return true;
}

Bool_t TrFitPar::check_hits() {
    if (is_check_) return is_check_;
    sort_hits();
   
    Bool_t passed = (nmes_cx_ > LMTN_CX && nmes_cy_ > LMTN_CY);
    if (passed) is_check_ = true;
   
    return is_check_;
}

SimpleTrFit& SimpleTrFit::operator=(const SimpleTrFit& rhs) {
    if (this != &rhs) {
        dynamic_cast<TrFitPar&>(*this) = dynamic_cast<const TrFitPar&>(rhs);
        succ_    = rhs.succ_;
        part_    = rhs.part_;
        ndof_    = rhs.ndof_;
        ndof_cx_ = rhs.ndof_cx_;
        ndof_cy_ = rhs.ndof_cy_;
        ndof_ib_ = rhs.ndof_ib_;
        nchi_    = rhs.nchi_;
        nchi_cx_ = rhs.nchi_cx_;
        nchi_cy_ = rhs.nchi_cy_;
        nchi_ib_ = rhs.nchi_ib_;
    }
    return *this;
}

SimpleTrFit::SimpleTrFit(const TrFitPar& fitPar) : TrFitPar(fitPar) {
    SimpleTrFit::clear();
    if (!check_hits()) return;
    ndof_cx_ = (nmes_cx_ > LMTN_CX) ? (nmes_cx_ - LMTN_CX) : 0;
    ndof_cy_ = (nmes_cy_ > LMTN_CY) ? (nmes_cy_ - LMTN_CY) : 0;
    ndof_ib_ = nmes_ib_ - (nmes_TOFt_ >= LMTN_TOF_T);
    ndof_    = (ndof_cx_ + ndof_cy_ + ndof_ib_);
    if (ndof_cx_            <= Numc::ONE<Short_t>) { SimpleTrFit::clear(); return; }
    if ((ndof_cy_+ndof_ib_) <= Numc::ONE<Short_t>) { SimpleTrFit::clear(); return; }

    succ_ = (analyticalFit() ? simpleFit() : false);
    if (!succ_) { SimpleTrFit::clear(); TrFitPar::clear(); }
}


void SimpleTrFit::clear() {
    noise_ctler_ = VirtualHitSt::NoiseController::OFF;
    sw_mscat_ = false;
    sw_eloss_ = false;
    
    succ_ = false;
    part_.reset(info_);
    part_.arg().reset(sw_mscat_, sw_eloss_);

    ndof_    = 0;
    ndof_cx_ = 0;
    ndof_cy_ = 0;
    ndof_ib_ = 0;
        
    nchi_    = 0;
    nchi_cx_ = 0;
    nchi_cy_ = 0;
    nchi_ib_ = 0;
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
    const Double_t LMT_ETA = 5.0;
    if (part_.eta_abs() > LMT_ETA)
        part_.set_eta(part_.eta_sign() * LMT_ETA);
    
    return true;
}


Bool_t SimpleTrFit::simpleFit() {
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
        
            if (resetTOF && Hit<HitStTOF>::IsSame(hit)) { // set reference
                HitStTOF::SetOffsetPath(ppst.path());
                HitStTOF::SetOffsetTime(ppst.time()-Hit<HitStTOF>::Cast(hit)->orgt());
                resetTOF = false;
            }
            hit->cal(ppst, noise_ctler_);

            SVecD<2>       rsC;
            SMtxD<2, DIMG> jbC;
            if (hit->scx()) rsC(0) += hit->nrmcx();
            if (hit->scy()) rsC(1) += hit->nrmcy();
            for (Short_t it = 0; it < DIMG; ++it) {
                if (hit->scx()) jbC(0, it) += hit->divcx() * ppjb(0, it);
                if (hit->scy()) jbC(1, it) += hit->divcy() * ppjb(1, it);
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

                if (hitTRK->sqx()) rsTRK(0) += hitTRK->nrmqx();
                if (hitTRK->sqy()) rsTRK(1) += hitTRK->nrmqy();
                    
                if (hitTRK->sqx()) jbTRK(0) += hitTRK->divqx_eta() * ppjb(4, 4);
                if (hitTRK->sqy()) jbTRK(1) += hitTRK->divqy_eta() * ppjb(4, 4);
                
                grdG(4)    += (jbTRK(0) * rsTRK(0) + jbTRK(1) * rsTRK(1));
                cvGG(4, 4) += (jbTRK(0) * jbTRK(0) + jbTRK(1) * jbTRK(1));
                if (hitTRK->sqx()) chi_ib += rsTRK(0) * rsTRK(0);
                if (hitTRK->sqy()) chi_ib += rsTRK(1) * rsTRK(1);
            }

            // TOF
            HitStTOF* hitTOF = Hit<HitStTOF>::Cast(hit);
            if (hitTOF != nullptr) {
                SVecD<2> rsTOF;
                SVecD<2> jbTOF;

                if (hitTOF->st()) rsTOF(0) += hitTOF->nrmt();
                if (hitTOF->sq()) rsTOF(1) += hitTOF->nrmq();
                
                if (hitTOF->st()) jbTOF(0) += hitTOF->divt_eta() * ppjb(4, 4);
                if (hitTOF->sq()) jbTOF(1) += hitTOF->divq_eta() * ppjb(4, 4);
                
                grdG(4)    += (jbTOF(0) * rsTOF(0) + jbTOF(1) * rsTOF(1));
                cvGG(4, 4) += (jbTOF(0) * jbTOF(0) + jbTOF(1) * jbTOF(1));
                if (hitTOF->st()) chi_ib += rsTOF(0) * rsTOF(0);
                if (hitTOF->sq()) chi_ib += rsTOF(1) * rsTOF(1);
            }
            
            // RICH
            HitStRICH* hitRICH = Hit<HitStRICH>::Cast(hit);
            if (hitRICH != nullptr) {
                Double_t rsRICH = Numc::ZERO<>;
                Double_t jbRICH = Numc::ZERO<>;

                if (hitRICH->sib()) rsRICH += hitRICH->nrmib();
                if (hitRICH->sib()) jbRICH += hitRICH->divib_eta() * ppjb(4, 4);

                grdG(4)    += (jbRICH * rsRICH);
                cvGG(4, 4) += (jbRICH * jbRICH);

                if (hitRICH->sib()) chi_ib += rsRICH * rsRICH;
            }

            // TRD
            HitStTRD* hitTRD = Hit<HitStTRD>::Cast(hit);
            if (hitTRD != nullptr) {
                Double_t rsTRD = Numc::ZERO<>;
                Double_t jbTRD = Numc::ZERO<>;

                if (hitTRD->sel()) rsTRD += hitTRD->nrmel();
                if (hitTRD->sel()) jbTRD += hitTRD->divel_eta() * ppjb(4, 4);

                grdG(4)    += (jbTRD * rsTRD);
                cvGG(4, 4) += (jbTRD * jbTRD);

                if (hitTRD->sel()) chi_ib += rsTRD * rsTRD;
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
    //if (!succ) CERR("FAIL. IT %2d %2d (RIG %14.8f MASS %14.8f CHI %14.8f) LAMBDA %14.8f\n", curIter, updIter, part_.rig(), part_.mass(), nchi_, lambda);
    //else       CERR("SUCC. IT %2d %2d (RIG %14.8f MASS %14.8f CHI %14.8f) LAMBDA %14.8f\n", curIter, updIter, part_.rig(), part_.mass(), nchi_, lambda);
    
    return succ;
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
    if (MuOpt::kFixed == mu_opt_) succ_ = (simpleFit() ? physicalFit(MuOpt::kFixed) : false);
    else                          succ_ = physicalMassFit();
    
    if (!succ_) { PhyTrFit::clear(); TrFitPar::clear(); }
}


Bool_t PhyTrFit::simpleFit() {
    SimpleTrFit simple(dynamic_cast<TrFitPar&>(*this));
    if (simple.status()) part_ = simple.part();
    part_.arg().reset(sw_mscat_, sw_eloss_);
    return simple.status();
}


Bool_t PhyTrFit::physicalFit(const MuOpt& mu_opt, Double_t fluc_eta, Double_t fluc_igb, Bool_t with_mu_est) {
    const Bool_t is_mu_free = (MuOpt::kFree == mu_opt);
    const Short_t DIMM = (is_mu_free ? 1 : 0);

    Double_t eta = part_.eta();
    const Bool_t is_fluc_eta = (MuOpt::kFree == mu_opt && Numc::Compare(fluc_eta) > 0);
    if (is_fluc_eta) {
        const Int_t niter = 5; Int_t iter = 0;
        do {
            Double_t rndm = Rndm::NormalGaussian();
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
            Double_t rndm = Rndm::NormalGaussian();
            if (std::fabs(rndm) > Numc::TWO<>) { iter++; continue; }
            igb = part_.igmbta() * (Numc::ONE<> + fluc_igb * rndm);
            iter++;
        } while ((iter < niter) && (igb < LMTL_INV_GB || igb > LMTU_INV_GB));
        if (iter >= niter) igb = part_.igmbta();
    }
    
    Short_t parIDigb = -1;
    std::vector<double> parameters({ part_.cx(), part_.cy(), part_.ux(), part_.uy(), eta });
    if (is_mu_free) { parIDigb = parameters.size(); parameters.push_back(igb); }

    if (nseg_ != 0) {
        std::vector<double> interaction_parameters(nseg_*DIML, Numc::ZERO<>);
        if (args_.size() != nseg_) args_ = std::move(std::vector<PhyArg>(nseg_, PhyArg(sw_mscat_, sw_eloss_)));
        else {
            for (Short_t is = 0; is < nseg_; ++is) {
                const PhyArg& arg = args_.at(is);
                interaction_parameters.at(is*DIML+0) = arg.tauu();
                interaction_parameters.at(is*DIML+1) = arg.rhou();
                interaction_parameters.at(is*DIML+2) = arg.taul();
                interaction_parameters.at(is*DIML+3) = arg.rhol();
            }
        }
        parameters.insert(parameters.end(), interaction_parameters.begin(), interaction_parameters.end());
    }
    if (nmes_TOFt_ >= LMTN_TOF_T) { TOFt_sft_ = Numc::ZERO<>; parameters.push_back(Numc::ZERO<>); } // TOF Shift Time
    else TOFt_sft_ = Numc::ZERO<>;

    ceres::CostFunction* cost_function = new VirtualPhyTrFit(dynamic_cast<TrFitPar&>(*this), part_, is_mu_free);

    ceres::Problem problem;
    problem.AddResidualBlock(cost_function, nullptr, parameters.data());
    problem.SetParameterLowerBound(parameters.data(), 2, -1.0);
    problem.SetParameterUpperBound(parameters.data(), 2,  1.0);
    problem.SetParameterLowerBound(parameters.data(), 3, -1.0);
    problem.SetParameterUpperBound(parameters.data(), 3,  1.0);
    if (is_mu_free) problem.SetParameterLowerBound(parameters.data(), parIDigb, LMTL_INV_GB);
    if (is_mu_free) problem.SetParameterUpperBound(parameters.data(), parIDigb, LMTU_INV_GB);
    
    ceres::Solver::Options options;
    //options.line_search_direction_type = ceres::LBFGS;

    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    if (!summary.IsSolutionUsable()) return false;
    
    Double_t partcz = part_.cz();
    Short_t  signuz = Numc::Compare(part_.uz());
    Double_t partmu = (is_mu_free ? std::fabs(parameters.at(parIDigb) / parameters.at(4)) : part_.mu());
    if (is_mu_free) part_.reset(partmu);
    else            part_.arg().clear();
    part_.set_state_with_uxy(parameters.at(0), parameters.at(1), partcz, parameters.at(2), parameters.at(3), signuz);
    part_.set_eta(parameters.at(4));

    for (Short_t is = 0; is < nseg_; ++is) {
        args_.at(is).set_mscat(
            parameters.at(DIMG+DIMM+is*DIML+0), 
            parameters.at(DIMG+DIMM+is*DIML+1), 
            parameters.at(DIMG+DIMM+is*DIML+2), 
            parameters.at(DIMG+DIMM+is*DIML+3));
    }
    
    if (nmes_TOFt_ >= LMTN_TOF_T) { TOFt_sft_ = parameters.at(DIMG+DIMM+nseg_*DIML); } // TOF Shift Time
    else TOFt_sft_ = Numc::ZERO<>;

    //std::cerr << summary.FullReport() << std::endl;
    Bool_t   rw_err_mu = (MuOpt::kFixed == mu_opt && with_mu_est && evolve(MuOpt::kFree) && Numc::Valid(err_.at(6)) && Numc::Compare(err_.at(6)) > 0);
    Double_t err_mu    = (rw_err_mu ? err_.at(6) : Numc::ZERO<>);
    
    Bool_t succ = evolve(mu_opt);
    if (rw_err_mu) err_.at(6) = err_mu;
    return succ;
}


Bool_t PhyTrFit::physicalMassFit() {
    Short_t chrg = std::abs(info_.chrg());
    if (chrg <= Numc::ZERO<Short_t> || chrg >= PartListMassQ.size()) return false;
    if (PartListMassQ.at(chrg).size() == 0) return false;

    class PartElem {
        public :
            PartElem() : succ(false), fluc_eta(0), fluc_igb(0), qlt(0) {}
            PartElem(const PhySt& _part, Double_t _fluc_eta, Double_t _fluc_ibta, Double_t _qlt) : succ(true), part(_part), fluc_eta(_fluc_eta), fluc_igb(_fluc_ibta), qlt(_qlt) {}
            Bool_t              succ;
            PhySt               part;
            Double_t            fluc_eta;
            Double_t            fluc_igb;
            Double_t            qlt;
    };
   
    // List of Particle Mass (Init)
    PartElem condElem;
    for (auto&& mass : PartListMassQ.at(chrg)) {
        args_.clear();
        info_.reset(chrg, mass);
        TOFt_sft_ = Numc::ZERO<>;
        if (!(simpleFit() ? physicalFit(MuOpt::kFixed, Numc::ZERO<>, Numc::ZERO<>, false) : false)) continue;
        
        Bool_t firstTime = (!condElem.succ);
        if (!firstTime && Numc::Compare(quality_.at(1), condElem.qlt) > 0) continue;
        if (!evolve(MuOpt::kFree)) continue;
        
        Double_t fluc_eta  = std::fabs(err_.at(4) / part_.eta());
        Double_t fluc_igb = (err_.at(5) / part_.igmbta());
        condElem = std::move(PartElem(part_, fluc_eta, fluc_igb, quality_.at(1)));
    }
    if (!condElem.succ) return false;

    for (Short_t iter = 1; iter <= LMTU_MU_ITER; ++iter) {
        args_.clear();
        TOFt_sft_ = Numc::ZERO<>;
        info_ = condElem.part.info();
        part_ = condElem.part;

        Double_t mass = part_.mass();
        Double_t igb  = part_.igmbta();
        if (!physicalFit(MuOpt::kFree, 
                         MU_FLUC * condElem.fluc_eta, 
                         MU_FLUC * condElem.fluc_igb)) return false;
        
        Double_t fluc_eta  = std::fabs(err_.at(4) / part_.eta());
        Double_t fluc_igb = (err_.at(5) / part_.igmbta());
        condElem = std::move(PartElem(part_, fluc_eta, fluc_igb, quality_.at(1)));
       
        Double_t convg_m = std::fabs((part_.mass()   - mass) / (part_.mass()   + mass));
        Double_t convg_b = std::fabs((part_.igmbta() - igb)  / (part_.igmbta() + igb)) / fluc_igb;
        if (!Numc::Valid(convg_m)) convg_m = Numc::ONE<>;
        if (!Numc::Valid(convg_b)) convg_b = Numc::ONE<>;

        Bool_t convg = (convg_m < CONVG_FLUC || convg_b < CONVG_FLUC);
        if (iter >= LMTL_MU_ITER && convg) break;
    }
    
    return true;
}


Bool_t PhyTrFit::evolve(const MuOpt& mu_opt) {
    const Bool_t is_mu_free = (MuOpt::kFree == mu_opt);
    const Short_t DIMM = (is_mu_free ? 1 : 0);
    
    // Number of Res and Par
    Short_t numOfRes = (nseq_ + nseg_*DIML);
    Short_t numOfPar = (DIMG + DIMM + nseg_*DIML);

    Short_t parIDigb  = -1;
    Short_t parIDtsft = -1;
    if (is_mu_free) { parIDigb = DIMG; }
    if (nmes_TOFt_ >= LMTN_TOF_T) { parIDtsft = numOfPar; numOfPar += 1; }
    
    // Seg
    Bool_t hasSeg = (nseg_ != 0);
   
    // Final State
    std::vector<PhySt> stts;

    // Reset TOF Time and Path
    Bool_t resetTOF = true;
    HitStTOF::SetOffsetTime(Numc::ZERO<>);
    HitStTOF::SetOffsetPath(Numc::ZERO<>);
    HitStTOF::SetTimeShiftCorr(true);
    
    // TOF time shift
    Double_t TOFt_sft = ((nmes_TOFt_ >= LMTN_TOF_T) ? TOFt_sft_ : Numc::ZERO<>);
   
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
    for (Short_t is = 0; is < nseg_; ++is) { // Interaction
        SVecD<5> inrm, idiv;
        args_.at(is).cal_nrm_and_div(inrm, idiv);
        
        chi_cx += inrm(0) * inrm(0); // tauu
        chi_cx += inrm(2) * inrm(2); // taul
        chi_cy += inrm(1) * inrm(1); // rhou
        chi_cy += inrm(3) * inrm(3); // rhol

        for (Short_t it = 0; it < DIML; ++it)
            jb(nseq_+is*DIML+it, DIMG+DIMM+is*DIML+it) += idiv(it);
    } // Interaction
    
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
        hit->cal(ppst, noise_ctler_);
        
        // Update Jacb
        jbGG = curjb.gg() * jbGG;
        if (hasLoc && hasCxy) jbGL.at(cnt_nseg) = curjb.gl();
        for (Short_t is = 0; is < cnt_nseg; ++is)
            jbGL.at(is) = curjb.gg() * jbGL.at(is);
        
        // State
        if (hasCxy) stts.push_back(ppst);

        // Coord
        if (hit->scx()) chi_cx += hit->nrmcx() * hit->nrmcx(); 
        if (hit->scy()) chi_cy += hit->nrmcy() * hit->nrmcy();
        for (Short_t it = 0; it < DIMG; ++it) {
            if (hit->scx()) jb(hit->seqIDcx(), it) += hit->divcx() * jbGG(0, it);
            if (hit->scy()) jb(hit->seqIDcy(), it) += hit->divcy() * jbGG(1, it);
        }
        if (hasLoc) {
            for (Short_t is = 0; is <= itnseg; ++is) {
                for (Short_t it = 0; it < DIML; ++it) {
                    if (hit->scx()) jb(hit->seqIDcx(), DIMG+DIMM+is*DIML+it) += hit->divcx() * jbGL.at(is)(0, it);
                    if (hit->scy()) jb(hit->seqIDcy(), DIMG+DIMM+is*DIML+it) += hit->divcy() * jbGL.at(is)(1, it);
                }
            }
        } // Local
        
         
        // TRK
        HitStTRK* hitTRK = Hit<HitStTRK>::Cast(hit);
        if (hitTRK != nullptr) {
            if (hitTRK->sqx()) chi_ib += hitTRK->nrmqx() * hitTRK->nrmqx();
            if (hitTRK->sqy()) chi_ib += hitTRK->nrmqy() * hitTRK->nrmqy();
            if (hitTRK->sqx()) {
                if (is_mu_free) jb(hitTRK->seqIDqx(), parIDigb) += hitTRK->divqx_igb() * jbGG(4, 4);
                else            jb(hitTRK->seqIDqx(),        4) += hitTRK->divqx_eta() * jbGG(4, 4);
            }
            if (hitTRK->sqy()) {
                if (is_mu_free) jb(hitTRK->seqIDqy(), parIDigb) += hitTRK->divqy_igb() * jbGG(4, 4);
                else            jb(hitTRK->seqIDqy(),        4) += hitTRK->divqy_eta() * jbGG(4, 4);
            }    
        }

        // TOF
        HitStTOF* hitTOF = Hit<HitStTOF>::Cast(hit);
        if (hitTOF != nullptr) {
            if (hitTOF->st()) chi_ib += hitTOF->nrmt() * hitTOF->nrmt();
            if (hitTOF->sq()) chi_ib += hitTOF->nrmq() * hitTOF->nrmq();
            if (hitTOF->st()) {
                if (is_mu_free) jb(hitTOF->seqIDt(), parIDigb) += hitTOF->divt_igb() * jbGG(4, 4);
                else            jb(hitTOF->seqIDt(),        4) += hitTOF->divt_eta() * jbGG(4, 4);
            }
            if (hitTOF->sq()) {
                if (is_mu_free) jb(hitTOF->seqIDq(), parIDigb) += hitTOF->divq_igb() * jbGG(4, 4);
                else            jb(hitTOF->seqIDq(),        4) += hitTOF->divq_eta() * jbGG(4, 4);
            }
            if (hitTOF->st() && (nmes_TOFt_ >= LMTN_TOF_T)) jb(hitTOF->seqIDt(), parIDtsft) += hitTOF->divt_sft(); // TOF time shift
        }
        
        // RICH
        HitStRICH* hitRICH = Hit<HitStRICH>::Cast(hit);
        if (hitRICH != nullptr) {
            if (hitRICH->sib()) chi_ib += hitRICH->nrmib() * hitRICH->nrmib();
            if (hitRICH->sib()) {
                if (is_mu_free) jb(hitRICH->seqIDib(), parIDigb) += hitRICH->divib_igb() * jbGG(4, 4);
                else            jb(hitRICH->seqIDib(),        4) += hitRICH->divib_eta() * jbGG(4, 4);
            }
        }
        
        // TRD
        HitStTRD* hitTRD = Hit<HitStTRD>::Cast(hit);
        if (hitTRD != nullptr) {
            if (hitTRD->sel()) chi_ib += hitTRD->nrmel() * hitTRD->nrmel();
            if (hitTRD->sel()) {
                if (is_mu_free) jb(hitTRD->seqIDel(), parIDigb) += hitTRD->divel_igb() * jbGG(4, 4);
                else            jb(hitTRD->seqIDel(),        4) += hitTRD->divel_eta() * jbGG(4, 4);
            }
        }
        
        if (hasCxy) {
            nearPpst = ppst;
            nearJbGG = jbGG;
            nearJbGL = jbGL;
            if (hasSeg) cnt_nseg++;
        }
        cnt_nhit++;
    }
    if (cnt_nhit != hits_.size()) return false;
    if (hasSeg && cnt_nseg != nseg_) return false;

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
    if (is_mu_free) {
        Double_t reEta = std::fabs(errs.at(4) / part_.eta());
        Double_t reIgb = (errs.at(parIDigb) / part_.igmbta());
        Double_t reMu  = std::sqrt((reEta * reEta) + (reIgb * reIgb));
        errMu = reMu * part_.mu();
        if (!Numc::Valid(errMu) || Numc::Compare(errMu) < 0) errMu = Numc::ZERO<>;
    }

    err_.fill(Numc::ZERO<>);
    err_.at(0) = errs.at(0);
    err_.at(1) = errs.at(1);
    err_.at(2) = errs.at(2);
    err_.at(3) = errs.at(3);
    err_.at(4) = errs.at(4);
    if (is_mu_free) err_.at(5) = errs.at(parIDigb);
    if (is_mu_free) err_.at(6) = errMu;

    stts_ = stts;
    info_ = part_.info();
 
    //CERR("MASS %14.8f RIG %14.8f QLT %14.8f %14.8f\n", part_.mass(), part_.rig(), quality_.at(0), quality_.at(1));
    return true;
}


bool VirtualPhyTrFit::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const {
    if (numOfRes_ <= 0 || numOfPar_ <= 0) return false;
    std::fill_n(residuals, numOfRes_, Numc::ZERO<>);
    Bool_t hasJacb = (jacobians != nullptr && jacobians[0] != nullptr);
    if (hasJacb) std::fill_n(jacobians[0], numOfRes_ * numOfPar_, Numc::ZERO<>);
    
    // Seg
    Bool_t hasSeg = (nseg_ != 0);

    // Reset TOF Time and Path
    Bool_t resetTOF = true;
    HitStTOF::SetOffsetTime(Numc::ZERO<>);
    HitStTOF::SetOffsetPath(Numc::ZERO<>);
    HitStTOF::SetTimeShiftCorr(true);

    // TOF time shift
    Double_t TOFt_sft = ((nmes_TOFt_ >= LMTN_TOF_T) ? parameters[0][parIDtsft_] : Numc::ZERO<>);

    // Particle Status
    PhySt ppst(part_);
    Double_t partmu = (is_mu_free_ ? std::fabs(parameters[0][parIDigb_] / parameters[0][4]) : part_.mu());
    if (is_mu_free_) ppst.reset(partmu);
    else             ppst.arg().clear();
    ppst.set_state_with_uxy(parameters[0][0], parameters[0][1], part_.cz(), parameters[0][2], parameters[0][3], Numc::Compare(part_.uz()));
    ppst.set_eta(parameters[0][4]);

    // Interaction
    std::vector<PhyArg> args(nseg_, PhyArg(sw_mscat_, sw_eloss_));
    for (Short_t is = 0; is < nseg_; ++is) {
        args.at(is).set_mscat(
            parameters[0][DIMG_+DIMM_+is*DIML_+0], 
            parameters[0][DIMG_+DIMM_+is*DIML_+1], 
            parameters[0][DIMG_+DIMM_+is*DIML_+2], 
            parameters[0][DIMG_+DIMM_+is*DIML_+3]);
    }

    // Matrix (Rs, Jb)
    PhyJb::SMtxDGG&& jbGG = SMtxId();
    std::vector<PhyJb::SMtxDGL> jbGL(nseg_);
    ceres::Vector rs = ceres::Vector::Zero(numOfRes_);
    ceres::Matrix jb = ceres::Matrix::Zero(numOfRes_, numOfPar_);

    for (Short_t is = 0; is < nseg_; ++is) { // Interaction
        SVecD<5> inrm, idiv;
        args.at(is).cal_nrm_and_div(inrm, idiv);

        for (Short_t it = 0; it < DIML_; ++it)
            rs(nseq_+is*DIML_+it) += inrm(it);

        if (hasJacb) {
            for (Short_t it = 0; it < DIML_; ++it)
                jb(nseq_+is*DIML_+it, DIMG_+DIMM_+is*DIML_+it) += idiv(it);
        } // hasJacb
    } // Interaction
    
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
        if (!PropMgnt::PropToZ(hit->cz(), ppst, nullptr, ((hasJacb)?&curjb:nullptr)))  break;
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
        hit->cal(ppst, noise_ctler_);

        // Update Jacb
        if (hasJacb) {
            jbGG = curjb.gg() * jbGG;
            if (hasLoc && hasCxy) jbGL.at(cnt_nseg) = curjb.gl();
            for (Short_t is = 0; is < cnt_nseg; ++is)
                jbGL.at(is) = curjb.gg() * jbGL.at(is);
        }
      
        // Coord
        if (hit->scx()) rs(hit->seqIDcx()) += hit->nrmcx();
        if (hit->scy()) rs(hit->seqIDcy()) += hit->nrmcy();
        if (hasJacb) {
            for (Short_t it = 0; it < DIMG_; ++it) {
                if (hit->scx()) jb(hit->seqIDcx(), it) += hit->divcx() * jbGG(0, it);
                if (hit->scy()) jb(hit->seqIDcy(), it) += hit->divcy() * jbGG(1, it);
            }
            if (hasLoc) {
                for (Short_t is = 0; is <= itnseg; ++is) {
                    for (Short_t it = 0; it < DIML_; ++it) {
                        if (hit->scx()) jb(hit->seqIDcx(), DIMG_+DIMM_+is*DIML_+it) += hit->divcx() * jbGL.at(is)(0, it);
                        if (hit->scy()) jb(hit->seqIDcy(), DIMG_+DIMM_+is*DIML_+it) += hit->divcy() * jbGL.at(is)(1, it);
                    }
                }
            } // Local
        } // hasJacb 

        // TRK
        HitStTRK* hitTRK = Hit<HitStTRK>::Cast(hit);
        if (hitTRK != nullptr) {
            if (hitTRK->sqx()) rs(hitTRK->seqIDqx()) += hitTRK->nrmqx();
            if (hitTRK->sqy()) rs(hitTRK->seqIDqy()) += hitTRK->nrmqy();
            if (hasJacb && hitTRK->sqx()) {
                if (is_mu_free_) jb(hitTRK->seqIDqx(), parIDigb_) += hitTRK->divqx_igb() * jbGG(4, 4);
                else             jb(hitTRK->seqIDqx(),         4) += hitTRK->divqx_eta() * jbGG(4, 4);
            } // hasJacb
            if (hasJacb && hitTRK->sqy()) {
                if (is_mu_free_) jb(hitTRK->seqIDqy(), parIDigb_) += hitTRK->divqy_igb() * jbGG(4, 4);
                else             jb(hitTRK->seqIDqy(),         4) += hitTRK->divqy_eta() * jbGG(4, 4);
            } // hasJacb
        }
        
        // TOF
        HitStTOF* hitTOF = Hit<HitStTOF>::Cast(hit);
        if (hitTOF != nullptr) {
            if (hitTOF->sq()) rs(hitTOF->seqIDq()) += hitTOF->nrmq();
            if (hitTOF->st()) rs(hitTOF->seqIDt()) += hitTOF->nrmt();
            if (hasJacb && hitTOF->st()) {
                if (is_mu_free_) jb(hitTOF->seqIDt(), parIDigb_) += hitTOF->divt_igb() * jbGG(4, 4);
                else             jb(hitTOF->seqIDt(),         4) += hitTOF->divt_eta() * jbGG(4, 4);
            } // hasJacb
            if (hasJacb && hitTOF->sq()) {
                if (is_mu_free_) jb(hitTOF->seqIDq(), parIDigb_) += hitTOF->divq_igb() * jbGG(4, 4);
                else             jb(hitTOF->seqIDq(),         4) += hitTOF->divq_eta() * jbGG(4, 4);
            } // hasJacb
            if (hasJacb && hitTOF->st() && (nmes_TOFt_ >= LMTN_TOF_T)) jb(hitTOF->seqIDt(), parIDtsft_) += hitTOF->divt_sft(); // TOF time shift
        }
        
        // RICH
        HitStRICH* hitRICH = Hit<HitStRICH>::Cast(hit);
        if (hitRICH != nullptr) {
            if (hitRICH->sib()) rs(hitRICH->seqIDib()) += hitRICH->nrmib();
            if (hasJacb && hitRICH->sib()) {
                if (is_mu_free_) jb(hitRICH->seqIDib(), parIDigb_) += hitRICH->divib_igb() * jbGG(4, 4);
                else             jb(hitRICH->seqIDib(),         4) += hitRICH->divib_eta() * jbGG(4, 4);
            } // hasJacb
        }
        
        // TRD
        HitStTRD* hitTRD = Hit<HitStTRD>::Cast(hit);
        if (hitTRD != nullptr) {
            if (hitTRD->sel()) rs(hitTRD->seqIDel()) += hitTRD->nrmel();
            if (hasJacb && hitTRD->sel()) {
                if (is_mu_free_) jb(hitTRD->seqIDel(), parIDigb_) += hitTRD->divel_igb() * jbGG(4, 4);
                else             jb(hitTRD->seqIDel(),         4) += hitTRD->divel_eta() * jbGG(4, 4);
            } // hasJacb
        }

        if (hasCxy) {
            nearPpst = ppst;
            nearJbGG = jbGG;
            nearJbGL = jbGL;
            if (hasSeg) cnt_nseg++;
        }
        cnt_nhit++;
    }
    if (cnt_nhit != hits_.size()) return false;
    if (hasSeg && cnt_nseg != nseg_) return false;
    
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


#endif // __TRACKLibs_PhyFit_C__
