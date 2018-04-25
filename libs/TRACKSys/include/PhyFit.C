#ifndef __TRACKLibs_PhyFit_C__
#define __TRACKLibs_PhyFit_C__


namespace TrackSys {
        
    
void TrFitPar::print() const {
    std::string outstr;
    outstr += "TrFitPar::Print()\n";
    for (auto&& hit : hits_TRK_) outstr += STR("TRK L%d COO %14.8f %14.8f %14.8f CERR %14.8f %14.8f\n", hit.lay(), hit.cx(), hit.cy(), hit.cz(), hit.ecx(), hit.ecy());
    for (auto&& hit : hits_TOF_) outstr += STR("TOF L%d COO %14.8f %14.8f %14.8f CERR %14.8f %14.8f\n", hit.lay(), hit.cx(), hit.cy(), hit.cz(), hit.ecx(), hit.ecy());
    outstr += "\n";
    COUT("%s", outstr.c_str());
}
        
TrFitPar& TrFitPar::operator=(const TrFitPar& rhs) {
    if (this != &rhs) {
        sw_mscat_      = rhs.sw_mscat_;
        sw_eloss_      = rhs.sw_eloss_;
        info_          = rhs.info_;
        ortt_          = rhs.ortt_;
        hits_TRK_      = rhs.hits_TRK_;
        hits_TOF_      = rhs.hits_TOF_;
        nseq_          = rhs.nseq_;
        nmes_          = rhs.nmes_;
        nmes_cx_       = rhs.nmes_cx_;
        nmes_cy_       = rhs.nmes_cy_;
        nmes_TRKqx_    = rhs.nmes_TRKqx_;
        nmes_TRKqy_    = rhs.nmes_TRKqy_;
        nmes_TOFq_     = rhs.nmes_TOFq_;
        nmes_TOFt_     = rhs.nmes_TOFt_;
        is_check_      = rhs.is_check_;
        
        hits_.clear();
        if (is_check_) {
            for (auto&& hit : hits_TRK_) hits_.push_back(&hit); 
            for (auto&& hit : hits_TOF_) hits_.push_back(&hit); 
            
            if (ortt_ == Orientation::kDownward) VirtualHitSt::Sort(hits_, VirtualHitSt::Orientation::kDownward);
            else                                 VirtualHitSt::Sort(hits_, VirtualHitSt::Orientation::kUpward);
        }
        else zero();
    }
   
    return *this;
}
    
TrFitPar::TrFitPar(const PartInfo& info, const Orientation& ortt, Bool_t sw_mscat, Bool_t sw_eloss) {
    clear();

    sw_mscat_ = sw_mscat;
    sw_eloss_ = sw_eloss;
    info_ = info;
    ortt_ = ortt;
}

void TrFitPar::zero() {
    hits_.clear();

    nseq_       = 0;
    nmes_       = 0;
    nmes_cx_    = 0;
    nmes_cy_    = 0;
    nmes_TRKqx_ = 0;
    nmes_TRKqy_ = 0;
    nmes_TOFq_  = 0;
    nmes_TOFt_  = 0;
    
    is_check_ = false;
}

void TrFitPar::clear() {
    sw_mscat_ = false;
    sw_eloss_ = false;
    
    info_ = PartInfo(PartType::Proton);
    ortt_ = Orientation::kDownward;
   
    hits_TRK_.clear();
    hits_TOF_.clear();

    zero();
}

Bool_t TrFitPar::sort_hits() {
    zero();

    if (ortt_ == Orientation::kDownward) Hit<HitStTRK>::Sort(hits_TRK_, VirtualHitSt::Orientation::kDownward);
    else                                 Hit<HitStTRK>::Sort(hits_TRK_, VirtualHitSt::Orientation::kUpward);
    if (ortt_ == Orientation::kDownward) Hit<HitStTOF>::Sort(hits_TOF_, VirtualHitSt::Orientation::kDownward);
    else                                 Hit<HitStTOF>::Sort(hits_TOF_, VirtualHitSt::Orientation::kUpward);
 
    for (auto&& hit : hits_TRK_) hits_.push_back(&hit); 
    for (auto&& hit : hits_TOF_) hits_.push_back(&hit); 
    
    if (ortt_ == Orientation::kDownward) VirtualHitSt::Sort(hits_, VirtualHitSt::Orientation::kDownward);
    else                                 VirtualHitSt::Sort(hits_, VirtualHitSt::Orientation::kUpward);
    
    Short_t nseq = 0;
    for (auto&& hit : hits_) {
        hit->set_type(info_);
        nseq += hit->set_seqID(nseq);
        
        if (hit->scx()) nmes_cx_++;
        if (hit->scy()) nmes_cy_++;
    }

    for (auto&& hit : hits_TRK_) {
        if (hit.sqx()) nmes_TRKqx_++;
        if (hit.sqy()) nmes_TRKqy_++;
    }
        
    for (auto&& hit : hits_TOF_) {
        if (hit.sq()) nmes_TOFq_++;
        if (hit.st()) nmes_TOFt_++;
    }

    nseq_ = nseq;
    nmes_ = (nmes_cx_ + nmes_cy_) + (nmes_TRKqx_ + nmes_TRKqy_) + (nmes_TOFq_ + nmes_TOFt_);

    return true;
}

Bool_t TrFitPar::check_hits() {
    if (is_check_) return is_check_;
    sort_hits();
   
    Bool_t passed = (nmes_cx_ >= LMTN_CX && nmes_cy_ >= LMTN_CY);
    if (passed) is_check_ = true;
   
    return is_check_;
}

SimpleTrFit& SimpleTrFit::operator=(const SimpleTrFit& rhs) {
    if (this != &rhs) {
        dynamic_cast<TrFitPar&>(*this) = dynamic_cast<const TrFitPar&>(rhs);
        succ_      = rhs.succ_;
        part_      = rhs.part_;
        ndof_      = rhs.ndof_;
        ndof_cx_   = rhs.ndof_cx_;
        ndof_cy_   = rhs.ndof_cy_;
        ndof_TRKq_ = rhs.ndof_TRKq_;
        ndof_TOFq_ = rhs.ndof_TOFq_;
        ndof_TOFt_ = rhs.ndof_TOFt_;
        nchi_      = rhs.nchi_;
        nchi_cx_   = rhs.nchi_cx_;
        nchi_cy_   = rhs.nchi_cy_;
        nchi_TRKq_ = rhs.nchi_TRKq_;
        nchi_TOFq_ = rhs.nchi_TOFq_;
        nchi_TOFt_ = rhs.nchi_TOFt_;
        if (!succ_) clear();
    }
    return *this;
}

SimpleTrFit::SimpleTrFit(const TrFitPar& fitPar) : TrFitPar(fitPar) {
    SimpleTrFit::clear();
    if (!check_hits()) return;
   
    ndof_cx_   = (nmes_cx_ >= LMTN_CX) ? (nmes_cx_ - Numc::TWO<Short_t>  ) : 0;
    ndof_cy_   = (nmes_cy_ >= LMTN_CY) ? (nmes_cy_ - Numc::THREE<Short_t>) : 0;
    ndof_TRKq_ = nmes_TRKqx_ + nmes_TRKqy_;
    ndof_TOFq_ = nmes_TOFq_;
    ndof_TOFt_ = (nmes_TOFt_ >= LMTN_TOF_T) ? (nmes_TOFt_ - Numc::ONE<Short_t>) : 0;
    ndof_      = (ndof_cx_ + ndof_cy_ + ndof_TRKq_ + ndof_TOFq_ + ndof_TOFt_);
    if (ndof_ == Numc::ZERO<Short_t>) return;

    succ_ = (analyticalFit() ? simpleFit() : false);
    
    if (succ_) info_ = part_.info();
    else { SimpleTrFit::clear(); COUT("SimpleTrFit::FAIL.\n"); }
}


void SimpleTrFit::clear() {
    sw_mscat_ = false;
    sw_eloss_ = false;
    
    succ_ = false;
    part_.reset(info_);
    part_.arg().reset(sw_mscat_, sw_eloss_);

    ndof_      = 0;
    ndof_cx_   = 0;
    ndof_cy_   = 0;
    ndof_TRKq_ = 0;
    ndof_TOFq_ = 0;
    ndof_TOFt_ = 0;
        
    nchi_      = 0;
    nchi_cx_   = 0;
    nchi_cy_   = 0;
    nchi_TRKq_ = 0;
    nchi_TOFq_ = 0;
    nchi_TOFt_ = 0;
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
        Double_t Lambda = PROP_FACT * part_.info().chrg_to_mass(); 
        
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
        Bool_t resetTOF = true;
        HitStTOF::SetOffsetTime(Numc::ZERO<>);
        HitStTOF::SetOffsetPath(Numc::ZERO<>);
        
        Double_t chi_cx   = 0;
        Double_t chi_cy   = 0;
        Double_t chi_TRKq = 0;
        Double_t chi_TOFq = 0;
        Double_t chi_TOFt = 0;
        
        SVecD<5>    grdG;
        SMtxSymD<5> cvGG;

        UInt_t cnt_nhit = 0;
        PhySt ppst(rltSt);
        SMtxD<5>&& ppjb = SMtxId();
        for (auto&& hit : hits_) {
            PhyJb curjb;
            if (!PropMgnt::PropToZ(hit->cz(), ppst, nullptr, &curjb)) break;
            ppjb = (curjb.gg()).Sub<SMtxD<5>>(0, 0) * ppjb;
        
            if (resetTOF && Hit<HitStTOF>::IsSame(hit)) { // set reference
                HitStTOF::SetOffsetPath(ppst.path());
                HitStTOF::SetOffsetTime(ppst.time()-Hit<HitStTOF>::Cast(hit)->t());
                resetTOF = false;
            }
            hit->cal(ppst);

            SVecD<2>    rsC;
            SMtxD<2, 5> jbC;
            if (hit->scx()) rsC(0) += hit->nrmcx();
            if (hit->scy()) rsC(1) += hit->nrmcy();
            for (UInt_t it = 0; it < 5; ++it) {
                if (hit->scx()) jbC(0, it) += hit->divcx() * ppjb(0, it);
                if (hit->scy()) jbC(1, it) += hit->divcy() * ppjb(1, it);
            }
           
            grdG += LA::Transpose(jbC) * rsC;
            cvGG += LA::SimilarityT(jbC, SMtxSymD<2>(SMtxId()));
            if (hit->scx()) chi_cx += rsC(0) * rsC(0); 
            if (hit->scy()) chi_cy += rsC(1) * rsC(1);
            
            HitStTRK* hitTRK = Hit<HitStTRK>::Cast(hit);
            if (hitTRK != nullptr) {
                SVecD<2>    rsTRK;
                SMtxD<2, 5> jbTRK;
                if (hitTRK->sqx()) rsTRK(0) += hitTRK->nrmqx();
                if (hitTRK->sqy()) rsTRK(1) += hitTRK->nrmqy();
                for (UInt_t it = 0; it < 5; ++it) {
                    if (hitTRK->sqx()) jbTRK(0, it) += hitTRK->divqx() * ppjb(4, it);
                    if (hitTRK->sqy()) jbTRK(1, it) += hitTRK->divqy() * ppjb(4, it);
                }
                grdG += LA::Transpose(jbTRK) * rsTRK;
                cvGG += LA::SimilarityT(jbTRK, SMtxSymD<2>(SMtxId()));
                if (hitTRK->sqx()) chi_TRKq += rsTRK(0) * rsTRK(0);
                if (hitTRK->sqy()) chi_TRKq += rsTRK(1) * rsTRK(1);
            }

            HitStTOF* hitTOF = Hit<HitStTOF>::Cast(hit);
            if (hitTOF != nullptr) {
                SVecD<2>    rsTOF;
                SMtxD<2, 5> jbTOF;
                if (hitTOF->sq()) rsTOF(0) += hitTOF->nrmq();
                if (hitTOF->st()) rsTOF(1) += hitTOF->nrmt();
                for (UInt_t it = 0; it < 5; ++it) {
                    if (hitTOF->sq()) jbTOF(0, it) += hitTOF->divq() * ppjb(4, it);
                    if (hitTOF->st()) jbTOF(1, it) += hitTOF->divt() * ppjb(4, it);
                }
                grdG += LA::Transpose(jbTOF) * rsTOF;
                cvGG += LA::SimilarityT(jbTOF, SMtxSymD<2>(SMtxId()));
                if (hitTOF->sq()) chi_TOFq += rsTOF(0) * rsTOF(0);
                if (hitTOF->st()) chi_TOFt += rsTOF(1) * rsTOF(1);
            }

            cnt_nhit++;
        }
        if (cnt_nhit != hits_.size()) break;
        Double_t chi  = (chi_cx + chi_cy + chi_TRKq + chi_TOFq + chi_TOFt);
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
                lambda     = std::max(lambda/LAMBDA_DN_FAC, LMTL_LAMBDA);
                nchi_cx_   = ((ndof_cx_   > 0) ? (chi_cx   / static_cast<Double_t>(ndof_cx_  )) : 0);
                nchi_cy_   = ((ndof_cy_   > 0) ? (chi_cy   / static_cast<Double_t>(ndof_cy_  )) : 0);
                nchi_TRKq_ = ((ndof_TRKq_ > 0) ? (chi_TRKq / static_cast<Double_t>(ndof_TRKq_)) : 0);
                nchi_TOFq_ = ((ndof_TOFq_ > 0) ? (chi_TOFq / static_cast<Double_t>(ndof_TOFq_)) : 0);
                nchi_TOFt_ = ((ndof_TOFt_ > 0) ? (chi_TOFt / static_cast<Double_t>(ndof_TOFt_)) : 0);
                nchi_      = nchi;
                part_      = rltSt;
                isUpdate   = true;
                updIter++;
            }
        }
        else { nchi_ = nchi; }

        SMtxSymD<5> lmCvGG(cvGG);
        SVecD<5>&&  diagCvGG = (lambda * cvGG.Diagonal());
        lmCvGG.SetDiagonal(SVecD<5>(lmCvGG.Diagonal() + diagCvGG));

        Bool_t isNomag = (Numc::EqualToZero(lmCvGG(4, 4)) && Numc::EqualToZero(grdG(4))); // Fast Check
        if (isNomag) {
            grdG(4) = Numc::ZERO<>;
            grdG(5) = Numc::ZERO<>;
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
            ((ortt_ == Orientation::kDownward) ? Numc::NEG<Short_t> : Numc::ONE<Short_t>)
        );
        rltSt.set_eta(rltSt.eta() - rslG(4));
        
        preSucc = curSucc;
        curSucc = (isSucc && updIter >= LMTL_ITER);
        succ    = (preSucc && curSucc);
        
        if (!succ) curIter++;
    }
    //if (!succ) COUT("FAIL. IT %2d %2d (RIG %14.8f MASS %14.8f CHI %14.8f) LAMBDA %14.8f\n", curIter, updIter, part_.rig(), part_.mass(), nchi_, lambda);
    //else       COUT("SUCC. IT %2d %2d (RIG %14.8f MASS %14.8f CHI %14.8f) LAMBDA %14.8f\n", curIter, updIter, part_.rig(), part_.mass(), nchi_, lambda);
    
    return succ;
}

PhyTrFit& PhyTrFit::operator=(const PhyTrFit& rhs) {
    if (this != &rhs) {
        dynamic_cast<TrFitPar&>(*this) = dynamic_cast<const TrFitPar&>(rhs);
        succ_      = rhs.succ_;
        part_      = rhs.part_;
        args_      = rhs.args_;
        nsegs_     = rhs.nsegs_;
        
        ndof_      = rhs.ndof_;
        ndof_cx_   = rhs.ndof_cx_;
        ndof_cy_   = rhs.ndof_cy_;
        ndof_TRKq_ = rhs.ndof_TRKq_;
        ndof_TOFq_ = rhs.ndof_TOFq_;
        ndof_TOFt_ = rhs.ndof_TOFt_;
        nchi_      = rhs.nchi_;
        nchi_cx_   = rhs.nchi_cx_;
        nchi_cy_   = rhs.nchi_cy_;
        nchi_TRKq_ = rhs.nchi_TRKq_;
        nchi_TOFq_ = rhs.nchi_TOFq_;
        nchi_TOFt_ = rhs.nchi_TOFt_;

        nrm_mstau_ = rhs.nrm_mstau_;
        nrm_msrho_ = rhs.nrm_msrho_;
        nrm_elion_ = rhs.nrm_elion_;
        
        if (succ_ && hits_.size() == stts_.size()) {
            for (UInt_t it = 0; it < stts_.size(); ++it) {
                stts_.at(it).first = hits_.at(it);
            }
        }
        else clear();
    }
    return *this;
}

PhyTrFit::PhyTrFit(const TrFitPar& fitPar, Bool_t is_fixed_mass) : TrFitPar(fitPar) {
    PhyTrFit::clear();
    if (!check_hits()) return;
    nsegs_ = (nhits() >= Numc::TWO<Short_t>) ? (nhits() - Numc::ONE<Short_t>) : 0;
    
    ndof_cx_   = (nmes_cx_ >= LMTN_CX) ? (nmes_cx_ - Numc::TWO<Short_t>  ) : 0;
    ndof_cy_   = (nmes_cy_ >= LMTN_CY) ? (nmes_cy_ - Numc::THREE<Short_t>) : 0;
    ndof_TRKq_ = nmes_TRKqx_ + nmes_TRKqy_;
    ndof_TOFq_ = nmes_TOFq_;
    ndof_TOFt_ = (nmes_TOFt_ >= LMTN_TOF_T) ? (nmes_TOFt_ - Numc::ONE<Short_t>) : 0;
    ndof_      = (ndof_cx_ + ndof_cy_ + ndof_TRKq_ + ndof_TOFq_ + ndof_TOFt_ - (is_fixed_mass?0:1));
    if (ndof_ == Numc::ZERO<Short_t>) return;
    Short_t absChrg = std::abs(info_.chrg());

    // Fixed Mass
    if (is_fixed_mass) succ_ = (simpleFit() ? physicalFit() : false);
    // Free Mass
    else if (absChrg > Numc::ZERO<Short_t> && absChrg < PartListMassQ.size()) {
        auto&& compless = [] (const Double_t& a, const Double_t& b) { return ((a > Numc::ZERO<>) ? ((b > Numc::ZERO<>) ? a < b : a > Numc::ZERO<>) : false); }; // return true, if a < b or only a
        
        // List of Mass
        const std::vector<Double_t>&     listMass = PartListMassQ.at(absChrg);
        std::vector<Bool_t>              refSucc(listMass.size(), false);
        std::vector<Double_t>            refNchi(listMass.size(), Numc::NEG<>);
        std::vector<PhySt>               refSt;
        std::vector<std::vector<PhyArg>> refArgs;
        for (UInt_t it = 0; it < listMass.size(); ++it) {
            info_.reset(absChrg, listMass.at(it)); 
            Bool_t succ = (simpleFit() ? physicalFit() : false);
            refSt.push_back(part_);
            refArgs.push_back(args_);
            if (!succ) continue;
            refSucc.at(it) = true;
            refNchi.at(it) = nchi_;
        }

        // Min NormChisq
        const Double_t LMT_NCHI = 0.05;
        Double_t minNchi = *std::min_element(refNchi.begin(), refNchi.end(), compless);
        
        // Parameters
        Double_t            nchi = Numc::NEG<>;
        PhySt               part = part_;
        std::vector<PhyArg> args = args_;
        
        Int_t seedIter = 0;
        const std::vector<Double_t> seedVec({ 0.3, 0.7, 0.7 });
        const std::vector<Double_t> qualityVec({ Numc::THREE<>/Numc::FIVE<>, Numc::FOUR<>/Numc::FIVE<>, Numc::FIVE<>/Numc::SIX<> });
        const std::vector<Double_t> flucVec({ 0.02, 0.05, 0.05 });
        while (Numc::Compare(minNchi) > 0 && seedIter < seedVec.size()) {
            CERR("ITER %d MIN_NCHI %14.8f\n", seedIter, minNchi);
            std::vector<Double_t> listVal;
            Short_t consistSelfCnt = Numc::ZERO<Int_t>;
            for (UInt_t it = 0; it < listMass.size(); ++it) {
                if (!refSucc.at(it)) continue;
                Double_t quality = (Numc::TWO<> * (minNchi + LMT_NCHI) / (minNchi + refNchi.at(it) + Numc::TWO<> * LMT_NCHI)); 
                if (quality < qualityVec.at(seedIter)) { refSucc.at(it) = false; refNchi.at(it) = Numc::NEG<>; continue; }

                part_ = refSt.at(it);
                args_ = refArgs.at(it);
                Bool_t succ = physicalFit(is_fixed_mass, seedVec.at(seedIter));
                if (!succ) { refSucc.at(it) = false; refNchi.at(it) = Numc::NEG<>; continue; }
                Double_t fluc = std::fabs(refNchi.at(it) - nchi_) / (refNchi.at(it) + nchi_);
                if (fluc < flucVec.at(seedIter)) consistSelfCnt++;
                
                refNchi.at(it) = nchi_;
                listVal.push_back(nchi_);
                Bool_t firstTime = (Numc::Compare(nchi) < 0);
                if (firstTime || Numc::Compare(nchi_, nchi) < 0) {
                    nchi = nchi_;
                    part = part_;
                    args = args_;
                }
            }
            if (listVal.size() == 0) break;
            auto minmax = std::minmax_element(listVal.begin(), listVal.end());
            if (seedIter == 0) minNchi = *minmax.first;
            else               minNchi = std::min(minNchi, *minmax.first);
            
            Bool_t consistSelf = (listVal.size() == consistSelfCnt);
            if (consistSelf) break;
            
            if (listVal.size() >= 2) {
                Bool_t consistGroup = true;
                Double_t normFluc = Numc::TWO<> / (*minmax.first + *minmax.second);
                for (UInt_t it = 0; it < listMass.size(); ++it) {
                    if (!refSucc.at(it)) continue;
                    Double_t fluc = std::fabs(normFluc * refNchi.at(it) - Numc::ONE<>);
                    if (fluc > flucVec.at(seedIter)) { consistGroup = false; break; }
                }
                if (consistGroup) break;
            }

            seedIter++;
        }

        // Final
        CERR("FINE\n");
        if (Numc::Compare(nchi) > 0) {
            part_ = part;
            args_ = args;
            succ_ = evolve();
        }
    }

    if (succ_) info_ = part_.info();
    else { PhyTrFit::clear(); COUT("PhyTrFit::FAIL.\n"); }
}


void PhyTrFit::clear() {
    succ_ = false;
    part_.reset(info_);
    part_.arg().reset(sw_mscat_, sw_eloss_);
    args_.clear();
    stts_.clear();
    nsegs_     = 0;

    ndof_      = 0;
    ndof_cx_   = 0;
    ndof_cy_   = 0;
    ndof_TRKq_ = 0;
    ndof_TOFq_ = 0;
    ndof_TOFt_ = 0;
        
    nchi_      = 0;
    nchi_cx_   = 0;
    nchi_cy_   = 0;
    nchi_TRKq_ = 0;
    nchi_TOFq_ = 0;
    nchi_TOFt_ = 0;

    nrm_mstau_ = 0;
    nrm_msrho_ = 0;
    nrm_elion_ = 0;
}


Bool_t PhyTrFit::simpleFit() {
    SimpleTrFit simple(dynamic_cast<TrFitPar&>(*this));
    if (simple.status()) part_ = simple.part();
    return simple.status();
}


Bool_t PhyTrFit::physicalFit(Bool_t is_fixed_mass, Double_t scl) {
    Bool_t resetArg = (args_.size() != nhits()-1);
    part_.arg().reset(sw_mscat_, sw_eloss_);
    if (resetArg) args_ = std::move(std::vector<PhyArg>(nhits()-1, PhyArg(sw_mscat_, sw_eloss_)));
    
    std::vector<double> parameters({ part_.cx(), part_.cy(), part_.ux(), part_.uy(), part_.eta(), part_.info().invu() });
    std::vector<double> interaction_parameters((nhits()-1)*DIML, Numc::ZERO<>);
    if (!resetArg) {
        Bool_t is_scl = (Numc::Compare(scl) > 0);
        for (Int_t it = 0; it < nhits()-1; ++it) {
            const PhyArg& arg = args_.at(it);
            interaction_parameters.at(it*DIML+0) = arg.tauu() + (is_scl ? (scl * (Numc::TWO<> * (Rndm::DecimalUniform() - Numc::HALF) + Rndm::NormalGaussian())) : Numc::ZERO<>);
            interaction_parameters.at(it*DIML+1) = arg.rhou() + (is_scl ? (scl * (Numc::TWO<> * (Rndm::DecimalUniform() - Numc::HALF) + Rndm::NormalGaussian())) : Numc::ZERO<>);
            interaction_parameters.at(it*DIML+2) = arg.taul() + (is_scl ? (scl * (Numc::TWO<> * (Rndm::DecimalUniform() - Numc::HALF) + Rndm::NormalGaussian())) : Numc::ZERO<>);
            interaction_parameters.at(it*DIML+3) = arg.rhol() + (is_scl ? (scl * (Numc::TWO<> * (Rndm::DecimalUniform() - Numc::HALF) + Rndm::NormalGaussian())) : Numc::ZERO<>);
        }
    }
    parameters.insert(parameters.end(), interaction_parameters.begin(), interaction_parameters.end());

    ceres::CostFunction* cost_function = new VirtualPhyTrFit(dynamic_cast<TrFitPar&>(*this), part_, is_fixed_mass);

    ceres::Problem problem;
    problem.AddResidualBlock(cost_function, nullptr, parameters.data());
    problem.SetParameterLowerBound(parameters.data(), 2, -1.0);
    problem.SetParameterUpperBound(parameters.data(), 2,  1.0);
    problem.SetParameterLowerBound(parameters.data(), 3, -1.0);
    problem.SetParameterUpperBound(parameters.data(), 3,  1.0);
    if (!is_fixed_mass && part_.chrg() < PartListMassQ.size()) {
        Double_t lmtl_invu = (Numc::ONE<> / PartListMassQ.at(part_.chrg()).back() ) / Numc::TEN<>;
        Double_t lmtu_invu = (Numc::ONE<> / PartListMassQ.at(part_.chrg()).front()) * Numc::TEN<>;
        problem.SetParameterLowerBound(parameters.data(), 5, lmtl_invu);
        problem.SetParameterUpperBound(parameters.data(), 5, lmtu_invu);
    }
    
    ceres::Solver::Options options;
    //options.line_search_direction_type = ceres::LBFGS;
    
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    if (!summary.IsSolutionUsable()) return false;

    part_.set_state_with_uxy(parameters.at(0), parameters.at(1), part_.cz(), parameters.at(2), parameters.at(3), Numc::Compare(part_.uz()));
    part_.set_eta(parameters.at(4));
    if (!is_fixed_mass)
        part_.recal(parameters.at(5));

    interaction_parameters = std::vector<double>(parameters.begin()+DIMG, parameters.end());
    for (Int_t it = 0; it < nhits()-1; ++it) {
        args_.at(it).set_mscat(
                interaction_parameters.at(it*DIML+0), 
                interaction_parameters.at(it*DIML+1), 
                interaction_parameters.at(it*DIML+2), 
                interaction_parameters.at(it*DIML+3));
        //args_.at(it).set_eloss(
        //        interaction_parameters.at(it*DIML+4));
    }

    //std::cout << summary.FullReport() << std::endl;
    Bool_t succ = evolve();
    return succ;
}


Bool_t PhyTrFit::evolve() {
    std::vector<std::pair<VirtualHitSt*, PhySt>> stts;
    
    Bool_t resetTOF = true;
    HitStTOF::SetOffsetTime(Numc::ZERO<>);
    HitStTOF::SetOffsetPath(Numc::ZERO<>);
    
    Double_t chi_cx   = 0;
    Double_t chi_cy   = 0;
    Double_t chi_TRKq = 0;
    Double_t chi_TOFq = 0;
    Double_t chi_TOFt = 0;
            
    Double_t chi_mstau = 0;
    Double_t chi_msrho = 0;
    Double_t chi_elion = 0;

    PhySt ppst(part_);
    UInt_t cnt_nhit = 0;
    for (auto&& hit : hits_) {
        // Interaction Local Parameters
        Bool_t hasLoc = (cnt_nhit != 0);
        
        // Propagate
        if (hasLoc) ppst.arg() = args_.at(cnt_nhit-1);
        if (!PropMgnt::PropToZ(hit->cz(), ppst)) break;
       
        // Particle State
        if (hasLoc) args_.at(cnt_nhit-1).setvar_mat(ppst.arg().mat(), ppst.arg().nrl(), ppst.arg().ela());
        PhyArg curArg = ppst.arg();
        ppst.symbk();

        // Hit Status: Setting TOF reference time and path
        if (resetTOF && Hit<HitStTOF>::IsSame(hit)) { // set reference
            HitStTOF::SetOffsetPath(ppst.path());
            HitStTOF::SetOffsetTime(ppst.time()-Hit<HitStTOF>::Cast(hit)->t());
            resetTOF = false;
        }
        hit->cal(ppst);

        // State Pair
        stts.push_back(std::make_pair(hit, ppst));

        // Coord
        if (hit->scx()) chi_cx += hit->nrmcx() * hit->nrmcx(); 
        if (hit->scy()) chi_cy += hit->nrmcy() * hit->nrmcy();
            
        // TRK
        HitStTRK* hitTRK = Hit<HitStTRK>::Cast(hit);
        if (hitTRK != nullptr) {
            if (hitTRK->sqx()) chi_TRKq += hitTRK->nrmqx() * hitTRK->nrmqx();
            if (hitTRK->sqy()) chi_TRKq += hitTRK->nrmqy() * hitTRK->nrmqy();
        }

        // TOF
        HitStTOF* hitTOF = Hit<HitStTOF>::Cast(hit);
        if (hitTOF != nullptr) {
            if (hitTOF->sq()) chi_TOFq += hitTOF->nrmq() * hitTOF->nrmq();
            if (hitTOF->st()) chi_TOFt += hitTOF->nrmt() * hitTOF->nrmt();
        }
        
        // Interaction Local Parameters
        if (hasLoc) {
            SVecD<5> inrm;
            curArg.cal_nrm(inrm);
            chi_mstau += inrm(0) * inrm(0);
            chi_msrho += inrm(1) * inrm(1);
            chi_mstau += inrm(2) * inrm(2);
            chi_msrho += inrm(3) * inrm(3);
            chi_elion += inrm(4) * inrm(4);
        }

        cnt_nhit++;
    }
    if (cnt_nhit != hits_.size()) return false;
    Double_t chi  = (chi_cx + chi_cy + chi_TRKq + chi_TOFq + chi_TOFt + chi_mstau + chi_msrho + chi_elion);
    Double_t nchi = (chi / static_cast<Double_t>(ndof_));

    if (stts.size() != 0) (stts.at(stts.size()-1).second).arg().clear();
    for (UInt_t it = 0; it < stts.size()-1; ++it) {
        PhySt&  stt = (stts.at(it).second);
        PhyArg& arg = args_.at(it);
        stt.arg().set_mscat(arg.tauu(), arg.rhou(), arg.taul(), arg.rhol());
        stt.arg().set_eloss(arg.elion(), arg.elbrm());
    }

    nchi_cx_   = ((ndof_cx_   > 0) ? (chi_cx   / static_cast<Double_t>(ndof_cx_  )) : 0);
    nchi_cy_   = ((ndof_cy_   > 0) ? (chi_cy   / static_cast<Double_t>(ndof_cy_  )) : 0);
    nchi_TRKq_ = ((ndof_TRKq_ > 0) ? (chi_TRKq / static_cast<Double_t>(ndof_TRKq_)) : 0);
    nchi_TOFq_ = ((ndof_TOFq_ > 0) ? (chi_TOFq / static_cast<Double_t>(ndof_TOFq_)) : 0);
    nchi_TOFt_ = ((ndof_TOFt_ > 0) ? (chi_TOFt / static_cast<Double_t>(ndof_TOFt_)) : 0);
    nchi_      = nchi;

    nrm_mstau_ = ((nsegs_ > 0) ? (chi_mstau / static_cast<Double_t>(nsegs_)) : 0);
    nrm_msrho_ = ((nsegs_ > 0) ? (chi_msrho / static_cast<Double_t>(nsegs_)) : 0);
    nrm_elion_ = ((nsegs_ > 0) ? (chi_elion / static_cast<Double_t>(nsegs_)) : 0);

    stts_ = std::move(stts);
    info_ = part_.info();
    
    CERR("PhyTrFit::SUCC (RIG %14.8f MASS %14.8f NCHI %14.8f)\n", part_.rig(), part_.mass(), nchi_);
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
    ppst.arg().clear();

    if (!is_fixed_mass_)
        ppst.recal(parameters[0][5]);

    UInt_t cnt_nhit = 0;
    PhyJb::SMtxDGG&& jbGG = SMtxId();
    std::vector<PhyJb::SMtxDGL> jbGL(hits_.size()-1);
    ceres::Vector rs = ceres::Vector::Zero(numOfRes_);
    ceres::Matrix jb = ceres::Matrix::Zero(numOfRes_, numOfPar_);
    for (auto&& hit : hits_) {
        // Interaction Local Parameters
        Bool_t hasLoc = (cnt_nhit != 0);
       
        // Propagate
        PhyJb curjb;
        if (!PropMgnt::PropToZ(hit->cz(), ppst, nullptr, ((hasJacb)?&curjb:nullptr)))  break;
        if (hasLoc) {
            ppst.arg().set_mscat(
                parameters[0][DIMG+(cnt_nhit-1)*DIML+0], 
                parameters[0][DIMG+(cnt_nhit-1)*DIML+1], 
                parameters[0][DIMG+(cnt_nhit-1)*DIML+2], 
                parameters[0][DIMG+(cnt_nhit-1)*DIML+3]);
            //ppst.arg().set_eloss(
            //    parameters[0][DIMG+(cnt_nhit-1)*DIML+4]);
        }
        PhyArg curArg = ppst.arg();
        ppst.symbk();

        // Hit Status: Setting TOF reference time and path
        if (resetTOF && Hit<HitStTOF>::IsSame(hit)) { // set reference
            HitStTOF::SetOffsetPath(ppst.path());
            HitStTOF::SetOffsetTime(ppst.time()-Hit<HitStTOF>::Cast(hit)->t());
            resetTOF = false;
        }
        hit->cal(ppst);

        // Update Jacb
        if (hasJacb) {
            jbGG = curjb.gg() * jbGG;
            if (hasLoc) { // Local
                jbGL.at(cnt_nhit-1) = curjb.gl();
                for (UInt_t it = 0; it < cnt_nhit-1; ++it)
                    jbGL.at(it) = curjb.gg() * jbGL.at(it);
            } // Local
        }

        // Interaction
        if (hasLoc) { // Local
            SVecD<5> inrm, idiv;
            curArg.cal_nrm_and_div(inrm, idiv);

            for (UInt_t it = 0; it < DIML; ++it)
                rs(nseq_+(cnt_nhit-1)*DIML+it) += inrm(it);

            if (hasJacb) {
                for (UInt_t it = 0; it < DIML; ++it)
                    jb(nseq_+(cnt_nhit-1)*DIML+it, DIMG+(cnt_nhit-1)*DIML+it) += idiv(it);
            } // hasJacb
        } // Local
      
        // Coord
        if (hit->scx()) rs(hit->seqIDcx()) += hit->nrmcx();
        if (hit->scy()) rs(hit->seqIDcy()) += hit->nrmcy();
        if (hasJacb) {
            for (UInt_t it = 0; it < DIMG-is_fixed_mass_; ++it) {
                if (hit->scx()) jb(hit->seqIDcx(), it) += hit->divcx() * jbGG(0, it);
                if (hit->scy()) jb(hit->seqIDcy(), it) += hit->divcy() * jbGG(1, it);
            }
            if (hasLoc) {
                for (UInt_t it = 0; it < cnt_nhit; ++it) {
                    for (UInt_t jl = 0; jl < DIML; ++jl) {
                        if (hit->scx()) jb(hit->seqIDcx(), DIMG+it*DIML+jl) += hit->divcx() * jbGL.at(it)(0, jl);
                        if (hit->scy()) jb(hit->seqIDcy(), DIMG+it*DIML+jl) += hit->divcy() * jbGL.at(it)(1, jl);
                    }
                }
            } // Local
        } // hasJacb 

        // TRK
        HitStTRK* hitTRK = Hit<HitStTRK>::Cast(hit);
        if (hitTRK != nullptr) {
            if (hitTRK->sqx()) rs(hitTRK->seqIDqx()) += hitTRK->nrmqx();
            if (hitTRK->sqy()) rs(hitTRK->seqIDqy()) += hitTRK->nrmqy();
            if (hasJacb) {
                for (UInt_t it = 0; it < DIMG-is_fixed_mass_; ++it) {
                    if (hitTRK->sqx()) jb(hitTRK->seqIDqx(), it) += hitTRK->divqx() * jbGG(4, it);
                    if (hitTRK->sqy()) jb(hitTRK->seqIDqy(), it) += hitTRK->divqy() * jbGG(4, it);
                }
                if (hasLoc) {
                    for (UInt_t it = 0; it < cnt_nhit; ++it) {
                        for (UInt_t jl = 0; jl < DIML; ++jl) {
                            if (hitTRK->sqx()) jb(hitTRK->seqIDqx(), DIMG+it*DIML+jl) += hitTRK->divqx() * jbGL.at(it)(4, jl);
                            if (hitTRK->sqy()) jb(hitTRK->seqIDqy(), DIMG+it*DIML+jl) += hitTRK->divqy() * jbGL.at(it)(4, jl);
                        }
                    }
                } // Local
            } // hasJacb
        }
        
        // TOF
        HitStTOF* hitTOF = Hit<HitStTOF>::Cast(hit);
        if (hitTOF != nullptr) {
            if (hitTOF->sq()) rs(hitTOF->seqIDq()) += hitTOF->nrmq();
            if (hitTOF->st()) rs(hitTOF->seqIDt()) += hitTOF->nrmt();
            if (hasJacb) {
                for (UInt_t it = 0; it < DIMG-is_fixed_mass_; ++it) {
                    if (hitTOF->sq()) jb(hitTOF->seqIDq(), it) += hitTOF->divq() * jbGG(4, it);
                    if (hitTOF->st()) jb(hitTOF->seqIDt(), it) += hitTOF->divt() * jbGG(4, it);
                }
                if (hasLoc) {
                    for (UInt_t it = 0; it < cnt_nhit; ++it) {
                        for (UInt_t jl = 0; jl < DIML; ++jl) {
                            if (hitTOF->sq()) jb(hitTOF->seqIDq(), DIMG+it*DIML+jl) += hitTOF->divq() * jbGL.at(it)(4, jl);
                            if (hitTOF->st()) jb(hitTOF->seqIDt(), DIMG+it*DIML+jl) += hitTOF->divt() * jbGL.at(it)(4, jl);
                        }
                    }
                } // Local
            } // hasJacb
        }

        cnt_nhit++;
    }
    if (cnt_nhit != hits_.size()) return false;
    
    for (Int_t it = 0; it < numOfRes_; ++it) {
        if ( !Numc::Valid(rs(it)) ) rs(it) = Numc::ZERO<>;
        residuals[it] = rs(it);
    }
    
    if (hasJacb) {
        for (Int_t it = 0; it < numOfRes_; ++it) {
        for (Int_t jt = 0; jt < numOfPar_; ++jt) {
            if ( !Numc::Valid(jb(it, jt)) ) jb(it, jt) = Numc::ZERO<>;
            jacobians[0][it * numOfPar_ + jt] = jb(it, jt); 
        }}
    }

    return true;
}
    

PhySt PhyTrFit::interpolate_to_z(Double_t zcoo) const {
    PhySt nullst = part_; nullst.reset(part_.info().type());
    if (!succ_ || stts_.size() == 0) return nullst;

    // Find Index
    Int_t idx = -1;
    if (ortt_ == Orientation::kDownward)
        for (UInt_t it = 0; it < stts_.size(); ++it)
            if (idx == -1 && Numc::Compare(zcoo, (stts_.at(it).first)->cz()) <= 0) idx = it;
    else // kUpward
        for (UInt_t it = 0; it < stts_.size(); ++it)
            if (idx == -1 && Numc::Compare(zcoo, (stts_.at(it).first)->cz()) >= 0) idx = it;

    // Set Particle
    PhySt ppst = ((idx == -1) ? (stts_.at(0).second) : (stts_.at(idx).second));
    if (idx == -1) ppst.arg().clear();
    PhyArg arg = ppst.arg();
        
    if (!PropMgnt::PropToZ(zcoo, ppst)) return nullst;
    ppst.symbk();

    ppst.arg().clear();
    ppst.arg().set_mscat(arg.tauu(), arg.rhou(), arg.taul(), arg.rhol());
    ppst.arg().set_eloss(arg.elion(), arg.elbrm());

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
            if (isat == -1 && Numc::Compare(zsat, (stts_.at(it).first)->cz()) <= 0) isat = it;
            if (iend == -1 && Numc::Compare(zend, (stts_.at(it).first)->cz()) <= 0) iend = it;
        }
    }
    else { // kUpward
        for (UInt_t it = 0; it < stts_.size(); ++it) {
            if (isat == -1 && Numc::Compare(zsat, (stts_.at(it).first)->cz()) >= 0) isat = it;
            if (iend == -1 && Numc::Compare(zend, (stts_.at(it).first)->cz()) >= 0) iend = it;
        }
    }

    // Material
    std::list<MatFld> mflds;
    for (Int_t it = isat; it <= iend; ++it) {
        PhySt ppst = (it == isat) ? interpolate_to_z(zsat) : (stts_.at(it).second);
        if (Numc::EqualToZero(ppst.mom())) return MatFld();

        MatFld fld;
        Double_t zcoo = (it == iend) ? zend : (stts_.at(it+1).first)->cz();
        if (!PropMgnt::PropToZ(zcoo, ppst, &fld)) return MatFld();
        mflds.push_back(fld);
    }
    
    MatFld&& mfld = MatFld::Merge(mflds);
    return mfld;
}


//PhyMassFit::PhyMassFit(const TrFitPar& fitPar, Short_t chrg) : TrFitPar(fitPar), phyTr_(nullptr) {
//    if (chrg <= Numc::ZERO<Short_t>) return;
//    if (chrg >= PartListMassQ.size()) return;
//    is_fixed_mass_ = false;
//
//    Double_t nchi = 0;
//    Bool_t   succ = false;
//    for (auto&& mass : PartListMassQ.at(chrg)) {
//        CERR("MASS %14.8f\n", mass);
//        info_.reset(chrg, mass);
//        PhyTrFit* phyTr = new PhyTrFit(static_cast<const TrFitPar&>(*this));
//        if (!phyTr->status()) { delete phyTr; continue; }
//        if (succ && phyTr->nchi() > nchi) { delete phyTr; continue; }
//        phyTr_ = phyTr;
//        nchi = phyTr->nchi();
//        succ = true;
//    }
//    if (!succ) clear();
//}


//bool VirtualPhyMassFit::operator() (const double* const x, double* residuals) const {
//    residuals[0] = Numc::ZERO<>;
//    if (!is_vary_mass()) return false;
//    if (Numc::Compare(x[0]) <= 0) return false;
//    
//    PartInfo::SetDefaultInvu(x[0]);
//
//    PhyTrFit trFit(fitPar_);
//    if (!trFit.status()) return false;
//
//    double nrm = std::sqrt(trFit.nchi());
//    if (!Numc::Valid(nrm) || Numc::Compare(nrm) < 0) return false;
//
//    residuals[0] = nrm;
//    return true;
//}
//
//
//PhyMassFit::PhyMassFit(const TrFitPar& fitPar, Short_t chrg, Double_t mass) : phyTr_(nullptr) {
//    if (fitPar.info().type() != PartType::Fixed) return;
//
//    std::vector<double> parameters{ PartInfo::ATOMIC_MASS/mass };
//
//    ceres::CostFunction* cost_function = new ceres::NumericDiffCostFunction<VirtualPhyMassFit, ceres::CENTRAL, 1, 1>(new VirtualPhyMassFit(fitPar, chrg));
//
//    ceres::Problem problem;
//    problem.AddResidualBlock(cost_function, nullptr, parameters.data());
//    
//    Double_t lmtl_invu = (1.0/56.0) * 0.1; // Fe
//    Double_t lmtu_invu = 2000.0 * 10.0; // El
//    problem.SetParameterLowerBound(parameters.data(), 0, lmtl_invu);
//    problem.SetParameterUpperBound(parameters.data(), 0, lmtu_invu);
//    
//    ceres::Solver::Options options;
//    //options.minimizer_type = ceres::LINE_SEARCH;
//    //options.line_search_direction_type = ceres::BFGS;
//
//    ceres::Solver::Summary summary;
//    ceres::Solve(options, &problem, &summary);
//    if (!summary.IsSolutionUsable()) return;
//
//    PartInfo::SetDefaultInvu(parameters.at(0));
//    phyTr_ = new PhyTrFit(fitPar);
//    if (!phyTr_->status()) clear();
//}


} // namespace TrackSys


#endif // __TRACKLibs_PhyFit_C__
