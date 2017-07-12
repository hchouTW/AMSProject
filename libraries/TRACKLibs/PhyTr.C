#ifndef __TRACKLibs_PhyTr_C__
#define __TRACKLibs_PhyTr_C__


namespace TrackSys {


PhyTr::PhyTr(const std::vector<HitSt>& hits, const PartType& type, const Orientation& ortt) : PhyTr() {
    if (!check(hits)) return;
    
    type_ = type;
    ortt_ = ortt;
    hits_ = hits;
    if (ortt == Orientation::kDownward) HitSt::Sort(hits_, HitSt::Orientation::kDownward);
    else                                HitSt::Sort(hits_, HitSt::Orientation::kUpward);
    part_.reset(type);

    if (!fit()) clear();
    else        succ_ = true;
}
        

void PhyTr::clear() {
    succ_ = false;
    ortt_ = Orientation::kDownward;
    hits_.clear();
    type_ = PartType::Proton;
    part_.reset(type_);
    phys_.clear();
    marg_.clear();
    nchi_ = 0.;
    ndf_ = 0;
}


Bool_t PhyTr::check(const std::vector<HitSt>& hits) {
    Short_t num_hit_x = 0;
    Short_t num_hit_y = 0;
    for (auto&& hit : hits) {
        num_hit_x += hit.sx();
        num_hit_y += hit.sy();
    }

    return (num_hit_x >= LMTL_HIT_X && num_hit_y >= LMTL_HIT_Y);
}


Bool_t PhyTr::fit() {
    MagMgnt::Load();
    MatGeoBoxAms::Load();

    //const Int_t ntimes = 1000;

    //COUT("\n==== Analysis ====\n");
    //MGClock::HrsStopwatch swa;
    //for (Int_t it = 0; it < ntimes; ++it)
    if (!fit_analysis()) return false;
    //swa.stop();
    //swa.print();
    //part_.print();

    //nchi_ = 0.;
    //ndf_ = 0;
    //COUT("\n==== Simple ====\n");
    //MGClock::HrsStopwatch sws;
    //for (Int_t it = 0; it < ntimes; ++it)
    if (!fit_simple()  ) return false;
    //sws.stop();
    //sws.print();
    //part_.print();
   
    //nchi_ = 0.;
    //ndf_ = 0;
    //COUT("\n==== Physics ====\n");
    //MGClock::HrsStopwatch swp;
    //for (Int_t it = 0; it < ntimes; ++it)
    //if (!fit_physics() ) return false;
    //swp.stop();
    //swp.print();
    //part_.print();
    
    return true; 
}


Bool_t PhyTr::fit_analysis() {
    // Linear Fit on X
    // Equation of Motion
    // X  = PX + UX * (Z - RefZ)
    // dZ = Z - RefZ
    // Chisq = (X - Xm)^2
    // | PX |   | Sum(1)     Sum(dZ)   |^-1   | Sum(Xm)    |
    // |    | = |                      |    * |            |
    // | UX |   | Sum(dZ)    Sum(dZ^2) |      | Sum(dZ*Xm) |
    Double_t prefit_pz = hits_.at(0).cz();
    Double_t prefit_px = MGMath::ZERO;
    Double_t prefit_ux = MGMath::ZERO;
    {
        SMtxSymD<2> mtx;
        SVecD<2>    res;
        for (auto&& hit : hits_) {
            if (!hit.sx()) continue;
            Double_t ex  = hit.ex(0.0);
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
        prefit_ux = rsl(1);
        
        for (auto&& hit : hits_) {
            if (hit.sx()) continue;
            Double_t dz = hit.cz() - prefit_pz;
            Double_t px = prefit_px + prefit_ux * dz;
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
        Double_t Lambda = PROP_FACT * part_.part().chrg_to_mass(); 

        std::vector<Double_t> stp(hits_.size(), MGMath::ZERO);
        std::vector<Double_t> crs(hits_.size(), MGMath::ZERO);

        for (Int_t ih = 1; ih < hits_.size(); ++ih) {
            SVecD<3>&& ref_l = (hits_.at(ih).coo() - hits_.at(ih-1).coo());
            SVecD<3>&& ref_u = LA::Unit(ref_l);
            Double_t   ref_s = LA::Mag(ref_l);
            SVecD<3>&& ref_m = (MGMath::HALF * (hits_.at(ih-1).coo() + hits_.at(ih).coo()));
            MagFld&&   mfld  = MagMgnt::Get(ref_m);
            Double_t   mucrs = Lambda * (ref_u(2) * mfld.x() - ref_u(0) * mfld.z());
            stp.at(ih) = ref_s;
            crs.at(ih) = mucrs;
        }
        
        SMtxSymD<3> mtx;
        SVecD<3>    res;
        Double_t    cur_Au = MGMath::ZERO;
        Double_t    cur_Ae = MGMath::ZERO;
        for (Int_t ih = 0; ih < hits_.size(); ++ih) {
            HitSt& hit = hits_.at(ih);
            Double_t ey  = hit.ey(0.0);
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
    Double_t prefit_uz = (MGMath::ONE - prefit_ux*prefit_ux - prefit_uy*prefit_uy);
    if (MGNumc::Compare(prefit_uz) <= 0) prefit_uz = MGMath::ZERO;
    else                                 prefit_uz = ((ortt_ == Orientation::kDownward) ? MGMath::NEG : MGMath::ONE) * std::sqrt(prefit_uz);

    part_.set_state_with_cos(
        prefit_px, prefit_py, prefit_pz,
        prefit_ux, prefit_uy, prefit_uz
    );
    part_.set_eta(prefit_ea);

    return true;
}


Bool_t PhyTr::fit_simple() {
    Bool_t succ = false;
    Int_t  curIter = 1;
    std::vector<MatFld> mflds;
    while (curIter <= LMTU_ITER && !succ) {
        Int_t ndfx = 0;
        Int_t ndfy = 0;
        Double_t chix = 0.;
        Double_t chiy = 0.;
        SVecD<5>    gradG;
        SMtxSymD<5> covGG;

        Int_t cnt_hits = 0;
        PhySt ppst(part_);
        PhyJb ppjb(PhyJb::Type::kIdentity);
        std::vector<Double_t> path_len;
        for (auto&& hit : hits_) {
            PhyJb curjb;
            SVecD<3> init_coo = ppst.coo();
            if (!PropMgnt::PropToZ(hit.cz(), ppst, MatArg(), &curjb)) break;
            //if (curIter == LMTL_ITER) mflds.push_back(MatMgnt::Get(init_coo, ppst.coo())); 
            path_len.push_back(curjb.path_len());
            ppjb.multiplied(curjb);
            
            SVecD<2> mres(ppst.cx() - hit.cx(), ppst.cy() - hit.cy());
            SVecD<2>&& merr = hit.err(mres);
            
            SMtxSymD<2> micov;
            micov(0, 0) = (hit.sx() ? (MGMath::ONE / merr(0) / merr(0)) : MGMath::ZERO);
            micov(1, 1) = (hit.sy() ? (MGMath::ONE / merr(1) / merr(1)) : MGMath::ZERO);

            SVecD<2> mgrad;
            mgrad(0) = (hit.sx() ? (micov(0, 0) * mres(0)) : MGMath::ZERO);
            mgrad(1) = (hit.sy() ? (micov(1, 1) * mres(1)) : MGMath::ZERO);
            
            PhyJb::SMtxDXYG&& subjb = ppjb.xyg();
            gradG += LA::Transpose(subjb) * mgrad;
            covGG += LA::SimilarityT(subjb, micov);

            if (hit.sx()) { ndfx++; chix += (mgrad(0) * mres(0)); } 
            if (hit.sy()) { ndfy++; chiy += (mgrad(1) * mres(1)); } 
            
            cnt_hits++;
            if (!hit.sx()) hit.set_dummy_x(ppst.cx());
        }
        if (cnt_hits != hits_.size()) return false;
        if (!covGG.Invert()) return false;
        SVecD<5>&& rslG = covGG * gradG;

        part_.set_state_with_uxy(
            part_.cx() - rslG(0),
            part_.cy() - rslG(1),
            part_.cz(),
            part_.dx() - rslG(2),
            part_.dy() - rslG(3),
            ((ortt_ == Orientation::kDownward) ? -1 : 1)
        );
        part_.set_eta(part_.eta() - rslG(4));
        
        Int_t ndf = ndfx + ndfy;
        Double_t nchi = ((chix + chiy) / static_cast<Double_t>(ndf));
        Double_t nchi_rat = std::fabs((nchi_ - nchi) / (nchi_ + nchi));
        
        ndf_ = ndf;
        nchi_ = nchi;

        succ = (curIter >= LMTL_ITER && MGNumc::Compare(nchi_rat, CONVG_TOLERANCE) < 0);
        curIter++;
    }
    
    return succ;
}


Bool_t PhyTr::fit_physics() {
    /*
    nchi_ = 0.;
    ndf_ = 0;
    marg_.resize(hits_.size(), MatArg(false, false));

    Bool_t succ = false;
    Int_t  curIter = 1;
    while (curIter <= LMTU_ITER && !succ) {
        Int_t ndfx = 0;
        Int_t ndfy = 0;
        Double_t chix = 0;
        Double_t chiy = 0;
        SVecD<5>    gradG;
        SMtxSymD<5> covGG;

        Int_t cnt_hits = 0;
        PhySt ppst(part_);
        PhyJb ppjb(PhyJb::Type::kIdentity);
        for (auto&& hit : hits_) {
            PhyJb curjb;
            if (!PropMgnt::PropToZ(hit.cz(), ppst, marg_.at(cnt_hits), &curjb)) break;
            ppjb.multiplied(curjb);
            PhyJb::SMtxDXYG&& subjb = ppjb.xyg();
            
            SVecD<2> mdif(ppst.cx() - hit.cx(), ppst.cy() - hit.cy());

            SMtxSymD<2> merr;
            merr(0, 0) = (hit.sx() ? (MGMath::ONE / hit.ex() / hit.ex()) : MGMath::ZERO);
            merr(1, 1) = (hit.sy() ? (MGMath::ONE / hit.ey() / hit.ey()) : MGMath::ZERO);

            SVecD<2> mres;
            mres(0) = (hit.sx() ? (mdif(0) * merr(0, 0)) : MGMath::ZERO);
            mres(1) = (hit.sy() ? (mdif(1) * merr(1, 1)) : MGMath::ZERO);
            
            gradG += LA::Transpose(subjb) * mres;
            covGG += LA::SimilarityT(subjb, merr);

            if (hit.sx()) { ndfx++; chix += (mdif(0) * mdif(0) * merr(0, 0)); } 
            if (hit.sy()) { ndfy++; chiy += (mdif(1) * mdif(1) * merr(1, 1)); } 
            
            cnt_hits++;
        }
        Int_t ndf = ndfx + ndfy;
        Double_t nchi = ((chix + chiy) / static_cast<Double_t>(ndf));
        Double_t nchi_rat = std::fabs((nchi_ - nchi) / (nchi_ + nchi));

        if (cnt_hits != hits_.size()) return false;
        if (!covGG.Invert()) return false;
        SVecD<5>&& rslG = covGG * gradG;

        part_.set_state_with_uxy(
            part_.cx() - rslG(0),
            part_.cy() - rslG(1),
            part_.cz(),
            part_.dx() - rslG(2),
            part_.dy() - rslG(3),
            ((ortt_ == Orientation::kDownward) ? -1 : 1)
        );
        part_.set_eta(part_.eta() - rslG(4));

        succ = (curIter >= LMTL_ITER && MGNumc::Compare(nchi_rat, CONVG_TOLERANCE) < 0);

        ndf_ = ndf;
        nchi_ = nchi;
        curIter++;
    }
    return succ;
    */
    return true;
}


} // namespace TrackSys


#endif // __TRACKLibs_PhyTr_C__
