#ifndef __TRACKLibs_PhyTr_C__
#define __TRACKLibs_PhyTr_C__


namespace TrackSys {


PhyTr::PhyTr(const std::vector<HitSt>& hits, const PartType& type, const Orientation& ortt) : PhyTr() {
    if (!check(hits)) return;
    
    type_ = type;
    ortt_ = ortt;
    hits_ = hits;
    for (auto&& hit : hits) {
        if (hit.sx()) nhit_(0)++;
        if (hit.sy()) nhit_(1)++;
    }
    if (ortt == Orientation::kDownward) HitSt::Sort(hits_, HitSt::Orientation::kDownward);
    else                                HitSt::Sort(hits_, HitSt::Orientation::kUpward);
    part_.reset(type);

    if (!fit()) clear();
    else        succ_ = true;
}
        

void PhyTr::clear() {
    succ_ = false;
    ortt_ = Orientation::kDownward;
    nhit_ = SVecI<2>();
    hits_.clear();
    type_ = PartType::Proton;
    part_.reset(type_);
    phys_.clear();
    marg_.clear();
    nchi_ = 0.;
    ndf_ = 0;
}


Bool_t PhyTr::check(const std::vector<HitSt>& hits) {
    Short_t nhit_x = 0;
    Short_t nhit_y = 0;
    for (auto&& hit : hits) {
        if (hit.sx()) nhit_x++;
        if (hit.sy()) nhit_y++;
    }

    return (nhit_x >= LMTL_NHIT_X && nhit_y >= LMTL_NHIT_Y);
}


Bool_t PhyTr::fit() {
    //MagMgnt::Load();
    //MatGeoBoxAms::Load();

    //const Int_t ntimes = 1000;

    COUT("\n==== Analysis ====\n");
    MGClock::HrsStopwatch swa;
    //for (Int_t it = 0; it < ntimes; ++it)
    if (!fit_analysis()) return false;
    swa.stop();
    swa.print();
    part_.print();

    //nchi_ = 0.;
    //ndf_ = 0;
    COUT("\n==== Simple ====\n");
    MGClock::HrsStopwatch sws;
    //for (Int_t it = 0; it < ntimes; ++it)
    if (!fit_simple()  ) return false;
    sws.stop();
    sws.print();
    part_.print();
    COUT("CHI %14.8f\n", nchi_);   
    
    COUT("\n==== Semi-Simple ====\n");
    MGClock::HrsStopwatch swss;
    //for (Int_t it = 0; it < ntimes; ++it)
    if (!fit_semi_simple()  ) return false;
    //if (!fit_simple()  ) return false;
    swss.stop();
    swss.print();
    part_.print();
    COUT("CHI %14.8f\n", nchi_);   
   
    //nchi_ = 0.;
    //ndf_ = 0;
    //COUT("\n==== Physics ====\n");
    //MGClock::HrsStopwatch swp;
    ////for (Int_t it = 0; it < ntimes; ++it)
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

        const Int_t nstp = 3;
        for (Int_t ih = 1; ih < hits_.size(); ++ih) {
            SVecD<3>&& ref_l = (hits_.at(ih).coo() - hits_.at(ih-1).coo());
            SVecD<3>&& ref_u = LA::Unit(ref_l);
            Double_t   ref_s = LA::Mag(ref_l);
            
            SVecD<3> mfldv;
            for (Int_t it = 0; it < nstp; ++it) {
                Double_t stp = ((static_cast<Double_t>(it) + MGMath::HALF) / static_cast<Double_t>(nstp));
                SVecD<3>&& ref_m = ((MGMath::ONE - stp) * hits_.at(ih-1).coo() + stp * hits_.at(ih).coo());
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
    Bool_t preSucc = false;
    Bool_t curSucc = false;
    Int_t  curIter = 1;
    while (curIter <= LMTU_ITER && !succ) {
        Int_t ndfx = 0;
        Int_t ndfy = 0;
        Double_t chix = 0.;
        Double_t chiy = 0.;
        SVecD<5>    gradG;
        SMtxSymD<5> covGG;

        Int_t cnt_nhit = 0;
        PhySt ppst(part_);
        PhyJb ppjb(PhyJb::Type::kIdentity);
        for (auto&& hit : hits_) {
            PhyJb curjb;
            if (!PropMgnt::PropToZ(hit.cz(), ppst, MatArg(), &curjb)) break;
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
            
            cnt_nhit++;
            if (!hit.sx()) hit.set_dummy_x(ppst.cx());
        }
        if (cnt_nhit != hits_.size()) return false;
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
        
        Int_t    ndf  = ndfx + ndfy;
        Double_t nchi = ((chix + chiy) / static_cast<Double_t>(ndf));
        Double_t nchi_rat = std::fabs((nchi - nchi_) / (nchi + nchi_ + CONVG_EPSILON));
        Bool_t   sign     = (MGNumc::Compare(nchi - nchi_, CONVG_EPSILON) < 0);
        
        ndf_ = ndf;
        nchi_ = nchi;

        curSucc = (curIter >= LMTL_ITER && (sign && MGNumc::Compare(nchi_rat, CONVG_TOLERANCE) < 0));
        
        succ = (preSucc && curSucc);
        preSucc = curSucc;
        
        curIter++;
    }
    
    return succ;
}


Bool_t PhyTr::fit_semi_simple() {
    std::vector<MatFld> mflds(hits_.size());
    {
        PhySt ppst(part_);
        Int_t cnt_nhit = 0;
        for (auto&& hit : hits_) {
            if (!PropMgnt::PropToZ(hit.cz(), ppst, MatArg(), nullptr, &mflds.at(cnt_nhit))) break;
            cnt_nhit++;
        }
        if (cnt_nhit != hits_.size()) return false;
    }
    
    Bool_t succ = false;
    Bool_t preSucc = false;
    Bool_t curSucc = false;
    Int_t  curIter = 1;
    while (curIter <= LMTU_ITER && !succ) {
        Int_t ndfx = 0;
        Int_t ndfy = 0;
        Double_t chix = 0.;
        Double_t chiy = 0.;
        SVecD<5>    gradG;
        SMtxSymD<5> covGG;

        Int_t cnt_nhit = 0;
        PhySt       ppst(part_);
        SMtxD<5>    ppjb = SMtxId5D;
        SMtxSymD<5> pper;
        for (auto&& hit : hits_) {
            MatMscatFld&& mscat = MatPhy::GetMscat(mflds.at(cnt_nhit), ppst);
            
            PhyJb curjb;
            if (!PropMgnt::PropToZ(hit.cz(), ppst, MatArg(), &curjb)) break;
            PhyJb::SMtxDXYG&& subjb = curjb.xyg() * ppjb;
            
            SVecD<2> mres(ppst.cx() - hit.cx(), ppst.cy() - hit.cy());
            SVecD<2>&& merr = hit.err(mres);
          
            SMtxSymD<5> pder = (LA::Similarity(curjb.gg(), pper) + mscat.cov());
            SMtxSymD<2>&& udS = pder.Sub<SMtxSymD<2>>(0, 0);
            udS(0, 0) += merr(0) * merr(0);
            udS(1, 1) += merr(1) * merr(1);
            
            Int_t udSfail = 0;
            SMtxSymD<2>&& iudS = udS.Inverse(udSfail);
            if (udSfail != 0) break;
            
            gradG += LA::Transpose(subjb) * iudS * mres;
            covGG += LA::SimilarityT(subjb, iudS);
            
            SMtxD<5, 2>&& udK   = pder.Sub<SMtxD<5, 2>>(0, 0) * iudS;
            SVecD<5>&&    rslst = udK * mres;
            SMtxSymD<5>&& rsler = LA::Similarity(udK, udS);
            
            pper = std::move(pder - rsler);
            ppst.set_state_with_uxy(
                ppst.cx() - rslst(0),
                ppst.cy() - rslst(1),
                ppst.cz(),
                ppst.dx() - rslst(2),
                ppst.dy() - rslst(3),
                ((ortt_ == Orientation::kDownward) ? -1 : 1)
            );
            ppst.set_eta(ppst.eta() - rslst(4));

            SMtxSymD<5> udjb(SMtxIdSym5D);
            udjb(0, 0) += udK(0, 0);
            udjb(0, 1) += udK(0, 1);
            udjb(1, 1) += udK(1, 1);
            ppjb = std::move(udjb * curjb.gg() * ppjb); 

            if (hit.sx()) { ndfx++; chix += iudS(0, 0) * mres(0) * mres(0); } 
            if (hit.sy()) { ndfy++; chiy += iudS(1, 1) * mres(1) * mres(1); } 
            
            cnt_nhit++;
            if (!hit.sx()) hit.set_dummy_x(ppst.cx());
        }
        if (cnt_nhit != hits_.size()) return false;
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
        
        Int_t    ndf  = ndfx + ndfy;
        Double_t nchi = ((chix + chiy) / static_cast<Double_t>(ndf));
        Double_t nchi_rat = std::fabs((nchi - nchi_) / (nchi + nchi_ + CONVG_EPSILON));
        Bool_t   sign     = (MGNumc::Compare(nchi - nchi_, CONVG_EPSILON) < 0);

        ndf_ = ndf;
        nchi_ = nchi;

        std::cout << Form("IT %d CHI %14.8f\n", curIter, nchi);

        curSucc = (curIter >= LMTL_ITER && (sign && MGNumc::Compare(nchi_rat, CONVG_TOLERANCE) < 0));

        succ = (preSucc && curSucc);
        preSucc = curSucc;

        curIter++;
    }
    
    return succ;
}


Bool_t PhyTr::fit_physics() {
    return true;
/*    
    ////////////////////////////////////////////////////////////////////////////////////////////
    Bool_t succ = false;
    Bool_t preSucc = false;
    Bool_t curSucc = false;
    Int_t  curIter = 1;
    while (curIter <= LMTU_ITER && !succ) {
        Int_t ndfx = 0;
        Int_t ndfy = 0;
        Double_t chix = 0.;
        Double_t chiy = 0.;
        SVecD<5>    gradG;
        SMtxSymD<5> covGG;

        Int_t cnt_nhit = 0;
        PhySt       ppst(part_);
        SMtxD<5>    ppjb = SMtxId5D;
        SMtxSymD<5> pper;
        for (auto&& hit : hits_) {
            MatMscatFld&& mscat = MatPhy::GetMscat(mflds.at(cnt_nhit), ppst);
            
            PhyJb curjb;
            if (!PropMgnt::PropToZ(hit.cz(), ppst, MatArg(), &curjb)) break;
            PhyJb::SMtxDXYG&& subjb = curjb.xyg() * ppjb;
            
            SVecD<2> mres(ppst.cx() - hit.cx(), ppst.cy() - hit.cy());
            SVecD<2>&& merr = hit.err(mres);
          
            SMtxSymD<5> pder = (LA::Similarity(curjb.gg(), pper) + mscat.cov());
            SMtxSymD<2>&& udS = pder.Sub<SMtxSymD<2>>(0, 0);
            udS(0, 0) += merr(0) * merr(0);
            udS(1, 1) += merr(1) * merr(1);
            
            Int_t udSfail = 0;
            SMtxSymD<2>&& iudS = udS.Inverse(udSfail);
            if (udSfail != 0) break;
            
            gradG += LA::Transpose(subjb) * iudS * mres;
            covGG += LA::SimilarityT(subjb, iudS);
            
            SMtxD<5, 2>&& udK   = pder.Sub<SMtxD<5, 2>>(0, 0) * iudS;
            SVecD<5>&&    rslst = udK * mres;
            SMtxSymD<5>&& rsler = LA::Similarity(udK, udS);
            
            pper = std::move(pder - rsler);
            ppst.set_state_with_uxy(
                ppst.cx() - rslst(0),
                ppst.cy() - rslst(1),
                ppst.cz(),
                ppst.dx() - rslst(2),
                ppst.dy() - rslst(3),
                ((ortt_ == Orientation::kDownward) ? -1 : 1)
            );
            ppst.set_eta(ppst.eta() - rslst(4));

            SMtxSymD<5> udjb(SMtxIdSym5D);
    /////////////////////////////////////////////////////////////////////////////////////
    */
    
    return true;
}


} // namespace TrackSys


#endif // __TRACKLibs_PhyTr_C__
