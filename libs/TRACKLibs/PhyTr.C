#ifndef __TRACKLibs_PhyTr_C__
#define __TRACKLibs_PhyTr_C__


namespace TrackSys {
        
    
Bool_t PhyTr::HitCheck(const std::vector<HitSt>& hits) {
    // Number of Hit Requirement
    const Int_t LMTL_NHIT_X = 3;
    const Int_t LMTL_NHIT_Y = 4;
    Short_t nx = 0, ny = 0;
    for (auto&& hit : hits) {
        if (hit.sx()) nx++;
        if (hit.sy()) ny++;
    }

    return (nx >= LMTL_NHIT_X && ny >= LMTL_NHIT_Y);
}


std::tuple<Short_t, std::vector<Short_t>, std::vector<Short_t>, std::vector<Short_t>, Short_t, Short_t> PhyTr::HitSort(std::vector<HitSt>& hits, Orientation ortt) {
    if (ortt == Orientation::kDownward) HitSt::Sort(hits, HitSt::Orientation::kDownward);
    else                                HitSt::Sort(hits, HitSt::Orientation::kUpward);
    
    Short_t hID = 0;
    Short_t nseq = 0;
    std::vector<Short_t> seqx;
    std::vector<Short_t> seqy;
    std::vector<Short_t> maps;
    Short_t              nhtx;
    Short_t              nhty;
    for (auto&& hit : hits) {
        hit.set_seqID(nseq);
        if (hit.sx()) { seqx.push_back(nseq); maps.push_back(0 * hits.size() + hID); } 
        else          { seqx.push_back(-1); }
        if (hit.sy()) {
            if (hit.sx()) seqy.push_back(nseq+1);
            else          seqy.push_back(nseq);
            maps.push_back(1 * hits.size() + hID);
        }
        else seqy.push_back(-1);

        hID++;
        nseq += (hit.sx() + hit.sy());
        nhtx += hit.sx();
        nhty += hit.sy();
    }
    return std::make_tuple(nseq, seqx, seqy, maps, nhtx, nhty);
}


void PhyTr::clear() {
    sw_mscat_ = false;
    sw_eloss_ = false;
    
    type_ = PartType::Proton;
    ortt_ = Orientation::kDownward;
    nseq_ = 0;
    seqx_.clear();
    seqy_.clear();
    maps_.clear();
    hits_.clear();
    nhtx_ = 0;
    nhty_ = 0;

    succ_ = false;
    part_.reset(type_);
    ndfx_ = 0;
    ndfy_ = 0;
    chix_ = 0.;
    chiy_ = 0.;
    nchi_ = 0.;
}
        

void PhyTr::print() const {
    std::cerr << "\n======================================\n";
    part_.print();
    for (auto&& hit : hits_) hit.print();
    std::cerr << "======================================\n";
}


PhyTr::PhyTr(const std::vector<HitSt>& hits, const PartType& type, const Orientation& ortt, Bool_t sw_mscat, Bool_t sw_eloss) {
    clear();
    if (!HitCheck(hits)) return;
    sw_mscat_ = sw_mscat;
    sw_eloss_ = sw_eloss;

    type_ = type;
    ortt_ = ortt;
    hits_ = hits;
    auto&& sort = HitSort(hits_);
    nseq_ = std::move(std::get<0>(sort)); 
    seqx_ = std::move(std::get<1>(sort)); 
    seqy_ = std::move(std::get<2>(sort)); 
    maps_ = std::move(std::get<3>(sort)); 
    nhtx_ = std::move(std::get<4>(sort)); 
    nhty_ = std::move(std::get<5>(sort)); 
    
    part_.reset(type);
}


Bool_t PhyTr::fit() {
    if (!fit_analysis()) return false;
    if (!fit_simple()) return false;

    if (sw_mscat_ || sw_eloss_) {
        //if (!fit_physics()) return false;
    }

    //const Int_t ntimes = 1000;
    //MGClock::HrsStopwatch swa;
    //for (Int_t it = 0; it < ntimes; ++it)
    //if (!fit_semi_simple()) return false;
    //swa.stop();
    //swa.print();
    //part_.print();
    
    return true; 
}


Bool_t PhyTr::fit_analysis() {
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


Bool_t PhyTr::fit_simple() {
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
            Double_t msl = (ppst.arg().mscat() ? ppst.arg().mscat_ll() : MGMath::ZERO);
            if (ppst.arg().mscat() && hit.sx()) mex = std::hypot(mex, msl);
            if (ppst.arg().mscat() && hit.sy()) mey = std::hypot(mey, msl);

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
        Double_t nchi = ((chi) / static_cast<Double_t>(nseq_ - 5));

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
                ndfx_    = (nhtx_ - 2);
                ndfy_    = (nhty_ - 3);
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
        
        curIter++;
    }
    //if (!succ) std::cout << Form("FAIL. IT %2d %2d (RIG %14.8f CHI %14.8f) LAMBDA %14.8f\n", curIter, updIter, part_.rig(), nchi_, lambda);
    //else       std::cout << Form("SUCC. IT %2d %2d (RIG %14.8f CHI %14.8f) LAMBDA %14.8f\n", curIter, updIter, part_.rig(), nchi_, lambda);
    
    return succ;
}




























Bool_t PhyTr::fit_physics() {
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
        Double_t chi = MGMath::ZERO;
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
            Double_t msl = (ppst.arg().mscat() ? ppst.arg().mscat_ll() : MGMath::ZERO);
            if (ppst.arg().mscat() && hit.sx()) mex = std::hypot(mex, msl);
            if (ppst.arg().mscat() && hit.sy()) mey = std::hypot(mey, msl);

            SMtxSymD<2> cvM;
            cvM(0, 0) = (hit.sx() ? (MGMath::ONE / mex / mex) : MGMath::ZERO);
            cvM(1, 1) = (hit.sy() ? (MGMath::ONE / mey / mey) : MGMath::ZERO);
            
            SVecD<2> rsM;
            rsM(0) = (hit.sx() ? cvM(0, 0) * (hit.cx() - ppst.cx()) : MGMath::ZERO);
            rsM(1) = (hit.sy() ? cvM(1, 1) * (hit.cy() - ppst.cy()) : MGMath::ZERO);
            
            PhyJb::SMtxDXYG&& subJbF = PhyJb::SubXYG(ppjb);
            grdG += LA::Transpose(subJbF) * rsM;
            cvGG += LA::SimilarityT(subJbF, cvM);

            if (hit.sx()) { chi += rsM(0) * (hit.cx() - ppst.cx()); }
            if (hit.sy()) { chi += rsM(1) * (hit.cy() - ppst.cy()); }

            cnt_nhit++;
            if (!hit.sx()) hit.set_dummy_x(ppst.cx());
        }
        if (cnt_nhit != hits_.size()) break;
        Double_t nchi = ((chi) / static_cast<Double_t>(nseq_ - 5));

        Bool_t isSucc   = false;
        Bool_t isUpdate = false;
        if (curIter != 0) {
            Double_t lmRho = (nchi_ - nchi) / curLmRhoDen;
            Bool_t   isLmt = (MGNumc::Compare(lambda, LMTU_LAMBDA) >= 0);
            Double_t convg = std::cbrt((MGMath::ONE + lambda) * (MGMath::ONE + lambda));
            isSucc = (MGNumc::Compare(((nchi_ - nchi) / (nchi_ + nchi + CONVG_TOLERANCE)) * convg, CONVG_TOLERANCE) <= 0);

            if (MGNumc::Compare(lmRho, CONVG_EPSILON) < 0) {
                lambda = std::min(lambda*LAMBDA_UP_FAC, LMTU_LAMBDA); 
                grdG   = curGrdG;
                cvGG   = curCvGG;
                rltSt  = part_;
                if (isSucc && isLmt) updIter++;
            }
            else {
                lambda   = std::max(lambda/LAMBDA_DN_FAC, LMTL_LAMBDA);
                ndfx_    = (nseq_ - 5);
                ndfy_    = (nseq_ - 5);
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
        
        curIter++;
    }
    //if (!succ) std::cout << Form("FAIL. IT %2d %2d (RIG %14.8f CHI %14.8f) LAMBDA %14.8f\n", curIter, updIter, part_.rig(), nchi_, lambda);
    //else       std::cout << Form("SUCC. IT %2d %2d (RIG %14.8f CHI %14.8f) LAMBDA %14.8f\n", curIter, updIter, part_.rig(), nchi_, lambda);
    
    return succ;
}


































//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
Bool_t PhyTr::fit_physics() {
    Bool_t succ = false;
    Bool_t preSucc = false;
    Bool_t curSucc = false;
    
    Int_t  curIter = 1;
    while (curIter <= LMTU_ITER && !succ) {
        Double_t chi = MGMath::ZERO;
        Int_t    ndf = (nseq_ - 5);
       
        TMtxD JbF(nseq_, DIM_G);
        TMtxD CvM(nseq_, nseq_);
        TVecD RsM(nseq_);

        Int_t nseq = 0;
        Int_t cnt_nhit = 0;
        PhySt ppst(part_);
        PhyJb::SMtxDGG&& ppjbF = SMtxId();
        Double_t         ppMsl = MGMath::ZERO;
        std::vector<Double_t> vMex(hits_.size());
        std::vector<Double_t> vMey(hits_.size());
        std::vector<Double_t> vMsl(hits_.size());
        for (auto&& hit : hits_) {
            PhyJb curjb;
            if (!PropMgnt::PropToZ(hit.cz(), ppst, nullptr, &curjb)) break;
            ppjbF = curjb.gg() * ppjbF;
            
            Double_t mex  = (hit.sx() ? hit.ex(hit.cx() - ppst.cx()) : MGMath::ZERO);
            Double_t mey  = (hit.sy() ? hit.ey(hit.cy() - ppst.cy()) : MGMath::ZERO);
            Double_t rsMx = (hit.sx() ? ((hit.cx() - ppst.cx()) / mex) : MGMath::ZERO);
            Double_t rsMy = (hit.sy() ? ((hit.cy() - ppst.cy()) / mey) : MGMath::ZERO);
            //ppMsl += (ppst.arg().mscat_ll() * ppst.arg().mscat_ll());
            ppMsl = (ppst.arg().mscat_ll() * ppst.arg().mscat_ll());
            
            vMex.at(cnt_nhit) = mex;
            vMey.at(cnt_nhit) = mey;
            vMsl.at(cnt_nhit) = ppMsl;

            if (hit.sx()) {
                for (Int_t elm = 0; elm < 5; ++elm)
                    JbF(nseq, elm) = (ppjbF(0, elm) / mex);
                CvM(nseq, nseq) = MGMath::ONE + (ppMsl/mex/mex);
                RsM(nseq) = rsMx;
                nseq++;
            }
            
            if (hit.sy()) {
                for (Int_t elm = 0; elm < 5; ++elm)
                    JbF(nseq, elm) = (ppjbF(1, elm) / mey);
                CvM(nseq, nseq) = MGMath::ONE + (ppMsl/mey/mey);
                RsM(nseq) = rsMy; 
                nseq++;
            }
        
            cnt_nhit++;
            if (!hit.sx()) hit.set_dummy_x(ppst.cx());
        }
        if (cnt_nhit != hits_.size()) break;

        for (Int_t im = 0; im < hits_.size(); ++im) {
            if (seqx_.at(im) >= 0) {
                for (Int_t is = seqx_.at(im)+1; is < maps_.size(); ++is) {
                    if ((maps_.at(is)/hits_.size()) != 0) continue;
                    Short_t  ih  = (maps_.at(is) % hits_.size());
                    Double_t msl = (vMsl.at(im) / vMex.at(im) / vMex.at(ih));
                    CvM(seqx_.at(im), is) += msl;
                    CvM(is, seqx_.at(im)) += msl;
                }
            }
            if (seqy_.at(im) >= 0) {
                for (Int_t is = seqy_.at(im)+1; is < maps_.size(); ++is) {
                    if ((maps_.at(is)/hits_.size()) != 1) continue;
                    Short_t  ih  = (maps_.at(is) % hits_.size());
                    Double_t msl = (vMsl.at(im) / vMey.at(im) / vMey.at(ih));
                    CvM(seqy_.at(im), is) += msl;
                    CvM(is, seqy_.at(im)) += msl;
                }
            }
        }
        Double_t detCvM = 0.0;
        CvM.Invert(&detCvM);

        TMtxD JbFTrans(TMtxD::kTransposed, JbF);
        TVecD GdF = JbFTrans * CvM * RsM;
        TMtxD CvF = JbFTrans * CvM * JbF;
        chi = CvM.Similarity(RsM);

        Double_t detCvF = 0.0;
        CvF.Invert(&detCvF);
        TVecD&& RslG = (CvF * GdF);

        part_.set_state_with_uxy(
            part_.cx() + RslG(0),
            part_.cy() + RslG(1),
            part_.cz(),
            part_.ux() + RslG(2),
            part_.uy() + RslG(3),
            ((ortt_ == Orientation::kDownward) ? -1 : 1)
        );
        part_.set_eta(part_.eta() + RslG(4));

        Double_t nchi = ((chi) / static_cast<Double_t>(ndf));
        Double_t nchi_rat = std::fabs((nchi - nchi_) / (nchi + nchi_ + CONVG_EPSILON));
        Bool_t   sign     = (MGNumc::Compare(nchi - nchi_, std::min(nchi, nchi_) * CONVG_EPSILON) < 0);
        
        nchi_ = nchi;
        ndf_  = ndf;

        curSucc = (curIter >= LMTL_ITER && (sign && MGNumc::Compare(nchi_rat, CONVG_TOLERANCE) < 0));
        
        succ = (preSucc && curSucc);
        preSucc = curSucc;
        
        //std::cout << Form("IT %d (RIG %14.8f CHI %14.8f)\n", curIter, part_.rig(), nchi_); // determinate problem
        
        curIter++;
    }
    if (!succ) std::cout << Form("FAIL. %d (MOM %14.8f CHI %14.8f)\n", curIter, part_.mom(), nchi_);
    
    return succ;
}
*/








//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////








/*

Bool_t PhyTr::fit_physics() {
    std::vector<PhyArg> args(hits_.size(), PhyArg(sw_mscat_, sw_eloss_));
    Int_t NPAR_G = 5;
    Int_t NPAR_L = 2 * (hits_.size());
    
    Bool_t succ = false;
    Bool_t preSucc = false;
    Bool_t curSucc = false;
    
    Int_t  curIter = 1;
    while (curIter <= LMTU_ITER && !succ) {
        Double_t chi = MGMath::ZERO;
        Int_t    ndf = -5;
        
        TMtxD JbFJ(nseq_, NPAR_G+NPAR_L);
        TVecD RsM(nseq_);
        TVecD RsL(NPAR_L);

        Int_t nseq = 0;
        Int_t cnt_nhit = 0;
        PhySt ppst(part_);
        PhyJb::SMtxDGG&&            ppjbF = SMtxId();
        std::vector<PhyJb::SMtxDGL> ppjbJ(hits_.size());
        for (auto&& hit : hits_) {
            PhyArg& arg = args.at(cnt_nhit);
            ppst.arg() = arg;
            PhyJb curjb;

            std::cout << Form("IT %2d HIT %2d CURSAT %14.8f %14.8f\n", 
                curIter,
                cnt_nhit,
                ppst.arg().tauu(),
                ppst.arg().rhou()
            );

            if (!PropMgnt::PropToZ(hit.cz(), ppst, nullptr, &curjb)) break;
            ppst.symbk();
            
            ppjbF              = curjb.gg() * ppjbF;
            ppjbJ.at(cnt_nhit) = curjb.gl();
            for (Int_t ih = 0; ih < cnt_nhit; ++ih)
                ppjbJ.at(ih) = curjb.gg() * ppjbJ.at(ih);
            
            Double_t mex = (hit.sx() ? hit.ex(hit.cx() - ppst.cx()) : MGMath::ZERO);
            Double_t mey = (hit.sy() ? hit.ey(hit.cy() - ppst.cy()) : MGMath::ZERO);
            
            SVecD<2> rsM;
            rsM(0) = (hit.sx() ? ((hit.cx() - ppst.cx()) / mex) : MGMath::ZERO);
            rsM(1) = (hit.sy() ? ((hit.cy() - ppst.cy()) / mey) : MGMath::ZERO);

            PhyJb::SMtxDXYG&&            subJbF = PhyJb::SubXYG(ppjbF);
            std::vector<PhyJb::SMtxDXYL> subJbL;
            for (Int_t ih = 0; ih <= cnt_nhit; ++ih)
                 subJbL.push_back(PhyJb::SubXYL(ppjbJ.at(ih)));
            
            if (hit.sx()) {
                for (Int_t elm = 0; elm < 5; ++elm)
                    JbFJ(nseq, elm) = (subJbF(0, elm) / mex);
                
                for (Int_t ih = 0; ih <= cnt_nhit; ++ih)
                    for (Int_t elm = 0; elm < 2; ++elm)
                        JbFJ(nseq, NPAR_G+2*ih+elm) = (subJbL.at(ih)(0, elm) / mex);
                
                RsM(nseq) = rsM(0); 
                nseq++;
            }
            
            if (hit.sy()) {
                for (Int_t elm = 0; elm < 5; ++elm)
                    JbFJ(nseq, elm) = (subJbF(1, elm) / mey);
                
                for (Int_t ih = 0; ih <= cnt_nhit; ++ih)
                    for (Int_t elm = 0; elm < 2; ++elm)
                        JbFJ(nseq, NPAR_G+2*ih+elm) = (subJbL.at(ih)(1, elm) / mey);
                
                RsM(nseq) = rsM(1); 
                nseq++;
            }
            RsL(2*cnt_nhit+0) = MGMath::NEG * arg.tauu();
            RsL(2*cnt_nhit+1) = MGMath::NEG * arg.rhou();

            if (hit.sx()) { ndf++; chi += rsM(0) * rsM(0); }
            if (hit.sy()) { ndf++; chi += rsM(1) * rsM(1); }
            chi += arg.tauu() * arg.tauu();
            chi += arg.rhou() * arg.rhou();
        
            cnt_nhit++;
            if (!hit.sx()) hit.set_dummy_x(ppst.cx());
        }
        if (cnt_nhit != hits_.size()) break;
        
        TMtxD JbFJTrans(TMtxD::kTransposed, JbFJ);
        TVecD GdFJ = JbFJTrans * RsM;
        for (Int_t il = 0; il < NPAR_L; ++il)
            GdFJ(NPAR_G+il) += RsL(il);
        TMtxD CvFJ = JbFJTrans * JbFJ;
        for (Int_t il = 0; il < NPAR_L; ++il)
            CvFJ(NPAR_G+il, NPAR_G+il) += MGMath::ONE;

        // testcode
        //Bool_t udloc = (curIter%5!=0);
        //if (udloc) {
            GdFJ.ResizeTo(NPAR_G+2, NPAR_G+NPAR_L-1);
            CvFJ.ResizeTo(NPAR_G+2, NPAR_G+NPAR_L-1, NPAR_G+2, NPAR_G+NPAR_L-1);
        //}
        //else {
        //    GdFJ.ResizeTo(5);
        //    CvFJ.ResizeTo(5, 5);
        //}

        Double_t detCvFJ = 0.0;
        CvFJ.Invert(&detCvFJ);
        //if (MGNumc::EqualToZero(detCvFJ)) break; // TODO fix
        
        TVecD&& RslGL = (CvFJ * GdFJ);
        
        // testcode
        Double_t learn = ((curIter+1)*0.1);
        //if (learn > 1.0) learn = 1.0;
        learn = 0.5;

        //if (!udloc) {
        //part_.set_state_with_uxy(
        //    part_.cx() + learn * RslGL(0),
        //    part_.cy() + learn * RslGL(1),
        //    part_.cz(),
        //    part_.ux() + learn * RslGL(2),
        //    part_.uy() + learn * RslGL(3),
        //    ((ortt_ == Orientation::kDownward) ? -1 : 1)
        //);
        //part_.set_eta(part_.eta() + learn * RslGL(4));
        //}

        //if (udloc)
        for (Int_t ih = 0; ih < hits_.size(); ++ih) {
            if (ih == 0) continue;
            PhyArg& arg = args.at(ih);
            arg.set_mscat(
                arg.tauu() + learn * RslGL(NPAR_G+2*ih+0),
                arg.rhou() + learn * RslGL(NPAR_G+2*ih+1)
            );
        }
        
        Double_t nchi = ((chi) / static_cast<Double_t>(ndf));
        Double_t nchi_rat = std::fabs((nchi - nchi_) / (nchi + nchi_ + CONVG_EPSILON));
        Bool_t   sign     = (MGNumc::Compare(nchi - nchi_, std::min(nchi, nchi_) * CONVG_EPSILON) < 0);
        
        nchi_ = nchi;
        ndf_  = ndf;

        curSucc = (curIter >= LMTM_ITER && (sign && MGNumc::Compare(nchi_rat, CONVG_TOLERANCE) < 0));
        
        succ = (preSucc && curSucc);
        preSucc = curSucc;
        
        std::cout << Form("IT %d (RIG %14.8f CHI %14.8f)\n", curIter, part_.rig(), nchi_); // determinate problem
        
        curIter++;
    }
    if (!succ) std::cout << Form("FAIL. %d (MOM %14.8f CHI %14.8f)\n", curIter, part_.mom(), nchi_);
    
    return succ;
}
*/






/////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
Bool_t PhyTr::fit_semi_simple() {
    Int_t NPAR_G = 5;
    Int_t NPAR_L = 4 * hits_.size();
    
    Bool_t succ = false;
    Bool_t preSucc = false;
    Bool_t curSucc = false;
    
    Int_t  curIter = 1;
    std::vector<PhyArg> args(hits_.size());
    while (curIter <= LMTU_ITER && !succ) {
        Double_t chi = MGMath::ZERO;
        Int_t    ndf = -5;
        
        TMtxD JbFJ(nseq_, NPAR_G+NPAR_L);
        TVecD RsM(nseq_);
        TVecD RsL(NPAR_L);

        Int_t nseq = 0;
        Int_t cnt_nhit = 0;
        PhySt ppst(part_);
        PhyJb::SMtxDGG&&            ppjbF = SMtxId();
        std::vector<PhyJb::SMtxDGL> ppjbJ(hits_.size());
        for (auto&& hit : hits_) {
            PhyJb curjb;
            
            ppst.arg() = args.at(cnt_nhit);
            //std::cout << Form("IT %2d HIT %2d CURSAT %14.8f %14.8f %14.8f %14.8f\n", 
            //    curIter,
            //    cnt_nhit,
            //    ppst.arg().tauu(),
            //    ppst.arg().rhou(),
            //    ppst.arg().taul(),
            //    ppst.arg().rhol()
            //);

            if (!PropMgnt::PropToZ(hit.cz(), ppst, nullptr, &curjb)) break;
            
            ppjbF              = curjb.gg() * ppjbF;
            ppjbJ.at(cnt_nhit) = curjb.gl();
            for (Int_t ih = 0; ih < cnt_nhit; ++ih)
                ppjbJ.at(ih) = curjb.gg() * ppjbJ.at(ih);
            
            Double_t mex = (hit.sx() ? hit.ex(hit.cx() - ppst.cx()) : MGMath::ZERO);
            Double_t mey = (hit.sy() ? hit.ey(hit.cy() - ppst.cy()) : MGMath::ZERO);
            
            SVecD<2> rsM;
            rsM(0) = (hit.sx() ? ((hit.cx() - ppst.cx()) / mex) : MGMath::ZERO);
            rsM(1) = (hit.sy() ? ((hit.cy() - ppst.cy()) / mey) : MGMath::ZERO);

            PhyJb::SMtxDXYG&&            subJbF = PhyJb::SubXYG(ppjbF);
            std::vector<PhyJb::SMtxDXYL> subJbL;
            for (Int_t ih = 0; ih <= cnt_nhit; ++ih)
                 subJbL.push_back(PhyJb::SubXYL(ppjbJ.at(ih)));
            
            if (hit.sx()) {
                for (Int_t elm = 0; elm < 5; ++elm)
                    JbFJ(nseq, elm) = (subJbF(0, elm) / mex);
                
                for (Int_t ih = 0; ih <= cnt_nhit; ++ih)
                    for (Int_t elm = 0; elm < 4; ++elm)
                        JbFJ(nseq, NPAR_G+4*ih+elm) = (subJbL.at(ih)(0, elm) / mex);
                
                RsM(nseq) = rsM(0); 
                nseq++;
            }
            
            if (hit.sy()) {
                for (Int_t elm = 0; elm < 5; ++elm)
                    JbFJ(nseq, elm) = (subJbF(1, elm) / mey);
                
                for (Int_t ih = 0; ih <= cnt_nhit; ++ih)
                    for (Int_t elm = 0; elm < 4; ++elm)
                        JbFJ(nseq, NPAR_G+4*ih+elm) = (subJbL.at(ih)(1, elm) / mey);
                
                RsM(nseq) = rsM(1); 
                nseq++;
            }
            RsL(4*cnt_nhit+0) = MGMath::NEG * ppst.arg().tauu();
            RsL(4*cnt_nhit+1) = MGMath::NEG * ppst.arg().rhou();
            RsL(4*cnt_nhit+2) = MGMath::NEG * ppst.arg().taul();
            RsL(4*cnt_nhit+3) = MGMath::NEG * ppst.arg().rhol();

            if (hit.sx()) { ndf++; chi += rsM(0) * rsM(0); }
            if (hit.sy()) { ndf++; chi += rsM(1) * rsM(1); }
            chi += ppst.arg().tauu() * ppst.arg().tauu();
            chi += ppst.arg().rhou() * ppst.arg().rhou();
            chi += ppst.arg().taul() * ppst.arg().taul();
            chi += ppst.arg().rhol() * ppst.arg().rhol();
        
            cnt_nhit++;
            if (!hit.sx()) hit.set_dummy_x(ppst.cx());
        }
        if (cnt_nhit != hits_.size()) break;
        
        TMtxD JbFJTrans(TMtxD::kTransposed, JbFJ);
        TVecD GdFJ = JbFJTrans * RsM;
        for (Int_t il = 0; il < NPAR_L; ++il)
            GdFJ(NPAR_G+il) += RsL(il);
        TMtxD CvFJ = JbFJTrans * JbFJ;
        for (Int_t il = 0; il < NPAR_L; ++il)
            CvFJ(NPAR_G+il, NPAR_G+il) += MGMath::ONE;

        // testcode
        //Bool_t udloc = (curIter%5!=0);
        //if (udloc) {
        //    GdFJ.ResizeTo(NPAR_G, NPAR_G+NPAR_L-1);
        //    CvFJ.ResizeTo(NPAR_G, NPAR_G+NPAR_L-1, NPAR_G, NPAR_G+NPAR_L-1);
        //}
        //else {
        //    GdFJ.ResizeTo(5);
        //    CvFJ.ResizeTo(5, 5);
        //}

        Double_t detCvFJ = 0.0;
        CvFJ.Invert(&detCvFJ);
        //if (MGNumc::EqualToZero(detCvFJ)) break; // TODO fix
        
        TVecD&& RslGL = (CvFJ * GdFJ);
        
        // testcode
        Double_t learn = ((curIter+1)*0.1);
        if (learn > 1.0) learn = 1.0;
        learn = 1.0;

        //if (!udloc) {
        part_.set_state_with_uxy(
            part_.cx() + learn * RslGL(0),
            part_.cy() + learn * RslGL(1),
            part_.cz(),
            part_.ux() + learn * RslGL(2),
            part_.uy() + learn * RslGL(3),
            ((ortt_ == Orientation::kDownward) ? -1 : 1)
        );
        part_.set_eta(part_.eta() + learn * RslGL(4));
        //}

        //if (udloc)
        for (Int_t ih = 0; ih < hits_.size(); ++ih) {
            PhyArg& arg = args.at(ih);
            arg.set_mscat(
                arg.tauu() + learn * RslGL(NPAR_G+4*ih+0),
                arg.rhou() + learn * RslGL(NPAR_G+4*ih+1),
                arg.taul() + learn * RslGL(NPAR_G+4*ih+2),
                arg.rhol() + learn * RslGL(NPAR_G+4*ih+3)
            );
            //std::cout << Form("IT %2d HIT %2d UPDATE VALUE %14.8f %14.8f %14.8f %14.8f\n",
            //    curIter,
            //    ih,
            //    RslGL(NPAR_G+4*ih+0),
            //    RslGL(NPAR_G+4*ih+1),
            //    RslGL(NPAR_G+4*ih+2),
            //    RslGL(NPAR_G+4*ih+3)
            //);
            //std::cout << Form("IT %2d HIT %2d UPDATE STATE %14.8f %14.8f %14.8f %14.8f\n", 
            //    curIter,
            //    ih,
            //    arg.tauu(),
            //    arg.rhou(),
            //    arg.taul(),
            //    arg.rhol()
            //);
        }
        
        Double_t nchi = ((chi) / static_cast<Double_t>(ndf));
        Double_t nchi_rat = std::fabs((nchi - nchi_) / (nchi + nchi_ + CONVG_EPSILON));
        Bool_t   sign     = (MGNumc::Compare(nchi - nchi_, std::min(nchi, nchi_) * CONVG_EPSILON) < 0);
        
        nchi_ = nchi;
        ndf_  = ndf;

        curSucc = (curIter >= LMTM_ITER && (sign && MGNumc::Compare(nchi_rat, CONVG_TOLERANCE) < 0));
        
        succ = (preSucc && curSucc);
        preSucc = curSucc;
        
        std::cout << Form("IT %d (RIG %14.8f CHI %14.8f)\n", curIter, part_.rig(), nchi_); // determinate problem
        
        curIter++;
    }
    if (!succ) std::cout << Form("FAIL. %d (MOM %14.8f CHI %14.8f)\n", curIter, part_.mom(), nchi_);
    
    return succ;
}
*/
/////////////////////////////////////////////////////////////////////////////////////////////////////////



/*
Bool_t PhyTr::fit_semi_simple() {
    Bool_t succ = false;
    Bool_t preSucc = false;
    Bool_t curSucc = false;
    Int_t  curIter = 1;
    //Double_t tuneIon = MGMath::ZERO;
    while (curIter <= LMTU_ITER && !succ) {
        std::vector<PhyJb::SMtxDGG> vJbF(hits_.size());
        std::vector<PhyJb::SMtxDGL> vJbJ(hits_.size());
        std::vector<SVecD<2>>       vRsM(hits_.size());
        std::vector<SVecD<2>>       vErM(hits_.size());

        PhySt ppst(part_);
        Int_t cnt_nhit = 0;
        for (auto&& hit : hits_) {
            PhyJb curjb;
            //ppst.arg().set_tune_ion(tuneIon);
            if (!PropMgnt::PropToZ(hit.cz(), ppst, nullptr, &curjb)) break;

            SVecD<2> rsM(hit.cx() - ppst.cx(), hit.cy() - ppst.cy());
            SVecD<2> erM((hit.sx() ? (hit.ex() * hit.ex()) : MGMath::ZERO), (hit.sy() ? (hit.ey() * hit.ey()) : MGMath::ZERO));

            vJbF.at(cnt_nhit) = std::move(curjb.gg());
            vJbJ.at(cnt_nhit) = std::move(curjb.gl());
            vRsM.at(cnt_nhit) = rsM;
            vErM.at(cnt_nhit) = erM;
        
            cnt_nhit++;
            if (!hit.sx()) hit.set_dummy_x(ppst.cx());
        }
        if (cnt_nhit != hits_.size()) break;
      
        PhyJb::SMtxDGG ppjbF = SMtxId(); 
        TMtxD JbF(nseq_, 5);
        TMtxD JbJ(nseq_, 4*hits_.size());
        TMtxD CvM(nseq_, nseq_);
        TVecD RsM(nseq_);
        for (Int_t it = 0; it < hits_.size(); ++it) {
            HitSt& hit = hits_.at(it);
            Short_t subID = hit.seqID();

            ppjbF *= vJbF.at(it);
            for (Int_t jt = 0; jt < it; ++jt)
                vJbJ.at(jt) = std::move(vJbF.at(it) * vJbJ.at(jt));
            
            SVecD<2>& rsM = vRsM.at(it);
            SVecD<2>& erM = vErM.at(it);
            
            if (hit.sx()) {
                for (Int_t elm = 0; elm < 5; ++elm) JbF(subID, elm) = ppjbF(0, elm);
                for (Int_t jt = 0; jt <= it; ++jt)
                    for (Int_t elm = 0; elm < 4; ++elm)
                        JbJ(subID, elm) = (vJbJ.at(jt))(0, elm);
                CvM(subID, subID) = erM(0);
                RsM(subID) = rsM(0);
                subID++;
            }
            if (hit.sy()) {
                for (Int_t elm = 0; elm < 5; ++elm) JbF(subID, elm) = ppjbF(1, elm);
                for (Int_t jt = 0; jt <= it; ++jt)
                    for (Int_t elm = 0; elm < 4; ++elm)
                        JbJ(subID, elm) = (vJbJ.at(jt))(1, elm);
                CvM(subID, subID) = erM(1);
                RsM(subID) = rsM(1);
                subID++;
            }
        }
        TMtxD JbFTrans(TMtxD::kTransposed, JbF);
        TMtxD JbJTrans(TMtxD::kTransposed, JbJ);
        
        // Interaction
        CvM += (JbJ * JbJTrans);

        Double_t detCvM = 0.0;
        CvM.Invert(&detCvM);
        //if (MGNumc::EqualToZero(detCvM)) break; // TODO fix
                
        Double_t chi = CvM.Similarity(RsM);

        TMtxD&& CvGG = (JbFTrans * CvM * JbF);
        TVecD&& grdG = (JbFTrans * CvM * RsM);

        Double_t detCvGG = 0.0;
        CvGG.Invert(&detCvGG);
        //if (MGNumc::EqualToZero(detCvGG)) break; // TODO fix

        TVecD&& rslG = (CvGG * grdG);

        part_.set_state_with_uxy(
            part_.cx() + rslG(0),
            part_.cy() + rslG(1),
            part_.cz(),
            part_.ux() + rslG(2),
            part_.uy() + rslG(3),
            ((ortt_ == Orientation::kDownward) ? -1 : 1)
        );
        part_.set_eta(part_.eta() + rslG(4));
        
        Int_t    ndf  = (nseq_ - 5);
        Double_t nchi = ((chi) / static_cast<Double_t>(ndf));
        Double_t nchi_rat = std::fabs((nchi - nchi_) / (nchi + nchi_ + CONVG_EPSILON));
        Bool_t   sign     = (MGNumc::Compare(nchi - nchi_, std::min(nchi, nchi_) * CONVG_EPSILON) < 0);
        
        ndf_ = ndf;
        nchi_ = nchi;

        curSucc = (curIter >= LMTM_ITER && (sign && MGNumc::Compare(nchi_rat, CONVG_TOLERANCE) < 0));
        
        succ = (preSucc && curSucc);
        preSucc = curSucc;
        
        //if (MGNumc::EqualToZero(detCvM)) std::cout << Form("IT %d (MOM %14.8f CHI %14.8f)\n", curIter, part_.mom(), nchi_); // determinate problem
        //std::cout << Form("IT %d (MOM %14.8f CHI %14.8f)\n", curIter, part_.mom(), nchi_); // determinate problem
        curIter++;
    }
    if (!succ) std::cout << Form("FAIL. %d (MOM %14.8f CHI %14.8f)\n", curIter, part_.mom(), nchi_);
    
    return succ;
}
*/











/*
Bool_t PhyTr::fit_physics() {
    Bool_t succ = false;
    Bool_t preSucc = false;
    Bool_t curSucc = false;
    Int_t  curIter = 1;
    std::vector<PhyArg> args(hits_.size()-1);
    while (curIter <= LMTU_ITER && !succ) {
        Int_t ndfx = 0;
        Int_t ndfy = 0;

        std::vector<SVecD<2>>       chiM(hits_.size());
        std::vector<SVecD<2>>       chiL(hits_.size()-1);

        std::vector<SVecD<2>>       resM(hits_.size());
        std::vector<SVecD<2>>       resL(hits_.size()-1);
        
        std::vector<SMtxSymD<2>>    covM(hits_.size());
        std::vector<SMtxSymD<2>>    covL(hits_.size()-1);
        
        std::vector<SVecD<2>>       gradM(hits_.size());
        std::vector<SVecD<2>>       gradL(hits_.size()-1);

        std::vector<PhyJb::SMtxDGG> jbG(hits_.size());
        std::vector<PhyJb::SMtxDGL> jbL(hits_.size()*hits_.size());

        PhyJb::SMtxDGG ppjb = SMtxIdSym5D;
        PhySt ppst(part_);
        ppst.arg().zero();
        Int_t cnt_nhit = 0;
        for (auto&& hit : hits_) {
            PhyJb curjb;
            if (cnt_nhit != 0) ppst.arg() = args.at(cnt_nhit-1);
            if (!PropMgnt::PropToZ(hit.cz(), ppst, nullptr, &curjb)) break;
            if (cnt_nhit != 0) ppst.symbk();
            ppjb = std::move(PhyJb::Multiply(curjb.gg(), ppjb));

            // Hit Measurement
            SVecD<2> mres(ppst.cx() - hit.cx(), ppst.cy() - hit.cy());
            (resM.at(cnt_nhit))(0) = (hit.sx() ? mres(0) : MGMath::ZERO);
            (resM.at(cnt_nhit))(1) = (hit.sy() ? mres(1) : MGMath::ZERO);

            SVecD<2> merr(hit.ex(mres(0)), hit.ey(mres(1)));
            (covM.at(cnt_nhit))(0, 0) = (hit.sx() ? (MGMath::ONE / merr(0) / merr(0)) : MGMath::ZERO);
            (covM.at(cnt_nhit))(1, 1) = (hit.sy() ? (MGMath::ONE / merr(1) / merr(1)) : MGMath::ZERO);
            
            (gradM.at(cnt_nhit))(0) = (hit.sx() ? ((covM.at(cnt_nhit))(0, 0) * mres(0)) : MGMath::ZERO);
            (gradM.at(cnt_nhit))(1) = (hit.sy() ? ((covM.at(cnt_nhit))(1, 1) * mres(1)) : MGMath::ZERO);

            (chiM.at(cnt_nhit))(0) = (hit.sx() ? (mres(0) * (gradM.at(cnt_nhit))(0)) : MGMath::ZERO);
            (chiM.at(cnt_nhit))(1) = (hit.sy() ? (mres(1) * (gradM.at(cnt_nhit))(1)) : MGMath::ZERO);
            
            jbG.at(cnt_nhit) = ppjb;
          
            // Physics Measurement
            if (cnt_nhit != 0) {
                SVecD<2> lres(args.at(cnt_nhit-1).tauu(), args.at(cnt_nhit-1).rhou());
                (resL.at(cnt_nhit-1))(0) = lres(0);
                (resL.at(cnt_nhit-1))(1) = lres(1);
                
                SVecD<2> lerr(args.at(cnt_nhit-1).etauu(), args.at(cnt_nhit-1).erhou());
                (covL.at(cnt_nhit-1))(0, 0) = (MGMath::ONE / lerr(0));
                (covL.at(cnt_nhit-1))(1, 1) = (MGMath::ONE / lerr(1));

                (gradL.at(cnt_nhit-1))(0) = ((covL.at(cnt_nhit-1))(0, 0) * lres(0));
                (gradL.at(cnt_nhit-1))(1) = ((covL.at(cnt_nhit-1))(1, 1) * lres(1));
                
                (chiL.at(cnt_nhit-1))(0) = (lres(0) * (gradL.at(cnt_nhit-1))(0));
                (chiL.at(cnt_nhit-1))(1) = (lres(1) * (gradL.at(cnt_nhit-1))(1));
            }

            if (cnt_nhit != 0) {
                Int_t currentM = cnt_nhit;
                Int_t currentL = (cnt_nhit-1);
                jbL.at(currentM*hits_.size()+currentL) = std::move(curjb.gl());
                for (Int_t itL = 0; itL < currentL; ++itL) {
                    Int_t preL = (currentM-1) * hits_.size() + itL;
                    Int_t curL = (currentM)   * hits_.size() + itL;
                    jbL.at(curL) = std::move(PhyJb::Multiply(curjb.gg(), jbL.at(preL)));
                }
            }
           
            cnt_nhit++;
            if (!hit.sx()) hit.set_dummy_x(ppst.cx());
            
            ndfx += hit.sx();
            ndfy += hit.sy();
        }
        if (cnt_nhit != hits_.size()) return false;
        
        std::vector<PhyJb::SMtxDXYG> jbXYG(hits_.size());
        std::vector<PhyJb::SMtxDXYL> jbXYL(hits_.size()*hits_.size());
        for (Int_t ih = 0; ih < hits_.size(); ++ih) jbXYG.at(ih) = std::move(PhyJb::SubXYG( jbG.at(ih) ));
        for (Int_t ih = 1; ih < hits_.size(); ++ih) {
            for (Int_t itL = 0; itL < ih; ++itL) {
                jbXYL.at(ih*hits_.size()+itL) = std::move(PhyJb::SubXYL( jbL.at(ih*hits_.size()+itL) ));
            }
        }

        ///////////////////
        TVecD grad(5+2*(cnt_nhit-1));
        {
            SVecD<5> gradG;
            for (Int_t ih = 0; ih < hits_.size(); ++ih) gradG += (LA::Transpose(jbXYG.at(ih)) * gradM.at(ih));
            for (Int_t ig = 0; ig < 5; ++ig) grad(ig) = gradG(ig);

            for (Int_t itL = 0; itL < hits_.size()-1; ++itL) {
                grad(5+2*itL+0) = gradL.at(itL)(0);
                grad(5+2*itL+1) = gradL.at(itL)(1);
            }
        }


        TMtxD cov(5+2*(cnt_nhit-1), 5+2*(cnt_nhit-1));
        {
            SMtxD<5> covGG;
            for (Int_t ih = 0; ih < hits_.size(); ++ih) covGG += LA::SimilarityT(jbXYG.at(ih), covM.at(ih));
            for (Int_t ig = 0; ig < 5; ++ig) {
                for (Int_t jg = 0; jg < 5; ++jg) {
                    cov(ig, jg) += covGG(ig, jg);
                }
            }
            
            for (Int_t itL = 0; itL < hits_.size()-1; ++itL) {
                cov(5+2*itL+0, 5+2*itL+0) += (covL.at(itL))(0, 0);
                cov(5+2*itL+1, 5+2*itL+1) += (covL.at(itL))(1, 1);
           
                // cross item  ll with m
                //for (Int_t jtL = 0; jtL < hits_.size()-1; ++jtL) {
                //    SMtxD<2>&& covLL = LA::Transpose(jbXYL.at(ih)) * covM.at(ih) * LA::Transpose(jbXYL.at(ih));
                //}
            }

            // cross item gl
        }
        Double_t det = 0;
        cov.Invert(&det);
        if (MGNumc::EqualToZero(det)) return false;

        TVecD rsl = cov * grad;
        
        part_.set_state_with_uxy(
            part_.cx() - rsl(0),
            part_.cy() - rsl(1),
            part_.cz(),
            part_.ux() - rsl(2),
            part_.uy() - rsl(3),
            ((ortt_ == Orientation::kDownward) ? -1 : 1)
        );
        part_.set_eta(part_.eta() - rsl(4));
        
        for (Int_t itL = 0; itL < hits_.size()-1; ++itL) {
            args.at(itL).set_mscat(rsl(5+2*itL+0), rsl(5+2*itL+1));
        }
        
        Int_t ndf = (ndfx + ndfy - 5);
        Double_t chix = 0, chiy = 0;
        Double_t chit = 0, chir = 0;
        for (Int_t ih = 0; ih < hits_.size(); ++ih) { chix += (chiM.at(ih))(0); chiy += (chiM.at(ih))(1); }
        for (Int_t itL = 0; itL < hits_.size()-1; ++itL) { chit += (chiL.at(itL))(0); chir += (chiL.at(itL))(1); }
        Double_t chi = (chix + chiy + chit + chir);
        Double_t nchi = ((chi) / static_cast<Double_t>(ndf));
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
*/
/////////////////////////////////////////////////////////////////////////////////////////////////


} // namespace TrackSys


#endif // __TRACKLibs_PhyTr_C__
