#ifndef __TRACKLibs_PhyTr_C__
#define __TRACKLibs_PhyTr_C__


namespace TrackSys {
        
    
void PhyTr::HitSort(std::vector<HitSt>& hits, Orientation ortt) {
    if (ortt == Orientation::kDownward) HitSt::Sort(hits, HitSt::Orientation::kDownward);
    else                                HitSt::Sort(hits, HitSt::Orientation::kUpward);
}

    
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


void PhyTr::clear() {
    succ_ = false;
    type_ = PartType::Proton;
    ortt_ = Orientation::kDownward;
    part_.reset(type_);
    hits_.clear();
    nchi_ = 0.;
    ndf_ = 0;
}
        

void PhyTr::print() const {
    std::cerr << "\n======================================\n";
    part_.print();
    for (auto&& hit : hits_) hit.print();
    std::cerr << "======================================\n";
}


PhyTr::PhyTr(const std::vector<HitSt>& hits, const PartType& type, const Orientation& ortt) {
    clear();
    if (!HitCheck(hits)) return;
    type_ = type;
    ortt_ = ortt;
    part_.reset(type);
    hits_ = hits;
    HitSort(hits_);
}


Bool_t PhyTr::fit() {
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
    ////for (Int_t it = 0; it < ntimes; ++it)
    //COUT("=========================================\n");
    //part_.print();
    //COUT("=========================================\n");
    if (!fit_simple()  ) return false;
    //sws.stop();
    //sws.print();
    //part_.print();
    //COUT("CHI %14.8f\n", nchi_);   
   
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
    Double_t tuneIon = MGMath::ZERO;
    while (curIter <= LMTU_ITER && !succ) {
        Int_t ndfx = 0;
        Int_t ndfy = 0;
        Double_t chix = 0.;
        Double_t chiy = 0.;
        SVecD<5>    gradG;
        SMtxSymD<5> covGG;

        Int_t cnt_nhit = 0;
        for (auto&& hit : hits_) {
            PhyJb curjb;
            PhySt curst(part_);
            curst.arg().set_tune_ion(tuneIon);
            if (!PropMgnt::PropToZ(hit.cz(), curst, nullptr, &curjb)) break;
            
            SVecD<2> resM(curst.cx() - hit.cx(), curst.cy() - hit.cy());

            SMtxSymD<2> icovM;
            icovM(0, 0) = (hit.sx() ? (hit.ex() * hit.ex()) : MGMath::ZERO);
            icovM(1, 1) = (hit.sy() ? (hit.ey() * hit.ey()) : MGMath::ZERO);
            
            if (cnt_nhit != 0) {
                SMtxD<5>&& covL = curjb.gl() * LA::Transpose(curjb.gl());
                icovM(0, 0) += covL(0, 0);
                icovM(0, 1) += covL(0, 1);
                icovM(1, 1) += covL(1, 1);
            }
            icovM.Invert();

            SVecD<2>&& gradM = (icovM * resM);
            
            PhyJb::SMtxDXYG&& subjb = PhyJb::SubXYG(curjb.gg());
            gradG += LA::Transpose(subjb) * gradM;
            covGG += LA::SimilarityT(subjb, icovM);

            if (hit.sx()) { ndfx++; chix += (gradM(0) * resM(0)); } 
            if (hit.sy()) { ndfy++; chiy += (gradM(1) * resM(1)); } 
            
            cnt_nhit++;
            if (!hit.sx()) hit.set_dummy_x(curst.cx());
        }
        //------------------//
        // testcode
        //CERR("END LOOP\n");
        if (cnt_nhit != hits_.size()) {
            if (curIter < LMTL_ITER) return false;
            else { tuneIon -= 0.1; curIter++; continue; }
        }

        //if (curIter >= LMTL_ITER && cnt_nhit == hits_.size()) tuneIon += 0.1;
        if (tuneIon < 0.) tuneIon = 0.;
        if (tuneIon > 1.) tuneIon = 1.;
        //CERR("ITER %d ION %14.8f RIG(%14.8f %14.8f) CHI %14.8f\n", curIter, tuneIon, part_.rig(), ppst.rig(), ((chix + chiy) / static_cast<Double_t>((ndfx + ndfy - 5))));
        //------------------//
                
        if (!covGG.Invert()) return false;
        SVecD<5>&& rslG = covGG * gradG;

        part_.set_state_with_uxy(
            part_.cx() - rslG(0),
            part_.cy() - rslG(1),
            part_.cz(),
            part_.ux() - rslG(2),
            part_.uy() - rslG(3),
            ((ortt_ == Orientation::kDownward) ? MGMath::NEG : MGMath::ONE)
        );
        part_.set_eta(part_.eta() - rslG(4));
        
        Int_t    ndf  = (ndfx + ndfy - 5);
        Double_t nchi = ((chix + chiy) / static_cast<Double_t>(ndf));
        Double_t nchi_rat = std::fabs((nchi - nchi_) / (nchi + nchi_ + CONVG_EPSILON));
        Bool_t   sign     = (MGNumc::Compare(nchi - nchi_, CONVG_EPSILON) < 0);
        
        ndf_ = ndf;
        nchi_ = nchi;

        curSucc = (curIter >= LMTM_ITER && (sign && MGNumc::Compare(nchi_rat, CONVG_TOLERANCE) < 0));
        
        succ = (preSucc && curSucc);
        preSucc = curSucc;
        
        curIter++;
    }
    
    return succ;
}


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

        PhyJb::SMtxDGG ppjb;
        PhySt ppst(part_);
        ppst.zero();
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

            SVecD<2>&& merr = std::move(hit.e(mres(0), mres(1)));
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
            ((ortt_ == Orientation::kDownward) ? MGMath::NEG : MGMath::ONE)
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
/////////////////////////////////////////////////////////////////////////////////////////////////


/*
Bool_t PhyTr::fit_semi_simple() {
    return true;
    std::vector<MatFld> mflds(hits_.size());
    {
        PhySt ppst(part_);
        Int_t cnt_nhit = 0;
        for (auto&& hit : hits_) {
            if (!PropMgnt::PropToZ(hit.cz(), ppst, nullptr, &mflds.at(cnt_nhit))) break;
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
            if (!PropMgnt::PropToZ(hit.cz(), ppst, &curjb)) break;
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
                ppst.ux() - rslst(2),
                ppst.uy() - rslst(3),
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
            part_.ux() - rslG(2),
            part_.uy() - rslG(3),
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
*/

/*
Bool_t PhyTr::fit_physics() {
    return true;
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
            if (!PropMgnt::PropToZ(hit.cz(), ppst, &curjb)) break;
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
                ppst.ux() - rslst(2),
                ppst.uy() - rslst(3),
                ((ortt_ == Orientation::kDownward) ? -1 : 1)
            );
            ppst.set_eta(ppst.eta() - rslst(4));

            SMtxSymD<5> udjb(SMtxIdSym5D);
    /////////////////////////////////////////////////////////////////////////////////////
    
    return true;
}
*/

} // namespace TrackSys


#endif // __TRACKLibs_PhyTr_C__
