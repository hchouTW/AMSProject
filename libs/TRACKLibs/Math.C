#ifndef __TRACKLibs_Math_C__
#define __TRACKLibs_Math_C__


namespace TrackSys {
        

MultiGauss::MultiGauss(Opt opt, long double sgm) : MultiGauss()  {
    multi_gauss_.push_back(std::make_pair(MGMath::ONE, sgm));
    bound_.first  = sgm;
    bound_.second = sgm;
    robust_ = opt;
}


MultiGauss::MultiGauss(Opt opt, long double wgt1, long double sgm1, long double wgt2, long double sgm2) : MultiGauss()  {
    long double norm = wgt1 + wgt2;
    multi_gauss_.push_back(std::make_pair(wgt1/norm, sgm1));
    multi_gauss_.push_back(std::make_pair(wgt2/norm, sgm2));
    bound_.first  = std::min(sgm1, sgm2);
    bound_.second = std::max(sgm1, sgm2);
    robust_ = opt;
}


MultiGauss::MultiGauss(Opt opt, long double wgt1, long double sgm1, long double wgt2, long double sgm2, long double wgt3, long double sgm3) : MultiGauss()  {
    long double norm = wgt1 + wgt2 + wgt3;
    multi_gauss_.push_back(std::make_pair(wgt1/norm, sgm1));
    multi_gauss_.push_back(std::make_pair(wgt2/norm, sgm2));
    multi_gauss_.push_back(std::make_pair(wgt3/norm, sgm3));
    bound_.first  = std::min(std::min(sgm1, sgm2), sgm3);
    bound_.second = std::max(std::max(sgm1, sgm2), sgm3);
    robust_ = opt;
}


MultiGauss::MultiGauss(Opt opt, long double wgt1, long double sgm1, long double wgt2, long double sgm2, long double wgt3, long double sgm3, long double wgt4, long double sgm4) : MultiGauss()  {
    long double norm = wgt1 + wgt2 + wgt3 + wgt4;
    multi_gauss_.push_back(std::make_pair(wgt1/norm, sgm1));
    multi_gauss_.push_back(std::make_pair(wgt2/norm, sgm2));
    multi_gauss_.push_back(std::make_pair(wgt3/norm, sgm3));
    multi_gauss_.push_back(std::make_pair(wgt4/norm, sgm4));
    bound_.first  = std::min(std::min(std::min(sgm1, sgm2), sgm3), sgm4);
    bound_.second = std::max(std::max(std::max(sgm1, sgm2), sgm3), sgm4);
    robust_ = opt;
}


MultiGauss::MultiGauss(Opt opt, long double wgt1, long double sgm1, long double wgt2, long double sgm2, long double wgt3, long double sgm3, long double wgt4, long double sgm4, long double wgt5, long double sgm5) : MultiGauss()  {
    long double norm = wgt1 + wgt2 + wgt3 + wgt4 + wgt5;
    multi_gauss_.push_back(std::make_pair(wgt1/norm, sgm1));
    multi_gauss_.push_back(std::make_pair(wgt2/norm, sgm2));
    multi_gauss_.push_back(std::make_pair(wgt3/norm, sgm3));
    multi_gauss_.push_back(std::make_pair(wgt4/norm, sgm4));
    multi_gauss_.push_back(std::make_pair(wgt5/norm, sgm5));
    bound_.first  = std::min(std::min(std::min(std::min(sgm1, sgm2), sgm3), sgm4), sgm5);
    bound_.second = std::max(std::max(std::max(std::max(sgm1, sgm2), sgm3), sgm4), sgm5);
    robust_ = opt;
    
    // find max sigma
    //long double max_sgm = std::max_element(
    //            multi_gauss_.begin(), multi_gauss_.end(), 
    //            [](const std::pair<long double, long double>& gauss1, const std::pair<long double, long double>& gauss2) { return (gauss1.second < gauss2.second); }
    //        )->second;
}


long double MultiGauss::efft_sgm(long double r) const {
    if (multi_gauss_.size() == 0) return MGMath::ZERO;

    long double sigma = 1.0L;
    long double absr = std::abs(r);
    if (multi_gauss_.size() == 1) sigma = multi_gauss_.at(0).second;
    else {
        long double ttl_wgt = MGMath::ZERO;
        long double inv_nrm = MGMath::ZERO;
        for (auto&& gauss : multi_gauss_) {
            long double res = (absr / gauss.second);
            long double nrm = (gauss.second / bound_.second);
            long double prb = (gauss.first / nrm) * std::exp(-MGMath::HALF * res * res);
            ttl_wgt += prb;
            inv_nrm += prb / (nrm * nrm);
        }
        
        sigma = bound_.second * ((MGNumc::Compare(ttl_wgt, LMTL_PROB_) <= 0 || MGNumc::EqualToZero(inv_nrm)) ? 1.0L : std::sqrt(ttl_wgt / inv_nrm));
        if (!MGNumc::Valid(sigma) || MGNumc::Compare(sigma, bound_.second) > 0 || MGNumc::Compare(sigma) <= 0) sigma = bound_.second;
    }
   
    // Robust Method (Modify-Cauchy)
    //if (robust_ == Opt::ROBUST) {
    //    long double nrmr   = absr / (ROBUST_SGM_ * sigma);
    //    long double cauchy = ((MGNumc::EqualToZero(nrmr)) ? 1.0L : (nrmr / std::sqrt(std::log1p(nrmr * nrmr))));
    //    long double corr   = std::pow(cauchy, ROBUST_POW_);
    //    if (MGNumc::Valid(corr) && MGNumc::Compare(corr, 1.0L) > 0) sigma *= corr;
    //}
   
    // Robust Method (Modify-Cauchy)
    if (robust_ == Opt::ROBUST) {
        long double sftnrmr = (absr / sigma) - ROBUST_SGM_;
        if (MGNumc::Compare(sftnrmr) > 0) {
            long double cauchy = (sftnrmr / std::sqrt(std::log1p(sftnrmr * sftnrmr)));
            long double modify_cauchy = ((!MGNumc::Valid(cauchy) || MGNumc::Compare(cauchy, 1.0L) < 0) ? 1.0L : std::sqrt(cauchy));
            if (MGNumc::Valid(modify_cauchy)) sigma *= modify_cauchy;
        }
    }

    return sigma;
}


long double MultiGauss::rndm() {
    if (multi_gauss_.size() == 0) return 0.0L;
    if (multi_gauss_.size() == 1) return (multi_gauss_.at(0).second * MGRndm::NormalGaussian());
    
    if (rand_func_) return rand_func_->GetRandom();
    else {
        std::string fmt = "([0]/[1])*TMath::Exp(-0.5*(x/[1])*(x/[1]))";
        for (unsigned int it = 1; it < multi_gauss_.size(); ++it) {
            int wgtID = 2 * it;
            int sgmID = 2 * it + 1;
            fmt += STR_FMT(" + ([%d]/[%d])*TMath::Exp(-0.5*(x/[%d])*(x/[%d]))", wgtID, sgmID, sgmID, sgmID);
        }
        rand_func_ = new TF1("rand_func", fmt.c_str(), -(MGMath::TWO*ROBUST_SGM_) * bound_.second, (MGMath::TWO*ROBUST_SGM_) * bound_.second);
        rand_func_->SetNpx(NPX_);
        
        int count = 0;
        for (auto&& gauss : multi_gauss_) {
            rand_func_->SetParameter(2*count+0, gauss.first);
            rand_func_->SetParameter(2*count+1, gauss.second);
            count++;
        }
        
        if (rndm_gen_ != nullptr) {
            gRandom->SetSeed(0);
            rndm_gen_ = gRandom;
        }

        return static_cast<long double>(rand_func_->GetRandom());
    }
    
    return 0.0L;
}


} // namesapce TrackSys


#endif // __TRACKLibs_Math_C__
