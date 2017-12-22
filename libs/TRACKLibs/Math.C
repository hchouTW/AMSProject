#ifndef __TRACKLibs_Math_C__
#define __TRACKLibs_Math_C__


namespace TrackSys {
        

MultiGauss::MultiGauss(Opt opt, Double_t sgm) : MultiGauss()  {
    multi_gauss_.push_back(std::make_pair(MGMath::ONE, sgm));
    bound_.first  = sgm;
    bound_.second = sgm;
    robust_ = opt;
}


MultiGauss::MultiGauss(Opt opt, Double_t wgt1, Double_t sgm1, Double_t wgt2, Double_t sgm2) : MultiGauss()  {
    Double_t norm = wgt1 + wgt2;
    multi_gauss_.push_back(std::make_pair(wgt1/norm, sgm1));
    multi_gauss_.push_back(std::make_pair(wgt2/norm, sgm2));
    bound_.first  = std::min(sgm1, sgm2);
    bound_.second = std::max(sgm1, sgm2);
    robust_ = opt;
}


MultiGauss::MultiGauss(Opt opt, Double_t wgt1, Double_t sgm1, Double_t wgt2, Double_t sgm2, Double_t wgt3, Double_t sgm3) : MultiGauss()  {
    Double_t norm = wgt1 + wgt2 + wgt3;
    multi_gauss_.push_back(std::make_pair(wgt1/norm, sgm1));
    multi_gauss_.push_back(std::make_pair(wgt2/norm, sgm2));
    multi_gauss_.push_back(std::make_pair(wgt3/norm, sgm3));
    bound_.first  = std::min(std::min(sgm1, sgm2), sgm3);
    bound_.second = std::max(std::max(sgm1, sgm2), sgm3);
    robust_ = opt;
}


MultiGauss::MultiGauss(Opt opt, Double_t wgt1, Double_t sgm1, Double_t wgt2, Double_t sgm2, Double_t wgt3, Double_t sgm3, Double_t wgt4, Double_t sgm4) : MultiGauss()  {
    Double_t norm = wgt1 + wgt2 + wgt3 + wgt4;
    multi_gauss_.push_back(std::make_pair(wgt1/norm, sgm1));
    multi_gauss_.push_back(std::make_pair(wgt2/norm, sgm2));
    multi_gauss_.push_back(std::make_pair(wgt3/norm, sgm3));
    multi_gauss_.push_back(std::make_pair(wgt4/norm, sgm4));
    bound_.first  = std::min(std::min(std::min(sgm1, sgm2), sgm3), sgm4);
    bound_.second = std::max(std::max(std::max(sgm1, sgm2), sgm3), sgm4);
    robust_ = opt;
}


MultiGauss::MultiGauss(Opt opt, Double_t wgt1, Double_t sgm1, Double_t wgt2, Double_t sgm2, Double_t wgt3, Double_t sgm3, Double_t wgt4, Double_t sgm4, Double_t wgt5, Double_t sgm5) : MultiGauss()  {
    Double_t norm = wgt1 + wgt2 + wgt3 + wgt4 + wgt5;
    multi_gauss_.push_back(std::make_pair(wgt1/norm, sgm1));
    multi_gauss_.push_back(std::make_pair(wgt2/norm, sgm2));
    multi_gauss_.push_back(std::make_pair(wgt3/norm, sgm3));
    multi_gauss_.push_back(std::make_pair(wgt4/norm, sgm4));
    multi_gauss_.push_back(std::make_pair(wgt5/norm, sgm5));
    bound_.first  = std::min(std::min(std::min(std::min(sgm1, sgm2), sgm3), sgm4), sgm5);
    bound_.second = std::max(std::max(std::max(std::max(sgm1, sgm2), sgm3), sgm4), sgm5);
    robust_ = opt;
    
    // find max sigma
    //Double_t max_sgm = std::max_element(
    //            multi_gauss_.begin(), multi_gauss_.end(), 
    //            [](const std::pair<Double_t, Double_t>& gauss1, const std::pair<Double_t, Double_t>& gauss2) { return (gauss1.second < gauss2.second); }
    //        )->second;
}


Double_t MultiGauss::efft_sgm(Double_t r) const {
    if (multi_gauss_.size() == 0) return MGMath::ZERO;

    long double sigma = MGMath::ZERO;
    if (multi_gauss_.size() == 1) sigma = multi_gauss_.at(0).second;
    else {
        long double ttl_wgt = MGMath::ZERO;
        long double inv_nrm = MGMath::ZERO;
        for (auto&& gauss : multi_gauss_) {
            long double res = (r / gauss.second);
            long double nrm = (gauss.second / bound_.second);
            long double prb = (gauss.first / nrm) * std::exp(-MGMath::HALF * res * res);
            ttl_wgt += prb;
            inv_nrm += prb / (nrm * nrm);
        }
        
        long double sgm = bound_.second * ((MGNumc::Compare(ttl_wgt, LMTL_PROB_) <= 0 || MGNumc::EqualToZero(inv_nrm)) ? MGMath::ONE : std::sqrt(ttl_wgt / inv_nrm));
        if (!MGNumc::Valid(sgm) || MGNumc::Compare(sgm, bound_.second) > 0 || MGNumc::Compare(sgm) <= 0) sgm = bound_.second;
        sigma = sgm;
    }
   
    // Robust Method
    if (robust_ == Opt::ROBUST) {
        long double nrmr = std::fabs(r) / (ROBUST_ * sigma);
        if (MGNumc::Compare(nrmr, 1.0L) > 0)
            sigma *= std::sqrt(1.0L + std::log(nrmr * nrmr));
    }

    return static_cast<Double_t>(sigma);
}


Double_t MultiGauss::rndm() {
    if (multi_gauss_.size() == 0) return MGMath::ZERO;
    if (multi_gauss_.size() == 1) return (multi_gauss_.at(0).second * MGRndm::NormalGaussian());
    
    if (rand_func_) return rand_func_->GetRandom();
    else {
        std::string fmt = "([0]/[1])*TMath::Exp(-0.5*(x/[1])*(x/[1]))";
        for (int it = 1; it < multi_gauss_.size(); ++it) {
            int wgtID = 2 * it;
            int sgmID = 2 * it + 1;
            fmt += STR_FMT(" + ([%d]/[%d])*TMath::Exp(-0.5*(x/[%d])*(x/[%d]))", wgtID, sgmID, sgmID, sgmID);
        }
        rand_func_ = new TF1("rand_func", fmt.c_str(), -(MGMath::TWO*ROBUST_) * bound_.second, (MGMath::TWO*ROBUST_) * bound_.second);
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

        return rand_func_->GetRandom();
    }
    
    return MGMath::ZERO;
}


} // namesapce TrackSys


#endif // __TRACKLibs_Math_C__
