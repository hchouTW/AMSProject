#ifndef __TRACKLibs_Math_C__
#define __TRACKLibs_Math_C__


namespace TrackSys {
        

MultiGauss::MultiGauss(Double_t sgm) : rand_func_(nullptr)  {
    multi_gauss_.push_back(std::make_pair(MGMath::ONE, sgm));
    bound_.first  = sgm;
    bound_.second = sgm;
}


MultiGauss::MultiGauss(Double_t wgt1, Double_t sgm1, Double_t wgt2, Double_t sgm2) : rand_func_(nullptr)  {
    Double_t norm = wgt1 + wgt2;
    multi_gauss_.push_back(std::make_pair(wgt1/norm, sgm1));
    multi_gauss_.push_back(std::make_pair(wgt2/norm, sgm2));
    bound_.first  = std::min(sgm1, sgm2);
    bound_.second = std::max(sgm1, sgm2);
}


MultiGauss::MultiGauss(Double_t wgt1, Double_t sgm1, Double_t wgt2, Double_t sgm2, Double_t wgt3, Double_t sgm3) : rand_func_(nullptr)  {
    Double_t norm = wgt1 + wgt2 + wgt3;
    multi_gauss_.push_back(std::make_pair(wgt1/norm, sgm1));
    multi_gauss_.push_back(std::make_pair(wgt2/norm, sgm2));
    multi_gauss_.push_back(std::make_pair(wgt3/norm, sgm3));
    bound_.first  = std::min(std::min(sgm1, sgm2), sgm3);
    bound_.second = std::max(std::max(sgm1, sgm2), sgm3);
}


Double_t MultiGauss::efft_sgm(Double_t r) const {
    if (multi_gauss_.size() == 0) return MGMath::ZERO;
    if (multi_gauss_.size() == 1) return multi_gauss_.at(0).second;

    // find max sigma
    //Double_t max_sgm = std::max_element(
    //            multi_gauss_.begin(), multi_gauss_.end(), 
    //            [](const std::pair<Double_t, Double_t>& gauss1, const std::pair<Double_t, Double_t>& gauss2) { return (gauss1.second < gauss2.second); }
    //        )->second;

    Double_t ttl_wgt = MGMath::ZERO;
    Double_t inv_nrm = MGMath::ZERO;
    for (auto&& gauss : multi_gauss_) {
        Double_t res = (r / gauss.second);
        Double_t nrm = (gauss.second / bound_.second);
        Double_t prb = (gauss.first / nrm) * std::exp(-MGMath::HALF * res * res);
        ttl_wgt += prb;
        inv_nrm += prb / (nrm * nrm);
    }
    
    Double_t sgm = bound_.second * ((MGNumc::Compare(ttl_wgt, LMTL_PROB) <= 0 || MGNumc::EqualToZero(inv_nrm)) ? MGMath::ONE : std::sqrt(ttl_wgt / inv_nrm));
    if (!MGNumc::Valid(sgm)) sgm = bound_.second;
    
    // Robust estimator
    if (std::fabs(r / sgm) > ROBUST) sgm = (ROBUST * std::fabs(r));

    return sgm;
}
        

Double_t MultiGauss::rndm() {
    if (multi_gauss_.size() == 0) return 0.0;
    if (multi_gauss_.size() == 1) return (multi_gauss_.at(0).second * MGRndm::NormalGaussian());
    
    if (rand_func_) return rand_func_->GetRandom();
    else {
        if (multi_gauss_.size() == 2) {
            rand_func_ = new TF1("rand_func", "([0]/[1])*TMath::Exp(-0.5*x*x/[1]/[1]) + ([2]/[3])*TMath::Exp(-0.5*x*x/[3]/[3])", -ROBUST*bound_.second, ROBUST*bound_.second);
            rand_func_->SetNpx(100000);

        }
        if (multi_gauss_.size() == 3) {
            rand_func_ = new TF1("rand_func", "([0]/[1])*TMath::Exp(-0.5*x*x/[1]/[1]) + ([2]/[3])*TMath::Exp(-0.5*x*x/[3]/[3]) + ([4]/[5])*TMath::Exp(-0.5*x*x/[5]/[5])", -ROBUST*bound_.second, ROBUST*bound_.second);
            rand_func_->SetNpx(100000);
        }
        if (rand_func_) {
            short count = 0;
            for (auto&& gauss : multi_gauss_) {
                rand_func_->SetParameter(2*count+0, gauss.first);
                rand_func_->SetParameter(2*count+1, gauss.second);
                count++;
            }
            return rand_func_->GetRandom();
        }
    }
    
    return 0.0;
}


} // namesapce TrackSys


#endif // __TRACKLibs_Math_C__
