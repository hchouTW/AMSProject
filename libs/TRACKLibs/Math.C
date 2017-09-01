#ifndef __TRACKLibs_Math_C__
#define __TRACKLibs_Math_C__


namespace TrackSys {
        

MultiGauss::MultiGauss(Double_t sgm) {
    multi_gauss_.push_back(std::make_pair(MGMath::ONE, sgm));
}


MultiGauss::MultiGauss(Double_t wgt1, Double_t sgm1, Double_t wgt2, Double_t sgm2) {
    Double_t norm = wgt1 + wgt2;
    multi_gauss_.push_back(std::make_pair(wgt1/norm, sgm1));
    multi_gauss_.push_back(std::make_pair(wgt2/norm, sgm2));
}


MultiGauss::MultiGauss(Double_t wgt1, Double_t sgm1, Double_t wgt2, Double_t sgm2, Double_t wgt3, Double_t sgm3) {
    Double_t norm = wgt1 + wgt2 + wgt3;
    multi_gauss_.push_back(std::make_pair(wgt1/norm, sgm1));
    multi_gauss_.push_back(std::make_pair(wgt2/norm, sgm2));
    multi_gauss_.push_back(std::make_pair(wgt3/norm, sgm3));
}


inline Double_t MultiGauss::efft_sgm(Double_t r) const {
    if (multi_gauss_.size() == 0) return MGMath::ZERO;
    if (multi_gauss_.size() == 1) return multi_gauss_.at(0).second;

    Double_t max_sgm = std::max_element(
                multi_gauss_.begin(), multi_gauss_.end(), 
                [](const std::pair<Double_t, Double_t>& gauss1, const std::pair<Double_t, Double_t>& gauss2) { return (gauss1.second < gauss2.second); }
            )->second;

    Double_t ttl_wgt = MGMath::ZERO;
    Double_t inv_nrm = MGMath::ZERO;
    for (auto&& gauss : multi_gauss_) {
        Double_t res = (r / gauss.second);
        Double_t nrm = (gauss.second / max_sgm);
        Double_t prb = (gauss.first / nrm) * std::exp(-MGMath::HALF * res * res);
        ttl_wgt += prb;
        inv_nrm += prb / (nrm * nrm);
    }
    
    Double_t sgm = max_sgm * ((MGNumc::Compare(ttl_wgt, LMTL_PROB) <= 0 || MGNumc::EqualToZero(inv_nrm)) ? MGMath::ONE : std::sqrt(ttl_wgt / inv_nrm));
    if (!MGNumc::Valid(sgm)) sgm = max_sgm;
    
    // Robust estimator
    if (std::fabs(r / sgm) > ROBUST) sgm = (ROBUST * std::fabs(r));

    return sgm;
}


} // namesapce TrackSys


#endif // __TRACKLibs_Math_C__
