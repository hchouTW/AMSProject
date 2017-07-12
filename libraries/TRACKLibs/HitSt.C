#ifndef __TRACKLibs_HitSt_C__
#define __TRACKLibs_HitSt_C__


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


inline Double_t MultiGauss::efft_sgm(Double_t x) const {
    if (multi_gauss_.size() == 0) return MGMath::ZERO;
    if (multi_gauss_.size() == 1) return multi_gauss_.at(0).second;

    Double_t ttl_wgt     = MGMath::ZERO;
    Double_t inv_sqr_sgm = MGMath::ZERO;
    for (auto&& gauss : multi_gauss_) {
        Double_t sqr_sgm = (gauss.second * gauss.second);
        Double_t prob = (gauss.first / gauss.second) * std::exp(-MGMath::HALF * x * x / sqr_sgm);
        ttl_wgt     += prob;
        inv_sqr_sgm += prob / sqr_sgm;
    }
    Double_t sgm = std::sqrt(ttl_wgt / inv_sqr_sgm);
    if (!MGNumc::Valid(sgm)) sgm = MGMath::ZERO;
    return sgm;
}


void HitSt::print() const {
    std::string printStr;
    printStr += STR_FMT("================= HitSt ==================\n");
    printStr += STR_FMT("Id   (%d)\n", id_);
    printStr += STR_FMT("Side (%d %d %d)\n", side_(0), side_(0), side_(0));
    printStr += STR_FMT("Coo  (%11.6f %11.6f %11.6f)\n", coo_(0), coo_(1), coo_(2));
    printStr += STR_FMT("Err  (%11.6f %11.6f %11.6f)\n", err_x_.efft_sgm(), err_y_.efft_sgm(), err_z_.efft_sgm());
    printStr += STR_FMT("==========================================\n");
    COUT(printStr);
}


} // namesapce TrackSys


#endif // __TRACKLibs_HitSt_C__
