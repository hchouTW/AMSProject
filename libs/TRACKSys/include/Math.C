#ifndef __TRACKLibs_Math_C__
#define __TRACKLibs_Math_C__


namespace TrackSys {
namespace Numc { 

// Compare
template <class IntType, typename std::enable_if<std::is_integral<IntType>::value, int>::type>
inline short Compare(IntType a, IntType b) {
	if (a == b)     return  0;
	else if (a > b) return  1;
	else            return -1;
}

// Compare
template <class RealType, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type>
inline short Compare(RealType a, RealType b) {
	RealType diff = std::fabs(a - b);
    if (!std::isfinite(diff)) return 0;
	if (diff < std::numeric_limits<RealType>::epsilon() * 5.0e3) return  0;
	else if (a > b)                                              return  1;
	else                                                         return -1;
}

} // namesapce Numc
} // namesapce TrackSys


namespace TrackSys {
namespace Rndm {

// Uniform Distributions
template <class IntType, typename std::enable_if<std::is_integral<IntType>::value, int>::type>
std::function<IntType()> Uniform(IntType a, IntType b) {
	std::uniform_int_distribution<IntType> distribution(a, b);
	std::function<IntType()>&& rngfunc = std::bind(distribution, std::ref(rndmEngMT64));
	return rngfunc;
}

template <class RealType, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type>
std::function<RealType()> Uniform(RealType a, RealType b) {
	std::uniform_real_distribution<RealType> distribution(a, b);
	std::function<RealType()>&& rngfunc = std::bind(distribution, std::ref(rndmEngMT64));
	return rngfunc;
}


// Normal Distributions
template <class RealType, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type>
std::function<RealType()> Gamma(RealType alpha, RealType beta) {
	std::gamma_distribution<RealType> distribution(alpha, beta);
	std::function<RealType()>&& rngfunc = std::bind(distribution, std::ref(rndmEngMT64));
	return rngfunc;
}	


// Gamma Distributions
template <class RealType, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type>
std::function<RealType()> Normal(RealType mean, RealType stddev) {
	std::normal_distribution<RealType> distribution(mean, stddev);
	std::function<RealType()>&& rngfunc = std::bind(distribution, std::ref(rndmEngMT64));
	return rngfunc;
}
	
} // namesapce Rndm
} // namesapce TrackSys



namespace TrackSys {

MultiGauss::MultiGauss(Opt opt, long double sgm) : MultiGauss()  {
    multi_gauss_.push_back(std::make_pair(Numc::ONE<long double>, sgm));
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
    if (multi_gauss_.size() == 0) return Numc::ZERO<long double>;

    long double sigma = Numc::ONE<long double>;
    long double absr = std::abs(r);
    if (multi_gauss_.size() == 1) sigma = multi_gauss_.at(0).second;
    else {
        // Note: inv_sgm_sqr = sum(prb * inv_sgm_sqr) / sum(prb)
        long double ttl_wgt = Numc::ZERO<long double>;
        long double inv_nrm = Numc::ZERO<long double>;
        for (auto&& gauss : multi_gauss_) {
            long double res = (absr / gauss.second);
            long double nrm = (gauss.second / bound_.second);
            long double prb = (gauss.first / nrm) * std::exp(-static_cast<long double>(Numc::HALF) * res * res);
            ttl_wgt += prb;
            inv_nrm += prb / (nrm * nrm);
        }
        
        sigma = bound_.second * ((Numc::Compare(ttl_wgt, LMTL_PROB_) <= 0 || Numc::EqualToZero(inv_nrm)) ? Numc::ONE<long double> : std::sqrt(ttl_wgt / inv_nrm));
        if (!Numc::Valid(sigma) || Numc::Compare(sigma, bound_.second) > 0 || Numc::Compare(sigma) <= 0) sigma = bound_.second;
    }
   
    // Robust Method (Modify-Cauchy)
    if (robust_ == Opt::ROBUST) {
        long double sftnrmr = (absr / sigma) - ROBUST_SGM_;
        if (Numc::Compare(sftnrmr) > 0) {
            long double cauchy = (sftnrmr / std::sqrt(std::log1p(sftnrmr * sftnrmr)));
            long double modify_cauchy = ((!Numc::Valid(cauchy) || Numc::Compare(cauchy, Numc::ONE<long double>) < 0) ? Numc::ONE<long double> : std::sqrt(cauchy));
            if (Numc::Valid(modify_cauchy)) sigma *= modify_cauchy;
        }
    }

    return sigma;
}

long double MultiGauss::rndm() {
    if (multi_gauss_.size() == 0) return Numc::ZERO<long double>;
    if (multi_gauss_.size() == 1) return (multi_gauss_.at(0).second * Rndm::NormalGaussian());
    
    if (rand_func_) return rand_func_->GetRandom();
    else {
        std::string fmt = "([0]/[1])*TMath::Exp(-0.5*(x/[1])*(x/[1]))";
        for (unsigned int it = 1; it < multi_gauss_.size(); ++it) {
            int wgtID = 2 * it;
            int sgmID = 2 * it + 1;
            fmt += STR(" + ([%d]/[%d])*TMath::Exp(-0.5*(x/[%d])*(x/[%d]))", wgtID, sgmID, sgmID, sgmID);
        }
        rand_func_ = new TF1("rand_func", fmt.c_str(), -(Numc::TWO<long double>*ROBUST_SGM_) * bound_.second, (Numc::TWO<long double>*ROBUST_SGM_) * bound_.second);
        rand_func_->SetNpx(NPX_);
        
        int count = 0;
        for (auto&& gauss : multi_gauss_) {
            rand_func_->SetParameter(2*count+0, gauss.first);
            rand_func_->SetParameter(2*count+1, gauss.second);
            count++;
        }
        
        if (rndm_gen_ == nullptr) {
            gRandom->SetSeed(0);
            rndm_gen_ = gRandom;
        }

        return static_cast<long double>(rand_func_->GetRandom());
    }
    return Numc::ZERO<long double>;
}


} // namesapce TrackSys


#endif // __TRACKLibs_Math_C__
