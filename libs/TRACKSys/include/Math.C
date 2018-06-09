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

// Chi-Square Distributions
template <class RealType, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type>
inline std::function<RealType()> ChiSquare(RealType k) {
    std::chi_squared_distribution<RealType> distribution(k);
	std::function<RealType()>&& rngfunc = std::bind(distribution, std::ref(rndmEngMT64));
	return rngfunc;
}
	
} // namesapce Rndm
} // namesapce TrackSys



namespace TrackSys {

MultiGaus::MultiGaus(Opt opt, long double sgm) : MultiGaus()  {
    multi_gaus_.push_back(std::make_pair(Numc::ONE<long double>, sgm));
    bound_.first  = sgm;
    bound_.second = sgm;
    robust_ = opt;
}

MultiGaus::MultiGaus(Opt opt, long double wgt1, long double sgm1, long double wgt2, long double sgm2) : MultiGaus()  {
    long double norm = wgt1 + wgt2;
    multi_gaus_.push_back(std::make_pair(wgt1/norm, sgm1));
    multi_gaus_.push_back(std::make_pair(wgt2/norm, sgm2));
    bound_.first  = std::min(sgm1, sgm2);
    bound_.second = std::max(sgm1, sgm2);
    robust_ = opt;
}

MultiGaus::MultiGaus(Opt opt, long double wgt1, long double sgm1, long double wgt2, long double sgm2, long double wgt3, long double sgm3) : MultiGaus()  {
    long double norm = wgt1 + wgt2 + wgt3;
    multi_gaus_.push_back(std::make_pair(wgt1/norm, sgm1));
    multi_gaus_.push_back(std::make_pair(wgt2/norm, sgm2));
    multi_gaus_.push_back(std::make_pair(wgt3/norm, sgm3));
    bound_.first  = std::min(std::min(sgm1, sgm2), sgm3);
    bound_.second = std::max(std::max(sgm1, sgm2), sgm3);
    robust_ = opt;
}

MultiGaus::MultiGaus(Opt opt, long double wgt1, long double sgm1, long double wgt2, long double sgm2, long double wgt3, long double sgm3, long double wgt4, long double sgm4) : MultiGaus()  {
    long double norm = wgt1 + wgt2 + wgt3 + wgt4;
    multi_gaus_.push_back(std::make_pair(wgt1/norm, sgm1));
    multi_gaus_.push_back(std::make_pair(wgt2/norm, sgm2));
    multi_gaus_.push_back(std::make_pair(wgt3/norm, sgm3));
    multi_gaus_.push_back(std::make_pair(wgt4/norm, sgm4));
    bound_.first  = std::min(std::min(std::min(sgm1, sgm2), sgm3), sgm4);
    bound_.second = std::max(std::max(std::max(sgm1, sgm2), sgm3), sgm4);
    robust_ = opt;
}

MultiGaus::MultiGaus(Opt opt, long double wgt1, long double sgm1, long double wgt2, long double sgm2, long double wgt3, long double sgm3, long double wgt4, long double sgm4, long double wgt5, long double sgm5) : MultiGaus()  {
    long double norm = wgt1 + wgt2 + wgt3 + wgt4 + wgt5;
    multi_gaus_.push_back(std::make_pair(wgt1/norm, sgm1));
    multi_gaus_.push_back(std::make_pair(wgt2/norm, sgm2));
    multi_gaus_.push_back(std::make_pair(wgt3/norm, sgm3));
    multi_gaus_.push_back(std::make_pair(wgt4/norm, sgm4));
    multi_gaus_.push_back(std::make_pair(wgt5/norm, sgm5));
    bound_.first  = std::min(std::min(std::min(std::min(sgm1, sgm2), sgm3), sgm4), sgm5);
    bound_.second = std::max(std::max(std::max(std::max(sgm1, sgm2), sgm3), sgm4), sgm5);
    robust_ = opt;
    
    // find max sigma
    //long double max_sgm = std::max_element(
    //            multi_gaus_.begin(), multi_gaus_.end(), 
    //            [](const std::pair<long double, long double>& gauss1, const std::pair<long double, long double>& gauss2) { return (gauss1.second < gauss2.second); }
    //        )->second;
}

long double MultiGaus::efft_sgm(long double r) const {
    if (multi_gaus_.size() == 0) return Numc::ZERO<long double>;

    long double sgm = Numc::ONE<long double>;
    if (multi_gaus_.size() == 1) sgm = multi_gaus_.at(0).second;
    else {
        // Note: inv_sgm_sqr = sum(prb * inv_sgm_sqr) / sum(prb)
        long double ttl_wgt = Numc::ZERO<long double>;
        long double inv_nrm = Numc::ZERO<long double>;
        for (auto&& gauss : multi_gaus_) {
            long double res = (r / gauss.second);
            long double nrm = (gauss.second / bound_.second);
            long double prb = (gauss.first / nrm) * std::exp(-static_cast<long double>(Numc::HALF) * res * res);
            ttl_wgt += prb;
            inv_nrm += prb / (nrm * nrm);
        }
        
        sgm = bound_.second * ((Numc::Compare(ttl_wgt, LMTL_PROB) <= 0 || Numc::EqualToZero(inv_nrm)) ? Numc::ONE<long double> : std::sqrt(ttl_wgt / inv_nrm));
        if (!Numc::Valid(sgm) || Numc::Compare(sgm, bound_.second) > 0 || Numc::Compare(sgm) <= 0) sgm = bound_.second;
    }
   
    // Robust Method (Modify-Cauchy)
    long double sigma = RobustSgm(r, sgm, robust_); 
    return sigma;
}

long double MultiGaus::rndm() {
    if (multi_gaus_.size() == 0) return Numc::ZERO<long double>;
    if (multi_gaus_.size() == 1) return (multi_gaus_.at(0).second * Rndm::NormalGaussian());
    
    if (rand_func_) return rand_func_->GetRandom();
    else {
        std::string fmt = "([0]/[1])*TMath::Exp(-0.5*(x/[1])*(x/[1]))";
        for (unsigned int it = 1; it < multi_gaus_.size(); ++it) {
            int wgtID = 2 * it;
            int sgmID = 2 * it + 1;
            fmt += STR(" + ([%d]/[%d])*TMath::Exp(-0.5*(x/[%d])*(x/[%d]))", wgtID, sgmID, sgmID, sgmID);
        }
        rand_func_ = new TF1("rand_func", fmt.c_str(), -(Numc::TWO<long double>*ROBUST_SGM) * bound_.second, (Numc::TWO<long double>*ROBUST_SGM) * bound_.second);
        rand_func_->SetNpx(NPX);
        
        int count = 0;
        for (auto&& gauss : multi_gaus_) {
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

long double MultiGaus::RobustSgm(long double res, long double sgm, Opt opt) {
    long double absr  = std::fabs(res);
    long double sigma = ((Numc::Compare(sgm) > 0) ? sgm : Numc::ONE<long double>);
    if (opt == Opt::NOROBUST) return sigma;

    // Robust Method (Modify-Cauchy)
    long double sftnrmr = (absr / sigma) - ROBUST_SGM;
    if (Numc::Compare(sftnrmr) > 0) {
        long double cauchy = (sftnrmr / std::sqrt(std::log1p(sftnrmr * sftnrmr)));
        long double modify_cauchy = ((!Numc::Valid(cauchy) || Numc::Compare(cauchy, Numc::ONE<long double>) < 0) ? Numc::ONE<long double> : std::sqrt(cauchy));
        if (Numc::Valid(modify_cauchy)) sigma *= modify_cauchy;
    }

    if (!Numc::Valid(sigma)) sigma = Numc::ONE<long double>;
    return sigma;
}

} // namesapce TrackSys


namespace TrackSys {

LandauGaus::LandauGaus(Opt opt, long double kpa, long double mpv, long double sgm, long double fluc) : robust_(Opt::NOROBUST), kpa_(Numc::ZERO<long double>), mpv_(Numc::ZERO<long double>), sgm_(Numc::ONE<long double>), fluc_(Numc::ZERO<long double>) {
    if (Numc::Compare(sgm) <= 0) return;
    robust_ = opt;
    if      (Numc::Compare(kpa, Numc::ZERO<long double>) <= 0) kpa_ = Numc::ZERO<long double>;
    else if (Numc::Compare(kpa, Numc::ONE<long double>)  >= 0) kpa_ = Numc::ONE<long double>;
    else kpa_ = kpa;
    mpv_ = mpv;
    sgm_ = sgm;
    fluc_ = ((Numc::Compare(fluc) <= 0) ? Numc::ZERO<long double> : fluc);
}

long double LandauGaus::operator() (long double x, long double wgt) const {
    long double norm = ((x - mpv_) / sgm_);
    if (!Numc::Valid(norm)) return Numc::ZERO<long double>;
    if (Numc::Compare(wgt) <= 0) return Numc::ZERO<long double>;

    long double landau = TMath::Landau(norm) / LANDAU0;
    long double gaus   = Numc::NEG<long double> * Numc::HALF * norm * norm;
    long double expcom = (Numc::ONE<long double> - kpa_) * landau + kpa_ * gaus;
    long double lg     = wgt * std::exp(lg);
    
    if (!Numc::Valid(lg) || Numc::Compare(lg) <= 0) lg = Numc::ZERO<long double>;
    return lg;
}

std::array<long double, 2> LandauGaus::minimizer(long double x) const {
    long double norm = ((x - mpv_) / sgm_);
    if (!Numc::Valid(norm))
        return std::array<long double, 2>({Numc::ZERO<long double>, Numc::ZERO<long double>});
   
    // Norm and Div
    long double nrmx = eval(norm); // norm x
    long double divx = div(norm) / sgm_;  // div x with sgm scale
    
    // Noise fluctuation
    if (!Numc::EqualToZero(fluc_)) {
        long double nsfluc = Numc::ONE<long double> / std::sqrt(Numc::ONE<long double> + (fluc_/sgm_) * (fluc_/sgm_));
        nrmx *= nsfluc; 
        divx *= nsfluc; 
    }
    
    // Robust Method (Modify-Cauchy)
    if (robust_ == Opt::ROBUST) {
        long double absnrmx = std::fabs(nrmx);
        long double sftnrmr = absnrmx - ROBUST_SGM;
        if (Numc::Compare(sftnrmr) > 0) {
            long double cauchy = (sftnrmr / std::sqrt(std::log1p(sftnrmr * sftnrmr)));
            long double modify_cauchy = ((!Numc::Valid(cauchy) || Numc::Compare(cauchy, Numc::ONE<long double>) < 0) ? Numc::ONE<long double> : std::sqrt(cauchy));
            if (Numc::Valid(modify_cauchy)) {
                nrmx /= modify_cauchy;
                divx /= modify_cauchy;
            }
        }
    }
    
    if (!Numc::Valid(nrmx) || !Numc::Valid(divx)) {
        nrmx = Numc::ZERO<long double>;
        divx = Numc::ZERO<long double>;
    }
    return std::array<long double, 2>({ nrmx, divx }); 
}

long double LandauGaus::eval(long double norm) const {
    short       sign   = Numc::Compare(norm);
    long double landau = (Numc::NEG<long double> * Numc::TWO<long double>) * std::log(TMath::Landau(norm) / LANDAU0);
    long double gaus   = norm * norm;
    long double ldgaus = (Numc::ONE<long double> - kpa_) * landau + kpa_ * gaus;
    long double nrmx   = sign * std::sqrt(ldgaus);
    if (!Numc::Valid(nrmx)) nrmx = Numc::ZERO<long double>;
    return nrmx;
}

long double LandauGaus::div(long double norm) const {
    long double normxlw = norm - DELTA;
    long double normxup = norm + DELTA;
    long double div = Numc::NEG<long double> * Numc::HALF * ((eval(normxup) - eval(normxlw)) / DELTA);
    if (!Numc::Valid(div)) div = Numc::ZERO<long double>;
    return div;
}

} // namesapce TrackSys



namespace TrackSys {

ErfGamma::ErfGamma(long double alpha, long double beta, long double eftm, long double efts) : alpha_(Numc::ONE<long double>), beta_(Numc::ONE<long double>), eftm_(Numc::ZERO<long double>), efts_(Numc::ONE<long double>) {
    if (Numc::Compare(alpha, Numc::ONE<long double>) <= 0) return;
    if (Numc::Compare(beta) <= 0) return;
    if (Numc::Compare(eftm) <= 0) return;
    if (Numc::Compare(efts) <= 0) return;
    alpha_  = alpha;
    beta_   = beta;
    eftm_   = eftm;
    efts_   = efts;
}

long double ErfGamma::operator() (long double x, long double wgt) const {
    if (Numc::Compare(wgt) <= 0) return Numc::ZERO<long double>;
    
    long double gm    = std::pow(x, alpha_) * std::exp(Numc::NEG<long double> * beta_ * x);
    long double efft  = (std::erf((x - eftm_) / efts_) + Numc::ONE<long double>);
    long double erfgm = wgt * efft * gm;
    
    if (!Numc::Valid(erfgm) || Numc::Compare(erfgm) <= 0) erfgm = Numc::ZERO<long double>;
    return erfgm;
}

} // namesapce TrackSys


#endif // __TRACKLibs_Math_C__
