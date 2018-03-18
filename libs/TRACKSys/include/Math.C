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

    long double sigma = Numc::ONE<long double>;
    long double absr = std::abs(r);
    if (multi_gaus_.size() == 1) sigma = multi_gaus_.at(0).second;
    else {
        // Note: inv_sgm_sqr = sum(prb * inv_sgm_sqr) / sum(prb)
        long double ttl_wgt = Numc::ZERO<long double>;
        long double inv_nrm = Numc::ZERO<long double>;
        for (auto&& gauss : multi_gaus_) {
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
        rand_func_ = new TF1("rand_func", fmt.c_str(), -(Numc::TWO<long double>*ROBUST_SGM_) * bound_.second, (Numc::TWO<long double>*ROBUST_SGM_) * bound_.second);
        rand_func_->SetNpx(NPX_);
        
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

std::array<long double, 2> LandauGaus::operator() (long double x) const {
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
        long double sftnrmr = absnrmx - ROBUST_SGM_;
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
    long double landau = (Numc::NEG<long double> * Numc::TWO<long double>) * std::log(TMath::Landau(norm) / LANDAU0_);
    long double gaus   = norm * norm;
    long double ldgaus = (Numc::ONE<long double> - kpa_) * landau + kpa_ * gaus;
    long double nrmx   = sign * std::sqrt(ldgaus);
    if (!Numc::Valid(nrmx)) nrmx = Numc::ZERO<long double>;
    return nrmx;
}

long double LandauGaus::div(long double norm) const {
    long double normxlw = norm - DELTA_;
    long double normxup = norm + DELTA_;
    long double div = Numc::NEG<long double> * Numc::HALF * ((eval(normxup) - eval(normxlw)) / DELTA_);
    if (!Numc::Valid(div)) div = Numc::ZERO<long double>;
    return div;
}

} // namesapce TrackSys


namespace TrackSys {

std::array<long double, 2> IonEloss::eval(long double x, long double eta) const {
    if (Numc::Compare(x) <= 0 || Numc::EqualToZero(eta))
        return std::array<long double, 2>({Numc::ZERO<long double>, Numc::ZERO<long double>});

    long double abseta = std::fabs(eta);
    long double ibsqr  = (Numc::ONE<long double> + abseta * abseta);
    
    // Robust Method (Reduce important index in high energy)
    long double lngm = std::log(std::sqrt(ibsqr) / abseta); // log(gamma)
    if (!Numc::Valid(lngm) || Numc::Compare(lngm) <= 0) lngm = Numc::ZERO<long double>;
    long double robust = Numc::ONE<long double> / (Numc::ONE<long double> + lngm * lngm);
  
    // PDF parameters
    long double kpa    = eval_kpa(abseta, ibsqr); 
    long double mpv    = eval_mpv(abseta, ibsqr); 
    long double sgm    = eval_sgm(abseta, ibsqr); 
    long double divmpv = eval_divmpv(abseta, ibsqr); 
  
    // Landau-Gaus with noise fluctuation 
    LandauGaus ldgaus(LandauGaus::Opt::ROBUST, kpa, mpv, sgm, fluc_);
    std::array<long double, 2>&& lg_par = ldgaus(x);
    
    long double res = robust * lg_par.at(0);          // res normx
    long double div = robust * lg_par.at(1) * divmpv; // div f/x * div x/eta
    if (!Numc::Valid(res) || !Numc::Valid(div)) { 
        res = Numc::ZERO<long double>;
        div = Numc::ZERO<long double>;
    }
    return std::array<long double, 2>({res, div});
}
        
long double IonEloss::eval_kpa(long double eta, long double ibsqr) const {
    long double kpa = kpa_.at(0) * std::pow(ibsqr,  kpa_.at(1)) * 
        (Numc::ONE<long double> +
         kpa_.at(2) * std::pow(eta, kpa_.at(3)) -
         std::log(kpa_.at(4) + std::pow(eta, kpa_.at(5)))
        );
    if (!Numc::Valid(kpa)) kpa = Numc::ZERO<long double>;
    else {
        if (Numc::Compare(kpa, Numc::ZERO<long double>) <= 0) kpa = Numc::ZERO<long double>;
        if (Numc::Compare(kpa, Numc::ONE<long double>)  >= 0) kpa = Numc::ONE<long double>;
    }
    return kpa;
}

long double IonEloss::eval_mpv(long double eta, long double ibsqr) const {
    long double mpv = mpv_.at(0) * std::pow(ibsqr,  mpv_.at(2)) * 
        (mpv_.at(1) - 
         std::pow(ibsqr, -mpv_.at(2)) - 
         std::log(mpv_.at(3) + std::pow(eta, mpv_.at(4)))
        );
    if (!Numc::Valid(mpv)) mpv = Numc::ZERO<long double>;
    return mpv;
}

long double IonEloss::eval_sgm(long double eta, long double ibsqr) const {
    long double sgm = sgm_.at(0) * std::pow(ibsqr,  sgm_.at(2)) * 
        (sgm_.at(1) - 
         std::pow(ibsqr, -sgm_.at(2)) - 
         std::log(sgm_.at(3) + std::pow(eta, sgm_.at(4)))
        );
    if (!Numc::Valid(sgm)) sgm = Numc::ZERO<long double>;
    return sgm;
}

long double IonEloss::eval_divmpv(long double eta, long double ibsqr) const {
    long double divbta = mpv_.at(2) * std::pow(ibsqr, mpv_.at(2)-Numc::ONE<long double>) * (Numc::TWO<long double> * eta);
    long double divlog = mpv_.at(4) * std::pow(eta, mpv_.at(4)-Numc::ONE<long double>) / (mpv_.at(3) + std::pow(eta, mpv_.at(4)));

    long double termA  = mpv_.at(1) * divbta;
    long double termB  = divbta * std::log(mpv_.at(3) + std::pow(eta, mpv_.at(4)));
    long double termC  = std::pow(ibsqr, mpv_.at(2)) * divlog;
    long double divmpv = mpv_.at(0) * (termA - termB - termC);

    if (!Numc::Valid(divmpv)) divmpv = Numc::ZERO<long double>;
    return divmpv;
}

} // namesapce TrackSys


#endif // __TRACKLibs_Math_C__
