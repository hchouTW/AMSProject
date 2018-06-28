#ifndef __TRACKLibs_Math_C__
#define __TRACKLibs_Math_C__


#include "Sys.h"
#include "Math.h"


namespace TrackSys {
namespace Numc { 
} // namesapce Numc
} // namesapce TrackSys


namespace TrackSys {
namespace Rndm {
} // namesapce Rndm
} // namesapce TrackSys



namespace TrackSys {

TRandom* MultiGaus::rndm_gen_ = nullptr;

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

std::array<long double, 2> LandauGaus::minimizer(long double x) const {
    long double norm = ((x - mpv_) / sgm_);
    if (!Numc::Valid(norm))
        return std::array<long double, 2>({Numc::ZERO<long double>, Numc::ZERO<long double>});
   
    // Norm and Div
    long double nrmx = eval(norm); // norm x
    long double divx = div(norm) * (Numc::NEG<long double> / sgm_);  // div x with sgm scale
    
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
    long double div = Numc::HALF * ((eval(normxup) - eval(normxlw)) / DELTA);
    if (!Numc::Valid(div)) div = Numc::ZERO<long double>;
    return div;
}

} // namesapce TrackSys


namespace TrackSys {

std::array<long double, 3> LgGeFunc::minimizer(long double x) const {
    long double rat    = get_ratio(x);
    long double wgt_lg = std::sqrt(Numc::ONE<long double> - rat);
    std::array<long double, 2> lg = lg_minimizer(x, wgt_lg);
    return std::array<long double, 3>({ lg.at(0), lg.at(1), rat }); 
}

long double LgGeFunc::get_ratio(long double x) const {
    long double lg  = get_lg(x) * (Numc::ONE<long double> - ratio_);
    long double ge  = get_ge(x) * (ratio_);
    long double sum = lg + ge;
    long double rat = (ge / sum);
    if (!Numc::Valid(rat) || Numc::Compare(rat) <= 0) rat = Numc::ZERO<long double>;
    return rat;
}

long double LgGeFunc::get_lg(long double x) const {
    long double norm   = ((x - lg_m_) / lg_s_);
    long double landau = (Numc::ONE<long double> - lg_k_) * std::log(TMath::Landau(norm) / LANDAU0);
    long double gaus   = lg_k_ * Numc::NEG<long double> * Numc::HALF * (norm * norm);
    long double lg     = (Numc::INV_SQRT_TWO * Numc::INV_SQRT_PI / lg_s_) * std::exp(landau + gaus);
    if (!Numc::Valid(lg) || Numc::Compare(lg) <= 0) lg = Numc::ZERO<long double>;
    return lg;
}

long double LgGeFunc::get_ge(long double x) const {
    long double gamma = (std::pow(ge_b_, ge_a_) / std::tgamma(ge_a_)) * std::pow(x, ge_a_ - Numc::ONE<long double>) * std::exp(Numc::NEG<long double> * ge_b_ * x);
    long double erf   = Numc::HALF * (Numc::ONE<long double> + std::erf((x - ge_m_) / ge_s_));
    long double ge    = gamma * erf;
    if (!Numc::Valid(ge) || Numc::Compare(ge) <= 0) ge = Numc::ZERO<long double>;
    return ge;
}

std::array<long double, 2> LgGeFunc::lg_minimizer(long double x, long double wgt) const {
    // Norm and Div
    long double norm = ((x - lg_m_) / lg_s_);
    long double nrmx = lg_eval(norm); // norm x
    long double divx = lg_div(norm) * (Numc::NEG<long double> / lg_s_);  // div x with sgm scale
    
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

    // testcode
    // Weight of (ION / (ION + TR));
    //if (Numc::Compare(wgt) >= 0) {
    //    nrmx *= wgt;
    //    divx *= wgt;
    //}
    
    if (!Numc::Valid(nrmx) || !Numc::Valid(divx)) {
        nrmx = Numc::ZERO<long double>;
        divx = Numc::ZERO<long double>;
    }
    return std::array<long double, 2>({ nrmx, divx }); 
}

long double LgGeFunc::lg_eval(long double norm) const {
    short       sign   = Numc::Compare(norm);
    long double landau = (Numc::NEG<long double> * Numc::TWO<long double>) * std::log(TMath::Landau(norm) / LANDAU0);
    long double gaus   = norm * norm;
    long double ldgaus = (Numc::ONE<long double> - lg_k_) * landau + lg_k_ * gaus;
    long double nrmx   = sign * std::sqrt(ldgaus);
    if (!Numc::Valid(nrmx)) nrmx = Numc::ZERO<long double>;
    return nrmx;
}

long double LgGeFunc::lg_div(long double norm) const {
    long double normxlw = norm - DELTA;
    long double normxup = norm + DELTA;
    long double div = Numc::HALF * ((lg_eval(normxup) - lg_eval(normxlw)) / DELTA);
    if (!Numc::Valid(div)) div = Numc::ZERO<long double>;
    return div;
}

} // namesapce TrackSys


#endif // __TRACKLibs_Math_C__
