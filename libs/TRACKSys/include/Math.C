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

long double Robust::chisq(long double chi) const {
    if (Opt::OFF == opt_) return (chi * chi);
    long double ss    = (chi / threshold_) * (chi / threshold_);
    long double sqrss = ss * ss;
    long double alpha = std::log1p(sqrss) / sqrss;
    if (!Numc::Valid(alpha) || Numc::EqualToZero(alpha)) alpha = Numc::ONE<long double>;
    long double chisq = (chi * chi) * std::pow(alpha, Numc::HALF * rate_);
    if (!Numc::Valid(chisq)) chisq = Numc::ZERO<long double>;
    return chisq;
}

long double Robust::rescale(long double chi) const {
    if (Opt::OFF == opt_) return Numc::ONE<long double>;
    long double ss    = (chi / threshold_) * (chi / threshold_);
    long double sqrss = ss * ss;
    long double alpha = std::log1p(sqrss) / sqrss;
    if (!Numc::Valid(alpha) || Numc::EqualToZero(alpha)) alpha = Numc::ONE<long double>;
    long double div = std::pow(alpha, Numc::HALF * rate_) * 
                      ((Numc::ONE<long double> - rate_) + 
                       (Numc::ONE<long double>/alpha/(Numc::ONE<long double> + sqrss)));
    if (!Numc::Valid(div) || Numc::Compare(div) <= 0) div = Numc::ZERO<long double>;
    long double scl = std::sqrt(div);
    return scl;
}

} // namesapce TrackSys


namespace TrackSys {

TRandom* MultiGaus::rndm_gen_ = nullptr;

MultiGaus::MultiGaus(Robust robust, long double sgm) : MultiGaus()  {
    multi_gaus_.push_back(std::make_pair(Numc::ONE<long double>, sgm));
    bound_.first  = sgm;
    bound_.second = sgm;
    robust_ = robust;
}

MultiGaus::MultiGaus(Robust robust, long double wgt1, long double sgm1, long double wgt2, long double sgm2) : MultiGaus()  {
    long double norm = wgt1 + wgt2;
    multi_gaus_.push_back(std::make_pair(wgt1/norm, sgm1));
    multi_gaus_.push_back(std::make_pair(wgt2/norm, sgm2));
    bound_.first  = std::min(sgm1, sgm2);
    bound_.second = std::max(sgm1, sgm2);
    robust_ = robust;
}

MultiGaus::MultiGaus(Robust robust, long double wgt1, long double sgm1, long double wgt2, long double sgm2, long double wgt3, long double sgm3) : MultiGaus()  {
    long double norm = wgt1 + wgt2 + wgt3;
    multi_gaus_.push_back(std::make_pair(wgt1/norm, sgm1));
    multi_gaus_.push_back(std::make_pair(wgt2/norm, sgm2));
    multi_gaus_.push_back(std::make_pair(wgt3/norm, sgm3));
    bound_.first  = std::min(std::min(sgm1, sgm2), sgm3);
    bound_.second = std::max(std::max(sgm1, sgm2), sgm3);
    robust_ = robust;
}

MultiGaus::MultiGaus(Robust robust, long double wgt1, long double sgm1, long double wgt2, long double sgm2, long double wgt3, long double sgm3, long double wgt4, long double sgm4) : MultiGaus()  {
    long double norm = wgt1 + wgt2 + wgt3 + wgt4;
    multi_gaus_.push_back(std::make_pair(wgt1/norm, sgm1));
    multi_gaus_.push_back(std::make_pair(wgt2/norm, sgm2));
    multi_gaus_.push_back(std::make_pair(wgt3/norm, sgm3));
    multi_gaus_.push_back(std::make_pair(wgt4/norm, sgm4));
    bound_.first  = std::min(std::min(std::min(sgm1, sgm2), sgm3), sgm4);
    bound_.second = std::max(std::max(std::max(sgm1, sgm2), sgm3), sgm4);
    robust_ = robust;
}

MultiGaus::MultiGaus(Robust robust, long double wgt1, long double sgm1, long double wgt2, long double sgm2, long double wgt3, long double sgm3, long double wgt4, long double sgm4, long double wgt5, long double sgm5) : MultiGaus()  {
    long double norm = wgt1 + wgt2 + wgt3 + wgt4 + wgt5;
    multi_gaus_.push_back(std::make_pair(wgt1/norm, sgm1));
    multi_gaus_.push_back(std::make_pair(wgt2/norm, sgm2));
    multi_gaus_.push_back(std::make_pair(wgt3/norm, sgm3));
    multi_gaus_.push_back(std::make_pair(wgt4/norm, sgm4));
    multi_gaus_.push_back(std::make_pair(wgt5/norm, sgm5));
    bound_.first  = std::min(std::min(std::min(std::min(sgm1, sgm2), sgm3), sgm4), sgm5);
    bound_.second = std::max(std::max(std::max(std::max(sgm1, sgm2), sgm3), sgm4), sgm5);
    robust_ = robust;
    
    // find max sigma
    //long double max_sgm = std::max_element(
    //            multi_gaus_.begin(), multi_gaus_.end(), 
    //            [](const std::pair<long double, long double>& gauss1, const std::pair<long double, long double>& gauss2) { return (gauss1.second < gauss2.second); }
    //        )->second;
}

std::array<long double, 2> MultiGaus::minimizer(long double r) const {
    std::array<long double, 2> obj({Numc::ZERO<long double>, Numc::ONE<long double>});
    if      (multi_gaus_.size() == 0) return obj;
    else if (multi_gaus_.size() == 1) {
        long double nrm = (r / multi_gaus_.at(0).second);
        long double div = (Numc::ONE<long double> / multi_gaus_.at(0).second);
        obj = std::array<long double, 2>({ nrm, div }); 
    }
    else {
        // Note: inv_sgm_sqr = sum(prb * inv_sgm_sqr) / sum(prb)
        long double sum_prb = Numc::ZERO<long double>;
        long double inv_sgm = Numc::ZERO<long double>;
        for (auto&& gaus : multi_gaus_) {
            long double res = (r / gaus.second);
            long double sgm = (gaus.second / bound_.second);
            long double prb = gaus.first * std::exp(-static_cast<long double>(Numc::HALF) * res * res);
            sum_prb += prb;
            inv_sgm += prb / (sgm * sgm);
        }
        
        short sign         = Numc::Compare(r);
        long double eftsgm = bound_.second * ((Numc::Compare(sum_prb, LMTL_PROB) <= 0 || Numc::EqualToZero(inv_sgm)) ? Numc::ONE<long double> : std::sqrt(sum_prb / inv_sgm));
        long double nrm    = ((Numc::Compare(sum_prb, LMTL_PROB) <= 0) ? (r / eftsgm) : sign * std::sqrt(Numc::NEG<long double> * Numc::TWO<long double> * std::log(sum_prb)));
        long double div    = Numc::ONE<long double> / eftsgm; 
        obj = std::array<long double, 2>({ nrm, div });
    }

    if (!Numc::Valid(obj.at(0)) || !Numc::Valid(obj.at(1))) 
        obj = std::array<long double, 2>({ Numc::ZERO<long double>, Numc::ONE<long double> });

    if (Robust::Opt::ON == robust_.opt()) {
        long double scl = robust_.rescale(obj.at(0));
        obj.at(0) *= scl;
        obj.at(1) *= scl;
    }
    return obj;
}

long double MultiGaus::rndm() {
    if (multi_gaus_.size() == 0) return Numc::ZERO<long double>;
    if (multi_gaus_.size() == 1) return (multi_gaus_.at(0).second * Rndm::NormalGaussian());
    
    if (rand_func_) return rand_func_->GetRandom();
    else {
        std::string fmt = "[0]*TMath::Exp(-0.5*(x/[1])*(x/[1]))";
        for (unsigned int it = 1; it < multi_gaus_.size(); ++it) {
            int wgtID = 2 * it;
            int sgmID = 2 * it + 1;
            fmt += STR(" + [%d]*TMath::Exp(-0.5*(x/[%d])*(x/[%d]))", wgtID, sgmID, sgmID);
        }
        rand_func_ = new TF1("rand_func", fmt.c_str(), -Numc::TWO<long double> * bound_.second, Numc::TWO<long double> * bound_.second);
        rand_func_->SetNpx(NPX);
        
        int count = 0;
        for (auto&& gaus : multi_gaus_) {
            rand_func_->SetParameter(2*count+0, gaus.first);
            rand_func_->SetParameter(2*count+1, gaus.second);
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

LandauGaus::LandauGaus(Robust robust, long double kpa, long double mpv, long double sgm, long double fluc) : robust_(robust), kpa_(Numc::ZERO<long double>), mpv_(Numc::ZERO<long double>), sgm_(Numc::ONE<long double>), fluc_(Numc::ZERO<long double>), corr_(Numc::ONE<long double>) {
    if (Numc::Compare(sgm) <= 0) return;
    if      (Numc::Compare(kpa, Numc::ZERO<long double>) <= 0) kpa_ = Numc::ZERO<long double>;
    else if (Numc::Compare(kpa, Numc::ONE<long double>)  >= 0) kpa_ = Numc::ONE<long double>;
    else kpa_ = kpa;
    mpv_ = mpv;
    sgm_ = sgm;
    fluc_ = ((Numc::Compare(fluc) <= 0) ? Numc::ZERO<long double> : fluc);
    corr_ = Numc::ONE<long double> / std::sqrt(Numc::ONE<long double> + (fluc_/sgm_) * (fluc_/sgm_));
}

std::array<long double, 2> LandauGaus::minimizer(long double x) const {
    long double norm = ((x - mpv_) / sgm_);
    if (!Numc::Valid(norm))
        return std::array<long double, 2>({Numc::ZERO<long double>, Numc::ZERO<long double>});
   
    // Norm and Div
    long double nrmx = (corr_ * eval(norm));                   // norm x
    long double divx = (corr_ / sgm_) * std::sqrt(icov(norm)); // div x
    
    if (Robust::Opt::ON == robust_.opt()) {
        long double scl = robust_.rescale(nrmx);
        nrmx *= scl;
        divx *= scl;
    }
    
    if (!Numc::Valid(nrmx) || !Numc::Valid(divx)) {
        nrmx = Numc::ZERO<long double>;
        divx = Numc::ONE<long double>;
    }
    return std::array<long double, 2>({ nrmx, divx }); 
}

long double LandauGaus::eval(long double norm) const {
    short       sign   = Numc::Compare(norm);
    long double landau = (Numc::NEG<long double> * Numc::TWO<long double>) * std::log(TMath::Landau(norm * WIDTH_SCL + LANDAU0_X) / LANDAU0);
    long double gaus   = norm * norm;
    long double ldgaus = (Numc::ONE<long double> - kpa_) * landau + kpa_ * gaus;
    long double nrmx   = sign * std::sqrt(ldgaus);
    if (!Numc::Valid(nrmx)) nrmx = Numc::ZERO<long double>;
    return nrmx;
}
        
long double LandauGaus::icov(long double norm) const {
    std::array<long double, 8> ldps({ 
            1.71213e-01, 1.05335e+00, 
            1.22843e-01, 2.44705e-01, 
            2.95471e-01, 5.75783e-01, 
            1.61325e-02, 6.61252e-02 });
    long double ldicov = 
        (ldps.at(0) * std::exp(-ldps.at(1) * norm) +
         ldps.at(2) * std::exp(-ldps.at(3) * norm) + 
         ldps.at(4) * std::exp(-ldps.at(5) * norm) + 
         ldps.at(6) * std::exp(-ldps.at(7) * norm));
    long double landau = (Numc::ONE<long double> - kpa_) * ldicov;
    long double gaus   = kpa_;
    long double icov   = (gaus + landau);
    if (!Numc::Valid(icov) || Numc::Compare(icov) <= 0) icov = Numc::ONE<long double>;
    return icov;
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
