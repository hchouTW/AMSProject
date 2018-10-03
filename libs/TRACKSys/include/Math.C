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


std::vector<long double> TrackSys::LandauNumc::data;
long double TrackSys::LandauNumc::EvalLn(long double x) {
    // Build data set
    long double newx = (LANDAU0_X + x * WIDTH_SCL);
    if (data.size() == 0) {
        for (long it = 0; it <= NSET; ++it)
            data.push_back( std::log(TMath::Landau(LANDAU0_X + (BOUNDL+it*STEP) * WIDTH_SCL) / LANDAU0) );
    }
    
    long double landln = Numc::ZERO<long double>;
    if (x <= BOUNDL || x >= BOUNDU) landln = std::log(TMath::Landau(newx) / LANDAU0);
    else {
        long iset = std::lrint( ((x - BOUNDL) / STEP) );
        landln = data.at(iset);
    }
    return landln;
}


namespace TrackSys {
        
std::array<long double, 4> Robust::minimizer(long double chi) const {
    std::array<long double, 4> mini { Numc::ONE<long double>, Numc::ONE<long double>, Numc::ONE<long double>, Numc::ONE<long double> };
    if (Opt::OFF == opt_) return mini;
     
    const long double ZERO  = Numc::ZERO<long double>;
    const long double ONE   = Numc::ONE<long double>;
    const long double TWO   = Numc::TWO<long double>;
    const long double THREE = Numc::THREE<long double>;
    const long double FOUR  = Numc::FOUR<long double>;
    const long double HALF  = Numc::ONE<long double> / Numc::TWO<long double>;
   
    if (Opt::ON == robust_opt_) {
        long double sh       = (chi / thres_) * (chi / thres_);
        long double sqrsh    = sh * sh;
        long double onesqrsh = (ONE + sqrsh);
        long double alpha    = std::log1p(sqrsh) / sqrsh;
        if (!Numc::Valid(alpha) || Numc::EqualToZero(sqrsh)) alpha = ONE;
        long double invalpha = (Numc::EqualToZero(alpha) ? ZERO : (ONE / alpha));

        long double crchi = std::pow(alpha, (HALF / FOUR) * rate_);
        if (!Numc::Valid(crchi) || Numc::Compare(crchi) <= 0) crchi = ONE;

        long double div1st = std::pow(alpha, rate_ / FOUR) * (
                             (ONE - HALF * rate_) + 
                             (HALF * rate_ * invalpha / onesqrsh)
                             );
        if (!Numc::Valid(div1st) || Numc::Compare(div1st) < 0) div1st = ONE;

        long double div2ndS = (rate_ * sqrsh / FOUR) * std::pow(alpha, ONE + rate_ / FOUR) * (
                              ((rate_ - FOUR) / (std::pow(alpha, THREE) * sqrsh * onesqrsh * onesqrsh)) +
                              ((rate_ - TWO) * invalpha / sqrsh) -
                              (TWO * ((rate_ - ONE) * sqrsh + (rate_ - THREE)) * (invalpha * invalpha) / (sqrsh * onesqrsh * onesqrsh))
                              );
        if (!Numc::Valid(div2ndS) || Numc::Compare(div2ndS) >= 0) div2ndS = ZERO;

        long double corr = std::sqrt(ONE + TWO * div2ndS / div1st);
        if (!Numc::Valid(corr) || Numc::Compare(corr) < 0) corr = ONE;

        long double sqrtdiv = std::sqrt(div1st);
        long double crnorm  = sqrtdiv / corr;
        long double crjacb  = sqrtdiv * corr;
        if (!Numc::Valid(crnorm)) crnorm = ONE;
        if (!Numc::Valid(crjacb)) crjacb = ONE;

        mini.at(0) = crchi;
        mini.at(1) = crnorm;
        mini.at(2) = crjacb;
    }

    if (Opt::ON == ghost_opt_) {
        long double sh    = (chi / ghost_thres_) * (chi / ghost_thres_);
        long double sixsh = std::pow(sh, Numc::SIX<long double>);
        long double alpha = std::log1p(sixsh) / sixsh;
        if (!Numc::Valid(alpha) || Numc::EqualToZero(sixsh)) alpha = ONE;
        alpha = std::pow(alpha, FOUR);
        
        long double crwgt = alpha * alpha;
        if (!Numc::Valid(crwgt)) crwgt = ONE;

        mini.at(0) *= alpha;
        mini.at(1) *= alpha;
        mini.at(2) *= alpha;
        mini.at(3)  = crwgt;
    }

    return mini;
}

} // namesapce TrackSys


namespace TrackSys {

TRandom* MultiGaus::rndm_gen_ = nullptr;

MultiGaus::MultiGaus(Robust robust, long double sgm) : MultiGaus()  {
    multi_gaus_.push_back(std::make_pair(Numc::ONE<long double>, sgm));
    bound_.first  = sgm;
    bound_.second = sgm;
    robust_ = robust;
    eftsgm_ = find_eftsgm();
}

MultiGaus::MultiGaus(Robust robust, long double wgt1, long double sgm1, long double wgt2, long double sgm2) : MultiGaus()  {
    long double norm = wgt1 + wgt2;
    multi_gaus_.push_back(std::make_pair(wgt1/norm, sgm1));
    multi_gaus_.push_back(std::make_pair(wgt2/norm, sgm2));
    bound_.first  = std::min(sgm1, sgm2);
    bound_.second = std::max(sgm1, sgm2);
    robust_ = robust;
    eftsgm_ = find_eftsgm();
}

MultiGaus::MultiGaus(Robust robust, long double wgt1, long double sgm1, long double wgt2, long double sgm2, long double wgt3, long double sgm3) : MultiGaus()  {
    long double norm = wgt1 + wgt2 + wgt3;
    multi_gaus_.push_back(std::make_pair(wgt1/norm, sgm1));
    multi_gaus_.push_back(std::make_pair(wgt2/norm, sgm2));
    multi_gaus_.push_back(std::make_pair(wgt3/norm, sgm3));
    bound_.first  = std::min(std::min(sgm1, sgm2), sgm3);
    bound_.second = std::max(std::max(sgm1, sgm2), sgm3);
    robust_ = robust;
    eftsgm_ = find_eftsgm();
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
    eftsgm_ = find_eftsgm();
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
    eftsgm_ = find_eftsgm();
}

MultiGaus::MultiGaus(Robust robust, long double wgt1, long double sgm1, long double wgt2, long double sgm2, long double wgt3, long double sgm3, long double wgt4, long double sgm4, long double wgt5, long double sgm5, long double wgt6, long double sgm6) : MultiGaus()  {
    long double norm = wgt1 + wgt2 + wgt3 + wgt4 + wgt5 + wgt6;
    multi_gaus_.push_back(std::make_pair(wgt1/norm, sgm1));
    multi_gaus_.push_back(std::make_pair(wgt2/norm, sgm2));
    multi_gaus_.push_back(std::make_pair(wgt3/norm, sgm3));
    multi_gaus_.push_back(std::make_pair(wgt4/norm, sgm4));
    multi_gaus_.push_back(std::make_pair(wgt5/norm, sgm5));
    multi_gaus_.push_back(std::make_pair(wgt6/norm, sgm6));
    bound_.first  = std::min(std::min(std::min(std::min(std::min(sgm1, sgm2), sgm3), sgm4), sgm5), sgm6);
    bound_.second = std::max(std::max(std::max(std::max(std::max(sgm1, sgm2), sgm3), sgm4), sgm5), sgm6);
    robust_ = robust;
    eftsgm_ = find_eftsgm();
}

MultiGaus::MultiGaus(Robust robust, long double wgt1, long double sgm1, long double wgt2, long double sgm2, long double wgt3, long double sgm3, long double wgt4, long double sgm4, long double wgt5, long double sgm5, long double wgt6, long double sgm6, long double wgt7, long double sgm7) : MultiGaus()  {
    long double norm = wgt1 + wgt2 + wgt3 + wgt4 + wgt5 + wgt6 + wgt7;
    multi_gaus_.push_back(std::make_pair(wgt1/norm, sgm1));
    multi_gaus_.push_back(std::make_pair(wgt2/norm, sgm2));
    multi_gaus_.push_back(std::make_pair(wgt3/norm, sgm3));
    multi_gaus_.push_back(std::make_pair(wgt4/norm, sgm4));
    multi_gaus_.push_back(std::make_pair(wgt5/norm, sgm5));
    multi_gaus_.push_back(std::make_pair(wgt6/norm, sgm6));
    multi_gaus_.push_back(std::make_pair(wgt7/norm, sgm7));
    bound_.first  = std::min(std::min(std::min(std::min(std::min(std::min(sgm1, sgm2), sgm3), sgm4), sgm5), sgm6), sgm7);
    bound_.second = std::max(std::max(std::max(std::max(std::max(std::max(sgm1, sgm2), sgm3), sgm4), sgm5), sgm6), sgm7);
    robust_ = robust;
    eftsgm_ = find_eftsgm();
}
        
MultiGaus::MultiGaus(Robust robust, long double wgt1, long double sgm1, long double wgt2, long double sgm2, long double wgt3, long double sgm3, long double wgt4, long double sgm4, long double wgt5, long double sgm5, long double wgt6, long double sgm6, long double wgt7, long double sgm7, long double wgt8, long double sgm8) : MultiGaus()  {
    long double norm = wgt1 + wgt2 + wgt3 + wgt4 + wgt5 + wgt6 + wgt7 + wgt8;
    multi_gaus_.push_back(std::make_pair(wgt1/norm, sgm1));
    multi_gaus_.push_back(std::make_pair(wgt2/norm, sgm2));
    multi_gaus_.push_back(std::make_pair(wgt3/norm, sgm3));
    multi_gaus_.push_back(std::make_pair(wgt4/norm, sgm4));
    multi_gaus_.push_back(std::make_pair(wgt5/norm, sgm5));
    multi_gaus_.push_back(std::make_pair(wgt6/norm, sgm6));
    multi_gaus_.push_back(std::make_pair(wgt7/norm, sgm7));
    multi_gaus_.push_back(std::make_pair(wgt8/norm, sgm8));
    bound_.first  = std::min(std::min(std::min(std::min(std::min(std::min(std::min(sgm1, sgm2), sgm3), sgm4), sgm5), sgm6), sgm7), sgm8);
    bound_.second = std::max(std::max(std::max(std::max(std::max(std::max(std::max(sgm1, sgm2), sgm3), sgm4), sgm5), sgm6), sgm7), sgm8);
    robust_ = robust;
    eftsgm_ = find_eftsgm();
}
        
MultiGaus::MultiGaus(Robust robust, long double wgt1, long double sgm1, long double wgt2, long double sgm2, long double wgt3, long double sgm3, long double wgt4, long double sgm4, long double wgt5, long double sgm5, long double wgt6, long double sgm6, long double wgt7, long double sgm7, long double wgt8, long double sgm8, long double wgt9, long double sgm9) : MultiGaus()  {
    long double norm = wgt1 + wgt2 + wgt3 + wgt4 + wgt5 + wgt6 + wgt7 + wgt8 + wgt9;
    multi_gaus_.push_back(std::make_pair(wgt1/norm, sgm1));
    multi_gaus_.push_back(std::make_pair(wgt2/norm, sgm2));
    multi_gaus_.push_back(std::make_pair(wgt3/norm, sgm3));
    multi_gaus_.push_back(std::make_pair(wgt4/norm, sgm4));
    multi_gaus_.push_back(std::make_pair(wgt5/norm, sgm5));
    multi_gaus_.push_back(std::make_pair(wgt6/norm, sgm6));
    multi_gaus_.push_back(std::make_pair(wgt7/norm, sgm7));
    multi_gaus_.push_back(std::make_pair(wgt8/norm, sgm8));
    multi_gaus_.push_back(std::make_pair(wgt9/norm, sgm9));
    bound_.first  = std::min(std::min(std::min(std::min(std::min(std::min(std::min(std::min(sgm1, sgm2), sgm3), sgm4), sgm5), sgm6), sgm7), sgm8), sgm9);
    bound_.second = std::max(std::max(std::max(std::max(std::max(std::max(std::max(std::max(sgm1, sgm2), sgm3), sgm4), sgm5), sgm6), sgm7), sgm8), sgm9);
    robust_ = robust;
    eftsgm_ = find_eftsgm();
}
        

long double MultiGaus::chi(long double r) const {
    long double chival = Numc::ZERO<long double>;
    if      (multi_gaus_.size() == 0) return Numc::ZERO<long double>;
    else if (multi_gaus_.size() == 1) return (r / multi_gaus_.at(0).second);
    else {
        long double sum_prb = Numc::ZERO<long double>;
        for (auto&& gaus : multi_gaus_) {
            long double res = (r / gaus.second);
            long double prb = gaus.first * std::exp(-static_cast<long double>(Numc::HALF) * res * res);
            sum_prb += prb;
        }
        long double nrm = Numc::Compare(r) * std::sqrt(Numc::NEG<long double> * Numc::TWO<long double> * std::log(sum_prb));
        if (!Numc::Valid(nrm)) nrm = Numc::ZERO<long double>;
        chival = nrm;
    }
    return chival;
}


std::array<long double, 4> MultiGaus::minimizer(long double r) const {
    std::array<long double, 4> obj { Numc::ZERO<long double>, Numc::ZERO<long double>, Numc::ONE<long double>, Numc::ONE<long double> };
    if      (multi_gaus_.size() == 0) return obj;
    else if (multi_gaus_.size() == 1) {
        long double nrm = (r / multi_gaus_.at(0).second);
        long double div = (Numc::ONE<long double> / multi_gaus_.at(0).second);
        obj = std::array<long double, 4>({ nrm, nrm, div, Numc::ONE<long double> }); 
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
        long double eftsgm = bound_.second * ((Numc::Compare(sum_prb, LMTL_PROB) <= 0 || Numc::Compare(inv_sgm) < 0) ? Numc::ONE<long double> : std::sqrt(sum_prb / inv_sgm));
        long double nrm    = ((Numc::Compare(sum_prb, LMTL_PROB) <= 0) ? (r / eftsgm) : sign * std::sqrt(Numc::NEG<long double> * Numc::TWO<long double> * std::log(sum_prb)));
        long double div    = Numc::ONE<long double> / eftsgm; 
        obj = std::array<long double, 4>({ nrm, nrm, div, Numc::ONE<long double> });
    }

    if (!Numc::Valid(obj.at(0)) || !Numc::Valid(obj.at(1)) || !Numc::Valid(obj.at(2))) 
        obj = std::array<long double, 4>({ Numc::ZERO<long double>, Numc::ZERO<long double>, Numc::ONE<long double>, Numc::ONE<long double> });

    if (Robust::Opt::ON == robust_.opt()) {
        std::array<long double, 4>&& rbmini = robust_.minimizer(obj.at(0));
        obj.at(0) *= rbmini.at(0);
        obj.at(1) *= rbmini.at(1);
        obj.at(2) *= rbmini.at(2);
        obj.at(3)  = rbmini.at(3);
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
        

long double MultiGaus::find_eftsgm() const {
    long double sgm = Numc::ONE<long double>;
    if      (multi_gaus_.size() == 0) return Numc::ONE<long double>;
    else if (multi_gaus_.size() == 1) return multi_gaus_.at(0).second;
    else {
        std::vector<long double> chivec;
        for (auto&& gaus : multi_gaus_) chivec.push_back( chi(gaus.second) );
        
        long double chil = Numc::ZERO<>;
        long double sgml = Numc::ZERO<>;
        for (int it = 0; it < chivec.size(); ++it)
            if (chivec.at(it) < Numc::ONE<long double>)
                { chil = chivec.at(it); sgml = multi_gaus_.at(it).second; }

        long double chiu = Numc::ZERO<>;
        long double sgmu = Numc::ZERO<>;
        for (int it = chivec.size()-1; it >= 0; --it)
            if (chivec.at(it) > Numc::ONE<long double>)
                { chiu = chivec.at(it); sgmu = multi_gaus_.at(it).second; }

        const int max_iter = 50;
        const long double lmt = 1.0e-04;
        if (std::fabs(chiu - chil) < lmt) sgm = 0.5L * (sgml + sgmu);
        else {
            for (int iter = 1; iter <= max_iter; ++iter) {
                if (std::fabs(chiu - chil) < lmt) { sgm = 0.5L * (sgml + sgmu); break; }
                long double wgtl = std::exp(-Numc::HALF * (chil * chil - Numc::ONE<long double>));
                long double wgtu = std::exp(-Numc::HALF * (chiu * chiu - Numc::ONE<long double>));
                long double sgmm = (wgtl * sgmu + wgtu * sgml) / (wgtl + wgtu);
                long double chim = chi(sgmm);

                sgm = sgmm;
                if (std::fabs(chim - Numc::ONE<long double>) < lmt) break;
                if (chim < Numc::ONE<long double>) { chil = chim; sgml = sgmm; }
                else                               { chiu = chim; sgmu = sgmm; }
            }
        }
    }
    if (!Numc::Valid(sgm)) sgm = multi_gaus_.at(0).second;
    return sgm;
}


} // namesapce TrackSys


namespace TrackSys {

LandauGaus::LandauGaus(Robust robust, long double kpa, long double mpv, long double sgm, long double mode, long double fluc) : robust_(robust), isfluc_(false), kpa_(Numc::ZERO<long double>), mpv_(Numc::ZERO<long double>), sgm_(Numc::ONE<long double>), mode_(Numc::ZERO<long double>), fluc_(Numc::ZERO<long double>), shft_(Numc::ZERO<long double>) {
    if (Numc::Compare(sgm) <= 0) return;
    if      (Numc::Compare(kpa, Numc::ZERO<long double>) <= 0) kpa_ = Numc::ZERO<long double>;
    else if (Numc::Compare(kpa, Numc::ONE<long double>)  >= 0) kpa_ = Numc::ONE<long double>;
    else kpa_ = kpa;
    mpv_ = mpv;
    sgm_ = sgm;

    isfluc_ = (Numc::Compare(fluc) > 0);
    if (isfluc_) {
        mode_ = mode;
        fluc_ = fluc;
        shft_ = ((mode_ - mpv_) / sgm_);
    }
    else mode_ = mpv_;
}

std::array<long double, 4> LandauGaus::minimizer(long double x) const {
    std::array<long double, 4> obj { Numc::ZERO<long double>, Numc::ZERO<long double>, Numc::ONE<long double>, Numc::ONE<long double> };
    long double norm = ((x - mpv_) / sgm_);
    if (!Numc::Valid(norm)) return obj;
   
    // Norm and Div
    long double nrmx = Numc::ZERO<long double>; // norm x
    long double divx = Numc::ONE<long double>;  // div x
    if (!isfluc_) {
        nrmx = eval_norm(norm);           
        divx = std::sqrt(eval_icov(norm));
    }
    else {
        std::array<long double, 2>&& conv = eval_conv(norm);
        nrmx = conv.at(0);
        divx = std::sqrt(conv.at(1));
    }
    if (!Numc::Valid(nrmx)) nrmx = Numc::ZERO<long double>;
    if (!Numc::Valid(divx)) divx = Numc::ONE<long double>;
    obj = std::array<long double, 4>({ nrmx, nrmx, divx, Numc::ONE<long double> });

    if (Robust::Opt::ON == robust_.opt()) {
        std::array<long double, 4>&& rbmini = robust_.minimizer(obj.at(0));
        obj.at(0) *= rbmini.at(0);
        obj.at(1) *= rbmini.at(1);
        obj.at(2) *= rbmini.at(2);
        obj.at(3)  = rbmini.at(3);
    }
    return obj; 
}

long double LandauGaus::eval_norm(long double norm) const {
    short       sign   = Numc::Compare(norm);
    long double landau = (-Numc::TWO<long double>) * LandauNumc::EvalLn(norm);
    long double ldgaus = (kpa_ * norm * norm) + (Numc::ONE<long double> - kpa_) * landau;
    long double nrmx   = sign * std::sqrt(ldgaus);
    if (!Numc::Valid(nrmx)) nrmx = Numc::ZERO<long double>;
    return nrmx;
}

long double LandauGaus::eval_icov(long double norm) const {
    long double ldinvc =
        (LAND_CONV[0] * std::exp(-LAND_CONV[1] * norm) +
         LAND_CONV[2] * std::exp(-LAND_CONV[3] * norm) + 
         LAND_CONV[4] * std::exp(-LAND_CONV[5] * norm) + 
         LAND_CONV[6] * std::exp(-LAND_CONV[7] * norm));
    long double icov = (kpa_ + (Numc::ONE<long double> - kpa_) * ldinvc);
    return icov;
}

std::array<long double, 2> LandauGaus::eval_conv(long double norm) const { // (nrm, icov)
    std::array<long double, 2> mini { Numc::ZERO<long double>, Numc::ONE<long double> };
    
    std::array<long double, GAUS_CONV_N>&& prob = convprob(norm);
    long double sumprob = std::accumulate(prob.begin(), prob.end(), Numc::ZERO<long double>);
    
    std::array<long double, GAUS_CONV_N>&& modeprob = convprob(shft_);
    long double summodeprob = std::accumulate(modeprob.begin(), modeprob.end(), Numc::ZERO<long double>);
    long double rate = (sumprob / summodeprob);
    
    long double nrmx = Numc::ZERO<long double>;
    if (Numc::Compare(rate, Numc::ONE<long double>) < 0) {
        nrmx = Numc::Compare(norm, shft_) * std::sqrt(-Numc::TWO<long double> * std::log(rate));
        if (!Numc::Valid(nrmx)) nrmx = Numc::ZERO<long double>;
    }

    long double icov = Numc::ZERO<long double>;
    for (int it = 0; it < GAUS_CONV_N; ++it) {
        long double newx = norm - GAUS_CONV_X[it] * (fluc_ / sgm_);
        long double invc = eval_icov(newx);
        icov += invc * (prob[it] / sumprob);
    }
    if (!Numc::Valid(icov) || Numc::Compare(icov) < 0) icov = Numc::ONE<long double>;

    mini.at(0) = nrmx;
    mini.at(1) = icov;

    return mini;
}

std::array<long double, LandauGaus::GAUS_CONV_N> LandauGaus::convprob(long double norm) const {
    std::array<long double, GAUS_CONV_N> prob; prob.fill(Numc::ZERO<long double>);
    for (int it = 0; it < GAUS_CONV_N; ++it) {
        long double newx = norm - GAUS_CONV_X[it] * (fluc_ / sgm_);
        long double land = (-Numc::TWO<long double>) * LandauNumc::EvalLn(newx);
        long double ldgs = kpa_ * (newx * newx) + (Numc::ONE<long double> - kpa_) * land;
        long double elem = std::exp(-Numc::HALF * ldgs) * GAUS_CONV_P[it];
        prob[it] = elem;
    }
    return prob;
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
