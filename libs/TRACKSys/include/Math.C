#ifndef __TRACKLibs_Math_C__
#define __TRACKLibs_Math_C__


#include "Sys.h"
#include "Math.h"


namespace TrackSys {
namespace Numc { 


Double_t NormQuality(Double_t nchi, Short_t ndof) {
    if (Numc::Compare(nchi) < 0 || ndof <= Numc::ZERO<Short_t>) return Numc::ZERO<>;
    Double_t chi = nchi * static_cast<Double_t>(ndof);
    if (Numc::EqualToZero(chi)) return Numc::ZERO<>;
    if (ndof <= Numc::TWO<Short_t>) return std::sqrt(nchi);
    Double_t qmin  = static_cast<Double_t>(ndof - Numc::TWO<Short_t>);
    Double_t sign  = static_cast<Double_t>(Numc::Compare(chi - qmin));
    Double_t qfunc = (chi - qmin) - qmin * std::log(chi / qmin);
    if (!Numc::Valid(qfunc)) return Numc::ZERO<>;
    Double_t xfunc = sign * std::sqrt(qfunc / static_cast<Double_t>(ndof));
    if (Numc::Valid(xfunc)) return xfunc;
    return Numc::ZERO<>;
}


Double_t NormQuality(Double_t nchi, Double_t ndof) {
    if (Numc::Compare(nchi) < 0 || Numc::Compare(ndof) <= 0) return Numc::ZERO<>;
    Double_t chi = nchi * ndof;
    if (Numc::EqualToZero(chi)) return Numc::ZERO<>;
    if (Numc::Compare(ndof, Numc::TWO<>) <= 0) return std::sqrt(nchi);
    Double_t qmin  = (ndof - Numc::TWO<>);
    Double_t sign  = static_cast<Double_t>(Numc::Compare(chi - qmin));
    Double_t qfunc = (chi - qmin) - qmin * std::log(chi / qmin);
    if (!Numc::Valid(qfunc)) return Numc::ZERO<>;
    Double_t xfunc = sign * std::sqrt(qfunc / ndof);
    if (Numc::Valid(xfunc)) return xfunc;
    return Numc::ZERO<>;
}
        

} // namesapce Numc
} // namesapce TrackSys


namespace TrackSys {
namespace Rndm {
} // namesapce Rndm
} // namesapce TrackSys


namespace TrackSys {
        
std::array<long double, 3> Robust::minimizer(long double chi) const {
    std::array<long double, 3> mini { Numc::ONE<long double>, Numc::ONE<long double>, Numc::ONE<long double> };
    if (Opt::OFF == opt_) return mini;
     
    const long double ZERO  = Numc::ZERO<long double>;
    const long double ONE   = Numc::ONE<long double>;
    const long double TWO   = Numc::TWO<long double>;
    const long double THREE = Numc::THREE<long double>;
    const long double FOUR  = Numc::FOUR<long double>;

    if (Opt::ON == robust_.opt) {
        long double abschi = std::fabs(chi / robust_.thres);
        long double rate   = robust_.rate * (Numc::EqualToZero(abschi) ? ZERO : Numc::HALF * (ONE + std::erf(robust_.thres * std::log(abschi))));
        if (!Numc::Valid(rate)) rate = ZERO;

        long double sh       = abschi * abschi;
        long double sqrsh    = sh * sh;
        long double onesqrsh = (ONE + sqrsh);
        long double alpha    = std::log1p(sqrsh) / sqrsh;
        if (!Numc::Valid(alpha) || Numc::EqualToZero(sqrsh)) alpha = ONE;
        long double invalpha = (Numc::EqualToZero(alpha) ? ZERO : (ONE / alpha));

        long double crchi = std::pow(alpha, (Numc::HALF / FOUR) * rate);
        if (!Numc::Valid(crchi) || Numc::Compare(crchi) <= 0) crchi = ONE;

        long double div1st = std::pow(alpha, rate / FOUR) * (
                             (ONE - Numc::HALF * rate) + 
                             (Numc::HALF * rate * invalpha / onesqrsh)
                             );
        if (!Numc::Valid(div1st) || Numc::Compare(div1st) < 0) div1st = ONE;

        long double div2ndS = (rate * sqrsh / FOUR) * std::pow(alpha, ONE + rate / FOUR) * (
                              ((rate - FOUR) / (std::pow(alpha, THREE) * sqrsh * onesqrsh * onesqrsh)) +
                              ((rate - TWO) * invalpha / sqrsh) -
                              (TWO * ((rate - ONE) * sqrsh + (rate - THREE)) * (invalpha * invalpha) / (sqrsh * onesqrsh * onesqrsh))
                              );
        if (!Numc::Valid(div2ndS) || Numc::Compare(div2ndS) >= 0) div2ndS = ZERO;

        long double corr = std::sqrt(ONE + TWO * div2ndS / div1st);
        if (!Numc::Valid(corr) || Numc::Compare(corr) < 0) corr = ONE;

        long double sqrtdiv = std::sqrt(div1st);
        long double crnorm  = sqrtdiv / corr;
        long double crjacb  = sqrtdiv * corr;
        if (!Numc::Valid(crnorm)) crnorm = ONE;
        if (!Numc::Valid(crjacb)) crjacb = ONE;

        mini.at(1) = crnorm;
        mini.at(2) = crjacb;
    }
    return mini;
}

} // namesapce TrackSys


namespace TrackSys {

TRandom* MultiGaus::rndm_gen_ = nullptr;
        
long double MultiGaus::Func(long double x, long double men, const std::vector<std::array<long double, 2>>& group) {
    if (group.size() == 0) return Numc::ZERO<long double>;

    long double sumw = Numc::ZERO<long double>;
    std::vector<std::array<long double, 2>> regroup;
    for (auto&& elem : group) {
        if (!Numc::Valid(elem[0]) || Numc::Compare(elem[0]) <= 0) continue;
        if (!Numc::Valid(elem[1]) || Numc::Compare(elem[1]) <= 0) continue;
        long double wgt = elem[0];
        long double res = (x - men) / elem[1];
        regroup.push_back(std::array<long double, 2>({ wgt, res }));
        sumw += wgt;
    }
    if (regroup.size() == 0) return Numc::ZERO<long double>;
    if (!Numc::Valid(sumw)) return Numc::ZERO<long double>;
    for (auto&& elem : regroup) elem[0] /= sumw;

    long double value = Numc::ZERO<long double>;
    for (auto&& elem : regroup) value += elem[0] * std::exp(-Numc::HALF * elem[1] * elem[1]);
    if (!Numc::Valid(value)) value = Numc::ZERO<long double>;
    return value;
}

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
            long double prb = gaus.first * std::exp(-Numc::HALF * res * res);
            sum_prb += prb;
        }
        long double nrm = Numc::Compare(r) * std::sqrt(Numc::NEG<long double> * Numc::TWO<long double> * std::log(sum_prb));
        if (!Numc::Valid(nrm)) nrm = Numc::ZERO<long double>;
        chival = nrm;
    }
    return chival;
}


std::array<long double, 3> MultiGaus::minimizer(long double r) const {
    std::array<long double, 3> obj { Numc::ZERO<long double>, Numc::ZERO<long double>, Numc::ZERO<long double> };
    if      (multi_gaus_.size() == 0) return obj;
    else if (multi_gaus_.size() == 1) {
        long double nrm = (r / multi_gaus_.at(0).second);
        long double div = (Numc::ONE<long double> / multi_gaus_.at(0).second);
        obj = std::array<long double, 3>({ nrm, nrm, div }); 
    }
    else {
        // Note: inv_sgm_sqr = sum(prb * inv_sgm_sqr) / sum(prb)
        long double sum_prb  = Numc::ZERO<long double>;
        long double sum_isqr = Numc::ZERO<long double>;
        for (auto&& gaus : multi_gaus_) {
            long double res = (r / gaus.second);
            long double sgm = (gaus.second / bound_.second);
            long double prb = gaus.first * std::exp(-Numc::HALF * res * res);
            sum_prb  += prb;
            sum_isqr += prb / (sgm * sgm);
        }
        short       sign    = Numc::Compare(r);
        long double sgmisqr = ((Numc::Compare(sum_prb, LMTL_PROB) <= 0 || Numc::Compare(sum_isqr) < 0) ? Numc::ONE<long double> : (sum_isqr / sum_prb));
        long double hessian = (sgmisqr / (bound_.second * bound_.second));

        long double nrm = ((Numc::Compare(sum_prb, LMTL_PROB) <= 0) ? (r * std::sqrt(hessian)) : sign * std::sqrt(Numc::NEG<long double> * Numc::TWO<long double> * std::log(sum_prb)));
        long double crr = ((sign == 0 || Numc::EqualToZero(nrm)) ? (Numc::ONE<long double> / std::sqrt(hessian)) : std::fabs(r / nrm));
        long double div = crr * hessian;

        obj = std::array<long double, 3>({ nrm, nrm, div });
    }

    if (!Numc::Valid(obj.at(0)) || !Numc::Valid(obj.at(1)) || !Numc::Valid(obj.at(2))) 
        obj = std::array<long double, 3>({ Numc::ZERO<long double>, Numc::ZERO<long double>, Numc::ZERO<long double> });

    if (Robust::Opt::ON == robust_.opt()) {
        std::array<long double, 3>&& rbmini = robust_.minimizer(obj.at(0));
        obj.at(0) *= rbmini.at(0);
        obj.at(1) *= rbmini.at(1);
        obj.at(2) *= rbmini.at(2);
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

long double LandauGaus::Func(long double x, long double kpa, long double mpv, long double sgm, long double fluc) {
    if (Numc::Compare(sgm) <= 0 || Numc::Compare(kpa) < 0 || Numc::Compare(kpa, Numc::ONE<long double>) > 0) return Numc::ZERO<long double>;
    bool is_fluc = (Numc::Compare(fluc) > 0);
    
    long double norm = (x - mpv) / sgm;
    if (!Numc::Valid(norm)) return Numc::ZERO<long double>;
    
    long double value = Numc::ZERO<long double>;
    if (!is_fluc) {
        long double gs = -Numc::HALF * norm * norm;
        long double ld = -LandauNumc::NegLn(norm);
        long double ldgs = kpa * gs + (Numc::ONE<long double> - kpa) * ld;
        value = std::exp(ldgs);
    }
    else {
        for (Int_t it = 0; it <= CONV_N; ++it) {
            long double ix = norm + (fluc / sgm) * CONV_X[it];
            long double gs = -Numc::HALF * ix * ix;
            long double ld = -LandauNumc::NegLn(ix);
            long double ldgs = kpa * gs + (Numc::ONE<long double> - kpa) * ld;
            value += std::exp(ldgs) * CONV_P[it];
        }
    }

    if (!Numc::Valid(value)) value = Numc::ZERO<long double>;
    return value;
}

LandauGaus::LandauGaus(Robust robust, long double kpa, long double mpv, long double sgm, long double mod, long double fluc) : robust_(robust), isfluc_(false), kpa_(Numc::ZERO<long double>), mpv_(Numc::ZERO<long double>), sgm_(Numc::ONE<long double>), mod_(Numc::ZERO<long double>), fluc_(Numc::ZERO<long double>), shft_(Numc::ZERO<long double>), nrmfluc_(Numc::ZERO<long double>) {
    if (Numc::Compare(sgm) <= 0) return;
    if      (Numc::Compare(kpa, Numc::ZERO<long double>) <= 0) kpa_ = Numc::ZERO<long double>;
    else if (Numc::Compare(kpa, Numc::ONE<long double>)  >= 0) kpa_ = Numc::ONE<long double>;
    else kpa_ = kpa;
    mpv_ = mpv;
    sgm_ = sgm;

    isfluc_ = (Numc::Compare(fluc) > 0);
    if (isfluc_) {
        mod_ = mod;
        fluc_ = fluc;
        shft_ = ((mod_ - mpv_) / sgm_);
        nrmfluc_ = (fluc_ / sgm_);
    }
    else mod_ = mpv_;
}

std::array<long double, 3> LandauGaus::minimizer(long double x) const {
    std::array<long double, 3> obj { Numc::ZERO<long double>, Numc::ZERO<long double>, Numc::ZERO<long double> };
    long double norm = ((x - mpv_) / sgm_);
    if (!Numc::Valid(norm)) return obj;
   
    // Norm and Div
    long double nrmx = Numc::ZERO<long double>; // norm x
    long double divx = Numc::ZERO<long double>;  // div x
    if (!isfluc_) {
        std::array<long double, 5>&& val = eval(norm);
        nrmx = val.at(0); 
        divx = val.at(1);
    }
    else {
        std::array<long double, 2>&& conv = eval_conv(norm);
        nrmx = conv.at(0);
        divx = conv.at(1);
    }
    if (!Numc::Valid(nrmx)) nrmx = Numc::ZERO<long double>;
    if (!Numc::Valid(divx)) divx = Numc::ZERO<long double>;
    obj = std::array<long double, 3>({ nrmx, nrmx, divx });

    if (Robust::Opt::ON == robust_.opt()) {
        std::array<long double, 3>&& rbmini = robust_.minimizer(obj.at(0));
        obj.at(0) *= rbmini.at(0);
        obj.at(1) *= rbmini.at(1);
        obj.at(2) *= rbmini.at(2);
    }
    return obj; 
}

std::array<long double, 5> LandauGaus::eval(long double norm) const { // (nrm, jacb, prb, dev, hes)
    short       sign   = Numc::Compare(norm);
    long double ldgaus = (kpa_ * norm * norm) + Numc::TWO<long double> * (Numc::ONE<long double> - kpa_) * LandauNumc::NegLn(norm);
    long double nrmx   = sign * std::sqrt(ldgaus);
    if (!Numc::Valid(nrmx)) nrmx = Numc::ZERO<long double>;

    long double dev = (kpa_ * norm + (Numc::ONE<long double> - kpa_) * LandauNumc::NegLnDev(norm));
    long double hes = (kpa_ + (Numc::ONE<long double> - kpa_) * LandauNumc::NegLnHes(norm));
    if (!Numc::Valid(dev) || !Numc::Valid(hes)) {
        dev = Numc::ZERO<long double>;
        hes = Numc::ZERO<long double>;
    }

    long double jacb = (Numc::ONE<long double> / nrmx) * dev;
    if (Numc::EqualToZero(norm) || Numc::EqualToZero(nrmx) || !Numc::Valid(jacb))
        jacb = std::sqrt(hes);

    long double prb = std::exp(-Numc::HALF * ldgaus);
    if (!Numc::Valid(nrmx) || !Numc::Valid(jacb) || !Numc::Valid(prb)) {
        nrmx = Numc::ZERO<long double>;
        jacb = Numc::ZERO<long double>;
        prb  = Numc::ZERO<long double>;
        dev  = Numc::ZERO<long double>;
        hes  = Numc::ZERO<long double>;
    }
    return std::array<long double, 5>({ nrmx, jacb, prb, dev, hes });
}

std::array<long double, 2> LandauGaus::eval_conv(long double norm) const { // (nrm, jacb)
    std::array<long double, 2> mini { Numc::ZERO<long double>, Numc::ZERO<long double> };

    std::array<long double, 2>&& mod = conv(shft_);
    std::array<long double, 2>&& val = conv(norm);

    short       sign = Numc::Compare(norm, shft_);
    long double nrm  = sign * std::sqrt(-Numc::TWO<long double> * std::log(val.at(0) / mod.at(0)));
    if (!Numc::Valid(nrm)) nrm = Numc::ZERO<long double>;
    
    long double jacb = (Numc::EqualToZero(nrm) ? mod.at(1) : (val.at(1) / nrm));
    if (!Numc::Valid(jacb) || Numc::Compare(jacb) <= 0) jacb = Numc::ZERO<long double>;

    mini.at(0) = nrm;
    mini.at(1) = jacb;
    return mini;
}

std::array<long double, 2> LandauGaus::conv(long double norm) const { // (nrm, jacb)
    long double sum_prb = Numc::ZERO<long double>;
    long double sum_dev = Numc::ZERO<long double>;
    long double sum_hes = Numc::ZERO<long double>;
    long double sum_jb  = Numc::ZERO<long double>;
    for (int it = 0; it < CONV_N; ++it) {
        long double ix = norm - CONV_X[it] * nrmfluc_;
        std::array<long double, 5>&& val = eval(ix);
        long double prb = val.at(2) * CONV_P[it];
        sum_prb += prb;
        sum_dev += prb * val.at(3) * val.at(3);
        sum_hes += prb * val.at(4);
        sum_jb  += prb * val.at(3);
    }
    long double dev = sum_dev / sum_prb;
    long double hes = sum_hes / sum_prb;
    long double jb  = sum_jb  / sum_prb;
    
    if (!Numc::Valid(sum_prb) || !Numc::Valid(dev) || !Numc::Valid(hes) || !Numc::Valid(jb)) {
        sum_prb = Numc::ZERO<long double>;
        dev     = Numc::ZERO<long double>;
        hes     = Numc::ZERO<long double>;
        jb      = Numc::ZERO<long double>;
    }
    long double var_dev = (dev - jb * jb);
    if (Numc::Compare(var_dev) <= 0) var_dev = Numc::ZERO<long double>;
    if (Numc::Compare(hes)     <= 0) hes     = Numc::ZERO<long double>;

    long double jacb = (Numc::EqualToZero(norm - shft_) ? std::sqrt(hes + var_dev) : jb);
    if (!Numc::Valid(jacb)) jacb = Numc::ZERO<long double>;

    return std::array<long double, 2>({ sum_prb, jacb });
}


long double SftLandauGaus::Func(long double x, long double kpa, long double mpv, long double sgm, long double sft) {
    if (Numc::Compare(sgm) <= 0 || Numc::Compare(kpa) < 0 || Numc::Compare(kpa, Numc::ONE<long double>) > 0) return Numc::ZERO<long double>;
    bool is_sft = !Numc::EqualToZero(sft);
    
    long double norm = (x - mpv) / sgm;
    if (!Numc::Valid(norm)) return Numc::ZERO<long double>;

    long double gs    = -Numc::HALF * (norm - sft) * (norm - sft);
    long double ld    = -LandauNumc::NegLn(norm);
    long double ldgs  = kpa * gs + (Numc::ONE<long double> - kpa) * ld;
    long double value = std::exp(ldgs);

    long double mod = Numc::ZERO<long double>;
    long double scl = Numc::ONE<long double>;
    if (is_sft && !Numc::EqualToZero(kpa)) {
        short iter = 1;
        long double newx = kpa * sft;
        while (iter <= LMT_ITER) {
            long double dev = kpa * (newx - sft) + (1.0 - kpa) * LandauNumc::NegLnDev(newx);
            long double hes = kpa + (1.0 - kpa) * LandauNumc::NegLnHes(newx);
            long double dlt = -dev / hes;
            if (Numc::Compare(std::fabs(dlt), LMT_CONV) <= 0) break;
            newx += dlt;
            iter++;
        }
        if (iter > LMT_ITER) return Numc::ZERO<long double>;
        
        long double newgs   = -Numc::HALF * (newx - sft) * (newx - sft);
        long double newld   = -LandauNumc::NegLn(newx);
        long double newldgs = kpa * newgs + (Numc::ONE<long double> - kpa) * newld;
        long double newval  = std::exp(newldgs);
        if (Numc::EqualToZero(newval)) return Numc::ZERO<long double>;
        
        mod = newx;
        scl = Numc::ONE<long double> / newval;
        value *= scl;
    }
    
    if (!Numc::Valid(value)) value = Numc::ZERO<long double>;
    return value;
}


long double GammaErf::Func(long double x, long double alp, long double bta, long double scl, long double tune) {
    if (Numc::Compare(x) <= 0) return Numc::ZERO<long double>;
    if (Numc::Compare(alp, Numc::ONE<long double>) <= 0 || Numc::Compare(bta) <= 0) return Numc::ZERO<long double>;
    if (Numc::Compare(scl) <= 0 || Numc::Compare(tune) <= 0 || Numc::Compare(tune, Numc::ONE<long double>) >= 0) return Numc::ZERO<long double>;
    long double alpm1 = alp - Numc::ONE<long double>;
    long double mpv   = alpm1 / bta;
    long double fact  = std::pow(mpv, -alpm1) * std::exp(alpm1);
    if (!Numc::Valid(mpv) || !Numc::Valid(fact)) return Numc::ZERO<long double>;

    long double gmmv = fact * std::pow(x, alpm1) * std::exp(-bta * x);
    long double erfv = Numc::HALF * (Numc::ONE<long double> + std::erf(scl * (x - tune * mpv)));
    long double value = gmmv * erfv;

    if (!Numc::Valid(value)) value = Numc::ZERO<long double>;
    return value;
}


} // namesapce TrackSys


#endif // __TRACKLibs_Math_C__
