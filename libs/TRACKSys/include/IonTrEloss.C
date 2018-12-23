#ifndef __TRACKLibs_IonTrEloss_C__
#define __TRACKLibs_IonTrEloss_C__


#include "Sys.h"
#include "Math.h"
#include "IonTrEloss.h"


namespace TrackSys {

long double IonTrEloss::FuncKpa(long double ibta, const std::array<long double, 8>& par) {
    if (Numc::Compare(ibta, Numc::ONE<long double>) <= 0) return Numc::ZERO<long double>;
    long double igbsqr = (ibta - Numc::ONE<long double>) * (ibta + Numc::ONE<long double>);
    long double logigb = std::log(igbsqr);

    long double kpa = Numc::ONE_TO_TWO * (
                      (Numc::ONE<long double> + std::erf(par[0] * logigb - par[1])) +
                      par[2] * std::erfc(par[3] * logigb + par[4]) +
                      par[5] * std::erfc(par[6] * logigb + par[7]));
    if (!Numc::Valid(kpa)) kpa = Numc::ZERO<long double>;
    else {
        if (Numc::Compare(kpa, Numc::ZERO<long double>) <= 0) kpa = Numc::ZERO<long double>;
        if (Numc::Compare(kpa, Numc::ONE<long double> ) >= 0) kpa = Numc::ONE<long double>;
    }
    return kpa;
}

long double IonTrEloss::FuncMpv(long double ibta, const std::array<long double, 8>& par) {
    if (Numc::Compare(ibta, Numc::ONE<long double>) <= 0) return Numc::ZERO<long double>;
    long double ibsqr  = ibta * ibta;
    long double igbsqr = (ibta - Numc::ONE<long double>) * (ibta + Numc::ONE<long double>);
    long double logigb = std::log(igbsqr);
   
    long double mpv = par[0] * std::pow(ibsqr, par[1]) + Numc::ONE_TO_TWO * (
                      par[2] * std::erfc(par[3] * logigb + par[4]) +
                      par[5] * std::erfc(par[6] * logigb + par[7]));
    if (!Numc::Valid(mpv) || Numc::Compare(mpv) <= 0) mpv = Numc::ZERO<long double>;
    return mpv;
}

long double IonTrEloss::FuncSgm(long double ibta, const std::array<long double, 8>& par) {
    if (Numc::Compare(ibta, Numc::ONE<long double>) <= 0) return Numc::ZERO<long double>;
    long double ibsqr  = ibta * ibta;
    long double igbsqr = (ibta - Numc::ONE<long double>) * (ibta + Numc::ONE<long double>);
    long double logigb = std::log(igbsqr);
 
    long double sgm = par[0] * std::pow(ibsqr, par[1]) + Numc::ONE_TO_TWO * (
                      par[2] * std::erfc(par[3] * logigb + par[4]) +
                      par[5] * std::erfc(par[6] * logigb + par[7]));
    if (!Numc::Valid(sgm) || Numc::Compare(sgm) <= 0) sgm = Numc::ZERO<long double>;
    return sgm;
}

std::array<long double, 3> IonTrEloss::minimizer(long double x, long double ibta, long double igb) const {
    if (Numc::Compare(x) <= 0 || Numc::Compare(ibta) <= 0 || Numc::Compare(igb) <= 0)
        return std::array<long double, 3>({ Numc::ZERO<long double>, Numc::ZERO<long double>, Numc::ZERO<long double> });
    long double ibsqr  = ibta * ibta;
    long double igbsqr = igb * igb;
    long double logigb = std::log(igbsqr);
    
    // PDF parameters
    long double kpa = get_kpa(igbsqr, logigb); 
    long double mpv = get_mpv(ibsqr, igbsqr, logigb); 
    long double sgm = get_sgm(ibsqr, igbsqr, logigb); 
    long double mod = (isfluc_ ? get_mod(ibsqr, igbsqr, logigb) : mpv);
    
    long double divIbta = Numc::NEG<long double> * (isfluc_ ? (get_divmod(ibta, ibsqr, igbsqr, logigb) / sgm) : (get_divmpv(ibta, ibsqr, igbsqr, logigb) / sgm));
    if (!Numc::Valid(divIbta)) divIbta = Numc::ZERO<long double>;

    // approximate Landau-Gaussian
    LandauGaus ldgaus(robust_, kpa, mpv, sgm, mod, fluc_);
    std::array<long double, 3>&& lg_par = ldgaus.minimizer(x);
    
    long double chi = lg_par.at(0);           // res chiz  (z)
    long double nrm = lg_par.at(1);           // res normz (z)
    long double div = lg_par.at(2) * divIbta; // div r/z * div z/ibta

    if (!Numc::Valid(chi) || !Numc::Valid(nrm) || !Numc::Valid(div)) { 
        chi = Numc::ZERO<long double>;
        nrm = Numc::ZERO<long double>;
        div = Numc::ZERO<long double>;
    }

    return std::array<long double, 3>({ chi, nrm, div });
}
        
long double IonTrEloss::get_kpa(long double igbsqr, long double logigb) const {
    long double kpa  = Numc::ONE_TO_TWO * (
                       (Numc::ONE<long double> + std::erf(kpa_[0] * logigb - kpa_[1])) +
                       kpa_[2] * std::erfc(kpa_[3] * logigb + kpa_[4]) +
                       kpa_[5] * std::erfc(kpa_[6] * logigb + kpa_[7]));
    if (!Numc::Valid(kpa)) kpa = Numc::ZERO<long double>;
    else {
        if (Numc::Compare(kpa, Numc::ZERO<long double>) <= 0) kpa = Numc::ZERO<long double>;
        if (Numc::Compare(kpa, Numc::ONE<long double>)  >= 0) kpa = Numc::ONE<long double>;
    }
    return kpa;
}

long double IonTrEloss::get_mpv(long double ibsqr, long double igbsqr, long double logigb) const {
    long double mpv = 
        (Numc::EqualToZero(mpv_[0]) ? Numc::ZERO<long double> : mpv_[0] * std::pow(ibsqr, mpv_[1])) + Numc::ONE_TO_TWO * (
        (Numc::EqualToZero(mpv_[2]) ? Numc::ZERO<long double> : mpv_[2] * std::erfc(mpv_[3] * logigb + mpv_[4])) +
        (Numc::EqualToZero(mpv_[5]) ? Numc::ZERO<long double> : mpv_[5] * std::erfc(mpv_[6] * logigb + mpv_[7])));
    if (!Numc::Valid(mpv)) mpv = Numc::ZERO<long double>;
    return mpv;
}

long double IonTrEloss::get_sgm(long double ibsqr, long double igbsqr, long double logigb) const {
    long double sgm = 
        sgm_[0] + 
        (Numc::EqualToZero(sgm_[0]) ? Numc::ZERO<long double> : sgm_[0] * std::pow(ibsqr, sgm_[1])) + Numc::ONE_TO_TWO * (
        (Numc::EqualToZero(sgm_[2]) ? Numc::ZERO<long double> : sgm_[2] * std::erfc(sgm_[3] * logigb + sgm_[4])) +
        (Numc::EqualToZero(sgm_[5]) ? Numc::ZERO<long double> : sgm_[5] * std::erfc(sgm_[6] * logigb + sgm_[7])));
    if (!Numc::Valid(sgm)) sgm = Numc::ZERO<long double>;
    return sgm;
}

long double IonTrEloss::get_mod(long double ibsqr, long double igbsqr, long double logigb) const {
    long double mod = 
        mod_[0] + 
        (Numc::EqualToZero(mod_[0]) ? Numc::ZERO<long double> : mod_[0] * std::pow(ibsqr, mod_[1])) + Numc::ONE_TO_TWO * (
        (Numc::EqualToZero(mod_[2]) ? Numc::ZERO<long double> : mod_[2] * std::erfc(mod_[3] * logigb + mod_[4])) +
        (Numc::EqualToZero(mod_[5]) ? Numc::ZERO<long double> : mod_[5] * std::erfc(mod_[6] * logigb + mod_[7])));
    if (!Numc::Valid(mod)) mod = Numc::ZERO<long double>;
    return mod;
}

long double IonTrEloss::get_divmpv(long double ibta, long double ibsqr, long double igbsqr, long double logigb) const {
    long double divbta = Numc::EqualToZero(mpv_[0]) ? Numc::ZERO<long double> : (mpv_[0] * mpv_[1] * std::pow(ibsqr, mpv_[1] - Numc::ONE<long double>));
    long double erfcv1 = Numc::EqualToZero(mpv_[2]) ? Numc::ZERO<long double> : (mpv_[3] * logigb + mpv_[4]);
    long double divtr1 = Numc::EqualToZero(mpv_[2]) ? Numc::ZERO<long double> : (-Numc::INV_SQRT_PI * mpv_[2] * mpv_[3] * (std::exp(-erfcv1 * erfcv1) / igbsqr));
    long double erfcv2 = Numc::EqualToZero(mpv_[5]) ? Numc::ZERO<long double> : (mpv_[6] * logigb + mpv_[7]);
    long double divtr2 = Numc::EqualToZero(mpv_[5]) ? Numc::ZERO<long double> : (-Numc::INV_SQRT_PI * mpv_[5] * mpv_[6] * (std::exp(-erfcv2 * erfcv2) / igbsqr));
    long double divmpv = (divbta + divtr1 + divtr2) * (Numc::TWO<long double> * ibta);
    if (!Numc::Valid(divmpv)) divmpv = Numc::ZERO<long double>;
    return divmpv;
}

long double IonTrEloss::get_divmod(long double ibta, long double ibsqr, long double igbsqr, long double logigb) const {
    long double divbta = Numc::EqualToZero(mod_[0]) ? Numc::ZERO<long double> : (mod_[0] * mod_[1] * std::pow(ibsqr, mod_[1] - Numc::ONE<long double>));
    long double erfcv1 = Numc::EqualToZero(mod_[2]) ? Numc::ZERO<long double> : (mod_[3] * logigb + mod_[4]);
    long double divtr1 = Numc::EqualToZero(mod_[2]) ? Numc::ZERO<long double> : (-Numc::INV_SQRT_PI * mod_[2] * mod_[3] * (std::exp(-erfcv1 * erfcv1) / igbsqr));
    long double erfcv2 = Numc::EqualToZero(mod_[5]) ? Numc::ZERO<long double> : (mod_[6] * logigb + mod_[7]);
    long double divtr2 = Numc::EqualToZero(mod_[5]) ? Numc::ZERO<long double> : (-Numc::INV_SQRT_PI * mod_[5] * mod_[6] * (std::exp(-erfcv2 * erfcv2) / igbsqr));
    long double divmod = (divbta + divtr1 + divtr2) * (Numc::TWO<long double> * ibta);
    if (!Numc::Valid(divmod)) divmod = Numc::ZERO<long double>;
    return divmod;
}

} // namesapce TrackSys


#endif // __TRACKLibs_IonTrEloss_C__
