#ifndef __TRACKLibs_IonTrEloss_C__
#define __TRACKLibs_IonTrEloss_C__


#include "Sys.h"
#include "Math.h"
#include "IonTrEloss.h"


namespace TrackSys {

long double IonTrEloss::FuncKpa(long double ibta, const std::array<long double, 7>& par) {
    if (Numc::Compare(ibta, Numc::ONE<long double>) < 0) return Numc::ZERO<long double>;
    long double ibsqr  = ibta * ibta;
    long double igbsqr = ibsqr - Numc::ONE<long double>;
    long double loggb = std::log(ibta * ibta - Numc::ONE<long double>);

    long double kpa = Numc::HALF * (
                      (Numc::ONE<long double> + std::erf(par[0] * loggb - par[1])) +
                      (par[2] - par[3] * std::log(LMT_IGBSQR + igbsqr)) +
                      par[4] * std::erfc(par[5] * loggb + par[6]));
    if (!Numc::Valid(kpa)) kpa = Numc::ZERO<long double>;
    else {
        if (Numc::Compare(kpa, Numc::ZERO<long double>) <= 0) kpa = Numc::ZERO<long double>;
        if (Numc::Compare(kpa, Numc::ONE<long double> ) >= 0) kpa = Numc::ONE<long double>;
    }
    return kpa;
}

long double IonTrEloss::FuncMpv(long double ibta, const std::array<long double, 10>& par) {
    if (Numc::Compare(ibta, Numc::ONE<long double>) < 0) return Numc::ZERO<long double>;
    long double ibsqr  = ibta * ibta;
    long double igbsqr = ibsqr - Numc::ONE<long double>;
    long double loggb  = std::log(igbsqr);
    
    long double mpv = par[0] + 
                      par[1] * std::pow(ibsqr, par[2]) - 
                      par[3] * std::log(LMT_IGBSQR + igbsqr) + 
                      par[4] * Numc::HALF * std::erfc(par[5] * loggb + par[6]) +
                      par[7] * Numc::HALF * std::erfc(par[8] * loggb + par[9]);
    if (!Numc::Valid(mpv)) mpv = Numc::ZERO<long double>;
    return mpv;
}

std::array<long double, 3> IonTrEloss::minimizer(long double x, long double ibta, long double igb) const {
    if (Numc::Compare(ibta, Numc::ONE<long double>) < 0 || Numc::Compare(igb) < 0)
        return std::array<long double, 3>({ Numc::ZERO<long double>, Numc::ZERO<long double>, Numc::ZERO<long double> });
    long double ibsqr  = ibta * ibta;
    long double igbsqr = igb * igb;
    long double loggb  = std::log(igbsqr);
    
    // PDF parameters
    long double kpa = get_kpa(igbsqr, loggb); 
    long double mpv = get_mpv(ibsqr, igbsqr, loggb); 
    long double sgm = get_sgm(ibsqr, igbsqr, loggb); 
    long double mod = (isfluc_ ? get_mod(ibsqr, igbsqr, loggb) : mpv);
    
    long double divIbta = Numc::NEG<long double> * (isfluc_ ? (get_divmod(ibta, igbsqr, loggb) / sgm) : (get_divmpv(ibta, igbsqr, loggb) / sgm));
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
        
long double IonTrEloss::get_kpa(long double igbsqr, long double loggb) const {
    long double kpa = Numc::HALF * (
                      (Numc::ONE<long double> + std::erf(kpa_[0] * loggb - kpa_[1])) + 
                      (kpa_[2] - kpa_[3] * std::log(LMT_IGBSQR + igbsqr)) +
                      kpa_[4] * std::erfc(kpa_[5] * loggb + kpa_[6]));
    if (!Numc::Valid(kpa)) kpa = Numc::ZERO<long double>;
    else {
        if (Numc::Compare(kpa, Numc::ZERO<long double>) <= 0) kpa = Numc::ZERO<long double>;
        if (Numc::Compare(kpa, Numc::ONE<long double>)  >= 0) kpa = Numc::ONE<long double>;
    }
    return kpa;
}

long double IonTrEloss::get_mpv(long double ibsqr, long double igbsqr, long double loggb) const {
    long double mpv = 
        mpv_[0] + 
        (Numc::EqualToZero(mpv_[1]) ? Numc::ZERO<long double> : mpv_[1] * std::pow(ibsqr, mpv_[2])) +
        (Numc::EqualToZero(mpv_[3]) ? Numc::ZERO<long double> : Numc::NEG<long double> * mpv_[3] * std::log(LMT_IGBSQR + igbsqr)) +
        (Numc::EqualToZero(mpv_[4]) ? Numc::ZERO<long double> : mpv_[4] * Numc::HALF * std::erfc(mpv_[5] * loggb + mpv_[6])) +
        (Numc::EqualToZero(mpv_[7]) ? Numc::ZERO<long double> : mpv_[7] * Numc::HALF * std::erfc(mpv_[8] * loggb + mpv_[9]));
    if (!Numc::Valid(mpv)) mpv = Numc::ZERO<long double>;
    return mpv;
}

long double IonTrEloss::get_sgm(long double ibsqr, long double igbsqr, long double loggb) const {
    long double sgm = 
        sgm_[0] + 
        (Numc::EqualToZero(sgm_[1]) ? Numc::ZERO<long double> : sgm_[1] * std::pow(ibsqr, sgm_[2])) +
        (Numc::EqualToZero(sgm_[3]) ? Numc::ZERO<long double> : Numc::NEG<long double> * sgm_[3] * std::log(LMT_IGBSQR + igbsqr)) +
        (Numc::EqualToZero(sgm_[4]) ? Numc::ZERO<long double> : sgm_[4] * Numc::HALF * std::erfc(sgm_[5] * loggb + sgm_[6])) +
        (Numc::EqualToZero(sgm_[7]) ? Numc::ZERO<long double> : sgm_[7] * Numc::HALF * std::erfc(sgm_[8] * loggb + sgm_[9]));
    if (!Numc::Valid(sgm)) sgm = Numc::ZERO<long double>;
    return sgm;
}

long double IonTrEloss::get_mod(long double ibsqr, long double igbsqr, long double loggb) const {
    long double mod = 
        mod_[0] + 
        (Numc::EqualToZero(mod_[1]) ? Numc::ZERO<long double> : mod_[1] * std::pow(ibsqr, mod_[2])) +
        (Numc::EqualToZero(mod_[3]) ? Numc::ZERO<long double> : Numc::NEG<long double> * mod_[3] * std::log(LMT_IGBSQR + igbsqr)) +
        (Numc::EqualToZero(mod_[4]) ? Numc::ZERO<long double> : mod_[4] * Numc::HALF * std::erfc(mod_[5] * loggb + mod_[6])) +
        (Numc::EqualToZero(mod_[7]) ? Numc::ZERO<long double> : mod_[7] * Numc::HALF * std::erfc(mod_[8] * loggb + mod_[9]));
    if (!Numc::Valid(mod)) mod = Numc::ZERO<long double>;
    return mod;
}

long double IonTrEloss::get_divmpv(long double ibta, long double igbsqr, long double loggb) const {
    long double divbta = Numc::EqualToZero(mpv_[1]) ? Numc::ZERO<long double> : mpv_[1] * mpv_[2] * std::pow(ibta * ibta, mpv_[2] - Numc::ONE<long double>);
    long double divlog = Numc::EqualToZero(mpv_[3]) ? Numc::ZERO<long double> : Numc::NEG<long double> * mpv_[3] / (LMT_IGBSQR + igbsqr);
    long double diverf1 = Numc::EqualToZero(mpv_[4]) ? Numc::ZERO<long double> : Numc::NEG<long double> * mpv_[4] * mpv_[5] * Numc::INV_SQRT_PI * std::exp(-(mpv_[5] * loggb + mpv_[6]) * (mpv_[5] * loggb + mpv_[6])) / igbsqr;
    long double diverf2 = Numc::EqualToZero(mpv_[7]) ? Numc::ZERO<long double> : Numc::NEG<long double> * mpv_[7] * mpv_[8] * Numc::INV_SQRT_PI * std::exp(-(mpv_[8] * loggb + mpv_[9]) * (mpv_[8] * loggb + mpv_[9])) / igbsqr;
    if (!Numc::Valid(divbta)) divbta = Numc::ZERO<long double>;
    if (!Numc::Valid(divlog)) divlog = Numc::ZERO<long double>;
    if (!Numc::Valid(diverf1)) diverf1 = Numc::ZERO<long double>;
    if (!Numc::Valid(diverf2)) diverf2 = Numc::ZERO<long double>;
    
    long double divmpv = (divbta + divlog + diverf1 + diverf2) * (Numc::TWO<long double> * ibta);
    if (!Numc::Valid(divmpv)) divmpv = Numc::ZERO<long double>;
    return divmpv;
}

long double IonTrEloss::get_divmod(long double ibta, long double igbsqr, long double loggb) const {
    long double divbta = Numc::EqualToZero(mod_[1]) ? Numc::ZERO<long double> : mod_[1] * mod_[2] * std::pow(ibta * ibta, mod_[2] - Numc::ONE<long double>);
    long double divlog = Numc::EqualToZero(mod_[3]) ? Numc::ZERO<long double> : Numc::NEG<long double> * mod_[3] / (LMT_IGBSQR + igbsqr);
    long double diverf1 = Numc::EqualToZero(mod_[4]) ? Numc::ZERO<long double> : Numc::NEG<long double> * mod_[4] * mod_[5] * Numc::INV_SQRT_PI * std::exp(-(mod_[5] * loggb + mod_[6]) * (mod_[5] * loggb + mod_[6])) / igbsqr;
    long double diverf2 = Numc::EqualToZero(mod_[7]) ? Numc::ZERO<long double> : Numc::NEG<long double> * mod_[7] * mod_[8] * Numc::INV_SQRT_PI * std::exp(-(mod_[8] * loggb + mod_[9]) * (mod_[8] * loggb + mod_[9])) / igbsqr;
    if (!Numc::Valid(divbta)) divbta = Numc::ZERO<long double>;
    if (!Numc::Valid(divlog)) divlog = Numc::ZERO<long double>;
    if (!Numc::Valid(diverf1)) diverf1 = Numc::ZERO<long double>;
    if (!Numc::Valid(diverf2)) diverf2 = Numc::ZERO<long double>;
    
    long double divmod = (divbta + divlog + diverf1 + diverf2) * (Numc::TWO<long double> * ibta);
    if (!Numc::Valid(divmod)) divmod = Numc::ZERO<long double>;
    return divmod;
}

} // namesapce TrackSys


#endif // __TRACKLibs_IonTrEloss_C__
