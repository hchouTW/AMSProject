#ifndef __TRACKLibs_IonEloss_C__
#define __TRACKLibs_IonEloss_C__


#include "Sys.h"
#include "Math.h"
#include "IonEloss.h"


namespace TrackSys {

long double IonEloss::FuncKpa(long double ibta, const std::array<long double, 5>& par) {
    if (Numc::Compare(ibta, Numc::ONE<long double>) < 0) return Numc::ZERO<long double>;
    long double loggb = std::log(ibta * ibta - Numc::ONE<long double>);
    
    long double kpa = Numc::HALF * (
        (Numc::ONE<long double> + std::erf(par[0] * loggb - par[1])) + 
        par[2] * std::erfc(par[3] * loggb + par[4]));
    if (!Numc::Valid(kpa)) kpa = Numc::ZERO<long double>;
    if (Numc::Compare(kpa, Numc::ZERO<long double>) <= 0) kpa = Numc::ZERO<long double>;
    if (Numc::Compare(kpa, Numc::ONE<long double> ) >= 0) kpa = Numc::ONE<long double>;
    return kpa;
}

long double IonEloss::FuncMpv(long double ibta, const std::array<long double, 4>& par) {
    if (Numc::Compare(ibta, Numc::ONE<long double>) < 0) return Numc::ZERO<long double>;
    long double ibsqr  = ibta * ibta;
    long double igbsqr = ibsqr - Numc::ONE<long double>;
 
    long double mpv = par[0] + par[1] * ibsqr - par[2] * std::log(par[3] + igbsqr);
    if (!Numc::Valid(mpv)) mpv = Numc::ZERO<long double>;
    return mpv;
}

std::array<long double, 3> IonEloss::minimizer(long double x, long double ibta, long double igb) const {
    if (Numc::Compare(ibta, Numc::ONE<long double>) < 0 || Numc::Compare(igb) < 0)
        return std::array<long double, 3>({ Numc::ZERO<long double>, Numc::ZERO<long double>, Numc::ZERO<long double> });
    long double ibsqr  = ibta * ibta;
    long double igbsqr = igb * igb;
    
    // PDF parameters
    long double kpa = get_kpa(igbsqr); 
    long double mpv = get_mpv(ibsqr, igbsqr); 
    long double sgm = get_sgm(ibsqr, igbsqr); 
    long double mod = (isfluc_ ? get_mod(ibsqr, igbsqr) : mpv);
    
    long double divIbta = Numc::NEG<long double> * (isfluc_ ? (get_divmod(ibta, igbsqr) / sgm) : (get_divmpv(ibta, igbsqr) / sgm));
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
        
long double IonEloss::get_kpa(long double igbsqr) const {
    long double loggb = std::log(igbsqr);
    long double kpa = Numc::HALF * ((Numc::ONE<long double> + std::erf(kpa_[0] * loggb - kpa_[1])) + kpa_[2] * std::erfc(kpa_[3] * loggb + kpa_[4]));
    if (!Numc::Valid(kpa)) kpa = Numc::ZERO<long double>;
    else {
        if (Numc::Compare(kpa, Numc::ZERO<long double>) <= 0) kpa = Numc::ZERO<long double>;
        if (Numc::Compare(kpa, Numc::ONE<long double>)  >= 0) kpa = Numc::ONE<long double>;
    }
    return kpa;
}

long double IonEloss::get_mpv(long double ibsqr, long double igbsqr) const {
    long double mpv = 
        mpv_[0] + 
        (Numc::EqualToZero(mpv_[1]) ? Numc::ZERO<long double> : mpv_[1] * ibsqr) +
        (Numc::EqualToZero(mpv_[2]) ? Numc::ZERO<long double> : Numc::NEG<long double> * mpv_[2] * std::log(mpv_[3] + igbsqr));
    if (!Numc::Valid(mpv)) mpv = Numc::ZERO<long double>;
    return mpv;
}

long double IonEloss::get_sgm(long double ibsqr, long double igbsqr) const {
    long double sgm = 
        sgm_[0] + 
        (Numc::EqualToZero(sgm_[1]) ? Numc::ZERO<long double> : sgm_[1] * ibsqr) +
        (Numc::EqualToZero(sgm_[2]) ? Numc::ZERO<long double> : Numc::NEG<long double> * sgm_[2] * std::log(sgm_[3] + igbsqr));
    if (!Numc::Valid(sgm)) sgm = Numc::ZERO<long double>;
    return sgm;
}

long double IonEloss::get_mod(long double ibsqr, long double igbsqr) const {
    long double mod = 
        mod_[0] + 
        (Numc::EqualToZero(mod_[1]) ? Numc::ZERO<long double> : mod_[1] * ibsqr) +
        (Numc::EqualToZero(mod_[2]) ? Numc::ZERO<long double> : Numc::NEG<long double> * mod_[2] * std::log(mod_[3] + igbsqr));
    if (!Numc::Valid(mod)) mod = Numc::ZERO<long double>;
    return mod;
}

long double IonEloss::get_divmpv(long double ibta, long double igbsqr) const {
    long double divbta = Numc::EqualToZero(mpv_[1]) ? Numc::ZERO<long double> : mpv_[1];
    long double divlog = Numc::EqualToZero(mpv_[2]) ? Numc::ZERO<long double> : Numc::NEG<long double> * mpv_[2] / (mpv_[3] + igbsqr);
    if (!Numc::Valid(divbta)) divbta = Numc::ZERO<long double>;
    if (!Numc::Valid(divlog)) divlog = Numc::ZERO<long double>;
    
    long double divmpv = (divbta + divlog) * (Numc::TWO<long double> * ibta);
    if (!Numc::Valid(divmpv)) divmpv = Numc::ZERO<long double>;
    return divmpv;
}

long double IonEloss::get_divmod(long double ibta, long double igbsqr) const {
    long double divbta = Numc::EqualToZero(mod_[1]) ? Numc::ZERO<long double> : mod_[1];
    long double divlog = Numc::EqualToZero(mod_[2]) ? Numc::ZERO<long double> : Numc::NEG<long double> * mod_[2] / (mod_[3] + igbsqr);
    if (!Numc::Valid(divbta)) divbta = Numc::ZERO<long double>;
    if (!Numc::Valid(divlog)) divlog = Numc::ZERO<long double>;
    
    long double divmod = (divbta + divlog) * (Numc::TWO<long double> * ibta);
    if (!Numc::Valid(divmod)) divmod = Numc::ZERO<long double>;
    return divmod;
}

} // namesapce TrackSys


#endif // __TRACKLibs_IonEloss_C__
