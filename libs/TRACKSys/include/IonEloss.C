#ifndef __TRACKLibs_IonEloss_C__
#define __TRACKLibs_IonEloss_C__


#include "Sys.h"
#include "Math.h"
#include "IonEloss.h"


namespace TrackSys {

std::array<long double, 2> IonEloss::minimizer(long double x, long double igmbta) const {
    if (Numc::Compare(x) <= 0 || Numc::EqualToZero(igmbta))
        return std::array<long double, 2>({ Numc::ZERO<long double>, Numc::ZERO<long double> });
    long double ibsqr = (Numc::ONE<long double> + igmbta * igmbta);
    
    // PDF parameters
    long double kpa  = get_kpa(igmbta, ibsqr); 
    long double mpv  = get_mpv(igmbta, ibsqr); 
    long double sgm  = get_sgm(igmbta, ibsqr); 
    long double mode = get_mode(igmbta, ibsqr); 
    
    long double divmpv = ((isfluc_) ? get_mode(igmbta, ibsqr) : get_divmpv(igmbta, ibsqr));
 
    // approximate Landau-Gaussian
    //LandauGaus ldgaus(Robust::Opt::ON, kpa, mpv, sgm, mode, fluc_);
    LandauGaus ldgaus(Robust::Opt::OFF, kpa, mpv, sgm, mode, fluc_); // testcode
    std::array<long double, 2>&& lg_par = ldgaus.minimizer(x);
    
    long double res = lg_par.at(0);                                   // res normx
    long double div = Numc::NEG<long double> * lg_par.at(1) * divmpv; // div r/x * div x/igmbta
    
    if (!Numc::Valid(res) || !Numc::Valid(div)) { 
        res = Numc::ZERO<long double>;
        div = Numc::ZERO<long double>;
    }
    return std::array<long double, 2>({ res, div });
}
        
long double IonEloss::get_kpa(long double igmbta, long double ibsqr) const {
    long double kpa = Numc::HALF * (Numc::ONE<long double> + std::erf(kpa_[0] * std::log1p(kpa_[1] * igmbta * igmbta) + kpa_[2]));
    if (!Numc::Valid(kpa)) kpa = Numc::ZERO<long double>;
    else {
        if (Numc::Compare(kpa, Numc::ZERO<long double>) <= 0) kpa = Numc::ZERO<long double>;
        if (Numc::Compare(kpa, Numc::ONE<long double>)  >= 0) kpa = Numc::ONE<long double>;
    }
    return kpa;
}

long double IonEloss::get_mpv(long double igmbta, long double ibsqr) const {
    long double mpv = mpv_[0] * std::pow(ibsqr,  mpv_[3]) * 
        (mpv_[1] - 
         mpv_[2] * std::pow(ibsqr, -mpv_[3]) - 
         std::log(mpv_[4] + std::pow(igmbta, mpv_[5]))
        );
    if (!Numc::Valid(mpv)) mpv = Numc::ZERO<long double>;
    return mpv;
}

long double IonEloss::get_sgm(long double igmbta, long double ibsqr) const {
    long double sgm = sgm_[0] * std::pow(ibsqr,  sgm_[3]) * 
        (sgm_[1] - 
         sgm_[2] * std::pow(ibsqr, -sgm_[3]) - 
         std::log(sgm_[4] + std::pow(igmbta, sgm_[5]))
        );
    if (!Numc::Valid(sgm)) sgm = Numc::ZERO<long double>;
    return sgm;
}

long double IonEloss::get_mode(long double igmbta, long double ibsqr) const {
    long double mode = mode_[0] * std::pow(ibsqr,  mode_[3]) * 
        (mode_[1] - 
         mode_[2] * std::pow(ibsqr, -mode_[3]) - 
         std::log(mode_[4] + std::pow(igmbta, mode_[5]))
        );
    if (!Numc::Valid(mode)) mode = Numc::ZERO<long double>;
    return mode;
}

long double IonEloss::get_divmpv(long double igmbta, long double ibsqr) const {
    long double divbta = mpv_[3] * std::pow(ibsqr, mpv_[3] - Numc::ONE<long double>) * (Numc::TWO<long double> * igmbta);
    long double divlog = mpv_[5] * std::pow(igmbta, mpv_[5] - Numc::ONE<long double>) / (mpv_[4] + std::pow(igmbta, mpv_[5]));

    long double termA  = mpv_[1] * divbta;
    long double termB  = divbta * std::log(mpv_[4] + std::pow(igmbta, mpv_[5]));
    long double termC  = std::pow(ibsqr, mpv_[3]) * divlog;
    long double divmpv = mpv_[0] * (termA - termB - termC);

    if (!Numc::Valid(divmpv)) divmpv = Numc::ZERO<long double>;
    return divmpv;
}

long double IonEloss::get_divmode(long double igmbta, long double ibsqr) const {
    long double divbta = mode_[3] * std::pow(ibsqr, mode_[3] - Numc::ONE<long double>) * (Numc::TWO<long double> * igmbta);
    long double divlog = mode_[5] * std::pow(igmbta, mode_[5] - Numc::ONE<long double>) / (mode_[4] + std::pow(igmbta, mode_[5]));

    long double termA  = mode_[1] * divbta;
    long double termB  = divbta * std::log(mode_[4] + std::pow(igmbta, mode_[5]));
    long double termC  = std::pow(ibsqr, mode_[3]) * divlog;
    long double divmode = mode_[0] * (termA - termB - termC);

    if (!Numc::Valid(divmode)) divmode = Numc::ZERO<long double>;
    return divmode;
}

} // namesapce TrackSys


#endif // __TRACKLibs_IonEloss_C__
