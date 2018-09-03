#ifndef __TRACKLibs_IonEloss_C__
#define __TRACKLibs_IonEloss_C__


#include "Sys.h"
#include "Math.h"
#include "IonEloss.h"


namespace TrackSys {

std::array<long double, 3> IonEloss::minimizer(long double x, long double igmbta) const {
    if (Numc::Compare(x) <= 0 || Numc::Compare(igmbta) <= 0)
        return std::array<long double, 3>({ Numc::ZERO<long double>, Numc::ZERO<long double>, Numc::ONE<long double> });
    long double ibsqr = (Numc::ONE<long double> + igmbta * igmbta);
    
    // PDF parameters
    long double kpa  = get_kpa(igmbta, ibsqr); 
    long double mpv  = get_mpv(igmbta, ibsqr); 
    long double sgm  = get_sgm(igmbta, ibsqr); 
    long double mode = get_mode(igmbta, ibsqr); 
    
    long double divmpv = ((isfluc_) ? get_divmode(igmbta, ibsqr) : get_divmpv(igmbta, ibsqr));
    long double divsgm = get_divsgm(igmbta, ibsqr);

    // div r/z
    long double normx = ((x - (isfluc_ ? mode : mpv)) / sgm);
    //long double divnx = (Numc::NEG<long double> / sgm) * (divmpv);
    long double divnx = (Numc::NEG<long double> / sgm) * (divmpv + normx * divsgm); // testcode

    // approximate Landau-Gaussian
    LandauGaus ldgaus(robust_, kpa, mpv, sgm, mode, fluc_);
    std::array<long double, 3>&& lg_par = ldgaus.minimizer(x);
    
    long double chi = lg_par.at(0);         // res chiz  (z)
    long double nrm = lg_par.at(1);         // res normz (z)
    long double div = lg_par.at(2) * divnx; // div r/z * div z/igmbta
    
    if (!Numc::Valid(chi) || !Numc::Valid(nrm) || !Numc::Valid(div)) { 
        chi = Numc::ZERO<long double>;
        nrm = Numc::ZERO<long double>;
        div = Numc::ONE<long double>;
    }
    return std::array<long double, 3>({ chi, nrm, div });
}
        
long double IonEloss::get_kpa(long double igmbta, long double ibsqr) const {
    long double kpa = kpa_[0] + (Numc::ONE<long double> - kpa_[0]) * Numc::HALF * (Numc::ONE<long double> + std::erf(kpa_[1] * std::log1p(kpa_[2] * std::pow(ibsqr, kpa_[3])) - kpa_[4]));
    if (!Numc::Valid(kpa)) kpa = Numc::ZERO<long double>;
    else {
        if (Numc::Compare(kpa, Numc::ZERO<long double>) <= 0) kpa = Numc::ZERO<long double>;
        if (Numc::Compare(kpa, Numc::ONE<long double>)  >= 0) kpa = Numc::ONE<long double>;
    }
    return kpa;
}

long double IonEloss::get_mpv(long double igmbta, long double ibsqr) const {
    long double mpv = mpv_[0] * std::pow(ibsqr,  mpv_[3]) * 
        (mpv_[1] + 
         mpv_[2] * std::pow(ibsqr, -mpv_[3]) - 
         std::log1p(std::pow(igmbta * igmbta, mpv_[4]))
        );
    if (!Numc::Valid(mpv)) mpv = Numc::ZERO<long double>;
    return mpv;
}

long double IonEloss::get_sgm(long double igmbta, long double ibsqr) const {
    long double sgm = sgm_[0] * std::pow(ibsqr,  sgm_[3]) * 
        (sgm_[1] + 
         sgm_[2] * std::pow(ibsqr, -sgm_[3]) - 
         std::log1p(std::pow(igmbta * igmbta, sgm_[4]))
        );
    if (!Numc::Valid(sgm)) sgm = Numc::ZERO<long double>;
    return sgm;
}

long double IonEloss::get_mode(long double igmbta, long double ibsqr) const {
    long double mode = mode_[0] * std::pow(ibsqr,  mode_[3]) * 
        (mode_[1] + 
         mode_[2] * std::pow(ibsqr, -mode_[3]) - 
         std::log1p(std::pow(igmbta * igmbta, mode_[4]))
        );
    if (!Numc::Valid(mode)) mode = Numc::ZERO<long double>;
    return mode;
}

long double IonEloss::get_divmpv(long double igmbta, long double ibsqr) const {
    long double igmbtasqr = igmbta * igmbta;
    long double divbta = mpv_[3] * std::pow(ibsqr, mpv_[3] - Numc::ONE<long double>) * (Numc::TWO<long double> * igmbta);
    long double divlog = (Numc::TWO<long double> * mpv_[4]) * std::pow(igmbtasqr, mpv_[4] - Numc::HALF) / (Numc::ONE<long double> + std::pow(igmbtasqr, mpv_[4]));

    long double termA  = mpv_[1] * divbta;
    long double termB  = divbta * std::log1p(std::pow(igmbtasqr, mpv_[4]));
    long double termC  = std::pow(ibsqr, mpv_[3]) * divlog;
    long double divmpv = mpv_[0] * (termA - termB - termC);

    if (!Numc::Valid(divmpv)) divmpv = Numc::ZERO<long double>;
    return divmpv;
}

long double IonEloss::get_divsgm(long double igmbta, long double ibsqr) const {
    long double igmbtasqr = igmbta * igmbta;
    long double divbta = sgm_[3] * std::pow(ibsqr, sgm_[3] - Numc::ONE<long double>) * (Numc::TWO<long double> * igmbta);
    long double divlog = (Numc::TWO<long double> * sgm_[4]) * std::pow(igmbtasqr, sgm_[4] - Numc::HALF) / (Numc::ONE<long double> + std::pow(igmbtasqr, sgm_[4]));

    long double termA  = sgm_[1] * divbta;
    long double termB  = divbta * std::log1p(std::pow(igmbtasqr, sgm_[4]));
    long double termC  = std::pow(ibsqr, sgm_[3]) * divlog;
    long double divsgm = sgm_[0] * (termA - termB - termC);

    if (!Numc::Valid(divsgm)) divsgm = Numc::ZERO<long double>;
    return divsgm;
}

long double IonEloss::get_divmode(long double igmbta, long double ibsqr) const {
    long double igmbtasqr = igmbta * igmbta;
    long double divbta = mode_[3] * std::pow(ibsqr, mode_[3] - Numc::ONE<long double>) * (Numc::TWO<long double> * igmbta);
    long double divlog = (Numc::TWO<long double> * mode_[4]) * std::pow(igmbtasqr, mode_[4] - Numc::HALF) / (Numc::ONE<long double> + std::pow(igmbtasqr, mode_[4]));

    long double termA  = mode_[1] * divbta;
    long double termB  = divbta * std::log1p(std::pow(igmbtasqr, mode_[4]));
    long double termC  = std::pow(ibsqr, mode_[3]) * divlog;
    long double divmode = mode_[0] * (termA - termB - termC);

    if (!Numc::Valid(divmode)) divmode = Numc::ZERO<long double>;
    return divmode;
}

} // namesapce TrackSys


#endif // __TRACKLibs_IonEloss_C__
