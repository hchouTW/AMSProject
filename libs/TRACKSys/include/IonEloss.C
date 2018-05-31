#ifndef __TRACKLibs_IonEloss_C__
#define __TRACKLibs_IonEloss_C__


namespace TrackSys {

std::array<long double, 2> IonEloss::eval(long double x, long double igmbta) const {
    if (Numc::Compare(x) <= 0 || Numc::EqualToZero(igmbta))
        return std::array<long double, 2>({Numc::ZERO<long double>, Numc::ZERO<long double>});
    long double ibsqr  = (Numc::ONE<long double> + igmbta * igmbta);
    
    // PDF parameters
    long double kpa    = eval_kpa(igmbta, ibsqr); 
    long double mpv    = eval_mpv(igmbta, ibsqr); 
    long double sgm    = eval_sgm(igmbta, ibsqr); 
    long double divmpv = eval_divmpv(igmbta, ibsqr); 
 
    // Landau-Gaus with noise fluctuation 
    LandauGaus ldgaus(LandauGaus::Opt::ROBUST, kpa, mpv, sgm, fluc_);
    std::array<long double, 2>&& lg_par = ldgaus(x);
    
    long double res = lg_par.at(0);          // res normx
    long double div = lg_par.at(1) * divmpv; // div r/x * div x/igmbta
    if (!Numc::Valid(res) || !Numc::Valid(div)) { 
        res = Numc::ZERO<long double>;
        div = Numc::ZERO<long double>;
    }
    return std::array<long double, 2>({res, div});
}
        
long double IonEloss::eval_kpa(long double igmbta, long double ibsqr) const {
    long double kpa = Numc::ONE<long double> - Numc::HALF *
        std::erfc(kpa_.at(0) * std::log(Numc::ONE<long double>+kpa_.at(1)*igmbta) + kpa_.at(2));
    if (!Numc::Valid(kpa)) kpa = Numc::ZERO<long double>;
    else {
        if (Numc::Compare(kpa, Numc::ZERO<long double>) <= 0) kpa = Numc::ZERO<long double>;
        if (Numc::Compare(kpa, Numc::ONE<long double>)  >= 0) kpa = Numc::ONE<long double>;
    }
    return kpa;
}

long double IonEloss::eval_mpv(long double igmbta, long double ibsqr) const {
    long double mpv = mpv_.at(0) * std::pow(ibsqr,  mpv_.at(2)) * 
        (mpv_.at(1) - 
         std::pow(ibsqr, -mpv_.at(2)) - 
         std::log(mpv_.at(3) + std::pow(igmbta, mpv_.at(4)))
        );
    if (!Numc::Valid(mpv)) mpv = Numc::ZERO<long double>;
    return mpv;
}

long double IonEloss::eval_sgm(long double igmbta, long double ibsqr) const {
    long double sgm = sgm_.at(0) * std::pow(ibsqr,  sgm_.at(2)) * 
        (sgm_.at(1) - 
         std::pow(ibsqr, -sgm_.at(2)) - 
         std::log(sgm_.at(3) + std::pow(igmbta, sgm_.at(4)))
        );
    if (!Numc::Valid(sgm)) sgm = Numc::ZERO<long double>;
    return sgm;
}

long double IonEloss::eval_divmpv(long double igmbta, long double ibsqr) const {
    long double divbta = mpv_.at(2) * std::pow(ibsqr, mpv_.at(2)-Numc::ONE<long double>) * (Numc::TWO<long double> * igmbta);
    long double divlog = mpv_.at(4) * std::pow(igmbta, mpv_.at(4)-Numc::ONE<long double>) / (mpv_.at(3) + std::pow(igmbta, mpv_.at(4)));

    long double termA  = mpv_.at(1) * divbta;
    long double termB  = divbta * std::log(mpv_.at(3) + std::pow(igmbta, mpv_.at(4)));
    long double termC  = std::pow(ibsqr, mpv_.at(2)) * divlog;
    long double divmpv = mpv_.at(0) * (termA - termB - termC);

    if (!Numc::Valid(divmpv)) divmpv = Numc::ZERO<long double>;
    
    //CERR("IGB %14.8f DIV %14.8f\n", static_cast<Double_t>(igmbta), static_cast<Double_t>(divmpv));
    return divmpv;
}

} // namesapce TrackSys


#endif // __TRACKLibs_IonEloss_C__
