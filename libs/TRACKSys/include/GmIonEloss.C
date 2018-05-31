#ifndef __TRACKLibs_GmIonEloss_C__
#define __TRACKLibs_GmIonEloss_C__


namespace TrackSys {

std::array<long double, 2> GmIonEloss::eval(long double x, long double igmbta) const {
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
        
long double GmIonEloss::eval_kpa(long double igmbta, long double ibsqr) const {
    long double kpa = 
        kpa_.at(0) * std::pow(ibsqr,  kpa_.at(2)) * 
        (kpa_.at(1) - 
         std::pow(ibsqr, -kpa_.at(2)) - 
         std::log(kpa_.at(3) + std::pow(igmbta, kpa_.at(4)))
        ) +
        kpa_.at(5) * std::erfc(kpa_.at(6) * std::log(igmbta) + kpa_.at(7));
    if (!Numc::Valid(kpa)) kpa = Numc::ZERO<long double>;
    return kpa;
}

long double GmIonEloss::eval_mpv(long double igmbta, long double ibsqr) const {
    long double mpv = 
        mpv_.at(0) * std::pow(ibsqr,  mpv_.at(2)) * 
        (mpv_.at(1) - 
         std::pow(ibsqr, -mpv_.at(2)) - 
         std::log(mpv_.at(3) + std::pow(igmbta, mpv_.at(4)))
        ) +
        mpv_.at(5) * std::erfc(mpv_.at(6) * std::log(igmbta) + mpv_.at(7));
    if (!Numc::Valid(mpv)) mpv = Numc::ZERO<long double>;
    return mpv;
}

long double GmIonEloss::eval_sgm(long double igmbta, long double ibsqr) const {
    long double sgm = 
        sgm_.at(0) * std::pow(ibsqr,  sgm_.at(2)) * 
        (sgm_.at(1) - 
         std::pow(ibsqr, -sgm_.at(2)) - 
         std::log(sgm_.at(3) + std::pow(igmbta, sgm_.at(4)))
        ) +
        sgm_.at(5) * std::erfc(sgm_.at(6) * std::log(igmbta) + sgm_.at(7));
    if (!Numc::Valid(sgm)) sgm = Numc::ZERO<long double>;
    return sgm;
}

long double GmIonEloss::eval_divmpv(long double igmbta, long double ibsqr) const {
    long double divbta = mpv_.at(2) * std::pow(ibsqr, mpv_.at(2)-Numc::ONE<long double>) * (Numc::TWO<long double> * igmbta);
    long double divlog = mpv_.at(4) * std::pow(igmbta, mpv_.at(4)-Numc::ONE<long double>) / (mpv_.at(3) + std::pow(igmbta, mpv_.at(4)));

    long double termA  = mpv_.at(1) * divbta;
    long double termB  = divbta * std::log(mpv_.at(3) + std::pow(igmbta, mpv_.at(4)));
    long double termC  = std::pow(ibsqr, mpv_.at(2)) * divlog;
    long double divion = mpv_.at(0) * (termA - termB - termC);
    if (!Numc::Valid(divion)) divion = Numc::ZERO<long double>;

    long double logg    = (mpv_.at(6) * std::log(igmbta) + mpv_.at(7));
    long double divlogg = (mpv_.at(6) / igmbta);
    long double diverfc = Numc::NEG<long double> * (Numc::TWO<long double> * Numc::INV_SQRT_PI) * std::exp(Numc::NEG<long double> * logg * logg);
    long double divgm   = mpv_.at(5) * diverfc * divlogg;
    if (!Numc::Valid(divgm)) divgm = Numc::ZERO<long double>;

    long double divmpv = (divion + divgm);
    if (!Numc::Valid(divmpv)) divmpv = Numc::ZERO<long double>;

    //CERR("IGMB %14.8f DIV %14.8f ION %14.8f GM %14.8f\n", static_cast<Double_t>(igmbta), static_cast<Double_t>(divmpv), static_cast<Double_t>(divion), static_cast<Double_t>(divgm));
    return divmpv;
}

} // namesapce TrackSys


#endif // __TRACKLibs_GmIonEloss_C__
