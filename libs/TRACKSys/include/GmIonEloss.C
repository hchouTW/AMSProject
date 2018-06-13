#ifndef __TRACKLibs_GmIonEloss_C__
#define __TRACKLibs_GmIonEloss_C__


namespace TrackSys {


std::array<long double, 3> GmIonEloss::eval(long double x, long double igmbta) const {
    if (Numc::Compare(x) <= 0 || Numc::EqualToZero(igmbta))
        return std::array<long double, 3>({ Numc::ZERO<long double>, Numc::ZERO<long double>, Numc::ZERO<long double> });
    long double ibsqr  = (Numc::ONE<long double> + igmbta * igmbta);
    
    // PDF parameters
    long double ratio     = get_ratio(igmbta);
    long double lg_kpa    = get_lg_k(igmbta); 
    long double lg_mpv    = get_lg_m(igmbta, ibsqr); 
    long double lg_sgm    = get_lg_s(igmbta, ibsqr); 
    long double ge_alpha  = get_ge_a();
    long double ge_beta   = get_ge_b(igmbta);
    long double ge_men    = get_ge_m();
    long double ge_sgm    = get_ge_s();
    long double lg_divmpv = get_lg_divm(igmbta, ibsqr); 
   
    // testcode
    //CERR("IGMBTA %14.8f DEDX %14.8f RAT %14.8f LG %14.8f %14.8f %14.8f GE %14.8f %14.8f %14.8f %14.8f\n", 
    //        static_cast<Double_t>(igmbta), 
    //        static_cast<Double_t>(x), 
    //        static_cast<Double_t>(ratio), 
    //        static_cast<Double_t>(lg_kpa), static_cast<Double_t>(lg_mpv), static_cast<Double_t>(lg_sgm), 
    //        static_cast<Double_t>(ge_alpha), static_cast<Double_t>(ge_beta), static_cast<Double_t>(ge_men), static_cast<Double_t>(ge_sgm));
 
    // Landau-Gaus with noise fluctuation 
    LgGeFunc lgge(LgGeFunc::Opt::ROBUST, ratio, lg_kpa, lg_mpv, lg_sgm, ge_alpha, ge_beta, ge_men, ge_sgm);
    std::array<long double, 3>&& lgge_par = lgge.minimizer(x);

    long double lg_res = lgge_par.at(0);             // res normx
    long double lg_div = lgge_par.at(1) * lg_divmpv; // div r/x * div x/igmbta
    if (!Numc::Valid(lg_res) || !Numc::Valid(lg_div)) { 
        lg_res = Numc::ZERO<long double>;
        lg_div = Numc::ZERO<long double>;
    }
    
    return std::array<long double, 3>({ lg_res, lg_div, ratio });
}

long double GmIonEloss::get_ratio(long double igmbta) const {
    long double rat = Numc::HALF * ratio_.at(0) * std::erfc(ratio_.at(1) * (std::log(igmbta) - ratio_.at(2)));
    if (!Numc::Valid(rat)) rat = Numc::ZERO<long double>;
    else {
        if (Numc::Compare(rat, Numc::ZERO<long double>) <= 0) rat = Numc::ZERO<long double>;
        if (Numc::Compare(rat, Numc::ONE<long double>)  >= 0) rat = Numc::ONE<long double>;
    }
    return rat;
}
        
long double GmIonEloss::get_lg_k(long double igmbta) const {
    long double kpa = lg_k_.at(0) + Numc::HALF * (Numc::ONE<long double> - lg_k_.at(0)) * (Numc::ONE<long double> + std::erf(igmbta - lg_k_.at(1)));
    if (!Numc::Valid(kpa)) kpa = Numc::ZERO<long double>;
    else {
        if (Numc::Compare(kpa, Numc::ZERO<long double>) <= 0) kpa = Numc::ZERO<long double>;
        if (Numc::Compare(kpa, Numc::ONE<long double>)  >= 0) kpa = Numc::ONE<long double>;
    }
    return kpa;
}
        
long double GmIonEloss::get_lg_m(long double igmbta, long double ibsqr) const {
    long double mpv = lg_m_.at(0) * std::pow(ibsqr, lg_m_.at(3)) * (
            lg_m_.at(1) - 
            lg_m_.at(2) * std::pow(ibsqr, -lg_m_.at(3)) - 
            std::log(lg_m_.at(4) + std::pow(igmbta, lg_m_.at(5))));
    if (!Numc::Valid(mpv) || Numc::Compare(mpv) <= 0) mpv = Numc::ZERO<long double>;
    return mpv;
}

long double GmIonEloss::get_lg_s(long double igmbta, long double ibsqr) const {
    long double sgm = lg_s_.at(0) * std::pow(ibsqr, lg_s_.at(3)) * (
            lg_s_.at(1) - 
            lg_s_.at(2) * std::pow(ibsqr, -lg_s_.at(3)) - 
            std::log(lg_s_.at(4) + std::pow(igmbta, lg_s_.at(5))));
    if (!Numc::Valid(sgm) || Numc::Compare(sgm) <= 0) sgm = Numc::ZERO<long double>;
    return sgm;
}

long double GmIonEloss::get_ge_b(long double igmbta) const { 
    long double beta = ge_b_.at(0) + Numc::HALF * (Numc::ONE<long double> - ge_b_.at(0)) * (Numc::ONE<long double> + std::erf(ge_b_.at(1) * (std::log(igmbta) - ge_b_.at(2))));
    if (!Numc::Valid(beta) || Numc::Compare(beta) <= 0) beta = Numc::ZERO<long double>;
    return beta; 
}

long double GmIonEloss::get_lg_divm(long double igmbta, long double ibsqr) const {
    long double divbta = lg_m_.at(3) * std::pow(ibsqr, lg_m_.at(3) - Numc::ONE<long double>) * (Numc::TWO<long double> * igmbta);
    long double divlog = lg_m_.at(5) * std::pow(igmbta, lg_m_.at(5) - Numc::ONE<long double>) / (lg_m_.at(4) + std::pow(igmbta, lg_m_.at(5)));

    long double termA  = lg_m_.at(1) * divbta;
    long double termB  = divbta * std::log(lg_m_.at(4) + std::pow(igmbta, lg_m_.at(5)));
    long double termC  = std::pow(ibsqr, lg_m_.at(3)) * divlog;
    long double divmpv = lg_m_.at(0) * (termA - termB - termC);

    if (!Numc::Valid(divmpv)) divmpv = Numc::ZERO<long double>;
    return divmpv;
}


} // namesapce TrackSys


#endif // __TRACKLibs_GmIonEloss_C__
