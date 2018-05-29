#ifndef __TRACKLibs_GmIonEloss_H__
#define __TRACKLibs_GmIonEloss_H__


namespace TrackSys {
//x := igmbta
//TF1* fmpv = new TF1("fmpv", "[0] * (1+x*x)^[2] * ([1] - (1+x*x)^(-[2]) - TMath::Log([3] + x^[4])) + [5] * (TMath::Erfc([6] * TMath::Log(x) + [7]))");
//fmpv->SetParameters(8.77636e-02, 1.47297e+01, 1.20992e+00, 2.15340e-07, 3.07490e+00, 2.14623e+00, 8.23813e-01, 5.99981e+00);
class GmIonEloss {
    public :
        GmIonEloss(const std::array<long double, 8>& kpa, const std::array<long double, 8>& mpv, const std::array<long double, 8>& sgm, long double fluc = Numc::ZERO<long double>) : kpa_(kpa), mpv_(mpv), sgm_(sgm), fluc_(fluc) { if (Numc::Compare(fluc_) <= 0) fluc_ = Numc::ZERO<long double>; }
        ~GmIonEloss() {}
        
        inline SVecD<2> operator() (long double x, long double igmbta) const { std::array<long double, 2>&& gmion = eval(x, igmbta); return SVecD<2>(gmion.at(0), gmion.at(1)); }
   
    protected :
        std::array<long double, 2> eval(long double x, long double igmbta) const;

        inline long double eval_kpa(long double igmbta, long double ibsqr) const;
        inline long double eval_mpv(long double igmbta, long double ibsqr) const;
        inline long double eval_sgm(long double igmbta, long double ibsqr) const;
        
        inline long double eval_divmpv(long double igmbta, long double ibsqr) const;

    private :
        std::array<long double, 8> kpa_;
        std::array<long double, 8> mpv_;
        std::array<long double, 8> sgm_;
        long double                fluc_;
};

} // namesapce TrackSys


#endif // __TRACKLibs_GmIonEloss_H__
