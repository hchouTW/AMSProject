#ifndef __TRACKLibs_IonTrEloss_H__
#define __TRACKLibs_IonTrEloss_H__


namespace TrackSys {
//x := igmbta
//TF1* fkpa = new TF1("fkpa", "0.5*(1+TMath::Erf([0]*log(x*x)-[1])) + ([2]-[3]*TMath::Log(1.0e-07+(x*x))) + [5]*0.5*TMath::Erfc([6]*log(x*x)+[7])");
//TF1* fmpv = new TF1("fmpv", "[0] + [1] * TMath::Power((1+x*x), [2]) - [3] * TMath::Log((x*x)) + [4]*0.5*TMath::Erfc([5] * log(x*x) + [6]) + [7]*0.5*TMath::Erfc([8] * log(x*x) + [9])");
class IonTrEloss {
    public :
        static long double FuncKpa(long double ibta, const std::array<long double, 7>& par);
        static long double FuncMpv(long double ibta, const std::array<long double, 10>& par);

    public :
        IonTrEloss(Robust robust, const std::array<long double, 7>& kpa, const std::array<long double, 10>& mpv, const std::array<long double, 10>& sgm, const std::array<long double, 10>& mod, long double fluc = Numc::ZERO<long double>) : robust_(robust), isfluc_(false), kpa_(kpa), mpv_(mpv), sgm_(sgm), mod_(mod), fluc_(fluc) { isfluc_ = (Numc::Compare(fluc_) > 0); if (!isfluc_) fluc_ = Numc::ZERO<long double>; }
        ~IonTrEloss() {}
        
        std::array<long double, 3> minimizer(long double x, long double ibta, long double igb) const;
   
    protected :
        long double get_kpa(long double igbsqr, long double loggb) const;
        long double get_mpv(long double ibsqr, long double igbsqr, long double loggb) const;
        long double get_sgm(long double ibsqr, long double igbsqr, long double loggb) const;
        long double get_mod(long double ibsqr, long double igbsqr, long double loggb) const;

        // derivative ibta
        long double get_divmpv(long double ibta, long double igbsqr, long double loggb) const;
        long double get_divmod(long double ibta, long double igbsqr, long double loggb) const;

    private :
        Robust robust_;

        bool isfluc_;
        std::array<long double, 7>  kpa_;
        std::array<long double, 10> mpv_;
        std::array<long double, 10> sgm_;
        std::array<long double, 10> mod_;
        long double                 fluc_;

    private :
        static constexpr long double LMT_IGBSQR = 1.0e-10;
};

} // namesapce TrackSys


#endif // __TRACKLibs_IonTrEloss_H__
