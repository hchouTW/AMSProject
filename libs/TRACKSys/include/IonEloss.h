#ifndef __TRACKLibs_IonEloss_H__
#define __TRACKLibs_IonEloss_H__


namespace TrackSys {
//x := igmbta
//TF1* fkpa = new TF1("fkpa", "0.5*(1+TMath::Erf([0]*log(x*x)-[1])) + [2]*0.5*TMath::Erfc([3]*log(x*x)+[4])");
//TF1* fmpv = new TF1("fmpv", "[0] + [1] * (1+x*x) - [2] * TMath::Log([3]+(x*x))");
class IonEloss {
    public :
        static long double FuncKpa(long double ibta, const std::array<long double, 5>& par);
        static long double FuncMpv(long double ibta, const std::array<long double, 4>& par);

    public :
        IonEloss(Robust robust, const std::array<long double, 5>& kpa, const std::array<long double, 4>& mpv, const std::array<long double, 4>& sgm, const std::array<long double, 4>& mod, long double fluc = Numc::ZERO<long double>) : robust_(robust), isfluc_(false), kpa_(kpa), mpv_(mpv), sgm_(sgm), mod_(mod), fluc_(fluc) { isfluc_ = (Numc::Compare(fluc_) > 0); if (!isfluc_) fluc_ = Numc::ZERO<long double>; }
        ~IonEloss() {}
        
        std::array<long double, 3> minimizer(long double x, long double ibta, long double igb) const;
   
    protected :
        long double get_kpa(long double igbsqr) const;
        long double get_mpv(long double ibsqr, long double igbsqr) const;
        long double get_sgm(long double ibsqr, long double igbsqr) const;
        long double get_mod(long double ibsqr, long double igbsqr) const;

        // derivative ibta
        long double get_divmpv(long double ibta, long double igbsqr) const;
        long double get_divmod(long double ibta, long double igbsqr) const;

    private :
        Robust robust_;

        bool isfluc_;
        std::array<long double, 5> kpa_;
        std::array<long double, 4> mpv_;
        std::array<long double, 4> sgm_;
        std::array<long double, 4> mod_;
        long double                fluc_;
};

} // namesapce TrackSys


#endif // __TRACKLibs_IonEloss_H__
