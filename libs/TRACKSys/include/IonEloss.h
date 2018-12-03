#ifndef __TRACKLibs_IonEloss_H__
#define __TRACKLibs_IonEloss_H__


namespace TrackSys {
//x := igmbta
//y := ibta = 1 + x^2
//TF1* fkpa = new TF1("fkpa", "0.5 * (1.0 + TMath::Erf([0] * TMath::Log(1+[1]*(1+x*x)^[2]) - [3]))");
//TF1* fmpv = new TF1("fmpv", "[0] + [1] * (1+x*x)^[2]  - [3] * TMath::Log([4]+(x*x)^[5])");
//TF1* fkpa = new TF1("fkpa", "0.5 * (1.0 + TMath::Erf([0] * TMath::Log(1+[1]*(y*y)^[2]) - [3]))");
//TF1* fmpv = new TF1("fmpv", "[0] + [1] * (y*y)^[2]  - [3] * TMath::Log([4]+(y*y-1)^[5])");
class IonEloss {
    public :
        IonEloss(Robust robust, const std::array<long double, 4>& kpa, const std::array<long double, 6>& mpv, const std::array<long double, 6>& sgm, const std::array<long double, 6>& mod, long double fluc = Numc::ZERO<long double>) : robust_(robust), isfluc_(false), kpa_(kpa), mpv_(mpv), sgm_(sgm), mod_(mod), fluc_(fluc) { isfluc_ = (Numc::Compare(fluc_) > 0); if (!isfluc_) fluc_ = Numc::ZERO<long double>; }
        ~IonEloss() {}
        
        std::array<long double, 3> minimizer(long double x, long double ibta, long double igb) const;
   
    protected :
        long double get_kpa(long double ibsqr) const;
        long double get_mpv(long double ibsqr, long double igbsqr) const;
        long double get_sgm(long double ibsqr, long double igbsqr) const;
        long double get_mod(long double ibsqr, long double igbsqr) const;

        long double get_divmpv(long double ibta, long double ibsqr, long double igbsqr) const;
        long double get_divmod(long double ibta, long double ibsqr, long double igbsqr) const;

    private :
        Robust robust_;

        bool isfluc_;
        std::array<long double, 4> kpa_;
        std::array<long double, 6> mpv_;
        std::array<long double, 6> sgm_;
        std::array<long double, 6> mod_;
        long double                fluc_;
};

} // namesapce TrackSys


#endif // __TRACKLibs_IonEloss_H__
