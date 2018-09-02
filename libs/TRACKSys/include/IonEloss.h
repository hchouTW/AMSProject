#ifndef __TRACKLibs_IonEloss_H__
#define __TRACKLibs_IonEloss_H__


namespace TrackSys {
//x := igmbta
//TF1* fkpa = new TF1("fkpa", "0.5 * (1.0 + TMath::Erf([0] * TMath::Log(1+[1]*(x*x)) + [2]))");
//fkpa->SetParameters(1.5, 1.6, -2.0);
//TF1* fmpv = new TF1("fmpv", "[0] * (1+x*x)^[3] * ([1] - [2]*(1+x*x)^(-[3]) - TMath::Log([4]+x^[5]))");
//fmpv->SetParameters(10, 6.5, 1.0, 1.0, 10.0, 1.0);
class IonEloss {
    public :
        IonEloss(Robust robust, const std::array<long double, 4>& kpa, const std::array<long double, 6>& mpv, const std::array<long double, 6>& sgm, const std::array<long double, 6>& mode, long double fluc = Numc::ZERO<long double>) : robust_(robust), isfluc_(false), kpa_(kpa), mpv_(mpv), sgm_(sgm), mode_(mode), fluc_(fluc) { isfluc_ = (Numc::Compare(fluc_) > 0); if (!isfluc_) fluc_ = Numc::ZERO<long double>; }
        ~IonEloss() {}
        
        std::array<long double, 3> minimizer(long double x, long double igmbta) const;
   
    protected :
        long double get_kpa(long double igmbta, long double ibsqr) const;
        long double get_mpv(long double igmbta, long double ibsqr) const;
        long double get_sgm(long double igmbta, long double ibsqr) const;
        long double get_mode(long double igmbta, long double ibsqr) const;

        long double get_divmpv(long double igmbta, long double ibsqr) const;
        long double get_divmode(long double igmbta, long double ibsqr) const;

    private :
        Robust robust_;

        bool isfluc_;
        std::array<long double, 4> kpa_;
        std::array<long double, 6> mpv_;
        std::array<long double, 6> sgm_;
        std::array<long double, 6> mode_;
        long double                fluc_;
};

} // namesapce TrackSys


#endif // __TRACKLibs_IonEloss_H__
