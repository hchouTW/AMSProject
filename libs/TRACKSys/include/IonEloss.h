#ifndef __TRACKLibs_IonEloss_H__
#define __TRACKLibs_IonEloss_H__


namespace TrackSys {
//x := igmbta
//TF1* fkpa = new TF1("fkpa", "1.0 - 0.5*TMath::Erfc([0]*TMath::Log(1.0+[1]*x)+[2])");
//fkpa->SetParameters(1.0, 0.3, 0.0);
//TF1* fmpv = new TF1("fmpv", "[0] * (1+x*x)^[2] * ([1] - (1+x*x)^(-[2]) - TMath::Log([3] + x^[4]))");
//fmpv->SetParameters(10, 6.5, 1.0, 10.0, 1.0);
class IonEloss {
    public :
        IonEloss(const std::array<long double, 3>& kpa, const std::array<long double, 5>& mpv, const std::array<long double, 5>& sgm, long double fluc = Numc::ZERO<long double>) : kpa_(kpa), mpv_(mpv), sgm_(sgm), fluc_(fluc) { if (Numc::Compare(fluc_) <= 0) fluc_ = Numc::ZERO<long double>; }
        ~IonEloss() {}
        
        inline SVecD<2> operator()(long double x, long double igmbta) const { std::array<long double, 2>&& ion = eval(x, igmbta); return SVecD<2>(ion.at(0), ion.at(1)); }
   
    protected :
        std::array<long double, 2> eval(long double x, long double igmbta) const;

        inline long double eval_kpa(long double igmbta, long double ibsqr) const;
        inline long double eval_mpv(long double igmbta, long double ibsqr) const;
        inline long double eval_sgm(long double igmbta, long double ibsqr) const;
        
        inline long double eval_divmpv(long double igmbta, long double ibsqr) const;

    private :
        std::array<long double, 3> kpa_;
        std::array<long double, 5> mpv_;
        std::array<long double, 5> sgm_;
        long double                fluc_;
};

} // namesapce TrackSys


#endif // __TRACKLibs_IonEloss_H__
