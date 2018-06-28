#ifndef __TRACKLibs_IonEloss_H__
#define __TRACKLibs_IonEloss_H__


namespace TrackSys {
//x := igmbta
//TF1* fkpa = new TF1("fkpa", "(1 - [0]) + [0] * 0.5 * (1.0 + TMath::Erf([1] * TMath::Log(1+[3]*(x*x)) + [2]))");
//fkpa->SetParameters(1.0, 1.0, 0.3, 1.0);
//TF1* fmpv = new TF1("fmpv", "[0] * (1+x*x)^[3] * ([1] - [2]*(1+x*x)^(-[3]) - TMath::Log([4]+x^[5]))");
//fmpv->SetParameters(10, 6.5, 1.0, 1.0, 10.0, 1.0);
    
class IonEloss {
    public :
        IonEloss(const std::array<long double, 4>& kpa, const std::array<long double, 6>& mpv, const std::array<long double, 6>& sgm, long double fluc = Numc::ZERO<long double>) : kpa_(kpa), mpv_(mpv), sgm_(sgm), fluc_(fluc) { if (Numc::Compare(fluc_) <= 0) fluc_ = Numc::ZERO<long double>; }
        ~IonEloss() {}
        
        inline SVecD<2> operator()(long double x, long double igmbta) const { std::array<long double, 2>&& ion = eval(x, igmbta); return SVecD<2>(ion.at(0), ion.at(1)); }
   
    protected :
        std::array<long double, 2> eval(long double x, long double igmbta) const;

        long double get_kpa(long double igmbta, long double ibsqr) const;
        long double get_mpv(long double igmbta, long double ibsqr) const;
        long double get_sgm(long double igmbta, long double ibsqr) const;
        
        long double get_divmpv(long double igmbta, long double ibsqr) const;

    private :
        std::array<long double, 4> kpa_;
        std::array<long double, 6> mpv_;
        std::array<long double, 6> sgm_;
        long double                fluc_;
};

} // namesapce TrackSys


#endif // __TRACKLibs_IonEloss_H__
