#ifndef __TRACKLibs_IonTrEloss_H__
#define __TRACKLibs_IonTrEloss_H__


namespace TrackSys {

class IonTrEloss {
    public :
        static long double FuncKpa(long double ibta, const std::array<long double, 8>& par);
        static long double FuncMpv(long double ibta, const std::array<long double, 8>& par);
        static long double FuncSgm(long double ibta, const std::array<long double, 8>& par);

    public :
        IonTrEloss(Robust robust, const std::array<long double, 8>& kpa, const std::array<long double, 8>& mpv, const std::array<long double, 8>& sgm, const std::array<long double, 8>& mod, long double fluc = Numc::ZERO<long double>) : robust_(robust), isfluc_(false), kpa_(kpa), mpv_(mpv), sgm_(sgm), mod_(mod), fluc_(fluc) { isfluc_ = (Numc::Compare(fluc_) > 0); if (!isfluc_) fluc_ = Numc::ZERO<long double>; }
        ~IonTrEloss() {}
        
        std::array<long double, 3> minimizer(long double x, long double ibta, long double igb) const;
   
    protected :
        long double get_kpa(long double igbsqr, long double logigb) const;
        long double get_mpv(long double ibsqr, long double igbsqr, long double logigb) const;
        long double get_sgm(long double ibsqr, long double igbsqr, long double logigb) const;
        long double get_mod(long double ibsqr, long double igbsqr, long double logigb) const;

        // derivative ibta
        long double get_divmpv(long double ibta, long double ibsqr, long double igbsqr, long double logigb) const;
        long double get_divmod(long double ibta, long double ibsqr, long double igbsqr, long double logigb) const;

    private :
        Robust robust_;

        bool isfluc_;
        std::array<long double, 8> kpa_;
        std::array<long double, 8> mpv_;
        std::array<long double, 8> sgm_;
        std::array<long double, 8> mod_;
        long double                fluc_;
};

} // namesapce TrackSys


#endif // __TRACKLibs_IonTrEloss_H__
