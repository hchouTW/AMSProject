#ifndef __TRACKLibs_GmIonEloss_H__
#define __TRACKLibs_GmIonEloss_H__


namespace TrackSys {

class GmIonEloss {
    public :
        GmIonEloss(
            const std::array<long double, 3>& ratio,
            const std::array<long double, 2>& lg_k, const std::array<long double, 6>& lg_m, const std::array<long double, 6>& lg_s,
            const std::array<long double, 1>& ge_a, const std::array<long double, 3>& ge_b, const std::array<long double, 1>& ge_m, const std::array<long double, 1>& ge_s
        ) : ratio_(ratio), lg_k_(lg_k), lg_m_(lg_m), lg_s_(lg_s), ge_a_(ge_a), ge_b_(ge_b), ge_m_(ge_m), ge_s_(ge_s) {}
        ~GmIonEloss() {}

        inline SVecD<3> operator()(long double x, long double igmbta) const { std::array<long double, 3>&& lgge = eval(x, igmbta); return SVecD<3>(lgge.at(0), lgge.at(1), lgge.at(2)); }
        
    protected :
        std::array<long double, 3> eval(long double x, long double igmbta) const;
        
        long double get_ratio(long double igmbta) const;
        
        long double get_lg_k(long double igmbta) const;
        long double get_lg_m(long double igmbta, long double ibsqr) const;
        long double get_lg_s(long double igmbta, long double ibsqr) const;

        long double get_ge_a() const { return ge_a_.at(0); }
        long double get_ge_b(long double igmbta) const; 
        long double get_ge_m() const { return ge_m_.at(0); }
        long double get_ge_s() const { return ge_s_.at(0); }
        
        long double get_lg_divm(long double igmbta, long double ibsqr) const;

    private :
        std::array<long double, 3> ratio_; // ratio := TR / (ION + TR)
                                                                      
        std::array<long double, 2> lg_k_;  // kpa
        std::array<long double, 6> lg_m_;  // mpv
        std::array<long double, 6> lg_s_;  // sgm
                                                                      
        std::array<long double, 1> ge_a_;  // alpha
        std::array<long double, 3> ge_b_;  // beta
        std::array<long double, 1> ge_m_;  // men
        std::array<long double, 1> ge_s_;  // sgm
};


} // namesapce TrackSys


#endif // __TRACKLibs_GmIonEloss_H__
