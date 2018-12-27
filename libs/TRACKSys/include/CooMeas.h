#ifndef __TRACKLibs_CooMeas_H__
#define __TRACKLibs_CooMeas_H__


namespace TrackSys {

class CooMeas {
    public :
        CooMeas(const MultiGaus& mgs) : mgs_(mgs), is_const_(true), isgm_(std::array<long double, 4>({0.0, 0.0, 0.0, 0.0})) {}
        CooMeas(const MultiGaus& mgs, const std::array<long double, 4>& isgm) : mgs_(mgs), is_const_(false), isgm_(isgm) {}
        ~CooMeas() {}
        
        std::array<long double, 3> minimizer(long double x, long double ibta = Numc::ONE<long double>, long double igb = Numc::ZERO<long double>) const;
        inline long double eftsgm(long double ibta = Numc::ONE<long double>, long double igb = Numc::ZERO<long double>) const { return (is_const_ ? mgs_.eftsgm() : (mgs_.eftsgm() / get_isgm(ibta, igb))); }

    protected :
        long double get_isgm(long double ibta, long double igb) const;

    private :
        MultiGaus mgs_;
        
        Bool_t is_const_;
        std::array<long double, 4> isgm_;
};

} // namesapce TrackSys


#endif // __TRACKLibs_CooMeas_H__
