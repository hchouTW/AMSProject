#ifndef __TRACKLibs_TmeMeas_H__
#define __TRACKLibs_TmeMeas_H__


namespace TrackSys {

class TmeMeas {
    public :
        TmeMeas(const MultiGaus& mgs) : mgs_(mgs), is_const_(true), isgm_(std::array<long double, 2>({0.0, 0.0})) {}
        TmeMeas(const MultiGaus& mgs, const std::array<long double, 2>& isgm) : mgs_(mgs), is_const_(false), isgm_(isgm) {}
        ~TmeMeas() {}
        
        std::array<long double, 3> minimizer(long double x, long double ibta = Numc::ONE<long double>, bool is_single = true) const;
        inline long double eftsgm(long double ibta = Numc::ONE<long double>) const { return (is_const_ ? mgs_.eftsgm() : (mgs_.eftsgm() / get_isgm(ibta))); }

    protected :
        long double get_isgm(long double ibta) const;

    private :
        MultiGaus mgs_;
        
        Bool_t is_const_;
        std::array<long double, 2> isgm_;
};

} // namesapce TrackSys


#endif // __TRACKLibs_TmeMeas_H__
