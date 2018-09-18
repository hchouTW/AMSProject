#ifndef __TRACKLibs_TmeMeas_H__
#define __TRACKLibs_TmeMeas_H__


namespace TrackSys {
//x := igmbta
class TmeMeas {
    public :
        TmeMeas(Robust robust, const std::array<long double, 6>& sgm) : robust_(robust), sgm_(sgm) {}
        ~TmeMeas() {}
        
        std::array<long double, 3> minimizer(long double x, long double igmbta) const;
   
    protected :
        long double get_sgm(long double igmbta, long double ibsqr) const;

    private :
        Robust robust_;
        std::array<long double, 6> sgm_;
};

} // namesapce TrackSys


#endif // __TRACKLibs_TmeMeas_H__
