#ifndef __TRACKLibs_TmeMeas_H__
#define __TRACKLibs_TmeMeas_H__


namespace TrackSys {
//x := igmbta
//TF1* ftme = new TF1("ftme", "[0] + [1] * TMath::Erfc([2] * (1+x*x)^[3] - [4])");
class TmeMeas {
    public :
        TmeMeas(Robust robust, const std::array<long double, 5>& sgm) : robust_(robust), sgm_(sgm) {}
        ~TmeMeas() {}
        
        std::array<long double, 4> minimizer(long double x, long double igmbta, bool is_single = true) const;
   
    protected :
        long double get_sgm(long double ibsqr) const;

    private :
        Robust robust_;
        std::array<long double, 5> sgm_;
};

} // namesapce TrackSys


#endif // __TRACKLibs_TmeMeas_H__
