#ifndef __TRACKLibs_TmeMeas_C__
#define __TRACKLibs_TmeMeas_C__


#include "Sys.h"
#include "Math.h"
#include "TmeMeas.h"


namespace TrackSys {


std::array<long double, 3> TmeMeas::minimizer(long double x, long double igmbta, bool is_single) const {
    if (Numc::Compare(igmbta) <= 0) return std::array<long double, 3>({ Numc::ZERO<long double>, Numc::ZERO<long double>, Numc::ONE<long double> });
    long double ibsqr = (Numc::ONE<long double> + igmbta * igmbta);
    
    // PDF parameters
    long double sgm = get_sgm(ibsqr);
    if (!is_single) sgm *= Numc::SQRT_TWO;

    // MultiGaus
    MultiGaus mgaus(robust_, sgm);
    return mgaus.minimizer(x);
}
        

long double TmeMeas::get_sgm(long double ibsqr) const {
    long double sgm = sgm_[0] + sgm_[1] * std::erfc(sgm_[2] * std::pow(ibsqr, sgm_[3]) - sgm_[4]);
    if (!Numc::Valid(sgm)) sgm = Numc::ZERO<long double>;
    return sgm;
}


} // namesapce TrackSys


#endif // __TRACKLibs_TmeMeas_C__
