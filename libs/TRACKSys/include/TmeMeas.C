#ifndef __TRACKLibs_TmeMeas_C__
#define __TRACKLibs_TmeMeas_C__


#include "Sys.h"
#include "Math.h"
#include "TmeMeas.h"


namespace TrackSys {


std::array<long double, 3> TmeMeas::minimizer(long double x, long double ibta, bool is_single) const {
    if (Numc::Compare(ibta, Numc::ONE<long double>) < 0)
        return std::array<long double, 3>({ Numc::ZERO<long double>, Numc::ZERO<long double>, Numc::ZERO<long double> });

    long double sclx = (is_const_ ? x : (x * get_isgm(ibta)));
    if (!is_single) sclx *= Numc::INV_SQRT_TWO;
    return mgs_.minimizer(sclx);
}
        

long double TmeMeas::get_isgm(long double ibta) const {
    long double isgm = isgm_[0] + isgm_[1] * ibta;
    if (!Numc::Valid(isgm)) isgm = Numc::ONE<long double>;
    return isgm;
}


} // namesapce TrackSys


#endif // __TRACKLibs_TmeMeas_C__
