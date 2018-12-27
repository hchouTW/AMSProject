#ifndef __TRACKLibs_CooMeas_C__
#define __TRACKLibs_CooMeas_C__


#include "Sys.h"
#include "Math.h"
#include "CooMeas.h"


namespace TrackSys {


std::array<long double, 3> CooMeas::minimizer(long double x, long double ibta, long double igb) const {
    if (Numc::Compare(ibta, Numc::ONE<long double>) < 0 || Numc::Compare(igb) < 0)
        return std::array<long double, 3>({ Numc::ZERO<long double>, Numc::ZERO<long double>, Numc::ZERO<long double> });

    long double sclx = (is_const_ ? x : (x * get_isgm(ibta, igb)));
    return mgs_.minimizer(sclx);
}
        

long double CooMeas::get_isgm(long double ibta, long double igb) const {
    long double isgm = isgm_[0] + isgm_[1] * ibta - isgm_[2] * std::log(isgm_[3] + igb);
    if (!Numc::Valid(isgm)) isgm = Numc::ONE<long double>;
    return isgm;
}


} // namesapce TrackSys


#endif // __TRACKLibs_CooMeas_C__
