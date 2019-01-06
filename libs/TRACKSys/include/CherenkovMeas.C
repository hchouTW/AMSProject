#ifndef __TRACKLibs_CherenkovMeas_C__
#define __TRACKLibs_CherenkovMeas_C__


#include "Sys.h"
#include "Math.h"
#include "CherenkovMeas.h"


namespace TrackSys {


std::array<long double, 3> CherenkovMeas::minimizer(long double x, long double ibta) const {
    if (Numc::Compare(ibta, Numc::ONE<long double>) < 0)
        return std::array<long double, 3>({ Numc::ZERO<long double>, Numc::ZERO<long double>, Numc::ZERO<long double> });

    bool rescl = false;
    long double eftsgm = Numc::ONE<long double>;
    if (Numc::Compare(rfr_, Numc::ONE<long double>) > 0) {
        long double thres = rfr_ - Numc::TWO<long double> * sgm_;
        if (Numc::Compare(ibta, thres) > 0) {
            long double extsgm = (ibta - thres) / sgm_;
            eftsgm = Numc::ONE<long double> + extsgm * extsgm;
            rescl = true;
        }
    }

    long double sclx = (rescl ? (x / eftsgm) : x);
    return mgs_.minimizer(sclx);
}


} // namesapce TrackSys


#endif // __TRACKLibs_CherenkovMeas_C__
