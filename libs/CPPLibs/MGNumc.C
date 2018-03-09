#ifndef __CPPLibs_MGNumc_C__
#define __CPPLibs_MGNumc_C__

#include "MGNumc.h"

namespace MGNumc {

template <class IntType, typename std::enable_if<std::is_integral<IntType>::value, int>::type>
inline short Compare(IntType a, IntType b) {
	if (a == b)     return  0;
	else if (a > b) return  1;
	else            return -1;
}


template <class RealType, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type>
inline short Compare(RealType a, RealType b) {
	RealType diff = std::fabs(a - b);
    if (!std::isfinite(diff)) return 0;
	if (diff < std::numeric_limits<RealType>::epsilon() * 5.0e3) return  0;
	else if (a > b)                                              return  1;
	else                                                         return -1;
}

} // namespace MGNumc

#endif // __CPPLibs_MGNumc_C__
