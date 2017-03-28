#ifndef __MGNumc_H__
#define __MGNumc_H__

#include <numeric>


namespace MGNumc {
	constexpr double PI    = 3.14159265358979323846;
	constexpr double EXP   = 2.71828182845904523536;
	constexpr double SQRT2 = 1.41421356237309504880;
	constexpr double SQRT3 = 1.73205080756887719318; 
	constexpr double LN2   = 0.69314718055994528623;
	constexpr double LN10  = 2.30258509299404590109;

	// TODO
	template <class IntType = long, typename std::enable_if<std::is_integral<IntType>::value, int>::type = 0>
	inline bool Valid(long long int var) { return (var == IntType(var)); }

	template <class RealType = double, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type = 0>
	inline bool Valid(RealType var) { int clf = std::fpclassify(var); return (clf == FP_NORMAL || clf == FP_ZERO); }
	
	template <class IntType = long, typename std::enable_if<std::is_integral<IntType>::value, int>::type = 0>
	inline short Compare(IntType a, IntType b = IntType(0));
	
	template <class RealType = double, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type = 0>
	inline short Compare(RealType a, RealType b = RealType(0.0));
	
	template <class IntType = long, typename std::enable_if<std::is_integral<IntType>::value, int>::type = 0>
	inline bool Equal(IntType a, IntType b) { return (Compare(a, b) == 0); }
	
	template <class RealType = double, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type = 0>
	inline bool Equal(RealType a, RealType b) { return (Compare(a, b) == 0); }
	
	template <class IntType = long, typename std::enable_if<std::is_integral<IntType>::value, int>::type = 0>
	inline bool EqualToZero(IntType a) { return (Compare(a) == 0); }
	
	template <class RealType = double, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type = 0>
	inline bool EqualToZero(RealType a) { return (Compare(a) == 0); }
	
	template <class IntType = long, typename std::enable_if<std::is_integral<IntType>::value, int>::type = 0>
	inline IntType Square(IntType val) { return (val * val); }
	
	template <class RealType = double, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type = 0>
	inline RealType Square(RealType val) { return (val * val); }
}


#endif // __MGNumc_H__
