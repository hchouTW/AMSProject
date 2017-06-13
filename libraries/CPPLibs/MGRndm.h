#ifndef __CPPLibs_MGRndm_H__
#define __CPPLibs_MGRndm_H__

#include <random>
#include <type_traits>
#include <functional>


namespace MGRndm {

// Random Engines
static long rndmSeed = std::chrono::system_clock::now().time_since_epoch().count();
static std::mt19937_64 rndmEngMT64(rndmSeed);

// Uniform Distributions
template <class IntType = long, typename std::enable_if<std::is_integral<IntType>::value, int>::type = 0>
inline std::function<IntType()> Uniform(IntType a = 0, IntType b = std::numeric_limits<IntType>::max());

template <class RealType = double, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type = 0>
inline std::function<RealType()> Uniform(RealType a = 0.0, RealType b = std::numeric_limits<RealType>::max());

// Bernoulli Distributions
inline std::function<bool()> Bernoulli(double p = 0.5);

// Bnomial Distributions
template <class IntType = long, typename std::enable_if<std::is_integral<IntType>::value, int>::type = 0>
inline std::function<IntType()> Bnomial(IntType t = 1, double p = 0.5);

// Negative Binomial Distributions
template <class IntType = long, typename std::enable_if<std::is_integral<IntType>::value, int>::type = 0>
inline std::function<IntType()> NegativeBinomial(IntType k = 1, double p = 0.5);

// Geometric Distributions
template <class IntType = long, typename std::enable_if<std::is_integral<IntType>::value, int>::type = 0>
inline std::function<IntType()> Geometric(double p = 0.5);

// Poisson Distributions
template <class IntType = long, typename std::enable_if<std::is_integral<IntType>::value, int>::type = 0>
inline std::function<IntType()> Poisson(double mean = 1.0);

// Exponential Distributions
template <class RealType = double, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type = 0>
inline std::function<RealType()> Exponential(RealType lambda = 1.0);

// Gamma Distributions
template <class RealType = double, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type = 0>
inline std::function<RealType()> Gamma(RealType alpha = 1.0, RealType beta = 1.0);

// Weibull Distributions
template <class RealType = double, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type = 0>
inline std::function<RealType()> Weibull(RealType a = 1.0, RealType b = 1.0);

// ExtremeValue Distributions
template <class RealType = double, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type = 0>
inline std::function<RealType()> ExtremeValue(RealType a = 0.0, RealType b = 1.0);

// Normal Distributions
template <class RealType = double, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type = 0>
inline std::function<RealType()> Normal(RealType mean = 0.0, RealType stddev = 1.0);

// Lognormal Distributions
template <class RealType = double, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type = 0>
inline std::function<RealType()> Lognormal(RealType m = 0.0, RealType s = 1.0);

// ChiSquared Distributions
template <class RealType = double, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type = 0>
inline std::function<RealType()> ChiSquared(RealType n = 1);

// Cauchy Distributions
template <class RealType = double, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type = 0>
inline std::function<RealType()> Cauchy(RealType a = 0.0, RealType b = 1.0);

// FisherF Distributions
template <class RealType = double, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type = 0>
inline std::function<RealType()> FisherF(RealType m = 1, RealType n = 1);

// StudentT Distributions
template <class RealType = double, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type = 0>
inline std::function<RealType()> StudentT(RealType n = 1);


// Special Distributions
static std::function<double()> DecimalUniform = Uniform<double>(0.0, 1.0);
static std::function<double()> NormalGaussian = Normal<double>(0.0, 1.0);

} // namespace MGRndm


#endif // __CPPLibs_MGRndm_H__
