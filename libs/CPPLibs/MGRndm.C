#ifndef __CPPLibs_MGRndm_C__
#define __CPPLibs_MGRndm_C__

#include "MGRndm.h"

namespace MGRndm {

template <class IntType, typename std::enable_if<std::is_integral<IntType>::value, int>::type>
std::function<IntType()> Uniform(IntType a, IntType b) {
	std::uniform_int_distribution<IntType> distribution(a, b);
	std::function<IntType()>&& rngfunc = std::bind(distribution, std::ref(rndmEngMT64));
	return rngfunc;
}


template <class RealType, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type>
std::function<RealType()> Uniform(RealType a, RealType b) {
	std::uniform_real_distribution<RealType> distribution(a, b);
	std::function<RealType()>&& rngfunc = std::bind(distribution, std::ref(rndmEngMT64));
	return rngfunc;
}


std::function<bool()> Bernoulli(double p) {
	std::bernoulli_distribution distribution(p);
	std::function<bool()>&& rngfunc = std::bind(distribution, std::ref(rndmEngMT64));
	return rngfunc;
}


template <class IntType, typename std::enable_if<std::is_integral<IntType>::value, int>::type>
std::function<IntType()> Bnomial(IntType t, double p) {
	std::binomial_distribution<IntType> distribution(t, p);
	std::function<IntType()>&& rngfunc = std::bind(distribution, std::ref(rndmEngMT64));
	return rngfunc;
}


template <class IntType, typename std::enable_if<std::is_integral<IntType>::value, int>::type>
std::function<IntType()> NegativeBinomial(IntType k, double p) {
	std::negative_binomial_distribution<IntType> distribution(k, p);
	std::function<IntType()>&& rngfunc = std::bind(distribution, std::ref(rndmEngMT64));
	return rngfunc;
}


template <class IntType, typename std::enable_if<std::is_integral<IntType>::value, int>::type>
std::function<IntType()> Geometric(double p) {
	std::geometric_distribution<IntType> distribution(p);
	std::function<IntType()>&& rngfunc = std::bind(distribution, std::ref(rndmEngMT64));
	return rngfunc;
}


template <class IntType, typename std::enable_if<std::is_integral<IntType>::value, int>::type>
std::function<IntType()> Poisson(double mean) {
	std::poisson_distribution<IntType> distribution(mean);
	std::function<IntType()>&& rngfunc = std::bind(distribution, std::ref(rndmEngMT64));
	return rngfunc;
}	


template <class RealType, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type>
std::function<RealType()> Exponential(RealType lambda) {
	std::exponential_distribution<RealType> distribution(lambda);
	std::function<RealType()>&& rngfunc = std::bind(distribution, std::ref(rndmEngMT64));
	return rngfunc;
}	


template <class RealType, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type>
std::function<RealType()> Gamma(RealType alpha, RealType beta) {
	std::gamma_distribution<RealType> distribution(alpha, beta);
	std::function<RealType()>&& rngfunc = std::bind(distribution, std::ref(rndmEngMT64));
	return rngfunc;
}	


template <class RealType, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type>
std::function<RealType()> Weibull(RealType a, RealType b) {
	std::weibull_distribution<RealType> distribution(a, b);
	std::function<RealType()>&& rngfunc = std::bind(distribution, std::ref(rndmEngMT64));
	return rngfunc;
}	


template <class RealType, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type>
std::function<RealType()> ExtremeValue(RealType a, RealType b) {
	std::extreme_value_distribution<RealType> distribution(a, b);
	std::function<RealType()>&& rngfunc = std::bind(distribution, std::ref(rndmEngMT64));
	return rngfunc;
}


template <class RealType, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type>
std::function<RealType()> Normal(RealType mean, RealType stddev) {
	std::normal_distribution<RealType> distribution(mean, stddev);
	std::function<RealType()>&& rngfunc = std::bind(distribution, std::ref(rndmEngMT64));
	return rngfunc;
}
	

template <class RealType, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type>
std::function<RealType()> Lognormal(RealType m, RealType s) {
	std::lognormal_distribution<RealType> distribution(m, s);
	std::function<RealType()>&& rngfunc = std::bind(distribution, std::ref(rndmEngMT64));
	return rngfunc;
}
	

template <class RealType, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type>
std::function<RealType()> ChiSquared(RealType n) {
	std::chi_squared_distribution<RealType> distribution(n);
	std::function<RealType()>&& rngfunc = std::bind(distribution, std::ref(rndmEngMT64));
	return rngfunc;
}
	

template <class RealType, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type>
std::function<RealType()> Cauchy(RealType a, RealType b) {
	std::cauchy_distribution<RealType> distribution(a, b);
	std::function<RealType()>&& rngfunc = std::bind(distribution, std::ref(rndmEngMT64));
	return rngfunc;
}
	

template <class RealType, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type>
std::function<RealType()> FisherF(RealType m, RealType n) {
	std::fisher_f_distribution<RealType> distribution(m, n);
	std::function<RealType()>&& rngfunc = std::bind(distribution, std::ref(rndmEngMT64));
	return rngfunc;
}
	

template <class RealType, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type>
std::function<RealType()> StudentT(RealType n) {
	std::student_t_distribution<RealType> distribution(n);
	std::function<RealType()>&& rngfunc = std::bind(distribution, std::ref(rndmEngMT64));
	return rngfunc;
}

} // namespace MGRndm

#endif // __CPPLibs_MGRndm_C__
