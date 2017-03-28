#ifndef __MGRndm_C__
#define __MGRndm_C__

#include "MGRndm.h"


template <class IntType = long, typename std::enable_if<std::is_integral<IntType>::value, int>::type = 0>
std::function<IntType()> MGRndm::Uniform(IntType a, IntType b) {
	std::uniform_int_distribution<IntType> distribution(a, b);
	std::function<IntType()>&& rngfunc = std::bind(distribution, std::ref(RndmEngMT64));
	return rngfunc;
}


template <class RealType = double, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type = 0>
std::function<RealType()> MGRndm::Uniform(RealType a, RealType b) {
	std::uniform_real_distribution<RealType> distribution(a, b);
	std::function<RealType()>&& rngfunc = std::bind(distribution, std::ref(RndmEngMT64));
	return rngfunc;
}


std::function<bool()> MGRndm::Bernoulli(double p) {
	std::bernoulli_distribution distribution(p);
	std::function<bool()>&& rngfunc = std::bind(distribution, std::ref(RndmEngMT64));
	return rngfunc;
}


template <class IntType = long, typename std::enable_if<std::is_integral<IntType>::value, int>::type = 0>
std::function<IntType()> MGRndm::Bnomial(IntType t, double p) {
	std::binomial_distribution<IntType> distribution(t, p);
	std::function<IntType()>&& rngfunc = std::bind(distribution, std::ref(RndmEngMT64));
	return rngfunc;
}


template <class IntType = long, typename std::enable_if<std::is_integral<IntType>::value, int>::type = 0>
std::function<IntType()> MGRndm::NegativeBinomial(IntType k, double p) {
	std::negative_binomial_distribution<IntType> distribution(k, p);
	std::function<IntType()>&& rngfunc = std::bind(distribution, std::ref(RndmEngMT64));
	return rngfunc;
}


template <class IntType = long, typename std::enable_if<std::is_integral<IntType>::value, int>::type = 0>
std::function<IntType()> MGRndm::Geometric(double p) {
	std::geometric_distribution<IntType> distribution(p);
	std::function<IntType()>&& rngfunc = std::bind(distribution, std::ref(RndmEngMT64));
	return rngfunc;
}


template <class IntType = long, typename std::enable_if<std::is_integral<IntType>::value, int>::type = 0>
std::function<IntType()> MGRndm::Poisson(double mean) {
	std::poisson_distribution<IntType> distribution(mean);
	std::function<IntType()>&& rngfunc = std::bind(distribution, std::ref(RndmEngMT64));
	return rngfunc;
}	


template <class RealType = double, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type = 0>
std::function<RealType()> MGRndm::Exponential(RealType lambda) {
	std::exponential_distribution<RealType> distribution(lambda);
	std::function<RealType()>&& rngfunc = std::bind(distribution, std::ref(RndmEngMT64));
	return rngfunc;
}	


template <class RealType = double, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type = 0>
std::function<RealType()> MGRndm::Gamma(RealType alpha, RealType beta) {
	std::gamma_distribution<RealType> distribution(alpha, beta);
	std::function<RealType()>&& rngfunc = std::bind(distribution, std::ref(RndmEngMT64));
	return rngfunc;
}	


template <class RealType = double, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type = 0>
std::function<RealType()> MGRndm::Weibull(RealType a, RealType b) {
	std::weibull_distribution<RealType> distribution(a, b);
	std::function<RealType()>&& rngfunc = std::bind(distribution, std::ref(RndmEngMT64));
	return rngfunc;
}	


template <class RealType = double, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type = 0>
std::function<RealType()> MGRndm::ExtremeValue(RealType a, RealType b) {
	std::extreme_value_distribution<RealType> distribution(a, b);
	std::function<RealType()>&& rngfunc = std::bind(distribution, std::ref(RndmEngMT64));
	return rngfunc;
}


template <class RealType = double, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type = 0>
std::function<RealType()> MGRndm::Normal(RealType mean, RealType stddev) {
	std::normal_distribution<RealType> distribution(mean, stddev);
	std::function<RealType()>&& rngfunc = std::bind(distribution, std::ref(RndmEngMT64));
	return rngfunc;
}
	

template <class RealType = double, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type = 0>
std::function<RealType()> MGRndm::Lognormal(RealType m, RealType s) {
	std::lognormal_distribution<RealType> distribution(m, s);
	std::function<RealType()>&& rngfunc = std::bind(distribution, std::ref(RndmEngMT64));
	return rngfunc;
}
	

template <class RealType = double, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type = 0>
std::function<RealType()> MGRndm::ChiSquared(RealType n) {
	std::chi_squared_distribution<RealType> distribution(n);
	std::function<RealType()>&& rngfunc = std::bind(distribution, std::ref(RndmEngMT64));
	return rngfunc;
}
	

template <class RealType = double, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type = 0>
std::function<RealType()> MGRndm::Cauchy(RealType a, RealType b) {
	std::cauchy_distribution<RealType> distribution(a, b);
	std::function<RealType()>&& rngfunc = std::bind(distribution, std::ref(RndmEngMT64));
	return rngfunc;
}
	

template <class RealType = double, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type = 0>
std::function<RealType()> MGRndm::FisherF(RealType m, RealType n) {
	std::fisher_f_distribution<RealType> distribution(m, n);
	std::function<RealType()>&& rngfunc = std::bind(distribution, std::ref(RndmEngMT64));
	return rngfunc;
}
	

template <class RealType = double, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type = 0>
std::function<RealType()> MGRndm::StudentT(RealType n) {
	std::student_t_distribution<RealType> distribution(n);
	std::function<RealType()>&& rngfunc = std::bind(distribution, std::ref(RndmEngMT64));
	return rngfunc;
}


#endif // __MGRndm_C__
