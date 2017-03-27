#ifndef __MgntRndm_H__
#define __MgntRndm_H__
#include <random>
#include <type_traits>
#include <functional>

namespace MgntRndm {
	// Random Engines
	static long RndmSeed = std::chrono::system_clock::now().time_since_epoch().count();
	static std::mt19937_64 RndmEngMT64(MgntRndm::RndmSeed);
	
	// Uniform Distributions
	template <class IntType = long, typename std::enable_if<std::is_integral<IntType>::value, int>::type = 0>
	inline std::function<IntType()> Uniform(IntType a = 0, IntType b = std::numeric_limits<IntType>::max()) {
		std::uniform_int_distribution<IntType> distribution(a, b);
		std::function<IntType()>&& rngfunc = std::bind(distribution, std::ref(MgntRndm::RndmEngMT64));
		return rngfunc;
	}
	
	template <class RealType = double, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type = 0>
	inline std::function<RealType()> Uniform(RealType a = 0.0, RealType b = std::numeric_limits<RealType>::max()) {
		std::uniform_real_distribution<RealType> distribution(a, b);
		std::function<RealType()>&& rngfunc = std::bind(distribution, std::ref(MgntRndm::RndmEngMT64));
		return rngfunc;
	}

	// Bernoulli Distributions
	inline std::function<bool()> Bernoulli(double p = 0.5) {
		std::bernoulli_distribution distribution(p);
		std::function<bool()>&& rngfunc = std::bind(distribution, std::ref(MgntRndm::RndmEngMT64));
		return rngfunc;
	}
	
	template <class IntType = long, typename std::enable_if<std::is_integral<IntType>::value, int>::type = 0>
	inline std::function<IntType()> Bnomial(IntType t = 1, double p = 0.5) {
		std::binomial_distribution<IntType> distribution(t, p);
		std::function<IntType()>&& rngfunc = std::bind(distribution, std::ref(MgntRndm::RndmEngMT64));
		return rngfunc;
	}

	template <class IntType = long, typename std::enable_if<std::is_integral<IntType>::value, int>::type = 0>
	inline std::function<IntType()> NegativeBinomial(IntType k = 1, double p = 0.5) {
		std::negative_binomial_distribution<IntType> distribution(k, p);
		std::function<IntType()>&& rngfunc = std::bind(distribution, std::ref(MgntRndm::RndmEngMT64));
		return rngfunc;
	}

	template <class IntType = long, typename std::enable_if<std::is_integral<IntType>::value, int>::type = 0>
	inline std::function<IntType()> Geometric(double p = 0.5) {
		std::geometric_distribution<IntType> distribution(p);
		std::function<IntType()>&& rngfunc = std::bind(distribution, std::ref(MgntRndm::RndmEngMT64));
		return rngfunc;
	}

	// Poisson Distributions
	template <class IntType = long, typename std::enable_if<std::is_integral<IntType>::value, int>::type = 0>
	inline std::function<IntType()> Poisson(double mean = 1.0) {
		std::poisson_distribution<IntType> distribution(mean);
		std::function<IntType()>&& rngfunc = std::bind(distribution, std::ref(MgntRndm::RndmEngMT64));
		return rngfunc;
	}	

	template <class RealType = double, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type = 0>
	inline std::function<RealType()> Exponential(RealType lambda = 1.0) {
		std::exponential_distribution<RealType> distribution(lambda);
		std::function<RealType()>&& rngfunc = std::bind(distribution, std::ref(MgntRndm::RndmEngMT64));
		return rngfunc;
	}	

	template <class RealType = double, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type = 0>
	inline std::function<RealType()> Gamma(RealType alpha = 1.0, RealType beta = 1.0) {
		std::gamma_distribution<RealType> distribution(alpha, beta);
		std::function<RealType()>&& rngfunc = std::bind(distribution, std::ref(MgntRndm::RndmEngMT64));
		return rngfunc;
	}	

	template <class RealType = double, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type = 0>
	inline std::function<RealType()> Weibull(RealType a = 1.0, RealType b = 1.0) {
		std::weibull_distribution<RealType> distribution(a, b);
		std::function<RealType()>&& rngfunc = std::bind(distribution, std::ref(MgntRndm::RndmEngMT64));
		return rngfunc;
	}	

	template <class RealType = double, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type = 0>
	inline std::function<RealType()> ExtremeValue(RealType a = 0.0, RealType b = 1.0) {
		std::extreme_value_distribution<RealType> distribution(a, b);
		std::function<RealType()>&& rngfunc = std::bind(distribution, std::ref(MgntRndm::RndmEngMT64));
		return rngfunc;
	}

	// Normal Distributions
	template <class RealType = double, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type = 0>
	inline std::function<RealType()> Normal(RealType mean = 0.0, RealType stddev = 1.0) {
		std::normal_distribution<RealType> distribution(mean, stddev);
		std::function<RealType()>&& rngfunc = std::bind(distribution, std::ref(MgntRndm::RndmEngMT64));
		return rngfunc;
	}
		
	template <class RealType = double, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type = 0>
	inline std::function<RealType()> Lognormal(RealType m = 0.0, RealType s = 1.0) {
		std::lognormal_distribution<RealType> distribution(m, s);
		std::function<RealType()>&& rngfunc = std::bind(distribution, std::ref(MgntRndm::RndmEngMT64));
		return rngfunc;
	}
		
	template <class RealType = double, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type = 0>
	inline std::function<RealType()> ChiSquared(RealType n = 1) {
		std::chi_squared_distribution<RealType> distribution(n);
		std::function<RealType()>&& rngfunc = std::bind(distribution, std::ref(MgntRndm::RndmEngMT64));
		return rngfunc;
	}
		
	template <class RealType = double, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type = 0>
	inline std::function<RealType()> Cauchy(RealType a = 0.0, RealType b = 1.0) {
		std::cauchy_distribution<RealType> distribution(a, b);
		std::function<RealType()>&& rngfunc = std::bind(distribution, std::ref(MgntRndm::RndmEngMT64));
		return rngfunc;
	}
		
	template <class RealType = double, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type = 0>
	inline std::function<RealType()> FisherF(RealType m = 1, RealType n = 1) {
		std::fisher_f_distribution<RealType> distribution(m, n);
		std::function<RealType()>&& rngfunc = std::bind(distribution, std::ref(MgntRndm::RndmEngMT64));
		return rngfunc;
	}
		
	template <class RealType = double, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type = 0>
	inline std::function<RealType()> StudentT(RealType n = 1) {
		std::student_t_distribution<RealType> distribution(n);
		std::function<RealType()>&& rngfunc = std::bind(distribution, std::ref(MgntRndm::RndmEngMT64));
		return rngfunc;
	}
}


#endif // __MgntRndm_H__
