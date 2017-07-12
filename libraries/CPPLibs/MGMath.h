#ifndef __CPPLibs_MGMath_H__
#define __CPPLibs_MGMath_H__

#include <cmath>
#include <TROOT.h> 

namespace MGMath {

constexpr Double_t NEG    = -1.0;
constexpr Double_t ZERO   =  0.0;
constexpr Double_t ONE    =  1.0;
constexpr Double_t TWO    =  2.0;
constexpr Double_t THREE  =  3.0;
constexpr Double_t FOUR   =  4.0;
constexpr Double_t FIVE   =  5.0;
constexpr Double_t SIX    =  6.0;
constexpr Double_t SEVEN  =  7.0;
constexpr Double_t EIGHT  =  8.0;
constexpr Double_t NINE   =  9.0;
constexpr Double_t TEN    = 10.0;
constexpr Double_t TWELVE = 12.0;
constexpr Double_t ONE_TO_TWO   = ONE / TWO  ;
constexpr Double_t ONE_TO_THREE = ONE / THREE;
constexpr Double_t ONE_TO_FOUR  = ONE / FOUR ;
constexpr Double_t ONE_TO_FIVE  = ONE / FIVE ;
constexpr Double_t ONE_TO_SIX   = ONE / SIX  ;
constexpr Double_t ONE_TO_SEVEN = ONE / SEVEN;
constexpr Double_t ONE_TO_EIGHT = ONE / EIGHT;
constexpr Double_t ONE_TO_NINE  = ONE / NINE ;
constexpr Double_t ONE_TO_TEN   = ONE / TEN  ;
constexpr Double_t HALF = ONE_TO_TWO;

const Double_t PI          = std::atan(ONE) * FOUR;
const Double_t SQR_PI      = PI * PI;
const Double_t SQRT_PI     = std::sqrt(PI);
const Double_t SQRT_TWO_PI = std::sqrt(TWO * PI);

const Double_t INV_PI          = ONE / PI;
const Double_t INV_SQR_PI      = ONE / SQR_PI;
const Double_t INV_SQRT_PI     = ONE / SQRT_PI;
const Double_t INV_SQRT_TWO_PI = ONE / SQRT_TWO_PI;

const Double_t SQRT_TWO    = std::sqrt(TWO);
const Double_t SQRT_THREE  = std::sqrt(THREE);
const Double_t SQRT_TWELVE = std::sqrt(TWELVE);

const Double_t INV_SQRT_TWO    = ONE / SQRT_TWO;
const Double_t INV_SQRT_THREE  = ONE / SQRT_THREE;
const Double_t INV_SQRT_TWELVE = ONE / SQRT_TWELVE;

const Double_t EXP     = std::exp(ONE);
const Double_t LOG_TWO = std::log(TWO);
const Double_t LOG_TEN = std::log(TEN);


inline double Quantile(double prob);
inline double Quantile(double prob, double sgm) { return (sgm * Quantile(prob)); }
inline double Quantile(double prob, double mean, double sgm) { return (mean + sgm * Quantile(prob)); }
inline double Quantile(double prob, double wgtA, double sgmA, double wgtB, double sgmB);


inline double Erf(double x) { return std::erf(x); }
inline double Erfc(double x) { return std::erfc(x); }

inline double InvErfc(double x) { return -INV_SQRT_TWO * Quantile(HALF * x); }
inline double InvErf(double x) { return INV_SQRT_TWO * Quantile(HALF * (x + ONE)); }


}


#endif // __CPPLibs_MGMath_H__
