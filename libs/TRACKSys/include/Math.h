#ifndef __TRACKLibs_Math_H__
#define __TRACKLibs_Math_H__

#include <TMath.h>

#include <array>
#include <random>
#include <type_traits>
#include <functional>
namespace TrackSys {
namespace Numc {

// Constant Number
template<typename T = Double_t> constexpr T NEG    = static_cast<T>(-1);
template<typename T = Double_t> constexpr T ZERO   = static_cast<T>(0);
template<typename T = Double_t> constexpr T ONE    = static_cast<T>(1);
template<typename T = Double_t> constexpr T TWO    = static_cast<T>(2);
template<typename T = Double_t> constexpr T THREE  = static_cast<T>(3);
template<typename T = Double_t> constexpr T FOUR   = static_cast<T>(4);
template<typename T = Double_t> constexpr T FIVE   = static_cast<T>(5);
template<typename T = Double_t> constexpr T SIX    = static_cast<T>(6);
template<typename T = Double_t> constexpr T SEVEN  = static_cast<T>(7);
template<typename T = Double_t> constexpr T EIGHT  = static_cast<T>(8);
template<typename T = Double_t> constexpr T NINE   = static_cast<T>(9);
template<typename T = Double_t> constexpr T TEN    = static_cast<T>(10);

template<typename T = Double_t> constexpr T HUNDRED  = static_cast<T>(100);
template<typename T = Double_t> constexpr T THOUSAND = static_cast<T>(1000);

constexpr Double_t HALF = 5.00000000000000000e-01;

constexpr Double_t ONE_TO_TWO   = 5.00000000000000000e-01;
constexpr Double_t ONE_TO_SIX   = 1.66666666666666657e-01;
constexpr Double_t ONE_TO_EIGHT = 1.25000000000000000e-01;

constexpr Double_t PI          = 3.14159265358979312e+00;
constexpr Double_t INV_PI      = 3.18309886183790691e-01;
constexpr Double_t SQRT_PI     = 1.77245385090551588e+00;
constexpr Double_t INV_SQRT_PI = 5.64189583547756279e-01;

constexpr Double_t SQRT_TWO    = 1.41421356237309515e+00;
constexpr Double_t SQRT_THREE  = 1.73205080756887719e+00;

constexpr Double_t INV_SQRT_TWO    = 7.07106781186547462e-01;
constexpr Double_t INV_SQRT_THREE  = 5.77350269189625842e-01;

constexpr Double_t LOG_TWO = 6.93147180559945286e-01;
constexpr Double_t LOG_TEN = 2.30258509299404590e+00;

// Valid
template <class IntType = long, typename std::enable_if<std::is_integral<IntType>::value, int>::type = 0>
inline bool Valid(long long int var) { return (var == IntType(var)); }

template <class RealType = double, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type = 0>
inline bool Valid(RealType var) { int clf = std::fpclassify(var); return (clf == FP_NORMAL || clf == FP_ZERO); }

// Commpare
template <class IntType = long, typename std::enable_if<std::is_integral<IntType>::value, int>::type = 0>
inline short Compare(IntType a, IntType b = IntType(0)) {
	if (a == b)     return  0;
	else if (a > b) return  1;
	else            return -1;
}

template <class RealType = double, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type = 0>
inline short Compare(RealType a, RealType b = RealType(0.0)) {
	RealType diff = std::fabs(a - b);
    if (!std::isfinite(diff)) return 0;
	if (diff < std::numeric_limits<RealType>::epsilon() * 5.0e3) return  0;
	else if (a > b)                                              return  1;
	else                                                         return -1;
}

// EqualToZero
template <class IntType = long, typename std::enable_if<std::is_integral<IntType>::value, int>::type = 0>
inline bool EqualToZero(IntType a) { return (Compare(a) == 0); }

template <class RealType = double, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type = 0>
inline bool EqualToZero(RealType a) { return (Compare(a) == 0); }
        

// Norm quality from chi-square
Double_t NormQuality(Double_t nchi, Short_t ndof);
Double_t NormQuality(Double_t nchi, Double_t ndof);


} // namesapce Numc
} // namesapce TrackSys


#include <chrono>
#include <TRandom3.h>
namespace TrackSys {
namespace Rndm {

static TRandom3 rndmEngTR3(0);
inline Double_t Gaus(Double_t mean = 0.0, Double_t sigma = 1.0) { return rndmEngTR3.Gaus(mean, sigma); }
inline Double_t Landau(Double_t mpv = 0.0, Double_t sigma = 1.0) { return rndmEngTR3.Landau(mpv, sigma); }

// Random Engines
static long rndmSeed = std::chrono::system_clock::now().time_since_epoch().count();
static std::mt19937_64 rndmEngMT64(rndmSeed);

// Uniform Distributions
template <class IntType = long, typename std::enable_if<std::is_integral<IntType>::value, int>::type = 0>
inline std::function<IntType()> Uniform(IntType a = 0, IntType b = std::numeric_limits<IntType>::max()) {
	std::uniform_int_distribution<IntType> distribution(a, b);
	std::function<IntType()>&& rngfunc = std::bind(distribution, std::ref(rndmEngMT64));
	return rngfunc;
}

template <class RealType = double, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type = 0>
inline std::function<RealType()> Uniform(RealType a = 0.0, RealType b = std::numeric_limits<RealType>::max()) {
	std::uniform_real_distribution<RealType> distribution(a, b);
	std::function<RealType()>&& rngfunc = std::bind(distribution, std::ref(rndmEngMT64));
	return rngfunc;
}

// Normal Distributions
template <class RealType = double, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type = 0>
inline std::function<RealType()> Normal(RealType mean = 0.0, RealType stddev = 1.0) {
	std::normal_distribution<RealType> distribution(mean, stddev);
	std::function<RealType()>&& rngfunc = std::bind(distribution, std::ref(rndmEngMT64));
	return rngfunc;
}

// Gamma Distributions
template <class RealType = double, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type = 0>
inline std::function<RealType()> Gamma(RealType alpha = 1.0, RealType beta = 1.0) {
	std::gamma_distribution<RealType> distribution(alpha, beta);
	std::function<RealType()>&& rngfunc = std::bind(distribution, std::ref(rndmEngMT64));
	return rngfunc;
}

// Chi-Square Distributions
template <class RealType = double, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type = 0>
inline std::function<RealType()> ChiSquare(RealType k = 1.0) {
    std::chi_squared_distribution<RealType> distribution(k);
	std::function<RealType()>&& rngfunc = std::bind(distribution, std::ref(rndmEngMT64));
	return rngfunc;
}

// Special Distributions
static std::function<double()> DecimalUniform = Uniform<double>(0.0, 1.0); // [0 ~ 1]
static std::function<double()> NormalGaussian = Normal<double>(0.0, 1.0);

} // namesapce Rndm
} // namesapce TrackSys


// SVector SMatrix Package
#include <Math/SVector.h>
#include <Math/SMatrix.h>
namespace TrackSys {

template <unsigned int D>
using SVecO = ROOT::Math::SVector<Bool_t, D>;
template <unsigned int D>
using SVecS = ROOT::Math::SVector<Short_t, D>;
template <unsigned int D>
using SVecI = ROOT::Math::SVector<Int_t, D>;
template <unsigned int D>
using SVecD = ROOT::Math::SVector<Double_t, D>;

template <unsigned int D1, unsigned int D2 = D1>
using SMtxD = ROOT::Math::SMatrix<Double_t, D1, D2, ROOT::Math::MatRepStd<Double_t, D1, D2>>;
template <unsigned int D>
using SMtxSymD = ROOT::Math::SMatrix<Double_t, D, D, ROOT::Math::MatRepSym<Double_t, D>>;

using SMtxId = ROOT::Math::SMatrixIdentity;

namespace LA { // Linear Algebra
	//// SVector Template Function 
	using ROOT::Math::Cross;
	using ROOT::Math::Dot;
	using ROOT::Math::Mag;
	using ROOT::Math::Unit;
	//
	//// SMatrix Template Function 
	using ROOT::Math::Transpose;
	using ROOT::Math::Similarity;   // U   M U^T
	using ROOT::Math::SimilarityT;  // U^T M U
} // namespace LA

} // namesapce TrackSys
	
        
namespace TrackSys {

class LandauNumc {
    public :
        LandauNumc() {}
        ~LandauNumc() {}

        static long double EvalLn(long double x);
        static long double DevLn(long double x);
        static long double IcovLn(long double x);

    private :
        static constexpr long double LANDAU0   =  1.80655634e-01;
        static constexpr long double LANDAU0_X = -2.22782980e-01;
        static constexpr long double WIDTH_SCL =  1.17741002e+00; // sqrt( 2*log(2) )

        static constexpr long double BOUNDL = -3.0;
        static constexpr long double BOUNDU = 22.0;
        static constexpr long        NSET = 10000;
        static constexpr long double STEP = 0.0025;

        static std::vector<long double> data;
        
        static const std::array<long double, 5> LAND_DEV;
        static const std::array<long double, 8> LAND_ICOV;
};

const std::array<long double, 5> LandauNumc::LAND_DEV 
    { 6.19538e-02, 4.98613e-01, 1.26216e-01, -5.59660e-01, 1.17070e+00 };

const std::array<long double, 8> LandauNumc::LAND_ICOV 
    { 1.71213e-01, 1.05335e+00, 1.22843e-01, 2.44705e-01, 2.95471e-01, 5.75783e-01, 1.61325e-02, 6.61252e-02 };


class ApproxLnX {
    public :
        ApproxLnX() {}
        ~ApproxLnX() {}

        static long double Eval(long double x);
        static long double Dev(long double x);
        
        inline static long double X0() { return LnX0; }
        inline static long double EvalWithX0(long double x) { return (Eval(x) - LnX0); }
        static std::array<long double, 2> EvalAndDevWithX0(long double x);

    private :
        static constexpr long double HALF_WIDTH = 0.025;
        static constexpr long double LnX0  = -4.68887945411395;
        static constexpr long double DevX0 = 55.45177444479562;
};

} // namesapce TrackSys


namespace TrackSys {
// Robust
// S := nrm * nrm
// H := thres * thres  (H>=1)
// R := rate * 0.5 * (1 + erf(thres * log( abs(nrm/thres)) ) )     (0<=rate<=1)
// X := (S/H)
// A := log(1+X*X) / (X*X)
// chisq(S)  := S * A^(R/4)
// div1st(S) := A^(R/4) * (
//                  (1-R/2) + 
//                  (R/2) * (1/A) * (1/(1+X*X))
//              )
// S*div2nd(S) := (R/4 * X*X) * A^(1+R/4) * (
//                    (R-4) * (1/A)^3 * (1/X)^2 * (1/(1+X*X))^2 +
//                    (R-2) * (1/A) * (1/X)^2 -
//                    2*((R-1)*X*X + (R-3)) * (1/A)^2 * (1/X)^2 * (1/(1+X*X))^2
//                )

class Robust {
    public :
        enum class Opt { OFF = 0, ON = 1 };

        struct Option { 
            Option(Opt _opt = Opt::OFF, long double _thres = 1.0L, long double _rate = 1.0L) : opt(_opt), thres(_thres), rate(_rate) {}
            Opt         opt;
            long double thres; 
            long double rate; 
        };

    public :
        Robust(const Option& robust = Option(Opt::OFF)) : opt_(Opt::OFF), robust_(robust) {
            if (Opt::OFF == robust_.opt) return;
            Bool_t check_robust = (Opt::ON == robust_.opt && Numc::Compare(robust_.thres, Numc::ONE<long double>) >= 0 && Numc::Compare(robust_.rate) > 0 && Numc::Compare(robust_.rate, Numc::ONE<long double>) <= 0);
            if (!check_robust) robust_ = Option(Opt::OFF);
            opt_ = (Opt::ON == robust_.opt) ? Opt::ON : Opt::OFF;
        }
        ~Robust() {}

    public :
        const Opt& opt() const { return opt_; }
        std::array<long double, 3> minimizer(long double chi) const;
        
    private :
        Opt    opt_;
        Option robust_;
};

} // namesapce TrackSys


#include <TF1.h>
namespace TrackSys {
// MultiGaus
// Gaussian Core f ~ [0]*exp(-0.5*x*x/[1]/[1])
class MultiGaus {
    public :
        MultiGaus() : robust_(Robust()), rand_func_(nullptr), eftsgm_(Numc::ONE<long double>) {}
        MultiGaus(Robust robust, long double sgm);
        MultiGaus(Robust robust, long double wgt1, long double sgm1, long double wgt2, long double sgm2);
        MultiGaus(Robust robust, long double wgt1, long double sgm1, long double wgt2, long double sgm2, long double wgt3, long double sgm3);
        MultiGaus(Robust robust, long double wgt1, long double sgm1, long double wgt2, long double sgm2, long double wgt3, long double sgm3, long double wgt4, long double sgm4);
        MultiGaus(Robust robust, long double wgt1, long double sgm1, long double wgt2, long double sgm2, long double wgt3, long double sgm3, long double wgt4, long double sgm4, long double wgt5, long double sgm5);
        MultiGaus(Robust robust, long double wgt1, long double sgm1, long double wgt2, long double sgm2, long double wgt3, long double sgm3, long double wgt4, long double sgm4, long double wgt5, long double sgm5, long double wgt6, long double sgm6);
        MultiGaus(Robust robust, long double wgt1, long double sgm1, long double wgt2, long double sgm2, long double wgt3, long double sgm3, long double wgt4, long double sgm4, long double wgt5, long double sgm5, long double wgt6, long double sgm6, long double wgt7, long double sgm7);
        MultiGaus(Robust robust, long double wgt1, long double sgm1, long double wgt2, long double sgm2, long double wgt3, long double sgm3, long double wgt4, long double sgm4, long double wgt5, long double sgm5, long double wgt6, long double sgm6, long double wgt7, long double sgm7, long double wgt8, long double sgm8);
        MultiGaus(Robust robust, long double wgt1, long double sgm1, long double wgt2, long double sgm2, long double wgt3, long double sgm3, long double wgt4, long double sgm4, long double wgt5, long double sgm5, long double wgt6, long double sgm6, long double wgt7, long double sgm7, long double wgt8, long double sgm8, long double wgt9, long double sgm9);
        ~MultiGaus() { if (rand_func_ != nullptr) { delete rand_func_; rand_func_ = nullptr; } }

        inline const Robust& robust() const { return robust_; }

        inline Int_t num() const { return multi_gaus_.size(); }
        inline const long double& wgt(Int_t i = 0) const { return multi_gaus_.at(i).first; }
        inline const long double& sgm(Int_t i = 0) const { return multi_gaus_.at(i).second; }
        
        inline const long double eftsgm() const { return eftsgm_; }
        long double chi(long double r) const;

        std::array<long double, 3> minimizer(long double r = 0.) const; 
        long double rndm();

    protected :
        long double find_eftsgm() const;

    private :
        std::pair<long double, long double> bound_;
        std::vector<std::pair<long double, long double>> multi_gaus_;
        long double eftsgm_;

    private :
        Robust robust_;
        TF1*   rand_func_;

    private :
        static TRandom* rndm_gen_;

        static constexpr Long64_t    NPX = 1000000;
        static constexpr long double LMTL_PROB = 1.0e-20;
};

} // namesapce TrackSys


namespace TrackSys {
//TF1* flg = new TF1("flg", "[0] * TMath::Exp([1] * (-0.5)*((x-[2])*(x-[2])/[3]/[3]) + (1-[1]) * TMath::Log(TMath::Landau(1.17741002*(x-[2])/[3]-2.22782980e-01)/1.80655634e-01))");
//flg->SetParameters(1.0, 0.1, 0.0, 1.0);
//flg->SetParLimits(1, 0.0, 1.0);
class LandauGaus {
    public :
        LandauGaus(Robust robust, long double kpa, long double mpv, long double sgm, long double mode = Numc::ZERO<long double>, long double fluc = Numc::ZERO<long double>);
        ~LandauGaus() {}

        std::array<long double, 3> minimizer(long double x) const;

        inline const long double& kpa() const { return kpa_; }
        inline const long double& mpv() const { return mpv_; }
        inline const long double& sgm() const { return sgm_; }
        inline const long double& mode() const { return mode_; }
        inline const long double& fluc() const { return fluc_; }
        inline const long double& shft() const { return shft_; }
        
    protected :
        long double eval_norm(long double norm) const;
        std::array<long double, 2> eval_conv(long double norm) const; // (nrm, icov)

    protected :
        bool isfluc_;
        long double kpa_;
        long double mpv_;
        long double sgm_;
        long double mode_;
        long double fluc_;
        long double shft_;
        
        long double nrmfluc_;

    private :
        Robust robust_;

    private :
        static constexpr int GAUS_CONV_N = 26;
        static const std::array<long double, GAUS_CONV_N> GAUS_CONV_X;
        static const std::array<long double, GAUS_CONV_N> GAUS_CONV_P;

    private :
        std::array<long double, GAUS_CONV_N+1> convprob(long double norm) const;
};

const std::array<long double, LandauGaus::GAUS_CONV_N> LandauGaus::GAUS_CONV_X 
    { -2.5, -2.3, -2.1, -1.9, -1.7, -1.5, -1.3, -1.1, -0.9, -0.7, 
      -0.5, -0.3, -0.1,  0.1,  0.3,  0.5,  0.7,  0.9,  1.1,  1.3, 
       1.5,  1.7,  1.9,  2.1,  2.3,  2.5 };

const std::array<long double, LandauGaus::GAUS_CONV_N> LandauGaus::GAUS_CONV_P 
    { 0.0439, 0.0710, 0.1103, 0.1645, 0.2357, 0.3247, 0.4296, 0.5461, 0.6670, 0.7827,  
      0.8825, 0.9560, 0.9950, 0.9950, 0.9560, 0.8825, 0.7827, 0.6670, 0.5461, 0.4296,  
      0.3247, 0.2357, 0.1645, 0.1103, 0.0710, 0.0439 };

} // namesapce TrackSys


namespace TrackSys {
//TF1* flggm = new TF1("flggm", "[0] *  (1.0/sqrt(2.0*TMath::Pi())/[3]) * TMath::Exp((1-[1]) * TMath::Log(TMath::Landau((x-[2])/[3])/TMath::Landau(0)) + [1] * (-0.5)*((x-[2])*(x-[2])/[3]/[3])) + [4] * (pow([6], [5]) / TMath::Gamma([5])) * TMath::Power(x,[5]-1) * TMath::Exp(-[6]*x) * (0.5 * (TMath::Erf((x-[7])/[8])+1))", 0.1, 100);
class LgGeFunc {
    public :
        enum class Opt { NOROBUST = 0, ROBUST = 1 };
    
    public :
        LgGeFunc(Opt opt, long double ratio, long double lg_k, long double lg_m, long double lg_s, long double ge_a, long double ge_b, long double ge_m, long double ge_s) : robust_(opt), ratio_(ratio), lg_k_(lg_k), lg_m_(lg_m), lg_s_(lg_s), ge_a_(ge_a), ge_b_(ge_b), ge_m_(ge_m), ge_s_(ge_s) {}
        ~LgGeFunc() {}

        std::array<long double, 3> minimizer(long double x) const;

    protected :
        // transition
        long double get_ratio(long double x) const;
        long double get_lg(long double x) const;
        long double get_ge(long double x) const;

        // ionization
        std::array<long double, 2> lg_minimizer(long double x, long double wgt = Numc::ONE<long double>) const;
        long double lg_eval(long double norm) const;
        long double lg_div(long double norm) const;

    private :
        long double ratio_; // ratio := TR / (ION + TR)

        long double lg_k_;  // kpa
        long double lg_m_;  // mpv
        long double lg_s_;  // sgm

        long double ge_a_;  // alpha
        long double ge_b_;  // beta
        long double ge_m_;  // men
        long double ge_s_;  // sgm
    
    private :
        Opt robust_;
    
    private :
        static constexpr long double LANDAU0    =  1.80655634e-01;
        static constexpr long double LANDAU0_X  = -2.22782980e-01;
        static constexpr long double DELTA      = 0.01;
        static constexpr long double ROBUST_SGM = 2.0;
};

} // namesapce TrackSys


#endif // __TRACKLibs_Math_H__
