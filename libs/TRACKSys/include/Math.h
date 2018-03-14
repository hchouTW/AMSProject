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
template<typename T = Double_t> constexpr T SIX    = static_cast<T>(6);
template<typename T = Double_t> constexpr T EIGHT  = static_cast<T>(8);
template<typename T = Double_t> constexpr T TEN    = static_cast<T>(10);

constexpr Double_t HALF = 5.00000000000000000e-01;

constexpr Double_t ONE_TO_TWO   = 5.00000000000000000e-01;
constexpr Double_t ONE_TO_SIX   = 1.66666666666666657e-01;
constexpr Double_t ONE_TO_EIGHT = 1.25000000000000000e-01;

constexpr Double_t PI     = 3.14159265358979312e+00;
constexpr Double_t INV_PI = 3.18309886183790691e-01;

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
inline short Compare(IntType a, IntType b = IntType(0));

template <class RealType = double, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type = 0>
inline short Compare(RealType a, RealType b = RealType(0.0));

// EqualToZero
template <class IntType = long, typename std::enable_if<std::is_integral<IntType>::value, int>::type = 0>
inline bool EqualToZero(IntType a) { return (Compare(a) == 0); }

template <class RealType = double, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type = 0>
inline bool EqualToZero(RealType a) { return (Compare(a) == 0); }

} // namesapce Numc
} // namesapce TrackSys


#include <chrono>
#include <TRandom3.h>
namespace TrackSys {
namespace Rndm {

TRandom3 rndmEngTR3(0);
inline Double_t Gaus(Double_t mean = 0.0, Double_t sigma = 1.0) { return rndmEngTR3.Gaus(mean, sigma); }
inline Double_t Landau(Double_t mpv = 0.0, Double_t sigma = 1.0) { return rndmEngTR3.Landau(mpv, sigma); }

// Random Engines
static long rndmSeed = std::chrono::system_clock::now().time_since_epoch().count();
static std::mt19937_64 rndmEngMT64(rndmSeed);

// Uniform Distributions
template <class IntType = long, typename std::enable_if<std::is_integral<IntType>::value, int>::type = 0>
inline std::function<IntType()> Uniform(IntType a = 0, IntType b = std::numeric_limits<IntType>::max());

template <class RealType = double, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type = 0>
inline std::function<RealType()> Uniform(RealType a = 0.0, RealType b = std::numeric_limits<RealType>::max());

// Normal Distributions
template <class RealType = double, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type = 0>
inline std::function<RealType()> Normal(RealType mean = 0.0, RealType stddev = 1.0);

// Gamma Distributions
template <class RealType = double, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type = 0>
inline std::function<RealType()> Gamma(RealType alpha = 1.0, RealType beta = 1.0);

// Special Distributions
static std::function<double()> DecimalUniform = Uniform<double>(0.0, 1.0);
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
	

#include <TF1.h>
namespace TrackSys {
// MultiGaus
// Gaussian Core f ~ ([0]/[1])*exp(-0.5*x*x/[1]/[1])
// Robust Method : efft_sigma = sigma * robust_factor if value/sigma > robust
class MultiGaus {
    public :
        enum class Opt { NOROBUST = 0, ROBUST = 1 };

    public :
        MultiGaus() : robust_(Opt::NOROBUST), rand_func_(nullptr) {}
        MultiGaus(Opt opt, long double sgm);
        MultiGaus(Opt opt, long double wgt1, long double sgm1, long double wgt2, long double sgm2);
        MultiGaus(Opt opt, long double wgt1, long double sgm1, long double wgt2, long double sgm2, long double wgt3, long double sgm3);
        MultiGaus(Opt opt, long double wgt1, long double sgm1, long double wgt2, long double sgm2, long double wgt3, long double sgm3, long double wgt4, long double sgm4);
        MultiGaus(Opt opt, long double wgt1, long double sgm1, long double wgt2, long double sgm2, long double wgt3, long double sgm3, long double wgt4, long double sgm4, long double wgt5, long double sgm5);
        ~MultiGaus() { robust_ = Opt::NOROBUST; if (rand_func_ != nullptr) { delete rand_func_; rand_func_ = nullptr; } }

        inline Int_t num() const { return multi_gaus_.size(); }
        inline const long double& wgt(Int_t i) const { return multi_gaus_.at(i).first; }
        inline const long double& sgm(Int_t i) const { return multi_gaus_.at(i).second; }

        inline long double efft_sgm(long double r = 0.) const; 

        inline long double rndm();

    private :
        std::pair<long double, long double> bound_;
        std::vector<std::pair<long double, long double>> multi_gaus_;

    private :
        Opt   robust_;
        TF1*  rand_func_;

    private :
        static TRandom* rndm_gen_;

        static constexpr Long64_t    NPX_ = 1000000;
        static constexpr long double LMTL_PROB_ = 1.0e-20;
        static constexpr long double ROBUST_SGM_ = 2.0;
};

TRandom* MultiGaus::rndm_gen_ = nullptr;
} // namesapce TrackSys


namespace TrackSys {
//TF1* flg = new TF1("flg", "[0] * TMath::Exp( (1-[1]) * TMath::Log(TMath::Landau((x-[2])/[3])/TMath::Landau(0)) + [1] * (-0.5)*((x-[2])*(x-[2])/[3]/[3]) )");
//flg->SetParameters(1.0, 0.1, 0.0, 1.0);
//flg->SetParLimits(1, 0.0, 1.0);
class LandauGaus {
    public :
        LandauGaus(long double kpa, long double mpv, long double sgm);
        ~LandauGaus() {}

        inline std::array<long double, 2> operator() (long double x) const { return std::array<long double, 2>({ eval(x), div(x) }); }

        inline const long double& kpa() const { return kpa_; }
        inline const long double& mpv() const { return mpv_; }
        inline const long double& sgm() const { return sgm_; }

    protected :
        long double eval(long double x) const;
        long double div(long double x) const;

    protected :
        long double kpa_;
        long double mpv_;
        long double sgm_;

    private :
        static constexpr long double LANDAU0_ = 1.78854160900000003e-01;
};

} // namesapce TrackSys


namespace TrackSys {
//TF1* fmpv = new TF1("fmpv", "[0] * (1+x*x)^[2] * ([1] - (1+x*x)^(-[2]) - TMath::Log([3] + abs(x)^[4]))");
//fmpv->SetParameters(10, 6.5, 1.0, 10.0, 1.0);
//TF1* fkpa = new TF1("fkpa", "[0] * (1+x*x)^[1] * (1 + [2]*abs(x)^[3] - TMath::Log([4] + abs(x)^[5]))");
//fkpa->SetParameters(10, 1.0, 1.5, 3.0, 1.0, 7.0);
class IonEloss {
    public :
        IonEloss(const std::array<long double, 6>& kpa, const std::array<long double, 5>& mpv, const std::array<long double, 5>& sgm, long double fluc = Numc::ZERO<long double>) : kpa_(kpa), mpv_(mpv), sgm_(sgm), fluc_(fluc) { if (Numc::Compare(fluc_) <= 0) fluc_ = Numc::ZERO<long double>; }
        ~IonEloss() {}
        
        inline SVecD<2> operator() (long double x, long double eta) const { std::array<long double, 2>&& ion = eval(x, eta); return SVecD<2>(ion.at(0), ion.at(1)); }
   
    protected :
        std::array<long double, 2> eval(long double x, long double eta) const;

        inline long double eval_kpa(long double eta, long double ibsqr) const;
        inline long double eval_mpv(long double eta, long double ibsqr) const;
        inline long double eval_sgm(long double eta, long double ibsqr) const;
        
        inline long double eval_divmpv(long double eta, long double ibsqr) const;

    private :
        std::array<long double, 6> kpa_;
        std::array<long double, 5> mpv_;
        std::array<long double, 5> sgm_;
        long double                fluc_;
};

} // namesapce TrackSys


#endif // __TRACKLibs_Math_H__
