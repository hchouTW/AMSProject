#ifndef __TRACKLibs_Math_H__
#define __TRACKLibs_Math_H__

#include <TRandom.h>
#include <TF1.h>


namespace TrackSys {


namespace Rndm {
	TRandom3 rndmEngTR3(0);
	inline Double_t Landau   (Double_t mpv = 0.0, Double_t sigma = 1.0)  { return rndmEngTR3.Landau(mpv, sigma); }
}

// SVector SMatrix Package
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
}
	

// MultiGauss
// Gaussian Core f ~ ([0]/[1])*exp(-0.5*x*x/[1]/[1])
// Robust Method : efft_sigma = sigma * robust_factor if value/sigma > robust
class MultiGauss {
    public :
        enum class Opt { NOROBUST = 0, ROBUST = 1 };

    public :
        MultiGauss() : robust_(Opt::NOROBUST), rand_func_(nullptr) {}
        MultiGauss(Opt opt, long double sgm);
        MultiGauss(Opt opt, long double wgt1, long double sgm1, long double wgt2, long double sgm2);
        MultiGauss(Opt opt, long double wgt1, long double sgm1, long double wgt2, long double sgm2, long double wgt3, long double sgm3);
        MultiGauss(Opt opt, long double wgt1, long double sgm1, long double wgt2, long double sgm2, long double wgt3, long double sgm3, long double wgt4, long double sgm4);
        MultiGauss(Opt opt, long double wgt1, long double sgm1, long double wgt2, long double sgm2, long double wgt3, long double sgm3, long double wgt4, long double sgm4, long double wgt5, long double sgm5);
        ~MultiGauss() { robust_ = Opt::NOROBUST; if (rand_func_ != nullptr) { delete rand_func_; rand_func_ = nullptr; } }

        inline Int_t num() const { return multi_gauss_.size(); }
        inline const long double& wgt(Int_t i) const { return multi_gauss_.at(i).first; }
        inline const long double& sgm(Int_t i) const { return multi_gauss_.at(i).second; }

        inline long double efft_sgm(long double r = 0.) const; 

        inline long double rndm();

    private :
        std::pair<long double, long double> bound_;
        std::vector<std::pair<long double, long double>> multi_gauss_;

    private :
        Opt   robust_;
        TF1*  rand_func_;

    private :
        static TRandom* rndm_gen_;

        static constexpr Long64_t    NPX_ = 1000000;
        static constexpr long double LMTL_PROB_ = 1.0e-20;
        static constexpr long double ROBUST_SGM_ = 2.0;
};

TRandom* MultiGauss::rndm_gen_ = nullptr;





class IonEloss {
    public :
        IonEloss(const std::array<Double_t, 5>& mpv, const std::array<Double_t, 5>& sgm) : mpv_par_(mpv), sgm_par_(sgm) {}
        ~IonEloss() {}

        SVecD<2> operator() (Double_t adc, Double_t eta, Double_t dzds = 1.0) {
            Double_t invbtasqr = (1.0 + eta * eta);
            Double_t abseta    = std::fabs(eta);
            if (MGNumc::Compare(abseta, 1.0) < 0) return SVecD<2>();

            Double_t mpv = mpv_par_.at(0) * std::pow(invbtasqr, mpv_par_.at(2)) * 
                (mpv_par_.at(1) - 
                 std::pow(invbtasqr, -mpv_par_.at(1)) - 
                 std::log(mpv_par_.at(3) + std::pow(abseta, mpv_par_.at(4))));

            Double_t sgm = sgm_par_.at(0) * std::pow(invbtasqr, sgm_par_.at(2)) * 
                (sgm_par_.at(1) - 
                 std::pow(invbtasqr, -sgm_par_.at(1)) - 
                 std::log(sgm_par_.at(3) + std::pow(abseta, sgm_par_.at(4))));

            Double_t difmpv = (mpv_par_.at(0) * mpv_par_.at(2) * abseta * std::pow(invbtasqr, mpv_par_.at(2)-1.0)) *
                (2.0 * mpv_par_.at(1) -
                 2.0 * std::log(mpv_par_.at(3) + std::pow(abseta, mpv_par_.at(4))) -
                 (mpv_par_.at(4)/mpv_par_.at(2)) * invbtasqr * std::pow(abseta, mpv_par_.at(4)-2.0) / (mpv_par_.at(3) + std::pow(abseta, mpv_par_.at(4))));
            
            Double_t difsgm = (sgm_par_.at(0) * sgm_par_.at(2) * abseta * std::pow(invbtasqr, sgm_par_.at(2)-1.0)) *
                (2.0 * sgm_par_.at(1) -
                 2.0 * std::log(sgm_par_.at(3) + std::pow(abseta, sgm_par_.at(4))) -
                 (sgm_par_.at(4)/sgm_par_.at(2)) * invbtasqr * std::pow(abseta, sgm_par_.at(4)-2.0) / (sgm_par_.at(3) + std::pow(abseta, sgm_par_.at(4))));

            if (!MGNumc::Valid(mpv) || MGNumc::Compare(mpv) <= 0) mpv = 0.0;
            if (!MGNumc::Valid(sgm) || MGNumc::Compare(sgm) <= 0) sgm = 0.0;
            if (!MGNumc::Valid(difmpv) || MGNumc::Compare(difmpv) <= 0) difmpv = 0.0;
            if (!MGNumc::Valid(difsgm) || MGNumc::Compare(difsgm) <= 0) difsgm = 0.0;

            Double_t scl = std::fabs(dzds);
            Double_t res = (scl * adc - mpv) / sgm;
            Double_t dif = -difmpv / sgm;

            if (!MGNumc::Valid(res)) res = 0.0;
            if (!MGNumc::Valid(dif)) dif = 0.0;

            return SVecD<2>(res, dif);
        }

    private :
        std::array<Double_t, 5> mpv_par_;
        std::array<Double_t, 5> sgm_par_;
};


} // namesapce TrackSys


#endif // __TRACKLibs_Math_H__
