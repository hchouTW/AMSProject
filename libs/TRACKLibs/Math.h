#ifndef __TRACKLibs_Math_H__
#define __TRACKLibs_Math_H__

#include <TRandom.h>
#include <TF1.h>


namespace TrackSys {


// MultiGauss
// Gaussian Core f ~ ([0]/[1])*exp(-0.5*x*x/[1]/[1])
// Robust Method : efft_sigma = sigma * robust_factor if value/sigma > robust
class MultiGauss {
    public :
        enum class Opt { NOROBUST = 0, ROBUST = 1 };

    public :
        MultiGauss() : robust_(Opt::NOROBUST), rand_func_(nullptr) {}
        MultiGauss(Opt opt, Double_t sgm);
        MultiGauss(Opt opt, Double_t wgt1, Double_t sgm1, Double_t wgt2, Double_t sgm2);
        MultiGauss(Opt opt, Double_t wgt1, Double_t sgm1, Double_t wgt2, Double_t sgm2, Double_t wgt3, Double_t sgm3);
        MultiGauss(Opt opt, Double_t wgt1, Double_t sgm1, Double_t wgt2, Double_t sgm2, Double_t wgt3, Double_t sgm3, Double_t wgt4, Double_t sgm4);
        MultiGauss(Opt opt, Double_t wgt1, Double_t sgm1, Double_t wgt2, Double_t sgm2, Double_t wgt3, Double_t sgm3, Double_t wgt4, Double_t sgm4, Double_t wgt5, Double_t sgm5);
        ~MultiGauss() { robust_ = Opt::NOROBUST; if (rand_func_ != nullptr) { delete rand_func_; rand_func_ = nullptr; } }

        inline Int_t num() const { return multi_gauss_.size(); }
        inline const long double& wgt(Int_t i) const { return multi_gauss_.at(i).first; }
        inline const long double& sgm(Int_t i) const { return multi_gauss_.at(i).second; }

        inline Double_t efft_sgm(Double_t r = 0.) const; 

        inline Double_t rndm();

    private :
        std::pair<long double, long double> bound_;
        std::vector<std::pair<long double, long double>> multi_gauss_;

    private :
        Opt   robust_;
        TF1*  rand_func_;

    private :
        static TRandom* rndm_gen_;

        static constexpr Long64_t    NPX_ = 100000;
        static constexpr long double LMTL_PROB_ = 1.0e-8;
        static constexpr long double ROBUST_ = 2.0;
};

TRandom* MultiGauss::rndm_gen_ = nullptr;


} // namesapce TrackSys


#endif // __TRACKLibs_Math_H__
