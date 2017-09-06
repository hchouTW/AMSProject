#ifndef __TRACKLibs_Math_H__
#define __TRACKLibs_Math_H__

// TODO
// include robust (efft_sigma is value/5 if value > 5 * maximum sigma)

namespace TrackSys {


class MultiGauss {
    public :
        MultiGauss(Double_t sgm);
        MultiGauss(Double_t wgt1, Double_t sgm1, Double_t wgt2, Double_t sgm2);
        MultiGauss(Double_t wgt1, Double_t sgm1, Double_t wgt2, Double_t sgm2, Double_t wgt3, Double_t sgm3);
        ~MultiGauss() { if (rand_func_ != nullptr) { delete rand_func_; rand_func_ = nullptr; } }

        inline Int_t num() const { return multi_gauss_.size(); }
        inline const Double_t& wgt(Int_t i) const { return multi_gauss_.at(i).first; }
        inline const Double_t& sgm(Int_t i) const { return multi_gauss_.at(i).second; }

        inline Double_t efft_sgm(Double_t r = 0.) const; 

        inline Double_t rndm();

    private :
        std::pair<Double_t, Double_t> bound_;
        std::vector<std::pair<Double_t, Double_t>> multi_gauss_;

    private :
        TF1 * rand_func_;

    private :
        static constexpr Double_t LMTL_PROB = 1.0e-6;
        static constexpr Double_t ROBUST = 5.0;
};


} // namesapce TrackSys


#endif // __TRACKLibs_Math_H__
