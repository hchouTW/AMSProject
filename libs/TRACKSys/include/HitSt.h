#ifndef __TRACKLibs_HitSt_H__
#define __TRACKLibs_HitSt_H__


namespace TrackSys {

class HitSt {
    public :
        HitSt(Bool_t sx = false, Bool_t sy = false, Int_t lay = 0) :
            seqID_(-1), seqIDcx_(-1), seqIDcy_(-1), seqIDex_(-1), seqIDey_(-1),
            lay_(lay), side_(sx, sy),
            nsr_(0, 0),
            coo_(0., 0., 0.),
            err_(DEFERR_X_, DEFERR_Y_),
            adc_(0., 0.),
            type_(PartType::Proton),
            pdf_cx_(&PDF_PR_CX_NN_), pdf_cy_(&PDF_PR_CY_NN_), 
            pdf_ex_(&PDF_PR_EX_), pdf_ey_(&PDF_PR_EY_) {}
        
        void print() const;
        
        inline Bool_t operator()() const { return (side_(0) || side_(1)); }
        
        inline Short_t set_seqID(Short_t id); 

        inline void set_coo(Double_t cx, Double_t cy, Double_t cz);
        inline void set_adc(Double_t ax, Double_t ay);
        inline void set_nsr(Int_t nx, Int_t ny);
        inline void set_err(const PartType& type = PartType::Proton);

        inline void set_dummy_x(Double_t cx) { if (!side_(0)) coo_(0) = cx; }

        inline const Short_t&  seqID() const { return seqID_; }
        inline const Short_t&  seqIDcx() const { return seqIDcx_; }
        inline const Short_t&  seqIDcy() const { return seqIDcy_; }
        inline const Short_t&  seqIDex() const { return seqIDex_; }
        inline const Short_t&  seqIDey() const { return seqIDey_; }
        
        inline const Short_t&  lay() const { return lay_; }
        inline const Bool_t&   sx()  const { return side_(0); }
        inline const Bool_t&   sy()  const { return side_(1); }

        inline const SVecD<3>& c()  const { return coo_; }
        inline const Double_t& cx() const { return coo_(0); }
        inline const Double_t& cy() const { return coo_(1); }
        inline const Double_t& cz() const { return coo_(2); }
        
        inline const Double_t& ex() const { return err_(0); }
        inline const Double_t& ey() const { return err_(1); }

        inline Double_t ex(Double_t r) const { return (side_(0) ? (err_(0) * pdf_cx_->efft_sgm( (r/err_(0)) )) : Numc::ZERO<>); }
        inline Double_t ey(Double_t r) const { return (side_(1) ? (err_(1) * pdf_cy_->efft_sgm( (r/err_(1)) )) : Numc::ZERO<>); }

        inline SVecD<2> ionx(Double_t eta, Double_t dzds = 1.0) const { return ((side_(0) && adc_(0) > 0) ? (*pdf_ex_)(adc_(0), eta, dzds) : SVecD<2>()); }
        inline SVecD<2> iony(Double_t eta, Double_t dzds = 1.0) const { return ((side_(1) && adc_(1) > 0) ? (*pdf_ey_)(adc_(1), eta, dzds) : SVecD<2>()); }

    private :
        Short_t    seqID_;
        Short_t    seqIDcx_;
        Short_t    seqIDcy_;
        Short_t    seqIDex_;
        Short_t    seqIDey_;

        Short_t    lay_;
        SVecO<2>   side_;
        SVecI<2>   nsr_;

        SVecD<3>   coo_;  // [cm]
        SVecD<2>   err_;  // [cm]
        SVecD<2>   adc_;

        PartType    type_;
        MultiGauss* pdf_cx_;
        MultiGauss* pdf_cy_;
        IonEloss*   pdf_ex_;
        IonEloss*   pdf_ey_;

    protected :
        static constexpr Double_t NORMFT_X_ =  1.80000e+1; // Normalized Factor X
        static constexpr Double_t NORMFT_Y_ =  8.00000e+0; // Normalized Factor Y
        static constexpr Double_t DEFERR_X_ =  NORMFT_X_ * 1.0e-4; // [cm]
        static constexpr Double_t DEFERR_Y_ =  NORMFT_Y_ * 1.0e-4; // [cm]

        static constexpr MultiGauss::Opt PDF_OPT_ = MultiGauss::Opt::ROBUST;

        static MultiGauss PDF_PR_CX_NN_;
        static MultiGauss PDF_PR_CX_N1_;
        static MultiGauss PDF_PR_CX_N2_;
        static MultiGauss PDF_PR_CX_N3_;
        
        static MultiGauss PDF_PR_CY_NN_;
        static MultiGauss PDF_PR_CY_N1_;
        static MultiGauss PDF_PR_CY_N2_;
        static MultiGauss PDF_PR_CY_N3_;
        static MultiGauss PDF_PR_CY_N4_;

        static IonEloss PDF_PR_EX_;
        static IonEloss PDF_PR_EY_;

    public :
        enum class Orientation {
            kDownward = 0, kUpward = 1
        };
    
        static void Sort(std::vector<HitSt>& hits, const Orientation& ortt = Orientation::kDownward) {
            if (hits.size() < 2) return;
            if (ortt == Orientation::kDownward) std::sort(hits.begin(), hits.end(), [](const HitSt& hit1, const HitSt& hit2) { return (hit1.cz() > hit2.cz()); } );
            else                                std::sort(hits.begin(), hits.end(), [](const HitSt& hit1, const HitSt& hit2) { return (hit1.cz() < hit2.cz()); } );
        }
};
        

MultiGauss HitSt::PDF_PR_CX_NN_(
    PDF_OPT_,
    3.31376712664997630e-01, 1.77875e+01/HitSt::NORMFT_X_,
    5.10255425401639595e-01, 2.65271e+01/HitSt::NORMFT_X_,
    1.58367861933362775e-01, 5.15837e+01/HitSt::NORMFT_X_
);

MultiGauss HitSt::PDF_PR_CX_N1_(
    PDF_OPT_,
    2.75608e+01/HitSt::NORMFT_X_
);

MultiGauss HitSt::PDF_PR_CX_N2_(
    PDF_OPT_,
    6.02616961110044480e-01, 1.82559e+01/HitSt::NORMFT_X_,
    3.97383038889955520e-01, 4.18164e+01/HitSt::NORMFT_X_
);

MultiGauss HitSt::PDF_PR_CX_N3_(
    PDF_OPT_,
    6.75185841999877190e-01, 1.84066e+01/HitSt::NORMFT_X_,
    3.24814158000122755e-01, 4.90605e+01/HitSt::NORMFT_X_
);

MultiGauss HitSt::PDF_PR_CY_NN_(
    PDF_OPT_,
    5.50759994181610257e-01, 7.90750e+00/HitSt::NORMFT_Y_,
    3.74341189078839287e-01, 1.52916e+01/HitSt::NORMFT_Y_,
    7.48988167395504278e-02, 3.63939e+01/HitSt::NORMFT_Y_
);

MultiGauss HitSt::PDF_PR_CY_N1_(
    PDF_OPT_,
    5.45247134760751928e-01, 9.87256e+00/HitSt::NORMFT_Y_,
    4.54752865239248016e-01, 1.65822e+01/HitSt::NORMFT_Y_
);

MultiGauss HitSt::PDF_PR_CY_N2_(
    PDF_OPT_,
    4.53177772355239150e-01, 7.32341e+00/HitSt::NORMFT_Y_,
    4.29177717623073274e-01, 1.19994e+01/HitSt::NORMFT_Y_,
    1.17644510021687618e-01, 2.17922e+01/HitSt::NORMFT_Y_
);

MultiGauss HitSt::PDF_PR_CY_N3_(
    PDF_OPT_,
    5.08190182558828307e-01, 7.74660e+00/HitSt::NORMFT_Y_,
    3.39669567581713017e-01, 1.56378e+01/HitSt::NORMFT_Y_,
    1.52140249859458648e-01, 3.18597e+01/HitSt::NORMFT_Y_
);

MultiGauss HitSt::PDF_PR_CY_N4_(
    PDF_OPT_,
    5.19847432537280163e-01, 7.84945e+00/HitSt::NORMFT_Y_,
    3.39808551148078952e-01, 1.72969e+01/HitSt::NORMFT_Y_,
    1.40344016314640940e-01, 3.85202e+01/HitSt::NORMFT_Y_
);

IonEloss HitSt::PDF_PR_EX_(
    { 2.69731e+00, 8.60498e+00, 1.12322e+00, 2.86260e-02, 1.64740e+00 },
    { 6.96247e-01, 9.55487e+00, 8.56826e-01, 7.95833e-02, 1.50747e+00 }
);

IonEloss HitSt::PDF_PR_EY_(
    { 3.37886e+00, 7.97910e+00, 1.04927e+00, 6.45436e-02, 1.20948e+00 },
    { 1.01843e+00, 7.12581e+00, 7.67523e-01, 1.15677e-01, 2.00335e+00 }
);


} // namesapce TrackSys


#endif // __TRACKLibs_HitSt_H__
