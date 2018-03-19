#ifndef __TRACKLibs_HitSt_H__
#define __TRACKLibs_HitSt_H__


namespace TrackSys {
/*
class VirtualHitSt {
    public :
        enum class Detector {
            TRK, TOF, TRD, RICH, ECAL
        };

    public :
        VirtualHitSt(Detector dec = Detector::TRK, Short_t lay = 0, Bool_t csx = false, Bool_t csy = false, Bool_t csz = true) : type_(PartType::Proton), dec_(dec), lay_(lay), coo_side_(csx, csy, csz), coo_(Numc::ZERO<>, Numc::ZERO<>, Numc::ZERO<>), cer_(Numc::ONE<>, Numc::ONE<>, Numc::ONE<>) {}
        ~VirtualHitSt() {}
       
        inline virtual void cal(const PhySt& part) = 0;
        inline virtual void set_type(PartType type) = 0;
        
        inline void set_coo(Double_t cx, Double_t cy, Double_t cz) { coo_ = std::move(SVecD<3>(cx, cy, cz)); }
        inline void set_cer(Double_t ex, Double_t ey, Double_t ez) { cer_ = std::move(SVecD<3>(ex, ey, ez)); }
        inline void set_dummy_x(Double_t cx) { if (!coo_side_(0)) coo_(0) = cx; }
        inline void set_dummy_y(Double_t cy) { if (!coo_side_(1)) coo_(1) = cy; }
        inline void set_dummy_z(Double_t cz) { if (!coo_side_(2)) coo_(2) = cz; }
        
        inline const PartType& type() const { return type_; }
        inline const Detector& dec()  const { return dec_; }
        inline const Short_t&  lay()  const { return lay_; }

        inline const SVecO<3>& cs()  const { return coo_side_; }
        inline const Bool_t&   csx() const { return coo_side_(0); }
        inline const Bool_t&   csy() const { return coo_side_(1); }
        inline const Bool_t&   csz() const { return coo_side_(2); }
        
        inline const SVecD<3>& c()  const { return coo_; }
        inline const Double_t& cx() const { return coo_(0); }
        inline const Double_t& cy() const { return coo_(1); }
        inline const Double_t& cz() const { return coo_(2); }
        
        inline const SVecD<3>& ce()  const { return cer_; }
        inline const Double_t& cex() const { return cer_(0); }
        inline const Double_t& cey() const { return cer_(1); }
        inline const Double_t& cez() const { return cer_(2); }

    protected :
        PartType type_; // particle type
        
        Detector dec_; // TRK TOF TRD RICH ECAL
        Short_t  lay_; // start from 0

        SVecO<3> side_; // (x, y, z)
        SVecD<3> coo_;  // [cm] coord
        SVecD<3> cer_;  // [cm] error
};


class HitStTRK : public VirtualHitSt {
    public :
        HitStTRK(Short_t lay = 0, Bool_t csx = false, Bool_t csy = false) : VirtualHitSt(VirtualHitSt::Detector::TRK, lay, csx, csy) { clear(); }
        ~HitStTRK() {}
        
        inline void cal(const PhySt& part);
        inline void set_type(PartType type);
        
        inline void set_nsr(Short_t nx, Short_t ny) {
            nsr_(0) = (Numc::Compare(nx)>0) ? nx : Numc::ZERO<Short_t>;
            nsr_(1) = (Numc::Compare(ny)>0) ? ny : Numc::ZERO<Short_t>;
        }
        inline void set_adc(Double_t ax, Double_t ay) {
            adc_side_(0) = (Numc::Compare(ax) > 0);
            adc_side_(1) = (Numc::Compare(ay) > 0);
            adc_(0) = (adc_side_(0) ? ax : Numc::ZERO<>);
            adc_(1) = (adc_side_(1) ? ay : Numc::ZERO<>);
        }
        
        inline const SVecD<2>& cnrm()  const { return cnrm_; }
        inline const Double_t& cnrmx() const { return cnrm_(0); }
        inline const Double_t& cnrmy() const { return cnrm_(1); }

        inline const SVecD<2>& cdiv()  const { return cdiv_; }
        inline const Double_t& cdivx() const { return cdiv_(0); }
        inline const Double_t& cdivy() const { return cdiv_(1); }

        inline const SVecD<2>& anrm()  const { return anrm_; }
        inline const Double_t& anrmx() const { return anrm_(0); }
        inline const Double_t& anrmy() const { return anrm_(1); }

        inline const SVecD<2>& adiv()  const { return adiv_; }
        inline const Double_t& adivx() const { return adiv_(0); }
        inline const Double_t& adivy() const { return adiv_(1); }

    protected :
        void clear();

    protected :
        SVecS<2> nsr_; // Number of strip
        
        SVecO<2> adc_side_;
        SVecD<2> adc_; // ADC

        SVecD<2> cnrm_; // coord norm
        SVecD<2> cdiv_; // coord div
        SVecD<2> anrm_; // adc nrom
        SVecD<2> adiv_; // adc div

        MultiGaus* pdf_cx_;
        MultiGaus* pdf_cy_;
        IonEloss*  pdf_ex_;
        IonEloss*  pdf_ey_;
    
    protected :
        static constexpr Double_t DEFERR_X_ = 1.80000e+1 * 1.0e-4; // [cm]
        static constexpr Double_t DEFERR_Y_ = 8.00000e+0 * 1.0e-4; // [cm]
        static constexpr Double_t DEFERR_Z_ = 1.50000e+2 * 1.0e-4; // [cm]
        
        static MultiGaus PDF_PR_CX_NN_;
        static MultiGaus PDF_PR_CX_N1_;
        static MultiGaus PDF_PR_CX_N2_;
        static MultiGaus PDF_PR_CX_N3_;
        
        static MultiGaus PDF_PR_CY_NN_;
        static MultiGaus PDF_PR_CY_N1_;
        static MultiGaus PDF_PR_CY_N2_;
        static MultiGaus PDF_PR_CY_N3_;
        static MultiGaus PDF_PR_CY_N4_;

        static IonEloss PDF_PR_EX_;
        static IonEloss PDF_PR_EY_;
};

*/




































////////////////////////////////////////////////////////////////////////////////////////
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

        inline SVecD<2> ionx(Double_t eta) const { return ((side_(0) && adc_(0) > 0) ? (*pdf_ex_)(adc_(0), eta) : SVecD<2>()); }
        inline SVecD<2> iony(Double_t eta) const { return ((side_(1) && adc_(1) > 0) ? (*pdf_ey_)(adc_(1), eta) : SVecD<2>()); }

    private :
        Short_t    seqID_;
        Short_t    seqIDcx_;
        Short_t    seqIDcy_;
        Short_t    seqIDex_;
        Short_t    seqIDey_;

        //Short_t    dec_;
        Short_t    lay_;
        SVecO<2>   side_;
        SVecI<2>   nsr_;

        SVecD<3>   coo_;  // [cm]
        SVecD<2>   err_;  // [cm]
        SVecD<2>   adc_;

        PartType    type_;
        MultiGaus* pdf_cx_;
        MultiGaus* pdf_cy_;
        IonEloss*   pdf_ex_;
        IonEloss*   pdf_ey_;

    protected :
        static constexpr Double_t NORMFT_X_ =  1.80000e+1; // Normalized Factor X
        static constexpr Double_t NORMFT_Y_ =  8.00000e+0; // Normalized Factor Y
        static constexpr Double_t DEFERR_X_ =  NORMFT_X_ * 1.0e-4; // [cm]
        static constexpr Double_t DEFERR_Y_ =  NORMFT_Y_ * 1.0e-4; // [cm]

        static constexpr MultiGaus::Opt PDF_OPT_ = MultiGaus::Opt::ROBUST;

        static MultiGaus PDF_PR_CX_NN_;
        static MultiGaus PDF_PR_CX_N1_;
        static MultiGaus PDF_PR_CX_N2_;
        static MultiGaus PDF_PR_CX_N3_;
        
        static MultiGaus PDF_PR_CY_NN_;
        static MultiGaus PDF_PR_CY_N1_;
        static MultiGaus PDF_PR_CY_N2_;
        static MultiGaus PDF_PR_CY_N3_;
        static MultiGaus PDF_PR_CY_N4_;

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
        

MultiGaus HitSt::PDF_PR_CX_NN_(
    PDF_OPT_,
    3.31376712664997630e-01, 1.77875e+01/HitSt::NORMFT_X_,
    5.10255425401639595e-01, 2.65271e+01/HitSt::NORMFT_X_,
    1.58367861933362775e-01, 5.15837e+01/HitSt::NORMFT_X_
);

MultiGaus HitSt::PDF_PR_CX_N1_(
    PDF_OPT_,
    2.75608e+01/HitSt::NORMFT_X_
);

MultiGaus HitSt::PDF_PR_CX_N2_(
    PDF_OPT_,
    6.02616961110044480e-01, 1.82559e+01/HitSt::NORMFT_X_,
    3.97383038889955520e-01, 4.18164e+01/HitSt::NORMFT_X_
);

MultiGaus HitSt::PDF_PR_CX_N3_(
    PDF_OPT_,
    6.75185841999877190e-01, 1.84066e+01/HitSt::NORMFT_X_,
    3.24814158000122755e-01, 4.90605e+01/HitSt::NORMFT_X_
);

MultiGaus HitSt::PDF_PR_CY_NN_(
    PDF_OPT_,
    5.50759994181610257e-01, 7.90750e+00/HitSt::NORMFT_Y_,
    3.74341189078839287e-01, 1.52916e+01/HitSt::NORMFT_Y_,
    7.48988167395504278e-02, 3.63939e+01/HitSt::NORMFT_Y_
);

MultiGaus HitSt::PDF_PR_CY_N1_(
    PDF_OPT_,
    5.45247134760751928e-01, 9.87256e+00/HitSt::NORMFT_Y_,
    4.54752865239248016e-01, 1.65822e+01/HitSt::NORMFT_Y_
);

MultiGaus HitSt::PDF_PR_CY_N2_(
    PDF_OPT_,
    4.53177772355239150e-01, 7.32341e+00/HitSt::NORMFT_Y_,
    4.29177717623073274e-01, 1.19994e+01/HitSt::NORMFT_Y_,
    1.17644510021687618e-01, 2.17922e+01/HitSt::NORMFT_Y_
);

MultiGaus HitSt::PDF_PR_CY_N3_(
    PDF_OPT_,
    5.08190182558828307e-01, 7.74660e+00/HitSt::NORMFT_Y_,
    3.39669567581713017e-01, 1.56378e+01/HitSt::NORMFT_Y_,
    1.52140249859458648e-01, 3.18597e+01/HitSt::NORMFT_Y_
);

MultiGaus HitSt::PDF_PR_CY_N4_(
    PDF_OPT_,
    5.19847432537280163e-01, 7.84945e+00/HitSt::NORMFT_Y_,
    3.39808551148078952e-01, 1.72969e+01/HitSt::NORMFT_Y_,
    1.40344016314640940e-01, 3.85202e+01/HitSt::NORMFT_Y_
);

IonEloss HitSt::PDF_PR_EX_(
    { 6.84708e-04, 1.22107e+00, 1.54109e+00, 2.83791e+00, 1.06167e+00, 7.87652e+00 }, // Kpa
    { 1.01510e+00, 2.06220e+01, 1.24078e+00, 4.82421e-04, 5.80771e+00 }, // Mpv
    { 6.21439e-02, 3.05480e+01, 1.37339e+00, 1.07762e-04, 8.70839e+00 }, // Sgm
    5.00000e+00 // Fluc
);

IonEloss HitSt::PDF_PR_EY_( // need to charge
    { 1.04571e-03, 1.56996e+00, 3.42208e+00, 1.60653e+00, 1.36027e+00, 1.00964e+01 }, // Kpa
    { 4.41904e+00, 6.24969e+00, 1.05197e+00, 1.39822e-01, 1.03750e+00 }, // Mpv
    { 7.13109e+00, 2.98634e+00, 5.79009e-01, 4.32315e+00, 7.58125e-01 }, // Sgm
    4.40000e+00 // Fluc
);


} // namesapce TrackSys


#endif // __TRACKLibs_HitSt_H__
