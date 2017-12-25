#ifndef __TRACKLibs_HitSt_H__
#define __TRACKLibs_HitSt_H__


namespace TrackSys {

class HitSt {
    public :
        HitSt(Bool_t sx = true, Bool_t sy = true) : seqID_(-1), coo_(0., 0., 0.), side_(sx, sy), nsr_(0, 0), err_(DEFERR_X_, DEFERR_Y_), pdf_x_(&PDF_X_NN_), pdf_y_(&PDF_Y_NN_) {}

        void print() const;
        
        inline void set_seqID(Short_t id) { seqID_ = id; }

        inline void set_coo(Double_t cx, Double_t cy, Double_t cz) { coo_(0) = cx; coo_(1) = cy; coo_(2) = cz; }
        
        inline void set_err(Int_t nx = 0, Int_t ny = 0) { 
            nsr_(0) = (nx > 0) ? nx : 0;
            nsr_(1) = (ny > 0) ? ny : 0;
            pdf_x_ = &PDF_X_NN_;
            pdf_y_ = &PDF_Y_NN_;
            if      (nsr_(1) == 1) pdf_y_ = &PDF_Y_N1_;
            else if (nsr_(1) == 2) pdf_y_ = &PDF_Y_N2_;
            else if (nsr_(1) == 3) pdf_y_ = &PDF_Y_N3_;
            else if (nsr_(1) >= 4) pdf_y_ = &PDF_Y_N4_;
        }
        
        inline void set_dummy_x(Double_t cx) { coo_(0) = cx; }

        inline const Short_t&  seqID() const { return seqID_; }
        
        inline const Double_t& cx() const { return coo_(0); }
        inline const Double_t& cy() const { return coo_(1); }
        inline const Double_t& cz() const { return coo_(2); }
        
        inline const Bool_t&   sx() const { return side_(0); }
        inline const Bool_t&   sy() const { return side_(1); }
        
        inline const Double_t& ex() const { return err_(0); }
        inline const Double_t& ey() const { return err_(1); }

        inline Double_t ex(Double_t r) const { return (side_(0) ? (err_(0) * pdf_x_->efft_sgm( (r/err_(0)) )) : MGMath::ZERO); }
        inline Double_t ey(Double_t r) const { return (side_(1) ? (err_(1) * pdf_y_->efft_sgm( (r/err_(1)) )) : MGMath::ZERO); }

        inline const SVecD<3>& c() const { return coo_; }
        inline const SVecO<2>& s() const { return side_; }
        inline const SVecD<2>& e() const { return err_; }

    private :
        Short_t    seqID_;
        SVecD<3>   coo_;  // [cm]
        
        SVecO<2>   side_;
        SVecI<2>   nsr_;
        SVecD<2>   err_;  // [cm]
   
        MultiGauss* pdf_x_;
        MultiGauss* pdf_y_;

    protected :
        static constexpr Double_t NORMFT_X_ =  2.00000e+1; // Normalized Factor X
        static constexpr Double_t NORMFT_Y_ =  7.00000e+0; // Normalized Factor Y
        static constexpr Double_t DEFERR_X_ =  NORMFT_X_ * 1.0e-4; // [cm]
        static constexpr Double_t DEFERR_Y_ =  NORMFT_Y_ * 1.0e-4; // [cm]

        static MultiGauss PDF_X_NN_;
        static MultiGauss PDF_Y_NN_;
        
        static MultiGauss PDF_Y_N1_;
        static MultiGauss PDF_Y_N2_;
        static MultiGauss PDF_Y_N3_;
        static MultiGauss PDF_Y_N4_;
        
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
        

MultiGauss HitSt::PDF_X_NN_(
    MultiGauss::Opt::ROBUST,
    5.48271239075105554e-01, 2.00679e+01/HitSt::NORMFT_X_,
    3.91208480516637092e-01, 3.35472e+01/HitSt::NORMFT_X_,
    6.05202804082573195e-02, 7.84590e+01/HitSt::NORMFT_X_
);

MultiGauss HitSt::PDF_Y_NN_(
    MultiGauss::Opt::ROBUST,
    5.19876182951942711e-01, 7.86516e+00/HitSt::NORMFT_Y_,
    3.87736752283213082e-01, 1.47139e+01/HitSt::NORMFT_Y_,
    9.23870647648442483e-02, 3.30725e+01/HitSt::NORMFT_Y_
);

MultiGauss HitSt::PDF_Y_N1_(
    MultiGauss::Opt::ROBUST,
    5.54810201488121213e-01, 1.00304e+01/HitSt::NORMFT_Y_,
    4.45189798511878787e-01, 1.69959e+01/HitSt::NORMFT_Y_
);

MultiGauss HitSt::PDF_Y_N2_(
    MultiGauss::Opt::ROBUST,
    4.60063022980839176e-01, 7.43073e+00/HitSt::NORMFT_Y_,
    4.24299470955614821e-01, 1.22331e+01/HitSt::NORMFT_Y_,
    1.15637506063545961e-01, 2.21414e+01/HitSt::NORMFT_Y_
);

MultiGauss HitSt::PDF_Y_N3_(
    MultiGauss::Opt::ROBUST,
    4.80200129367643502e-01, 7.71034e+00/HitSt::NORMFT_Y_,
    3.37921187601758155e-01, 1.47931e+01/HitSt::NORMFT_Y_,
    1.81878683030598343e-01, 3.01339e+01/HitSt::NORMFT_Y_
);

MultiGauss HitSt::PDF_Y_N4_(
    MultiGauss::Opt::ROBUST,
    4.88650682054746988e-01, 7.79726e+00/HitSt::NORMFT_Y_,
    3.32706018621177213e-01, 1.60200e+01/HitSt::NORMFT_Y_,
    1.78643299324075799e-01, 3.46379e+01/HitSt::NORMFT_Y_
);


} // namesapce TrackSys


#endif // __TRACKLibs_HitSt_H__
