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
    5.23888719314210882e-01, 1.94095e+01/HitSt::NORMFT_X_,
    3.98349799428234230e-01, 3.23686e+01/HitSt::NORMFT_X_,
    6.64332057378534818e-02, 6.75576e+01/HitSt::NORMFT_X_,
    1.13282755197015063e-02, 1.07105e+02/HitSt::NORMFT_X_
);

MultiGauss HitSt::PDF_Y_NN_(
    MultiGauss::Opt::ROBUST,
    3.20904952818954037e-01, 7.10154e+00/HitSt::NORMFT_Y_,
    4.08240588544126815e-01, 1.08899e+01/HitSt::NORMFT_Y_,
    2.19072329640065483e-01, 1.90135e+01/HitSt::NORMFT_Y_,
    5.17821289968536519e-02, 4.39795e+01/HitSt::NORMFT_Y_
);

MultiGauss HitSt::PDF_Y_N1_(
    MultiGauss::Opt::ROBUST,
    5.05132735024615420e-01, 9.78539e+00/HitSt::NORMFT_Y_,
    4.94867264975384524e-01, 1.62914e+01/HitSt::NORMFT_Y_
);

MultiGauss HitSt::PDF_Y_N2_(
    MultiGauss::Opt::ROBUST,
    4.50592395256901423e-01, 7.38835e+00/HitSt::NORMFT_Y_,
    4.31419140679972724e-01, 1.21601e+01/HitSt::NORMFT_Y_,
    1.17988464063125839e-01, 2.22308e+01/HitSt::NORMFT_Y_
);

MultiGauss HitSt::PDF_Y_N3_(
    MultiGauss::Opt::ROBUST,
    3.62401409806186903e-01, 7.23893e+00/HitSt::NORMFT_Y_,
    3.31112856289004998e-01, 1.15761e+01/HitSt::NORMFT_Y_,
    2.41202035075139293e-01, 2.21391e+01/HitSt::NORMFT_Y_,
    6.52836988296688209e-02, 4.21914e+01/HitSt::NORMFT_Y_
);

MultiGauss HitSt::PDF_Y_N4_(
    MultiGauss::Opt::ROBUST,
    3.35396173686006438e-01, 7.14830e+00/HitSt::NORMFT_Y_,
    3.15241468764665189e-01, 1.12436e+01/HitSt::NORMFT_Y_,
    2.64802023251396956e-01, 2.24157e+01/HitSt::NORMFT_Y_,
    8.45603342979313888e-02, 5.22629e+01/HitSt::NORMFT_Y_
);


} // namesapce TrackSys


#endif // __TRACKLibs_HitSt_H__
