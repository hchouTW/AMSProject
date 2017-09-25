#ifndef __TRACKLibs_HitSt_H__
#define __TRACKLibs_HitSt_H__


namespace TrackSys {

class HitSt {
    public :
        HitSt() : side_(true, true, true), coo_(0., 0., 0.), err_(DEFERR_X_, DEFERR_Y_, DEFERR_Z_) {}
        
        HitSt(Double_t cx, Double_t cy, Double_t cz) : HitSt() { set_coo(cx, cy, cz); }

        void print() const;
        
        inline void set_coo(Double_t cx, Double_t cy, Double_t cz) { coo_(0) = cx, coo_(1) = cy, coo_(2) = cz; }
        inline void set_err(Double_t ex, Double_t ey, Double_t ez) { err_(0) = ex, err_(1) = ey, err_(2) = ez; }
        
        inline void set_dummy_x(Double_t cx) { coo_(0) = cx; }

        inline const Bool_t&   sx() const { return side_(0); }
        inline const Bool_t&   sy() const { return side_(1); }
        inline const Bool_t&   sz() const { return side_(2); }
        
        inline const Double_t& cx() const { return coo_(0); }
        inline const Double_t& cy() const { return coo_(1); }
        inline const Double_t& cz() const { return coo_(2); }
        
        inline const Double_t& ex() const { return err_(0); }
        inline const Double_t& ey() const { return err_(1); }
        inline const Double_t& ez() const { return err_(2); }

        inline Double_t ex(Double_t r) const { return (side_(0) ? (err_(0) * pdf_meas_x_.efft_sgm( (r/err_(0)) )) : MGMath::ZERO); }
        inline Double_t ey(Double_t r) const { return (side_(1) ? (err_(1) * pdf_meas_y_.efft_sgm( (r/err_(1)) )) : MGMath::ZERO); }
        inline Double_t ez(Double_t r) const { return (side_(2) ? (err_(2) * pdf_meas_z_.efft_sgm( (r/err_(2)) )) : MGMath::ZERO); }

        inline const SVecO<3>& s() const { return side_; }
        inline const SVecD<3>& c() const { return coo_; }
        inline const SVecD<3>& e() const { return err_; }
        
        inline SVecD<2> e(Double_t rx, Double_t ry)              const { return SVecD<2>(ex(rx), ey(ry)); }
        inline SVecD<3> e(Double_t rx, Double_t ry, Double_t rz) const { return SVecD<3>(ex(rx), ey(ry), ez(rz)); }
        
    private :
        SVecO<3>   side_;
        SVecD<3>   coo_;  // [cm]
        SVecD<3>   err_;  // [cm]
    
    protected :
        static constexpr Double_t DEFERR_X_ =  24.0e-4;
        static constexpr Double_t DEFERR_Y_ =  10.0e-4;
        static constexpr Double_t DEFERR_Z_ = 300.0e-4;
        static MultiGauss pdf_meas_x_;
        static MultiGauss pdf_meas_y_;
        static MultiGauss pdf_meas_z_;
        
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
        
MultiGauss HitSt::pdf_meas_x_(1.0);
MultiGauss HitSt::pdf_meas_y_(1.0);
MultiGauss HitSt::pdf_meas_z_;


} // namesapce TrackSys


#endif // __TRACKLibs_HitSt_H__
