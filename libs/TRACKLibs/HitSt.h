#ifndef __TRACKLibs_HitSt_H__
#define __TRACKLibs_HitSt_H__


namespace TrackSys {

class HitSt {
    public :
        HitSt(Bool_t sx = true, Bool_t sy = true, Bool_t sz = true) : seqID_(-1), side_(sx, sy, sz), coo_(0., 0., 0.), err_(DEFERR_X_, DEFERR_Y_, DEFERR_Z_) {}
        
        HitSt(Double_t cx, Double_t cy, Double_t cz) : HitSt() { set_coo(cx, cy, cz); }

        void print() const;
        
        inline void set_seqID(Short_t id) { seqID_ = id; }

        inline void set_coo(Double_t cx, Double_t cy, Double_t cz) { coo_(0) = cx, coo_(1) = cy, coo_(2) = cz; }
        inline void set_err(Double_t ex, Double_t ey, Double_t ez) { err_(0) = ex, err_(1) = ey, err_(2) = ez; }
        
        inline void set_dummy_x(Double_t cx) { coo_(0) = cx; }

        inline const Short_t&  seqID() const { return seqID_; }
        
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

        inline Bool_t is_norm_x(Double_t r) { return pdf_meas_x_.is_norm( (r/err_(0)) ); }
        inline Bool_t is_norm_y(Double_t r) { return pdf_meas_y_.is_norm( (r/err_(1)) ); }
        inline Bool_t is_norm_z(Double_t r) { return pdf_meas_z_.is_norm( (r/err_(2)) ); }
        
    private :
        Short_t    seqID_;

        SVecO<3>   side_;
        SVecD<3>   coo_;  // [cm]
        SVecD<3>   err_;  // [cm]
    
    protected :
        //static constexpr Double_t DEFERR_X_ =  1.04562e-3;
        //static constexpr Double_t DEFERR_Y_ =  7.65345e-4;
        static constexpr Double_t DEFERR_Z_ =  3.00000e-2;
        static constexpr Double_t DEFERR_X_ =  1.12212e-3;
        static constexpr Double_t DEFERR_Y_ =  5.33089e-4;
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
        
//MultiGauss HitSt::pdf_meas_x_(3.11933447340336616e-02, 1.0, 6.05897933730214833e-01, 1.87175073162334304e+00, 3.11933447340336616e-02, 5.60371836804957812e+00);
//MultiGauss HitSt::pdf_meas_y_(3.17504515061161174e-01, 1.0, 5.04722176029590797e-01, 2.11550346575727266e+00, 1.77773308909248029e-01, 4.22574133234031724e+00);
MultiGauss HitSt::pdf_meas_z_;

MultiGauss HitSt::pdf_meas_x_(7.72558811989003946e-02, 1.0, 5.65389654829016375e-01, 1.81143727943535460e+00, 3.57354463972083258e-01, 5.44595052222578691e+00);
MultiGauss HitSt::pdf_meas_y_(1.44671932077881177e-01, 1.0, 4.74432270171227743e-01, 1.80172166373719955e+00, 3.37116091570108967e-01, 3.25272140299274604e+00, 4.37797061807821616e-02, 8.08217764763482194e+00);


} // namesapce TrackSys


#endif // __TRACKLibs_HitSt_H__
