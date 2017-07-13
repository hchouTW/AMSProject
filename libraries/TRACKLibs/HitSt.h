#ifndef __TRACKLibs_HitSt_H__
#define __TRACKLibs_HitSt_H__


namespace TrackSys {


class MultiGauss {
    public :
        MultiGauss() {}
        MultiGauss(Double_t sgm);
        MultiGauss(Double_t wgt1, Double_t sgm1, Double_t wgt2, Double_t sgm2);
        MultiGauss(Double_t wgt1, Double_t sgm1, Double_t wgt2, Double_t sgm2, Double_t wgt3, Double_t sgm3);
        ~MultiGauss() {}

        inline Int_t num() const { return multi_gauss_.size(); }
        inline const Double_t& wgt(Int_t i) const { return multi_gauss_.at(i).first; }
        inline const Double_t& sgm(Int_t i) const { return multi_gauss_.at(i).second; }

        inline Double_t efft_sgm(Double_t r = 0.) const; 

    private :
        std::vector<std::pair<Double_t, Double_t>> multi_gauss_;

    private :
        static constexpr Double_t LMTL_PROB = 1.0e-6;
};


class HitSt {
    public :
        enum class Orientation {
            kDownward = 0, kUpward = 1
        };
    
    public :
        HitSt(Int_t id = -1, Bool_t sx = true, Bool_t sy = true, Bool_t sz = true) : id_(id), side_(sx, sy, sz), coo_(0., 0., 0.), err_x_(DEFAULT_ERR_X), err_y_(DEFAULT_ERR_Y), err_z_(DEFAULT_ERR_Z) {}
        ~HitSt() {}

        inline void set_id(Int_t id) { id_ = id; }
        inline void set_coo(Double_t cx, Double_t cy, Double_t cz) { coo_(0) = cx, coo_(1) = cy, coo_(2) = cz; }
        inline void set_err(const MultiGauss& ex = DEFAULT_ERR_X, const MultiGauss& ey = DEFAULT_ERR_Y, const MultiGauss& ez = DEFAULT_ERR_Z) { err_x_ = ex; err_y_ = ey; err_z_ = ez; }
        inline void set_coo_and_err(Double_t cx, Double_t cy, Double_t cz, const MultiGauss& ex = DEFAULT_ERR_X, const MultiGauss& ey = DEFAULT_ERR_Y, const MultiGauss& ez = DEFAULT_ERR_Z) { set_coo(cx, cy, cz); set_err(ex, ey, ez); }

        inline void set_dummy_x(Double_t cx) { coo_(0) = cx; }

        void print() const;

        inline const Int_t& id() const { return id_; }
        
        inline const SVecO<3>& side() const { return side_; }
        inline const SVecD<3>& coo() const { return coo_; }

        inline SVecD<3> err(const SVecD<3>& res = SVecD<3>()) const { return SVecD<3>(err_x_.efft_sgm(res(0)), err_y_.efft_sgm(res(1)), err_z_.efft_sgm(res(2))); }
        inline SVecD<2> err(const SVecD<2>& res) const { return SVecD<2>(err_x_.efft_sgm(res(0)), err_y_.efft_sgm(res(1))); }
        
        inline const Bool_t&   sx() const { return side_(0); }
        inline const Bool_t&   sy() const { return side_(1); }
        inline const Bool_t&   sz() const { return side_(2); }

        inline const Double_t& cx() const { return coo_(0); }
        inline const Double_t& cy() const { return coo_(1); }
        inline const Double_t& cz() const { return coo_(2); }
        
        inline const MultiGauss& ex() const { return err_x_; }
        inline const MultiGauss& ey() const { return err_y_; }
        inline const MultiGauss& ez() const { return err_z_; }

        inline Double_t ex(Double_t r) const { return err_x_.efft_sgm(r); }
        inline Double_t ey(Double_t r) const { return err_y_.efft_sgm(r); }
        inline Double_t ez(Double_t r) const { return err_z_.efft_sgm(r); }

    private :
        Int_t      id_;
        SVecO<3>   side_;
        SVecD<3>   coo_;     // [cm]
        MultiGauss err_x_;   // [cm]
        MultiGauss err_y_;   // [cm]
        MultiGauss err_z_;   // [cm]
    
    public :
        static void SetDefaultErrX(const MultiGauss& ex = DEFAULT_ERR_X) { DEFAULT_ERR_X = ex; }
        static void SetDefaultErrY(const MultiGauss& ey = DEFAULT_ERR_Y) { DEFAULT_ERR_Y = ey; }
        static void SetDefaultErrZ(const MultiGauss& ez = DEFAULT_ERR_Z) { DEFAULT_ERR_Z = ez; }

        static const MultiGauss& GetDefaultErrX() { return DEFAULT_ERR_X; }
        static const MultiGauss& GetDefaultErrY() { return DEFAULT_ERR_Y; }
        static const MultiGauss& GetDefaultErrZ() { return DEFAULT_ERR_Z; }
        
        static Double_t GetDefaultErrX(Double_t r) { return DEFAULT_ERR_X.efft_sgm(r); }
        static Double_t GetDefaultErrY(Double_t r) { return DEFAULT_ERR_Y.efft_sgm(r); }
        static Double_t GetDefaultErrZ(Double_t r) { return DEFAULT_ERR_Z.efft_sgm(r); }

    protected :
        static MultiGauss DEFAULT_ERR_X;
        static MultiGauss DEFAULT_ERR_Y;
        static MultiGauss DEFAULT_ERR_Z;
        
    public :
        static void Sort(std::vector<HitSt>& hits, const Orientation& ortt = Orientation::kDownward) {
            if (hits.size() < 2) return;
            if (ortt == Orientation::kDownward) std::sort(hits.begin(), hits.end(), [](const HitSt& hit1, const HitSt& hit2) { return (hit1.cz() > hit2.cz()); } );
            else                                std::sort(hits.begin(), hits.end(), [](const HitSt& hit1, const HitSt& hit2) { return (hit1.cz() < hit2.cz()); } );
        }
};
        
MultiGauss HitSt::DEFAULT_ERR_X(24.0e-4);
MultiGauss HitSt::DEFAULT_ERR_Y(10.0e-4);
MultiGauss HitSt::DEFAULT_ERR_Z(300.0e-4);


} // namesapce TrackSys


#endif // __TRACKLibs_HitSt_H__
