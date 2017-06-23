#ifndef __TRACKLibs_HitSt_H__
#define __TRACKLibs_HitSt_H__

namespace TrackSys {

class HitSt {
    public :
        HitSt(Int_t id = -1, Bool_t sx = true, Bool_t sy = true, Bool_t sz = true) : id_(id), side_(sx, sy, sz) {}
        ~HitSt() {}

        void set_coo(Double_t cx, Double_t cy, Double_t cz) { coo_(0) = cx, coo_(1) = cy, coo_(2) = cz; }
        void set_err(Double_t ex = DEFAULT_ERR(0), Double_t ey = DEFAULT_ERR(1), Double_t ez = DEFAULT_ERR(2)) { err_(0) = ex; err_(1) = ey; err_(2) = ez; }
        void set_coo_and_err(Double_t cx, Double_t cy, Double_t cz, Double_t ex = DEFAULT_ERR(0), Double_t ey = DEFAULT_ERR(1), Double_t ez = DEFAULT_ERR(2)) { set_coo(cx, cy, cz); set_err(ex, ey, ez); }

        inline const Int_t& id() const { return id_; }
        
        inline const SVecO<3>& side() const { return side_; }
        inline const SVecD<3>& coo() const { return coo_; }
        inline const SVecD<3>& err() const { return err_; }
        
        inline const Bool_t&   sx() const { return side_(0); }
        inline const Bool_t&   sy() const { return side_(1); }
        inline const Bool_t&   sz() const { return side_(2); }

        inline const Double_t& cx() const { return coo_(0); }
        inline const Double_t& cy() const { return coo_(1); }
        inline const Double_t& cz() const { return coo_(2); }

        inline const Double_t& ex() const { return err_(0); }
        inline const Double_t& ey() const { return err_(1); }
        inline const Double_t& ez() const { return err_(2); }

    private :
        Int_t    id_;
        SVecO<3> side_;
        SVecD<3> coo_;
        SVecD<3> err_;
    
    public :
        static void SetDefaultErr(Double_t ex = DEFAULT_ERR(0), Double_t ey = DEFAULT_ERR(1), Double_t ez = DEFAULT_ERR(2)) { DEFAULT_ERR(0) = ex; DEFAULT_ERR(1) = ey; DEFAULT_ERR(2) = ez; }

    protected :
        static SVecD<3> DEFAULT_ERR;
        
    public :
        static void Sort(std::vector<HitSt>& hits, Bool_t dirType = true) {
            if (hits.size() < 2) return;
            if (dirType) std::sort(hits.begin(), hits.end(), [](const HitSt& hit1, const HitSt& hit2) { return (hit1.cz() < hit2.cz()); } );
            else         std::sort(hits.begin(), hits.end(), [](const HitSt& hit1, const HitSt& hit2) { return (hit1.cz() > hit2.cz()); } );
        }
};
        
SVecD<3> HitSt::DEFAULT_ERR(24.0e-4, 10.0e-4, 300.0e-4);

}

#endif // __TRACKLibs_HitSt_H__
