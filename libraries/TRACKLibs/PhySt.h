#ifndef __TRACKLibs_PhySt_H__
#define __TRACKLibs_PhySt_H__


namespace TrackSys {


class PhySt {
    public :
        PhySt(const PartType& type = PartType::Proton) { reset(type); }
        PhySt(const PartInfo& part) { reset(part.type()); }
        ~PhySt() {}

        void reset(const PartType& type = PartType::Proton);

        void set_state_with_cos(Double_t cx, Double_t cy, Double_t cz, Double_t dx = 0., Double_t dy = 0., Double_t dz = -1.);
        void set_state_with_tan(Double_t cx, Double_t cy, Double_t cz, Double_t tx = 0., Double_t ty = 0., Double_t dz = -1.);
        void set_state_with_uxy(Double_t cx, Double_t cy, Double_t cz, Double_t ux = 0., Double_t uy = 0., Short_t signz = -1);

        void set_mom(Double_t mom, Double_t sign = 0.);
        
        void set_eta(Double_t eta);
        
        void set_irig(Double_t irig);

        void set_rig(Double_t rig);

        void print() const;

        inline const PartInfo& part() const { return part_; }

        inline const Double_t& mom()   const { return mom_; }
        inline const Double_t& eng()   const { return eng_; } 
        inline const Double_t& bta()   const { return bta_; }
        inline const Double_t& gmbta() const { return gmbta_; }
        inline const Double_t& eta()   const { return eta_; }
       
        inline Short_t  eta_sign() const { return (MGNumc::Compare(eta_)); }
        inline Double_t eta_abs()  const { return std::fabs(eta_); }

        inline const Double_t& irig()  const { return irig_; }
        inline Double_t rig() const { return (MGNumc::EqualToZero(irig_) ? MGMath::ZERO : MGMath::ONE / irig_); }

        inline const SVecD<3>& coo() const { return coo_; } 
        inline const SVecD<3>& dir() const { return dir_; }

        inline const Double_t& cx() const { return coo_(0); } 
        inline const Double_t& cy() const { return coo_(1); } 
        inline const Double_t& cz() const { return coo_(2); } 
        
        inline const Double_t& dx() const { return dir_(0); } 
        inline const Double_t& dy() const { return dir_(1); } 
        inline const Double_t& dz() const { return dir_(2); } 
        
        inline Double_t tx() const { return ((MGNumc::EqualToZero(dir_(2))) ? 0. : dir_(0)/dir_(2)); } 
        inline Double_t ty() const { return ((MGNumc::EqualToZero(dir_(2))) ? 0. : dir_(1)/dir_(2)); } 
   
    private :
        PartInfo part_;

        Double_t mom_;
        Double_t eng_;
        Double_t bta_;
        
        Double_t gmbta_;
        Double_t eta_;
        
        Double_t irig_;
        
        SVecD<3> coo_;
        SVecD<3> dir_;
};


} // namespace TrackSys


#endif // __TRACKLibs_PhySt_H__
