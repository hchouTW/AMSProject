#ifndef __TRACKLibs_PhySt_H__
#define __TRACKLibs_PhySt_H__

#include <tuple>

namespace TrackSys {
        
class PhySt {
    public :
        PhySt(PartType type = PartType::Proton) : part_(type), mom_(0), eng_(0), bta_(0), ieta_(0), irig_(0), coo_(0, 0, 0), dir_(0, 0, -1) { eng_ = part_.mass(); if (part_.is_massless()) bta_ = 1.;}
        PhySt(PartInfo part) : part_(part), mom_(0), eng_(0), bta_(0), ieta_(0), irig_(0), coo_(0, 0, 0), dir_(0, 0, -1) { eng_ = part_.mass(); if (part_.is_massless()) bta_ = 1.; }
        ~PhySt() {}

        void set_state_with_cos(Double_t cx, Double_t cy, Double_t cz, Double_t dx, Double_t dy, Double_t dz);
        void set_state_with_tan(Double_t cx, Double_t cy, Double_t cz, Double_t tx, Double_t ty, Double_t dz);

        void set_mom(Double_t mom, Double_t sign = 0.);
        
        void set_ieta(Double_t ieta);
        
        void set_eta(Double_t eta);
        
        void set_irig(Double_t irig);

        void set_rig(Double_t rig);

        inline const PartInfo& part() const { return part_; }

        inline const Double_t& mom() const { return mom_; }
        inline const Double_t& eng() const { return eng_; } 
        inline const Double_t& bta() const { return bta_; }
        inline const Double_t& gmbta() const { return gmbta_; }
        inline const Double_t& ieta() const { return ieta_; }
        inline const Double_t& irig() const { return irig_; }

        inline Double_t eta() const { return (MGNumc::EqualToZero(ieta_) ? MGMath::ZERO : MGMath::ONE / ieta_); }
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
        Double_t ieta_;
        Double_t irig_;
        
        SVecD<3> coo_;
        SVecD<3> dir_;
};

}


#endif // __TRACKLibs_PhySt_H__
