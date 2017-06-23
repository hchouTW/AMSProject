#ifndef __TRACKLibs_PhySt_H__
#define __TRACKLibs_PhySt_H__

#include <tuple>

namespace TrackSys {
        
class PhySt {
    public :
        PhySt(PartType type = PartType::Proton) : part_(type), mom_(0), eng_(0), bta_(0) {}
        PhySt(PartInfo part) : part_(part), mom_(0), eng_(0), bta_(0) {}
        ~PhySt() {}

        //void set_space_with_cos();
        //void set_space_with_tan();

        //void set_mom();
        //void set_eng();
        //void set_beta();

        //void set_eta();
        //void set_ieta();
        //void set_rig();
        //void set_irig();

        inline const PartInfo& part() const { return part_; }

        inline const Double_t& mom() const { return mom_; }
        inline const Double_t& eng() const { return eng_; } 
        inline const Double_t& bta() const { return bta_; }

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

        //Double_t ieta_;
        //Double_t irig_;

        SVecD<3> coo_;
        SVecD<3> dir_;
};

}


#endif // __TRACKLibs_PhySt_H__
