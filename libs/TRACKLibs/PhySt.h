#ifndef __TRACKLibs_PhySt_H__
#define __TRACKLibs_PhySt_H__


namespace TrackSys {


class PhyArg {
    public :
        static void SetOpt(Bool_t opt_mscat = true, Bool_t opt_eloss = true) { opt_mscat_ = opt_mscat; opt_eloss_ = opt_eloss; }
        static void SetOptMscat(Bool_t opt_mscat = true) { opt_mscat_ = opt_mscat; }
        static void SetOptEloss(Bool_t opt_eloss = true) { opt_eloss_ = opt_eloss; }

        static const Bool_t& OptMscat() { return opt_mscat_; }
        static const Bool_t& OptEloss() { return opt_eloss_; }
   
    private :
        static Bool_t opt_mscat_;
        static Bool_t opt_eloss_;

    public :
        PhyArg(Bool_t sw_mscat = OptMscat(), Bool_t sw_eloss = OptEloss()) { reset(sw_mscat, sw_eloss); }
        ~PhyArg() {}

        inline void reset(Bool_t sw_mscat = false, Bool_t sw_eloss = false) { mat_ = (sw_mscat || sw_eloss); sw_mscat_ = sw_mscat; sw_eloss_ = sw_eloss; mscatu_tau_ = 0.0; mscatu_rho_ = 0.0; mscatc_tau_ = 0.0; mscatc_rho_ = 0.0; eloss_ion_ = 0.0; eloss_brm_ = 0.0; } 
        
        inline void set_mscat(Double_t tauu = 0.0, Double_t rhou = 0.0, Double_t tauc = 0.0, Double_t rhoc = 0.0) { if (sw_mscat_) {mscatu_tau_ = tauu; mscatu_rho_ = rhou; mscatc_tau_ = tauc; mscatc_rho_ = rhoc; } }
        
        inline void set_eloss(Double_t ion = 0.0, Double_t brm = 0.0) { if (sw_eloss_) { eloss_ion_ = ((ion<-4.0)?0.0:ion); eloss_brm_ = ((brm<0.0)?0.0:brm); } }
        
        inline void set(Double_t tauu = 0.0, Double_t rhou = 0.0, Double_t tauc = 0.0, Double_t rhoc = 0.0, Double_t ion = 0.0, Double_t brm = 0.0) { set_mscat(tauu, rhou, tauc, rhoc); set_eloss(ion, brm); }

        void rndm_mscatu() { mscatu_tau_ = 0.0; mscatu_rho_ = 0.0; if (sw_mscat_) { mscatu_tau_ = pdf_mscat_.rndm(); mscatu_rho_ = pdf_mscat_.rndm(); } }
        void rndm_mscatc() { mscatc_tau_ = 0.0; mscatc_rho_ = 0.0; if (sw_mscat_) { mscatc_tau_ = pdf_mscat_.rndm(); mscatc_rho_ = pdf_mscat_.rndm(); } }
        void rndm_mscat() { rndm_mscatu(); rndm_mscatc(); }

        void rndm_eloss_ion(Double_t kpa = 0.0, Double_t mos = 15.0);
        void rndm_eloss_brm(Double_t nrl = 0.0) { if (sw_eloss_) { Double_t bremslen = nrl / MGMath::LOG_TWO; eloss_brm_=((bremslen<=0.0)?0.0:MGRndm::Gamma(bremslen,1.0/bremslen)()); } }
        void rndm_eloss(Double_t kpa = 0.0, Double_t mos = 15.0, Double_t nrl = 0.0) { rndm_eloss_ion(kpa, mos); rndm_eloss_brm(nrl); }

        void rndm(Double_t kpa = 0.0, Double_t mos = 15.0, Double_t nrl = 0.0) { rndm_mscatu(); rndm_mscatc(); rndm_eloss_ion(kpa, mos); rndm_eloss_brm(nrl); }

        inline const Bool_t& operator() () const { return mat_; }

        inline const Bool_t& mscat() const { return sw_mscat_; }
        inline const Bool_t& eloss() const { return sw_eloss_; }

        inline const Double_t& tauu() const { return mscatu_tau_; }
        inline const Double_t& rhou() const { return mscatu_rho_; }
        inline const Double_t& tauc() const { return mscatc_tau_; }
        inline const Double_t& rhoc() const { return mscatc_rho_; }
        
        inline const Double_t& ion() const { return eloss_ion_; }
        inline const Double_t& brm() const { return eloss_brm_; }
        
        //inline static MultiGauss& pdf_mscat() { return pdf_mscat_; }

    private :
        Bool_t   mat_;
        Bool_t   sw_mscat_;
        Bool_t   sw_eloss_;
        Double_t mscatu_tau_;
        Double_t mscatu_rho_;
        Double_t mscatc_tau_;
        Double_t mscatc_rho_;
        Double_t eloss_ion_;
        Double_t eloss_brm_;
        
    protected :
        static MultiGauss pdf_mscat_;

    protected :
        static constexpr Int_t    NPX_      = 400;
        static constexpr Double_t LMT_SGM_  = 4.0;
        static constexpr Double_t STEP_KPA_ = 0.02;
        static constexpr Double_t STEP_MOS_ = 0.2;
        static std::map<std::pair<Int_t, Int_t>, TF1 *> pdf_eloss_ion_;
};

Bool_t PhyArg::opt_mscat_ = true;
Bool_t PhyArg::opt_eloss_ = true;

MultiGauss PhyArg::pdf_mscat_(8.931154e-01, 1.000000e+00, 9.211506e-02, 1.851060e+00, 1.167630e-02, 4.478370e+00, 3.093243e-03, 2.390957e+01);

std::map<std::pair<Int_t, Int_t>, TF1 *> PhyArg::pdf_eloss_ion_;


class VirtualPhySt {
    public :
        VirtualPhySt() { reset(); }
        ~VirtualPhySt() {}

        inline void reset() { len_ = 0.0; nrl_ = 0.0; tau_ = std::move(SVecD<3>(1.0, 0.0, 0.0)); rho_ = std::move(SVecD<3>(0.0, -1.0, 0.0)); mscatu_ = 0.0; mscatcu_ = 0.0; mscatcl_ = 0.0; eloss_ion_kpa_ = 0.0; eloss_ion_mos_ = 15.0; eloss_ion_ = 0.0; eloss_brm_ = 0.0; }

        inline void set_len(Double_t len) { len_ = len; }
        inline void set_nrl(Double_t nrl) { nrl_ = nrl; }
        
        inline const Double_t& len() const { return len_; }
        inline const Double_t& nrl() const { return nrl_; }

        inline void set_orth(const SVecD<3>& tau, const SVecD<3>& rho) { tau_ = tau; rho_ = rho; }
        inline const SVecD<3>& tau() const { return tau_; }
        inline const SVecD<3>& rho() const { return rho_; }

        inline void set_mscatu(Double_t mscatu) { mscatu_ = mscatu; }
        inline void set_mscatc(Double_t mscatcu, Double_t mscatcl) { mscatcu_ = mscatcu; mscatcl_ = mscatcl; } 

        inline const Double_t& mscatu()  const { return mscatu_; }
        inline const Double_t& mscatcu() const { return mscatcu_; }
        inline const Double_t& mscatcl() const { return mscatcl_; }

        inline SVecD<3> symbk_mscatu(Double_t tauu = 0.0, Double_t rhou = 0.0) { return ((tauu*mscatu_)*tau_ + (rhou*mscatu_)*rho_); }
        inline SVecD<3> symbk_mscatc(Double_t tauu = 0.0, Double_t rhou = 0.0, Double_t tauc = 0.0, Double_t rhoc = 0.0) { return ((tauu*mscatcu_+tauc*mscatcl_)*tau_ + (rhou*mscatcu_+rhoc*mscatcl_)*rho_); }

        inline void set_eloss_ion(Double_t eloss_ion, Double_t kpa, Double_t mos) { eloss_ion_ = eloss_ion; eloss_ion_kpa_ = kpa; eloss_ion_mos_ = mos; }
        inline void set_eloss_brm(Double_t eloss_brm) { eloss_brm_ = eloss_brm; }
        
        inline const Double_t& eloss_ion_kpa() const { return eloss_ion_kpa_; }
        inline const Double_t& eloss_ion_mos() const { return eloss_ion_mos_; }

        inline const Double_t& eloss_ion() const { return eloss_ion_; }
        inline const Double_t& eloss_brm() const { return eloss_brm_; }
        
        inline Double_t symbk_eloss(Double_t ion = 0.0, Double_t brm = 0.0) { return (ion * eloss_ion_ + brm * eloss_brm_); }

    protected :
        //inline Double_t proj_xy(const SVecD<3>& vec) { return (MGMath::ONE - vec(2) * vec(2)); }

    private :
        Double_t len_;
        Double_t nrl_;

        SVecD<3> tau_;
        SVecD<3> rho_;

        Double_t mscatu_;
        Double_t mscatcu_;
        Double_t mscatcl_;

        Double_t eloss_ion_kpa_;
        Double_t eloss_ion_mos_;

        Double_t eloss_ion_;
        Double_t eloss_brm_;
};


class PhySt {
    public :
        PhySt(const PartType& type = PartType::Proton, Bool_t sw_mscat = PhyArg::OptMscat(), Bool_t sw_eloss = PhyArg::OptEloss()) : arg_(sw_mscat, sw_eloss) { reset(type); }
        PhySt(const PartInfo& part, Bool_t sw_mscat = PhyArg::OptMscat(), Bool_t sw_eloss = PhyArg::OptEloss()) : arg_(sw_mscat, sw_eloss) { reset(part.type()); }
        ~PhySt() {}

        void reset(const PartType& type = PartType::Proton, Bool_t sw_mscat = PhyArg::OptMscat(), Bool_t sw_eloss = PhyArg::OptEloss());
        
        void set_state_with_cos(Double_t cx, Double_t cy, Double_t cz, Double_t ux = 0., Double_t uy = 0., Double_t uz = -1.);
        void set_state_with_tan(Double_t cx, Double_t cy, Double_t cz, Double_t tx = 0., Double_t ty = 0., Double_t uz = -1.);
        void set_state_with_uxy(Double_t cx, Double_t cy, Double_t cz, Double_t ux = 0., Double_t uy = 0., Short_t signz = -1);
        
        void set_state(Double_t cx, Double_t cy, Double_t cz, Double_t mx, Double_t my, Double_t mz);
        
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
      
        inline Double_t gm() const { return ((MGNumc::EqualToZero(bta_)) ? MGMath::ONE : (gmbta_/bta_)); }

        inline Short_t  eta_sign() const { return (MGNumc::Compare(eta_)); }
        inline Double_t eta_abs()  const { return std::fabs(eta_); }

        inline const Double_t& irig()  const { return irig_; }
        inline Double_t rig() const { return (MGNumc::EqualToZero(irig_) ? MGMath::ZERO : MGMath::ONE / irig_); }

        inline const SVecD<3>& c()  const { return coo_; } 
        inline const Double_t& cx() const { return coo_(0); } 
        inline const Double_t& cy() const { return coo_(1); } 
        inline const Double_t& cz() const { return coo_(2); } 
        
        inline const SVecD<3>& u()  const { return dir_; }
        inline const Double_t& ux() const { return dir_(0); } 
        inline const Double_t& uy() const { return dir_(1); } 
        inline const Double_t& uz() const { return dir_(2); } 
        
        inline Double_t tx() const { return ((MGNumc::EqualToZero(dir_(2))) ? 0. : dir_(0)/dir_(2)); } 
        inline Double_t ty() const { return ((MGNumc::EqualToZero(dir_(2))) ? 0. : dir_(1)/dir_(2)); } 
  
        inline PhyArg&       arg() { return arg_; }
        inline VirtualPhySt& vst() { return vst_; }
        
        inline const Bool_t& field() const { return arg_(); }

        void symbk(Bool_t is_rndm = false);

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

        PhyArg       arg_;
        VirtualPhySt vst_;
};
        

} // namespace TrackSys


#endif // __TRACKLibs_PhySt_H__
