#ifndef __TRACKLibs_Prop_H__
#define __TRACKLibs_Prop_H__

namespace TrackSys {


//---- OrthCoord ----//
// Orthogonal Coordinate based on a seed vector
// Org is particle direction
// Axis Org := Tau x Rho
class OrthCoord {
    public :
        static const SVecD<3> AXIS_X;
        static const SVecD<3> AXIS_Y;
        static const SVecD<3> AXIS_Z;

    public :
        OrthCoord() : org_(0, 0, 1), tau_(1, 0, 0), rho_(0, 1, 0) {}
        OrthCoord(const SVecD<3>& org, const SVecD<3>& seed = AXIS_X) : OrthCoord() { reset(org, seed); }
        ~OrthCoord() {}

        inline void reset(const SVecD<3>& org, const SVecD<3>& seed = AXIS_X);

        inline const SVecD<3>& org() const { return org_; }
        inline const SVecD<3>& tau() const { return tau_; }
        inline const SVecD<3>& rho() const { return rho_; }
        
        inline const Double_t& ox() const { return org_(0); }
        inline const Double_t& oy() const { return org_(1); }
        inline const Double_t& oz() const { return org_(2); }
        
        inline const Double_t& tx() const { return tau_(0); }
        inline const Double_t& ty() const { return tau_(1); }
        inline const Double_t& tz() const { return tau_(2); }

        inline const Double_t& rx() const { return rho_(0); }
        inline const Double_t& ry() const { return rho_(1); }
        inline const Double_t& rz() const { return rho_(2); }

    private :
        SVecD<3> org_;
        SVecD<3> tau_;
        SVecD<3> rho_;
};
        
const SVecD<3> OrthCoord::AXIS_X(1, 0, 0);
const SVecD<3> OrthCoord::AXIS_Y(0, 1, 0);
const SVecD<3> OrthCoord::AXIS_Z(0, 0, 1);


class MotionFunc {
    public :
        MotionFunc(PhySt& part, const MatPhyFld* mphy = nullptr);
        ~MotionFunc() {}
        
        inline const OrthCoord& orth() const { return orth_; }
        
        inline const SVecD<3>& c() const { return zeta_c_; }
        inline const SVecD<3>& u() const { return zeta_u_; }
        inline const Double_t& e() const { return zeta_e_; }
        
        inline const Double_t& cx() const { return zeta_c_(0); }
        inline const Double_t& cy() const { return zeta_c_(1); }
        inline const Double_t& cz() const { return zeta_c_(2); }

        inline const Double_t& ux() const { return zeta_u_(0); }
        inline const Double_t& uy() const { return zeta_u_(1); }
        inline const Double_t& uz() const { return zeta_u_(2); }

    private :
        static constexpr Double_t PROP_FACT = 2.99792458e-04;
        
        OrthCoord orth_;
        
        SVecD<3> zeta_c_;
        SVecD<3> zeta_u_;
        Double_t zeta_e_;
};


class TransferFunc {
    public :
        TransferFunc(PhySt& part, const MatPhyFld* mphy = nullptr);
        ~TransferFunc() {}
        
        inline const SVecD<3>&    cu() const { return kappa_cu_; }
        inline const Double_t&    cu(Int_t i) const { return kappa_cu_(i); }
        
        inline const SMtxD<3, 3>& uu() const { return kappa_uu_; }
        inline const Double_t&    uu(Int_t i, Int_t j) const { return kappa_uu_(i, j); }
        
        inline const SVecD<3>&    ue() const { return kappa_ue_; }
        inline const Double_t&    ue(Int_t i) const { return kappa_ue_(i); }
        
        inline const Double_t&    ee() const { return kappa_ee_; }
    
    private :
        static constexpr Double_t PROP_FACT = 2.99792458e-04;
        
        SVecD<3>    kappa_cu_; // d(dp/ds) / du
        
        SMtxD<3, 3> kappa_uu_; // d(du/ds) / du
        SVecD<3>    kappa_ue_; // d(du/ds) / de
        
        Double_t    kappa_ee_; // d(de/ds) / de
};


class PhyJb {
    public :
        static constexpr Int_t DIM_G = 5;
        static constexpr Int_t DIM_L = 2;

        using SMtxDGG = SMtxD<DIM_G, DIM_G>;
        using SMtxDGL = SMtxD<DIM_G, DIM_L>;
        
        using SMtxDXYG = SMtxD<2, DIM_G>;
        using SMtxDXYL = SMtxD<2, DIM_L>;
    
    public :
        PhyJb() { init(); }
        ~PhyJb() {}

        inline void init();
        
        inline void set(PhySt& part);

        inline void multiplied(PhyJb& phyJb);

        inline void print() const;

        inline const Bool_t& field() const { return field_; }

        inline SMtxDGG& gg() { return jb_gg_; }
        inline SMtxDGL& gl() { return jb_gl_; }
        
        inline Double_t& gg(Int_t i, Int_t j) { return jb_gg_(i, j); }
        inline Double_t& gl(Int_t i, Int_t j) { return jb_gl_(i, j); }

        inline SMtxDXYG xyg() { return jb_gg_.Sub<SMtxDXYG>(0, 0); }
        inline SMtxDXYL xyl() { return jb_gl_.Sub<SMtxDXYL>(0, 0); }

    private :
        Bool_t    field_;
        SMtxDGG   jb_gg_;
        SMtxDGL   jb_gl_;
    
    private :
        static constexpr Short_t X = 0;
        static constexpr Short_t Y = 1;
        static constexpr Short_t E = 2;
        static constexpr Short_t JPX = 0;
        static constexpr Short_t JPY = 1;
        static constexpr Short_t JUX = 2;
        static constexpr Short_t JUY = 3;
        static constexpr Short_t JEA = 4;
        static constexpr Short_t JTAU = 0;
        static constexpr Short_t JRHO = 1;
        static constexpr Short_t JION = 2;
        static constexpr Short_t JBRM = 3;
};


class TransferPhyJb {
    public :
        TransferPhyJb(const TransferFunc& tf, PhyJb& jb);
        ~TransferPhyJb() {}

        inline const SMtxD<2, 2>& uu() const { return uu_; }
        inline const Double_t&    uu(Int_t i, Int_t j) const { return uu_(i, j); }
        
        inline const SVecD<2>& ue() const { return ue_; }
        inline const Double_t& ue(Int_t i) const { return ue_(i); }
        
        inline const Double_t& ee() const { return ee_; }

    private :
        SMtxD<2, 2> uu_;
        SVecD<2>    ue_;
        Double_t    ee_;
    
    private :
        static constexpr Short_t X = 0;
        static constexpr Short_t Y = 1;
        static constexpr Short_t E = 2;
        static constexpr Short_t JPX = 0;
        static constexpr Short_t JPY = 1;
        static constexpr Short_t JUX = 2;
        static constexpr Short_t JUY = 3;
        static constexpr Short_t JEA = 4;
};


class PropPhyCal {
    public :
        PropPhyCal(PhySt& part, Double_t sign = 1.) { init(); sw_mscat_ = part.arg().mscat(); sw_eloss_ = part.arg().eloss(); eta0_abs_ = part.eta_abs(); sign_ = ((MGNumc::Compare(sign)>=0)?1:-1); }
        ~PropPhyCal() {}

        void init(); 
        void normalized(const MatFld& mfld, const PhySt& part);

        void push(const PhySt& part, const MatPhyFld& mpfld, const SVecD<3>& tau, const SVecD<3>& rho);

        //Bool_t operator() () const { return (sw_mscat_ || sw_eloss_); }

        const Short_t& sign() const { return sign_; }
        
        const Double_t& len() const { return len_; }
        const Double_t& nrl() const { return nrl_; }

        const SVecD<3>& tau() const { return tau_; }
        const SVecD<3>& rho() const { return rho_; }
        
        const Double_t& mscatu()  const { return mscatu_; }
        const Double_t& mscatcu() const { return mscatcu_; }
        const Double_t& mscatcl() const { return mscatcl_; }

        const Double_t& eloss_ion_kpa() const { return eloss_ion_kpa_; }
        const Double_t& eloss_ion_mpv() const { return eloss_ion_mpv_; }
        const Double_t& eloss_ion_sgm() const { return eloss_ion_sgm_; }

        const Double_t& eloss_brm_men() const { return eloss_brm_men_; }

        inline void set_virtualPhySt(PhySt& part) const;

    private :
        Bool_t   sw_mscat_;
        Bool_t   sw_eloss_;
        Double_t eta0_abs_;
       
        Bool_t   field_;
        Short_t  sign_;
        Double_t len_;
        Double_t nrl_;

        std::vector<Bool_t>   vec_vac_;
        std::vector<Double_t> vec_len_;
        std::vector<Double_t> vec_eft_;
        std::vector<Double_t> vec_invloc_;
        std::vector<Double_t> vec_invlocsqr_;
        std::vector<Double_t> vec_mscat_;

        SVecD<3> tau_;
        SVecD<3> rho_;
        
        Double_t mscatu_;
        Double_t mscatcu_;
        Double_t mscatcl_;

        std::vector<Double_t> vec_ieta_kpa_;
        std::vector<Double_t> vec_ieta_sgm_;

        Double_t eloss_ion_kpa_;
        Double_t eloss_ion_mpv_;
        Double_t eloss_ion_sgm_;

        Double_t eloss_brm_men_;
};


class PropMgnt {
    public :
        PropMgnt() {}
        ~PropMgnt() {}
    
    public :
        enum class Method {
            kEuler = 1, kEulerHeun = 2, kRungeKuttaNystrom = 4
        };

        static void   SetMethod(Method method = Method::kRungeKuttaNystrom) { method_ = method; }
        static Bool_t CheckMethod(Method method) { return (method == method_); }

    private :
        static Method method_;

    public :
#ifdef __HAS_AMS_OFFICE_LIBS__
        static Bool_t PropToZ_AMSLibs(const Double_t zcoo, PhySt& part);
#endif // __HAS_AMS_OFFICE_LIBS__
        
        static Bool_t Prop(const Double_t step, PhySt& part, MatFld* mfld = nullptr, PhyJb* phyJb = nullptr);
        static Bool_t PropToZ(const Double_t zcoo, PhySt& part, MatFld* mfld = nullptr, PhyJb* phyJb = nullptr);
        
        static Bool_t PropWithMC(const Double_t step, PhySt& part, MatFld* mfld = nullptr);
        static Bool_t PropToZWithMC(const Double_t zcoo, PhySt& part, MatFld* mfld = nullptr);

    protected :
        // Step Length
        // step < 0, backward trace
        // step > 0, forward trace
        static Double_t GetPropStep(PhySt& part, Short_t ward);
        static Double_t GetStep(PhySt& part, Double_t resStep);
        static Double_t GetStepToZ(PhySt& part, Double_t resStepZ);
        
        static Bool_t PropWithEuler(const Double_t step, PhySt& part, const MatFld& mfld, PropPhyCal& ppcal, PhyJb* phyJb = nullptr); 
        static Bool_t PropWithEulerHeun(const Double_t step, PhySt& part, const MatFld& mfld, PropPhyCal& ppcal, PhyJb* phyJb = nullptr); 
        static Bool_t PropWithRungeKuttaNystrom(const Double_t step, PhySt& part, const MatFld& mfld, PropPhyCal& ppcal, PhyJb* phyJb = nullptr);
        
    private :
        static constexpr Double_t PROP_FACT  = 2.99792458e-04;
        static constexpr Double_t LMTL_CURVE = 2.0e-6; // (du/ds threshold)
        static constexpr Double_t TUNE_STEP  = 1.0e-3; // (du threshold)
        static constexpr Double_t PROP_STEP  = 30.0;   // (ds threshold)
        static constexpr Double_t LMTU_STEP  = 50.0;   // (ds threshold)
        static constexpr Double_t LMTL_STEP  =  8.0;   // (ds threshold)
        static constexpr Double_t TUNE_MAT   =  0.1;  // (number radiation length threshold)
        
        static constexpr Long64_t LMTU_ITER  = 100;
        static constexpr Double_t CONV_STEP  = 1.0e-4; // [cm]

    private :
        static constexpr Short_t X = 0;
        static constexpr Short_t Y = 1;
        static constexpr Short_t E = 2;
        static constexpr Short_t JPX = 0;
        static constexpr Short_t JPY = 1;
        static constexpr Short_t JUX = 2;
        static constexpr Short_t JUY = 3;
        static constexpr Short_t JEA = 4;
        static constexpr Short_t JTAU = 0;
        static constexpr Short_t JRHO = 1;
        static constexpr Short_t JION = 2;
        static constexpr Short_t JBRM = 3;
};

PropMgnt::Method PropMgnt::method_ = PropMgnt::Method::kRungeKuttaNystrom;
//PropMgnt::Method PropMgnt::method_ = PropMgnt::Method::kEulerHeun;
//PropMgnt::Method PropMgnt::method_ = PropMgnt::Method::kEuler;


} // namespace TrackSys

#endif // __TRACKLibs_Prop_H__
