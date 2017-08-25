#ifndef __TRACKLibs_Prop_H__
#define __TRACKLibs_Prop_H__

namespace TrackSys {


//---- OrthCoord ----//
// Orthogonal Coordinate based on a seed vector
// Org is particle direction
// Axis Org := Tau x Rho
class OrthCoord {
    public :
        OrthCoord() : org_(0, 0, 1), tau_(1, 0, 0), rho_(0, 1, 0) {}
        OrthCoord(const SVecD<3>& org, const SVecD<3>& seed = AXIS_X);
        ~OrthCoord() {}
        
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
        static const SVecD<3> AXIS_X;
        static const SVecD<3> AXIS_Y;
        static const SVecD<3> AXIS_Z;

        SVecD<3> org_;
        SVecD<3> tau_;
        SVecD<3> rho_;
};
        
const SVecD<3> OrthCoord::AXIS_X(1, 0, 0);
const SVecD<3> OrthCoord::AXIS_Y(0, 1, 0);
const SVecD<3> OrthCoord::AXIS_Z(0, 0, 1);


class MotionFunc {
    public :
        MotionFunc(const PhySt& part);
        MotionFunc(const PhySt& part, const MatArg& marg, const MatPhyFld& mphy);
        ~MotionFunc() {}
        
        inline const SVecD<3>& p() const { return zeta_p_; }
        inline const SVecD<3>& u() const { return zeta_u_; }
        inline const Double_t& e() const { return zeta_e_; }
        
        inline const Double_t& px() const { return zeta_p_(0); }
        inline const Double_t& py() const { return zeta_p_(1); }
        inline const Double_t& pz() const { return zeta_p_(2); }

        inline const Double_t& ux() const { return zeta_u_(0); }
        inline const Double_t& uy() const { return zeta_u_(1); }
        inline const Double_t& uz() const { return zeta_u_(2); }
        
    private :
        static constexpr Double_t PROP_FACT = 2.99792458e-04;
        
        SVecD<3> zeta_p_;
        SVecD<3> zeta_u_;
        Double_t zeta_e_;
};


class TransferFunc {
    public :
        TransferFunc(const PhySt& part);
        TransferFunc(const PhySt& part, const MatArg& marg, const MatPhyFld& mphy);
        ~TransferFunc() {}
        
        inline const SVecD<3>&    pu() const { return kappa_pu_; }
        inline const Double_t&    pu(Int_t i) const { return kappa_pu_(i); }
        
        inline const SMtxD<3, 3>& uu() const { return kappa_uu_; }
        inline const Double_t&    uu(Int_t i, Int_t j) const { return kappa_uu_(i, j); }
        
        inline const SVecD<3>&    ue() const { return kappa_ue_; }
        inline const Double_t&    ue(Int_t i) const { return kappa_ue_(i); }
        
        inline const SVecD<3>&    ut() const { return kappa_ut_; }
        inline const Double_t&    ut(Int_t i) const { return kappa_ut_(i); }
        
        inline const SVecD<3>&    ur() const { return kappa_ur_; }
        inline const Double_t&    ur(Int_t i) const { return kappa_ur_(i); }
        
        inline const Double_t&    ee() const { return kappa_ee_; }
        inline const Double_t&    ei() const { return kappa_ei_; }
        inline const Double_t&    eb() const { return kappa_eb_; }
    
    private :
        static constexpr Double_t PROP_FACT = 2.99792458e-04;
        
        SVecD<3>    kappa_pu_; // d(dp/ds) / du
        
        SMtxD<3, 3> kappa_uu_; // d(du/ds) / du
        SVecD<3>    kappa_ue_; // d(du/ds) / de
        
        SVecD<3>    kappa_ut_; // d(du/ds) / dtau
        SVecD<3>    kappa_ur_; // d(du/ds) / drho
        
        Double_t    kappa_ee_; // d(de/ds) / de
        Double_t    kappa_ei_; // d(de/ds) / dion
        Double_t    kappa_eb_; // d(de/ds) / dbrm
};


class PhyJb {
    public :
        enum class Type {
            kZero = 0, kIdentity = 1
        };

    public :
        static constexpr Int_t DIM_G = 5;
        static constexpr Int_t DIM_L = 4;

        using SMtxDGG = SMtxD<DIM_G, DIM_G>;
        using SMtxDGL = SMtxD<DIM_G, DIM_L>;
        
        using SMtxDXYG = SMtxD<2, DIM_G>;
        using SMtxDXYL = SMtxD<2, DIM_L>;
    
    public :
        PhyJb() : mat_(false), num_rad_len_(0.0) {}
        PhyJb(const Type& type) : mat_(false), num_rad_len_(0.0) { if (type == Type::kIdentity) jb_gg_ = std::move(SMtxId()); }
        ~PhyJb() {}

        inline void init(const Type& type);

        inline void set_mat(Bool_t mat, Double_t num_rad_len = 0.);

        inline void multiplied(PhyJb& phyJb);

        inline const Bool_t& mat() const { return mat_; }

        inline const Double_t& num_rad_len() const { return num_rad_len_; }

        inline SMtxDGG& gg() { return jb_gg_; }
        inline SMtxDGL& gl() { return jb_gl_; }
        
        inline Double_t& gg(Int_t i, Int_t j) { return jb_gg_(i, j); }
        inline Double_t& gl(Int_t i, Int_t j) { return jb_gl_(i, j); }

        inline SMtxDXYG xyg() { return jb_gg_.Sub<SMtxDXYG>(0, 0); }
        inline SMtxDXYL xyl() { return jb_gl_.Sub<SMtxDXYL>(0, 0); }

    private :
        Bool_t    mat_;
        Double_t  num_rad_len_;
        SMtxDGG   jb_gg_;
        SMtxDGL   jb_gl_;
};


class TransferPhyJb {
    public :
        TransferPhyJb(Bool_t mat, const TransferFunc& tf, PhyJb& jb);
        ~TransferPhyJb() {}

        inline const SMtxD<2, 2>& uu() const { return uu_; }
        inline const Double_t&    uu(Int_t i, Int_t j) const { return uu_(i, j); }
        
        inline const SVecD<2>& ue() const { return ue_; }
        inline const Double_t& ue(Int_t i) const { return ue_(i); }
        
        inline const SVecD<2>& ut() const { return ut_; }
        inline const Double_t& ut(Int_t i) const { return ut_(i); }

        inline const SVecD<2>& ur() const { return ur_; }
        inline const Double_t& ur(Int_t i) const { return ur_(i); }

        inline const Double_t& ee() const { return ee_; }
        inline const Double_t& ei() const { return ei_; }
        inline const Double_t& eb() const { return eb_; }

    private :
        Bool_t      mat_;
        SMtxD<2, 2> uu_;
        SVecD<2>    ue_;
        SVecD<2>    ut_;
        SVecD<2>    ur_;
        Double_t    ee_;
        Double_t    ei_;
        Double_t    eb_;
    
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
        
        static Bool_t Prop(const Double_t step, PhySt& part, const MatArg& marg = MatArg(), PhyJb* phyJb = nullptr, MatFld* mfld = nullptr);
        static Bool_t PropToZ(const Double_t zcoo, PhySt& part, const MatArg& marg = MatArg(), PhyJb* phyJb = nullptr, MatFld* mfld = nullptr);
        
        static Bool_t PropWithMC(const Double_t step, PhySt& part, const MatArg& marg = MatArg(true, true));
        static Bool_t PropToZWithMC(const Double_t zcoo, PhySt& part, const MatArg& marg = MatArg(true, true));

    protected :
        // Step Length
        // step < 0, backward trace
        // step > 0, forward trace
        static Double_t GetPropStep(const PhySt& part, Short_t ward, Bool_t mat = false);
        static Double_t GetStep(const PhySt& part, Double_t resStep, Bool_t mat = false);
        static Double_t GetStepToZ(const PhySt& part, Double_t resStepZ, Bool_t mat = false);

        static Bool_t PropWithEuler(const Double_t step, PhySt& part, const MatArg& marg, const MatFld& mfld, PhyJb* phyJb = nullptr); 
        static Bool_t PropWithEulerHeun(const Double_t step, PhySt& part, const MatArg& marg, const MatFld& mfld, PhyJb* phyJb = nullptr); 
        static Bool_t PropWithRungeKuttaNystrom(const Double_t step, PhySt& part, const MatArg& marg, const MatFld& mfld, PhyJb* phyJb = nullptr); 
    
    private :
        static constexpr Double_t PROP_FACT = 2.99792458e-04;
        static constexpr Double_t LMTL_CURVE = 1.0e-6; // (du/ds threshold)
        static constexpr Double_t TUNE_STEP  = 2.0e-4; // (du threshold)
        static constexpr Double_t LMTU_STEP  = 40.0;   // (ds threshold)
        static constexpr Double_t LMTL_STEP  =  5.0;   // (ds threshold)
        static constexpr Double_t TUNE_MAT   =  0.2;   // (number radiation length threshold)

        static constexpr Long64_t LMTU_ITER = 100;
        static constexpr Double_t CONV_STEP = 1.0e-4; // [cm]

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
