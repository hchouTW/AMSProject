#ifndef __TRACKLibs_PhyFit_H__
#define __TRACKLibs_PhyFit_H__


#include "ceres/ceres.h" // Ceres-Solver


namespace TrackSys {


class TrFitPar {
    public :
        enum class Orientation {
            kDownward = 0, kUpward = 1
        };

    public :
        TrFitPar& operator=(const TrFitPar& rhs);
        TrFitPar(const TrFitPar& fitPar) { *this = fitPar; }
        
        TrFitPar(const PartInfo& info = PartInfo(PartType::Proton), const Orientation& ortt = Orientation::kDownward, Bool_t sw_mscat = PhyArg::OptMscat(), Bool_t sw_eloss = PhyArg::OptEloss());
        ~TrFitPar() { TrFitPar::clear(); }

    public :
        inline Bool_t check() { return check_hits(); }

        inline void add_hit(const HitStTRK& hit) { hits_TRK_.push_back(hit); zero(); }
        inline void add_hit(const HitStTOF& hit) { hits_TOF_.push_back(hit); zero(); }
        
        inline void add_hit(const std::vector<HitStTRK>& hits) { hits_TRK_.insert(hits_TRK_.end(), hits.begin(), hits.end()); zero(); }
        inline void add_hit(const std::vector<HitStTOF>& hits) { hits_TOF_.insert(hits_TOF_.end(), hits.begin(), hits.end()); zero(); }
       
        inline const PartInfo&    info() const { return info_; }
        inline const Orientation& ortt() const { return ortt_; }

        inline const Short_t& nseq() const { return nseq_; }
        inline const Short_t& nseg() const { return nseg_; }

        inline const Short_t nhit() const { return hits_.size(); }
        inline const std::vector<VirtualHitSt*>& hits() const { return hits_; }
        inline const VirtualHitSt* hit(Int_t idx) const { return ((idx<0 || idx>=hits_.size()) ? nullptr : hits_.at(idx)); }

    protected :
        void zero();
        void clear();
        
        Bool_t sort_hits();
        Bool_t check_hits();

    protected :
        Bool_t      sw_mscat_;
        Bool_t      sw_eloss_;
        PartInfo    info_;
        Orientation ortt_;

        std::vector<VirtualHitSt*> hits_;
        std::vector<HitStTRK>      hits_TRK_;
        std::vector<HitStTOF>      hits_TOF_;

        Short_t nseq_;
        Short_t nseg_;

        Short_t nmes_;
        Short_t nmes_cx_;
        Short_t nmes_cy_;
        Short_t nmes_TRKqx_;
        Short_t nmes_TRKqy_;
        Short_t nmes_TOFq_;
        Short_t nmes_TOFt_;

    private :
        Bool_t  is_check_;

    protected :
        // Number of Hit Requirement
        static constexpr Short_t LMTN_CX = 3;
        static constexpr Short_t LMTN_CY = 4;
        static constexpr Short_t LMTN_TOF_T = 2;
};


class SimpleTrFit : public TrFitPar {
    public :
        SimpleTrFit& operator=(const SimpleTrFit& rhs);
        SimpleTrFit(const SimpleTrFit& trFit) { *this = trFit; }
        
        SimpleTrFit(const TrFitPar& fitPar); 
        ~SimpleTrFit() { SimpleTrFit::clear(); }
        
    public :
        inline const Bool_t& status() const { return succ_; }
        inline const PhySt& part() const { return part_; }
        
        inline const Short_t& ndof()      const { return ndof_; }
        inline const Short_t& ndof_cx()   const { return ndof_cx_; }
        inline const Short_t& ndof_cy()   const { return ndof_cy_; }
        inline const Short_t& ndof_TRKq() const { return ndof_TRKq_; }
        inline const Short_t& ndof_TOFq() const { return ndof_TOFq_; }
        inline const Short_t& ndof_TOFt() const { return ndof_TOFt_; }
        
        inline const Double_t& nchi()      const { return nchi_; }
        inline const Double_t& nchi_cx()   const { return nchi_cx_; }
        inline const Double_t& nchi_cy()   const { return nchi_cy_; }
        inline const Double_t& nchi_TRKq() const { return nchi_TRKq_; }
        inline const Double_t& nchi_TOFq() const { return nchi_TOFq_; }
        inline const Double_t& nchi_TOFt() const { return nchi_TOFt_; }

    protected :
        void clear();

        Bool_t analyticalFit();
        Bool_t simpleFit();

    protected :
        Bool_t succ_;
        PhySt  part_;

        Short_t ndof_;
        Short_t ndof_cx_;
        Short_t ndof_cy_;
        Short_t ndof_TRKq_;
        Short_t ndof_TOFq_;
        Short_t ndof_TOFt_;
        
        Double_t nchi_;
        Double_t nchi_cx_;
        Double_t nchi_cy_;
        Double_t nchi_TRKq_;
        Double_t nchi_TOFq_;
        Double_t nchi_TOFt_;

    protected :
        static constexpr Short_t DIMG = 5;
        
        // Minimization (Levenberg-Marquardt Method)
        static constexpr Short_t LMTL_ITER = 3;
        static constexpr Short_t LMTU_ITER = 25;
        
        static constexpr Double_t LAMBDA0 = 1.0e-2;
        static constexpr Double_t LAMBDA_DN_FAC = 5.0;
        static constexpr Double_t LAMBDA_UP_FAC = 7.0;
        static constexpr Double_t LMTL_LAMBDA = 1.0e-4;
        static constexpr Double_t LMTU_LAMBDA = 1.0e+3;
        
        static constexpr Double_t CONVG_EPSILON   = 3.0e-3;
        static constexpr Double_t CONVG_TOLERANCE = 7.0e-3;
};


class VirtualPhyTrFit : protected TrFitPar, public ceres::CostFunction {
    public :
        VirtualPhyTrFit(const TrFitPar& fitPar, const PhySt& part, Bool_t is_mass_fixed = true, Double_t lmtdiv_invu = Numc::ZERO<>) : 
            TrFitPar(fitPar), part_(part), is_mass_fixed_(is_mass_fixed),
            is_invu_bound_(Numc::Compare(lmtdiv_invu)>0), lmtdiv_invu_(is_invu_bound_?lmtdiv_invu:Numc::ZERO<>),
            DIMG_(Numc::SIX<Short_t>-is_mass_fixed), DIML_(Numc::FOUR<Short_t>), 
            numOfRes_(0), numOfPar_(0) { if (check_hits()) setvar(nseq_+nseg_*DIML_+is_invu_bound_, DIMG_+nseg_*DIML_); }
    
    public :
        virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;

    protected :
        inline void setvar(const Short_t num_of_residual = 0, const Short_t num_of_parameter = 0) {
            mutable_parameter_block_sizes()->clear(); 
            if (num_of_residual  > 0) { numOfRes_ = num_of_residual;  set_num_residuals(num_of_residual); }
            if (num_of_parameter > 0) { numOfPar_ = num_of_parameter; mutable_parameter_block_sizes()->push_back(num_of_parameter); }
        }
    
    protected :
        const Short_t DIMG_;
        const Short_t DIML_;
       
        const Bool_t   is_mass_fixed_;
        const Bool_t   is_invu_bound_;
        const Double_t lmtdiv_invu_;
        const PhySt    part_;
        
        Short_t numOfRes_;
        Short_t numOfPar_;
};


class PhyTrFit : public TrFitPar {
    public :
        enum class MassOpt {
            kFixed = 0, kFree = 1
        };
        
    protected :
        static constexpr Short_t DIMG_ = 6;
        static constexpr Short_t DIML_ = 4;

    public :
        PhyTrFit& operator=(const PhyTrFit& rhs);
        PhyTrFit(const PhyTrFit& trFit) { *this = trFit; }
        
        PhyTrFit(const TrFitPar& fitPar, const MassOpt& mass_opt = MassOpt::kFixed);
        ~PhyTrFit() { PhyTrFit::clear(); }
        
    public :
        inline const Bool_t& status() const { return succ_; }
        inline const PhySt& part() const { return part_; }
        
        inline const Double_t& quality() const { return quality_; }

        inline const Short_t& ndof()      const { return ndof_; }
        inline const Short_t& ndof_cx()   const { return ndof_cx_; }
        inline const Short_t& ndof_cy()   const { return ndof_cy_; }
        inline const Short_t& ndof_TRKq() const { return ndof_TRKq_; }
        inline const Short_t& ndof_TOFq() const { return ndof_TOFq_; }
        inline const Short_t& ndof_TOFt() const { return ndof_TOFt_; }
        
        inline const Double_t& nchi()      const { return nchi_; }
        inline const Double_t& nchi_cx()   const { return nchi_cx_; }
        inline const Double_t& nchi_cy()   const { return nchi_cy_; }
        inline const Double_t& nchi_TRKq() const { return nchi_TRKq_; }
        inline const Double_t& nchi_TOFq() const { return nchi_TOFq_; }
        inline const Double_t& nchi_TOFt() const { return nchi_TOFt_; }
        
        inline const Double_t& nrm_mstau() const { return nrm_mstau_; }
        inline const Double_t& nrm_msrho() const { return nrm_msrho_; }
        inline const Double_t& nrm_elion() const { return nrm_elion_; }
    
        inline const SVecD<DIMG_>& errG()      const { return errG_; }
        inline const Double_t&     errG_cx()   const { return errG_(0); }
        inline const Double_t&     errG_cy()   const { return errG_(1); }
        inline const Double_t&     errG_ux()   const { return errG_(2); }
        inline const Double_t&     errG_uy()   const { return errG_(3); }
        inline const Double_t&     errG_eta()  const { return errG_(4); }
        inline const Double_t&     errG_invu() const { return errG_(5); }

        inline const std::vector<SVecD<DIML_>>& errL() const { errL_; }
        inline const SVecD<DIML_>& errL(Int_t it)      const { return errL_.at(it); }
        inline const Double_t&     errL_tauu(Int_t it) const { return errL_.at(it)(0); }
        inline const Double_t&     errL_rhou(Int_t it) const { return errL_.at(it)(1); }
        inline const Double_t&     errL_taul(Int_t it) const { return errL_.at(it)(2); }
        inline const Double_t&     errL_rhol(Int_t it) const { return errL_.at(it)(3); }

    protected :
        void clear();

        Bool_t simpleFit();
        Bool_t physicalFit(const MassOpt& mass_opt = MassOpt::kFixed, Double_t sgm_invu = 0);
        Bool_t physicalMassFit();

        Bool_t evolve();

    protected :
        Bool_t              succ_;
        PhySt               part_;
        std::vector<PhyArg> args_; 

        Double_t lmtl_invu_;
        Double_t lmtu_invu_;
        Double_t lmtdiv_invu_;

        Double_t quality_;

        Short_t ndof_;
        Short_t ndof_cx_;
        Short_t ndof_cy_;
        Short_t ndof_TRKq_;
        Short_t ndof_TOFq_;
        Short_t ndof_TOFt_;
        
        Double_t nchi_;
        Double_t nchi_cx_;
        Double_t nchi_cy_;
        Double_t nchi_TRKq_;
        Double_t nchi_TOFq_;
        Double_t nchi_TOFt_;

        Double_t nrm_mstau_;
        Double_t nrm_msrho_;
        Double_t nrm_elion_;

        SVecD<DIMG_>              errG_; // (cx, cy, ux, uy, eta, invu)
        std::vector<SVecD<DIML_>> errL_; // (tauu, rhou, taul, rhol)

    public :
        PhySt interpolate_to_z(Double_t zcoo = 0) const;
        MatFld get_mat(Double_t zbd1 = 0, Double_t zbd2 = 0) const;

    protected :
        std::vector<PhySt> stts_;

    private :
        static Double_t NormQuality(Double_t nchi, Short_t ndof = Numc::THREE<Short_t>) {
            if (ndof < Numc::THREE<Short_t>) return Numc::ZERO<>;
            Double_t chi   = nchi * static_cast<Double_t>(ndof);
            if (Numc::EqualToZero(chi)) return Numc::ZERO<>;
            Double_t qmin  = static_cast<Double_t>(ndof - Numc::TWO<Short_t>);
            Double_t sign  = static_cast<Double_t>(Numc::Compare(chi - qmin));
            Double_t qfunc = (chi - qmin) - qmin * std::log(chi / qmin);
            if (!Numc::Valid(qfunc)) return Numc::ZERO<>;
            Double_t xfunc = sign * std::sqrt(qfunc / static_cast<Double_t>(ndof));
            if (Numc::Valid(xfunc)) return xfunc;
            return Numc::ZERO<>;
        }
};


} // namespace TrackSys


#endif // __TRACKLibs_PhyFit_H__
