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
        TrFitPar(const PartType& type = PartType::Proton, const Orientation& ortt = Orientation::kDownward, Bool_t sw_mscat = PhyArg::OptMscat(), Bool_t sw_eloss = PhyArg::OptEloss());
        ~TrFitPar() { TrFitPar::clear(); }

        void print() const;

        inline void addHit(HitStTRK& hit) { hits_TRK_.push_back(hit); hits_.clear(); is_check_ = false; }
        inline void addHit(HitStTOF& hit) { hits_TOF_.push_back(hit); hits_.clear(); is_check_ = false; }
        
        inline Bool_t check() { return check_hits(); }

        inline const Short_t& nseq() const { return nseq_; }
        
        inline const Short_t nhits() const { return hits_.size(); }

    protected :
        void zero();
        void clear();
        
        Bool_t sort_hits();
        Bool_t check_hits();
        Bool_t recheck_hits() { is_check_ = false; return check_hits(); }

    protected :
        Bool_t      sw_mscat_;
        Bool_t      sw_eloss_;
        PartType    type_;
        Orientation ortt_;

        std::vector<VirtualHitSt*> hits_;
        std::vector<HitStTRK>      hits_TRK_;
        std::vector<HitStTOF>      hits_TOF_;

        Short_t nseq_;
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


class SimpleTrFit : protected TrFitPar {
    public :
        SimpleTrFit(TrFitPar& fitPar); 
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
        // Minimization (Levenberg-Marquardt Method)
        static constexpr Int_t LMTL_ITER = 3;
        static constexpr Int_t LMTU_ITER = 25;
        
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
        VirtualPhyTrFit(TrFitPar& fitPar, PhySt& part) : TrFitPar(fitPar), numOfRes_(0), numOfPar_(0), part_(part) { if (recheck_hits()) setvar(nseq_+(nhits()-1)*PhyJb::DIM_L, PhyJb::DIM_G+(nhits()-1)*PhyJb::DIM_L); }
        ~VirtualPhyTrFit() { VirtualPhyTrFit::clear(); }
    
    public :
        virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;
    
    protected :
        inline void setvar(Int_t num_of_residual = 0, Int_t num_of_parameter = 0) {
            numOfRes_ = 0;
            numOfPar_ = 0;
            mutable_parameter_block_sizes()->clear(); 
            if (num_of_residual  > 0) { numOfRes_ = num_of_residual;  set_num_residuals(num_of_residual); }
            if (num_of_parameter > 0) { numOfPar_ = num_of_parameter; mutable_parameter_block_sizes()->push_back(num_of_parameter); }
        }

        inline void clear() { part_.reset(type_); part_.arg().reset(sw_mscat_, sw_eloss_); setvar(); }
    
    protected :
        Int_t numOfRes_;
        Int_t numOfPar_;
        PhySt part_;
};


class PhyTrFit : protected TrFitPar {
    public :
        PhyTrFit(TrFitPar& fitPar);
        ~PhyTrFit() { PhyTrFit::clear(); }
        
    public :
        inline const Bool_t& status() const { return succ_; }
        inline const PhySt& part() const { return part_; }

        inline const Int_t nargs() const { return args_.size(); }
        inline const std::vector<PhyArg>& args() const { return args_; }
        inline const PhyArg* args(UInt_t it) const { return ((succ_ && it<args_.size()) ? &args_.at(it) : nullptr); }

        inline const Int_t nstts() const { return stts_.size(); }
        inline const std::vector<PhySt>& stts() const { return stts_; }
        inline const PhySt* stts(Int_t lay) const { auto&& stt = map_stts_.find(lay); return ((succ_ && stt!=map_stts_.end()) ? stt->second : nullptr); }
        
        inline const Int_t nhits() const { return hits_.size(); }
        inline const std::vector<VirtualHitSt*>& hits() const { return hits_; }
        inline const VirtualHitSt* hits(Int_t lay) const { auto&& hit = map_hits_.find(lay); return ((succ_ && hit!=map_hits_.end()) ? hit->second : nullptr); }

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
        
        inline const Short_t&  nsegs()     const { return nsegs_; }
        inline const Double_t& nrm_mscat() const { return nrm_mscat_; }
        inline const Double_t& nrm_eloss() const { return nrm_eloss_; }

    protected :
        void clear();

        Bool_t simpleFit();
        Bool_t physicalFit();

        Bool_t evolve();

    protected :
        Bool_t              succ_;
        PhySt               part_;
        std::vector<PhyArg> args_; 
        std::vector<PhySt>  stts_; 
        
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

        Short_t  nsegs_;
        Double_t nrm_mscat_;
        Double_t nrm_eloss_;

    protected :
        std::map<Int_t, VirtualHitSt*> map_hits_;
        std::map<Int_t, PhySt*>        map_stts_;

    protected :
        // (z, hit, part, args)
        //std::vector<std::tuple<Double_t, VirtualHitSt*, PhySt, SVecD<5>>> stts;
        // Interpolate form -inf to inf; return part;
        // Material return nrl, ela
};


} // namespace TrackSys


#endif // __TRACKLibs_PhyFit_H__
