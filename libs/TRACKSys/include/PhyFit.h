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
        void print() const;
        
        inline Bool_t check() { return check_hits(); }

        inline void add_hit(const HitStTRK& hit) { hits_TRK_.push_back(hit); zero(); }
        inline void add_hit(const HitStTOF& hit) { hits_TOF_.push_back(hit); zero(); }
        
        inline void add_hit(const std::vector<HitStTRK>& hits) { hits_TRK_.insert(hits_TRK_.end(), hits.begin(), hits.end()); zero(); }
        inline void add_hit(const std::vector<HitStTOF>& hits) { hits_TOF_.insert(hits_TOF_.end(), hits.begin(), hits.end()); zero(); }
       
        inline void set_info(const PartInfo& info = PartInfo(PartType::Proton)) { info_ = info; zero(); }
        inline void set_ortt(const Orientation& ortt = Orientation::kDownward) { ortt_ = ortt; zero(); }

        inline const PartInfo&    info() const { return info_; }
        inline const Orientation& ortt() const { return ortt_; }

        inline const Short_t& nseq() const { return nseq_; }
        
        inline const Short_t nhits() const { return hits_.size(); }
        inline const std::vector<VirtualHitSt*>& hits() const { return hits_; }
        inline const VirtualHitSt* hits(Int_t idx) const { return ((idx<0 || idx>=hits_.size()) ? nullptr : hits_.at(idx)); }

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

    protected :
        static constexpr Short_t DIMG = 6;
        static constexpr Short_t DIML = 4;
};


class SimpleTrFit : protected TrFitPar {
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
        VirtualPhyTrFit(const TrFitPar& fitPar, const PhySt& part, Bool_t is_mass_fixed = true) : TrFitPar(fitPar), is_mass_fixed_(is_mass_fixed), numOfRes_(0), numOfPar_(0), part_(part) { if (check_hits()) setvar(nseq_+(nhits()-1)*DIML, DIMG+(nhits()-1)*DIML); }
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

        inline void clear() { info_ = part_.info(); is_mass_fixed_ = true; part_.arg().reset(sw_mscat_, sw_eloss_); setvar(); part_.reset(part_.info()); }
    
    protected :
        Bool_t is_mass_fixed_;
        Int_t numOfRes_;
        Int_t numOfPar_;
        PhySt part_;
};


class PhyTrFit : protected TrFitPar {
    public :
        enum class MassOpt {
            kFixed = 0, kFree = 1
        };

    public :
        PhyTrFit& operator=(const PhyTrFit& rhs);
        PhyTrFit(const PhyTrFit& trFit) { *this = trFit; }
        
        PhyTrFit(const TrFitPar& fitPar, const MassOpt& massOpt = MassOpt::kFixed);
        ~PhyTrFit() { PhyTrFit::clear(); }
        
    public :
        inline const Bool_t& status() const { return succ_; }
        inline const PhySt& part() const { return part_; }

        inline const Int_t nstts() const { return stts_.size(); }
        inline const std::vector<std::pair<VirtualHitSt*, PhySt>>& stts() const { return stts_; }
        inline const std::pair<VirtualHitSt*, PhySt>* stts(Int_t idx) const { return ((idx<0 || idx>=stts_.size()) ? nullptr : &stts_.at(idx)); }

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
        inline const Double_t& nrm_mstau() const { return nrm_mstau_; }
        inline const Double_t& nrm_msrho() const { return nrm_msrho_; }
        inline const Double_t& nrm_elion() const { return nrm_elion_; }

    public :
        PhySt interpolate_to_z(Double_t zcoo = 0) const;
        MatFld get_mat(Double_t zbd1 = 0, Double_t zbd2 = 0) const;

    protected :
        void clear();

        Bool_t simpleFit();
        Bool_t physicalFit(const MassOpt& massOpt = MassOpt::kFixed, Double_t scl = 0);

        Bool_t evolve();

    protected :
        Bool_t              succ_;
        PhySt               part_;
        std::vector<PhyArg> args_; 
        
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
        Double_t nrm_mstau_;
        Double_t nrm_msrho_;
        Double_t nrm_elion_;

    protected :
        std::vector<std::pair<VirtualHitSt*, PhySt>> stts_;
};


} // namespace TrackSys


#endif // __TRACKLibs_PhyFit_H__
