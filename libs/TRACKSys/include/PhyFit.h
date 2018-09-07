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
        
        TrFitPar(const PartInfo& info = PartInfo(PartType::Proton), const Orientation& ortt = Orientation::kDownward, const Bool_t& sw_mscat = PhyArg::OptMscat(), const Bool_t& sw_eloss = PhyArg::OptEloss());
        ~TrFitPar() { TrFitPar::clear(); }

    public :
        inline Bool_t check() { return check_hits(); }

        inline void add_hit(const HitStTRK&  hit) { hits_TRK_.push_back(hit); zero(); }
        inline void add_hit(const HitStTOF&  hit) { hits_TOF_.push_back(hit); zero(); }
        inline void add_hit(const HitStRICH& hit) { hits_RICH_.push_back(hit); zero(); }
        inline void add_hit(const HitStTRD&  hit) { hits_TRD_.push_back(hit); zero(); }
        
        inline void add_hit(const std::vector<HitStTRK>& hits) { hits_TRK_.insert(hits_TRK_.end(), hits.begin(), hits.end()); zero(); }
        inline void add_hit(const std::vector<HitStTOF>& hits) { hits_TOF_.insert(hits_TOF_.end(), hits.begin(), hits.end()); zero(); }
       
        inline const PartInfo&    info() const { return info_; }
        inline const Orientation& ortt() const { return ortt_; }

        inline const Short_t& nseq() const { return nseq_; }
        inline const Short_t& nseg() const { return nseg_; }

        inline Short_t nhit() const { return hits_.size(); }
        inline const std::vector<VirtualHitSt*>& hits() const { return hits_; }
        inline const VirtualHitSt* hit(Int_t idx) const { return ((idx<0 || idx>=hits_.size()) ? nullptr : hits_.at(idx)); }
        
        inline const std::vector<HitStTRK>&  hitsTRK()  const { return hits_TRK_; }
        inline const std::vector<HitStTOF>&  hitsTOF()  const { return hits_TOF_; }
        inline const std::vector<HitStRICH>& hitsRICH() const { return hits_RICH_; }
        inline const std::vector<HitStTRD>&  hitsTRD()  const { return hits_TRD_; }

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
        std::vector<HitStRICH>     hits_RICH_;
        std::vector<HitStTRD>      hits_TRD_;

        Short_t nseq_;
        Short_t nseg_;

        Short_t nmes_;
        Short_t nmes_cx_;
        Short_t nmes_cy_;
        Short_t nmes_ib_;
        Short_t nmes_TRKqx_;
        Short_t nmes_TRKqy_;
        Short_t nmes_TOFq_;
        Short_t nmes_TOFt_;
        Short_t nmes_RICHib_;
        Short_t nmes_TRDel_;

    private :
        Bool_t  is_check_;

    protected :
        // Number of Hit Requirement
        static constexpr Short_t LMTN_CX = 2;
        static constexpr Short_t LMTN_CY = 3;
        static constexpr Short_t LMTN_TOF_T = 1;
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
        
        inline const Short_t& ndof()    const { return ndof_; }
        inline const Short_t& ndof_cx() const { return ndof_cx_; }
        inline const Short_t& ndof_cy() const { return ndof_cy_; }
        inline const Short_t& ndof_ib() const { return ndof_ib_; }
        
        inline const Double_t& nchi()    const { return nchi_; }
        inline const Double_t& nchi_cx() const { return nchi_cx_; }
        inline const Double_t& nchi_cy() const { return nchi_cy_; }
        inline const Double_t& nchi_ib() const { return nchi_ib_; }

    protected :
        void clear();

        Bool_t analyticalFit();
        Bool_t simpleFit();

    protected :
        Bool_t succ_;
        PhySt  part_;

        Short_t ndof_;
        Short_t ndof_cx_; // (cx)
        Short_t ndof_cy_; // (cy)
        Short_t ndof_ib_; // (TRKq + TOFt + TOFq + RICHib + TRDel)
        
        Double_t nchi_;
        Double_t nchi_cx_;
        Double_t nchi_cy_;
        Double_t nchi_ib_;

    protected :
        static constexpr Short_t DIMG = 5;

        // Minimization (Analytical Method)
        static constexpr Double_t LMTU_ETA = 5.0;
        
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
        VirtualPhyTrFit(const TrFitPar& fitPar, const PhySt& part, Bool_t opt_mu = false) : 
            TrFitPar(fitPar), part_(part), opt_mu_(opt_mu), opt_int_(sw_mscat_), opt_tsft_(nmes_TOFt_>=LMTN_TOF_T),
            DIMG_(5+opt_mu_+opt_tsft_), DIML_(4), 
            numOfRes_(0), numOfParGlb_(0), numOfParInt_(0),
            parIDeta_(-1), parIDigb_(-1), parIDtsft_(-1)
            { if (check_hits()) { setvar(nseq_, (opt_int_?nseg_:0)); } }
    
    public :
        virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;

    protected :
        inline void setvar(const Short_t nseq = 0, const Short_t nseg = 0) {
            if (nseq <= 0 || nseg < 0) return;
            numOfRes_    = nseq + nseg * DIML_;
            numOfParGlb_ = DIMG_;
            numOfParInt_ = nseg * DIML_;

            parIDeta_ = 4;
            if (opt_mu_)   parIDigb_  = 5;
            if (opt_tsft_) parIDtsft_ = DIMG_ - 1;

            set_num_residuals(numOfRes_);
            mutable_parameter_block_sizes()->clear(); 
            mutable_parameter_block_sizes()->push_back(numOfParGlb_);
            if (opt_int_) mutable_parameter_block_sizes()->push_back(numOfParInt_);
        }
    
    protected :
        const Bool_t opt_mu_;
        const Bool_t opt_int_;
        const Bool_t opt_tsft_;
        const PhySt  part_;
        
        const Short_t DIMG_;
        const Short_t DIML_;
        
        Short_t numOfRes_;
        Short_t numOfParGlb_;
        Short_t numOfParInt_;

        Short_t parIDeta_;
        Short_t parIDigb_;
        Short_t parIDtsft_;
};


class PhyTrFit : public TrFitPar {
    public :
        enum class MuOpt { kFixed = 0, kFree = 1 };
        
    protected :
        static constexpr Short_t DIMG = 5;
        static constexpr Short_t DIML = 4;

    public :
        PhyTrFit& operator=(const PhyTrFit& rhs);
        PhyTrFit(const PhyTrFit& trFit) { *this = trFit; }
        
        PhyTrFit(const TrFitPar& fitPar, const MuOpt& mu_opt = MuOpt::kFixed);
        ~PhyTrFit() { PhyTrFit::clear(); }
        
    public :
        inline const MuOpt&  muOpt()  const { return mu_opt_; }
        inline const Bool_t& status() const { return succ_; }
        inline const PhySt&  part()   const { return part_; }
        
        inline const Short_t&  ndof(Int_t it)    const { return ndof_.at(it); }
        inline const Double_t& nchi(Int_t it)    const { return nchi_.at(it); }
        inline const Double_t& quality(Int_t it) const { return quality_.at(it); }

        inline const Short_t& ndof_cx() const { return ndof_cx_; }
        inline const Short_t& ndof_cy() const { return ndof_cy_; }
        inline const Short_t& ndof_ib() const { return ndof_ib_; }
        
        inline const Double_t& nchi_cx() const { return nchi_cx_; }
        inline const Double_t& nchi_cy() const { return nchi_cy_; }
        inline const Double_t& nchi_ib() const { return nchi_ib_; }
    
        inline const Double_t& err_cx()  const { return err_[0]; }
        inline const Double_t& err_cy()  const { return err_[1]; }
        inline const Double_t& err_ux()  const { return err_[2]; }
        inline const Double_t& err_uy()  const { return err_[3]; }
        inline const Double_t& err_eta() const { return err_[4]; }
        inline const Double_t& err_igb() const { return err_[5]; }
        inline const Double_t& err_mu()  const { return err_[6]; }
        
        inline Double_t err_irig() const { return (err_[4] / PartInfo::ATOMIC_MASS); }
        inline Double_t err_mass() const { return (err_[6] * PartInfo::ATOMIC_MASS); }

    protected :
        void clear();
       
        inline void resetPhyArg(Bool_t sw_mscat, Bool_t sw_eloss) { sw_mscat_ = sw_mscat; sw_eloss_ = sw_eloss; part_.arg().reset(sw_mscat_, sw_eloss_); }
        Bool_t survivalTestAndModify();

        Bool_t simpleFit();
        Bool_t physicalFit(const MuOpt& mu_opt = MuOpt::kFixed, Double_t fluc_eta = Numc::ZERO<>, Double_t fluc_igb = Numc::ZERO<>, Bool_t with_mu_est = true);
        Bool_t physicalTrFit();
        Bool_t physicalMuFit();

        Bool_t evolve(const MuOpt& mu_opt = MuOpt::kFixed);

    protected :
        MuOpt               mu_opt_;
        Bool_t              succ_;
        PhySt               part_;
        std::vector<PhyArg> args_; 
        Double_t            TOFt_sft_; // TOF time shift [cm]

        std::array<Short_t,  2> ndof_;
        std::array<Double_t, 2> nchi_;
        std::array<Double_t, 2> quality_;

        Short_t ndof_tt_; // (total)
        Short_t ndof_cx_; // (cx + mstau)
        Short_t ndof_cy_; // (cy + msrho)
        Short_t ndof_ib_; // (TRKq + TOFt + TOFq + RICHib + TRDel)

        Double_t nchi_tt_;
        Double_t nchi_cx_;
        Double_t nchi_cy_;
        Double_t nchi_ib_;
        
        std::array<Double_t, 7> err_; // (cx, cy, ux, uy, eta, igb, mu)

    public :
        PhySt interpolate_to_z(Double_t zcoo = 0) const;
        MatFld get_mat(Double_t zbd1 = 0, Double_t zbd2 = 0) const;

    protected :
        std::vector<PhySt> stts_;

    private :
        static constexpr Short_t  SURVIVAL_LMTN = 12;
        static constexpr Double_t SURVIVAL_FACT = 0.85;
        static constexpr Double_t SURVIVAL_BETA = 0.25;
        static constexpr Double_t SURVIVAL_MASS = 0.000510999;
        static constexpr Double_t LMTL_INV_GB   = 1.0e-12;
        static constexpr Double_t LMTU_INV_GB   = 1.0e+3;
        static constexpr Short_t  LMT_MU_ITER   = 2;
        static constexpr Double_t MU_FLUC_BASE  = 3.00e-1;
        static constexpr Double_t MU_FLUC       = 5.00e-3;

        static Double_t NormQuality(Double_t nchi, Short_t ndof);
};


} // namespace TrackSys


#endif // __TRACKLibs_PhyFit_H__
