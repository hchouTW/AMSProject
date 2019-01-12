#ifndef __TRACKLibs_PhyMuFit_H__
#define __TRACKLibs_PhyMuFit_H__


#include "ceres/ceres.h" // Ceres-Solver


namespace TrackSys {


class VirtualPhyMuFit : protected TrFitPar, public ceres::CostFunction {
    public :
        VirtualPhyMuFit(const TrFitPar& fitPar, const PhySt& part) : 
            TrFitPar(fitPar), part_(part), opt_loc_(sw_mscat_), opt_tsft_(nmes_TOFt_>LMTN_TOF_T),
            DIMG_((PhyJb::DIMG+1) + (opt_tsft_?1:0)), 
            numOfRes_(0), numOfParGlb_(0), numOfParLoc_(0)
            { if (check_hits()) setvar(nseq_, nseg_); }
    
    public :
        virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;

    protected :
        inline void setvar(const Short_t nseq = 0, const Short_t nseg = 0) {
            if (nseq <= 0 || nseg < 0) return;
            Short_t nloc = (opt_loc_ ? nseg * PhyJb::DIML : 0);
            numOfRes_    = nseq + nloc;
            numOfParGlb_ = DIMG_;
            numOfParLoc_ = nloc;

            set_num_residuals(numOfRes_);
            mutable_parameter_block_sizes()->clear(); 
            mutable_parameter_block_sizes()->push_back(numOfParGlb_);
            if (opt_loc_) mutable_parameter_block_sizes()->push_back(numOfParLoc_);
        }
    
    protected :
        const Bool_t  opt_loc_;
        const Bool_t  opt_tsft_;
        const Short_t DIMG_;
        const PhySt   part_;
        
        Short_t numOfRes_;
        Short_t numOfParGlb_;
        Short_t numOfParLoc_;
        
    private :
        static constexpr Short_t parIDeta  = 4;
        static constexpr Short_t parIDibta = 5;
        static constexpr Short_t parIDtsft = 6;
};


class PhyMuFit : public TrFitPar {
    public :
        PhyMuFit& operator=(const PhyMuFit& rhs);
        PhyMuFit(const PhyMuFit& muFit) { *this = muFit; }
        
        PhyMuFit(const TrFitPar& fitPar);
        ~PhyMuFit() { PhyMuFit::clear(); }
        
    public :
        inline const Bool_t&   status() const { return succ_; }
        inline const PhySt&    part() const { return part_; }
        inline const Double_t& tsft() const { return tsft_; }
        
        inline const std::vector<PhyArg>& args() const { return args_; }
        
        inline const Short_t&  ndof(Int_t it)    const { return ndof_.at(it); }
        inline const Double_t& nchi(Int_t it)    const { return nchi_.at(it); }
        inline const Double_t& quality(Int_t it) const { return quality_.at(it); }

        PhySt interpolate_to_z(Double_t zcoo = 0) const;
        MatFld get_mat(Double_t zbd1 = 0, Double_t zbd2 = 0) const;

        // Check mass with electron-like or high beta events
        inline Bool_t is_like_el(Double_t fact = Numc::FIFTY<>) const { return (succ_ ? false : (Numc::Compare(part_.mu(), fact * EL_MU) <= 0)); }

        // Expert : do this after check status
        inline const TrFitPar& get() const { return dynamic_cast<const TrFitPar&>(*this); }           // fitPar for track fitting
        inline       PhyTrFit  fit() const { return PhyTrFit(dynamic_cast<const TrFitPar&>(*this)); } // physics track fitting
        
        inline const Sys::HrsStopwatch& timer() const { return timer_; }

    protected :
        void clear();

        Bool_t simpleScan(const TrFitPar& fitPar);
        Bool_t physicalFit();
        Bool_t evolve();

    protected :
        Bool_t              succ_;
        PhySt               part_;
        Double_t            tsft_; // time shift [cm]
        std::vector<PhyArg> args_; 

        // [0] (cx + mstau)
        // [1] (cy + msrho)
        // [2] (TRKq + TOFt + TOFq + RICHib + TRDel)
        std::array<Short_t,  3> ndof_;
        std::array<Double_t, 3> nchi_;
        std::array<Double_t, 3> quality_;

        Short_t  ndof_tt_; // (total)
        Double_t nchi_tt_; // (total)
    
    protected :
        std::vector<PhySt> stts_;
    
    private :
        static constexpr Short_t parIDeta  = 4;
        static constexpr Short_t parIDibta = 5;
        static constexpr Short_t parIDtsft = 6;
    
    protected :
        Sys::HrsStopwatch timer_;
};


} // namespace TrackSys


#endif // __TRACKLibs_PhyMuFit_H__
