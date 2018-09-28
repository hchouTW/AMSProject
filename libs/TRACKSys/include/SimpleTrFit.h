#ifndef __TRACKLibs_SimpleTrFit_H__
#define __TRACKLibs_SimpleTrFit_H__


#include "ceres/ceres.h" // Ceres-Solver


namespace TrackSys {


class VirtualSimpleTrCooFit : public TrFitPar, public ceres::FirstOrderFunction {
    public :
        VirtualSimpleTrCooFit(const TrFitPar& fitPar, const PhySt& part) :
            TrFitPar(fitPar), part_(part), opt_loc_(sw_mscat_),
            numOfRes_(0), numOfPar_(0), numOfDof_(0)
            { if (check_hits()) setvar(onlyc_nseq_, PhyJb::DIMG); } 
    
    public :
        virtual bool Evaluate(const double* parameters, double* cost, double* gradient) const;
        virtual int NumParameters() const { return numOfPar_; }
        
    protected :
        inline void setvar(const Short_t nres = 0, const Short_t npar = 0) {
            if (nres <= 0 || npar <= 0) return;
            numOfRes_ = nres; 
            numOfPar_ = npar;
            numOfDof_ = (numOfRes_ - numOfPar_);
        }

    protected :
        const Bool_t opt_loc_;
        const PhySt  part_;
        
        Short_t numOfRes_;
        Short_t numOfPar_;
        Short_t numOfDof_;
};


class VirtualSimpleTrFit : public TrFitPar, public ceres::FirstOrderFunction {
    public :
        VirtualSimpleTrFit(const TrFitPar& fitPar, const PhySt& part) :
            TrFitPar(fitPar), part_(part), opt_loc_(sw_mscat_), opt_tsft_(nmes_TOFt_>=LMTN_TOF_T),
            numOfRes_(0), numOfPar_(0), numOfDof_(0)
            { if (check_hits()) setvar(nseq_, PhyJb::DIMG); } 
    
    public :
        virtual bool Evaluate(const double* parameters, double* cost, double* gradient) const;
        virtual int NumParameters() const { return numOfPar_; }
        
    protected :
        inline void setvar(const Short_t nres = 0, const Short_t npar = 0) {
            if (nres <= 0 || npar <= 0) return;
            numOfRes_ = nres; 
            numOfPar_ = npar;
            numOfDof_ = (numOfRes_ - numOfPar_);
        }

    protected :
        const Bool_t opt_loc_;
        const Bool_t opt_tsft_;
        const PhySt  part_;
        
        Short_t numOfRes_;
        Short_t numOfPar_;
        Short_t numOfDof_;
};


class VirtualSimpleTrLocFit : public TrFitPar, public ceres::CostFunction {
    public :
        VirtualSimpleTrLocFit(const TrFitPar& fitPar, const PhySt& part) :
            TrFitPar(fitPar), part_(part),
            numOfRes_(0), numOfPar_(0), numOfDof_(0)
            { if (check_hits()) setvar(onlyc_nseq_ + nseg_ * PhyJb::DIML, nseg_ * PhyJb::DIML); } 
    
    public :
        virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;
        
    protected :
        inline void setvar(const Short_t nres = 0, const Short_t npar = 0) {
            if (nres <= 0 || npar <= 0) return;
            numOfRes_ = nres; 
            numOfPar_ = npar;
            numOfDof_ = (numOfRes_ - numOfPar_);
            
            set_num_residuals(numOfRes_);
            mutable_parameter_block_sizes()->clear(); 
            mutable_parameter_block_sizes()->push_back(numOfPar_);
        }
    
    protected :
        const PhySt part_;

        Short_t numOfRes_;
        Short_t numOfPar_;
        Short_t numOfDof_;
};


class SimpleTrFit : public TrFitPar {
    public :
        SimpleTrFit& operator=(const SimpleTrFit& rhs);
        SimpleTrFit(const SimpleTrFit& trFit) { *this = trFit; }
        
        SimpleTrFit(const TrFitPar& fitPar, Bool_t withLocal = false); 
        ~SimpleTrFit() { SimpleTrFit::clear(); }
        
    public :
        inline const Bool_t&   status() const { return succ_; }
        inline const PhySt&    part() const { return part_; }
        
        inline const std::vector<PhyArg>& args() const { return args_; }

        inline const Short_t&  ndof(Int_t it)    const { return ndof_.at(it); }
        inline const Double_t& nchi(Int_t it)    const { return nchi_.at(it); }
        inline const Double_t& quality(Int_t it) const { return quality_.at(it); }
        
        inline const Short_t& ndof_cx() const { return ndof_cx_; }
        inline const Short_t& ndof_cy() const { return ndof_cy_; }
        inline const Short_t& ndof_ib() const { return ndof_ib_; }
        
        inline const Double_t& nchi_cx() const { return nchi_cx_; }
        inline const Double_t& nchi_cy() const { return nchi_cy_; }
        inline const Double_t& nchi_ib() const { return nchi_ib_; }

        PhySt interpolate_to_z(Double_t zcoo = 0) const;

    protected :
        void clear();

        Bool_t analyticalFit();
        Bool_t simpleFit();
        Bool_t advancedSimpleCooFit();
        Bool_t advancedSimpleFit();
        Bool_t localSimpleFit();
        Bool_t evolve();

    protected :
        Bool_t              succ_;
        PhySt               part_;
        std::vector<PhyArg> args_; 
        
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
    
    protected :
        std::vector<PhySt> stts_;

    protected :
        static constexpr Short_t DIMG = PhyJb::DIMG;
        static constexpr Short_t DIML = PhyJb::DIML;
        
        // Minimization (Analytical Method)
        static constexpr Double_t LMTU_ETA = 3.0;
        
        // Minimization (Levenberg-Marquardt Method)
        static constexpr Short_t LMTL_ITER = 3;
        static constexpr Short_t LMTU_ITER = 20;
        
        static constexpr Double_t LAMBDA0 = 1.0e-2;
        static constexpr Double_t LAMBDA_DN_FAC = 5.0;
        static constexpr Double_t LAMBDA_UP_FAC = 7.0;
        static constexpr Double_t LMTL_LAMBDA = 1.0e-4;
        static constexpr Double_t LMTU_LAMBDA = 1.0e+3;
        
        static constexpr Double_t CONVG_EPSILON   = 5.0e-2;
        static constexpr Double_t CONVG_TOLERANCE = 5.0e-2;
};


} // namespace TrackSys


#endif // __TRACKLibs_SimpleTrFit_H__
