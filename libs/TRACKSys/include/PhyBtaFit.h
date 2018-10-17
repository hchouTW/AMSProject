#ifndef __TRACKLibs_PhyBtaFit_H__
#define __TRACKLibs_PhyBtaFit_H__


#include "ceres/ceres.h" // Ceres-Solver


namespace TrackSys {


class VirtualPhyBtaFit : protected TrFitPar, public ceres::CostFunction {
    public :
        VirtualPhyBtaFit(const TrFitPar& fitPar, const PhySt& part) : 
            TrFitPar(fitPar), part_(part), opt_tsft_(nmes_TOFt_>=LMTN_TOF_T),
            numOfRes_(0), numOfPar_(0)
            { if (check_hits()) setvar(nseq_); }
    
    public :
        virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;

    protected :
        inline void setvar(const Short_t nseq = 0) {
            if (nseq <= 0) return;
            numOfRes_ = nseq;
            numOfPar_ = 1 + opt_tsft_;

            set_num_residuals(numOfRes_);
            mutable_parameter_block_sizes()->clear(); 
            mutable_parameter_block_sizes()->push_back(numOfPar_);
        }
    
    protected :
        const Bool_t opt_tsft_;
        const PhySt  part_;
        
        Short_t numOfRes_;
        Short_t numOfPar_;

    private :
        static constexpr Short_t parIDigb  = 0;
        static constexpr Short_t parIDtsft = 1;
};


class PhyBtaFit : public TrFitPar {
    public :
        TrFitPar BulidFitPar(const TrFitPar& fitPar); // only for bta fit
        PhySt BulidRefSt(const TrFitPar& fitPar, Double_t refz = 0.0); // only coord info
        PhySt BulidRefSt(const PhySt& refSt, Double_t refz = 0.0);

    public :
        PhyBtaFit& operator=(const PhyBtaFit& rhs);
        PhyBtaFit(const PhyBtaFit& trFit) { *this = trFit; }
        
        PhyBtaFit(const TrFitPar& fitPar);
        PhyBtaFit(const TrFitPar& fitPar, const PhySt& refSt);
        ~PhyBtaFit() { PhyBtaFit::clear(); }
        
    public :
        inline const Bool_t&   status() const { return succ_; }
        inline const PhySt&    part() const { return part_; }
        inline const Double_t& tsft() const { return tsft_; }

        inline const Short_t&  ndof() const { return ndof_; }
        inline const Double_t& nchi() const { return nchi_; }
        inline const Double_t& quality() const { return quality_; }
       
        PhySt interpolate_to_z(Double_t zcoo = 0) const;
        MatFld get_mat(Double_t zbd1 = 0, Double_t zbd2 = 0) const;

    protected :
        void   clear();
        Bool_t physicalFit();
        Bool_t evolve();

    protected :
        Bool_t   succ_;
        PhySt    part_;
        Double_t tsft_; // time shift [cm]

        Short_t  ndof_;
        Double_t nchi_;
        Double_t quality_;
    
    private :
        static constexpr Double_t LMT_IGB   = 1.0e-10;
        static constexpr Double_t CONV_IGB  = 1.0e-08;
        static constexpr Short_t  parIDigb  = 0;
        static constexpr Short_t  parIDtsft = 1;
};


} // namespace TrackSys


#endif // __TRACKLibs_PhyBtaFit_H__
