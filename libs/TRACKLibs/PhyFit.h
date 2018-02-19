#ifndef __TRACKLibs_PhyFit_H__
#define __TRACKLibs_PhyFit_H__

// Ceres Solver
#ifdef __CeresSolver__
#include "ceres/ceres.h"
#endif

namespace TrackSys {


class TrFitPar {
    public :
        enum class Orientation {
            kDownward = 0, kUpward = 1
        };

    public :
        TrFitPar(const PartType& type = PartType::Proton, const Orientation& ortt = Orientation::kDownward, Bool_t sw_mscat = PhyArg::OptMscat(), Bool_t sw_eloss = PhyArg::OptEloss());
        ~TrFitPar() { clear(); }
   
        inline void addHit(HitSt& hit) { hits_.push_back(hit); is_check_ = false; }
        inline void delHit(Int_t it) { if (it >= 0 && it < static_cast<Int_t>(hits_.size())) { hits_.erase(hits_.begin()+it); is_check_ = false; } }
        
        inline Short_t numOfHit() const { return hits_.size(); }
        inline Short_t numOfSeq() const { return (nhtx_ + nhty_); }

    protected :
        Bool_t checkHit();
        void clear();

    protected :
        Bool_t      sw_mscat_;
        Bool_t      sw_eloss_;
        PartType    type_;
        Orientation ortt_;

        std::vector<HitSt> hits_;
        Short_t            nhtx_;
        Short_t            nhty_;

    private :
        Bool_t is_check_;

    protected :
        // Number of Hit Requirement
        static constexpr Short_t LMTL_NHIT_X = 3;
        static constexpr Short_t LMTL_NHIT_Y = 4;
};


class SimpleTrFit : public TrFitPar {
    public :
        SimpleTrFit(TrFitPar& fitPar); 
        ~SimpleTrFit() { SimpleTrFit::clear(); TrFitPar::clear(); }
        
    public :
        inline const Bool_t& status() { return succ_; }
        inline const PhySt& part() const { return part_; }
        
        inline const Int_t& ndof() const { return ndof_; }
        inline const Double_t& nchi() const { return nchi_; }
        inline Double_t nchix() const { return succ_?(chix_/ndfx_):0; }
        inline Double_t nchiy() const { return succ_?(chiy_/ndfy_):0; }

    protected :
        void clear();

        Bool_t analyticalFit();
        Bool_t simpleFit();

    protected :
        Bool_t   succ_;
        PhySt    part_;
        
        Double_t ndfx_;
        Double_t ndfy_;
        Double_t chix_;
        Double_t chiy_;

        Int_t    ndof_;
        Double_t nchi_;

    protected :
        // Minimization (Levenberg-Marquardt Method)
        static constexpr Int_t    LMTL_ITER = 3;
        static constexpr Int_t    LMTU_ITER = 25;
        
        static constexpr Double_t LAMBDA0 = 1.0e-2;
        static constexpr Double_t LAMBDA_DN_FAC = 5.0;
        static constexpr Double_t LAMBDA_UP_FAC = 7.0;
        static constexpr Double_t LMTL_LAMBDA = 1.0e-4;
        static constexpr Double_t LMTU_LAMBDA = 1.0e+3;
        
        static constexpr Double_t CONVG_EPSILON   = 3.0e-3;
        static constexpr Double_t CONVG_TOLERANCE = 7.0e-3;
};


#ifdef __CeresSolver__
using EMtxXD = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
using EVecXD = Eigen::Matrix<double, Eigen::Dynamic, 1>;

class VirtualPhyTrFit : public TrFitPar, public ceres::CostFunction {
    public :
        VirtualPhyTrFit(TrFitPar& fitPar, PhySt& part) : TrFitPar(fitPar), part_(part) { checkHit(); set_num_residuals(numOfSeq()); mutable_parameter_block_sizes()->push_back(5); }
        ~VirtualPhyTrFit() { VirtualPhyTrFit::clear(); TrFitPar::clear(); }
    
    public :
        virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;
    
    protected :
        void clear() { part_.reset(type_); part_.arg().reset(sw_mscat_, sw_eloss_); }
    
    protected :
        PhySt part_;
};


class PhyTrFit : public TrFitPar {
    public :
        PhyTrFit(TrFitPar& fitPar);
        ~PhyTrFit() { PhyTrFit::clear(); TrFitPar::clear(); }
        
    public :
        inline const Bool_t& status() { return succ_; }
        inline const PhySt& part() const { return part_; }
        
        inline const Int_t& ndof() const { return ndof_; }
        inline const Double_t& nchi() const { return nchi_; }
        inline const Double_t& nchix() const { return nchi_; }
        inline const Double_t& nchiy() const { return nchi_; }

    protected :
        void clear();

        Bool_t simpleFit();
        Bool_t physicalFit();
    
    protected :
        Bool_t   succ_;
        PhySt    part_;
        
        Int_t    ndof_;
        Double_t nchi_;
};
#endif


} // namespace TrackSys


#endif // __TRACKLibs_PhyFit_H__
