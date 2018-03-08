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
        inline Short_t numOfSeq() const { return 2 * (nhtx_ + nhty_); }

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


class SimpleTrFit : protected TrFitPar {
    public :
        SimpleTrFit(TrFitPar& fitPar); 
        ~SimpleTrFit() { SimpleTrFit::clear(); TrFitPar::clear(); }
        
    public :
        inline const Bool_t& status() const { return succ_; }
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
class VirtualPhyTrFit : protected TrFitPar, public ceres::CostFunction {
    public :
        VirtualPhyTrFit(TrFitPar& fitPar, PhySt& part) : TrFitPar(fitPar), part_(part), numOfRes_(0), numOfPar_(0) { if (checkHit()) setvar(numOfSeq()+(numOfHit()-1)*PhyJb::DIM_L, PhyJb::DIM_G+(numOfHit()-1)*PhyJb::DIM_L); }
        ~VirtualPhyTrFit() { VirtualPhyTrFit::clear(); TrFitPar::clear(); }
    
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
        ~PhyTrFit() { PhyTrFit::clear(); TrFitPar::clear(); }
        
    public :
        inline const Bool_t& status() const { return succ_; }
        inline const PhySt& part() const { return part_; }

        inline const Int_t nargs() const { return args_.size(); }
        inline const std::vector<PhyArg>& args() const { return args_; }
        inline const PhyArg* args(Int_t it) const { return ((succ_ && it>=0 && it<args_.size()) ? &args_.at(it) : nullptr); }

        inline const Int_t nstts() const { return stts_.size(); }
        inline const std::vector<PhySt>& stts() const { return stts_; }
        inline const PhySt* stts(Int_t lay) const { auto&& stt = map_stts_.find(lay); return ((succ_ && stt!=map_stts_.end()) ? stt->second : nullptr); }
        
        inline const Int_t nhits() const { return hits_.size(); }
        inline const std::vector<HitSt>& hits() const { return hits_; }
        inline const HitSt* hits(Int_t lay) const { auto&& hit = map_hits_.find(lay); return ((succ_ && hit!=map_hits_.end()) ? hit->second : nullptr); }

        inline const Int_t& ndof() const { return ndof_; }
        inline const Double_t& nchi() const { return nchi_; }
        inline Double_t nchix() const { return succ_?((chix_+chit_)/ndfx_):0; }
        inline Double_t nchiy() const { return succ_?((chiy_+chir_)/ndfy_):0; }
        //inline Double_t nchit() const { return succ_?(chit_/ndfx_):0; }
        //inline Double_t nchir() const { return succ_?(chir_/ndfy_):0; }

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
        
        Double_t ndfx_;
        Double_t ndfy_;
        Double_t chix_;
        Double_t chiy_;
        Double_t chit_;
        Double_t chir_;
        
        Int_t    ndof_;
        Double_t nchi_;

    protected :
        std::map<Int_t, HitSt*> map_hits_;
        std::map<Int_t, PhySt*> map_stts_;
};
#endif


} // namespace TrackSys


#endif // __TRACKLibs_PhyFit_H__
