#ifndef __TRACKLibs_PhyFit_H__
#define __TRACKLibs_PhyFit_H__


#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"

// Ceres Solver
#ifdef __CeresSolver__
#include "ceres/ceres.h"
#include "glog/logging.h"
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
   
        inline void addHit(HitSt& hit) { hits_.push_back(hit); }
        inline Short_t numOfHit() const { return hits_.size(); }
    
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
class VirtualPhyTrFit : public TrFitPar, public ceres::FirstOrderFunction {
    public :
        VirtualPhyTrFit(TrFitPar& fitPar);
        ~VirtualPhyTrFit() { VirtualPhyTrFit::clear(); TrFitPar::clear(); }
    
    public :
        virtual bool Evaluate(const double* parameters, double* cost, double* gradient) const;
        virtual int NumParameters() const { return 5; }
        
    public :
        inline const Bool_t& status() { return succ_; }
        inline const PhySt& part() const { return part_; }
        
        inline const Int_t& ndof() const { return ndof_; }
        inline const Double_t& nchi() const { return nchi_; }
    
    protected :
        void clear();
    
    protected :
        Bool_t   succ_;
        PhySt    part_;
        
        Int_t    ndof_;
        Double_t nchi_;
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
/*
    // Used to Minimizer
    protected :
        //double Chisq(const double* x);

        void FdF(const double *x, double &est, double *grad);
    private :
        IBaseFunctionMultiDim * Clone() const { return nullptr; }
        unsigned int NDim() const { return 5; }
        double DoEval(const double* x) const { return 0; }
        double DoDerivative(const double * x, unsigned int icoord ) const { return 0; }

    protected :
        PhySt               mnz_part_;
        std::vector<PhyArg> mnz_args_;
*/
};
#endif



} // namespace TrackSys


#endif // __TRACKLibs_PhyFit_H__
