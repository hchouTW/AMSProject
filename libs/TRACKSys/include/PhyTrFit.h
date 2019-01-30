#ifndef __TRACKLibs_PhyTrFit_H__
#define __TRACKLibs_PhyTrFit_H__


#include "ceres/ceres.h" // Ceres-Solver


namespace TrackSys {


class VirtualPhyTrFit : protected TrFitPar, public ceres::CostFunction {
    public :
        VirtualPhyTrFit(const TrFitPar& fitPar, const PhySt& part) : 
            TrFitPar(fitPar), part_(part), opt_loc_(sw_mscat_), opt_tsft_(nmes_TOFt_>LMTN_TOF_T),
            DIMG_(PhyJb::DIMG + (opt_tsft_?1:0)), 
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
        static constexpr Short_t parIDtsft = 5;
};


class PhyTrFit : public TrFitPar {
    public :
        class LocScat;

    public :
        PhyTrFit& operator=(const PhyTrFit& rhs);
        PhyTrFit(const PhyTrFit& trFit) { *this = trFit; }
        
        PhyTrFit(const TrFitPar& fitPar);
        ~PhyTrFit() { PhyTrFit::clear(); }
        
    public :
        inline const Bool_t&   status() const { return succ_; }
        inline const PhySt&    part() const { return part_; }
        inline const Double_t& tsft() const { return tsft_; }
        
        inline const std::vector<PhyArg>&  args()  const { return args_; }
        inline const std::vector<PhySt>&   stts()  const { return stts_; }
        inline const std::vector<LocScat>& lscat() const { return lscat_; }
        
        inline const Short_t&  ndof(Int_t it)    const { return ndof_.at(it); }
        inline const Double_t& nchi(Int_t it)    const { return nchi_.at(it); }
        inline const Double_t& quality(Int_t it) const { return quality_.at(it); }
        
        inline const Short_t&  ndof_all() const { return ndof_all_; }
        inline const Double_t& nchi_all() const { return nchi_all_; }
        inline const Double_t& quality_all() const { return quality_all_; }

        inline const Short_t& ndof_cx() const { return ndof_cx_; }
        inline const Short_t& ndof_cy() const { return ndof_cy_; }
        inline const Short_t& ndof_ib() const { return ndof_ib_; }
        
        inline const Double_t& nchi_cx() const { return nchi_cx_; }
        inline const Double_t& nchi_cy() const { return nchi_cy_; }
        inline const Double_t& nchi_ib() const { return nchi_ib_; }
        
        PhySt interpolate_to_z(Double_t zcoo = 0) const;
        MatFld get_mat(Double_t zbd1 = 0, Double_t zbd2 = 0) const;

    protected :
        void clear();

        Bool_t simpleFit();
        Bool_t physicalFit();
        Bool_t evolve();
    
    protected :
        Bool_t              succ_;
        PhySt               part_;
        Double_t            tsft_; // time shift [cm]
        std::vector<PhyArg> args_; 

        std::array<Short_t,  2> ndof_;
        std::array<Double_t, 2> nchi_;
        std::array<Double_t, 2> quality_;

        Short_t  ndof_all_;
        Double_t nchi_all_;
        Double_t quality_all_;

        Short_t ndof_tt_; // (total)
        Short_t ndof_cx_; // (cx + mstau)
        Short_t ndof_cy_; // (cy + msrho)
        Short_t ndof_ib_; // (TRKq + TOFt + TOFq + RICHib + TRDel)

        Double_t nchi_tt_;
        Double_t nchi_cx_;
        Double_t nchi_cy_;
        Double_t nchi_ib_;
        
    protected :
        std::vector<PhySt>   stts_;
        std::vector<LocScat> lscat_;

    private :
        static constexpr Short_t parIDeta  = 4;
        static constexpr Short_t parIDtsft = 5;
};
        

class PhyTrFit::LocScat {
    public :
        LocScat(Double_t cx = 0, Double_t cy = 0, Double_t cz = 0, Bool_t scx = false, Bool_t scy = false, Double_t chix = 0, Double_t chiy = 0, Double_t chit = 0, Double_t chir = 0) : coo_({cx, cy, cz}), side_({scx, scy}) {
            chic_ = std::array<Double_t, 2>({ (scx?chix:Numc::ZERO<>), (scy?chiy:Numc::ZERO<>) });
            chis_ = std::array<Double_t, 2>({ (scx?chit:Numc::ZERO<>), (scy?chir:Numc::ZERO<>) });
        }
        ~LocScat() {}

        inline const Bool_t&  scx() const { return side_[0]; }
        inline const Bool_t&  scy() const { return side_[1]; }

        inline const Double_t& cx() const { return coo_[0]; }
        inline const Double_t& cy() const { return coo_[1]; }
        inline const Double_t& cz() const { return coo_[2]; }

        inline const Double_t& chic(Int_t it) const { return chic_.at(it); }
        inline const Double_t& chis(Int_t it) const { return chis_.at(it); }

    private :
        std::array<Double_t, 3> coo_;
        std::array<Bool_t, 2>   side_;
        std::array<Double_t, 2> chic_;
        std::array<Double_t, 2> chis_;
};


} // namespace TrackSys


#endif // __TRACKLibs_PhyTrFit_H__
