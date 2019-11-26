#ifndef __TRACKLibs_SimpleMuScan_H__
#define __TRACKLibs_SimpleMuScan_H__


#include "ceres/ceres.h" // Ceres-Solver


namespace TrackSys {


class SimpleMuScan {
    protected :
        class MuScanObj;

    public :
        SimpleMuScan& operator=(const SimpleMuScan& rhs);
        SimpleMuScan(const SimpleMuScan& muScan) { *this = muScan; }
        
        SimpleMuScan(const TrFitPar& fitPar);
        ~SimpleMuScan() { SimpleMuScan::clear(); }
        
    public :
        inline const Bool_t&   status() const { return succ_; }
        inline const PhySt&    part() const { return part_; }
        inline const Double_t& ibta() const { return ibta_; }
        inline const Double_t& tsft() const { return tsft_; }
        
        inline const std::vector<PhyArg>& args() const { return args_; }
        
        inline const Short_t&  ndof(Int_t it)    const { return ndof_.at(it); }
        inline const Double_t& nchi(Int_t it)    const { return nchi_.at(it); }
        inline const Double_t& quality(Int_t it) const { return quality_.at(it); }
        
        inline const Short_t&  ndof_all()    const { return ndof_all_; }
        inline const Double_t& nchi_all()    const { return nchi_all_; }
        inline const Double_t& quality_all() const { return quality_all_; }
        
        inline const Sys::HrsStopwatch& timer() const { return timer_; }

    protected :
        void clear();
        MuScanObj scan(const TrFitPar& fitPar);

    protected :
        Bool_t              succ_;
        PhySt               part_;
        Double_t            ibta_; // beta from fit (may be not physical 1/beta)
        Double_t            tsft_; // time shift [cm]
        std::vector<PhyArg> args_; 
        
        std::array<Short_t,  2> ndof_;
        std::array<Double_t, 2> nchi_;
        std::array<Double_t, 2> quality_;
        
        Short_t  ndof_all_;
        Double_t nchi_all_;
        Double_t quality_all_;

    protected :
        static constexpr Double_t EL_MASS = 0.000510999;  // electron mass
        static const std::vector<std::vector<Double_t>> LIST_MASS_Q;
    
    protected :
        Sys::HrsStopwatch timer_;

    private :
        // Limit of 1/beta
        static constexpr Double_t LMTL_IBTA_APPROX_LIGHT = 1.00000000000001;
};


// List of Particle
const std::vector<std::vector<Double_t>> SimpleMuScan::LIST_MASS_Q({
    { PartInfo(PartType::Photon).mass() }, // Q0
    { PartInfo(PartType::Tritium).mass(), PartInfo(PartType::Deuterium).mass(), PartInfo(PartType::Proton).mass(), PartInfo(PartType::KaonPlus).mass(), PartInfo(PartType::PionPlus).mass() }, // Q1
    { PartInfo(PartType::Helium4).mass(), PartInfo(PartType::Helium3).mass() }, // Q2
    { PartInfo(PartType::Lithium7).mass(), PartInfo(PartType::Lithium6).mass() }, // Q3
    { PartInfo(PartType::Beryllium10).mass(), PartInfo(PartType::Beryllium9).mass(), PartInfo(PartType::Beryllium7).mass() }, // Q4
    { PartInfo(PartType::Boron11).mass(), PartInfo(PartType::Boron10).mass() }, // Q5
    { PartInfo(PartType::Carbon14).mass(), PartInfo(PartType::Carbon13).mass(), PartInfo(PartType::Carbon12).mass() }, // Q6
    { PartInfo(PartType::Nitrogen15).mass(), PartInfo(PartType::Nitrogen14).mass(), PartInfo(PartType::Nitrogen13).mass() }, // Q7
    { PartInfo(PartType::Oxygen18).mass(), PartInfo(PartType::Oxygen17).mass(), PartInfo(PartType::Oxygen16).mass() } // Q8
});


class SimpleMuScan::MuScanObj {
    public :
        MuScanObj(Short_t chrg = 0, Double_t mass = 0, Double_t qltr = 0, Double_t qltb = 0, Double_t mom = 0, Double_t ibta = 0, Double_t sqrm = 0, Double_t paribta = 0) : chrg_(chrg), mass_(mass), qltr_(qltr), qltb_(qltb), mom_(mom), ibta_(ibta), sqrm_(sqrm), paribta_(paribta) {}
        ~MuScanObj() {}

        inline void reset(Short_t chrg, Double_t mass) { if (chrg != 0 && Numc::Compare(mass) > 0) { chrg_ = chrg; mass_ = mass; qltr_ = 0; qltb_ = 0; mom_ = 0; ibta_ = 0; sqrm_ = 0; paribta_ = 0; } }
        
        inline const Short_t&  chrg() const { return chrg_; }
        inline const Double_t& mass() const { return mass_; }
        inline const Double_t& qltr() const { return qltr_; }
        inline const Double_t& qltb() const { return qltb_; }

        inline const Double_t& mom()  const { return mom_; }
        inline const Double_t& ibta() const { return ibta_; }
        
        inline const Double_t& sqrm() const { return sqrm_; }
        
        inline const Double_t& paribta() const { return paribta_; }

    private :
        Short_t  chrg_;
        Double_t mass_;
        Double_t qltr_;
        Double_t qltb_;

        Double_t mom_;
        Double_t ibta_;
        
        Double_t sqrm_; // sqrt_mass := mom^2 * (ibta^2 - 1) (may be not physical mass)
        
        Double_t paribta_; // beta from fitting
};


} // namespace TrackSys


#endif // __TRACKLibs_SimpleMuScan_H__
