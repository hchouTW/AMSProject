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
        inline const Double_t& tsft() const { return tsft_; }
        
        inline const std::vector<PhyArg>& args() const { return args_; }
        
        inline const Short_t&  ndof(Int_t it)    const { return ndof_.at(it); }
        inline const Double_t& nchi(Int_t it)    const { return nchi_.at(it); }
        inline const Double_t& quality(Int_t it) const { return quality_.at(it); }
        
        inline const Sys::HrsStopwatch& timer() const { return timer_; }

    protected :
        void clear();
        MuScanObj scan(const TrFitPar& fitPar, Double_t mass, Bool_t selopt = false);

    protected :
        Bool_t              succ_;
        PhySt               part_;
        Double_t            tsft_; // time shift [cm]
        std::vector<PhyArg> args_; 
        
        std::array<Short_t,  2> ndof_;
        std::array<Double_t, 2> nchi_;
        std::array<Double_t, 2> quality_;

    protected :
        static constexpr Short_t  LMT_SCAN = 3;
        static constexpr Double_t CONVG_DM = 0.005; // |preM - aftM|
        static constexpr Double_t CONVG_RM = 0.030; // |preM - aftM| / (preM + aftM)

        static constexpr Double_t LMT_MASS = 0.000510999; // electron mass
        static constexpr Double_t LMT_QLTR = 4.0;
        static constexpr Double_t LMT_QLTB = 3.7;
        static constexpr Double_t LMT_QLT  = 3.5;

        static const std::vector<std::vector<Double_t>> LIST_MASS_Q;
    
    protected :
        Sys::HrsStopwatch timer_;
};


// List of Particle
const std::vector<std::vector<Double_t>> SimpleMuScan::LIST_MASS_Q({
    { PartInfo(PartType::Photon).mass() }, // Q0
    { PartInfo(PartType::PionPlus).mass(), PartInfo(PartType::KaonPlus).mass(), PartInfo(PartType::Proton).mass(), PartInfo(PartType::Deuterium).mass() }, // Q1
    { PartInfo(PartType::Helium3).mass(), PartInfo(PartType::Helium4).mass() }, // Q2
    { PartInfo(PartType::Lithium6).mass(), PartInfo(PartType::Lithium7).mass() }, // Q3
    { PartInfo(PartType::Beryllium7).mass(), PartInfo(PartType::Beryllium9).mass(), PartInfo(PartType::Beryllium10).mass() }, // Q4
    { PartInfo(PartType::Boron10).mass(), PartInfo(PartType::Boron11).mass() }, // Q5
    { PartInfo(PartType::Carbon12).mass(), PartInfo(PartType::Carbon13).mass(), PartInfo(PartType::Carbon14).mass() }, // Q6
    { PartInfo(PartType::Nitrogen13).mass(), PartInfo(PartType::Nitrogen14).mass(), PartInfo(PartType::Nitrogen15).mass() }, // Q7
    { PartInfo(PartType::Oxygen16).mass(), PartInfo(PartType::Oxygen17).mass(), PartInfo(PartType::Oxygen18).mass() } // Q8
});


class SimpleMuScan::MuScanObj {
    public :
        MuScanObj(Short_t chrg = 0, Double_t mass = 0, Double_t qltr = 0, Double_t qltb = 0, Double_t qlt = 0) : chrg_(chrg), mass_(mass), qltr_(qltr), qltb_(qltb), qlt_(qlt) {}
        ~MuScanObj() {}

        inline const Short_t&  chrg() const { return chrg_; }
        inline const Double_t& mass() const { return mass_; }
        inline const Double_t& qltr() const { return qltr_; }
        inline const Double_t& qltb() const { return qltb_; }
        inline const Double_t& qlt()  const { return qlt_; }

    private :
        Short_t  chrg_;
        Double_t mass_;
        Double_t qltr_;
        Double_t qltb_;
        Double_t qlt_;
};


} // namespace TrackSys


#endif // __TRACKLibs_SimpleMuScan_H__
