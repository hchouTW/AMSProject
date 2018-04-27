#ifndef __TRACKLibs_PartInfo_H__
#define __TRACKLibs_PartInfo_H__


namespace TrackSys {


enum class PartType {
    Fixed,
    Photon, 
    Electron, Positron, 
    Muon, 
    PionPlus, PionMinus,
    KaonPlus, KaonMinus,
    Proton, AntiProton,
    Deuterium, AntiDeuterium,
    Helium3, AntiHelium3, Helium4, AntiHelium4,
    Lithium6, Lithium7,
    Beryllium7, Beryllium9, Beryllium10,
    Boron10, Boron11,
    Carbon12, Carbon13, Carbon14,
    Nitrogen13, Nitrogen14, Nitrogen15,
    Oxygen16, Oxygen17, Oxygen18
};


class PartInfo {
    public :
        PartInfo(const PartType& type = PartType::Proton) : type_(PartType::Fixed), name_(""), chrg_(0), mass_(0), invu_(0), is_chrgless_(true), is_massless_(true), chrg_to_mass_(0), chrg_to_atomic_mass_(0) { reset(type); }
        PartInfo(Short_t chrg, Double_t mass) { reset(chrg, mass); }
        ~PartInfo() {}
        
        inline void reset(const PartType& type);
        inline void reset(Short_t chrg, Double_t mass) { reset(PartType::Fixed, "", chrg, mass); }
        inline void reset(Double_t invu) { reset(PartType::Fixed, "", chrg_, ((Numc::Compare(invu)<=0) ? Numc::ZERO<> : ATOMIC_MASS/invu) ); }

        void print() const;

        inline Bool_t is_std() const { return (type_ != PartType::Fixed); }

        inline const PartType&     type() const { return type_; }
        inline const std::string&  name() const { return name_; }
        inline const Short_t&      chrg() const { return chrg_; } // [1]
        inline const Double_t&     mass() const { return mass_; } // [GeV]
        inline const Double_t&     invu() const { return invu_; } // [1]

        inline const Bool_t& is_chrgless() const { return is_chrgless_; }
        inline const Bool_t& is_massless() const { return is_massless_; }
        
        inline const Double_t& chrg_to_mass() const { return chrg_to_mass_; } // [1/GeV]
        inline const Double_t& chrg_to_atomic_mass() const { return chrg_to_atomic_mass_; } // [1/GeV]

    protected :
        inline void reset(const PartType& type, const std::string& name, Short_t chrg, Double_t mass);

    private :
        PartType    type_;
        std::string name_;
        Short_t     chrg_;
        Double_t    mass_; // [GeV]
        Double_t    invu_; // [1]

        Bool_t      is_chrgless_;
        Bool_t      is_massless_;
        
        Double_t    chrg_to_mass_;
        Double_t    chrg_to_atomic_mass_;

    public :
        static inline void SetDefault(const PartType& type) { PartInfo info(type); DefaultChrg_ = info.chrg(); DefaultMass_ = info.mass(); }
        static inline void SetDefault(Short_t chrg, Double_t mass) { DefaultChrg_ = chrg; DefaultMass_ = (Numc::Compare(mass)<=0 ? Numc::ZERO<> : mass); }
        static inline void SetDefaultChrg(Short_t chrg) { DefaultChrg_ = chrg; }
        static inline void SetDefaultMass(Double_t mass) { DefaultMass_ = (Numc::Compare(mass)<=0 ? Numc::ZERO<> : mass); }
        static inline void SetDefaultInvu(Double_t invu) { DefaultMass_ = (Numc::Compare(invu)<=0 ? Numc::ZERO<> : (ATOMIC_MASS/invu)); }
        static inline const Short_t&  DefaultChrg() { return DefaultChrg_; }
        static inline const Double_t& DefaultMass() { return DefaultMass_; }

    private :
        static Short_t  DefaultChrg_;
        static Double_t DefaultMass_;
    
    public :
        // Atomic mass unit u = 0.931494095 GeV/c^2
        static constexpr Double_t ATOMIC_MASS = 0.931494095;
};


// Self Particle
Short_t  PartInfo::DefaultChrg_ = 0;
Double_t PartInfo::DefaultMass_ = 0;


// List of Particle
const std::vector<std::vector<Double_t>> PartListMassQ({
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

const Short_t PartListMaxQ = PartListMassQ.size()-1;


} // namespace TrackSys


#endif // __TRACKLibs_PartInfo_H__
