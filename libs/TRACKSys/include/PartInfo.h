#ifndef __TRACKLibs_PartInfo_H__
#define __TRACKLibs_PartInfo_H__


namespace TrackSys {


enum class PartType {
    Fixed,
    Q0, Q1, Q2, Q3, Q4, Q5, Q6, Q7, Q8,
    Photon, 
    Electron, Positron, 
    Muon, 
    PionPlus, PionMinus,
    KaonPlus, KaonMinus,
    Proton, AntiProton,
    Deuterium, AntiDeuterium,
    Tritium, AntiTritium,
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
        PartInfo(const PartInfo& info) { *this = info; }
        PartInfo(const PartType& type = PartType::Proton) : type_(PartType::Fixed), name_(""), chrg_(0), mass_(0), mu_(0), is_chrgless_(true), is_massless_(true), chrg_to_mass_(0), chrg_to_atomic_mass_(0) { reset(type); }
        PartInfo(Short_t chrg, Double_t mass) { reset(chrg, mass); }
        ~PartInfo() {}
        
        void reset(const PartType& type);
        inline void reset(Short_t chrg, Double_t mass) { reset(PartType::Fixed, "", chrg, mass); }
        inline void reset(Double_t mu) { reset(PartType::Fixed, "", chrg_, ((Numc::Compare(mu)<=0) ? Numc::ZERO<> : mu * ATOMIC_MASS) ); }

        void print() const;

        inline Bool_t is_std() const { return (type_ != PartType::Fixed); }

        inline const PartType&     type() const { return type_; }
        inline const std::string&  name() const { return name_; }
        inline const Short_t&      chrg() const { return chrg_; } // [1]
        inline const Double_t&     mass() const { return mass_; } // [GeV]
        inline const Double_t&     mu()   const { return mu_; }   // [1]

        inline const Bool_t& is_chrgless() const { return is_chrgless_; }
        inline const Bool_t& is_massless() const { return is_massless_; }
        
        inline const Double_t& chrg_to_mass() const { return chrg_to_mass_; } // [1/GeV]
        inline const Double_t& chrg_to_atomic_mass() const { return chrg_to_atomic_mass_; } // [1/GeV]

    protected :
        void reset(const PartType& type, const std::string& name, Short_t chrg, Double_t mass);

    private :
        PartType    type_;
        std::string name_;
        Short_t     chrg_;
        Double_t    mass_; // [GeV]
        Double_t    mu_;   // Mass / AtomicMass [1]

        Bool_t      is_chrgless_;
        Bool_t      is_massless_;
        
        Double_t    chrg_to_mass_;
        Double_t    chrg_to_atomic_mass_;

    public :
        static inline void SetDefault(const PartType& type) { PartInfo info(type); DefaultChrg_ = info.chrg(); DefaultMass_ = info.mass(); }
        static inline void SetDefault(Short_t chrg, Double_t mass) { DefaultChrg_ = chrg; DefaultMass_ = (Numc::Compare(mass)<=0 ? Numc::ZERO<> : mass); }
        static inline void SetDefaultChrg(Short_t chrg) { DefaultChrg_ = chrg; }
        static inline void SetDefaultMass(Double_t mass) { DefaultMass_ = (Numc::Compare(mass)<=0 ? Numc::ZERO<> : mass); }
        static inline void SetDefaultMu(Double_t mu) { DefaultMass_ = (Numc::Compare(mu)<=0 ? Numc::ZERO<> : (mu * ATOMIC_MASS)); }
        static inline const Short_t&  DefaultChrg() { return DefaultChrg_; }
        static inline const Double_t& DefaultMass() { return DefaultMass_; }

    private :
        static Short_t  DefaultChrg_;
        static Double_t DefaultMass_;
    
    public :
        // Atomic mass unit u = 0.931494095 GeV/c^2
        static constexpr Double_t ATOMIC_MASS = 0.931494095;
};


} // namespace TrackSys


#endif // __TRACKLibs_PartInfo_H__
