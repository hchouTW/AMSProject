#ifndef __TRACKLibs_PartInfo_H__
#define __TRACKLibs_PartInfo_H__


namespace TrackSys {


enum class PartType {
    Self,
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
        PartInfo(const PartType& type = PartType::Proton);
        ~PartInfo() {}

        void print() const;

        inline const PartType&     type() const { return type_; }
        inline const std::string&  name() const { return name_; }
        inline const Short_t&      chrg() const { return chrg_; } // [1]
        inline const Double_t&     mass() const { return mass_; } // [GeV]

        inline const Bool_t& is_chrgless() const { return is_chrgless_; }
        inline const Bool_t& is_massless() const { return is_massless_; }
        
        inline const Double_t& mass_to_chrg() const { return mass_to_chrg_; } // [GeV]
        inline const Double_t& chrg_to_mass() const { return chrg_to_mass_; } // [1/GeV]
       
        // testcode
        //inline Double_t umass() const { return (mass_ / ATOMIC_MASS); } // [1]
        //inline Double_t umass_to_chrg() const { return (mass_to_chrg_ / ATOMIC_MASS); } // [1]
        //inline Double_t chrg_to_umass() const { return (chrg_to_mass_ * ATOMIC_MASS); } // [1]

    private :
        PartType    type_;
        std::string name_;
        Short_t     chrg_;
        Double_t    mass_; // [GeV]
        
        Bool_t      is_chrgless_;
        Bool_t      is_massless_;
        
        Double_t    mass_to_chrg_;
        Double_t    chrg_to_mass_;

    public :
        static inline void SetSelf(const PartType& type = PartType::Proton) { PartInfo info(type); SelfMass_ = info.mass(); SelfChrg_ = info.chrg(); SelfName_ = info.name(); }
        static inline void SetSelf(Double_t mass, Short_t chrg = 1, const std::string& name = "Self") { SelfMass_ = mass; SelfChrg_ = chrg; SelfName_ = name; }
        static inline const std::string& SelfName() { return SelfName_; }
        static inline const Short_t&     SelfChrg() { return SelfChrg_; }
        static inline const Double_t&    SelfMass() { return SelfMass_; }

    private :
        static std::string SelfName_;
        static Short_t     SelfChrg_;
        static Double_t    SelfMass_;
    
    public :
        // Atomic mass unit u = 0.931494095 GeV/c^2
        static constexpr Double_t ATOMIC_MASS = 0.931494095;
};

// List of Particle
static const PartInfo PIElectron(PartType::Electron);
static const PartInfo PIProton(PartType::Proton);
static const PartInfo PIDeuterium(PartType::Deuterium);
static const PartInfo PIHelium3(PartType::Helium3);
static const PartInfo PIHelium4(PartType::Helium4);

// Self Particle
std::string PartInfo::SelfName_ = "Self";
Short_t     PartInfo::SelfChrg_ = PIProton.chrg();
Double_t    PartInfo::SelfMass_ = PIProton.mass();

} // namespace TrackSys


#endif // __TRACKLibs_PartInfo_H__
