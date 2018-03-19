#ifndef __TRACKLibs_PartInfo_C__
#define __TRACKLibs_PartInfo_C__


namespace TrackSys {


PartInfo::PartInfo(const PartType& type) {
    // Atomic mass unit u = 0.931494095 GeV/c^2
    type_ = type;
    name_ = "";
    chrg_ = 0;
    mass_ = 0.;

    switch (type) {
        case PartType::None          : name_ = "None"         ; chrg_ =  0; mass_ =  0.000000000; break;
        case PartType::Photon        : name_ = "Photon"       ; chrg_ =  0; mass_ =  0.000000000; break;
        case PartType::Electron      : name_ = "Electron"     ; chrg_ = -1; mass_ =  0.000510999; break;
        case PartType::Positron      : name_ = "Positron"     ; chrg_ =  1; mass_ =  0.000510999; break;
        case PartType::Muon          : name_ = "Muon"         ; chrg_ = -1; mass_ =  0.105658367; break;
        case PartType::PionPlus      : name_ = "PionPlus"     ; chrg_ =  1; mass_ =  0.139570180; break;
        case PartType::PionMinus     : name_ = "PionMinus"    ; chrg_ = -1; mass_ =  0.139570180; break;
        case PartType::KaonPlus      : name_ = "KaonPlus"     ; chrg_ =  1; mass_ =  0.493677000; break;
        case PartType::KaonMinus     : name_ = "KaonMinus"    ; chrg_ = -1; mass_ =  0.493677000; break;
        case PartType::Proton        : name_ = "Proton"       ; chrg_ =  1; mass_ =  0.938272297; break;
        case PartType::AntiProton    : name_ = "AntiProton"   ; chrg_ = -1; mass_ =  0.938272297; break;
        case PartType::Deuterium     : name_ = "Deuterium"    ; chrg_ =  1; mass_ =  1.876123915; break;
        case PartType::AntiDeuterium : name_ = "AntiDeuterium"; chrg_ = -1; mass_ =  1.876123915; break;
        case PartType::Helium3       : name_ = "Helium3"      ; chrg_ =  2; mass_ =  2.809413500; break;
        case PartType::AntiHelium3   : name_ = "AntiHelium3"  ; chrg_ = -2; mass_ =  2.809413500; break;
        case PartType::Helium4       : name_ = "Helium4"      ; chrg_ =  2; mass_ =  3.727379240; break;
        case PartType::AntiHelium4   : name_ = "AntiHelium4"  ; chrg_ = -2; mass_ =  3.727379240; break;
        case PartType::Lithium6      : name_ = "Lithium6"     ; chrg_ =  3; mass_ =  5.603051363; break;
        case PartType::Lithium7      : name_ = "Lithium7"     ; chrg_ =  3; mass_ =  6.535366807; break;
        case PartType::Beryllium7    : name_ = "Beryllium7"   ; chrg_ =  4; mass_ =  6.536228700; break;
        case PartType::Beryllium9    : name_ = "Beryllium9"   ; chrg_ =  4; mass_ =  8.394794503; break;
        case PartType::Beryllium10   : name_ = "Beryllium10"  ; chrg_ =  4; mass_ =  9.327547622; break;
        case PartType::Boron10       : name_ = "Boron10"      ; chrg_ =  5; mass_ =  9.326991682; break;
        case PartType::Boron11       : name_ = "Boron11"      ; chrg_ =  5; mass_ = 10.255102976; break;
        case PartType::Carbon12      : name_ = "Carbon12"     ; chrg_ =  6; mass_ = 11.177929140; break;
        case PartType::Carbon13      : name_ = "Carbon13"     ; chrg_ =  6; mass_ = 12.112548247; break;
        case PartType::Carbon14      : name_ = "Carbon14"     ; chrg_ =  6; mass_ = 13.043937223; break;
        case PartType::Nitrogen13    : name_ = "Nitrogen13"   ; chrg_ =  7; mass_ = 12.114768715; break;
        case PartType::Nitrogen14    : name_ = "Nitrogen14"   ; chrg_ =  7; mass_ = 13.043780747; break;
        case PartType::Nitrogen15    : name_ = "Nitrogen15"   ; chrg_ =  7; mass_ = 13.972512863; break;
        case PartType::Oxygen16      : name_ = "Oxygen16"     ; chrg_ =  8; mass_ = 14.899168518; break;
        case PartType::Oxygen17      : name_ = "Oxygen17"     ; chrg_ =  8; mass_ = 15.834590801; break;
        case PartType::Oxygen18      : name_ = "Oxygen18"     ; chrg_ =  8; mass_ = 16.766112187; break;
        default : CERR("(ERROR) PartInfo : It is not in list of PartType.\n");
    }
    
    is_chrgless_ = Numc::EqualToZero(chrg_);
    is_massless_ = Numc::EqualToZero(mass_);
    mass_to_chrg_ = ((is_chrgless_ || is_massless_) ? 0.0 : std::fabs(mass_ / chrg_));
    chrg_to_mass_ = ((is_chrgless_ || is_massless_) ? 0.0 : std::fabs(chrg_ / mass_));
}


void PartInfo::print() const {
    std::string printStr;
    printStr += STR("==== %-15s ====\n", name_.c_str());
    printStr += STR("Chrg %3d\n"   , chrg_);
    printStr += STR("Mass %10.6f\n", mass_);
    printStr += STR("M/Q  %10.6f\n", mass_to_chrg_);
    printStr += STR("=========================\n");
    COUT(printStr.c_str());
}


} // namespace TrackSys


#endif // __TRACKLibs_PartInfo_C__
