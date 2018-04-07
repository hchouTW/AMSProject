#ifndef __TRACKLibs_PartInfo_C__
#define __TRACKLibs_PartInfo_C__


namespace TrackSys {


void PartInfo::reset(const PartType& type) {
    // Atomic mass unit u = 0.931494095 GeV/c^2
    switch (type) {
        case PartType::Self          : reset(PartType::Self, "", SelfChrg_, SelfMass_); break;
        case PartType::Photon        : reset(PartType::Photon       , "Photon"       ,  0,  0.000000000); break;
        case PartType::Electron      : reset(PartType::Electron     , "Electron"     , -1,  0.000510999); break;
        case PartType::Positron      : reset(PartType::Positron     , "Positron"     ,  1,  0.000510999); break;
        case PartType::Muon          : reset(PartType::Muon         , "Muon"         , -1,  0.105658367); break;
        case PartType::PionPlus      : reset(PartType::PionPlus     , "PionPlus"     ,  1,  0.139570180); break;
        case PartType::PionMinus     : reset(PartType::PionMinus    , "PionMinus"    , -1,  0.139570180); break;
        case PartType::KaonPlus      : reset(PartType::KaonPlus     , "KaonPlus"     ,  1,  0.493677000); break;
        case PartType::KaonMinus     : reset(PartType::KaonMinus    , "KaonMinus"    , -1,  0.493677000); break;
        case PartType::Proton        : reset(PartType::Proton       , "Proton"       ,  1,  0.938272297); break;
        case PartType::AntiProton    : reset(PartType::AntiProton   , "AntiProton"   , -1,  0.938272297); break;
        case PartType::Deuterium     : reset(PartType::Deuterium    , "Deuterium"    ,  1,  1.876123915); break;
        case PartType::AntiDeuterium : reset(PartType::AntiDeuterium, "AntiDeuterium", -1,  1.876123915); break;
        case PartType::Helium3       : reset(PartType::Helium3      , "Helium3"      ,  2,  2.809413500); break;
        case PartType::AntiHelium3   : reset(PartType::AntiHelium3  , "AntiHelium3"  , -2,  2.809413500); break;
        case PartType::Helium4       : reset(PartType::Helium4      , "Helium4"      ,  2,  3.727379240); break;
        case PartType::AntiHelium4   : reset(PartType::AntiHelium4  , "AntiHelium4"  , -2,  3.727379240); break;
        case PartType::Lithium6      : reset(PartType::Lithium6     , "Lithium6"     ,  3,  5.603051363); break;
        case PartType::Lithium7      : reset(PartType::Lithium7     , "Lithium7"     ,  3,  6.535366807); break;
        case PartType::Beryllium7    : reset(PartType::Beryllium7   , "Beryllium7"   ,  4,  6.536228700); break;
        case PartType::Beryllium9    : reset(PartType::Beryllium9   , "Beryllium9"   ,  4,  8.394794503); break;
        case PartType::Beryllium10   : reset(PartType::Beryllium10  , "Beryllium10"  ,  4,  9.327547622); break;
        case PartType::Boron10       : reset(PartType::Boron10      , "Boron10"      ,  5,  9.326991682); break;
        case PartType::Boron11       : reset(PartType::Boron11      , "Boron11"      ,  5, 10.255102976); break;
        case PartType::Carbon12      : reset(PartType::Carbon12     , "Carbon12"     ,  6, 11.177929140); break;
        case PartType::Carbon13      : reset(PartType::Carbon13     , "Carbon13"     ,  6, 12.112548247); break;
        case PartType::Carbon14      : reset(PartType::Carbon14     , "Carbon14"     ,  6, 13.043937223); break;
        case PartType::Nitrogen13    : reset(PartType::Nitrogen13   , "Nitrogen13"   ,  7, 12.114768715); break;
        case PartType::Nitrogen14    : reset(PartType::Nitrogen14   , "Nitrogen14"   ,  7, 13.043780747); break;
        case PartType::Nitrogen15    : reset(PartType::Nitrogen15   , "Nitrogen15"   ,  7, 13.972512863); break;
        case PartType::Oxygen16      : reset(PartType::Oxygen16     , "Oxygen16"     ,  8, 14.899168518); break;
        case PartType::Oxygen17      : reset(PartType::Oxygen17     , "Oxygen17"     ,  8, 15.834590801); break;
        case PartType::Oxygen18      : reset(PartType::Oxygen18     , "Oxygen18"     ,  8, 16.766112187); break;
        default : CERR("(ERROR) PartInfo : It is not in list of PartType.\n");
    }
}

void PartInfo::reset(const PartType& type, const std::string& name, Short_t chrg, Double_t mass) {
    type_ = type;
    name_ = name;
    chrg_ = chrg;
    mass_ = mass;
    is_chrgless_ = Numc::EqualToZero(chrg_);
    is_massless_ = Numc::EqualToZero(mass_);
    mass_to_chrg_ = ((is_chrgless_ || is_massless_) ? Numc::ZERO<> : std::fabs(mass_ / chrg_));
    chrg_to_mass_ = ((is_chrgless_ || is_massless_) ? Numc::ZERO<> : std::fabs(chrg_ / mass_));
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
