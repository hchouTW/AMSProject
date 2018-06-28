#ifndef __TRACKLibs_MatEnvAms_H__
#define __TRACKLibs_MatEnvAms_H__

#include <TFile.h>


namespace TrackSys {

namespace MatAms {
    // TRK
    static constexpr Double_t                             TRL1_STP = Numc::ONE<> / Numc::SIX<>;
    static constexpr std::array<Long64_t, MATGEOBOX_NDIM> TRL1_N   {  65,  45,            1 };
    static constexpr std::array<Double_t, MATGEOBOX_NDIM> TRL1_MIN { -65, -45, 158.90500032 };
    static constexpr std::array<Double_t, MATGEOBOX_NDIM> TRL1_MAX {  65,  45, 158.93500031 };
    
    static constexpr Double_t                             TRL2_STP = Numc::ONE<> / Numc::SIX<>;
    static constexpr std::array<Long64_t, MATGEOBOX_NDIM> TRL2_N   {  65,  45,           1 };
    static constexpr std::array<Double_t, MATGEOBOX_NDIM> TRL2_MIN { -65, -45, 53.04499947 };
    static constexpr std::array<Double_t, MATGEOBOX_NDIM> TRL2_MAX {  65,  45, 53.07499947 };
    
    static constexpr Double_t                             TRL3_STP = Numc::ONE<> / Numc::SIX<>;
    static constexpr std::array<Long64_t, MATGEOBOX_NDIM> TRL3_N   {  55,  55,           1 };
    static constexpr std::array<Double_t, MATGEOBOX_NDIM> TRL3_MIN { -55, -55, 29.21300064 };
    static constexpr std::array<Double_t, MATGEOBOX_NDIM> TRL3_MAX {  55,  55, 29.24300064 };
    
    static constexpr Double_t                             TRL4_STP = Numc::ONE<> / Numc::SIX<>;
    static constexpr std::array<Long64_t, MATGEOBOX_NDIM> TRL4_N   {  55,  55,           1 };
    static constexpr std::array<Double_t, MATGEOBOX_NDIM> TRL4_MIN { -55, -55, 25.19700061 };
    static constexpr std::array<Double_t, MATGEOBOX_NDIM> TRL4_MAX {  55,  55, 25.22700061 };
    
    static constexpr Double_t                             TRL5_STP = Numc::ONE<> / Numc::SIX<>;
    static constexpr std::array<Long64_t, MATGEOBOX_NDIM> TRL5_N   {  55,  55,          1 };
    static constexpr std::array<Double_t, MATGEOBOX_NDIM> TRL5_MIN { -55, -55, 1.68299995 };
    static constexpr std::array<Double_t, MATGEOBOX_NDIM> TRL5_MAX {  55,  55, 1.71299995 };
    
    static constexpr Double_t                             TRL6_STP = Numc::ONE<> / Numc::SIX<>;
    static constexpr std::array<Long64_t, MATGEOBOX_NDIM> TRL6_N   {  55,  55,           1 };
    static constexpr std::array<Double_t, MATGEOBOX_NDIM> TRL6_MIN { -55, -55, -2.33300008 };
    static constexpr std::array<Double_t, MATGEOBOX_NDIM> TRL6_MAX {  55,  55, -2.30300008 };
    
    static constexpr Double_t                             TRL7_STP = Numc::ONE<> / Numc::SIX<>;
    static constexpr std::array<Long64_t, MATGEOBOX_NDIM> TRL7_N   {  55,  55,            1 };
    static constexpr std::array<Double_t, MATGEOBOX_NDIM> TRL7_MIN { -55, -55, -25.22700061 };
    static constexpr std::array<Double_t, MATGEOBOX_NDIM> TRL7_MAX {  55,  55, -25.19700061 };
    
    static constexpr Double_t                             TRL8_STP = Numc::ONE<> / Numc::SIX<>;
    static constexpr std::array<Long64_t, MATGEOBOX_NDIM> TRL8_N   {  55,  55,            1 };
    static constexpr std::array<Double_t, MATGEOBOX_NDIM> TRL8_MIN { -55, -55, -29.24300064 };
    static constexpr std::array<Double_t, MATGEOBOX_NDIM> TRL8_MAX {  55,  55, -29.21300064 };
    
    static constexpr Double_t                             TRL9_STP = Numc::ONE<> / Numc::SIX<>;
    static constexpr std::array<Long64_t, MATGEOBOX_NDIM> TRL9_N   {  50,  30,             1 };
    static constexpr std::array<Double_t, MATGEOBOX_NDIM> TRL9_MIN { -50, -30, -135.89699270 };
    static constexpr std::array<Double_t, MATGEOBOX_NDIM> TRL9_MAX {  50,  30, -135.86699270 };
    
    // TR Inner
    static constexpr Double_t                             TR34_STP = Numc::ONE<> / Numc::FOUR<>;
    static constexpr std::array<Long64_t, MATGEOBOX_NDIM> TR34_N   {   370,   370,           8 };
    static constexpr std::array<Double_t, MATGEOBOX_NDIM> TR34_MIN { -55.5, -55.5, TRL4_MAX[2] };
    static constexpr std::array<Double_t, MATGEOBOX_NDIM> TR34_MAX {  55.5,  55.5, TRL3_MIN[2] };
    
    static constexpr Double_t                             TR56_STP = Numc::ONE<> / Numc::FOUR<>;
    static constexpr std::array<Long64_t, MATGEOBOX_NDIM> TR56_N   {   370,   370,           8 };
    static constexpr std::array<Double_t, MATGEOBOX_NDIM> TR56_MIN { -55.5, -55.5, TRL6_MAX[2] };
    static constexpr std::array<Double_t, MATGEOBOX_NDIM> TR56_MAX {  55.5,  55.5, TRL5_MIN[2] };
    
    static constexpr Double_t                             TR78_STP = Numc::ONE<> / Numc::FOUR<>;
    static constexpr std::array<Long64_t, MATGEOBOX_NDIM> TR78_N   {   370,   370,           8 };
    static constexpr std::array<Double_t, MATGEOBOX_NDIM> TR78_MIN { -55.5, -55.5, TRL8_MAX[2] };
    static constexpr std::array<Double_t, MATGEOBOX_NDIM> TR78_MAX {  55.5,  55.5, TRL7_MIN[2] };
    
    static constexpr Double_t                             TRS1_STP = Numc::ONE<> / Numc::THREE<>;
    static constexpr std::array<Long64_t, MATGEOBOX_NDIM> TRS1_N   {  65,  45,            1 };
    static constexpr std::array<Double_t, MATGEOBOX_NDIM> TRS1_MIN { -65, -45, 158.79500019 };
    static constexpr std::array<Double_t, MATGEOBOX_NDIM> TRS1_MAX {  65,  45, 158.80500019 };

    static constexpr Double_t                             TRS2_STP = Numc::ONE<> / Numc::THREE<>;
    static constexpr std::array<Long64_t, MATGEOBOX_NDIM> TRS2_N   {  65,  45,           1 };
    static constexpr std::array<Double_t, MATGEOBOX_NDIM> TRS2_MIN { -65, -45, 52.93499934 };
    static constexpr std::array<Double_t, MATGEOBOX_NDIM> TRS2_MAX {  65,  45, 52.94499934 };
    
    static constexpr Double_t                             TRS3_STP = Numc::ONE<> / Numc::THREE<>;
    static constexpr std::array<Long64_t, MATGEOBOX_NDIM> TRS3_N   {  55,  55,           1 };
    static constexpr std::array<Double_t, MATGEOBOX_NDIM> TRS3_MIN { -55, -55, 29.34300065 };
    static constexpr std::array<Double_t, MATGEOBOX_NDIM> TRS3_MAX {  55,  55, 29.35300065 };
    
    static constexpr Double_t                             TRS4_STP = Numc::ONE<> / Numc::THREE<>;
    static constexpr std::array<Long64_t, MATGEOBOX_NDIM> TRS4_N   {  55,  55,           1 };
    static constexpr std::array<Double_t, MATGEOBOX_NDIM> TRS4_MIN { -55, -55, 25.08700072 };
    static constexpr std::array<Double_t, MATGEOBOX_NDIM> TRS4_MAX {  55,  55, 25.09700072 };
    
    static constexpr Double_t                             TRS5_STP = Numc::ONE<> / Numc::THREE<>;
    static constexpr std::array<Long64_t, MATGEOBOX_NDIM> TRS5_N   {  55,  55,           1 };
    static constexpr std::array<Double_t, MATGEOBOX_NDIM> TRS5_MIN { -55, -55,  1.81299996 };
    static constexpr std::array<Double_t, MATGEOBOX_NDIM> TRS5_MAX {  55,  55,  1.82299996 };

    static constexpr Double_t                             TRS6_STP = Numc::ONE<> / Numc::THREE<>;
    static constexpr std::array<Long64_t, MATGEOBOX_NDIM> TRS6_N   {  55,  55,           1 };
    static constexpr std::array<Double_t, MATGEOBOX_NDIM> TRS6_MIN { -55, -55, -2.44299996 };
    static constexpr std::array<Double_t, MATGEOBOX_NDIM> TRS6_MAX {  55,  55, -2.43299996 };
    
    static constexpr Double_t                             TRS7_STP = Numc::ONE<> / Numc::THREE<>;
    static constexpr std::array<Long64_t, MATGEOBOX_NDIM> TRS7_N   {  55,  55,            1 };
    static constexpr std::array<Double_t, MATGEOBOX_NDIM> TRS7_MIN { -55, -55, -25.09700072 };
    static constexpr std::array<Double_t, MATGEOBOX_NDIM> TRS7_MAX {  55,  55, -25.08700072 };
    
    static constexpr Double_t                             TRS8_STP = Numc::ONE<> / Numc::THREE<>;
    static constexpr std::array<Long64_t, MATGEOBOX_NDIM> TRS8_N   {  55,  55,            1 };
    static constexpr std::array<Double_t, MATGEOBOX_NDIM> TRS8_MIN { -55, -55, -29.35300065 };
    static constexpr std::array<Double_t, MATGEOBOX_NDIM> TRS8_MAX {  55,  55, -29.34300065 };
    
    static constexpr Double_t                             TRS9_STP = Numc::ONE<> / Numc::THREE<>;
    static constexpr std::array<Long64_t, MATGEOBOX_NDIM> TRS9_N   {  55,  55,           3 };
    static constexpr std::array<Double_t, MATGEOBOX_NDIM> TRS9_MIN { -55, -55, TRL9_MAX[2] };
    static constexpr std::array<Double_t, MATGEOBOX_NDIM> TRS9_MAX {  55,  55,     -134.75 };

    // RAD
    static constexpr Double_t                             RAD_STP = Numc::ONE<> / Numc::THREE<>;
    static constexpr std::array<Long64_t, MATGEOBOX_NDIM> RAD_N    {  120,  120,          10 };
    static constexpr std::array<Double_t, MATGEOBOX_NDIM> RAD_MIN  { -120, -120, TRL1_MAX[2] };
    static constexpr std::array<Double_t, MATGEOBOX_NDIM> RAD_MAX  {  120,  120,      173.00 };
    
    // TRD
    static constexpr Double_t                             TRD_STP = Numc::ONE<> / Numc::FOUR<>;
    static constexpr std::array<Long64_t, MATGEOBOX_NDIM> TRD_N    {   220,  220,    20 };
    static constexpr std::array<Double_t, MATGEOBOX_NDIM> TRD_MIN  {  -110, -110,  79.0 };
    static constexpr std::array<Double_t, MATGEOBOX_NDIM> TRD_MAX  {   110,  110, 155.5 };
    
    // SUPU
    static constexpr Double_t                             SUPU_STP = Numc::ONE<> / Numc::THREE<>;
    static constexpr std::array<Long64_t, MATGEOBOX_NDIM> SUPU_N   {  130, 130,     5 };
    static constexpr std::array<Double_t, MATGEOBOX_NDIM> SUPU_MIN {  -65, -65, 66.74 };
    static constexpr std::array<Double_t, MATGEOBOX_NDIM> SUPU_MAX {   65,  65, 76.70 };

    // TOFU
    static constexpr Double_t                             TOFU_STP = Numc::ONE<> / Numc::THREE<>;
    static constexpr std::array<Long64_t, MATGEOBOX_NDIM> TOFU_N   { 260, 260,    10 };
    static constexpr std::array<Double_t, MATGEOBOX_NDIM> TOFU_MIN { -65, -65, 60.24 };
    static constexpr std::array<Double_t, MATGEOBOX_NDIM> TOFU_MAX {  65,  65, 66.71 };
    
    // SPIU
    static constexpr Double_t                             SPIU_STP = Numc::ONE<> / Numc::THREE<>;
    static constexpr std::array<Long64_t, MATGEOBOX_NDIM> SPIU_N   {  130, 130,     5 };
    static constexpr std::array<Double_t, MATGEOBOX_NDIM> SPIU_MIN {  -65, -65, 53.12 };
    static constexpr std::array<Double_t, MATGEOBOX_NDIM> SPIU_MAX {   65,  65, 59.71 };
    
    // SPIL
    static constexpr Double_t                             SPIL_STP = Numc::ONE<> / Numc::THREE<>;
    static constexpr std::array<Long64_t, MATGEOBOX_NDIM> SPIL_N   {  130, 130,      5 };
    static constexpr std::array<Double_t, MATGEOBOX_NDIM> SPIL_MIN {  -65, -65, -58.21 };
    static constexpr std::array<Double_t, MATGEOBOX_NDIM> SPIL_MAX {   65,  65, -54.07 };
    
    // TOFL
    static constexpr Double_t                             TOFL_STP = Numc::ONE<> / Numc::THREE<>;
    static constexpr std::array<Long64_t, MATGEOBOX_NDIM> TOFL_N   { 260, 260,     10 };
    static constexpr std::array<Double_t, MATGEOBOX_NDIM> TOFL_MIN { -65, -65, -66.71 };
    static constexpr std::array<Double_t, MATGEOBOX_NDIM> TOFL_MAX {  65,  65, -60.24 };
    
    // SUPL
    static constexpr Double_t                             SUPL_STP = Numc::ONE<> / Numc::THREE<>;
    static constexpr std::array<Long64_t, MATGEOBOX_NDIM> SUPL_N   {  130, 130,      5 };
    static constexpr std::array<Double_t, MATGEOBOX_NDIM> SUPL_MIN {  -65, -65, -71.67 };
    static constexpr std::array<Double_t, MATGEOBOX_NDIM> SUPL_MAX {   65,  65, -66.82 };
    
    // RICH
    static constexpr Double_t                             RICH_STP = Numc::ONE<> / Numc::THREE<>;
    static constexpr std::array<Long64_t, MATGEOBOX_NDIM> RICH_N   {  65,  65,     6 };
    static constexpr std::array<Double_t, MATGEOBOX_NDIM> RICH_MIN { -65, -65, -75.9 };
    static constexpr std::array<Double_t, MATGEOBOX_NDIM> RICH_MAX {  65,  65, -72.8 };
    
    // PMT
    static constexpr Double_t                             PMT_STP = Numc::ONE<> / Numc::THREE<>;
    static constexpr std::array<Long64_t, MATGEOBOX_NDIM> PMT_N    {  65,  65,      5 };
    static constexpr std::array<Double_t, MATGEOBOX_NDIM> PMT_MIN  { -65, -65, -132.0 };
    static constexpr std::array<Double_t, MATGEOBOX_NDIM> PMT_MAX  {  65,  65, -122.8 };
    
    // ECAL
    static constexpr Double_t                             ECAL_STP = Numc::ONE<> / Numc::THREE<>;
    static constexpr std::array<Long64_t, MATGEOBOX_NDIM> ECAL_N   {  42,  42,          10 };
    static constexpr std::array<Double_t, MATGEOBOX_NDIM> ECAL_MIN { -42, -42,      -163.6 };
    static constexpr std::array<Double_t, MATGEOBOX_NDIM> ECAL_MAX {  42,  42, TRL9_MIN[2] };
}


class MatGeoBoxAms {
    public :
        MatGeoBoxAms() {}
        ~MatGeoBoxAms() {}

        static Bool_t CreateMatGeoBoxFromG4MatTree(const std::string& dpath = ".", const std::string& fpath = "g4mscan.root");

        static Bool_t Load();
        
        inline static std::list<MatGeoBoxReader*>& Reader() { return reader_; }

    private :
        static Bool_t is_load_;
        static std::list<MatGeoBoxReader*> reader_;

        static MatGeoBoxReader reader_TRL1_;
        static MatGeoBoxReader reader_TRL2_;
        static MatGeoBoxReader reader_TRL3_;
        static MatGeoBoxReader reader_TRL4_;
        static MatGeoBoxReader reader_TRL5_;
        static MatGeoBoxReader reader_TRL6_;
        static MatGeoBoxReader reader_TRL7_;
        static MatGeoBoxReader reader_TRL8_;
        static MatGeoBoxReader reader_TRL9_;
        
        static MatGeoBoxReader reader_TR34_;
        static MatGeoBoxReader reader_TR56_;
        static MatGeoBoxReader reader_TR78_;
        
        static MatGeoBoxReader reader_TRS1_;
        static MatGeoBoxReader reader_TRS2_;
        static MatGeoBoxReader reader_TRS3_;
        static MatGeoBoxReader reader_TRS4_;
        static MatGeoBoxReader reader_TRS5_;
        static MatGeoBoxReader reader_TRS6_;
        static MatGeoBoxReader reader_TRS7_;
        static MatGeoBoxReader reader_TRS8_;
        static MatGeoBoxReader reader_TRS9_;
        
        static MatGeoBoxReader reader_RAD_;
        static MatGeoBoxReader reader_TRD_;
        static MatGeoBoxReader reader_SUPU_;
        static MatGeoBoxReader reader_SPIU_;
        static MatGeoBoxReader reader_TOFU_;
        static MatGeoBoxReader reader_TOFL_;
        static MatGeoBoxReader reader_SPIL_;
        static MatGeoBoxReader reader_SUPL_;
        static MatGeoBoxReader reader_RICH_;
        static MatGeoBoxReader reader_PMT_;
        static MatGeoBoxReader reader_ECAL_;
};


} // namespace TrackSys


#endif // __TRACKLibs_MatEnvAms_H__
