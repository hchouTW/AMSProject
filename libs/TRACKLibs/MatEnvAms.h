#ifdef __HAS_AMS_OFFICE_LIBS__


#ifndef __TRACKLibs_MatEnvAms_H__
#define __TRACKLibs_MatEnvAms_H__


namespace TrackSys {


namespace MatAms {
    constexpr Long64_t DIM = 3;

    // TRACKER
    constexpr Float_t INNER_TRACKER_RADIUS = 53.0;
    
    constexpr std::array<Long64_t, DIM> TRL1_N   {  260, 180,       4 };
    constexpr std::array<Double_t, DIM> TRL1_MIN {  -65, -45,  158.90 };
    constexpr std::array<Double_t, DIM> TRL1_MAX {   65,  45,  158.94 };

    constexpr std::array<Long64_t, DIM> TRL2_N   {  260, 180,       4 };
    constexpr std::array<Double_t, DIM> TRL2_MIN {  -65, -45,   53.04 };
    constexpr std::array<Double_t, DIM> TRL2_MAX {   65,  45,   53.08 };

    constexpr std::array<Long64_t, DIM> TRL3_N   {  220, 220,       4 };
    constexpr std::array<Double_t, DIM> TRL3_MIN {  -55, -55,   29.21 };
    constexpr std::array<Double_t, DIM> TRL3_MAX {   55,  55,   29.25 };

    constexpr std::array<Long64_t, DIM> TRL4_N   {  220, 220,       4 };
    constexpr std::array<Double_t, DIM> TRL4_MIN {  -55, -55,   25.19 };
    constexpr std::array<Double_t, DIM> TRL4_MAX {   55,  55,   25.23 };

    constexpr std::array<Long64_t, DIM> TRL5_N   {  220, 180,       4 };
    constexpr std::array<Double_t, DIM> TRL5_MIN {  -55, -45,    1.68 };
    constexpr std::array<Double_t, DIM> TRL5_MAX {   55,  45,    1.72 };

    constexpr std::array<Long64_t, DIM> TRL6_N   {  220, 180,       4 };
    constexpr std::array<Double_t, DIM> TRL6_MIN {  -55, -45,   -2.34 };
    constexpr std::array<Double_t, DIM> TRL6_MAX {   55,  45,   -2.30 };

    constexpr std::array<Long64_t, DIM> TRL7_N   {  220, 220,       4 };
    constexpr std::array<Double_t, DIM> TRL7_MIN {  -55, -55,  -25.23 };
    constexpr std::array<Double_t, DIM> TRL7_MAX {   55,  55,  -25.19 };

    constexpr std::array<Long64_t, DIM> TRL8_N   {  220, 220,       4 };
    constexpr std::array<Double_t, DIM> TRL8_MIN {  -55, -55,  -29.25 };
    constexpr std::array<Double_t, DIM> TRL8_MAX {   55,  55,  -29.21 };
    
    constexpr std::array<Long64_t, DIM> TRL9_N   {  200, 120,       4 };
    constexpr std::array<Double_t, DIM> TRL9_MIN {  -50, -30, -135.90 };
    constexpr std::array<Double_t, DIM> TRL9_MAX {   50,  30, -135.86 };
    
    // TRD
    constexpr Float_t TRDL_Z = 85.0;
    constexpr Float_t TRDU_Z = 145.0;
    constexpr Float_t TRDL_RADIUS = 65.0;
    constexpr Float_t TRDU_RADIUS = 105.0;
    constexpr Float_t TRD_SLOPE = (TRDU_RADIUS - TRDL_RADIUS) / (TRDU_Z - TRDL_Z);
    constexpr Float_t TRD_FACTOR = 1.5;
    
    constexpr std::array<Long64_t, DIM> TRDS_N   {   440,  440,    12 };
    constexpr std::array<Double_t, DIM> TRDS_MIN {  -110, -110, 144.0 };
    constexpr std::array<Double_t, DIM> TRDS_MAX {   110,  110, 156.0 };

    constexpr std::array<Long64_t, DIM> TRDU_N   {   440,  440,    12 };
    constexpr std::array<Double_t, DIM> TRDU_MIN {  -110, -110, 132.0 };
    constexpr std::array<Double_t, DIM> TRDU_MAX {   110,  110, 144.0 };
    
    constexpr std::array<Long64_t, DIM> TRDM_N   {   400,  400,    18 };
    constexpr std::array<Double_t, DIM> TRDM_MIN {  -100, -100, 114.0 };
    constexpr std::array<Double_t, DIM> TRDM_MAX {   100,  100, 132.0 };
    
    constexpr std::array<Long64_t, DIM> TRDI_N   {   400,  400,    18 };
    constexpr std::array<Double_t, DIM> TRDI_MIN {  -100, -100,  96.0 };
    constexpr std::array<Double_t, DIM> TRDI_MAX {   100,  100, 114.0 };
    
    constexpr std::array<Long64_t, DIM> TRDL_N   {  360, 360,    12 };
    constexpr std::array<Double_t, DIM> TRDL_MIN {  -90, -90,  84.0 };
    constexpr std::array<Double_t, DIM> TRDL_MAX {   90,  90,  96.0 };

    // TOF
    constexpr std::array<Long64_t, DIM> TOFU_N   { 260, 260,   13 };
    constexpr std::array<Double_t, DIM> TOFU_MIN { -65, -65, 60.0 };
    constexpr std::array<Double_t, DIM> TOFU_MAX {  65,  65, 66.5 };
    
    constexpr std::array<Long64_t, DIM> TOFL_N   { 260, 260,    13 };
    constexpr std::array<Double_t, DIM> TOFL_MIN { -65, -65, -66.5 };
    constexpr std::array<Double_t, DIM> TOFL_MAX {  65,  65, -59.5 };
    
    // RICH
    constexpr Float_t RICH_BOUND_INNER = 18;
    constexpr Float_t RICH_BOUND_OUTER = 65;

    constexpr std::array<Long64_t, DIM> NAF_N   {  72,  72,     7 };
    constexpr std::array<Double_t, DIM> NAF_MIN { -18, -18, -76.0 };
    constexpr std::array<Double_t, DIM> NAF_MAX {  18,  18, -72.5 };
    
    constexpr std::array<Long64_t, DIM> AGL_N   { 260, 260,     7 };
    constexpr std::array<Double_t, DIM> AGL_MIN { -65, -65, -76.0 };
    constexpr std::array<Double_t, DIM> AGL_MAX {  65,  65, -72.5 };
    
    // ECAL
    constexpr std::array<Long64_t, DIM> ECAL_N   {  168, 168,     22 };
    constexpr std::array<Double_t, DIM> ECAL_MIN {  -42, -42, -164.0 };
    constexpr std::array<Double_t, DIM> ECAL_MAX {   42,  42, -142.0 };

    // RAD
    constexpr std::array<Long64_t, DIM> RAD_N   {   480,  480,    10 };
    constexpr std::array<Double_t, DIM> RAD_MIN {  -120, -120, 163.0 };
    constexpr std::array<Double_t, DIM> RAD_MAX {   120,  120, 173.0 };
    
    // PMT
    constexpr std::array<Long64_t, DIM> PMT_N   {  260, 260,     10 };
    constexpr std::array<Double_t, DIM> PMT_MIN {  -65, -65, -132.0 };
    constexpr std::array<Double_t, DIM> PMT_MAX {   65,  65, -122.0 };
    
    // SUPPORT U
    constexpr std::array<Long64_t, DIM> SUPU1_N   {  280, 200,      2 };
    constexpr std::array<Double_t, DIM> SUPU1_MIN {  -70, -50, 158.98 };
    constexpr std::array<Double_t, DIM> SUPU1_MAX {   70,  50, 160.00 };

    constexpr std::array<Long64_t, DIM> SUPU2_N   {  280, 200,      1 };
    constexpr std::array<Double_t, DIM> SUPU2_MIN {  -70, -50, 158.79 };
    constexpr std::array<Double_t, DIM> SUPU2_MAX {   70,  50, 158.81 };

    constexpr std::array<Long64_t, DIM> SUPU3_N   {  220, 220,    7 };
    constexpr std::array<Double_t, DIM> SUPU3_MIN {  -55, -55, 66.5 };
    constexpr std::array<Double_t, DIM> SUPU3_MAX {   55,  55, 84.0 };
    
    // SUPPORT M
    constexpr std::array<Long64_t, DIM> SUPM1_N   {  220, 220,     6 };
    constexpr std::array<Double_t, DIM> SUPM1_MIN {  -55, -55, 28.60 };
    constexpr std::array<Double_t, DIM> SUPM1_MAX {   55,  55, 29.20 };
    
    constexpr std::array<Long64_t, DIM> SUPM2_N   {  220, 220,     6 };
    constexpr std::array<Double_t, DIM> SUPM2_MIN {  -55, -55, 25.25 };
    constexpr std::array<Double_t, DIM> SUPM2_MAX {   55,  55, 25.85 };
    
    constexpr std::array<Long64_t, DIM> SUPM3_N   {  220, 220,     6 };
    constexpr std::array<Double_t, DIM> SUPM3_MIN {  -55, -55,  1.05 };
    constexpr std::array<Double_t, DIM> SUPM3_MAX {   55,  55,  1.65 };
    
    constexpr std::array<Long64_t, DIM> SUPM4_N   {  220, 220,     6 };
    constexpr std::array<Double_t, DIM> SUPM4_MIN {  -55, -55, -2.30 };
    constexpr std::array<Double_t, DIM> SUPM4_MAX {   55,  55, -1.70 };
    
    constexpr std::array<Long64_t, DIM> SUPM5_N   {  220, 220,      6 };
    constexpr std::array<Double_t, DIM> SUPM5_MIN {  -55, -55, -25.85 };
    constexpr std::array<Double_t, DIM> SUPM5_MAX {   55,  55, -25.25 };
    
    constexpr std::array<Long64_t, DIM> SUPM6_N   {  220, 220,      6 };
    constexpr std::array<Double_t, DIM> SUPM6_MIN {  -55, -55, -29.20 };
    constexpr std::array<Double_t, DIM> SUPM6_MAX {   55,  55, -28.60 };
    
    // SUPPORT T
    constexpr std::array<Long64_t, DIM> SUPT1_N   {  260, 180,     2 };
    constexpr std::array<Double_t, DIM> SUPT1_MIN {  -65, -45, 52.93 };
    constexpr std::array<Double_t, DIM> SUPT1_MAX {   65,  45, 52.95 };
    
    constexpr std::array<Long64_t, DIM> SUPT2_N   {  220, 220,     2 };
    constexpr std::array<Double_t, DIM> SUPT2_MIN {  -55, -55, 29.34 };
    constexpr std::array<Double_t, DIM> SUPT2_MAX {   55,  55, 29.36 };
    
    constexpr std::array<Long64_t, DIM> SUPT3_N   {  220, 220,     2 };
    constexpr std::array<Double_t, DIM> SUPT3_MIN {  -55, -55, 25.08 };
    constexpr std::array<Double_t, DIM> SUPT3_MAX {   55,  55, 25.10 };
    
    constexpr std::array<Long64_t, DIM> SUPT4_N   {  220, 220,     2 };
    constexpr std::array<Double_t, DIM> SUPT4_MIN {  -55, -55,  1.81 };
    constexpr std::array<Double_t, DIM> SUPT4_MAX {   55,  55,  1.83 };
    
    constexpr std::array<Long64_t, DIM> SUPT5_N   {  220, 220,     2 };
    constexpr std::array<Double_t, DIM> SUPT5_MIN {  -55, -55, -2.45 };
    constexpr std::array<Double_t, DIM> SUPT5_MAX {   55,  55, -2.43 };
    
    constexpr std::array<Long64_t, DIM> SUPT6_N   {  220, 220,      2 };
    constexpr std::array<Double_t, DIM> SUPT6_MIN {  -55, -55, -25.10 };
    constexpr std::array<Double_t, DIM> SUPT6_MAX {   55,  55, -25.08 };
    
    constexpr std::array<Long64_t, DIM> SUPT7_N   {  220, 220,      2 };
    constexpr std::array<Double_t, DIM> SUPT7_MIN {  -55, -55, -29.36 };
    constexpr std::array<Double_t, DIM> SUPT7_MAX {   55,  55, -29.34 };

    // SUPPORT L
    constexpr std::array<Long64_t, DIM> SUPL1_N   {   280,  280,     2 };
    constexpr std::array<Double_t, DIM> SUPL1_MIN {   -70,  -70, -54.2 };
    constexpr std::array<Double_t, DIM> SUPL1_MAX {    70,   70, -54.0 };

    constexpr std::array<Long64_t, DIM> SUPL2_N   {   240,  240,     1 };
    constexpr std::array<Double_t, DIM> SUPL2_MIN {   -60,  -60, -60.0 };
    constexpr std::array<Double_t, DIM> SUPL2_MAX {    60,   60, -59.5 };

    constexpr std::array<Long64_t, DIM> SUPL3_N   {   240,  240,     6 };
    constexpr std::array<Double_t, DIM> SUPL3_MIN {   -60,  -60, -72.5 };
    constexpr std::array<Double_t, DIM> SUPL3_MAX {    60,   60, -66.5 };
    
    constexpr std::array<Long64_t, DIM> SUPL4_N   {   240,  160,       5 };
    constexpr std::array<Double_t, DIM> SUPL4_MIN {   -60,  -40, -135.85 };
    constexpr std::array<Double_t, DIM> SUPL4_MAX {    60,   40, -135.10 };

    constexpr std::array<Long64_t, DIM> SUPL5_N   {   240,  160,      3 };
    constexpr std::array<Double_t, DIM> SUPL5_MIN {   -60,  -40, -138.9 };
    constexpr std::array<Double_t, DIM> SUPL5_MAX {    60,   40, -135.9 };
}


class MatGeoBoxAms {
    public :
        MatGeoBoxAms() {}
        ~MatGeoBoxAms() {}

        static Bool_t CreateMatGeoBoxFromG4MatTree();

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
        
        static MatGeoBoxReader reader_TRDS_;
        static MatGeoBoxReader reader_TRDU_;
        static MatGeoBoxReader reader_TRDM_;
        static MatGeoBoxReader reader_TRDI_;
        static MatGeoBoxReader reader_TRDL_;
        
        static MatGeoBoxReader reader_TOFU_;
        static MatGeoBoxReader reader_TOFL_;
        
        static MatGeoBoxReader reader_NAF_;
        static MatGeoBoxReader reader_AGL_;
        
        static MatGeoBoxReader reader_ECAL_;
        
        static MatGeoBoxReader reader_RAD_;
        
        static MatGeoBoxReader reader_PMT_;
        
        static MatGeoBoxReader reader_SUPU1_;
        static MatGeoBoxReader reader_SUPU2_;
        static MatGeoBoxReader reader_SUPU3_;
        
        static MatGeoBoxReader reader_SUPM1_;
        static MatGeoBoxReader reader_SUPM2_;
        static MatGeoBoxReader reader_SUPM3_;
        static MatGeoBoxReader reader_SUPM4_;
        static MatGeoBoxReader reader_SUPM5_;
        static MatGeoBoxReader reader_SUPM6_;
        
        static MatGeoBoxReader reader_SUPT1_;
        static MatGeoBoxReader reader_SUPT2_;
        static MatGeoBoxReader reader_SUPT3_;
        static MatGeoBoxReader reader_SUPT4_;
        static MatGeoBoxReader reader_SUPT5_;
        static MatGeoBoxReader reader_SUPT6_;
        static MatGeoBoxReader reader_SUPT7_;

        static MatGeoBoxReader reader_SUPL1_;
        static MatGeoBoxReader reader_SUPL2_;
        static MatGeoBoxReader reader_SUPL3_;
        static MatGeoBoxReader reader_SUPL4_;
        static MatGeoBoxReader reader_SUPL5_;
};

Bool_t MatGeoBoxAms::is_load_ = false;
std::list<MatGeoBoxReader*> MatGeoBoxAms::reader_;

MatGeoBoxReader MatGeoBoxAms::reader_TRL1_;
MatGeoBoxReader MatGeoBoxAms::reader_TRL2_;
MatGeoBoxReader MatGeoBoxAms::reader_TRL3_;
MatGeoBoxReader MatGeoBoxAms::reader_TRL4_;
MatGeoBoxReader MatGeoBoxAms::reader_TRL5_;
MatGeoBoxReader MatGeoBoxAms::reader_TRL6_;
MatGeoBoxReader MatGeoBoxAms::reader_TRL7_;
MatGeoBoxReader MatGeoBoxAms::reader_TRL8_;
MatGeoBoxReader MatGeoBoxAms::reader_TRL9_;

MatGeoBoxReader MatGeoBoxAms::reader_TRDS_;
MatGeoBoxReader MatGeoBoxAms::reader_TRDU_;
MatGeoBoxReader MatGeoBoxAms::reader_TRDM_;
MatGeoBoxReader MatGeoBoxAms::reader_TRDI_;
MatGeoBoxReader MatGeoBoxAms::reader_TRDL_;

MatGeoBoxReader MatGeoBoxAms::reader_TOFU_;
MatGeoBoxReader MatGeoBoxAms::reader_TOFL_;

MatGeoBoxReader MatGeoBoxAms::reader_NAF_;
MatGeoBoxReader MatGeoBoxAms::reader_AGL_;

MatGeoBoxReader MatGeoBoxAms::reader_ECAL_;

MatGeoBoxReader MatGeoBoxAms::reader_RAD_;

MatGeoBoxReader MatGeoBoxAms::reader_PMT_;

MatGeoBoxReader MatGeoBoxAms::reader_SUPU1_;
MatGeoBoxReader MatGeoBoxAms::reader_SUPU2_;
MatGeoBoxReader MatGeoBoxAms::reader_SUPU3_;

MatGeoBoxReader MatGeoBoxAms::reader_SUPM1_;
MatGeoBoxReader MatGeoBoxAms::reader_SUPM2_;
MatGeoBoxReader MatGeoBoxAms::reader_SUPM3_;
MatGeoBoxReader MatGeoBoxAms::reader_SUPM4_;
MatGeoBoxReader MatGeoBoxAms::reader_SUPM5_;
MatGeoBoxReader MatGeoBoxAms::reader_SUPM6_;

MatGeoBoxReader MatGeoBoxAms::reader_SUPT1_;
MatGeoBoxReader MatGeoBoxAms::reader_SUPT2_;
MatGeoBoxReader MatGeoBoxAms::reader_SUPT3_;
MatGeoBoxReader MatGeoBoxAms::reader_SUPT4_;
MatGeoBoxReader MatGeoBoxAms::reader_SUPT5_;
MatGeoBoxReader MatGeoBoxAms::reader_SUPT6_;
MatGeoBoxReader MatGeoBoxAms::reader_SUPT7_;

MatGeoBoxReader MatGeoBoxAms::reader_SUPL1_;
MatGeoBoxReader MatGeoBoxAms::reader_SUPL2_;
MatGeoBoxReader MatGeoBoxAms::reader_SUPL3_;
MatGeoBoxReader MatGeoBoxAms::reader_SUPL4_;
MatGeoBoxReader MatGeoBoxAms::reader_SUPL5_;

} // namespace TrackSys


#endif // __TRACKLibs_MatEnvAms_H__


#endif // __HAS_AMS_OFFICE_LIBS__
