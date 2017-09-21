#ifdef __HAS_AMS_OFFICE_LIBS__


#ifndef __TRACKLibs_MatEnvAms_H__
#define __TRACKLibs_MatEnvAms_H__


namespace TrackSys {


namespace MatAms {
    constexpr Long64_t DIM = 3;
    constexpr Float_t TRDL_Z = 85.0;
    constexpr Float_t TRDU_Z = 145.0;
    constexpr Float_t TRDL_RADIUS = 75.0;
    constexpr Float_t TRDU_RADIUS = 100.0;
    constexpr Float_t TRD_SLOPE = (TRDU_RADIUS - TRDL_RADIUS) / (TRDU_Z - TRDL_Z);
    constexpr Float_t TRACKER_RADIUS = 70.0;
    constexpr Float_t MAGNETIC_RADIUS = 53.0;
    constexpr Float_t RICH_BOUND = 17;
    
    constexpr std::array<Long64_t, DIM> RAD_N   {   480,  480,    16 };
    constexpr std::array<Double_t, DIM> RAD_MIN {  -120, -120, 165.0 };
    constexpr std::array<Double_t, DIM> RAD_MAX {   120,  120, 173.0 };
    
    constexpr std::array<Long64_t, DIM> TRL1_N   {  280, 280,    14 };
    constexpr std::array<Double_t, DIM> TRL1_MIN {  -70, -70, 158.0 };
    constexpr std::array<Double_t, DIM> TRL1_MAX {   70,  70, 165.0 };

    constexpr std::array<Long64_t, DIM> UTRD_N   {  440,  440,    20 };
    constexpr std::array<Double_t, DIM> UTRD_MIN { -110, -110, 146.0 };
    constexpr std::array<Double_t, DIM> UTRD_MAX {  110,  110, 156.0 };

    constexpr std::array<Long64_t, DIM> TRD_N   {   400,  400,   120 };
    constexpr std::array<Double_t, DIM> TRD_MIN {  -100, -100,  85.0 };
    constexpr std::array<Double_t, DIM> TRD_MAX {   100,  100, 145.0 };
    
    constexpr std::array<Long64_t, DIM> LTRD_N   {  320, 320,   10 };
    constexpr std::array<Double_t, DIM> LTRD_MIN {  -80, -80, 79.0 };
    constexpr std::array<Double_t, DIM> LTRD_MAX {   80,  80, 84.0 };

    constexpr std::array<Long64_t, DIM> UTOF_N   { 260, 260,   34 };
    constexpr std::array<Double_t, DIM> UTOF_MIN { -65, -65, 60.0 };
    constexpr std::array<Double_t, DIM> UTOF_MAX {  65,  65, 77.0 };
    
    constexpr std::array<Long64_t, DIM> UITR_N   {  280, 280,    4 };
    constexpr std::array<Double_t, DIM> UITR_MIN {  -70, -70, 58.0 };
    constexpr std::array<Double_t, DIM> UITR_MAX {   70,  70, 60.0 };
    
    constexpr std::array<Long64_t, DIM> TRS1_N   { 212, 212,   10 };
    constexpr std::array<Double_t, DIM> TRS1_MIN { -53, -53, 25.0 };
    constexpr std::array<Double_t, DIM> TRS1_MAX {  53,  53, 30.0 };
    
    constexpr std::array<Long64_t, DIM> TRS2_N   { 212, 212,   10 };
    constexpr std::array<Double_t, DIM> TRS2_MIN { -53, -53, -3.0 };
    constexpr std::array<Double_t, DIM> TRS2_MAX {  53,  53,  2.0 };
    
    constexpr std::array<Long64_t, DIM> TRS3_N   { 212, 212,    10 };
    constexpr std::array<Double_t, DIM> TRS3_MIN { -53, -53, -30.0 };
    constexpr std::array<Double_t, DIM> TRS3_MAX {  53,  53, -25.0 };

    constexpr std::array<Long64_t, DIM> LITR_N   {  280, 280,    14 };
    constexpr std::array<Double_t, DIM> LITR_MIN {  -70, -70, -60.0 };
    constexpr std::array<Double_t, DIM> LITR_MAX {   70,  70, -54.0 };
    
    constexpr std::array<Long64_t, DIM> LTOF_N   {  260, 260,    24 };
    constexpr std::array<Double_t, DIM> LTOF_MIN {  -65, -65, -72.0 };
    constexpr std::array<Double_t, DIM> LTOF_MAX {   65,  65, -60.0 };
    
    constexpr std::array<Long64_t, DIM> NAF_N   {  68,  68,     6 };
    constexpr std::array<Double_t, DIM> NAF_MIN { -17, -17, -76.0 };
    constexpr std::array<Double_t, DIM> NAF_MAX {  17,  17, -73.0 };
    
    constexpr std::array<Long64_t, DIM> AGL_N   { 600, 600,     6 };
    constexpr std::array<Double_t, DIM> AGL_MIN { -70, -70, -76.0 };
    constexpr std::array<Double_t, DIM> AGL_MAX {  70,  70, -73.0 };
    
    constexpr std::array<Long64_t, DIM> PMT_N   {  268, 268,     20 };
    constexpr std::array<Double_t, DIM> PMT_MIN {  -67, -67, -132.0 };
    constexpr std::array<Double_t, DIM> PMT_MAX {   67,  67, -122.0 };

    constexpr std::array<Long64_t, DIM> TRL9_N   {  240, 160,      6 };
    constexpr std::array<Double_t, DIM> TRL9_MIN {  -60, -40, -137.5 };
    constexpr std::array<Double_t, DIM> TRL9_MAX {   60,  40, -134.5 };

    constexpr std::array<Long64_t, DIM> ECAL_N   {  168, 168,     52 };
    constexpr std::array<Double_t, DIM> ECAL_MIN {  -42, -42, -164.0 };
    constexpr std::array<Double_t, DIM> ECAL_MAX {   42,  42, -138.0 };
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
        static MatGeoBoxReader reader_AMS02RAD_ ;
        static MatGeoBoxReader reader_AMS02TRL1_;
        static MatGeoBoxReader reader_AMS02UTRD_;
        static MatGeoBoxReader reader_AMS02TRD_ ;
        static MatGeoBoxReader reader_AMS02LTRD_;
        static MatGeoBoxReader reader_AMS02UTOF_;
        static MatGeoBoxReader reader_AMS02UITR_;
        static MatGeoBoxReader reader_AMS02TRS1_;
        static MatGeoBoxReader reader_AMS02TRS2_;
        static MatGeoBoxReader reader_AMS02TRS3_;
        static MatGeoBoxReader reader_AMS02LITR_;
        static MatGeoBoxReader reader_AMS02LTOF_;
        static MatGeoBoxReader reader_AMS02NAF_;
        static MatGeoBoxReader reader_AMS02AGL_;
        static MatGeoBoxReader reader_AMS02PMT_ ;
        static MatGeoBoxReader reader_AMS02TRL9_;
        static MatGeoBoxReader reader_AMS02ECAL_;
};

Bool_t MatGeoBoxAms::is_load_ = false;
std::list<MatGeoBoxReader*> MatGeoBoxAms::reader_;
MatGeoBoxReader MatGeoBoxAms::reader_AMS02RAD_ ;
MatGeoBoxReader MatGeoBoxAms::reader_AMS02TRL1_;
MatGeoBoxReader MatGeoBoxAms::reader_AMS02UTRD_;
MatGeoBoxReader MatGeoBoxAms::reader_AMS02TRD_ ;
MatGeoBoxReader MatGeoBoxAms::reader_AMS02LTRD_;
MatGeoBoxReader MatGeoBoxAms::reader_AMS02UTOF_;
MatGeoBoxReader MatGeoBoxAms::reader_AMS02UITR_;
MatGeoBoxReader MatGeoBoxAms::reader_AMS02TRS1_;
MatGeoBoxReader MatGeoBoxAms::reader_AMS02TRS2_;
MatGeoBoxReader MatGeoBoxAms::reader_AMS02TRS3_;
MatGeoBoxReader MatGeoBoxAms::reader_AMS02LITR_;
MatGeoBoxReader MatGeoBoxAms::reader_AMS02LTOF_;
MatGeoBoxReader MatGeoBoxAms::reader_AMS02NAF_;
MatGeoBoxReader MatGeoBoxAms::reader_AMS02AGL_;
MatGeoBoxReader MatGeoBoxAms::reader_AMS02PMT_ ;
MatGeoBoxReader MatGeoBoxAms::reader_AMS02TRL9_;
MatGeoBoxReader MatGeoBoxAms::reader_AMS02ECAL_;


} // namespace TrackSys


#endif // __TRACKLibs_MatEnvAms_H__


#endif // __HAS_AMS_OFFICE_LIBS__
