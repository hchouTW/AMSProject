#ifdef __HAS_TESTPROP__


#ifndef __TRACKLibs_MatEnvTestProp_H__
#define __TRACKLibs_MatEnvTestProp_H__


namespace TrackSys {


namespace MatTestProp {
    constexpr Long64_t DIM = 3;
    
    constexpr std::array<Long64_t, DIM> TRL1_N   {   400,  400,       8 };
    constexpr std::array<Double_t, DIM> TRL1_MIN {  -100, -100,  57.985 };
    constexpr std::array<Double_t, DIM> TRL1_MAX {   100,  100,  58.015 };
    
    constexpr std::array<Long64_t, DIM> TRL2_N   {   400,  400,       8 };
    constexpr std::array<Double_t, DIM> TRL2_MIN {  -100, -100,  49.985 };
    constexpr std::array<Double_t, DIM> TRL2_MAX {   100,  100,  50.015 };
    
    constexpr std::array<Long64_t, DIM> TRL3_N   {   400,  400,      8 };
    constexpr std::array<Double_t, DIM> TRL3_MIN {  -100, -100, -0.015 };
    constexpr std::array<Double_t, DIM> TRL3_MAX {   100,  100,  0.015 };
    
    constexpr std::array<Long64_t, DIM> MATC_N   {   400,  400,    20 };
    constexpr std::array<Double_t, DIM> MATC_MIN {  -100, -100,  50.5 };
    constexpr std::array<Double_t, DIM> MATC_MAX {   100,  100,  57.5 };
}


class MatGeoBoxTestProp {
    public :
        MatGeoBoxTestProp() {}
        ~MatGeoBoxTestProp() {}

        static Bool_t CreateMatGeoBox();

        static Bool_t Load();

        inline static std::list<MatGeoBoxReader*>& Reader() { return reader_; }

    private :
        static Bool_t is_load_;
        static std::list<MatGeoBoxReader*> reader_;
        static MatGeoBoxReader reader_TRL1_;
        static MatGeoBoxReader reader_TRL2_;
        static MatGeoBoxReader reader_TRL3_;
        static MatGeoBoxReader reader_MATC_;
};

Bool_t MatGeoBoxTestProp::is_load_ = false;
std::list<MatGeoBoxReader*> MatGeoBoxTestProp::reader_;
MatGeoBoxReader MatGeoBoxTestProp::reader_TRL1_;
MatGeoBoxReader MatGeoBoxTestProp::reader_TRL2_;
MatGeoBoxReader MatGeoBoxTestProp::reader_TRL3_;
MatGeoBoxReader MatGeoBoxTestProp::reader_MATC_;


} // namespace TrackSys


#endif // __TRACKLibs_MatEnvTestProp_H__


#endif // __HAS_TESTROP__
