#ifdef __HAS_TESTFIT__


#ifndef __TRACKLibs_MatEnvTestFit_H__
#define __TRACKLibs_MatEnvTestFit_H__


namespace TrackSys {


namespace MatTestFit {
    constexpr Long64_t DIM = 3;
    
    constexpr std::array<Long64_t, DIM> TRL01_N   {   400,  400,        2 };
    constexpr std::array<Double_t, DIM> TRL01_MIN {  -100, -100,   99.985 };
    constexpr std::array<Double_t, DIM> TRL01_MAX {   100,  100,  100.015 };
    
    constexpr std::array<Long64_t, DIM> TRL02_N   {   400,  400,       2 };
    constexpr std::array<Double_t, DIM> TRL02_MIN {  -100, -100,  69.985 };
    constexpr std::array<Double_t, DIM> TRL02_MAX {   100,  100,  70.015 };
    
    constexpr std::array<Long64_t, DIM> TRL03_N   {   400,  400,       2 };
    constexpr std::array<Double_t, DIM> TRL03_MIN {  -100, -100,  59.985 };
    constexpr std::array<Double_t, DIM> TRL03_MAX {   100,  100,  60.015 };
    
    constexpr std::array<Long64_t, DIM> TRL04_N   {   400,  400,       2 };
    constexpr std::array<Double_t, DIM> TRL04_MIN {  -100, -100,  26.985 };
    constexpr std::array<Double_t, DIM> TRL04_MAX {   100,  100,  27.015 };
    
    constexpr std::array<Long64_t, DIM> TRL05_N   {   400,  400,       2 };
    constexpr std::array<Double_t, DIM> TRL05_MIN {  -100, -100,  23.985 };
    constexpr std::array<Double_t, DIM> TRL05_MAX {   100,  100,  24.015 };
    
    constexpr std::array<Long64_t, DIM> TRL06_N   {   400,  400,      2 };
    constexpr std::array<Double_t, DIM> TRL06_MIN {  -100, -100,  1.985 };
    constexpr std::array<Double_t, DIM> TRL06_MAX {   100,  100,  2.015 };
    
    constexpr std::array<Long64_t, DIM> TRL07_N   {   400,  400,      2 };
    constexpr std::array<Double_t, DIM> TRL07_MIN {  -100, -100, -2.015 };
    constexpr std::array<Double_t, DIM> TRL07_MAX {   100,  100, -1.985 };
    
    constexpr std::array<Long64_t, DIM> TRL08_N   {   400,  400,       2 };
    constexpr std::array<Double_t, DIM> TRL08_MIN {  -100, -100, -23.015 };
    constexpr std::array<Double_t, DIM> TRL08_MAX {   100,  100, -22.985 };
    
    constexpr std::array<Long64_t, DIM> TRL09_N   {   400,  400,       2 };
    constexpr std::array<Double_t, DIM> TRL09_MIN {  -100, -100, -27.015 };
    constexpr std::array<Double_t, DIM> TRL09_MAX {   100,  100, -26.985 };
    
    constexpr std::array<Long64_t, DIM> TRL10_N   {   400,  400,       2 };
    constexpr std::array<Double_t, DIM> TRL10_MIN {  -100, -100, -60.015 };
    constexpr std::array<Double_t, DIM> TRL10_MAX {   100,  100, -59.985 };
    
    constexpr std::array<Long64_t, DIM> TRL11_N   {   400,  400,       2 };
    constexpr std::array<Double_t, DIM> TRL11_MIN {  -100, -100, -70.015 };
    constexpr std::array<Double_t, DIM> TRL11_MAX {   100,  100, -69.985 };
    
    constexpr std::array<Long64_t, DIM> TRL12_N   {   400,  400,        2 };
    constexpr std::array<Double_t, DIM> TRL12_MIN {  -100, -100, -100.015 };
    constexpr std::array<Double_t, DIM> TRL12_MAX {   100,  100,  -99.985 };
    
    constexpr std::array<Long64_t, DIM> MATC01_N   {   400,  400,    16 };
    constexpr std::array<Double_t, DIM> MATC01_MIN {  -100, -100,  61.0 };
    constexpr std::array<Double_t, DIM> MATC01_MAX {   100,  100,  69.0 };
    
    constexpr std::array<Long64_t, DIM> MATC02_N   {   400,  400,    16 };
    constexpr std::array<Double_t, DIM> MATC02_MIN {  -100, -100, -69.0 };
    constexpr std::array<Double_t, DIM> MATC02_MAX {   100,  100, -61.0 };
    
    constexpr std::array<Long64_t, DIM> MATAL01_N   {   400,  400,       2 };
    constexpr std::array<Double_t, DIM> MATAL01_MIN {  -100, -100,  26.875 };
    constexpr std::array<Double_t, DIM> MATAL01_MAX {   100,  100,  26.925 };
    
    constexpr std::array<Long64_t, DIM> MATAL02_N   {   400,  400,       2 };
    constexpr std::array<Double_t, DIM> MATAL02_MIN {  -100, -100,  23.075 };
    constexpr std::array<Double_t, DIM> MATAL02_MAX {   100,  100,  23.125 };
    
    constexpr std::array<Long64_t, DIM> MATAL03_N   {   400,  400,      2 };
    constexpr std::array<Double_t, DIM> MATAL03_MIN {  -100, -100,  1.875 };
    constexpr std::array<Double_t, DIM> MATAL03_MAX {   100,  100,  1.925 };
    
    constexpr std::array<Long64_t, DIM> MATAL04_N   {   400,  400,      2 };
    constexpr std::array<Double_t, DIM> MATAL04_MIN {  -100, -100, -1.925 };
    constexpr std::array<Double_t, DIM> MATAL04_MAX {   100,  100, -1.875 };
    
    constexpr std::array<Long64_t, DIM> MATAL05_N   {   400,  400,       2 };
    constexpr std::array<Double_t, DIM> MATAL05_MIN {  -100, -100, -23.125 };
    constexpr std::array<Double_t, DIM> MATAL05_MAX {   100,  100, -23.075 };
    
    constexpr std::array<Long64_t, DIM> MATAL06_N   {   400,  400,       2 };
    constexpr std::array<Double_t, DIM> MATAL06_MIN {  -100, -100, -26.925 };
    constexpr std::array<Double_t, DIM> MATAL06_MAX {   100,  100, -26.875 };
}


class MatGeoBoxTestFit {
    public :
        MatGeoBoxTestFit() {}
        ~MatGeoBoxTestFit() {}

        static Bool_t CreateMatGeoBox();

        static Bool_t Load();

        inline static std::list<MatGeoBoxReader*>& Reader() { return reader_; }

    private :
        static Bool_t is_load_;
        static std::list<MatGeoBoxReader*> reader_;
        static MatGeoBoxReader reader_TRL01_;
        static MatGeoBoxReader reader_TRL02_;
        static MatGeoBoxReader reader_TRL03_;
        static MatGeoBoxReader reader_TRL04_;
        static MatGeoBoxReader reader_TRL05_;
        static MatGeoBoxReader reader_TRL06_;
        static MatGeoBoxReader reader_TRL07_;
        static MatGeoBoxReader reader_TRL08_;
        static MatGeoBoxReader reader_TRL09_;
        static MatGeoBoxReader reader_TRL10_;
        static MatGeoBoxReader reader_TRL11_;
        static MatGeoBoxReader reader_TRL12_;
        static MatGeoBoxReader reader_MATC01_;
        static MatGeoBoxReader reader_MATC02_;
        static MatGeoBoxReader reader_MATAL01_;
        static MatGeoBoxReader reader_MATAL02_;
        static MatGeoBoxReader reader_MATAL03_;
        static MatGeoBoxReader reader_MATAL04_;
        static MatGeoBoxReader reader_MATAL05_;
        static MatGeoBoxReader reader_MATAL06_;
};

Bool_t MatGeoBoxTestFit::is_load_ = false;
std::list<MatGeoBoxReader*> MatGeoBoxTestFit::reader_;
MatGeoBoxReader MatGeoBoxTestFit::reader_TRL01_;
MatGeoBoxReader MatGeoBoxTestFit::reader_TRL02_;
MatGeoBoxReader MatGeoBoxTestFit::reader_TRL03_;
MatGeoBoxReader MatGeoBoxTestFit::reader_TRL04_;
MatGeoBoxReader MatGeoBoxTestFit::reader_TRL05_;
MatGeoBoxReader MatGeoBoxTestFit::reader_TRL06_;
MatGeoBoxReader MatGeoBoxTestFit::reader_TRL07_;
MatGeoBoxReader MatGeoBoxTestFit::reader_TRL08_;
MatGeoBoxReader MatGeoBoxTestFit::reader_TRL09_;
MatGeoBoxReader MatGeoBoxTestFit::reader_TRL10_;
MatGeoBoxReader MatGeoBoxTestFit::reader_TRL11_;
MatGeoBoxReader MatGeoBoxTestFit::reader_TRL12_;
MatGeoBoxReader MatGeoBoxTestFit::reader_MATC01_;
MatGeoBoxReader MatGeoBoxTestFit::reader_MATC02_;
MatGeoBoxReader MatGeoBoxTestFit::reader_MATAL01_;
MatGeoBoxReader MatGeoBoxTestFit::reader_MATAL02_;
MatGeoBoxReader MatGeoBoxTestFit::reader_MATAL03_;
MatGeoBoxReader MatGeoBoxTestFit::reader_MATAL04_;
MatGeoBoxReader MatGeoBoxTestFit::reader_MATAL05_;
MatGeoBoxReader MatGeoBoxTestFit::reader_MATAL06_;


} // namespace TrackSys


#endif // __TRACKLibs_MatEnvTestFit_H__


#endif // __HAS_TESTFIT__
