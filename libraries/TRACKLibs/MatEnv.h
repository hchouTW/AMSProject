#ifndef __TRACKLibs_MatEnv_H__
#define __TRACKLibs_MatEnv_H__

namespace TrackSys {

namespace MatProperty {
    constexpr Int_t NUM_ELM = 9;
    const     std::array<std::string, NUM_ELM> NAME { "Hydrogen(H)", "Carbon(C)", "Nitrogen(N)", "Oxygen(O)", "Fluorine(F)", "Sodium(Na)", "Aluminum(Al)", "Silicon(Si)", "Lead(Pb)" };
    constexpr std::array<Int_t, NUM_ELM>       CHRG { 1, 6, 7, 8, 9, 11, 13, 14, 82 };
    constexpr std::array<Double_t, NUM_ELM>    MASS { 1.007947, 12.01078, 14.00672, 15.99943, 18.99840325, 22.989769282, 26.98153868, 28.08553, 207.21 }; // [g mol^-1]

    constexpr std::array<Double_t, NUM_ELM> RAD_LEN { 63.04, 42.70, 37.99, 34.24, 32.93, 27.74, 24.01, 21.82, 6.37 }; // [g cm^-2]
    
    constexpr std::array<Double_t, NUM_ELM> MEAN_EXENG        {  1.92e-05,  8.10e-05,  8.20e-05,  9.50e-05,  1.15e-04,  1.49e-04,  1.66e-04,  1.73e-04,  8.23e-04 }; // From NIST, ESTART
    constexpr std::array<Double_t, NUM_ELM> NEG_LN_MEAN_EXENG { 1.086e+01, 9.421e+00, 9.409e+00, 9.262e+00, 9.071e+00, 8.812e+00, 8.704e+00, 8.662e+00, 7.103e+00 }; // From NIST, ESTART

    constexpr std::array<Double_t, NUM_ELM> DEN_EFF_CORR_C    { 3.2632,  2.9925, 10.5400, 10.7004, 10.9653, 5.0526, 4.2395, 4.4351, 6.2018 }; // Form Geant4
    constexpr std::array<Double_t, NUM_ELM> DEN_EFF_CORR_X0   { 0.4759, -0.0351,  1.7378,  1.7541,  1.8433, 0.2880, 0.1708, 0.2014, 0.3776 }; // Form Geant4
    constexpr std::array<Double_t, NUM_ELM> DEN_EFF_CORR_X1   { 1.9215,  2.4860,  4.1323,  4.3213,  4.4096, 3.1962, 3.0127, 2.8715, 3.8073 }; // Form Geant4
    constexpr std::array<Double_t, NUM_ELM> DEN_EFF_CORR_A    { 0.1348,  0.2024,  0.1535,  0.1178,  0.1108, 0.0777, 0.0802, 0.1492, 0.0936 }; // Form Geant4
    constexpr std::array<Double_t, NUM_ELM> DEN_EFF_CORR_M    { 5.6249,  3.0036,  3.2125,  3.2913,  3.2962, 3.6452, 3.6345, 3.2546, 3.1608 }; // Form Geant4
    constexpr std::array<Double_t, NUM_ELM> DEN_EFF_CORR_DLT0 {   0.13,     0.1,    0.19,    0.11,    0.11,   0.08,   0.12,   0.14,   0.14 }; // Form Geant4
    constexpr std::array<Double_t, NUM_ELM> DEN_EFF_CORR_DLTM {  0.021,   0.038,  0.086,  0.101,     0.121,  0.098,  0.061,  0.059,  0.019 }; // Form Geant4
}


class MatFld {
    public :
        MatFld() : mat_(false), inv_rad_len_(0.) {}
        //MagFld(Bool_t mat, const std::array<Bool_t, MatProperty::NUM_ELM>& elm, const std::array<Double_t, MatProperty::NUM_ELM>& den, Double_t inv_rad_len = 0.0);
        ~MatFld() {}

        inline const Bool_t& operator() () const { return mat_; }
        inline const std::array<Double_t, MatProperty::NUM_ELM>& elm() const { return elm_; }
        inline const std::array<Double_t, MatProperty::NUM_ELM>& den() const { return den_; }
        inline const Double_t& inv_rad_len() const { return inv_rad_len_; }

    protected :
        inline void clear() { mat = false; elm_.fill(false); den_.fill(0.); inv_rad_len_ = 0.; real_len_ = 0.; efft_len_ = 0.; efft_ = 0.; }

    private :
        Bool_t                                     mat_;
        std::array<Bool_t, MatProperty::NUM_ELM>   elm_;
        std::array<Double_t, MatProperty::NUM_ELM> den_;
        Double_t                                   inv_rad_len_;
        Double_t                                   real_len_;
        Double_t                                   efft_len_;
        Double_t                                   efft_;
};


//======================================================================//
// MatGeoBox - Data Format                                              //
//----------------------------------------------------------------------//
// Header :                                                             //
// Long64_t N[3]                                                        //
// Double_t MIN[3]                                                      //
// Double_t MAX[3]                                                      //
// Bool_t   ELM[NUM_ELM]                                                //
// Double_t DEN[NUM_ELM]                                                //
// Double_t INV_RAD_LEN
// Content :                                                            //
// Bool_t   MAT[index]                                                  //
// ## Index = ( xi*YN*ZN + yi*ZN + zi )                                 //
// ##         xi[0 XN-1] yi[0 YN-1] zi[0 ZN-1]                          //
//======================================================================//
#include <sys/mman.h>
#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>


struct MatGeoBox {
    Long64_t n[3];
    Double_t min[3];
    Double_t max[3];
    Bool_t   elm[MatProperty::NUM_ELM];
    Double_t den[MatProperty::NUM_ELM];
    Double_t inv_rad_len;
    Bool_t   mat;
};


class MatGeoBoxCreator {
    public :
        MatGeoBoxCreator(Long64_t xn, Double_t xmin, Double_t xmax, Long64_t yn, Double_t ymin, Double_t ymax, Long64_t zn, Double_t zmin, Double_t zmax, const std::string& file_path = "MatGeoBox.bin");
        ~MatGeoBoxCreator() { save_and_close(); }

        void fill(Bool_t elm[MatProperty::NUM_ELM], Float_t den[MatProperty::NUM_ELM], Bool_t calculated = true); // from min to max coord

        inline Bool_t is_open() { return is_open_; }
       
        inline Bool_t is_in_box(Float_t coo[3]);

        void save_and_close();

    protected :
        inline void clear() { is_open_ = false; file_path_ = ""; file_des_ = -1; file_len_ = 0; file_ptr_ = reinterpret_cast<void*>(-1); max_len_ = 0; cur_len_= 0; geo_box_ = nullptr; elm_.fill(false); cnt_ = 0; den_.fill(0); inv_rad_len_ = 0.; }

    private :
        static constexpr Long64_t DIM_ = 3;
        Bool_t                    is_open_;
        std::string               file_path_;
        Int_t                     file_des_;
        Int_t                     file_len_;
        void*                     file_ptr_;
        Long64_t                  max_len_;
        Long64_t                  cur_len_;
        MatGeoBox*                geo_box_;
        Long64_t                                   cnt_;
        std::array<Bool_t,   MatProperty::NUM_ELM> elm_;
        std::array<Double_t, MatProperty::NUM_ELM> den_;
        Double_t                                   inv_rad_len_;
};



class MatGeoBoxReader {
    public :
        MatGeoBoxReader() { clear(); }
        MatGeoBoxReader(const std::string& file_path) { clear(); load(file_path); }
        ~MatGeoBoxReader() { clear(); }

        inline Bool_t exist() { return is_load_; }
        //inline MatFld get(const SVecD<3>& coo) { return MatFld( do_trilinear_interpolation( get_index(coo) ) ); }
        
        Bool_t load(const std::string& file_path);

        inline Bool_t is_cross(const SVecD<3>& vcoo, const SVecD<3>& wcoo) { return ((!is_load_) ? false : is_cross_box_coord(get_box_coord(vcoo), get_box_coord(wcoo))); }

    protected :
        inline void clear() { is_load_ = false; file_path_ = ""; file_ptr_ = reinterpret_cast<void*>(-1); mat_ptr_ = nullptr; n_.fill(0); min_.fill(0.); max_.fill(0.); len_.fill(0.); dlt_.fill(0.); fact_.fill(0); elm_.fill(false); den_.fill(0.); inv_rad_len_ = 0.; }

        inline SVecD<3> get_box_coord(const SVecD<3>& coo);
        inline SVecD<3> get_loc_coord(const SVecD<3>& coo);
        
        inline std::tuple<Long64_t, Bool_t> get_index(const SVecD<3>& loc);

        inline Bool_t is_cross_box_coord(const SVecD<3>& vbox, const SVecD<3>& wbox);

    private :
        static constexpr Long64_t                  DIM_ = 3;
        static constexpr Double_t                  STEP_LEN_ = 0.15;
        Bool_t                                     is_load_;
        std::string                                file_path_;
        void*                                      file_ptr_;
        Bool_t*                                    mat_ptr_;
        std::array<Long64_t, DIM_>                 n_;
        std::array<Double_t, DIM_>                 min_;
        std::array<Double_t, DIM_>                 max_;
        std::array<Double_t, DIM_>                 len_;
        std::array<Double_t, DIM_>                 dlt_;
        std::array<Long64_t, 2>                    fact_;
        std::array<Bool_t, MatProperty::NUM_ELM>   elm_;
        std::array<Double_t, MatProperty::NUM_ELM> den_;
        Double_t                                   inv_rad_len_;
};


#ifdef __HAS_AMS_OFFICE_LIBS__
namespace MatAms {
    constexpr Long64_t DIM = 3;
    
    constexpr Float_t MAGNETIC_RADIUS = 53.0;

    constexpr std::array<Long64_t, DIM> TRL1_N   {  260, 200,     1 };
    constexpr std::array<Double_t, DIM> TRL1_MAX {   65,  50, 159.0 };
    constexpr std::array<Double_t, DIM> TRL1_MIN {  -65, -50, 158.5 };
    
    constexpr std::array<Long64_t, DIM> TRD_N   {  400,  400, 120 };
    constexpr std::array<Double_t, DIM> TRD_MAX {  100,  100, 145 };
    constexpr std::array<Double_t, DIM> TRD_MIN { -100, -100,  85 };
    
    constexpr std::array<Long64_t, DIM> UTOF_N   { 320, 360, 14 };
    constexpr std::array<Double_t, DIM> UTOF_MAX {  80,  90, 67 };
    constexpr std::array<Double_t, DIM> UTOF_MIN { -80, -90, 60 };
    
    constexpr std::array<Long64_t, DIM> TRL2_N   { 260, 160,    1 };
    constexpr std::array<Double_t, DIM> TRL2_MAX {  65,  40, 53.0 };
    constexpr std::array<Double_t, DIM> TRL2_MIN { -65, -40, 52.5 };
    
    constexpr std::array<Long64_t, DIM> TRL3_N   { 240, 240,    1 };
    constexpr std::array<Double_t, DIM> TRL3_MAX {  60,  60, 29.5 };
    constexpr std::array<Double_t, DIM> TRL3_MIN { -60, -60, 29.0 };
    
    constexpr std::array<Long64_t, DIM> TRL4_N   { 240, 240,    1 };
    constexpr std::array<Double_t, DIM> TRL4_MAX {  60,  60, 25.5 };
    constexpr std::array<Double_t, DIM> TRL4_MIN { -60, -60, 25.0 };
    
    constexpr std::array<Long64_t, DIM> TRL5_N   { 240, 240,   1 };
    constexpr std::array<Double_t, DIM> TRL5_MAX {  60,  60, 2.0 };
    constexpr std::array<Double_t, DIM> TRL5_MIN { -60, -60, 1.5 };
    
    constexpr std::array<Long64_t, DIM> TRL6_N   { 240, 240,    1 };
    constexpr std::array<Double_t, DIM> TRL6_MAX {  60,  60, -2.0 };
    constexpr std::array<Double_t, DIM> TRL6_MIN { -60, -60, -2.5 };
    
    constexpr std::array<Long64_t, DIM> TRL7_N   { 240, 240,     1 };
    constexpr std::array<Double_t, DIM> TRL7_MAX {  60,  60, -25.0 };
    constexpr std::array<Double_t, DIM> TRL7_MIN { -60, -60, -25.5 };
    
    constexpr std::array<Long64_t, DIM> TRL8_N   { 240, 240,     1 };
    constexpr std::array<Double_t, DIM> TRL8_MAX {  60,  60, -29.0 };
    constexpr std::array<Double_t, DIM> TRL8_MIN { -60, -60, -29.5 };
    
    constexpr std::array<Long64_t, DIM> LTOF_N   { 320, 360,  14 };
    constexpr std::array<Double_t, DIM> LTOF_MAX {  80,  90, -60 };
    constexpr std::array<Double_t, DIM> LTOF_MIN { -80, -90, -67 };
    
    constexpr std::array<Long64_t, DIM> RICH_N   { 240, 240,   2 };
    constexpr std::array<Double_t, DIM> RICH_MAX {  60,  60, -75 };
    constexpr std::array<Double_t, DIM> RICH_MIN { -60, -60, -76 };
    
    constexpr std::array<Long64_t, DIM> PMT_N   { 280, 280,    8 };
    constexpr std::array<Double_t, DIM> PMT_MAX {  70,  70, -122 };
    constexpr std::array<Double_t, DIM> PMT_MIN { -70, -70, -126 };
    
    constexpr std::array<Long64_t, DIM> TRL9_N   {  200, 120,      1 };
    constexpr std::array<Double_t, DIM> TRL9_MAX {   50,  30, -135.5 };
    constexpr std::array<Double_t, DIM> TRL9_MIN {  -50, -30, -136.0 };
    
    constexpr std::array<Long64_t, DIM> ECAL_N   { 160, 160,   36 };
    constexpr std::array<Double_t, DIM> ECAL_MAX {  40,  40, -142 };
    constexpr std::array<Double_t, DIM> ECAL_MIN { -40, -40, -160 };
}

class MatGeoBoxAms {
    public :
        MatGeoBoxAms() {}
        ~MatGeoBoxAms() {}

        static Bool_t CreateMatGeoBoxFromG4MatTree();

        static Bool_t Load();

    private :
        static Bool_t is_load_;
        static std::vector<MatGeoBoxReader*> reader_;
        static MatGeoBoxReader reader_AMS02TRL1_;
        static MatGeoBoxReader reader_AMS02TRD_;
        static MatGeoBoxReader reader_AMS02UTOF_;
        static MatGeoBoxReader reader_AMS02TRL2_;
        static MatGeoBoxReader reader_AMS02TRL3_;
        static MatGeoBoxReader reader_AMS02TRL4_;
        static MatGeoBoxReader reader_AMS02TRL5_;
        static MatGeoBoxReader reader_AMS02TRL6_;
        static MatGeoBoxReader reader_AMS02TRL7_;
        static MatGeoBoxReader reader_AMS02TRL8_;
        static MatGeoBoxReader reader_AMS02LTOF_;
        static MatGeoBoxReader reader_AMS02RICH_;
        static MatGeoBoxReader reader_AMS02PMT_;
        static MatGeoBoxReader reader_AMS02TRL9_;
        static MatGeoBoxReader reader_AMS02ECAL_;
};

Bool_t MatGeoBoxAms::is_load_ = false;
std::vector<MatGeoBoxReader*> MatGeoBoxAms::reader_;
MatGeoBoxReader MatGeoBoxAms::reader_AMS02TRL1_;
MatGeoBoxReader MatGeoBoxAms::reader_AMS02TRD_;
MatGeoBoxReader MatGeoBoxAms::reader_AMS02UTOF_;
MatGeoBoxReader MatGeoBoxAms::reader_AMS02TRL2_;
MatGeoBoxReader MatGeoBoxAms::reader_AMS02TRL3_;
MatGeoBoxReader MatGeoBoxAms::reader_AMS02TRL4_;
MatGeoBoxReader MatGeoBoxAms::reader_AMS02TRL5_;
MatGeoBoxReader MatGeoBoxAms::reader_AMS02TRL6_;
MatGeoBoxReader MatGeoBoxAms::reader_AMS02TRL7_;
MatGeoBoxReader MatGeoBoxAms::reader_AMS02TRL8_;
MatGeoBoxReader MatGeoBoxAms::reader_AMS02LTOF_;
MatGeoBoxReader MatGeoBoxAms::reader_AMS02RICH_;
MatGeoBoxReader MatGeoBoxAms::reader_AMS02PMT_;
MatGeoBoxReader MatGeoBoxAms::reader_AMS02TRL9_;
MatGeoBoxReader MatGeoBoxAms::reader_AMS02ECAL_;


#endif // __HAS_AMS_OFFICE_LIBS__





























/*

class MatPhyBox : public MatGeoBox {
    public :
        MatPhyBox(PhySt, MatGeoBox)

    protected :
        static Double_t GetInverseRadiationLength();
        static Double_t GetMultipleScatteringLength();
        static Double_t GetMultipleScatteringDirection();
        static Double_t GetDensityEffectCorrection();
        static Double_t GetIonizationEnergyLoss();
        static Double_t GetBremsstrahlungEnergyLoss();
};


class MatPhyPar {
    public :
    private :
        Bool_t   fVacuum;    // vacuum
        Double_t fNumRadLen; // number of radiation length [1]
        Double_t fInvRadLen; // inverse radiation length [cm^-1]
        Double_t fMscatL;    // multiple-scattering length [1]
        Double_t fMscatD;    // multiple-scattering direction [cm^-1]
        Double_t fEnglsISGM; // ionization-energy-loss SGM [cm^-1]
        Double_t fEnglsIMPV; // ionization-energy-loss MPV [cm^-1]
        Double_t fEnglsBMEN; // bremsstrahlung-energy-loss Mean [cm^-1]
};

class MatPhyMgnt {
    public :
        static void LoadAMS02Env();

    private :
        std::set<MayPhyBox> ...;
};
*/


} // namespace TrackSys



#endif // __TRACKLibs_MatEnv_H__

