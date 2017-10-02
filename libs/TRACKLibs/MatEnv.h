#ifndef __TRACKLibs_MatEnv_H__
#define __TRACKLibs_MatEnv_H__

namespace TrackSys {

namespace MatProperty {
    constexpr Int_t NUM_ELM = 9;
    const     std::array<const std::string, NUM_ELM> NAME { "Hydrogen(H)", "Carbon(C)", "Nitrogen(N)", "Oxygen(O)", "Fluorine(F)", "Sodium(Na)", "Aluminum(Al)", "Silicon(Si)", "Lead(Pb)" };
    constexpr std::array<Int_t, NUM_ELM>             CHRG { 1, 6, 7, 8, 9, 11, 13, 14, 82 };
    constexpr std::array<Double_t, NUM_ELM>          MASS { 1.007947, 12.01078, 14.00672, 15.99943, 18.99840325, 22.989769282, 26.98153868, 28.08553, 207.21 }; // [g mol^-1]

    constexpr std::array<Double_t, NUM_ELM> RAD_LEN { 63.04, 42.70, 37.99, 34.24, 32.93, 27.74, 24.01, 21.82, 6.37 }; // [g cm^-2]

    // Mean Excitation Energy I = 16eV * Z^0.9 [eV]
    constexpr std::array<Double_t, NUM_ELM> MEAN_EXENG        {  1.92e-05,  8.10e-05,  8.20e-05,  9.50e-05,  1.15e-04,  1.49e-04,  1.66e-04,  1.73e-04,  8.23e-04 }; // [MeV] From NIST, ESTART
    constexpr std::array<Double_t, NUM_ELM> NEG_LN_MEAN_EXENG { 1.086e+01, 9.421e+00, 9.409e+00, 9.262e+00, 9.071e+00, 8.812e+00, 8.704e+00, 8.662e+00, 7.103e+00 }; // [MeV] From NIST, ESTART

    // Density Effect Correction
    constexpr std::array<Double_t, NUM_ELM> DEN_EFF_CORR_C    { 3.2632,  2.9925, 10.5400, 10.7004, 10.9653, 5.0526, 4.2395, 4.4351, 6.2018 }; // Form Geant4
    constexpr std::array<Double_t, NUM_ELM> DEN_EFF_CORR_X0   { 0.4759, -0.0351,  1.7378,  1.7541,  1.8433, 0.2880, 0.1708, 0.2014, 0.3776 }; // Form Geant4
    constexpr std::array<Double_t, NUM_ELM> DEN_EFF_CORR_X1   { 1.9215,  2.4860,  4.1323,  4.3213,  4.4096, 3.1962, 3.0127, 2.8715, 3.8073 }; // Form Geant4
    constexpr std::array<Double_t, NUM_ELM> DEN_EFF_CORR_A    { 0.1348,  0.2024,  0.1535,  0.1178,  0.1108, 0.0777, 0.0802, 0.1492, 0.0936 }; // Form Geant4
    constexpr std::array<Double_t, NUM_ELM> DEN_EFF_CORR_K    { 5.6249,  3.0036,  3.2125,  3.2913,  3.2962, 3.6452, 3.6345, 3.2546, 3.1608 }; // Form Geant4
    constexpr std::array<Double_t, NUM_ELM> DEN_EFF_CORR_DLT0 {   0.13,     0.1,    0.19,    0.11,    0.11,   0.08,   0.12,   0.14,   0.14 }; // Form Geant4
    constexpr std::array<Double_t, NUM_ELM> DEN_EFF_CORR_DLTM {  0.021,   0.038,  0.086,  0.101,     0.121,  0.098,  0.061,  0.059,  0.019 }; // Form Geant4
}


class MatFld {
    public :
        MatFld() { clear(); }
        MatFld(Double_t real_len) { clear(); real_len_ = real_len; }
        MatFld(Bool_t mat, const std::array<Bool_t, MatProperty::NUM_ELM>& elm, const std::array<Double_t, MatProperty::NUM_ELM>& den, Double_t inv_rad_len = 0.0, Double_t real_len = 0.0, Double_t efft_len = 0.0, Double_t efft = 0.0, Double_t loc = 0.0, Double_t locsqr = 0.0) : mat_(mat), elm_(elm), den_(den), inv_rad_len_(inv_rad_len), real_len_(real_len), efft_len_(efft_len), efft_(efft), loc_(loc), locsqr_(locsqr) {}
        ~MatFld() {}

        void print() const;

        inline const Bool_t& operator() () const { return mat_; }
        inline const std::array<Bool_t, MatProperty::NUM_ELM>&   elm() const { return elm_; }
        inline const std::array<Double_t, MatProperty::NUM_ELM>& den() const { return den_; }
        inline const Double_t& inv_rad_len() const { return inv_rad_len_; }
        inline const Double_t& real_len() const { return real_len_; }
        inline const Double_t& efft_len() const { return efft_len_; }
        inline const Double_t& efft() const { return efft_; }
        inline const Double_t& loc() const { return loc_; }
        inline const Double_t& locsqr() const { return locsqr_; }

        inline Double_t num_rad_len() const { return (mat_ ? (inv_rad_len_ * efft_len_) : 0.); }

    protected :
        inline void clear() { mat_ = false; elm_.fill(false); den_.fill(0.); inv_rad_len_ = 0.; real_len_ = 0.; efft_len_ = 0.; efft_ = 0.; loc_ = 0.; locsqr_ = 0.; }

    private :
        Bool_t                                     mat_;
        std::array<Bool_t, MatProperty::NUM_ELM>   elm_;
        std::array<Double_t, MatProperty::NUM_ELM> den_;
        Double_t                                   inv_rad_len_;
        Double_t                                   real_len_;
        Double_t                                   efft_len_;
        Double_t                                   efft_;
        Double_t                                   loc_;
        Double_t                                   locsqr_;

    public :
        static MatFld Merge(const std::list<MatFld>& mflds);
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
// Double_t INV_RAD_LEN                                                 //
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

        void fill(Float_t coo[3], Bool_t elm[MatProperty::NUM_ELM], Float_t mol[MatProperty::NUM_ELM], Bool_t calculated = true);

        inline Bool_t is_open() { return is_open_; }
       
        void save_and_close();

        void save_and_close(Bool_t elm[MatProperty::NUM_ELM], Float_t den[MatProperty::NUM_ELM]);

    protected :
        inline void clear() { is_open_ = false; file_path_ = ""; file_des_ = -1; file_len_ = 0; file_ptr_ = reinterpret_cast<void*>(-1); max_len_ = 0; geo_box_ = nullptr; vol_ = 0.; dlt_.fill(0); fact_.fill(0); cnt_ = 0; elm_.fill(false); den_.fill(0); inv_rad_len_ = 0.; }

    private :
        static constexpr Long64_t  DIM_ = 3;
        Bool_t                     is_open_;
        std::string                file_path_;
        Int_t                      file_des_;
        Int_t                      file_len_;
        void*                      file_ptr_;
        Long64_t                   max_len_;
        MatGeoBox*                 geo_box_;
        Double_t                   vol_;
        std::array<Double_t, DIM_> dlt_;
        std::array<Long64_t, 2>    fact_;
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

        void print() const;

        Bool_t load(const std::string& file_path);
        
        inline Bool_t is_in_box(const SVecD<3>& coo);

        inline Bool_t is_cross(const SVecD<3>& vcoo, const SVecD<3>& wcoo);
        
        MatFld get(const SVecD<3>& coo);
        MatFld get(const SVecD<3>& vcoo, const SVecD<3>& wcoo, Bool_t is_std = true);
        
    protected :
        inline void clear() { is_load_ = false; file_path_ = ""; file_ptr_ = reinterpret_cast<void*>(-1); mat_ptr_ = nullptr; max_len_ = 0; n_.fill(0); min_.fill(0.); max_.fill(0.); len_.fill(0.); dlt_.fill(0.); fact_.fill(0); elm_.fill(false); den_.fill(0.); inv_rad_len_ = 0.; }

    private :
        static constexpr Long64_t                  DIM_ = 3;
        static constexpr Double_t                  STD_STEP_LEN_ = MGMath::ONE_TO_THREE;
        static constexpr Double_t                  FST_STEP_LEN_ = MGMath::ONE_TO_THREE + MGMath::ONE;
        Bool_t                                     is_load_;
        std::string                                file_path_;
        void*                                      file_ptr_;
        Bool_t*                                    mat_ptr_;
        Long64_t                                   max_len_;
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


class MatMgnt {
    public :
        MatMgnt() {}
        ~MatMgnt() {}

        static Bool_t Load();

        static MatFld Get(const SVecD<3>& coo);
        static MatFld Get(const SVecD<3>& vcoo, const SVecD<3>& wcoo, Bool_t is_std = true);
        
        static MatFld Get(Double_t stp_len, const PhySt& part, Bool_t is_std = true);

    protected :
        static Bool_t is_load_;
        static std::list<MatGeoBoxReader*> * reader_;
};

Bool_t MatMgnt::is_load_ = false;
std::list<MatGeoBoxReader*> * MatMgnt::reader_ = nullptr;


class MatPhyFld {
    public :
        MatPhyFld() { clear(); }
        MatPhyFld(Double_t len) { clear(); len_ = len; }
        MatPhyFld(Bool_t mat, Double_t len = 0., Double_t efft = 0, Double_t loc = 0., Double_t locsqr = 0., Double_t inv_rad_len = 0., Double_t num_rad_len = 0., Double_t mult_scat_sgm = 0., Double_t eloss_ion_kpa = 0., Double_t eloss_ion_mpv = 0., Double_t eloss_ion_sgm = 0., Double_t eloss_brm_men = 0.) : mat_(mat), len_(len), efft_(efft), loc_(loc), locsqr_(locsqr), inv_rad_len_(inv_rad_len), num_rad_len_(num_rad_len), mult_scat_sgm_(mult_scat_sgm), eloss_ion_kpa_(eloss_ion_kpa), eloss_ion_mpv_(eloss_ion_mpv), eloss_ion_sgm_(eloss_ion_sgm), eloss_brm_men_(eloss_brm_men) {}
        ~MatPhyFld() {}

        inline const Bool_t& operator() () const { return mat_; }
        inline const Double_t& len() const { return len_; }
        inline const Double_t& efft() const { return efft_; }
        inline const Double_t& loc() const { return loc_; }
        inline const Double_t& locsqr() const { return locsqr_; }
        inline const Double_t& inv_rad_len() const { return inv_rad_len_; }
        inline const Double_t& num_rad_len() const { return num_rad_len_; }
        inline const Double_t& mult_scat_sgm() const { return mult_scat_sgm_; }
        inline const Double_t& eloss_ion_kpa() const { return eloss_ion_kpa_; }
        inline const Double_t& eloss_ion_mpv() const { return eloss_ion_mpv_; }
        inline const Double_t& eloss_ion_sgm() const { return eloss_ion_sgm_; }
        inline const Double_t& eloss_brm_men() const { return eloss_brm_men_; }

    protected :
        inline void clear() { mat_ = false; len_ = 0.; efft_ - 0.; loc_ = 0.; locsqr_ = 0.; inv_rad_len_ = 0.; num_rad_len_ = 0.; mult_scat_sgm_ = 0.; eloss_ion_kpa_ = 0.; eloss_ion_mpv_ = 0.; eloss_ion_sgm_ = 0.; eloss_brm_men_ = 0.; }

    private :
        Bool_t   mat_;           // has matter?
        Double_t len_;
        Double_t efft_;
        Double_t loc_;
        Double_t locsqr_;
        Double_t inv_rad_len_;   // inverse radiation length [cm^-1]
        Double_t num_rad_len_;   // number of radiation length [1]
        Double_t mult_scat_sgm_; // multiple-scattering length [1]
        Double_t eloss_ion_kpa_; // ionization-energy-loss KPA [1]
        Double_t eloss_ion_mpv_; // ionization-energy-loss MPV [1]
        Double_t eloss_ion_sgm_; // ionization-energy-loss SGM [1]
        Double_t eloss_brm_men_; // bremsstrahlung-energy-loss Mean [1]
};


class MatPhy {
    public :
        MatPhy() {}
        ~MatPhy() {}

        static Double_t GetNumRadLen(const Double_t stp_len, const PhySt& part, Bool_t is_std = true);

        static MatPhyFld Get(const Double_t stp_len, const PhySt& part, Bool_t is_std = true);
        
        static MatPhyFld Get(const MatFld& mfld, const PhySt& part);

    protected :
        static std::array<Double_t, MatProperty::NUM_ELM> GetDensityEffectCorrection(const MatFld& mfld, const PhySt& part);
        
        static Double_t GetRadiationLength(const MatFld& mfld, const PhySt& part);
        static Double_t GetMultipleScattering(const MatFld& mfld, const PhySt& part);
        static std::tuple<Double_t, Double_t, Double_t>  GetIonizationEnergyLoss(const MatFld& mfld, const PhySt& part);
        static Double_t GetBremsstrahlungEnergyLoss(const MatFld& mfld, const PhySt& part);

    private :
        // Coulomb Multiple Scattering, the Highland-Lynch-Dahl equation
        // Sigma_plane_angle = (RydbergConstant / abs(beta * rigidity) *
        //                     sqrt( radiationLength ) *
        //                     (1. + 0.038 * log(radiationLength)) )
        static constexpr Double_t RYDBERG_CONST = 0.0136; // [GeV]
        static constexpr Double_t NRL_CORR_FACT = 0.0380; // [1]
        static constexpr Double_t LMTL_NUM_RAD_LEN = 1.0e-3;
        static constexpr Double_t LMTU_NUM_RAD_LEN = 100.;
        
        // Energy Loss from ionization, the Bethe-Bloch equation
        static constexpr Double_t BETHE_BLOCH_K = 0.307075; // [MeV mol^-1 cm^2]
        static constexpr Double_t LANDAU_ELOSS_CORR = 0.2;
        static constexpr Double_t MASS_EL_IN_MEV = 0.510999; // [MeV]
        static constexpr Double_t MASS_EL_IN_GEV = 0.000510999; // [GeV]
        
        // Beta Limit (0.2)
        static constexpr Double_t LMT_BTA           = 0.2;
        static constexpr Double_t LMT_SQR_BTA       = 0.04;
        static constexpr Double_t LMT_INV_SQR_BTA   = 2.500000e+01;
        static constexpr Double_t LMT_GMBTA         = 2.041241e-01;
        static constexpr Double_t LMT_SQR_GMBTA     = 4.166666e-02;
        static constexpr Double_t LMT_INV_GMBTA     = 4.898979e+00;
        static constexpr Double_t LMT_INV_SQR_GMBTA = 2.400000e+01;
        static constexpr Double_t LMT_GM            = 1.020620e+00;
       
        // Unit
        static constexpr Double_t MEV_TO_GEV = 1.0e-3;
        static constexpr Double_t GEV_TO_MEV = 1.0e+3;
};


} // namespace TrackSys


#endif // __TRACKLibs_MatEnv_H__
