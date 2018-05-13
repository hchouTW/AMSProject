#ifndef __TRACKLibs_MatEnv_H__
#define __TRACKLibs_MatEnv_H__

#include <sys/mman.h>
#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <TTree.h>


namespace TrackSys {


class MatFld {
    public :
        MatFld() { clear(); }
        MatFld(Double_t rlen) { clear(); rlen_ = rlen; }
        
        MatFld(Bool_t mat, Double_t irl = 0, Double_t eld = 0, Double_t lme = 0, Double_t dec = 0, Double_t rlen = 0, Double_t elen = 0, Double_t loc1 = 0, Double_t loc2 = 0) : mat_(mat), irl_(irl), eld_(eld), lme_(lme), dec_(dec), rlen_(rlen), elen_(elen), loc1_(loc1), loc2_(loc2) {}
        
        ~MatFld() {}

        void print() const;

        inline const Bool_t& operator() () const { return mat_; }
        inline const Double_t& irl() const { return irl_; }
        inline const Double_t& eld() const { return eld_; }
        inline const Double_t& lme() const { return lme_; }
        inline const Double_t& dec() const { return dec_; }
        
        inline const Double_t& rlen() const { return rlen_; }
        inline const Double_t& elen() const { return elen_; }
        inline const Double_t& loc1() const { return loc1_; }
        inline const Double_t& loc2() const { return loc2_; }

        inline Double_t nrl() const { return (mat_ ? (irl_ * elen_) : Numc::ZERO<>); }
        inline Double_t ela() const { return (mat_ ? (eld_ * elen_) : Numc::ZERO<>); }

    protected :
        inline void clear() {
            mat_ = false; irl_ = 0; eld_ = 0; lme_ = 0; dec_ = 0;
            rlen_ = 0; elen_ = 0; loc1_ = 0; loc2_ = 0;
        }

    private :
        Bool_t   mat_;
        Double_t irl_;
        Double_t eld_;
        Double_t lme_;
        Double_t dec_;

        Double_t rlen_;
        Double_t elen_;
        Double_t loc1_;
        Double_t loc2_;

    public :
        static MatFld Merge(const std::list<MatFld>& mflds);
};


//======================================================================//
// MatGeoBox - Data Format                                              //
//----------------------------------------------------------------------//
// Header :                                                             //
// Long64_t n[3]       [1]        n-cell                                //
// Double_t min[3]     [cm]       boundary                              //
// Double_t max[3]     [cm]       boundary                              //
// Double_t stp        [1]        n-cell per step                       //
// Content :                                                            //
// ## Index = ( xi*YN*ZN + yi*ZN + zi )                                 //
// ##         xi[0 XN-1] yi[0 YN-1] zi[0 ZN-1]                          //
// Bool_t   mat[index]                                                  //
// Double_t var[npar*index+elm]                                         //
// elm_0: [1/cm]     inverse radiation length                           //
// elm_1: [mole/cm3] electron density                                   //
// elm_2: [MeV]      log mean excitation energy                         //
// elm_3: [1]        density effect correction C                        //
// elm_4: [1]        density effect correction M                        //
// elm_5: [1]        density effect correction A                        //
// elm_6: [1]        density effect correction X0                       //
// elm_7: [1]        density effect correction X1                       //
//======================================================================//
class G4MatStep {
    public :
        G4MatStep(TTree* tree = nullptr) : x(0), y(0), area(0), nstp(0), min(nullptr), max(nullptr), irl(nullptr), eld(nullptr), lme(nullptr), dcC(nullptr), dcM(nullptr), dcA(nullptr), dcX0(nullptr), dcX1(nullptr) { SetBranchAddress(tree); }
        ~G4MatStep() {}

        void SetBranchAddress(TTree* tree) {
            if (tree == nullptr) return;
            tree->SetBranchAddress("x",    &x);
            tree->SetBranchAddress("y",    &y);
            tree->SetBranchAddress("area", &area);
            tree->SetBranchAddress("nstp", &nstp);
            tree->SetBranchAddress("min",  &min);
            tree->SetBranchAddress("max",  &max);
            tree->SetBranchAddress("irl",  &irl);
            tree->SetBranchAddress("eld",  &eld);
            tree->SetBranchAddress("lme",  &lme);
            
            tree->SetBranchAddress("dcC",  &dcC);
            tree->SetBranchAddress("dcM",  &dcM);
            tree->SetBranchAddress("dcA",  &dcA);
            tree->SetBranchAddress("dcX0", &dcX0);
            tree->SetBranchAddress("dcX1", &dcX1);
        }

    public :
        Double_t x;
        Double_t y;
        Double_t area;
        Long64_t nstp;
        std::vector<Double_t>* min;
        std::vector<Double_t>* max;
        std::vector<Double_t>* irl;
        std::vector<Double_t>* eld;
        std::vector<Double_t>* lme;
        
        std::vector<Double_t>* dcC;
        std::vector<Double_t>* dcM;
        std::vector<Double_t>* dcA;
        std::vector<Double_t>* dcX0;
        std::vector<Double_t>* dcX1;
};


constexpr Long64_t MATGEOBOX_NDIM = 3;
constexpr Long64_t MATGEOBOX_NPAR = 8;
constexpr Long64_t MATVAR_IRL = 0;
constexpr Long64_t MATVAR_ELD = 1;
constexpr Long64_t MATVAR_LME = 2;
constexpr Long64_t MATVAR_C   = 3;
constexpr Long64_t MATVAR_M   = 4;
constexpr Long64_t MATVAR_A   = 5;
constexpr Long64_t MATVAR_X0  = 6;
constexpr Long64_t MATVAR_X1  = 7;

struct MatGeoBoxInf {
    Long64_t  n[3];
    Double_t  min[3];
    Double_t  max[3];
    Double_t  stp;
    Bool_t    mat;
};

struct MatGeoBoxVar {
    Double_t  var;
};


class MatGeoBoxCreator {
    public :
        MatGeoBoxCreator(const std::array<Long64_t, 3>& n, const std::array<Double_t, 3>& min, const std::array<Double_t, 3>& max, Double_t stp = Numc::HALF, const std::string& file_path = "MatGeoBox", const std::string& dir_path = ".");
        ~MatGeoBoxCreator() { save_and_close(); }

        void fill(const G4MatStep& g4mat);

        inline Bool_t is_open() { return is_open_; }
       
        void save_and_close();

    protected :
        inline void clear() { is_open_ = false; fdes_inf_ = -1; flen_inf_ = 0; fptr_inf_ = reinterpret_cast<void*>(-1); gbox_inf_ = nullptr; fdes_var_ = -1; flen_var_ = 0; fptr_var_ = reinterpret_cast<void*>(-1); gbox_var_ = nullptr; mat_ptr_ = nullptr; var_ptr_ = nullptr; max_len_ = 0; area_ = 0; dlt_.fill(0); fact_.fill(0); }

    private :
        Bool_t                               is_open_;
        
        Int_t                                fdes_inf_;
        Int_t                                flen_inf_;
        void*                                fptr_inf_;
        MatGeoBoxInf*                        gbox_inf_;
        
        Int_t                                fdes_var_;
        Int_t                                flen_var_;
        void*                                fptr_var_;
        MatGeoBoxVar*                        gbox_var_;
        
        Bool_t*                              mat_ptr_;
        Double_t*                            var_ptr_;
        Long64_t                             max_len_;
        Double_t                             area_;
        std::array<Double_t, MATGEOBOX_NDIM> dlt_;
        std::array<Long64_t, 2>              fact_;
};



class MatGeoBoxReader {
    public :
        MatGeoBoxReader() { clear(); }
        MatGeoBoxReader(const std::string& fname, const std::string& dpath) { clear(); load(fname, dpath); }
        ~MatGeoBoxReader() { clear(); }

        inline const Bool_t& exist() const { return is_load_; }

        void print() const;

        Bool_t load(const std::string& fname, const std::string& dpath);
        
        inline Bool_t is_in_box(const SVecD<3>& coo);

        inline Bool_t is_cross(const SVecD<3>& vcoo, const SVecD<3>& wcoo);
        
        MatFld get(const SVecD<3>& coo, Double_t log10gb = -10);
        MatFld get(const SVecD<3>& vcoo, const SVecD<3>& wcoo, Double_t log10gb = -10, Bool_t is_std = true);
        
    protected :
        inline void clear() { is_load_ = false; max_len_ = 0; n_.fill(0); min_.fill(0.); max_.fill(0.); len_.fill(0.); dlt_.fill(0.); fact_.fill(0); stp_ = 0; tmp_dec_.first = -1; tmp_dec_.second = 0.; mat_.clear(); var_.clear(); }
        
        inline Double_t get_density_effect_correction(Long64_t idx = -1, Double_t log10gb = -10);

    private :
        static constexpr Double_t STD_STEP_LEN = Numc::ONE<>;
        static constexpr Double_t FST_STEP_LEN = Numc::ONE<> + Numc::HALF;

        Bool_t                               is_load_;
        Long64_t                             max_len_;
        std::array<Long64_t, MATGEOBOX_NDIM> n_;
        std::array<Double_t, MATGEOBOX_NDIM> min_;
        std::array<Double_t, MATGEOBOX_NDIM> max_;
        std::array<Double_t, MATGEOBOX_NDIM> len_;
        std::array<Double_t, MATGEOBOX_NDIM> dlt_;
        std::array<Long64_t, 2>              fact_;
        Double_t                             stp_;

        std::vector<Bool_t>                               mat_;
        std::vector<std::array<Double_t, MATGEOBOX_NPAR>> var_;

    private :
        std::pair<Long64_t, Double_t> tmp_dec_; // for speed up
        inline void reset_tmp_dec() { tmp_dec_.first = -1; tmp_dec_.second = 0.; }
};


class MatMgnt {
    public :
        MatMgnt() {}
        ~MatMgnt() {}

        static Bool_t Load();

        static MatFld Get(const SVecD<3>& coo, Double_t log10gb = -10);
        static MatFld Get(const SVecD<3>& vcoo, const SVecD<3>& wcoo, Double_t log10gb = -10, Bool_t is_std = true);
        
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
        MatPhyFld(Bool_t mat, Double_t mscat_sgm = 0, Double_t elion_mpv = 0, Double_t elion_sgm = 0, Double_t elion_men = 0, Double_t elbrm_men = 0) : mat_(mat), mscat_sgm_(mscat_sgm), elion_mpv_(elion_mpv), elion_sgm_(elion_sgm), elion_men_(elion_men), elbrm_men_(elbrm_men) {}
        ~MatPhyFld() {}

        inline const Bool_t& operator() () const { return mat_; }
        inline const Double_t& mscat_sgm() const { return mscat_sgm_; }
        inline const Double_t& elion_mpv() const { return elion_mpv_; }
        inline const Double_t& elion_sgm() const { return elion_sgm_; }
        inline const Double_t& elion_men() const { return elion_men_; }
        inline const Double_t& elbrm_men() const { return elbrm_men_; }
        
    protected :
        inline void clear() { mat_ = false; mscat_sgm_ = 0; elion_mpv_ = 0; elion_sgm_ = 0; elion_men_ = 0; elbrm_men_ = 0; }

    private :
        Bool_t   mat_;       // has matter?
        Double_t mscat_sgm_; // multiple-scattering [1]
        Double_t elion_mpv_; // ionization-energy-loss MPV [1]
        Double_t elion_sgm_; // ionization-energy-loss SGM [1]
        Double_t elion_men_; // ionization-energy-loss MEN [1]
        Double_t elbrm_men_; // bremsstrahlung-energy-loss MEN [1]
};


class MatPhy {
    public :
        MatPhy() {}
        ~MatPhy() {}

        static MatPhyFld Get(const Double_t stp_len, PhySt& part, Bool_t is_std = true);
        static MatPhyFld Get(const MatFld& mfld, PhySt& part);

    protected :
        static Double_t GetMultipleScattering(const MatFld& mfld, PhySt& part);
        static std::tuple<Double_t, Double_t, Double_t> GetIonizationEnergyLoss(const MatFld& mfld, PhySt& part);
        static Double_t GetBremsstrahlungEnergyLoss(const MatFld& mfld, PhySt& part);
    
    // Expert only
    public :
        static void SetCorrFactor(const MatFld* mfld = nullptr, PhySt* part = nullptr, Bool_t sw_mscat = true, Bool_t sw_eloss = true) {
            Bool_t sw = ((mfld != nullptr) && (*mfld)()) && (sw_mscat || sw_eloss); 
            if (sw) { corr_sw_mscat_ = sw_mscat; corr_sw_eloss_ = sw_eloss; corr_use_elion_mpv_ = false; corr_mfld_ = *mfld; }
            else    { corr_sw_mscat_ = false;    corr_sw_eloss_ = false;    corr_use_elion_mpv_ = false; corr_mfld_ = std::move(MatFld()); }
            if (corr_sw_eloss_ && corr_mfld_() && part != nullptr) {
                MatPhyFld&& mpfld = MatPhy::Get(corr_mfld_, *part);
                corr_use_elion_mpv_ = (Numc::Compare(mpfld.elion_mpv(), mpfld.elion_men()) < 0);
            }
        }

        static const Bool_t& UseElionMpv() { return corr_use_elion_mpv_; }

    private :
        static Bool_t corr_sw_mscat_;
        static Bool_t corr_sw_eloss_;
        static Bool_t corr_use_elion_mpv_;
        static MatFld corr_mfld_;

    private :
        // Coulomb Multiple Scattering, the Highland-Lynch-Dahl equation
        // Sigma_plane_angle = (RydbergConstant / abs(beta * rigidity) *
        //                     sqrt( radiationLength ) *
        //                     (1. + 0.038 * log(radiationLength)) )
        static constexpr Double_t RYDBERG_CONST  = 0.0136; // [GeV]
        static constexpr Double_t NRL_CORR_FACT  = 0.0380; // [1]
        static constexpr Double_t NRL_CORR_FACT1 = 0.1050; // [1]
        static constexpr Double_t NRL_CORR_FACT2 = 0.0035; // [1]
        
        // Energy Loss from ionization, the Bethe-Bloch equation
        static constexpr Double_t BETHE_BLOCH_K = 0.307075; // [MeV mol^-1 cm^2]
        static constexpr Double_t LANDAU_ELOSS_CORR = 0.2;
        static constexpr Double_t MASS_EL_IN_MEV = 0.510999; // [MeV]
        static constexpr Double_t MASS_EL_IN_GEV = 0.000510999; // [GeV]
        
        // Beta Limit (0.2)
        static constexpr Double_t LMT_BTA           = 0.2;
        static constexpr Double_t LMT_INV_BTA       = 5.0;
        static constexpr Double_t LMT_SQR_BTA       = 0.04;
        static constexpr Double_t LMT_GMBTA         = 2.041241e-01;
        static constexpr Double_t LMT_INV_GMBTA     = 4.898981e+00;
        static constexpr Double_t LMT_SQR_GMBTA     = 4.166666e-02;
       
        // Unit
        static constexpr Double_t MEV_TO_GEV = 1.0e-3; // (MeV/GeV)
        static constexpr Double_t GEV_TO_MEV = 1.0e+3; // (GeV/MeV)
};
        

Bool_t   MatPhy::corr_sw_mscat_      = false;
Bool_t   MatPhy::corr_sw_eloss_      = false;
Bool_t   MatPhy::corr_use_elion_mpv_ = false;
MatFld   MatPhy::corr_mfld_;


} // namespace TrackSys


#endif // __TRACKLibs_MatEnv_H__
