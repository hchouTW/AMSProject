#ifndef __TRACKLibs_MagEnv_H__
#define __TRACKLibs_MagEnv_H__

namespace TrackSys {

class MagFld {
    public :
        MagFld() {}
        MagFld(const SVecD<3>& mag) : mag_(mag) {}
        MagFld(Float_t mag[3]) : mag_(mag[0], mag[1], mag[2]) {}
        MagFld(Double_t mag[3]) : mag_(mag[0], mag[1], mag[2]) {}
        MagFld(Double_t x, Double_t y, Double_t z) : mag_(x, y, z) {}
        ~MagFld() {}

        inline const SVecD<3>& operator() () const { return mag_; }
        inline const Double_t& x() const { return mag_(0); }
        inline const Double_t& y() const { return mag_(1); }
        inline const Double_t& z() const { return mag_(2); }

    private :
        SVecD<3> mag_;
};


//======================================================================//
// MagGeoBox - Data Format                                              //
//----------------------------------------------------------------------//
// Header :                                                             //
// Long64_t N[3]                                                        //
// Double_t MIN[3]                                                      //
// Double_t MAX[3]                                                      //
// Content :                                                            //
// Double_t MAG[(index)*3 + di]                                         //
// ## Index = ( xi*YN*ZN + yi*ZN + zi )                                 //
// ##         xi[0 XN-1] yi[0 YN-1] zi[0 ZN-1] di[0 2]                  //
//======================================================================//
#include <sys/mman.h>
#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>


struct MagGeoBox {
    Long64_t n[3];
    Double_t min[3];
    Double_t max[3];
    Double_t mag;
};


class MagGeoBoxCreator {
    public :
        MagGeoBoxCreator(Long64_t xn, Double_t xmin, Double_t xmax, Long64_t yn, Double_t ymin, Double_t ymax, Long64_t zn, Double_t zmin, Double_t zmax, const std::string& file_path = "MagGeoBox.bin");
        ~MagGeoBoxCreator() { save_and_close(); }

        void fill(Double_t mx = 0., Double_t my = 0., Double_t mz = 0.); // from min to max coord

        inline Bool_t is_open() { return is_open_; }
        void save_and_close();

    protected :
        inline void clear() { is_open_ = false; file_path_ = ""; file_des_ = -1; file_len_ = 0; file_ptr_ = reinterpret_cast<void*>(-1); max_len_ = 0; cur_len_= 0; geo_box_ = nullptr; }

    private :
        static constexpr Long64_t DIM_ = 3;
        Bool_t                    is_open_;
        std::string               file_path_;
        Int_t                     file_des_;
        Int_t                     file_len_;
        void*                     file_ptr_;
        Long64_t                  max_len_;
        Long64_t                  cur_len_;
        MagGeoBox*                geo_box_;
};


class MagGeoBoxReader {
    public :
        MagGeoBoxReader() { clear(); }
        MagGeoBoxReader(const std::string& file_path) { clear(); load(file_path); }
        ~MagGeoBoxReader() { clear(); }

        inline Bool_t exist() { return is_load_; }
        inline MagFld get(const SVecD<3>& coo) { return ((!is_load_) ? MagFld() : MagFld( do_trilinear_interpolation( get_index(coo) ) )); }
        
        Bool_t load(const std::string& file_path);

    protected :
        inline void clear() { is_load_ = false; file_path_ = ""; file_ptr_ = reinterpret_cast<void*>(-1); mag_ptr_ = nullptr;  n_.fill(0); min_.fill(0.); max_.fill(0.); dlt_.fill(0.); fact_.fill(0); }

        using Index = std::tuple<Long64_t, Double_t, Double_t, Double_t>;
        Index get_index(const SVecD<3>& coo);
        SVecD<3> do_trilinear_interpolation(const Index& index);

    private :
        static constexpr Long64_t DIM_ = 3;
        Bool_t                    is_load_;
        std::string               file_path_;
        void*                     file_ptr_;
        Double_t*                 mag_ptr_;
        std::array<Long64_t, 3>   n_;
        std::array<Double_t, 3>   min_;
        std::array<Double_t, 3>   max_;
        std::array<Double_t, 3>   dlt_;
        std::array<Long64_t, 2>   fact_;
};


#ifdef __HAS_AMS_OFFICE_LIBS__
class MagGeoBoxAms {
    public :
        MagGeoBoxAms() {}
        ~MagGeoBoxAms() {}

        static Bool_t Load();

        inline static MagFld Get(const SVecD<3>& coo);

        static void Output(const std::string& file_path = "MagGeoBox.bin");

    private :
        static Bool_t    is_load_;
        static MagField* mag_field_;
};

Bool_t MagGeoBoxAms::is_load_ = false;
MagField* MagGeoBoxAms::mag_field_ = nullptr;


class MagFuncAms {
    public :
        MagFuncAms() {}
        ~MagFuncAms() {}

        inline static MagFld Get(const SVecD<3>& coo) { return MagFld(GetMagx(coo[2]), MGMath::ZERO, MGMath::ZERO); }

        inline static Double_t GetMagx(Double_t cooz);
        inline static Double_t GetMagxInt1st(Double_t cooz);
        inline static Double_t GetMagxInt2nd(Double_t cooz);

    private :
        static const std::array<Double_t, 2> PAR_MAG;
        static const std::array<Double_t, 2> PAR_SGM;
};

const std::array<Double_t, 2> MagFuncAms::PAR_MAG = { 123.21342990, 20.80309443 };
const std::array<Double_t, 2> MagFuncAms::PAR_SGM = {  37.01735081, 79.25552614 };
#endif // __HAS_AMS_OFFICE_LIBS__


enum class MagType {
#ifdef __HAS_AMS_OFFICE_LIBS__
    kGeoBoxAms,
    kFuncAms,
#endif // __HAS_AMS_OFFICE_LIBS__
    kGeoBox
};


class MagMgnt {
    public :
        static Bool_t Load();

        static MagFld Get(const SVecD<3>& coo,
#ifdef __HAS_AMS_OFFICE_LIBS__
                MagType type = MagType::kGeoBoxAms
#elif
                MagType type = MagType::kGeoBox
#endif // __HAS_AMS_OFFICE_LIBS__
            );

#ifdef __HAS_AMS_OFFICE_LIBS__
        static Double_t GetMagx(Double_t cooz) { return MagFuncAms::GetMagx(cooz); }
        static Double_t GetMagxInt1st(Double_t cooz) { return MagFuncAms::GetMagxInt1st(cooz); }
        static Double_t GetMagxInt2nd(Double_t cooz) { return MagFuncAms::GetMagxInt2nd(cooz); }
#endif // __HAS_AMS_OFFICE_LIBS__

    protected :
        static Bool_t is_load_;
        static MagGeoBoxReader geo_box_reader_;
};

Bool_t          MagMgnt::is_load_ = false;
MagGeoBoxReader MagMgnt::geo_box_reader_;


} // namespace TrackSys


#endif // __TRACKLibs_MagEnv_H__

