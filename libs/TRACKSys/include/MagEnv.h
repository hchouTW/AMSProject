#ifndef __TRACKLibs_MagEnv_H__
#define __TRACKLibs_MagEnv_H__


namespace TrackSys {


class MagFld {
    public :
        MagFld() {}
        MagFld(const SVecD<3>& mag) : mag_(mag) {}
        MagFld(Double_t x, Double_t y, Double_t z) : mag_(x, y, z) {}
        MagFld(Float_t mag[3]) : mag_(mag[0], mag[1], mag[2]) {}
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
// Float_t  MIN[3]                                                      //
// Float_t  MAX[3]                                                      //
// Content :                                                            //
// Float_t MAG[(index)*3 + di]                                          //
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
    Float_t  min[3];
    Float_t  max[3];
    Float_t  mag;
};


class MagGeoBoxCreator {
    public :
        MagGeoBoxCreator(const Long64_t n[3], const Float_t min[3], const Float_t max[3], const std::string& fpath = "MagGeoBox.bin");
        ~MagGeoBoxCreator() { if (is_open_) save_and_close(); }

        void fill(Long64_t idx, Float_t bx = 0., Float_t by = 0., Float_t bz = 0.);

        inline Bool_t is_open() { return is_open_; }
        void save_and_close();
        
        void save_and_close(Float_t bx, Float_t by, Float_t bz);

    protected :
        inline void clear() { is_open_ = false; fdes_ = -1; flen_ = 0; fptr_ = reinterpret_cast<void*>(-1); max_len_ = 0; geo_box_ = nullptr; }

    private :
        static constexpr Long64_t DIM_ = 3;
        Bool_t                    is_open_;
        Int_t                     fdes_;
        Int_t                     flen_;
        void*                     fptr_;
        Long64_t                  max_len_;
        MagGeoBox*                geo_box_;
};


class MagGeoBoxReader {
    public :
        MagGeoBoxReader() { clear(); }
        MagGeoBoxReader(const std::string& fpath) { clear(); load(fpath); }
        ~MagGeoBoxReader() { clear(); }

        inline Bool_t exist() { return is_load_; }
        inline MagFld get(const SVecD<3>& coo);
        
        Bool_t load(const std::string& fpath);

    protected :
        inline void clear() { is_load_ = false; fpath_ = ""; fptr_ = reinterpret_cast<void*>(-1); mag_ptr_ = nullptr;  n_.fill(0); min_.fill(0.); max_.fill(0.); dlt_.fill(0.); fact_.fill(0); }

    private :
        static constexpr Long64_t   DIM_ = 3;
        Bool_t                      is_load_;
        std::string                 fpath_;
        void*                       fptr_;
        Float_t*                    mag_ptr_;
        std::array<Long64_t, DIM_>  n_;
        std::array<Float_t,  DIM_>  min_;
        std::array<Float_t,  DIM_>  max_;
        std::array<Float_t,  DIM_>  dlt_;
        std::array<Long64_t, 2>     fact_;
};


class MagMgnt {
    public :
        static Bool_t Load();
        inline static MagFld Get(const SVecD<3>& coo) { if (!Load()) return MagFld(); return geo_box_reader_.get(coo); }

    protected :
        static Bool_t is_load_;
        static MagGeoBoxReader geo_box_reader_;
};

Bool_t          MagMgnt::is_load_ = false;
MagGeoBoxReader MagMgnt::geo_box_reader_;


} // namespace TrackSys


#endif // __TRACKLibs_MagEnv_H__

