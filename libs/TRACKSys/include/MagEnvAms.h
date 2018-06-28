#if defined(_PGTRACK_) || defined(__ROOTSHAREDLIBRARY__)
#ifndef __TRACKLibs_MagEnvAms_H__
#define __TRACKLibs_MagEnvAms_H__

#include <MagField.h>


namespace TrackSys {


class MagGeoBoxAms {
    public :
        MagGeoBoxAms() {}
        ~MagGeoBoxAms() {}

        static Bool_t Load();

        static MagFld Get(const SVecD<3>& coo);

        static void Output(const std::string& file_path = "AMS02Mag.bin");

    private :
        static Bool_t    is_load_;
        static MagField* mag_field_;
};

} // namespace TrackSys


#endif // __TRACKLibs_MagEnvAms_H__
#endif // _PGTRACK_ __ROOTSHAREDLIBRARY__ 
