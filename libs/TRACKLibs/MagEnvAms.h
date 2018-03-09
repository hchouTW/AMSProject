#ifndef __TRACKLibs_MagEnvAms_H__
#define __TRACKLibs_MagEnvAms_H__


#ifdef __HAS_AMS_OFFICE_LIBS__
#include <MagField.h>

namespace TrackSys {


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

} // namespace TrackSys


#endif // __TRACKLibs_MagEnvAms_H__

