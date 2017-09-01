#ifndef __TRACKLibs_HitSt_C__
#define __TRACKLibs_HitSt_C__


namespace TrackSys {
        

void HitSt::print() const {
    std::string printStr;
    printStr += STR_FMT("================= HitSt ==================\n");
    printStr += STR_FMT("Id   (%d)\n", id_);
    printStr += STR_FMT("Side (%d %d %d)\n", side_(0), side_(0), side_(0));
    printStr += STR_FMT("Coo  (%11.6f %11.6f %11.6f)\n", coo_(0), coo_(1), coo_(2));
    printStr += STR_FMT("Err  (%11.6f %11.6f %11.6f)\n", err_x_.efft_sgm(), err_y_.efft_sgm(), err_z_.efft_sgm());
    printStr += STR_FMT("==========================================\n");
    COUT(printStr);
}


} // namesapce TrackSys


#endif // __TRACKLibs_HitSt_C__
