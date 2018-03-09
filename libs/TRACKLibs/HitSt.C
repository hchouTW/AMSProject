#ifndef __TRACKLibs_HitSt_C__
#define __TRACKLibs_HitSt_C__


namespace TrackSys {
        

void HitSt::print() const {
    std::string printStr;
    printStr += STR("================= HitSt ==================\n");
    printStr += STR("Lay  (%d)\n", lay_);
    printStr += STR("Side (%d %d %d)\n", side_(0), side_(1));
    printStr += STR("Coo  (%11.6f %11.6f %11.6f)\n", coo_(0), coo_(1), coo_(2));
    printStr += STR("Err  (%11.6f %11.6f)\n", err_(0), err_(1));
    printStr += STR("==========================================\n");
    COUT(printStr.c_str());
}


} // namesapce TrackSys


#endif // __TRACKLibs_HitSt_C__
