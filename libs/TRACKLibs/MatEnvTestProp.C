#ifdef __HAS_TESTPROP__


#ifndef __TRACKLibs_MatEnvTestProp_C__
#define __TRACKLibs_MatEnvTestProp_C__


namespace TrackSys {


// Set to MatMgnt::Load()
Bool_t MatMgnt::Load() {
    if (is_load_ && reader_ != nullptr) return true;
    is_load_ = false;
    reader_ = nullptr;

    is_load_ = MatGeoBoxTestProp::Load();
    if (is_load_) reader_ = &MatGeoBoxTestProp::Reader();
    
    return is_load_;
};


Bool_t MatGeoBoxTestProp::CreateMatGeoBox() {
    std::string dir_path = "/data3/hchou/AMSData/MatTestProp";
    
    MatGeoBoxCreator creator_TRL1(
            MatTestProp::TRL1_N.at(0), MatTestProp::TRL1_MIN.at(0), MatTestProp::TRL1_MAX.at(0),
            MatTestProp::TRL1_N.at(1), MatTestProp::TRL1_MIN.at(1), MatTestProp::TRL1_MAX.at(1),
            MatTestProp::TRL1_N.at(2), MatTestProp::TRL1_MIN.at(2), MatTestProp::TRL1_MAX.at(2),
            CSTR_FMT("%s/TRL1.bin", dir_path.c_str())
        );
    
    MatGeoBoxCreator creator_TRL2(
            MatTestProp::TRL2_N.at(0), MatTestProp::TRL2_MIN.at(0), MatTestProp::TRL2_MAX.at(0),
            MatTestProp::TRL2_N.at(1), MatTestProp::TRL2_MIN.at(1), MatTestProp::TRL2_MAX.at(1),
            MatTestProp::TRL2_N.at(2), MatTestProp::TRL2_MIN.at(2), MatTestProp::TRL2_MAX.at(2),
            CSTR_FMT("%s/TRL2.bin", dir_path.c_str())
        );
    
    MatGeoBoxCreator creator_TRL3(
            MatTestProp::TRL3_N.at(0), MatTestProp::TRL3_MIN.at(0), MatTestProp::TRL3_MAX.at(0),
            MatTestProp::TRL3_N.at(1), MatTestProp::TRL3_MIN.at(1), MatTestProp::TRL3_MAX.at(1),
            MatTestProp::TRL3_N.at(2), MatTestProp::TRL3_MIN.at(2), MatTestProp::TRL3_MAX.at(2),
            CSTR_FMT("%s/TRL3.bin", dir_path.c_str())
        );
    
    MatGeoBoxCreator creator_MATC(
            MatTestProp::MATC_N.at(0), MatTestProp::MATC_MIN.at(0), MatTestProp::MATC_MAX.at(0),
            MatTestProp::MATC_N.at(1), MatTestProp::MATC_MIN.at(1), MatTestProp::MATC_MAX.at(1),
            MatTestProp::MATC_N.at(2), MatTestProp::MATC_MIN.at(2), MatTestProp::MATC_MAX.at(2),
            CSTR_FMT("%s/MATC.bin", dir_path.c_str())
        );
    
    Bool_t  TR_ELM[9] = { 0, 0, 0, 0, 0, 0, 0,     1, 0 };
    Float_t TR_DEN[9] = { 0, 0, 0, 0, 0, 0, 0, 0.083, 0 };
    creator_TRL1.save_and_close(TR_ELM, TR_DEN);
    creator_TRL2.save_and_close(TR_ELM, TR_DEN);
    creator_TRL3.save_and_close(TR_ELM, TR_DEN);
    
    Bool_t  C_ELM[9] = { 0,    1, 0, 0, 0, 0, 0, 0, 0 };
    Float_t C_DEN[9] = { 0, 0.01, 0, 0, 0, 0, 0, 0, 0 };
    creator_MATC.save_and_close(C_ELM, C_DEN);

    return true;
}


Bool_t MatGeoBoxTestProp::Load() {
    if (is_load_) return is_load_;
    std::string g4mat_dir_path = "/data3/hchou/AMSData/MatTestProp";

    reader_TRL1_.load(STR_FMT("%s/TRL1.bin", g4mat_dir_path.c_str()));
    reader_TRL2_.load(STR_FMT("%s/TRL2.bin", g4mat_dir_path.c_str()));
    reader_TRL3_.load(STR_FMT("%s/TRL3.bin", g4mat_dir_path.c_str()));
    reader_MATC_.load(STR_FMT("%s/MATC.bin", g4mat_dir_path.c_str()));
    
    reader_.push_back(&reader_TRL1_);
    reader_.push_back(&reader_TRL2_);
    reader_.push_back(&reader_TRL3_);
    reader_.push_back(&reader_MATC_);
   
    is_load_ = true;
    for (auto&& reader : reader_) {
        if (!reader->exist()) {
            is_load_ = false;
            break;
        }
    }

    return is_load_;
}


} // namespace TrackSys


#endif // __TRACKLibs_MatEnvTestProp_C__


#endif // __HAS_TESTPROP__
