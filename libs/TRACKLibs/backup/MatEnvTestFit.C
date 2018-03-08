#ifdef __HAS_TESTFIT__


#ifndef __TRACKLibs_MatEnvTestFit_C__
#define __TRACKLibs_MatEnvTestFit_C__


namespace TrackSys {


// Set to MatMgnt::Load()
Bool_t MatMgnt::Load() {
    if (is_load_ && reader_ != nullptr) return true;
    is_load_ = false;
    reader_ = nullptr;

    is_load_ = MatGeoBoxTestFit::Load();
    if (is_load_) reader_ = &MatGeoBoxTestFit::Reader();
    
    return is_load_;
};


Bool_t MatGeoBoxTestFit::CreateMatGeoBox() {
    std::string dir_path = "/data3/hchou/AMSData/MatTestFit";
    
    MatGeoBoxCreator creator_TRL01(
            MatTestFit::TRL01_N.at(0), MatTestFit::TRL01_MIN.at(0), MatTestFit::TRL01_MAX.at(0),
            MatTestFit::TRL01_N.at(1), MatTestFit::TRL01_MIN.at(1), MatTestFit::TRL01_MAX.at(1),
            MatTestFit::TRL01_N.at(2), MatTestFit::TRL01_MIN.at(2), MatTestFit::TRL01_MAX.at(2),
            CSTR_FMT("%s/TRL01.bin", dir_path.c_str())
        );
    
    MatGeoBoxCreator creator_TRL02(
            MatTestFit::TRL02_N.at(0), MatTestFit::TRL02_MIN.at(0), MatTestFit::TRL02_MAX.at(0),
            MatTestFit::TRL02_N.at(1), MatTestFit::TRL02_MIN.at(1), MatTestFit::TRL02_MAX.at(1),
            MatTestFit::TRL02_N.at(2), MatTestFit::TRL02_MIN.at(2), MatTestFit::TRL02_MAX.at(2),
            CSTR_FMT("%s/TRL02.bin", dir_path.c_str())
        );
    
    MatGeoBoxCreator creator_TRL03(
            MatTestFit::TRL03_N.at(0), MatTestFit::TRL03_MIN.at(0), MatTestFit::TRL03_MAX.at(0),
            MatTestFit::TRL03_N.at(1), MatTestFit::TRL03_MIN.at(1), MatTestFit::TRL03_MAX.at(1),
            MatTestFit::TRL03_N.at(2), MatTestFit::TRL03_MIN.at(2), MatTestFit::TRL03_MAX.at(2),
            CSTR_FMT("%s/TRL03.bin", dir_path.c_str())
        );
    
    MatGeoBoxCreator creator_TRL04(
            MatTestFit::TRL04_N.at(0), MatTestFit::TRL04_MIN.at(0), MatTestFit::TRL04_MAX.at(0),
            MatTestFit::TRL04_N.at(1), MatTestFit::TRL04_MIN.at(1), MatTestFit::TRL04_MAX.at(1),
            MatTestFit::TRL04_N.at(2), MatTestFit::TRL04_MIN.at(2), MatTestFit::TRL04_MAX.at(2),
            CSTR_FMT("%s/TRL04.bin", dir_path.c_str())
        );
    
    MatGeoBoxCreator creator_TRL05(
            MatTestFit::TRL05_N.at(0), MatTestFit::TRL05_MIN.at(0), MatTestFit::TRL05_MAX.at(0),
            MatTestFit::TRL05_N.at(1), MatTestFit::TRL05_MIN.at(1), MatTestFit::TRL05_MAX.at(1),
            MatTestFit::TRL05_N.at(2), MatTestFit::TRL05_MIN.at(2), MatTestFit::TRL05_MAX.at(2),
            CSTR_FMT("%s/TRL05.bin", dir_path.c_str())
        );
    
    MatGeoBoxCreator creator_TRL06(
            MatTestFit::TRL06_N.at(0), MatTestFit::TRL06_MIN.at(0), MatTestFit::TRL06_MAX.at(0),
            MatTestFit::TRL06_N.at(1), MatTestFit::TRL06_MIN.at(1), MatTestFit::TRL06_MAX.at(1),
            MatTestFit::TRL06_N.at(2), MatTestFit::TRL06_MIN.at(2), MatTestFit::TRL06_MAX.at(2),
            CSTR_FMT("%s/TRL06.bin", dir_path.c_str())
        );
    
    MatGeoBoxCreator creator_TRL07(
            MatTestFit::TRL07_N.at(0), MatTestFit::TRL07_MIN.at(0), MatTestFit::TRL07_MAX.at(0),
            MatTestFit::TRL07_N.at(1), MatTestFit::TRL07_MIN.at(1), MatTestFit::TRL07_MAX.at(1),
            MatTestFit::TRL07_N.at(2), MatTestFit::TRL07_MIN.at(2), MatTestFit::TRL07_MAX.at(2),
            CSTR_FMT("%s/TRL07.bin", dir_path.c_str())
        );
    
    MatGeoBoxCreator creator_TRL08(
            MatTestFit::TRL08_N.at(0), MatTestFit::TRL08_MIN.at(0), MatTestFit::TRL08_MAX.at(0),
            MatTestFit::TRL08_N.at(1), MatTestFit::TRL08_MIN.at(1), MatTestFit::TRL08_MAX.at(1),
            MatTestFit::TRL08_N.at(2), MatTestFit::TRL08_MIN.at(2), MatTestFit::TRL08_MAX.at(2),
            CSTR_FMT("%s/TRL08.bin", dir_path.c_str())
        );
    
    MatGeoBoxCreator creator_TRL09(
            MatTestFit::TRL09_N.at(0), MatTestFit::TRL09_MIN.at(0), MatTestFit::TRL09_MAX.at(0),
            MatTestFit::TRL09_N.at(1), MatTestFit::TRL09_MIN.at(1), MatTestFit::TRL09_MAX.at(1),
            MatTestFit::TRL09_N.at(2), MatTestFit::TRL09_MIN.at(2), MatTestFit::TRL09_MAX.at(2),
            CSTR_FMT("%s/TRL09.bin", dir_path.c_str())
        );
    
    MatGeoBoxCreator creator_TRL10(
            MatTestFit::TRL10_N.at(0), MatTestFit::TRL10_MIN.at(0), MatTestFit::TRL10_MAX.at(0),
            MatTestFit::TRL10_N.at(1), MatTestFit::TRL10_MIN.at(1), MatTestFit::TRL10_MAX.at(1),
            MatTestFit::TRL10_N.at(2), MatTestFit::TRL10_MIN.at(2), MatTestFit::TRL10_MAX.at(2),
            CSTR_FMT("%s/TRL10.bin", dir_path.c_str())
        );
    
    MatGeoBoxCreator creator_TRL11(
            MatTestFit::TRL11_N.at(0), MatTestFit::TRL11_MIN.at(0), MatTestFit::TRL11_MAX.at(0),
            MatTestFit::TRL11_N.at(1), MatTestFit::TRL11_MIN.at(1), MatTestFit::TRL11_MAX.at(1),
            MatTestFit::TRL11_N.at(2), MatTestFit::TRL11_MIN.at(2), MatTestFit::TRL11_MAX.at(2),
            CSTR_FMT("%s/TRL11.bin", dir_path.c_str())
        );
    
    MatGeoBoxCreator creator_TRL12(
            MatTestFit::TRL12_N.at(0), MatTestFit::TRL12_MIN.at(0), MatTestFit::TRL12_MAX.at(0),
            MatTestFit::TRL12_N.at(1), MatTestFit::TRL12_MIN.at(1), MatTestFit::TRL12_MAX.at(1),
            MatTestFit::TRL12_N.at(2), MatTestFit::TRL12_MIN.at(2), MatTestFit::TRL12_MAX.at(2),
            CSTR_FMT("%s/TRL12.bin", dir_path.c_str())
        );
    
    MatGeoBoxCreator creator_MATC01(
            MatTestFit::MATC01_N.at(0), MatTestFit::MATC01_MIN.at(0), MatTestFit::MATC01_MAX.at(0),
            MatTestFit::MATC01_N.at(1), MatTestFit::MATC01_MIN.at(1), MatTestFit::MATC01_MAX.at(1),
            MatTestFit::MATC01_N.at(2), MatTestFit::MATC01_MIN.at(2), MatTestFit::MATC01_MAX.at(2),
            CSTR_FMT("%s/MATC01.bin", dir_path.c_str())
        );
    
    MatGeoBoxCreator creator_MATC02(
            MatTestFit::MATC02_N.at(0), MatTestFit::MATC02_MIN.at(0), MatTestFit::MATC02_MAX.at(0),
            MatTestFit::MATC02_N.at(1), MatTestFit::MATC02_MIN.at(1), MatTestFit::MATC02_MAX.at(1),
            MatTestFit::MATC02_N.at(2), MatTestFit::MATC02_MIN.at(2), MatTestFit::MATC02_MAX.at(2),
            CSTR_FMT("%s/MATC02.bin", dir_path.c_str())
        );
    
    MatGeoBoxCreator creator_MATAL01(
            MatTestFit::MATAL01_N.at(0), MatTestFit::MATAL01_MIN.at(0), MatTestFit::MATAL01_MAX.at(0),
            MatTestFit::MATAL01_N.at(1), MatTestFit::MATAL01_MIN.at(1), MatTestFit::MATAL01_MAX.at(1),
            MatTestFit::MATAL01_N.at(2), MatTestFit::MATAL01_MIN.at(2), MatTestFit::MATAL01_MAX.at(2),
            CSTR_FMT("%s/MATAL01.bin", dir_path.c_str())
        );
    
    MatGeoBoxCreator creator_MATAL02(
            MatTestFit::MATAL02_N.at(0), MatTestFit::MATAL02_MIN.at(0), MatTestFit::MATAL02_MAX.at(0),
            MatTestFit::MATAL02_N.at(1), MatTestFit::MATAL02_MIN.at(1), MatTestFit::MATAL02_MAX.at(1),
            MatTestFit::MATAL02_N.at(2), MatTestFit::MATAL02_MIN.at(2), MatTestFit::MATAL02_MAX.at(2),
            CSTR_FMT("%s/MATAL02.bin", dir_path.c_str())
        );
    
    MatGeoBoxCreator creator_MATAL03(
            MatTestFit::MATAL03_N.at(0), MatTestFit::MATAL03_MIN.at(0), MatTestFit::MATAL03_MAX.at(0),
            MatTestFit::MATAL03_N.at(1), MatTestFit::MATAL03_MIN.at(1), MatTestFit::MATAL03_MAX.at(1),
            MatTestFit::MATAL03_N.at(2), MatTestFit::MATAL03_MIN.at(2), MatTestFit::MATAL03_MAX.at(2),
            CSTR_FMT("%s/MATAL03.bin", dir_path.c_str())
        );
    
    MatGeoBoxCreator creator_MATAL04(
            MatTestFit::MATAL04_N.at(0), MatTestFit::MATAL04_MIN.at(0), MatTestFit::MATAL04_MAX.at(0),
            MatTestFit::MATAL04_N.at(1), MatTestFit::MATAL04_MIN.at(1), MatTestFit::MATAL04_MAX.at(1),
            MatTestFit::MATAL04_N.at(2), MatTestFit::MATAL04_MIN.at(2), MatTestFit::MATAL04_MAX.at(2),
            CSTR_FMT("%s/MATAL04.bin", dir_path.c_str())
        );
    
    MatGeoBoxCreator creator_MATAL05(
            MatTestFit::MATAL05_N.at(0), MatTestFit::MATAL05_MIN.at(0), MatTestFit::MATAL05_MAX.at(0),
            MatTestFit::MATAL05_N.at(1), MatTestFit::MATAL05_MIN.at(1), MatTestFit::MATAL05_MAX.at(1),
            MatTestFit::MATAL05_N.at(2), MatTestFit::MATAL05_MIN.at(2), MatTestFit::MATAL05_MAX.at(2),
            CSTR_FMT("%s/MATAL05.bin", dir_path.c_str())
        );
    
    MatGeoBoxCreator creator_MATAL06(
            MatTestFit::MATAL06_N.at(0), MatTestFit::MATAL06_MIN.at(0), MatTestFit::MATAL06_MAX.at(0),
            MatTestFit::MATAL06_N.at(1), MatTestFit::MATAL06_MIN.at(1), MatTestFit::MATAL06_MAX.at(1),
            MatTestFit::MATAL06_N.at(2), MatTestFit::MATAL06_MIN.at(2), MatTestFit::MATAL06_MAX.at(2),
            CSTR_FMT("%s/MATAL06.bin", dir_path.c_str())
        );
    
    Bool_t   TR_ELM[9] = { 0, 0, 0, 0, 0, 0, 0,     1, 0 };
    Double_t TR_DEN[9] = { 0, 0, 0, 0, 0, 0, 0, 0.083, 0 };
    creator_TRL01.save_and_close(TR_ELM, TR_DEN);
    creator_TRL02.save_and_close(TR_ELM, TR_DEN);
    creator_TRL03.save_and_close(TR_ELM, TR_DEN);
    creator_TRL04.save_and_close(TR_ELM, TR_DEN);
    creator_TRL05.save_and_close(TR_ELM, TR_DEN);
    creator_TRL06.save_and_close(TR_ELM, TR_DEN);
    creator_TRL07.save_and_close(TR_ELM, TR_DEN);
    creator_TRL08.save_and_close(TR_ELM, TR_DEN);
    creator_TRL09.save_and_close(TR_ELM, TR_DEN);
    creator_TRL10.save_and_close(TR_ELM, TR_DEN);
    creator_TRL11.save_and_close(TR_ELM, TR_DEN);
    creator_TRL12.save_and_close(TR_ELM, TR_DEN);
    
    Bool_t   C_ELM[9] = { 0,    1, 0, 0, 0, 0, 0, 0, 0 };
    Double_t C_DEN[9] = { 0, 0.08, 0, 0, 0, 0, 0, 0, 0 };
    creator_MATC01.save_and_close(C_ELM, C_DEN);
    creator_MATC02.save_and_close(C_ELM, C_DEN);
    
    Bool_t   AL_ELM[9] = { 0, 0, 0, 0, 0, 0,   1, 0, 0 };
    Double_t AL_DEN[9] = { 0, 0, 0, 0, 0, 0, 0.1, 0, 0 };
    creator_MATAL01.save_and_close(AL_ELM, AL_DEN);
    creator_MATAL02.save_and_close(AL_ELM, AL_DEN);
    creator_MATAL03.save_and_close(AL_ELM, AL_DEN);
    creator_MATAL04.save_and_close(AL_ELM, AL_DEN);
    creator_MATAL05.save_and_close(AL_ELM, AL_DEN);
    creator_MATAL06.save_and_close(AL_ELM, AL_DEN);

    return true;
}


Bool_t MatGeoBoxTestFit::Load() {
    if (is_load_) return is_load_;
    std::string g4mat_dir_path = "/data3/hchou/AMSData/MatTestFit";

    reader_TRL01_.load(STR_FMT("%s/TRL01.bin", g4mat_dir_path.c_str()));
    reader_TRL02_.load(STR_FMT("%s/TRL02.bin", g4mat_dir_path.c_str()));
    reader_TRL03_.load(STR_FMT("%s/TRL03.bin", g4mat_dir_path.c_str()));
    reader_TRL04_.load(STR_FMT("%s/TRL04.bin", g4mat_dir_path.c_str()));
    reader_TRL05_.load(STR_FMT("%s/TRL05.bin", g4mat_dir_path.c_str()));
    reader_TRL06_.load(STR_FMT("%s/TRL06.bin", g4mat_dir_path.c_str()));
    reader_TRL07_.load(STR_FMT("%s/TRL07.bin", g4mat_dir_path.c_str()));
    reader_TRL08_.load(STR_FMT("%s/TRL08.bin", g4mat_dir_path.c_str()));
    reader_TRL09_.load(STR_FMT("%s/TRL09.bin", g4mat_dir_path.c_str()));
    reader_TRL10_.load(STR_FMT("%s/TRL10.bin", g4mat_dir_path.c_str()));
    reader_TRL11_.load(STR_FMT("%s/TRL11.bin", g4mat_dir_path.c_str()));
    reader_TRL12_.load(STR_FMT("%s/TRL12.bin", g4mat_dir_path.c_str()));
    reader_MATC01_.load(STR_FMT("%s/MATC01.bin", g4mat_dir_path.c_str()));
    reader_MATC02_.load(STR_FMT("%s/MATC02.bin", g4mat_dir_path.c_str()));
    reader_MATAL01_.load(STR_FMT("%s/MATAL01.bin", g4mat_dir_path.c_str()));
    reader_MATAL02_.load(STR_FMT("%s/MATAL02.bin", g4mat_dir_path.c_str()));
    reader_MATAL03_.load(STR_FMT("%s/MATAL03.bin", g4mat_dir_path.c_str()));
    reader_MATAL04_.load(STR_FMT("%s/MATAL04.bin", g4mat_dir_path.c_str()));
    reader_MATAL05_.load(STR_FMT("%s/MATAL05.bin", g4mat_dir_path.c_str()));
    reader_MATAL06_.load(STR_FMT("%s/MATAL06.bin", g4mat_dir_path.c_str()));
    
    reader_.push_back(&reader_TRL01_);
    reader_.push_back(&reader_TRL02_);
    reader_.push_back(&reader_TRL03_);
    reader_.push_back(&reader_TRL04_);
    reader_.push_back(&reader_TRL05_);
    reader_.push_back(&reader_TRL06_);
    reader_.push_back(&reader_TRL07_);
    reader_.push_back(&reader_TRL08_);
    reader_.push_back(&reader_TRL09_);
    reader_.push_back(&reader_TRL10_);
    reader_.push_back(&reader_TRL11_);
    reader_.push_back(&reader_TRL12_);
    reader_.push_back(&reader_MATC01_);
    reader_.push_back(&reader_MATC02_);
    reader_.push_back(&reader_MATAL01_);
    reader_.push_back(&reader_MATAL02_);
    reader_.push_back(&reader_MATAL03_);
    reader_.push_back(&reader_MATAL04_);
    reader_.push_back(&reader_MATAL05_);
    reader_.push_back(&reader_MATAL06_);
   
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


#endif // __TRACKLibs_MatEnvTestFit_C__


#endif // __HAS_TESTFIT__
