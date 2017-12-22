#ifdef __HAS_AMS_OFFICE_LIBS__


#ifndef __TRACKLibs_MatEnvAms_C__
#define __TRACKLibs_MatEnvAms_C__


namespace TrackSys {


// Set to MatMgnt::Load()
Bool_t MatMgnt::Load() {
    if (is_load_ && reader_ != nullptr) return true;
    is_load_ = false;
    reader_ = nullptr;

    is_load_ = MatGeoBoxAms::Load();
    if (is_load_) reader_ = &MatGeoBoxAms::Reader();
    
    return is_load_;
};


Bool_t MatGeoBoxAms::CreateMatGeoBoxFromG4MatTree() {
    std::string dpath = "/data1/hchou/material15"; // at NCU
   
    std::string flepath = "/data1/hchou/newstr_material/scan14/g4mscan.root";
    //TFile* ifle = TFile::Open(CSTR_FMT("%s/g4mscan.root", dpath.c_str()));
    TFile* ifle = TFile::Open(flepath.c_str());
    if (ifle == nullptr || ifle->IsZombie()) return false;
    TTree* tree = (TTree*)ifle->Get("tree"); 
    G4MatStep g4mat(tree);
   
    // TRD
    MatGeoBoxCreator creator_TRL1(&MatAms::TRL1_N.front(), &MatAms::TRL1_MIN.front(), &MatAms::TRL1_MAX.front(), MatAms::TRL1_STP, "AMS02TRL1", dpath.c_str());
    MatGeoBoxCreator creator_TRL2(&MatAms::TRL2_N.front(), &MatAms::TRL2_MIN.front(), &MatAms::TRL2_MAX.front(), MatAms::TRL2_STP, "AMS02TRL2", dpath.c_str());
    MatGeoBoxCreator creator_TRL3(&MatAms::TRL3_N.front(), &MatAms::TRL3_MIN.front(), &MatAms::TRL3_MAX.front(), MatAms::TRL3_STP, "AMS02TRL3", dpath.c_str());
    MatGeoBoxCreator creator_TRL4(&MatAms::TRL4_N.front(), &MatAms::TRL4_MIN.front(), &MatAms::TRL4_MAX.front(), MatAms::TRL4_STP, "AMS02TRL4", dpath.c_str());
    MatGeoBoxCreator creator_TRL5(&MatAms::TRL5_N.front(), &MatAms::TRL5_MIN.front(), &MatAms::TRL5_MAX.front(), MatAms::TRL5_STP, "AMS02TRL5", dpath.c_str());
    MatGeoBoxCreator creator_TRL6(&MatAms::TRL6_N.front(), &MatAms::TRL6_MIN.front(), &MatAms::TRL6_MAX.front(), MatAms::TRL6_STP, "AMS02TRL6", dpath.c_str());
    MatGeoBoxCreator creator_TRL7(&MatAms::TRL7_N.front(), &MatAms::TRL7_MIN.front(), &MatAms::TRL7_MAX.front(), MatAms::TRL7_STP, "AMS02TRL7", dpath.c_str());
    MatGeoBoxCreator creator_TRL8(&MatAms::TRL8_N.front(), &MatAms::TRL8_MIN.front(), &MatAms::TRL8_MAX.front(), MatAms::TRL8_STP, "AMS02TRL8", dpath.c_str());
    MatGeoBoxCreator creator_TRL9(&MatAms::TRL9_N.front(), &MatAms::TRL9_MIN.front(), &MatAms::TRL9_MAX.front(), MatAms::TRL9_STP, "AMS02TRL9", dpath.c_str());
    
    MatGeoBoxCreator creator_TR34(&MatAms::TR34_N.front(), &MatAms::TR34_MIN.front(), &MatAms::TR34_MAX.front(), MatAms::TR34_STP, "AMS02TR34", dpath.c_str());
    MatGeoBoxCreator creator_TR56(&MatAms::TR56_N.front(), &MatAms::TR56_MIN.front(), &MatAms::TR56_MAX.front(), MatAms::TR56_STP, "AMS02TR56", dpath.c_str());
    MatGeoBoxCreator creator_TR78(&MatAms::TR78_N.front(), &MatAms::TR78_MIN.front(), &MatAms::TR78_MAX.front(), MatAms::TR78_STP, "AMS02TR78", dpath.c_str());
    
    MatGeoBoxCreator creator_TRS1(&MatAms::TRS1_N.front(), &MatAms::TRS1_MIN.front(), &MatAms::TRS1_MAX.front(), MatAms::TRS1_STP, "AMS02TRS1", dpath.c_str());
    MatGeoBoxCreator creator_TRS2(&MatAms::TRS2_N.front(), &MatAms::TRS2_MIN.front(), &MatAms::TRS2_MAX.front(), MatAms::TRS2_STP, "AMS02TRS2", dpath.c_str());
    MatGeoBoxCreator creator_TRS3(&MatAms::TRS3_N.front(), &MatAms::TRS3_MIN.front(), &MatAms::TRS3_MAX.front(), MatAms::TRS3_STP, "AMS02TRS3", dpath.c_str());
    MatGeoBoxCreator creator_TRS4(&MatAms::TRS4_N.front(), &MatAms::TRS4_MIN.front(), &MatAms::TRS4_MAX.front(), MatAms::TRS4_STP, "AMS02TRS4", dpath.c_str());
    MatGeoBoxCreator creator_TRS5(&MatAms::TRS5_N.front(), &MatAms::TRS5_MIN.front(), &MatAms::TRS5_MAX.front(), MatAms::TRS5_STP, "AMS02TRS5", dpath.c_str());
    MatGeoBoxCreator creator_TRS6(&MatAms::TRS6_N.front(), &MatAms::TRS6_MIN.front(), &MatAms::TRS6_MAX.front(), MatAms::TRS6_STP, "AMS02TRS6", dpath.c_str());
    MatGeoBoxCreator creator_TRS7(&MatAms::TRS7_N.front(), &MatAms::TRS7_MIN.front(), &MatAms::TRS7_MAX.front(), MatAms::TRS7_STP, "AMS02TRS7", dpath.c_str());
    MatGeoBoxCreator creator_TRS8(&MatAms::TRS8_N.front(), &MatAms::TRS8_MIN.front(), &MatAms::TRS8_MAX.front(), MatAms::TRS8_STP, "AMS02TRS8", dpath.c_str());
    MatGeoBoxCreator creator_TRS9(&MatAms::TRS9_N.front(), &MatAms::TRS9_MIN.front(), &MatAms::TRS9_MAX.front(), MatAms::TRS9_STP, "AMS02TRS9", dpath.c_str());

    MatGeoBoxCreator creator_RAD (&MatAms::RAD_N.front() , &MatAms::RAD_MIN.front() , &MatAms::RAD_MAX.front() , MatAms::RAD_STP , "AMS02RAD" , dpath.c_str());
    MatGeoBoxCreator creator_TRD (&MatAms::TRD_N.front() , &MatAms::TRD_MIN.front() , &MatAms::TRD_MAX.front() , MatAms::TRD_STP , "AMS02TRD" , dpath.c_str());
    MatGeoBoxCreator creator_SUPU(&MatAms::SUPU_N.front(), &MatAms::SUPU_MIN.front(), &MatAms::SUPU_MAX.front(), MatAms::SUPU_STP, "AMS02SUPU", dpath.c_str());
    MatGeoBoxCreator creator_SPIU(&MatAms::SPIU_N.front(), &MatAms::SPIU_MIN.front(), &MatAms::SPIU_MAX.front(), MatAms::SPIU_STP, "AMS02SPIU", dpath.c_str());
    MatGeoBoxCreator creator_TOFU(&MatAms::TOFU_N.front(), &MatAms::TOFU_MIN.front(), &MatAms::TOFU_MAX.front(), MatAms::TOFU_STP, "AMS02TOFU", dpath.c_str());
    MatGeoBoxCreator creator_TOFL(&MatAms::TOFL_N.front(), &MatAms::TOFL_MIN.front(), &MatAms::TOFL_MAX.front(), MatAms::TOFL_STP, "AMS02TOFL", dpath.c_str());
    MatGeoBoxCreator creator_SPIL(&MatAms::SPIL_N.front(), &MatAms::SPIL_MIN.front(), &MatAms::SPIL_MAX.front(), MatAms::SPIL_STP, "AMS02SPIL", dpath.c_str());
    MatGeoBoxCreator creator_SUPL(&MatAms::SUPL_N.front(), &MatAms::SUPL_MIN.front(), &MatAms::SUPL_MAX.front(), MatAms::SUPL_STP, "AMS02SUPL", dpath.c_str());
    MatGeoBoxCreator creator_RICH(&MatAms::RICH_N.front(), &MatAms::RICH_MIN.front(), &MatAms::RICH_MAX.front(), MatAms::RICH_STP, "AMS02RICH", dpath.c_str());
    MatGeoBoxCreator creator_PMT (&MatAms::PMT_N.front() , &MatAms::PMT_MIN.front() , &MatAms::PMT_MAX.front() , MatAms::PMT_STP , "AMS02PMT" , dpath.c_str());
    MatGeoBoxCreator creator_ECAL(&MatAms::ECAL_N.front(), &MatAms::ECAL_MIN.front(), &MatAms::ECAL_MAX.front(), MatAms::ECAL_STP, "AMS02ECAL", dpath.c_str());

    COUT("====================================================================\n");
    COUT("===========  MatGeoBoxAms::CreateMatGeoBoxFromG4MatTree  ===========\n");
    COUT("====================================================================\n");
    Long64_t nentries = tree->GetEntries();
    Long64_t printRat = (nentries / 100);
    for (Long64_t entry = 0; entry < nentries; ++entry) {
        if ((entry%printRat) == 0 || (entry == nentries-1)) {
            COUT("Entry %lld/%lld  (%6.2f %)\n", entry+1, nentries, 100. * static_cast<Double_t>(entry+1)/static_cast<Double_t>(nentries));
        }
        tree->GetEntry(entry);
        
        creator_TRL1.fill(g4mat);
        creator_TRL2.fill(g4mat);
        creator_TRL3.fill(g4mat);
        creator_TRL4.fill(g4mat);
        creator_TRL5.fill(g4mat);
        creator_TRL6.fill(g4mat);
        creator_TRL7.fill(g4mat);
        creator_TRL8.fill(g4mat);
        creator_TRL9.fill(g4mat);
        
        creator_TR34.fill(g4mat);
        creator_TR56.fill(g4mat);
        creator_TR78.fill(g4mat);
        
        creator_TRS1.fill(g4mat);
        creator_TRS2.fill(g4mat);
        creator_TRS3.fill(g4mat);
        creator_TRS4.fill(g4mat);
        creator_TRS5.fill(g4mat);
        creator_TRS6.fill(g4mat);
        creator_TRS7.fill(g4mat);
        creator_TRS8.fill(g4mat);
        creator_TRS9.fill(g4mat);

        creator_RAD.fill(g4mat);
        creator_TRD.fill(g4mat);
        creator_SUPU.fill(g4mat);
        creator_SPIU.fill(g4mat);
        creator_TOFU.fill(g4mat);
        creator_TOFL.fill(g4mat);
        creator_SPIL.fill(g4mat);
        creator_SUPL.fill(g4mat);
        creator_RICH.fill(g4mat);
        creator_PMT.fill(g4mat);
        creator_ECAL.fill(g4mat);
    }
    COUT("====================================================================\n");
    
    creator_TRL1.save_and_close();
    creator_TRL2.save_and_close();
    creator_TRL3.save_and_close();
    creator_TRL4.save_and_close();
    creator_TRL5.save_and_close();
    creator_TRL6.save_and_close();
    creator_TRL7.save_and_close();
    creator_TRL8.save_and_close();
    creator_TRL9.save_and_close();
    
    creator_TR34.save_and_close();
    creator_TR56.save_and_close();
    creator_TR78.save_and_close();
    
    creator_TRS1.save_and_close();
    creator_TRS2.save_and_close();
    creator_TRS3.save_and_close();
    creator_TRS4.save_and_close();
    creator_TRS5.save_and_close();
    creator_TRS6.save_and_close();
    creator_TRS7.save_and_close();
    creator_TRS8.save_and_close();
    creator_TRS9.save_and_close();

    creator_RAD.save_and_close();
    creator_TRD.save_and_close();
    creator_SUPU.save_and_close();
    creator_SPIU.save_and_close();
    creator_TOFU.save_and_close();
    creator_TOFL.save_and_close();
    creator_SPIL.save_and_close();
    creator_SUPL.save_and_close();
    creator_RICH.save_and_close();
    creator_PMT.save_and_close();
    creator_ECAL.save_and_close();

    ifle->Close();

    return true;
}


Bool_t MatGeoBoxAms::Load() {
    if (is_load_) return is_load_;
    std::string dpath = "/data1/hchou/material15"; // at NCU

    reader_TRL1_.load("AMS02TRL1" , dpath);
    reader_TRL2_.load("AMS02TRL2" , dpath);
    reader_TRL3_.load("AMS02TRL3" , dpath);
    reader_TRL4_.load("AMS02TRL4" , dpath);
    reader_TRL5_.load("AMS02TRL5" , dpath);
    reader_TRL6_.load("AMS02TRL6" , dpath);
    reader_TRL7_.load("AMS02TRL7" , dpath);
    reader_TRL8_.load("AMS02TRL8" , dpath);
    reader_TRL9_.load("AMS02TRL9" , dpath);
    
    reader_TR34_.load("AMS02TR34" , dpath);
    reader_TR56_.load("AMS02TR56" , dpath);
    reader_TR78_.load("AMS02TR78" , dpath);
    
    reader_TRS1_.load("AMS02TRS1" , dpath);
    reader_TRS2_.load("AMS02TRS2" , dpath);
    reader_TRS3_.load("AMS02TRS3" , dpath);
    reader_TRS4_.load("AMS02TRS4" , dpath);
    reader_TRS5_.load("AMS02TRS5" , dpath);
    reader_TRS6_.load("AMS02TRS6" , dpath);
    reader_TRS7_.load("AMS02TRS7" , dpath);
    reader_TRS8_.load("AMS02TRS8" , dpath);
    reader_TRS9_.load("AMS02TRS9" , dpath);
    
    reader_RAD_.load("AMS02RAD" , dpath);
    reader_TRD_.load("AMS02TRD" , dpath);
    reader_SUPU_.load("AMS02SUPU" , dpath);
    reader_SPIU_.load("AMS02SPIU" , dpath);
    reader_TOFU_.load("AMS02TOFU" , dpath);
    reader_TOFL_.load("AMS02TOFL" , dpath);
    reader_SPIL_.load("AMS02SPIL" , dpath);
    reader_SUPL_.load("AMS02SUPL" , dpath);
    reader_RICH_.load("AMS02RICH" , dpath);
    reader_PMT_.load("AMS02PMT" , dpath);
    reader_ECAL_.load("AMS02ECAL" , dpath);
   
    reader_.clear();
    
    reader_.push_back(&reader_TRL1_);
    reader_.push_back(&reader_TRL2_);
    reader_.push_back(&reader_TRL3_);
    reader_.push_back(&reader_TRL4_);
    reader_.push_back(&reader_TRL5_);
    reader_.push_back(&reader_TRL6_);
    reader_.push_back(&reader_TRL7_);
    reader_.push_back(&reader_TRL8_);
    reader_.push_back(&reader_TRL9_);
    
    reader_.push_back(&reader_TR34_);
    reader_.push_back(&reader_TR56_);
    reader_.push_back(&reader_TR78_);
    
    reader_.push_back(&reader_TRS1_);
    reader_.push_back(&reader_TRS2_);
    reader_.push_back(&reader_TRS3_);
    reader_.push_back(&reader_TRS4_);
    reader_.push_back(&reader_TRS5_);
    reader_.push_back(&reader_TRS6_);
    reader_.push_back(&reader_TRS7_);
    reader_.push_back(&reader_TRS8_);
    reader_.push_back(&reader_TRS9_);

    reader_.push_back(&reader_RAD_);
    reader_.push_back(&reader_TRD_);
    reader_.push_back(&reader_SUPU_);
    reader_.push_back(&reader_SPIU_);
    reader_.push_back(&reader_TOFU_);
    reader_.push_back(&reader_TOFL_);
    reader_.push_back(&reader_SPIL_);
    reader_.push_back(&reader_SUPL_);
    reader_.push_back(&reader_RICH_);
    reader_.push_back(&reader_PMT_);
    reader_.push_back(&reader_ECAL_);
   
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


#endif // __TRACKLibs_MatEnvAms_C__


#endif // __HAS_AMS_OFFICE_LIBS__
