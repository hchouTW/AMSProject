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
    std::string g4mat_file_path = "/afs/cern.ch/work/h/hchou/public/DATABASE/detector/g4mat_AMS02.root";
    TFile * root_file = TFile::Open(g4mat_file_path.c_str());
    if (root_file == nullptr || root_file->IsZombie()) return false;

    std::string dir_path = "/afs/cern.ch/work/h/hchou/public/DATABASE/detector/material";
    
    MatGeoBoxCreator creator_AMS02RAD(
            MatAms::RAD_N.at(0), MatAms::RAD_MIN.at(0), MatAms::RAD_MAX.at(0),
            MatAms::RAD_N.at(1), MatAms::RAD_MIN.at(1), MatAms::RAD_MAX.at(1),
            MatAms::RAD_N.at(2), MatAms::RAD_MIN.at(2), MatAms::RAD_MAX.at(2),
            CSTR_FMT("%s/AMS02RAD.bin", dir_path.c_str())
        );
    
    MatGeoBoxCreator creator_AMS02TRL1(
            MatAms::TRL1_N.at(0), MatAms::TRL1_MIN.at(0), MatAms::TRL1_MAX.at(0),
            MatAms::TRL1_N.at(1), MatAms::TRL1_MIN.at(1), MatAms::TRL1_MAX.at(1),
            MatAms::TRL1_N.at(2), MatAms::TRL1_MIN.at(2), MatAms::TRL1_MAX.at(2),
            CSTR_FMT("%s/AMS02TRL1.bin", dir_path.c_str())
        );
    
    MatGeoBoxCreator creator_AMS02UTRD(
            MatAms::UTRD_N.at(0), MatAms::UTRD_MIN.at(0), MatAms::UTRD_MAX.at(0),
            MatAms::UTRD_N.at(1), MatAms::UTRD_MIN.at(1), MatAms::UTRD_MAX.at(1),
            MatAms::UTRD_N.at(2), MatAms::UTRD_MIN.at(2), MatAms::UTRD_MAX.at(2),
            CSTR_FMT("%s/AMS02UTRD.bin", dir_path.c_str())
        );
    
    MatGeoBoxCreator creator_AMS02TRD(
            MatAms::TRD_N.at(0), MatAms::TRD_MIN.at(0), MatAms::TRD_MAX.at(0),
            MatAms::TRD_N.at(1), MatAms::TRD_MIN.at(1), MatAms::TRD_MAX.at(1),
            MatAms::TRD_N.at(2), MatAms::TRD_MIN.at(2), MatAms::TRD_MAX.at(2),
            CSTR_FMT("%s/AMS02TRD.bin", dir_path.c_str())
        );
    
    MatGeoBoxCreator creator_AMS02LTRD(
            MatAms::LTRD_N.at(0), MatAms::LTRD_MIN.at(0), MatAms::LTRD_MAX.at(0),
            MatAms::LTRD_N.at(1), MatAms::LTRD_MIN.at(1), MatAms::LTRD_MAX.at(1),
            MatAms::LTRD_N.at(2), MatAms::LTRD_MIN.at(2), MatAms::LTRD_MAX.at(2),
            CSTR_FMT("%s/AMS02LTRD.bin", dir_path.c_str())
        );
    
    MatGeoBoxCreator creator_AMS02UTOF(
            MatAms::UTOF_N.at(0), MatAms::UTOF_MIN.at(0), MatAms::UTOF_MAX.at(0),
            MatAms::UTOF_N.at(1), MatAms::UTOF_MIN.at(1), MatAms::UTOF_MAX.at(1),
            MatAms::UTOF_N.at(2), MatAms::UTOF_MIN.at(2), MatAms::UTOF_MAX.at(2),
            CSTR_FMT("%s/AMS02UTOF.bin", dir_path.c_str())
        );
    
    MatGeoBoxCreator creator_AMS02UITR(
            MatAms::UITR_N.at(0), MatAms::UITR_MIN.at(0), MatAms::UITR_MAX.at(0),
            MatAms::UITR_N.at(1), MatAms::UITR_MIN.at(1), MatAms::UITR_MAX.at(1),
            MatAms::UITR_N.at(2), MatAms::UITR_MIN.at(2), MatAms::UITR_MAX.at(2),
            CSTR_FMT("%s/AMS02UITR.bin", dir_path.c_str())
        );
    
    MatGeoBoxCreator creator_AMS02TRS1(
            MatAms::TRS1_N.at(0), MatAms::TRS1_MIN.at(0), MatAms::TRS1_MAX.at(0),
            MatAms::TRS1_N.at(1), MatAms::TRS1_MIN.at(1), MatAms::TRS1_MAX.at(1),
            MatAms::TRS1_N.at(2), MatAms::TRS1_MIN.at(2), MatAms::TRS1_MAX.at(2),
            CSTR_FMT("%s/AMS02TRS1.bin", dir_path.c_str())
        );
    
    MatGeoBoxCreator creator_AMS02TRS2(
            MatAms::TRS2_N.at(0), MatAms::TRS2_MIN.at(0), MatAms::TRS2_MAX.at(0),
            MatAms::TRS2_N.at(1), MatAms::TRS2_MIN.at(1), MatAms::TRS2_MAX.at(1),
            MatAms::TRS2_N.at(2), MatAms::TRS2_MIN.at(2), MatAms::TRS2_MAX.at(2),
            CSTR_FMT("%s/AMS02TRS2.bin", dir_path.c_str())
        );
    
    MatGeoBoxCreator creator_AMS02TRS3(
            MatAms::TRS3_N.at(0), MatAms::TRS3_MIN.at(0), MatAms::TRS3_MAX.at(0),
            MatAms::TRS3_N.at(1), MatAms::TRS3_MIN.at(1), MatAms::TRS3_MAX.at(1),
            MatAms::TRS3_N.at(2), MatAms::TRS3_MIN.at(2), MatAms::TRS3_MAX.at(2),
            CSTR_FMT("%s/AMS02TRS3.bin", dir_path.c_str())
        );
    
    MatGeoBoxCreator creator_AMS02LITR(
            MatAms::LITR_N.at(0), MatAms::LITR_MIN.at(0), MatAms::LITR_MAX.at(0),
            MatAms::LITR_N.at(1), MatAms::LITR_MIN.at(1), MatAms::LITR_MAX.at(1),
            MatAms::LITR_N.at(2), MatAms::LITR_MIN.at(2), MatAms::LITR_MAX.at(2),
            CSTR_FMT("%s/AMS02LITR.bin", dir_path.c_str())
        );
    
    MatGeoBoxCreator creator_AMS02LTOF(
            MatAms::LTOF_N.at(0), MatAms::LTOF_MIN.at(0), MatAms::LTOF_MAX.at(0),
            MatAms::LTOF_N.at(1), MatAms::LTOF_MIN.at(1), MatAms::LTOF_MAX.at(1),
            MatAms::LTOF_N.at(2), MatAms::LTOF_MIN.at(2), MatAms::LTOF_MAX.at(2),
            CSTR_FMT("%s/AMS02LTOF.bin", dir_path.c_str())
        );
    
    MatGeoBoxCreator creator_AMS02NAF(
            MatAms::NAF_N.at(0), MatAms::NAF_MIN.at(0), MatAms::NAF_MAX.at(0),
            MatAms::NAF_N.at(1), MatAms::NAF_MIN.at(1), MatAms::NAF_MAX.at(1),
            MatAms::NAF_N.at(2), MatAms::NAF_MIN.at(2), MatAms::NAF_MAX.at(2),
            CSTR_FMT("%s/AMS02NAF.bin", dir_path.c_str())
        );
    
    MatGeoBoxCreator creator_AMS02AGL(
            MatAms::AGL_N.at(0), MatAms::AGL_MIN.at(0), MatAms::AGL_MAX.at(0),
            MatAms::AGL_N.at(1), MatAms::AGL_MIN.at(1), MatAms::AGL_MAX.at(1),
            MatAms::AGL_N.at(2), MatAms::AGL_MIN.at(2), MatAms::AGL_MAX.at(2),
            CSTR_FMT("%s/AMS02AGL.bin", dir_path.c_str())
        );
    
    MatGeoBoxCreator creator_AMS02PMT(
            MatAms::PMT_N.at(0), MatAms::PMT_MIN.at(0), MatAms::PMT_MAX.at(0),
            MatAms::PMT_N.at(1), MatAms::PMT_MIN.at(1), MatAms::PMT_MAX.at(1),
            MatAms::PMT_N.at(2), MatAms::PMT_MIN.at(2), MatAms::PMT_MAX.at(2),
            CSTR_FMT("%s/AMS02PMT.bin", dir_path.c_str())
        );

    MatGeoBoxCreator creator_AMS02TRL9(
            MatAms::TRL9_N.at(0), MatAms::TRL9_MIN.at(0), MatAms::TRL9_MAX.at(0),
            MatAms::TRL9_N.at(1), MatAms::TRL9_MIN.at(1), MatAms::TRL9_MAX.at(1),
            MatAms::TRL9_N.at(2), MatAms::TRL9_MIN.at(2), MatAms::TRL9_MAX.at(2),
            CSTR_FMT("%s/AMS02TRL9.bin", dir_path.c_str())
        );

    MatGeoBoxCreator creator_AMS02ECAL(
            MatAms::ECAL_N.at(0), MatAms::ECAL_MIN.at(0), MatAms::ECAL_MAX.at(0),
            MatAms::ECAL_N.at(1), MatAms::ECAL_MIN.at(1), MatAms::ECAL_MAX.at(1),
            MatAms::ECAL_N.at(2), MatAms::ECAL_MIN.at(2), MatAms::ECAL_MAX.at(2),
            CSTR_FMT("%s/AMS02ECAL.bin", dir_path.c_str())
        );
    
    root_file->cd();

    Long64_t n[3];
    Float_t  min[3];
    Float_t  max[3];
    TTree * tree_inf = (TTree*) root_file->Get("g4mat_inf");
    tree_inf->SetBranchAddress("n",   &n);
    tree_inf->SetBranchAddress("min",  min);
    tree_inf->SetBranchAddress("max",  max);
    tree_inf->GetEntry(0);

    Bool_t  mat;
    Int_t   idx[3];
    Float_t coo[3];
    Bool_t  elm[9];
    Float_t den[9];
    TTree * tree_elm = (TTree*) root_file->Get("g4mat_elm");
    tree_elm->SetBranchAddress("mat", &mat);
    tree_elm->SetBranchAddress("idx",  idx);
    tree_elm->SetBranchAddress("coo",  coo);
    tree_elm->SetBranchAddress("elm",  elm);
    tree_elm->SetBranchAddress("den",  den);

    for (Long64_t entry = 0; entry < tree_elm->GetEntries(); ++entry) {
        tree_elm->GetEntry(entry);
        if (!mat) continue;
        Float_t radius = std::sqrt(coo[0]*coo[0] + coo[1]*coo[1]);
        
        Float_t trd_radius = (MatAms::TRDL_RADIUS + (coo[2] - MatAms::TRDL_Z) * MatAms::TRD_SLOPE);
        
        Bool_t is_trd      = (MGNumc::Compare(radius, trd_radius) < 0);
        Bool_t is_tracker  = (MGNumc::Compare(radius, MatAms::TRACKER_RADIUS) < 0);
        Bool_t is_mag_hole = (MGNumc::Compare(radius, MatAms::MAGNETIC_RADIUS) < 0);

        Bool_t is_naf = (std::max(std::fabs(coo[0]), std::fabs(coo[1])) < MatAms::RICH_BOUND);
        Bool_t is_agl = !is_naf;

        creator_AMS02RAD .fill(coo, elm, den);
        creator_AMS02TRL1.fill(coo, elm, den, is_tracker);
        creator_AMS02UTRD.fill(coo, elm, den, is_trd);
        creator_AMS02TRD .fill(coo, elm, den, is_trd);
        creator_AMS02LTRD.fill(coo, elm, den, is_trd);
        creator_AMS02UTOF.fill(coo, elm, den);
        creator_AMS02UITR.fill(coo, elm, den, is_tracker);
        creator_AMS02TRS1.fill(coo, elm, den, is_mag_hole);
        creator_AMS02TRS2.fill(coo, elm, den, is_mag_hole);
        creator_AMS02TRS3.fill(coo, elm, den, is_mag_hole);
        creator_AMS02LITR.fill(coo, elm, den, is_tracker);
        creator_AMS02LTOF.fill(coo, elm, den);
        creator_AMS02NAF .fill(coo, elm, den, is_naf);
        creator_AMS02AGL .fill(coo, elm, den, is_agl);
        creator_AMS02PMT .fill(coo, elm, den);
        creator_AMS02TRL9.fill(coo, elm, den);
        creator_AMS02ECAL.fill(coo, elm, den);
    }
    
    creator_AMS02RAD .save_and_close();
    creator_AMS02TRL1.save_and_close();
    creator_AMS02UTRD.save_and_close();
    creator_AMS02TRD .save_and_close();
    creator_AMS02LTRD.save_and_close();
    creator_AMS02UTOF.save_and_close();
    creator_AMS02UITR.save_and_close();
    creator_AMS02TRS1.save_and_close();
    creator_AMS02TRS2.save_and_close();
    creator_AMS02TRS3.save_and_close();
    creator_AMS02LITR.save_and_close();
    creator_AMS02LTOF.save_and_close();
    creator_AMS02NAF .save_and_close();
    creator_AMS02AGL .save_and_close();
    creator_AMS02PMT .save_and_close();
    creator_AMS02TRL9.save_and_close();
    creator_AMS02ECAL.save_and_close();

    root_file->Close();

    return true;
}


Bool_t MatGeoBoxAms::Load() {
    if (is_load_) return is_load_;
    std::string g4mat_dir_path = "/afs/cern.ch/work/h/hchou/public/DATABASE/detector/material";

    reader_AMS02RAD_ .load(STR_FMT("%s/AMS02RAD.bin" , g4mat_dir_path.c_str()));
    reader_AMS02TRL1_.load(STR_FMT("%s/AMS02TRL1.bin", g4mat_dir_path.c_str()));
    reader_AMS02UTRD_.load(STR_FMT("%s/AMS02UTRD.bin", g4mat_dir_path.c_str()));
    reader_AMS02TRD_ .load(STR_FMT("%s/AMS02TRD.bin" , g4mat_dir_path.c_str()));
    reader_AMS02LTRD_.load(STR_FMT("%s/AMS02LTRD.bin", g4mat_dir_path.c_str()));
    reader_AMS02UTOF_.load(STR_FMT("%s/AMS02UTOF.bin", g4mat_dir_path.c_str()));
    reader_AMS02UITR_.load(STR_FMT("%s/AMS02UITR.bin", g4mat_dir_path.c_str()));
    reader_AMS02TRS1_.load(STR_FMT("%s/AMS02TRS1.bin", g4mat_dir_path.c_str()));
    reader_AMS02TRS2_.load(STR_FMT("%s/AMS02TRS2.bin", g4mat_dir_path.c_str()));
    reader_AMS02TRS3_.load(STR_FMT("%s/AMS02TRS3.bin", g4mat_dir_path.c_str()));
    reader_AMS02LITR_.load(STR_FMT("%s/AMS02LITR.bin", g4mat_dir_path.c_str()));
    reader_AMS02LTOF_.load(STR_FMT("%s/AMS02LTOF.bin", g4mat_dir_path.c_str()));
    reader_AMS02NAF_ .load(STR_FMT("%s/AMS02NAF.bin",  g4mat_dir_path.c_str()));
    reader_AMS02AGL_ .load(STR_FMT("%s/AMS02AGL.bin",  g4mat_dir_path.c_str()));
    reader_AMS02PMT_ .load(STR_FMT("%s/AMS02PMT.bin" , g4mat_dir_path.c_str()));
    reader_AMS02TRL9_.load(STR_FMT("%s/AMS02TRL9.bin", g4mat_dir_path.c_str()));
    reader_AMS02ECAL_.load(STR_FMT("%s/AMS02ECAL.bin", g4mat_dir_path.c_str()));

    reader_.clear();
    reader_.push_back(&reader_AMS02RAD_ );
    reader_.push_back(&reader_AMS02TRL1_);
    reader_.push_back(&reader_AMS02UTRD_);
    reader_.push_back(&reader_AMS02TRD_ );
    reader_.push_back(&reader_AMS02LTRD_);
    reader_.push_back(&reader_AMS02UTOF_);
    reader_.push_back(&reader_AMS02UITR_);
    reader_.push_back(&reader_AMS02TRS1_);
    reader_.push_back(&reader_AMS02TRS2_);
    reader_.push_back(&reader_AMS02TRS3_);
    reader_.push_back(&reader_AMS02LITR_);
    reader_.push_back(&reader_AMS02LTOF_);
    reader_.push_back(&reader_AMS02NAF_ );
    reader_.push_back(&reader_AMS02AGL_ );
    reader_.push_back(&reader_AMS02PMT_ );
    reader_.push_back(&reader_AMS02TRL9_);
    reader_.push_back(&reader_AMS02ECAL_);
   
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
