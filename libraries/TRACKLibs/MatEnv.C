#ifndef __TRACKLibs_MatEnv_C__
#define __TRACKLibs_MatEnv_C__

namespace TrackSys {

MatGeoBoxCreator::MatGeoBoxCreator(Long64_t xn, Double_t xmin, Double_t xmax, Long64_t yn, Double_t ymin, Double_t ymax, Long64_t zn, Double_t zmin, Double_t zmax, const std::string& file_path) {
    clear();
    if (xn < 1 || yn < 1 || zn < 1) { MGSys::ShowError("MatGeoBoxCreator::MatGeoBoxCreator() : Size failure."); return; }
    if (MGNumc::Compare(xmin, xmax) >= 0 || MGNumc::Compare(ymin, ymax) >= 0 || MGNumc::Compare(zmin, zmax) >= 0) { MGSys::ShowError("MatGeoBoxCreator::MatGeoBoxCreator() : Range failure."); return; }
    if (file_path.size() == 0) { MGSys::ShowError("MatGeoBoxCreator::MatGeoBoxCreator() : File path not found."); return; }
    
    Int_t file_len = ((sizeof(Long64_t) + sizeof(Double_t) + sizeof(Double_t)) * DIM_ + (sizeof(Bool_t) + sizeof(Double_t)) * MatProperty::NUM_ELM + sizeof(Double_t) + (xn * yn * zn) * sizeof(Bool_t));
    Int_t file_des = open(file_path.c_str(), O_CREAT | O_RDWR, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
    if (file_des < 0) { MGSys::ShowError("MatGeoBoxCreator::MatGeoBoxCreator() : File not opened."); return; }

    off_t sret = lseek(file_des, file_len, SEEK_SET);
    write(file_des, "\0", 1);

    void* file_ptr = mmap(nullptr, file_len, PROT_READ | PROT_WRITE, MAP_SHARED, file_des, 0);
    if (file_ptr == reinterpret_cast<void*>(-1)) { MGSys::ShowError("MatGeoBoxCreator::MatGeoBoxCreator() : mmap() failure."); close(file_des); return; }
    
    is_open_   = true;
    file_path_ = file_path;
    file_des_  = file_des;
    file_len_  = file_len;
    file_ptr_  = file_ptr;
    max_len_   = (xn * yn * zn);
    cur_len_   = 0;
    geo_box_   = reinterpret_cast<MatGeoBox*>(file_ptr);
    elm_.fill(false);
    cnt_ = 0;
    den_.fill(0);
    inv_rad_len_ = 0.;

    geo_box_->n[0] = xn;
    geo_box_->n[1] = yn;
    geo_box_->n[2] = zn;
    geo_box_->min[0] = xmin;
    geo_box_->min[1] = ymin;
    geo_box_->min[2] = zmin;
    geo_box_->max[0] = xmax;
    geo_box_->max[1] = ymax;
    geo_box_->max[2] = zmax;
    std::fill_n(geo_box_->elm, 9, false);
    std::fill_n(geo_box_->den, 9, 0.0);
    geo_box_->inv_rad_len = 0.;
    std::fill_n(static_cast<Bool_t*>(&(geo_box_->mat)), max_len_, false);
}


void MatGeoBoxCreator::fill(Bool_t elm[MatProperty::NUM_ELM], Float_t den[MatProperty::NUM_ELM], Bool_t calculated) {
    if (!is_open_) return;
    if (cur_len_ >= max_len_) { MGSys::ShowError("MatGeoBoxCreator::fill() : Out range."); return; }

    Bool_t has_mat = false;
    if (calculated) {
        for (Int_t it = 0; it < MatProperty::NUM_ELM; ++it) {
            if (!elm[it]) continue;
            has_mat = true;
            elm_.at(it) = true;
            den_.at(it) += den[it];
        }
        if (has_mat) cnt_++;
    }

    Bool_t& mat = *(static_cast<Bool_t*>(&(geo_box_->mat) + cur_len_));
    mat = has_mat;
    cur_len_++;
}


Bool_t MatGeoBoxCreator::is_in_box(Float_t coo[3]) {
    if (!is_open_) return false;
    if (MGNumc::Compare(coo[2], static_cast<Float_t>(geo_box_->min[2])) <= 0 || MGNumc::Compare(coo[2], static_cast<Float_t>(geo_box_->max[2])) >= 0) return false;
    if (MGNumc::Compare(coo[1], static_cast<Float_t>(geo_box_->min[1])) <= 0 || MGNumc::Compare(coo[1], static_cast<Float_t>(geo_box_->max[1])) >= 0) return false;
    if (MGNumc::Compare(coo[0], static_cast<Float_t>(geo_box_->min[0])) <= 0 || MGNumc::Compare(coo[0], static_cast<Float_t>(geo_box_->max[0])) >= 0) return false;
    return true;
}


void MatGeoBoxCreator::save_and_close() {
    if (!is_open_) { clear(); return; }
    if (cur_len_ < max_len_) {
        MGSys::ShowWarning(STR_FMT("MatGeoBoxCreator::close() : fill not finished. (%lld/%lld)", cur_len_, max_len_));
        Bool_t  elm[MatProperty::NUM_ELM];
        Float_t den[MatProperty::NUM_ELM];
        std::fill_n(elm, 9, false);
        std::fill_n(den, 9, 0.0);
        for (Long64_t it = cur_len_; it < max_len_; ++it) fill(elm, den);
    }
    
    if (cnt_ > 0) {
        Double_t inv_rad_len = 0.;
        for (Int_t it = 0; it < MatProperty::NUM_ELM; ++it) {
            if (!elm_.at(it)) continue;
            den_.at(it) = (den_.at(it) / static_cast<Double_t>(cnt_));
            geo_box_->elm[it] = elm_.at(it);
            geo_box_->den[it] = den_.at(it);
        
            inv_rad_len += (den_.at(it) * MatProperty::MASS[it] / MatProperty::RAD_LEN[it]);
        }
        geo_box_->inv_rad_len = inv_rad_len;
    }

    munmap(file_ptr_, file_len_);
    close(file_des_);

    clear();
}


Bool_t MatGeoBoxReader::load(const std::string& file_path) {
    if (is_load_) {
        if (file_path_ == file_path) return is_load_;
        else {
            MGSys::ShowWarning("MatGeoBoxReader::Load() : re-load different file.");
            clear();
        }
    }

    Int_t file_des = open(file_path.c_str(), O_RDONLY);
    Int_t file_len = lseek(file_des, 0, SEEK_END); 
    if (file_des < 0) {
        MGSys::ShowError(STR_FMT("MatGeoBoxReader::Load() : Matnetic field map not found (%s)", file_path.c_str()));
        is_load_ = false;
        return is_load_;
    }
    
    void* file_ptr = mmap(nullptr, file_len, PROT_READ, MAP_SHARED, file_des, 0);
    if (file_ptr == reinterpret_cast<void*>(-1)) {
        MGSys::ShowError("MatGeoBoxReader::Load() : mmap() failure.");
        close(file_des);
        is_load_ = false;
        return is_load_;
    }
    file_ptr_ = file_ptr;

    MatGeoBox* geo_box = reinterpret_cast<MatGeoBox*>(file_ptr_);
    for (Int_t ig = 0; ig < DIM_; ++ig) {
        n_.at(ig)   = geo_box->n[ig]; 
        min_.at(ig) = geo_box->min[ig];
        max_.at(ig) = geo_box->max[ig];
        len_.at(ig) = (max_.at(ig) - min_.at(ig));
        dlt_.at(ig) = (len_.at(ig) / static_cast<Double_t>(n_.at(ig)));
        elm_.at(ig) = geo_box->elm[ig];
        den_.at(ig) = geo_box->den[ig];
    }
    fact_.at(0) = n_.at(1) * n_.at(2);
    fact_.at(1) = n_.at(2);
    inv_rad_len_ = geo_box->inv_rad_len;
    mat_ptr_ = &(geo_box->mat);
    
    is_load_ = true;
    file_path_ = file_path;
    COUT("MatGeoBoxReader::Load() : Open file (%s)\n", file_path.c_str());
    return is_load_;
}

        
SVecD<3> MatGeoBoxReader::get_box_coord(const SVecD<3>& coo) {
    if (!is_load_) return SVecD<3>();
    Double_t xbox = (coo[0] - min_.at(0)) / len_.at(0);
    Double_t ybox = (coo[1] - min_.at(1)) / len_.at(1);
    Double_t zbox = (coo[2] - min_.at(2)) / len_.at(2);
    return SVecD<3>(xbox, ybox, zbox);
}

SVecD<3> MatGeoBoxReader::get_loc_coord(const SVecD<3>& coo) {
    if (!is_load_) return SVecD<3>();
    Double_t xloc = (coo[0] - min_.at(0)) / dlt_.at(0);
    Double_t yloc = (coo[1] - min_.at(1)) / dlt_.at(1);
    Double_t zloc = (coo[2] - min_.at(2)) / dlt_.at(2);
    return SVecD<3>(xloc, yloc, zloc);
}


std::tuple<Long64_t, Bool_t> MatGeoBoxReader::get_index(const SVecD<3>& loc) {
    Bool_t xin = (MGNumc::Compare(loc(0), MGMath::ZERO) > 0 && MGNumc::Compare(loc(0), static_cast<Double_t>(n_.at(0))) < 0);
    Bool_t yin = (MGNumc::Compare(loc(1), MGMath::ZERO) > 0 && MGNumc::Compare(loc(1), static_cast<Double_t>(n_.at(1))) < 0);
    Bool_t zin = (MGNumc::Compare(loc(2), MGMath::ZERO) > 0 && MGNumc::Compare(loc(2), static_cast<Double_t>(n_.at(2))) < 0);
    Bool_t inbox = (xin && yin && zin);
    if (!inbox) return std::make_tuple(-1, false);
    Long64_t xi = static_cast<Long64_t>(loc(0));
    Long64_t yi = static_cast<Long64_t>(loc(1));
    Long64_t zi = static_cast<Long64_t>(loc(2));
    Long64_t index = (xi * fact_.at(0) + yi * fact_.at(1) + zi);
    Bool_t   mat = mat_ptr_[index];
    return std::make_tuple(index, mat);
}
        

Bool_t MatGeoBoxReader::is_cross_box_coord(const SVecD<3>& vbox, const SVecD<3>& wbox) {
    if (MGNumc::Compare(vbox(2), MGMath::ZERO) <= 0 && MGNumc::Compare(wbox(2), MGMath::ZERO) <= 0) return false;
    if (MGNumc::Compare(vbox(2), MGMath::ONE ) >= 0 && MGNumc::Compare(wbox(2), MGMath::ONE ) >= 0) return false;
    if (MGNumc::Compare(vbox(1), MGMath::ZERO) <= 0 && MGNumc::Compare(wbox(1), MGMath::ZERO) <= 0) return false;
    if (MGNumc::Compare(vbox(1), MGMath::ONE ) >= 0 && MGNumc::Compare(wbox(1), MGMath::ONE ) >= 0) return false;
    if (MGNumc::Compare(vbox(0), MGMath::ZERO) <= 0 && MGNumc::Compare(wbox(0), MGMath::ZERO) <= 0) return false;
    if (MGNumc::Compare(vbox(0), MGMath::ONE ) >= 0 && MGNumc::Compare(wbox(0), MGMath::ONE ) >= 0) return false;
    return true;
}


Bool_t MatGeoBoxAms::CreateMatGeoBoxFromG4MatTree() {
    std::string g4mat_file_path = "/afs/cern.ch/work/h/hchou/public/DATABASE/detector/g4mat.root";
    TFile * root_file = TFile::Open(g4mat_file_path.c_str());
    if (root_file == nullptr || root_file->IsZombie()) return false;
    
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

    std::string dir_path = "/afs/cern.ch/work/h/hchou/public/DATABASE/detector/material";
    
    MatGeoBoxCreator creator_AMS02TRL1(
            MatAms::TRL1_N.at(0), MatAms::TRL1_MIN.at(0), MatAms::TRL1_MAX.at(0),
            MatAms::TRL1_N.at(1), MatAms::TRL1_MIN.at(1), MatAms::TRL1_MAX.at(1),
            MatAms::TRL1_N.at(2), MatAms::TRL1_MIN.at(2), MatAms::TRL1_MAX.at(2),
            CSTR_FMT("%s/AMS02TRL1.bin", dir_path.c_str())
        );
    
    MatGeoBoxCreator creator_AMS02TRD(
            MatAms::TRD_N.at(0), MatAms::TRD_MIN.at(0), MatAms::TRD_MAX.at(0),
            MatAms::TRD_N.at(1), MatAms::TRD_MIN.at(1), MatAms::TRD_MAX.at(1),
            MatAms::TRD_N.at(2), MatAms::TRD_MIN.at(2), MatAms::TRD_MAX.at(2),
            CSTR_FMT("%s/AMS02TRD.bin", dir_path.c_str())
        );
    
    MatGeoBoxCreator creator_AMS02UTOF(
            MatAms::UTOF_N.at(0), MatAms::UTOF_MIN.at(0), MatAms::UTOF_MAX.at(0),
            MatAms::UTOF_N.at(1), MatAms::UTOF_MIN.at(1), MatAms::UTOF_MAX.at(1),
            MatAms::UTOF_N.at(2), MatAms::UTOF_MIN.at(2), MatAms::UTOF_MAX.at(2),
            CSTR_FMT("%s/AMS02UTOF.bin", dir_path.c_str())
        );
    
    MatGeoBoxCreator creator_AMS02TRL2(
            MatAms::TRL2_N.at(0), MatAms::TRL2_MIN.at(0), MatAms::TRL2_MAX.at(0),
            MatAms::TRL2_N.at(1), MatAms::TRL2_MIN.at(1), MatAms::TRL2_MAX.at(1),
            MatAms::TRL2_N.at(2), MatAms::TRL2_MIN.at(2), MatAms::TRL2_MAX.at(2),
            CSTR_FMT("%s/AMS02TRL2.bin", dir_path.c_str())
        );
    
    MatGeoBoxCreator creator_AMS02TRL3(
            MatAms::TRL3_N.at(0), MatAms::TRL3_MIN.at(0), MatAms::TRL3_MAX.at(0),
            MatAms::TRL3_N.at(1), MatAms::TRL3_MIN.at(1), MatAms::TRL3_MAX.at(1),
            MatAms::TRL3_N.at(2), MatAms::TRL3_MIN.at(2), MatAms::TRL3_MAX.at(2),
            CSTR_FMT("%s/AMS02TRL3.bin", dir_path.c_str())
        );
    
    MatGeoBoxCreator creator_AMS02TRL4(
            MatAms::TRL4_N.at(0), MatAms::TRL4_MIN.at(0), MatAms::TRL4_MAX.at(0),
            MatAms::TRL4_N.at(1), MatAms::TRL4_MIN.at(1), MatAms::TRL4_MAX.at(1),
            MatAms::TRL4_N.at(2), MatAms::TRL4_MIN.at(2), MatAms::TRL4_MAX.at(2),
            CSTR_FMT("%s/AMS02TRL4.bin", dir_path.c_str())
        );
    
    MatGeoBoxCreator creator_AMS02TRL5(
            MatAms::TRL5_N.at(0), MatAms::TRL5_MIN.at(0), MatAms::TRL5_MAX.at(0),
            MatAms::TRL5_N.at(1), MatAms::TRL5_MIN.at(1), MatAms::TRL5_MAX.at(1),
            MatAms::TRL5_N.at(2), MatAms::TRL5_MIN.at(2), MatAms::TRL5_MAX.at(2),
            CSTR_FMT("%s/AMS02TRL5.bin", dir_path.c_str())
        );
    
    MatGeoBoxCreator creator_AMS02TRL6(
            MatAms::TRL6_N.at(0), MatAms::TRL6_MIN.at(0), MatAms::TRL6_MAX.at(0),
            MatAms::TRL6_N.at(1), MatAms::TRL6_MIN.at(1), MatAms::TRL6_MAX.at(1),
            MatAms::TRL6_N.at(2), MatAms::TRL6_MIN.at(2), MatAms::TRL6_MAX.at(2),
            CSTR_FMT("%s/AMS02TRL6.bin", dir_path.c_str())
        );
    
    MatGeoBoxCreator creator_AMS02TRL7(
            MatAms::TRL7_N.at(0), MatAms::TRL7_MIN.at(0), MatAms::TRL7_MAX.at(0),
            MatAms::TRL7_N.at(1), MatAms::TRL7_MIN.at(1), MatAms::TRL7_MAX.at(1),
            MatAms::TRL7_N.at(2), MatAms::TRL7_MIN.at(2), MatAms::TRL7_MAX.at(2),
            CSTR_FMT("%s/AMS02TRL7.bin", dir_path.c_str())
        );
    
    MatGeoBoxCreator creator_AMS02TRL8(
            MatAms::TRL8_N.at(0), MatAms::TRL8_MIN.at(0), MatAms::TRL8_MAX.at(0),
            MatAms::TRL8_N.at(1), MatAms::TRL8_MIN.at(1), MatAms::TRL8_MAX.at(1),
            MatAms::TRL8_N.at(2), MatAms::TRL8_MIN.at(2), MatAms::TRL8_MAX.at(2),
            CSTR_FMT("%s/AMS02TRL8.bin", dir_path.c_str())
        );
    
    MatGeoBoxCreator creator_AMS02LTOF(
            MatAms::LTOF_N.at(0), MatAms::LTOF_MIN.at(0), MatAms::LTOF_MAX.at(0),
            MatAms::LTOF_N.at(1), MatAms::LTOF_MIN.at(1), MatAms::LTOF_MAX.at(1),
            MatAms::LTOF_N.at(2), MatAms::LTOF_MIN.at(2), MatAms::LTOF_MAX.at(2),
            CSTR_FMT("%s/AMS02LTOF.bin", dir_path.c_str())
        );

    MatGeoBoxCreator creator_AMS02RICH(
            MatAms::RICH_N.at(0), MatAms::RICH_MIN.at(0), MatAms::RICH_MAX.at(0),
            MatAms::RICH_N.at(1), MatAms::RICH_MIN.at(1), MatAms::RICH_MAX.at(1),
            MatAms::RICH_N.at(2), MatAms::RICH_MIN.at(2), MatAms::RICH_MAX.at(2),
            CSTR_FMT("%s/AMS02RICH.bin", dir_path.c_str())
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

    for (Long64_t entry = 0; entry < tree_elm->GetEntries(); ++entry) {
        tree_elm->GetEntry(entry);
        Bool_t is_in_magnetic = (MGNumc::Compare(std::sqrt(coo[0] * coo[0] + coo[1] * coo[1]), MatAms::MAGNETIC_RADIUS) < 0);

        if (creator_AMS02TRL1.is_in_box(coo)) creator_AMS02TRL1.fill(elm, den);
        if (creator_AMS02TRD .is_in_box(coo)) creator_AMS02TRD .fill(elm, den);
        if (creator_AMS02UTOF.is_in_box(coo)) creator_AMS02UTOF.fill(elm, den);
        if (creator_AMS02TRL2.is_in_box(coo)) creator_AMS02TRL2.fill(elm, den);
        if (creator_AMS02TRL3.is_in_box(coo)) creator_AMS02TRL3.fill(elm, den, is_in_magnetic);
        if (creator_AMS02TRL4.is_in_box(coo)) creator_AMS02TRL4.fill(elm, den, is_in_magnetic);
        if (creator_AMS02TRL5.is_in_box(coo)) creator_AMS02TRL5.fill(elm, den, is_in_magnetic);
        if (creator_AMS02TRL6.is_in_box(coo)) creator_AMS02TRL6.fill(elm, den, is_in_magnetic);
        if (creator_AMS02TRL7.is_in_box(coo)) creator_AMS02TRL7.fill(elm, den, is_in_magnetic);
        if (creator_AMS02TRL8.is_in_box(coo)) creator_AMS02TRL8.fill(elm, den, is_in_magnetic);
        if (creator_AMS02LTOF.is_in_box(coo)) creator_AMS02LTOF.fill(elm, den);
        if (creator_AMS02RICH.is_in_box(coo)) creator_AMS02RICH.fill(elm, den);
        if (creator_AMS02PMT .is_in_box(coo)) creator_AMS02PMT .fill(elm, den);
        if (creator_AMS02TRL9.is_in_box(coo)) creator_AMS02TRL9.fill(elm, den);
        if (creator_AMS02ECAL.is_in_box(coo)) creator_AMS02ECAL.fill(elm, den);
    }
    
    creator_AMS02TRL1.save_and_close();
    creator_AMS02TRD .save_and_close();
    creator_AMS02UTOF.save_and_close();
    creator_AMS02TRL2.save_and_close();
    creator_AMS02TRL3.save_and_close();
    creator_AMS02TRL4.save_and_close();
    creator_AMS02TRL5.save_and_close();
    creator_AMS02TRL6.save_and_close();
    creator_AMS02TRL7.save_and_close();
    creator_AMS02TRL8.save_and_close();
    creator_AMS02LTOF.save_and_close();
    creator_AMS02RICH.save_and_close();
    creator_AMS02PMT .save_and_close();
    creator_AMS02TRL9.save_and_close();
    creator_AMS02ECAL.save_and_close();

    root_file->Close();

    return true;
}


Bool_t MatGeoBoxAms::Load() {
    if (is_load_) return is_load_;
    std::string g4mat_dir_path = "/afs/cern.ch/work/h/hchou/public/DATABASE/detector/material";

    reader_AMS02TRL1_.load(STR_FMT("%s/AMS02TRL1.bin",  g4mat_dir_path.c_str()));
    reader_AMS02TRD_ .load(STR_FMT("%s/AMS02TRD.bin",   g4mat_dir_path.c_str()));
    reader_AMS02UTOF_.load(STR_FMT("%s/AMS02UTOF.bin" , g4mat_dir_path.c_str()));
    reader_AMS02TRL2_.load(STR_FMT("%s/AMS02TRL2.bin" , g4mat_dir_path.c_str()));
    reader_AMS02TRL3_.load(STR_FMT("%s/AMS02TRL3.bin" , g4mat_dir_path.c_str()));
    reader_AMS02TRL4_.load(STR_FMT("%s/AMS02TRL4.bin" , g4mat_dir_path.c_str()));
    reader_AMS02TRL5_.load(STR_FMT("%s/AMS02TRL5.bin" , g4mat_dir_path.c_str()));
    reader_AMS02TRL6_.load(STR_FMT("%s/AMS02TRL6.bin" , g4mat_dir_path.c_str()));
    reader_AMS02TRL7_.load(STR_FMT("%s/AMS02TRL7.bin" , g4mat_dir_path.c_str()));
    reader_AMS02TRL8_.load(STR_FMT("%s/AMS02TRL8.bin" , g4mat_dir_path.c_str()));
    reader_AMS02LTOF_.load(STR_FMT("%s/AMS02LTOF.bin" , g4mat_dir_path.c_str()));
    reader_AMS02RICH_.load(STR_FMT("%s/AMS02RICH.bin" , g4mat_dir_path.c_str()));
    reader_AMS02PMT_ .load(STR_FMT("%s/AMS02PMT.bin",   g4mat_dir_path.c_str()));
    reader_AMS02TRL9_.load(STR_FMT("%s/AMS02TRL9.bin" , g4mat_dir_path.c_str()));
    reader_AMS02ECAL_.load(STR_FMT("%s/AMS02ECAL.bin" , g4mat_dir_path.c_str()));

    reader_.clear();
    reader_.push_back(&reader_AMS02TRL1_);
    reader_.push_back(&reader_AMS02TRD_ );
    reader_.push_back(&reader_AMS02UTOF_);
    reader_.push_back(&reader_AMS02TRL2_);
    reader_.push_back(&reader_AMS02TRL3_);
    reader_.push_back(&reader_AMS02TRL4_);
    reader_.push_back(&reader_AMS02TRL5_);
    reader_.push_back(&reader_AMS02TRL6_);
    reader_.push_back(&reader_AMS02TRL7_);
    reader_.push_back(&reader_AMS02TRL8_);
    reader_.push_back(&reader_AMS02LTOF_);
    reader_.push_back(&reader_AMS02RICH_);
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



#endif // __TRACKLibs_MatEnv_C__

