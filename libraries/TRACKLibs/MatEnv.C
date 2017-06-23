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
    cnt_.fill(0);
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
}


void MatGeoBoxCreator::fill(Bool_t elm[MatProperty::NUM_ELM], Double_t den[MatProperty::NUM_ELM]) {
    if (!is_open_) return;
    if (cur_len_ >= max_len_) { MGSys::ShowError("MatGeoBoxCreator::fill() : Out range."); return; }

    Bool_t has_mat = false;
    for (Int_t it = 0; it < MatProperty::NUM_ELM; ++it) {
        if (elm[it]) elm_.at(it) = true;
        else continue;
        den_.at(it) += den[it];
        cnt_.at(it)++;
        has_mat = true;
    }

    Bool_t& mat = *(static_cast<Bool_t*>(&(geo_box_->mat) + cur_len_));
    mat = has_mat;
    cur_len_++;
}


void MatGeoBoxCreator::save_and_close() {
    if (!is_open_) { clear(); return; }
    if (cur_len_ < max_len_) {
        MGSys::ShowWarning(STR_FMT("MatGeoBoxCreator::close() : fill not finished. (%lld/%lld)", cur_len_, max_len_));
        Bool_t   elm[MatProperty::NUM_ELM];
        Double_t den[MatProperty::NUM_ELM];
        std::fill_n(elm, 9, false);
        std::fill_n(den, 9, 0.0);
        for (Long64_t it = cur_len_; it < max_len_; ++it) fill(elm, den);
    }
    
    Long64_t mat_cnt = 0;
    for (Long64_t it = 0; it < max_len_; ++it) {
        Bool_t mat = *(static_cast<Bool_t*>(&(geo_box_->mat) + it));
        if (!mat) continue;
        mat_cnt++;
    }

    Double_t inv_rad_len = 0.;
    for (Int_t it = 0; it < MatProperty::NUM_ELM; ++it) {
        if (!elm_.at(it)) continue;
        den_.at(it) = (den_.at(it) / static_cast<Double_t>(mat_cnt));
        geo_box_->elm[it] = elm_.at(it);
        geo_box_->den[it] = den_.at(it);
    
        inv_rad_len += (den_.at(it) * MatProperty::MASS[it] / MatProperty::RAD_LEN[it]);
    }
    if (mat_cnt > 0) geo_box_->inv_rad_len = inv_rad_len;

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


} // namespace TrackSys



#endif // __TRACKLibs_MatEnv_C__

