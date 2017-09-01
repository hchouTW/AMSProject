#ifndef __TRACKLibs_MatEnv_C__
#define __TRACKLibs_MatEnv_C__


namespace TrackSys {


void MatFld::print() const {
    const Double_t scl = 1.0e+4;
    std::string printStr;
    printStr += STR_FMT("========================= MatFld =========================\n");
    printStr += STR_FMT("Mat       %-d\n", mat_);
    printStr += STR_FMT("InvRadLen %-8.6f\n", inv_rad_len_);
    printStr += STR_FMT("NumRadLen %-8.5f\n", num_rad_len());
    printStr += STR_FMT("RealLen   %-7.2f\n", real_len_);
    printStr += STR_FMT("EfftLen   %-7.2f\n", efft_len_);
    printStr += STR_FMT("Efft      %-6.4f\n", efft_);
    printStr += STR_FMT("Loc       %-6.4f\n", loc_);
    printStr += STR_FMT("Elm       H(%d) C(%d) N(%d) O(%d) F(%d) Na(%d) Al(%d) Si(%d) Pb(%d)\n", elm_.at(0), elm_.at(1), elm_.at(2), elm_.at(3), elm_.at(4), elm_.at(5), elm_.at(6), elm_.at(7), elm_.at(8));
    printStr += STR_FMT("Den [10^4 mol cm^-3]    H (%6.2f)  C (%6.2f)  N (%6.2f)\n", den_.at(0)*scl, den_.at(1)*scl, den_.at(2)*scl);
    printStr += STR_FMT("Den [10^4 mol cm^-3]    O (%6.2f)  F (%6.2f)  Na(%6.2f)\n", den_.at(3)*scl, den_.at(4)*scl, den_.at(5)*scl);
    printStr += STR_FMT("Den [10^4 mol cm^-3]    Al(%6.2f)  Si(%6.2f)  Pb(%6.2f)\n", den_.at(6)*scl, den_.at(7)*scl, den_.at(8)*scl);
    printStr += STR_FMT("==========================================================\n");
    COUT(printStr);
}
        

MatFld MatFld::Merge(const std::list<MatFld>& mflds) {
    if      (mflds.size() == 0) return MatFld();
    else if (mflds.size() == 1) return mflds.front();

    Bool_t                                     mat = false;
    std::array<Bool_t, MatProperty::NUM_ELM>   elm; elm.fill(false);
    std::array<Double_t, MatProperty::NUM_ELM> den; den.fill(0.);
    Double_t                                   inv_rad_len = 0.;
    Double_t                                   real_len = 0.;
    Double_t                                   efft_len = 0.;
    Double_t                                   efft = 0.;
    Double_t                                   loc = 0.;
    
    for (auto&& mfld : mflds) {
        if (!mfld()) { real_len += mfld.real_len(); continue; }
        
        for (Int_t it = 0; it < MatProperty::NUM_ELM; ++it) {
            if (!mfld.elm().at(it)) continue;
            elm.at(it) = true;
            den.at(it) = (mfld.efft_len() * mfld.den().at(it));
        }
        Double_t num_rad_len = (mfld.efft_len() * mfld.inv_rad_len());
        loc         += (real_len + mfld.loc() * mfld.real_len()) * num_rad_len;
        inv_rad_len += num_rad_len;
        efft_len += mfld.efft_len();
        real_len += mfld.real_len();

        mat = true;
    }

    if (mat && !MGNumc::EqualToZero(efft_len)) {
        for (Int_t it = 0; it < MatProperty::NUM_ELM; ++it) {
            if (!elm.at(it)) continue;
            den.at(it) = (den.at(it) / efft_len);
        }
        loc         = (loc / real_len / inv_rad_len);
        inv_rad_len = (inv_rad_len / efft_len);

        efft = (efft_len / real_len);
        
        return MatFld(mat, elm, den, inv_rad_len, real_len, efft_len, efft, loc);
    }
    else return MatFld(real_len);
}


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
    geo_box_   = reinterpret_cast<MatGeoBox*>(file_ptr);
    dlt_.fill(0);
    fact_.fill(0);
    cnt_ = 0;
    elm_.fill(false);
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
    dlt_.at(0) = (xmax - xmin) / static_cast<Double_t>(xn);
    dlt_.at(1) = (ymax - ymin) / static_cast<Double_t>(yn);
    dlt_.at(2) = (zmax - zmin) / static_cast<Double_t>(zn);
    fact_.at(0) = yn * zn;
    fact_.at(1) = zn;

    std::fill_n(geo_box_->elm, MatProperty::NUM_ELM, false);
    std::fill_n(geo_box_->den, MatProperty::NUM_ELM, 0.0);
    geo_box_->inv_rad_len = 0.;
    std::fill_n(static_cast<Bool_t*>(&(geo_box_->mat)), max_len_, false);
}


void MatGeoBoxCreator::fill(Float_t coo[3], Bool_t elm[MatProperty::NUM_ELM], Float_t den[MatProperty::NUM_ELM], Bool_t calculated) {
    if (!is_open_) return;
    if (!calculated) return;
    
    Long64_t zi = static_cast<Long64_t>(((static_cast<Double_t>(coo[2]) - geo_box_->min[2]) / dlt_.at(2)));
    if (zi < 0 || zi >= geo_box_->n[2]) return;
    
    Long64_t yi = static_cast<Long64_t>(((static_cast<Double_t>(coo[1]) - geo_box_->min[1]) / dlt_.at(1)));
    if (yi < 0 || yi >= geo_box_->n[1]) return;
    
    Long64_t xi = static_cast<Long64_t>(((static_cast<Double_t>(coo[0]) - geo_box_->min[0]) / dlt_.at(0)));
    if (xi < 0 || xi >= geo_box_->n[0]) return;
    
    Long64_t idx = (xi * fact_.at(0) + yi * fact_.at(1) + zi);

    Bool_t has_mat = false;
    for (Int_t it = 0; it < MatProperty::NUM_ELM; ++it) {
        if (!elm[it]) continue;
        has_mat = true;
        elm_.at(it) = true;
        den_.at(it) += static_cast<Double_t>(den[it]);
    }
    
    if (has_mat) {
        Bool_t& mat = *(static_cast<Bool_t*>(&(geo_box_->mat) + idx));
        mat = true;
        cnt_++;
    }
}


void MatGeoBoxCreator::save_and_close() {
    if (!is_open_) { clear(); return; }
    
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
        

void MatGeoBoxCreator::save_and_close(Bool_t elm[MatProperty::NUM_ELM], Float_t den[MatProperty::NUM_ELM]) {
    if (!is_open_) { clear(); return; }
    
    Bool_t has_mat = false;
    Double_t inv_rad_len = 0.;
    for (Int_t it = 0; it < MatProperty::NUM_ELM; ++it) {
        if (!elm[it]) continue;
        has_mat = true;
        elm_.at(it) = true;
        den_.at(it) = static_cast<Double_t>(den[it]);
        
        geo_box_->elm[it] = elm_.at(it);
        geo_box_->den[it] = den_.at(it);
        
        inv_rad_len += (den_.at(it) * MatProperty::MASS[it] / MatProperty::RAD_LEN[it]);
    }

    if (has_mat) {
        std::fill_n(static_cast<Bool_t*>(&(geo_box_->mat)), max_len_, true);
        geo_box_->inv_rad_len = inv_rad_len;
    }
    
    munmap(file_ptr_, file_len_);
    close(file_des_);

    clear();
}


void MatGeoBoxReader::print() const {
    if (!is_load_) return;
    std::string printStr;
    printStr += STR_FMT("===================== MatGeoBoxReader ====================\n");
    printStr += STR_FMT("BOX X     (%3d %7.2f %7.2f)\n", n_.at(0), min_.at(0), max_.at(0));
    printStr += STR_FMT("BOX Y     (%3d %7.2f %7.2f)\n", n_.at(1), min_.at(1), max_.at(1));
    printStr += STR_FMT("BOX Z     (%3d %7.2f %7.2f)\n", n_.at(2), min_.at(2), max_.at(2));
    printStr += STR_FMT("InvRadLen %-8.6f\n", inv_rad_len_);
    printStr += STR_FMT("Elm       H(%d) C(%d) N(%d) O(%d) F(%d) Na(%d) Al(%d) Si(%d) Pb(%d)\n", elm_.at(0), elm_.at(1), elm_.at(2), elm_.at(3), elm_.at(4), elm_.at(5), elm_.at(6), elm_.at(7), elm_.at(8));
    printStr += STR_FMT("==========================================================\n");
    COUT(printStr);
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
    }
    for (Int_t ie = 0; ie < MatProperty::NUM_ELM; ++ie) {
        elm_.at(ie) = geo_box->elm[ie];
        den_.at(ie) = geo_box->den[ie];
    }
    fact_.at(0) = n_.at(1) * n_.at(2);
    fact_.at(1) = n_.at(2);
    inv_rad_len_ = geo_box->inv_rad_len;
    mat_ptr_ = &(geo_box->mat);
    max_len_ = n_.at(0) * n_.at(1) * n_.at(2);
    
    is_load_ = true;
    file_path_ = file_path;
    COUT("MatGeoBoxReader::Load() : Open file (%s)\n", file_path.c_str());
    return is_load_;
}


Bool_t MatGeoBoxReader::is_in_box(const SVecD<3>& coo) {
    if (!is_load_) return false;
    if (MGNumc::Compare(coo(2), min_.at(2)) <= 0 || MGNumc::Compare(coo(2), max_.at(2)) >= 0) return false;
    if (MGNumc::Compare(coo(1), min_.at(1)) <= 0 || MGNumc::Compare(coo(1), max_.at(1)) >= 0) return false;
    if (MGNumc::Compare(coo(0), min_.at(0)) <= 0 || MGNumc::Compare(coo(0), max_.at(0)) >= 0) return false;
    return true;
}


Bool_t MatGeoBoxReader::is_cross(const SVecD<3>& vcoo, const SVecD<3>& wcoo) {
    if (!is_load_) return false;
    if (MGNumc::Compare(vcoo(2), min_.at(2)) <= 0 && MGNumc::Compare(wcoo(2), min_.at(2)) <= 0) return false;
    if (MGNumc::Compare(vcoo(2), max_.at(2)) >= 0 && MGNumc::Compare(wcoo(2), max_.at(2)) >= 0) return false;
    if (MGNumc::Compare(vcoo(1), min_.at(1)) <= 0 && MGNumc::Compare(wcoo(1), min_.at(1)) <= 0) return false;
    if (MGNumc::Compare(vcoo(1), max_.at(1)) >= 0 && MGNumc::Compare(wcoo(1), max_.at(1)) >= 0) return false;
    if (MGNumc::Compare(vcoo(0), min_.at(0)) <= 0 && MGNumc::Compare(wcoo(0), min_.at(0)) <= 0) return false;
    if (MGNumc::Compare(vcoo(0), max_.at(0)) >= 0 && MGNumc::Compare(wcoo(0), max_.at(0)) >= 0) return false;
    return true;
}
        

MatFld MatGeoBoxReader::get(const SVecD<3>& coo) {
    if (!is_load_) return MatFld();

    Long64_t zi = static_cast<Long64_t>(((coo(2) - min_.at(2)) / dlt_.at(2)));
    if (zi < 0 || zi >= n_.at(2)) return MatFld();
    
    Long64_t yi = static_cast<Long64_t>(((coo(1) - min_.at(1)) / dlt_.at(1)));
    if (yi < 0 || yi >= n_.at(1)) return MatFld();

    Long64_t xi = static_cast<Long64_t>(((coo(0) - min_.at(0)) / dlt_.at(0)));
    if (xi < 0 || xi >= n_.at(0)) return MatFld();

    Long64_t idx = (xi * fact_.at(0) + yi * fact_.at(1) + zi);
    Bool_t mat = mat_ptr_[idx];

    if (mat) return MatFld(mat, elm_, den_, inv_rad_len_);
    else     return MatFld();
}


MatFld MatGeoBoxReader::get(const SVecD<3>& vcoo, const SVecD<3>& wcoo, Bool_t is_std) {
    if (!is_load_) return MatFld();
    
    Double_t vzloc = std::move((vcoo(2) - min_.at(2)) / dlt_.at(2));
    Double_t wzloc = std::move((wcoo(2) - min_.at(2)) / dlt_.at(2));
    Long64_t vzi = static_cast<Long64_t>(vzloc);
    Long64_t wzi = static_cast<Long64_t>(wzloc);
    if ((vzi < 0 && wzi < 0) || (vzi >= n_.at(2) && wzi >= n_.at(2))) return MatFld();
    
    Double_t vyloc = std::move((vcoo(1) - min_.at(1)) / dlt_.at(1));
    Double_t wyloc = std::move((wcoo(1) - min_.at(1)) / dlt_.at(1));
    Long64_t vyi = static_cast<Long64_t>(vyloc);
    Long64_t wyi = static_cast<Long64_t>(wyloc);
    if ((vyi < 0 && wyi < 0) || (vyi >= n_.at(1) && wyi >= n_.at(1))) return MatFld();
    
    Double_t vxloc = std::move((vcoo(0) - min_.at(0)) / dlt_.at(0));
    Double_t wxloc = std::move((wcoo(0) - min_.at(0)) / dlt_.at(0));
    Long64_t vxi = static_cast<Long64_t>(vxloc);
    Long64_t wxi = static_cast<Long64_t>(wxloc);
    if ((vxi < 0 && wxi < 0) || (vxi >= n_.at(0) && wxi >= n_.at(0))) return MatFld();

    Double_t real_len = LA::Mag((wcoo - vcoo));

    if (MGNumc::EqualToZero(real_len)) {
        Long64_t idx = (vxi * fact_.at(0) + vyi * fact_.at(1) + vzi);
        Bool_t mat = mat_ptr_[idx];
        if (mat) return MatFld(mat, elm_, den_, inv_rad_len_);
        else     return MatFld();
    }
   
    SVecD<3> vwvec((wxloc - vxloc), (wyloc - vyloc), (wzloc - vzloc));
    
    Double_t   vwlen  = LA::Mag(vwvec);
    Long64_t   nstp   = static_cast<Long64_t>(vwlen / (is_std ? STD_STEP_LEN_ : FST_STEP_LEN_)) + 2;
    SVecD<3>&& unit   = (vwvec / static_cast<Double_t>(nstp));

    SVecD<3> itloc((vxloc + MGMath::HALF * unit(0)), (vyloc + MGMath::HALF * unit(1)), (vzloc + MGMath::HALF * unit(2)));

    Long64_t suml = 0;
    Long64_t itcnt = 0;
    for (Long64_t it = 0; it < nstp; ++it, itloc += unit) {
        Long64_t zi = static_cast<Long64_t>(itloc(2));
        if (zi < 0 || zi >= n_.at(2)) continue;
        Long64_t yi = static_cast<Long64_t>(itloc(1));
        if (yi < 0 || yi >= n_.at(1)) continue;
        Long64_t xi = static_cast<Long64_t>(itloc(0));
        if (xi < 0 || xi >= n_.at(0)) continue;

        Long64_t idx = (xi * fact_.at(0) + yi * fact_.at(1) + zi);
        if (mat_ptr_[idx]) {
            suml += it;
            itcnt++;
        }
    }

    if (itcnt != 0) {
        Double_t efft = static_cast<Double_t>(itcnt) / static_cast<Double_t>(nstp);
        Double_t efft_len = efft * real_len;
        Double_t loc = ((MGMath::HALF + static_cast<Double_t>(suml) / static_cast<Double_t>(itcnt)) / static_cast<Double_t>(nstp));

        return MatFld(true, elm_, den_, inv_rad_len_, real_len, efft_len, efft, loc);
    }
    else return MatFld(real_len);
}


#ifdef __HAS_TESTING__
Bool_t MatGeoBoxTest::CreateMatGeoBox() {
    std::string dir_path = "/data3/hchou/AMSData/MatTest";
    
    MatGeoBoxCreator creator_TRL1(
            MatTest::TRL1_N.at(0), MatTest::TRL1_MIN.at(0), MatTest::TRL1_MAX.at(0),
            MatTest::TRL1_N.at(1), MatTest::TRL1_MIN.at(1), MatTest::TRL1_MAX.at(1),
            MatTest::TRL1_N.at(2), MatTest::TRL1_MIN.at(2), MatTest::TRL1_MAX.at(2),
            CSTR_FMT("%s/TRL1.bin", dir_path.c_str())
        );
    
    MatGeoBoxCreator creator_TRL2(
            MatTest::TRL2_N.at(0), MatTest::TRL2_MIN.at(0), MatTest::TRL2_MAX.at(0),
            MatTest::TRL2_N.at(1), MatTest::TRL2_MIN.at(1), MatTest::TRL2_MAX.at(1),
            MatTest::TRL2_N.at(2), MatTest::TRL2_MIN.at(2), MatTest::TRL2_MAX.at(2),
            CSTR_FMT("%s/TRL2.bin", dir_path.c_str())
        );
    
    bool  tr_elm[9] = { 0, 0, 0, 0, 0, 0, 0, 1, 0 };
    float tr_den[9] = { 0, 0, 0, 0, 0, 0, 0, 0.083, 0 };
    creator_TRL1.save_and_close(tr_elm, tr_den);
    creator_TRL2.save_and_close(tr_elm, tr_den);

    return true;
}


Bool_t MatGeoBoxTest::Load() {
    if (is_load_) return is_load_;
    std::string g4mat_dir_path = "/data3/hchou/AMSData/MatTest";

    reader_TRL1_.load(STR_FMT("%s/TRL1.bin", g4mat_dir_path.c_str()));
    reader_TRL2_.load(STR_FMT("%s/TRL2.bin", g4mat_dir_path.c_str()));
    
    reader_.push_back(&reader_TRL1_);
    reader_.push_back(&reader_TRL2_);
   
    is_load_ = true;
    for (auto&& reader : reader_) {
        if (!reader->exist()) {
            is_load_ = false;
            break;
        }
    }

    return is_load_;
}
#endif // __HAS_TESTING__


#ifdef __HAS_AMS_OFFICE_LIBS__
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
#endif // __HAS_AMS_OFFICE_LIBS__


Bool_t MatMgnt::Load() {
    if (is_load_ && reader_ != nullptr) return true;
    is_load_ = false;
    reader_ = nullptr;

#ifdef __HAS_TESTING__
    is_load_ = MatGeoBoxTest::Load();
    if (is_load_) reader_ = &MatGeoBoxTest::Reader();
#endif // __HAS_TESTING__

#ifdef __HAS_AMS_OFFICE_LIBS__
    is_load_ = MatGeoBoxAms::Load();
    if (is_load_) reader_ = &MatGeoBoxAms::Reader();
#endif // __HAS_AMS_OFFICE_LIBS__
    
    return is_load_;
};


MatFld MatMgnt::Get(const SVecD<3>& coo) {
    if (!Load()) return MatFld();

    Bool_t mat = false;
    std::array<Bool_t,   MatProperty::NUM_ELM> elm; elm.fill(false);
    std::array<Double_t, MatProperty::NUM_ELM> den; den.fill(0.);
    Double_t                                   inv_rad_len = 0.0;

    for (auto&& reader : *reader_) {
        if (!reader->is_in_box(coo)) continue;
        MatFld&& mfld = reader->get(coo);
        if (!mfld()) continue;

        for (Int_t it = 0; it < MatProperty::NUM_ELM; ++it) {
            if (!mfld.elm().at(it)) continue;
            elm.at(it) = true;
            den.at(it) += mfld.den().at(it);
        }
        inv_rad_len += mfld.inv_rad_len();
        mat = true;
    }

    if (mat) return MatFld(mat, elm, den, inv_rad_len);
    else     return MatFld();
}
    

MatFld MatMgnt::Get(const SVecD<3>& vcoo, const SVecD<3>& wcoo, Bool_t is_std) {
    if (!Load()) return MatFld();
    
    Double_t real_len = LA::Mag(wcoo - vcoo);
    if (MGNumc::EqualToZero(real_len)) return Get(vcoo);
    
    Bool_t                                     mat = false;
    std::array<Bool_t,   MatProperty::NUM_ELM> elm; elm.fill(false);
    std::array<Double_t, MatProperty::NUM_ELM> den; den.fill(0.);
    Double_t                                   inv_rad_len = 0.0;
    Double_t                                   efft_len = 0.0;
    Double_t                                   efft = 0.0;
    Double_t                                   loc = 0.0;

    for (auto&& reader : *reader_) {
        if (!reader->is_cross(vcoo, wcoo)) continue;
        MatFld&& mfld = reader->get(vcoo, wcoo, is_std);
        if (!mfld()) continue;
        
        for (Int_t it = 0; it < MatProperty::NUM_ELM; ++it) {
            if (!mfld.elm().at(it)) continue;
            elm.at(it) = true;
            den.at(it) += (mfld.efft_len() * mfld.den().at(it));
        }
        Double_t num_rad_len = (mfld.efft_len() * mfld.inv_rad_len());
        loc         += mfld.loc() * num_rad_len;
        inv_rad_len += num_rad_len;
        efft_len    += mfld.efft_len();
        
        mat = true;
    }

    if (mat && !MGNumc::EqualToZero(efft_len)) {
        for (Int_t it = 0; it < MatProperty::NUM_ELM; ++it) {
            if (!elm.at(it)) continue;
            den.at(it) = (den.at(it) / efft_len);
        }
        loc         = (loc / inv_rad_len);
        inv_rad_len = (inv_rad_len / efft_len);
        efft = (efft_len / real_len);

        return MatFld(mat, elm, den, inv_rad_len, real_len, efft_len, efft, loc);
    }
    else return MatFld(real_len);
}


MatFld MatMgnt::Get(Double_t stp_len, const PhySt& part, Bool_t is_std) {
    const SVecD<3>&  vcoo = part.coo();
    SVecD<3>&&       wcoo = part.coo() + stp_len * part.dir();

    return Get(vcoo, wcoo, is_std);
}


void MatArg::rndm(const MatFld& mfld) {
    if (!mfld()) {
        tau_mscat_ = MGMath::ZERO;
        rho_mscat_ = MGMath::ZERO;
        ion_eloss_ = MGMath::ZERO;
        brm_eloss_ = MGMath::ZERO;
        return;
    }

    if (sw_mscat_) {
        tau_mscat_ = MGRndm::NormalGaussian();
        rho_mscat_ = MGRndm::NormalGaussian();
    }
    if (sw_eloss_) {
        Int_t cnt = 0;
        const Int_t MAX_CNT = 3;
        const Double_t ION_LIMIT = -10.0;
        do {
            ion_eloss_ = Rndm::Landau();
            cnt++;
        } while(MGNumc::Compare(ion_eloss_, ION_LIMIT) < 0 && cnt <= MAX_CNT);
        
        Double_t bremslen = mfld.num_rad_len() / MGMath::LOG_TWO;
        brm_eloss_ = (MGNumc::Compare(bremslen) <= 0) ? MGMath::ZERO : MGRndm::Gamma(bremslen, (MGMath::ONE / bremslen))();
    }
}


void MatArg::rndm(const MatPhyFld& mphy) {
    if (!mphy()) {
        tau_mscat_ = MGMath::ZERO;
        rho_mscat_ = MGMath::ZERO;
        ion_eloss_ = MGMath::ZERO;
        brm_eloss_ = MGMath::ZERO;
        return;
    }

    if (sw_mscat_) {
        tau_mscat_ = MGRndm::NormalGaussian();
        rho_mscat_ = MGRndm::NormalGaussian();
    }
    if (sw_eloss_) {
        Int_t cnt = 0;
        const Int_t MAX_CNT = 3;
        Double_t ionlimit = MGMath::NEG * mphy.ion_eloss_mpv() / mphy.ion_eloss_sgm();
        do {
            ion_eloss_ = Rndm::Landau();
            cnt++;
        } while(MGNumc::Compare(ion_eloss_, ionlimit) < 0 && cnt <= MAX_CNT);
        
        Double_t bremslen = mphy.num_rad_len() / MGMath::LOG_TWO;
        brm_eloss_ = (MGNumc::Compare(bremslen) <= 0) ? MGMath::ZERO : MGRndm::Gamma(bremslen, (MGMath::ONE / bremslen))();
    }
}


Double_t MatPhy::GetNumRadLen(const Double_t stp_len, const PhySt& part, Bool_t is_std) {
    const SVecD<3>&  vcoo = part.coo();
    SVecD<3>&&       wcoo = part.coo() + stp_len * part.dir();

    MatFld&& mfld = MatMgnt::Get(vcoo, wcoo, is_std);
   
    return mfld.num_rad_len();
}


MatPhyFld MatPhy::Get(const Double_t stp_len, const PhySt& part, const MatArg& marg, Bool_t is_std) {
    if (!marg()) return MatPhyFld();
    if (part.part().is_chrgless() || part.part().is_massless()) return MatPhyFld();
    if (MGNumc::EqualToZero(stp_len)) return MatPhyFld();
    if (MGNumc::EqualToZero(part.mom())) return MatPhyFld();

    const SVecD<3>&  vcoo = part.coo();
    SVecD<3>&&       wcoo = part.coo() + stp_len * part.dir();

    MatFld&& mfld = MatMgnt::Get(vcoo, wcoo, is_std);
    
    if (!mfld()) return MatPhyFld();

    Double_t mult_scat_sgm = GetMultipleScattering(mfld, part);
    std::pair<Double_t, Double_t>&& ion_eloss = GetIonizationEnergyLoss(mfld, part);
    Double_t brm_eloss_men = GetBremsstrahlungEnergyLoss(mfld, part);

    return MatPhyFld(mfld(), mfld.inv_rad_len(), mfld.num_rad_len(), mult_scat_sgm, ion_eloss.first, ion_eloss.second, brm_eloss_men);
}


MatPhyFld MatPhy::Get(const MatFld& mfld, const PhySt& part, const MatArg& marg) {
    if (!mfld() || !marg()) return MatPhyFld();
    if (MGNumc::EqualToZero(mfld.efft_len())) return MatPhyFld();
    if (part.part().is_chrgless() || part.part().is_massless()) return MatPhyFld();
    if (MGNumc::EqualToZero(part.mom())) return MatPhyFld();
    
    Double_t mult_scat_sgm = GetMultipleScattering(mfld, part);
    std::pair<Double_t, Double_t>&& ion_eloss = GetIonizationEnergyLoss(mfld, part);
    Double_t brm_eloss_men = GetBremsstrahlungEnergyLoss(mfld, part);

    return MatPhyFld(mfld(), mfld.inv_rad_len(), mfld.num_rad_len(), mult_scat_sgm, ion_eloss.first, ion_eloss.second, brm_eloss_men);
}
        

void MatMscatFld::print() const {
    std::string printStr;
    printStr += STR_FMT("===== MatMscatFld =====\n");
    printStr += STR_FMT("Mat     %-d\n", mat_);
    printStr += STR_FMT("SgmLoc  %-9.6f\n", sgm_loc_);
    printStr += STR_FMT("SgmTha  %-9.6f\n", sgm_tha_);
    printStr += STR_FMT("CovLoc  %-9.6f %-9.6f\n", std::sqrt(cov_(0, 0)), std::sqrt(cov_(1, 1)));
    printStr += STR_FMT("CovTha  %-9.6f %-9.6f\n", std::sqrt(cov_(2, 2)), std::sqrt(cov_(3, 3)));
    printStr += STR_FMT("=======================\n");
    COUT(printStr);
}


MatMscatFld MatPhy::GetMscat(const MatFld& mfld, const PhySt& part) {
    if (!mfld()) return MatMscatFld();
    if (MGNumc::EqualToZero(mfld.efft_len())) return MatMscatFld();
    if (part.part().is_chrgless() || part.part().is_massless()) return MatMscatFld();
    if (MGNumc::EqualToZero(part.mom())) return MatMscatFld();
    
    Double_t num_rad_len = mfld.num_rad_len();
    if (MGNumc::Compare(num_rad_len, LMTL_NUM_RAD_LEN) < 0) num_rad_len = LMTL_NUM_RAD_LEN;
    if (MGNumc::Compare(num_rad_len, LMTU_NUM_RAD_LEN) > 0) num_rad_len = LMTU_NUM_RAD_LEN;
    
    Bool_t is_over_lmt = (MGNumc::Compare(part.bta(), LMT_BTA) > 0);
    Double_t eta       = ((is_over_lmt) ? part.eta_abs() : LMT_INV_GMBTA);
    Double_t eta_part  = (part.part().chrg_to_mass() * eta) * std::sqrt(eta * eta + MGMath::ONE);
    Double_t sgm_tha   = RYDBERG_CONST * eta_part * std::sqrt(num_rad_len) * (MGMath::ONE + NRL_CORR_FACT * std::log(num_rad_len));
    Double_t sgm_loc   = MGMath::INV_SQRT_TWELVE * mfld.efft_len() * sgm_tha;
    Double_t res_len   = (MGMath::ONE - mfld.loc()) * mfld.real_len();

    Double_t mscat_tha_tha = sgm_tha * sgm_tha;
    Double_t mscat_loc_loc = (sgm_loc * sgm_loc + res_len * res_len * mscat_tha_tha);
    Double_t mscat_loc_tha = res_len * mscat_tha_tha;

    Double_t dxx = (MGMath::ONE - part.dy() * part.dy());
    Double_t dyy = (MGMath::ONE - part.dx() * part.dx());
    Double_t dxy = (part.dx() * part.dy());

    SMtxSymD<5> mscat;
    mscat(0, 0) = mscat_loc_loc * dxx;
    mscat(1, 1) = mscat_loc_loc * dyy;
    mscat(0, 1) = mscat_loc_loc * dxy;
    mscat(2, 2) = mscat_tha_tha * dxx;
    mscat(3, 3) = mscat_tha_tha * dyy;
    mscat(2, 3) = mscat_tha_tha * dxy;
    mscat(0, 2) = mscat_loc_tha;
    mscat(1, 3) = mscat_loc_tha;

    return MatMscatFld(sgm_loc, sgm_tha, mscat);
}


std::array<Double_t, MatProperty::NUM_ELM> MatPhy::GetDensityEffectCorrection(const MatFld& mfld, const PhySt& part) {
    Double_t gmbta     = ((MGNumc::Compare(part.bta(), LMT_BTA) > 0) ? part.gmbta() : LMT_GMBTA);
    Double_t log_gmbta = std::log10(gmbta);

    std::array<Double_t, MatProperty::NUM_ELM> dlt; dlt.fill(0.);
    for (Int_t it = 0; it < MatProperty::NUM_ELM; ++it) {
       if (!mfld.elm().at(it)) continue;
       if (log_gmbta < MatProperty::DEN_EFF_CORR_X0.at(it)) continue; // if nonconductors
       Double_t dif_log_gmbta = MatProperty::DEN_EFF_CORR_X1.at(it) - log_gmbta;
       dlt.at(it) = MGMath::TWO * MGMath::LOG_TEN * log_gmbta - MatProperty::DEN_EFF_CORR_C.at(it);
       if (MGNumc::Compare(dif_log_gmbta) <= 0) continue;
       dlt.at(it) += MatProperty::DEN_EFF_CORR_A.at(it) * std::pow(dif_log_gmbta, MatProperty::DEN_EFF_CORR_M.at(it));
       if (MGNumc::Compare(dlt.at(it)) <= 0) dlt.at(it) = MGMath::ZERO;
    }

    return dlt;
}


Double_t MatPhy::GetMultipleScattering(const MatFld& mfld, const PhySt& part) {
    Double_t num_rad_len = mfld.num_rad_len();
    if (MGNumc::Compare(num_rad_len, LMTL_NUM_RAD_LEN) < 0) num_rad_len = LMTL_NUM_RAD_LEN;
    if (MGNumc::Compare(num_rad_len, LMTU_NUM_RAD_LEN) > 0) num_rad_len = LMTU_NUM_RAD_LEN;
    Bool_t is_over_lmt = (MGNumc::Compare(part.bta(), LMT_BTA) > 0);
    Double_t eta = ((is_over_lmt) ? part.eta_abs() : LMT_INV_GMBTA);
    Double_t eta_part = std::sqrt(eta * eta + MGMath::ONE);
    Double_t mscat_sgm = RYDBERG_CONST * part.part().chrg_to_mass() * eta_part * std::sqrt(num_rad_len) * (MGMath::ONE + NRL_CORR_FACT * std::log(num_rad_len)) * (eta / mfld.real_len());
    return mscat_sgm;
}


std::pair<Double_t, Double_t> MatPhy::GetIonizationEnergyLoss(const MatFld& mfld, const PhySt& part) {
    Bool_t is_over_lmt = (MGNumc::Compare(part.bta(), LMT_BTA) > 0);
    Double_t sqr_gmbta = ((is_over_lmt) ? (part.gmbta() * part.gmbta()) : LMT_SQR_GMBTA);
    Double_t sqr_bta   = ((is_over_lmt) ? (part.bta() * part.bta()) : LMT_SQR_BTA);
    Double_t sqr_chrg  = part.part().chrg() * part.part().chrg();
    Double_t mass      = part.part().mass() * GEV_TO_MEV;
        
    std::array<Double_t, MatProperty::NUM_ELM>&& dlt = GetDensityEffectCorrection(mfld, part);
    
    Double_t avg_dlt = MGMath::ZERO;
    Double_t m_exeng = MGMath::ZERO;
    Double_t ttl_wgt = MGMath::ZERO;
    for (Int_t it = 0; it < MatProperty::NUM_ELM; ++it) {
       if (!mfld.elm().at(it)) continue;
       Double_t wgt = MatProperty::CHRG.at(it) * mfld.den().at(it);
       m_exeng += wgt * MatProperty::NEG_LN_MEAN_EXENG.at(it);
       avg_dlt += wgt * dlt.at(it);
       ttl_wgt += wgt;
    }
    avg_dlt /= ttl_wgt;
    m_exeng /= ttl_wgt;

    if (!MGNumc::Valid(ttl_wgt) || MGNumc::Compare(ttl_wgt) <= 0) ttl_wgt = MGMath::ZERO;
    if (!MGNumc::Valid(avg_dlt) || MGNumc::Compare(avg_dlt) <= 0) avg_dlt = MGMath::ZERO;
    if (!MGNumc::Valid(m_exeng) || MGNumc::Compare(m_exeng) <= 0) m_exeng = MGMath::ZERO;
  
    Double_t el_abs   = (ttl_wgt * mfld.efft_len());
    Double_t eng_fact = (MGMath::HALF * BETHE_BLOCH_K * el_abs * sqr_chrg / sqr_bta); 
    Double_t eta_part = (std::sqrt(sqr_gmbta + MGMath::ONE) / sqr_gmbta);
    Double_t sgm_fact = (eng_fact * eta_part / mass / mfld.real_len());

    Double_t exeng    = std::exp(-m_exeng);
    Double_t max_keng = std::log(MGMath::TWO * MASS_EL_IN_MEV * sqr_gmbta / exeng);
    Double_t ion_keng = std::log(eng_fact / exeng);
    
    if (!MGNumc::Valid(max_keng) || MGNumc::Compare(max_keng) <= 0) max_keng = MGMath::ZERO;
    if (!MGNumc::Valid(ion_keng) || MGNumc::Compare(ion_keng) <= 0) ion_keng = MGMath::ZERO;

    Double_t ion_eloss_sgm = sgm_fact;
    Double_t ion_eloss_mpv = sgm_fact * (max_keng + ion_keng + LANDAU_ELOSS_CORR - sqr_bta - avg_dlt);
    
    if (!MGNumc::Valid(ion_eloss_sgm) || MGNumc::Compare(ion_eloss_sgm) <= 0) ion_eloss_sgm = MGMath::ZERO;
    if (!MGNumc::Valid(ion_eloss_mpv) || MGNumc::Compare(ion_eloss_mpv) <= 0) ion_eloss_mpv = MGMath::ZERO;

    return std::make_pair(ion_eloss_sgm, ion_eloss_mpv);
}


Double_t MatPhy::GetBremsstrahlungEnergyLoss(const MatFld& mfld, const PhySt& part) {
    Bool_t   is_over_lmt  = (MGNumc::Compare(part.bta(), LMT_BTA) > 0);
    Double_t inv_sqr_bta  = ((is_over_lmt) ? (MGMath::ONE / part.bta() / part.bta()) : LMT_INV_SQR_BTA);
    Double_t sqr_chrg_rat = (part.part().chrg() * part.part().chrg());
    Double_t sqr_mass_rat = (MASS_EL_IN_GEV * MASS_EL_IN_GEV / part.part().mass() / part.part().mass());
    Double_t brm_eloss_men = sqr_chrg_rat * sqr_mass_rat * inv_sqr_bta * (mfld.num_rad_len() / MGMath::LOG_TWO) / mfld.real_len();
    if (!MGNumc::Valid(brm_eloss_men) || MGNumc::Compare(brm_eloss_men) <= 0) brm_eloss_men = MGMath::ZERO;
    return brm_eloss_men;
}


} // namespace TrackSys



#endif // __TRACKLibs_MatEnv_C__

