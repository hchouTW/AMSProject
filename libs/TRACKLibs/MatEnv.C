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
    printStr += STR_FMT("LocSqr    %-6.4f\n", locsqr_);
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
    Double_t                                   locsqr = 0.;
    
    for (auto&& mfld : mflds) {
        if (!mfld()) { real_len += mfld.real_len(); continue; }
        
        for (Int_t it = 0; it < MatProperty::NUM_ELM; ++it) {
            if (!mfld.elm().at(it)) continue;
            elm.at(it)  = true;
            den.at(it) += (mfld.efft_len() * mfld.den().at(it));
        }
        Double_t loclen    = mfld.loc() * mfld.real_len();
        Double_t locsqrlen = mfld.locsqr() * mfld.real_len() * mfld.real_len();
        Double_t num_rad_len = mfld.num_rad_len();
        loc         += (real_len + loclen) * num_rad_len;
        locsqr      += (real_len * real_len + locsqrlen + MGMath::TWO * real_len * loclen) * num_rad_len;
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
        loc         = ((loc / (real_len)) / inv_rad_len);
        locsqr      = ((locsqr / (real_len * real_len)) / inv_rad_len);
        inv_rad_len = (inv_rad_len / efft_len);

        efft = (efft_len / real_len);
        
        return MatFld(mat, elm, den, inv_rad_len, real_len, efft_len, efft, loc, locsqr);
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
    vol_ = 0.;
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
    vol_ = (dlt_.at(0) * dlt_.at(1) * dlt_.at(2));

    std::fill_n(geo_box_->elm, MatProperty::NUM_ELM, false);
    std::fill_n(geo_box_->den, MatProperty::NUM_ELM, 0.0);
    geo_box_->inv_rad_len = 0.;
    std::fill_n(static_cast<Bool_t*>(&(geo_box_->mat)), max_len_, false);
}


void MatGeoBoxCreator::fill(Double_t coo[3], Bool_t elm[MatProperty::NUM_ELM], Double_t mol[MatProperty::NUM_ELM], Bool_t calculated) {
    if (!is_open_) return;
    if (!calculated) return;
    
    Long64_t zi = static_cast<Long64_t>(std::floor((coo[2] - geo_box_->min[2]) / dlt_.at(2)));
    if (zi < 0 || zi >= geo_box_->n[2]) return;
    
    Long64_t yi = static_cast<Long64_t>(std::floor((coo[1] - geo_box_->min[1]) / dlt_.at(1)));
    if (yi < 0 || yi >= geo_box_->n[1]) return;
    
    Long64_t xi = static_cast<Long64_t>(std::floor((coo[0] - geo_box_->min[0]) / dlt_.at(0)));
    if (xi < 0 || xi >= geo_box_->n[0]) return;
    
    Long64_t idx = (xi * fact_.at(0) + yi * fact_.at(1) + zi);
    Bool_t&  mat = *(static_cast<Bool_t*>(&(geo_box_->mat) + idx));

    Bool_t has_mat = false;
    for (Int_t it = 0; it < MatProperty::NUM_ELM; ++it) {
        if (!elm[it]) continue;
        has_mat = true;
        elm_.at(it) = true;
        den_.at(it) += mol[it];
    }
    
    if (!mat && has_mat) {
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
            den_.at(it) = (den_.at(it) / (vol_ * static_cast<Double_t>(cnt_)));
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
        

void MatGeoBoxCreator::save_and_close(Bool_t elm[MatProperty::NUM_ELM], Double_t den[MatProperty::NUM_ELM]) {
    if (!is_open_) { clear(); return; }
    
    Bool_t has_mat = false;
    Double_t inv_rad_len = 0.;
    for (Int_t it = 0; it < MatProperty::NUM_ELM; ++it) {
        if (!elm[it]) continue;
        has_mat = true;
        elm_.at(it) = true;
        den_.at(it) = den[it];
        
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
        MGSys::ShowError(STR_FMT("MatGeoBoxReader::Load() : Mat field map not found (%s)", file_path.c_str()));
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

    Long64_t zi = static_cast<Long64_t>(std::floor((coo(2) - min_.at(2)) / dlt_.at(2)));
    if (zi < 0 || zi >= n_.at(2)) return MatFld();
    
    Long64_t yi = static_cast<Long64_t>(std::floor((coo(1) - min_.at(1)) / dlt_.at(1)));
    if (yi < 0 || yi >= n_.at(1)) return MatFld();

    Long64_t xi = static_cast<Long64_t>(std::floor((coo(0) - min_.at(0)) / dlt_.at(0)));
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
    Long64_t vzi = static_cast<Long64_t>(std::floor(vzloc));
    Long64_t wzi = static_cast<Long64_t>(std::floor(wzloc));
    if ((vzi < 0 && wzi < 0) || (vzi >= n_.at(2) && wzi >= n_.at(2))) return MatFld();
    
    Double_t vyloc = std::move((vcoo(1) - min_.at(1)) / dlt_.at(1));
    Double_t wyloc = std::move((wcoo(1) - min_.at(1)) / dlt_.at(1));
    Long64_t vyi = static_cast<Long64_t>(std::floor(vyloc));
    Long64_t wyi = static_cast<Long64_t>(std::floor(wyloc));
    if ((vyi < 0 && wyi < 0) || (vyi >= n_.at(1) && wyi >= n_.at(1))) return MatFld();
    
    Double_t vxloc = std::move((vcoo(0) - min_.at(0)) / dlt_.at(0));
    Double_t wxloc = std::move((wcoo(0) - min_.at(0)) / dlt_.at(0));
    Long64_t vxi = static_cast<Long64_t>(std::floor(vxloc));
    Long64_t wxi = static_cast<Long64_t>(std::floor(wxloc));
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
    Long64_t   nstp   = static_cast<Long64_t>(std::floor(vwlen / (is_std ? STD_STEP_LEN_ : FST_STEP_LEN_))) + 2;
    SVecD<3>&& unit   = (vwvec / static_cast<Double_t>(nstp));
    Double_t   ulen   = LA::Mag(unit);

    SVecD<3> itloc((vxloc + MGMath::HALF * unit(0)), (vyloc + MGMath::HALF * unit(1)), (vzloc + MGMath::HALF * unit(2)));
    Long64_t itsat = 0;
    Long64_t itend = nstp;

    //==== faster method (tuned by axis-Z)
    Short_t uz_sign = MGNumc::Compare(unit(2));
    if (uz_sign != 0) {
        Double_t sat = ((uz_sign == 1) ? 0. : n_.at(2));
        Double_t end = ((uz_sign == 1) ? n_.at(2) : 0.);
        Long64_t satID = static_cast<Long64_t>(std::floor((sat - itloc(2)) / unit(2)));
        Long64_t endID = static_cast<Long64_t>(std::floor((end - itloc(2)) / unit(2))) + 1;
        if (MGNumc::Valid(satID) && MGNumc::Valid(endID)) {
            if (satID > 0 && satID < nstp) itsat = satID;
            if (endID > 0 && endID < nstp) itend = endID;
            itloc += (itsat * unit);
        }
    }
    //====

    Long64_t suml = 0;
    Long64_t sumll = 0;
    Long64_t itcnt = 0;
    for (Long64_t it = itsat; it < itend; ++it, itloc += unit) {
        Long64_t zi = static_cast<Long64_t>(std::floor(itloc(2)));
        if (zi < 0 || zi >= n_.at(2)) continue;
        Long64_t yi = static_cast<Long64_t>(std::floor(itloc(1)));
        if (yi < 0 || yi >= n_.at(1)) continue;
        Long64_t xi = static_cast<Long64_t>(std::floor(itloc(0)));
        if (xi < 0 || xi >= n_.at(0)) continue;

        Long64_t idx = (xi * fact_.at(0) + yi * fact_.at(1) + zi);
        if (mat_ptr_[idx]) {
            suml += it;
            sumll += it*it;
            itcnt++;
        }
    }

    if (itcnt != 0) {
        Double_t efft = static_cast<Double_t>(itcnt) / static_cast<Double_t>(nstp);
        Double_t efft_len = efft * real_len;
        Double_t sftloc = (MGMath::HALF / static_cast<Double_t>(nstp));
        Double_t sftsqr = (sftloc * sftloc);
        Double_t avgloc = (static_cast<Double_t>(suml)  / static_cast<Double_t>(itcnt) / static_cast<Double_t>(nstp));
        Double_t avgsqr = (static_cast<Double_t>(sumll) / static_cast<Double_t>(itcnt) / static_cast<Double_t>(nstp * nstp));
        Double_t loc    = (avgloc + sftloc);
        Double_t locsqr = (avgsqr + MGMath::TWO * sftloc * avgloc + sftsqr);

        return MatFld(true, elm_, den_, inv_rad_len_, real_len, efft_len, efft, loc, locsqr);
    }
    else return MatFld(real_len);
}


MatFld MatMgnt::Get(const SVecD<3>& coo) {
    if (!Load()) return MatFld();

    Bool_t mat = false;
    Double_t real_len = 0.0;
    std::array<Bool_t,   MatProperty::NUM_ELM> elm; elm.fill(false);
    std::array<Double_t, MatProperty::NUM_ELM> den; den.fill(0.);
    Double_t                                   inv_rad_len = 0.0;

    for (auto&& reader : *reader_) {
        if (!reader->is_in_box(coo)) continue;
        MatFld&& mfld = reader->get(coo);
        if (!mfld()) continue;

        for (Int_t it = 0; it < MatProperty::NUM_ELM; ++it) {
            if (!mfld.elm().at(it)) continue;
            elm.at(it)  = true;
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
    Double_t                                   locsqr = 0.0;

    for (auto&& reader : *reader_) {
        if (!reader->is_cross(vcoo, wcoo)) continue;
        MatFld&& mfld = reader->get(vcoo, wcoo, is_std);
        if (!mfld()) continue;
        
        for (Int_t it = 0; it < MatProperty::NUM_ELM; ++it) {
            if (!mfld.elm().at(it)) continue;
            elm.at(it)  = true;
            den.at(it) += (mfld.efft_len() * mfld.den().at(it));
        }
        Double_t num_rad_len = mfld.num_rad_len();
        loc         += mfld.loc() * num_rad_len;
        locsqr      += mfld.locsqr() * num_rad_len;
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
        locsqr      = (locsqr / inv_rad_len);
        inv_rad_len = (inv_rad_len / efft_len);
        efft = (efft_len / real_len);

        return MatFld(mat, elm, den, inv_rad_len, real_len, efft_len, efft, loc, locsqr);
    }
    else return MatFld(real_len);
}


MatFld MatMgnt::Get(Double_t stp_len, const PhySt& part, Bool_t is_std) {
    const SVecD<3>&  vcoo = part.c();
    SVecD<3>&&       wcoo = part.c() + stp_len * part.u();

    return Get(vcoo, wcoo, is_std);
}


Double_t MatPhy::GetNumRadLen(const Double_t stp_len, const PhySt& part, Bool_t is_std) {
    const SVecD<3>&  vcoo = part.c();
    SVecD<3>&&       wcoo = part.c() + stp_len * part.u();

    MatFld&& mfld = MatMgnt::Get(vcoo, wcoo, is_std);
   
    return mfld.num_rad_len();
}


MatPhyFld MatPhy::Get(const Double_t stp_len, PhySt& part, Bool_t is_std) {
    Double_t len = std::fabs(stp_len);
    if (!part.field()) return MatPhyFld(len);
    if (part.info().is_chrgless() || part.info().is_massless()) return MatPhyFld(len);
    if (MGNumc::EqualToZero(stp_len)) return MatPhyFld(len);
    if (MGNumc::EqualToZero(part.mom())) return MatPhyFld(len);

    const SVecD<3>&  vcoo = part.c();
    SVecD<3>&&       wcoo = part.c() + stp_len * part.u();

    MatFld&& mfld = MatMgnt::Get(vcoo, wcoo, is_std);
    
    if (!mfld()) return MatPhyFld(mfld.real_len());

    Double_t mult_scat_sgm = GetMultipleScattering(mfld, part);
    std::tuple<Double_t, Double_t, Double_t, Double_t>&& ion_eloss = GetIonizationEnergyLoss(mfld, part);
    Double_t eloss_brm_men = GetBremsstrahlungEnergyLoss(mfld, part);

    return MatPhyFld(mfld(), mfld.real_len(), mfld.efft(), mfld.loc(), mfld.locsqr(), mfld.inv_rad_len(), mfld.num_rad_len(), mult_scat_sgm, std::get<0>(ion_eloss), std::get<1>(ion_eloss), std::get<2>(ion_eloss), std::get<3>(ion_eloss), eloss_brm_men);
}


MatPhyFld MatPhy::Get(const MatFld& mfld, PhySt& part) {
    Double_t len = mfld.real_len();
    if (!mfld() || !part.field()) return MatPhyFld(len);
    if (MGNumc::EqualToZero(mfld.efft_len())) return MatPhyFld(len);
    if (part.info().is_chrgless() || part.info().is_massless()) return MatPhyFld(len);
    if (MGNumc::EqualToZero(part.mom())) return MatPhyFld(len);
    
    Double_t mult_scat_sgm = GetMultipleScattering(mfld, part);
    std::tuple<Double_t, Double_t, Double_t, Double_t>&& ion_eloss = GetIonizationEnergyLoss(mfld, part);
    Double_t eloss_brm_men = GetBremsstrahlungEnergyLoss(mfld, part);
    
    return MatPhyFld(mfld(), mfld.real_len(), mfld.efft(), mfld.loc(), mfld.locsqr(), mfld.inv_rad_len(), mfld.num_rad_len(), mult_scat_sgm, std::get<0>(ion_eloss), std::get<1>(ion_eloss), std::get<2>(ion_eloss), std::get<3>(ion_eloss), eloss_brm_men);
}
        

std::array<Double_t, MatProperty::NUM_ELM> MatPhy::GetDensityEffectCorrection(const MatFld& mfld, PhySt& part) {
    Double_t gmbta     = ((MGNumc::Compare(part.bta(), LMT_BTA) > 0) ? part.gmbta() : LMT_GMBTA);
    Double_t log_gmbta = std::log10(gmbta);

    std::array<Double_t, MatProperty::NUM_ELM> dlt; dlt.fill(0.);
    for (Int_t it = 0; it < MatProperty::NUM_ELM; ++it) {
       if (!mfld.elm().at(it)) continue;
       if (log_gmbta < MatProperty::DEN_EFF_CORR_X0.at(it)) continue; // if nonconductors
       Double_t dif_log_gmbta = MatProperty::DEN_EFF_CORR_X1.at(it) - log_gmbta;
       dlt.at(it) = MGMath::TWO * MGMath::LOG_TEN * log_gmbta - MatProperty::DEN_EFF_CORR_C.at(it);
       if (MGNumc::Compare(dif_log_gmbta) <= 0) continue;
       dlt.at(it) += MatProperty::DEN_EFF_CORR_A.at(it) * std::pow(dif_log_gmbta, MatProperty::DEN_EFF_CORR_K.at(it));
       if (MGNumc::Compare(dlt.at(it)) <= 0) dlt.at(it) = MGMath::ZERO;
    }

    return dlt;
}


Double_t MatPhy::GetMultipleScattering(const MatFld& mfld, PhySt& part) {
    if (!part.arg().mscat()) return MGMath::ZERO;
    
    Double_t num_rad_len = mfld.num_rad_len();
    if (MGNumc::Compare(num_rad_len, LMTL_NUM_RAD_LEN) < 0) num_rad_len = LMTL_NUM_RAD_LEN;
    if (MGNumc::Compare(num_rad_len, LMTU_NUM_RAD_LEN) > 0) num_rad_len = LMTU_NUM_RAD_LEN;
    Bool_t is_over_lmt = (MGNumc::Compare(part.bta(), LMT_BTA) > 0);
    Double_t bta = ((is_over_lmt) ? part.bta() : LMT_BTA);
    Double_t eta = ((is_over_lmt) ? part.eta_abs() : LMT_INV_GMBTA);
    Double_t eta_part = (eta / bta);
    Double_t log_nrl = std::log(num_rad_len);
    
    // Highland-Lynch-Dahl formula
    Double_t mscat_sgm = RYDBERG_CONST * part.info().chrg_to_mass() * eta_part * std::sqrt(num_rad_len) * (MGMath::ONE + NRL_CORR_FACT * log_nrl);
    
    // Modified Highland-Lynch-Dahl formula
    //Double_t mscat_sgm = RYDBERG_CONST * part.info().chrg_to_mass() * eta_part * std::sqrt(num_rad_len * (MGMath::ONE + NRL_CORR_FACT1 * log_nrl + NRL_CORR_FACT2 * log_nrl * log_nrl));
    
    // Tune by Hsin-Yi Chou
    //const Double_t pars[3] = { 8.47474822486500079e-01, 6.44047741589626535e-03, -2.74914 }; 
    //const Double_t tune_sgm = pars[0] * (MGMath::ONE + pars[1] * TMath::Power(bta, pars[2]));
    //mscat_sgm *= tune_sgm;
    
    //----------------- AMS
    //Double_t mscat_sgm = RYDBERG_CONST * part.info().chrg_to_mass() * eta_part * std::sqrt(num_rad_len);
    //----------------- AMS
    
    if (!MGNumc::Valid(mscat_sgm) || MGNumc::Compare(mscat_sgm) <= 0) mscat_sgm = MGMath::ZERO;

    return mscat_sgm;
}


std::tuple<Double_t, Double_t, Double_t, Double_t> MatPhy::GetIonizationEnergyLoss(const MatFld& mfld, PhySt& part) {
    if (!part.arg().eloss()) return std::make_tuple(MGMath::ZERO, MGMath::ZERO, MGMath::ZERO, MGMath::ZERO);
    
    Bool_t is_over_lmt   = (MGNumc::Compare(part.bta(), LMT_BTA) > 0);
    Double_t gmbta       = ((is_over_lmt) ? part.gmbta() : LMT_GMBTA);
    Double_t sqr_gmbta   = ((is_over_lmt) ? (gmbta * gmbta) : LMT_SQR_GMBTA);
    Double_t bta         = ((is_over_lmt) ? part.bta() : LMT_BTA);
    Double_t sqr_bta     = ((is_over_lmt) ? (part.bta() * part.bta()) : LMT_SQR_BTA);
    Double_t gm          = ((is_over_lmt) ? part.gm() : LMT_GM);
    Double_t sqr_chrg    = part.chrg() * part.chrg();
    Double_t mass_in_GeV = part.mass();
    Double_t mass_in_MeV = part.mass() * GEV_TO_MEV;
    Double_t rel_mass    = (MASS_EL_IN_GEV / mass_in_GeV);
        
    // Density Effect Correction
    std::array<Double_t, MatProperty::NUM_ELM>&& dlt = GetDensityEffectCorrection(mfld, part);
    
    // Calculate Weighted Average
    Double_t avg_dlt = MGMath::ZERO;
    Double_t m_exeng = MGMath::ZERO;
    Double_t ttl_den = MGMath::ZERO;
    for (Int_t it = 0; it < MatProperty::NUM_ELM; ++it) {
       if (!mfld.elm().at(it)) continue;
       Double_t den = MatProperty::CHRG.at(it) * mfld.den().at(it);
       m_exeng += den * MatProperty::NEG_LN_MEAN_EXENG.at(it);
       avg_dlt += den * dlt.at(it);
       ttl_den += den;
    }
    avg_dlt /= ttl_den;
    m_exeng /= ttl_den;

    if (!MGNumc::Valid(ttl_den) || MGNumc::Compare(ttl_den) <= 0) ttl_den = MGMath::ZERO;
    if (!MGNumc::Valid(avg_dlt) || MGNumc::Compare(avg_dlt) <= 0) avg_dlt = MGMath::ZERO;
    if (!MGNumc::Valid(m_exeng) || MGNumc::Compare(m_exeng) <= 0) m_exeng = MGMath::ZERO;
    Double_t exeng = std::exp(-m_exeng); // [MeV]
 
    // Calculate Matterial Quality
    Double_t electron_cloud_abundance = (ttl_den * mfld.efft_len());
    Double_t Bethe_Bloch_fact = (MGMath::HALF * BETHE_BLOCH_K * electron_cloud_abundance * sqr_chrg / sqr_bta);
  
    // Calculate Eta Trans
    Double_t eta_trans = (std::sqrt(sqr_gmbta + MGMath::ONE) / sqr_gmbta);
    
    Double_t Bethe_Bloch_eta_trans = (Bethe_Bloch_fact * eta_trans / mass_in_MeV);
    
    // Calculate Sigma
    Double_t eloss_ion_sgm = Bethe_Bloch_eta_trans;
    
    // Calculate Peak
    Double_t trans_eng     = MGMath::TWO * MASS_EL_IN_MEV * sqr_gmbta;
    Double_t max_keng_part = std::log(trans_eng / exeng);
    Double_t ion_keng_part = std::log(Bethe_Bloch_fact / exeng);
    if (!MGNumc::Valid(max_keng_part) || MGNumc::Compare(max_keng_part) <= 0) max_keng_part = MGMath::ZERO;
    if (!MGNumc::Valid(ion_keng_part) || MGNumc::Compare(ion_keng_part) <= 0) ion_keng_part = MGMath::ZERO;
    Double_t eloss_ion_mpv = Bethe_Bloch_eta_trans * (max_keng_part + ion_keng_part + LANDAU_ELOSS_CORR - sqr_bta - avg_dlt);
    
    // Calculate Mean
    Double_t mass_rat = (MASS_EL_IN_GEV / mass_in_GeV);
    Double_t max_trans_eng_part = std::log((trans_eng / exeng) * (trans_eng / exeng) / (MGMath::ONE + MGMath::TWO * gm * mass_rat + mass_rat * mass_rat));
    Double_t eloss_ion_men = Bethe_Bloch_eta_trans * (max_trans_eng_part - MGMath::TWO * sqr_bta - avg_dlt);
  
    // Calculate Kappa
    const Double_t eloss_ion_kpa0 = MGMath::SQRT_TWO;
    Double_t eloss_ion_kpa = eloss_ion_kpa0 / (MGMath::HALF * (sqr_bta - MGMath::ONE) + MGMath::ONE/sqr_bta);

    // Tune by Hsin-Yi Chou
    //const Double_t tune_kpa = 1.05379361892339474e+00;
    //eloss_ion_kpa *= tune_kpa;

    //const Double_t tune_mpv_pars[2] = { 0.903547, 0.0045 };
    //const Double_t tune_mpv = tune_mpv_pars[0] / (MGMath::ONE + tune_mpv_pars[1] * (std::pow(sqr_gmbta, -MGMath::SQRT_TWO) + MGMath::TWO * (sqr_bta - MGMath::ONE)));
    //eloss_ion_mpv *= tune_mpv;

    //----------------- AMS
    eloss_ion_kpa  = MGMath::ONE;
    eloss_ion_mpv  = eloss_ion_men;
    eloss_ion_sgm  *= MGMath::TWO;
    //----------------- AMS

    if (!MGNumc::Valid(eloss_ion_kpa) || MGNumc::Compare(eloss_ion_kpa) <= 0) eloss_ion_kpa = MGMath::ZERO;
    if (!MGNumc::Valid(eloss_ion_mpv) || MGNumc::Compare(eloss_ion_mpv) <= 0) eloss_ion_mpv = MGMath::ZERO;
    if (!MGNumc::Valid(eloss_ion_sgm) || MGNumc::Compare(eloss_ion_sgm) <= 0) eloss_ion_sgm = MGMath::ZERO;
    if (!MGNumc::Valid(eloss_ion_men) || MGNumc::Compare(eloss_ion_men) <= 0) eloss_ion_men = MGMath::ZERO;

    return std::make_tuple(eloss_ion_kpa, eloss_ion_mpv, eloss_ion_sgm, eloss_ion_men);
}


Double_t MatPhy::GetBremsstrahlungEnergyLoss(const MatFld& mfld, PhySt& part) {
    if (!part.arg().eloss()) return MGMath::ZERO;
    // testcode
    return 0.0; // TODO: Include Bremsstrahlung

    Bool_t   is_over_lmt  = (MGNumc::Compare(part.bta(), LMT_BTA) > 0);
    Double_t sqr_gmbta    = ((is_over_lmt) ? (part.gmbta() * part.gmbta()) : LMT_SQR_GMBTA);
    Double_t sqr_chrg_rat = (part.chrg() * part.chrg());
    Double_t sqr_mass_rat = (MASS_EL_IN_GEV * MASS_EL_IN_GEV / part.mass() / part.mass());
    Double_t eng          = std::sqrt(sqr_gmbta + MGMath::ONE);
    Double_t ke_part      = (eng - MGMath::ONE);
    Double_t eta_trans    = (eng / sqr_gmbta);
    Double_t eloss_brm_men = sqr_chrg_rat * sqr_mass_rat * (ke_part * eta_trans) * (mfld.num_rad_len() / MGMath::LOG_TWO);
    
    if (!MGNumc::Valid(eloss_brm_men) || MGNumc::Compare(eloss_brm_men) <= 0) eloss_brm_men = MGMath::ZERO;

    return eloss_brm_men;
}


} // namespace TrackSys


#endif // __TRACKLibs_MatEnv_C__
