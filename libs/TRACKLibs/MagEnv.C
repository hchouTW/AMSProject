#ifndef __TRACKLibs_MagEnv_C__
#define __TRACKLibs_MagEnv_C__


namespace TrackSys {


void MagFld::print() const {
    std::string printStr = STR_FMT("MAG (%8.5f %8.5f %8.5f)\n", mag_(0), mag_(1), mag_(2));
    COUT(printStr);
}


MagGeoBoxCreator::MagGeoBoxCreator(Long64_t xn, Float_t xmin, Float_t xmax, Long64_t yn, Float_t ymin, Float_t ymax, Long64_t zn, Float_t zmin, Float_t zmax, const std::string& file_path) {
    clear();
    if (xn < 2 || yn < 2 || zn < 2) { MGSys::ShowError("MagGeoBoxCreator::MagGeoBoxCreator() : Size failure."); return; }
    if (MGNumc::Compare(xmin, xmax) >= 0 || MGNumc::Compare(ymin, ymax) >= 0 || MGNumc::Compare(zmin, zmax) >= 0) { MGSys::ShowError("MagGeoBoxCreator::MagGeoBoxCreator() : Range failure."); return; }
    if (file_path.size() == 0) { MGSys::ShowError("MagGeoBoxCreator::MagGeoBoxCreator() : File path not found."); return; }
    
    Int_t file_len = ((sizeof(Long64_t) + sizeof(Float_t) + sizeof(Float_t)) + (xn * yn * zn) * sizeof(Float_t)) * DIM_;
    Int_t file_des = open(file_path.c_str(), O_CREAT | O_RDWR, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
    if (file_des < 0) { MGSys::ShowError("MagGeoBoxCreator::MagGeoBoxCreator() : File not opened."); return; }

    off_t sret = lseek(file_des, file_len, SEEK_SET);
    write(file_des, "\0", 1);

    void* file_ptr = mmap(nullptr, file_len, PROT_READ | PROT_WRITE, MAP_SHARED, file_des, 0);
    if (file_ptr == reinterpret_cast<void*>(-1)) { MGSys::ShowError("MagGeoBoxCreator::MagGeoBoxCreator() : mmap() failure."); close(file_des); return; }
    
    is_open_   = true;
    file_path_ = file_path;
    file_des_  = file_des;
    file_len_  = file_len;
    file_ptr_  = file_ptr;
    max_len_   = (xn * yn * zn);
    geo_box_   = reinterpret_cast<MagGeoBox*>(file_ptr);

    geo_box_->n[0] = xn;
    geo_box_->n[1] = yn;
    geo_box_->n[2] = zn;
    geo_box_->min[0] = xmin;
    geo_box_->min[1] = ymin;
    geo_box_->min[2] = zmin;
    geo_box_->max[0] = xmax;
    geo_box_->max[1] = ymax;
    geo_box_->max[2] = zmax;
    std::fill_n(static_cast<Float_t*>(&(geo_box_->mag)), max_len_ * DIM_, 0.0);
}


void MagGeoBoxCreator::fill(Long64_t idx, Float_t bx, Float_t by, Float_t bz) {
    if (!is_open_) return;
    if (idx < 0 || idx >= max_len_) return;

    Float_t * mag = static_cast<Float_t*>(&(geo_box_->mag) + idx * DIM_);
    mag[0] = bx;
    mag[1] = by;
    mag[2] = bz;
}


void MagGeoBoxCreator::save_and_close() {
    if (!is_open_) { clear(); return; }
        
    munmap(file_ptr_, file_len_);
    close(file_des_);

    clear();
}
        

void MagGeoBoxCreator::save_and_close(Float_t bx, Float_t by, Float_t bz) {
    if (!is_open_) return;

    for (Long64_t idx = 0; idx < max_len_; ++idx) {
        Float_t * mag = static_cast<Float_t*>(&(geo_box_->mag) + idx * DIM_);
        mag[0] = bx;
        mag[1] = by;
        mag[2] = bz;
    }

    munmap(file_ptr_, file_len_);
    close(file_des_);

    clear();
}


Bool_t MagGeoBoxReader::load(const std::string& file_path) {
    if (is_load_) {
        if (file_path_ == file_path) return is_load_;
        else {
            MGSys::ShowWarning("MagGeoBoxReader::Load() : re-load different file.");
            clear();
        }
    }

    Int_t file_des = open(file_path.c_str(), O_RDONLY);
    Int_t file_len = lseek(file_des, 0, SEEK_END); 
    if (file_des < 0) {
        MGSys::ShowError(STR_FMT("MagGeoBoxReader::Load() : Magnetic field map not found (%s)", file_path.c_str()));
        is_load_ = false;
        return is_load_;
    }
    
    void* file_ptr = mmap(nullptr, file_len, PROT_READ, MAP_SHARED, file_des, 0);
    if (file_ptr == reinterpret_cast<void*>(-1)) {
        MGSys::ShowError("MagGeoBoxReader::Load() : mmap() failure.");
        close(file_des);
        is_load_ = false;
        return is_load_;
    }
    file_ptr_ = file_ptr;

    MagGeoBox* geo_box = reinterpret_cast<MagGeoBox*>(file_ptr_);
    for (Long64_t ig = 0; ig < DIM_; ++ig) {
        n_.at(ig)   = geo_box->n[ig]; 
        min_.at(ig) = geo_box->min[ig];
        max_.at(ig) = geo_box->max[ig];
        dlt_.at(ig) = (max_.at(ig) - min_.at(ig)) / static_cast<Float_t>(n_.at(ig)-1);
    }
    fact_.at(0) = n_.at(1) * n_.at(2);
    fact_.at(1) = n_.at(2);
    mag_ptr_ = &(geo_box->mag);
    
    is_load_ = true;
    file_path_ = file_path;
    COUT("MagGeoBoxReader::Load() : Open file (%s)\n", file_path.c_str());
    return is_load_;
}


MagFld MagGeoBoxReader::get(const SVecD<3>& coo) {
    if (!is_load_) return MagFld();
    
    // Get Local Coord
    Float_t xloc = (coo[0] - min_.at(0)) / dlt_.at(0);
    Float_t yloc = (coo[1] - min_.at(1)) / dlt_.at(1);
    Float_t zloc = (coo[2] - min_.at(2)) / dlt_.at(2);

    Long64_t xi = static_cast<Long64_t>(xloc);
    Long64_t yi = static_cast<Long64_t>(yloc);
    Long64_t zi = static_cast<Long64_t>(zloc);

    if (xi < 0 || xi >= (n_.at(0)-1) || 
        yi < 0 || yi >= (n_.at(1)-1) || 
        zi < 0 || zi >= (n_.at(2)-1)) return MagFld();

    // Do Trilinear Interpolation
    Long64_t idx = (xi * fact_.at(0) + yi * fact_.at(1) + zi);
    Float_t xu = xloc - static_cast<Float_t>(xi);
    Float_t yu = yloc - static_cast<Float_t>(yi);
    Float_t zu = zloc - static_cast<Float_t>(zi);
    Float_t xl = (MGMath::ONE - xu);
    Float_t yl = (MGMath::ONE - yu);
    Float_t zl = (MGMath::ONE - zu);

    Float_t * cell = (mag_ptr_ + (idx * DIM_));
    Float_t * clll = (cell                                          );
    Float_t * cllu = (cell + (                          + 1) * DIM_ );
    Float_t * clul = (cell + (            + fact_.at(1)    ) * DIM_ );
    Float_t * cluu = (cell + (            + fact_.at(1) + 1) * DIM_ );
    Float_t * cull = (cell + (fact_.at(0)                  ) * DIM_ );
    Float_t * culu = (cell + (fact_.at(0)               + 1) * DIM_ );
    Float_t * cuul = (cell + (fact_.at(0) + fact_.at(1)    ) * DIM_ );
    Float_t * cuuu = (cell + (fact_.at(0) + fact_.at(1) + 1) * DIM_ );

    Float_t cll[3] { (clll[0] * xl + cull[0] * xu),  (clll[1] * xl + cull[1] * xu), (clll[2] * xl + cull[2] * xu) };
    Float_t clu[3] { (cllu[0] * xl + culu[0] * xu),  (cllu[1] * xl + culu[1] * xu), (cllu[2] * xl + culu[2] * xu) };
    Float_t cul[3] { (clul[0] * xl + cuul[0] * xu),  (clul[1] * xl + cuul[1] * xu), (clul[2] * xl + cuul[2] * xu) };
    Float_t cuu[3] { (cluu[0] * xl + cuuu[0] * xu),  (cluu[1] * xl + cuuu[1] * xu), (cluu[2] * xl + cuuu[2] * xu) };
    
    Float_t cl[3] { (cll[0] * yl + cul[0] * yu),  (cll[1] * yl + cul[1] * yu), (cll[2] * yl + cul[2] * yu) };
    Float_t cu[3] { (clu[0] * yl + cuu[0] * yu),  (clu[1] * yl + cuu[1] * yu), (clu[2] * yl + cuu[2] * yu) };

    Float_t c[3] { (cl[0] * zl + cu[0] * zu), (cl[1] * zl + cu[1] * zu), (cl[2] * zl + cu[2] * zu) };

    return MagFld(c[0], c[1], c[2]);
}


#ifdef __HAS_AMS_OFFICE_LIBS__
Bool_t MagGeoBoxAms::Load() {
    if (is_load_) return is_load_;
    is_load_ = false;
    
    mag_field_ = MagField::GetPtr();
    if (mag_field_ == nullptr) return is_load_; 
    
    if (mag_field_->GetMap()) is_load_ = true;
    else {
        std::string filePath = STR_FMT("%s/v5.00/MagneticFieldMapPermanent_NEW_FULL.bin", MGSys::GetEnv("AMSDataDir").c_str());
        if ((mag_field_->Read(filePath.c_str())) < 0) {
            CERR("Magnetic Field map not found : %s\n", filePath.c_str());
            mag_field_ = nullptr;
        }
        else {
            COUT("MagGeoBoxAms::Load() Open file : %s\n", filePath.c_str());
            mag_field_->SetMagstat(1);
            mag_field_->SetScale(1);
            is_load_ = true;
        }
    }

    return is_load_;
}


MagFld MagGeoBoxAms::Get(const SVecD<3>& coo) {
    if (!Load()) return MagFld();
    Float_t incoo[3] = { static_cast<Float_t>(coo[0]), static_cast<Float_t>(coo[1]), static_cast<Float_t>(coo[2]) };
    Float_t outmag[3] = { 0, 0, 0 };
    mag_field_->GuFld(incoo, outmag);
    return MagFld(outmag);
}
        

void MagGeoBoxAms::Output(const std::string& file_path) {
    if (!Load()) return;
    const Long64_t n = 201;
    const Float_t  min = -200.0;
    const Float_t  max =  200.0;
    const Float_t  dlt =    2.0;
    MagGeoBoxCreator creator(
            n, min, max,
            n, min, max,
            n, min, max,
            file_path
        );
    if (!creator.is_open()) return;
    

    COUT("\n");
    COUT("MagGeoBoxAms::Output() : start.\n");
    COUT("MagGeoBoxAms::Output() : X (%lld %8.2f %8.2f)\n", n, min, max);
    COUT("MagGeoBoxAms::Output() : Y (%lld %8.2f %8.2f)\n", n, min, max);
    COUT("MagGeoBoxAms::Output() : Z (%lld %8.2f %8.2f)\n", n, min, max);

    Long64_t fact[2] { (n * n), (n) };
    for (Long64_t xi = 0; xi < n; ++xi) {
        for (Long64_t yi = 0; yi < n; ++yi) {
            for (Long64_t zi = 0; zi < n; ++zi) {
                Long64_t idx = xi * fact[0] + yi * fact[1] + zi;
                SVecD<3> coo((min + static_cast<Double_t>(xi) * dlt), (min + static_cast<Double_t>(yi) * dlt), (min + static_cast<Double_t>(zi) * dlt));
                MagFld&& mag = MagGeoBoxAms::Get(coo);
                creator.fill(idx, static_cast<Float_t>(mag.x()), static_cast<Float_t>(mag.y()), static_cast<Float_t>(mag.z()));
            }
        }
    }
    COUT("MagGeoBoxAms::Output() : save.\n");
    creator.save_and_close();
    
    COUT("MagGeoBoxAms::Output() : finish.\n");
    COUT("\n");
}


Double_t MagFuncAms::GetMagx(Double_t cooz) {
    Double_t norm0 = MGMath::INV_SQRT_TWO * cooz / PAR_SGM.at(0);
    Double_t norm1 = MGMath::INV_SQRT_TWO * cooz / PAR_SGM.at(1);
    Double_t magx0 = PAR_MAG.at(0) * (MGMath::INV_SQRT_TWO_PI / PAR_SGM.at(0)) * std::exp(MGMath::NEG * norm0 * norm0);
    Double_t magx1 = PAR_MAG.at(1) * (MGMath::INV_SQRT_TWO_PI / PAR_SGM.at(1)) * std::exp(MGMath::NEG * norm1 * norm1);
    return (magx0 + magx1);
}


Double_t MagFuncAms::GetMagxInt1st(Double_t cooz) {
    Double_t norm0 = MGMath::INV_SQRT_TWO * cooz / PAR_SGM.at(0);
    Double_t norm1 = MGMath::INV_SQRT_TWO * cooz / PAR_SGM.at(1);
    Double_t magxint0 = PAR_MAG.at(0) * MGMath::HALF * std::erf(norm0);
    Double_t magxint1 = PAR_MAG.at(1) * MGMath::HALF * std::erf(norm1);
    return (magxint0 + magxint1);
}


Double_t MagFuncAms::GetMagxInt2nd(Double_t cooz) {
    Double_t fact0 = MGMath::SQRT_TWO * PAR_SGM.at(0);
    Double_t fact1 = MGMath::SQRT_TWO * PAR_SGM.at(1);
    Double_t norm0 = cooz / fact0;
    Double_t norm1 = cooz / fact1;
    Double_t magxint0 = PAR_MAG.at(0) * (MGMath::HALF * fact0) * (norm0 * std::erf(norm0) + ((std::exp(MGMath::NEG * norm0 * norm0) - MGMath::ONE) * MGMath::SQRT_PI));
    Double_t magxint1 = PAR_MAG.at(1) * (MGMath::HALF * fact1) * (norm1 * std::erf(norm1) + ((std::exp(MGMath::NEG * norm1 * norm1) - MGMath::ONE) * MGMath::SQRT_PI));
    return (magxint0 + magxint1);
}
#endif // __HAS_AMS_OFFICE_LIBS__


Bool_t MagMgnt::Load() {
    if (is_load_) return is_load_;

#ifdef __HAS_AMS_OFFICE_LIBS__
    is_load_ = MagGeoBoxAms::Load();
#else
    if (!geo_box_reader_.exist()) {
        std::string file_path = "/data3/hchou/AMSData/MagDB/MagTest.bin";
        geo_box_reader_.load(file_path);
    }
    is_load_ = (geo_box_reader_.exist());
#endif // __HAS_AMS_OFFICE_LIBS__

    return is_load_;
}


MagFld MagMgnt::Get(const SVecD<3>& coo) {
    if (!Load()) return MagFld();

#ifdef __HAS_AMS_OFFICE_LIBS__
    return MagGeoBoxAms::Get(coo);
#endif // __HAS_AMS_OFFICE_LIBS__

#ifdef __HAS_TESTPROP__
    return MagFld();
#endif // __HAS_TESTPROP__

#ifdef __HAS_TESTFIT__
    return geo_box_reader_.get(coo);
#endif // __HAS_TESTFIT__
}


} // namespace TrackSys


#endif // __TRACKLibs_MagEnv_C__
