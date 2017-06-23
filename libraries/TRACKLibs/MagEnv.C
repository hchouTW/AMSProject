#ifndef __TRACKLibs_MagEnv_C__
#define __TRACKLibs_MagEnv_C__

namespace TrackSys {

MagGeoBoxCreator::MagGeoBoxCreator(Long64_t xn, Double_t xmin, Double_t xmax, Long64_t yn, Double_t ymin, Double_t ymax, Long64_t zn, Double_t zmin, Double_t zmax, const std::string& file_path) {
    clear();
    if (xn < 2 || yn < 2 || zn < 2) { MGSys::ShowError("MagGeoBoxCreator::MagGeoBoxCreator() : Size failure."); return; }
    if (MGNumc::Compare(xmin, xmax) >= 0 || MGNumc::Compare(ymin, ymax) >= 0 || MGNumc::Compare(zmin, zmax) >= 0) { MGSys::ShowError("MagGeoBoxCreator::MagGeoBoxCreator() : Range failure."); return; }
    if (file_path.size() == 0) { MGSys::ShowError("MagGeoBoxCreator::MagGeoBoxCreator() : File path not found."); return; }
    
    Int_t file_len = ((sizeof(Long64_t) + sizeof(Double_t)) + (xn * yn * zn) * sizeof(Double_t)) * DIM_;
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
    cur_len_   = 0;
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
}


void MagGeoBoxCreator::fill(Double_t mx, Double_t my, Double_t mz) {
    if (!is_open_) return;
    if (cur_len_ >= max_len_) { MGSys::ShowError("MagGeoBoxCreator::fill() : Out range."); return; }

    Double_t * mag = static_cast<Double_t*>(&(geo_box_->mag) + cur_len_ * DIM_);
    mag[0] = mx;
    mag[1] = my;
    mag[2] = mz;
    cur_len_++;
}


void MagGeoBoxCreator::save_and_close() {
    if (!is_open_) { clear(); return; }
    if (cur_len_ < max_len_) {
        MGSys::ShowWarning(STR_FMT("MagGeoBoxCreator::close() : fill not finished. (%lld/%lld)", cur_len_, max_len_));
        for (Long64_t it = cur_len_; it < max_len_; ++it) fill();
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
    for (Int_t ig = 0; ig < DIM_; ++ig) {
        n_.at(ig)   = geo_box->n[ig]; 
        min_.at(ig) = geo_box->min[ig];
        max_.at(ig) = geo_box->max[ig];
        dlt_.at(ig) = (max_.at(ig) - min_.at(ig)) / static_cast<Double_t>(n_.at(ig)-1);
    }
    fact_.at(0) = n_.at(1) * n_.at(2);
    fact_.at(1) = n_.at(2);
    mag_ptr_ = &(geo_box->mag);
    
    is_load_ = true;
    file_path_ = file_path;
    COUT("MagGeoBoxReader::Load() : Open file (%s)\n", file_path.c_str());
    return is_load_;
}


MagGeoBoxReader::Index MagGeoBoxReader::get_index(const SVecD<3>& coo) {
    Double_t xloc = (coo[0] - min_.at(0)) / dlt_.at(0);
    Double_t yloc = (coo[1] - min_.at(1)) / dlt_.at(1);
    Double_t zloc = (coo[2] - min_.at(2)) / dlt_.at(2);
    Bool_t xin = (MGNumc::Compare(xloc, MGMath::ZERO) > 0 && MGNumc::Compare(xloc, static_cast<Double_t>(n_.at(0)-1)) < 0);
    Bool_t yin = (MGNumc::Compare(yloc, MGMath::ZERO) > 0 && MGNumc::Compare(yloc, static_cast<Double_t>(n_.at(1)-1)) < 0);
    Bool_t zin = (MGNumc::Compare(zloc, MGMath::ZERO) > 0 && MGNumc::Compare(zloc, static_cast<Double_t>(n_.at(2)-1)) < 0);
    Bool_t inbox = (xin && yin && zin);
    if (!inbox) return std::make_tuple(-1, 0.0, 0.0, 0.0);;
    Long64_t xi = static_cast<Long64_t>(xloc);
    Long64_t yi = static_cast<Long64_t>(yloc);
    Long64_t zi = static_cast<Long64_t>(zloc);
    Long64_t index = (xi * fact_.at(0) + yi * fact_.at(1) + zi);
    Double_t xc = xloc - static_cast<Double_t>(xi);
    Double_t yc = yloc - static_cast<Double_t>(yi);
    Double_t zc = zloc - static_cast<Double_t>(zi);
    return std::make_tuple(index, xc, yc, zc);
}


SVecD<3> MagGeoBoxReader::do_trilinear_interpolation(const MagGeoBoxReader::Index& index) {
    if (std::get<0>(index) < 0) return SVecD<3>();
    Long64_t idx = std::get<0>(index);
    Double_t xl = (MGMath::ONE - std::get<1>(index));
    Double_t yl = (MGMath::ONE - std::get<2>(index));
    Double_t zl = (MGMath::ONE - std::get<3>(index));
    Double_t xu = std::get<1>(index);
    Double_t yu = std::get<2>(index);
    Double_t zu = std::get<3>(index);

    Double_t * clll = (mag_ptr_ + (idx                                ) * DIM_ );
    Double_t * cllu = (mag_ptr_ + (idx                             + 1) * DIM_ );
    Double_t * clul = (mag_ptr_ + (idx               + fact_.at(1)    ) * DIM_ );
    Double_t * cluu = (mag_ptr_ + (idx               + fact_.at(1) + 1) * DIM_ );
    Double_t * cull = (mag_ptr_ + (idx + fact_.at(0)                  ) * DIM_ );
    Double_t * culu = (mag_ptr_ + (idx + fact_.at(0)               + 1) * DIM_ );
    Double_t * cuul = (mag_ptr_ + (idx + fact_.at(0) + fact_.at(1)    ) * DIM_ );
    Double_t * cuuu = (mag_ptr_ + (idx + fact_.at(0) + fact_.at(1) + 1) * DIM_ );

    Double_t cll[3] = { (clll[0] * xl + cull[0] * xu),  (clll[1] * xl + cull[1] * xu), (clll[2] * xl + cull[2] * xu) };
    Double_t clu[3] = { (cllu[0] * xl + culu[0] * xu),  (cllu[1] * xl + culu[1] * xu), (cllu[2] * xl + culu[2] * xu) };
    Double_t cul[3] = { (clul[0] * xl + cuul[0] * xu),  (clul[1] * xl + cuul[1] * xu), (clul[2] * xl + cuul[2] * xu) };
    Double_t cuu[3] = { (cluu[0] * xl + cuuu[0] * xu),  (cluu[1] * xl + cuuu[1] * xu), (cluu[2] * xl + cuuu[2] * xu) };
    
    Double_t cl[3] = { (cll[0] * yl + cul[0] * yu),  (cll[1] * yl + cul[1] * yu), (cll[2] * yl + cul[2] * yu) };
    Double_t cu[3] = { (clu[0] * yl + cuu[0] * yu),  (clu[1] * yl + cuu[1] * yu), (clu[2] * yl + cuu[2] * yu) };

    Double_t c[3] = { (cl[0] * zl + cu[0] * zu),  (cl[1] * zl + cu[1] * zu), (cl[2] * zl + cu[2] * zu) };

    return SVecD<3>(c[0], c[1], c[2]);
}


#ifdef __HAS_AMS_OFFICE_LIBS__
Bool_t MagGeoBoxAms::Load() {
    if (is_load_) return is_load_;
    mag_field_ = MagField::GetPtr();
    if (!mag_field_->GetMap()) {
        bool iswork = true;
        static int magerr = 0;
        if (!mag_field_->GetMap() && !magerr) {
            std::string filePath = STR_FMT("%s/v5.00/MagneticFieldMapPermanent_NEW_FULL.bin", MGSys::GetEnv("AMSDataDir").c_str());
            if ((mag_field_->Read(filePath.c_str())) < 0) {
                CERR("Magnetic Field map not found : %s\n", filePath.c_str());
                magerr = -1;
                iswork = false;
            }
            else {
                COUT("MagGeoBoxAms::Load() Open file : %s\n", filePath.c_str());
            }
            mag_field_->SetMagstat(1);
            mag_field_->SetScale(1);
        }
        if (!iswork) mag_field_ = 0;
        else         is_load_ = true;
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
    const Double_t min = -200.0;
    const Double_t max =  200.0;
    const Double_t dlt =    2.0;
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

    Long64_t len = n * n * n;
    for (Long64_t xi = 0; xi < n; ++xi) {
        for (Long64_t yi = 0; yi < n; ++yi) {
            for (Long64_t zi = 0; zi < n; ++zi) {
                Long64_t index = (xi * n * n + yi * n + zi);
                SVecD<3> coo(
                        (min + static_cast<Double_t>(xi) * dlt),
                        (min + static_cast<Double_t>(yi) * dlt),
                        (min + static_cast<Double_t>(zi) * dlt)
                    );
                MagFld&& mag = MagGeoBoxAms::Get(coo);
                creator.fill(mag.x(), mag.y(), mag.z());
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
    Double_t magxint0 = PAR_MAG.at(0) * MGMath::ONE_TO_TWO * std::erf(norm0);
    Double_t magxint1 = PAR_MAG.at(1) * MGMath::ONE_TO_TWO * std::erf(norm1);
    return (magxint0 + magxint1);
}


Double_t MagFuncAms::GetMagxInt2nd(Double_t cooz) {
    Double_t fact0 = MGMath::SQRT_TWO * PAR_SGM.at(0);
    Double_t fact1 = MGMath::SQRT_TWO * PAR_SGM.at(1);
    Double_t norm0 = cooz / fact0;
    Double_t norm1 = cooz / fact1;
    Double_t magxint0 = PAR_MAG.at(0) * (MGMath::ONE_TO_TWO * fact0) * (norm0 * std::erf(norm0) + ((std::exp(MGMath::NEG * norm0 * norm0) - MGMath::ONE) * MGMath::SQRT_PI));
    Double_t magxint1 = PAR_MAG.at(1) * (MGMath::ONE_TO_TWO * fact1) * (norm1 * std::erf(norm1) + ((std::exp(MGMath::NEG * norm1 * norm1) - MGMath::ONE) * MGMath::SQRT_PI));
    return (magxint0 + magxint1);
}
#endif // __HAS_AMS_OFFICE_LIBS__


Bool_t MagMgnt::Load() {
    if (is_load_) return is_load_;
#ifdef __HAS_AMS_OFFICE_LIBS__
    Bool_t is_load_ams = MagGeoBoxAms::Load();
#endif

    if (!geo_box_reader_.exist()) {
        std::string file_path = "/afs/cern.ch/work/h/hchou/public/DATABASE/detector/MagGeoBox_AMS02.bin";
        geo_box_reader_.load(file_path);
    }

#ifdef __HAS_AMS_OFFICE_LIBS__
    is_load_ = (is_load_ams && geo_box_reader_.exist());
#elif
    is_load_ = (geo_box_reader_.exist());
#endif

    return is_load_;
}


MagFld MagMgnt::Get(const SVecD<3>& coo, MagType type) {
    if (!Load()) return MagFld();

    switch (type) {
#ifdef __HAS_AMS_OFFICE_LIBS__
        case MagType::kGeoBoxAms : return MagGeoBoxAms::Get(coo); break;
        case MagType::kFuncAms : return MagFuncAms::Get(coo); break;
#elif
        case MagType::kGeoBox : return geo_box_reader_.get(coo); break;
#endif // __HAS_AMS_OFFICE_LIBS__
        default : break;
    }
    return MagFld(); 
}


} // namespace TrackSys


#endif // __TRACKLibs_MagEnv_C__
