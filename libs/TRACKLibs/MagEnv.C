#ifndef __TRACKLibs_MagEnv_C__
#define __TRACKLibs_MagEnv_C__


namespace TrackSys {


MagGeoBoxCreator::MagGeoBoxCreator(const Long64_t n[3], const Float_t min[3], const Float_t max[3], const std::string& fpath) {
    clear();
    if (n[0] < 2 || n[1] < 2 || n[2] < 2) { 
        Sys::ShowWarning("MagGeoBoxCreator::MagGeoBoxCreator() : Size failure."); return; 
    }
    if (Numc::Compare(min[0], max[0]) >= 0 || Numc::Compare(min[1], max[1]) >= 0 || Numc::Compare(min[2], max[2]) >= 0) { 
        Sys::ShowWarning("MagGeoBoxCreator::MagGeoBoxCreator() : Range failure."); return; 
    }
    if (fpath.size() == 0) { 
        Sys::ShowWarning("MagGeoBoxCreator::MagGeoBoxCreator() : File path not found."); return; 
    }
    
    Long64_t flen = ((sizeof(Long64_t) + sizeof(Float_t) + sizeof(Float_t)) + (n[0] * n[1] * n[2]) * sizeof(Float_t)) * DIM_;
    Long64_t fdes = open(fpath.c_str(), O_CREAT | O_RDWR, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
    if (fdes < 0) { 
        Sys::ShowWarning("MagGeoBoxCreator::MagGeoBoxCreator() : File not opened."); return; 
    }
    lseek(fdes, flen, SEEK_SET); 
    
    void* fptr = mmap(nullptr, flen, PROT_READ | PROT_WRITE, MAP_SHARED, fdes, 0);
    if (fptr == reinterpret_cast<void*>(-1)) { 
        Sys::ShowWarning("MagGeoBoxCreator::MagGeoBoxCreator() : mmap() failure."); close(fdes); return; 
    }
    
    is_open_ = true;
    fdes_    = fdes;
    flen_    = flen;
    fptr_    = fptr;
    max_len_ = (n[0] * n[1] * n[2]);
    geo_box_ = reinterpret_cast<MagGeoBox*>(fptr_);

    geo_box_->n[0] = n[0];
    geo_box_->n[1] = n[1];
    geo_box_->n[2] = n[2];
    geo_box_->min[0] = min[0];
    geo_box_->min[1] = min[1];
    geo_box_->min[2] = min[2];
    geo_box_->max[0] = max[0];
    geo_box_->max[1] = max[1];
    geo_box_->max[2] = max[2];

    std::fill_n(static_cast<Float_t*>(&(geo_box_->mag)), max_len_ * DIM_, Numc::ZERO<Float_t>);
}


void MagGeoBoxCreator::fill(Long64_t idx, Float_t bx, Float_t by, Float_t bz) {
    if (!is_open_) return;
    if (idx < 0 || idx >= max_len_) return;

    Float_t* mag = static_cast<Float_t*>(&(geo_box_->mag) + idx * DIM_);
    mag[0] = bx;
    mag[1] = by;
    mag[2] = bz;
}


void MagGeoBoxCreator::save_and_close() {
    if (!is_open_) { clear(); return; }
    msync(fptr_, flen_, MS_SYNC);
    munmap(fptr_, flen_);
    close(fdes_);
    clear();
}
        

void MagGeoBoxCreator::save_and_close(Float_t bx, Float_t by, Float_t bz) {
    if (!is_open_) return;

    for (Long64_t idx = 0; idx < max_len_; ++idx) {
        Float_t* mag = static_cast<Float_t*>(&(geo_box_->mag) + idx * DIM_);
        mag[0] = bx;
        mag[1] = by;
        mag[2] = bz;
    }

    munmap(fptr_, flen_);
    close(fdes_);
    clear();
}


Bool_t MagGeoBoxReader::load(const std::string& fpath) {
    if (is_load_) {
        if (fpath_ == fpath) return is_load_;
        else {
            Sys::ShowWarning("MagGeoBoxReader::Load() : re-load different file.");
            clear();
        }
    }

    Int_t fdes = open(fpath.c_str(), O_RDONLY);
    Int_t flen = lseek(fdes, 0, SEEK_END); 
    if (fdes < 0) {
        Sys::ShowWarning(STR("MagGeoBoxReader::Load() : Magnetic field map not found (%s)", fpath.c_str()));
        is_load_ = false;
        return is_load_;
    }
    
    void* fptr = mmap(nullptr, flen, PROT_READ, MAP_SHARED, fdes, 0);
    if (fptr == reinterpret_cast<void*>(-1)) {
        Sys::ShowWarning("MagGeoBoxReader::Load() : mmap() failure.");
        close(fdes);
        is_load_ = false;
        return is_load_;
    }
    fptr_ = fptr;

    MagGeoBox* geo_box = reinterpret_cast<MagGeoBox*>(fptr_);
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
    fpath_ = fpath;
    COUT("MagGeoBoxReader::Load() : Open file (%s)\n", fpath.c_str());
    return is_load_;
}


MagFld MagGeoBoxReader::get(const SVecD<3>& coo) {
    if (!is_load_) return MagFld();
    
    // Get Local Coord
    Float_t xloc = (coo[0] - min_.at(0)) / dlt_.at(0);
    Float_t yloc = (coo[1] - min_.at(1)) / dlt_.at(1);
    Float_t zloc = (coo[2] - min_.at(2)) / dlt_.at(2);

    Long64_t xi = static_cast<Long64_t>(std::floor(xloc));
    Long64_t yi = static_cast<Long64_t>(std::floor(yloc));
    Long64_t zi = static_cast<Long64_t>(std::floor(zloc));

    if (xi < 0 || xi >= (n_.at(0)-1) || 
        yi < 0 || yi >= (n_.at(1)-1) || 
        zi < 0 || zi >= (n_.at(2)-1)) return MagFld();

    // Do Trilinear Interpolation
    Long64_t idx = (xi * fact_.at(0) + yi * fact_.at(1) + zi);
    Float_t xu = xloc - static_cast<Float_t>(xi);
    Float_t yu = yloc - static_cast<Float_t>(yi);
    Float_t zu = zloc - static_cast<Float_t>(zi);
    Float_t xl = (Numc::ONE<> - xu);
    Float_t yl = (Numc::ONE<> - yu);
    Float_t zl = (Numc::ONE<> - zu);

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

    return MagFld(c);
}


Bool_t MagMgnt::Load() {
    if (is_load_) return is_load_;

    if (!geo_box_reader_.exist() && Sys::IsEnv("TRACKSys_MagBox")) {
        //std::string fpath = "/data3/hchou/AMSData/MagDB/MagTest.bin";
        std::string fpath = Sys::GetEnv("TRACKSys_MagBox");
        geo_box_reader_.load(fpath);
    }
    is_load_ = (geo_box_reader_.exist());

    return is_load_;
}


} // namespace TrackSys


#endif // __TRACKLibs_MagEnv_C__
