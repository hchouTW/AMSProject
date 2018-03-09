#if defined(_PGTRACK_) || defined(__ROOTSHAREDLIBRARY__)
#ifndef __TRACKLibs_MagEnvAms_C__
#define __TRACKLibs_MagEnvAms_C__


namespace TrackSys {


Bool_t MagGeoBoxAms::Load() {
    if (is_load_) return is_load_;
    is_load_ = false;
    
    mag_field_ = MagField::GetPtr();
    if (mag_field_ == nullptr) return is_load_; 
    
    if (mag_field_->GetMap()) is_load_ = true;
    else {
        std::string fdir = Sys::GetEnv("AMSDataDir");
        std::string fpath = STR("%s/v5.00/MagneticFieldMapPermanent_NEW_FULL.bin", fdir.c_str());
        if (fdir == "" || (mag_field_->Read(fpath.c_str())) < 0) {
            CERR("Magnetic Field map not found : %s\n", fpath.c_str());
            mag_field_ = nullptr;
        }
        else {
            COUT("MagGeoBoxAms::Load() Had file : %s\n", fpath.c_str());
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
        

void MagGeoBoxAms::Output(const std::string& fpath) {
    if (!Load()) return;
    const Long64_t n[3]   = {    201,    201,    201 };
    const Double_t min[3] = { -200.0, -200.0, -200.0 };
    const Double_t max[3] = {  200.0,  200.0,  200.0 };
    const Double_t dlt    = 2.0;
    MagGeoBoxCreator creator(n, min, max, fpath);
    if (!creator.is_open()) return;

    COUT("\n");
    COUT("MagGeoBoxAms::Output() : Start.\n");
    COUT("MagGeoBoxAms::Output() : X (%lld %8.2f %8.2f)\n", n[0], min[0], max[0]);
    COUT("MagGeoBoxAms::Output() : Y (%lld %8.2f %8.2f)\n", n[1], min[1], max[1]);
    COUT("MagGeoBoxAms::Output() : Z (%lld %8.2f %8.2f)\n", n[2], min[2], max[2]);

    Long64_t fact[2] { (n[2] * n[1]), (n[2]) };
    for (Long64_t xi = 0; xi < n[0]; ++xi) {
        for (Long64_t yi = 0; yi < n[1]; ++yi) {
            for (Long64_t zi = 0; zi < n[2]; ++zi) {
                Long64_t idx = xi * fact[0] + yi * fact[1] + zi;
                SVecD<3> coo((min[0] + static_cast<Double_t>(xi) * dlt), (min[1] + static_cast<Double_t>(yi) * dlt), (min[2] + static_cast<Double_t>(zi) * dlt));
                MagFld&& mag = MagGeoBoxAms::Get(coo);
                creator.fill(idx, static_cast<Float_t>(mag.x()), static_cast<Float_t>(mag.y()), static_cast<Float_t>(mag.z()));
            }
        }
    }
    COUT("MagGeoBoxAms::Output() : Save.\n");
    creator.save_and_close();
    
    COUT("MagGeoBoxAms::Output() : Finish.\n");
    COUT("\n");
}

} // namespace TrackSys


#endif // __TRACKLibs_MagEnvAms_C__
#endif // _PGTRACK_ __ROOTSHAREDLIBRARY__ 
