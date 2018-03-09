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

} // namespace TrackSys


#endif // __TRACKLibs_MagEnvAms_C__
