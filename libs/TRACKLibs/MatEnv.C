#ifndef __TRACKLibs_MatEnv_C__
#define __TRACKLibs_MatEnv_C__


namespace TrackSys {


void MatFld::print() const {
    std::string printStr;
    printStr += STR("========================= MatFld =========================\n");
    printStr += STR("Mat       %-d\n", mat_);
    printStr += STR("InvRadLen %-8.5f\n", irl_);
    printStr += STR("ECloudDen %-8.5f\n", eld_);
    printStr += STR("LogMenExc %-8.5f\n", lme_);
    printStr += STR("DenEftCor %-8.5f\n", dec_);
    printStr += STR("NumRadLen %-8.5f\n", nrl());
    printStr += STR("ECloudAbs %-8.5f\n", ela());
    printStr += STR("RealLen   %-7.2f\n", rlen_);
    printStr += STR("EfftLen   %-7.2f\n", elen_);
    printStr += STR("Loc       %-6.4f\n", loc1_);
    printStr += STR("LocSqr    %-6.4f\n", loc2_);
    printStr += STR("==========================================================\n");
    COUT(printStr.c_str());
}
        

MatFld MatFld::Merge(const std::list<MatFld>& mflds) {
    if      (mflds.size() == 0) return MatFld();
    else if (mflds.size() == 1) return mflds.front();

    Bool_t   mat  = false;
    Double_t irl  = 0.;
    Double_t eld  = 0.;
    Double_t lme  = 0.;
    Double_t dec  = 0.;
    Double_t rlen = 0.;
    Double_t elen = 0.;
    Double_t loc1 = 0.;
    Double_t loc2 = 0.;
    
    for (auto&& mfld : mflds) {
        if (!mfld()) { rlen += mfld.rlen(); continue; }
        
        Double_t nrl = mfld.nrl();
        Double_t ela = mfld.ela();
        irl  += nrl;
        eld  += ela;
        lme  += mfld.lme() * ela;
        dec  += mfld.dec() * ela;
        
        Double_t loc1len = mfld.loc1() * mfld.rlen();
        Double_t loc2len = mfld.loc2() * mfld.rlen() * mfld.rlen();
        loc1 += (rlen + loc1len) * nrl;
        loc2 += (rlen * rlen + Numc::TWO<> * rlen * loc1len + loc2len) * nrl;
        elen += mfld.elen();
        rlen += mfld.rlen();

        mat = true;
    }

    if (mat && !Numc::EqualToZero(elen)) {
        loc1 = (loc1 / irl) / (rlen);
        loc2 = (loc2 / irl) / (rlen * rlen);
        lme  = (lme / eld);
        dec  = (dec / eld);
        irl  = (irl / elen);
        eld  = (eld / elen);
        
        return MatFld(mat, irl, eld, lme, dec, rlen, elen, loc1, loc2);
    }
    else return MatFld(rlen);
}


MatGeoBoxCreator::MatGeoBoxCreator(const Long64_t n[3], const Double_t min[3], const Double_t max[3], Double_t stp, const std::string& fname, const std::string& dpath) {
    clear();
    if (n[0] < 1 || n[1] < 1 || n[2] < 1) { MGSys::ShowError("MatGeoBoxCreator::MatGeoBoxCreator() : Size failure."); return; }
    if (Numc::Compare(min[0], max[0]) >= 0 || Numc::Compare(min[1], max[1]) >= 0 || Numc::Compare(min[2], max[2]) >= 0) { MGSys::ShowError("MatGeoBoxCreator::MatGeoBoxCreator() : Range failure."); return; }
    if (Numc::Compare(stp) <= 0) { MGSys::ShowError("MatGeoBoxCreator::MatGeoBoxCreator() : Step is negative or zero."); return; }
  
    // Inf
    Long64_t flen_inf = (MATGEOBOX_NDIM*(sizeof(Long64_t)+sizeof(Double_t)+sizeof(Double_t)) + sizeof(Double_t) + (n[0]*n[1]*n[2])*sizeof(Bool_t));
    Long64_t fdes_inf = open(CSTR("%s/%s.inf", dpath.c_str(), fname.c_str()), O_CREAT | O_RDWR, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
    if (fdes_inf < 0) { MGSys::ShowError("MatGeoBoxCreator::MatGeoBoxCreator() : File not opened."); return; }
    write(fdes_inf, "\0", 1);

    void* fptr_inf = mmap(nullptr, flen_inf, PROT_READ | PROT_WRITE, MAP_SHARED, fdes_inf, 0);
    if (fptr_inf == reinterpret_cast<void*>(-1)) { MGSys::ShowError("MatGeoBoxCreator::MatGeoBoxCreator() : mmap() failure."); close(fdes_inf); return; }

    // Var
    Long64_t flen_var = (MATGEOBOX_NPAR * (n[0]*n[1]*n[2]) * sizeof(Double_t));
    Long64_t fdes_var = open(CSTR("%s/%s.var", dpath.c_str(), fname.c_str()), O_CREAT | O_RDWR, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
    if (fdes_var < 0) { MGSys::ShowError("MatGeoBoxCreator::MatGeoBoxCreator() : File not opened."); return; }
    write(fdes_var, "\0", 1);

    void* fptr_var = mmap(nullptr, flen_var, PROT_READ | PROT_WRITE, MAP_SHARED, fdes_var, 0);
    if (fptr_var == reinterpret_cast<void*>(-1)) { MGSys::ShowError("MatGeoBoxCreator::MatGeoBoxCreator() : mmap() failure."); close(fdes_var); return; }
    
    // Create
    is_open_  = true;
    
    fdes_inf_ = fdes_inf;
    flen_inf_ = flen_inf;
    fptr_inf_ = fptr_inf;
    gbox_inf_ = reinterpret_cast<MatGeoBoxInf*>(fptr_inf);
    
    fdes_var_ = fdes_var;
    flen_var_ = flen_var;
    fptr_var_ = fptr_var;
    gbox_var_ = reinterpret_cast<MatGeoBoxVar*>(fptr_var);

    max_len_  = (n[0] * n[1] * n[2]);
    
    dlt_.fill(0);
    dlt_.at(0) = (max[0] - min[0]) / static_cast<Double_t>(n[0]);
    dlt_.at(1) = (max[1] - min[1]) / static_cast<Double_t>(n[1]);
    dlt_.at(2) = (max[2] - min[2]) / static_cast<Double_t>(n[2]);
    
    area_ = dlt_.at(0) * dlt_.at(1);
    
    fact_.fill(0);
    fact_.at(0) = n[1] * n[2];
    fact_.at(1) = n[2];

    gbox_inf_->n[0] = n[0];
    gbox_inf_->n[1] = n[1];
    gbox_inf_->n[2] = n[2];
    gbox_inf_->min[0] = min[0];
    gbox_inf_->min[1] = min[1];
    gbox_inf_->min[2] = min[2];
    gbox_inf_->max[0] = max[0];
    gbox_inf_->max[1] = max[1];
    gbox_inf_->max[2] = max[2];
    gbox_inf_->stp    = stp;
    std::fill_n(&(gbox_inf_->mat), max_len_, false);
   
    std::fill_n(&(gbox_var_->var), MATGEOBOX_NPAR * max_len_, 0.);

    mat_ptr_ = &(gbox_inf_->mat);
    var_ptr_ = &(gbox_var_->var);
}


void MatGeoBoxCreator::fill(const G4MatStep& g4mat) {
    if (!is_open_) return;
    if (g4mat.nstp <= 0) return;
    
    Long64_t xi = static_cast<Long64_t>(std::floor((g4mat.x - gbox_inf_->min[0]) / dlt_.at(0)));
    if (xi < 0 || xi >= gbox_inf_->n[0]) return;
    Long64_t yi = static_cast<Long64_t>(std::floor((g4mat.y - gbox_inf_->min[1]) / dlt_.at(1)));
    if (yi < 0 || yi >= gbox_inf_->n[1]) return;
   
    for (Long64_t istp = 0; istp < g4mat.nstp; ++istp) {
        Double_t g4min = g4mat.min->at(istp);
        Double_t g4max = g4mat.max->at(istp);
        if (g4min > gbox_inf_->max[2]) continue;
        if (g4max < gbox_inf_->min[2]) continue;

        Long64_t zimin = static_cast<Long64_t>(std::floor((g4min - gbox_inf_->min[2]) / dlt_.at(2)));
        Long64_t zimax = static_cast<Long64_t>(std::floor((g4max - gbox_inf_->min[2]) / dlt_.at(2)));
        if (zimin <= 0              ) zimin = 0;    
        if (zimax >= gbox_inf_->n[2]) zimax = gbox_inf_->n[2]-1;    

        Double_t scale = (g4mat.area / area_);
        Long64_t xyidx = (xi * fact_.at(0) + yi * fact_.at(1));
        for (Long64_t zi = zimin; zi <= zimax; ++zi) {
            Long64_t idx = (xyidx + zi);
            Double_t lbv = gbox_inf_->min[2] + dlt_.at(2) * (zi);
            Double_t ubv = gbox_inf_->min[2] + dlt_.at(2) * (zi+1);
            Double_t dbv = (ubv - lbv);

            if (g4min > lbv) lbv = g4min; 
            if (g4max < ubv) ubv = g4max;
            Double_t len = (ubv - lbv);
            Double_t eff = (len / dbv);
            Double_t scl = scale * eff;
            if (Numc::Compare(eff) <= 0) continue;

            Double_t irl  = scl * g4mat.irl->at(istp);
            Double_t eld  = scl * g4mat.eld->at(istp);
            Double_t lme  =       g4mat.lme->at(istp);
            Double_t dcC  =       g4mat.dcC->at(istp);
            Double_t dcM  =       g4mat.dcM->at(istp);
            Double_t dcA  =       g4mat.dcA->at(istp);
            Double_t dcX0 =       g4mat.dcX0->at(istp);
            Double_t dcX1 =       g4mat.dcX1->at(istp);

            mat_ptr_[idx] = true;
            var_ptr_[idx*MATGEOBOX_NPAR+MATVAR_IRL] += irl;
            var_ptr_[idx*MATGEOBOX_NPAR+MATVAR_ELD] += eld;
            var_ptr_[idx*MATGEOBOX_NPAR+MATVAR_LME] += lme * eld;
            var_ptr_[idx*MATGEOBOX_NPAR+MATVAR_C]   += dcC * eld;
            var_ptr_[idx*MATGEOBOX_NPAR+MATVAR_M]   += dcM * eld;
            var_ptr_[idx*MATGEOBOX_NPAR+MATVAR_A]   += dcA * eld;
            var_ptr_[idx*MATGEOBOX_NPAR+MATVAR_X0]  += dcX0 * eld;
            var_ptr_[idx*MATGEOBOX_NPAR+MATVAR_X1]  += dcX1 * eld;
        }
    }
}


void MatGeoBoxCreator::save_and_close() {
    if (!is_open_) { clear(); return; }
    
    for (Long64_t xi = 0; xi < gbox_inf_->n[0]; ++xi) {
    for (Long64_t yi = 0; yi < gbox_inf_->n[1]; ++yi) {
    for (Long64_t zi = 0; zi < gbox_inf_->n[2]; ++zi) {
        Long64_t idx = (xi * fact_.at(0) + yi * fact_.at(1) + zi);
        if (!mat_ptr_[idx]) continue;
        var_ptr_[idx*MATGEOBOX_NPAR+MATVAR_LME] /= var_ptr_[idx*MATGEOBOX_NPAR+MATVAR_ELD];
        var_ptr_[idx*MATGEOBOX_NPAR+MATVAR_C]   /= var_ptr_[idx*MATGEOBOX_NPAR+MATVAR_ELD];
        var_ptr_[idx*MATGEOBOX_NPAR+MATVAR_M]   /= var_ptr_[idx*MATGEOBOX_NPAR+MATVAR_ELD];
        var_ptr_[idx*MATGEOBOX_NPAR+MATVAR_A]   /= var_ptr_[idx*MATGEOBOX_NPAR+MATVAR_ELD];
        var_ptr_[idx*MATGEOBOX_NPAR+MATVAR_X0]  /= var_ptr_[idx*MATGEOBOX_NPAR+MATVAR_ELD];
        var_ptr_[idx*MATGEOBOX_NPAR+MATVAR_X1]  /= var_ptr_[idx*MATGEOBOX_NPAR+MATVAR_ELD];
    }}}

    munmap(fptr_inf_, flen_inf_);
    close(fdes_inf_);
    
    munmap(fptr_var_, flen_var_);
    close(fdes_var_);

    clear();
}
        

void MatGeoBoxReader::print() const {
    if (!is_load_) return;
    std::string printStr;
    printStr += STR("===================== MatGeoBoxReader ====================\n");
    printStr += STR("BOX X     (%3d %7.2f %7.2f)\n", n_.at(0), min_.at(0), max_.at(0));
    printStr += STR("BOX Y     (%3d %7.2f %7.2f)\n", n_.at(1), min_.at(1), max_.at(1));
    printStr += STR("BOX Z     (%3d %7.2f %7.2f)\n", n_.at(2), min_.at(2), max_.at(2));
    printStr += STR("==========================================================\n");
    COUT(printStr.c_str());
}


Bool_t MatGeoBoxReader::load(const std::string& fname, const std::string& dpath) {
    if (is_load_) return is_load_;
    //MGSys::ShowWarning("MatGeoBoxReader::Load() : re-load different file.");

    // Inf
    Int_t fdes_inf = open(CSTR("%s/%s.inf", dpath.c_str(), fname.c_str()), O_RDONLY);
    Int_t flen_inf = lseek(fdes_inf, 0, SEEK_END); 
    if (fdes_inf < 0) {
        MGSys::ShowError(STR("MatGeoBoxReader::Load() : Mat field map not found (%s)", fname.c_str()));
        is_load_ = false;
        return is_load_;
    }
    
    void* fptr_inf = mmap(nullptr, flen_inf, PROT_READ, MAP_SHARED, fdes_inf, 0);
    if (fptr_inf == reinterpret_cast<void*>(-1)) {
        MGSys::ShowError("MatGeoBoxReader::Load() : mmap() failure.");
        close(fdes_inf);
        is_load_ = false;
        return is_load_;
    }

    // Var
    Int_t fdes_var = open(CSTR("%s/%s.var", dpath.c_str(), fname.c_str()), O_RDONLY);
    Int_t flen_var = lseek(fdes_var, 0, SEEK_END); 
    if (fdes_var < 0) {
        MGSys::ShowError(STR("MatGeoBoxReader::Load() : Mat field map not found (%s)", fname.c_str()));
        is_load_ = false;
        return is_load_;
    }
    
    void* fptr_var = mmap(nullptr, flen_var, PROT_READ, MAP_SHARED, fdes_var, 0);
    if (fptr_var == reinterpret_cast<void*>(-1)) {
        MGSys::ShowError("MatGeoBoxReader::Load() : mmap() failure.");
        close(fdes_var);
        is_load_ = false;
        return is_load_;
    }
  
    // Inf
    MatGeoBoxInf* gbox_inf = reinterpret_cast<MatGeoBoxInf*>(fptr_inf);
    mat_ptr_ = &(gbox_inf->mat);
    
    for (Int_t ig = 0; ig < MATGEOBOX_NDIM; ++ig) {
        n_.at(ig)   = gbox_inf->n[ig]; 
        min_.at(ig) = gbox_inf->min[ig];
        max_.at(ig) = gbox_inf->max[ig];
        len_.at(ig) = (max_.at(ig) - min_.at(ig));
        dlt_.at(ig) = (len_.at(ig) / static_cast<Double_t>(n_.at(ig)));
    }
    max_len_ = n_.at(0) * n_.at(1) * n_.at(2);
    fact_.at(0) = n_.at(1) * n_.at(2);
    fact_.at(1) = n_.at(2);
    stp_ = gbox_inf->stp;

    // Var
    MatGeoBoxVar* gbox_var = reinterpret_cast<MatGeoBoxVar*>(fptr_var);
    var_ptr_ = &(gbox_var->var);

    is_load_ = true;
    COUT("MatGeoBoxReader::Load() : Open file (%s)\n", fname.c_str());
    return is_load_;
}


Double_t MatGeoBoxReader::get_density_effect_correction(Long64_t idx, Double_t log10gb) {
    if (idx < 0 || idx >= max_len_) return Numc::ZERO<>;
    if (!mat_ptr_[idx]) return Numc::ZERO<>;
    if (idx == tmp_dec_.first) return tmp_dec_.second;
    Double_t dec = Numc::ZERO<>;
    
    Double_t& X0 = var_ptr_[idx*MATGEOBOX_NPAR+MATVAR_X0];
    if (log10gb >= X0) {
        Double_t& C = var_ptr_[idx*MATGEOBOX_NPAR+MATVAR_C];
        dec += Numc::TWO<> * Numc::LOG_TEN<> * log10gb - C;
        
        Double_t& X1 = var_ptr_[idx*MATGEOBOX_NPAR+MATVAR_X1];
        if (log10gb < X1) { 
            Double_t& A = var_ptr_[idx*MATGEOBOX_NPAR+MATVAR_A];
            Double_t& M = var_ptr_[idx*MATGEOBOX_NPAR+MATVAR_M];
            dec += A * std::pow(X1 - log10gb, M); 
        }
    }
    
    if (!Numc::Valid(dec)) dec = Numc::ZERO<>;
    tmp_dec_.first  = idx;
    tmp_dec_.second = dec; 
    return dec;  
}

Bool_t MatGeoBoxReader::is_in_box(const SVecD<3>& coo) {
    if (!is_load_) return false;
    if (Numc::Compare(coo(2), min_.at(2)) <= 0 || Numc::Compare(coo(2), max_.at(2)) >= 0) return false;
    if (Numc::Compare(coo(1), min_.at(1)) <= 0 || Numc::Compare(coo(1), max_.at(1)) >= 0) return false;
    if (Numc::Compare(coo(0), min_.at(0)) <= 0 || Numc::Compare(coo(0), max_.at(0)) >= 0) return false;
    return true;
}


Bool_t MatGeoBoxReader::is_cross(const SVecD<3>& vcoo, const SVecD<3>& wcoo) {
    if (!is_load_) return false;
    if (Numc::Compare(vcoo(2), min_.at(2)) <= 0 && Numc::Compare(wcoo(2), min_.at(2)) <= 0) return false;
    if (Numc::Compare(vcoo(2), max_.at(2)) >= 0 && Numc::Compare(wcoo(2), max_.at(2)) >= 0) return false;
    if (Numc::Compare(vcoo(1), min_.at(1)) <= 0 && Numc::Compare(wcoo(1), min_.at(1)) <= 0) return false;
    if (Numc::Compare(vcoo(1), max_.at(1)) >= 0 && Numc::Compare(wcoo(1), max_.at(1)) >= 0) return false;
    if (Numc::Compare(vcoo(0), min_.at(0)) <= 0 && Numc::Compare(wcoo(0), min_.at(0)) <= 0) return false;
    if (Numc::Compare(vcoo(0), max_.at(0)) >= 0 && Numc::Compare(wcoo(0), max_.at(0)) >= 0) return false;
    return true;
}
        

MatFld MatGeoBoxReader::get(const SVecD<3>& coo, Double_t log10gb) {
    if (!is_load_) return MatFld();

    Long64_t zi = static_cast<Long64_t>(std::floor((coo(2) - min_.at(2)) / dlt_.at(2)));
    if (zi < 0 || zi >= n_.at(2)) return MatFld();
    
    Long64_t yi = static_cast<Long64_t>(std::floor((coo(1) - min_.at(1)) / dlt_.at(1)));
    if (yi < 0 || yi >= n_.at(1)) return MatFld();

    Long64_t xi = static_cast<Long64_t>(std::floor((coo(0) - min_.at(0)) / dlt_.at(0)));
    if (xi < 0 || xi >= n_.at(0)) return MatFld();

    Long64_t idx = (xi * fact_.at(0) + yi * fact_.at(1) + zi);
    Bool_t   mat = mat_ptr_[idx];
    Double_t irl = var_ptr_[idx*MATGEOBOX_NPAR+MATVAR_IRL];
    Double_t eld = var_ptr_[idx*MATGEOBOX_NPAR+MATVAR_ELD];
    Double_t lme = var_ptr_[idx*MATGEOBOX_NPAR+MATVAR_LME];
    Double_t dec = get_density_effect_correction(idx, log10gb);

    reset_tmp_dec();
    if (mat) return MatFld(mat, irl, eld, lme, dec);
    else     return MatFld();
}


MatFld MatGeoBoxReader::get(const SVecD<3>& vcoo, const SVecD<3>& wcoo, Double_t log10gb, Bool_t is_std) {
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

    Double_t rlen = LA::Mag((wcoo - vcoo));
    if (Numc::EqualToZero(rlen)) {
        Long64_t idx = (vxi * fact_.at(0) + vyi * fact_.at(1) + vzi);
        Bool_t   mat = mat_ptr_[idx];
        Double_t irl = var_ptr_[idx*MATGEOBOX_NPAR+MATVAR_IRL];
        Double_t eld = var_ptr_[idx*MATGEOBOX_NPAR+MATVAR_ELD];
        Double_t lme = var_ptr_[idx*MATGEOBOX_NPAR+MATVAR_LME];
        Double_t dec = get_density_effect_correction(idx, log10gb);
        
        reset_tmp_dec();
        if (mat) return MatFld(mat, irl, eld, lme, dec);
        else     return MatFld();
    }
   
    SVecD<3> vwvec((wxloc - vxloc), (wyloc - vyloc), (wzloc - vzloc));
    
    Double_t   vwlen = LA::Mag(vwvec);
    Long64_t   nstp  = static_cast<Long64_t>(std::floor((vwlen / stp_) / (is_std ? STD_STEP_LEN_ : FST_STEP_LEN_))) + 2;
    SVecD<3>&& unit  = (vwvec / static_cast<Double_t>(nstp));

    SVecD<3> itloc((vxloc + Numc::HALF<> * unit(0)), (vyloc + Numc::HALF<> * unit(1)), (vzloc + Numc::HALF<> * unit(2)));
    Long64_t itsat = 0;
    Long64_t itend = nstp;

    //==== faster method (tuned by axis-Z)
    Short_t uz_sign = Numc::Compare(unit(2));
    if (uz_sign != 0) {
        Double_t sat = ((uz_sign == 1) ? 0. : n_.at(2));
        Double_t end = ((uz_sign == 1) ? n_.at(2) : 0.);
        Long64_t satID = static_cast<Long64_t>(std::floor((sat - itloc(2)) / unit(2)));
        Long64_t endID = static_cast<Long64_t>(std::floor((end - itloc(2)) / unit(2))) + 1;
        if (Numc::Valid(satID) && Numc::Valid(endID)) {
            if (satID > 0 && satID < nstp) itsat = satID;
            if (endID > 0 && endID < nstp) itend = endID;
            itloc += (itsat * unit);
        }
    }
    //====

    Long64_t itcnt    = 0;
    Double_t sum_irl  = 0;
    Double_t sum_eld  = 0;
    Double_t sum_lme  = 0;
    Double_t sum_dec  = 0;
    Double_t sum_loc1 = 0;
    Double_t sum_loc2 = 0;
    for (Long64_t it = itsat; it < itend; ++it, itloc += unit) {
        Long64_t zi = static_cast<Long64_t>(std::floor(itloc(2)));
        if (zi < 0 || zi >= n_.at(2)) continue;
        Long64_t yi = static_cast<Long64_t>(std::floor(itloc(1)));
        if (yi < 0 || yi >= n_.at(1)) continue;
        Long64_t xi = static_cast<Long64_t>(std::floor(itloc(0)));
        if (xi < 0 || xi >= n_.at(0)) continue;

        Long64_t idx = (xi * fact_.at(0) + yi * fact_.at(1) + zi);
        Bool_t   mat = mat_ptr_[idx];
        if (mat) {
            Double_t itrat = ((Numc::HALF<> + static_cast<Double_t>(it)) / static_cast<Double_t>(nstp));
            Double_t irl = var_ptr_[idx*MATGEOBOX_NPAR+MATVAR_IRL];
            Double_t eld = var_ptr_[idx*MATGEOBOX_NPAR+MATVAR_ELD];
            Double_t lme = var_ptr_[idx*MATGEOBOX_NPAR+MATVAR_LME];
            Double_t dec = get_density_effect_correction(idx, log10gb);
    
            sum_irl  += irl;
            sum_eld  += eld;
            sum_lme  += lme * eld;
            sum_dec  += dec * eld;
            sum_loc1 += irl * (itrat);
            sum_loc2 += irl * (itrat * itrat);
            itcnt++;
        }
    }
    reset_tmp_dec();

    if (itcnt != 0) {
        Double_t lme  = (sum_lme / sum_eld);
        Double_t dec  = (sum_dec / sum_eld);
        Double_t irl  = (sum_irl / static_cast<Double_t>(itcnt));
        Double_t eld  = (sum_eld / static_cast<Double_t>(itcnt));
        Double_t efft = (static_cast<Double_t>(itcnt) / static_cast<Double_t>(nstp));
        Double_t elen = efft * rlen;
        Double_t loc1 = (sum_loc1 / sum_irl);
        Double_t loc2 = (sum_loc2 / sum_irl);

        return MatFld(true, irl, eld, lme, dec, rlen, elen, loc1, loc2);
    }
    else return MatFld(rlen);
}


MatFld MatMgnt::Get(const SVecD<3>& coo, Double_t log10gb) {
    if (!Load()) return MatFld();

    Bool_t   mat  = false;
    Double_t irl  = 0.0;
    Double_t eld  = 0.0;
    Double_t lme  = 0.0;
    Double_t dec  = 0.0;

    for (auto&& reader : *reader_) {
        if (!reader->is_in_box(coo)) continue;
        MatFld&& mfld = reader->get(coo, log10gb);
        if (!mfld()) continue;

        irl  += mfld.irl();
        eld  += mfld.eld();
        lme  += mfld.lme() * mfld.eld();
        dec  += mfld.dec() * mfld.eld();
        mat = true;
    }
    if (mat) { 
        lme  = (lme / eld); 
        dec  = (dec / eld); 
    }

    if (mat) return MatFld(mat, irl, eld, lme, dec);
    else     return MatFld();
}
    

MatFld MatMgnt::Get(const SVecD<3>& vcoo, const SVecD<3>& wcoo, Double_t log10gb, Bool_t is_std) {
    if (!Load()) return MatFld();
    
    Double_t rlen = LA::Mag(wcoo - vcoo);
    if (Numc::EqualToZero(rlen)) return Get(vcoo);
    
    Bool_t   mat  = false;
    Double_t irl  = 0.0;
    Double_t eld  = 0.0;
    Double_t lme  = 0.0;
    Double_t dec  = 0.0;
    Double_t elen = 0.0;
    Double_t loc1 = 0.0;
    Double_t loc2 = 0.0;
        
    for (auto&& reader : *reader_) {
        if (!reader->is_cross(vcoo, wcoo)) continue;
        MatFld&& mfld = reader->get(vcoo, wcoo, log10gb, is_std);
        if (!mfld()) continue;
        
        Double_t nrl = mfld.nrl();
        Double_t ela = mfld.ela();
        loc1 += mfld.loc1() * nrl;
        loc2 += mfld.loc2() * nrl;
        irl  += nrl;
        eld  += ela;
        lme  += mfld.lme() * ela;
        dec  += mfld.dec() * ela;
        elen += mfld.elen();
        
        mat = true;
    }

    if (mat && !Numc::EqualToZero(elen)) {
        loc1 = (loc1 / irl);
        loc2 = (loc2 / irl);
        lme  = (lme / eld);
        dec  = (dec / eld);
        irl  = (irl / elen);
        eld  = (eld / elen);

        return MatFld(mat, irl, eld, lme, dec, rlen, elen, loc1, loc2);
    }
    else return MatFld(rlen);
}


MatFld MatMgnt::Get(Double_t stp_len, const PhySt& part, Bool_t is_std) {
    const SVecD<3>&  vcoo = part.c();
    SVecD<3>&&       wcoo = part.c() + stp_len * part.u();
    Double_t log10gb = std::log10(part.gmbta());
    if (!Numc::Valid(log10gb)) log10gb = Numc::ZERO<>;

    return Get(vcoo, wcoo, log10gb, is_std);
}


MatPhyFld MatPhy::Get(const Double_t stp_len, PhySt& part, Bool_t is_std) {
    if (!part.field()) return MatPhyFld();
    if (part.info().is_chrgless() || part.info().is_massless()) return MatPhyFld();
    if (Numc::EqualToZero(stp_len)) return MatPhyFld();
    if (Numc::EqualToZero(part.mom())) return MatPhyFld();

    const SVecD<3>&  vcoo = part.c();
    SVecD<3>&&       wcoo = part.c() + stp_len * part.u();
    Double_t log10gb = std::log10(part.gmbta());
    if (!Numc::Valid(log10gb)) log10gb = Numc::ZERO<>;

    MatFld&& mfld = MatMgnt::Get(vcoo, wcoo, log10gb, is_std);
    
    if (!mfld()) return MatPhyFld();

    Double_t mscat_sgm = GetMultipleScattering(mfld, part);
    std::tuple<Double_t, Double_t, Double_t>&& ion_eloss = GetIonizationEnergyLoss(mfld, part);
    Double_t eloss_brm_men = GetBremsstrahlungEnergyLoss(mfld, part);

    return MatPhyFld(mfld(), mscat_sgm, std::get<0>(ion_eloss), std::get<1>(ion_eloss), std::get<2>(ion_eloss), eloss_brm_men);
}


MatPhyFld MatPhy::Get(const MatFld& mfld, PhySt& part) {
    if (!mfld() || !part.field()) return MatPhyFld();
    if (Numc::EqualToZero(mfld.elen())) return MatPhyFld();
    if (part.info().is_chrgless() || part.info().is_massless()) return MatPhyFld();
    if (Numc::EqualToZero(part.mom())) return MatPhyFld();
    
    Double_t mscat_sgm = GetMultipleScattering(mfld, part);
    std::tuple<Double_t, Double_t, Double_t>&& ion_eloss = GetIonizationEnergyLoss(mfld, part);
    Double_t eloss_brm_men = GetBremsstrahlungEnergyLoss(mfld, part);
    
    return MatPhyFld(mfld(), mscat_sgm, std::get<0>(ion_eloss), std::get<1>(ion_eloss), std::get<2>(ion_eloss), eloss_brm_men);
}
        

Double_t MatPhy::GetMultipleScattering(const MatFld& mfld, PhySt& part) {
    if (!part.arg().mscat()) return Numc::ZERO<>;
    
    Bool_t is_over_lmt = (Numc::Compare(part.bta(), LMT_BTA) > 0);
    Double_t eta_part = ((is_over_lmt) ? (part.eta_abs() / part.bta()) : (LMT_INV_GMBTA / LMT_BTA));

    Double_t nrl     = mfld.nrl();
    Double_t sqr_nrl = std::sqrt(nrl);
    Double_t log_nrl = ((corr_sw_mscat_) ? std::log(corr_mfld_.nrl()) : std::log(nrl));
    
    // Highland-Lynch-Dahl formula
    //Double_t mscat_sgm = RYDBERG_CONST * part.info().chrg_to_mass() * eta_part * sqr_nrl * (Numc::ONE<> + NRL_CORR_FACT * log_nrl);
    
    // Modified Highland-Lynch-Dahl formula
    Double_t corr_fact = (1.02246 + 0.0282457 * TMath::Erfc(3.38323 * (part.bta() - 0.691661))); // testcode
    Double_t mscat_sgm = corr_fact * RYDBERG_CONST * part.info().chrg_to_mass() * eta_part * sqr_nrl * std::sqrt(Numc::ONE<> + NRL_CORR_FACT1 * log_nrl + NRL_CORR_FACT2 * log_nrl * log_nrl);
   
    if (!Numc::Valid(mscat_sgm) || Numc::Compare(mscat_sgm) <= 0) mscat_sgm = Numc::ZERO<>;
    return mscat_sgm;
}


std::tuple<Double_t, Double_t, Double_t> MatPhy::GetIonizationEnergyLoss(const MatFld& mfld, PhySt& part) {
    if (!part.arg().eloss()) return std::make_tuple(Numc::ZERO<>, Numc::ZERO<>, Numc::ZERO<>);

    Bool_t is_over_lmt   = (Numc::Compare(part.bta(), LMT_BTA) > 0);
    Double_t sqr_gmbta   = ((is_over_lmt) ? (part.gmbta() * part.gmbta()) : LMT_SQR_GMBTA);
    Double_t sqr_bta     = ((is_over_lmt) ? (part.bta() * part.bta()) : LMT_SQR_BTA);
    Double_t gm          = ((is_over_lmt) ? part.gm() : LMT_GM);
    Double_t sqr_chrg    = part.chrg() * part.chrg();
    Double_t mass_in_GeV = part.mass();
    Double_t mass_in_MeV = part.mass() * GEV_TO_MEV;

    // Calculate Matterial Quality and Eta Trans
    Double_t log_mean_exc_eng  = mfld.lme(); // log[MeV]
    Double_t elcloud_abundance = mfld.ela(); // [mol cm^-2]
    Double_t density_corr      = mfld.dec(); // [1]
    Double_t Bethe_Bloch       = (Numc::HALF<> * BETHE_BLOCH_K * elcloud_abundance * sqr_chrg / sqr_bta); // [MeV]

    // Calculate Sigma
    Double_t eta_trans = (std::sqrt(sqr_gmbta + Numc::ONE<>) / sqr_gmbta); // ke to eta
    Double_t elion_sgm = ((Bethe_Bloch / mass_in_MeV) * eta_trans); // [1]
    
    // Calculate Peak
    Double_t corr_fact  = ((corr_sw_eloss_) ? (corr_mfld_.ela() / elcloud_abundance) : Numc::ONE<>);
    Double_t trans_eng  = Numc::TWO<> * MASS_EL_IN_MEV * sqr_gmbta; // [MeV]
    Double_t maxke_part = std::log(trans_eng) - log_mean_exc_eng;   // [1]
    Double_t ionke_part = std::log(Bethe_Bloch * corr_fact) - log_mean_exc_eng; // [1]
    Double_t elion_mpv  = elion_sgm * (maxke_part + ionke_part + LANDAU_ELOSS_CORR - sqr_bta - density_corr); //[1]
    
    // Calculate Mean
    Double_t mass_rat     = (MASS_EL_IN_GEV / mass_in_GeV);
    Double_t mass_rel     = (Numc::ONE<> + mass_rat * (Numc::TWO<> * gm + mass_rat));
    Double_t transke_part = Numc::TWO<> * maxke_part - std::log(mass_rel);
    Double_t elion_men    = elion_sgm * (transke_part - Numc::TWO<> * sqr_bta - density_corr); // [1]

    if (!Numc::Valid(elion_mpv) || Numc::Compare(elion_mpv) <= 0) elion_mpv = Numc::ZERO<>;
    if (!Numc::Valid(elion_sgm) || Numc::Compare(elion_sgm) <= 0) elion_sgm = Numc::ZERO<>;
    if (!Numc::Valid(elion_men) || Numc::Compare(elion_men) <= 0) elion_men = Numc::ZERO<>;

    return std::make_tuple(elion_mpv, elion_sgm, elion_men);
}


Double_t MatPhy::GetBremsstrahlungEnergyLoss(const MatFld& mfld, PhySt& part) {
    if (!part.arg().eloss()) return Numc::ZERO<>;
    Double_t chrgmass_sqr = (MASS_EL_IN_GEV * part.info().chrg_to_mass()) * (MASS_EL_IN_GEV * part.info().chrg_to_mass());

    Bool_t   is_over_lmt = (Numc::Compare(part.bta(), LMT_BTA) > 0);
    Double_t sqr_gmbta   = ((is_over_lmt) ? (part.gmbta() * part.gmbta()) : LMT_SQR_GMBTA);
    Double_t eng         = std::sqrt(sqr_gmbta + Numc::ONE<>);
    Double_t ke_part     = (eng - Numc::ONE<>);
    Double_t eta_trans   = (eng / sqr_gmbta);

    Double_t elbrm_men = chrgmass_sqr * (ke_part * eta_trans) * (mfld.nrl() / Numc::LOG_TWO<>);
    
    if (!Numc::Valid(elbrm_men) || Numc::Compare(elbrm_men) <= 0) elbrm_men = Numc::ZERO<>;
    return elbrm_men;
}


} // namespace TrackSys


#endif // __TRACKLibs_MatEnv_C__
