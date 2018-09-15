#ifndef __TRACKLibs_MatEnv_C__
#define __TRACKLibs_MatEnv_C__


#include "Sys.h"
#include "Math.h"
#include "PartInfo.h"
#include "PhySt.h"
#include "MatEnv.h"


namespace TrackSys {


void MatFld::print() const {
    std::string printStr;
    printStr += STR("================= MatFld =================\n");
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
    printStr += STR("==========================================\n");
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


MatGeoBoxCreator::MatGeoBoxCreator(const std::array<Long64_t, 3>& n, const std::array<Double_t, 3>& min, const std::array<Double_t, 3>& max, Double_t stp, const std::string& fname, const std::string& dpath) {
    clear();
    if (n[0] < 1 || n[1] < 1 || n[2] < 1) { 
        Sys::ShowWarning("MatGeoBoxCreator::MatGeoBoxCreator() : Size failure."); return; 
    }
    if (Numc::Compare(min[0], max[0]) >= 0 || Numc::Compare(min[1], max[1]) >= 0 || Numc::Compare(min[2], max[2]) >= 0) { 
        Sys::ShowWarning("MatGeoBoxCreator::MatGeoBoxCreator() : Range failure."); return; 
    }
    if (Numc::Compare(stp) <= 0) { 
        Sys::ShowWarning("MatGeoBoxCreator::MatGeoBoxCreator() : Step is negative or zero."); return; 
    }
  
    // Inf
    Long64_t flen_inf = (MATGEOBOX_NDIM*(sizeof(Long64_t)+sizeof(Double_t)+sizeof(Double_t)) + sizeof(Double_t) + (n[0]*n[1]*n[2])*sizeof(Bool_t));
    Long64_t fdes_inf = open(CSTR("%s/%s.inf", dpath.c_str(), fname.c_str()), O_RDWR | O_CREAT | O_TRUNC, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
    if (fdes_inf < 0) { Sys::ShowWarning("MatGeoBoxCreator::MatGeoBoxCreator() : File not opened."); return; }
    if (lseek(fdes_inf, flen_inf, SEEK_SET) == -1) { 
        Sys::ShowWarning("MatGeoBoxCreator::MatGeoBoxCreator() : Error calling lseek() to 'stretch' the file"); close(fdes_inf); return; 
    }
    if (write(fdes_inf, "", 1) == -1) {
        Sys::ShowWarning("MatGeoBoxCreator::MatGeoBoxCreator() : Error writing last byte of the file"); close(fdes_inf); return; 
    }

    void* fptr_inf = mmap(nullptr, flen_inf, PROT_READ | PROT_WRITE, MAP_SHARED, fdes_inf, 0);
    if (fptr_inf == reinterpret_cast<void*>(-1)) { 
        Sys::ShowWarning("MatGeoBoxCreator::MatGeoBoxCreator() : mmap() failure."); close(fdes_inf); return; 
    }

    // Var
    Long64_t flen_var = (MATGEOBOX_NPAR * (n[0]*n[1]*n[2]) * sizeof(Double_t));
    Long64_t fdes_var = open(CSTR("%s/%s.var", dpath.c_str(), fname.c_str()), O_RDWR | O_CREAT | O_TRUNC, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
    if (fdes_var < 0) { Sys::ShowWarning("MatGeoBoxCreator::MatGeoBoxCreator() : File not opened."); return; }
    if (lseek(fdes_var, flen_var, SEEK_SET) == -1) {
        Sys::ShowWarning("MatGeoBoxCreator::MatGeoBoxCreator() : Error calling lseek() to 'stretch' the file"); close(fdes_var); return; 
    } 
    if (write(fdes_var, "", 1) == -1) {
        Sys::ShowWarning("MatGeoBoxCreator::MatGeoBoxCreator() : Error writing last byte of the file"); close(fdes_var); return; 
    }

    void* fptr_var = mmap(nullptr, flen_var, PROT_READ | PROT_WRITE, MAP_SHARED, fdes_var, 0);
    if (fptr_var == reinterpret_cast<void*>(-1)) { 
        Sys::ShowWarning("MatGeoBoxCreator::MatGeoBoxCreator() : mmap() failure."); close(fdes_var); return; 
    }
    
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
    dlt_[0] = (max[0] - min[0]) / static_cast<Double_t>(n[0]);
    dlt_[1] = (max[1] - min[1]) / static_cast<Double_t>(n[1]);
    dlt_[2] = (max[2] - min[2]) / static_cast<Double_t>(n[2]);
    
    area_ = dlt_[0] * dlt_[1];
    
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
    
    Long64_t xi = static_cast<Long64_t>(std::floor((g4mat.x - gbox_inf_->min[0]) / dlt_[0]));
    if (xi < 0 || xi >= gbox_inf_->n[0]) return;
    Long64_t yi = static_cast<Long64_t>(std::floor((g4mat.y - gbox_inf_->min[1]) / dlt_[1]));
    if (yi < 0 || yi >= gbox_inf_->n[1]) return;
   
    for (Long64_t istp = 0; istp < g4mat.nstp; ++istp) {
        Double_t g4min = g4mat.min->at(istp);
        Double_t g4max = g4mat.max->at(istp);
        if (g4min > gbox_inf_->max[2]) continue;
        if (g4max < gbox_inf_->min[2]) continue;

        Long64_t zimin = static_cast<Long64_t>(std::floor((g4min - gbox_inf_->min[2]) / dlt_[2]));
        Long64_t zimax = static_cast<Long64_t>(std::floor((g4max - gbox_inf_->min[2]) / dlt_[2]));
        if (zimin <= 0              ) zimin = 0;    
        if (zimax >= gbox_inf_->n[2]) zimax = gbox_inf_->n[2]-1;    

        Double_t scale = (g4mat.area / area_);
        Long64_t xyidx = (xi * fact_.at(0) + yi * fact_.at(1));
        for (Long64_t zi = zimin; zi <= zimax; ++zi) {
            Long64_t idx = (xyidx + zi);
            Double_t lbv = gbox_inf_->min[2] + dlt_[2] * (zi);
            Double_t ubv = gbox_inf_->min[2] + dlt_[2] * (zi+1);
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

    if (msync(fptr_inf_, flen_inf_, MS_SYNC) == -1) {
        Sys::ShowWarning("MatGeoBoxCreator::save_and_close() : Could not sync the file to disk"); return; 
    }
    if (munmap(fptr_inf_, flen_inf_) == -1) {
        Sys::ShowWarning("MatGeoBoxCreator::save_and_close() : Error un-mmapping the file"); close(fdes_inf_); return; 
    }
    close(fdes_inf_);
    
    if (msync(fptr_var_, flen_var_, MS_SYNC) == -1) {
        Sys::ShowWarning("MatGeoBoxCreator::save_and_close() : Could not sync the file to disk"); return; 
    }
    if (munmap(fptr_var_, flen_var_) == -1) {
        Sys::ShowWarning("MatGeoBoxCreator::save_and_close() : Error un-mmapping the file"); close(fdes_var_); return; 
    }
    close(fdes_var_);

    clear();
}
        

void MatGeoBoxReader::print() const {
    if (!is_load_) return;
    std::string printStr;
    printStr += STR("===================== MatGeoBoxReader ====================\n");
    printStr += STR("BOX X     (%3lld %7.2f %7.2f)\n", n_.at(0), min_[0], max_[0]);
    printStr += STR("BOX Y     (%3lld %7.2f %7.2f)\n", n_.at(1), min_[1], max_[1]);
    printStr += STR("BOX Z     (%3lld %7.2f %7.2f)\n", n_.at(2), min_[2], max_[2]);
    printStr += STR("==========================================================\n");
    COUT(printStr.c_str());
}


Bool_t MatGeoBoxReader::load(const std::string& fname, const std::string& dpath) {
    if (is_load_) return is_load_;
    is_load_ = false;

    // Inf
    Long64_t fdes_inf = open(CSTR("%s/%s.inf", dpath.c_str(), fname.c_str()), O_RDONLY);
    if (fdes_inf < 0) {
        Sys::ShowWarningExit(STR("MatGeoBoxReader::Load() : Mat field map not found (%s/%s.inf)", dpath.c_str(), fname.c_str()));
        return is_load_;
    }
    Long64_t flen_inf = lseek(fdes_inf, 0, SEEK_END); 
    
    void* fptr_inf = mmap(nullptr, flen_inf, PROT_READ, MAP_PRIVATE, fdes_inf, 0);
    if (fptr_inf == reinterpret_cast<void*>(-1)) {
        Sys::ShowWarningExit("MatGeoBoxReader::Load() : mmap() failure.");
        close(fdes_inf);
        return is_load_;
    }

    // Var
    Long64_t fdes_var = open(CSTR("%s/%s.var", dpath.c_str(), fname.c_str()), O_RDONLY);
    if (fdes_var < 0) {
        Sys::ShowWarningExit(STR("MatGeoBoxReader::Load() : Mat field map not found (%s/%s.var)", dpath.c_str(), fname.c_str()));
        return is_load_;
    }
    Long64_t flen_var = lseek(fdes_var, 0, SEEK_END); 
    
    void* fptr_var = mmap(nullptr, flen_var, PROT_READ, MAP_PRIVATE, fdes_var, 0);
    if (fptr_var == reinterpret_cast<void*>(-1)) {
        Sys::ShowWarningExit("MatGeoBoxReader::Load() : mmap() failure.");
        close(fdes_var);
        return is_load_;
    }
  
    // Inf
    MatGeoBoxInf* gbox_inf = reinterpret_cast<MatGeoBoxInf*>(fptr_inf);
    Bool_t* mat_ptr = &(gbox_inf->mat);
    
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
    Double_t* var_ptr = &(gbox_var->var);
        
    // Load To Memory
    mat_ = std::vector<Bool_t>(max_len_, false);
    for (Long64_t it = 0; it < mat_.size(); ++it)
        mat_[it] = mat_ptr[it];

    var_ = std::vector<std::array<Double_t, MATGEOBOX_NPAR>>(max_len_, std::array<Double_t, MATGEOBOX_NPAR>());
    for (Long64_t it = 0; it < var_.size(); ++it) {
        var_.at(it).fill(Numc::ZERO<>);
        for (Long64_t ipar = 0; ipar < MATGEOBOX_NPAR; ++ipar)
            var_[it][ipar] = var_ptr[it*MATGEOBOX_NPAR+ipar];
    }
    
    // Release Memory and File
    if (munmap(fptr_inf, flen_inf) == -1) {
        Sys::ShowWarningExit("MatGeoBoxReader::Load() : Error un-mmapping the file");
        if (fdes_inf >= 0) close(fdes_inf);
        return is_load_;
    }
    if (fdes_inf >= 0) close(fdes_inf);
    
    if (munmap(fptr_var, flen_var) == -1) {
        Sys::ShowWarningExit("MatGeoBoxReader::Load() : Error un-mmapping the file");
        if (fdes_var >= 0) close(fdes_var);
        return is_load_;
    }
    if (fdes_var >= 0) close(fdes_var);

    is_load_ = true;
    COUT("MatGeoBoxReader::Load() : Open file (%s/%s)\n", dpath.c_str(), fname.c_str());
    return is_load_;
}


MatFld MatGeoBoxReader::get(const SVecD<3>& coo, Double_t log10gb) {
    if (!is_load_) return MatFld();

    Double_t zloc = (coo(2) - min_[2]) / dlt_[2];
    Long64_t zi = static_cast<Long64_t>(std::floor(zloc));
    if (!Numc::Valid(zloc) || !Numc::Valid(zi) || zi < 0 || zi >= n_.at(2)) return MatFld();
    
    Double_t yloc = (coo(1) - min_[1]) / dlt_[1];
    Long64_t yi = static_cast<Long64_t>(std::floor(yloc));
    if (!Numc::Valid(yloc) || !Numc::Valid(yi) || yi < 0 || yi >= n_.at(1)) return MatFld();

    Double_t xloc = (coo(0) - min_[0]) / dlt_[0];
    Long64_t xi = static_cast<Long64_t>(std::floor(xloc));
    if (!Numc::Valid(xloc) || !Numc::Valid(xi) || xi < 0 || xi >= n_.at(0)) return MatFld();

    Long64_t idx = (xi * fact_.at(0) + yi * fact_.at(1) + zi);
    Bool_t   mat = mat_[idx];
    Double_t irl = var_[idx][MATVAR_IRL];
    Double_t eld = var_[idx][MATVAR_ELD];
    Double_t lme = var_[idx][MATVAR_LME];
    Double_t dec = get_density_effect_correction(idx, log10gb);

    reset_tmp_dec();
    if (mat) return MatFld(mat, irl, eld, lme, dec);
    else     return MatFld();
}


MatFld MatGeoBoxReader::get(const SVecD<3>& vcoo, const SVecD<3>& wcoo, Double_t log10gb) {
    if (!is_load_) return MatFld();
    
    Double_t vzloc = (vcoo(2) - min_[2]) / dlt_[2];
    Double_t wzloc = (wcoo(2) - min_[2]) / dlt_[2];
    Long64_t vzi = static_cast<Long64_t>(std::floor(vzloc));
    Long64_t wzi = static_cast<Long64_t>(std::floor(wzloc));
    if (!Numc::Valid(vzloc) || !Numc::Valid(vzi)) return MatFld();
    if (!Numc::Valid(wzloc) || !Numc::Valid(wzi)) return MatFld();
    if ((vzi < 0 && wzi < 0) || (vzi >= n_.at(2) && wzi >= n_.at(2))) return MatFld();
    
    Double_t vyloc = (vcoo(1) - min_[1]) / dlt_[1];
    Double_t wyloc = (wcoo(1) - min_[1]) / dlt_[1];
    Long64_t vyi = static_cast<Long64_t>(std::floor(vyloc));
    Long64_t wyi = static_cast<Long64_t>(std::floor(wyloc));
    if (!Numc::Valid(vyloc) || !Numc::Valid(vyi)) return MatFld();
    if (!Numc::Valid(wyloc) || !Numc::Valid(wyi)) return MatFld();
    if ((vyi < 0 && wyi < 0) || (vyi >= n_.at(1) && wyi >= n_.at(1))) return MatFld();
    
    Double_t vxloc = (vcoo(0) - min_[0]) / dlt_[0];
    Double_t wxloc = (wcoo(0) - min_[0]) / dlt_[0];
    Long64_t vxi = static_cast<Long64_t>(std::floor(vxloc));
    Long64_t wxi = static_cast<Long64_t>(std::floor(wxloc));
    if (!Numc::Valid(vxloc) || !Numc::Valid(vxi)) return MatFld();
    if (!Numc::Valid(wxloc) || !Numc::Valid(wxi)) return MatFld();
    if ((vxi < 0 && wxi < 0) || (vxi >= n_.at(0) && wxi >= n_.at(0))) return MatFld();

    Double_t rlen = LA::Mag((wcoo - vcoo));
    if (Numc::EqualToZero(rlen)) {
        Long64_t idx = (vxi * fact_.at(0) + vyi * fact_.at(1) + vzi);
        Bool_t   mat = mat_[idx];
        Double_t irl = var_[idx][MATVAR_IRL];
        Double_t eld = var_[idx][MATVAR_ELD];
        Double_t lme = var_[idx][MATVAR_LME];
        Double_t dec = get_density_effect_correction(idx, log10gb);
        
        reset_tmp_dec();
        if (mat) return MatFld(mat, irl, eld, lme, dec);
        else     return MatFld();
    }
   
    SVecD<3> vwvec((wxloc - vxloc), (wyloc - vyloc), (wzloc - vzloc));
    
    Double_t   vwlen = LA::Mag(vwvec);
    Long64_t   nstp  = static_cast<Long64_t>(std::floor((vwlen / stp_) / STD_STEP_LEN)) + 2;
    SVecD<3>&& unit  = (vwvec / static_cast<Double_t>(nstp));

    SVecD<3> itloc((vxloc + Numc::HALF * unit(0)), (vyloc + Numc::HALF * unit(1)), (vzloc + Numc::HALF * unit(2)));
    Long64_t itsat = 0;
    Long64_t itend = nstp;

    //==== faster method (tuned by axis-ZYX)
    Short_t iu    = -1;
    Short_t signu =  0;
    const Double_t LMTU = 5.0e-2;
    if      (Numc::Compare(std::fabs(unit(2)), LMTU) > 0) { iu = 2; signu = Numc::Compare(unit(2)); }
    else if (Numc::Compare(std::fabs(unit(1)), LMTU) > 0) { iu = 1; signu = Numc::Compare(unit(1)); }
    else if (Numc::Compare(std::fabs(unit(0)), LMTU) > 0) { iu = 0; signu = Numc::Compare(unit(0)); }
    if (iu >= 0 && signu != 0) {
        Double_t sat = ((signu == 1) ? Numc::ZERO<> :    n_.at(iu));
        Double_t end = ((signu == 1) ?    n_.at(iu) : Numc::ZERO<>);
        sat = std::floor((sat - itloc(iu)) / unit(iu));
        end = std::floor((end - itloc(iu)) / unit(iu));
        if (Numc::Valid(sat) && Numc::Valid(end)) {
            Long64_t satID = static_cast<Long64_t>(sat);
            Long64_t endID = static_cast<Long64_t>(end) + Numc::ONE<Long64_t>;
            if (Numc::Valid(satID) && Numc::Valid(endID)) {
                if (satID > 0 && satID < nstp) itsat = satID;
                if (endID > 0 && endID < nstp) itend = endID;
                itloc += (itsat * unit);
            }
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
        Bool_t   mat = mat_[idx];
        if (mat) {
            Double_t itrat = ((Numc::HALF + static_cast<Double_t>(it)) / static_cast<Double_t>(nstp));
            Double_t irl = var_[idx][MATVAR_IRL];
            Double_t eld = var_[idx][MATVAR_ELD];
            Double_t lme = var_[idx][MATVAR_LME];
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


Bool_t MatMgnt::is_load_ = false;
std::list<MatGeoBoxReader*> * MatMgnt::reader_ = nullptr;


MatFld MatMgnt::Get(const SVecD<3>& coo, Double_t log10gb) {
    if (!Load()) return MatFld();
    if (!(Numc::Valid(coo(0)) && Numc::Valid(coo(1)) && Numc::Valid(coo(2)))) return MatFld();
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
    

MatFld MatMgnt::Get(const SVecD<3>& vcoo, const SVecD<3>& wcoo, Double_t log10gb) {
    if (!Load()) return MatFld();
    if (!(Numc::Valid(vcoo(0)) && Numc::Valid(vcoo(1)) && Numc::Valid(vcoo(2)))) return MatFld();
    if (!(Numc::Valid(wcoo(0)) && Numc::Valid(wcoo(1)) && Numc::Valid(wcoo(2)))) return MatFld();
    
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
        MatFld&& mfld = reader->get(vcoo, wcoo, log10gb);
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


MatFld MatMgnt::Get(Double_t stp_len, const PhySt& part) {
    const SVecD<3>&  vcoo = part.c();
    SVecD<3>&&       wcoo = part.c() + stp_len * part.u();
    Double_t log10gb = std::log10(part.gmbta());
    if (!Numc::Valid(log10gb)) log10gb = Numc::ZERO<>;

    return Get(vcoo, wcoo, log10gb);
}


MatFld MatPhy::corr_mfld_;
Bool_t MatPhy::corr_sw_mscat_  = false;
Bool_t MatPhy::corr_sw_eloss_  = false;


MatPhyFld MatPhy::Get(const Double_t stp_len, PhySt& part) {
    if (!part.field()) return MatPhyFld();
    if (part.info().is_chrgless() || part.info().is_massless()) return MatPhyFld();
    if (Numc::EqualToZero(stp_len)) return MatPhyFld();
    if (Numc::EqualToZero(part.mom())) return MatPhyFld();

    const SVecD<3>&  vcoo = part.c();
    SVecD<3>&&       wcoo = part.c() + stp_len * part.u();
    Double_t log10gb = std::log10(part.gmbta());
    if (!Numc::Valid(log10gb)) log10gb = Numc::ZERO<>;

    MatFld&& mfld = MatMgnt::Get(vcoo, wcoo, log10gb);
    
    if (!mfld()) return MatPhyFld();

    Double_t&& mscat_sgm = GetMultipleScattering(mfld, part);
    std::tuple<Double_t, Double_t, Double_t, Double_t>&& ion_eloss = GetIonizationEnergyLoss(mfld, part);
    Double_t eloss_brm_men = GetBremsstrahlungEnergyLoss(mfld, part);

    return MatPhyFld(mfld(), mscat_sgm, std::get<0>(ion_eloss), std::get<1>(ion_eloss), std::get<2>(ion_eloss), std::get<3>(ion_eloss), eloss_brm_men);
}


MatPhyFld MatPhy::Get(const MatFld& mfld, PhySt& part) {
    if (!mfld() || !part.field()) return MatPhyFld();
    if (Numc::EqualToZero(mfld.elen())) return MatPhyFld();
    if (part.info().is_chrgless() || part.info().is_massless()) return MatPhyFld();
    if (Numc::EqualToZero(part.mom())) return MatPhyFld();
    
    Double_t&& mscat_sgm = GetMultipleScattering(mfld, part);
    std::tuple<Double_t, Double_t, Double_t, Double_t>&& ion_eloss = GetIonizationEnergyLoss(mfld, part);
    Double_t eloss_brm_men = GetBremsstrahlungEnergyLoss(mfld, part);
    
    return MatPhyFld(mfld(), mscat_sgm, std::get<0>(ion_eloss), std::get<1>(ion_eloss), std::get<2>(ion_eloss), std::get<3>(ion_eloss), eloss_brm_men);
}
        

Double_t MatPhy::GetMultipleScattering(const MatFld& mfld, PhySt& part) {
    if (!part.arg().mscat()) return Numc::ZERO<>;
    
    Bool_t is_over_lmt = (Numc::Compare(part.bta(), LMT_BTA) > 0);
    Double_t ibta = ((is_over_lmt) ? part.ibta() : LMT_INV_BTA);
    Double_t eta  = ((is_over_lmt) ? part.eta_abs() : (LMT_INV_GMBTA / part.mu()));

    Double_t nrl      = mfld.nrl();
    Double_t sqrt_nrl = std::sqrt(nrl);

    // Highland-Lynch-Dahl formula
    Double_t log_nrl   = ((corr_sw_mscat_) ? std::log(corr_mfld_.nrl()) : std::log(nrl));
    Double_t mscat_crr = (Numc::ONE<> + NRL_CORR_FACT * log_nrl);
    Double_t mscat_fat = RYDBERG_CONST * part.info().chrg_to_atomic_mass() * sqrt_nrl * ((Numc::Valid(mscat_crr) && Numc::Compare(mscat_crr)>0) ? mscat_crr : Numc::ZERO<>);
    
    // Modified Highland-Lynch-Dahl formula
    //Double_t log_nrl   = ((corr_sw_mscat_) ? std::log(corr_mfld_.nrl()) : std::log(nrl));
    //Double_t mscat_crr = std::sqrt(Numc::ONE<> + NRL_CORR_FACT1 * log_nrl + NRL_CORR_FACT2 * log_nrl * log_nrl);
    //Double_t mscat_fat = RYDBERG_CONST * part.info().chrg_to_atomic_mass() * sqrt_nrl * (Numc::Valid(mscat_crr) ? mscat_crr : Numc::ZERO<>);
   
    // Lynch & Dahl formula (PDG 2018)
    //Double_t nrl_crr   = (ibta * ibta * part.chrg() * part.chrg());
    //Double_t log_nrl   = ((corr_sw_mscat_) ? std::log(corr_mfld_.nrl() * nrl_crr) : std::log(nrl * nrl_crr));
    //Double_t mscat_crr = (Numc::ONE<> + NRL_CORR_FACT * log_nrl);
    //Double_t mscat_fat = RYDBERG_CONST * part.info().chrg_to_atomic_mass() * sqrt_nrl * ((Numc::Valid(mscat_crr) && Numc::Compare(mscat_crr)>0) ? mscat_crr : Numc::ZERO<>);

    Double_t mscat_sgm = mscat_fat * (eta * ibta);
    
    // Correction (Tune)
    //mscat_sgm *= 1.01001619337256598e+00;
   
    if (!Numc::Valid(mscat_sgm) || Numc::Compare(mscat_sgm) <= 0) mscat_sgm = Numc::ZERO<>;
    return mscat_sgm;
}


std::tuple<Double_t, Double_t, Double_t, Double_t> MatPhy::GetIonizationEnergyLoss(const MatFld& mfld, PhySt& part) {
    if (!part.arg().eloss()) return std::make_tuple(Numc::ZERO<>, Numc::ZERO<>, Numc::ZERO<>, Numc::ZERO<>);

    Bool_t is_over_lmt   = (Numc::Compare(part.bta(), LMT_BTA) > 0);
    Double_t sqr_gmbta   = ((is_over_lmt) ? (part.gmbta() * part.gmbta()) : LMT_SQR_GMBTA);
    Double_t sqr_bta     = ((is_over_lmt) ? (part.bta() * part.bta()) : LMT_SQR_BTA);
    Double_t gm          = ((is_over_lmt) ? part.gm() : (LMT_GMBTA / LMT_BTA));
    Double_t sqr_chrg    = part.chrg() * part.chrg();
    
    // Trans KE[MeV] to Eta[1]
    Double_t eta_trans = (std::sqrt(Numc::ONE<> + sqr_gmbta) / (sqr_gmbta * part.info().mass() * GEV_TO_MEV)); // [1/MeV]
    if (!Numc::Valid(eta_trans)) eta_trans = Numc::ZERO<>;
    
    // Calculate Matterial Quality
    Double_t log_mean_exc_eng  = mfld.lme(); // log[MeV]
    Double_t elcloud_abundance = mfld.ela(); // [mol cm^-2]
    Double_t density_corr      = mfld.dec(); // [1]
    Double_t Bethe_Bloch       = (Numc::HALF * BETHE_BLOCH_K * elcloud_abundance * sqr_chrg / sqr_bta); // [MeV]
    
    // Calculate Maximum Trans KE
    Double_t mass_rat       = (MASS_EL_IN_GEV / part.mass()); // [1]
    Double_t mass_rel       = (Numc::ONE<> + mass_rat * (Numc::TWO<> * gm + mass_rat)); // [1]
    Double_t elpair_keng    = (Numc::TWO<> * MASS_EL_IN_MEV * sqr_gmbta); // [MeV]
    Double_t max_trans_keng = (elpair_keng / mass_rel); // [MeV]
    
    // Calculate Sigma
    // TODO: (GEANT3 manual W5013) (From GenFit package)
    Double_t elion_sgm = (eta_trans * Bethe_Bloch); // [1]
    
    // Calculate Global Material Quality
    // Calculate Ncl (Number of men eloss to maximum trans eng)
    Double_t global_ela = (corr_sw_eloss_ ? corr_mfld_.ela() : mfld.ela()); // Elcloud Abundance [mol cm^-2]
    Double_t global_BB  = (Numc::HALF * BETHE_BLOCH_K * global_ela * sqr_chrg / sqr_bta); // Bethe Bloch [MeV]

    // Calculate Mean
    Double_t elke_part = (std::log(elpair_keng) - log_mean_exc_eng); // [1]
    Double_t mtke_part = (std::log(max_trans_keng) - log_mean_exc_eng); // [1]
    Double_t elion_men = elion_sgm * (elke_part + mtke_part - Numc::TWO<> * sqr_bta - density_corr); // [1]
    
    // Calculate Mpv
    Double_t bbke_part = (std::log(global_BB) - log_mean_exc_eng); // [1]
    Double_t elion_mpv = elion_sgm * (elke_part + bbke_part + LANDAU_ELOSS_CORR - sqr_bta - density_corr); //[1]
   
    // Correction (Tune)
    //elion_men *= 9.72577696526508273e-01;
    //elion_mpv *= 1.02514814814814826e+00;

    // Els := min(Mean, Mpv)
    Double_t elion_els = std::min(elion_men, elion_mpv);

    // Correction (Tune) 
    //elion_els *= 1.0;

    if (!Numc::Valid(elion_mpv) || Numc::Compare(elion_mpv) <= 0) elion_mpv = Numc::ZERO<>;
    if (!Numc::Valid(elion_sgm) || Numc::Compare(elion_sgm) <= 0) elion_sgm = Numc::ZERO<>;
    if (!Numc::Valid(elion_men) || Numc::Compare(elion_men) <= 0) elion_men = Numc::ZERO<>;
    if (!Numc::Valid(elion_els) || Numc::Compare(elion_els) <= 0) elion_els = Numc::ZERO<>;
   
    return std::make_tuple(elion_mpv, elion_sgm, elion_men, elion_els);
}


// TODO : need to modify to new struction
Double_t MatPhy::GetBremsstrahlungEnergyLoss(const MatFld& mfld, PhySt& part) {
    if (!part.arg().eloss()) return Numc::ZERO<>;
    Double_t chrgmass_sqr = (MASS_EL_IN_GEV * part.info().chrg_to_mass()) * (MASS_EL_IN_GEV * part.info().chrg_to_mass());

    Bool_t   is_over_lmt = (Numc::Compare(part.bta(), LMT_BTA) > 0);
    Double_t sqr_gmbta   = ((is_over_lmt) ? (part.gmbta() * part.gmbta()) : LMT_SQR_GMBTA);
    Double_t eng         = std::sqrt(sqr_gmbta + Numc::ONE<>);
    Double_t ke_part     = (eng - Numc::ONE<>);
    Double_t eta_trans   = (eng / sqr_gmbta);

    Double_t elbrm_men = chrgmass_sqr * (ke_part * eta_trans) * (mfld.nrl() / Numc::LOG_TWO);
    
    if (!Numc::Valid(elbrm_men) || Numc::Compare(elbrm_men) <= 0) elbrm_men = Numc::ZERO<>;
    return elbrm_men;
}


void MatPhy::SetCorrFactor(const MatFld* mfld, PhySt* part, Bool_t sw_mscat, Bool_t sw_eloss) {
    Bool_t sw = ((mfld != nullptr) && (*mfld)()) && (sw_mscat || sw_eloss); 
    if (sw) { corr_sw_mscat_ = sw_mscat; corr_sw_eloss_ = sw_eloss; corr_mfld_ = *mfld; }
    else    { corr_sw_mscat_ = false;    corr_sw_eloss_ = false;    corr_mfld_ = std::move(MatFld()); }
}


} // namespace TrackSys


#endif // __TRACKLibs_MatEnv_C__
