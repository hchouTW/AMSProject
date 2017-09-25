#ifndef __TRACKLibs_PhySt_C__
#define __TRACKLibs_PhySt_C__


namespace TrackSys {

        
void PhyArg::rndm_eloss_ion(Double_t kpa, Double_t mos) {
    if (sw_eloss_) {
        if (MGNumc::EqualToZero(kpa) || MGNumc::EqualToZero(mos)) return;
        if (pdf_eloss_ion_.size() == 0) gRandom->SetSeed(0);

        std::string ion_var = "([1]/(x*[0]+[1])/[0]/[0])";
        std::string ion_fmt = STR_FMT("(x<-%f)?0.0:TMath::Power(%s, %s)/TMath::Gamma(%s)*TMath::Exp(-%s*([0]*x+TMath::Exp(-[0]*x)))", std::min(LMT_SGM_, MGMath::HALF*mos), ion_var.c_str(), ion_var.c_str(), ion_var.c_str(), ion_var.c_str());

        Int_t idx_kpa = static_cast<Int_t>(std::rint(kpa / STEP_KPA_));
        Int_t idx_mos = static_cast<Int_t>(std::rint(mos / STEP_MOS_));
        
        TF1 * func = nullptr;
        std::pair<Int_t, Int_t> idx(idx_kpa, idx_mos);
        std::map<std::pair<Int_t, Int_t>, TF1*>::iterator it = pdf_eloss_ion_.find(idx);
        if (it == pdf_eloss_ion_.end()) {
            Double_t ref_kpa = idx_kpa * STEP_KPA_;
            Double_t ref_mos = idx_mos * STEP_MOS_;
            func = new TF1(CSTR_FMT("fElossIon%06d%06d", idx_kpa, idx_mos), ion_fmt.c_str(), -LMT_SGM_, LMT_SGM_ + MGMath::EIGHT * LMT_SGM_ * ref_kpa);
            func->SetParameters(ref_kpa, ref_mos);
            func->SetNpx(NPX_);

            pdf_eloss_ion_[idx] = func; 
        }
        else func = it->second;

        eloss_ion_ = func->GetRandom();
    }
}


void PhySt::reset(const PartType& type, Bool_t sw_mscat, Bool_t sw_eloss) {
    part_ = std::move(PartInfo(type));
    mom_ = 0;
    eng_ = 0;
    bta_ = 0;
    eta_ = part_.mass();
    irig_ = 0;
    coo_ = std::move(SVecD<3>(0, 0, 0));
    dir_ = std::move(SVecD<3>(0, 0, -1));
    arg_.reset(sw_mscat, sw_eloss);
    vst_.reset();
}
        

void PhySt::set_state_with_cos(Double_t cx, Double_t cy, Double_t cz, Double_t ux, Double_t uy, Double_t uz) {
    Double_t norm = std::sqrt(ux * ux + uy * uy + uz * uz);
    if (MGNumc::EqualToZero(norm)) return;
    coo_(0) = cx;
    coo_(1) = cy;
    coo_(2) = cz;
    dir_(0) = ux / norm;
    dir_(1) = uy / norm;
    dir_(2) = uz / norm;
}


void PhySt::set_state_with_tan(Double_t cx, Double_t cy, Double_t cz, Double_t tx, Double_t ty, Double_t uz) {
    Short_t tz = MGNumc::Compare(uz);
    Double_t norm_uz = static_cast<Double_t>(tz) / std::sqrt(tx * tx + ty * ty + MGMath::ONE);
    if (tz == 0) return;
    coo_(0) = cx;
    coo_(1) = cy;
    coo_(2) = cz;
    dir_(0) = tx * norm_uz;
    dir_(1) = ty * norm_uz;
    dir_(2) = norm_uz;
}


void PhySt::set_state_with_uxy(Double_t cx, Double_t cy, Double_t cz, Double_t ux, Double_t uy, Short_t signz) {
    Short_t sign = MGNumc::Compare(signz);
    if (sign == 0) return;
    Double_t uz = (MGMath::ONE - ux * ux - uy * uy);
    if (MGNumc::Compare(uz) <= 0) uz = MGMath::ZERO;
    else                          uz = sign * std::sqrt(uz);
    Double_t norm = std::sqrt(ux * ux + uy * uy + uz * uz);
    if (MGNumc::EqualToZero(norm)) return;
    coo_(0) = cx;
    coo_(1) = cy;
    coo_(2) = cz;
    dir_(0) = ux / norm;
    dir_(1) = uy / norm;
    dir_(2) = uz / norm;
}
        

void PhySt::set_state(Double_t cx, Double_t cy, Double_t cz, Double_t mx, Double_t my, Double_t mz) {
    Double_t norm = std::sqrt(mx * mx + my * my + mz * mz);
    if (MGNumc::EqualToZero(norm)) {
        mom_   = MGMath::ZERO;
        eng_   = part_.mass();
        bta_   = MGMath::ZERO;
        gmbta_ = MGMath::ZERO;
        eta_   = MGMath::ZERO;
        irig_  = MGMath::ZERO;
    }
    else {
        coo_(0) = cx;
        coo_(1) = cy;
        coo_(2) = cz;
        dir_(0) = mx / norm;
        dir_(1) = my / norm;
        dir_(2) = mz / norm;
        set_mom(norm);
    }
}
 

void PhySt::set_mom(Double_t mom, Double_t sign) {
    Short_t mom_sign = MGNumc::Compare(mom);
    Short_t eta_sign = MGNumc::Compare(sign);
    if (mom_sign  < 0) return;
    if (eta_sign == 0) eta_sign = ((MGNumc::Compare(part_.chrg()) >= 0) ? 1 : -1);

    if (mom_sign == 0) {
        mom_   = MGMath::ZERO;
        eng_   = part_.mass();
        bta_   = MGMath::ZERO;
        gmbta_ = MGMath::ZERO;
        eta_   = MGMath::ZERO;
        irig_  = MGMath::ZERO;
    }
    else {
        mom_   = mom;
        eng_   = std::sqrt(mom_ * mom_ + part_.mass() * part_.mass());
        bta_   = (part_.is_massless() ? MGMath::ONE : (MGMath::ONE / std::sqrt(part_.mass() * part_.mass() / mom_ / mom_ + MGMath::ONE)));
        gmbta_ = (part_.is_massless() ? mom_ : (mom_ / part_.mass()));
        eta_   = static_cast<Double_t>(eta_sign) / gmbta_;
        irig_  = part_.chrg_to_mass() * eta_;
    }
}


void PhySt::set_eta(Double_t eta) {
    Short_t eta_sign = MGNumc::Compare(eta);
    if (eta_sign == 0) {
        mom_   = MGMath::ZERO;
        eng_   = part_.mass();
        bta_   = MGMath::ZERO;
        gmbta_ = MGMath::ZERO;
        eta_   = MGMath::ZERO;
        irig_  = MGMath::ZERO;
    }
    else {
        eta_   = eta;
        gmbta_ = (MGMath::ONE / std::fabs(eta_));
        irig_  = part_.chrg_to_mass() * eta_;
        mom_   = (part_.is_massless() ? gmbta_ : (part_.mass() * gmbta_));
        eng_   = std::sqrt(mom_ * mom_ + part_.mass() * part_.mass());
        bta_   = (part_.is_massless() ? MGMath::ONE : (MGMath::ONE / std::sqrt(part_.mass() * part_.mass() / mom_ / mom_ + MGMath::ONE)));
    }
}


void PhySt::set_irig(Double_t irig) {
    Short_t rig_sign = MGNumc::Compare(irig);
    if (part_.is_chrgless()) return;
    if (part_.is_massless()) return;
    if (rig_sign == 0) {
        mom_ = MGMath::ZERO;
        eng_   = part_.mass();
        bta_   = MGMath::ZERO;
        gmbta_ = MGMath::ZERO;
        eta_   = MGMath::ZERO;
        irig_  = MGMath::ZERO;
    }
    else {
        irig_  = irig;
        eta_   = part_.mass_to_chrg() * irig;
        gmbta_ = (MGMath::ONE / std::fabs(eta_));
        mom_   = (part_.is_massless() ? gmbta_ : (part_.mass() * gmbta_));
        eng_   = std::sqrt(mom_ * mom_ + part_.mass() * part_.mass());
        bta_   = (part_.is_massless() ? MGMath::ONE : (MGMath::ONE / std::sqrt(part_.mass() * part_.mass() / mom_ / mom_ + MGMath::ONE)));
    }
}


void PhySt::set_rig(Double_t rig) {
    Short_t rig_sign = MGNumc::Compare(rig);
    Double_t irig = ((rig_sign==0) ? MGMath::ZERO : (MGMath::ONE / rig));
    set_irig(irig);
}


void PhySt::print() const {
    std::string printStr;
    printStr += STR_FMT("============ %-15s ============\n", part_.name().c_str());
    printStr += STR_FMT("Bta %14.8f\n", bta_);
    printStr += STR_FMT("Mom %14.8f\n", mom_);
    printStr += STR_FMT("Eta %14.8f\n", eta_);
    printStr += STR_FMT("Rig %14.8f\n", rig());
    printStr += STR_FMT("Coo (%11.6f %11.6f %11.6f)\n", coo_(0), coo_(1), coo_(2));
    printStr += STR_FMT("Dir (%11.8f %11.8f %11.8f)\n", dir_(0), dir_(1), dir_(2));
    if (arg_()) {
        printStr += STR_FMT("Mscat    Tauu %6.2f  Rhou %6.2f\n", arg_.tauu(), arg_.rhou());
        printStr += STR_FMT("Mscat    Tauc %6.2f  Rhoc %6.2f\n", arg_.tauc(), arg_.rhoc());
        printStr += STR_FMT("Eloss    Ion  %6.2f  Brm  %6.2f\n", arg_.ion(),  arg_.brm());
    }
    printStr += STR_FMT("=========================================\n");
    COUT(printStr);
}
        

void PhySt::symbk(Bool_t is_rndm) {
    if (!arg_()) { zero(); return; }
    if (is_rndm) { arg_.rndm(vst_.eloss_ion_kpa(), vst_.eloss_ion_mos(), vst_.nrl()); }

    if (arg_.mscat()) {
        coo_ = std::move(coo_ + vst_.symbk_mscatc(arg_.tauu(), arg_.rhou(), arg_.tauc(), arg_.rhoc()));
        dir_ = std::move(LA::Unit(dir_ + vst_.symbk_mscatu(arg_.tauu(), arg_.rhou())));
    }
    if (arg_.eloss()) {
        Short_t org_sign = eta_sign();
        set_eta(eta_ * (MGMath::ONE + vst_.symbk_eloss(arg_.ion(), arg_.brm())));
        Short_t sym_sign = eta_sign();
        if (org_sign != sym_sign) set_eta(MGMath::ZERO);
    }
    
    zero();
}


} // namespace TrackSys


#endif // __TRACKLibs_PhySt_C__
