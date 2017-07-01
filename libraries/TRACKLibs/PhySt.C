#ifndef __TRACKLibs_PhySt_C__
#define __TRACKLibs_PhySt_C__

namespace TrackSys {

void PhySt::set_state_with_cos(Double_t cx, Double_t cy, Double_t cz, Double_t dx, Double_t dy, Double_t dz) {
    Double_t norm = std::sqrt(dx * dx + dy * dy + dz * dz);
    if (MGNumc::EqualToZero(norm)) return;
    coo_(0) = cx;
    coo_(1) = cy;
    coo_(2) = cz;
    dir_(0) = dx / norm;
    dir_(1) = dy / norm;
    dir_(2) = dz / norm;
}


void PhySt::set_state_with_tan(Double_t cx, Double_t cy, Double_t cz, Double_t tx, Double_t ty, Double_t dz) {
    Short_t tz = MGNumc::Compare(dz);
    Double_t norm_dz = static_cast<Double_t>(tz) / std::sqrt(tx * tx + ty * ty + MGMath::ONE);
    if (tz == 0) return;
    coo_(0) = cx;
    coo_(1) = cy;
    coo_(2) = cz;
    dir_(0) = tx * norm_dz;
    dir_(1) = ty * norm_dz;
    dir_(2) = norm_dz;
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
        eta_  = MGMath::ZERO;
        irig_  = MGMath::ZERO;
    }
    else {
        mom_   = mom;
        eng_   = std::sqrt(mom_ * mom_ + part_.mass() * part_.mass());
        bta_   = (part_.is_massless() ? MGMath::ONE : (MGMath::ONE / std::sqrt(part_.mass() * part_.mass() / mom_ / mom_ + MGMath::ONE)));
        gmbta_ = (part_.is_massless() ? mom_ : (mom_ / part_.mass()));
        eta_   = static_cast<Double_t>(eta_sign) / gmbta_;
        irig_  = std::fabs(part_.chrg_to_mass()) * eta_;
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
        irig_  = std::fabs(part_.chrg_to_mass()) * eta_;
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
        eta_   = std::fabs(part_.mass_to_chrg()) * irig;
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


void PhySt::print() {
    std::string printStr;
    printStr += STR_FMT("============ %-15s ============\n", part_.name().c_str());
    printStr += STR_FMT("Bta %14.8f\n", bta_);
    printStr += STR_FMT("Mom %14.8f\n", mom_);
    printStr += STR_FMT("Eta %14.8f\n", eta_);
    printStr += STR_FMT("Rig %14.8f\n", rig());
    printStr += STR_FMT("Coo (%11.6f %11.6f %11.6f)\n", coo_(0), coo_(1), coo_(2));
    printStr += STR_FMT("Dir (%11.8f %11.8f %11.8f)\n", dir_(0), dir_(1), dir_(2));
    printStr += STR_FMT("=========================================\n");
    COUT(printStr);
}



}

#endif // __TRACKLibs_PhySt_C__
