#ifndef __TRACKLibs_PhySt_C__
#define __TRACKLibs_PhySt_C__


namespace TrackSys {


void PhyArg::zero() {
    mat_ = false;
    tme_ = 0.;
    len_ = 0.;
    nrl_ = 0.;
    ela_ = 0.;
    
    sign_ = 1;
    orth_tau_ = SVecD<3>(1, 0, 0);
    orth_rho_ = SVecD<3>(0, 1, 0);

    mscat_uu_ = 0.;
    mscat_ul_ = 0.;
    mscat_ll_ = 0.;

    elion_mpv_ = 0.;
    elion_sgm_ = 0.;
    elbrm_men_ = 0.;
}


void PhyArg::reset(Bool_t sw_mscat, Bool_t sw_eloss) {
    sw_mscat_ = sw_mscat;
    sw_eloss_ = sw_eloss;

    field_ = (sw_mscat || sw_eloss);

    tauu_ = 0.;
    rhou_ = 0.;
    taul_ = 0.;
    rhol_ = 0.;
    elion_ = 0.;
    elbrm_ = 0.;

    zero();
}


void PhySt::reset(const PartType& type) {
    info_ = std::move(PartInfo(type));
    mom_ = 0;
    eng_ = 0;
    ke_  = 0;
    bta_ = 0;
    eta_ = info_.mass();
    irig_ = 0;
    coo_ = std::move(SVecD<3>(0, 0, 0));
    dir_ = std::move(SVecD<3>(0, 0, -1));
    arg_.clear();
    
    path_ = Numc::ZERO<>;
    time_ = Numc::ZERO<>;
}
        

void PhySt::set_state_with_cos(Double_t cx, Double_t cy, Double_t cz, Double_t ux, Double_t uy, Double_t uz) {
    Double_t norm = std::sqrt(ux * ux + uy * uy + uz * uz);
    if (Numc::EqualToZero(norm)) return;
    coo_(0) = cx;
    coo_(1) = cy;
    coo_(2) = cz;
    dir_(0) = ux / norm;
    dir_(1) = uy / norm;
    dir_(2) = uz / norm;
}


void PhySt::set_state_with_tan(Double_t cx, Double_t cy, Double_t cz, Double_t tx, Double_t ty, Double_t uz) {
    Short_t tz = Numc::Compare(uz);
    if (tz == 0) return;
    Double_t norm_uz = static_cast<Double_t>(tz) / std::sqrt(tx * tx + ty * ty + Numc::ONE<>);
    coo_(0) = cx;
    coo_(1) = cy;
    coo_(2) = cz;
    dir_(0) = tx * norm_uz;
    dir_(1) = ty * norm_uz;
    dir_(2) = norm_uz;
}


void PhySt::set_state_with_uxy(Double_t cx, Double_t cy, Double_t cz, Double_t ux, Double_t uy, Short_t signz) {
    Short_t sign = Numc::Compare(signz);
    if (sign == 0) return;
    Double_t uz = (Numc::ONE<> - ux * ux - uy * uy);
    if (Numc::Compare(uz) <= 0) uz = Numc::ZERO<>;
    else                          uz = sign * std::sqrt(uz);
    Double_t norm = std::sqrt(ux * ux + uy * uy + uz * uz);
    if (Numc::EqualToZero(norm)) return;
    coo_(0) = cx;
    coo_(1) = cy;
    coo_(2) = cz;
    dir_(0) = ux / norm;
    dir_(1) = uy / norm;
    dir_(2) = uz / norm;
}
        

void PhySt::set_state(Double_t cx, Double_t cy, Double_t cz, Double_t mx, Double_t my, Double_t mz) {
    Double_t norm = std::sqrt(mx * mx + my * my + mz * mz);
    if (Numc::EqualToZero(norm)) {
        mom_   = Numc::ZERO<>;
        eng_   = info_.mass();
        ke_    = Numc::ZERO<>;
        bta_   = Numc::ZERO<>;
        gmbta_ = Numc::ZERO<>;
        eta_   = Numc::ZERO<>;
        irig_  = Numc::ZERO<>;
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
    Short_t mom_sign = Numc::Compare(mom);
    Short_t eta_sign = Numc::Compare(sign);
    if (mom_sign  < 0) return;
    if (eta_sign == 0) eta_sign = ((Numc::Compare(info_.chrg()) >= 0) ? 1 : -1);

    if (mom_sign == 0) {
        mom_   = Numc::ZERO<>;
        eng_   = info_.mass();
        ke_    = Numc::ZERO<>;
        bta_   = Numc::ZERO<>;
        gmbta_ = Numc::ZERO<>;
        eta_   = Numc::ZERO<>;
        irig_  = Numc::ZERO<>;
    }
    else {
        mom_   = mom;
        eng_   = std::sqrt(mom_ * mom_ + info_.mass() * info_.mass());
        //eng_   = std::sqrt(static_cast<long double>(mom_) * static_cast<long double>(mom_) + static_cast<long double>(info_.mass()) * static_cast<long double>(info_.mass()));
        //eng_   = std::hypot(mom_, info_.mass());
        ke_    = (eng_ - info_.mass());
        bta_   = (info_.is_massless() ? Numc::ONE<> : (Numc::ONE<> / std::sqrt(info_.mass() * info_.mass() / mom_ / mom_ + Numc::ONE<>)));
        gmbta_ = (info_.is_massless() ? mom_ : (mom_ / info_.mass()));
        eta_   = static_cast<Double_t>(eta_sign) / gmbta_;
        irig_  = info_.chrg_to_mass() * eta_;
    }
}


void PhySt::set_eta(Double_t eta) {
    Short_t eta_sign = Numc::Compare(eta);
    if (eta_sign == 0) {
        mom_   = Numc::ZERO<>;
        eng_   = info_.mass();
        ke_    = Numc::ZERO<>;
        bta_   = Numc::ZERO<>;
        gmbta_ = Numc::ZERO<>;
        eta_   = Numc::ZERO<>;
        irig_  = Numc::ZERO<>;
    }
    else {
        eta_   = eta;
        gmbta_ = (Numc::ONE<> / std::fabs(eta_));
        irig_  = info_.chrg_to_mass() * eta_;
        mom_   = (info_.is_massless() ? gmbta_ : (info_.mass() * gmbta_));
        eng_   = (info_.is_massless() ? gmbta_ : info_.mass() * std::sqrt(gmbta_ * gmbta_ + Numc::ONE<>));
        ke_    = (eng_ - info_.mass());
        bta_   = (info_.is_massless() ? Numc::ONE<> : (Numc::ONE<> / std::sqrt(info_.mass() * info_.mass() / mom_ / mom_ + Numc::ONE<>)));
    }
}


void PhySt::set_irig(Double_t irig) {
    Short_t rig_sign = Numc::Compare(irig);
    if (info_.is_chrgless()) return;
    if (info_.is_massless()) return;
    if (rig_sign == 0) {
        mom_   = Numc::ZERO<>;
        eng_   = info_.mass();
        ke_    = Numc::ZERO<>;
        bta_   = Numc::ZERO<>;
        gmbta_ = Numc::ZERO<>;
        eta_   = Numc::ZERO<>;
        irig_  = Numc::ZERO<>;
    }
    else {
        irig_  = irig;
        eta_   = info_.mass_to_chrg() * irig;
        gmbta_ = (Numc::ONE<> / std::fabs(eta_));
        mom_   = (info_.is_massless() ? gmbta_ : (info_.mass() * gmbta_));
        eng_   = (info_.is_massless() ? gmbta_ : info_.mass() * std::sqrt(gmbta_ * gmbta_ + Numc::ONE<>));
        ke_    = (eng_ - info_.mass());
        bta_   = (info_.is_massless() ? Numc::ONE<> : (Numc::ONE<> / std::sqrt(eta_ * eta_ + Numc::ONE<>)));
    }
}


void PhySt::set_rig(Double_t rig) {
    Short_t rig_sign = Numc::Compare(rig);
    Double_t irig = ((rig_sign==0) ? Numc::ZERO<> : (Numc::ONE<> / rig));
    set_irig(irig);
}


void PhySt::print() const {
    std::string printStr;
    printStr += STR("============ %-15s ============\n", info_.name().c_str());
    printStr += STR("Bta %14.8f\n", bta_);
    printStr += STR("Mom %14.8f\n", mom_);
    printStr += STR("Eng %14.8f\n", eng_);
    printStr += STR("KE  %14.8f\n", ke_);
    printStr += STR("Eta %14.8f\n", eta_);
    printStr += STR("Rig %14.8f\n", rig());
    printStr += STR("Coo (%11.6f %11.6f %11.6f)\n", coo_(0), coo_(1), coo_(2));
    printStr += STR("Dir (%11.8f %11.8f %11.8f)\n", dir_(0), dir_(1), dir_(2));
    if (arg_()) {
        printStr += STR("Mscat    Tauu %6.2f  Rhou %6.2f\n", arg_.tauu(),  arg_.rhou());
        printStr += STR("Mscat    Taul %6.2f  Rhol %6.2f\n", arg_.taul(),  arg_.rhol());
        printStr += STR("Eloss    Ion  %6.2f  Brm  %6.2f\n", arg_.elion(), arg_.elbrm());
    }
    printStr += STR("=========================================\n");
    COUT(printStr.c_str());
}
        

void PhySt::symbk(Bool_t is_rndm) {
    if (!arg_()) { arg_.clear(); return; }
    if (is_rndm) arg_.rndm();
    
    if (arg_.mscat()) {
        SVecD<3>&& mscatl = arg_.symbk_mscatl();
        mscatl(2) = Numc::ZERO<>; // set dz to zero
        coo_ = std::move(coo_ + mscatl);
        
        //coo_ = std::move(coo_ + arg_.symbk_mscatl());
        dir_ = std::move(LA::Unit(dir_ + arg_.symbk_mscatu()));
    }
    if (arg_.eloss()) {
        Short_t org_sign = eta_sign();
        set_eta(eta_ * (Numc::ONE<> + arg_.symbk_eloss()));
        Short_t sym_sign = eta_sign();
        if (org_sign != sym_sign) set_eta(Numc::ZERO<>);
    }

    arg_.clear();
}


} // namespace TrackSys


#endif // __TRACKLibs_PhySt_C__
