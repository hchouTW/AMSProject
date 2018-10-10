#ifndef __TRACKLibs_PhySt_C__
#define __TRACKLibs_PhySt_C__


#include "Sys.h"
#include "Math.h"
#include "PartInfo.h"
#include "PhySt.h"


namespace TrackSys {

Bool_t PhyArg::opt_mscat_ = true;
Bool_t PhyArg::opt_eloss_ = true;

MultiGaus PhyArg::pdf_mscatu_(
    Robust(), 1.0
);

MultiGaus PhyArg::pdf_mscatl_(
    Robust(), 1.0 
);

MultiGaus PhyArg::pdf_elion_(
    Robust(), 1.0
);


void PhyArg::cal_chi(SVecD<6>& chi) const {
    chi = std::move(SVecD<6>());
    if (!field_) return;
    if (sw_mscat_) {
        chi(0) = (pdf_mscatu_.minimizer(Numc::NEG<> * tauu_).at(0)); 
        chi(1) = (pdf_mscatu_.minimizer(Numc::NEG<> * rhou_).at(0)); 
        chi(2) = (pdf_mscatl_.minimizer(Numc::NEG<> * taul_).at(0)); 
        chi(3) = (pdf_mscatl_.minimizer(Numc::NEG<> * rhol_).at(0));
    }
    if (sw_eloss_) {
        chi(4) = (pdf_elion_.minimizer(Numc::NEG<> * elion_).at(0));
    }
}


void PhyArg::cal_nrm(SVecD<6>& nrm) const {
    nrm = std::move(SVecD<6>());
    if (!field_) return;
    if (sw_mscat_) {
        nrm(0) = (pdf_mscatu_.minimizer(Numc::NEG<> * tauu_).at(1)); 
        nrm(1) = (pdf_mscatu_.minimizer(Numc::NEG<> * rhou_).at(1)); 
        nrm(2) = (pdf_mscatl_.minimizer(Numc::NEG<> * taul_).at(1)); 
        nrm(3) = (pdf_mscatl_.minimizer(Numc::NEG<> * rhol_).at(1));
    }
    if (sw_eloss_) {
        nrm(4) = (pdf_elion_.minimizer(Numc::NEG<> * elion_).at(1));
    }
}


void PhyArg::cal_chi_and_div(SVecD<6>& chi, SVecD<6>& div) const {
    chi = std::move(SVecD<6>());
    div = std::move(SVecD<6>());
    if (!field_) return;
    if (sw_mscat_) {
        std::array<long double, 4>&& minitauu = pdf_mscatu_.minimizer(Numc::NEG<> * tauu_);
        std::array<long double, 4>&& minirhou = pdf_mscatu_.minimizer(Numc::NEG<> * rhou_);
        std::array<long double, 4>&& minitaul = pdf_mscatl_.minimizer(Numc::NEG<> * taul_);
        std::array<long double, 4>&& minirhol = pdf_mscatl_.minimizer(Numc::NEG<> * rhol_);
        div(0) = (Numc::NEG<> * minitauu.at(2));
        div(1) = (Numc::NEG<> * minirhou.at(2));
        div(2) = (Numc::NEG<> * minitaul.at(2));
        div(3) = (Numc::NEG<> * minirhol.at(2));
        chi(0) = (minitauu.at(0)); 
        chi(1) = (minirhou.at(0)); 
        chi(2) = (minitaul.at(0)); 
        chi(3) = (minirhol.at(0)); 
    }
    if (sw_eloss_) {
        std::array<long double, 4>&& minielion = pdf_elion_.minimizer(Numc::NEG<> * elion_);
        div(4) = (Numc::NEG<> * minielion.at(2));
        chi(4) = (minielion.at(0));
    }
}


void PhyArg::cal_nrm_and_div(SVecD<6>& nrm, SVecD<6>& div) const {
    nrm = std::move(SVecD<6>());
    div = std::move(SVecD<6>());
    if (!field_) return;
    if (sw_mscat_) {
        std::array<long double, 4>&& minitauu = pdf_mscatu_.minimizer(Numc::NEG<> * tauu_);
        std::array<long double, 4>&& minirhou = pdf_mscatu_.minimizer(Numc::NEG<> * rhou_);
        std::array<long double, 4>&& minitaul = pdf_mscatl_.minimizer(Numc::NEG<> * taul_);
        std::array<long double, 4>&& minirhol = pdf_mscatl_.minimizer(Numc::NEG<> * rhol_);
        div(0) = (Numc::NEG<> * minitauu.at(2));
        div(1) = (Numc::NEG<> * minirhou.at(2));
        div(2) = (Numc::NEG<> * minitaul.at(2));
        div(3) = (Numc::NEG<> * minirhol.at(2));
        nrm(0) = (minitauu.at(1)); 
        nrm(1) = (minirhou.at(1)); 
        nrm(2) = (minitaul.at(1)); 
        nrm(3) = (minirhol.at(1)); 
    }
    if (sw_eloss_) {
        std::array<long double, 4>&& minielion = pdf_elion_.minimizer(Numc::NEG<> * elion_);
        div(4) = (Numc::NEG<> * minielion.at(2));
        nrm(4) = (minielion.at(1));
    }
}


void PhyArg::cal_chi_and_nrm_and_div(SVecD<6>& chi, SVecD<6>& nrm, SVecD<6>& div) const {
    chi = std::move(SVecD<6>());
    nrm = std::move(SVecD<6>());
    div = std::move(SVecD<6>());
    if (!field_) return;
    if (sw_mscat_) {
        std::array<long double, 4>&& minitauu = pdf_mscatu_.minimizer(Numc::NEG<> * tauu_);
        std::array<long double, 4>&& minirhou = pdf_mscatu_.minimizer(Numc::NEG<> * rhou_);
        std::array<long double, 4>&& minitaul = pdf_mscatl_.minimizer(Numc::NEG<> * taul_);
        std::array<long double, 4>&& minirhol = pdf_mscatl_.minimizer(Numc::NEG<> * rhol_);
        div(0) = (Numc::NEG<> * minitauu.at(2));
        div(1) = (Numc::NEG<> * minirhou.at(2));
        div(2) = (Numc::NEG<> * minitaul.at(2));
        div(3) = (Numc::NEG<> * minirhol.at(2));
        nrm(0) = (minitauu.at(1)); 
        nrm(1) = (minirhou.at(1)); 
        nrm(2) = (minitaul.at(1)); 
        nrm(3) = (minirhol.at(1)); 
        chi(0) = (minitauu.at(0)); 
        chi(1) = (minirhou.at(0)); 
        chi(2) = (minitaul.at(0)); 
        chi(3) = (minirhol.at(0)); 
    }
    if (sw_eloss_) {
        std::array<long double, 4>&& minielion = pdf_elion_.minimizer(Numc::NEG<> * elion_);
        div(4) = (Numc::NEG<> * minielion.at(2));
        nrm(4) = (minielion.at(1));
        chi(4) = (minielion.at(0));
    }
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
    printStr += STR("-----------------------------------------\n");
    printStr += STR("Path %14.8f\n", path_);
    printStr += STR("Time %14.8f\n", time_);
    printStr += STR("=========================================\n");
    COUT(printStr.c_str());
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
    else                        uz = sign * std::sqrt(uz);
    Double_t norm = std::sqrt(ux * ux + uy * uy + uz * uz);
    if (Numc::EqualToZero(norm)) return;
    coo_(0) = cx;
    coo_(1) = cy;
    coo_(2) = cz;
    dir_(0) = ux / norm;
    dir_(1) = uy / norm;
    dir_(2) = uz / norm;
}
        

void PhySt::set_state_with_mom(Double_t cx, Double_t cy, Double_t cz, Double_t mx, Double_t my, Double_t mz) {
    Double_t norm = std::sqrt(mx * mx + my * my + mz * mz);
    if (Numc::EqualToZero(norm)) zero();
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
 

void PhySt::set_mom(Double_t mom, Short_t sign) {
    Short_t mom_sign = Numc::Compare(mom);
    Short_t eta_sign = Numc::Compare(sign);
    if (mom_sign  < 0) return;
    if (eta_sign == 0) eta_sign = ((Numc::Compare(info_.chrg()) >= 0) ? 1 : -1);

    if (mom_sign == 0) zero();
    else {
        mom_   = mom;
        eng_   = std::hypot(mom_, info_.mass());
        ke_    = (eng_ - info_.mass());
        eta_   = static_cast<Double_t>(eta_sign) * PartInfo::ATOMIC_MASS / mom_;
        gmbta_ = (info_.is_massless() ? mom_ : (mom_ / info_.mass()));
        bta_   = (info_.is_massless() ? Numc::ONE<> : (Numc::ONE<> / std::hypot(Numc::ONE<>, Numc::ONE<> / gmbta_)));
        irig_  = info_.chrg_to_atomic_mass() * eta_;
    }
}


void PhySt::set_eta(Double_t eta) {
    Short_t eta_sign = Numc::Compare(eta);
    if (eta_sign == 0) zero();
    else {
        eta_   = eta;
        irig_  = info_.chrg_to_atomic_mass() * eta_;
        gmbta_ = Numc::ONE<> / (info_.is_massless() ? std::fabs(eta_) : std::fabs(info_.mu() * eta_));
        bta_   = (info_.is_massless() ? Numc::ONE<> : (Numc::ONE<> / std::hypot(Numc::ONE<>, Numc::ONE<> / gmbta_)));
        mom_   = (info_.is_massless() ? gmbta_ : (info_.mass() * gmbta_));
        eng_   = (info_.is_massless() ? gmbta_ : std::hypot(info_.mass(), mom_));
        ke_    = (eng_ - info_.mass());
    }
}


void PhySt::set_irig(Double_t irig) {
    Short_t irig_sign = Numc::Compare(irig);
    if (info_.is_chrgless() || info_.is_massless()) { zero(); return; }
    if (irig_sign == 0) zero();
    else {
        irig_  = irig;
        eta_   = irig / info_.chrg_to_atomic_mass();
        gmbta_ = Numc::ONE<> / std::fabs(info_.mu() * eta_);
        bta_   = (Numc::ONE<> / std::hypot(Numc::ONE<>, Numc::ONE<> / gmbta_));
        mom_   = (info_.mass() * gmbta_);
        eng_   = std::hypot(info_.mass(), mom_);
        ke_    = (eng_ - info_.mass());
    }
}


void PhySt::set_rig(Double_t rig) {
    Short_t rig_sign = Numc::Compare(rig);
    Double_t irig = ((rig_sign==0) ? Numc::ZERO<> : (Numc::ONE<> / rig));
    set_irig(irig);
}


void PhySt::symbk(Bool_t is_rndm) {
    if (!(arg_.field() && arg_.mat())) { arg_.clear(); return; }
    if (is_rndm) arg_.rndm();
    
    if (arg_.mscat()) {
        SVecD<3>&& mscatu = arg_.symbk_mscatu();
        SVecD<3>&& mscatl = arg_.symbk_mscatl();
        mscatl(2) = Numc::ZERO<>; // proj to z-plane: set dz to zero
        dir_ = std::move(LA::Unit(dir_ + mscatu));
        coo_ = std::move(coo_ + mscatl);
    }
    if (arg_.eloss()) {
        Short_t org_sign = eta_sign();
        set_eta(eta_ * (Numc::ONE<> + arg_.symbk_elion() + arg_.symbk_elbrm()));
        Short_t sym_sign = eta_sign();
        if (org_sign != sym_sign) set_eta(Numc::ZERO<>);
    }

    arg_.clear();
}


} // namespace TrackSys


#endif // __TRACKLibs_PhySt_C__
