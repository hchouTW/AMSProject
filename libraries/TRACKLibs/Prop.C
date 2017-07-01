#ifndef __TRACKLibs_Prop_C__
#define __TRACKLibs_Prop_C__

namespace TrackSys {


OrthCoord::OrthCoord(const SVecD<3>& org, const SVecD<3>& seed) : OrthCoord() {
    Double_t org_mag = LA::Mag(org);
    if (MGNumc::EqualToZero(org_mag)) return;
    SVecD<3>&& uorg = org / org_mag;
    
    SVecD<3> tag = seed; 
    Double_t tag_mag = LA::Mag(tag);
    if (!MGNumc::EqualToZero(tag_mag)) tag /= tag_mag;
    else {
        tag = std::move(AXIS_X);
        Double_t dotx = std::fabs(LA::Dot(uorg, tag));
        if (MGNumc::Compare(dotx, MGMath::ONE) == 0) {
            tag = std::move(AXIS_Y);
            Double_t doty = std::fabs(LA::Dot(uorg, tag));
            if (MGNumc::Compare(doty, MGMath::ONE) == 0) {
                tag = std::move(AXIS_Z);
                Double_t dotz = std::fabs(LA::Dot(uorg, tag));
                if (MGNumc::Compare(doty, MGMath::ONE) == 0) return;
            }
        }
    }

    Double_t&&   dot = LA::Dot(uorg, tag);
    SVecD<3>&& cross = LA::Cross(uorg, tag);
    Double_t&&  norm = LA::Mag(cross);

    org_ = std::move(uorg);
    tau_ = std::move((tag - dot * uorg) / norm);
    rho_ = std::move(cross / norm);
}


MotionFunc::MotionFunc(const PhySt& part) {
    Double_t Lambda = PROP_FACT * std::fabs(part.part().chrg_to_mass());
    
    MagFld&& mag = MagMgnt::Get(part.coo());
    SVecD<3>&& crsub = LA::Cross(part.dir(), mag());
    
    zeta_p_ = std::move(part.dir());
    zeta_u_ = std::move((Lambda * part.eta()) * crsub);
    zeta_e_ = 0.0;
}


MotionFunc::MotionFunc(const PhySt& part, const MatArg& arg, const MatPhyFld& mat) {
    Double_t Lambda = PROP_FACT * std::fabs(part.part().chrg_to_mass());
    
    MagFld&& mag = MagMgnt::Get(part.coo());
    SVecD<3>&& crsub = LA::Cross(part.dir(), mag());
    
    zeta_p_ = std::move(part.dir());
    zeta_u_ = std::move((Lambda * part.eta()) * crsub);
    zeta_e_ = 0.0;

    if (arg() && mat()) {
        if (arg.mscat()) {
            OrthCoord orth(part.dir(), mag());
            SVecD<3>&& tau_mscat = (arg.tau() * mat.mult_scat_sgm()) * orth.tau();
            SVecD<3>&& rho_mscat = (arg.rho() * mat.mult_scat_sgm()) * orth.rho();
            zeta_u_ += (tau_mscat + rho_mscat);
        }
        if (arg.eloss()) {
            Double_t ion = (arg.ion() * mat.ion_eloss_sgm() + mat.ion_eloss_mpv());
            Double_t brm = (arg.brm() * mat.brm_eloss_men());
            zeta_e_ += ((ion + brm) * part.eta_sign());
        }
    }
}


TransferFunc::TransferFunc(const PhySt& part) {
    Double_t Lambda = PROP_FACT * std::fabs(part.part().chrg_to_mass());
    
    MagFld&& mag = MagMgnt::Get(part.coo());
    SVecD<3>&& crsub = LA::Cross(part.dir(), mag());
    
    kappa_pu_(0) = MGMath::ONE;
    kappa_pu_(1) = MGMath::ONE;
    kappa_pu_(2) = MGMath::ONE;

    kappa_uu_(0, 1) =  mag.z();
    kappa_uu_(0, 2) = -mag.y();
    kappa_uu_(1, 0) = -mag.z();
    kappa_uu_(1, 2) =  mag.x();
    kappa_uu_(2, 0) =  mag.y();
    kappa_uu_(2, 1) = -mag.x();
    kappa_uu_ = std::move((Lambda * part.eta()) * kappa_uu_);
 
    kappa_ue_ = std::move(Lambda * crsub);
    
    kappa_ei_ = 0.0;
    kappa_eb_ = 0.0;
}


TransferFunc::TransferFunc(const PhySt& part, const MatArg& arg, const MatPhyFld& mat) {
    Double_t Lambda = PROP_FACT * std::fabs(part.part().chrg_to_mass());
    
    MagFld&& mag = MagMgnt::Get(part.coo());
    SVecD<3>&& crsub = LA::Cross(part.dir(), mag());
    
    kappa_pu_(0) = MGMath::ONE;
    kappa_pu_(1) = MGMath::ONE;
    kappa_pu_(2) = MGMath::ONE;

    kappa_uu_(0, 1) =  mag.z();
    kappa_uu_(0, 2) = -mag.y();
    kappa_uu_(1, 0) = -mag.z();
    kappa_uu_(1, 2) =  mag.x();
    kappa_uu_(2, 0) =  mag.y();
    kappa_uu_(2, 1) = -mag.x();
    kappa_uu_ = std::move((Lambda * part.eta()) * kappa_uu_);
 
    kappa_ue_ = std::move(Lambda * crsub);

    kappa_ei_ = 0.0;
    kappa_eb_ = 0.0;
    
    if (arg() && mat()) {
        if (arg.mscat()) {
            OrthCoord orth(part.dir(), mag());
            kappa_ut_ = std::move(mat.mult_scat_sgm() * orth.tau());
            kappa_ur_ = std::move(mat.mult_scat_sgm() * orth.rho());
        }
        if (arg.eloss()) {
            kappa_ei_ = mat.ion_eloss_sgm() * part.eta_sign();
            kappa_eb_ = mat.brm_eloss_men() * part.eta_sign();
        }
    }
}
 

void PhyJb::init(Bool_t is_identity) {
    mat_ = false;
    num_rad_len_ = 0.0; 
    if (is_identity) jb_gg_ = std::move(SMtxId()); 
    else             jb_gg_ = std::move(SMtxD<5, 5>());
    jb_gl_ = std::move(SMtxD<5, 4>()); 
}
 

// Step Length
// ward < 0, backward trace
// ward > 0, forward trace
Double_t PropMgnt::GetPropStep(const PhySt& part, Short_t ward, Bool_t mat) {
    Double_t sign = static_cast<Double_t>(MGNumc::Compare(ward));

    // Current
    Double_t cur_mag = LA::Mag(MagMgnt::Get(part.coo())());
    Double_t curve = std::fabs(PROP_FACT * part.irig() * cur_mag);
    if (MGNumc::Compare(curve, LMTL_CURVE) < 0) curve = LMTL_CURVE;
    Double_t pred_step = TUNE_STEP / curve;
    if (MGNumc::Compare(pred_step, LMTU_STEP) > 0) pred_step = LMTU_STEP;
    if (MGNumc::Compare(pred_step, LMTL_STEP) < 0) pred_step = LMTL_STEP;
    
    // Predict
    SVecD<3>&& pred_coo = part.coo() + (sign * pred_step) * part.dir();
    Double_t pred_mag = LA::Mag(MagMgnt::Get(pred_coo)());
    curve = std::fabs(PROP_FACT * part.irig() * MGMath::HALF * (cur_mag + pred_mag));
    if (MGNumc::Compare(curve, LMTL_CURVE) < 0) curve = LMTL_CURVE;
    pred_step = TUNE_STEP / curve;
    if (MGNumc::Compare(pred_step, LMTU_STEP) > 0) pred_step = LMTU_STEP;
    if (MGNumc::Compare(pred_step, LMTL_STEP) < 0) pred_step = LMTL_STEP;

    Double_t step = sign * pred_step;

    // TODO num_rad_len is faliure
    if (mat) {
        Double_t num_rad_len = MatPhy::GetNumRadLen(step, part);
        //step = std::fabs(step / (MGMath::ONE + (num_rad_len / TUNE_MAT)));
        std::cout << Form("STP %14.8f, NRL %14.8f\n", step, num_rad_len);
    }

    return step; 
}


Double_t PropMgnt::GetStep(const PhySt& part, Double_t resStep, Bool_t mat) {
    Short_t  sign = MGNumc::Compare(resStep);
    Double_t len  = GetPropStep(part, sign, mat);
    Double_t res  = std::fabs(resStep);
    
    Double_t length = MGMath::ZERO;
    if      (res < 1.2 * len) length = res;
    else if (res < 1.7 * len) length = 0.5 * len;
    else                      length = len;
    Double_t step = static_cast<Double_t>(sign) * length;

    std::cout << Form("Z %8.2f step %8.2f\n", part.cz(), step);

    return step;
}


Double_t PropMgnt::GetStepToZ(const PhySt& part, Double_t resStepZ, Bool_t mat) {
    Short_t  signz = MGNumc::Compare(resStepZ);
    Short_t  signs = (MGNumc::Compare(part.dz() * resStepZ) >= 0 ? 1 : -1);
    Double_t lens  = GetPropStep(part, signs, mat);

    MotionFunc mnfunc(part);
    Double_t lenz  = std::fabs((static_cast<Double_t>(signs) * lens) * mnfunc.pz() + MGMath::HALF * (lens * lens) * mnfunc.uz());
    Double_t resz  = std::fabs(resStepZ);
    
    Double_t lengthz = MGMath::ZERO;
    Double_t lengths = MGMath::ZERO;
    if      (resz < 1.2 * lenz) { lengthz = resz;       lengths = lens * (resz / lenz); }
    else if (resz < 1.7 * lenz) { lengthz = 0.5 * lenz; lengths = 0.5 * lens;           }
    else                        { lengthz = lenz;       lengths = lens;                 }
    lengthz *= static_cast<Double_t>(signz);
    if (!MGNumc::Valid(lengths)) lengths = lens;
    lengths *= static_cast<Double_t>(signs);

    // step solver
    Double_t step = MGMath::ZERO;
    if (MGNumc::EqualToZero(mnfunc.uz())) {
        step = lengthz / mnfunc.pz();
        if (!MGNumc::Valid(step))
            step = lengths;
    }
    else {
        Double_t discriminant = (mnfunc.pz() * mnfunc.pz() + MGMath::TWO * mnfunc.uz() * lengthz);
        if (MGNumc::Compare(discriminant) < 0) step = lengths;
        else {
            discriminant = std::sqrt(discriminant);
            Double_t solveA = ((MGMath::NEG * mnfunc.pz() + discriminant) / mnfunc.uz());
            Double_t solveB = ((MGMath::NEG * mnfunc.pz() - discriminant) / mnfunc.uz());
            Short_t  signA = MGNumc::Compare(solveA);
            Short_t  signB = MGNumc::Compare(solveB);
            Bool_t   isSolA = (signs == signA);
            Bool_t   isSolB = (signs == signA);
            if (isSolA && isSolB) step = (signs>0) ? std::min(solveA, solveB) : std::max(solveA, solveB); 
            else if (isSolA)      step = solveA;
            else if (isSolB)      step = solveB;
            else                  step = lengths;
        }
    }
    if (!MGNumc::Valid(step)) step = lengths;
    
    std::cout << Form("Z %8.2f step %8.2f\n", part.cz(), step);

    return step;
}


Bool_t PropMgnt::Prop(const Double_t step, PhySt& part, const MatArg& marg, PhyJb* phyJb) {
    Bool_t withJb = (phyJb != nullptr);
    if (withJb) phyJb->init();

    Long64_t iter     = 1;
    Bool_t   is_succ  = false;
    Double_t int_step = MGMath::ZERO;
    while (iter <= LMTU_ITER && !is_succ) {
        Double_t res_step = step - int_step;
        Double_t cur_step = GetStep(part, res_step, marg());

        Bool_t valid = false;
        PhyJb * curJb = (!withJb) ? nullptr : (new PhyJb());
        switch (method_) {
            case Method::kEuler             : valid = PropWithEuler(cur_step, part, marg, phyJb); break;
            case Method::kEulerHeun         : break;
            case Method::kRungeKuttaNystrom : break;
            default : break;
        }
        if (!valid) break;

        // TODO jacb multiplied
        //if (WithJacb) { phyJb->multiplied(*subPhyJb); delete subPhyJb; subPhyJb = nullptr; }
        if (withJb) {
            delete curJb;
        }

        iter++;
        int_step += cur_step;
        is_succ = (MGNumc::Compare(std::fabs(step - int_step), CONV_STEP) < 0);
    }

    return is_succ;
}
        
Bool_t PropMgnt::PropToZ(const Double_t zcoo, PhySt& part, const MatArg& marg, PhyJb* phyJb) {
    Bool_t withJb = (phyJb != nullptr);
    if (withJb) phyJb->init();

    Long64_t iter     = 1;
    Bool_t   is_succ  = false;
    Double_t int_step = MGMath::ZERO;
    while (iter <= LMTU_ITER && !is_succ) {
        Double_t res_stepz = zcoo - part.cz();
        Double_t cur_step  = GetStepToZ(part, res_stepz, marg());

        Bool_t valid = false;
        PhyJb * curJb = (!withJb) ? nullptr : (new PhyJb());
        switch (method_) {
            case Method::kEuler             : valid = PropWithEuler(cur_step, part, marg, phyJb); break;
            case Method::kEulerHeun         : break;
            case Method::kRungeKuttaNystrom : break;
            default : break;
        }
        if (!valid) break;

        // TODO jacb multiplied
        //if (WithJacb) { phyJb->multiplied(*subPhyJb); delete subPhyJb; subPhyJb = nullptr; }
        if (withJb) {
            delete curJb;
        }

        iter++;
        int_step += cur_step;
        is_succ = (MGNumc::Compare(std::fabs(zcoo - part.cz()), CONV_STEP) < 0);
    }

    return is_succ;
}


Bool_t PropMgnt::PropWithEuler(const Double_t step, PhySt& part, const MatArg& marg, PhyJb* phyJb) {
    Bool_t        withJb = (phyJb != nullptr);
    MatPhyFld&&     mfld = (marg() ? MatPhy::Get(step, part, marg) : MatPhyFld());
    Bool_t           mat = (marg() && mfld());
    Short_t     eta_sign = part.eta_sign();
    
    Double_t half_ss = MGMath::HALF * step * step;

    PhySt st0 = part;
    MotionFunc mn0(st0, marg, mfld);
    
    part.set_state_with_cos(
        st0.cx() + step * mn0.px() + half_ss * mn0.ux(),
        st0.cy() + step * mn0.py() + half_ss * mn0.uy(),
        st0.cz() + step * mn0.pz() + half_ss * mn0.uz(),
        st0.dx() + step * mn0.ux(),
        st0.dy() + step * mn0.uy(),
        st0.dz() + step * mn0.uz()
    );
    if (mat) {
        Double_t eta    = (st0.eta() + step * mn0.e());
        Bool_t   is_mch = (MGNumc::Compare(eta) == eta_sign);
        if (is_mch) part.set_eta(eta);
        else        return false;
    }

    if (withJb) {
        TransferFunc tf0(st0, marg, mfld);

        phyJb->init();
        
        phyJb->gg(JPX, JUX) += step * tf0.pu(X);
        phyJb->gg(JPY, JUY) += step * tf0.pu(Y);
        for (Int_t it = JUX; it <= JUY; ++it) {
            Int_t dm = it - JUX;
            phyJb->gg(JPX, it) += half_ss * tf0.uu(X, dm);
            phyJb->gg(JPY, it) += half_ss * tf0.uu(Y, dm);
            phyJb->gg(JUX, it) += step * tf0.uu(X, dm);
            phyJb->gg(JUY, it) += step * tf0.uu(Y, dm);
            phyJb->gg(JEA, it) += step * tf0.ue(dm);
        }

        if (mat) phyJb->num_rad_len() += mfld.num_rad_len();
        if (mat && marg.mscat()) {
            phyJb->gl(JPX, JTAU) += half_ss * tf0.ut(X);
            phyJb->gl(JPY, JTAU) += half_ss * tf0.ut(Y);
            phyJb->gl(JUX, JTAU) += step * tf0.ut(X);
            phyJb->gl(JUY, JTAU) += step * tf0.ut(Y);
            
            phyJb->gl(JPX, JRHO) += half_ss * tf0.ur(X);
            phyJb->gl(JPY, JRHO) += half_ss * tf0.ur(Y);
            phyJb->gl(JUX, JRHO) += step * tf0.ur(X);
            phyJb->gl(JUY, JRHO) += step * tf0.ur(Y);
        }
        if (mat && marg.eloss()) {
            phyJb->gl(JEA, JION) += step * tf0.ei();
            phyJb->gl(JEA, JBRM) += step * tf0.eb();
        }
    }

    return true;
}


} // namespace TrackSys

#endif // __TRACKLibs_Prop_C__
