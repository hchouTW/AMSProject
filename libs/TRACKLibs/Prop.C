#ifndef __TRACKLibs_Prop_C__
#define __TRACKLibs_Prop_C__

namespace TrackSys {


void OrthCoord::reset(const SVecD<3>& org, const SVecD<3>& seed) {
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


MotionFunc::MotionFunc(PhySt& part) {
    Double_t Lambda = PROP_FACT * part.part().chrg_to_mass();
    
    MagFld&& mag = MagMgnt::Get(part.c());
    SVecD<3>&& crsub = LA::Cross(part.u(), mag());
    
    zeta_c_ = std::move(part.u());
    zeta_u_ = std::move((Lambda * part.eta()) * crsub);
    zeta_e_ = MGMath::ZERO;
}


MotionFunc::MotionFunc(PhySt& part, const MatPhyFld& mphy) {
    Double_t Lambda = PROP_FACT * part.part().chrg_to_mass();
    
    MagFld&& mag = MagMgnt::Get(part.c());
    SVecD<3>&& crsub = LA::Cross(part.u(), mag());
    
    zeta_c_ = std::move(part.u());
    zeta_u_ = std::move((Lambda * part.eta()) * crsub);
    
    if (part.field()) orth_.reset(part.u(), mag());
    
    if (part.arg().eloss() && mphy()) zeta_e_ = mphy.eloss_ion_mpv();
    else                              zeta_e_ = MGMath::ZERO;
}


TransferFunc::TransferFunc(PhySt& part) {
    Double_t Lambda = PROP_FACT * part.part().chrg_to_mass();
    
    MagFld&& mag = MagMgnt::Get(part.c());
    SVecD<3>&& crsub = LA::Cross(part.u(), mag());
    
    kappa_cu_(0) = MGMath::ONE;
    kappa_cu_(1) = MGMath::ONE;
    kappa_cu_(2) = MGMath::ONE;

    kappa_uu_(0, 1) =  mag.z();
    kappa_uu_(0, 2) = -mag.y();
    kappa_uu_(1, 0) = -mag.z();
    kappa_uu_(1, 2) =  mag.x();
    kappa_uu_(2, 0) =  mag.y();
    kappa_uu_(2, 1) = -mag.x();
    kappa_uu_ = std::move((Lambda * part.eta()) * kappa_uu_);
 
    kappa_ue_ = std::move(Lambda * crsub);
    
    kappa_ee_ = 0.0;
    kappa_ei_ = 0.0;
    kappa_eb_ = 0.0;
}


TransferFunc::TransferFunc(PhySt& part, const MatPhyFld& mphy) {
    Double_t Lambda = PROP_FACT * part.part().chrg_to_mass();
    
    MagFld&& mag = MagMgnt::Get(part.c());
    SVecD<3>&& crsub = LA::Cross(part.u(), mag());
    
    kappa_cu_(0) = MGMath::ONE;
    kappa_cu_(1) = MGMath::ONE;
    kappa_cu_(2) = MGMath::ONE;

    kappa_uu_(0, 1) =  mag.z();
    kappa_uu_(0, 2) = -mag.y();
    kappa_uu_(1, 0) = -mag.z();
    kappa_uu_(1, 2) =  mag.x();
    kappa_uu_(2, 0) =  mag.y();
    kappa_uu_(2, 1) = -mag.x();
    kappa_uu_ = std::move((Lambda * part.eta()) * kappa_uu_);
 
    kappa_ue_ = std::move(Lambda * crsub);

    kappa_ee_ = 0.0;
    kappa_ei_ = 0.0;
    kappa_eb_ = 0.0;
    
    if (part.field() && mphy()) {
        if (part.arg().mscat()) {
            OrthCoord orth(part.u(), mag());
            kappa_ut_ = std::move(mphy.mult_scat_sgm() * orth.tau());
            kappa_ur_ = std::move(mphy.mult_scat_sgm() * orth.rho());
        }
        if (part.arg().eloss()) {
            Double_t ion = (part.arg().ion() * mphy.eloss_ion_sgm() + mphy.eloss_ion_mpv());
            Double_t brm = (part.arg().brm() * mphy.eloss_brm_men());
            kappa_ee_ = (ion + brm);
            kappa_ei_ = part.eta() * mphy.eloss_ion_sgm();
            kappa_eb_ = part.eta() * mphy.eloss_brm_men();
        }
    }
}
 

void PhyJb::init(const Type& type) {
    mat_ = false;
    num_rad_len_ = 0.0; 
    if (type == Type::kIdentity) jb_gg_ = std::move(SMtxId()); 
    else                         jb_gg_ = std::move(SMtxD<5, 5>());
    jb_gl_ = std::move(SMtxD<5, 4>()); 
}
 

void PhyJb::set_mat(Bool_t mat, Double_t num_rad_len) {
    num_rad_len_ = (mat) ? num_rad_len : 0.0;
    mat_ = mat;
}


void PhyJb::multiplied(PhyJb& phyJb) {
    if      (mat_ && phyJb.mat()) jb_gl_ = std::move(phyJb.gg() * jb_gl_ + phyJb.gl());
    else if (mat_)                jb_gl_ = std::move(phyJb.gg() * jb_gl_);
    else if (phyJb.mat())         jb_gl_ = std::move(phyJb.gl());

    jb_gg_ = std::move(phyJb.gg() * jb_gg_);

    mat_ = (mat_ || phyJb.mat());
    if (mat_) num_rad_len_ = num_rad_len_ + phyJb.num_rad_len();
}
        

TransferPhyJb::TransferPhyJb(Bool_t mat, const TransferFunc& tf, PhyJb& jb) {
    uu_(X, X) = tf.uu(X, X) * jb.gg(JUX, JUX) + tf.uu(X, Y) * jb.gg(JUY, JUX) + tf.ue(X) * jb.gg(JEA, JUX);
    uu_(X, Y) = tf.uu(X, X) * jb.gg(JUX, JUY) + tf.uu(X, Y) * jb.gg(JUY, JUY) + tf.ue(X) * jb.gg(JEA, JUY);
    
    uu_(Y, X) = tf.uu(Y, X) * jb.gg(JUX, JUX) + tf.uu(Y, Y) * jb.gg(JUY, JUX) + tf.ue(Y) * jb.gg(JEA, JUX);
    uu_(Y, Y) = tf.uu(Y, X) * jb.gg(JUX, JUY) + tf.uu(Y, Y) * jb.gg(JUY, JUY) + tf.ue(Y) * jb.gg(JEA, JUY);

    ue_(X) = tf.uu(X, X) * jb.gg(JUX, JEA) + tf.uu(X, Y) * jb.gg(JUY, JEA) + tf.ue(X) * jb.gg(JEA, JEA);
    ue_(Y) = tf.uu(Y, X) * jb.gg(JUX, JEA) + tf.uu(Y, Y) * jb.gg(JUY, JEA) + tf.ue(Y) * jb.gg(JEA, JEA);

    mat_ = mat;
    ee_ = 0.;
    ei_ = 0.;
    eb_ = 0.;
    if (mat) {
        ut_(X) = tf.ut(X) + tf.uu(X, X) * jb.gl(JUX, JTAU) + tf.uu(X, Y) * jb.gl(JUY, JTAU) + tf.ue(X) * jb.gl(JEA, JTAU);
        ut_(Y) = tf.ut(Y) + tf.uu(Y, X) * jb.gl(JUX, JTAU) + tf.uu(Y, Y) * jb.gl(JUY, JTAU) + tf.ue(Y) * jb.gl(JEA, JTAU);
    
        ur_(X) = tf.ur(X) + tf.uu(X, X) * jb.gl(JUX, JRHO) + tf.uu(X, Y) * jb.gl(JUY, JRHO) + tf.ue(X) * jb.gl(JEA, JRHO);
        ur_(Y) = tf.ur(Y) + tf.uu(Y, X) * jb.gl(JUX, JRHO) + tf.uu(Y, Y) * jb.gl(JUY, JRHO) + tf.ue(Y) * jb.gl(JEA, JRHO);

        ee_ = tf.ee() * jb.gg(JEA, JEA);
        ei_ = tf.ee() * jb.gl(JEA, JION) + tf.ei();
        eb_ = tf.ee() * jb.gl(JEA, JBRM) + tf.eb();
    }
}


void PropPhyCal::init(Double_t eta) {
    sw_mscat_      = false;
    sw_eloss_      = false;
    eta_           = eta;
    len_           = MGMath::ZERO;
    nrl_           = MGMath::ZERO;
    tau_           = std::move(SVecD<3>(1.0,  0.0, 0.0));
    rho_           = std::move(SVecD<3>(0.0, -1.0, 0.0));
    mscatu_        = MGMath::ZERO;
    mscatcu_       = MGMath::ZERO;
    mscatcl_       = MGMath::ZERO;
    eloss_ion_kpa_ = MGMath::ZERO;
    eloss_ion_mpv_ = MGMath::ZERO;
    eloss_ion_sgm_ = MGMath::ZERO;
    eloss_ion_     = MGMath::ZERO;
    eloss_brm_     = MGMath::ZERO;
    vec_vac_.clear();
    vec_len_.clear();
    vec_eft_.clear();
    vec_invloc_.clear();
    vec_invlocsqr_.clear();
    vec_mscat_.clear();
}


void PropPhyCal::normalized(PhySt& part, MatFld& mfld) {
    if (sw_mscat_ && !MGNumc::EqualToZero(mscatu_)) {
        tau_.Unit();
        rho_.Unit();
        mscatu_  = std::sqrt(mscatu_);
        mscatcl_ = (mscatcl_ / mscatu_) * MGMath::INV_SQRT_TWELVE;
        mscatcu_ = MGMath::ZERO;

        Double_t real_len = 0.0;
        Double_t efft_len = 0.0;
        for (Int_t it = vec_vac_.size()-1; it >= 0; --it) {
            if (it != vec_vac_.size()-1) {
                real_len += vec_len_.at(it+1);
                efft_len += vec_eft_.at(it+1);
            }
            if (vec_vac_.at(it)) continue;
            Double_t len = vec_len_.at(it);
            Double_t eft = vec_eft_.at(it);
            Double_t est_efft = (efft_len + eft) / (real_len + len);
            Double_t est_eftl = (MGMath::ONE - (MGMath::ONE - MGMath::INV_SQRT_THREE) * est_efft);
            Double_t est_len2 = (real_len*real_len + MGMath::TWO * real_len * vec_invloc_.at(it) + vec_invlocsqr_.at(it));
            Double_t est_wgt = (est_eftl * est_eftl) * est_len2 * vec_mscat_.at(it); // version 1
            //Double_t est_wgt = est_eftl * std::sqrt(est_len2) * vec_mscat_.at(it); // version 2
            mscatcu_ += est_wgt;
        }
        mscatcu_ = std::sqrt(mscatcu_); // version 1
        //mscatcu_ = (mscatcu_ / mscatu_); // version 2

        vec_vac_.clear();
        vec_len_.clear();
        vec_eft_.clear();
        vec_invloc_.clear();
        vec_invlocsqr_.clear();
        vec_mscat_.clear();
    }
    //if (sw_eloss_ && !MGNumc::EqualToZero(eloss_ion_)) {
    if (sw_eloss_) {
        PhySt mid = part;
        //Double_t efft_eta = (MGNumc::Compare(eta_)) * std::sqrt(part.eta() * eta_);
        //mid.set_eta(efft_eta);
        MatPhyFld&& mpfld = MatPhy::Get(mfld, mid);


        //eloss_ion_mpv_ = mpfld.eloss_ion_mpv();
        eloss_ion_mpv_ = std::fabs(part.eta() - eta_);
        eloss_ion_kpa_ = mpfld.eloss_ion_kpa();
        Double_t mos = mpfld.eloss_ion_mpv() / mpfld.eloss_ion_sgm();
        eloss_ion_sgm_ = eloss_ion_mpv_ / mos;
        eloss_ion_ = eloss_ion_kpa_ * eloss_ion_sgm_;

        //std::cout << Form("MPV %14.8f %14.8f  RAT %14.8f\n", eloss_ion_mpv_, mpfld.eloss_ion_mpv(), eloss_ion_mpv_/mpfld.eloss_ion_mpv());

        //std::cout << Form("ETA %14.8f %14.8f KPA %14.8f MPV %14.8f SGM %14.8f MPV/SGM %14.8f SIGX %14.8f\n", eta_, eta, eloss_ion_kpa_, eloss_ion_mpv_, eloss_ion_sgm_, eloss_ion_mpv_/eloss_ion_sgm_, eloss_ion_);
    }
    if (sw_eloss_ && !MGNumc::EqualToZero(eloss_brm_)) {
        MatPhyFld&& mpfld = MatPhy::Get(mfld, part);
        eloss_brm_ = mpfld.eloss_brm_men();
    }
}


void PropPhyCal::push(const MatPhyFld& mpfld, const SVecD<3>& tau, const SVecD<3>& rho) {
    len_ += mpfld.len();
    nrl_ += mpfld.num_rad_len();
    if (sw_mscat_) {
        Double_t len = mpfld.len();
        Double_t mscat = mpfld.mult_scat_sgm(); 
        Bool_t vac = MGNumc::EqualToZero(mscat);
        Double_t wgt = mscat * mscat;
        vec_vac_.push_back(vac);
        vec_len_.push_back(len);
        if (vac) {
            vec_eft_.push_back(0.0);
            vec_invloc_.push_back(0.0);
            vec_invlocsqr_.push_back(0.0);
            vec_mscat_.push_back(0.0);
        }
        else {
            Double_t invloc = (MGMath::ONE - mpfld.loc()) * len;
            Double_t invlocsqr = (MGMath::ONE - MGMath::TWO * mpfld.loc() + mpfld.locsqr()) * (len * len);
            vec_eft_.push_back(len * mpfld.efft());
            vec_invloc_.push_back(invloc);
            vec_invlocsqr_.push_back(invlocsqr);
            vec_mscat_.push_back(wgt);
        }
        if (vac) return;
        mscatu_  += wgt;
        tau_     += wgt * tau;
        rho_     += wgt * rho;
        mscatcl_ += wgt * (len * mpfld.efft());

        //std::cout << Form("LEN %14.8f MSCAT %14.8f LOC %14.8f LOCSQR %14.8f\n", len, mscat, mpfld.loc(), mpfld.locsqr());
    }
}


void PropPhyCal::set_virtualPhySt(PhySt& part) const {
    part.vst().set_len(len_);
    part.vst().set_nrl(nrl_);
    part.vst().set_orth(tau_, rho_);
    part.vst().set_mscatu(mscatu_);
    part.vst().set_mscatc(mscatcu_, mscatcl_);
    part.vst().set_eloss_ion(eloss_ion_, eloss_ion_kpa_, eloss_ion_mpv_/eloss_ion_sgm_);
    part.vst().set_eloss_brm(eloss_brm_);
}


#ifdef __HAS_AMS_OFFICE_LIBS__
#include <TrProp.h>

Bool_t PropMgnt::PropToZ_AMSLibs(const Double_t zcoo, PhySt& part) {
    Short_t sign = MGNumc::Compare(part.dz());
    if (sign == 0) return false;
    if (MGNumc::Compare(std::fabs(zcoo - part.cz()), CONV_STEP) < 0) return true;

    AMSPoint pos(part.cx(), part.cy(), part.cz());
    AMSDir   dir(part.ux(), part.uy(), part.uz());
    Double_t rig = static_cast<Double_t>(-sign) * part.rig();

    TrProp trProp(pos, dir, rig);
    Double_t len = trProp.Propagate(zcoo);
    if (len < 0) return false;
   
    part.set_state_with_cos(
        trProp.GetP0x(),
        trProp.GetP0y(),
        trProp.GetP0z(),
        trProp.GetD0x() * sign,
        trProp.GetD0y() * sign,
        trProp.GetD0z() * sign
    );
    return true;
}
#endif
 

// Step Length
// ward < 0, backward trace
// ward > 0, forward trace
Double_t PropMgnt::GetPropStep(PhySt& part, Short_t ward) {
    Double_t sign = static_cast<Double_t>(MGNumc::Compare(ward));

    // Current
    Double_t cur_mag = LA::Mag(MagMgnt::Get(part.c())());
    Double_t curve = std::fabs(PROP_FACT * part.irig() * cur_mag);
    if (MGNumc::Compare(curve, LMTL_CURVE) < 0) curve = LMTL_CURVE;
    Double_t pred_step = TUNE_STEP / curve;
    if (MGNumc::Compare(pred_step, LMTU_STEP) > 0) pred_step = LMTU_STEP;
    if (MGNumc::Compare(pred_step, LMTL_STEP) < 0) pred_step = LMTL_STEP;

    // Predict
    SVecD<3>&& pred_coo = part.c() + (sign * pred_step) * part.u();
    Double_t pred_mag = LA::Mag(MagMgnt::Get(pred_coo)());
    curve = std::fabs(PROP_FACT * part.irig() * MGMath::HALF * (cur_mag + pred_mag));
    if (MGNumc::Compare(curve, LMTL_CURVE) < 0) curve = LMTL_CURVE;
    pred_step = TUNE_STEP / curve;
    if (MGNumc::Compare(pred_step, LMTU_STEP) > 0) pred_step = LMTU_STEP;
    if (MGNumc::Compare(pred_step, LMTL_STEP) < 0) pred_step = LMTL_STEP;

    if (MGNumc::EqualToZero(cur_mag) && MGNumc::EqualToZero(pred_mag)) pred_step = PROP_STEP;

    Double_t step = pred_step;
    if (part.field()) {
        Double_t beta_corr = std::sqrt(part.bta());
        Double_t num_rad_len = MatPhy::GetNumRadLen((sign * step), part, false);
        step = std::fabs(step / (MGMath::ONE + (num_rad_len / TUNE_MAT))) * beta_corr;
    }

    return step; 
}


Double_t PropMgnt::GetStep(PhySt& part, Double_t resStep) {
    Short_t  sign = MGNumc::Compare(resStep);
    Double_t len  = GetPropStep(part, sign);
    Double_t res  = std::fabs(resStep);
    
    Double_t length = MGMath::ZERO;
    if      (res < 1.2 * len) length = res;
    else if (res < 1.7 * len) length = 0.5 * len;
    else                      length = len;
    Double_t step = static_cast<Double_t>(sign) * length;

    return step;
}


Double_t PropMgnt::GetStepToZ(PhySt& part, Double_t resStepZ) {
    Short_t  signz = MGNumc::Compare(resStepZ);
    Short_t  signs = (MGNumc::Compare(part.uz() * resStepZ) >= 0 ? 1 : -1);
    Double_t lens  = GetPropStep(part, signs);

    MotionFunc mnfunc(part);
    Double_t lenz  = std::fabs((static_cast<Double_t>(signs) * lens) * mnfunc.cz() + MGMath::HALF * (lens * lens) * mnfunc.uz());
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
        step = lengthz / mnfunc.cz();
        if (!MGNumc::Valid(step))
            step = lengths;
    }
    else {
        Double_t discriminant = (mnfunc.cz() * mnfunc.cz() + MGMath::TWO * mnfunc.uz() * lengthz);
        if (MGNumc::Compare(discriminant) < 0) step = lengths;
        else {
            discriminant = std::sqrt(discriminant);
            Double_t solveA = ((MGMath::NEG * mnfunc.cz() + discriminant) / mnfunc.uz());
            Double_t solveB = ((MGMath::NEG * mnfunc.cz() - discriminant) / mnfunc.uz());
            Short_t  signA = MGNumc::Compare(solveA);
            Short_t  signB = MGNumc::Compare(solveB);
            Bool_t   isSolA = (signs == signA);
            Bool_t   isSolB = (signs == signB);
            if (isSolA && isSolB) step = (signs>0) ? std::min(solveA, solveB) : std::max(solveA, solveB); 
            else if (isSolA)      step = solveA;
            else if (isSolB)      step = solveB;
            else                  step = lengths;
        }
    }
    if (!MGNumc::Valid(step)) step = lengths;
    
    return step;
}
        

Bool_t PropMgnt::Prop(const Double_t step, PhySt& part, PhyJb* phyJb, MatFld* mfld) {
    Bool_t withJb = (phyJb != nullptr);
    Bool_t withMf = (mfld != nullptr);
    if (withJb) phyJb->init(PhyJb::Type::kIdentity);
    std::list<MatFld> mflds;

    Long64_t iter     = 1;
    Bool_t   is_succ  = false;
    Double_t int_step = MGMath::ZERO;
    while (iter <= LMTU_ITER && !is_succ) {
        Double_t res_step = step - int_step;
        Double_t cur_step = GetStep(part, res_step);

        Bool_t valid = false;
        PhyJb * curJb = ((withJb) ? (new PhyJb()) : nullptr);
        MatFld&& curMfld = ((part.field() || withMf) ? MatMgnt::Get(cur_step, part) : MatFld(std::fabs(cur_step)));
        switch (method_) {
            case Method::kRungeKuttaNystrom : valid = PropWithRungeKuttaNystrom(cur_step, part, curMfld, curJb); break;
            case Method::kEulerHeun         : valid = PropWithEulerHeun(cur_step, part, curMfld, curJb); break;
            case Method::kEuler             : valid = PropWithEuler(cur_step, part, curMfld, curJb); break;
            default : break;
        }
        if (!valid) { delete curJb; break; }

        if (withJb) {
            phyJb->multiplied(*curJb);
            delete curJb;
        }

        if (withMf) mflds.push_back(curMfld);

        iter++;
        int_step += cur_step;
        is_succ = (MGNumc::Compare(std::fabs(step - int_step), CONV_STEP) < 0);
    }
    if (withMf) (*mfld) = std::move(MatFld::Merge(mflds));

    return is_succ;
}


Bool_t PropMgnt::PropToZ(const Double_t zcoo, PhySt& part, PhyJb* phyJb, MatFld* mfld) {
    Bool_t withJb = (phyJb != nullptr);
    Bool_t withMf = (mfld != nullptr);
    if (withJb) phyJb->init(PhyJb::Type::kIdentity);
    std::list<MatFld> mflds;

    Long64_t iter     = 1;
    Bool_t   is_succ  = false;
    Double_t int_step = MGMath::ZERO;
    while (iter <= LMTU_ITER && !is_succ) {
        Double_t res_stepz = zcoo - part.cz();
        Double_t cur_step  = GetStepToZ(part, res_stepz);

        Bool_t valid = false;
        PhyJb * curJb = ((withJb) ? (new PhyJb()) : nullptr);
        MatFld&& curMfld = ((part.field() || withMf) ? MatMgnt::Get(cur_step, part) : MatFld(std::fabs(cur_step)));
        switch (method_) {
            case Method::kRungeKuttaNystrom : valid = PropWithRungeKuttaNystrom(cur_step, part, curMfld, curJb); break;
            case Method::kEulerHeun         : valid = PropWithEulerHeun(cur_step, part, curMfld, curJb); break;
            case Method::kEuler             : valid = PropWithEuler(cur_step, part, curMfld, curJb); break;
            default : break;
        }
        if (!valid) { delete curJb; break; }

        if (withJb) {
            phyJb->multiplied(*curJb);
            delete curJb;
        }
        
        if (withMf) mflds.push_back(curMfld);

        iter++;
        int_step += cur_step;
        is_succ = (MGNumc::Compare(std::fabs(zcoo - part.cz()), CONV_STEP) < 0);
    }
    if (withMf) (*mfld) = std::move(MatFld::Merge(mflds));

    return is_succ;
}
        

Bool_t PropMgnt::PropWithMC(const Double_t step, PhySt& part) {
    Long64_t iter     = 1;
    Bool_t   is_succ  = false;
    Double_t int_step = MGMath::ZERO;
    
    while (iter <= LMTU_ITER && !is_succ) {
        Double_t res_step = step - int_step;
        Double_t cur_step = GetStep(part, res_step);
    
        Bool_t valid = false;
        MatFld&& mfld = (part.field() ? MatMgnt::Get(cur_step, part) : MatFld(std::fabs(cur_step)));
        part.arg().rndm(mfld.num_rad_len());
        
        switch (method_) {
            case Method::kRungeKuttaNystrom : valid = PropWithRungeKuttaNystrom(cur_step, part, mfld); break;
            case Method::kEulerHeun         : valid = PropWithEulerHeun(cur_step, part, mfld); break;
            case Method::kEuler             : valid = PropWithEuler(cur_step, part, mfld); break;
            default : break;
        }
        if (!valid) break;

        iter++;
        int_step += cur_step;
        is_succ = (MGNumc::Compare(std::fabs(step - int_step), CONV_STEP) < 0);
    }
    
    return is_succ;
}


/*
Bool_t PropMgnt::PropToZWithMC(const Double_t zcoo, PhySt& part) {
    Long64_t iter     = 1;
    Bool_t   is_succ  = false;
    Double_t int_step = MGMath::ZERO;
    
    while (iter <= LMTU_ITER && !is_succ) {
        Double_t res_stepz = zcoo - part.cz();
        Double_t cur_step  = GetStepToZ(part, res_stepz);
        
        Bool_t valid = false;
        MatFld&& mfld = (part.field() ? MatMgnt::Get(cur_step, part) : MatFld(std::fabs(cur_step)));
        part.arg().rndm(mfld.num_rad_len());
        
        switch (method_) {
            case Method::kRungeKuttaNystrom : valid = PropWithRungeKuttaNystrom(cur_step, part, mfld); break;
            case Method::kEulerHeun         : valid = PropWithEulerHeun(cur_step, part, mfld); break;
            case Method::kEuler             : valid = PropWithEuler(cur_step, part, mfld); break;
            default : break;
        }
        if (!valid) break;

        iter++;
        int_step += cur_step;
        is_succ = (MGNumc::Compare(std::fabs(zcoo - part.cz()), CONV_STEP) < 0);
    }

    return is_succ;
}
*/


// testcode
Bool_t PropMgnt::PropToZWithMC(const Double_t zcoo, PhySt& part) {
    Long64_t iter     = 1;
    Bool_t   is_succ  = false;
    Double_t int_step = MGMath::ZERO;
    std::list<MatFld> mflds;

    PropPhyCal ppcal(part.eta(), part.arg().mscat(), part.arg().eloss());
    while (iter <= LMTU_ITER && !is_succ) {
        Double_t res_stepz = zcoo - part.cz();
        Double_t cur_step  = GetStepToZ(part, res_stepz);
        
        Bool_t valid = false;
        MatFld&& mfld = (part.field() ? MatMgnt::Get(cur_step, part) : MatFld(std::fabs(cur_step)));
        MatPhyFld&& mpfld = MatPhy::Get(mfld, part);
        part.arg().rndm(mpfld.eloss_ion_kpa(), mpfld.eloss_ion_mpv()/mpfld.eloss_ion_sgm(), mfld.num_rad_len());

        valid = PropWithEuler2(cur_step, part, mfld, ppcal);
        if (!valid) break;

        iter++;
        int_step += cur_step;
        is_succ = (MGNumc::Compare(std::fabs(zcoo - part.cz()), CONV_STEP) < 0);
        
        mflds.push_back(mfld);
    }
    MatFld&& fld = std::move(MatFld::Merge(mflds));
    ppcal.normalized(part, fld);

    ppcal.set_virtualPhySt(part);
    part.symbk(true);

    return is_succ;
}





//////////////////////////////////////////////////
// testcode
Bool_t PropMgnt::PropWithEuler2(const Double_t step, PhySt& part, const MatFld& mfld, PropPhyCal& ppcal, PhyJb* phyJb) {
    Bool_t        withJb = (phyJb != nullptr);
    Bool_t           mat = (part.field() && mfld());
    Short_t     eta_sign = part.eta_sign();
    
    Double_t sstep = step * step;
    Double_t ss1o2 = MGMath::HALF * sstep;

    PhySt st0 = part;
    MatPhyFld&& mp0 = MatPhy::Get(mfld, st0);
    MotionFunc mn0(st0, mp0);

    part.set_state_with_cos(
        st0.cx() + step * mn0.cx() + ss1o2 * mn0.ux(),
        st0.cy() + step * mn0.cy() + ss1o2 * mn0.uy(),
        st0.cz() + step * mn0.cz() + ss1o2 * mn0.uz(),
        st0.ux() + step * mn0.ux(),
        st0.uy() + step * mn0.uy(),
        st0.uz() + step * mn0.uz()
    );
    if (mat) {
        Double_t eta    = st0.eta() + mn0.e();
        Bool_t   is_mch = (MGNumc::Compare(eta) == eta_sign);
        if (is_mch) part.set_eta(eta);
        else        return false;
    }
    
    ppcal.push(mp0, mn0.orth().tau(), mn0.orth().rho());
    Double_t sgm = st0.eta() * mp0.eloss_ion_kpa() * mp0.eloss_ion_sgm(); 

    return true;
}
/////////////////////////////////////////





Bool_t PropMgnt::PropWithEuler(const Double_t step, PhySt& part, const MatFld& mfld, PhyJb* phyJb) {
    Bool_t        withJb = (phyJb != nullptr);
    Bool_t           mat = (part.field() && mfld());
    Short_t     eta_sign = part.eta_sign();
    
    Double_t sstep = step * step;
    Double_t ss1o2 = MGMath::HALF * sstep;

    PhySt st0 = part;
    MatPhyFld&& mp0 = MatPhy::Get(mfld, st0);
    MotionFunc mn0(st0, mp0);
    
    part.set_state_with_cos(
        st0.cx() + step * mn0.cx() + ss1o2 * mn0.ux(),
        st0.cy() + step * mn0.cy() + ss1o2 * mn0.uy(),
        st0.cz() + step * mn0.cz() + ss1o2 * mn0.uz(),
        st0.ux() + step * mn0.ux(),
        st0.uy() + step * mn0.uy(),
        st0.uz() + step * mn0.uz()
    );
    if (mat) {
        Double_t eta    = st0.eta() + step * mn0.e();
        Bool_t   is_mch = (MGNumc::Compare(eta) == eta_sign);
        if (is_mch) part.set_eta(eta);
        else        return false;
    }
  
    /*
    if (withJb) {
        TransferFunc tf0(st0, marg, mp0);

        phyJb->init(PhyJb::Type::kIdentity);
        phyJb->set_mat(mat, mfld.num_rad_len());
        
        phyJb->gg(JPX, JUX) += step * tf0.cu(X);
        phyJb->gg(JPY, JUY) += step * tf0.cu(Y);
        
        phyJb->gg(JPX, JUX) += ss1o2 * tf0.uu(X, X);
        phyJb->gg(JPX, JUY) += ss1o2 * tf0.uu(X, Y);
        phyJb->gg(JPX, JEA) += ss1o2 * tf0.ue(X);
        
        phyJb->gg(JPY, JUX) += ss1o2 * tf0.uu(Y, X);
        phyJb->gg(JPY, JUY) += ss1o2 * tf0.uu(Y, Y);
        phyJb->gg(JPY, JEA) += ss1o2 * tf0.ue(Y);
        
        phyJb->gg(JUX, JUX) += step * tf0.uu(X, X);
        phyJb->gg(JUX, JUY) += step * tf0.uu(X, Y);
        phyJb->gg(JUX, JEA) += step * tf0.ue(X);
        
        phyJb->gg(JUY, JUX) += step * tf0.uu(Y, X);
        phyJb->gg(JUY, JUY) += step * tf0.uu(Y, Y);
        phyJb->gg(JUY, JEA) += step * tf0.ue(Y);

        if (mat && marg.mscat()) {
            phyJb->gl(JPX, JTAU) += ss1o2 * tf0.ut(X);
            phyJb->gl(JPX, JRHO) += ss1o2 * tf0.ur(X);
            
            phyJb->gl(JPY, JTAU) += ss1o2 * tf0.ut(Y);
            phyJb->gl(JPY, JRHO) += ss1o2 * tf0.ur(Y);
            
            phyJb->gl(JUX, JTAU) += step * tf0.ut(X);
            phyJb->gl(JUX, JRHO) += step * tf0.ur(X);
            
            phyJb->gl(JUY, JTAU) += step * tf0.ut(Y);
            phyJb->gl(JUY, JRHO) += step * tf0.ur(Y);
        }
        if (mat && marg.eloss()) {
            phyJb->gg(JEA, JEA) += step * tf0.ee();
            
            phyJb->gl(JEA, JION) += step * tf0.ei();
            phyJb->gl(JEA, JBRM) += step * tf0.eb();
        }
    }
    */

    return true;
}


Bool_t PropMgnt::PropWithEulerHeun(const Double_t step, PhySt& part, const MatFld& mfld, PhyJb* phyJb) {
    Bool_t        withJb = (phyJb != nullptr);
    Bool_t           mat = (part.field() && mfld());
    Short_t     eta_sign = part.eta_sign();
   
    Double_t s1o2  = MGMath::HALF * step;
    Double_t sstep = step * step;
    Double_t ss1o2 = MGMath::HALF * sstep;
    Double_t ss1o6 = MGMath::ONE_TO_SIX * sstep;

    PhySt st0 = part;
    MatPhyFld&& mp0 = MatPhy::Get(mfld, st0);
    MotionFunc mn0(st0, mp0);
    
    PhySt st1 = part;
    st1.set_state_with_cos(
        st0.cx() + step * mn0.cx() + ss1o2 * mn0.ux(),
        st0.cy() + step * mn0.cy() + ss1o2 * mn0.uy(),
        st0.cz() + step * mn0.cz() + ss1o2 * mn0.uz(),
        st0.ux() + step * mn0.ux(),
        st0.uy() + step * mn0.uy(),
        st0.uz() + step * mn0.uz()
    );
    if (mat) {
        Double_t eta    = st0.eta() + step * mn0.e();
        Bool_t   is_mch = (MGNumc::Compare(eta) == eta_sign);
        if (is_mch) st1.set_eta(eta);
        else        return false;
    }
    MatPhyFld&& mp1 = MatPhy::Get(mfld, st1);
    MotionFunc mn1(st1, mp1);
   
    part.set_state_with_cos(
        st0.cx() + step * mn0.cx() + ss1o6 * (MGMath::TWO * mn0.ux() + mn1.ux()),
        st0.cy() + step * mn0.cy() + ss1o6 * (MGMath::TWO * mn0.uy() + mn1.uy()),
        st0.cz() + step * mn0.cz() + ss1o6 * (MGMath::TWO * mn0.uz() + mn1.uz()),
        st0.ux() + s1o2 * (mn0.ux() + mn1.ux()),
        st0.uy() + s1o2 * (mn0.uy() + mn1.uy()),
        st0.uz() + s1o2 * (mn0.uz() + mn1.uz())
    );
    if (mat) {
        Double_t eta    = st0.eta() + s1o2 * (mn0.e() + mn1.e());
        Bool_t   is_mch = (MGNumc::Compare(eta) == eta_sign);
        if (is_mch) part.set_eta(eta);
        else        return false;
    }
 
    /*
    if (withJb) {
        TransferFunc tf0(st0, marg, mp0);
        TransferFunc tf1(st1, marg, mp1);

        PhyJb jb1(PhyJb::Type::kIdentity);
        
        jb1.gg(JPX, JUX) += step * tf0.cu(X);
        jb1.gg(JPY, JUY) += step * tf0.cu(Y);
        
        jb1.gg(JPX, JUX) += ss1o2 * tf0.uu(X, X);
        jb1.gg(JPX, JUY) += ss1o2 * tf0.uu(X, Y);
        jb1.gg(JPX, JEA) += ss1o2 * tf0.ue(X);
        
        jb1.gg(JPY, JUX) += ss1o2 * tf0.uu(Y, X);
        jb1.gg(JPY, JUY) += ss1o2 * tf0.uu(Y, Y);
        jb1.gg(JPY, JEA) += ss1o2 * tf0.ue(Y);
        
        jb1.gg(JUX, JUX) += step * tf0.uu(X, X);
        jb1.gg(JUX, JUY) += step * tf0.uu(X, Y);
        jb1.gg(JUX, JEA) += step * tf0.ue(X);
        
        jb1.gg(JUY, JUX) += step * tf0.uu(Y, X);
        jb1.gg(JUY, JUY) += step * tf0.uu(Y, Y);
        jb1.gg(JUY, JEA) += step * tf0.ue(Y);

        if (mat && marg.mscat()) {
            jb1.gl(JPX, JTAU) += ss1o2 * tf0.ut(X);
            jb1.gl(JPX, JRHO) += ss1o2 * tf0.ur(X);
            
            jb1.gl(JPY, JTAU) += ss1o2 * tf0.ut(Y);
            jb1.gl(JPY, JRHO) += ss1o2 * tf0.ur(Y);
            
            jb1.gl(JUX, JTAU) += step * tf0.ut(X);
            jb1.gl(JUX, JRHO) += step * tf0.ur(X);
            
            jb1.gl(JUY, JTAU) += step * tf0.ut(Y);
            jb1.gl(JUY, JRHO) += step * tf0.ur(Y);
        }
        if (mat && marg.eloss()) {
            jb1.gg(JEA, JEA) += step * tf0.ee();
            
            jb1.gl(JEA, JION) += step * tf0.ei();
            jb1.gl(JEA, JBRM) += step * tf0.eb();
        }

        TransferPhyJb tj1(mat, tf1, jb1);

        phyJb->init(PhyJb::Type::kIdentity);
        phyJb->set_mat(mat, mfld.num_rad_len());
        
        phyJb->gg(JPX, JUX) += step * tf0.cu(X);
        phyJb->gg(JPY, JUY) += step * tf0.cu(Y);
        
        phyJb->gg(JPX, JUX) += ss1o6 * (MGMath::TWO * tf0.uu(X, X) + tj1.uu(X, X));
        phyJb->gg(JPX, JUY) += ss1o6 * (MGMath::TWO * tf0.uu(X, Y) + tj1.uu(X, Y));
        phyJb->gg(JPX, JEA) += ss1o6 * (MGMath::TWO * tf0.ue(X)    + tj1.ue(X)   );
        
        phyJb->gg(JPY, JUX) += ss1o6 * (MGMath::TWO * tf0.uu(Y, X) + tj1.uu(Y, X));
        phyJb->gg(JPY, JUY) += ss1o6 * (MGMath::TWO * tf0.uu(Y, Y) + tj1.uu(Y, Y));
        phyJb->gg(JPY, JEA) += ss1o6 * (MGMath::TWO * tf0.ue(Y)    + tj1.ue(Y)   );
        
        phyJb->gg(JUX, JUX) += s1o2 * (tf0.uu(X, X) + tj1.uu(X, X));
        phyJb->gg(JUX, JUY) += s1o2 * (tf0.uu(X, Y) + tj1.uu(X, Y));
        phyJb->gg(JUX, JEA) += s1o2 * (tf0.ue(X)    + tj1.ue(X)   );
        
        phyJb->gg(JUY, JUX) += s1o2 * (tf0.uu(Y, X) + tj1.uu(Y, X));
        phyJb->gg(JUY, JUY) += s1o2 * (tf0.uu(Y, Y) + tj1.uu(Y, Y));
        phyJb->gg(JUY, JEA) += s1o2 * (tf0.ue(Y)    + tj1.ue(Y)   );

        if (mat && marg.mscat()) {
            phyJb->gl(JPX, JTAU) += ss1o6 * (MGMath::TWO * tf0.ut(X) + tj1.ut(X));
            phyJb->gl(JPX, JRHO) += ss1o6 * (MGMath::TWO * tf0.ur(X) + tj1.ur(X));
            
            phyJb->gl(JPY, JTAU) += ss1o6 * (MGMath::TWO * tf0.ut(Y) + tj1.ut(Y));
            phyJb->gl(JPY, JRHO) += ss1o6 * (MGMath::TWO * tf0.ur(Y) + tj1.ur(Y));
            
            phyJb->gl(JUX, JTAU) += s1o2 * (tf0.ut(X) + tj1.ut(X));
            phyJb->gl(JUX, JRHO) += s1o2 * (tf0.ur(X) + tj1.ur(X));
            
            phyJb->gl(JUY, JTAU) += s1o2 * (tf0.ut(Y) + tj1.ut(Y));
            phyJb->gl(JUY, JRHO) += s1o2 * (tf0.ur(Y) + tj1.ur(Y));
        }
        if (mat && marg.eloss()) {
            phyJb->gg(JEA, JEA) += s1o2 * (tf0.ee() + tj1.ee());
            
            phyJb->gl(JEA, JION) += s1o2 * (tf0.ei() + tj1.ei());
            phyJb->gl(JEA, JBRM) += s1o2 * (tf0.eb() + tj1.eb());
        }
    }
    */

    return true;
}


Bool_t PropMgnt::PropWithRungeKuttaNystrom(const Double_t step, PhySt& part, const MatFld& mfld, PhyJb* phyJb) {
    Bool_t        withJb = (phyJb != nullptr);
    Bool_t           mat = (part.field() && mfld());
    Short_t     eta_sign = part.eta_sign();

    Double_t s1o2  = MGMath::HALF * step;
    Double_t s1o6  = MGMath::ONE_TO_SIX * step;
    Double_t sstep = step * step;
    Double_t ss1o2 = MGMath::HALF * sstep;
    Double_t ss1o6 = MGMath::ONE_TO_SIX * sstep;
    Double_t ss1o8 = MGMath::ONE_TO_EIGHT * sstep;

    PhySt st0 = part;
    MatPhyFld&& mp0 = MatPhy::Get(mfld, st0);
    MotionFunc mn0(st0, mp0);

    PhySt st1 = part;
    st1.set_state_with_cos(
        st0.cx() + s1o2 * mn0.cx() + ss1o8 * mn0.ux(),
        st0.cy() + s1o2 * mn0.cy() + ss1o8 * mn0.uy(),
        st0.cz() + s1o2 * mn0.cz() + ss1o8 * mn0.uz(),
        st0.ux() + s1o2 * mn0.ux(),
        st0.uy() + s1o2 * mn0.uy(),
        st0.uz() + s1o2 * mn0.uz()
    );
    if (mat) {
        Double_t eta    = st0.eta() + s1o2 * mn0.e();
        Bool_t   is_mch = (MGNumc::Compare(eta) == eta_sign);
        if (is_mch) st1.set_eta(eta);
        else        return false;
    }
    MatPhyFld&& mp1 = MatPhy::Get(mfld, st1);
    MotionFunc mn1(st1, mp1);
    
    PhySt st2 = part;
    st2.set_state_with_cos(
        st0.cx() + s1o2 * mn0.cx() + ss1o8 * mn0.ux(),
        st0.cy() + s1o2 * mn0.cy() + ss1o8 * mn0.uy(),
        st0.cz() + s1o2 * mn0.cz() + ss1o8 * mn0.uz(),
        st0.ux() + s1o2 * mn1.ux(),
        st0.uy() + s1o2 * mn1.uy(),
        st0.uz() + s1o2 * mn1.uz()
    );
    if (mat) {
        Double_t eta    = st0.eta() + s1o2 * mn1.e();
        Bool_t   is_mch = (MGNumc::Compare(eta) == eta_sign);
        if (is_mch) st2.set_eta(eta);
        else        return false;
    }
    MatPhyFld&& mp2 = MatPhy::Get(mfld, st2);
    MotionFunc mn2(st2, mp2);
    
    PhySt st3 = part;
    st3.set_state_with_cos(
        st0.cx() + step * mn0.cx() + ss1o2 * mn2.ux(),
        st0.cy() + step * mn0.cy() + ss1o2 * mn2.uy(),
        st0.cz() + step * mn0.cz() + ss1o2 * mn2.uz(),
        st0.ux() + step * mn2.ux(),
        st0.uy() + step * mn2.uy(),
        st0.uz() + step * mn2.uz()
    );
    if (mat) {
        Double_t eta    = st0.eta() + step * mn2.e();
        Bool_t   is_mch = (MGNumc::Compare(eta) == eta_sign);
        if (is_mch) st3.set_eta(eta);
        else        return false;
    }
    MatPhyFld&& mp3 = MatPhy::Get(mfld, st3);
    MotionFunc mn3(st3, mp3);
    
    part.set_state_with_cos(
        st0.cx() + step * mn0.cx() + ss1o6 * (mn0.ux() + mn1.ux() + mn2.ux()),
        st0.cy() + step * mn0.cy() + ss1o6 * (mn0.uy() + mn1.uy() + mn2.uy()),
        st0.cz() + step * mn0.cz() + ss1o6 * (mn0.uz() + mn1.uz() + mn2.uz()),
        st0.ux() + s1o6 * (mn0.ux() + MGMath::TWO * mn1.ux() + MGMath::TWO * mn2.ux() + mn3.ux()),
        st0.uy() + s1o6 * (mn0.uy() + MGMath::TWO * mn1.uy() + MGMath::TWO * mn2.uy() + mn3.uy()),
        st0.uz() + s1o6 * (mn0.uz() + MGMath::TWO * mn1.uz() + MGMath::TWO * mn2.uz() + mn3.uz())
    );
    if (mat) {
        Double_t eta    = st0.eta() + s1o6 * (mn0.e() + MGMath::TWO * mn1.e() + MGMath::TWO * mn2.e() + mn3.e());
        Bool_t   is_mch = (MGNumc::Compare(eta) == eta_sign);
        if (is_mch) part.set_eta(eta);
        else        return false;
    }
    
    /*
    if (withJb) {
        TransferFunc tf0(st0, marg, mp0);
        TransferFunc tf1(st1, marg, mp1);
        TransferFunc tf2(st2, marg, mp2);
        TransferFunc tf3(st3, marg, mp3);

        PhyJb jb1(PhyJb::Type::kIdentity);
        
        jb1.gg(JPX, JUX) += s1o2 * tf0.cu(X);
        jb1.gg(JPY, JUY) += s1o2 * tf0.cu(Y);
        
        jb1.gg(JPX, JUX) += ss1o8 * tf0.uu(X, X);
        jb1.gg(JPX, JUY) += ss1o8 * tf0.uu(X, Y);
        jb1.gg(JPX, JEA) += ss1o8 * tf0.ue(X);
        
        jb1.gg(JPY, JUX) += ss1o8 * tf0.uu(Y, X);
        jb1.gg(JPY, JUY) += ss1o8 * tf0.uu(Y, Y);
        jb1.gg(JPY, JEA) += ss1o8 * tf0.ue(Y);
        
        jb1.gg(JUX, JUX) += s1o2 * tf0.uu(X, X);
        jb1.gg(JUX, JUY) += s1o2 * tf0.uu(X, Y);
        jb1.gg(JUX, JEA) += s1o2 * tf0.ue(X);
        
        jb1.gg(JUY, JUX) += s1o2 * tf0.uu(Y, X);
        jb1.gg(JUY, JUY) += s1o2 * tf0.uu(Y, Y);
        jb1.gg(JUY, JEA) += s1o2 * tf0.ue(Y);

        if (mat && marg.mscat()) {
            jb1.gl(JPX, JTAU) += ss1o8 * tf0.ut(X);
            jb1.gl(JPX, JRHO) += ss1o8 * tf0.ur(X);
            
            jb1.gl(JPY, JTAU) += ss1o8 * tf0.ut(Y);
            jb1.gl(JPY, JRHO) += ss1o8 * tf0.ur(Y);
            
            jb1.gl(JUX, JTAU) += s1o2 * tf0.ut(X);
            jb1.gl(JUX, JRHO) += s1o2 * tf0.ur(X);
            
            jb1.gl(JUY, JTAU) += s1o2 * tf0.ut(Y);
            jb1.gl(JUY, JRHO) += s1o2 * tf0.ur(Y);
        }
        if (mat && marg.eloss()) {
            jb1.gg(JEA, JEA) += s1o2 * tf0.ee();
            
            jb1.gl(JEA, JION) += s1o2 * tf0.ei();
            jb1.gl(JEA, JBRM) += s1o2 * tf0.eb();
        }
        
        TransferPhyJb tj1(mat, tf1, jb1);
        
        PhyJb jb2(PhyJb::Type::kIdentity);
        
        jb2.gg(JPX, JUX) += s1o2 * tf0.cu(X);
        jb2.gg(JPY, JUY) += s1o2 * tf0.cu(Y);
        
        jb2.gg(JPX, JUX) += ss1o8 * tf0.uu(X, X);
        jb2.gg(JPX, JUY) += ss1o8 * tf0.uu(X, Y);
        jb2.gg(JPX, JEA) += ss1o8 * tf0.ue(X);
        
        jb2.gg(JPY, JUX) += ss1o8 * tf0.uu(Y, X);
        jb2.gg(JPY, JUY) += ss1o8 * tf0.uu(Y, Y);
        jb2.gg(JPY, JEA) += ss1o8 * tf0.ue(Y);
        
        jb2.gg(JUX, JUX) += s1o2 * tj1.uu(X, X);
        jb2.gg(JUX, JUY) += s1o2 * tj1.uu(X, Y);
        jb2.gg(JUX, JEA) += s1o2 * tj1.ue(X);
        
        jb2.gg(JUY, JUX) += s1o2 * tj1.uu(Y, X);
        jb2.gg(JUY, JUY) += s1o2 * tj1.uu(Y, Y);
        jb2.gg(JUY, JEA) += s1o2 * tj1.ue(Y);

        if (mat && marg.mscat()) {
            jb2.gl(JPX, JTAU) += ss1o8 * tf0.ut(X);
            jb2.gl(JPX, JRHO) += ss1o8 * tf0.ur(X);
            
            jb2.gl(JPY, JTAU) += ss1o8 * tf0.ut(Y);
            jb2.gl(JPY, JRHO) += ss1o8 * tf0.ur(Y);
            
            jb2.gl(JUX, JTAU) += s1o2 * tj1.ut(X);
            jb2.gl(JUX, JRHO) += s1o2 * tj1.ur(X);
            
            jb2.gl(JUY, JTAU) += s1o2 * tj1.ut(Y);
            jb2.gl(JUY, JRHO) += s1o2 * tj1.ur(Y);
        }
        if (mat && marg.eloss()) {
            jb2.gg(JEA, JEA) += s1o2 * tj1.ee();
            
            jb2.gl(JEA, JION) += s1o2 * tj1.ei();
            jb2.gl(JEA, JBRM) += s1o2 * tj1.eb();
        }
        
        TransferPhyJb tj2(mat, tf2, jb2);
        
        PhyJb jb3(PhyJb::Type::kIdentity);
        
        jb3.gg(JPX, JUX) += step * tf0.cu(X);
        jb3.gg(JPY, JUY) += step * tf0.cu(Y);
        
        jb3.gg(JPX, JUX) += ss1o2 * tf2.uu(X, X);
        jb3.gg(JPX, JUY) += ss1o2 * tf2.uu(X, Y);
        jb3.gg(JPX, JEA) += ss1o2 * tf2.ue(X);
        
        jb3.gg(JPY, JUX) += ss1o2 * tf2.uu(Y, X);
        jb3.gg(JPY, JUY) += ss1o2 * tf2.uu(Y, Y);
        jb3.gg(JPY, JEA) += ss1o2 * tf2.ue(Y);
        
        jb3.gg(JUX, JUX) += step * tj2.uu(X, X);
        jb3.gg(JUX, JUY) += step * tj2.uu(X, Y);
        jb3.gg(JUX, JEA) += step * tj2.ue(X);
        
        jb3.gg(JUY, JUX) += step * tj2.uu(Y, X);
        jb3.gg(JUY, JUY) += step * tj2.uu(Y, Y);
        jb3.gg(JUY, JEA) += step * tj2.ue(Y);

        if (mat && marg.mscat()) {
            jb3.gl(JPX, JTAU) += ss1o2 * tf2.ut(X);
            jb3.gl(JPX, JRHO) += ss1o2 * tf2.ur(X);
            
            jb3.gl(JPY, JTAU) += ss1o2 * tf2.ut(Y);
            jb3.gl(JPY, JRHO) += ss1o2 * tf2.ur(Y);
            
            jb3.gl(JUX, JTAU) += step * tj2.ut(X);
            jb3.gl(JUX, JRHO) += step * tj2.ur(X);
            
            jb3.gl(JUY, JTAU) += step * tj2.ut(Y);
            jb3.gl(JUY, JRHO) += step * tj2.ur(Y);
        }
        if (mat && marg.eloss()) {
            jb3.gg(JEA, JEA) += step * tj2.ee();

            jb3.gl(JEA, JION) += step * tj2.ei();
            jb3.gl(JEA, JBRM) += step * tj2.eb();
        }
        
        TransferPhyJb tj3(mat, tf3, jb3);
        
        phyJb->init(PhyJb::Type::kIdentity);
        phyJb->set_mat(mat, mfld.num_rad_len());
        
        phyJb->gg(JPX, JUX) += step * tf0.cu(X);
        phyJb->gg(JPY, JUY) += step * tf0.cu(Y);
        
        phyJb->gg(JPX, JUX) += ss1o6 * (tf0.uu(X, X) + tj1.uu(X, X) + tj2.uu(X, X));
        phyJb->gg(JPX, JUY) += ss1o6 * (tf0.uu(X, Y) + tj1.uu(X, Y) + tj2.uu(X, Y));
        phyJb->gg(JPX, JEA) += ss1o6 * (tf0.ue(X)    + tj1.ue(X)    + tj2.ue(X)   );
        
        phyJb->gg(JPY, JUX) += ss1o6 * (tf0.uu(Y, X) + tj1.uu(Y, X) + tj2.uu(Y, X));
        phyJb->gg(JPY, JUY) += ss1o6 * (tf0.uu(Y, Y) + tj1.uu(Y, Y) + tj2.uu(Y, Y));
        phyJb->gg(JPY, JEA) += ss1o6 * (tf0.ue(Y)    + tj1.ue(Y)    + tj2.ue(Y)   );
        
        phyJb->gg(JUX, JUX) += s1o6 * (tf0.uu(X, X) + MGMath::TWO * tj1.uu(X, X) + MGMath::TWO * tj2.uu(X, X) + tj3.uu(X, X));
        phyJb->gg(JUX, JUY) += s1o6 * (tf0.uu(X, Y) + MGMath::TWO * tj1.uu(X, Y) + MGMath::TWO * tj2.uu(X, Y) + tj3.uu(X, Y));
        phyJb->gg(JUX, JEA) += s1o6 * (tf0.ue(X)    + MGMath::TWO * tj1.ue(X)    + MGMath::TWO * tj2.ue(X)    + tj3.ue(X)   );
        
        phyJb->gg(JUY, JUX) += s1o6 * (tf0.uu(Y, X) + MGMath::TWO * tj1.uu(Y, X) + MGMath::TWO * tj2.uu(Y, X) + tj3.uu(Y, X));
        phyJb->gg(JUY, JUY) += s1o6 * (tf0.uu(Y, Y) + MGMath::TWO * tj1.uu(Y, Y) + MGMath::TWO * tj2.uu(Y, Y) + tj3.uu(Y, Y));
        phyJb->gg(JUY, JEA) += s1o6 * (tf0.ue(Y)    + MGMath::TWO * tj1.ue(Y)    + MGMath::TWO * tj2.ue(Y)    + tj3.ue(Y)   );

        if (mat && marg.mscat()) {
            phyJb->gl(JPX, JTAU) += ss1o6 * (tf0.ut(X) + tj1.ut(X) + tj2.ut(X));
            phyJb->gl(JPX, JRHO) += ss1o6 * (tf0.ur(X) + tj1.ur(X) + tj2.ur(X));
            
            phyJb->gl(JPY, JTAU) += ss1o6 * (tf0.ut(Y) + tj1.ut(Y) + tj2.ut(Y));
            phyJb->gl(JPY, JRHO) += ss1o6 * (tf0.ur(Y) + tj1.ur(Y) + tj2.ur(Y));
            
            phyJb->gl(JUX, JTAU) += s1o6 * (tf0.ut(X) + MGMath::TWO * tj1.ut(X) + MGMath::TWO * tj2.ut(X) + tj3.ut(X));
            phyJb->gl(JUX, JRHO) += s1o6 * (tf0.ur(X) + MGMath::TWO * tj1.ur(X) + MGMath::TWO * tj2.ur(X) + tj3.ur(X));
            
            phyJb->gl(JUY, JTAU) += s1o6 * (tf0.ut(Y) + MGMath::TWO * tj1.ut(Y) + MGMath::TWO * tj2.ut(Y) + tj3.ut(Y));
            phyJb->gl(JUY, JRHO) += s1o6 * (tf0.ur(Y) + MGMath::TWO * tj1.ur(Y) + MGMath::TWO * tj2.ur(Y) + tj3.ur(Y));
        }
        if (mat && marg.eloss()) {
            phyJb->gg(JEA, JEA) += s1o6 * (tf0.ee() + MGMath::TWO * tj1.ee() + MGMath::TWO * tj2.ee() + tj3.ee());
            
            phyJb->gl(JEA, JION) += s1o6 * (tf0.ei() + MGMath::TWO * tj1.ei() + MGMath::TWO * tj2.ei() + tj3.ei());
            phyJb->gl(JEA, JBRM) += s1o6 * (tf0.eb() + MGMath::TWO * tj1.eb() + MGMath::TWO * tj2.eb() + tj3.eb());
        }
    }
    */
    
    return true;
}


} // namespace TrackSys


#endif // __TRACKLibs_Prop_C__
