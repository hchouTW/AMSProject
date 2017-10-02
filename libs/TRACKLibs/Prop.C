#ifndef __TRACKLibs_Prop_C__
#define __TRACKLibs_Prop_C__


#ifdef __HAS_AMS_OFFICE_LIBS__
#include <TrFit.h>
#endif // __HAS_AMS_OFFICE_LIBS__


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


MotionFunc::MotionFunc(PhySt& part, const MatPhyFld* mphy) {
    Bool_t field = (part.field() && mphy != nullptr && (*mphy)());
    Double_t Lambda = PROP_FACT * part.info().chrg_to_mass();
    
    MagFld&& mag = MagMgnt::Get(part.c());
    SVecD<3>&& crsub = LA::Cross(part.u(), mag());
    
    zeta_c_ = std::move(part.u());
    zeta_u_ = std::move((Lambda * part.eta()) * crsub);
    
    if (field) orth_.reset(part.u(), mag());
    
    if (field && part.arg().eloss()) zeta_e_ = part.arg().tune_ion() * mphy->eloss_ion_mpv();
    else                             zeta_e_ = MGMath::ZERO;
}


TransferFunc::TransferFunc(PhySt& part, const MatPhyFld* mphy) {
    Bool_t field = (part.field() && mphy != nullptr && (*mphy)());
    Double_t Lambda = PROP_FACT * part.info().chrg_to_mass();
    
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

    if (field && part.arg().eloss()) kappa_ee_ = part.arg().tune_ion() * mphy->eloss_ion_mpv();
    else                             kappa_ee_ = MGMath::ZERO;
}
 

void PhyJb::init() {
    field_ = false;
    jb_gg_ = std::move(SMtxId()); 
    jb_gl_ = std::move(SMtxDGL()); 
}
 

void PhyJb::set(PhySt& part) {
    if (!part.field()) return;
    VirtualPhySt& vst = part.vst();
    field_ = (part.field() && vst());
    if (!field_) return;
    
    if (part.arg().mscat()) {
        jb_gl_(JPX, JTAU) = vst.mscatcu() * vst.tau(X);
        jb_gl_(JPY, JTAU) = vst.mscatcu() * vst.tau(Y);
        jb_gl_(JUX, JTAU) = vst.mscatu()  * vst.tau(X);
        jb_gl_(JUY, JTAU) = vst.mscatu()  * vst.tau(Y);
        jb_gl_(JEA, JTAU) = MGMath::ZERO;
        
        jb_gl_(JPX, JRHO) = vst.mscatcu() * vst.rho(X);
        jb_gl_(JPY, JRHO) = vst.mscatcu() * vst.rho(Y);
        jb_gl_(JUX, JRHO) = vst.mscatu()  * vst.rho(X);
        jb_gl_(JUY, JRHO) = vst.mscatu()  * vst.rho(Y);
        jb_gl_(JEA, JRHO) = MGMath::ZERO;
    }
}


void PhyJb::multiplied(PhyJb& phyJb) {
    if (field_) jb_gl_ = std::move(phyJb.gg() * jb_gl_);
    jb_gg_ = std::move(phyJb.gg() * jb_gg_);
}
        

void PhyJb::print() const {
    std::string printStr;
    printStr += STR_FMT("============================================= PhyJb =============================================\n");
    printStr += STR_FMT("%12.7f %12.7f %12.7f %12.7f %12.7f        %12.7f %12.7f\n", jb_gg_(JPX, JPX), jb_gg_(JPX, JPY), jb_gg_(JPX, JUX), jb_gg_(JPX, JUY), jb_gg_(JPX, JEA), jb_gl_(JPX, JTAU), jb_gl_(JPX, JRHO));
    printStr += STR_FMT("%12.7f %12.7f %12.7f %12.7f %12.7f        %12.7f %12.7f\n", jb_gg_(JPY, JPX), jb_gg_(JPY, JPY), jb_gg_(JPY, JUX), jb_gg_(JPY, JUY), jb_gg_(JPY, JEA), jb_gl_(JPY, JTAU), jb_gl_(JPY, JRHO));
    printStr += STR_FMT("%12.7f %12.7f %12.7f %12.7f %12.7f        %12.7f %12.7f\n", jb_gg_(JUX, JPX), jb_gg_(JUX, JPY), jb_gg_(JUX, JUX), jb_gg_(JUX, JUY), jb_gg_(JUX, JEA), jb_gl_(JUX, JTAU), jb_gl_(JUX, JRHO));
    printStr += STR_FMT("%12.7f %12.7f %12.7f %12.7f %12.7f        %12.7f %12.7f\n", jb_gg_(JUY, JPX), jb_gg_(JUY, JPY), jb_gg_(JUY, JUX), jb_gg_(JUY, JUY), jb_gg_(JUY, JEA), jb_gl_(JUY, JTAU), jb_gl_(JUY, JRHO));
    printStr += STR_FMT("%12.7f %12.7f %12.7f %12.7f %12.7f        %12.7f %12.7f\n", jb_gg_(JEA, JPX), jb_gg_(JEA, JPY), jb_gg_(JEA, JUX), jb_gg_(JEA, JUY), jb_gg_(JEA, JEA), jb_gl_(JEA, JTAU), jb_gl_(JEA, JRHO));
    printStr += STR_FMT("=================================================================================================\n");
    COUT(printStr); 
}


TransferPhyJb::TransferPhyJb(const TransferFunc& tf, PhyJb& jb) {
    uu_(X, X) = tf.uu(X, X) * jb.gg(JUX, JUX) + tf.uu(X, Y) * jb.gg(JUY, JUX) + tf.ue(X) * jb.gg(JEA, JUX);
    uu_(X, Y) = tf.uu(X, X) * jb.gg(JUX, JUY) + tf.uu(X, Y) * jb.gg(JUY, JUY) + tf.ue(X) * jb.gg(JEA, JUY);
    
    uu_(Y, X) = tf.uu(Y, X) * jb.gg(JUX, JUX) + tf.uu(Y, Y) * jb.gg(JUY, JUX) + tf.ue(Y) * jb.gg(JEA, JUX);
    uu_(Y, Y) = tf.uu(Y, X) * jb.gg(JUX, JUY) + tf.uu(Y, Y) * jb.gg(JUY, JUY) + tf.ue(Y) * jb.gg(JEA, JUY);

    ue_(X) = tf.uu(X, X) * jb.gg(JUX, JEA) + tf.uu(X, Y) * jb.gg(JUY, JEA) + tf.ue(X) * jb.gg(JEA, JEA);
    ue_(Y) = tf.uu(Y, X) * jb.gg(JUX, JEA) + tf.uu(Y, Y) * jb.gg(JUY, JEA) + tf.ue(Y) * jb.gg(JEA, JEA);

    ee_ = tf.ee() * jb.gg(JEA, JEA);
}


void PropPhyCal::init() {
    sw_mscat_      = false;
    sw_eloss_      = false;
    eta0_abs_      = MGMath::ZERO;
    field_         = false;
    sign_          = 1; 
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
    eloss_brm_men_ = MGMath::ZERO;
    vec_vac_.clear();
    vec_len_.clear();
    vec_eft_.clear();
    vec_invloc_.clear();
    vec_invlocsqr_.clear();
    vec_mscat_.clear();
    vec_ieta_kpa_.clear();
    vec_ieta_sgm_.clear();
}


void PropPhyCal::normalized(const MatFld& mfld, const PhySt& part) {
    if (sw_mscat_ && !MGNumc::EqualToZero(mscatu_)) {
        tau_.Unit();
        rho_.Unit();
        mscatu_  = std::sqrt(mscatu_);
        //mscatcl_ = (mscatcl_ / mscatu_) * MGMath::INV_SQRT_TWELVE; // testcode
        mscatcu_ = MGMath::ZERO;
    
        //----------------- AMS
        mscatcl_ = (mscatcl_ * mscatu_) * MGMath::INV_SQRT_THREE;
        //mscatcl_ = (mscatcl_ * mscatu_) * MGMath::INV_SQRT_TWELVE;
        //----------------- AMS

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
            //Double_t est_eftl = (MGMath::ONE - (MGMath::ONE - MGMath::INV_SQRT_THREE) * est_efft);
            Double_t est_eftl = (MGMath::ONE); // testcode
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

        if (!MGNumc::Valid(mscatu_)  || MGNumc::Compare(mscatu_)  <= 0) mscatu_  = MGMath::ZERO;
        if (!MGNumc::Valid(mscatcu_) || MGNumc::Compare(mscatcu_) <= 0) mscatcu_ = MGMath::ZERO;
        if (!MGNumc::Valid(mscatcl_) || MGNumc::Compare(mscatcl_) <= 0) mscatcl_ = MGMath::ZERO;
        field_ = true;
    }
    if (sw_eloss_) {
        MatPhyFld&& mpfld = MatPhy::Get(mfld, part);
        if (mpfld()) {
            Double_t ttlwgt = MGMath::ZERO;
            for (auto&& wgt : vec_ieta_sgm_) ttlwgt += wgt;
            Double_t eloss_ieta_sgm = std::sqrt(ttlwgt);
            if (MGNumc::Compare(eloss_ieta_sgm) > 0) {
                Double_t ttlkpa = MGMath::ZERO;
                for (auto&& kpa : vec_ieta_kpa_) ttlkpa += kpa;
                Double_t eloss_eta_kpa = (ttlkpa / ttlwgt); 
                
                Double_t eloss_eta_sgm = MGMath::HALF * (std::sqrt(MGMath::ONE/ttlwgt + MGMath::FOUR * part.eta() * part.eta()) - (MGMath::ONE / eloss_ieta_sgm));
                Double_t eloss_ion_sgm = (eloss_eta_sgm / eloss_eta_kpa / part.eta_abs());

                eloss_ion_kpa_ = eloss_eta_kpa;
                eloss_ion_sgm_ = eloss_ion_sgm;
                eloss_ion_mpv_ = (std::fabs(part.eta_abs() - eta0_abs_) / std::sqrt(part.eta_abs() * eta0_abs_));
                
                //------------ AMS
                eloss_ion_kpa_ = mpfld.eloss_ion_kpa();
                eloss_ion_sgm_ = mpfld.eloss_ion_sgm();
                eloss_ion_mpv_ = (std::fabs(part.eta_abs() - eta0_abs_) / std::sqrt(part.eta_abs() * eta0_abs_));
                //------------ AMS
                
                if (!MGNumc::Valid(eloss_ion_kpa_) || MGNumc::Compare(eloss_ion_kpa_) <= 0) eloss_ion_kpa_ = MGMath::ZERO;
                if (!MGNumc::Valid(eloss_ion_sgm_) || MGNumc::Compare(eloss_ion_sgm_) <= 0) eloss_ion_sgm_ = MGMath::ZERO;
                if (!MGNumc::Valid(eloss_ion_mpv_) || MGNumc::Compare(eloss_ion_mpv_) <= 0) eloss_ion_mpv_ = MGMath::ZERO;
            }

            eloss_brm_men_ = mpfld.eloss_brm_men();
            if (!MGNumc::Valid(eloss_brm_men_) || MGNumc::Compare(eloss_brm_men_) <= 0) eloss_brm_men_ = MGMath::ZERO;
            field_ = true;
        }
    }
}


void PropPhyCal::push(PhySt& part, const MatPhyFld& mpfld, const SVecD<3>& tau, const SVecD<3>& rho) {
    len_ += mpfld.len();
    if (mpfld()) nrl_ += mpfld.num_rad_len();
    if (sw_mscat_) {
        Bool_t vac = !mpfld();
        Double_t len = mpfld.len();
        Double_t mscat = mpfld.mult_scat_sgm(); 
        Double_t wgt = mscat * mscat;
        vec_vac_.push_back(vac);
        vec_len_.push_back(len);
        if (vac) {
            vec_eft_.push_back(MGMath::ZERO);
            vec_invloc_.push_back(MGMath::ZERO);
            vec_invlocsqr_.push_back(MGMath::ZERO);
            vec_mscat_.push_back(MGMath::ZERO);
        }
        else {
            Double_t invloc    = (MGMath::ONE - mpfld.loc());
            Double_t invlocsqr = (MGMath::ONE - MGMath::TWO * mpfld.loc() + mpfld.locsqr());
            
            if (MGNumc::Compare(invloc)    <= 0) invloc    = MGMath::ZERO;
            if (MGNumc::Compare(invlocsqr) <= 0) invlocsqr = MGMath::ZERO;

            vec_eft_.push_back(len * mpfld.efft());
            vec_invloc_.push_back((invloc * len));
            vec_invlocsqr_.push_back((invlocsqr * len * len));
            vec_mscat_.push_back(wgt);
        }
        if (vac) return;
        mscatu_  += wgt;
        tau_     += wgt * tau;
        rho_     += wgt * rho;
        //mscatcl_ += wgt * (len * mpfld.efft()); // testcode
        //----------------- AMS
        mscatcl_ += (len * mpfld.efft()); // testcode
        //----------------- AMS
    }
    if (sw_eloss_ && mpfld()) {
        Double_t ion_eloss = (part.arg().tune_ion() * mpfld.eloss_ion_kpa() * mpfld.eloss_ion_sgm() / (MGMath::ONE + sign_ * part.arg().tune_ion() * mpfld.eloss_ion_mpv()));
        Double_t ion_sgm = ((ion_eloss / (MGMath::ONE - ion_eloss * ion_eloss)) / part.eta_abs());
        Double_t wgt = ion_sgm * ion_sgm;
        vec_ieta_kpa_.push_back(mpfld.eloss_ion_kpa() * wgt);
        vec_ieta_sgm_.push_back(wgt);
    }
}


void PropPhyCal::set_virtualPhySt(PhySt& part) const {
    part.vst().reset(field_);
    part.vst().set_sign(sign_);
    part.vst().set_len(len_);
    if (part.vst()()) {
        part.vst().set_nrl(nrl_);
        if (part.arg().mscat()) {
            part.vst().set_orth(tau_, rho_);
            part.vst().set_mscatu(mscatu_);
            part.vst().set_mscatc(mscatcu_, mscatcl_);
        }
        if (part.arg().eloss()) {
            part.vst().set_eloss_ion(eloss_ion_kpa_, eloss_ion_mpv_, eloss_ion_sgm_);
            part.vst().set_eloss_brm(eloss_brm_men_);
        }
    }
}


#ifdef __HAS_AMS_OFFICE_LIBS__
Bool_t PropMgnt::PropToZ_AMSLibs(const Double_t zcoo, PhySt& part) {
    Short_t sign = MGNumc::Compare(part.uz());
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
#endif // __HAS_AMS_OFFICE_LIBS__


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
    if (sign == 0) return MGMath::ZERO;

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
    if (signz == 0) return MGMath::ZERO;

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
        

Bool_t PropMgnt::Prop(const Double_t step, PhySt& part, MatFld* mfld, PhyJb* phyJb) {
    Bool_t withMf = (mfld != nullptr);
    Bool_t withJb = (phyJb != nullptr);
    if (withJb) phyJb->init();

    Long64_t iter     = 1;
    Bool_t   is_succ  = false;
    Double_t int_step = MGMath::ZERO;
    
    std::list<MatFld> mflds;
    PropPhyCal ppcal(part, step);
    while (iter <= LMTU_ITER && !is_succ) {
        Double_t res_step = step - int_step;
        Double_t cur_step = GetStep(part, res_step);

        Bool_t valid = false;
        PhyJb* curJb = ((withJb) ? (new PhyJb()) : nullptr);
        MatFld&& curMfld = (part.field() ? MatMgnt::Get(cur_step, part) : MatFld(std::fabs(cur_step)));
        
        switch (method_) {
            case Method::kRungeKuttaNystrom : valid = PropWithRungeKuttaNystrom(cur_step, part, curMfld, ppcal, curJb); break;
            case Method::kEulerHeun         : valid = PropWithEulerHeun(cur_step, part, curMfld, ppcal, curJb); break;
            case Method::kEuler             : valid = PropWithEuler(cur_step, part, curMfld, ppcal, curJb); break;
            default : break;
        }
        if (!valid) { 
            if (withJb && curJb != nullptr) delete curJb;
            curJb = nullptr;
            break;
        }

        mflds.push_back(curMfld);
        if (withJb) {
            phyJb->multiplied(*curJb);
            delete curJb;
        }

        iter++;
        int_step += cur_step;
        is_succ = (MGNumc::Compare(std::fabs(step - int_step), CONV_STEP) < 0);
    }
    MatFld&& mgfld = std::move(MatFld::Merge(mflds));
    if (withMf) *mfld = mgfld;
    
    ppcal.normalized(mgfld, part);
    ppcal.set_virtualPhySt(part);
    if (withJb) phyJb->set(part);

    return is_succ;
}


Bool_t PropMgnt::PropToZ(const Double_t zcoo, PhySt& part, MatFld* mfld, PhyJb* phyJb) {
    Bool_t withMf = (mfld != nullptr);
    Bool_t withJb = (phyJb != nullptr);
    if (withJb) phyJb->init();
    
    Long64_t iter     = 1;
    Bool_t   is_succ  = false;
    Double_t int_step = MGMath::ZERO;
    
    std::list<MatFld> mflds;
    PropPhyCal ppcal(part, part.uz()*(zcoo-part.cz()));
    while (iter <= LMTU_ITER && !is_succ) {
        Double_t res_stepz = zcoo - part.cz();
        Double_t cur_step  = GetStepToZ(part, res_stepz);

        Bool_t valid = false;
        PhyJb* curJb = ((withJb) ? (new PhyJb()) : nullptr);
        MatFld&& curMfld = (part.field() ? MatMgnt::Get(cur_step, part) : MatFld(std::fabs(cur_step)));

        switch (method_) {
            case Method::kRungeKuttaNystrom : valid = PropWithRungeKuttaNystrom(cur_step, part, curMfld, ppcal, curJb); break;
            case Method::kEulerHeun         : valid = PropWithEulerHeun(cur_step, part, curMfld, ppcal, curJb); break;
            case Method::kEuler             : valid = PropWithEuler(cur_step, part, curMfld, ppcal, curJb); break;
            default : break;
        }
        if (!valid) { 
            if (withJb && curJb != nullptr) delete curJb;
            curJb = nullptr;
            break;
        }

        mflds.push_back(curMfld);
        if (withJb) {
            phyJb->multiplied(*curJb);
            delete curJb;
        }

        iter++;
        int_step += cur_step;
        is_succ = (MGNumc::Compare(std::fabs(zcoo - part.cz()), CONV_STEP) < 0);
    }
    MatFld&& mgfld = std::move(MatFld::Merge(mflds));
    if (withMf) *mfld = mgfld;
    
    ppcal.normalized(mgfld, part);
    ppcal.set_virtualPhySt(part);
    if (withJb) phyJb->set(part);
   
    return is_succ;
}
        

Bool_t PropMgnt::PropWithMC(const Double_t step, PhySt& part, MatFld* mfld) {
    Bool_t is_succ = PropMgnt::Prop(step, part, mfld, nullptr);
    if (is_succ) part.symbk(true);
    else         part.zero();
    return is_succ;
}


Bool_t PropMgnt::PropToZWithMC(const Double_t zcoo, PhySt& part, MatFld* mfld) {
    Bool_t is_succ = PropMgnt::PropToZ(zcoo, part, mfld, nullptr);
    if (is_succ) part.symbk(true);
    else         part.zero();
    return is_succ;
}


Bool_t PropMgnt::PropWithEuler(const Double_t step, PhySt& part, const MatFld& mfld, PropPhyCal& ppcal, PhyJb* phyJb) {
    Short_t    prop_sign = MGNumc::Compare(step);
    Short_t     eta_sign = part.eta_sign();
    Bool_t        withMf = (part.field() && mfld());
    Bool_t     withEloss = (part.arg().eloss() && mfld());
    Bool_t        withJb = (phyJb != nullptr);
    
    Double_t sstep = step * step;
    Double_t ss1o2 = sstep * MGMath::HALF;
    
    Double_t step_ps = prop_sign;

    PhySt st0 = part;
    MatPhyFld&& mp0 = MatPhy::Get(mfld, st0);
    MotionFunc mn0(st0, &mp0);
    
    part.set_state_with_cos(
        st0.cx() + step * mn0.cx() + ss1o2 * mn0.ux(),
        st0.cy() + step * mn0.cy() + ss1o2 * mn0.uy(),
        st0.cz() + step * mn0.cz() + ss1o2 * mn0.uz(),
        st0.ux() + step * mn0.ux(),
        st0.uy() + step * mn0.uy(),
        st0.uz() + step * mn0.uz()
    );
    if (withEloss) {
        Double_t eta    = st0.eta() * (MGMath::ONE + step_ps * mn0.e());
        Bool_t   is_mch = (MGNumc::Compare(eta) == eta_sign);
        if (is_mch) part.set_eta(eta);
        else        return false;
    }
    
    ppcal.push(part, mp0, mn0.orth().tau(), mn0.orth().rho());
 
   
    if (withJb) {
        TransferFunc tf0(st0, &mp0);

        phyJb->init();
        
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

        if (withEloss) phyJb->gg(JEA, JEA) += step_ps * tf0.ee();
    }

    return true;
}


Bool_t PropMgnt::PropWithEulerHeun(const Double_t step, PhySt& part, const MatFld& mfld, PropPhyCal& ppcal, PhyJb* phyJb) {
    Short_t    prop_sign = MGNumc::Compare(step);
    Short_t     eta_sign = part.eta_sign();
    Bool_t        withMf = (part.field() && mfld());
    Bool_t     withEloss = (part.arg().eloss() && mfld());
    Bool_t        withJb = (phyJb != nullptr);
   
    Double_t s1o2  = step * MGMath::HALF;
    Double_t sstep = step * step;
    Double_t ss1o2 = sstep * MGMath::HALF;
    Double_t ss1o6 = sstep * MGMath::ONE_TO_SIX;
    
    Double_t step_ps = prop_sign;
    Double_t s1o2_ps = prop_sign * MGMath::HALF;

    PhySt st0 = part;
    MatPhyFld&& mp0 = MatPhy::Get(mfld, st0);
    MotionFunc mn0(st0, &mp0);
    
    PhySt st1 = part;
    st1.set_state_with_cos(
        st0.cx() + step * mn0.cx() + ss1o2 * mn0.ux(),
        st0.cy() + step * mn0.cy() + ss1o2 * mn0.uy(),
        st0.cz() + step * mn0.cz() + ss1o2 * mn0.uz(),
        st0.ux() + step * mn0.ux(),
        st0.uy() + step * mn0.uy(),
        st0.uz() + step * mn0.uz()
    );
    if (withEloss) {
        Double_t eta    = st0.eta() * (MGMath::ONE + step_ps * mn0.e());
        Bool_t   is_mch = (MGNumc::Compare(eta) == eta_sign);
        if (is_mch) st1.set_eta(eta);
        else        return false;
    }
    MatPhyFld&& mp1 = MatPhy::Get(mfld, st1);
    MotionFunc mn1(st1, &mp1);
   
    part.set_state_with_cos(
        st0.cx() + step * mn0.cx() + ss1o6 * (MGMath::TWO * mn0.ux() + mn1.ux()),
        st0.cy() + step * mn0.cy() + ss1o6 * (MGMath::TWO * mn0.uy() + mn1.uy()),
        st0.cz() + step * mn0.cz() + ss1o6 * (MGMath::TWO * mn0.uz() + mn1.uz()),
        st0.ux() + s1o2 * (mn0.ux() + mn1.ux()),
        st0.uy() + s1o2 * (mn0.uy() + mn1.uy()),
        st0.uz() + s1o2 * (mn0.uz() + mn1.uz())
    );
    if (withEloss) {
        Double_t eta    = st0.eta() * (MGMath::ONE + s1o2_ps * (mn0.e() + mn1.e()));
        Bool_t   is_mch = (MGNumc::Compare(eta) == eta_sign);
        if (is_mch) part.set_eta(eta);
        else        return false;
    }
    
    ppcal.push(part, mp0, mn0.orth().tau(), mn0.orth().rho());

    
    if (withJb) {
        TransferFunc tf0(st0, &mp0);
        TransferFunc tf1(st1, &mp1);

        PhyJb jb1;
        
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

        if (withEloss) jb1.gg(JEA, JEA) += step_ps * tf0.ee();

        TransferPhyJb tj1(tf1, jb1);

        phyJb->init();
        
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
        
        if (withEloss) phyJb->gg(JEA, JEA) += s1o2_ps * (tf0.ee() + tj1.ee());
    }

    return true;
}


Bool_t PropMgnt::PropWithRungeKuttaNystrom(const Double_t step, PhySt& part, const MatFld& mfld, PropPhyCal& ppcal, PhyJb* phyJb) {
    Short_t    prop_sign = MGNumc::Compare(step);
    Short_t     eta_sign = part.eta_sign();
    Bool_t        withMf = (part.field() && mfld());
    Bool_t     withEloss = (part.arg().eloss() && mfld());
    Bool_t        withJb = (phyJb != nullptr);

    Double_t s1o2  = step * MGMath::HALF;
    Double_t s1o6  = step * MGMath::ONE_TO_SIX;
    Double_t sstep = step * step;
    Double_t ss1o2 = sstep * MGMath::HALF;
    Double_t ss1o6 = sstep * MGMath::ONE_TO_SIX;
    Double_t ss1o8 = sstep * MGMath::ONE_TO_EIGHT;
    
    Double_t step_ps = prop_sign;
    Double_t s1o2_ps = prop_sign * MGMath::HALF;
    Double_t s1o6_ps = prop_sign * MGMath::ONE_TO_SIX;

    PhySt st0 = part;
    MatPhyFld&& mp0 = MatPhy::Get(mfld, st0);
    MotionFunc mn0(st0, &mp0);

    PhySt st1 = part;
    st1.set_state_with_cos(
        st0.cx() + s1o2 * mn0.cx() + ss1o8 * mn0.ux(),
        st0.cy() + s1o2 * mn0.cy() + ss1o8 * mn0.uy(),
        st0.cz() + s1o2 * mn0.cz() + ss1o8 * mn0.uz(),
        st0.ux() + s1o2 * mn0.ux(),
        st0.uy() + s1o2 * mn0.uy(),
        st0.uz() + s1o2 * mn0.uz()
    );
    if (withEloss) {
        Double_t eta    = st0.eta() * (MGMath::ONE + s1o2_ps * mn0.e());
        Bool_t   is_mch = (MGNumc::Compare(eta) == eta_sign);
        if (is_mch) st1.set_eta(eta);
        else        return false;
    }
    MatPhyFld&& mp1 = MatPhy::Get(mfld, st1);
    MotionFunc mn1(st1, &mp1);
    
    PhySt st2 = part;
    st2.set_state_with_cos(
        st0.cx() + s1o2 * mn0.cx() + ss1o8 * mn0.ux(),
        st0.cy() + s1o2 * mn0.cy() + ss1o8 * mn0.uy(),
        st0.cz() + s1o2 * mn0.cz() + ss1o8 * mn0.uz(),
        st0.ux() + s1o2 * mn1.ux(),
        st0.uy() + s1o2 * mn1.uy(),
        st0.uz() + s1o2 * mn1.uz()
    );
    if (withEloss) {
        Double_t eta    = st0.eta() * (MGMath::ONE + s1o2_ps * mn1.e());
        Bool_t   is_mch = (MGNumc::Compare(eta) == eta_sign);
        if (is_mch) st2.set_eta(eta);
        else        return false;
    }
    MatPhyFld&& mp2 = MatPhy::Get(mfld, st2);
    MotionFunc mn2(st2, &mp2);
    
    PhySt st3 = part;
    st3.set_state_with_cos(
        st0.cx() + step * mn0.cx() + ss1o2 * mn2.ux(),
        st0.cy() + step * mn0.cy() + ss1o2 * mn2.uy(),
        st0.cz() + step * mn0.cz() + ss1o2 * mn2.uz(),
        st0.ux() + step * mn2.ux(),
        st0.uy() + step * mn2.uy(),
        st0.uz() + step * mn2.uz()
    );
    if (withEloss) {
        Double_t eta    = st0.eta() * (MGMath::ONE + step_ps * mn2.e());
        Bool_t   is_mch = (MGNumc::Compare(eta) == eta_sign);
        if (is_mch) st3.set_eta(eta);
        else        return false;
    }
    MatPhyFld&& mp3 = MatPhy::Get(mfld, st3);
    MotionFunc mn3(st3, &mp3);
    
    part.set_state_with_cos(
        st0.cx() + step * mn0.cx() + ss1o6 * (mn0.ux() + mn1.ux() + mn2.ux()),
        st0.cy() + step * mn0.cy() + ss1o6 * (mn0.uy() + mn1.uy() + mn2.uy()),
        st0.cz() + step * mn0.cz() + ss1o6 * (mn0.uz() + mn1.uz() + mn2.uz()),
        st0.ux() + s1o6 * (mn0.ux() + MGMath::TWO * mn1.ux() + MGMath::TWO * mn2.ux() + mn3.ux()),
        st0.uy() + s1o6 * (mn0.uy() + MGMath::TWO * mn1.uy() + MGMath::TWO * mn2.uy() + mn3.uy()),
        st0.uz() + s1o6 * (mn0.uz() + MGMath::TWO * mn1.uz() + MGMath::TWO * mn2.uz() + mn3.uz())
    );
    if (withEloss) {
        Double_t eta    = st0.eta() * (MGMath::ONE + s1o6_ps * (mn0.e() + MGMath::TWO * mn1.e() + MGMath::TWO * mn2.e() + mn3.e()));
        Bool_t   is_mch = (MGNumc::Compare(eta) == eta_sign);
        if (is_mch) part.set_eta(eta);
        else        return false;
    }
    
    ppcal.push(part, mp0, mn0.orth().tau(), mn0.orth().rho());
    
    
    if (withJb) {
        TransferFunc tf0(st0, &mp0);
        TransferFunc tf1(st1, &mp1);
        TransferFunc tf2(st2, &mp2);
        TransferFunc tf3(st3, &mp3);

        PhyJb jb1;
        
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

        if (withEloss) jb1.gg(JEA, JEA) += s1o2_ps * tf0.ee();
        
        TransferPhyJb tj1(tf1, jb1);
        
        PhyJb jb2;
        
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

        if (withEloss) jb2.gg(JEA, JEA) += s1o2_ps * tj1.ee();
        
        TransferPhyJb tj2(tf2, jb2);
        
        PhyJb jb3;
        
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

        if (withEloss) jb3.gg(JEA, JEA) += step_ps * tj2.ee();
        
        TransferPhyJb tj3(tf3, jb3);
        
        phyJb->init();
        
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

        if (withEloss) phyJb->gg(JEA, JEA) += s1o6_ps * (tf0.ee() + MGMath::TWO * tj1.ee() + MGMath::TWO * tj2.ee() + tj3.ee());
    }
    
    return true;
}


} // namespace TrackSys


#endif // __TRACKLibs_Prop_C__
