#ifndef __TRACKLibs_HitSt_C__
#define __TRACKLibs_HitSt_C__


namespace TrackSys {

// VirtualHitSt
Double_t VirtualHitSt::DoNoiseControllerLU(Double_t norm, Double_t threshold) {
    Double_t thres  = ((Numc::Compare(threshold) < 0) ? NOISE_THRESHOLD_DEFAULT : threshold);
    Double_t sqrnrm = (norm * norm);
    Double_t sclnrm = std::log1p(sqrnrm * sqrnrm);
    Double_t controller = Numc::HALF * (Numc::ONE<> + std::erf(thres - sclnrm));
    if (!Numc::Valid(controller)) controller = Numc::ONE<>;
    return controller;
}


Double_t VirtualHitSt::DoNoiseControllerL(Double_t norm, Double_t threshold) {
    Double_t thres  = ((Numc::Compare(threshold) < 0) ? NOISE_THRESHOLD_DEFAULT : threshold);
    Short_t  sign   = Numc::Compare(norm);
    Double_t sqrnrm = (norm * norm);
    Double_t sclnrm = sign * std::log1p(sqrnrm * sqrnrm);
    Double_t controller = Numc::HALF * (Numc::ONE<> + std::erf(thres + sclnrm));
    if (!Numc::Valid(controller)) controller = Numc::ONE<>;
    return controller;
}


Double_t VirtualHitSt::DoNoiseControllerU(Double_t norm, Double_t threshold) {
    Double_t thres  = ((Numc::Compare(threshold) < 0) ? NOISE_THRESHOLD_DEFAULT : threshold);
    Short_t  sign   = Numc::Compare(norm);
    Double_t sqrnrm = (norm * norm);
    Double_t sclnrm = sign * std::log1p(sqrnrm * sqrnrm);
    Double_t controller = Numc::HALF * (Numc::ONE<> + std::erf(thres - sclnrm));
    if (!Numc::Valid(controller)) controller = Numc::ONE<>;
    return controller;
}


Double_t VirtualHitSt::DoNoiseSlowControllerLU(Double_t norm, Double_t threshold) {
    Double_t thres  = ((Numc::Compare(threshold) < 0) ? NOISE_THRESHOLD_DEFAULT : threshold);
    Double_t sclnrm = std::log1p(norm * norm);
    Double_t controller = Numc::HALF * (Numc::ONE<> + std::erf(thres - sclnrm));
    if (!Numc::Valid(controller)) controller = Numc::ONE<>;
    return controller;
}


Double_t VirtualHitSt::DoNoiseSlowControllerL(Double_t norm, Double_t threshold) {
    Double_t thres  = ((Numc::Compare(threshold) < 0) ? NOISE_THRESHOLD_DEFAULT : threshold);
    Short_t  sign   = Numc::Compare(norm);
    Double_t sclnrm = sign * std::log1p(norm * norm);
    Double_t controller = Numc::HALF * (Numc::ONE<> + std::erf(thres + sclnrm));
    if (!Numc::Valid(controller)) controller = Numc::ONE<>;
    return controller;
}


Double_t VirtualHitSt::DoNoiseSlowControllerU(Double_t norm, Double_t threshold) {
    Double_t thres  = ((Numc::Compare(threshold) < 0) ? NOISE_THRESHOLD_DEFAULT : threshold);
    Short_t  sign   = Numc::Compare(norm);
    Double_t sclnrm = sign * std::log1p(norm * norm);
    Double_t controller = Numc::HALF * (Numc::ONE<> + std::erf(thres - sclnrm));
    if (!Numc::Valid(controller)) controller = Numc::ONE<>;
    return controller;
}


VirtualHitSt::VirtualHitSt(Detector dec, Short_t lay, Bool_t scx, Bool_t scy, Bool_t scz) {
    clear();
    dec_ = dec;
    lay_ = lay;
    side_coo_ = std::move(SVecO<3>(scx, scy, scz));
}
    
void VirtualHitSt::clear() {
    seqID_   = -1;
    seqIDcx_ = -1;
    seqIDcy_ = -1;

    type_ = PartType::Proton;
    dec_  = Detector::NONE;
    lay_  = 0;
    
    side_coo_ = std::move(SVecO<3>());
    coo_      = std::move(SVecD<3>());
    erc_      = std::move(SVecD<2>(Numc::ONE<>, Numc::ONE<>));
    
    nrmc_.fill(Numc::ZERO<>);
    divc_.fill(Numc::ZERO<>);
}
        

// HitStTRK
void HitStTRK::clear() {
    seqIDcx_ = -1;
    seqIDcy_ = -1;
    seqIDqx_ = -1;
    seqIDqy_ = -1;

    nsr_.fill(Numc::ZERO<Short_t>);

    side_q_.fill(false);
    q_.fill(Numc::ZERO<>);
    
    nrmq_.fill(Numc::ZERO<>);
    divq_.fill(Numc::ZERO<>);

    pdf_cx_ = nullptr;
    pdf_cy_ = nullptr;
    pdf_qx_ = nullptr;
    pdf_qy_ = nullptr;

    set_type();
    if (pdf_cx_ != nullptr && pdf_cy_ != nullptr)
        erc_ = std::move(SVecD<2>(pdf_cx_->efft_sgm(), pdf_cy_->efft_sgm()));
}

Short_t HitStTRK::set_seqID(Short_t seqID) {
    if (seqID < 0) { seqID_ = -1; return 0; }

    Short_t iter = 0;
    if (side_coo_(0)) { seqIDcx_ = seqID + iter; iter++; } else seqIDcx_ = -1;
    if (side_coo_(1)) { seqIDcy_ = seqID + iter; iter++; } else seqIDcy_ = -1;
    if (side_q_[0])   { seqIDqx_ = seqID + iter; iter++; } else seqIDqx_ = -1;
    if (side_q_[1])   { seqIDqy_ = seqID + iter; iter++; } else seqIDqy_ = -1;
    if (iter != 0) seqID_ = seqID; else seqID_ = -1;
    return iter;
}

void HitStTRK::cal(const PhySt& part, const NoiseController& ctler) {
    if (!set_type(part.info())) return;

    nrmc_.fill(Numc::ZERO<>);
    divc_.fill(Numc::ZERO<>);
    SVecD<3>&& crs = (coo_ - part.c());
    if (side_coo_(0)) {
        erc_(0)  = pdf_cx_->efft_sgm(crs(0));
        nrmc_[0] = (crs(0) / erc_(0));
        divc_[0] = (Numc::NEG<> / erc_(0));
        if (NoiseController::ON == ctler)
            divc_[0] *= DoNoiseSlowControllerLU(nrmc_[0], NOISE_THRESHOLD_COORD);
        if (!Numc::Valid(nrmc_[0]) || !Numc::Valid(divc_[0])) {
            nrmc_[0] = Numc::ZERO<>;
            divc_[0] = Numc::ZERO<>;
        }
    }
    if (side_coo_(1)) {
        erc_(1) = pdf_cy_->efft_sgm(crs(1));
        nrmc_[1] = (crs(1) / erc_(1));
        divc_[1] = (Numc::NEG<> / erc_(1));
        if (NoiseController::ON == ctler)
            divc_[1] *= DoNoiseSlowControllerLU(nrmc_[1], NOISE_THRESHOLD_COORD); 
        if (!Numc::Valid(nrmc_[1]) || !Numc::Valid(divc_[1])) {
            nrmc_[1] = Numc::ZERO<>;
            divc_[1] = Numc::ZERO<>;
        }
    }

    nrmq_.fill(Numc::ZERO<>);
    divq_.fill(Numc::ZERO<>);
    if (side_q_[0]) {
        SVecD<2>&& ionx = (*pdf_qx_)(q_[0]*q_[0], part.igmbta());
        nrmq_[0] = ionx(0);
        divq_[0] = ionx(1) * (part.mu() * part.eta_sign());
        divq_[1] = ionx(1);
        
        Double_t ctlerq = DoNoiseControllerL(nrmq_[0], NOISE_THRESHOLD_DEDX);
        divq_[0] *= ctlerq;
        divq_[1] *= ctlerq;
        
        if (!Numc::Valid(nrmq_[0]) || !Numc::Valid(divq_[0]) || !Numc::Valid(divq_[1])) {
            nrmq_[0] = Numc::ZERO<>;
            divq_[0] = Numc::ZERO<>;
            divq_[1] = Numc::ZERO<>;
        }
    }
    if (side_q_[1]) {
        SVecD<2>&& iony = (*pdf_qy_)(q_[1]*q_[1], part.igmbta());
        nrmq_[1] = iony(0);
        divq_[2] = iony(1) * (part.mu() * part.eta_sign());
        divq_[3] = iony(1);
        
        Double_t ctlerq = DoNoiseControllerL(nrmq_[1], NOISE_THRESHOLD_DEDX);
        divq_[2] *= ctlerq;
        divq_[3] *= ctlerq;
        
        if (!Numc::Valid(nrmq_[1]) || !Numc::Valid(divq_[2]) || !Numc::Valid(divq_[3])) {
            nrmq_[1] = Numc::ZERO<>;
            divq_[2] = Numc::ZERO<>;
            divq_[3] = Numc::ZERO<>;
        }
    }

    set_dummy_x(part.cx());
    set_dummy_y(part.cy());
}

Bool_t HitStTRK::set_type(const PartInfo& info) {
    if ((info.is_std() && type_ == info.type()) && (pdf_cx_ && pdf_cy_ && pdf_qx_ && pdf_qy_)) return true;

    switch (info.chrg()) {
        case 1 : case -1 :
        {
            pdf_cx_ = &PDF_Q01_CX_NN_;
            if      (nsr_[0] == 1) pdf_cx_ = &PDF_Q01_CX_N1_;
            else if (nsr_[0] == 2) pdf_cx_ = &PDF_Q01_CX_N2_;
            else if (nsr_[0] >= 3) pdf_cx_ = &PDF_Q01_CX_N3_;
            pdf_cy_ = &PDF_Q01_CY_NN_;
            if      (nsr_[1] == 1) pdf_cy_ = &PDF_Q01_CY_N1_;
            else if (nsr_[1] == 2) pdf_cy_ = &PDF_Q01_CY_N2_;
            else if (nsr_[1] == 3) pdf_cy_ = &PDF_Q01_CY_N3_;
            else if (nsr_[1] >= 4) pdf_cy_ = &PDF_Q01_CY_N4_;
            pdf_qx_ = &PDF_Q01_QX_;
            pdf_qy_ = &PDF_Q01_QY_;
            type_ = info.type();
            break;
        }
        default :
            CERR("HitStTRK::set_type() NO PartType Setting.\n");
            return false;
    }
    return true;
}

MultiGaus HitStTRK::PDF_Q01_CX_NN_(
    MultiGaus::Opt::ROBUST,
    3.31376712664997630e-01, 1.77875e-03,
    5.10255425401639595e-01, 2.65271e-03,
    1.58367861933362775e-01, 5.15837e-03
);

MultiGaus HitStTRK::PDF_Q01_CX_N1_(
    MultiGaus::Opt::ROBUST,
    2.75608e-03
);

MultiGaus HitStTRK::PDF_Q01_CX_N2_(
    MultiGaus::Opt::ROBUST,
    6.02616961110044480e-01, 1.82559e-03,
    3.97383038889955520e-01, 4.18164e-03
);

MultiGaus HitStTRK::PDF_Q01_CX_N3_(
    MultiGaus::Opt::ROBUST,
    6.75185841999877190e-01, 1.84066e-03,
    3.24814158000122755e-01, 4.90605e-03
);

MultiGaus HitStTRK::PDF_Q01_CY_NN_(
    MultiGaus::Opt::ROBUST,
    5.50759994181610257e-01, 7.90750e-04,
    3.74341189078839287e-01, 1.52916e-03,
    7.48988167395504278e-02, 3.63939e-03
);

MultiGaus HitStTRK::PDF_Q01_CY_N1_(
    MultiGaus::Opt::ROBUST,
    5.45247134760751928e-01, 9.87256e-04,
    4.54752865239248016e-01, 1.65822e-03
);

MultiGaus HitStTRK::PDF_Q01_CY_N2_(
    MultiGaus::Opt::ROBUST,
    4.53177772355239150e-01, 7.32341e-04,
    4.29177717623073274e-01, 1.19994e-03,
    1.17644510021687618e-01, 2.17922e-03
);

MultiGaus HitStTRK::PDF_Q01_CY_N3_(
    MultiGaus::Opt::ROBUST,
    5.08190182558828307e-01, 7.74660e-04,
    3.39669567581713017e-01, 1.56378e-03,
    1.52140249859458648e-01, 3.18597e-03
);

MultiGaus HitStTRK::PDF_Q01_CY_N4_(
    MultiGaus::Opt::ROBUST,
    5.19847432537280163e-01, 7.84945e-04,
    3.39808551148078952e-01, 1.72969e-03,
    1.40344016314640940e-01, 3.85202e-03
);

IonEloss HitStTRK::PDF_Q01_QX_(
    { 9.97328e-01, 1.26479e+00, -1.33731e+01, 2.21958e+03 }, // Kpa
    { 5.35798e-02, 1.31568e+01, -6.12405e+00, 1.89315e+00, 1.12457e-03, 2.64436e+00 }, // Mpv
    { 2.49410e-03, 6.80308e+01, -5.45630e+00, 1.54873e+00, 2.68543e-10, 5.02291e+00 }, // Sgm
    0.266246 // Fluc
);

IonEloss HitStTRK::PDF_Q01_QY_(
    { 9.96926e-01, 1.33625e+00, -1.32089e+01, 1.23685e+03 }, // Kpa
    { 7.07719e-02, 1.09349e+01, -2.19142e+00, 1.96544e+00, 1.79601e-03, 1.92512e+00 }, // Mpv
    { 1.94316e-03, 1.03082e+02, -3.56452e+01, 1.62118e+00, 3.08950e-15, 8.63832e+00 }, // Sgm
    0.152518 // Fluc
);


// HitStTOF
void HitStTOF::clear() {
    seqIDcx_ = -1;
    seqIDcy_ = -1;
    seqIDt_  = -1;
    seqIDq_  = -1;
    
    side_t_ = false;
    orgt_   = Numc::ZERO<>;
    sftt_   = Numc::ZERO<>;

    side_q_ = false;
    q_      = Numc::ZERO<>;
    
    nrmt_     = Numc::ZERO<>;
    divt_sft_ = Numc::ZERO<>;
    divt_.fill(Numc::ZERO<>);
    
    nrmq_ = Numc::ZERO<>;
    divq_.fill(Numc::ZERO<>);

    pdf_t_ = nullptr;
    pdf_q_ = nullptr;

    set_type();
}

Short_t HitStTOF::set_seqID(Short_t seqID) {
    if (seqID < 0) { seqID_ = -1; return 0; }

    Short_t iter = 0;
    if (side_coo_(0)) { seqIDcx_ = seqID + iter; iter++; } else seqIDcx_ = -1;
    if (side_coo_(1)) { seqIDcy_ = seqID + iter; iter++; } else seqIDcy_ = -1;
    if (side_t_     ) { seqIDt_  = seqID + iter; iter++; } else seqIDt_  = -1;
    if (side_q_     ) { seqIDq_  = seqID + iter; iter++; } else seqIDq_  = -1;
    if (iter != 0) seqID_ = seqID; else seqID_ = -1;
    return iter;
}

void HitStTOF::cal(const PhySt& part, const NoiseController& ctler) {
    if (!set_type(part.info())) return;

    nrmt_     = Numc::ZERO<>;
    divt_sft_ = Numc::ZERO<>;
    divt_.fill(Numc::ZERO<>);
    if (side_t_) {
        // 1/bta := (1+igmbta*igmbta)^(1/2)
        // t     := (delta_S / bta)
        // res_t := (meas_t - part_t)
        sftt_ = orgt_ + OFFSET_T_; // update shift time
        Double_t ds = std::fabs(part.path() - OFFSET_S_);
        Double_t dt = sftt_ - part.time();
        if (Numc::Compare(ds) >= 0) {
            Double_t sgm = ( 3.40971e+00 - 1.32705e+00 * std::erf( std::log(part.igmbta()) - 5.78158e-02 ) ); // ONE Layer [cm]
            if (!TShiftCorr_) sgm *= Numc::SQRT_TWO; // TWO Layer [cm]
            Double_t ter = MultiGaus::RobustSgm(dt, sgm, pdf_t_->opt());
            nrmt_     = (dt / ter);
            divt_sft_ = (Numc::ONE<> / ter);
            Double_t divt = (Numc::NEG<> * ds / ter) * (part.bta() * part.igmbta());
            divt_[0] = divt * (part.mu() * part.eta_sign());
            divt_[1] = divt;
            if (NoiseController::ON == ctler) {
                Double_t ctlert = DoNoiseSlowControllerLU(nrmt_, NOISE_THRESHOLD_TIME);
                divt_sft_ *= ctlert;
                divt_[0] *= ctlert;
                divt_[1] *= ctlert;
            }
            if (!Numc::Valid(nrmt_) || !Numc::Valid(divt_sft_) || !Numc::Valid(divt_[0]) || !Numc::Valid(divt_[1])) {
                nrmt_     = Numc::ZERO<>;
                divt_sft_ = Numc::ZERO<>;
                divt_.fill(Numc::ZERO<>);
            }
        }
    }
    
    nrmq_ = Numc::ZERO<>;
    divq_.fill(Numc::ZERO<>);
    if (side_q_) {
        SVecD<2>&& ion = (*pdf_q_)(q_*q_, part.igmbta());
        nrmq_ = ion(0);
        divq_[0] = ion(1) * (part.mu() * part.eta_sign());
        divq_[1] = ion(1);
        
        Double_t ctlerq = DoNoiseControllerL(nrmq_, NOISE_THRESHOLD_DEDX);
        divq_[0] *= ctlerq;
        divq_[1] *= ctlerq;
        
        if (!Numc::Valid(nrmq_) || !Numc::Valid(divq_[0]) || !Numc::Valid(divq_[1])) {
            nrmq_ = Numc::ZERO<>;
            divq_.fill(Numc::ZERO<>);
        }
    }

    set_dummy_x(part.cx());
    set_dummy_y(part.cy());
}

Bool_t HitStTOF::set_type(const PartInfo& info) {
    if ((info.is_std() && type_ == info.type()) && (pdf_t_ && pdf_q_)) return true;
    
    switch (info.chrg()) {
        case 1 : case -1 :
        {
            pdf_t_ = &PDF_Q01_T_;
            pdf_q_ = &PDF_Q01_Q_;
            type_ = info.type();
            break;
        }
        default :
            CERR("HitStTOF::set_type() NO PartType Setting.\n");
            return false;
    }
    return true;
}

MultiGaus HitStTOF::PDF_Q01_T_(
    MultiGaus::Opt::ROBUST,
    6.698790e+00 // TWO Time Fluc
);

IonEloss HitStTOF::PDF_Q01_Q_(
    { 1.00000e+00, 1.23333e+00, -2.31113e+00, 1.97051e+00 }, // Kpa
    { 4.90401e-01, 2.21732e+00, 4.49282e-03, 1.33820e+00, 1.11419e+00, 1.21069e+00 }, // Mpv
    { 1.21818e+00, 9.31334e+00, 9.55853e+00, 1.42066e-01, 7.40717e-01, 1.98209e+00 }, // Sgm
    0.0829427 // Fluc
);


// HitStRICH
void HitStRICH::clear() {
    seqIDcx_ = -1;
    seqIDcy_ = -1;
    seqIDib_ = -1;
   
    rad_ = Radiator::AGL;

    side_ib_ = false;
    ib_      = Numc::ZERO<>;

    nrmib_ = Numc::ZERO<>;
    divib_.fill(Numc::ZERO<>);

    pdf_ib_ = nullptr;

    set_type();
}

Short_t HitStRICH::set_seqID(Short_t seqID) {
    if (seqID < 0) { seqID_ = -1; return 0; }

    Short_t iter = 0;
    if (side_coo_(0)) { seqIDcx_ = seqID + iter; iter++; } else seqIDcx_ = -1;
    if (side_coo_(1)) { seqIDcy_ = seqID + iter; iter++; } else seqIDcy_ = -1;
    if (side_ib_    ) { seqIDib_ = seqID + iter; iter++; } else seqIDib_ = -1;
    if (iter != 0) seqID_ = seqID; else seqID_ = -1;
    return iter;
}

void HitStRICH::cal(const PhySt& part, const NoiseController& ctler) {
    if (!set_type(part.info())) return;

    nrmib_ = Numc::ZERO<>;
    divib_.fill(Numc::ZERO<>);
    if (side_ib_) {
        // 1/bta := (1+igmbta*igmbta)^(1/2)
        Double_t dib = ib_ - part.ibta();
        Double_t sgm = pdf_ib_->efft_sgm(dib);
        nrmib_ = (dib / sgm);
        Double_t divib = (Numc::NEG<> / sgm) * (part.bta() * part.igmbta());
        divib_[0] = divib * (part.mu() * part.eta_sign());
        divib_[1] = divib;
        if (NoiseController::ON == ctler) {
            Double_t ctlerib = DoNoiseSlowControllerLU(nrmib_, NOISE_THRESHOLD_BETA);
            divib_[0] *= ctlerib;
            divib_[1] *= ctlerib;
        }
        if (!Numc::Valid(nrmib_) || !Numc::Valid(divib)) {
            nrmib_ = Numc::ZERO<>;
            divib_.fill(Numc::ZERO<>);
        }
    }
    
    set_dummy_x(part.cx());
    set_dummy_y(part.cy());
}

Bool_t HitStRICH::set_type(const PartInfo& info) {
    if ((info.is_std() && type_ == info.type()) && (pdf_ib_)) return true;
    
    switch (info.chrg()) {
        case 1 : case -1 :
        {
            switch (rad_) {
                case Radiator::AGL : pdf_ib_ = &PDF_AGL_Q01_IB_; break;
                case Radiator::NAF : pdf_ib_ = &PDF_NAF_Q01_IB_; break;
            }
            type_ = info.type();
            break;
        }
        default :
            CERR("HitStRICH::set_type() NO PartType Setting.\n");
            return false;
    }
    return true;
}

MultiGaus HitStRICH::PDF_AGL_Q01_IB_(
    MultiGaus::Opt::ROBUST,
    4.64372169946583258e-01, 9.46452e-04,
    5.35627830053416742e-01, 1.55197e-03

);

MultiGaus HitStRICH::PDF_NAF_Q01_IB_(
    MultiGaus::Opt::ROBUST,
    8.30132503718388537e-01, 3.20352e-03,
    1.69867496281611519e-01, 5.35715e-03
);


void HitStTRD::clear() {
    seqIDcx_ = -1;
    seqIDcy_ = -1;
    seqIDel_ = -1;
    
    side_el_ = false;
    el_      = Numc::ZERO<>;

    nrmel_ = Numc::ZERO<>;
    divel_.fill(Numc::ZERO<>);

    pdf_el_ = nullptr;

    set_type();
}

Short_t HitStTRD::set_seqID(Short_t seqID) {
    if (seqID < 0) { seqID_ = -1; return 0; }

    Short_t iter = 0;
    if (side_coo_(0)) { seqIDcx_ = seqID + iter; iter++; } else seqIDcx_ = -1;
    if (side_coo_(1)) { seqIDcy_ = seqID + iter; iter++; } else seqIDcy_ = -1;
    if (side_el_    ) { seqIDel_ = seqID + iter; iter++; } else seqIDel_ = -1;
    if (iter != 0) seqID_ = seqID; else seqID_ = -1;
    return iter;
}

void HitStTRD::cal(const PhySt& part, const NoiseController& ctler) {
    if (!set_type(part.info())) return;
    
    nrmel_ = Numc::ZERO<>;
    divel_.fill(Numc::ZERO<>);
    if (side_el_) {
        SVecD<3>&& lgge = (*pdf_el_)(el_, part.igmbta());
        nrmel_ = lgge(0);
        divel_[0] = lgge(1) * part.mu() * part.eta_sign();
        divel_[1] = lgge(1);
        if (NoiseController::ON == ctler) {
            Double_t ctlerel = DoNoiseControllerL(nrmel_, NOISE_THRESHOLD_DEDX);
            divel_[0] *= ctlerel;
            divel_[1] *= ctlerel;
        }
        if (!Numc::Valid(nrmel_) || !Numc::Valid(divel_[0]) || !Numc::Valid(divel_[1])) {
            nrmel_ = Numc::ZERO<>;
            divel_.fill(Numc::ZERO<>);
        }
    }

    set_dummy_x(part.cx());
    set_dummy_y(part.cy());
}

Bool_t HitStTRD::set_type(const PartInfo& info) {
    if ((info.is_std() && type_ == info.type()) && (pdf_el_)) return true;
    
    switch (info.chrg()) {
        case 1 : case -1 :
        {
            pdf_el_ = &PDF_Q01_EL_;
            type_ = info.type();
            break;
        }
        default :
            CERR("HitStTRD::set_type() NO PartType Setting.\n");
            return false;
    }
    return true;
}

GmIonEloss HitStTRD::PDF_Q01_EL_(
    { 6.87250e-01, 8.49158e-01, -6.50955e+00 }, // Ratio
    { 2.17486e-04, 4.27812e+00 }, // Kpa
    { 9.62033e+00, 1.26248e+00, 8.36500e-03, 1.24422e+00, 2.22255e+00, 8.06270e-02 }, // Mpv
    { 2.13993e-01, 2.78925e+00, 1.65428e+00, 6.55129e-01, 1.17441e-01, 3.61303e-01 }, // Sgm
    { 2.11427 }, // Alpha
    { 1.83481e-01, 5.07634e-01, -5.30224e+00 }, // Beta
    { 6.52658434 }, // Erf Man
    { 1.79160092 }  // Erf Sgm
);


} // namesapce TrackSys


#endif // __TRACKLibs_HitSt_C__
