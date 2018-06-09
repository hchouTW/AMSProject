#ifndef __TRACKLibs_HitSt_C__
#define __TRACKLibs_HitSt_C__


namespace TrackSys {

// VirtualHitSt
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

void HitStTRK::cal(const PhySt& part) {
    if (!set_type(part.info())) return;

    nrmc_.fill(Numc::ZERO<>);
    divc_.fill(Numc::ZERO<>);
    SVecD<3>&& crs = (coo_ - part.c());
    if (side_coo_(0)) {
        erc_(0)  = pdf_cx_->efft_sgm(crs(0));
        nrmc_[0] = crs(0) / erc_(0);
        divc_[0] = Numc::NEG<> / erc_(0);
    }
    if (side_coo_(1)) {
        erc_(1) = pdf_cy_->efft_sgm(crs(1));
        nrmc_[1] = crs(1) / erc_(1);
        divc_[1] = Numc::NEG<> / erc_(1);
    }

    nrmq_.fill(Numc::ZERO<>);
    divq_.fill(Numc::ZERO<>);
    if (side_q_[0]) {
        SVecD<2>&& ionx = (*pdf_qx_)(q_[0]*q_[0], part.igmbta());
        nrmq_[0] = ionx(0);
        divq_[0] = ionx(1) * part.mu() * part.eta_sign();
        divq_[1] = ionx(1) * part.gm();
    }
    if (side_q_[1]) {
        SVecD<2>&& iony = (*pdf_qy_)(q_[1]*q_[1], part.igmbta());
        nrmq_[1] = iony(0);
        divq_[2] = iony(1) * part.mu() * part.eta_sign();
        divq_[3] = iony(1) * part.gm();
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

//IonEloss HitStTRK::PDF_Q01_QX_( // TEMP
//    { 5.96811e-01, 2.03022e+00, -2.63964e+00 }, // Kpa
//    { 2.04435e-01, 4.72938e+00, 1.74079e+00, 4.02238e-02, 1.71754e+00 }, // Mpv
//    { 2.82318e-01, 2.51241e+00, 9.71154e-01, 1.82888e+00, 4.67111e-01 }, // Sgm
//    0.266246 // Fluc
//);
//
//IonEloss HitStTRK::PDF_Q01_QY_( // TEMP
//    { 3.17864e+04,1.27740e-05, -2.30014e+00 }, // Kpa
//    { 1.61937e-01, 5.76721e+00, 1.86636e+00, 2.56620e-02, 1.51650e+00 }, // Mpv
//    { 1.39973e-01, 3.00167e+00, 1.24903e+00, 5.41780e-01, 4.81160e-01 }, // Sgm
//    0.152518 // Fluc
//);


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

void HitStTOF::cal(const PhySt& part) {
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
            Double_t divt = (Numc::NEG<> * ds / ter);
            divt_[0] = divt * (part.bta() * part.eta()) * (part.mu() * part.mu());
            divt_[1] = divt;
            if (!Numc::Valid(nrmt_) || !Numc::Valid(divt_sft_) || !Numc::Valid(divt)) {
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
        divq_[0] = ion(1) * part.mu() * part.eta_sign();
        divq_[1] = ion(1) * part.gm();
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

//IonEloss HitStTOF::PDF_Q01_Q_( // NEW2
//    { 8.56923e+00, 2.03064e-01, -2.46093e+00 }, // Kpa
//    { 3.19209e-01, 3.78021e+00, 1.03394e+00, 6.42416e-01, 1.27161e+00 }, // Mpv
//    { 5.62652e-02, 1.84529e+00, 1.41309e+00, 7.72896e-01, 1.52086e+00 }, // Sgm
//    0.0829427 // Fluc
//);


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

void HitStRICH::cal(const PhySt& part) {
    if (!set_type(part.info())) return;

    nrmib_ = Numc::ZERO<>;
    divib_.fill(Numc::ZERO<>);
    if (side_ib_) {
        // 1/bta := (1+igmbta*igmbta)^(1/2)
        Double_t dib = ib_ - part.ibta();
        Double_t sgm = pdf_ib_->efft_sgm(dib);
        nrmib_ = (dib / sgm);
        Double_t divib = (Numc::NEG<> / sgm);
        divib_[0] = divib * (part.bta() * part.eta()) * (part.mu() * part.mu());
        divib_[1] = divib;
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
    4.50204264345890337e-01, 9.12779e-04,
    5.49795735654109663e-01, 1.51231e-03
);

MultiGaus HitStRICH::PDF_NAF_Q01_IB_(
    MultiGaus::Opt::ROBUST,
    9.28354466225493113e-01, 3.24032e-03,
    7.16455337745068865e-02, 6.91529e-03
);

/*
void HitStTRD::clear() {
    seqIDcx_ = -1;
    seqIDcy_ = -1;
    seqIDel_ = -1;
    
    side_el_ = false;
    ndofel_  = Numc::ZERO<Short_t>;
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

void HitStTRD::cal(const PhySt& part) {
    if (!set_type(part.info())) return;
    
    nrmel_ = Numc::ZERO<>;
    divel_.fill(Numc::ZERO<>);
    if (side_el_) {
        SVecD<2>&& gmion = (*pdf_el_)(el_, part.igmbta());
        nrmel_ = gmion(0);
        divel_[0] = gmion(1) * part.mu() * part.eta_sign();
        divel_[1] = gmion(1) * part.gm();
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
    { 3.42286e-01, 1.49108e+00, 3.79622e-01, 6.51468e-01, 1.45158e-01, 1.29064e-01, 1.16096e+00, 8.47144e+00 }, // Kpa
    { 1.13611e-01, 1.13937e+01, 1.23836e+00, 2.67476e-06, 2.51050e+00, 2.14009e+00, 8.59030e-01, 6.26050e+00 }, // Mpv
    { 7.30837e-02, 2.85852e+00, 1.18204e+00, 2.29138e-02, 5.94741e-01, 5.14559e-01, 9.92186e-01, 7.07455e+00 }, // Sgm
    0.170 // Fluc
);
*/

} // namesapce TrackSys


#endif // __TRACKLibs_HitSt_C__
