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
    
    nrmc_ = std::move(SVecD<2>());
    divc_ = std::move(SVecD<2>());
}
        

// HitStTRK
void HitStTRK::clear() {
    seqIDcx_ = -1;
    seqIDcy_ = -1;
    seqIDqx_ = -1;
    seqIDqy_ = -1;

    nsr_    = std::move(SVecS<2>());
    side_q_ = std::move(SVecO<2>());
    q_      = std::move(SVecD<2>());
    
    nrmq_ = std::move(SVecD<2>());
    divq_ = std::move(SVecD<2>());

    set_type();
    if (pdf_cx_ != nullptr && pdf_cy_ != nullptr)
        erc_ = std::move(SVecD<2>(pdf_cx_->efft_sgm(), pdf_cy_->efft_sgm()));
}

Short_t HitStTRK::set_seqID(Short_t seqID) {
    if (seqID < 0) { seqID_ = -1; return 0; }

    Short_t iter = 0;
    if (side_coo_(0)) { seqIDcx_ = seqID + iter; iter++; } else seqIDcx_ = -1;
    if (side_coo_(1)) { seqIDcy_ = seqID + iter; iter++; } else seqIDcy_ = -1;
    if (side_q_(0))   { seqIDqx_ = seqID + iter; iter++; } else seqIDqx_ = -1;
    if (side_q_(1))   { seqIDqy_ = seqID + iter; iter++; } else seqIDqy_ = -1;
    if (iter != 0) seqID_ = seqID; else seqID_ = -1;
    return iter;
}

void HitStTRK::cal(const PhySt& part) {
    set_type(part.info());

    nrmc_ = std::move(SVecD<2>());
    divc_ = std::move(SVecD<2>());
    SVecD<3>&& crs = (coo_ - part.c());
    if (side_coo_(0)) {
        erc_(0)  = pdf_cx_->efft_sgm(crs(0));
        nrmc_(0) = crs(0) / erc_(0);
        divc_(0) = Numc::NEG<> / erc_(0);
    }
    if (side_coo_(1)) {
        erc_(1) = pdf_cy_->efft_sgm(crs(1));
        nrmc_(1) = crs(1) / erc_(1);
        divc_(1) = Numc::NEG<> / erc_(1);
    }

    nrmq_ = std::move(SVecD<2>());
    divq_ = std::move(SVecD<2>());
    if (side_q_(0)) {
        SVecD<2>&& ionx = (*pdf_qx_)(q_(0)*q_(0), part.eta());
        nrmq_(0) = ionx(0);
        divq_(0) = ionx(1);
    }
    if (side_q_(1)) {
        SVecD<2>&& iony = (*pdf_qy_)(q_(1)*q_(1), part.eta());
        nrmq_(1) = iony(0);
        divq_(1) = iony(1);
    }

    set_dummy_x(part.cx());
    set_dummy_y(part.cy());
}

void HitStTRK::set_type(const PartInfo& info) {
    if ((info.is_std() && type_ == info.type()) && (pdf_cx_ && pdf_cy_ && pdf_qx_ && pdf_qy_)) return;

    switch (info.chrg()) {
        case 1 : case -1 :
        {
            pdf_cx_ = &PDF_Q01_CX_NN_;
            if      (nsr_(0) == 1) pdf_cx_ = &PDF_Q01_CX_N1_;
            else if (nsr_(0) == 2) pdf_cx_ = &PDF_Q01_CX_N2_;
            else if (nsr_(0) >= 3) pdf_cx_ = &PDF_Q01_CX_N3_;
            pdf_cy_ = &PDF_Q01_CY_NN_;
            if      (nsr_(1) == 1) pdf_cy_ = &PDF_Q01_CY_N1_;
            else if (nsr_(1) == 2) pdf_cy_ = &PDF_Q01_CY_N2_;
            else if (nsr_(1) == 3) pdf_cy_ = &PDF_Q01_CY_N3_;
            else if (nsr_(1) >= 4) pdf_cy_ = &PDF_Q01_CY_N4_;
            pdf_qx_ = &PDF_Q01_QX_;
            pdf_qy_ = &PDF_Q01_QY_;
            type_ = info.type();
            break;
        }
        default :
            CERR("HitStTRK::set_type() NO PartType Setting.\n");
            break;
    }
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

//IonEloss HitStTRK::PDF_Q01_QX_( // OLD
//    { 6.05261e-03, 2.12431e-01, 5.25164e-01, 2.81963e+00, 2.40051e+00, 4.40867e+00 }, // Kpa
//    { 1.06985e-01, 9.37601e+00, 8.91247e-01, 6.73657e-02, 1.83183e+00 }, // Mpv
//    { 1.34009e-05, 4.93244e+03, 7.61516e-01, 0.00000e+00, 3.29793e+01 }, // Sgm
//    0.178087 // Fluc
//);
//
//IonEloss HitStTRK::PDF_Q01_QY_( // OLD
//    { 8.52848e-03, 3.24435e+00, 2.06864e-01, 1.60633e+00, 2.09756e+00, 1.53539e+00 }, // Kpa
//    { 1.40446e-01, 7.48716e+00, 9.13645e-01, 1.31821e-01, 1.07348e+00 }, // Mpv
//    { 2.46606e-01, 2.89505e+00, 2.79787e-01, 3.49813e+00, 5.30240e-02 }, // Sgm
//    0.166633 // Fluc
//);

IonEloss HitStTRK::PDF_Q01_QX_( // TEMP
    { 1.42180e+01, 1.03536e-01, -2.12863e+00 }, // Kpa
    { 1.10206e-02, 8.62096e+01, 9.49522e-01, 3.59949e-04, 2.70646e+00 }, // Mpv
    { 1.89554e-01, 2.80196e+00, 6.28283e-01, 4.37731e+00, 1.25767e+00 }, // Sgm
    0.0829427 // Fluc
);

IonEloss HitStTRK::PDF_Q01_QY_( // TEMP
    { 1.42180e+01, 1.03536e-01, -2.12863e+00 }, // Kpa
    { 1.10206e-02, 8.62096e+01, 9.49522e-01, 3.59949e-04, 2.70646e+00 }, // Mpv
    { 1.89554e-01, 2.80196e+00, 6.28283e-01, 4.37731e+00, 1.25767e+00 }, // Sgm
    0.0829427 // Fluc
);


// HitStTOF
void HitStTOF::clear() {
    seqIDcx_ = -1;
    seqIDcx_ = -1;
    seqIDt_  = -1;
    seqIDq_  = -1;
    
    side_t_ = false;
    t_      = Numc::ZERO<>;

    side_q_ = false;
    q_      = Numc::ZERO<>;
    
    nrmt_ = Numc::ZERO<>;
    divt_ = Numc::ZERO<>;
    
    nrmq_ = Numc::ZERO<>;
    divq_ = Numc::ZERO<>;

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
    set_type(part.info());

    nrmt_ = Numc::ZERO<>;
    divt_ = Numc::ZERO<>;
    if (side_t_) {
        // t := (delta_S / beta)
        // res_t = (meas_t - part_t);
        Double_t ds = std::fabs(part.path() - OFFSET_S_);
        Double_t dt = (t_ + OFFSET_T_) - part.time();
        if (!Numc::EqualToZero(ds)) {
            Double_t sgm = ( 4.82206e+00 + 1.87673e+00 * std::erf( -(std::log(std::fabs(part.eta()))-5.78158e-02) ) ); // [cm]
            Double_t ter = MultiGaus::RobustSgm(dt, sgm, pdf_t_->opt());
            nrmt_ = (dt / ter);
            divt_ = (Numc::NEG<> * ds / ter) * (part.bta() * part.eta());
            if (!Numc::Valid(nrmt_) || !Numc::Valid(divt_)) {
                nrmt_ = Numc::ZERO<>;
                divt_ = Numc::ZERO<>;
            }
        }
    }
    
    nrmq_ = Numc::ZERO<>;
    divq_ = Numc::ZERO<>;
    if (side_q_) {
        SVecD<2>&& ion = (*pdf_q_)(q_*q_, part.eta());
        nrmq_ = ion(0);
        divq_ = ion(1);
    }

    set_dummy_x(part.cx());
    set_dummy_y(part.cy());
}

void HitStTOF::set_type(const PartInfo& info) {
    if ((info.is_std() && type_ == info.type()) && (pdf_t_ && pdf_q_)) return;
    
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
            break;
    }
}

MultiGaus HitStTOF::PDF_Q01_T_(
    MultiGaus::Opt::ROBUST,
    6.698790e+00 // TWO Time Fluc
);

//IonEloss HitStTOF::PDF_Q01_Q_( // OLD
//    { 3.08657e-02, 6.90186e-01, 1.81945e+00, 3.04335e+00, 1.45315e+00, 0.00000e+00 }, // Kpa
//    { 3.39070e-01, 4.05151e+00, 5.83036e-01, 1.05523e+00, 1.39493e+00 }, // Mpv
//    { 1.36404e-02, 1.69478e+00, 1.50198e+00, 1.44517e-01, 1.96318e+00 }, // Sgm
//    3.9000000e-02 // Fluc
//);

IonEloss HitStTOF::PDF_Q01_Q_( // NEW
    { 1.42180e+01, 1.03536e-01, -2.12863e+00 }, // Kpa
    { 1.10206e-02, 8.62096e+01, 9.49522e-01, 3.59949e-04, 2.70646e+00 }, // Mpv
    { 1.89554e-01, 2.80196e+00, 6.28283e-01, 4.37731e+00, 1.25767e+00 }, // Sgm
    0.0829427 // Fluc
);


} // namesapce TrackSys


#endif // __TRACKLibs_HitSt_C__
