#ifndef __TRACKLibs_HitSt_C__
#define __TRACKLibs_HitSt_C__


namespace TrackSys {
        
    
void VirtualHitSt::clear() {
    type_ = PartType::Proton;
    dec_  = Detector::NONE;
    lay_  = 0;
    
    coo_side_ = std::move(SVecO<3>(false, false, true));
    coo_      = std::move(SVecD<3>(Numc::ZERO<>, Numc::ZERO<>, Numc::ZERO<>));
    cer_      = std::move(SVecD<3>(Numc::ONE<>, Numc::ONE<>, Numc::ONE<>));
}

void HitStTRK::clear() {
    nsr_      = std::move(SVecS<2>());
    adc_side_ = std::move(SVecO<2>());
    adc_      = std::move(SVecD<2>());

    cnrm_ = std::move(SVecD<2>());
    cdiv_ = std::move(SVecD<2>());
    
    anrm_ = std::move(SVecD<2>());
    adiv_ = std::move(SVecD<2>());

    set_type(PartType::Proton);
    cer_ = std::move(SVecD<3>(DEF_CER_X_, DEF_CER_Y_, DEF_CER_Z_));
}

void HitStTRK::cal(const PhySt& part) {
    set_type(part.info().type());

    cnrm_ = std::move(SVecD<2>());
    cdiv_ = std::move(SVecD<2>());
    SVecD<3>&& crs = (coo_ - part.c());
    if (coo_side_(0)) {
        cer_(0)  = pdf_cx_->efft_sgm(crs(0));
        cnrm_(0) = crs(0) / cer_(0);
        cdiv_(0) = Numc::NEG<> / cer_(0);
    }
    if (coo_side_(1)) {
        cer_(1) = pdf_cy_->efft_sgm(crs(1));
        cnrm_(1) = crs(1) / cer_(1);
        cdiv_(1) = Numc::NEG<> / cer_(1);
    }

    anrm_ = std::move(SVecD<2>());
    adiv_ = std::move(SVecD<2>());
    if (adc_side_(0)) {
        SVecD<2>&& ionx = (*pdf_ax_)(adc_(0), part.eta());
        anrm_(0) = ionx(0);
        adiv_(0) = ionx(1);
    }
    if (adc_side_(1)) {
        SVecD<2>&& iony = (*pdf_ay_)(adc_(1), part.eta());
        anrm_(1) = iony(0);
        adiv_(1) = iony(1);
    }

    set_dummy_x(part.cx());
    set_dummy_y(part.cy());
}

void HitStTRK::set_type(PartType type) {
    if (type_ == type) return;
    switch (type) {
        case PartType::Proton :
        {
            pdf_cx_ = &PDF_PR_CX_NN_;
            if      (nsr_(0) == 1) pdf_cx_ = &PDF_PR_CX_N1_;
            else if (nsr_(0) == 2) pdf_cx_ = &PDF_PR_CX_N2_;
            else if (nsr_(0) >= 3) pdf_cx_ = &PDF_PR_CX_N3_;
            pdf_cy_ = &PDF_PR_CY_NN_;
            if      (nsr_(1) == 1) pdf_cy_ = &PDF_PR_CY_N1_;
            else if (nsr_(1) == 2) pdf_cy_ = &PDF_PR_CY_N2_;
            else if (nsr_(1) == 3) pdf_cy_ = &PDF_PR_CY_N3_;
            else if (nsr_(1) >= 4) pdf_cy_ = &PDF_PR_CY_N4_;
            pdf_ax_ = &PDF_PR_AX_;
            pdf_ay_ = &PDF_PR_AY_;
            type_ = type;
            break;
        }
        default :
            CERR("HitStTRK::set_type() NO PartType Setting.\n");
            break;
    }
}

MultiGaus HitStTRK::PDF_PR_CX_NN_(
    MultiGaus::Opt::ROBUST,
    3.31376712664997630e-01, 1.77875e-03,
    5.10255425401639595e-01, 2.65271e-03,
    1.58367861933362775e-01, 5.15837e-03
);

MultiGaus HitStTRK::PDF_PR_CX_N1_(
    MultiGaus::Opt::ROBUST,
    2.75608e-03
);

MultiGaus HitStTRK::PDF_PR_CX_N2_(
    MultiGaus::Opt::ROBUST,
    6.02616961110044480e-01, 1.82559e-03,
    3.97383038889955520e-01, 4.18164e-03
);

MultiGaus HitStTRK::PDF_PR_CX_N3_(
    MultiGaus::Opt::ROBUST,
    6.75185841999877190e-01, 1.84066e-03,
    3.24814158000122755e-01, 4.90605e-03
);

MultiGaus HitStTRK::PDF_PR_CY_NN_(
    MultiGaus::Opt::ROBUST,
    5.50759994181610257e-01, 7.90750e-04,
    3.74341189078839287e-01, 1.52916e-03,
    7.48988167395504278e-02, 3.63939e-03
);

MultiGaus HitStTRK::PDF_PR_CY_N1_(
    MultiGaus::Opt::ROBUST,
    5.45247134760751928e-01, 9.87256e-04,
    4.54752865239248016e-01, 1.65822e-03
);

MultiGaus HitStTRK::PDF_PR_CY_N2_(
    MultiGaus::Opt::ROBUST,
    4.53177772355239150e-01, 7.32341e-04,
    4.29177717623073274e-01, 1.19994e-03,
    1.17644510021687618e-01, 2.17922e-03
);

MultiGaus HitStTRK::PDF_PR_CY_N3_(
    MultiGaus::Opt::ROBUST,
    5.08190182558828307e-01, 7.74660e-04,
    3.39669567581713017e-01, 1.56378e-03,
    1.52140249859458648e-01, 3.18597e-03
);

MultiGaus HitStTRK::PDF_PR_CY_N4_(
    MultiGaus::Opt::ROBUST,
    5.19847432537280163e-01, 7.84945e-04,
    3.39808551148078952e-01, 1.72969e-03,
    1.40344016314640940e-01, 3.85202e-03
);

IonEloss HitStTRK::PDF_PR_AX_(
    { 6.84708e-04, 1.22107e+00, 1.54109e+00, 2.83791e+00, 1.06167e+00, 7.87652e+00 }, // Kpa
    { 1.01510e+00, 2.06220e+01, 1.24078e+00, 4.82421e-04, 5.80771e+00 }, // Mpv
    { 6.21439e-02, 3.05480e+01, 1.37339e+00, 1.07762e-04, 8.70839e+00 }, // Sgm
    5.00000e+00 // Fluc
);

IonEloss HitStTRK::PDF_PR_AY_(
    { 1.04571e-03, 1.56996e+00, 3.42208e+00, 1.60653e+00, 1.36027e+00, 1.00964e+01 }, // Kpa
    { 4.41904e+00, 6.24969e+00, 1.05197e+00, 1.39822e-01, 1.03750e+00 }, // Mpv
    { 7.13109e+00, 2.98634e+00, 5.79009e-01, 4.32315e+00, 7.58125e-01 }, // Sgm
    4.40000e+00 // Fluc
);


void HitStTOF::clear() {
    q_side_ = false;
    q_      = Numc::ZERO<>;

    cnrm_ = std::move(SVecD<2>());
    cdiv_ = std::move(SVecD<2>());
    
    qnrm_ = Numc::ZERO<>;
    qdiv_ = Numc::ZERO<>;

    set_type(PartType::Proton);
    cer_ = std::move(SVecD<3>(DEF_CER_X_, DEF_CER_Y_, DEF_CER_Z_));
}

void HitStTOF::cal(const PhySt& part) {
    set_type(part.info().type());

    cnrm_ = std::move(SVecD<2>());
    cdiv_ = std::move(SVecD<2>());
    SVecD<3>&& crs = (coo_ - part.c());
    if (coo_side_(0)) {
        cer_(0)  = pdf_c_->efft_sgm(crs(0));
        cnrm_(0) = crs(0) / cer_(0);
        cdiv_(0) = Numc::NEG<> / cer_(0);
    }
    if (coo_side_(1)) {
        cer_(1) = pdf_c_->efft_sgm(crs(1));
        cnrm_(1) = crs(1) / cer_(1);
        cdiv_(1) = Numc::NEG<> / cer_(1);
    }

    qnrm_ = Numc::ZERO<>;
    qdiv_ = Numc::ZERO<>;
    if (q_side_) {
        SVecD<2>&& ion = (*pdf_q_)(q_, part.eta());
        qnrm_ = ion(0);
        qdiv_ = ion(1);
    }

    set_dummy_x(part.cx());
    set_dummy_y(part.cy());
}

void HitStTOF::set_type(PartType type) {
    if (type_ == type) return;
    switch (type) {
        case PartType::Proton :
        {
            pdf_c_ = &PDF_PR_C_;
            pdf_q_ = &PDF_PR_Q_;
            type_ = type;
            break;
        }
        default :
            CERR("HitStTOF::set_type() NO PartType Setting.\n");
            break;
    }
}

MultiGaus HitStTOF::PDF_PR_C_(
    MultiGaus::Opt::ROBUST,
    2.75608e-03
);

IonEloss HitStTOF::PDF_PR_Q_(
    { 6.84708e-04, 1.22107e+00, 1.54109e+00, 2.83791e+00, 1.06167e+00, 7.87652e+00 }, // Kpa
    { 1.01510e+00, 2.06220e+01, 1.24078e+00, 4.82421e-04, 5.80771e+00 }, // Mpv
    { 6.21439e-02, 3.05480e+01, 1.37339e+00, 1.07762e-04, 8.70839e+00 }, // Sgm
    5.00000e+00 // Fluc
);












































////////////////////////////////////////////////////////////////////////////////////////

void HitSt::print() const {
    std::string printStr;
    printStr += STR("================= HitSt ==================\n");
    printStr += STR("Lay  (%d)\n", lay_);
    printStr += STR("Side (%d %d)\n", side_(0), side_(1));
    printStr += STR("Coo  (%11.6f %11.6f %11.6f)\n", coo_(0), coo_(1), coo_(2));
    printStr += STR("Err  (%11.6f %11.6f)\n", err_(0), err_(1));
    printStr += STR("==========================================\n");
    COUT(printStr.c_str());
}
        

void HitSt::set_coo(Double_t cx, Double_t cy, Double_t cz) {
    coo_(0) = cx;
    coo_(1) = cy;
    coo_(2) = cz;
}
    

void HitSt::set_adc(Double_t ax, Double_t ay) {
    adc_(0) = (side_(0) && ax > 0) ? ax : 0.;
    adc_(1) = (side_(1) && ay > 0) ? ay : 0.; 
}


void HitSt::set_nsr(Int_t nx, Int_t ny) { 
    nsr_(0) = (side_(0) && nx > 0) ? nx : 0;
    nsr_(1) = (side_(1) && ny > 0) ? ny : 0; 
}


Short_t HitSt::set_seqID(Short_t id) {
    Short_t iter = 0;
    if (side_(0)) { seqIDcx_ = id + iter; iter++; } else seqIDcx_ = -1;
    if (side_(1)) { seqIDcy_ = id + iter; iter++; } else seqIDcy_ = -1;
    if (side_(0) && adc_(0) > 0) { seqIDex_ = id + iter; iter++; } else seqIDex_ = -1;
    if (side_(1) && adc_(1) > 0) { seqIDey_ = id + iter; iter++; } else seqIDey_ = -1;
    if (iter != 0) seqID_ = id; else seqID_ = -1;
    return iter;
}


void HitSt::set_err(const PartType& type) {
    if (type == PartType::Proton) {
        type_ = type;

        pdf_cx_ = &PDF_PR_CX_NN_;
        if      (nsr_(0) == 1) pdf_cx_ = &PDF_PR_CX_N1_;
        else if (nsr_(0) == 2) pdf_cx_ = &PDF_PR_CX_N2_;
        else if (nsr_(0) >= 3) pdf_cx_ = &PDF_PR_CX_N3_;
        
        pdf_cy_ = &PDF_PR_CY_NN_;
        if      (nsr_(1) == 1) pdf_cy_ = &PDF_PR_CY_N1_;
        else if (nsr_(1) == 2) pdf_cy_ = &PDF_PR_CY_N2_;
        else if (nsr_(1) == 3) pdf_cy_ = &PDF_PR_CY_N3_;
        else if (nsr_(1) >= 4) pdf_cy_ = &PDF_PR_CY_N4_;
        
        pdf_ex_ = &PDF_PR_EX_;
        pdf_ey_ = &PDF_PR_EY_;
    }
}


} // namesapce TrackSys


#endif // __TRACKLibs_HitSt_C__
