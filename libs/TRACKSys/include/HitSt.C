#ifndef __TRACKLibs_HitSt_C__
#define __TRACKLibs_HitSt_C__


namespace TrackSys {
 /*       
void HitStTRK::clear() {
    nsr_      = std::move(SVecS<2>());
    adc_side_ = std::move(SVecO<2>());
    adc_      = std::move(SVecD<2>());

    cnrm_ = std::move(SVecD<2>());
    cdiv_ = std::move(SVecD<2>());
    
    anrm_ = std::move(SVecD<2>());
    adiv_ = std::move(SVecD<2>());

    set_type(PartType::Proton);
}

void HitStTRK::cal(const PhySt& part) {
    SVecD<3>&& res = (coo_ - part.c());
}

void HitStTRK::set_type(PartType type) {
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
            pdf_ex_ = &PDF_PR_EX_;
            pdf_ey_ = &PDF_PR_EY_;
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
    3.31376712664997630e-01, 1.77875e+01/HitSt::NORMFT_X_,
    5.10255425401639595e-01, 2.65271e+01/HitSt::NORMFT_X_,
    1.58367861933362775e-01, 5.15837e+01/HitSt::NORMFT_X_
);

*/




































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
