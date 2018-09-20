#ifndef __TRACKLibs_HitSt_C__
#define __TRACKLibs_HitSt_C__


#include "Sys.h"
#include "Math.h"
#include "TmeMeas.h"
#include "IonEloss.h"
#include "GmIonEloss.h"
#include "PartInfo.h"
#include "PhySt.h"
#include "HitSt.h"


namespace TrackSys {


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
    
    onlyc_seqID_   = -1;
    onlyc_seqIDcx_ = -1;
    onlyc_seqIDcy_ = -1;

    onlycx_seqID_ = -1;
    onlycy_seqID_ = -1;
    
    type_ = PartType::Proton;
    dec_  = Detector::NONE;
    lay_  = 0;
    
    side_coo_ = std::move(SVecO<3>());
    coo_      = std::move(SVecD<3>());
    erc_      = std::move(SVecD<2>(Numc::ONE<>, Numc::ONE<>));
    
    chic_.fill(Numc::ZERO<>);
    nrmc_.fill(Numc::ZERO<>);
    divc_.fill(Numc::ZERO<>);
}


Short_t VirtualHitSt::set_onlycx_seqID(Short_t onlycx_seqID) {
    if (onlycx_seqID < 0) { onlycx_seqID_ = -1; return 0; }
    if (side_coo_(0)) onlycx_seqID_ = onlycx_seqID;
    else              onlycx_seqID_ = -1;
    return (side_coo_(0)?1:0);
}


Short_t VirtualHitSt::set_onlycy_seqID(Short_t onlycy_seqID) {
    if (onlycy_seqID < 0) { onlycy_seqID_ = -1; return 0; }
    if (side_coo_(1)) onlycy_seqID_ = onlycy_seqID;
    else              onlycy_seqID_ = -1;
    return (side_coo_(1)?1:0);
}


Short_t VirtualHitSt::set_onlyc_seqID(Short_t onlyc_seqID) {
    if (onlyc_seqID < 0) { onlyc_seqID_ = -1; return 0; }
    
    Short_t iter = 0;
    if (side_coo_(0)) { onlyc_seqIDcx_ = onlyc_seqID + iter; iter++; } else onlyc_seqIDcx_ = -1;
    if (side_coo_(1)) { onlyc_seqIDcy_ = onlyc_seqID + iter; iter++; } else onlyc_seqIDcy_ = -1;
    if (iter != 0) onlyc_seqID_ = onlyc_seqID; else onlyc_seqID_ = -1;
    return iter;
}
 

// HitStTRK
void HitStTRK::clear() {
    seqIDcx_ = -1;
    seqIDcy_ = -1;
    seqIDqx_ = -1;
    seqIDqy_ = -1;

    side_q_.fill(false);
    q_.fill(Numc::ZERO<>);
    
    chiq_.fill(Numc::ZERO<>);
    nrmq_.fill(Numc::ZERO<>);
    divq_.fill(Numc::ZERO<>);

    pdf_cx_ = nullptr;
    pdf_cy_ = nullptr;
    pdf_qx_ = nullptr;
    pdf_qy_ = nullptr;

    set_type();
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

    chic_.fill(Numc::ZERO<>);
    nrmc_.fill(Numc::ZERO<>);
    divc_.fill(Numc::ZERO<>);
    SVecD<3>&& crs = (coo_ - part.c());
    if (side_coo_(0) && pdf_cx_ != nullptr) {
        std::array<long double, 3> minix = pdf_cx_->minimizer(crs(0));
        chic_[0] = minix.at(0);
        nrmc_[0] = minix.at(1);
        divc_[0] = Numc::NEG<> * minix.at(2);
        if (!Numc::Valid(chic_[0]) || !Numc::Valid(nrmc_[0]) || !Numc::Valid(divc_[0])) {
            chic_[0] = Numc::ZERO<>;
            nrmc_[0] = Numc::ZERO<>;
            divc_[0] = Numc::ZERO<>;
        }
    }
    if (side_coo_(1) && pdf_cy_ != nullptr) {
        std::array<long double, 3> miniy = pdf_cy_->minimizer(crs(1));
        chic_[1] = miniy.at(0);
        nrmc_[1] = miniy.at(1);
        divc_[1] = Numc::NEG<> * miniy.at(2);
        if (!Numc::Valid(chic_[1]) || !Numc::Valid(nrmc_[1]) || !Numc::Valid(divc_[1])) {
            chic_[1] = Numc::ZERO<>;
            nrmc_[1] = Numc::ZERO<>;
            divc_[1] = Numc::ZERO<>;
        }
    }

    chiq_.fill(Numc::ZERO<>);
    nrmq_.fill(Numc::ZERO<>);
    divq_.fill(Numc::ZERO<>);
    if (side_q_[0] && pdf_qx_ != nullptr) {
        std::array<long double, 3>&& ionx = pdf_qx_->minimizer(q_[0]*q_[0], part.igmbta());
        chiq_[0] = ionx.at(0);
        nrmq_[0] = ionx.at(1);
        divq_[0] = ionx.at(2) * (part.mu() * part.eta_sign());
        divq_[1] = ionx.at(2);
        
        if (!Numc::Valid(chiq_[0]) || !Numc::Valid(nrmq_[0]) || !Numc::Valid(divq_[0]) || !Numc::Valid(divq_[1])) {
            chiq_[0] = Numc::ZERO<>;
            nrmq_[0] = Numc::ZERO<>;
            divq_[0] = Numc::ZERO<>;
            divq_[1] = Numc::ZERO<>;
        }
    }
    if (side_q_[1] && pdf_qy_ != nullptr) {
        std::array<long double, 3>&& iony = pdf_qy_->minimizer(q_[1]*q_[1], part.igmbta());
        chiq_[1] = iony.at(0);
        nrmq_[1] = iony.at(1);
        divq_[2] = iony.at(2) * (part.mu() * part.eta_sign());
        divq_[3] = iony.at(2);
        
        if (!Numc::Valid(chiq_[1]) || !Numc::Valid(nrmq_[1]) || !Numc::Valid(divq_[2]) || !Numc::Valid(divq_[3])) {
            chiq_[1] = Numc::ZERO<>;
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
            pdf_cx_ = &PDF_Q01_CX_;
            pdf_cy_ = &PDF_Q01_CY_;
            pdf_qx_ = &PDF_Q01_QX_;
            pdf_qy_ = &PDF_Q01_QY_;
            type_ = info.type();
            break;
        }
        case 2 : case -2 :
        {
            pdf_cx_ = &PDF_Q02_CX_;
            pdf_cy_ = &PDF_Q02_CY_;
            pdf_qx_ = &PDF_Q02_QX_;
            pdf_qy_ = &PDF_Q02_QY_;
            type_ = info.type();
            break;
        }
        case 6 : case -6 :
        {
            pdf_cx_ = &PDF_Q06_CX_;
            pdf_cy_ = &PDF_Q06_CY_;
            type_ = info.type();
            break;
        }
        default :
            CERR("HitStTRK::set_type() NO PartType Setting.\n");
            return false;
    }
    if (pdf_cx_ != nullptr && pdf_cy_ != nullptr) {
        erc_ = std::move(SVecD<2>(
                   pdf_cx_->eftsgm(),
                   pdf_cy_->eftsgm()
               ));
    }

    return true;
}

MultiGaus HitStTRK::PDF_Q01_CX_(
    Robust(Robust::Opt::ON, 3.5L),
    7.41501407849083582e-01, 2.00206e-03,
    2.44716738877455764e-01, 3.55350e-03,
    1.37818532734607074e-02, 7.81871e-03
);

MultiGaus HitStTRK::PDF_Q01_CY_(
    Robust(Robust::Opt::ON, 3.5L),
    6.95999248725403641e-01, 7.74321e-04,
    2.65910693360992156e-01, 1.43279e-03,
    3.67491744210614660e-02, 2.86097e-03,
    1.34088349254265398e-03, 7.66690e-03
);

IonEloss HitStTRK::PDF_Q01_QX_(
    Robust(Robust::Opt::ON, 3.0L),
    { 3.63156e+04, 7.49003e-07, 1.67437e+00, 1.99507e+00 }, // Kpa
    { 1.24004e-01, 5.08937e+00, 3.07125e+00, 1.96594e+00, 3.91458e-02, 6.07854e-01 }, // Mpv
    { 2.71789e-03, 5.70229e+01, 2.81955e+01, 1.72472e+00, 1.40243e-09, 2.17644e+00 }, // Sgm
    { 1.19970e-01, 4.96139e+00, 4.21256e+00, 2.01933e+00, 3.41026e-02, 7.01322e-01 }, // Mode
    0.284416 // Fluc
);

IonEloss HitStTRK::PDF_Q01_QY_(
    Robust(Robust::Opt::ON, 3.0L),
    { 1.81137e+04, 9.99578e-07, 1.92633e+00, 1.98827e+00 }, // Kpa
    { 9.91754e-02, 7.06198e+00, 1.79925e+00, 2.01905e+00, 1.11139e-02, 6.94616e-01 }, // Mpv
    { 3.72390e-03, 5.13365e+01, 3.66595e+01, 1.76450e+00, 1.83094e-08, 2.00564e+00 }, // Sgm
    { 9.80744e-02, 7.00622e+00, 2.22096e+00, 2.03910e+00, 1.00794e-02, 7.57248e-01 }, // Mode
    0.175361 // Fluc
);

MultiGaus HitStTRK::PDF_Q02_CX_(
    Robust(Robust::Opt::ON, 3.5L),
    8.22819915794514078e-01, 1.58166e-03,
    1.58040635379077365e-01, 4.38282e-04,
    1.73994701031708332e-02, 3.41656e-03,
    1.73997872323774366e-03, 8.34167e-03
);

MultiGaus HitStTRK::PDF_Q02_CY_(
    Robust(Robust::Opt::ON, 3.5L),
    7.23281437547637074e-01, 3.04592e-04,
    2.44942684365078739e-01, 7.18113e-04,
    2.87918299651523865e-02, 1.42243e-03,
    2.98404812213171663e-03, 3.56694e-03
);

IonEloss HitStTRK::PDF_Q02_QX_(
    Robust(Robust::Opt::ON, 3.0L),
    { 3.63156e+04, 7.49003e-07, 1.67437e+00, 1.99507e+00 }, // Kpa
    { 1.24004e-01, 5.08937e+00, 3.07125e+00, 1.96594e+00, 3.91458e-02, 6.07854e-01 }, // Mpv
    { 2.71789e-03, 5.70229e+01, 2.81955e+01, 1.72472e+00, 1.40243e-09, 2.17644e+00 }, // Sgm
    { 1.19970e-01, 4.96139e+00, 4.21256e+00, 2.01933e+00, 3.41026e-02, 7.01322e-01 }, // Mode
    0.284416 // Fluc
);

IonEloss HitStTRK::PDF_Q02_QY_(
    Robust(Robust::Opt::ON, 3.0L),
    { 1.81137e+04, 9.99578e-07, 1.92633e+00, 1.98827e+00 }, // Kpa
    { 9.91754e-02, 7.06198e+00, 1.79925e+00, 2.01905e+00, 1.11139e-02, 6.94616e-01 }, // Mpv
    { 3.72390e-03, 5.13365e+01, 3.66595e+01, 1.76450e+00, 1.83094e-08, 2.00564e+00 }, // Sgm
    { 9.80744e-02, 7.00622e+00, 2.22096e+00, 2.03910e+00, 1.00794e-02, 7.57248e-01 }, // Mode
    0.175361 // Fluc
);

MultiGaus HitStTRK::PDF_Q06_CX_(
    Robust(Robust::Opt::ON, 3.5L),
    2.08789803988027378e-01, 3.48111e-04,
    5.99859624100181121e-01, 1.42959e-03,
    1.91350571911791417e-01, 1.42959e-03
);

MultiGaus HitStTRK::PDF_Q06_CY_(
    Robust(Robust::Opt::ON, 3.5L),
    4.65507932150528658e-01, 1.67264e-04,
    3.14393779518838512e-01, 3.18042e-04,
    2.12846546538308951e-01, 6.65752e-04,
    7.25174179232377942e-03, 1.13413e-03
);


// HitStTOF
Double_t HitStTOF::OFFSET_T_ = Numc::ZERO<>;
Double_t HitStTOF::OFFSET_S_ = Numc::ZERO<>;

void HitStTOF::clear() {
    seqIDcx_ = -1;
    seqIDcy_ = -1;
    seqIDt_  = -1;
    seqIDq_  = -1;
    
    side_t_ = false;
    orgt_   = Numc::ZERO<>;
    tsft_   = Numc::ZERO<>;

    side_q_ = false;
    q_      = Numc::ZERO<>;
    
    chit_    = Numc::ZERO<>;
    nrmt_    = Numc::ZERO<>;
    divtsft_ = Numc::ZERO<>;
    divt_.fill(Numc::ZERO<>);
    
    chiq_ = Numc::ZERO<>;
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

    chit_    = Numc::ZERO<>;
    nrmt_    = Numc::ZERO<>;
    divtsft_ = Numc::ZERO<>;
    divt_.fill(Numc::ZERO<>);
    if (side_t_) {
        // 1/bta := (1+igmbta*igmbta)^(1/2)
        // t     := (delta_S / bta)
        // res_t := (meas_t - part_t)
        tsft_ = orgt_ + OFFSET_T_; // update shift time
        Double_t ds = std::fabs(part.path() - OFFSET_S_);
        Double_t dt = tsft_ - part.time();
        if (Numc::Compare(ds) >= 0) {
            //Double_t sgm = ( 3.40971e+00 - 1.32705e+00 * std::erf( std::log(part.igmbta()) - 5.78158e-02 ) ); // ONE Layer [cm] // testcode
            //MultiGaus pdf_t(pdf_t_->robust(), sgm); // testcode
            std::array<long double, 3> minit = pdf_t_->minimizer(dt);
            chit_    = minit.at(0);
            nrmt_    = minit.at(1);
            divtsft_ = minit.at(2);
            Double_t divt = (Numc::NEG<> * minit.at(2) * ds) * (part.bta() * part.igmbta());
            divt_[0] = divt * (part.mu() * part.eta_sign());
            divt_[1] = divt;
            
            if (!Numc::Valid(chit_) || !Numc::Valid(nrmt_) || !Numc::Valid(divtsft_) || !Numc::Valid(divt_[0]) || !Numc::Valid(divt_[1])) {
                chit_    = Numc::ZERO<>;
                nrmt_    = Numc::ZERO<>;
                divtsft_ = Numc::ZERO<>;
                divt_.fill(Numc::ZERO<>);
            }
        }
    }
    
    chiq_ = Numc::ZERO<>;
    nrmq_ = Numc::ZERO<>;
    divq_.fill(Numc::ZERO<>);
    if (side_q_) {
        std::array<long double, 3>&& ion = pdf_q_->minimizer(q_*q_, part.igmbta());
        chiq_    = ion.at(0);
        nrmq_    = ion.at(1);
        divq_[0] = ion.at(2) * (part.mu() * part.eta_sign());
        divq_[1] = ion.at(2);
        
        if (!Numc::Valid(chiq_) || !Numc::Valid(nrmq_) || !Numc::Valid(divq_[0]) || !Numc::Valid(divq_[1])) {
            chiq_ = Numc::ZERO<>;
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
        case 2 : case -2 :
        {
            pdf_t_ = &PDF_Q02_T_;
            pdf_q_ = &PDF_Q02_Q_;
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
    Robust(Robust::Opt::ON, 3.0L),
    4.73675983474463180e+00
);

IonEloss HitStTOF::PDF_Q01_Q_(
    Robust(Robust::Opt::ON, 3.0L),
    { 7.23786e-02, 8.94961e+01, 2.54620e+01, 2.42002e+00 }, // Kpa
    { 4.47107e-01, 2.23454e+00, 6.68218e-02, 1.37587e+00, 1.03677e+00, 6.22281e-01 }, // Mpv
    { 5.25754e-03, 1.90548e+01, 2.83244e-01, 1.00000e+00, 8.98326e+01, 6.23076e+00 }, // Sgm
    { 4.38287e-01, 2.21155e+00, 2.25673e-01, 1.39373e+00, 1.06060e+00, 6.17605e-01 }, // Mode
    0.0772006 // Fluc
);

MultiGaus HitStTOF::PDF_Q02_T_(
    Robust(Robust::Opt::ON, 3.0L),
    2.26894302626795774e+00
);

IonEloss HitStTOF::PDF_Q02_Q_(
    Robust(Robust::Opt::ON, 3.0L),
    { 7.23786e-02, 8.94961e+01, 2.54620e+01, 2.42002e+00 }, // Kpa
    { 4.47107e-01, 2.23454e+00, 6.68218e-02, 1.37587e+00, 1.03677e+00, 6.22281e-01 }, // Mpv
    { 5.25754e-03, 1.90548e+01, 2.83244e-01, 1.00000e+00, 8.98326e+01, 6.23076e+00 }, // Sgm
    { 4.38287e-01, 2.21155e+00, 2.25673e-01, 1.39373e+00, 1.06060e+00, 6.17605e-01 }, // Mode
    0.0772006 // Fluc
);


// HitStRICH
void HitStRICH::clear() {
    seqIDcx_ = -1;
    seqIDcy_ = -1;
    seqIDib_ = -1;
   
    rad_ = Radiator::AGL;

    side_ib_ = false;
    ib_      = Numc::ZERO<>;

    chiib_ = Numc::ZERO<>;
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

    chiib_ = Numc::ZERO<>;
    nrmib_ = Numc::ZERO<>;
    divib_.fill(Numc::ZERO<>);
    if (side_ib_) {
        // 1/bta := (1+igmbta*igmbta)^(1/2)
        Double_t dib = ib_ - part.ibta();
        std::array<long double, 3> miniib = pdf_ib_->minimizer(dib);
        chiib_ = miniib.at(0);
        nrmib_ = miniib.at(1);
        Double_t divib = (Numc::NEG<> * miniib.at(2)) * (part.bta() * part.igmbta());
        divib_[0] = divib * (part.mu() * part.eta_sign());
        divib_[1] = divib;
        
        if (!Numc::Valid(chiib_) || !Numc::Valid(nrmib_) || !Numc::Valid(divib)) {
            chiib_ = Numc::ZERO<>;
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
        case 2 : case -2 :
        {
            switch (rad_) {
                case Radiator::AGL : pdf_ib_ = &PDF_AGL_Q02_IB_; break;
                case Radiator::NAF : pdf_ib_ = &PDF_NAF_Q02_IB_; break;
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
    Robust(Robust::Opt::ON, 3.0L),
    8.84159092979881156e-01, 9.43386e-04,
    1.15840907020118802e-01, 1.94402e-03
);

MultiGaus HitStRICH::PDF_NAF_Q01_IB_(
    Robust(Robust::Opt::ON, 3.0L),
    8.53098145164126631e-01, 3.26879e-03,
    1.46901854835873258e-01, 5.62675e-03
);

MultiGaus HitStRICH::PDF_AGL_Q02_IB_(
    Robust(Robust::Opt::ON, 3.0L),
    8.00614903314544657e-01, 6.24044e-04,
    1.99385096685455399e-01, 1.11872e-03
);

MultiGaus HitStRICH::PDF_NAF_Q02_IB_(
    Robust(Robust::Opt::ON, 3.0L),
    8.37610079454359169e-01, 1.98967e-03,
    1.62389920545640720e-01, 3.57482e-03
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

void HitStTRD::cal(const PhySt& part) {
    if (!set_type(part.info())) return;
    
    nrmel_ = Numc::ZERO<>;
    divel_.fill(Numc::ZERO<>);
    if (side_el_) {
        SVecD<3>&& lgge = (*pdf_el_)(el_, part.igmbta());
        nrmel_ = lgge(0);
        divel_[0] = lgge(1) * part.mu() * part.eta_sign();
        divel_[1] = lgge(1);
        
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
