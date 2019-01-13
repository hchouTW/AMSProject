#ifndef __TRACKLibs_HitSt_C__
#define __TRACKLibs_HitSt_C__


#include "Sys.h"
#include "Math.h"
#include "CooMeas.h"
#include "TmeMeas.h"
#include "CherenkovMeas.h"
#include "IonEloss.h"
#include "IonTrEloss.h"
#include "PartInfo.h"
#include "PhySt.h"
#include "HitSt.h"


namespace TrackSys {


VirtualHitSt::VirtualHitSt(Detector dec, Short_t lay, Bool_t scx, Bool_t scy, Bool_t scz) {
    clear();
    dec_ = dec;
    lay_ = lay;
    side_c_ = std::move(SVecO<3>(scx, scy, scz));
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
    
    side_c_ = std::move(SVecO<3>());
    coo_    = std::move(SVecD<3>());
    erc_    = std::move(SVecD<2>(Numc::ONE<>, Numc::ONE<>));
    
    chic_.fill(Numc::ZERO<>);
    nrmc_.fill(Numc::ZERO<>);
    divc_.fill(Numc::ZERO<>);
}


Short_t VirtualHitSt::set_onlycx_seqID(Short_t onlycx_seqID) {
    if (onlycx_seqID < 0) { onlycx_seqID_ = -1; return 0; }
    if (side_c_(0)) onlycx_seqID_ = onlycx_seqID;
    else            onlycx_seqID_ = -1;
    return (side_c_(0)?1:0);
}


Short_t VirtualHitSt::set_onlycy_seqID(Short_t onlycy_seqID) {
    if (onlycy_seqID < 0) { onlycy_seqID_ = -1; return 0; }
    if (side_c_(1)) onlycy_seqID_ = onlycy_seqID;
    else            onlycy_seqID_ = -1;
    return (side_c_(1)?1:0);
}


Short_t VirtualHitSt::set_onlyc_seqID(Short_t onlyc_seqID) {
    if (onlyc_seqID < 0) { onlyc_seqID_ = -1; return 0; }
    
    Short_t iter = 0;
    if (side_c_(0)) { onlyc_seqIDcx_ = onlyc_seqID + iter; iter++; } else onlyc_seqIDcx_ = -1;
    if (side_c_(1)) { onlyc_seqIDcy_ = onlyc_seqID + iter; iter++; } else onlyc_seqIDcy_ = -1;
    if (iter != 0) onlyc_seqID_ = onlyc_seqID; else onlyc_seqID_ = -1;
    return iter;
}
 

// HitStTRK
void HitStTRK::clear() {
    seqIDcx_ = -1;
    seqIDcy_ = -1;
    seqIDq_  = -1;

    side_q_ = false;
    q_      = Numc::ZERO<>;
    qx_     = Numc::ZERO<>;
    qy_     = Numc::ZERO<>;
    chiq_   = Numc::ZERO<>;
    nrmq_   = Numc::ZERO<>;
    divq_.fill(Numc::ZERO<>);

    pdf_cx_ = nullptr;
    pdf_cy_ = nullptr;
    pdf_q_  = nullptr;

    set_type();
}

Short_t HitStTRK::set_seqID(Short_t seqID) {
    if (seqID < 0) { seqID_ = -1; return 0; }

    Short_t iter = 0;
    if (side_c_[0]) { seqIDcx_ = seqID + iter; iter++; } else seqIDcx_ = -1;
    if (side_c_[1]) { seqIDcy_ = seqID + iter; iter++; } else seqIDcy_ = -1;
    if (side_q_   ) { seqIDq_  = seqID + iter; iter++; } else seqIDq_  = -1;
    if (iter != 0) seqID_ = seqID; else seqID_ = -1;
    return iter;
}

void HitStTRK::cal(const PhySt& part) {
    if (!set_type(part.info())) return;

    chic_.fill(Numc::ZERO<>);
    nrmc_.fill(Numc::ZERO<>);
    divc_.fill(Numc::ZERO<>);
    SVecD<3>&& crs = (coo_ - part.c());
    if (side_c_[0] && pdf_cx_ != nullptr) {
        std::array<long double, 3> minix = pdf_cx_->minimizer(crs(0), part.ibta(), part.igb());
        chic_[0] = minix[0];
        nrmc_[0] = minix[1];
        divc_[0] = Numc::NEG<> * minix[2];
        if (!Numc::Valid(chic_[0]) || !Numc::Valid(nrmc_[0]) || !Numc::Valid(divc_[0])) {
            chic_[0] = Numc::ZERO<>;
            nrmc_[0] = Numc::ZERO<>;
            divc_[0] = Numc::ZERO<>;
        }
    }
    if (side_c_[1] && pdf_cy_ != nullptr) {
        std::array<long double, 3> miniy = pdf_cy_->minimizer(crs(1), part.ibta(), part.igb());
        chic_[1] = miniy[0];
        nrmc_[1] = miniy[1];
        divc_[1] = Numc::NEG<> * miniy[2];
        if (!Numc::Valid(chic_[1]) || !Numc::Valid(nrmc_[1]) || !Numc::Valid(divc_[1])) {
            chic_[1] = Numc::ZERO<>;
            nrmc_[1] = Numc::ZERO<>;
            divc_[1] = Numc::ZERO<>;
        }
    }
    
    chiq_ = Numc::ZERO<>;
    nrmq_ = Numc::ZERO<>;
    divq_.fill(Numc::ZERO<>);
    if (side_q_ && pdf_q_ != nullptr) {
        std::array<long double, 3>&& ion = pdf_q_->minimizer(q_*q_, part.ibta(), part.igb());
        chiq_    = ion[0];
        nrmq_    = ion[1];
        divq_[0] = ion[2];
        divq_[1] = ion[2] * (part.bta() * part.eta()) * (part.mu() * part.mu());
        
        if (!Numc::Valid(chiq_) || !Numc::Valid(nrmq_) || !Numc::Valid(divq_[0]) || !Numc::Valid(divq_[1])) {
            chiq_    = Numc::ZERO<>;
            nrmq_    = Numc::ZERO<>;
            divq_[0] = Numc::ZERO<>;
            divq_[1] = Numc::ZERO<>;
        }
    }

    set_dummy_x(part.cx());
    set_dummy_y(part.cy());
}

Bool_t HitStTRK::set_type(const PartInfo& info) {
    if ((info.is_std() && type_ == info.type()) && (pdf_cx_ && pdf_cy_ && pdf_q_)) return true;

    Short_t absq = std::abs(info.chrg());
    if (absq == 1) {
        pdf_cx_ = &PDF_Q01_CX_;
        pdf_cy_ = &PDF_Q01_CY_;
        pdf_q_  = &PDF_Q01_QXY_;
        type_ = info.type();
    } else if (absq >= 2) {
        pdf_cx_ = &PDF_Q02_CX_;
        pdf_cy_ = &PDF_Q02_CY_;
        pdf_q_  = &PDF_Q02_QXY_;
        type_ = info.type();
    } else {
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


CooMeas HitStTRK::PDF_Q01_CX_(
    MultiGaus(Robust::Option(Robust::Opt::ON, 4.5L, 0.5L),
    7.61961330645686608e-01, 2.03263e-03,
    2.26788284954245273e-01, 3.67874e-03,
    1.12503844000680815e-02, 8.39460e-03),
    { 4.20214e-01, 5.10576e-01, 1.44097e-02, 8.20524e-03 }
);

CooMeas HitStTRK::PDF_Q01_CY_(
    MultiGaus(Robust::Option(Robust::Opt::ON, 4.5L, 0.5L),
    6.33751302183348963e-01, 7.52618e-04,
    3.04067161529244734e-01, 1.30097e-03,
    5.86416168641642407e-02, 2.46420e-03,
    3.53991942324217079e-03, 5.65686e-03),
    { -8.15267e-01, 1.61188e+00, 1.32014e-01, 2.14243e-01 }
);

CooMeas HitStTRK::PDF_Q02_CX_(
    MultiGaus(Robust::Option(Robust::Opt::ON, 4.5L, 0.5L),
    8.23483402633133132e-01, 1.58076e-03,
    1.56886848767209708e-01, 4.40970e-04,
    1.78728108648472869e-02, 3.43942e-03,
    1.75693773480982223e-03, 8.60774e-03),
    { 7.83985e-01, 1.48644e-01, 1.43967e-02, 9.28284e-03 }
);

CooMeas HitStTRK::PDF_Q02_CY_(
    MultiGaus(Robust::Option(Robust::Opt::ON, 4.5L, 0.5L),
    6.41741353237298306e-01, 2.89126e-04,
    2.13916925500147043e-01, 5.29917e-04,
    1.32054183043477047e-01, 8.97440e-04,
    1.02257539432194320e-02, 1.92346e-03,
    2.06178427585825040e-03, 3.83662e-03), 
    { 3.91096e-01, 5.42372e-01, 1.60491e-02, 1.58355e-02 }
);


IonEloss HitStTRK::PDF_Q01_QXY_(
    Robust(Robust::Option(Robust::Opt::ON, 4.0L, 1.0L)),
    { 3.08975e-01, 1.35977e+00, 7.18478e-03, 1.98220e-01, 8.06396e-01 }, // Kpa
    { 5.28872e-02, 6.90213e-01, 5.71953e-02, 3.86236e-02 }, // Mpv
    { 2.03321e-02, 7.09628e-02, 0.0, 0.0 }, // Sgm
    { 9.06801e-02, 6.83408e-01, 5.93934e-02, 4.14200e-02 }, // Mod
    0.089323 // Fluc
);

IonEloss HitStTRK::PDF_Q02_QXY_(
    Robust(Robust::Option(Robust::Opt::ON, 4.0L, 1.0L)),
    { 5.37172e-01, 8.19362e-01, 5.63577e-03, 6.49063e-01, 1.87820e+00 }, // Kpa
    { 5.52093e-02, 3.06309e+00, 1.70756e-01, 1.82024e-02 }, // Mpv
    { 0.0, 3.27483e-01, 0.0, 0.0 }, // Sgm
    { 6.31500e-01, 2.76368e+00, 8.55058e-02, 4.55644e-03 }, // Mod
    0.174984 // Fluc
);


// HitStTOF
Bool_t   HitStTOF::USE_TSHF_ = false;
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

    chit_    = Numc::ZERO<>;
    nrmt_    = Numc::ZERO<>;
    divtsft_ = Numc::ZERO<>;
    divt_.fill(Numc::ZERO<>);
    
    side_q_ = false;
    q_      = Numc::ZERO<>;

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
    if (side_t_) { seqIDt_  = seqID + iter; iter++; } else seqIDt_  = -1;
    if (side_q_) { seqIDq_  = seqID + iter; iter++; } else seqIDq_  = -1;
    if (iter != 0) seqID_ = seqID; else seqID_ = -1;
    return iter;
}

void HitStTOF::cal(const PhySt& part) {
    if (!set_type(part.info())) return;

    chit_    = Numc::ZERO<>;
    nrmt_    = Numc::ZERO<>;
    divtsft_ = Numc::ZERO<>;
    divt_.fill(Numc::ZERO<>);
    if (side_t_ && pdf_t_ != nullptr) {
        // 1/bta := (1+igb*igb)^(1/2)
        // t     := (delta_S / bta)
        // res_t := (meas_t - part_t)
        tsft_ = orgt_ + OFFSET_T_; // update shift time
        Double_t ds = std::fabs(part.path() - OFFSET_S_);
        Double_t dt = tsft_ - part.time();
        if (Numc::Compare(ds) >= 0) {
            std::array<long double, 3> minit = pdf_t_->minimizer(dt, part.ibta(), USE_TSHF_);
            chit_    = minit[0];
            nrmt_    = minit[1];
            divtsft_ = minit[2];
            Double_t divt = (Numc::NEG<> * minit[2] * ds); // d(t) / d(ibta)
            divt_[0] = divt;
            divt_[1] = divt * (part.bta() * part.eta()) * (part.mu() * part.mu());


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
    if (side_q_ && pdf_q_ != nullptr) {
        std::array<long double, 3>&& ion = pdf_q_->minimizer(q_*q_, part.ibta(), part.igb());
        chiq_    = ion[0];
        nrmq_    = ion[1];
        divq_[0] = ion[2];
        divq_[1] = ion[2] * (part.bta() * part.eta()) * (part.mu() * part.mu());
    
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
    
    Short_t absq = std::abs(info.chrg());
    if (absq == 1) {
        pdf_t_ = &PDF_Q01_T_;
        pdf_q_ = &PDF_Q01_Q_;
        type_ = info.type();
    } else if (absq >= 2) {
        pdf_t_ = &PDF_Q02_T_;
        pdf_q_ = &PDF_Q02_Q_;
        type_ = info.type();
    } else {
        CERR("HitStTOF::set_type() NO PartType Setting.\n");
        return false;
    }
    return true;
}


TmeMeas HitStTOF::PDF_Q01_T_(
    MultiGaus(Robust::Option(Robust::Opt::ON, 4.5L, 0.5L),
    7.66722158352867389e-01, 3.49983e+00,
    2.33277841647132611e-01, 5.28969e+00), 
    { 2.05267e-01, 7.94733e-01 }
);

TmeMeas HitStTOF::PDF_Q02_T_(
    MultiGaus(Robust::Option(Robust::Opt::ON, 4.5L, 0.5L),
    6.57749858134140331e-01, 2.19754e+00,
    3.42250141865859614e-01, 2.97941e+00), 
    { 7.72450e-01, 2.27550e-01 }
);


IonEloss HitStTOF::PDF_Q01_Q_(
    Robust(Robust::Option(Robust::Opt::ON, 4.0L, 0.5L)),
    { 8.77160e-01, 7.84314e-01, 5.45507e-03, 2.84523e+01, 1.61855e+01 }, // Kpa
    { 6.63128e-02, 8.60609e-01, 1.57473e-02, 3.90757e-03 }, // Mpv
    { 4.48021e-03, 7.18446e-02, 0.0, 0.0 }, // Sgm
    { 1.10112e-01, 8.47921e-01, 1.60048e-02, 4.16204e-03 }, // Mod
    0.0810673 // Fluc
);

IonEloss HitStTOF::PDF_Q02_Q_(
    Robust(Robust::Option(Robust::Opt::ON, 4.0L, 0.5L)),
    { 4.75812e-01, 2.92518e-01, 8.86139e-02, 0.00000e+00, 4.92484e-02 }, // Kpa
    { 2.01919e+00, 2.32107e+00, 1.23154e-02, 4.92361e-05 }, // Mpv
    { -4.90115e-02, 2.06891e-01, 9.39886e-02, 1.13678e-01 }, // Sgm
    { 2.18181e+00, 2.22977e+00, 7.31254e-03, 3.37621e-06 }, // Mod
    0.1588111 // Fluc
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
    if (side_ib_) { seqIDib_ = seqID + iter; iter++; } else seqIDib_ = -1;
    if (iter != 0) seqID_ = seqID; else seqID_ = -1;
    return iter;
}

void HitStRICH::cal(const PhySt& part) {
    if (!set_type(part.info())) return;

    chiib_ = Numc::ZERO<>;
    nrmib_ = Numc::ZERO<>;
    divib_.fill(Numc::ZERO<>);
    if (side_ib_ && pdf_ib_ != nullptr) {
        // 1/bta := (1+igb*igb)^(1/2)
        Double_t dib = ib_ - part.ibta();
        std::array<long double, 3> miniib = pdf_ib_->minimizer(dib, part.ibta());
        chiib_ = miniib[0];
        nrmib_ = miniib[1];
        Double_t divib = (Numc::NEG<> * miniib[2]);
        divib_[0] = divib;
        divib_[1] = divib * (part.bta() * part.eta()) * (part.mu() * part.mu());
        
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
    
    Short_t absq = std::abs(info.chrg());
    if (absq == 1) {
        switch (rad_) {
            case Radiator::AGL : pdf_ib_ = &PDF_AGL_Q01_IB_; break;
            case Radiator::NAF : pdf_ib_ = &PDF_NAF_Q01_IB_; break;
        }
        type_ = info.type();
    } else if (absq >= 2) {
        switch (rad_) {
            case Radiator::AGL : pdf_ib_ = &PDF_AGL_Q02_IB_; break;
            case Radiator::NAF : pdf_ib_ = &PDF_NAF_Q02_IB_; break;
        }
        type_ = info.type();
    } else {
        CERR("HitStRICH::set_type() NO PartType Setting.\n");
        return false;
    }
    return true;
}

CherenkovMeas HitStRICH::PDF_AGL_Q01_IB_(
    MultiGaus(Robust::Option(Robust::Opt::ON, 3.5L, 0.5L),
    6.84734464661977515e-01, 9.59147e-04,
    3.15265535338022374e-01, 1.64991e-03),
    RFR_INDEX_AGL
);

CherenkovMeas HitStRICH::PDF_AGL_Q02_IB_(
    MultiGaus(Robust::Option(Robust::Opt::ON, 3.5L, 0.5L),
    7.75514988289395357e-01, 6.15751e-04,
    2.24485011710604615e-01, 1.07079e-03),
    RFR_INDEX_NAF
);


CherenkovMeas HitStRICH::PDF_NAF_Q01_IB_(
    MultiGaus(Robust::Option(Robust::Opt::ON, 3.5L, 0.5L),
    8.53153985032869877e-01, 3.16209e-03,
    1.46846014967130206e-01, 5.38069e-03),
    RFR_INDEX_AGL
);

CherenkovMeas HitStRICH::PDF_NAF_Q02_IB_(
    MultiGaus(Robust::Option(Robust::Opt::ON, 3.5L, 0.5L),
    7.62865861184466976e-01, 1.90407e-03,
    2.37134138815532997e-01, 3.11423e-03),
    RFR_INDEX_NAF
);


void HitStTRD::clear() {
    seqIDcx_  = -1;
    seqIDcy_  = -1;
    seqIDelm_ = -1;
    seqIDels_ = -1;
    
    side_elm_ = false;
    elm_      = Numc::ZERO<>;
    chielm_   = Numc::ZERO<>;
    nrmelm_   = Numc::ZERO<>;
    divelm_.fill(Numc::ZERO<>);
    
    side_els_ = false;
    els_      = Numc::ZERO<>;
    chiels_   = Numc::ZERO<>;
    nrmels_   = Numc::ZERO<>;
    divels_.fill(Numc::ZERO<>);

    pdf_elm_ = nullptr;
    pdf_els_ = nullptr;

    set_type();
}

Short_t HitStTRD::set_seqID(Short_t seqID) {
    if (seqID < 0) { seqID_ = -1; return 0; }

    Short_t iter = 0;
    if (side_elm_) { seqIDelm_ = seqID + iter; iter++; } else seqIDelm_ = -1;
    if (side_els_) { seqIDels_ = seqID + iter; iter++; } else seqIDels_ = -1;
    if (iter != 0) seqID_ = seqID; else seqID_ = -1;
    return iter;
}

void HitStTRD::cal(const PhySt& part) {
    if (!set_type(part.info())) return;
    
    chielm_ = Numc::ZERO<>;
    nrmelm_ = Numc::ZERO<>;
    divelm_.fill(Numc::ZERO<>);
    if (side_elm_) {
        std::array<long double, 3>&& iontr = pdf_elm_->minimizer(elm_, part.ibta(), part.igb());
        chielm_    = iontr[0];
        nrmelm_    = iontr[1];
        divelm_[0] = iontr[2];
        divelm_[1] = iontr[2] * (part.bta() * part.eta()) * (part.mu() * part.mu());

        if (!Numc::Valid(chielm_) || !Numc::Valid(nrmelm_) || !Numc::Valid(divelm_[0]) || !Numc::Valid(divelm_[1])) {
            chielm_ = Numc::ZERO<>;
            nrmelm_ = Numc::ZERO<>;
            divelm_.fill(Numc::ZERO<>);
        }
    }

    chiels_ = Numc::ZERO<>;
    nrmels_ = Numc::ZERO<>;
    divels_.fill(Numc::ZERO<>);
    if (side_els_) {
        std::array<long double, 3>&& iontr = pdf_els_->minimizer(els_, part.ibta(), part.igb());
        chiels_    = iontr[0];
        nrmels_    = iontr[1];
        divels_[0] = iontr[2];
        divels_[1] = iontr[2] * (part.bta() * part.eta()) * (part.mu() * part.mu());

        if (!Numc::Valid(chiels_) || !Numc::Valid(nrmels_) || !Numc::Valid(divels_[0]) || !Numc::Valid(divels_[1])) {
            chiels_ = Numc::ZERO<>;
            nrmels_ = Numc::ZERO<>;
            divels_.fill(Numc::ZERO<>);
        }
    }

    set_dummy_x(part.cx());
    set_dummy_y(part.cy());
}

Bool_t HitStTRD::set_type(const PartInfo& info) {
    if ((info.is_std() && type_ == info.type()) && (pdf_elm_ && pdf_els_)) return true;
    
    Short_t absq = std::abs(info.chrg());
    if (absq == 1) {
        pdf_elm_ = &PDF_Q01_EL_MEN_;
        pdf_els_ = &PDF_Q01_EL_SGM_;
        type_ = info.type();
    } else if (absq >= 2) {
        pdf_elm_ = &PDF_Q02_EL_MEN_;
        pdf_els_ = &PDF_Q02_EL_SGM_;
        type_ = info.type();
    } else {
        CERR("HitStTRD::set_type() NO PartType Setting.\n");
        return false;
    }
    return true;
}

IonTrEloss HitStTRD::PDF_Q01_EL_MEN_(
    Robust(Robust::Option(Robust::Opt::ON, 3.5L, 0.5L)),
    { 4.32812e-01, 1.17304e+00, 6.69675e-02, 1.57961e-01, 1.17835e+00, 3.33984e-01, 5.29717e-01 }, // Kpa
    { -1.62672e-01, 1.74195e+00, 1.36595e-01, 4.56864e+00, 4.88916e-01, 7.08228e+00 }, // Mpv
    { -1.84424e-02, 3.21822e-01, 1.82253e-02, 1.30338e+00, 4.80543e-01, 6.78145e+00 }, // Sgm
    { -7.85233e-02, 1.71828e+00, 1.34960e-01, 4.53189e+00, 4.91589e-01, 7.12373e+00 }, // Mode
    0.20 // Fluc
);

IonTrEloss HitStTRD::PDF_Q01_EL_SGM_(
    Robust(Robust::Option(Robust::Opt::ON, 3.5L, 0.5L)),
    { 0.0, 5.56791e+00, 3.01595e-01, 3.37129e-01, 4.17765e+00, 0.0, 0.0 }, // Kpa
    { 1.10478e-01, 8.66312e-01, 5.61390e-02, 4.35985e+00, 4.58503e-01, 6.43397e+00 }, // Mpv
    { 5.68283e-02, 3.28905e-01, 1.53413e-02, 1.05370e+00, 4.08763e-01, 5.45840e+00 }, // Sgm
    { 1.12279e-01, 8.64621e-01, 5.61388e-02, 4.35974e+00, 4.58504e-01, 6.43308e+00 }, // Mode
    0.00 // Fluc
);

IonTrEloss HitStTRD::PDF_Q02_EL_MEN_(
    Robust(Robust::Option(Robust::Opt::ON, 3.5L, 0.5L)),
    { 5.02636e-01, 6.29782e-01, 3.96850e+00, 2.29894e-02, 1.54215e+00, 4.31925e-01, 5.88217e-01 }, // Kpa
    { 5.17292e-01, 6.72552e+00, 4.93536e-01, 1.73039e+01, 4.75069e-01, 6.86345e+00 }, // Mpv
    { 2.55399e-01, 7.67649e-01, 3.17313e-02, 2.83073e+00, 4.64441e-01, 6.52800e+00 }, // Sgm
    { 8.63297e-01, 6.62279e+00, 4.89280e-01, 1.71202e+01, 4.78403e-01, 6.91607e+00 }, // Mode
    0.80 // Fluc
);

IonTrEloss HitStTRD::PDF_Q02_EL_SGM_(
    Robust(Robust::Option(Robust::Opt::ON, 3.5L, 0.5L)),
    { 7.47084e-01, 1.17180e+00, 4.69442e-01, 3.83918e-01, 5.18781e+00, 0.0, 0.0}, // Kpa
    { 1.49998e+00, 2.02300e+00, 5.81187e-02, 7.36798e+00, 4.20838e-01, 5.92585e+00 }, // Mpv
    { 8.24817e-01, 5.10777e-01, 0.00000e+00, 1.69321e+00, 4.01223e-01, 5.71100e+00 }, // Sgm
    { 1.50404e+00, 2.01905e+00, 5.81185e-02, 7.36787e+00, 4.20838e-01, 5.92503e+00 }, // Mode
    0.00 // Fluc
);


} // namesapce TrackSys


#endif // __TRACKLibs_HitSt_C__
