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
        divq_[2] = ion[2] * (part.bta() * part.mu()) * (part.eta() * part.eta());
        
        if (!Numc::Valid(chiq_) || !Numc::Valid(nrmq_) || !Numc::Valid(divq_[0]) || !Numc::Valid(divq_[1]) || !Numc::Valid(divq_[2])) {
            chiq_ = Numc::ZERO<>;
            nrmq_ = Numc::ZERO<>;
            divq_.fill(Numc::ZERO<>);
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
Bool_t   HitStTOF::USE_TSHF_   = false;
Double_t HitStTOF::OFFSET_T_   = Numc::ZERO<>;
Double_t HitStTOF::OFFSET_S_   = Numc::ZERO<>;
Double_t HitStTOF::BASELINE_T_ = Numc::ZERO<>;

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
            divt_[2] = divt * (part.bta() * part.mu()) * (part.eta() * part.eta());

            if (!Numc::Valid(chit_) || !Numc::Valid(nrmt_) || !Numc::Valid(divtsft_) || !Numc::Valid(divt_[0]) || !Numc::Valid(divt_[1]) || !Numc::Valid(divt_[2])) {
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
        divq_[2] = ion[2] * (part.bta() * part.mu()) * (part.eta() * part.eta());
    
        if (!Numc::Valid(chiq_) || !Numc::Valid(nrmq_) || !Numc::Valid(divq_[0]) || !Numc::Valid(divq_[1]) || !Numc::Valid(divq_[2])) {
            chiq_ = Numc::ZERO<>;
            nrmq_ = Numc::ZERO<>;
            divq_.fill(Numc::ZERO<>);
        }
    }
    
    set_dummy_x(part.cx());
    set_dummy_y(part.cy());
}
        
void HitStTOF::recal_t(const PhySt& part, double bta) {
    bool over_light_speed = (bta >= Numc::ONE<>);
    chit_    = Numc::ZERO<>;
    nrmt_    = Numc::ZERO<>;
    divtsft_ = Numc::ZERO<>;
    divt_.fill(Numc::ZERO<>);
    if (side_t_ && pdf_t_ != nullptr) {
        // correct time when velocity > 1
        Double_t crrt = (over_light_speed) ? (orgt_ - BASELINE_T_) * (bta - Numc::ONE<>) : Numc::ZERO<>;
        // 1/bta := (1+igb*igb)^(1/2)
        // t     := (delta_S / bta)
        // res_t := (meas_t - part_t)
        tsft_ = (orgt_ + crrt) + OFFSET_T_; // update shift time
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
            divt_[2] = divt * (part.bta() * part.mu()) * (part.eta() * part.eta());

            if (!Numc::Valid(chit_) || !Numc::Valid(nrmt_) || !Numc::Valid(divtsft_) || !Numc::Valid(divt_[0]) || !Numc::Valid(divt_[1]) || !Numc::Valid(divt_[2])) {
                chit_    = Numc::ZERO<>;
                nrmt_    = Numc::ZERO<>;
                divtsft_ = Numc::ZERO<>;
                divt_.fill(Numc::ZERO<>);
            }
        }
    }
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
        divib_[2] = divib * (part.bta() * part.mu()) * (part.eta() * part.eta());
        
        if (!Numc::Valid(chiib_) || !Numc::Valid(nrmib_) || !Numc::Valid(divib_[0]) || !Numc::Valid(divib_[1]) || !Numc::Valid(divib_[2])) {
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
    seqIDcx_ = -1;
    seqIDcy_ = -1;
    seqIDel_ = -1;
    
    side_el_ = false;
    el_      = Numc::ZERO<>;
    chiel_   = Numc::ZERO<>;
    nrmel_   = Numc::ZERO<>;
    divel_.fill(Numc::ZERO<>);

    pdf_el_ = nullptr;

    set_type();
}

Short_t HitStTRD::set_seqID(Short_t seqID) {
    if (seqID < 0) { seqID_ = -1; return 0; }

    Short_t iter = 0;
    if (side_el_) { seqIDel_ = seqID + iter; iter++; } else seqIDel_ = -1;
    if (iter != 0) seqID_ = seqID; else seqID_ = -1;
    return iter;
}

void HitStTRD::cal(const PhySt& part) {
    if (!set_type(part.info())) return;
    
    chiel_ = Numc::ZERO<>;
    nrmel_ = Numc::ZERO<>;
    divel_.fill(Numc::ZERO<>);
    if (side_el_) {
        std::array<long double, 3>&& iontr = pdf_el_->minimizer(el_, part.ibta(), part.igb());
        chiel_    = iontr[0];
        nrmel_    = iontr[1];
        divel_[0] = iontr[2];
        divel_[1] = iontr[2] * (part.bta() * part.eta()) * (part.mu() * part.mu());
        divel_[2] = iontr[2] * (part.bta() * part.mu()) * (part.eta() * part.eta());

        if (!Numc::Valid(chiel_) || !Numc::Valid(nrmel_) || !Numc::Valid(divel_[0]) || !Numc::Valid(divel_[1]) || !Numc::Valid(divel_[2])) {
            chiel_ = Numc::ZERO<>;
            nrmel_ = Numc::ZERO<>;
            divel_.fill(Numc::ZERO<>);
        }
    }

    set_dummy_x(part.cx());
    set_dummy_y(part.cy());
}

Bool_t HitStTRD::set_type(const PartInfo& info) {
    if ((info.is_std() && type_ == info.type()) && pdf_el_) return true;
    
    Short_t absq = std::abs(info.chrg());
    if (absq == 1) {
        pdf_el_ = &PDF_Q01_EL_;
        type_ = info.type();
    } else if (absq >= 2) {
        pdf_el_ = &PDF_Q02_EL_;
        type_ = info.type();
    } else {
        CERR("HitStTRD::set_type() NO PartType Setting.\n");
        return false;
    }
    return true;
}

IonTrEloss HitStTRD::PDF_Q01_EL_(
    Robust(Robust::Option(Robust::Opt::ON, 3.5L, 0.5L)),
    { 4.63229e-01, 5.60790e-01, 1.56910e-01, 9.58806e-03, 3.30968e-01, 5.10023e-01, 6.75037e+00 }, // Kpa
    { -4.31028e+00, 1.13090e+01, 7.10537e-01, 5.15120e-01, 1.54264e+01, 4.83010e-01, 6.98184e+00 }, // Mpv
    { -1.70979e+00, 2.87568e+00, 3.43449e-01, 2.87868e-02, 2.50949e+00, 4.57862e-01, 6.49599e+00 }, // Sgm
    { -4.29505e+00, 1.12943e+01, 7.10707e-01, 5.15123e-01, 1.54264e+01, 4.83011e-01, 6.98133e+00 }, // Mode
    0.0 // Fluc
);

IonTrEloss HitStTRD::PDF_Q02_EL_(
    Robust(Robust::Option(Robust::Opt::ON, 3.5L, 0.5L)),
    { 4.64173e-01, 5.76225e-01, 1.95796e-01, 1.24523e-02, 2.52173e-01, 6.21123e-01, 7.99860e+00 }, // Kpa
    { -4.36915e+00, 1.14534e+01, 7.03697e-01, 5.27722e-01, 4.75415e+00, 7.34352e-01, 9.79897e+00, 1.08747e+01, 5.18116e-01, 7.82150e+00 }, // Mpv
    { -8.58152e-01, 2.10421e+00, 4.02927e-01, 2.70797e-02, 6.10239e-01, 7.70214e-01, 1.00874e+01, 1.69976e+00, 4.55862e-01, 6.77782e+00 }, // Sgm
    { -4.36915e+00, 1.14534e+01, 7.03697e-01, 5.27722e-01, 4.75415e+00, 7.34352e-01, 9.79897e+00, 1.08747e+01, 5.18116e-01, 7.82150e+00 }, // Mode
    0.0 // Fluc
);

} // namesapce TrackSys


#endif // __TRACKLibs_HitSt_C__
