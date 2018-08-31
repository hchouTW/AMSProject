#ifndef __TRACKLibs_HitSt_C__
#define __TRACKLibs_HitSt_C__


#include "Sys.h"
#include "Math.h"
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
        erc_ = std::move(SVecD<2>(pdf_cx_->sgm(0), pdf_cy_->sgm(0)));
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
        std::array<long double, 2> minix = pdf_cx_->minimizer(crs(0));
        nrmc_[0] = minix.at(0);
        divc_[0] = Numc::NEG<> * minix.at(1);
        if (!Numc::Valid(nrmc_[0]) || !Numc::Valid(divc_[0])) {
            nrmc_[0] = Numc::ZERO<>;
            divc_[0] = Numc::ZERO<>;
        }
    }
    if (side_coo_(1)) {
        std::array<long double, 2> miniy = pdf_cy_->minimizer(crs(1));
        nrmc_[1] = miniy.at(0);
        divc_[1] = Numc::NEG<> * miniy.at(1);
        if (!Numc::Valid(nrmc_[1]) || !Numc::Valid(divc_[1])) {
            nrmc_[1] = Numc::ZERO<>;
            divc_[1] = Numc::ZERO<>;
        }
    }

    nrmq_.fill(Numc::ZERO<>);
    divq_.fill(Numc::ZERO<>);
    if (side_q_[0]) {
        std::array<long double, 2>&& ionx = pdf_qx_->minimizer(q_[0]*q_[0], part.igmbta());
        nrmq_[0] = ionx.at(0);
        divq_[0] = ionx.at(1) * (part.mu() * part.eta_sign());
        divq_[1] = ionx.at(1);
        
        if (!Numc::Valid(nrmq_[0]) || !Numc::Valid(divq_[0]) || !Numc::Valid(divq_[1])) {
            nrmq_[0] = Numc::ZERO<>;
            divq_[0] = Numc::ZERO<>;
            divq_[1] = Numc::ZERO<>;
        }
    }
    if (side_q_[1]) {
        std::array<long double, 2>&& iony = pdf_qy_->minimizer(q_[1]*q_[1], part.igmbta());
        nrmq_[1] = iony.at(0);
        divq_[2] = iony.at(1) * (part.mu() * part.eta_sign());
        divq_[3] = iony.at(1);
        
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
            pdf_cx_ = &PDF_Q01_CX_;
            pdf_cy_ = &PDF_Q01_CY_;
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

MultiGaus HitStTRK::PDF_Q01_CX_(
    Robust::Opt::ON,
    7.82641585556508090e-01, 1.96972e-03,
    2.17358414443491854e-01, 5.23170e-03
);

MultiGaus HitStTRK::PDF_Q01_CY_(
    Robust::Opt::ON,
    6.49340441291077486e-01, 7.62258e-04,
    2.95655485170381982e-01, 1.34376e-03,
    5.23290530119385547e-02, 2.57560e-03,
    2.67502052660201174e-03, 6.15151e-03
);

IonEloss HitStTRK::PDF_Q01_QX_(
    { 1.00000e+00, 1.07629e+04, -1.23171e+01 }, // Kpa
    { 6.37951e-02, 1.00878e+01, -6.12065e+00, 1.92281e+00, 3.26620e-03, 1.62024e+00 }, // Mpv
    { 6.79662e-05, 1.09175e+03, -3.53937e+03, 2.38184e+00, 8.08407e-38, 5.56362e+02 }, // Sgm
    { 5.16261e-02, 1.22170e+01, -8.76369e+00, 1.93889e+00, 9.67711e-04, 2.12383e+00 }, // Mode
    0.200000 // Fluc
);

IonEloss HitStTRK::PDF_Q01_QY_(
    { 1.00000e+00, 1.38677e+01, -5.63572e+00 }, // Kpa
    { 1.13078e-01, 6.63146e+00, -8.65647e-01, 1.97895e+00, 1.44020e-02, 1.33740e+00 }, // Mpv
    { 7.16082e-03, 3.30044e+01, -1.37265e+01, 1.63458e+00, 5.48096e-06, 2.97555e+00 }, // Sgm
    { 2.93377e-02, 1.83070e+01, -2.21313e+01, 2.29597e+00, 4.55659e+00, 8.05005e+00 }, // Mode
    0.110000 // Fluc
);


// HitStTOF
Double_t HitStTOF::OFFSET_T_ = Numc::ZERO<>;
Double_t HitStTOF::OFFSET_S_ = Numc::ZERO<>;
Bool_t HitStTOF::TShiftCorr_ = true;

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
            MultiGaus pdf_t(pdf_t_->robust(), sgm);
            std::array<long double, 2> minit = pdf_t.minimizer(dt);
            nrmt_     = minit.at(0);
            divt_sft_ = minit.at(1);
            Double_t divt = (Numc::NEG<> * minit.at(1) * ds) * (part.bta() * part.igmbta());
            divt_[0] = divt * (part.mu() * part.eta_sign());
            divt_[1] = divt;
            
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
        std::array<long double, 2>&& ion = pdf_q_->minimizer(q_*q_, part.igmbta());
        nrmq_    = ion.at(0);
        divq_[0] = ion.at(1) * (part.mu() * part.eta_sign());
        divq_[1] = ion.at(1);
        
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
    Robust::Opt::ON,
    6.698790e+00 // TWO Time Fluc
);

IonEloss HitStTOF::PDF_Q01_Q_(
    { 1.39506e+00, 1.80190e+00, -2.23104e+00 }, // Kpa
    { 4.41922e-01, 2.24614e+00, 1.17808e-02, 1.36798e+00, 9.43086e-01, 1.26089e+00 }, // Mpv
    { 5.77739e+00, 1.63125e+02, 1.63097e+02, 6.20364e-03, 1.01299e+00, 1.97666e+00 }, // Sgm
    //{ 5.97153e+00, -8.51391e-01, 5.05831e-01, -9.96197e+00, -4.13979e-01, 1.87597e+00, 2.74693e-01 }, // Shft
    { 2.54924e-01, 3.03910e+00, 2.94049e-01, 1.40558e+00, 2.76002e-01, 1.78231e+00 }, // Mode
    0.0751353 // Fluc
);

//IonEloss HitStTOF::PDF_Q01_Q_(
//    { 1.00000e+00, 1.23333e+00, -2.31113e+00, 1.97051e+00 }, // Kpa
//    { 4.90401e-01, 2.21732e+00, 4.49282e-03, 1.33820e+00, 1.11419e+00, 1.21069e+00 }, // Mpv
//    { 1.21818e+00, 9.31334e+00, 9.55853e+00, 1.42066e-01, 7.40717e-01, 1.98209e+00 }, // Sgm
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
        std::array<long double, 2> miniib = pdf_ib_->minimizer(dib);
        nrmib_ = miniib.at(0);
        Double_t divib = (Numc::NEG<> * miniib.at(1)) * (part.bta() * part.igmbta());
        divib_[0] = divib * (part.mu() * part.eta_sign());
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
    Robust::Opt::ON,
    8.84159092979881156e-01, 9.43386e-04,
    1.15840907020118802e-01, 1.94402e-03
);

MultiGaus HitStRICH::PDF_NAF_Q01_IB_(
    Robust::Opt::ON,
    8.53098145164126631e-01, 3.26879e-03,
    1.46901854835873258e-01, 5.62675e-03
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
