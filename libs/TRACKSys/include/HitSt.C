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
    coo_      = std::move(SVecD<3>());
    erc_      = std::move(SVecD<2>(Numc::ONE<>, Numc::ONE<>));
    
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
    isInnTr_ = true;
    
    nsr_.fill(Numc::ZERO<Short_t>);

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
    if (side_c_(0)) { seqIDcx_ = seqID + iter; iter++; } else seqIDcx_ = -1;
    if (side_c_(1)) { seqIDcy_ = seqID + iter; iter++; } else seqIDcy_ = -1;
    if (side_q_[0]) { seqIDqx_ = seqID + iter; iter++; } else seqIDqx_ = -1;
    if (side_q_[1]) { seqIDqy_ = seqID + iter; iter++; } else seqIDqy_ = -1;
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
    if (side_c_[1] && pdf_cy_ != nullptr) {
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
        
        // testcode
        if (chiq_[0] < LMTL_CHI_Q) { nrmq_[0] = 0.0; divq_[0] = 0.0; divq_[1] = 0.0; }

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
        
        // testcode
        if (chiq_[1] < LMTL_CHI_Q) { nrmq_[1] = 0.0; divq_[2] = 0.0; divq_[3] = 0.0; }
        
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
            if (isInnTr_) { 
                pdf_cx_ = &PDF_Q01_CX_INN_;
                pdf_cy_ = &PDF_Q01_CY_INN_;
                if      (nsr_[0] == 1) pdf_cx_ = &PDF_Q01_CX_INN_S1_;
                else if (nsr_[0] == 2) pdf_cx_ = &PDF_Q01_CX_INN_S2_;
                else if (nsr_[0] == 3) pdf_cx_ = &PDF_Q01_CX_INN_S3_;
                else if (nsr_[0] == 4) pdf_cx_ = &PDF_Q01_CX_INN_S4_;
                else if (nsr_[0] == 5) pdf_cx_ = &PDF_Q01_CX_INN_S5_;
                if      (nsr_[1] == 1) pdf_cy_ = &PDF_Q01_CY_INN_S1_;
                else if (nsr_[1] == 2) pdf_cy_ = &PDF_Q01_CY_INN_S2_;
                else if (nsr_[1] == 3) pdf_cy_ = &PDF_Q01_CY_INN_S3_;
                else if (nsr_[1] == 4) pdf_cy_ = &PDF_Q01_CY_INN_S4_;
                else if (nsr_[1] == 5) pdf_cy_ = &PDF_Q01_CY_INN_S5_;
            }
            else {
                pdf_cx_ = &PDF_Q01_CX_EXT_;
                pdf_cy_ = &PDF_Q01_CY_EXT_;
                if      (nsr_[0] == 1) pdf_cx_ = &PDF_Q01_CX_EXT_S1_;
                else if (nsr_[0] == 2) pdf_cx_ = &PDF_Q01_CX_EXT_S2_;
                else if (nsr_[0] == 3) pdf_cx_ = &PDF_Q01_CX_EXT_S3_;
                else if (nsr_[0] == 4) pdf_cx_ = &PDF_Q01_CX_EXT_S4_;
                else if (nsr_[0] == 5) pdf_cx_ = &PDF_Q01_CX_EXT_S5_;
                if      (nsr_[1] == 1) pdf_cy_ = &PDF_Q01_CY_EXT_S1_;
                else if (nsr_[1] == 2) pdf_cy_ = &PDF_Q01_CY_EXT_S2_;
                else if (nsr_[1] == 3) pdf_cy_ = &PDF_Q01_CY_EXT_S3_;
                else if (nsr_[1] == 4) pdf_cy_ = &PDF_Q01_CY_EXT_S4_;
                else if (nsr_[1] == 5) pdf_cy_ = &PDF_Q01_CY_EXT_S5_;
            }
            pdf_qx_ = &PDF_Q01_QX_;
            pdf_qy_ = &PDF_Q01_QY_;
            type_ = info.type();
            break;
        }
        case 2 : case -2 :
        {
            if (isInnTr_) { 
                pdf_cx_ = &PDF_Q02_CX_INN_;
                pdf_cy_ = &PDF_Q02_CY_INN_;
            }
            else {
                pdf_cx_ = &PDF_Q02_CX_EXT_;
                pdf_cy_ = &PDF_Q02_CY_EXT_;
            }
            pdf_qx_ = &PDF_Q02_QX_;
            pdf_qy_ = &PDF_Q02_QY_;
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


MultiGaus HitStTRK::PDF_Q01_CX_INN_(
    Robust(Robust::Opt::ON, 3.5L),
    7.61961330645686608e-01, 2.03263e-03,
    2.26788284954245273e-01, 3.67874e-03,
    1.12503844000680815e-02, 8.39460e-03
);

MultiGaus HitStTRK::PDF_Q01_CY_INN_(
    Robust(Robust::Opt::ON, 3.5L),
    6.33751302183348963e-01, 7.52618e-04,
    3.04067161529244734e-01, 1.30097e-03,
    5.86416168641642407e-02, 2.46420e-03,
    3.53991942324217079e-03, 5.65686e-03
);

MultiGaus HitStTRK::PDF_Q01_CX_EXT_(
    Robust(Robust::Opt::ON, 3.5L),
    5.36810408155873331e-01, 1.84480e-03,
    4.44451633843737703e-01, 3.22181e-03,
    1.87379580003889309e-02, 7.59401e-03
);

MultiGaus HitStTRK::PDF_Q01_CY_EXT_(
    Robust(Robust::Opt::ON, 3.5L),
    5.95352696369181866e-01, 1.06724e-03,
    3.24409486956099280e-01, 1.55332e-03,
    7.34103231635499370e-02, 2.69880e-03,
    6.82749351116899734e-03, 5.26917e-03
);

MultiGaus HitStTRK::PDF_Q01_CX_INN_S1_(
    Robust(Robust::Opt::ON, 3.5L),
    9.76797404354876053e-01, 2.64251e-03,
    2.32025956451238566e-02, 6.95598e-03
);

MultiGaus HitStTRK::PDF_Q01_CY_INN_S1_(
    Robust(Robust::Opt::ON, 3.5L),
    7.33529943942763829e-01, 1.12865e-03,
    2.48219756537832709e-01, 1.81506e-03,
    1.82502995194035143e-02, 3.63912e-03
);

MultiGaus HitStTRK::PDF_Q01_CX_INN_S2_(
    Robust(Robust::Opt::ON, 3.5L),
    6.54967856835809936e-01, 1.79355e-03,
    3.39622662994680469e-01, 3.33380e-03,
    5.40948016950957081e-03, 6.99129e-03
);

MultiGaus HitStTRK::PDF_Q01_CY_INN_S2_(
    Robust(Robust::Opt::ON, 3.5L),
    7.25079987394906200e-01, 8.79921e-04,
    2.50580069821453977e-01, 1.44768e-03,
    2.43399427836398438e-02, 2.85266e-03
);

MultiGaus HitStTRK::PDF_Q01_CX_INN_S3_(
    Robust(Robust::Opt::ON, 3.5L),
    7.78944110488754338e-01, 1.88904e-03,
    2.06679046809400080e-01, 4.06013e-03,
    1.43768427018455668e-02, 8.54558e-03
);

MultiGaus HitStTRK::PDF_Q01_CY_INN_S3_(
    Robust(Robust::Opt::ON, 3.5L),
    7.56069971935815599e-01, 9.09088e-04,
    2.08680431553891349e-01, 1.75711e-03,
    3.52495965102930103e-02, 3.88510e-03
);

MultiGaus HitStTRK::PDF_Q01_CX_INN_S4_(
    Robust(Robust::Opt::ON, 3.5L),
    7.93036362622009605e-01, 1.87295e-03,
    1.87045606202200482e-01, 4.67200e-03,
    1.99180311757898333e-02, 1.07436e-02
);

MultiGaus HitStTRK::PDF_Q01_CY_INN_S4_(
    Robust(Robust::Opt::ON, 3.5L),
    5.79134848510354350e-01, 8.48334e-04,
    3.26263023537911367e-01, 1.51258e-03,
    8.20015606890014681e-02, 3.43712e-03,
    1.26005672627326291e-02, 6.82896e-03
);

MultiGaus HitStTRK::PDF_Q01_CX_INN_S5_(
    Robust(Robust::Opt::ON, 3.5L),
    7.94651298159730457e-01, 1.89048e-03,
    1.81898822964718127e-01, 4.69307e-03,
    2.34498788755514269e-02, 1.31538e-02
);

MultiGaus HitStTRK::PDF_Q01_CY_INN_S5_(
    Robust(Robust::Opt::ON, 3.5L),
    5.92314252634837390e-01, 8.50086e-04,
    3.10451855879154470e-01, 1.53744e-03,
    7.01804729695321983e-02, 3.41980e-03,
    2.70534185164758474e-02, 6.94972e-03
);

MultiGaus HitStTRK::PDF_Q01_CX_EXT_S1_(
    Robust(Robust::Opt::ON, 3.5L),
    9.83408942163766953e-01, 2.77597e-03,
    1.65910578362329572e-02, 7.03008e-03
);

MultiGaus HitStTRK::PDF_Q01_CY_EXT_S1_(
    Robust(Robust::Opt::ON, 3.5L),
    7.01774331427257092e-01, 1.47035e-03,
    2.75795035425676571e-01, 2.20973e-03,
    2.24306331470663783e-02, 3.69646e-03
);

MultiGaus HitStTRK::PDF_Q01_CX_EXT_S2_(
    Robust(Robust::Opt::ON, 3.5L),
    5.42316416038590265e-01, 1.70502e-03,
    4.52793442991823669e-01, 3.24588e-03,
    4.89014096958611307e-03, 7.66383e-03
);

MultiGaus HitStTRK::PDF_Q01_CY_EXT_S2_(
    Robust(Robust::Opt::ON, 3.5L),
    8.02263996603330476e-01, 1.29785e-03,
    1.72678496642768659e-01, 1.99864e-03,
    2.50575067539009731e-02, 3.42195e-03
);

MultiGaus HitStTRK::PDF_Q01_CX_EXT_S3_(
    Robust(Robust::Opt::ON, 3.5L),
    6.60095044685769738e-01, 1.81411e-03,
    3.23478449251639544e-01, 3.70079e-03,
    1.64265060625908600e-02, 8.47580e-03
);

MultiGaus HitStTRK::PDF_Q01_CY_EXT_S3_(
    Robust(Robust::Opt::ON, 3.5L),
    8.27508524457650996e-01, 1.32903e-03,
    1.33528726995349289e-01, 2.44747e-03,
    3.89627485469997911e-02, 4.47588e-03
);

MultiGaus HitStTRK::PDF_Q01_CX_EXT_S4_(
    Robust(Robust::Opt::ON, 3.5L),
    7.01805552229827101e-01, 1.84649e-03,
    2.71688854582317341e-01, 4.18918e-03,
    2.65055931878556555e-02, 9.88614e-03
);

MultiGaus HitStTRK::PDF_Q01_CY_EXT_S4_(
    Robust(Robust::Opt::ON, 3.5L),
    7.93966871796025653e-01, 1.35410e-03,
    1.40004403665006466e-01, 2.68690e-03,
    6.60287245389678812e-02, 5.22093e-03
);

MultiGaus HitStTRK::PDF_Q01_CX_EXT_S5_(
    Robust(Robust::Opt::ON, 3.5L),
    6.93444716388343041e-01, 1.85022e-03,
    2.68734283508543670e-01, 4.10178e-03,
    3.78210001031133031e-02, 1.08402e-02
);

MultiGaus HitStTRK::PDF_Q01_CY_EXT_S5_(
    Robust(Robust::Opt::ON, 3.5L),
    7.93076536980707236e-01, 1.34906e-03,
    1.36321533051417448e-01, 2.70994e-03,
    7.06019299678752460e-02, 5.94222e-03
);

MultiGaus HitStTRK::PDF_Q02_CX_INN_(
    Robust(Robust::Opt::ON, 3.5L),
    8.23483402633133132e-01, 1.58076e-03,
    1.56886848767209708e-01, 4.40970e-04,
    1.78728108648472869e-02, 3.43942e-03,
    1.75693773480982223e-03, 8.60774e-03
);

MultiGaus HitStTRK::PDF_Q02_CY_INN_(
    Robust(Robust::Opt::ON, 3.5L),
    6.41741353237298306e-01, 2.89126e-04,
    2.13916925500147043e-01, 5.29917e-04,
    1.32054183043477047e-01, 8.97440e-04,
    1.02257539432194320e-02, 1.92346e-03,
    2.06178427585825040e-03, 3.83662e-03
);

MultiGaus HitStTRK::PDF_Q02_CX_EXT_(
    Robust(Robust::Opt::ON, 3.5L),
    8.46914506107380061e-01, 1.67402e-03,
    9.50578907857529132e-02, 2.56901e-03,
    5.45306637637586744e-02, 5.20314e-04,
    3.49693934310836519e-03, 7.57400e-03
);

MultiGaus HitStTRK::PDF_Q02_CY_EXT_(
    Robust(Robust::Opt::ON, 3.5L),
    5.57047413429423655e-01, 7.27972e-04,
    3.93570943518477712e-01, 9.99026e-04,
    4.19827792386289200e-02, 1.60732e-03,
    6.09344578181839416e-03, 3.13326e-03,
    1.30541803165129056e-03, 4.54254e-03
);

IonEloss HitStTRK::PDF_Q01_QX_(
    Robust(Robust::Opt::ON, 3.0L),
    { 3.83633e+02, 1.22278e-03, 6.08278e-01, 2.28694e+00 }, // Kpa
    { 2.60420e-02, 2.03425e+01, 6.97964e+00, 1.14792e+00, 1.51381e-03, 1.00000e+00 }, // Mpv
    { 1.28682e-03, 4.30291e+01, 2.62414e+01, 1.33026e+00, 1.19680e-02, 6.87456e-01 }, // Sgm
    { 2.58153e-02, 1.99418e+01, 9.93034e+00, 1.15273e+00, 1.65802e-03, 9.52642e-01 }, // Mode
    0.135 // Fluc
);

IonEloss HitStTRK::PDF_Q01_QY_(
    Robust(Robust::Opt::ON, 3.0L),
    { 9.94116e+00, 1.91792e+00, 1.20427e-01, 1.27033e+01 }, // Kpa
    { 3.57557e-02, 2.13704e+01,  1.49338e-01, 1.05680e+00, 1.86159e-03, 1.00000e+00 }, // Mpv
    { 2.36976e-02, 4.77391e+00, -1.02411e+00, 1.19441e+00, 5.82170e-01, 5.01835e-01 }, // Sgm
    { 3.49041e-02, 2.10310e+01,  3.00574e+00, 1.06691e+00, 1.93229e-03, 9.57865e-01 }, // Mode
    0.15 // Fluc
);

IonEloss HitStTRK::PDF_Q02_QX_(
    Robust(Robust::Opt::ON, 3.0L),
    { 4.09744e+04, 2.01209e-03, 5.41564e-03, 8.34717e+01 }, // Kpa
    { 1.21035e-01, 2.44086e+01,  1.70501e+00, 1.09843e+00, 1.37206e-03, 1.00000e+00 }, // Mpv
    { 3.72130e+00, 3.47179e+00, -6.79454e-02, 8.93500e-01, 2.58386e+01, 3.74295e-01 }, // Sgm
    { 1.21051e-01, 2.44113e+01,  1.72232e+00, 1.09842e+00, 1.35359e-03, 1.00686e+00 }, // Mode
    0.08 // Fluc
);

IonEloss HitStTRK::PDF_Q02_QY_(
    Robust(Robust::Opt::ON, 3.0L),
    { 3.80213e+00, 5.81308e+03, 4.34709e-01, 3.47169e+01 }, // Kpa
    { 1.49870e-01, 2.17453e+01, -6.95377e-01, 1.03360e+00, 2.69577e-03, 1.00000e+00 }, // Mpv
    { 5.43837e-01, 2.02673e+00, -1.29748e+00, 1.23882e+00, 1.18455e+00, 6.55918e-01 }, // Sgm
    { 1.42602e-01, 2.07106e+01,  3.17443e+00, 1.07317e+00, 2.76366e-03, 9.47687e-01 }, // Mode
    0.49 // Fluc
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
    if (side_c_(0)) { seqIDcx_ = seqID + iter; iter++; } else seqIDcx_ = -1;
    if (side_c_(1)) { seqIDcy_ = seqID + iter; iter++; } else seqIDcy_ = -1;
    if (side_t_   ) { seqIDt_  = seqID + iter; iter++; } else seqIDt_  = -1;
    if (side_q_   ) { seqIDq_  = seqID + iter; iter++; } else seqIDq_  = -1;
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
        // 1/bta := (1+igmbta*igmbta)^(1/2)
        // t     := (delta_S / bta)
        // res_t := (meas_t - part_t)
        tsft_ = orgt_ + OFFSET_T_; // update shift time
        Double_t ds = std::fabs(part.path() - OFFSET_S_);
        Double_t dt = tsft_ - part.time();
        if (Numc::Compare(ds) >= 0) {
            std::array<long double, 3> minit = pdf_t_->minimizer(dt, part.igmbta(), USE_TSHF_);
            chit_    = minit.at(0);
            nrmt_    = minit.at(1);
            if (USE_TSHF_) divtsft_ = minit.at(2);
            Double_t divt = (Numc::NEG<> * minit.at(2) * ds);
            divt_[0] = divt * (part.bta() * part.eta()) * (part.mu() * part.mu());
            divt_[1] = divt * (part.bta() * part.igmbta());
            
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
        std::array<long double, 3>&& ion = pdf_q_->minimizer(q_*q_, part.igmbta());
        chiq_    = ion.at(0);
        nrmq_    = ion.at(1);
        divq_[0] = ion.at(2) * (part.mu() * part.eta_sign());
        divq_[1] = ion.at(2);
    
        // testcode
        if (chiq_ < LMTL_CHI_Q) { nrmq_ = Numc::ZERO<>; divq_[0] = Numc::ZERO<>; divq_[1] = Numc::ZERO<>; }

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

TmeMeas HitStTOF::PDF_Q01_T_(
    Robust(Robust::Opt::ON, 3.0L),
    { 2.38464e+00, 2.15741e+00, 7.72451e+01, 1.03508e-02, 7.73432e+01 } // Sgm
);

IonEloss HitStTOF::PDF_Q01_Q_(
    Robust(Robust::Opt::ON, 3.0L),
    { 2.94768e+01, 2.21358e+01, 6.06181e-02, 9.45926e+01 }, // Kpa
    { 3.81512e-01, 2.60771e+00, -2.27962e-01, 1.27420e+00, 7.54779e-01, 6.57663e-01 }, // Mpv
    { 5.53873e-01, 8.14537e+00, -8.41010e+00, 1.87171e-01, 6.68142e-01, 9.80451e-01 }, // Sgm
    { 3.77430e-01, 2.61789e+00, -4.19394e-02, 1.26715e+00, 8.19513e-01, 6.40270e-01 }, // Mode
    0.082 // Fluc
);

TmeMeas HitStTOF::PDF_Q02_T_(
    Robust(Robust::Opt::ON, 3.0L),
    { 2.13989e+00, 2.36987e-01, 7.71676e+01, 1.72728e-02, 7.76650e+01 } // Sgm
);

IonEloss HitStTOF::PDF_Q02_Q_(
    Robust(Robust::Opt::ON, 3.0L),
    { 3.17895e+00, 1.64932e+04, 5.05553e-01, 3.19234e+01 }, // Kpa
    { 1.00385e-02, 2.71829e+02,  1.53998e+02, 9.24822e-01, 2.50440e-10, 3.01499e+00 }, // Mpv
    { 4.41616e+00, 6.46553e+01, -6.48173e+01, 1.84870e-02, 7.73260e-01, 9.88732e-01 }, // Sgm
    { 1.00385e-02, 2.71829e+02,  1.53998e+02, 9.24822e-01, 2.50440e-10, 3.01499e+00 }, // Mode
    0.00 // Fluc
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
    if (side_c_(0)) { seqIDcx_ = seqID + iter; iter++; } else seqIDcx_ = -1;
    if (side_c_(1)) { seqIDcy_ = seqID + iter; iter++; } else seqIDcy_ = -1;
    if (side_ib_  ) { seqIDib_ = seqID + iter; iter++; } else seqIDib_ = -1;
    if (iter != 0) seqID_ = seqID; else seqID_ = -1;
    return iter;
}

void HitStRICH::cal(const PhySt& part) {
    if (!set_type(part.info())) return;

    chiib_ = Numc::ZERO<>;
    nrmib_ = Numc::ZERO<>;
    divib_.fill(Numc::ZERO<>);
    if (side_ib_ && pdf_ib_ != nullptr) {
        // 1/bta := (1+igmbta*igmbta)^(1/2)
        Double_t dib = ib_ - part.ibta();
        std::array<long double, 3> miniib = pdf_ib_->minimizer(dib);
        chiib_ = miniib.at(0);
        nrmib_ = miniib.at(1);
        Double_t divib = (Numc::NEG<> * miniib.at(2));
        divib_[0] = divib * (part.bta() * part.eta()) * (part.mu() * part.mu());
        divib_[1] = divib * (part.bta() * part.igmbta());
        
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
    6.84734464661977515e-01, 9.59147e-04,
    3.15265535338022374e-01, 1.64991e-03
);

MultiGaus HitStRICH::PDF_NAF_Q01_IB_(
    Robust(Robust::Opt::ON, 3.0L),
    8.53153985032869877e-01, 3.16209e-03,
    1.46846014967130206e-01, 5.38069e-03
);

MultiGaus HitStRICH::PDF_AGL_Q02_IB_(
    Robust(Robust::Opt::ON, 3.0L),
    7.75514988289395357e-01, 6.15751e-04,
    2.24485011710604615e-01, 1.07079e-03
);

MultiGaus HitStRICH::PDF_NAF_Q02_IB_(
    Robust(Robust::Opt::ON, 3.0L),
    7.62865861184466976e-01, 1.90407e-03,
    2.37134138815532997e-01, 3.11423e-03
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
    if (side_c_(0)) { seqIDcx_ = seqID + iter; iter++; } else seqIDcx_ = -1;
    if (side_c_(1)) { seqIDcy_ = seqID + iter; iter++; } else seqIDcy_ = -1;
    if (side_el_  ) { seqIDel_ = seqID + iter; iter++; } else seqIDel_ = -1;
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
