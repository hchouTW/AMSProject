#ifndef __TRACKLibs_HitSt_H__
#define __TRACKLibs_HitSt_H__
#include <typeinfo>
#include <typeindex>


namespace TrackSys {

class VirtualHitSt {
    public :
        enum class Detector {
            NONE, TRK, TOF, TRD, RICH, ECAL
        };

    public :
        VirtualHitSt(Detector dec = Detector::NONE, Short_t lay = 0, Bool_t scx = false, Bool_t scy = false, Bool_t scz = true);
        ~VirtualHitSt() { clear(); }

        virtual Short_t set_seqID(Short_t seqID) = 0; 
        virtual Short_t set_onlycx_seqID(Short_t onlycx_seqID); 
        virtual Short_t set_onlycy_seqID(Short_t onlycy_seqID); 
        virtual Short_t set_onlyc_seqID(Short_t onlyc_seqID); 

        virtual void cal(const PhySt& part) = 0;
        virtual Bool_t set_type(const PartInfo& info = PartInfo(PartType::Proton)) = 0;
        
        inline void set_coo(Double_t cx, Double_t cy, Double_t cz) { 
            coo_ = std::move(SVecD<3>(cx, cy, cz)); 
            gstc_.fill(Numc::ZERO<>); 
            chic_.fill(Numc::ZERO<>); 
            nrmc_.fill(Numc::ZERO<>); 
            divc_.fill(Numc::ZERO<>); 
        }

        inline void set_dummy_x(Double_t cx) { if (!side_c_(0)) coo_(0) = cx; }
        inline void set_dummy_y(Double_t cy) { if (!side_c_(1)) coo_(1) = cy; }
        
        inline const Short_t& seqID()   const { return seqID_; }
        inline const Short_t& seqIDcx() const { return seqIDcx_; }
        inline const Short_t& seqIDcy() const { return seqIDcy_; }
        
        inline const Short_t& onlyc_seqID()   const { return onlyc_seqID_; }
        inline const Short_t& onlyc_seqIDcx() const { return onlyc_seqIDcx_; }
        inline const Short_t& onlyc_seqIDcy() const { return onlyc_seqIDcy_; }
        
        inline const Short_t& onlycx_seqID() const { return onlycx_seqID_; }
        inline const Short_t& onlycy_seqID() const { return onlycy_seqID_; }
        
        inline const PartType& type() const { return type_; }
        inline const Detector& dec()  const { return dec_; }
        inline const Short_t&  lay()  const { return lay_; }

        inline const Bool_t& scx() const { return side_c_(0); }
        inline const Bool_t& scy() const { return side_c_(1); }
        inline const Bool_t& scz() const { return side_c_(2); }
        
        inline const SVecD<3>& c()  const { return coo_; }
        inline const Double_t& cx() const { return coo_(0); }
        inline const Double_t& cy() const { return coo_(1); }
        inline const Double_t& cz() const { return coo_(2); }
        
        inline const Double_t& ecx() const { return erc_(0); }
        inline const Double_t& ecy() const { return erc_(1); }
        
        inline const Double_t& gstcx() const { return gstc_[0]; }
        inline const Double_t& gstcy() const { return gstc_[1]; }
        
        inline const Double_t& chicx() const { return chic_[0]; }
        inline const Double_t& chicy() const { return chic_[1]; }
        
        inline const Double_t& nrmcx() const { return nrmc_[0]; }
        inline const Double_t& nrmcy() const { return nrmc_[1]; }

        inline const Double_t& divcx() const { return divc_[0]; }
        inline const Double_t& divcy() const { return divc_[1]; }

    protected :
        void clear();

    protected :
        Short_t seqID_;
        Short_t seqIDcx_;
        Short_t seqIDcy_;
        
        Short_t onlyc_seqID_;
        Short_t onlyc_seqIDcx_;
        Short_t onlyc_seqIDcy_;
        
        Short_t onlycx_seqID_;
        Short_t onlycy_seqID_;

        PartType type_; // particle type
        
        Detector dec_; // TRK TOF TRD RICH ECAL
        Short_t  lay_; // start from 0

        SVecO<3> side_c_; // (x, y, z)
        SVecD<3> coo_;    // [cm] coord
        SVecD<2> erc_;    // [cm] error
        
        std::array<Double_t, 2> gstc_; // coord ghost
        std::array<Double_t, 2> chic_; // coord chi
        std::array<Double_t, 2> nrmc_; // coord norm
        std::array<Double_t, 2> divc_; // coord div
    
    public :
        enum class Orientation {
            kDownward = 0, kUpward = 1
        };

        inline static void Sort(std::vector<VirtualHitSt*>& hits, const Orientation& ortt = Orientation::kDownward) {
            if (hits.size() < 2) return;
            if (ortt == Orientation::kDownward) std::sort(hits.begin(), hits.end(), [](const VirtualHitSt* hit1, const VirtualHitSt* hit2) { return (hit1->cz() > hit2->cz()); } );
            else                                std::sort(hits.begin(), hits.end(), [](const VirtualHitSt* hit1, const VirtualHitSt* hit2) { return (hit1->cz() < hit2->cz()); } );
        }
};


class HitStTRK : public VirtualHitSt {
    public :
        static constexpr VirtualHitSt::Detector DEC = VirtualHitSt::Detector::TRK;

    public :
        HitStTRK(Bool_t scx = false, Bool_t scy = false, Short_t lay = 0, Bool_t isInnTr = true) : VirtualHitSt(DEC, lay, scx, scy) { clear(); isInnTr_ = isInnTr; }
        ~HitStTRK() { clear(); }
        
        Short_t set_seqID(Short_t seqID); 
        
        void cal(const PhySt& part);
        Bool_t set_type(const PartInfo& info = PartInfo(PartType::Proton));
        
        inline void set_nsr(Short_t nx, Short_t ny) {
            nsr_.at(0) = (side_c_[0] && nx > 0) ? nx : 0;
            nsr_.at(1) = (side_c_[1] && ny > 0) ? ny : 0;
        }

        inline void set_q(Double_t qx, Double_t qy, Short_t chrg = 0) {
            Short_t chrgz = std::abs(chrg);
            side_q_[0] = (Numc::Compare(qx) > (chrgz * THRES_Q));
            side_q_[1] = (Numc::Compare(qy) > (chrgz * THRES_Q));
            q_[0] = (side_q_[0] ? qx : Numc::ZERO<>);
            q_[1] = (side_q_[1] ? qy : Numc::ZERO<>);
            gstq_.fill(Numc::ZERO<>);
            chiq_.fill(Numc::ZERO<>);
            nrmq_.fill(Numc::ZERO<>);
            divq_.fill(Numc::ZERO<>);
        }
        
        inline const Short_t& seqIDqx() const { return seqIDqx_; }
        inline const Short_t& seqIDqy() const { return seqIDqy_; }
        
        inline const Short_t& nsrx() const { return nsr_[0]; }
        inline const Short_t& nsry() const { return nsr_[1]; }

        inline const Bool_t& sqx() const { return side_q_[0]; }
        inline const Bool_t& sqy() const { return side_q_[1]; }
        
        inline const Double_t& gstqx() const { return gstq_[0]; }
        inline const Double_t& gstqy() const { return gstq_[1]; }
        
        inline const Double_t& chiqx() const { return chiq_[0]; }
        inline const Double_t& chiqy() const { return chiq_[1]; }

        inline const Double_t& nrmqx() const { return nrmq_[0]; }
        inline const Double_t& nrmqy() const { return nrmq_[1]; }

        inline const Double_t& divqx_eta() const { return divq_[0]; }
        inline const Double_t& divqx_igb() const { return divq_[1]; }

        inline const Double_t& divqy_eta() const { return divq_[2]; }
        inline const Double_t& divqy_igb() const { return divq_[3]; }

    protected :
        void clear();

    protected :
        Bool_t isInnTr_; // is inner tracker
        
        std::array<Short_t, 2> nsr_; // num of strip (x, y)
        
        Short_t seqIDqx_;
        Short_t seqIDqy_;

        std::array<Bool_t, 2>   side_q_;
        std::array<Double_t, 2> q_; // ADC

        std::array<Double_t, 2> gstq_; // q ghost (x, y)
        std::array<Double_t, 2> chiq_; // q chi (x, y)
        std::array<Double_t, 2> nrmq_; // q nrom (x, y)
        std::array<Double_t, 4> divq_; // q div (igmbta) (x, y) [eta, igb]

        MultiGaus* pdf_cx_;
        MultiGaus* pdf_cy_;
        IonEloss*  pdf_qx_;
        IonEloss*  pdf_qy_;
    
    protected :
        static constexpr Double_t THRES_Q = 0.8;

        static MultiGaus PDF_Q01_CX_INN_;
        static MultiGaus PDF_Q01_CY_INN_;
        static MultiGaus PDF_Q01_CX_EXT_;
        static MultiGaus PDF_Q01_CY_EXT_;
        
        static MultiGaus PDF_Q01_CX_INN_S1_;
        static MultiGaus PDF_Q01_CY_INN_S1_;
        static MultiGaus PDF_Q01_CX_INN_S2_;
        static MultiGaus PDF_Q01_CY_INN_S2_;
        static MultiGaus PDF_Q01_CX_INN_S3_;
        static MultiGaus PDF_Q01_CY_INN_S3_;
        static MultiGaus PDF_Q01_CX_INN_S4_;
        static MultiGaus PDF_Q01_CY_INN_S4_;
        static MultiGaus PDF_Q01_CX_INN_S5_;
        static MultiGaus PDF_Q01_CY_INN_S5_;

        static MultiGaus PDF_Q01_CX_EXT_S1_;
        static MultiGaus PDF_Q01_CY_EXT_S1_;
        static MultiGaus PDF_Q01_CX_EXT_S2_;
        static MultiGaus PDF_Q01_CY_EXT_S2_;
        static MultiGaus PDF_Q01_CX_EXT_S3_;
        static MultiGaus PDF_Q01_CY_EXT_S3_;
        static MultiGaus PDF_Q01_CX_EXT_S4_;
        static MultiGaus PDF_Q01_CY_EXT_S4_;
        static MultiGaus PDF_Q01_CX_EXT_S5_;
        static MultiGaus PDF_Q01_CY_EXT_S5_;
        
        static MultiGaus PDF_Q02_CX_INN_;
        static MultiGaus PDF_Q02_CY_INN_;
        static MultiGaus PDF_Q02_CX_EXT_;
        static MultiGaus PDF_Q02_CY_EXT_;

        static IonEloss PDF_Q01_QX_;
        static IonEloss PDF_Q01_QY_;
        static IonEloss PDF_Q02_QX_;
        static IonEloss PDF_Q02_QY_;
};


class HitStTOF : public VirtualHitSt {
    public :
        static constexpr VirtualHitSt::Detector DEC = VirtualHitSt::Detector::TOF;
    
    public :
        HitStTOF(Short_t lay = 0) : VirtualHitSt(DEC, lay, false, false) { clear(); }
        ~HitStTOF() { clear(); }
        
        Short_t set_seqID(Short_t seqID); 
        
        void cal(const PhySt& part);
        Bool_t set_type(const PartInfo& info = PartInfo(PartType::Proton));

        inline void set_t(Double_t t) {
            side_t_  = (Numc::Compare(t) >= 0);
            orgt_    = (side_t_ ? t : Numc::ZERO<>);
            tsft_    = (side_t_ ? (orgt_ + OFFSET_T_) : Numc::ZERO<>);
            gstt_    = Numc::ZERO<>;
            chit_    = Numc::ZERO<>;
            nrmt_    = Numc::ZERO<>;
            divtsft_ = Numc::ZERO<>;
            divt_.fill(Numc::ZERO<>);
        }
        
        inline void set_q(Double_t q, Short_t chrg = 0) {
            Short_t chrgz = std::abs(chrg);
            side_q_ = (Numc::Compare(q) > (chrgz * THRES_Q));
            q_      = (side_q_ ? q : Numc::ZERO<>);
            gstq_   = Numc::ZERO<>;
            chiq_   = Numc::ZERO<>;
            nrmq_   = Numc::ZERO<>;
            divq_.fill(Numc::ZERO<>);
        }

        inline const Double_t& orgt() const { return orgt_; }
        inline const Double_t& tsft() const { return tsft_; }

        inline const Short_t&  seqIDt() const { return seqIDt_; }
        inline const Short_t&  seqIDq() const { return seqIDq_; }
        
        inline const Bool_t& st() const { return side_t_; }
        inline const Bool_t& sq() const { return side_q_; }

        inline const Double_t& gstt() const { return gstt_; }
        inline const Double_t& chit() const { return chit_; }
        inline const Double_t& nrmt() const { return nrmt_; }
        inline const Double_t& divtsft() const { return divtsft_; }
        inline const Double_t& divt_eta() const { return divt_[0]; }
        inline const Double_t& divt_igb() const { return divt_[1]; }
        
        inline const Double_t& gstq() const { return gstq_; }
        inline const Double_t& chiq() const { return chiq_; }
        inline const Double_t& nrmq() const { return nrmq_; }
        inline const Double_t& divq_eta() const { return divq_[0]; }
        inline const Double_t& divq_igb() const { return divq_[1]; }

    protected :
        void clear();

    protected :
        Short_t seqIDt_;
        Short_t seqIDq_;
        
        Bool_t   side_t_;
        Double_t orgt_; // T [cm]
        Double_t tsft_; // T [cm]

        Bool_t   side_q_;
        Double_t q_; // Q
        
        Double_t gstt_; // T ghost
        Double_t chit_; // T chi
        Double_t nrmt_; // T nrom
        Double_t divtsft_; // T(shift)
        std::array<Double_t, 2> divt_; // T div (igmbta) [eta, igb]

        Double_t gstq_; // Q ghost
        Double_t chiq_; // Q chi
        Double_t nrmq_; // Q nrom
        std::array<Double_t, 2> divq_; // Q div (igmbta) [eta, igb]

        TmeMeas*  pdf_t_;
        IonEloss* pdf_q_;
    
    protected :
        static constexpr Double_t THRES_Q = 0.6;
        
        static TmeMeas  PDF_Q01_T_;
        static IonEloss PDF_Q01_Q_;
        static TmeMeas  PDF_Q02_T_;
        static IonEloss PDF_Q02_Q_;

    public :
        static constexpr Double_t TRANS_NS_TO_CM = 2.99792458e+01; // [ns] -> [cm]
        static void SetOffsetPathTime(Double_t offset_s = Numc::ZERO<>, Double_t offset_t = Numc::ZERO<>, Bool_t use_tshf = false) { 
            OFFSET_S_ = offset_s; OFFSET_T_ = offset_t; USE_TSHF_ = use_tshf; 
        }
        static const Bool_t&   USED_TSHF()  { return USE_TSHF_; }
        static const Double_t& OffsetTime() { return OFFSET_T_; }
        static const Double_t& OffsetPath() { return OFFSET_S_; }

    protected :
        static Bool_t   USE_TSHF_;
        static Double_t OFFSET_S_; // move TOF to particle path (FIRST-TOF)
        static Double_t OFFSET_T_; // move TOF to particle time (FIRST-TOF | PART-TOF)
};


class HitStRICH : public VirtualHitSt {
    public :
        static constexpr VirtualHitSt::Detector DEC = VirtualHitSt::Detector::RICH;
        enum class Radiator { AGL, NAF };

    public :
        HitStRICH(const Radiator& rad = Radiator::AGL) : VirtualHitSt(DEC, 0, false, false) { clear(); rad_ = rad; }
        ~HitStRICH() { clear(); }
        
        Short_t set_seqID(Short_t seqID); 
        
        void cal(const PhySt& part);
        Bool_t set_type(const PartInfo& info = PartInfo(PartType::Proton));

        inline void set_ib(Double_t ib) {
            side_ib_ = (Numc::Compare(ib) > 0);
            ib_      = (side_ib_ ? ib : Numc::ZERO<>);
            gstib_   = Numc::ZERO<>;
            chiib_   = Numc::ZERO<>;
            nrmib_   = Numc::ZERO<>;
            divib_.fill(Numc::ZERO<>);
        }
        
        inline const Radiator& rad() const { return rad_; }
        
        inline const Double_t& ib() const { return ib_; }

        inline const Short_t& seqIDib() const { return seqIDib_; }
        
        inline const Bool_t& sib() const { return side_ib_; }

        inline const Double_t& gstib() const { return gstib_; }
        inline const Double_t& chiib() const { return chiib_; }
        inline const Double_t& nrmib() const { return nrmib_; }
        inline const Double_t& divib_eta() const { return divib_[0]; }
        inline const Double_t& divib_igb() const { return divib_[1]; }
        
    protected :
        void clear();

    protected :
        Radiator rad_;

        Short_t seqIDib_;
        
        Bool_t   side_ib_;
        Double_t ib_; // 1/Beta

        Double_t gstib_; // 1/Beta ghost
        Double_t chiib_; // 1/Beta chi
        Double_t nrmib_; // 1/Beta nrom
        std::array<Double_t, 2> divib_; // 1/Beta div (igmbta) [eta, igb]

        MultiGaus* pdf_ib_;
    
    protected :
        static MultiGaus PDF_AGL_Q01_IB_;
        static MultiGaus PDF_NAF_Q01_IB_;
        static MultiGaus PDF_AGL_Q02_IB_;
        static MultiGaus PDF_NAF_Q02_IB_;
};


class HitStTRD : public VirtualHitSt {
    public :
        static constexpr VirtualHitSt::Detector DEC = VirtualHitSt::Detector::TRD;

    public :
        HitStTRD(Short_t lay = 0) : VirtualHitSt(DEC, lay, false, false) { clear(); }
        ~HitStTRD() { clear(); }
        
        Short_t set_seqID(Short_t seqID); 
        
        void cal(const PhySt& part);
        Bool_t set_type(const PartInfo& info = PartInfo(PartType::Proton));

        inline void set_el(Double_t el) {
            side_el_ = (Numc::Compare(el) > 0);
            el_      = (side_el_ ? el : Numc::ZERO<>);
            nrmel_   = Numc::ZERO<>;
            divel_.fill(Numc::ZERO<>);
        }
        
        inline const Double_t& el() const { return el_; }

        inline const Short_t& seqIDel() const { return seqIDel_; }
        
        inline const Bool_t& sel() const { return side_el_; }

        inline const Double_t& nrmel() const { return nrmel_; }
        inline const Double_t& divel_eta() const { return divel_[0]; }
        inline const Double_t& divel_igb() const { return divel_[1]; }
        
    protected :
        void clear();

    protected :
        Short_t seqIDel_;
        
        Bool_t   side_el_;
        Double_t el_; // energy loss dE/dx

        Double_t nrmel_; // dE/dx nrom
        std::array<Double_t, 2> divel_; // dE/dx div (igmbta) [eta, igb]

        GmIonEloss* pdf_el_;
    
    protected :
        static GmIonEloss PDF_Q01_EL_;
};


template<class HitStType>
class Hit {
    public :
        inline static Bool_t     IsSame(VirtualHitSt* hit) { return (hit->dec() == HitStType::DEC); }
        inline static HitStType* Cast(VirtualHitSt* hit)   { return (IsSame(hit) ? dynamic_cast<HitStType*>(hit) : nullptr); }
        
        inline static void Sort(std::vector<HitStType>& hits, const VirtualHitSt::Orientation& ortt = VirtualHitSt::Orientation::kDownward) {
            if (hits.size() < 2) return;
            if (ortt == VirtualHitSt::Orientation::kDownward) std::sort(hits.begin(), hits.end(), [](const HitStType& hit1, const HitStType& hit2) { return (hit1.cz() > hit2.cz()); } );
            else                                              std::sort(hits.begin(), hits.end(), [](const HitStType& hit1, const HitStType& hit2) { return (hit1.cz() < hit2.cz()); } );
        }
};

template class Hit<HitStTRK>;
template class Hit<HitStTOF>;
template class Hit<HitStRICH>;
template class Hit<HitStTRD>;


} // namesapce TrackSys


#endif // __TRACKLibs_HitSt_H__
