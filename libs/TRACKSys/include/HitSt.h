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
        HitStTRK(Bool_t scx = false, Bool_t scy = false, Short_t lay = 0) : VirtualHitSt(DEC, lay, scx, scy) { clear(); }
        ~HitStTRK() { clear(); }
        
        Short_t set_seqID(Short_t seqID); 
        
        void cal(const PhySt& part);
        Bool_t set_type(const PartInfo& info = PartInfo(PartType::Proton));
        
        inline void set_q(Double_t qxy, Double_t qx, Double_t qy) {
            Bool_t side_qxy = (qxy > THRES_Q);
            Bool_t side_qx  = (qx  > THRES_Q);
            Bool_t side_qy  = (qy  > THRES_Q);
            side_q_ = (side_qxy && side_qx && side_qy);
            q_  = (side_q_ ? qxy : Numc::ZERO<>);
            qx_ = (side_q_ ? qx  : Numc::ZERO<>);
            qy_ = (side_q_ ? qy  : Numc::ZERO<>);

            chiq_ = Numc::ZERO<>;
            nrmq_ = Numc::ZERO<>;
            divq_.fill(Numc::ZERO<>);
        }
        
        inline const Short_t& seqIDq() const { return seqIDq_; }
        
        inline const Bool_t&   sq()        const { return side_q_; }
        inline const Double_t& q()         const { return q_; }
        inline const Double_t& qx()        const { return qx_; }
        inline const Double_t& qy()        const { return qy_; }
        inline const Double_t& chiq()      const { return chiq_; }
        inline const Double_t& nrmq()      const { return nrmq_; }
        inline const Double_t& divq_ibta() const { return divq_[0]; }
        inline const Double_t& divq_eta()  const { return divq_[1]; }
        inline const Double_t& divq_mu()   const { return divq_[2]; }

    protected :
        void clear();

    protected :
        Short_t seqIDq_;

        Bool_t                  side_q_;
        Double_t                q_;
        Double_t                qx_;
        Double_t                qy_;
        Double_t                chiq_; // q chi
        Double_t                nrmq_; // q nrom
        std::array<Double_t, 3> divq_; // q div (ibta) [ibta, eta, mu]
        
        CooMeas*  pdf_cx_;
        CooMeas*  pdf_cy_;
        IonEloss* pdf_q_;
    
    protected :
        static constexpr Double_t THRES_Q = 0.775;

        static CooMeas PDF_Q01_CX_;
        static CooMeas PDF_Q01_CY_;
        
        static CooMeas PDF_Q02_CX_;
        static CooMeas PDF_Q02_CY_;
        
        static IonEloss PDF_Q01_QXY_;
        static IonEloss PDF_Q02_QXY_;
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
            chit_    = Numc::ZERO<>;
            nrmt_    = Numc::ZERO<>;
            divtsft_ = Numc::ZERO<>;
            divt_.fill(Numc::ZERO<>);
        }
        
        inline void set_q(Double_t q) {
            side_q_ = (q > THRES_Q);
            q_      = (side_q_ ? q : Numc::ZERO<>);
            chiq_   = Numc::ZERO<>;
            nrmq_   = Numc::ZERO<>;
            divq_.fill(Numc::ZERO<>);
        }

        inline const Double_t& orgt() const { return orgt_; }
        inline const Double_t& tsft() const { return tsft_; }

        inline const Short_t& seqIDt() const { return seqIDt_; }
        inline const Short_t& seqIDq() const { return seqIDq_; }

        inline const Bool_t&   st()        const { return side_t_; }
        inline const Double_t& chit()      const { return chit_; }
        inline const Double_t& nrmt()      const { return nrmt_; }
        inline const Double_t& divtsft()   const { return divtsft_; }
        inline const Double_t& divt_ibta() const { return divt_[0]; }
        inline const Double_t& divt_eta()  const { return divt_[1]; }
        inline const Double_t& divt_mu()   const { return divt_[2]; }
        
        inline const Bool_t&   sq()        const { return side_q_; }
        inline const Double_t& chiq()      const { return chiq_; }
        inline const Double_t& nrmq()      const { return nrmq_; }
        inline const Double_t& divq_ibta() const { return divq_[0]; }
        inline const Double_t& divq_eta()  const { return divq_[1]; }
        inline const Double_t& divq_mu()   const { return divq_[2]; }

    protected :
        void clear();

    protected :
        Short_t seqIDt_;
        Short_t seqIDq_;
        
        Bool_t   side_t_;
        Double_t orgt_; // T [cm]
        Double_t tsft_; // T [cm]
        
        Double_t chit_; // T chi
        Double_t nrmt_; // T nrom
        Double_t divtsft_; // T(shift)
        std::array<Double_t, 3> divt_; // T div (ibta) [ibta, eta, mu]
        
        Bool_t   side_q_;
        Double_t q_; // Q

        Double_t chiq_; // Q chi
        Double_t nrmq_; // Q nrom
        std::array<Double_t, 3> divq_; // Q div (ibta) [ibta, eta, mu]

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
            chiib_   = Numc::ZERO<>;
            nrmib_   = Numc::ZERO<>;
            divib_.fill(Numc::ZERO<>);
        }
        
        inline const Radiator& rad() const { return rad_; }
        
        inline const Double_t& ib() const { return ib_; }

        inline const Short_t& seqIDib() const { return seqIDib_; }
        
        inline const Bool_t& sib() const { return side_ib_; }

        inline const Double_t& chiib() const { return chiib_; }
        inline const Double_t& nrmib() const { return nrmib_; }
        inline const Double_t& divib_ibta() const { return divib_[0]; }
        inline const Double_t& divib_eta()  const { return divib_[1]; }
        inline const Double_t& divib_mu()   const { return divib_[2]; }
        
    protected :
        void clear();

    protected :
        Radiator rad_;

        Short_t seqIDib_;
        
        Bool_t   side_ib_;
        Double_t ib_; // 1/Beta

        Double_t chiib_; // 1/Beta chi
        Double_t nrmib_; // 1/Beta nrom
        std::array<Double_t, 3> divib_; // 1/Beta div (ibta) [ibta, eta, mu]

        CherenkovMeas* pdf_ib_;
    
    protected :
        static constexpr long double RFR_INDEX_AGL = 1.0529;
        static constexpr long double RFR_INDEX_NAF = 1.33;

        static CherenkovMeas PDF_AGL_Q01_IB_;
        static CherenkovMeas PDF_NAF_Q01_IB_;
        static CherenkovMeas PDF_AGL_Q02_IB_;
        static CherenkovMeas PDF_NAF_Q02_IB_;
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
            chiel_   = Numc::ZERO<>;
            nrmel_   = Numc::ZERO<>;
            divel_.fill(Numc::ZERO<>);
        }
        
        inline const Double_t& el() const { return el_; }

        inline const Short_t& seqIDel() const { return seqIDel_; }
        
        inline const Bool_t&   sel() const { return side_el_; }
        inline const Double_t& chiel() const { return chiel_; }
        inline const Double_t& nrmel() const { return nrmel_; }
        inline const Double_t& divel_ibta() const { return divel_[0]; }
        inline const Double_t& divel_eta()  const { return divel_[1]; }
        inline const Double_t& divel_mu()   const { return divel_[2]; }
        
    protected :
        void clear();

    protected :
        Short_t seqIDel_;
        
        Bool_t   side_el_;
        Double_t el_; // energy loss dE/dx mean
        Double_t chiel_; // dE/dx chi
        Double_t nrmel_; // dE/dx nrom
        std::array<Double_t, 3> divel_; // dE/dx div (ibta) [ibta, eta, mu]
        
        IonTrEloss* pdf_el_;
    
    protected :
        static IonTrEloss PDF_Q01_EL_;
        static IonTrEloss PDF_Q02_EL_;
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
