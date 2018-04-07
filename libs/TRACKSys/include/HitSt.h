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

        virtual void cal(const PhySt& part) = 0;
        virtual void set_type(const PartInfo& info = PartInfo(PartType::Proton)) = 0;
        
        inline void set_coo(Double_t cx, Double_t cy, Double_t cz) { coo_ = std::move(SVecD<3>(cx, cy, cz)); }
        inline void set_dummy_x(Double_t cx) { if (!side_coo_(0)) coo_(0) = cx; }
        inline void set_dummy_y(Double_t cy) { if (!side_coo_(1)) coo_(1) = cy; }
        
        inline const Short_t&  seqID() const { return seqID_; }
        inline const Short_t&  seqIDcx() const { return seqIDcx_; }
        inline const Short_t&  seqIDcy() const { return seqIDcy_; }
        
        inline const PartType& type() const { return type_; }
        inline const Detector& dec()  const { return dec_; }
        inline const Short_t&  lay()  const { return lay_; }

        inline const SVecO<3>& sc()  const { return side_coo_; }
        inline const Bool_t&   scx() const { return side_coo_(0); }
        inline const Bool_t&   scy() const { return side_coo_(1); }
        inline const Bool_t&   scz() const { return side_coo_(2); }
        
        inline const SVecD<3>& c()  const { return coo_; }
        inline const Double_t& cx() const { return coo_(0); }
        inline const Double_t& cy() const { return coo_(1); }
        inline const Double_t& cz() const { return coo_(2); }
        
        inline const SVecD<2>& ec()  const { return erc_; }
        inline const Double_t& ecx() const { return erc_(0); }
        inline const Double_t& ecy() const { return erc_(1); }
        
        inline const SVecD<2>& nrmc()  const { return nrmc_; }
        inline const Double_t& nrmcx() const { return nrmc_(0); }
        inline const Double_t& nrmcy() const { return nrmc_(1); }

        inline const SVecD<2>& divc()  const { return divc_; }
        inline const Double_t& divcx() const { return divc_(0); }
        inline const Double_t& divcy() const { return divc_(1); }

    protected :
        void clear();

    protected :
        Short_t seqID_;
        Short_t seqIDcx_;
        Short_t seqIDcy_;

        PartType type_; // particle type
        
        Detector dec_; // TRK TOF TRD RICH ECAL
        Short_t  lay_; // start from 0

        SVecO<3> side_coo_; // (x, y, z)
        SVecD<3> coo_;      // [cm] coord
        SVecD<2> erc_;      // [cm] error
        
        SVecD<2> nrmc_; // coord norm
        SVecD<2> divc_; // coord div
    
    public :
        enum class Orientation {
            kDownward = 0, kUpward = 1
        };

        static void Sort(std::vector<VirtualHitSt*>& hits, const Orientation& ortt = Orientation::kDownward) {
            if (hits.size() < 2) return;
            if (ortt == Orientation::kDownward) std::sort(hits.begin(), hits.end(), [](const VirtualHitSt* hit1, const VirtualHitSt* hit2) { return (hit1->cz() > hit2->cz()); } );
            else                                std::sort(hits.begin(), hits.end(), [](const VirtualHitSt* hit1, const VirtualHitSt* hit2) { return (hit1->cz() < hit2->cz()); } );
        }
};


class HitStTRK : public VirtualHitSt {
    public :
        static constexpr VirtualHitSt::Detector DEC = VirtualHitSt::Detector::TRK;

    public :
        HitStTRK(Bool_t scx = false, Bool_t scy = false, Short_t lay = 0) : VirtualHitSt(VirtualHitSt::Detector::TRK, lay, scx, scy), pdf_cx_(nullptr), pdf_cy_(nullptr), pdf_qx_(nullptr), pdf_qy_(nullptr) { clear(); }
        ~HitStTRK() { clear(); }
        
        Short_t set_seqID(Short_t seqID); 
        
        void cal(const PhySt& part);
        void set_type(const PartInfo& info = PartInfo(PartType::Proton));
        
        inline void set_nsr(Short_t nx, Short_t ny) {
            nsr_(0) = ((Numc::Compare(nx)>0) ? nx : Numc::ZERO<Short_t>);
            nsr_(1) = ((Numc::Compare(ny)>0) ? ny : Numc::ZERO<Short_t>);
        }
        // TODO: // rebuild template
        //inline void set_q(Double_t qx, Double_t qy) {
        //    side_q_(0) = (Numc::Compare(qx) > 0);
        //    side_q_(1) = (Numc::Compare(qy) > 0);
        //    q_(0) = (side_q_(0) ? qx : Numc::ZERO<>);
        //    q_(1) = (side_q_(1) ? qy : Numc::ZERO<>);
        //}
        
        inline const Short_t&  seqIDqx() const { return seqIDqx_; }
        inline const Short_t&  seqIDqy() const { return seqIDqy_; }

        inline const Bool_t& sqx() const { return side_q_(0); }
        inline const Bool_t& sqy() const { return side_q_(1); }

        inline const SVecD<2>& nrmq()  const { return nrmq_; }
        inline const Double_t& nrmqx() const { return nrmq_(0); }
        inline const Double_t& nrmqy() const { return nrmq_(1); }

        inline const SVecD<2>& divq()  const { return divq_; }
        inline const Double_t& divqx() const { return divq_(0); }
        inline const Double_t& divqy() const { return divq_(1); }

    protected :
        void clear();

    protected :
        Short_t seqIDqx_;
        Short_t seqIDqy_;
        
        SVecS<2> nsr_; // Number of strip
        
        SVecO<2> side_q_;
        SVecD<2> q_; // ADC

        SVecD<2> nrmq_; // q nrom
        SVecD<2> divq_; // q div

        MultiGaus* pdf_cx_;
        MultiGaus* pdf_cy_;
        IonEloss*  pdf_qx_;
        IonEloss*  pdf_qy_;
    
    protected :
        static MultiGaus PDF_Q01_CX_NN_;
        static MultiGaus PDF_Q01_CX_N1_;
        static MultiGaus PDF_Q01_CX_N2_;
        static MultiGaus PDF_Q01_CX_N3_;
        
        static MultiGaus PDF_Q01_CY_NN_;
        static MultiGaus PDF_Q01_CY_N1_;
        static MultiGaus PDF_Q01_CY_N2_;
        static MultiGaus PDF_Q01_CY_N3_;
        static MultiGaus PDF_Q01_CY_N4_;

        static IonEloss PDF_Q01_QX_;
        static IonEloss PDF_Q01_QY_;
};


class HitStTOF : public VirtualHitSt {
    public :
        static constexpr VirtualHitSt::Detector DEC = VirtualHitSt::Detector::TOF;
    
    public :
        HitStTOF(Bool_t scx = false, Bool_t scy = false, Short_t lay = 0) : VirtualHitSt(VirtualHitSt::Detector::TOF, lay, scx, scy), pdf_c_(nullptr), pdf_q_(nullptr), pdf_t_(nullptr) { clear(); }
        ~HitStTOF() { clear(); }
        
        Short_t set_seqID(Short_t seqID); 
        
        void cal(const PhySt& part);
        void set_type(const PartInfo& info = PartInfo(PartType::Proton));
        
        inline void set_q(Double_t q) {
            side_q_ = (Numc::Compare(q) > 0);
            q_      = (side_q_ ? q : Numc::ZERO<>);
        }
        
        inline void set_t(Double_t t) {
            side_t_ = (Numc::Compare(t) >= 0);
            t_      = (side_t_ ? t : Numc::ZERO<>);
        }

        inline Double_t t() const { return (t_ + OFFSET_T_); }

        inline const Short_t&  seqIDq() const { return seqIDq_; }
        inline const Short_t&  seqIDt() const { return seqIDt_; }
        
        inline const Bool_t& sq() const { return side_q_; }
        inline const Bool_t& st() const { return side_t_; }

        inline const Double_t& nrmq() const { return nrmq_; }
        inline const Double_t& divq() const { return divq_; }

        inline const Double_t& nrmt() const { return nrmt_; }
        inline const Double_t& divt() const { return divt_; }

    protected :
        void clear();

    protected :
        Short_t seqIDq_;
        Short_t seqIDt_;
        
        Bool_t   side_q_;
        Double_t q_; // Q
        
        Bool_t   side_t_;
        Double_t t_; // T [cm]

        Double_t nrmq_; // Q nrom
        Double_t divq_; // Q div
        
        Double_t nrmt_; // T nrom
        Double_t divt_; // T div

        MultiGaus* pdf_c_;
        IonEloss*  pdf_q_;
        MultiGaus* pdf_t_;
    
    protected :
        static MultiGaus PDF_Q01_C_;
        static IonEloss  PDF_Q01_Q_;
        static MultiGaus PDF_Q01_T_;

    public :
        static constexpr Double_t TRANS_NS_TO_CM = 2.99792458e+01; // [ns] -> [cm]
        static void SetOffsetTime(Double_t offset_t = Numc::ZERO<>) { OFFSET_T_ = offset_t; }
        static void SetOffsetPath(Double_t offset_s = Numc::ZERO<>) { OFFSET_S_ = offset_s; }
        static const Double_t& OffsetTime() { return OFFSET_T_; }
        static const Double_t& OffsetPath() { return OFFSET_S_; }
    
    protected :
        static Double_t OFFSET_T_; // move TOF to particle time (FIRST-TOF PART-TOF)
        static Double_t OFFSET_S_; // move TOF to particle path (FIRST-TOF)
};

Double_t HitStTOF::OFFSET_T_ = Numc::ZERO<>;
Double_t HitStTOF::OFFSET_S_ = Numc::ZERO<>;


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


} // namesapce TrackSys


#endif // __TRACKLibs_HitSt_H__
