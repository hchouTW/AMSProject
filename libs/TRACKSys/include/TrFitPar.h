#ifndef __TRACKLibs_TrFitPar_H__
#define __TRACKLibs_TrFitPar_H__


namespace TrackSys {


class TrFitPar {
    public :
        enum class Orientation {
            kDownward = 0, kUpward = 1
        };

    public :
        TrFitPar& operator=(const TrFitPar& rhs);
        TrFitPar(const TrFitPar& fitPar) { *this = fitPar; }
        
        TrFitPar(const PartInfo& info = PartInfo(PartType::Proton), const Orientation& ortt = Orientation::kDownward, const Bool_t& sw_mscat = PhyArg::OptMscat(), const Bool_t& sw_eloss = PhyArg::OptEloss());
        ~TrFitPar() { TrFitPar::clear(); }

    public :
        inline Bool_t check() { return check_hits(); }

        inline void add_hit(const HitStTRK&  hit) { hits_TRK_.push_back(hit);  zero(); }
        inline void add_hit(const HitStTOF&  hit) { hits_TOF_.push_back(hit);  zero(); }
        inline void add_hit(const HitStRICH& hit) { hits_RICH_.push_back(hit); zero(); }
        inline void add_hit(const HitStTRD&  hit) { hits_TRD_.push_back(hit);  zero(); }
        
        inline void add_hit(const std::vector<HitStTRK>&  hits) { hits_TRK_.insert(hits_TRK_.end(), hits.begin(), hits.end()); zero(); }
        inline void add_hit(const std::vector<HitStTOF>&  hits) { hits_TOF_.insert(hits_TOF_.end(), hits.begin(), hits.end()); zero(); }
        inline void add_hit(const std::vector<HitStRICH>& hits) { hits_RICH_.insert(hits_RICH_.end(), hits.begin(), hits.end()); zero(); }
        inline void add_hit(const std::vector<HitStTRD>&  hits) { hits_TRD_.insert(hits_TRD_.end(), hits.begin(), hits.end()); zero(); }
      
        inline const Bool_t& sw_mscat() const { return sw_mscat_; }
        inline const Bool_t& sw_eloss() const { return sw_eloss_; }
        inline const PartInfo&    info() const { return info_; }
        inline const Orientation& ortt() const { return ortt_; }

        inline const Short_t& onlycx_nseq() const { return onlycx_nseq_; }
        inline const Short_t& onlycy_nseq() const { return onlycy_nseq_; }
        inline const Short_t& onlyc_nseq() const { return onlyc_nseq_; }
        inline const Short_t& nseq() const { return nseq_; }
        inline const Short_t& nseg() const { return nseg_; }

        inline Short_t nhit() const { return hits_.size(); }
        inline const std::vector<VirtualHitSt*>& hits() const { return hits_; }
        inline const VirtualHitSt* hit(Int_t idx) const { return ((idx<0 || idx>=hits_.size()) ? nullptr : hits_.at(idx)); }
        
        inline const std::vector<HitStTRK>&  hitsTRK()  const { return hits_TRK_; }
        inline const std::vector<HitStTOF>&  hitsTOF()  const { return hits_TOF_; }
        inline const std::vector<HitStRICH>& hitsRICH() const { return hits_RICH_; }
        inline const std::vector<HitStTRD>&  hitsTRD()  const { return hits_TRD_; }

    protected :
        void zero();
        void clear();
        
        Bool_t sort_hits();
        Bool_t check_hits();

        Bool_t survival_test_and_modify(PhySt& part, Bool_t with_eloss = false);

    protected :
        Bool_t      sw_mscat_;
        Bool_t      sw_eloss_;
        PartInfo    info_;
        Orientation ortt_;

        std::vector<VirtualHitSt*> hits_;
        std::vector<HitStTRK>      hits_TRK_;
        std::vector<HitStTOF>      hits_TOF_;
        std::vector<HitStRICH>     hits_RICH_;
        std::vector<HitStTRD>      hits_TRD_;

        Short_t onlycx_nseq_;
        Short_t onlycy_nseq_;
        Short_t onlyc_nseq_;
        Short_t nseq_;
        Short_t nseg_;

        Short_t nmes_;
        Short_t nmes_cx_;
        Short_t nmes_cy_;
        Short_t nmes_ib_;
        Short_t nmes_TRKq_;
        Short_t nmes_TOFq_;
        Short_t nmes_TOFt_;
        Short_t nmes_RICHib_;
        Short_t nmes_TRDel_;

    private :
        Bool_t  is_check_;

    protected :
        // Number of Hit Requirement
        static constexpr Short_t LMTN_CX = 2;
        static constexpr Short_t LMTN_CY = 3;
        static constexpr Short_t LMTN_TOF_T = 1;
        
        // Survival Test
        static constexpr Short_t  SURVIVAL_LMTN = 15;
        static constexpr Double_t SURVIVAL_FACT = 0.85;
        static constexpr Double_t SURVIVAL_BETA = 0.30;
        
        // Limit of 1/beta
        static constexpr Double_t LMTL_IBTA = 1.000000001;
        static constexpr Double_t LMTU_IBTA = 100.;
};


} // namespace TrackSys


#endif // __TRACKLibs_TrFitPar_H__
