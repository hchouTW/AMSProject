#ifndef __TRACKLibs_PhyTr_H__
#define __TRACKLibs_PhyTr_H__


namespace TrackSys {


class PhyTr {
    public :
        enum class Orientation {
            kDownward = 0, kUpward = 1
        };

    public :
        PhyTr() { clear(); }
        PhyTr(const std::vector<HitSt>& hits, const PartType& type = PartType::Proton, const Orientation& ortt = Orientation::kDownward);
        ~PhyTr() {}

        PhySt prop_to_z(const Double_t zcoo) { return PhySt(); }

        inline const Bool_t& exist() const { return succ_; }
        inline const PhySt& part() const { return part_; }
        
        inline const std::vector<HitSt>& hits() const { return hits_; }
        inline const std::vector<PhySt>& phys() const { return phys_; }
        inline const std::vector<MatArg>& marg() const { return marg_; }

    protected :
        void clear();

        Bool_t check(const std::vector<HitSt>& hits);

        Bool_t fit();
        Bool_t fit_analysis();
        Bool_t fit_simple();
        Bool_t fit_physics();

    private :
        Bool_t              succ_;
        Orientation         ortt_;
        std::vector<HitSt>  hits_;
        PartType            type_;
        PhySt               part_;
        
        std::vector<MatArg> marg_;
        std::vector<PhySt>  phys_;
        Double_t nchi_;
        Double_t ndf_;

    private:
        // Number of Hit Requirement
        static constexpr Int_t LMTL_HIT_X = 3;
        static constexpr Int_t LMTL_HIT_Y = 4;

        // Minimization
        static constexpr Int_t    LMTL_ITER = 3;
        static constexpr Int_t    LMTU_ITER = 10;
        static constexpr Double_t CONVG_EPSILON   = 1.0e-3;
        static constexpr Double_t CONVG_TOLERANCE = 5.0e-2;
};


} // namespace TrackSys


#endif // __TRACKLibs_PhyTr_H__
