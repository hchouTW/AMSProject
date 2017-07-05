#ifndef __TRACKLibs_PhyTr_H__
#define __TRACKLibs_PhyTr_H__


namespace TrackSys {


class PhyTr {
    public :
        enum class Orientation {
            kUpward = 0, kDownward = 1
        };

    public :
        PhyTr(const std::vector<HitSt>& hits, const PartType& type = PartType::Proton, const Orientation& ortt = Orientation::kDownward) {}
        ~PhyTr() {}

        inline const Bool_t& exist() const { return succ_; }
        inline const PhySt& part() const { return part_; }

        inline const std::vector<HitSt>& hits() const { return hits_; }
        inline const std::vector<PhySt>& phys() const { return phys_; }
        inline const std::vector<MatArg>& marg() const { return marg_; }

    protected :
        Bool_t fit() { return true; }

        Bool_t fit_analyis() { return true; }
        Bool_t fit_simple() { return true; }
        Bool_t fit_physics() { return true; }

    private :
        Bool_t              succ_;
        PartType            type_;
        Orientation         ortt_;
        PhySt               part_;
        std::vector<HitSt>  hits_;
        std::vector<PhySt>  phys_;
        std::vector<MatArg> marg_;

        Double_t chi_;
        Double_t ndf_;

    private:
        // Number of Hit Requirement
        static constexpr Int_t LMTL_HIT_X = 3;
        static constexpr Int_t LMTL_HIT_Y = 4;

        // Minimization
};


} // namespace TrackSys


#endif // __TRACKLibs_PhyTr_H__
