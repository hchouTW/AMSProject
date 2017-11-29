#ifndef __TRACKLibs_PhyTr_H__
#define __TRACKLibs_PhyTr_H__


namespace TrackSys {

class PhyTr {
    public :
        enum class Orientation {
            kDownward = 0, kUpward = 1
        };

        inline static Short_t HitSort(std::vector<HitSt>& hits, Orientation ortt = Orientation::kDownward);

        inline static Bool_t HitCheck(const std::vector<HitSt>& hits);
        
    public :
        PhyTr(const std::vector<HitSt>& hits, const PartType& type = PartType::Proton, const Orientation& ortt = Orientation::kDownward);
        ~PhyTr() { clear(); }

        void print() const;

        Bool_t fit();

        inline const PhySt& part() const { return part_; }

    protected :
        void   clear();
        Bool_t fit_analysis();
        Bool_t fit_simple();
        Bool_t fit_semi_simple();
        Bool_t fit_physics();
        
    private :
        Bool_t              succ_;
        PartType            type_;
        Orientation         ortt_;
        PhySt               part_;
        std::vector<HitSt>  hits_;
        Short_t             nseq_;
        Double_t nchi_;
        Double_t ndf_;

    protected :
        // Minimization
        static constexpr Int_t    LMTL_ITER = 3;
        static constexpr Int_t    LMTM_ITER = 5;
        static constexpr Int_t    LMTU_ITER = 100;
        static constexpr Double_t CONVG_EPSILON   = 5.0e-2;
        static constexpr Double_t CONVG_TOLERANCE = 1.0e-2;
       
        // Dimension
        static constexpr Int_t DIM_G = 5;
        static constexpr Int_t DIM_L = 4;
};








/*
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
        
        inline const SVecI<2>& nhit() const { return nhit_; }
        inline const std::vector<HitSt>& hits() const { return hits_; }
        inline const std::vector<PhySt>& phys() const { return phys_; }

    protected :
        void clear();

        Bool_t check(const std::vector<HitSt>& hits);

        Bool_t fit();
        Bool_t fit_analysis();
        Bool_t fit_simple();
        Bool_t fit_semi_simple();
        Bool_t fit_physics();

    private :
        Bool_t              succ_;
        Orientation         ortt_;
        SVecI<2>            nhit_;
        std::vector<HitSt>  hits_;
        PartType            type_;
        PhySt               part_;
        
        std::vector<PhySt>  phys_;
        Double_t nchi_;
        Double_t ndf_;

    private:
        // Number of Hit Requirement
        static constexpr Int_t LMTL_NHIT_X = 3;
        static constexpr Int_t LMTL_NHIT_Y = 4;

        // Minimization
        static constexpr Int_t    LMTL_ITER = 3;
        static constexpr Int_t    LMTU_ITER = 20 * 100;
        static constexpr Double_t CONVG_EPSILON   = 1.0e-3;
        static constexpr Double_t CONVG_TOLERANCE = 1.0e-2;
       
        // Dimension
        static constexpr Int_t DIM_G = 5;
        static constexpr Int_t DIM_L = 4;
};
*/

} // namespace TrackSys


#endif // __TRACKLibs_PhyTr_H__
