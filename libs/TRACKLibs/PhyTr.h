#ifndef __TRACKLibs_PhyTr_H__
#define __TRACKLibs_PhyTr_H__


namespace TrackSys {


class PhyTr {
    public :
        enum class Orientation {
            kDownward = 0, kUpward = 1
        };

        inline static Bool_t HitCheck(const std::vector<HitSt>& hits);
        inline static std::tuple<Short_t, std::vector<Short_t>, std::vector<Short_t>, std::vector<Short_t>> HitSort(std::vector<HitSt>& hits, Orientation ortt = Orientation::kDownward);
        
    public :
        PhyTr(const std::vector<HitSt>& hits, const PartType& type = PartType::Proton, const Orientation& ortt = Orientation::kDownward, Bool_t sw_mscat = PhyArg::OptMscat(), Bool_t sw_eloss = PhyArg::OptEloss());
        ~PhyTr() { clear(); }

        void print() const;

        Bool_t fit();

        inline const PhySt& part() const { return part_; }
        inline const Double_t& nchi() const { return nchi_; }

    protected :
        void clear();
        
        // Fitting
        Bool_t fit_analysis();
        Bool_t fit_simple();
        Bool_t fit_physics();
        
    private :
        Bool_t sw_mscat_;
        Bool_t sw_eloss_;
        
        PartType             type_;
        Orientation          ortt_;
        Short_t              nseq_; // tot seq
        std::vector<Short_t> seqx_; // x -> seq
        std::vector<Short_t> seqy_; // y -> seq
        std::vector<Short_t> maps_; // seq -> (x=0 or y=1) * nhit + hitID
        std::vector<HitSt>   hits_;

        // Fitting
        Bool_t   succ_;
        PhySt    part_;
        Double_t nchi_;
        Double_t ndf_;

    protected :
        // Dimension
        static constexpr Int_t DIM_G = 5;
        static constexpr Int_t DIM_L = 4;

        // Minimization (Levenberg-Marquardt Method)
        static constexpr Int_t    LMTL_ITER = 3;
        static constexpr Int_t    LMTU_ITER = 30;
        
        static constexpr Double_t LAMBDA0 = 1.0e-2;
        static constexpr Double_t LAMBDA_DN_FAC = 7.0;
        static constexpr Double_t LAMBDA_UP_FAC = 9.0;
        static constexpr Double_t LMTL_LAMBDA = 1.0e-4;
        static constexpr Double_t LMTU_LAMBDA = 1.0e+3;
        
        static constexpr Double_t CONVG_EPSILON   = 1.0e-3;
        static constexpr Double_t CONVG_TOLERANCE = 1.0e-2;
};


} // namespace TrackSys


#endif // __TRACKLibs_PhyTr_H__
