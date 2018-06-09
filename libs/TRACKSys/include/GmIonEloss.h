#ifndef __TRACKLibs_GmIonEloss_H__
#define __TRACKLibs_GmIonEloss_H__


namespace TrackSys {

class GmIonEloss {
    public :
        GmIonEloss(
            const std::array<long double, 8>& lg_wgt, const std::array<long double, 8>& lg_kpa, const std::array<long double, 8>& lg_mpv, const std::array<long double, 8>& lg_sgm,
            const std::array<long double, 8>& gm_wgt, const std::array<long double, 8>& gm_alpha, const std::array<long double, 8>& gm_beta, const std::array<long double, 8>& gm_eftm, const std::array<long double, 8>& gm_efts
        ) : lg_wgt_(lg_wgt), lg_kpa_(lg_kpa), lg_mpv_(lg_mpv), lg_sgm_(lg_sgm), gm_wgt_(gm_wgt), gm_alpha_(gm_alpha), gm_beta_(gm_beta), gm_eftm_(gm_eftm), gm_efts_(gm_efts) {}
        ~GmIonEloss() {}
        
    protected :

    private :
        std::array<long double, 8> lg_wgt_;
        std::array<long double, 8> lg_kpa_;
        std::array<long double, 8> lg_mpv_;
        std::array<long double, 8> lg_sgm_;

        std::array<long double, 8> gm_wgt_;
        std::array<long double, 8> gm_alpha_;
        std::array<long double, 8> gm_beta_;
        std::array<long double, 8> gm_eftm_;
        std::array<long double, 8> gm_efts_;
};

} // namesapce TrackSys


#endif // __TRACKLibs_GmIonEloss_H__
