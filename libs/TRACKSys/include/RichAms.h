#if defined(_PGTRACK_) || defined(__ROOTSHAREDLIBRARY__)
#ifndef __TRACKLibs_RichAms_H__
#define __TRACKLibs_RichAms_H__


#include "root.h"
#include "point.h"
#include "richidOff.h"
#include "richdbcOff.h"
#include "richrecOff.h"
#include "RichBeta.h"


namespace TrackSys {


class RichHitAms {
    public :
        RichHitAms() { clear(); }
        RichHitAms(RichHitR* hit, double dbta, double rbta, double dist);
        ~RichHitAms() {}

        inline const bool& status() const { return status_; }

        inline const short& chann() const { return chann_; }
        inline const short& pmtid() const { return pmtid_; }
        
        inline const short&  type() const { return type_; }
        inline const double& dbta() const { return dbta_; }
        inline const double& rbta() const { return rbta_; }
        inline const double& npe()  const { return npe_; }
        
        inline const double& cx() const { return cx_; }
        inline const double& cy() const { return cy_; }
        inline const double& cz() const { return cz_; }
        
        inline const double& dist() const { return dist_; }
        
        CherenkovHit trans() const;

    protected :
        void clear();

    protected :
        bool   status_;
        short  chann_;
        short  pmtid_;
        short  type_;
        double dbta_;
        double rbta_;
        double npe_;
        double cx_;
        double cy_;
        double cz_;
        double dist_;

    protected :
        RichHitR* hit_;
};

struct RichHitAms_sort {
    bool operator() (const RichHitAms& hit1, const RichHitAms& hit2) {
        if      (hit1.chann() < hit2.chann()) return true;
        else if (hit1.chann() > hit2.chann()) return false;
        return false;
    }
};


class RichAms {
    public :
        inline static void LoadEnv(bool is_iss = false) {
            if (is_iss) {
                RichRingR::reloadRunTag = true;
                RichRingR::useTemperatureCorrections = true;
                RichRingR::loadChargeUniformityCorrection = true;
		        RichRingR::setBetaCorrection(RichRingR::fullUniformityCorrection);
            }
            else {
                RichRingR::setBetaCorrection(RichRingR::noCorrection);
            }
        }

    public :
        RichAms() { clear(); }
        RichAms(AMSEventR* event, TrTrackR* trtk = nullptr);
        ~RichAms() {}
        
        inline const bool& status() const { return status_; }

        inline const short& kind() const { return kind_; }
        inline const short& tile() const { return tile_; }
       
        inline bool is_agl() const { return (kind_ == KIND_AGL); }
        inline bool is_naf() const { return (kind_ == KIND_NAF); }

        inline const double& index() const { return index_; }
        inline const double& dist()  const { return dist_; }

        inline const AMSPoint& dirp() const { return dirp_; }
        inline const AMSDir&   dird() const { return dird_; }
        
        inline const AMSPoint& refp() const { return refp_; }
        inline const AMSDir&   refd() const { return refd_; }
        
        inline const AMSPoint& radp() const { return radp_; }
        inline const AMSDir&   radd() const { return radd_; }
        inline const AMSPoint& pmtp() const { return pmtp_; }
        
        inline const double& bta_crr() const { return bta_crr_; }

        inline const std::vector<RichHitAms>&   rhhits() const { return rhhits_; }
        inline const std::vector<CherenkovHit>& chhits() const { return chhits_; }
        
        inline const Sys::HrsStopwatch& timer() const { return timer_; }

    protected :
        void clear();
        bool build(RichOffline::TrTrack track);
    
    protected :
        Sys::HrsStopwatch timer_;

    protected :
        bool     status_;
        short    kind_;
        short    tile_;
        double   index_;
        double   dist_;
        AMSPoint dirp_;
        AMSDir   dird_;
        AMSPoint refp_;
        AMSDir   refd_;

        AMSPoint radp_; // charged particle in RAD plane (AMS coord)
        AMSPoint radd_; // charged particle in RAD plane (AMS coord)
        AMSPoint pmtp_; // charged particle in PMT plane (AMS coord)

        double bta_crr_;

        std::vector<RichHitAms>   rhhits_;
        std::vector<CherenkovHit> chhits_;

    protected :
        AMSEventR* event_;
        TrTrackR*  trtk_;

    private :
        static constexpr short KIND_EMPTY = 0;
        static constexpr short KIND_AGL   = 1;
        static constexpr short KIND_NAF   = 2;

        static constexpr double PMT_CZ = -121.89;
};


} // namespace TrackSys


#endif // __TRACKLibs_RichAms_H__
#endif // _PGTRACK_ __ROOTSHAREDLIBRARY__ 
