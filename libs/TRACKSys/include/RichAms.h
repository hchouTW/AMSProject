#if defined(_PGTRACK_) || defined(__ROOTSHAREDLIBRARY__)
#ifndef __TRACKLibs_RichAms_H__
#define __TRACKLibs_RichAms_H__


#include "root.h"
#include "point.h"
#include "richidOff.h"
#include "richradidOff.h"
#include "richtrrecOff.h"
#include "richdbcOff.h"
#include "richrecOff.h"
#include "RichBeta.h"


namespace TrackSys {


class RichHitAms {
    public :
        RichHitAms() { clear(); }
        RichHitAms(RichHitR* hit, double dbta, double rbtaA, double rbtaB, double dist);
        ~RichHitAms() {}

        inline const bool& status() const { return status_; }

        inline const short& chann() const { return chann_; }
        inline const short& pmtid() const { return pmtid_; }
        inline const short& pixel() const { return pixel_; }
        inline const short& locid(int it)  const { return ((it == 0) ? locid_[0] : locid_[1]); }
        
        inline const short&  type()  const { return type_; }
        inline const double& dbta()  const { return dbta_; }
        inline const double& rbtaA() const { return rbtaA_; }
        inline const double& rbtaB() const { return rbtaB_; }
        inline const double& npe()   const { return npe_; }
        
        inline const double& cx() const { return cx_; }
        inline const double& cy() const { return cy_; }
        inline const double& cz() const { return cz_; }
        
        inline const double& dist() const { return dist_; }
       
        inline const bool& cross() const { return cross_; }

        CherenkovHit trans() const;

    protected :
        void clear();

    protected :
        bool   status_;
        short  chann_;
        short  pmtid_;
        short  pixel_;
        short  locid_[2];

        short  type_;
        double dbta_;
        double rbtaA_;
        double rbtaB_;
        double npe_;
        double cx_;
        double cy_;
        double cz_;
        double dist_;

        bool   cross_;

    protected :
        RichHitR* hit_;
        
    private :
        static constexpr short NUM_CHANN_IN_1D_PMT = 4;  // 4
        static constexpr short NUM_CHANN_IN_2D_PMT = 16; // 4 x 4
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

        CherenkovFit fit() const;
       
        std::array<double, 4> cal_trace(double cbta = 1.0, const CherenkovCloud* cloud = nullptr) const;

        inline const bool& status() const { return status_; }

        inline const short& kind() const { return kind_; }
        inline const short& tile() const { return tile_; }
       
        inline bool is_agl() const { return (kind_ == KIND_AGL); }
        inline bool is_naf() const { return (kind_ == KIND_NAF); }

        inline const double& index()  const { return index_; }
        inline const double& dist()   const { return dist_; }

        inline const AMSPoint& dirp() const { return dirp_; }
        inline const AMSDir&   dird() const { return dird_; }
        
        inline const AMSPoint& refp() const { return refp_; }
        inline const AMSDir&   refd() const { return refd_; }
        
        inline const AMSPoint& radp() const { return radp_; }
        inline const AMSDir&   radd() const { return radd_; }
        inline const AMSPoint& pmtp() const { return pmtp_; }
        inline const AMSPoint& rayp() const { return rayp_; }
        inline const AMSDir&   rayd() const { return rayd_; }
        
        inline const bool& is_good_in_geometry() const { return is_good_in_geometry_; }
        inline const bool& is_bad_tile() const { return is_bad_tile_; }
        
        inline const bool& is_in_pmt_plane() const { return is_in_pmt_plane_; }

        inline const double& npe_col() const { return npe_col_; }
        inline const double& bta_crr() const { return bta_crr_; }

        inline const std::vector<RichHitAms>&   rhhits() const { return rhhits_; }
        inline const std::vector<CherenkovHit>& chhits() const { return chhits_; }

        inline const Sys::HrsStopwatch& timer() const { return timer_; }

    protected :
        void clear();
        bool build(RichOffline::TrTrack track);
    
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
        AMSDir   radd_; // charged particle in RAD plane (AMS coord)
        AMSPoint pmtp_; // charged particle in PMT plane (AMS coord)
        AMSPoint rayp_; // ray of charged particle in PMT plane (AMS coord)
        AMSDir   rayd_; // ray of charged particle in PMT plane (AMS coord)

        bool is_good_in_geometry_;
        bool is_bad_tile_;

        bool is_in_pmt_plane_;

        double npe_col_;
        double bta_crr_;

        std::vector<RichHitAms>   rhhits_;
        std::vector<CherenkovHit> chhits_;

    protected :
        AMSEventR* event_;
        TrTrackR*  trtk_;
    
    protected :
        Sys::HrsStopwatch timer_;

    private :
        static constexpr short KIND_EMPTY = 0;
        static constexpr short KIND_AGL   = 1;
        static constexpr short KIND_NAF   = 2;

        static constexpr double WIDTH_CELL = 0.85;
        static constexpr double WIDTH_PMT  = 3.40;
        
        static constexpr double PMT_CZ = -121.89;

        static constexpr double EXTERNAL_RAD_RADIUS = 58.0;
        static constexpr double MIRROR_TOP_RADIUS = 60.10;
        static constexpr double MIRROR_BTM_RADIUS = 67.00;
        static constexpr double MIRROR_HEIGHT     = 46.32;
        static constexpr std::array<double, 2> PMT_HOLE { 31.9, 32.15 };
        
        //static constexpr std::array<double, 2> RAD_HEIGHT { 2.5, 0.5 }; // AGL, NAF
        static constexpr std::array<double, 2> RAD_BOUNDARY { 19.0, 17.0 }; // AGL, NAF
        //static constexpr std::array<short, 7>  BAD_TILE_INDEX { 3, 7, 12, 20, 87, 100, 108 };

    // trace
    private :
        static constexpr short TRACE_NSET = 5;
        static constexpr short TRACE_NPHI = 780;

    // interface of clustering, fitting
    private :
        static constexpr double AGL_BETA_WIDTH = 2.0e-03;
        static constexpr double NAF_BETA_WIDTH = 7.0e-03;
        static MultiGaus AGL_BETA_PDF;
        static MultiGaus NAF_BETA_PDF;
};

MultiGaus RichAms::AGL_BETA_PDF(Robust::Option(Robust::Opt::ON, 3.0L, 0.5L), 7.53836973171249758e-01, 1.87228e-03, 2.46163026828750242e-01, 3.72749e-03);
MultiGaus RichAms::NAF_BETA_PDF(Robust::Option(Robust::Opt::ON, 3.0L, 0.5L), 7.00000e-03);


} // namespace TrackSys


#endif // __TRACKLibs_RichAms_H__
#endif // _PGTRACK_ __ROOTSHAREDLIBRARY__ 
