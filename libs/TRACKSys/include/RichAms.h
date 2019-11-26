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


class RichObjAms {
    public :
        RichObjAms() { init(); }
        ~RichObjAms() {}

        inline void init() {
            status = false;
            kind = -1;
            tile = -1;

            index = 0;
            dist = 0;

            is_good_geom = false;
            is_bad_tile = false;
            
            bta_crr = 1.0;

            std::fill_n(radp, 3, 0);
            std::fill_n(radd, 3, 0);
            std::fill_n(pmtp, 3, 0);

            nstn = 0;
            ncld = 0;
            ntmr = 0;
            ngst = 0; 

            nhit_ttl = 0;
            nhit_stn = 0;
            nhit_cld = 0;
            nhit_tmr = 0;
            nhit_gst = 0;
            nhit_oth = 0;
            
            npe_ttl = 0;
            npe_stn = 0;
            npe_cld = 0;
            npe_tmr = 0;
            npe_gst = 0;
            npe_oth = 0;

            stn_status = false;
            stn_nhit = 0;
            stn_npmt = 0;
            stn_cx = 0;
            stn_cy = 0;
            stn_dist = 0;
            stn_npe = 0;
            stn_cnt = 0;
            stn_nchi = 0;
            stn_chic = 0;

            cld_status = false;
            cld_nhit = 0;
            cld_npmt = 0;
            cld_border = 0;
            cld_trace = 0;
            cld_accuracy = 0;
            cld_uniform = 0;
            cld_crrch = 0;
            cld_beta = 0;
            cld_cbta = 0;
            cld_npe = 0;
            cld_expnpe = 0;
            cld_cnt = 0;
            cld_nchi = 0;
            cld_misjudge = 0;

            cldhit_beta.clear();
            cldhit_npe.clear();
            cldhit_cx.clear();
            cldhit_cy.clear();

            nhit = 0;
            hit_chann.clear();
            hit_pmtid.clear();
            hit_cls.clear();
            hit_mode.clear();
            hit_beta.clear();
            hit_npe.clear();
            hit_cx.clear();
            hit_cy.clear();

            cpuTime = 0;
        }

    public :
        bool  status;
        short kind;
        short tile;
        float index;
        float dist;

        bool is_good_geom;
        bool is_bad_tile;
        
        float bta_crr;

        float radp[3];
        float radd[3];
        float pmtp[3];

        short nstn;
        short ncld;
        short ntmr;
        short ngst;

        short nhit_ttl;
        short nhit_stn;
        short nhit_cld;
        short nhit_tmr;
        short nhit_gst;
        short nhit_oth;
        
        float npe_ttl;
        float npe_stn;
        float npe_cld;
        float npe_tmr;
        float npe_gst;
        float npe_oth;

        bool  stn_status;
        short stn_nhit;
        short stn_npmt;
        float stn_cx;
        float stn_cy;
        float stn_dist;
        float stn_npe;
        float stn_cnt;
        float stn_nchi;
        float stn_chic;

        bool  cld_status;
        short cld_nhit;
        short cld_npmt;
        float cld_border;
        float cld_trace;
        float cld_accuracy;
        float cld_uniform;
        float cld_crrch;
        float cld_beta;
        float cld_cbta;
        float cld_npe;
        float cld_expnpe;
        float cld_cnt;
        float cld_nchi;
        float cld_misjudge;
        
        std::vector<float> cldhit_beta;
        std::vector<float> cldhit_npe;
        std::vector<float> cldhit_cx;
        std::vector<float> cldhit_cy;

        short              nhit;
        std::vector<short> hit_chann;
        std::vector<short> hit_pmtid;
        std::vector<short> hit_cls;
        std::vector<short> hit_mode;
        std::vector<float> hit_beta;
        std::vector<float> hit_npe;
        std::vector<float> hit_cx;
        std::vector<float> hit_cy;

        float cpuTime;
};


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

        static std::vector<std::array<double, 2>> RayTrace(const std::array<double, 6>& part, double cbta, double index, double height, short tile);

    public :
        RichAms() { clear(); }
        RichAms(AMSEventR* event, TrTrackR* trtk = nullptr);
        ~RichAms() {}

        CherenkovFit fit();
        RichObjAms get_obj_by_fit(int zin = 1);

        std::array<double, 5> cal_trace(double cbta = 1.0, const CherenkovCloud* cloud = nullptr) const;

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
        
        inline const bool& is_good_geom() const { return is_good_geom_; }
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

        bool is_good_geom_;
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
