#ifndef __AmsRich_H__
#define __AmsRich_H__

#include "root.h"
#include "point.h"
#include "richidOff.h"
#include "richradidOff.h"
#include "richtrrecOff.h"
#include "richdbcOff.h"
#include "richrecOff.h"
#include "RichBeta.h"


class AmsRichHit {
    public :
        AmsRichHit() { clear(); }
        AmsRichHit(RichHitR* hit, double dbeta, double rbetaA, double rbetaB);
        ~AmsRichHit() {}

        inline const bool& status() const { return status_; }

        inline const short& chann() const { return chann_; }
        inline const short& pmtid() const { return pmtid_; }
        inline const short& pixel() const { return pixel_; }
        
        inline const short&  type()   const { return type_; }
        inline const double& dbeta()  const { return dbeta_; }
        inline const double& rbetaA() const { return rbetaA_; }
        inline const double& rbetaB() const { return rbetaB_; }
        inline const double& npe()    const { return npe_; }
        
        inline const double& cx() const { return cx_; }
        inline const double& cy() const { return cy_; }
        inline const double& cz() const { return cz_; }
        
    protected :
        void clear();

    protected :
        bool   status_;
        short  chann_;
        short  pmtid_;
        short  pixel_;

        short  type_;
        double dbeta_;
        double rbetaA_;
        double rbetaB_;
        double npe_;
        double cx_;
        double cy_;
        double cz_;
        double dist_;

    protected :
        RichHitR* hit_;
        
    private :
        static constexpr short NUM_CHANN_IN_1D_PMT = 4;  // 4
        static constexpr short NUM_CHANN_IN_2D_PMT = 16; // 4 x 4
};

struct AmsRichHit_sort {
    bool operator() (const AmsRichHit& hit1, const AmsRichHit& hit2) {
        if      (hit1.chann() < hit2.chann()) return true;
        else if (hit1.chann() > hit2.chann()) return false;
        return false;
    }
};


class AmsRich {
    public :
        AmsRich() { clear(); }
        AmsRich(AMSEventR* event);
        ~AmsRich() {}

        inline const bool& status() const { return status_; }

        inline const short& kind() const { return kind_; }
        inline const short& tile() const { return tile_; }
       
        inline bool is_agl() const { return (kind_ == KIND_AGL); }
        inline bool is_naf() const { return (kind_ == KIND_NAF); }

        inline const double& index()  const { return index_; }
        inline const double& dist()   const { return dist_; }

        inline const double& locx() const { return locx_; }
        inline const double& locy() const { return locy_; }
        inline const double& loctha() const { return loctha_; }
        inline const double& locphi() const { return locphi_; }

        inline const AMSPoint& dirp() const { return dirp_; }
        inline const AMSDir&   dird() const { return dird_; }
        
        inline const AMSPoint& refp() const { return refp_; }
        inline const AMSDir&   refd() const { return refd_; }
        
        inline const AMSPoint& radp() const { return radp_; }
        inline const AMSDir&   radd() const { return radd_; }
        inline const AMSPoint& pmtp() const { return pmtp_; }
        
        inline const bool& is_good_geom() const { return is_good_geom_; }
        inline const bool& is_bad_tile() const { return is_bad_tile_; }
        
        inline const double& beta_crr() const { return beta_crr_; }

        inline const std::vector<AmsRichHit>& hits() const { return hits_; }

    protected :
        void clear();
        bool build(RichOffline::TrTrack track);
    
    protected :
        bool     status_;
        short    kind_;
        short    tile_;
        double   index_;
        double   dist_;
        double   locx_;
        double   locy_;
        double   loctha_;
        double   locphi_;

        AMSPoint dirp_;
        AMSDir   dird_;
        AMSPoint refp_;
        AMSDir   refd_;

        AMSPoint radp_; // charged particle in RAD plane (AMS coord)
        AMSDir   radd_; // charged particle in RAD plane (AMS coord)
        AMSPoint pmtp_; // charged particle in PMT plane (AMS coord)

        bool is_good_geom_;
        bool is_bad_tile_;

        double beta_crr_;

        std::vector<AmsRichHit> hits_;

    protected :
        AMSEventR* event_;
        TrTrackR*  trtk_;
    
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
        
        static constexpr std::array<double, 2> RAD_HEIGHT { 2.5, 0.5 }; // AGL, NAF
        static constexpr std::array<double, 2> RAD_BOUNDARY { 19.0, 17.0 }; // AGL, NAF
        static constexpr std::array<short, 7>  BAD_TILE_INDEX { 3, 7, 12, 20, 87, 100, 108 };
};


#endif // __AmsRich_H__
