#ifndef __TRACKLibs_CherenkovMeas_H__
#define __TRACKLibs_CherenkovMeas_H__


namespace TrackSys {

class CherenkovHit {
    public :
        enum Cluster { stone = 0, cloud = 1, tumor = 2, ghost = 3, emery = 4, other = 5 };

    public :
        CherenkovHit() { clear(); }
        ~CherenkovHit() { clear(); }
        
        CherenkovHit(short chann, short pmtid = -1, double dbta = -1.0, double rbtaA = -1.0, double rbtaB = -1.0, double npe = -1.0, double cx = 0.0, double cy = 0.0) : CherenkovHit() {
            chann_ = chann; pmtid_ = pmtid;
            dbta_ = dbta; rbtaA_ = rbtaA; rbtaB_ = rbtaB;
            npe_ = npe;
            cx_ = cx; cy_ = cy;
            status_ = (chann_ >= 0 && pmtid_ >= 0 && (dbta_ > 0.0 || rbtaA_ > 0.0 || rbtaB_ > 0.0) && npe_ > 0.0);
            if (!status_) clear(); 
            else {
                set_lmt_bta();
                search_closest_beta(1.0);
                wgt_ = Numc::ONE<>;
                cnt_ = Numc::ONE<>;
            }
        }
        
        const double& search_closest_beta(double bta);

        // expert only
        inline void set_cluster(short cls) { cluster_ = ((status_ && cls >= static_cast<int>(Cluster::stone) && cls <= static_cast<int>(Cluster::other)) ? cls : -1); }

        // expert only
        inline void set_wgt(double wgt) { wgt_ = ((status_ && Numc::Compare(wgt) > 0) ? wgt : Numc::ZERO<>); }
        inline void set_cnt(double cnt) { cnt_ = ((status_ && Numc::Compare(cnt) > 0) ? cnt : Numc::ZERO<>); }

        // expert only
        void set_lmt_bta(double lmtl_bta = 0.0, double lmtu_bta = 0.0);

     public :
        inline const bool&  status() const { return status_; }
        inline const short& chann()  const { return chann_; }
        inline const short& pmtid()  const { return pmtid_; }
       
        inline const short& cluster() const { return cluster_; }

        inline const short&  mode() const { return mode_; }
        inline const double& beta() const { return beta_; }

        inline const short&  type()  const { return type_; }
        inline const double& dbta()  const { return dbta_; }
        inline const double& rbtaA() const { return rbtaA_; }
        inline const double& rbtaB() const { return rbtaB_; }

        inline const bool hasDb()  const { return ((type_&1) == 1); } 
        inline const bool hasRbA() const { return ((type_&2) == 2); } 
        inline const bool hasRbB() const { return ((type_&4) == 4); } 

        inline const double& npe()   const { return npe_; }
        inline const double& cx()    const { return cx_; }
        inline const double& cy()    const { return cy_; }

        inline const double& wgt() const { return wgt_; }
        inline const double& cnt() const { return cnt_; }
    
    protected :
        inline void clear() {
            status_ = false; chann_ = -1; pmtid_ = -1;
            cluster_ = -1;
            mode_ = -1; beta_ = -1.0;
            dbta_ = -1.0; rbtaA_ = -1.0; rbtaB_ = -1.0;
            npe_ = -1.0; 
            cx_ = 0.0; cy_ = 0.0;
            wgt_ = 0.0; cnt_ = 0.0;
        }

    protected :
        bool   status_;
        short  chann_;
        short  pmtid_;

        short  cluster_; // 0 (stone), 1(cloud), 2(tumor), 3(ghost), 4(emery), 5(other)

        short  mode_; // -1(null), 0(d), 1(rA), 2(rB) 
        double beta_; // raw beta

        short  type_; // 0(null), 1(d), 2(rA), 3(drA), 4(rB), 5(drB), 6(rArB), 7(drArB)
        double dbta_;
        double rbtaA_;
        double rbtaB_;

        double npe_;
        double cx_;
        double cy_;
        
        double wgt_;
        double cnt_;
};

struct CherenkovHit_sort {
    bool operator() (const CherenkovHit& cls1, const CherenkovHit& cls2) {
        if      (cls1.cluster() < cls2.cluster()) return true;
        else if (cls1.cluster() > cls2.cluster()) return false;
        else {
            if      (cls1.chann() < cls2.chann()) return true;
            else if (cls1.chann() > cls2.chann()) return false;
        }
        return false;
    }
};


class CherenkovStone {
    public :
        CherenkovStone() { clear(); }
        ~CherenkovStone() { clear(); }
        
        CherenkovStone(const std::vector<CherenkovHit>& hits, short nhit, short npmt,
                       double cx, double cy, double npe,
                       double cnt, double nchi, double chic,
                       double dist) : CherenkovStone() { 
            hits_ = hits; 
            status_ = true;
            nhit_ = nhit; npmt_ = npmt;
            cx_ = cx; cy_ = cy; npe_ = npe;
            cnt_ = cnt; nchi_ = nchi; chic_ = chic;
            dist_ = dist;
        }
        
        inline const bool& status() const { return status_; }
        
        inline const short& nhit() const { return nhit_; }
        inline const short& npmt() const { return npmt_; }
        
        inline const double& cx()  const { return cx_; }
        inline const double& cy()  const { return cy_; }
        inline const double& npe() const { return npe_; }

        inline const double& cnt()  const { return cnt_; }
        inline const double& nchi() const { return nchi_; }
        inline const double& chic() const { return chic_; }
       
        inline const double& dist() const { return dist_; }

        inline const std::vector<CherenkovHit>& hits() const { return hits_; }
        
    protected :
        inline void clear() {
            hits_.clear(); 
            status_ = false;
            nhit_ = 0; npmt_ = 0;
            cx_ = 0; cy_ = 0; npe_ = 0;
            cnt_ = 0; nchi_ = 0; chic_ = 0;
            dist_ = 0;
        }

    protected :
        bool status_;
        
        short nhit_;
        short npmt_;
        
        double cx_;
        double cy_;
        double npe_;

        double cnt_;
        double nchi_;
        double chic_;

        double dist_;

    protected :
        std::vector<CherenkovHit> hits_;
};

struct CherenkovStone_sort {
    bool operator() (const CherenkovStone& stn1, const CherenkovStone& stn2) {
        if      (stn1.cnt() > stn2.cnt()) return true;
        else if (stn1.cnt() < stn2.cnt()) return false;
        else {
            if      (stn1.npe() >= stn2.npe()) return true;
            else if (stn1.npe() <  stn2.npe()) return false;
        }
        return false;
    }
};


class CherenkovCloud {
    public :
        CherenkovCloud() { clear(); }
        ~CherenkovCloud() { clear(); }
        
        CherenkovCloud(const std::vector<CherenkovHit>& hits,
                       short nhit, short npmt, short nhit_dir, short nhit_rfl,
                       double beta, double cbta, double npe,
                       double cnt, double nchi,
                       double misjudge) : CherenkovCloud() { 
            hits_ = hits; 
            status_ = true;
            nhit_ = nhit; npmt_ = npmt;
            nhit_dir_ = nhit_dir; nhit_rfl_ = nhit_rfl;
            beta_ = beta; cbta_ = cbta; npe_ = npe;
            cnt_ = cnt; nchi_ = nchi;
            misjudge_ = misjudge;
        }

        inline const bool& status() const { return status_; }
        
        inline const short& nhit() const { return nhit_; }
        inline const short& npmt() const { return npmt_; }
        
        inline const short& nhit_dir() const { return nhit_dir_; }
        inline const short& nhit_rfl() const { return nhit_rfl_; }
        
        inline const double& beta() const { return beta_; }
        inline const double& cbta() const { return cbta_; }
        inline const double& npe()  const { return npe_; }

        inline const double& cnt()     const { return cnt_; }
        inline const double& nchi()    const { return nchi_; }
        
        inline const double& misjudge() const { return misjudge_; }
        
        inline const std::vector<CherenkovHit>& hits() const { return hits_; }

        void set_misjudge(double mj) { misjudge_ = (mj > 0.0) ? mj : 0.0; }

    protected :
        inline void clear() {
            hits_.clear(); 
            status_ = false;
            nhit_ = 0; npmt_ = 0;
            nhit_dir_ = 0; nhit_rfl_ = 0;
            beta_ = 0; cbta_ = 0; npe_ = 0;
            cnt_ = 0; nchi_ = 0;
            misjudge_ = 0;
        }

    protected :
        bool status_;
        
        short nhit_;
        short npmt_;
       
        short nhit_dir_;
        short nhit_rfl_;

        double beta_;
        double cbta_;
        double npe_;

        double cnt_;
        double nchi_;

        double misjudge_;
        
    protected :
        std::vector<CherenkovHit> hits_;
};

struct CherenkovCloud_sort {
    bool operator() (const CherenkovCloud& cld1, const CherenkovCloud& cld2) {
        if      (cld1.cnt() > cld2.cnt()) return true;
        else if (cld1.cnt() < cld2.cnt()) return false;
        else {
            if      (cld1.npe() >= cld2.npe()) return true;
            else if (cld1.npe() <  cld2.npe()) return false;
        }
        return false;
    }
};


class CherenkovTumor {
    public :
        CherenkovTumor() { clear(); }
        ~CherenkovTumor() { clear(); }
        
        CherenkovTumor(const std::vector<CherenkovHit>& hits, short nhit, short npmt, 
                       double beta, double cbta, double npe) : CherenkovTumor() { 
            hits_ = hits; 
            status_ = true;
            nhit_ = nhit; npmt_ = npmt;
            beta_ = beta; cbta_ = cbta; npe_ = npe;
        }
        
        inline const bool& status() const { return status_; }
        
        inline const short& nhit() const { return nhit_; }
        inline const short& npmt() const { return npmt_; }
        
        inline const double& beta() const { return beta_; }
        inline const double& cbta() const { return cbta_; }
        inline const double& npe()  const { return npe_; }

        inline const std::vector<CherenkovHit>& hits() const { return hits_; }
    
    protected :
        inline void clear() {
            hits_.clear(); 
            status_ = false;
            nhit_ = 0; npmt_ = 0;
            beta_ = 0; cbta_ = 0; npe_ = 0;
        }

    protected :
        bool  status_;

        short nhit_;
        short npmt_;
        
        double beta_;
        double cbta_;
        double npe_;
    
    protected :
        std::vector<CherenkovHit> hits_;
};

struct CherenkovTumor_sort {
    bool operator() (const CherenkovTumor& tmr1, const CherenkovTumor& tmr2) {
        if      (tmr1.nhit() > tmr2.nhit()) return true;
        else if (tmr1.nhit() < tmr2.nhit()) return false;
        else {
            if      (tmr1.npmt() > tmr2.npmt()) return true;
            else if (tmr1.npmt() < tmr2.npmt()) return false;
            else {
                if      (tmr1.npe() >= tmr2.npe()) return true;
                else if (tmr1.npe() <  tmr2.npe()) return false;
            }
        }
        return false;
    }
};


class CherenkovGhost {
    public :
        CherenkovGhost() { clear(); }
        ~CherenkovGhost() { clear(); }
        
        CherenkovGhost(const std::vector<CherenkovHit>& hits, short nhit, short npmt, 
                       double beta, double cbta, double npe) : CherenkovGhost() { 
            hits_ = hits;
            status_ = true;
            nhit_ = nhit; npmt_ = npmt;
            beta_ = beta; cbta_ = cbta; npe_ = npe;
        }
        
        inline const bool& status() const { return status_; }
        
        inline const short& nhit() const { return nhit_; }
        inline const short& npmt() const { return npmt_; }
        
        inline const double& beta() const { return beta_; }
        inline const double& cbta() const { return cbta_; }
        inline const double& npe()  const { return npe_; }

        inline const std::vector<CherenkovHit>& hits() const { return hits_; }

    protected :
        inline void clear() {
            hits_.clear(); 
            status_ = false;
            nhit_ = 0; npmt_ = 0;
            beta_ = 0; cbta_ = 0; npe_ = 0;
        }

    protected :
        bool  status_;

        short nhit_;
        short npmt_;
        
        double beta_;
        double cbta_;
        double npe_;
    
    protected :
        std::vector<CherenkovHit> hits_;
};

struct CherenkovGhost_sort {
    bool operator() (const CherenkovGhost& gst1, const CherenkovGhost& gst2) {
        if      (gst1.nhit() > gst2.nhit()) return true;
        else if (gst1.nhit() < gst2.nhit()) return false;
        else {
            if      (gst1.npmt() > gst2.npmt()) return true;
            else if (gst1.npmt() < gst2.npmt()) return false;
            else {
                if      (gst1.npe() >= gst2.npe()) return true;
                else if (gst1.npe() <  gst2.npe()) return false;
            }
        }
        return false;
    }
};


class CherenkovFit {
    public :
        CherenkovFit() { clear(); }
        ~CherenkovFit() { clear(); }

        CherenkovFit(const std::vector<CherenkovHit>& args_hits, const std::array<double, 2>& pmtc, double rfr_index, double width_bta, const MultiGaus& pdf_bta, double bta_crr = 1.0);
   
        inline const bool& status() const { return succ_; }
        
        inline const std::vector<CherenkovStone>& stns() const { return stns_; }
        inline const std::vector<CherenkovCloud>& clds() const { return clds_; }
        inline const std::vector<CherenkovTumor>& tmrs() const { return tmrs_; }
        inline const std::vector<CherenkovGhost>& gsts() const { return gsts_; }
        inline const std::vector<CherenkovHit>&   hits() const { return hits_; }
        
        inline const short& nhit_total() const { return nhit_total_; }
        inline const short& nhit_stone() const { return nhit_stone_; }
        inline const short& nhit_cloud() const { return nhit_cloud_; }
        inline const short& nhit_tumor() const { return nhit_tumor_; }
        inline const short& nhit_ghost() const { return nhit_ghost_; }
        inline const short& nhit_emery() const { return nhit_emery_; }
        inline const short& nhit_other() const { return nhit_other_; }
        inline const short& nhit_other_inn() const { return nhit_other_inn_; }
        inline const short& nhit_other_out() const { return nhit_other_out_; }
        
        inline const double& npe_total() const { return npe_total_; }
        inline const double& npe_stone() const { return npe_stone_; }
        inline const double& npe_cloud() const { return npe_cloud_; }
        inline const double& npe_tumor() const { return npe_tumor_; }
        inline const double& npe_ghost() const { return npe_ghost_; }
        inline const double& npe_emery() const { return npe_emery_; }
        inline const double& npe_other() const { return npe_other_; }
        inline const double& npe_other_inn() const { return npe_other_inn_; }
        inline const double& npe_other_out() const { return npe_other_out_; }
        
        inline const Sys::HrsStopwatch& timer() const { return timer_; }

    protected :
        void clear();
        bool check();

        inline double cal_dist_to_pmtc(double cx, double cy) { return std::hypot(cx - pmtc_[0], cy - pmtc_[1]); }
        inline bool is_within_pmtc(double cx, double cy) { return (std::hypot(cx - pmtc_[0], cy - pmtc_[1]) < (WIDTH_CELL * Numc::SIX<>)); } // testcode   org = 1.0
        inline bool is_within_detectable(double cx, double cy) { return (std::hypot(cx, cy) < (MIRROR_BTM_RADIUS - WIDTH_CELL)); }

        std::vector<std::vector<CherenkovHit>> make_group_table(const std::vector<CherenkovHit>& args_hits, bool opt_stone = false, double width_stone = 1.0, bool opt_cloud = false, double width_cloud = 1.0);

        std::tuple<std::vector<CherenkovStone>, std::vector<CherenkovCloud>, std::vector<CherenkovTumor>, std::vector<CherenkovGhost>, std::vector<CherenkovHit>> build(const std::vector<CherenkovHit>& args_hits);

        std::vector<CherenkovStone> build_stone(const std::vector<CherenkovHit>& args_hits, bool weighted = false);
        std::vector<CherenkovStone> fit_stone(std::vector<CherenkovHit>& hits);
        std::vector<std::array<double, 3>> clustering_stone(std::vector<CherenkovHit>& hits);
        
        std::vector<CherenkovCloud> build_cloud(const std::vector<CherenkovHit>& args_hits, bool weighted = false);
        CherenkovCloud refit_cloud(const CherenkovCloud& cand_cld);
        std::vector<CherenkovCloud> fit_cloud(std::vector<CherenkovHit>& hits);
        std::vector<std::array<double, 2>> clustering_cloud(std::vector<CherenkovHit>& hits);
        
        std::vector<CherenkovTumor> build_tumor(const std::vector<CherenkovHit>& args_hits, bool weighted = false);
        CherenkovTumor fit_tumor(std::vector<CherenkovHit>& hits);
        
        std::vector<CherenkovGhost> build_ghost(const std::vector<CherenkovHit>& args_hits, bool weighted = false);
        CherenkovGhost fit_ghost(std::vector<CherenkovHit>& hits);
    
    protected :
        bool succ_;
        
        std::vector<CherenkovStone> stns_;
        std::vector<CherenkovCloud> clds_;
        std::vector<CherenkovTumor> tmrs_;
        std::vector<CherenkovGhost> gsts_;
        std::vector<CherenkovHit>   hits_;

        short nhit_total_;
        short nhit_stone_;
        short nhit_cloud_;
        short nhit_tumor_;
        short nhit_ghost_;
        short nhit_emery_;
        short nhit_other_;
        short nhit_other_inn_;
        short nhit_other_out_;

        double npe_total_;
        double npe_stone_;
        double npe_cloud_;
        double npe_tumor_;
        double npe_ghost_;
        double npe_emery_;
        double npe_other_;
        double npe_other_inn_;
        double npe_other_out_;

    protected :
        std::array<double, 2> pmtc_; // particle on the PMTs plane
        
        double    rfr_index_;
        double    width_bta_;
        MultiGaus scan_bta_;
        MultiGaus pdf_bta_;
        double    bta_crr_;
        double    lmtl_bta_;
        double    lmtu_bta_;

    protected :
        Sys::HrsStopwatch timer_;

    private :
        static constexpr double MIRROR_TOP_RADIUS = 60.10;
        static constexpr double MIRROR_BTM_RADIUS = 67.00;
        static constexpr double WIDTH_CELL = 0.85;
        static constexpr double WIDTH_PMT  = 3.40;
        
        static constexpr double WIDTH_CORE_COS = 1.4e-03; // (half width of PMT in the cos-core)
        static constexpr double WIDTH_CORE_COO = 0.85;    // (width of CELL in the coo-core)
        static constexpr double WIDTH_COO      = 0.66;

        static constexpr short  LMTMIN_GROUP_HIT = 3;
        static constexpr double LMTMAX_GROUP_SGM = 2.5;
        
        static constexpr short  LMTMIN_STONE_PMT_HITS_L = 3;
        static constexpr short  LMTMIN_STONE_PMT_HITS_H = 5;
        static constexpr short  LMTMIN_STONE_HITS = 6;
      
        static constexpr short  LMTMIN_CLOUD_SETS = 2;
        static constexpr short  LMTMIN_CLOUD_HITS = 3;
        static constexpr short  LMTMIN_CLOUD_PMTS = 3;
        static constexpr short  LMTMIN_CLOUD_SPCS = 3;
        
        static constexpr short  LMTMIN_TUMOR_HITS = 2;
        static constexpr short  LMTMIN_TUMOR_PMTS = 2;
        
        static constexpr short  LMTMIN_GHOST_HITS = 2;
        static constexpr short  LMTMIN_GHOST_PMTS = 2;

    private :
        static constexpr short  LMTL_ITER   = 6;
        static constexpr short  LMTM_ITER   = 12;
        static constexpr short  LMTU_ITER   = 20;
        static constexpr short  LMTMAX_ITER = 200;
       
        static constexpr double                SMOOTH_RATE  = 0.25;
        static constexpr std::array<double, 2> SMOOTH_BOUND { 3.0, 6.0 };

        static constexpr double CONVG_EPSILON    = 1.0e-06;
        static constexpr double CONVG_TOLERANCE  = 1.0e-06;
        static constexpr double CONVG_CLOSED     = 0.7;

        static constexpr double CONVG_PROB_SGM30 = 1.110900e-02; // sigma ~ 3.0
        static constexpr double CONVG_PROB_SGM40 = 3.354626e-04; // sigma ~ 4.0
        static constexpr double CONVG_PROB_SGM50 = 3.726653e-06; // sigma ~ 5.0
        static constexpr double CONVG_PROB_SGM60 = 1.522998e-08; // sigma ~ 6.0
        static constexpr double CONVG_PROB_SGM70 = 2.289735e-11; // sigma ~ 7.0
};


class CherenkovMeas {
    public :
        CherenkovMeas(const MultiGaus& mgs) : mgs_(mgs), rfr_(Numc::ONE<long double>) { sgm_ = mgs_.eftsgm(); }
        CherenkovMeas(const MultiGaus& mgs, long double rfr) : mgs_(mgs) { rfr_ = (Numc::Compare(rfr, Numc::ONE<long double>) > 0) ? rfr : Numc::ONE<long double>; sgm_ = mgs_.eftsgm(); }
        ~CherenkovMeas() {}
        
        std::array<long double, 3> minimizer(long double x, long double ibta = Numc::ONE<long double>) const;

    private :
        MultiGaus mgs_;

        long double rfr_;
        long double sgm_;
};

} // namesapce TrackSys


#endif // __TRACKLibs_CherenkovMeas_H__
