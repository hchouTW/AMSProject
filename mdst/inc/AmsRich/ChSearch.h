#ifndef __ChSearch_H__
#define __ChSearch_H__

class ChHit {
    public :
        enum Cluster { stone = 0, cloud = 1, ghost = 2, other_inn = 3, other_out = 4 };

    public :
        ChHit() { clear(); }
        ~ChHit() { clear(); }
        
        ChHit(short chann, short pmtid = -1, double dbta = -1.0, double rbtaA = -1.0, double rbtaB = -1.0, double npe = -1.0, double lx = 0.0, double ly = 0.0) : ChHit() {
            chann_ = chann; pmtid_ = pmtid;
            dbta_ = dbta; rbtaA_ = rbtaA; rbtaB_ = rbtaB;
            npe_ = npe;
            lx_ = lx; ly_ = ly;
            status_ = (chann_ >= 0 && pmtid_ >= 0 && (dbta_ > 0.0 || rbtaA_ > 0.0 || rbtaB_ > 0.0) && npe_ > 0.0);
            if (!status_) clear(); 
            else {
                set_lmt_bta();
                search_closest_beta(1.0);
            }
        }
        
        const double& search_closest_beta(double bta);

        // expert only
        inline void set_cluster(short cls) { cluster_ = ((status_ && cls >= static_cast<int>(Cluster::stone) && cls <= static_cast<int>(Cluster::other_out)) ? cls : -1); }

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

        inline bool hasDb()  const { return ((type_&1) == 1); } 
        inline bool hasRbA() const { return ((type_&2) == 2); } 
        inline bool hasRbB() const { return ((type_&4) == 4); } 

        inline const double& npe() const { return npe_; }
        inline const double& lx()  const { return lx_; }
        inline const double& ly()  const { return ly_; }
    
    protected :
        inline void clear() {
            status_ = false; chann_ = -1; pmtid_ = -1;
            cluster_ = -1;
            mode_ = -1; beta_ = -1.0;
            type_ = 0;
            dbta_ = -1.0; rbtaA_ = -1.0; rbtaB_ = -1.0;
            npe_ = -1.0; 
            lx_ = 0.0; ly_ = 0.0;
        }

    protected :
        bool   status_;
        short  chann_;
        short  pmtid_;

        short  cluster_; // 0 (stone), 1(cloud), 2(other_inn), 3(other_out)

        short  mode_; // -1(null), 0(d), 1(rA), 2(rB) 
        double beta_; // raw beta

        short  type_; // 0(null), 1(d), 2(rA), 3(drA), 4(rB), 5(drB), 6(rArB), 7(drArB)
        double dbta_;
        double rbtaA_;
        double rbtaB_;

        double npe_;
        double lx_;
        double ly_;
};

struct ChHit_sort {
    bool operator() (const ChHit& cls1, const ChHit& cls2) {
        if      (cls1.cluster() < cls2.cluster()) return true;
        else if (cls1.cluster() > cls2.cluster()) return false;
        else {
            if      (cls1.chann() < cls2.chann()) return true;
            else if (cls1.chann() > cls2.chann()) return false;
        }
        return false;
    }
};


class ChStone {
    public :
        ChStone() { clear(); }
        ~ChStone() { clear(); }
        
        ChStone(const std::vector<ChHit>& hits, short nhit, short npmt,
                double lx, double ly, double npe,
                double nchi, double dist) : ChStone() { 
            hits_ = hits; 
            status_ = true;
            nhit_ = nhit; npmt_ = npmt;
            lx_ = lx; ly_ = ly; npe_ = npe;
            nchi_ = nchi;
            dist_ = dist;
        }
        
        inline const bool& status() const { return status_; }
        
        inline const short& nhit() const { return nhit_; }
        inline const short& npmt() const { return npmt_; }
        
        inline const double& lx()  const { return lx_; }
        inline const double& ly()  const { return ly_; }
        inline const double& npe() const { return npe_; }

        inline const double& nchi() const { return nchi_; }
        inline const double& dist() const { return dist_; }

        inline const std::vector<ChHit>& hits() const { return hits_; }
        
    protected :
        inline void clear() {
            hits_.clear(); 
            status_ = false;
            nhit_ = 0; npmt_ = 0;
            lx_ = 0; ly_ = 0; npe_ = 0;
            nchi_ = 0;
            dist_ = 0;
        }

    protected :
        bool status_;
        
        short nhit_;
        short npmt_;
        
        double lx_;
        double ly_;
        double npe_;

        double nchi_;
        double dist_;

    protected :
        std::vector<ChHit> hits_;
};

struct ChStone_sort {
    bool operator() (const ChStone& stn1, const ChStone& stn2) {
        if      (stn1.nhit() > stn2.nhit()) return true;
        else if (stn1.nhit() < stn2.nhit()) return false;
        else {
            if      (stn1.npe() >= stn2.npe()) return true;
            else if (stn1.npe() <  stn2.npe()) return false;
        }
        return false;
    }
};


class ChCloud {
    public :
        ChCloud() { clear(); }
        ~ChCloud() { clear(); }
        
        ChCloud(const std::vector<ChHit>& hits,
                short nhit, short npmt, short nhit_dir, short nhit_rfl, short nhit_ght,
                double beta, double cbta, double npe,
                double nchi) : ChCloud() { 
            hits_ = hits; 
            status_ = true;
            nhit_ = nhit; npmt_ = npmt;
            nhit_dir_ = nhit_dir; nhit_rfl_ = nhit_rfl; nhit_ght_ = nhit_ght;
            beta_ = beta; cbta_ = cbta; npe_ = npe;
            nchi_ = nchi;
        }

        inline const bool& status() const { return status_; }
        
        inline const short& nhit() const { return nhit_; }
        inline const short& npmt() const { return npmt_; }
        
        inline const short& nhit_dir() const { return nhit_dir_; }
        inline const short& nhit_rfl() const { return nhit_rfl_; }
        inline const short& nhit_ght() const { return nhit_ght_; }
        
        inline const double& beta() const { return beta_; }
        inline const double& cbta() const { return cbta_; }
        inline const double& npe()  const { return npe_; }

        inline const double& nchi() const { return nchi_; }
        
        inline const std::vector<ChHit>& hits() const { return hits_; }

    protected :
        inline void clear() {
            hits_.clear(); 
            status_ = false;
            nhit_ = 0; npmt_ = 0;
            nhit_dir_ = 0; nhit_rfl_ = 0; nhit_ght_ = 0;
            beta_ = 0; cbta_ = 0; npe_ = 0;
            nchi_ = 0;
        }

    protected :
        bool status_;
        
        short nhit_;
        short npmt_;
       
        short nhit_dir_;
        short nhit_rfl_;
        short nhit_ght_;

        double beta_;
        double cbta_;
        double npe_;

        double nchi_;
        
    protected :
        std::vector<ChHit> hits_;
};

struct ChCloud_sort {
    bool operator() (const ChCloud& cld1, const ChCloud& cld2) {
        if      (cld1.nhit() > cld2.nhit()) return true;
        else if (cld1.nhit() < cld2.nhit()) return false;
        else {
            if      (cld1.npe() >= cld2.npe()) return true;
            else if (cld1.npe() <  cld2.npe()) return false;
        }
        return false;
    }
};


class ChFit {
    public :
        ChFit() { clear(); }
        ~ChFit() { clear(); }

        ChFit(const std::vector<ChHit>& args_hits, const std::array<double, 2>& pmtc, double rfr_index, double width_bta, double bta_crr = 1.0);
   
        inline const bool& status() const { return succ_; }
       
        inline const std::vector<ChStone>& stns() const { return stns_; }
        inline const std::vector<ChCloud>& clds() const { return clds_; }
        inline const std::vector<ChCloud>& ghts() const { return ghts_; }
        inline const std::vector<ChHit>&   hits() const { return hits_; }
        
        const short& nhit_total() const { return nhit_total_; }
        const short& nhit_stone() const { return nhit_stone_; } 
        const short& nhit_cloud() const { return nhit_cloud_; } 
        const short& nhit_ghost() const { return nhit_ghost_; } 
        const short& nhit_other_inn() const { return nhit_other_inn_; }
        const short& nhit_other_out() const { return nhit_other_out_; }

        const double& npe_total() const { return npe_total_; }
        const double& npe_stone() const { return npe_stone_; }
        const double& npe_cloud() const { return npe_cloud_; }
        const double& npe_ghost() const { return npe_ghost_; }
        const double& npe_other_inn() const { return npe_other_inn_; }
        const double& npe_other_out() const { return npe_other_out_; }

    protected :
        void clear();
        bool check();

        inline double cal_dist_to_pmtc(double lx, double ly) { return std::hypot(lx - pmtc_[0], ly - pmtc_[1]); }
        inline bool is_within_pmtc(double lx, double ly) { return (std::hypot(lx - pmtc_[0], ly - pmtc_[1]) < WIDTH_PMT); }
        inline bool is_within_detectable(double lx, double ly) { return (std::hypot(lx, ly) < (MIRROR_BTM_RADIUS - WIDTH_CELL)); }

        std::tuple<std::vector<ChStone>, std::vector<ChCloud>, std::vector<ChCloud>, std::vector<ChHit>> build(const std::vector<ChHit>& args_hits);
        std::vector<ChStone> build_stone(const std::vector<ChHit>& args_hits);
        std::vector<ChCloud> build_cloud(const std::vector<ChHit>& args_hits);
        std::tuple<std::vector<ChCloud>, std::vector<ChCloud>> build_cloud_and_ghost(const std::vector<ChCloud>& args_clds);

    protected :
        bool succ_;
        
        std::vector<ChStone> stns_;
        std::vector<ChCloud> clds_;
        std::vector<ChCloud> ghts_;
        std::vector<ChHit>   hits_;

        short nhit_total_;
        short nhit_stone_;
        short nhit_cloud_;
        short nhit_ghost_;
        short nhit_other_inn_;
        short nhit_other_out_;

        double npe_total_;
        double npe_stone_;
        double npe_cloud_;
        double npe_ghost_;
        double npe_other_inn_;
        double npe_other_out_;

    protected :
        std::array<double, 2> pmtc_; // particle on the PMTs plane
        
        double rfr_index_;
        double width_bta_;
        double bta_crr_;
        
        double lmtl_bta_;
        double lmtu_bta_;

    private :
        static constexpr double MIRROR_TOP_RADIUS = 60.10;
        static constexpr double MIRROR_BTM_RADIUS = 67.00;
        static constexpr double WIDTH_CELL = 0.85;
        static constexpr double WIDTH_PMT  = 3.40;
        
        static constexpr double WIDTH_CORE_COS = 1.4e-03; // (half width of PMT in the cos-core)
        
        static constexpr short LMTMIN_STONE_HITS_L = 4;
        static constexpr short LMTMIN_STONE_HITS_H = 6;
        
        static constexpr short LMTMIN_CLOUD_HITS = 3;
        static constexpr short LMTMIN_CLOUD_PMTS = 3;
        static constexpr short LMTMIN_CLOUD_ALONE_HITS = 2;
    
        static constexpr short LMTMAX_ITER = 150;
       
        static constexpr short  SMOOTH_NUM = 3;
        static constexpr double SMOOTH_BOUND = 3.0;
        static constexpr std::array<double, SMOOTH_NUM> SMOOTH_WIDTH = { 5.0, 4.0, 3.0 };

    public :
        static constexpr double AGL_BETA_WIDTH = 2.0e-03;
        static constexpr double NAF_BETA_WIDTH = 7.0e-03;
};


#endif // __ChSearch_H__
