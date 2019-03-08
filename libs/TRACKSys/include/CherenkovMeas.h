#ifndef __TRACKLibs_CherenkovMeas_H__
#define __TRACKLibs_CherenkovMeas_H__


namespace TrackSys {

class CherenkovHit {
    public :
        CherenkovHit(int chann = -1, double dbta = -1.0, double rbta = -1.0, double npe = -1.0, double cx = 0.0, double cy = 0.0) : status_(false), chann_(chann), pmtid_(chann/16), mode_(-1), beta_(-1.0), type_(0), dbta_(dbta), rbta_(rbta), npe_(npe), cx_(cx), cy_(cy), cnt_(0.0), wgt_(0.0) {
            status_ = (chann_ >= 0 && pmtid_ >= 0 && (dbta_ > 0.0 || rbta_ > 0.0) && npe_ > 0.0);
            if (!status_) clear(); 
            else {
                if (dbta_ > 0.0) type_ += 1;
                if (rbta_ > 0.0) type_ += 2;
                search_closed_beta(1.0);
                cnt_ = Numc::ONE<>;
                wgt_ = Numc::ONE<>;
            }
        }
        ~CherenkovHit() {}
        
        inline const double& search_closed_beta(double bta) {
            mode_ = -1; beta_ = -1.0;
            if      (type_ == 1) { mode_ = 0; beta_ = dbta_; }
            else if (type_ == 2) { mode_ = 1; beta_ = rbta_; }
            else if (type_ == 3) {
                mode_ = (std::fabs(dbta_ - bta) <= std::fabs(rbta_ - bta)) ? 0 : 1;
                beta_ = (mode_ == 0 ? dbta_ : rbta_);
            }
            return beta_;
        }

        // expert only
        inline void set_cnt(double cnt) { if (Numc::Compare(cnt) > 0) cnt_ = cnt; else cnt_ = Numc::ZERO<>; }
        inline void set_wgt(double wgt) { if (Numc::Compare(wgt) > 0) wgt_ = wgt; else wgt_ = Numc::ZERO<>; }

     public :
        inline const bool& status() const { return status_; }
        inline const int&  chann()  const { return chann_; }
        inline const int&  pmtid()  const { return pmtid_; }
        
        inline const short&  mode() const { return mode_; }
        inline const double& beta() const { return beta_; }

        inline const short&  type() const { return type_; }
        inline const double& dbta() const { return dbta_; }
        inline const double& rbta() const { return rbta_; }
        inline const double& npe()  const { return npe_; }
        inline const double& cx() const { return cx_; }
        inline const double& cy() const { return cy_; }

        inline const double& cnt() const { return cnt_; }
        inline const double& wgt() const { return wgt_; }
    
    protected :
        inline void clear() {
            status_ = false; chann_ = -1; pmtid_ = -1;
            mode_ = -1; beta_ = -1.0;
            dbta_ = -1.0; rbta_ = -1.0; 
            npe_ = -1.0; 
            cx_ = 0.0; cy_ = 0.0;
            cnt_ = 0.0; wgt_ = 0.0; 
        }

    protected :
        bool   status_;
        int    chann_;
        short  pmtid_;
        short  mode_; // -1(null), 0(d), 1(r) 
        double beta_; // raw beta

        short  type_; // 0(null), 1(d), 2(r), 3(dr)
        double dbta_;
        double rbta_;
        double npe_;
        double cx_;
        double cy_;
        
        double cnt_;
        double wgt_;
};


struct CherenkovHit_sort {
    bool operator() (const CherenkovHit& cls1, const CherenkovHit& cls2) {
        if      (cls1.chann() < cls2.chann()) return true;
        else if (cls1.chann() > cls2.chann()) return false;
        return false;
    }
};


class CherenkovStone {
    public :
        CherenkovStone() { clear(); }
        ~CherenkovStone() { clear(); }
        
        CherenkovStone(const std::vector<CherenkovHit>& hits, double cx, double cy, double dist, double npe,
                       double cnt, double nchi, double quality) : CherenkovStone() { 
            hits_ = hits; 
            status_ = true; cx_ = cx; cy_ = cy; dist_ = dist; npe_ = npe;
            cnt_ = cnt; nchi_ = nchi; quality_ = quality;
        }
        
        inline const bool& status() const { return status_; }
        inline const double& cx()   const { return cx_; }
        inline const double& cy()   const { return cy_; }
        inline const double& dist() const { return dist_; }
        inline const double& npe()  const { return npe_; }

        inline const double& cnt()     const { return cnt_; }
        inline const double& nchi()    const { return nchi_; }
        inline const double& quality() const { return quality_; }
        
        inline const std::vector<CherenkovHit>& hits() const { return hits_; }
        
    protected :
        inline void clear() {
            hits_.clear(); 
            status_ = false; cx_ = 0; cy_ = 0; dist_ = 0; npe_ = 0;
            cnt_ = 0; nchi_ = 0; quality_ = 0;
        }

    protected :
        bool   status_;
        double cx_;
        double cy_;
        double dist_;
        double npe_;

        double cnt_;
        double nchi_;
        double quality_;
    
    protected :
        std::vector<CherenkovHit> hits_;
};


struct CherenkovStone_sort {
    bool operator() (const CherenkovStone& stn1, const CherenkovStone& stn2) {
        if      (stn1.cnt() > stn2.cnt()) return true;
        else if (stn1.cnt() < stn2.cnt()) return false;
        return false;
    }
};



class CherenkovCloud {
    public :
        CherenkovCloud() { clear(); }
        ~CherenkovCloud() { clear(); }
        
        CherenkovCloud(const std::vector<CherenkovHit>& hits, double beta, double cbta, double npe, short nhit, short npmt,
                       double cnt, double wgt, double all, double nchi, double quality) : CherenkovCloud() { 
            hits_ = hits; 
            status_ = true; beta_ = beta; cbta_ = cbta; npe_ = npe;
            nhit_ = nhit; npmt_ = npmt;
            cnt_ = cnt; wgt_ = wgt; all_ = all; nchi_ = nchi; quality_ = quality;
        }

        inline const bool& status() const { return status_; }
        inline const double& beta() const { return beta_; }
        inline const double& cbta() const { return cbta_; }
        inline const double& npe()  const { return npe_; }
        
        inline const short& nhit() const { return nhit_; }
        inline const short& npmt() const { return npmt_; }

        inline const double& cnt()     const { return cnt_; }
        inline const double& wgt()     const { return wgt_; }
        inline const double& all()     const { return all_; }
        inline const double& nchi()    const { return nchi_; }
        inline const double& quality() const { return quality_; }
        
        inline const std::vector<CherenkovHit>& hits() const { return hits_; }

    protected :
        inline void clear() {
            hits_.clear(); 
            status_ = false; beta_ = 0; cbta_ = 0; npe_ = 0;
            nhit_ = 0; npmt_ = 0;
            cnt_ = 0; wgt_ = 0; all_ = 0; nchi_ = 0; quality_ = 0; 
        }

    protected :
        bool   status_;
        double beta_;
        double cbta_;
        double npe_;

        short nhit_;
        short npmt_;

        double cnt_;
        double wgt_;
        double all_;
        double nchi_;
        double quality_;

    protected :
        std::vector<CherenkovHit> hits_;
};


struct CherenkovCloud_sort {
    bool operator() (const CherenkovCloud& cls1, const CherenkovCloud& cls2) {
        if      (cls1.all() > cls2.all()) return true;
        else if (cls1.all() < cls2.all()) return false;
        return false;
    }
};


class CherenkovFit {
    public :
        CherenkovFit(const std::vector<CherenkovHit>& args_hits, const std::array<double, 2>& args_core, double width_bta, const std::array<double, 4>& args_bta, const std::array<double, 4>& args_npe, double bta_crr = 1.0);
        ~CherenkovFit() { clear(); }
   
        inline const bool& status() const { return succ_; }
        inline const std::vector<CherenkovStone>& stns() const { return stns_; }
        inline const std::vector<CherenkovCloud>& clds() const { return clds_; }
    
        inline const Sys::HrsStopwatch& timer() const { return timer_; }

    protected :
        void clear();
        bool check();
       
        std::vector<CherenkovStone> fit_stone(std::vector<CherenkovHit>& hits);
        std::vector<std::array<double, 3>> clustering_stone(std::vector<CherenkovHit>& hits);
        std::vector<std::array<double, 3>> clustering_evolve_stone(std::vector<CherenkovHit>& hits);
        
        CherenkovCloud refit_cloud(const CherenkovCloud& cand_cld);

        std::vector<CherenkovCloud> fit_cloud(std::vector<CherenkovHit>& hits);
        std::vector<std::array<double, 2>> clustering_cloud(std::vector<CherenkovHit>& hits);
        std::vector<std::array<double, 2>> clustering_evolve_cloud(std::vector<CherenkovHit>& hits);
    
    protected :
        bool succ_;
        std::vector<CherenkovStone> stns_;
        std::vector<CherenkovCloud> clds_;

    protected :
        std::array<double, 2> args_core_; // cx, cy
        std::array<double, 4> args_bta_; // wgt1, sgm1, wgt2, sgm2, wgt3, sgm3
        std::array<double, 3> args_npe_sig_; // kpa, mpv, sgm
        double                args_npe_nos_;
        double                width_bta_;
        MultiGaus             pdf_bta_;
        MultiGaus             pdf_bta_acc_;
        double                bta_crr_;

    protected :
        Sys::HrsStopwatch timer_;

    private :
        static constexpr double WIDTH_CELL = 0.85;
        static constexpr double WIDTH_PMT  = 3.40;
        static constexpr double WIDTH_BETA = 2.0e-03;

        static constexpr short  LMTMIN_STONE_PMT_HITS_L = 3;
        static constexpr short  LMTMIN_STONE_PMT_HITS_H = 5;
        static constexpr short  LMTMIN_STONE_HITS = 6;
        
        static constexpr short  LMTMIN_CLOUD_PMTS = 3;
        static constexpr short  LMTMIN_CLOUD_HITS = 3;

    private :
        static constexpr short  LMTL_ITER   = 6;
        static constexpr short  LMTM_ITER   = 12;
        static constexpr short  LMTU_ITER   = 20;
        static constexpr short  LMTMAX_ITER = 200;
        
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
