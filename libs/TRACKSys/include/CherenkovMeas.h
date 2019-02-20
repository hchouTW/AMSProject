#ifndef __TRACKLibs_CherenkovMeas_H__
#define __TRACKLibs_CherenkovMeas_H__


namespace TrackSys {

class CherenkovHit {
    public :
        CherenkovHit(double dbta = -1.0, double rbta = -1.0, double dist = -1.0, double npe = -1.0) : status_(false), mode_(-1), beta_(-1.0), type_(0), dbta_(dbta), rbta_(rbta), dist_(dist), npe_(npe) {
            status_ = ((dbta_ > 0.0 || rbta_ > 0.0) && dist_ > 0.0 && npe_ > 0.0);
            if (!status_) { dbta_ = -1.0; rbta_ = -1.0; dist_ = -1.0; npe_ = -1.0; }
            else {
                if (dbta_ > 0.0) type_ += 1;
                if (rbta_ > 0.0) type_ += 2;
                search_closed_beta(1.0);
            }
        }
        ~CherenkovHit() {}
        
        inline const bool& status() const { return status_; }
        
        inline const short&  mode() const { return mode_; }
        inline const double& beta() const { return beta_; }

        inline const short&  type() const { return type_; }
        inline const double& dbta() const { return dbta_; }
        inline const double& rbta() const { return rbta_; }
        inline const double& dist() const { return dist_; }
        inline const double& npe()  const { return npe_; }
        
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

    protected :
        bool   status_;
        short  mode_; // -1(null), 0(d), 1(r) 
        double beta_;

        short  type_; // 0(null), 1(d), 2(r), 3(dr)
        double dbta_;
        double rbta_;
        double dist_;
        double npe_;
};


class CherenkovCls {
    public :
        CherenkovCls() { clear(); }
        ~CherenkovCls() { clear(); }
        
        CherenkovCls(const std::vector<CherenkovHit>& hits, double beta, double cnt, double nos, double eta, short ndof, double nchi, double quality, double compact_c, double compact_b, double compact_s) : CherenkovCls() { 
            hits_ = hits; 
            status_ = true; beta_ = beta;
            cnt_ = cnt; nos_ = nos; eta_ = eta;
            ndof_ = ndof; nchi_ = nchi; quality_ = quality; 
            compact_c_ = compact_c; compact_b_ = compact_b; compact_s_ = compact_s; 
        }

        inline const bool& status() const { return status_; }

        inline const double& beta() const { return beta_; }
        inline const double& cnt()  const { return cnt_; }
        inline const double& nos()  const { return nos_; }
        inline const double& eta()  const { return eta_; }
        
        inline const short&  ndof()    const { return ndof_; }
        inline const double& nchi()    const { return nchi_; }
        inline const double& quality() const { return quality_; }
        
        inline const double& compact_c() const { return compact_c_; }
        inline const double& compact_b() const { return compact_b_; }
        inline const double& compact_s() const { return compact_s_; }
        
        inline const std::vector<CherenkovHit>& hits() const { return hits_; }

    protected :
        inline void clear() {
            hits_.clear(); 
            status_ = false; beta_ = 0;
            cnt_ = 0; nos_ = 0; eta_ = 0;
            ndof_ = 0; nchi_ = 0; quality_ = 0; 
            compact_c_ = 0; compact_b_ = 0; compact_s_ = 0; 
        }

    protected :
        bool   status_;
        double beta_;
        double cnt_;
        double nos_;
        double eta_;

        short  ndof_;
        double nchi_;
        double quality_;

        double compact_c_; // compact of cluster
        double compact_b_; // compact of beta
        double compact_s_; // compact of signal(npe)

    protected :
        std::vector<CherenkovHit> hits_;
};


struct CherenkovCls_sort {
    bool operator() (const CherenkovCls& cls1, const CherenkovCls& cls2) {
        if      (cls1.cnt() > cls2.cnt()) return true;
        else if (cls1.cnt() < cls2.cnt()) return false;
        return false;
    }
};


class CherenkovFit {
    public :
        CherenkovFit(const std::vector<CherenkovHit>& args_hits, const std::array<double, 2>& scan_bta, const std::array<double, 4>& scan_npe, const std::array<double, 5>& args_bta);
        ~CherenkovFit() { clear(); }
   
        inline const bool& status() const { return succ_; }
        inline const std::vector<CherenkovHit>& hits() const { return hits_; }
        inline const std::vector<CherenkovCls>& clss() const { return clss_; }
        
        inline const Sys::HrsStopwatch& timer() const { return timer_; }

    protected :
        void clear();
        bool check();

        CherenkovCls physicalFit(const std::pair<double, std::vector<CherenkovHit>>& param);
        
        std::vector<std::pair<double, std::vector<CherenkovHit>>> clustering(std::vector<CherenkovHit>& hits);
        std::array<double, 2> clustering_evolve(std::vector<CherenkovHit>& hits, double sbta);
    
    protected :
        bool succ_;
        std::vector<CherenkovHit> hits_;
        std::vector<CherenkovCls> clss_;

    protected :
        std::array<double, 2> scan_bta_; // sgm, noise
        std::array<double, 4> scan_npe_; // kpa, mpv, sgm, noise
        std::array<double, 5> args_bta_; // wgt1, sgm1, wgt2, sgm2, noise
        
        MultiGaus mgscan_;
        MultiGaus mgfit_;

    protected :
        Sys::HrsStopwatch timer_;

    private :
        static constexpr short  LMTL_ITER = 4;
        static constexpr short  LMTU_ITER = 50;
        static constexpr double LMTL_DIST = 2.5;
        static constexpr double CONVG_EPSILON    = 5.0e-07;
        static constexpr double CONVG_TOLERANCE  = 5.0e-07;
        static constexpr double CONVG_CLUSTERING = 2.0e-05;
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
