#ifndef __TRACKLibs_CherenkovMeas_H__
#define __TRACKLibs_CherenkovMeas_H__


namespace TrackSys {

class CherenkovHit {
    public :
        CherenkovHit(short idx = -1, double dbta = -1.0, double rbta = -1.0, double dist = -1.0, double npe = -1.0) : status_(false), index_(idx), mode_(-1), beta_(-1.0), type_(0), dbta_(dbta), rbta_(rbta), dist_(dist), npe_(npe) {
            if (dbta_ && dist_ < LMT_DIST_FOR_D) dbta_ = 0.0;
            if (rbta_ && dist_ < LMT_DIST_FOR_R) rbta_ = 0.0;
            status_ = (index_ >= 0 && (dbta_ > 0.0 || rbta_ > 0.0) && npe_ > 0.0);
            if (!status_) { index_ = -1; dbta_ = -1.0; rbta_ = -1.0; dist_ = -1.0; npe_ = -1.0; }
            else {
                if (dbta_ > 0.0) type_ += 1;
                if (rbta_ > 0.0) type_ += 2;
                search_closed_beta(1.0);
            }
        }
        ~CherenkovHit() {}
        
        inline const bool&  status() const { return status_; }
        inline const short& index()  const { return index_; }
        
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
        short  index_;
        short  mode_; // -1(null), 0(d), 1(r) 
        double beta_;

        short  type_; // 0(null), 1(d), 2(r), 3(dr)
        double dbta_;
        double rbta_;
        double dist_;
        double npe_;
        
    private :
        static constexpr double LMT_DIST_FOR_D = 5.0;
        static constexpr double LMT_DIST_FOR_R = 3.0;
};


class CherenkovCls {
    public :
        CherenkovCls() { clear(); }
        ~CherenkovCls() { clear(); }
        
        CherenkovCls(const std::vector<CherenkovHit>& hits, double beta, 
                     short ndof, double nchi, double quality, 
                     const std::array<double, 4>& sig, const std::array<double, 4>& nos, const std::array<double, 4>& compact,
                     double signal_s2s4, double signal_s3s4, double noise_s2s4, double noise_s3s4, double sn_s2s4, double sn_s3s4) : CherenkovCls() { 
            hits_ = hits; 
            status_ = true; beta_ = beta;
            ndof_ = ndof; nchi_ = nchi; quality_ = quality; 
            sig_ = sig; nos_ = nos; compact_ = compact;
            signal_s2s4_ = signal_s2s4; noise_s2s4_ = noise_s2s4;
            signal_s3s4_ = signal_s3s4; noise_s3s4_ = noise_s3s4;
            sn_s2s4_ = sn_s2s4; sn_s3s4_ = sn_s3s4;
        }

        inline const bool& status() const { return status_; }
        inline const double& beta() const { return beta_; }

        inline const short&  ndof()    const { return ndof_; }
        inline const double& nchi()    const { return nchi_; }
        inline const double& quality() const { return quality_; }
        
        inline const double& sig(int it)     const { return sig_.at(it); }
        inline const double& nos(int it)     const { return nos_.at(it); }
        inline const double& compact(int it) const { return compact_.at(it); }
        
        inline const std::array<double, 4>& sig()     const { return sig_; }
        inline const std::array<double, 4>& nos()     const { return nos_; }
        inline const std::array<double, 4>& compact() const { return compact_; }
        
        inline const double& signal_s2s4() const { return signal_s2s4_; }
        inline const double& signal_s3s4() const { return signal_s3s4_; }
        
        inline const double& noise_s2s4() const { return noise_s2s4_; }
        inline const double& noise_s3s4() const { return noise_s3s4_; }
        
        inline const double& sn_s2s4() const { return sn_s2s4_; }
        inline const double& sn_s3s4() const { return sn_s3s4_; }
        
        inline const std::vector<CherenkovHit>& hits() const { return hits_; }

    protected :
        inline void clear() {
            hits_.clear(); 
            status_ = false; beta_ = 0;
            ndof_ = 0; nchi_ = 0; quality_ = 0; 
            sig_.fill(0); nos_.fill(0); compact_.fill(0);
            signal_s2s4_ = 0; noise_s2s4_ = 0;
            signal_s3s4_ = 0; noise_s3s4_ = 0;
            sn_s2s4_ = 0; sn_s3s4_ = 0;
        }

    protected :
        bool   status_;
        double beta_;
        
        short  ndof_;
        double nchi_;
        double quality_;
        
        std::array<double, 4> sig_;
        std::array<double, 4> nos_;
        std::array<double, 4> compact_; // compact of cluster

        double signal_s2s4_;
        double signal_s3s4_;
        
        double noise_s2s4_;
        double noise_s3s4_;

        double sn_s2s4_;
        double sn_s3s4_;

    protected :
        std::vector<CherenkovHit> hits_;
};


struct CherenkovCls_sort {
    bool operator() (const CherenkovCls& cls1, const CherenkovCls& cls2) {
        if      (cls1.sig(0) > cls2.sig(0)) return true;
        else if (cls1.sig(0) < cls2.sig(0)) return false;
        return false;
    }
};


class CherenkovFit {
    public :
        CherenkovFit(const std::vector<CherenkovHit>& args_hits, const std::array<double, 4>& args_bta, const std::array<double, 2>& scan_bta, const std::array<double, 4>& scan_npe);
        ~CherenkovFit() { clear(); }
   
        inline const bool& status() const { return succ_; }
        inline const std::vector<CherenkovHit>& hits() const { return hits_; }
        inline const std::vector<CherenkovCls>& clss() const { return clss_; }
        
        inline const Sys::HrsStopwatch& timer() const { return timer_; }

    protected :
        void clear();
        bool check();

        bool fit(const std::vector<CherenkovHit>& hits);

        CherenkovCls physicalFit(const std::pair<double, std::vector<CherenkovHit>>& param);
        
        std::vector<std::pair<double, std::vector<CherenkovHit>>> clustering(std::vector<CherenkovHit>& hits);
        std::array<double, 2> clustering_evolve(std::vector<CherenkovHit>& hits, double sbta);
    
    protected :
        bool succ_;
        std::vector<CherenkovHit> hits_;
        std::vector<CherenkovCls> clss_;

    protected :
        std::array<double, 4> args_bta_; // wgt1, sgm1, wgt2, sgm2
        MultiGaus pdf_bta_;
        
        double scan_bta_sig_;
        double scan_bta_nos_;
        MultiGaus pdf_scan_bta_;

        std::array<double, 3> scan_npe_sig_; // kpa, mpv, sgm
        double                scan_npe_nos_;

        double convg_epsilon_;
        double convg_tolerance_;
        double convg_closed_;

    protected :
        Sys::HrsStopwatch timer_;

    private :
        static constexpr short  LMTL_ITER = 3;
        static constexpr short  LMTU_ITER = 50;
        static constexpr double CONVG_EPSILON   = 1.0e-05;
        static constexpr double CONVG_TOLERANCE = 1.0e-05;
        static constexpr double CONVG_CLOSED    = 0.3;

        static constexpr double PROB_SIGMA2 = 0.135335280000;
        static constexpr double PROB_SIGMA3 = 0.011108997000;
        static constexpr double PROB_SIGMA4 = 0.000335462630;
        static constexpr double PROB_SIGMA5 = 0.000003726653;
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
