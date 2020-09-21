#ifndef __HistFit1D_H__
#define __HistFit1D_H__

#include "ceres/ceres.h" // Ceres-Solver
#include "TH1D.h"

namespace HistFit {

class Point1D {
    public :
        Point1D() : val(0.0), err(0.0) {}
        Point1D(double value, double error = 0.0) : val(value), err(error) {}

    public :
        double val;
        double err;
};

class Data1D {
    public :
        Data1D(const std::string name = "", const std::string title = "") : name_(name), title_(title), nbins_(0), sum_(0.0) {}
        Data1D(const std::string name, const std::string title, const std::vector<Point1D> data);

        const std::string& name()  const { return name_; }
        const std::string& title() const { return title_; }

        const int&    nbins() const { return nbins_; }
        const double& sum()   const { return sum_; }

        const std::vector<Point1D>& data() const { return data_; }
        const std::vector<Point1D>& pdf()  const { return pdf_; }

        const Point1D& data(int i) const { return data_.at(i); }
        const Point1D& pdf(int i)  const { return pdf_.at(i); }
        
        std::vector<Point1D> weight_pdf(double w) const {
            std::vector<Point1D> wpdf;
            for (int i = 0; i < pdf_.size(); ++i) wpdf.push_back( Point1D(w * pdf_.at(i).val, w * pdf_.at(i).err) );
            return wpdf;
        }

        Point1D weight_pdf(double w, int i) const { return Point1D(w * pdf_.at(i).val, w * pdf_.at(i).err); }

        void push_back(double val, double err = 0.0) { data_.push_back(Point1D(val, err)); }
        void push_back(const Point1D& pnt) { data_.push_back(pnt); }

        bool build();
        std::vector<Data1D> build_fluc(int nsample, int method) const;
        Data1D build_fluc(int method = 0) const { return build_fluc(1, method).at(0); }

    protected :
        std::string name_;
        std::string title_;
        std::vector<Point1D> data_;

    protected :
        int    nbins_;
        double sum_;
        std::vector<Point1D> pdf_;
};

class DataFit1D {
    public :
        DataFit1D(const Data1D& smp, const std::vector<Data1D>& tmps);

        const Data1D&              smp()  const { return smp_; }
        const std::vector<Data1D>& tmps() const { return tmps_; }
        
        const Data1D& tmps(int i) const { return tmps_.at(i); }

        const double&              nsmp() const { return nsmp_; }
        const std::vector<double>& wgts() const { return wgts_; }
        const std::vector<double>& errs() const { return errs_; }
        
        const double& wgts(int i) const { return wgts_.at(i); }
        const double& errs(int i) const { return errs_.at(i); }
        
        const Data1D&              chi_fits() const { return chi_fits_; }
        const Data1D&              sum_tmps() const { return sum_tmps_; }
        const std::vector<Data1D>& wgt_tmps() const { return wgt_tmps_; }
        
        const Data1D& wgt_tmps(int i) const { return wgt_tmps_.at(i); }

        const int& nbins() const { return nbins_; }
        const int& ntmps() const { return ntmps_; }

        const bool&   status() const { return status_; }
        const double& nchi()   const { return nchi_; }
        const int&    ndof()   const { return ndof_; }

    protected :
        void clear();
        bool build(const Data1D& smp, const std::vector<Data1D>& tmps);
        bool fit();
        bool fine();

    protected :
        Data1D              smp_;
        std::vector<Data1D> tmps_;

        double              nsmp_;
        std::vector<double> wgts_;
        std::vector<double> errs_;
    
        Data1D              chi_fits_;
        Data1D              sum_tmps_;
        std::vector<Data1D> wgt_tmps_;

        int nbins_;
        int ntmps_;

        bool   status_;
        double nchi_;
        int    ndof_;

        std::vector<bool>   bins_;
        Data1D              reduce_smp_;
        std::vector<Data1D> reduce_tmps_;
        std::vector<int>    reduce_pars_;
};

class VirtualDataFit1D : public ceres::CostFunction  {
    public :
        VirtualDataFit1D(const Data1D& smp, const std::vector<Data1D>& tmps) {
            smp_  = smp;
            tmps_ = tmps;
            
            nbins_ = smp_.nbins();
            ntmps_ = tmps_.size();

            set_num_residuals(nbins_);
            mutable_parameter_block_sizes()->clear(); 
            mutable_parameter_block_sizes()->push_back(ntmps_);
        }

    public :
        virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;

    protected :
        Data1D              smp_;
        std::vector<Data1D> tmps_;
        
        int nbins_;
        int ntmps_;
};

class Axis1D {
    public :
        Axis1D() : status_(false), ref_sbin_(-1), ref_ebin_(-1), ref_nbins_(0), nbins_(0) {}
        
        Axis1D(const std::string& name_x, const std::string& name_y, const std::vector<double>& bins) : Axis1D() {
            if (bins.size() == 0) return;
            status_ = true;
            name_x_ = name_x;
            name_y_ = name_y;
            ref_nbins_ = bins.size() - 1;
            ref_bins_  = bins;
            ref_sbin_ = 1;
            ref_ebin_ = ref_nbins_;
            nbins_    = ref_nbins_;
            bins_     = ref_bins_;
        }
        Axis1D(const std::string& name_x, const std::string& name_y, const TAxis* axis, int sbin = -1, int ebin = -1) : Axis1D() {
            if (axis == nullptr || axis->GetNbins() == 0) return;
            int ref_nbins = axis->GetNbins();
            std::vector<double> ref_bins;
            for (int ib = 1; ib <= axis->GetNbins()+1; ++ib) {
                ref_bins.push_back(axis->GetBinLowEdge(ib));
            }
            int ref_sbin = (sbin <= 0 || sbin > axis->GetNbins()) ?                1 : sbin;
            int ref_ebin = (ebin <= 0 || ebin > axis->GetNbins()) ? axis->GetNbins() : ebin;
            if (ref_sbin > ref_ebin) return;

            int nbins = ref_ebin - ref_sbin + 1;
            std::vector<double> bins;
            for (int it = ref_sbin-1; it <= ref_ebin; ++it) {
                bins.push_back(ref_bins.at(it));
            }

            status_ = true;
            name_x_ = name_x;
            name_y_ = name_y;
            ref_sbin_ = ref_sbin;
            ref_ebin_ = ref_ebin;
            ref_nbins_ = ref_nbins;
            ref_bins_  = ref_bins;
            nbins_ = nbins;
            bins_  = bins;
        }
        
        const bool& status() const { return status_; }

        const std::string& name_x() const { return name_x_; }
        const std::string& name_y() const { return name_y_; }

        const int&                 ref_sbin() const { return ref_sbin_; }
        const int&                 ref_ebin() const { return ref_ebin_; }
        const int&                 ref_nbins() const { return ref_nbins_; }
        const std::vector<double>& ref_bins()  const { return ref_bins_; }
        
        const int&                 nbins() const { return nbins_; }
        const std::vector<double>& bins()  const { return bins_; }

    protected :
        bool status_;

        std::string name_x_;
        std::string name_y_;

        int                 ref_sbin_;
        int                 ref_ebin_;
        int                 ref_nbins_;
        std::vector<double> ref_bins_;

        int                 nbins_;
        std::vector<double> bins_;
};

class Hist1D {
    public :
        Hist1D(const std::string& name = "hist", const std::string& title = "") : status_(false), name_(name), title_(title) {}
        
        Hist1D(const std::string& name, const std::string& title, const Data1D& data, const Axis1D& axis, bool sw_limit = false) : Hist1D(name, title) {
            if (data.nbins() != axis.nbins()) return;
            if (data.nbins() == 0) return;
            if (axis.nbins() == 0) return;
            status_ = true;
            axis_   = axis;
            data_   = Data1D(name_, title_);
            for (int ib = 0; ib < axis.nbins(); ++ib) {
                Point1D pnt = data.data(ib);
                if (sw_limit && pnt.val < 1.0e-3) { pnt.val = 0.0; pnt.err = 1.0; }
                data_.push_back(pnt);
            }
            data_.build();
        }
        
        Hist1D(const std::string& name, const std::string& title, const TH1D* hist, const Axis1D& axis, bool sw_limit = false) : Hist1D(name, title) {
            if (hist == nullptr || !axis.status() || hist->GetXaxis()->GetNbins() != axis.ref_nbins()) return;
            if (hist->GetXaxis()->GetNbins() == 0) return;
            if (axis.ref_nbins() == 0) return;
            status_ = true;
            axis_   = axis;
            data_   = Data1D(name_, title_);
            for (int ib = axis.ref_sbin(); ib <= axis.ref_ebin(); ++ib) {
                Point1D pnt(hist->GetBinContent(ib), hist->GetBinError(ib));
                if (sw_limit && pnt.val < 1.0e-5) { pnt.val = 0.0; pnt.err = 1.0e-5; }
                data_.push_back(pnt);
            }
            data_.build();
        }

        const bool& status() const { return status_; }

        const std::string& name()  const { return name_; }
        const std::string& title() const { return title_; }

        const Axis1D& axis() const { return axis_; }
        const Data1D& data() const { return data_; }

        const std::shared_ptr<TH1D>& operator()() const { return hist_; }
        const TH1D* get() const { return hist_.get(); }

        void build_hist() {
            if (!status_) return;
            hist_ = std::make_shared<TH1D>(name_.c_str(), title_.c_str(), axis_.nbins(), &(axis_.bins().at(0)));
            hist_.get()->GetXaxis()->SetTitle(axis_.name_x().c_str());
            hist_.get()->GetYaxis()->SetTitle(axis_.name_y().c_str());
            for (int ib = 1; ib <= axis_.nbins(); ++ib) {
                Point1D pnt = data_.data(ib-1);
                hist_.get()->SetBinContent(ib, pnt.val);
                hist_.get()->SetBinError  (ib, ((pnt.val == 0.0) ? 0.0 : pnt.err));
            }
        }

        void delete_hist() {
            if (!status_) return;
            if (!hist_ || hist_.get() == nullptr) return;
            hist_.get()->Delete();
            hist_.reset();
        }

    protected :
        bool status_;

        std::string name_;
        std::string title_;

        Axis1D                axis_;
        Data1D                data_;
        std::shared_ptr<TH1D> hist_;
};

class HistFit1D {
    public :
        HistFit1D(const TH1D* smp, const std::vector<TH1D*>& tmps, const Axis1D& axis, const std::string& prefix_name = "FIT1D_", bool with_fluc = false, bool with_build_hist = false);
        HistFit1D(const Hist1D& smp, const std::vector<Hist1D>& tmps, const Axis1D& axis, const std::string& prefix_name = "FIT1D_", bool with_fluc = false, bool with_build_hist = false);

        const bool& with_fluc() const { return with_fluc_; }

        const Hist1D&              ref_smp()  const { return hist_ref_smp_; }
        const Hist1D&              sum_tmps() const { return hist_sum_tmps_; }
        const std::vector<Hist1D>& wgt_tmps() const { return hist_wgt_tmps_; }

        const Hist1D& wgt_tmps(int i) const { return hist_wgt_tmps_.at(i); }

        const double&               nsmp() const { return nsmp_; }
        const std::vector<double>&  wgts() const { return wgts_; }
        const std::vector<double>&  errs() const { return errs_; }
        const std::vector<Point1D>& fluc() const { return fluc_; }
        
        const double&  wgts(int i) const { return wgts_.at(i); }
        const double&  errs(int i) const { return errs_.at(i); }
        const Point1D& fluc(int i) const { return fluc_.at(i); }
        
        const std::vector<std::vector<Point1D>>& fluc_datas() const { return fluc_datas_; }
        const std::vector<Point1D>& fluc_datas(int i) const { return fluc_datas_.at(i); }
        
        const std::vector<Hist1D>& fluc_hists() const { return fluc_hists_; }
        const Hist1D& fluc_hists(int i) const { return fluc_hists_.at(i); }
        
        const std::vector<double>& fluc_nchi() const { return fluc_nchi_; }
        const double& fluc_nchi(int i) const { return fluc_nchi_.at(i); }
       
        const bool&   status() const { return status_; }
        const int&    ntmps()  const { return ntmps_; }
        const double& nchi()   const { return nchi_; }
        const int&    ndof()   const { return ndof_; }

        bool build_hist() {
            if (!status_) return false;
            hist_ref_smp_.build_hist();
            hist_sum_tmps_.build_hist();
            for (auto&& hist : hist_wgt_tmps_) hist.build_hist();
            for (auto&& hist : fluc_hists_) hist.build_hist();
            with_build_hist_ = true;
            return true;
        }

        bool delete_hist() {
            if (!status_) return false;
            hist_ref_smp_.delete_hist();
            hist_sum_tmps_.delete_hist();
            for (auto&& hist : hist_wgt_tmps_) hist.delete_hist();
            for (auto&& hist : fluc_hists_) hist.delete_hist();
            with_build_hist_ = false;
            return true;
        }

    protected :
        bool init(const TH1D* smp, const std::vector<TH1D*>& tmps, const Axis1D& axis, const std::string& prefix_name, bool with_fluc, bool with_build_hist);
        bool init(const Hist1D& smp, const std::vector<Hist1D>& tmps, const Axis1D& axis, const std::string& prefix_name, bool with_fluc, bool with_build_hist);
        
        void clear();
        bool fit();
        bool flucfit();

    protected :
        bool with_fluc_;
        bool with_build_hist_;
        std::string prefix_name_;
        Axis1D axis_;

        Hist1D              hist_smp_;
        std::vector<Hist1D> hist_tmps_;

        Hist1D              hist_ref_smp_;
        Hist1D              hist_sum_tmps_;
        std::vector<Hist1D> hist_wgt_tmps_;

        double               nsmp_;
        std::vector<double>  wgts_;
        std::vector<double>  errs_;
        std::vector<Point1D> fluc_;


        std::vector<std::vector<Point1D>> fluc_datas_;
        std::vector<Hist1D> fluc_hists_;
        std::vector<double> fluc_nchi_;

        bool   status_;
        int    ntmps_;
        double nchi_;
        int    ndof_;
};


class CeresSolverGaussianResidual {
    public :
        CeresSolverGaussianResidual(const std::vector<std::array<Point1D, 2>>& datas) : datas_(datas) {}
        
        template <typename T>
        bool operator()(T const* const* parameters, T* residuals) const {
            if (datas_.size() == 0) return false;
            std::fill_n(residuals, datas_.size(), T(0));
            const T& parN = parameters[0][0];
            const T& parM = parameters[0][1];
            const T& parS = parameters[0][2];
            const T NEGHALF = T(-0.5);
            const T INV2PI  = T(0.3989422804);
            for (int it = 0; it < datas_.size(); ++it) {
                const Point1D& pntx = datas_.at(it).at(0);
                const Point1D& pnty = datas_.at(it).at(1);
                T x = (T(pntx.val) - parM) / parS;
                T y = INV2PI * (parN / parS) * exp(NEGHALF * x * x);
                residuals[it] = (T(pnty.val) - y) / T(pnty.err);
            }
            return true;
        }
        
    private :
        const std::vector<std::array<Point1D, 2>> datas_;
};

class VirtualFlucFit1D {
    public :
        // Factory method.
        static ceres::CostFunction* Create(const std::vector<std::array<Point1D, 2>>& datas, const std::vector<double>& params) {
            // Stride. Number of derivatives to calculate. See below for details:
            // http://ceres-solver.org/nnls_modeling.html#dynamicautodiffcostfunction
            const int kStride = 1;
            
            const int kNumResiduals  = datas.size();
            const int kNumParameters = 3;
            if (datas.size()  <= kNumParameters) return nullptr;
            if (params.size() != kNumParameters) return nullptr;
            
            ceres::DynamicAutoDiffCostFunction<CeresSolverGaussianResidual, kStride>* cost_function =
                new ceres::DynamicAutoDiffCostFunction<CeresSolverGaussianResidual, kStride>(
                        new CeresSolverGaussianResidual(datas));
            cost_function->AddParameterBlock(kNumParameters);
            cost_function->SetNumResiduals(kNumResiduals);

            return cost_function;
        }
};

} // namespace HistFit

#endif // __HistFit1D_H__
