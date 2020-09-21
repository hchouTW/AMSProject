#ifndef __DataFit1D_H__
#define __DataFit1D_H__

#include "ceres/ceres.h" // Ceres-Solver

class ResultFit1D {
    public :
        ResultFit1D() { clear(); }

        ResultFit1D(
            const std::string& name, int ndim, int ntmp,
            const std::vector<std::vector<std::array<double, 2>>>& tmps, 
            double num_smp, const std::vector<std::array<double, 2>>& smp, 
            double num_ref = 0.0, const std::vector<std::array<double, 2>>& ref = std::vector<std::array<double, 2>>());
        
        void clear();

        void set_fitting_result(int ndof, double nchi, const std::vector<double>& wgt_tmps, const std::vector<double>& wgt_errs);
        void set_histogram_result(bool hist_build = false, TAxis* hist_axis = nullptr, const std::string& hist_xtitle = "", const std::string& hist_ytitle = "");

        inline int    ndof() { return ndof_; }
        inline double nchi() { return nchi_; }
       
        inline int data_ndim() { return data_ndim_; }
        inline int data_ntmp() { return data_ntmp_; }

        inline const std::vector<std::vector<std::array<double, 2>>>& data_tmps() const { return data_tmps_; }
        inline const std::vector<std::array<double, 2>>& data_smp() const { return data_smp_; }
        inline const std::vector<std::array<double, 2>>& data_ref() const { return data_ref_; }

        inline const std::vector<double>& wgt_tmps() const { return wgt_tmps_; }
        inline const std::vector<double>& wgt_errs() const { return wgt_errs_; }
    
        inline const std::vector<double>& num_tmps() const { return num_tmps_; }
        inline const std::vector<double>& num_errs() const { return num_errs_; }
        inline const double& num_sum() const { return num_sum_; }
        inline const double& num_smp() const { return num_smp_; }
        inline const double& num_ref() const { return num_ref_; }

        inline std::vector<std::shared_ptr<TH1D>>& htmps() { return hist_tmps_; }
        inline std::shared_ptr<TH1D>& hsum() { return hist_sum_; }
        inline std::shared_ptr<TH1D>& hsmp() { return hist_smp_; }
        inline std::shared_ptr<TH1D>& href() { return hist_ref_; }

    protected :
        bool  hist_build_;
        TAxis hist_axis_;

        std::string data_name_;
        int         data_ndim_;
        int         data_ntmp_;
        
        std::vector<std::vector<std::array<double, 2>>> data_tmps_;
        std::vector<std::array<double, 2>> data_smp_;
        std::vector<std::array<double, 2>> data_ref_;

        std::vector<double> wgt_tmps_;
        std::vector<double> wgt_errs_;

        std::vector<double> num_tmps_;
        std::vector<double> num_errs_;
        double num_sum_;
        double num_smp_;
        double num_ref_;

        std::vector<std::vector<std::array<double, 2>>> tmps_;
        std::vector<std::array<double, 2>> sum_;
        std::vector<std::array<double, 2>> smp_;
        std::vector<std::array<double, 2>> ref_;

        std::vector<std::shared_ptr<TH1D>> hist_tmps_;
        std::shared_ptr<TH1D> hist_sum_;
        std::shared_ptr<TH1D> hist_smp_;
        std::shared_ptr<TH1D> hist_ref_;
        
        int    ndof_;
        double nchi_;
};


class DataFit1D {
    public :
        DataFit1D(const std::string& name, std::vector<TH1D*> htmps, TH1D* hsmp, TH1D* href = nullptr, bool hist_build = false, const std::string& hist_xtitle = "", const std::string& hist_ytitle = "");
       
        inline const int& data_ndim() const { return data_ndim_; }
        inline const int& data_ntmp() const { return data_ntmp_; }

        inline const double& data_num_smp()      const { return data_num_smp_; }
        inline const double& data_num_ref()      const { return data_num_ref_; }
        inline const double& data_num_tmp(int i) const { return data_num_tmps_.at(i); }

        inline const ResultFit1D& result() const { return result_; }

    protected :
        void clear();

        std::tuple<double, std::vector<std::array<double, 2>>> build_data(TH1D* hist);
        bool fit();
    
    private :
        bool hist_build_;
        TAxis hist_axis_;
        std::string hist_xtitle_;
        std::string hist_ytitle_;

        std::string data_name_;
        int         data_ndim_;
        int         data_ntmp_;

        std::vector<TH1D*> htmps_;
        TH1D*              hsmp_;
        TH1D*              href_;

        std::vector<double> data_num_tmps_;
        double              data_num_smp_;
        double              data_num_ref_;

        std::vector<std::vector<std::array<double, 2>>> data_tmps_;
        std::vector<std::array<double, 2>>              data_smp_;
        std::vector<std::array<double, 2>>              data_ref_;

        ResultFit1D result_;
};


class VirtualDataFit1D : public ceres::CostFunction {
    public :
        VirtualDataFit1D(int ndim, int ntmp, const std::vector<std::vector<std::array<double, 2>>>& tmps, const std::vector<std::array<double, 2>>& smp) {
            ndim_ = ndim;
            ntmp_ = ntmp;
            tmps_ = tmps;
            smp_  = smp;

            nres_ = ndim_;
            npar_ = tmps_.size();
            ndof_ = -npar_;
        
            has_          = std::vector<bool>(ndim_, false);
            has_tmps_     = std::vector<std::vector<bool>>(ntmp_, std::vector<bool>(ndim_, false));
            has_tmps_all_ = std::vector<bool>(ndim_, false);
            has_smp_      = std::vector<bool>(ndim_, false);
            
            for (int ib = 0; ib < ndim_; ++ib) {
                for (int it = 0; it < ntmp_; ++it) {
                    has_tmps_.at(it).at(ib) = (tmps_.at(it).at(ib)[0] > 0.0);
                    if (has_tmps_.at(it).at(ib)) has_tmps_all_.at(ib) = true;
                }
                has_smp_.at(ib) = (smp_.at(ib)[0] > 0.0);
                has_.at(ib) = (has_tmps_all_.at(ib) || has_smp_.at(ib));
                if (!has_.at(ib)) continue;
                ndof_++;
            }

            set_num_residuals(nres_);
            mutable_parameter_block_sizes()->clear(); 
            mutable_parameter_block_sizes()->push_back(npar_);
        }
    
    public :
        virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;

        inline const int& nres() const { return nres_; }
        inline const int& npar() const { return npar_; }
        inline const int& ndof() const { return ndof_; }

    protected :
        int ndim_;
        int ntmp_;

        std::vector<bool>              has_;
        std::vector<std::vector<bool>> has_tmps_;
        std::vector<bool>              has_tmps_all_;
        std::vector<bool>              has_smp_;
        
        std::vector<std::vector<std::array<double, 2>>> tmps_;
        std::vector<std::array<double, 2>>              smp_;

        int nres_;
        int npar_;
        int ndof_;
};

#endif // __DataFit1D_H__
