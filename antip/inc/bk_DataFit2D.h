#ifndef __DataFit2D_H__
#define __DataFit2D_H__

#include "ceres/ceres.h" // Ceres-Solver

class ResultFit2D {
    public :
        ResultFit2D() { clear(); }

        ResultFit2D(
            const std::string& name, int ndimx, int ndimy, int ntmp,
            const std::vector<std::vector<std::array<double, 2>>>& tmps, 
            double num_smp, const std::vector<std::array<double, 2>>& smp, 
            double num_ref = 0.0, const std::vector<std::array<double, 2>>& ref = std::vector<std::array<double, 2>>());
        
        void clear();

        void set_fitting_result(int ndof, double nchi, const std::vector<double>& wgt_tmps, const std::vector<double>& wgt_errs);
        void set_histogram_result(bool hist_build = false, TAxis* hist_xaxis = nullptr, TAxis* hist_yaxis = nullptr, const std::string& hist_xtitle = "", const std::string& hist_ytitle = "", const std::string& hist_ztitle = "");

        inline int    ndof() { return ndof_; }
        inline double nchi() { return nchi_; }
       
        inline int data_ndimx() { return data_ndimx_; }
        inline int data_ndimy() { return data_ndimy_; }
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

        inline std::vector<std::shared_ptr<TH2D>>& htmps() { return hist_tmps_; }
        inline std::shared_ptr<TH2D>& hsum() { return hist_sum_; }
        inline std::shared_ptr<TH2D>& hsmp() { return hist_smp_; }
        inline std::shared_ptr<TH2D>& href() { return hist_ref_; }
        
        inline std::vector<std::shared_ptr<TH1D>>& htmps_x() { return hist_tmps_x_; }
        inline std::shared_ptr<TH1D>& hsum_x() { return hist_sum_x_; }
        inline std::shared_ptr<TH1D>& hsmp_x() { return hist_smp_x_; }
        inline std::shared_ptr<TH1D>& href_x() { return hist_ref_x_; }
        
        inline std::vector<std::shared_ptr<TH1D>>& htmps_y() { return hist_tmps_y_; }
        inline std::shared_ptr<TH1D>& hsum_y() { return hist_sum_y_; }
        inline std::shared_ptr<TH1D>& hsmp_y() { return hist_smp_y_; }
        inline std::shared_ptr<TH1D>& href_y() { return hist_ref_y_; }

    protected :
        bool  hist_build_;
        TAxis hist_xaxis_;
        TAxis hist_yaxis_;

        std::string data_name_;
        int         data_ndimx_;
        int         data_ndimy_;
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
        
        std::vector<std::vector<std::array<double, 2>>> tmps_x_;
        std::vector<std::array<double, 2>> sum_x_;
        std::vector<std::array<double, 2>> smp_x_;
        std::vector<std::array<double, 2>> ref_x_;
        
        std::vector<std::vector<std::array<double, 2>>> tmps_y_;
        std::vector<std::array<double, 2>> sum_y_;
        std::vector<std::array<double, 2>> smp_y_;
        std::vector<std::array<double, 2>> ref_y_;

        std::vector<std::shared_ptr<TH2D>> hist_tmps_;
        std::shared_ptr<TH2D> hist_sum_;
        std::shared_ptr<TH2D> hist_smp_;
        std::shared_ptr<TH2D> hist_ref_;
        
        std::vector<std::shared_ptr<TH1D>> hist_tmps_x_;
        std::shared_ptr<TH1D> hist_sum_x_;
        std::shared_ptr<TH1D> hist_smp_x_;
        std::shared_ptr<TH1D> hist_ref_x_;
        
        std::vector<std::shared_ptr<TH1D>> hist_tmps_y_;
        std::shared_ptr<TH1D> hist_sum_y_;
        std::shared_ptr<TH1D> hist_smp_y_;
        std::shared_ptr<TH1D> hist_ref_y_;
        
        int    ndof_;
        double nchi_;
};


class DataFit2D {
    public :
        DataFit2D(const std::string& name, std::vector<TH2D*> htmps, TH2D* hsmp, TH2D* href = nullptr, bool hist_build = false, const std::string& hist_xtitle = "", const std::string& hist_ytitle = "", const std::string& hist_ztitle = "");
       
        inline const int& data_ndimx() const { return data_ndimx_; }
        inline const int& data_ndimy() const { return data_ndimy_; }
        inline const int& data_ntmp() const { return data_ntmp_; }

        inline const double& data_num_smp()      const { return data_num_smp_; }
        inline const double& data_num_ref()      const { return data_num_ref_; }
        inline const double& data_num_tmp(int i) const { return data_num_tmps_.at(i); }

        inline const ResultFit2D& result() const { return result_; }

    protected :
        void clear();

        std::tuple<double, std::vector<std::array<double, 2>>> build_data(TH2D* hist);
        bool fit();
    
    private :
        bool hist_build_;
        TAxis hist_xaxis_;
        TAxis hist_yaxis_;
        std::string hist_xtitle_;
        std::string hist_ytitle_;
        std::string hist_ztitle_;

        std::string data_name_;
        int         data_ndimx_;
        int         data_ndimy_;
        int         data_ntmp_;

        std::vector<TH2D*> htmps_;
        TH2D*              hsmp_;
        TH2D*              href_;

        std::vector<double> data_num_tmps_;
        double              data_num_smp_;
        double              data_num_ref_;

        std::vector<std::vector<std::array<double, 2>>> data_tmps_;
        std::vector<std::array<double, 2>>              data_smp_;
        std::vector<std::array<double, 2>>              data_ref_;

        ResultFit2D result_;
};


class VirtualDataFit2D : public ceres::CostFunction {
    public :
        VirtualDataFit2D(int ndimx, int ndimy, int ntmp, const std::vector<std::vector<std::array<double, 2>>>& tmps, const std::vector<std::array<double, 2>>& smp) {
            ndimx_ = ndimx;
            ndimy_ = ndimy;
            ntmp_ = ntmp;
            tmps_ = tmps;
            smp_  = smp;

            nres_ = ndimy_*ndimx_;
            npar_ = tmps_.size();
            ndof_ = -npar_;
        
            has_      = std::vector<bool>(ndimy_*ndimx_, false);
            has_tmps_ = std::vector<std::vector<bool>>(ntmp_, std::vector<bool>(ndimy_*ndimx_, false));
            has_smp_  = std::vector<bool>(ndimy_*ndimx_, false);
            
            for (int iby = 0; iby < ndimy_; ++iby) {
            for (int ibx = 0; ibx < ndimx_; ++ibx) {
                int idx = iby * ndimx_ + ibx;
                bool has = false;
                for (int it = 0; it < ntmp_; ++it) {
                    has_tmps_.at(it).at(idx) = (tmps_.at(it).at(idx)[0] > 0.0);
                    if (has_tmps_.at(it).at(idx)) has = true;
                }
                has_smp_.at(idx) = (smp_.at(idx)[0] > 0.0);
                if (has_smp_.at(idx)) has = true;
                if (!has) continue;
                has_.at(idx) = true;
                ndof_++;
            }}

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
        int ndimx_;
        int ndimy_;
        int ntmp_;

        std::vector<bool>              has_;
        std::vector<std::vector<bool>> has_tmps_;
        std::vector<bool>              has_smp_;
        
        std::vector<std::vector<std::array<double, 2>>> tmps_;
        std::vector<std::array<double, 2>>              smp_;

        int nres_;
        int npar_;
        int ndof_;
};

#endif // __DataFit2D_H__
