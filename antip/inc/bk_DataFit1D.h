#ifndef __DataFit1D_H__
#define __DataFit1D_H__

#include "ceres/ceres.h" // Ceres-Solver

class ResultFit1D {
    public :
        ResultFit1D() { clear(); }

        ResultFit1D(
            const std::string& name, int ndim,
            const std::vector<std::array<double, 2>>& smp, 
            const std::vector<std::array<double, 2>>& sig, 
            const std::vector<std::array<double, 2>>& bkg,
            double num_ref = 0.0, const std::vector<std::array<double, 2>>& ref = std::vector<std::array<double, 2>>());
        
        void clear();
        void evolve(
            double num_smp, double wgt_sig, double wgt_bkg, double wgt_syst_err = 0, 
            bool build_hist = false, TAxis* axis = nullptr);

        inline int    ndof() { return ndof_; }
        inline double nchi() { return nchi_; }
       
        inline int data_ndim() { return data_ndim_; }
        inline const std::vector<std::array<double, 2>>& data_smp() const { return data_smp_; }
        inline const std::vector<std::array<double, 2>>& data_sig() const { return data_sig_; }
        inline const std::vector<std::array<double, 2>>& data_bkg() const { return data_bkg_; }
        inline const std::vector<std::array<double, 2>>& data_ref() const { return data_ref_; }

        inline const double& wgt_sig() const { return wgt_sig_; }
        inline const double& wgt_bkg() const { return wgt_bkg_; }
        inline const double& wgt_err() const { return wgt_err_; }
        inline const double& wgt_stat_err() const { return wgt_stat_err_; }
        inline const double& wgt_syst_err() const { return wgt_syst_err_; }
        
        inline const double& num_ref() const { return num_ref_; }
        inline const double& num_smp() const { return num_smp_; }
        inline const double& num_sig() const { return num_sig_; }
        inline const double& num_bkg() const { return num_bkg_; }
        inline const double& num_err() const { return num_err_; }
        inline const double& num_stat_err() const { return num_stat_err_; }
        inline const double& num_syst_err() const { return num_syst_err_; }

        inline std::shared_ptr<TH1D>& href() { return hist_ref_; }
        inline std::shared_ptr<TH1D>& hsmp() { return hist_smp_; }
        inline std::shared_ptr<TH1D>& hsum() { return hist_sum_; }
        inline std::shared_ptr<TH1D>& hsig() { return hist_sig_; }
        inline std::shared_ptr<TH1D>& hbkg() { return hist_bkg_; }

    protected :
        std::string name_;

        int data_ndim_;
        std::vector<std::array<double, 2>> data_smp_;
        std::vector<std::array<double, 2>> data_sig_;
        std::vector<std::array<double, 2>> data_bkg_;
        std::vector<std::array<double, 2>> data_ref_;

        int    ndof_;
        double nchi_;
        
        double wgt_sig_;
        double wgt_bkg_;
        double wgt_err_;
        double wgt_stat_err_;
        double wgt_syst_err_;

        double num_ref_;
        double num_smp_;
        double num_sig_;
        double num_bkg_;
        double num_err_;
        double num_stat_err_;
        double num_syst_err_;

        std::vector<std::array<double, 2>> ref_;
        std::vector<std::array<double, 2>> smp_;
        std::vector<std::array<double, 2>> sum_;
        std::vector<std::array<double, 2>> sig_;
        std::vector<std::array<double, 2>> bkg_;
        
        std::shared_ptr<TH1D> hist_ref_;
        std::shared_ptr<TH1D> hist_smp_;
        std::shared_ptr<TH1D> hist_sum_;
        std::shared_ptr<TH1D> hist_sig_;
        std::shared_ptr<TH1D> hist_bkg_;
};


class DataFit1D {
    public :
        DataFit1D(const std::string& name, TH1D* hsmp, TH1D* hsig, TH1D* hbkg, TH1D* href = nullptr);
       
        inline const int& data_ndim() const { return data_ndim_; }

        inline const double& data_num_smp() const { return data_num_smp_; }
        inline const double& data_num_sig() const { return data_num_sig_; }
        inline const double& data_num_bkg() const { return data_num_bkg_; }
        inline const double& data_num_ref() const { return data_num_ref_; }

        inline ResultFit1D& result()      { return result_; }
        inline ResultFit1D& result_syst() { return result_syst_; }

    protected :
        void clear();

        std::tuple<double, std::vector<std::array<double, 2>>> build_data(TH1D* hist);
        ResultFit1D fit(const std::string& name, const std::vector<std::array<double, 2>>& smp, const std::vector<std::array<double, 2>>& sig, const std::vector<std::array<double, 2>>& bkg, bool build_hist = false);
        bool evolve();
    
    private :
        std::string name_;

        TH1D* hsmp_;
        TH1D* hsig_;
        TH1D* hbkg_;
        TH1D* href_;

        int data_ndim_;
        
        double data_num_smp_;
        double data_num_sig_;
        double data_num_bkg_;
        double data_num_ref_;

        std::vector<std::array<double, 2>> data_smp_;
        std::vector<std::array<double, 2>> data_sig_;
        std::vector<std::array<double, 2>> data_bkg_;
        std::vector<std::array<double, 2>> data_ref_;

        ResultFit1D result_;
        ResultFit1D result_syst_;
};


class VirtualDataFit1D : public ceres::CostFunction {
    public :
        VirtualDataFit1D(int ndim, const std::vector<std::array<double, 2>>& smp, const std::vector<std::array<double, 2>>& sig, const std::vector<std::array<double, 2>>& bkg) {
            ndim_ = ndim;
            smp_ = smp;
            sig_ = sig;
            bkg_ = bkg;
            
            nres_ = ndim_;
            npar_ = 1;
            ndof_ = -npar_;
        
            has_ = std::vector<bool>(ndim_, false);
            has_smp_ = std::vector<bool>(ndim_, false);
            has_sig_ = std::vector<bool>(ndim_, false);
            has_bkg_ = std::vector<bool>(ndim_, false);
            for (int it = 0; it < ndim; ++it) {
                if (smp_[it][0] <= 0.0 && sig_[it][0] <= 0.0 && bkg_[it][0] <= 0.0) continue;
                if (smp_[it][0] > 0.0) has_smp_[it] = true;
                if (sig_[it][0] > 0.0) has_sig_[it] = true;
                if (bkg_[it][0] > 0.0) has_bkg_[it] = true;
                has_[it] = true;
                ndof_++;
            }

            set_num_residuals(nres_);
            mutable_parameter_block_sizes()->clear(); 
            mutable_parameter_block_sizes()->push_back(npar_);
        }
    
    public :
        virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;

    protected :
        int ndim_;
        std::vector<bool> has_;
        std::vector<bool> has_smp_;
        std::vector<bool> has_sig_;
        std::vector<bool> has_bkg_;
        std::vector<std::array<double, 2>> smp_;
        std::vector<std::array<double, 2>> sig_;
        std::vector<std::array<double, 2>> bkg_;

        int nres_;
        int npar_;
        int ndof_;
};

#endif // __DataFit1D_H__
