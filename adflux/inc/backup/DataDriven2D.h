#ifndef __DataDriven2D_H__
#define __DataDriven2D_H__

#include "ceres/ceres.h" // Ceres-Solver

class DataDriven2D {
    public :
        DataDriven2D(TH2D* hsmp, TH2D* hsig);

        inline TH1D* hsmp_x() { return hfit_smp_x_; }
        inline TH1D* hsmp_y() { return hfit_smp_y_; }
        
        inline TH1D* hsig_x() { return hfit_sig_x_; }
        inline TH1D* hsig_y() { return hfit_sig_y_; }
        
        inline TH1D* hbkg_x() { return hfit_bkg_x_; }
        inline TH1D* hbkg_y() { return hfit_bkg_y_; }

        inline THStack* hx() { return hfit_x_; }
        inline THStack* hy() { return hfit_y_; }

    protected :
        void clear();

        std::tuple<std::vector<std::vector<std::array<double, 2>>>, std::vector<std::vector<std::array<double, 2>>>> build_data(TH2D* hist);
        bool fit();
        bool evolve();

    private :
        TH2D* hsmp_;
        TH2D* hsig_;

        std::vector<std::vector<std::array<double, 2>>> data_smp_;
        std::vector<std::vector<std::array<double, 2>>> data_sig_;
        std::vector<std::vector<std::array<double, 2>>> data_bkg_;

        std::vector<std::vector<std::array<double, 2>>> data_int_smp_;
        std::vector<std::vector<std::array<double, 2>>> data_int_sig_;
        std::vector<std::vector<std::array<double, 2>>> data_int_bkg_;

        int param_nx_;
        int param_ny_;
        double param_wgt_sig_;
        double param_wgt_bkg_;
        std::vector<double> param_bkg_eftx_;
        std::vector<double> param_bkg_efty_;

        std::vector<std::array<double, 2>> fit_smp_x_;
        std::vector<std::array<double, 2>> fit_smp_y_;
        std::vector<std::array<double, 2>> fit_sig_x_;
        std::vector<std::array<double, 2>> fit_sig_y_;
        std::vector<std::array<double, 2>> fit_bkg_x_;
        std::vector<std::array<double, 2>> fit_bkg_y_;

        TH1D* hfit_smp_x_;
        TH1D* hfit_smp_y_;
        TH1D* hfit_sig_x_;
        TH1D* hfit_sig_y_;
        TH1D* hfit_bkg_x_;
        TH1D* hfit_bkg_y_;

        THStack* hfit_x_;
        THStack* hfit_y_;
};

class VirtualDataDriven2D : public ceres::CostFunction {
    public :
        VirtualDataDriven2D(int nx, int ny, const std::vector<std::vector<std::array<double, 2>>>& smp, const std::vector<std::vector<std::array<double, 2>>>& sig) {
            nx_ = nx;
            ny_ = ny;
            smp_ = smp;
            sig_ = sig;
            
            nres_ = nx_ * ny_ + (nx_ - 2) * (ny_ - 2);
            npar_ = nx_ + ny_ - 1;
            ndof_ = nres_ - npar_;
            
            set_num_residuals(nres_);
            mutable_parameter_block_sizes()->clear(); 
            mutable_parameter_block_sizes()->push_back(1);
            mutable_parameter_block_sizes()->push_back(nx_-1);
            mutable_parameter_block_sizes()->push_back(ny_-1);
        }
    
    public :
        virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;

    protected :
        int nx_;
        int ny_;
        std::vector<std::vector<std::array<double, 2>>> smp_;
        std::vector<std::vector<std::array<double, 2>>> sig_;

        int nres_;
        int npar_;
        int ndof_;
};

#endif // __DataDriven2D_H__
