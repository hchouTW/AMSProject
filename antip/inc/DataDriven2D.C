#ifndef __DataDriven2D_C__
#define __DataDriven2D_C__

#include "DataDriven2D.h"

DataDriven2D::DataDriven2D(TH2D* hsmp, TH2D* hsig) {
    clear();
    if (hsmp == nullptr || hsig == nullptr) return;
    
    hsmp_ = hsmp;
    auto&& hsmp_data = build_data(hsmp);
    data_smp_     = std::get<0>(hsmp_data);
    data_int_smp_ = std::get<1>(hsmp_data);
    
    hsig_ = hsig;
    auto&& hsig_data = build_data(hsig);
    data_sig_     = std::get<0>(hsig_data);
    data_int_sig_ = std::get<1>(hsig_data);

    data_bkg_     = std::vector<std::vector<std::array<double, 2>>>(data_sig_    .size(), std::vector<std::array<double, 2>>(data_sig_    .at(0).size()));
    data_int_bkg_ = std::vector<std::vector<std::array<double, 2>>>(data_int_sig_.size(), std::vector<std::array<double, 2>>(data_int_sig_.at(0).size()));

    param_bkg_eftx_ = std::vector<double>(data_bkg_.size());
    param_bkg_efty_ = std::vector<double>(data_bkg_.at(0).size());

    param_wgt_sig_ = 0.0;
    param_wgt_bkg_ = 0.0;
    param_nx_ = data_bkg_.size();
    param_ny_ = data_bkg_.at(0).size();

    fit();
}

void DataDriven2D::clear() {
    hsmp_ = nullptr;
    hsig_ = nullptr;
    param_nx_ = 0;
    param_ny_ = 0;
    param_wgt_sig_ = 0;
    param_wgt_bkg_ = 0;
    hfit_smp_x_ = nullptr;
    hfit_smp_y_ = nullptr;
    hfit_sig_x_ = nullptr;
    hfit_sig_y_ = nullptr;
    hfit_bkg_x_ = nullptr;
    hfit_bkg_y_ = nullptr;
    hfit_x_ = nullptr;
    hfit_y_ = nullptr;
}

std::tuple<std::vector<std::vector<std::array<double, 2>>>, std::vector<std::vector<std::array<double, 2>>>> DataDriven2D::build_data(TH2D* hist) {
    std::tuple<std::vector<std::vector<std::array<double, 2>>>, std::vector<std::vector<std::array<double, 2>>>> datas;
    if (hist == nullptr) return datas;

    std::vector<std::vector<std::array<double, 2>>> data    (hist->GetXaxis()->GetNbins(), std::vector<std::array<double, 2>>(hist->GetYaxis()->GetNbins()));
    std::vector<std::vector<std::array<double, 2>>> data_int(hist->GetXaxis()->GetNbins(), std::vector<std::array<double, 2>>(hist->GetYaxis()->GetNbins()));

    for (int ib = 1; ib <= hist->GetXaxis()->GetNbins(); ++ib) {
    for (int jb = 1; jb <= hist->GetYaxis()->GetNbins(); ++jb) {
        double val = hist->GetBinContent(ib, jb);
        double err = hist->GetBinError  (ib, jb);
        data[ib-1][jb-1][0] = val;
        data[ib-1][jb-1][1] = err;
    }}
    
    for (int it = 0; it < data.size(); ++it) {
    for (int jt = 0; jt < data.at(0).size(); ++jt) {
        double val = 0.0;
        double err = 0.0;
        for (int ib = 0; ib <= it; ++ib) {
        for (int jb = 0; jb <= jt; ++jb) {
            val += data[ib][jb][0];
            err += data[ib][jb][1] * data[ib][jb][1];
        }}
        err = std::sqrt(err);
        data_int[it][jt][0] = val;
        data_int[it][jt][1] = err;
    }}

    double total = (data_int.back().back())[0];
    double error = std::sqrt(1.0 / total);
    for (auto&& yvec : data) {
    for (auto&& arr : yvec) {
        arr[0] /= total;
        arr[1] /= total;
        if (arr[1] == 0.0) arr[1] = error;
    }}
    for (auto&& yvec : data_int) {
    for (auto&& arr : yvec) {
        arr[0] /= total;
        arr[1] /= total;
        if (arr[1] == 0.0) arr[1] = error;
    }}

    std::get<0>(datas) = data;
    std::get<1>(datas) = data_int;
    return datas;
}

bool DataDriven2D::fit() {
    param_wgt_bkg_ = 0.5;

    //for (int it = 1; it < param_bkg_eftx_.size(); ++it) {
    //    param_bkg_eftx_[it] = 1.0 / param_bkg_eftx_.size();
    //}
    //for (int it = 1; it < param_bkg_efty_.size(); ++it) {
    //    param_bkg_efty_[it] = 1.0 / param_bkg_efty_.size();
    //}

    for (int it = 0; it < param_bkg_eftx_.size(); ++it) {
        param_bkg_eftx_[it] = data_int_smp_.at(it).back()[0];
    } 

    for (int it = 0; it < param_bkg_efty_.size(); ++it) {
        param_bkg_efty_[it] = data_int_smp_.back().at(it)[0];
    } 
    
    // CeresSolver: Cost Function
    ceres::CostFunction* cost_function = new VirtualDataDriven2D(param_nx_, param_ny_, data_int_smp_, data_int_sig_);

    // CeresSolver: Problem
    ceres::Problem problem;
    problem.AddResidualBlock(cost_function, nullptr, &param_wgt_bkg_, param_bkg_eftx_.data(), param_bkg_efty_.data());

    problem.SetParameterLowerBound(&param_wgt_bkg_, 0, 0.0);
    problem.SetParameterUpperBound(&param_wgt_bkg_, 0, 1.0);

    for (int it = 0; it < param_bkg_eftx_.size()-1; ++it) {
        problem.SetParameterLowerBound(param_bkg_eftx_.data(), it, 0.0);
        problem.SetParameterUpperBound(param_bkg_eftx_.data(), it, 1.0);
    }
    
    for (int it = 0; it < param_bkg_efty_.size()-1; ++it) {
        problem.SetParameterLowerBound(param_bkg_efty_.data(), it, 0.0);
        problem.SetParameterUpperBound(param_bkg_efty_.data(), it, 1.0);
    }

    // CeresSolver: Options
    ceres::Solver::Options options;
    options.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT;
    options.max_num_iterations = 100;

    // CeresSolver: Summary
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    if (!summary.IsSolutionUsable()) return false;
    //if (ceres::NO_CONVERGENCE == summary.termination_type) return false;
    std::cerr << summary.FullReport() << std::endl;
    
    param_wgt_sig_ = 1.0 - param_wgt_bkg_;

    for (int it = 0; it < param_bkg_eftx_.size(); ++it) {
    for (int jt = 0; jt < param_bkg_efty_.size(); ++jt) {
        double eft = param_bkg_eftx_[it] * param_bkg_efty_[jt];
        data_int_bkg_[it][jt][0] = eft;
        
        double elem_eftx = (it == 0) ? param_bkg_eftx_[it] : (param_bkg_eftx_[it] - param_bkg_eftx_[it-1]);
        double elem_efty = (jt == 0) ? param_bkg_efty_[jt] : (param_bkg_efty_[jt] - param_bkg_efty_[jt-1]);
        double elem_eft  = elem_eftx * elem_efty;
        data_bkg_[it][jt][0] = elem_eft;
    }}

    return evolve();
}
        
bool DataDriven2D::evolve() {
    fit_smp_x_ = std::vector<std::array<double, 2>>(param_nx_, { 0.0, 0.0 });
    fit_smp_y_ = std::vector<std::array<double, 2>>(param_ny_, { 0.0, 0.0 });
    for (int ii = 0; ii < param_nx_; ++ii) {
        fit_smp_x_.at(ii)[0] = (ii == 0) ? data_int_smp_.at(0).back()[0] : (data_int_smp_.at(ii).back()[0] - data_int_smp_.at(ii-1).back()[0]);
    }
    for (int ii = 0; ii < param_ny_; ++ii) {
        fit_smp_y_.at(ii)[0] = (ii == 0) ? data_int_smp_.back().at(0)[0] : (data_int_smp_.back().at(ii)[0] - data_int_smp_.back().at(ii-1)[0]);
    }
    
    fit_sig_x_ = std::vector<std::array<double, 2>>(param_nx_, { 0.0, 0.0 });
    fit_sig_y_ = std::vector<std::array<double, 2>>(param_ny_, { 0.0, 0.0 });
    for (int ii = 0; ii < param_nx_; ++ii) {
        fit_sig_x_.at(ii)[0] = param_wgt_sig_ * ((ii == 0) ? data_int_sig_.at(0).back()[0] : (data_int_sig_.at(ii).back()[0] - data_int_sig_.at(ii-1).back()[0]));
    }
    for (int ii = 0; ii < param_ny_; ++ii) {
        fit_sig_y_.at(ii)[0] = param_wgt_sig_ * ((ii == 0) ? data_int_sig_.back().at(0)[0] : (data_int_sig_.back().at(ii)[0] - data_int_sig_.back().at(ii-1)[0]));
    }
    
    fit_bkg_x_ = std::vector<std::array<double, 2>>(param_nx_, { 0.0, 0.0 });
    fit_bkg_y_ = std::vector<std::array<double, 2>>(param_ny_, { 0.0, 0.0 });
    for (int ii = 0; ii < param_nx_; ++ii) {
        fit_bkg_x_.at(ii)[0] = param_wgt_bkg_ * ((ii == 0) ? data_int_bkg_.at(0).back()[0] : (data_int_bkg_.at(ii).back()[0] - data_int_bkg_.at(ii-1).back()[0]));
    }
    for (int ii = 0; ii < param_ny_; ++ii) {
        fit_bkg_y_.at(ii)[0] = param_wgt_bkg_ * ((ii == 0) ? data_int_bkg_.back().at(0)[0] : (data_int_bkg_.back().at(ii)[0] - data_int_bkg_.back().at(ii-1)[0]));
    }

    TAxis xaxis = *(hsmp_->GetXaxis());
    TAxis yaxis = *(hsmp_->GetYaxis());

    hfit_smp_x_ = new TH1D("hfit_smp_x", "", xaxis.GetNbins(), xaxis.GetXbins()->GetArray());
    for (int ii = 0; ii < param_nx_; ++ii) {
        hfit_smp_x_->SetBinContent(ii+1, fit_smp_x_.at(ii)[0]);
        hfit_smp_x_->SetBinError  (ii+1, fit_smp_x_.at(ii)[1]);
    }

    hfit_smp_y_ = new TH1D("hfit_smp_y", "", yaxis.GetNbins(), yaxis.GetXbins()->GetArray());
    for (int ii = 0; ii < param_ny_; ++ii) {
        hfit_smp_y_->SetBinContent(ii+1, fit_smp_y_.at(ii)[0]);
        hfit_smp_y_->SetBinError  (ii+1, fit_smp_y_.at(ii)[1]);
    }
    
    hfit_sig_x_ = new TH1D("hfit_sig_x", "", xaxis.GetNbins(), xaxis.GetXbins()->GetArray());
    for (int ii = 0; ii < param_nx_; ++ii) {
        hfit_sig_x_->SetBinContent(ii+1, fit_sig_x_.at(ii)[0]);
        hfit_sig_x_->SetBinError  (ii+1, fit_sig_x_.at(ii)[1]);
    }

    hfit_sig_y_ = new TH1D("hfit_sig_y", "", yaxis.GetNbins(), yaxis.GetXbins()->GetArray());
    for (int ii = 0; ii < param_ny_; ++ii) {
        hfit_sig_y_->SetBinContent(ii+1, fit_sig_y_.at(ii)[0]);
        hfit_sig_y_->SetBinError  (ii+1, fit_sig_y_.at(ii)[1]);
    }

    hfit_bkg_x_ = new TH1D("hfit_bkg_x", "", xaxis.GetNbins(), xaxis.GetXbins()->GetArray());
    for (int ii = 0; ii < param_nx_; ++ii) {
        hfit_bkg_x_->SetBinContent(ii+1, fit_bkg_x_.at(ii)[0]);
        hfit_bkg_x_->SetBinError  (ii+1, fit_bkg_x_.at(ii)[1]);
    }

    hfit_bkg_y_ = new TH1D("hfit_bkg_y", "", yaxis.GetNbins(), yaxis.GetXbins()->GetArray());
    for (int ii = 0; ii < param_ny_; ++ii) {
        hfit_bkg_y_->SetBinContent(ii+1, fit_bkg_y_.at(ii)[0]);
        hfit_bkg_y_->SetBinError  (ii+1, fit_bkg_y_.at(ii)[1]);
    }

    hfit_x_ = new THStack();
    hfit_x_->SetTitle("hfit_x");
    hfit_x_->Add(hfit_smp_x_);
    hfit_x_->Add(hfit_sig_x_);
    hfit_x_->Add(hfit_bkg_x_);

    hfit_y_ = new THStack(); 
    hfit_y_->SetTitle("hfit_y");
    hfit_y_->Add(hfit_smp_y_);
    hfit_y_->Add(hfit_sig_y_);
    hfit_y_->Add(hfit_bkg_y_);

    return true;
}

bool VirtualDataDriven2D::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const {
    if (nx_ <= 1 || ny_ <= 1) return false;
    std::fill_n(residuals, nres_, 0.0);
    
    Bool_t hasJacb = (jacobians != nullptr && jacobians[0] != nullptr);
    if (hasJacb) std::fill_n(jacobians[0], nres_ *       1, 0.0);
    if (hasJacb) std::fill_n(jacobians[1], nres_ * (nx_-1), 0.0);
    if (hasJacb) std::fill_n(jacobians[2], nres_ * (ny_-1), 0.0);
   
    double wgt_bkg = parameters[0][0];

    std::vector<double> eftx(nx_, 0.0);
    std::vector<double> efty(ny_, 0.0);
    for (int it = 0; it < nx_-1; ++it) eftx[it] = parameters[1][it];
    for (int it = 0; it < ny_-1; ++it) efty[it] = parameters[2][it];
    eftx[nx_-1] = 1.0;
    efty[ny_-1] = 1.0;

    std::vector<double> res(nres_, 0.0);

    std::vector<double> jacb_wgt(nres_, 0.0);
    std::vector<std::vector<double>> jacb_eftx(nres_, std::vector<double>(nx_-1, 0.0));
    std::vector<std::vector<double>> jacb_efty(nres_, std::vector<double>(ny_-1, 0.0));

    for (int it = 0; it < nx_; ++it) {
    for (int jt = 0; jt < ny_; ++jt) {
        int idx = it*ny_+jt;
      
        double res_err = std::sqrt(smp_[it][jt][1] * smp_[it][jt][1] + sig_[it][jt][1] * sig_[it][jt][1]) * std::sqrt(ndof_);
        double res_ttl = (smp_[it][jt][0] - (1.0 - wgt_bkg) * sig_[it][jt][0] - wgt_bkg * eftx[it] * efty[jt]);

        double jacb_corr_sig = std::exp(-sig_[it][jt][1] / sig_[it][jt][0]);
        if (!std::isfinite(jacb_corr_sig)) jacb_corr_sig = 0.0;

        res[idx] += res_ttl / res_err;
        if (hasJacb) {
            jacb_wgt[idx] += (sig_[it][jt][0] * jacb_corr_sig - eftx[it] * efty[jt]) / res_err;
            if (it != nx_-1) jacb_eftx[idx][it] += (-1.0 * wgt_bkg * efty[jt]) / res_err;
            if (jt != ny_-1) jacb_efty[idx][jt] += (-1.0 * wgt_bkg * eftx[it]) / res_err;
        }
    }}
    
    for (int it = 1; it < nx_-1; ++it) {
    for (int jt = 1; jt < ny_-1; ++jt) {
        int idx = nx_*ny_ + (it-1)*(ny_-2)+(jt-1);
        
        double error = std::sqrt(smp_[it][jt][1] * smp_[it][jt][1] + sig_[it][jt][1] * sig_[it][jt][1]) * std::sqrt(ndof_);

        double center = eftx[it] * efty[jt];
        std::array<double, 4> corner({
            eftx[it-1] * efty[jt-1], eftx[it-1] * efty[jt+1],
            eftx[it+1] * efty[jt-1], eftx[it+1] * efty[jt+1]
        });
        std::array<double, 4> side({
            eftx[it]   * efty[jt-1], eftx[it]   * efty[jt+1],
            eftx[it-1] * efty[jt]  , eftx[it+1] * efty[jt]  
        });

        res[idx] += (center + (-1.0/12.0) * (corner[0] + corner[1] + corner[2] + corner[3]) + (-1.0/6.0) * (side[0] + side[1] + side[2] + side[3])) / error;
        
        if (hasJacb) {
            if (it   != nx_-1) jacb_eftx[idx][it]   += (efty[jt] + (-1.0/6.0) * (efty[jt-1] + efty[jt+1])) / error;
            if (it-1 != nx_-1) jacb_eftx[idx][it-1] += ((-1.0/12.0) * (efty[jt-1] + efty[jt+1]) + (-1.0/6.0) * (efty[jt])) / error;
            if (it+1 != nx_-1) jacb_eftx[idx][it+1] += ((-1.0/12.0) * (efty[jt-1] + efty[jt+1]) + (-1.0/6.0) * (efty[jt])) / error;
            
            if (jt   != ny_-1) jacb_efty[idx][jt]   += (eftx[jt] + (-1.0/6.0) * (eftx[jt-1] + eftx[jt+1])) / error;
            if (jt-1 != ny_-1) jacb_efty[idx][jt-1] += ((-1.0/12.0) * (eftx[jt-1] + eftx[jt+1]) + (-1.0/6.0) * (eftx[jt])) / error;
            if (jt+1 != ny_-1) jacb_efty[idx][jt+1] += ((-1.0/12.0) * (eftx[jt-1] + eftx[jt+1]) + (-1.0/6.0) * (eftx[jt])) / error;
        }
    }}

    for (int ii = 0; ii < nres_; ++ii) {
        residuals[ii] = res[ii];
    }

    if (hasJacb) {
        for (int it = 0; it < nx_; ++it) { 
        for (int jt = 0; jt < ny_; ++jt) {
            int idx = it*ny_+jt;
            jacobians[0][idx] = jacb_wgt[idx];
            for (int ii = 0; ii < nx_-1; ++ii) jacobians[1][idx*(nx_-1)+ii] = jacb_eftx[idx][ii];
            for (int ii = 0; ii < ny_-1; ++ii) jacobians[2][idx*(ny_-1)+ii] = jacb_efty[idx][ii];
        }}
        for (int it = 1; it < nx_-1; ++it) { 
        for (int jt = 1; jt < ny_-1; ++jt) {
            int idx = nx_*ny_ + (it-1)*(ny_-2)+(jt-1);
            for (int ii = 0; ii < nx_-1; ++ii) jacobians[1][idx*(nx_-1)+ii] = jacb_eftx[idx][ii];
            for (int ii = 0; ii < ny_-1; ++ii) jacobians[2][idx*(ny_-1)+ii] = jacb_efty[idx][ii];
        }}
    }

    return true;
}

#endif // __DataDriven2D_C__
