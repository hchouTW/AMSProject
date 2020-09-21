#ifndef __DataFit1D_C__
#define __DataFit1D_C__

#include "DataFit1D.h"


ResultFit1D::ResultFit1D(
    const std::string& name, int ndim, int ntmp,
    const std::vector<std::vector<std::array<double, 2>>>& tmps, 
    double num_smp, const std::vector<std::array<double, 2>>& smp, 
    double num_ref, const std::vector<std::array<double, 2>>& ref) {
    clear();
    data_name_ = name;
    data_ndim_ = ndim;
    data_ntmp_ = ntmp;
    data_tmps_ = tmps;
    
    num_smp_  = num_smp;
    data_smp_ = smp;

    if (num_ref > 0.0 && ref.size() == data_ndim_) {
        num_ref_ = num_ref;
        data_ref_ = ref;
    }
}

void ResultFit1D::clear() {
    hist_build_ = false;
    hist_axis_  = TAxis();

    data_name_ = "";
    data_ndim_ = 0;
    data_ntmp_ = 0;

    data_tmps_.clear();
    data_smp_.clear();
    data_ref_.clear();

    wgt_tmps_.clear();
    wgt_errs_.clear();
    
    num_tmps_.clear();
    num_errs_.clear();
    num_sum_ = 0;
    num_smp_ = 0;
    num_ref_ = 0;

    tmps_.clear();
    sum_.clear();
    smp_.clear();
    ref_.clear();

    ndof_ = 0;
    nchi_ = 0;
}

void ResultFit1D::set_fitting_result(int ndof, double nchi, const std::vector<double>& wgt_tmps, const std::vector<double>& wgt_errs) {
    if (wgt_tmps.size() != data_ntmp_) return;
    if (wgt_errs.size() != data_ntmp_) return;

    ndof_ = ndof;
    nchi_ = nchi;

    wgt_tmps_ = wgt_tmps;
    wgt_errs_ = wgt_errs;
    
    num_sum_ = 0.0;
    num_tmps_ = std::vector<double>(data_ntmp_, 0.0);
    num_errs_ = std::vector<double>(data_ntmp_, 0.0);
    for (int itmp = 0; itmp < data_ntmp_; ++itmp) {
        num_tmps_.at(itmp) = wgt_tmps_.at(itmp) * num_smp_;
        num_errs_.at(itmp) = wgt_errs_.at(itmp) * num_smp_;
        num_sum_ += num_tmps_.at(itmp);
    }
    
    tmps_.clear();
    for (int itmp = 0; itmp < data_ntmp_; ++itmp) {
        tmps_.push_back(std::vector<std::array<double, 2>>(data_ndim_, { 0.0, 0.0 }));
        for (int it = 0; it < data_ndim_; ++it) {
            tmps_.back()[it][0] = num_tmps_.at(itmp) * data_tmps_.at(itmp).at(it)[0];
            tmps_.back()[it][1] = num_tmps_.at(itmp) * data_tmps_.at(itmp).at(it)[1];
            if (tmps_.back().at(it)[0] <= 0.0) tmps_.back().at(it)[1] = 0.0;
        }
    }

    sum_ = std::vector<std::array<double, 2>>(data_ndim_, { 0.0, 0.0 });
    for (int it = 0; it < data_ndim_; ++it) {
        double sum_val = 0.0;
        double sum_err = 0.0;
        for (int itmp = 0; itmp < data_ntmp_; ++itmp) {
            sum_val += tmps_.at(itmp).at(it)[0];
            sum_err += tmps_.at(itmp).at(it)[1] * tmps_.at(itmp).at(it)[1];
        }
        sum_err = std::sqrt(sum_err);
        sum_.at(it)[0] = sum_val;
        sum_.at(it)[1] = sum_err;
        
        if (sum_.at(it)[0] <= 0.0) sum_.at(it)[1] = 0.0;
    }
    
    smp_ = std::vector<std::array<double, 2>>(data_ndim_, { 0.0, 0.0 });
    for (int it = 0; it < data_ndim_; ++it) {
        smp_.at(it)[0] = num_smp_ * data_smp_.at(it)[0];
        smp_.at(it)[1] = num_smp_ * data_smp_.at(it)[1];
        if (smp_.at(it)[0] <= 0.0) smp_.at(it)[1] = 0.0;
    }

    if (num_ref_ > 0.0) {
        ref_ = std::vector<std::array<double, 2>>(data_ndim_, { 0.0, 0.0 });
        for (int it = 0; it < data_ndim_; ++it) {
            ref_.at(it)[0] = num_ref_ * data_ref_.at(it)[0];
            ref_.at(it)[1] = num_ref_ * data_ref_.at(it)[1];
            if (ref_.at(it)[0] <= 0.0) ref_.at(it)[1] = 0.0;
        }
    }
}
        

void ResultFit1D::set_histogram_result(bool hist_build, TAxis* hist_axis, const std::string& hist_xtitle, const std::string& hist_ytitle) {
    hist_tmps_.clear();
    if (hist_build && hist_axis != nullptr) {
        hist_build_ = hist_build;
        hist_axis_  = (*hist_axis);
    }
    if (!hist_build_) return;
    
    for (int itmp = 0; itmp < data_ntmp_; ++itmp) {
        std::shared_ptr<TH1D> hist_tmp_ = std::make_shared<TH1D>(Form("hFit1D_tmp%d_%s", itmp, data_name_.c_str()), "", hist_axis_.GetNbins(), hist_axis_.GetXbins()->GetArray());
        for (int it = 0; it < data_ndim_; ++it) {
            hist_tmp_->SetBinContent(it+1, tmps_.at(itmp).at(it)[0]);
            hist_tmp_->SetBinError  (it+1, tmps_.at(itmp).at(it)[1]);
        }
        hist_tmp_->GetXaxis()->SetTitle(hist_xtitle.c_str());
        hist_tmp_->GetYaxis()->SetTitle(hist_ytitle.c_str());
        hist_tmps_.push_back(hist_tmp_);
    }
    
    hist_sum_ = std::make_shared<TH1D>(Form("hFit1D_sum_%s", data_name_.c_str()), "", hist_axis_.GetNbins(), hist_axis_.GetXbins()->GetArray());
    for (int it = 0; it < data_ndim_; ++it) {
        hist_sum_->SetBinContent(it+1, sum_.at(it)[0]);
        hist_sum_->SetBinError  (it+1, sum_.at(it)[1]);
    }
    hist_sum_->GetXaxis()->SetTitle(hist_xtitle.c_str());
    hist_sum_->GetYaxis()->SetTitle(hist_ytitle.c_str());
    
    hist_smp_ = std::make_shared<TH1D>(Form("hFit1D_smp_%s", data_name_.c_str()), "", hist_axis_.GetNbins(), hist_axis_.GetXbins()->GetArray());
    for (int it = 0; it < data_ndim_; ++it) {
        hist_smp_->SetBinContent(it+1, smp_.at(it)[0]);
        hist_smp_->SetBinError  (it+1, smp_.at(it)[1]);
    }
    hist_smp_->GetXaxis()->SetTitle(hist_xtitle.c_str());
    hist_smp_->GetYaxis()->SetTitle(hist_ytitle.c_str());
    
    if (num_ref_ > 0.0) {
        hist_ref_ = std::make_shared<TH1D>(Form("hFit1D_ref_%s", data_name_.c_str()), "", hist_axis_.GetNbins(), hist_axis_.GetXbins()->GetArray());
        for (int it = 0; it < data_ndim_; ++it) {
            hist_ref_->SetBinContent(it+1, ref_.at(it)[0]);
            hist_ref_->SetBinError  (it+1, ref_.at(it)[1]);
        }
        hist_ref_->GetXaxis()->SetTitle(hist_xtitle.c_str());
        hist_ref_->GetYaxis()->SetTitle(hist_ytitle.c_str());
    }
}

DataFit1D::DataFit1D(const std::string& name, std::vector<TH1D*> htmps, TH1D* hsmp, TH1D* href, bool hist_build, const std::string& hist_xtitle, const std::string& hist_ytitle) {
    clear();

    if (htmps.size() == 0) return;
    for (auto&& htmp : htmps) { if (htmp == nullptr) return; }
    if (hsmp == nullptr) return;

    hist_xtitle_ = hist_xtitle;
    hist_ytitle_ = hist_ytitle;

    // Init
    for (auto&& htmp : htmps) {
        htmps_.push_back(htmp);
        auto&& htmp_data = build_data(htmps_.back());
        data_num_tmps_.push_back(std::get<0>(htmp_data));
        data_tmps_    .push_back(std::get<1>(htmp_data));
    }

    hsmp_ = hsmp;
    auto&& hsmp_data = build_data(hsmp_);
    data_num_smp_ = std::get<0>(hsmp_data);
    data_smp_     = std::get<1>(hsmp_data);

    if (href != nullptr) {
        href_ = href;
        auto&& href_data = build_data(href_);
        data_num_ref_ = std::get<0>(href_data);
        data_ref_     = std::get<1>(href_data);
    }

    if (hist_build) {
        hist_build_ = hist_build;
        hist_axis_  = (*(hsmp_->GetXaxis()));
    }
    
    data_name_ = name;
    data_ndim_ = data_smp_.size();
    data_ntmp_ = data_tmps_.size();

    if (data_ndim_ == 0) { clear(); return; }
    if (data_ntmp_ == 0) { clear(); return; }

    // Check
    int count_not_empty_tmps = 0;
    for (int it = 0; it < data_ntmp_; ++it) {
        if (data_tmps_.at(it).size() != data_ndim_) { clear(); return; }
        if (data_num_tmps_.at(it) > 0.0) count_not_empty_tmps++;
    }
    if (count_not_empty_tmps == 0) { clear(); return; }
    
    if (data_num_smp_ <= 0.0) { clear(); return; }
    if (data_smp_.size() != data_ndim_) { clear(); return; }
    
    if (href_ != nullptr && data_num_ref_ <= 0.0) { clear(); return; }
    if (href_ != nullptr && data_ref_.size() != data_ndim_) { clear(); return; }

    if (!fit()) clear();
}

void DataFit1D::clear() {
    hist_build_ = false;
    hist_axis_ = TAxis();
    hist_xtitle_ = "";
    hist_ytitle_ = "";

    data_name_ = "";
    data_ndim_ = 0;
    data_ntmp_ = 0;

    htmps_.clear();
    hsmp_ = nullptr;
    href_ = nullptr;
   
    data_num_tmps_.clear();
    data_num_smp_ = 0.0;
    data_num_ref_ = 0.0;

    data_tmps_.clear();
    data_smp_.clear();
    data_ref_.clear();

    result_ = ResultFit1D();
}

std::tuple<double, std::vector<std::array<double, 2>>> DataFit1D::build_data(TH1D* hist) {
    std::tuple<double, std::vector<std::array<double, 2>>> datas;
    std::get<0>(datas) = 0.0;

    if (hist == nullptr) return datas;
    std::get<1>(datas) = std::vector<std::array<double, 2>>(hist->GetXaxis()->GetNbins(), {0.0, 0.0});
    std::vector<std::array<double, 2>> data(hist->GetXaxis()->GetNbins(), {0.0, 0.0});

    double total = 0.0;
    for (int ib = 1; ib <= hist->GetXaxis()->GetNbins(); ++ib) {
        if (hist->GetBinContent(ib) <= 0.0) continue;
        data[ib-1][0] = hist->GetBinContent(ib);
        data[ib-1][1] = hist->GetBinError  (ib);
        total += data[ib-1][0];
    }
    //if (total <= 1.001) return datas;

    double error = (total > 0.0) ? 1.0 / total : 0.0;
    for (auto&& arr : data) {
        arr[0] /= total;
        arr[1] /= total;
        if (arr[0] <= 0.0) arr[0] = 0.0;
        if (arr[1] <= 0.0) arr[1] = error;
    }

    std::get<0>(datas) = total;
    std::get<1>(datas) = data;
    return datas;
}

bool DataFit1D::fit() {
    std::vector<double> params(data_ntmp_, 1.0/static_cast<double>(data_ntmp_));
    for (int it = 0; it < params.size(); ++it) {
        if (data_num_tmps_.at(it) <= 0.0) params.at(it) = 0.0;
    }

    // CeresSolver: Cost Function
    ceres::CostFunction* cost_function = new VirtualDataFit1D(data_ndim_, data_ntmp_, data_tmps_, data_smp_);

    // CeresSolver: Loss Function
    ceres::LossFunction* loss_function = new ceres::HuberLoss(3.0);

    // CeresSolver: Problem
    ceres::Problem problem;
    problem.AddResidualBlock(cost_function, loss_function, params.data());

    for (int it = 0; it < params.size(); ++it) {
        problem.SetParameterLowerBound(params.data(), it, 0.0);
        problem.SetParameterUpperBound(params.data(), it, 1.5);
        if (data_num_tmps_.at(it) <= 0.0) {
            problem.SetParameterLowerBound(params.data(), it, 0.0);
            problem.SetParameterUpperBound(params.data(), it, 0.0);
        }
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
    //std::cerr << summary.FullReport() << std::endl;
        
    for (int it = 0; it < params.size(); ++it) {
        if (!std::isfinite(params.at(it))) return false;
        if (data_num_tmps_.at(it) <= 0.0) params.at(it) = 0.0;
    }

    ceres::Covariance::Options covariance_options;
    ceres::Covariance covariance(covariance_options);

    std::vector<std::pair<const double*, const double*> > covariance_blocks;
    covariance_blocks.push_back(std::make_pair(params.data(), params.data()));

    CHECK(covariance.Compute(covariance_blocks, &problem));

    std::vector<double> covariance_params(params.size() * params.size(), 0.0);
    covariance.GetCovarianceBlock(params.data(), params.data(), covariance_params.data());

    std::vector<double> errors(params.size(), 0.0);
    for (int it = 0; it < params.size(); ++it) {
        errors[it] = std::sqrt(covariance_params[it*params.size()+it]);
        if (!std::isfinite(errors[it])) errors[it] = 0.0;
    }

    int    ndof = static_cast<VirtualDataFit1D*>(cost_function)->ndof();
    double nchi = 2.0 * summary.final_cost / static_cast<double>(ndof);

    ResultFit1D result(data_name_, data_ndim_, data_ntmp_, data_tmps_, data_num_smp_, data_smp_, data_num_ref_, data_ref_);
    result.set_fitting_result(ndof, nchi, params, errors);
    result.set_histogram_result(hist_build_, &hist_axis_, hist_xtitle_, hist_ytitle_);

    result_ = result;
    return true;
}
        
bool VirtualDataFit1D::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const {
    if (ndof_ < 1) return false;
    std::fill_n(residuals, nres_, 0.0);
    
    Bool_t hasJacb = (jacobians != nullptr && jacobians[0] != nullptr);
    if (hasJacb) std::fill_n(jacobians[0], nres_ * npar_, 0.0);
  
    std::vector<double> wgt_tmps(npar_, 0.0);
    for (int it = 0; it < npar_; ++it)
        wgt_tmps[it] = parameters[0][it];

    std::vector<double> res (nres_, 0.0);
    std::vector<double> jacb(nres_ * npar_, 0.0);

    for (int it = 0; it < nres_; ++it) {
        if (!has_[it]) continue;
        
        double sum_tmps = 0.0;
        for (int itmp = 0; itmp < ntmp_; ++itmp) sum_tmps += wgt_tmps[itmp] * tmps_.at(itmp).at(it)[0];

        res[it] += (smp_[it][0] - sum_tmps) / smp_[it][1];
        if (hasJacb) {
            for (int itmp = 0; itmp < ntmp_; ++itmp) {
                if (!has_tmps_[itmp][it]) continue;
                jacb[it*npar_+itmp] += (-tmps_[itmp][it][0] / smp_[it][1]);
            }
        }
    }

    for (int it = 0; it < nres_; ++it) {
        if (!std::isfinite(res[it])) continue;
        residuals[it] = res[it];
    }

    if (hasJacb) {
        for (int it = 0; it < nres_ * npar_; ++it) { 
            if (!std::isfinite(jacb[it])) continue;
            jacobians[0][it] = jacb[it];
        }
    }

    return true;
}

#endif // __DataFit1D_C__
