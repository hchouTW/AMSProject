#ifndef __DataFit1D_C__
#define __DataFit1D_C__

#include "DataFit1D.h"

ResultFit1D::ResultFit1D(
    const std::string& name, int ndim,
    const std::vector<std::array<double, 2>>& smp, 
    const std::vector<std::array<double, 2>>& sig, 
    const std::vector<std::array<double, 2>>& bkg,
    double num_ref, const std::vector<std::array<double, 2>>& ref) {
    clear();
    name_ = name;
    data_ndim_ = ndim;
    data_smp_ = smp;
    data_sig_ = sig;
    data_bkg_ = bkg;
    if (num_ref > 0.0 && ref.size() == data_ndim_) {
        num_ref_  = num_ref;
        data_ref_ = ref;
    }
    else {
        num_ref_  = 0.0;
        data_ref_ = std::vector<std::array<double, 2>>(data_ndim_, { 0.0, 0.0 });
    }
}

void ResultFit1D::clear() {
    name_ = "";
    data_ndim_ = 0;
    data_smp_.clear();
    data_sig_.clear();
    data_bkg_.clear();
    data_ref_.clear();
    ndof_ = 0;
    nchi_ = 0;
    wgt_sig_ = 0;
    wgt_bkg_ = 0;
    wgt_err_ = 0;
    wgt_stat_err_ = 0;
    wgt_syst_err_ = 0;
    num_ref_ = 0;
    num_smp_ = 0;
    num_sig_ = 0;
    num_bkg_ = 0;
    num_err_ = 0;
    num_stat_err_ = 0;
    num_syst_err_ = 0;
    smp_.clear();
    sum_.clear();
    sig_.clear();
    bkg_.clear();
    ref_.clear();
}

void ResultFit1D::evolve(
    double num_smp, double wgt_sig, double wgt_bkg, double wgt_syst_err,
    bool build_hist, TAxis* axis) {

    std::vector<bool> has     = std::vector<bool>(data_ndim_, false);
    std::vector<bool> has_smp = std::vector<bool>(data_ndim_, false);
    std::vector<bool> has_sig = std::vector<bool>(data_ndim_, false);
    std::vector<bool> has_bkg = std::vector<bool>(data_ndim_, false);
    
    int ndof = -1;
    for (int it = 0; it < data_ndim_; ++it) {
        if (data_smp_[it][0] <= 0.0 && data_sig_[it][0] <= 0.0 && data_bkg_[it][0] <= 0.0) continue;
        if (data_smp_[it][0] > 0.0) has_smp[it] = true;
        if (data_sig_[it][0] > 0.0) has_sig[it] = true;
        if (data_bkg_[it][0] > 0.0) has_bkg[it] = true;
        has[it] = true;
        ndof++;
    }
   
    std::vector<double> res (data_ndim_, 0.0);
    std::vector<double> jacb(data_ndim_, 0.0);
    for (int it = 0; it < data_ndim_; ++it) {
        if (!has[it]) continue;
        
        double dif = (data_smp_[it][0] - wgt_sig * data_sig_[it][0] - wgt_bkg * data_bkg_[it][0]);

        double err_smp = data_smp_[it][1];
        double err_sig = has_sig[it] ? data_sig_[it][1] : 0.0;
        double err_bkg = has_bkg[it] ? data_bkg_[it][1] : 0.0;
        double err     = std::sqrt(err_smp * err_smp + err_sig * err_sig + err_bkg * err_bkg);

        res[it]  += (dif / err);
        jacb[it] += (-data_sig_[it][0] / err);
        jacb[it] += (+data_bkg_[it][0] / err);
    }

    double chi2 = 0.0;
    for (auto&& elem : res) chi2 += elem * elem;
    double nchi = chi2 / static_cast<double>(ndof);

    double conv = 0.0;
    for (auto&& elem : jacb) conv += elem * elem;
    double wgt_stat_err = (!std::isfinite(conv) || conv == 0.0) ? 0.0 : (1.0 / std::sqrt(conv));
 
    double wgt_err = std::sqrt(wgt_stat_err * wgt_stat_err + wgt_syst_err * wgt_syst_err);

    ndof_ = ndof;
    nchi_ = nchi;

    wgt_sig_ = wgt_sig;
    wgt_bkg_ = wgt_bkg;
    wgt_err_ = wgt_err;
    wgt_stat_err_ = wgt_stat_err;
    wgt_syst_err_ = wgt_syst_err;

    num_smp_ = num_smp;
    num_sig_ = num_smp * wgt_sig;
    num_bkg_ = num_smp * wgt_bkg;
    num_err_ = num_smp * wgt_err;
    num_stat_err_ = num_smp * wgt_stat_err;
    num_syst_err_ = num_smp * wgt_syst_err;
    
    smp_ = std::vector<std::array<double, 2>>(data_ndim_, { 0.0, 0.0 });
    sum_ = std::vector<std::array<double, 2>>(data_ndim_, { 0.0, 0.0 });
    sig_ = std::vector<std::array<double, 2>>(data_ndim_, { 0.0, 0.0 });
    bkg_ = std::vector<std::array<double, 2>>(data_ndim_, { 0.0, 0.0 });
    ref_ = std::vector<std::array<double, 2>>(data_ndim_, { 0.0, 0.0 });
    for (int it = 0; it < data_ndim_; ++it) {
        smp_.at(it)[0] = num_smp_ * data_smp_.at(it)[0];
        smp_.at(it)[1] = num_smp_ * data_smp_.at(it)[1];
        sig_.at(it)[0] = num_sig_ * data_sig_.at(it)[0];
        sig_.at(it)[1] = num_sig_ * data_sig_.at(it)[1];
        bkg_.at(it)[0] = num_bkg_ * data_bkg_.at(it)[0];
        bkg_.at(it)[1] = num_bkg_ * data_bkg_.at(it)[1];
        if (num_ref_ != 0.0) ref_.at(it)[0] = num_ref_ * data_ref_.at(it)[0];
        if (num_ref_ != 0.0) ref_.at(it)[1] = num_ref_ * data_ref_.at(it)[1];

        sum_.at(it)[0] = sig_.at(it)[0] + bkg_.at(it)[0];
        sum_.at(it)[1] = std::sqrt(sig_.at(it)[1] * sig_.at(it)[1] + bkg_.at(it)[1] * bkg_.at(it)[1]);

        if (smp_.at(it)[0] <= 0.0) smp_.at(it)[1] = 0.0;
        if (sig_.at(it)[0] <= 0.0) sig_.at(it)[1] = 0.0;
        if (bkg_.at(it)[0] <= 0.0) bkg_.at(it)[1] = 0.0;
        if (sum_.at(it)[0] <= 0.0) sum_.at(it)[1] = 0.0;
        if (ref_.at(it)[0] <= 0.0) ref_.at(it)[1] = 0.0;
    }
   
    if (build_hist && axis != nullptr) {
        hist_smp_ = std::make_shared<TH1D>(Form("hFit1D_smp_%s", name_.c_str()), "", axis->GetNbins(), axis->GetXbins()->GetArray());
        hist_sum_ = std::make_shared<TH1D>(Form("hFit1D_sum_%s", name_.c_str()), "", axis->GetNbins(), axis->GetXbins()->GetArray());
        hist_sig_ = std::make_shared<TH1D>(Form("hFit1D_sig_%s", name_.c_str()), "", axis->GetNbins(), axis->GetXbins()->GetArray());
        hist_bkg_ = std::make_shared<TH1D>(Form("hFit1D_bkg_%s", name_.c_str()), "", axis->GetNbins(), axis->GetXbins()->GetArray());
        hist_ref_ = std::make_shared<TH1D>(Form("hFit1D_ref_%s", name_.c_str()), "", axis->GetNbins(), axis->GetXbins()->GetArray());
        
        for (int it = 0; it < data_ndim_; ++it) {
            hist_smp_->SetBinContent(it+1, smp_.at(it)[0]);
            hist_smp_->SetBinError  (it+1, smp_.at(it)[1]);
            hist_sum_->SetBinContent(it+1, sum_.at(it)[0]);
            hist_sum_->SetBinError  (it+1, sum_.at(it)[1]);
            hist_sig_->SetBinContent(it+1, sig_.at(it)[0]);
            hist_sig_->SetBinError  (it+1, sig_.at(it)[1]);
            hist_bkg_->SetBinContent(it+1, bkg_.at(it)[0]);
            hist_bkg_->SetBinError  (it+1, bkg_.at(it)[1]);
            hist_ref_->SetBinContent(it+1, ref_.at(it)[0]);
            hist_ref_->SetBinError  (it+1, ref_.at(it)[1]);
        }
    }
}

DataFit1D::DataFit1D(const std::string& name, TH1D* hsmp, TH1D* hsig, TH1D* hbkg, TH1D* href) {
    clear();
    if (hsmp == nullptr || hsig == nullptr || hbkg == nullptr) return;
   
    name_ = name;

    hsmp_ = hsmp;
    auto&& hsmp_data = build_data(hsmp_);
    data_num_smp_ = std::get<0>(hsmp_data);
    data_smp_     = std::get<1>(hsmp_data);
    
    hsig_ = hsig;
    auto&& hsig_data = build_data(hsig_);
    data_num_sig_ = std::get<0>(hsig_data);
    data_sig_     = std::get<1>(hsig_data);

    hbkg_ = hbkg;
    auto&& hbkg_data = build_data(hbkg_);
    data_num_bkg_ = std::get<0>(hbkg_data);
    data_bkg_     = std::get<1>(hbkg_data);

    if (href != nullptr) {
        href_ = href;
        auto&& href_data = build_data(href_);
        data_num_ref_ = std::get<0>(href_data);
        data_ref_     = std::get<1>(href_data);
    }

    data_ndim_ = data_smp_.size();

    if (!evolve()) clear();
}

void DataFit1D::clear() {
    name_ = "";
    hsmp_ = nullptr;
    hsig_ = nullptr;
    hbkg_ = nullptr;
    href_ = nullptr;
    
    data_ndim_ = 0;
    
    data_num_smp_ = 0.0;
    data_num_sig_ = 0.0;
    data_num_bkg_ = 0.0;
    data_num_ref_ = 0.0;

    data_smp_.clear();
    data_sig_.clear();
    data_bkg_.clear();
    data_ref_.clear();

    result_ = ResultFit1D();
}

std::tuple<double, std::vector<std::array<double, 2>>> DataFit1D::build_data(TH1D* hist) {
    std::tuple<double, std::vector<std::array<double, 2>>> datas;
    std::get<0>(datas) = 0.0;

    if (hist == nullptr) return datas;
    std::vector<std::array<double, 2>> data(hist->GetXaxis()->GetNbins(), {0.0, 0.0});

    double sum_val = 0.0;
    for (int ib = 1; ib <= hist->GetXaxis()->GetNbins(); ++ib) {
        if (hist->GetBinContent(ib) <= 0.0) continue;
        data[ib-1][0] = hist->GetBinContent(ib);
        data[ib-1][1] = hist->GetBinError  (ib);
        sum_val += data[ib-1][0];
    }
   
    double total = sum_val;
    double error = std::sqrt(1.0 / total);
    if (total > 0.0) {
        for (auto&& arr : data) {
            arr[0] /= total;
            arr[1] /= total;
            if (arr[0] <= 0.0) arr[0] = 0.0;
            if (arr[1] <= 0.0) arr[1] = 1.0 / total;
        }
    }

    std::get<0>(datas) = total;
    std::get<1>(datas) = data;
    return datas;
}

ResultFit1D DataFit1D::fit(
    const std::string& name,
    const std::vector<std::array<double, 2>>& smp, 
    const std::vector<std::array<double, 2>>& sig, 
    const std::vector<std::array<double, 2>>& bkg,
    bool build_hist) {
    if (data_num_smp_ == 0) return ResultFit1D();
    
    std::vector<double> params({ 0.5 });

    // CeresSolver: Cost Function
    ceres::CostFunction* cost_function = new VirtualDataFit1D(data_ndim_, smp, sig, bkg);

    // CeresSolver: Loss Function
    ceres::LossFunction* loss_function = new ceres::HuberLoss(1.0);

    // CeresSolver: Problem
    ceres::Problem problem;
    problem.AddResidualBlock(cost_function, loss_function, params.data());

    problem.SetParameterLowerBound(params.data(), 0, 0.0);
    problem.SetParameterUpperBound(params.data(), 0, 1.0);

    // CeresSolver: Options
    ceres::Solver::Options options;
    options.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT;
    options.max_num_iterations = 100;

    // CeresSolver: Summary
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    if (!summary.IsSolutionUsable()) return ResultFit1D();
    //if (ceres::NO_CONVERGENCE == summary.termination_type) return false;
    //std::cerr << summary.FullReport() << std::endl;
    
    if (!std::isfinite(params[0])) return ResultFit1D();

    double wgt_sig = params[0];
    double wgt_bkg = 1.0 - wgt_sig;



    ceres::Covariance::Options options2;
    ceres::Covariance covariance(options2);

    std::vector<std::pair<const double*, const double*> > covariance_blocks;
    covariance_blocks.push_back(std::make_pair(params.data(), params.data()));

    CHECK(covariance.Compute(covariance_blocks, &problem));

    double covariance_params[1 * 1];
    covariance.GetCovarianceBlock(params.data(), params.data(), covariance_params);



    ResultFit1D result(name_, data_ndim_, smp, sig, bkg, data_num_ref_, data_ref_);
    result.evolve(data_num_smp_, wgt_sig, wgt_bkg, sqrt(covariance_params[0]), build_hist, hsmp_->GetXaxis());
    
    //std::cerr << Form("COV ERR %14.8f %14.8f    NUM %14.8f\n", result.wgt_stat_err(), sqrt(covariance_params[0]), result.wgt_sig());

    return result;
}
        
bool DataFit1D::evolve() {
    ResultFit1D&& result = fit(name_, data_smp_, data_sig_, data_bkg_, false);
    if (result.ndof() == 0) return false;

    static long rndmSeed = std::chrono::system_clock::now().time_since_epoch().count();
    static std::mt19937_64 rndmEngMT64(rndmSeed);
    std::normal_distribution<double> distribution(0.0, 1.0);
    std::function<double()>&& normal = std::bind(distribution, std::ref(rndmEngMT64));

    int counts = 0;
    int counts_while = 0;
    std::vector<std::array<double, 3>> wgt_syst_fluc;
    while (counts < result.ndof() * 100 && counts_while < result.ndof() * 1000) {
        counts_while++;

        auto&& data_smp = data_smp_; 
        auto&& data_sig = data_sig_; 
        auto&& data_bkg = data_bkg_;
        
        //double sum_num_smp = 0.0;
        //for (auto&& elem : data_smp) {
        //    if (elem[0] <= 0.0) continue;
        //    elem[0] = elem[0] + normal() * elem[1];
        //    if (elem[0] <= 0.0) elem[0] = 0.0;
        //    sum_num_smp += elem[0];
        //}
        //if (sum_num_smp > 0.0) {
        //    for (auto&& elem : data_smp) {
        //        elem[0] /= sum_num_smp;
        //        elem[1] /= sum_num_smp;
        //    }
        //}
        
        double sum_num_sig = 0.0;
        for (auto&& elem : data_sig) {
            if (elem[0] <= 0.0) continue;
            elem[0] = elem[0] + normal() * elem[1];
            if (elem[0] <= 0.0) elem[0] = 0.0;
            sum_num_sig += elem[0];
        }
        if (sum_num_sig > 0.0) {
            for (auto&& elem : data_sig) {
                elem[0] /= sum_num_sig;
                elem[1] /= sum_num_sig;
            }
        }
        
        double sum_num_bkg = 0.0;
        for (auto&& elem : data_bkg) {
            if (elem[0] <= 0.0) continue;
            elem[0] = elem[0] + normal() * elem[1];
            if (elem[0] <= 0.0) elem[0] = 0.0;
            sum_num_bkg += elem[0];
        }
        if (sum_num_bkg > 0.0) {
            for (auto&& elem : data_bkg) {
                elem[0] /= sum_num_bkg;
                elem[1] /= sum_num_bkg;
            }
        }
    
        auto&& sub_result = fit("", data_smp, data_sig, data_bkg, false);
        if (sub_result.ndof() == 0) continue;
        
        wgt_syst_fluc.push_back(std::array<double, 3>({ sub_result.nchi(), sub_result.wgt_sig(), sub_result.wgt_stat_err() }));
        counts++;
    }
    if (wgt_syst_fluc.size() < 50 * result.ndof()) return false;

    double sum_ww = 0.0;
    double sum_x1 = 0.0;
    double sum_x2 = 0.0;
    for (auto&& elem : wgt_syst_fluc) {
        double wgt = 1.0 / elem[2] / elem[2];
        sum_ww += wgt;
        sum_x1 += wgt * elem[1];
        sum_x2 += wgt * elem[1] * elem[1];
    }
    sum_x1 /= sum_ww;
    sum_x2 /= sum_ww;
    double syst_err_corr = std::sqrt(static_cast<double>(wgt_syst_fluc.size()) / (static_cast<double>(wgt_syst_fluc.size()) - 1.0));
    double wgt_syst_err = syst_err_corr * std::sqrt(sum_x2 - sum_x1 * sum_x1);

    ResultFit1D rlt(name_, result.data_ndim(), result.data_smp(), result.data_sig(), result.data_bkg(), result.num_ref(), result.data_ref());
    //rlt.evolve(result.num_smp(), result.wgt_sig(), result.wgt_bkg(), wgt_syst_err, true, hsmp_->GetXaxis());
    rlt.evolve(result.num_smp(), result.wgt_sig(), result.wgt_bkg(), result.wgt_syst_err(), true, hsmp_->GetXaxis());
    result_ = rlt;
    
    return true;
}
        
bool VirtualDataFit1D::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const {
    if (ndof_ < 1) return false;
    std::fill_n(residuals, nres_, 0.0);
    
    Bool_t hasJacb = (jacobians != nullptr && jacobians[0] != nullptr);
    if (hasJacb) std::fill_n(jacobians[0], nres_ * npar_, 0.0);
   
    double wgt_sig = parameters[0][0];
    double wgt_bkg = 1.0 - wgt_sig;

    std::vector<double> res (nres_, 0.0);
    std::vector<double> jacb(nres_ * npar_, 0.0);

    //double scale = std::sqrt(1.0 / static_cast<double>(ndof_));
    double scale = 1.0;
    for (int it = 0; it < nres_; ++it) {
        if (!has_[it]) continue;
        
        double dif = (smp_[it][0] - wgt_sig * sig_[it][0] - wgt_bkg * bkg_[it][0]);

        double err_smp = smp_[it][1];
        double err_sig = has_sig_[it] ? sig_[it][1] : 0.0;
        double err_bkg = has_bkg_[it] ? bkg_[it][1] : 0.0;
        //double err     = std::sqrt(err_smp * err_smp + err_sig * err_sig + err_bkg * err_bkg);
        double err     = std::sqrt(err_smp * err_smp + wgt_sig * wgt_sig * err_sig * err_sig + wgt_bkg * wgt_bkg * err_bkg * err_bkg);

        res[it] += scale * (dif / err);
        if (hasJacb) {
            jacb[it*npar_+0] += scale * (-sig_[it][0] / err);
            jacb[it*npar_+0] += scale * (+bkg_[it][0] / err);
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
