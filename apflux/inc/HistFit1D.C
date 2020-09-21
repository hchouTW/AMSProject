#ifndef __HistFit1D_C__
#define __HistFit1D_C__

#include "HistFit1D.h"

namespace HistFit {

Data1D::Data1D(const std::string name, const std::string title, const std::vector<Point1D> data) : Data1D(name, title) {
    if (data.size() == 0) return;
    for (int i = 0; i < data.size(); ++i) data_.push_back(data.at(i));
    build();
}

bool Data1D::build() {
    nbins_ = 0;
    sum_ = 0.0;
    pdf_.clear();
    
    if (data_.size() == 0) return false;
    nbins_ = data_.size();
    
    for (int i = 0; i < data_.size(); ++i) sum_ += data_.at(i).val;
    if (!std::isfinite(sum_)) sum_ = 0.0;
    
    if (sum_ == 0.0) {
        data_ = std::vector<Point1D>(nbins_);
        pdf_  = std::vector<Point1D>(nbins_);
    }
    else {
        for (int i = 0; i < data_.size(); ++i) { 
            pdf_.push_back(Point1D(data_.at(i).val / sum_, data_.at(i).err / sum_));
        }
    }
    return true;
}

std::vector<Data1D> Data1D::build_fluc(int nsample, int method) const {
    std::vector<Data1D> fluc_data;
    if (nsample <= 0) return fluc_data;

    long total = static_cast<long>(sum_);
    if (nbins_ == 0 || total < 1) return std::vector<Data1D>(nsample, Data1D(name_, title_, std::vector<Point1D>(nbins_)));

    std::vector<std::function<long()>> funcs;
    std::mt19937_64 generatorMT64(std::chrono::system_clock::now().time_since_epoch().count());

    if (method == 0) { // binomial
        for (int ib = 0; ib < nbins_; ++ib) {
            std::binomial_distribution<long> distribution(total, pdf_.at(ib).val);
            std::function<long()>&& func = std::bind(distribution, std::ref(generatorMT64));
            funcs.push_back(func);
        }
    }
    else { // poisson
        for (int ib = 0; ib < nbins_; ++ib) {
            std::poisson_distribution<long> distribution(data_.at(ib).val);
            std::function<long()>&& func = std::bind(distribution, std::ref(generatorMT64));
            funcs.push_back(func);
        }
    }
    
    for (int isample = 0; isample < nsample; ++isample) {
        double prob = 0.0;
        Data1D data(name_, title_);
        for (int ib = 0; ib < nbins_; ++ib) {
            data.push_back(static_cast<double>(funcs.at(ib)()), data_.at(ib).err);
            double res = (data.data(ib).val - data_.at(ib).val) / data_.at(ib).err;
            if (!std::isfinite(res)) continue;
            prob += (res * res);
        }
        prob = std::exp(-0.5 * prob / static_cast<double>(nbins_));
        data.build();
        fluc_data.push_back(data);
    }

    return fluc_data;
}

DataFit1D::DataFit1D(const Data1D& smp, const std::vector<Data1D>& tmps) {
    if (!build(smp, tmps)) { clear(); return; }
    if (!fit()) { clear(); return; }
    if (!fine()) { clear(); return; }
}

void DataFit1D::clear() {
    smp_  = Data1D();
    tmps_.clear();

    nsmp_ = 0;
    wgts_ = std::vector<double>();
    errs_ = std::vector<double>();
    
    chi_fits_ = Data1D();
    sum_tmps_ = Data1D();
    wgt_tmps_.clear();
    
    nbins_ = 0;
    ntmps_ = 0;
    
    status_ = false;
    nchi_ = 0;
    ndof_ = 0;
    
    bins_.clear();
    reduce_smp_ = Data1D();
    reduce_tmps_.clear();
    reduce_pars_.clear();
}

bool DataFit1D::build(const Data1D& smp, const std::vector<Data1D>& tmps) {
    const double LIMIT = 1.0e-6;
    if (smp.nbins() == 0 || tmps.size() == 0) return false;
    for (int it = 0; it < tmps.size(); ++it) {
        if (smp.nbins() != tmps.at(it).nbins()) return false;
    }
    clear();

    smp_  = smp;
    tmps_ = tmps;
       
    smp_.build();
    for (int it = 0; it < ntmps_; ++it) tmps_.at(it).build();

    nsmp_ = smp_.sum();
    wgts_ = std::vector<double>(tmps_.size(), 0.0);
    errs_ = std::vector<double>(tmps_.size(), 0.0);
    
    nbins_ = smp_.nbins();
    ntmps_ = tmps_.size();

    std::vector<bool> status_bins(nbins_, false);
    for (int ib = 0; ib < nbins_; ++ib) {
        bool status_smp  = (smp_.data(ib).val > LIMIT);
        bool status_tmps = false;
        for (int it = 0; it < ntmps_; ++it) {
            status_tmps = (status_tmps || (tmps_.at(it).data(ib).val > LIMIT));
        }
        bool status_bin = (status_smp && status_tmps);
        status_bins.at(ib) = status_bin;
    }
    int nbins_fitting = std::accumulate(status_bins.begin(), status_bins.end(), 0);
    int ndof = nbins_fitting - ntmps_;

    if (ndof < 1) return false;
    ndof_ = ndof;
    bins_ = status_bins;
   
    reduce_smp_ = Data1D(smp_.name(), smp_.title());
    for (int ib = 0; ib < nbins_; ++ib) {
        if (!status_bins.at(ib)) continue;
        reduce_smp_.push_back(smp_.data(ib));
    }
    reduce_smp_.build();

    for (int it = 0; it < ntmps_; ++it) {
        Data1D reduce_tmp(tmps_.at(it).name(), tmps_.at(it).title());
        for (int ib = 0; ib < nbins_; ++ib) {
            if (!status_bins.at(ib)) continue;
            reduce_tmp.push_back(tmps_.at(it).data(ib));
        }
        reduce_tmp.build();
        if (reduce_tmp.sum() < LIMIT) continue;
        reduce_tmps_.push_back(reduce_tmp);
        reduce_pars_.push_back(it);
    }
    if (reduce_pars_.size() == 0) return false;
    ndof_ += (tmps_.size() - reduce_tmps_.size());

    return true;
}

bool DataFit1D::fit() {
    const double LIMIT = 1.0e-8;
    std::vector<double> params(reduce_pars_.size(), nsmp_ / static_cast<double>(reduce_pars_.size()));

    // CeresSolver: Cost Function
    ceres::CostFunction* cost_function = new VirtualDataFit1D(reduce_smp_, reduce_tmps_);

    // CeresSolver: Loss Function
    //ceres::LossFunction* loss_function = new ceres::HuberLoss(3.0);

    // CeresSolver: Problem
    ceres::Problem problem;
    //problem.AddResidualBlock(cost_function, loss_function, params.data());
    problem.AddResidualBlock(cost_function, nullptr, params.data());

    double bound = (nsmp_ < LIMIT) ? 0.0 : std::sqrt(nsmp_);
    for (int it = 0; it < params.size(); ++it) {
        problem.SetParameterLowerBound(params.data(), it, 0.0);
        problem.SetParameterUpperBound(params.data(), it, nsmp_ + 5.0 * bound);
    }

    // CeresSolver: Options
    ceres::Solver::Options options;
    options.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT;
    options.max_num_iterations = 100;

    // CeresSolver: Summary
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    if (!summary.IsSolutionUsable()) return false;
    if (ceres::NO_CONVERGENCE == summary.termination_type) return false;
    //std::cerr << summary.FullReport() << std::endl;
    
    double sum_wgts = 0.0;
    std::vector<double> tmps_sum(ntmps_, 0.0);
    for (int it = 0; it < params.size(); ++it) {
        int idx = reduce_pars_.at(it);
        double tmp_sum = reduce_tmps_.at(it).sum();
        if (!std::isfinite(params.at(it))) continue;
        tmps_sum.at(idx) = tmp_sum;
        wgts_.at(idx) = params.at(it);
        sum_wgts += wgts_.at(idx);
    }
    if (sum_wgts < LIMIT) return false;

    for (int it = 0; it < wgts_.size(); ++it) {
        errs_.at(it) = std::sqrt(wgts_.at(it));
    }

    for (int it = 0; it < ntmps_; ++it) {
        double wgt_crr = wgts_.at(it) * (tmps_.at(it).sum() / tmps_sum.at(it));
        if (!std::isfinite(wgt_crr) || wgt_crr <= 0.0) wgt_crr = 0.0;
        wgt_tmps_.push_back(
            Data1D(tmps_.at(it).name(), tmps_.at(it).title(), tmps_.at(it).weight_pdf(wgt_crr))
        );
    }

    sum_tmps_ = Data1D("Summation Template", "Summation Template");
    for (int ib = 0; ib < nbins_; ++ib) {
        Point1D pnt;
        for (int it = 0; it < ntmps_; ++it) {
            Point1D ipnt = wgt_tmps_.at(it).data(ib);
            pnt.val += ipnt.val;
            pnt.err += (ipnt.val == 0.0) ? 0.0 : (ipnt.err * ipnt.err);
        }
        pnt.err = std::sqrt(pnt.err);
        sum_tmps_.push_back(pnt);
    }
    sum_tmps_.build();

    chi_fits_ = Data1D("Chi of Fitting", "Chi of Fitting");
    for (int ib = 0; ib < nbins_; ++ib) {
        Point1D pnt_smp = smp_.data(ib);
        Point1D pnt_tmp = sum_tmps_.data(ib);
        double err = std::hypot(pnt_smp.err, pnt_tmp.err);
        double chi = (pnt_tmp.val - pnt_smp.val) / err;
        chi_fits_.push_back( (bins_.at(ib) ? chi : 0.0) );
    }

    double nchi = 0.0;
    for (int ib = 0; ib < nbins_; ++ib) {
        if (!bins_.at(ib)) continue;
        nchi += chi_fits_.data(ib).val * chi_fits_.data(ib).val;
    }
    nchi = std::sqrt(nchi / static_cast<double>(ndof_));
    if (!std::isfinite(nchi)) return false;
    nchi_ = nchi;
   
    ceres::Covariance::Options covariance_options;
    ceres::Covariance covariance(covariance_options);

    std::vector<std::pair<const double*, const double*> > covariance_blocks;
    covariance_blocks.push_back(std::make_pair(params.data(), params.data()));

    CHECK(covariance.Compute(covariance_blocks, &problem));

    std::vector<double> covariance_params(params.size() * params.size(), 0.0);
    covariance.GetCovarianceBlock(params.data(), params.data(), covariance_params.data());

    std::vector<double> cov_errs(errs_.size(), 0.0);
    for (int it = 0; it < params.size(); ++it) {
        int idx = reduce_pars_.at(it);
        cov_errs.at(idx) = std::sqrt(covariance_params.at(it*params.size()+it));
        if (!std::isfinite(cov_errs.at(idx))) cov_errs.at(idx) = 0.0;
    }
    double ceres_solver_nchi = 2.0 * summary.final_cost / static_cast<double>(ndof_);
    if (!std::isfinite(ceres_solver_nchi) || ceres_solver_nchi <= 0.0) return false;

    // Errors from covariance
    for (int it = 0; it < errs_.size(); ++it) {
        errs_.at(it) = (cov_errs.at(it) > errs_.at(it)) ? cov_errs.at(it) : errs_.at(it);
        if (errs_.at(it) == 0.0) errs_.at(it) = 1.0;
    }

    status_ = true;
    return true;
}

bool DataFit1D::fine() {
    if (!status_) return false;

    return true;
}

bool VirtualDataFit1D::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const {
    std::fill_n(residuals, nbins_, 0.0);
    if (ntmps_ == 0) return false;
    
    Bool_t hasJacb = (jacobians != nullptr && jacobians[0] != nullptr);
    if (hasJacb) std::fill_n(jacobians[0], nbins_ * ntmps_, 0.0);
 
    std::vector<double> wgts(ntmps_, 0.0);
    for (int it = 0; it < ntmps_; ++it) wgts.at(it) = parameters[0][it];

    std::vector<double> resd(nbins_, 0.0);
    std::vector<double> jacb(nbins_ * ntmps_, 0.0);

    for (int ib = 0; ib < nbins_; ++ib) {
        Point1D data_smp = smp_.data(ib);
        
        Point1D data_tmps;
        for (int it = 0; it < ntmps_; ++it) {
            Point1D data_tmp = tmps_.at(it).weight_pdf(wgts.at(it), ib);
            data_tmps.val += data_tmp.val;
            data_tmps.err += (data_tmp.val == 0.0) ? 0.0 : (data_tmp.err * data_tmp.err);
        }
        data_tmps.err = std::sqrt(data_tmps.err);

        double error = std::hypot(data_smp.err, data_tmps.err);
        double value = (data_tmps.val - data_smp.val) / error;

        resd.at(ib) = value;
        for (int it = 0; it < ntmps_; ++it) {
            jacb.at(ib * ntmps_ + it) += (tmps_.at(it).pdf(ib).val / error);
        }
    }
    
    for (int ib = 0; ib < nbins_; ++ib) {
        if (!std::isfinite(resd.at(ib))) continue;
        residuals[ib] = resd.at(ib);
    }
    if (hasJacb) {
        for (int it = 0; it < nbins_ * ntmps_; ++it) { 
            if (!std::isfinite(jacb.at(it))) continue;
            jacobians[0][it] = jacb.at(it);
        }
    }

    return true;
}
        
HistFit1D::HistFit1D(const TH1D* smp, const std::vector<TH1D*>& tmps, const Axis1D& axis, const std::string& prefix_name, bool with_fluc, bool with_build_hist) {
    clear();
    if (!init(smp, tmps, axis, prefix_name, with_fluc, with_build_hist)) { clear(); return; }
    if (!fit()) { clear(); return; }
    if (status_ && with_build_hist_) build_hist();
}
        
HistFit1D::HistFit1D(const Hist1D& smp, const std::vector<Hist1D>& tmps, const Axis1D& axis, const std::string& prefix_name, bool with_fluc, bool with_build_hist) {
    clear();
    if (!init(smp, tmps, axis, prefix_name, with_fluc, with_build_hist)) { clear(); return; }
    if (!fit()) { clear(); return; }
    if (status_ && with_build_hist_) build_hist();
}
        
bool HistFit1D::init(const TH1D* smp, const std::vector<TH1D*>& tmps, const Axis1D& axis, const std::string& prefix_name, bool with_fluc, bool with_build_hist) {
    if (smp == nullptr || tmps.size() == 0) return false;
    for (const auto& tmp : tmps) if (tmp == nullptr) return false;
   
    hist_smp_ = Hist1D(smp->GetName(), smp->GetTitle(), smp, axis, true);
    hist_tmps_.clear();
    for (const auto& tmp : tmps) {
        hist_tmps_.push_back(
            Hist1D(tmp->GetName(), tmp->GetTitle(), tmp, axis, true)
        );
    }
    if (!hist_smp_.status()) return false;
    for (const auto& hist_tmp : hist_tmps_) if (!hist_tmp.status()) return false;

    with_fluc_ = with_fluc;
    with_build_hist_ = with_build_hist;
    prefix_name_ = prefix_name;
    ntmps_ = hist_tmps_.size();
    axis_ = axis;

    return true;
}
        
bool HistFit1D::init(const Hist1D& smp, const std::vector<Hist1D>& tmps, const Axis1D& axis, const std::string& prefix_name, bool with_fluc, bool with_build_hist) {
    if (!smp.status() || tmps.size() == 0) return false;
    for (const auto& tmp : tmps) if (!tmp.status()) return false;
    
    if (smp.axis().nbins() != axis.nbins()) return false;
    for (const auto& tmp : tmps) if (tmp.axis().nbins() != axis.nbins()) return false;

    hist_smp_ = Hist1D(smp.name(), smp.title(), smp.data(), smp.axis(), true);
    hist_tmps_.clear();
    for (const auto& tmp : tmps) {
        hist_tmps_.push_back(
            Hist1D(tmp.name(), tmp.title(), tmp.data(), tmp.axis(), true)
        );
    }

    if (!hist_smp_.status()) return false;
    for (const auto& hist_tmp : hist_tmps_) if (!hist_tmp.status()) return false;

    with_fluc_ = with_fluc;
    with_build_hist_ = with_build_hist;
    prefix_name_ = prefix_name;
    ntmps_ = hist_tmps_.size();
    axis_ = axis;

    return true;
}

void HistFit1D::clear() {
    with_fluc_ = false;
    with_build_hist_ = false;
    prefix_name_ = "FIT1D_";
    axis_ = Axis1D();

    hist_smp_ = Hist1D();
    hist_tmps_.clear();
    
    hist_ref_smp_  = Hist1D();
    hist_sum_tmps_ = Hist1D();
    hist_wgt_tmps_.clear();

    nsmp_ = 0;
    wgts_.clear();
    errs_.clear();
    fluc_.clear();

    fluc_datas_.clear();
    fluc_hists_.clear();
    fluc_nchi_.clear();
    
    status_ = false;
    ntmps_ = 0;
    nchi_ = 0;
    ndof_ = 0;
}
        
bool HistFit1D::fit() {
    Data1D data_smp = hist_smp_.data();
    std::vector<Data1D> data_tmps;
    for (const auto& hist_tmp : hist_tmps_) {
        data_tmps.push_back(hist_tmp.data());
    }

    DataFit1D fit1D(data_smp, data_tmps);
    if (!fit1D.status()) return false;
    
    hist_ref_smp_ = Hist1D(
        Form("%sSMP", prefix_name_.c_str()), "Sample",
        fit1D.smp(), axis_);
    //hist_ref_smp_.build_hist();
    
    hist_sum_tmps_ = Hist1D(
        Form("%sSUMTMPS", prefix_name_.c_str()), "Summation Template", 
        fit1D.sum_tmps(), axis_);
    //hist_sum_tmps_.build_hist();
   
    for (int it = 0; it < fit1D.ntmps(); ++it) {
        hist_wgt_tmps_.push_back(Hist1D(
            Form("%sWGTTMPS_%02d", prefix_name_.c_str(), it), Form("Template %2d", it), 
            fit1D.wgt_tmps(it), axis_));
        //hist_wgt_tmps_.back().build_hist();
    }
    
    nchi_ = fit1D.nchi();
    ndof_ = fit1D.ndof();
    
    nsmp_ = fit1D.nsmp();
    wgts_ = fit1D.wgts();
    errs_ = fit1D.errs();
    
    fluc_       = std::vector<Point1D>(ntmps_);
    fluc_datas_ = std::vector<std::vector<Point1D>>(ntmps_, std::vector<Point1D>());
    fluc_hists_ = std::vector<Hist1D>(ntmps_);
    fluc_nchi_  = std::vector<double>(ntmps_, 0.0);
    //std::cerr << "LINE  " << __LINE__ << std::endl;
    if (with_fluc_ && !flucfit()) return false;
    //std::cerr << "LINE  " << __LINE__ << std::endl;

    status_ = true;
    return true;
}

bool HistFit1D::flucfit() {
    const int MAXITER = 30000;
    const int NSAMPLE = 10000;
    
    Data1D data_smp = hist_smp_.data();
    std::vector<Data1D> data_tmps;
    for (const auto& hist_tmp : hist_tmps_) {
        data_tmps.push_back(hist_tmp.data());
    }
    
    int iter = 0;
    int isample = 0;
    std::vector<std::vector<Point1D>> fluc_datas(ntmps_, std::vector<Point1D>());
    while (isample < NSAMPLE && iter < MAXITER) {
        std::vector<Data1D> fluc_data_tmps;
        for (int it = 0; it < ntmps_; ++it) {
            fluc_data_tmps.push_back(data_tmps.at(it).build_fluc());
        }

        DataFit1D fluc_fit1D(data_smp, fluc_data_tmps);
        if (!fluc_fit1D.status()) { iter++; continue; }
        for (int it = 0; it < ntmps_; ++it) {
            Point1D pnt(fluc_fit1D.wgts(it), fluc_fit1D.errs(it));
            fluc_datas.at(it).push_back(pnt);
        }
        
        isample++;
        iter++;
    }
    if (isample < static_cast<int>(0.3 * NSAMPLE)) return false;

    std::vector<double> fluc_nchi;
    std::vector<Point1D> fluc_params;
    for (int itmp = 0; itmp < ntmps_; ++itmp) {
        if (wgts_.at(itmp) == 0.0 || errs_.at(itmp) == 0.0) {
            fluc_params.push_back(Point1D(0.0, 0.0));
            continue;
        }

        double sum_v1 = 0.0;
        for (const auto& pnt : fluc_datas.at(itmp)) sum_v1 += pnt.val;
        double mean = sum_v1 / static_cast<double>(fluc_datas.at(itmp).size());
        if (!std::isfinite(mean)) mean = 0.0;
        
        double sum_v2 = 0.0;
        for (const auto& pnt : fluc_datas.at(itmp)) sum_v2 += (pnt.val - mean) * (pnt.val - mean);
        double sigma = std::sqrt(sum_v2 / static_cast<double>(fluc_datas.at(itmp).size()));
        if (!std::isfinite(sigma)) sigma = 0.0;

        std::vector<double> params({ static_cast<double>(fluc_datas.at(itmp).size()), mean, sigma });

        if (params.at(0) == 0.0 || params.at(1) == 0.0 || params.at(2) == 0.0) {
            fluc_params.push_back(Point1D(0.0, 0.0));
            continue;
        }

        // Build histogram
        constexpr int axis_half_nbin = 70;
        std::vector<double> axis_bins(2*axis_half_nbin+1, 0.0);
        double axis_bin_width = 7.0 * params.at(2) / static_cast<double>(axis_half_nbin);
        for (int ib = -axis_half_nbin; ib <= axis_half_nbin; ++ib) {
            axis_bins.at(axis_half_nbin + ib) = params.at(1) + ib * axis_bin_width;
        }

        std::vector<int> content_bins(2*axis_half_nbin+2, 0);
        for (const auto& pnt : fluc_datas.at(itmp)) {
            int bin = 1 + static_cast<int>(std::floor((pnt.val - axis_bins.at(0)) / axis_bin_width));
            if (bin < 0) bin = 0;
            if (bin > 2*axis_half_nbin+1) bin = 2*axis_half_nbin+1;
            content_bins.at(bin)++;
        }

        std::vector<std::array<Point1D, 2>> hist_datas;
        for (int ib = 1; ib <= 2*axis_half_nbin; ++ib) {
            double x = 0.5 * (axis_bins.at(ib-1) + axis_bins.at(ib));
            double y = static_cast<double>(content_bins.at(ib));
            double ex = 0.5 * axis_bin_width;
            double ey = std::sqrt(y);
            if (!std::isfinite(ey)) continue;
            if (y == 0.0 || ey == 0.0) continue;
            std::array<Point1D, 2> point({ Point1D(x, ex), Point1D(y, ey) });
            hist_datas.push_back(point);
        }

        // CeresSolver: Cost Function
        ceres::CostFunction* cost_function = VirtualFlucFit1D::Create(hist_datas, params);
        if (cost_function == nullptr) {
            fluc_params.push_back(Point1D(0.0, 0.0));
            continue;
        }
        
        // CeresSolver: Loss Function
        //ceres::LossFunction* loss_function = new ceres::HuberLoss(3.0);

        // CeresSolver: Problem
        ceres::Problem problem;
        //problem.AddResidualBlock(cost_function, loss_function, params.data());
        problem.AddResidualBlock(cost_function, nullptr, params.data());
        problem.SetParameterLowerBound(params.data(), 0, 0.0);
        problem.SetParameterLowerBound(params.data(), 2, 0.0);
            
        // CeresSolver: Options
        ceres::Solver::Options options;
        options.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT;
        options.max_num_iterations = 100;

        // CeresSolver: Summary
        ceres::Solver::Summary summary;
        ceres::Solve(options, &problem, &summary);
        if (!summary.IsSolutionUsable()) continue;
        if (ceres::NO_CONVERGENCE == summary.termination_type) continue;
        //std::cerr << summary.FullReport() << std::endl;
        
        double nchi = 2.0 * summary.final_cost / static_cast<double>(hist_datas.size() - params.size());
        if (!std::isfinite(nchi)) continue;
        fluc_nchi.push_back(nchi);

        Point1D fluc_param(params.at(1), params.at(2));
        fluc_params.push_back(fluc_param);
    }
    if (fluc_params.size() != ntmps_) return false;
    
    std::vector<Hist1D> fluc_hists;
    for (int itmp = 0; itmp < ntmps_; ++itmp) {
        const std::vector<Point1D>& data = fluc_datas.at(itmp);
        const Point1D& param = fluc_params.at(itmp);

        constexpr int axis_half_nbin = 70;
        std::vector<double> axis_bins(2*axis_half_nbin+1, 0.0);
        double axis_bin_width = 7.0 * param.err / static_cast<double>(axis_half_nbin);
        for (int ib = -axis_half_nbin; ib <= axis_half_nbin; ++ib) {
            axis_bins.at(axis_half_nbin + ib) = param.val + ib * axis_bin_width;
        }
        Axis1D axis1D(Form("Number [Temp %2d]", itmp), "Events/Bin", axis_bins); 

        std::vector<int> content_bins(2*axis_half_nbin+2, 0);
        for (const auto& pnt : data) {
            int bin = 1 + static_cast<int>(std::floor((pnt.val - axis_bins.at(0)) / axis_bin_width));
            if (bin < 0) bin = 0;
            if (bin > 2*axis_half_nbin+1) bin = 2*axis_half_nbin+1;
            content_bins.at(bin)++;
        }

        Data1D data1D(Form("%sFLUCS_%02d", prefix_name_.c_str(), itmp), Form("Fluc Template %2d", itmp));
        for (int ib = 1; ib <= 2*axis_half_nbin; ++ib) {
            data1D.push_back(content_bins.at(ib), std::sqrt(static_cast<double>(content_bins.at(ib))));
        }
        data1D.build();

        Hist1D hist1D(data1D.name(), data1D.title(), data1D, axis1D);
        //hist1D.build_hist();

        fluc_hists.push_back(hist1D);
    }
    
    fluc_       = fluc_params;
    fluc_datas_ = fluc_datas;
    fluc_hists_ = fluc_hists;
    fluc_nchi_  = fluc_nchi;

    return true;
}

} // namespace HistFit

#endif // __HistFit1D_C__
