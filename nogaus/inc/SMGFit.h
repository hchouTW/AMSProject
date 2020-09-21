class SMGFit {
    protected :
        std::vector<double> hits_;
        TrSys::MultiGaus    smg_;

        bool   status_;
        double param_;
        double error_;
        double nchi_;

    public :
        SMGFit(const std::vector<double>& hits, const TrSys::MultiGaus& smg);

        const bool& status() const { return status_; }

        const double& param() const { return param_; }
        const double& error() const { return error_; }
        
        const double& nchi() const { return nchi_; }
};

class VirtualSMGFit : public ceres::CostFunction {
    protected :
        const std::vector<double> hits_;
        const TrSys::MultiGaus    smg_;

        int nres_;
        int npar_;

    public :
        VirtualSMGFit(const std::vector<double>& hits, const TrSys::MultiGaus& smg) : hits_(hits), smg_(smg), nres_(0), npar_(0) {
            nres_ = hits_.size();
            npar_ = 1;

            set_num_residuals(nres_);
            mutable_parameter_block_sizes()->clear(); 
            mutable_parameter_block_sizes()->push_back(npar_);
        }
    
    public :
        virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;
};
    

SMGFit::SMGFit(const std::vector<double>& hits, const TrSys::MultiGaus& smg) : hits_(hits), smg_(smg) {
    // init params
    status_ = false;
    param_ = 0.0;
    error_ = 0.0;
    nchi_ = 0.0;

    // fit parameters
    std::vector<double> params({ std::accumulate(hits_.begin(), hits_.end(), 0.0) / hits_.size() });

    // CeresSolver: Cost Function
    ceres::CostFunction* cost_function = new VirtualSMGFit(hits_, smg_);

    // CeresSolver: Problem
    ceres::Problem problem;
    problem.AddResidualBlock(cost_function, nullptr, params.data());
    //problem.SetParameterLowerBound(params.data(), 0, -10.0);
    //problem.SetParameterUpperBound(params.data(), 0,  10.0);
   
    // CeresSolver: Options
    ceres::Solver::Options options;
    options.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT;
    options.max_num_iterations = 100;

    // CeresSolver: Summary
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    if (!summary.IsSolutionUsable()) return;
    if (ceres::NO_CONVERGENCE == summary.termination_type) return;
    //std::cerr << summary.FullReport() << std::endl;
    
    ceres::Covariance::Options covariance_options;
    ceres::Covariance covariance(covariance_options);

    std::vector<std::pair<const double*, const double*> > covariance_blocks;
    covariance_blocks.push_back(std::make_pair(params.data(), params.data()));

    CHECK(covariance.Compute(covariance_blocks, &problem));

    std::vector<double> covariance_params(params.size() * params.size(), 0.0);
    covariance.GetCovarianceBlock(params.data(), params.data(), covariance_params.data());

    std::vector<double> errors(params.size(), 0.0);
    for (int it = 0; it < params.size(); ++it) {
        errors.at(it) = std::sqrt(covariance_params.at(it*params.size()+it));
        if (!std::isfinite(errors.at(it))) errors.at(it) = 0.0;
    }
    double nchi = 2.0 * summary.final_cost / static_cast<double>(hits_.size()-1);
    if (!std::isfinite(nchi) || nchi < 0.0) return;

    // Result
    status_ = true;
    param_ = params[0];
    error_ = errors[0];
    nchi_ = nchi;
}


bool VirtualSMGFit::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const {
    std::fill_n(residuals, nres_, 0.0);

    bool hasJacb = (jacobians != nullptr && jacobians[0] != nullptr);
    if (hasJacb) std::fill_n(jacobians[0], nres_ * npar_, 0.0);
   
    double param = parameters[0][0];

    TrSys::Vec<double>&& rs = TrSys::Vec<double>::Zero(nres_);        // Residual
    TrSys::Mat<double>&& jb = TrSys::Mat<double>::Zero(nres_, npar_); // Jacb

    int cnt_nhit = 0;
    for (auto&& hit : hits_) {
        std::array<double, 3>&& minimizer = smg_.minimizer(hit - param);

        rs(cnt_nhit) += minimizer[0];
        if (hasJacb) jb(cnt_nhit, 0) += minimizer[2];

        cnt_nhit++;
    }
    if (cnt_nhit != hits_.size()) return false;
    
    for (int it = 0; it < nres_; ++it) {
        if (!std::isfinite(rs(it))) rs(it) = 0.0;
        residuals[it] = rs(it);
    }
    
    if (hasJacb) {
        for (int it = 0; it < nres_; ++it) {
        for (int jt = 0; jt < npar_; ++jt) {
            if (!std::isfinite(jb(it, jt))) jb(it, jt) = 0.0;
            jacobians[0][it * npar_ + jt] = jb(it, jt); 
        }}
    }

    return true;
}
