#ifndef __DataFit2D_C__
#define __DataFit2D_C__

#include "DataFit2D.h"


ResultFit2D::ResultFit2D(
    const std::string& name, int ndimx, int ndimy, int ntmp,
    const std::vector<std::vector<std::array<double, 2>>>& tmps, 
    double num_smp, const std::vector<std::array<double, 2>>& smp, 
    double num_ref, const std::vector<std::array<double, 2>>& ref) {
    clear();
    data_name_ = name;
    data_ndimx_ = ndimx;
    data_ndimy_ = ndimy;
    data_ntmp_ = ntmp;
    data_tmps_ = tmps;
    
    num_smp_  = num_smp;
    data_smp_ = smp;

    if (num_ref > 0.0 && ref.size() == data_ndimy_*data_ndimx_) {
        num_ref_  = num_ref;
        data_ref_ = ref;
    }
}

void ResultFit2D::clear() {
    hist_build_ = false;
    hist_xaxis_ = TAxis();
    hist_yaxis_ = TAxis();

    data_name_ = "";
    data_ndimx_ = 0;
    data_ndimy_ = 0;
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
    
    tmps_x_.clear();
    sum_x_.clear();
    smp_x_.clear();
    ref_x_.clear();
    
    tmps_y_.clear();
    sum_y_.clear();
    smp_y_.clear();
    ref_y_.clear();
}

void ResultFit2D::set_fitting_result(int ndof, double nchi, const std::vector<double>& wgt_tmps, const std::vector<double>& wgt_errs) {
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
        tmps_.push_back(std::vector<std::array<double, 2>>(data_ndimy_*data_ndimx_, { 0.0, 0.0 }));
        for (int jt = 0; jt < data_ndimy_; ++jt) {
        for (int it = 0; it < data_ndimx_; ++it) {
            int idx = jt * data_ndimx_ + it;
            tmps_.back()[idx][0] = num_tmps_.at(itmp) * data_tmps_.at(itmp).at(idx)[0];
            tmps_.back()[idx][1] = num_tmps_.at(itmp) * data_tmps_.at(itmp).at(idx)[1];
            if (tmps_.back().at(idx)[0] <= 0.0) tmps_.back().at(idx)[1] = 0.0;
        }}
    }

    sum_ = std::vector<std::array<double, 2>>(data_ndimy_*data_ndimx_, { 0.0, 0.0 });
    for (int jt = 0; jt < data_ndimy_; ++jt) {
    for (int it = 0; it < data_ndimx_; ++it) {
        int idx = jt * data_ndimx_ + it;
        double sum_val = 0.0;
        double sum_err = 0.0;
        for (int itmp = 0; itmp < data_ntmp_; ++itmp) {
            sum_val += tmps_.at(itmp).at(idx)[0];
            sum_err += tmps_.at(itmp).at(idx)[1] * tmps_.at(itmp).at(idx)[1];
        }
        sum_err = std::sqrt(sum_err);
        sum_.at(idx)[0] = sum_val;
        sum_.at(idx)[1] = sum_err;
        if (sum_.at(idx)[0] <= 0.0) sum_.at(idx)[1] = 0.0;
    }}
    
    smp_ = std::vector<std::array<double, 2>>(data_ndimy_*data_ndimx_, { 0.0, 0.0 });
    for (int jt = 0; jt < data_ndimy_; ++jt) {
    for (int it = 0; it < data_ndimx_; ++it) {
        int idx = jt * data_ndimx_ + it;
        smp_.at(idx)[0] = num_smp_ * data_smp_.at(idx)[0];
        smp_.at(idx)[1] = num_smp_ * data_smp_.at(idx)[1];
        if (smp_.at(idx)[0] <= 0.0) smp_.at(idx)[1] = 0.0;
    }}

    if (num_ref_ > 0.0) {
        ref_ = std::vector<std::array<double, 2>>(data_ndimy_*data_ndimx_, { 0.0, 0.0 });
        for (int jt = 0; jt < data_ndimy_; ++jt) {
        for (int it = 0; it < data_ndimx_; ++it) {
            int idx = jt * data_ndimx_ + it;
            ref_.at(idx)[0] = num_ref_ * data_ref_.at(idx)[0];
            ref_.at(idx)[1] = num_ref_ * data_ref_.at(idx)[1];
            if (ref_.at(idx)[0] <= 0.0) ref_.at(idx)[1] = 0.0;
        }}
    }
    
    tmps_x_.clear();
    for (int itmp = 0; itmp < data_ntmp_; ++itmp) {
        tmps_x_.push_back(std::vector<std::array<double, 2>>(data_ndimx_, { 0.0, 0.0 }));
        for (int it = 0; it < data_ndimx_; ++it) {
            double sum_val = 0.0;
            double sum_err = 0.0;
            for (int jt = 0; jt < data_ndimy_; ++jt) {
                int idx = jt * data_ndimx_ + it;
                sum_val += tmps_.at(itmp).at(idx)[0];
                sum_err += tmps_.at(itmp).at(idx)[1] * tmps_.at(itmp).at(idx)[1];
            }
            tmps_x_.back().at(it)[0] = sum_val;
            tmps_x_.back().at(it)[1] = std::sqrt(sum_err);
        }
    }
    
    sum_x_ = std::vector<std::array<double, 2>>(data_ndimx_, { 0.0, 0.0 });
    for (int it = 0; it < data_ndimx_; ++it) {
        double sum_val = 0.0;
        double sum_err = 0.0;
        for (int jt = 0; jt < data_ndimy_; ++jt) {
            int idx = jt * data_ndimx_ + it;
            sum_val += sum_.at(idx)[0];
            sum_err += sum_.at(idx)[1] * sum_.at(idx)[1];
        }
        sum_x_.at(it)[0] = sum_val;
        sum_x_.at(it)[1] = std::sqrt(sum_err);
    }

    smp_x_ = std::vector<std::array<double, 2>>(data_ndimx_, { 0.0, 0.0 });
    for (int it = 0; it < data_ndimx_; ++it) {
        double sum_val = 0.0;
        double sum_err = 0.0;
        for (int jt = 0; jt < data_ndimy_; ++jt) {
            int idx = jt * data_ndimx_ + it;
            sum_val += smp_.at(idx)[0];
            sum_err += smp_.at(idx)[1] * smp_.at(idx)[1];
        }
        smp_x_.at(it)[0] = sum_val;
        smp_x_.at(it)[1] = std::sqrt(sum_err);
    }
    
    if (num_ref_ > 0.0) {
        ref_x_ = std::vector<std::array<double, 2>>(data_ndimx_, { 0.0, 0.0 });
        for (int it = 0; it < data_ndimx_; ++it) {
            double sum_val = 0.0;
            double sum_err = 0.0;
            for (int jt = 0; jt < data_ndimy_; ++jt) {
                int idx = jt * data_ndimx_ + it;
                sum_val += ref_.at(idx)[0];
                sum_err += ref_.at(idx)[1] * ref_.at(idx)[1];
            }
            ref_x_.at(it)[0] = sum_val;
            ref_x_.at(it)[1] = std::sqrt(sum_err);
        }
    }
    
    tmps_y_.clear();
    for (int itmp = 0; itmp < data_ntmp_; ++itmp) {
        tmps_y_.push_back(std::vector<std::array<double, 2>>(data_ndimy_, { 0.0, 0.0 }));
        for (int jt = 0; jt < data_ndimy_; ++jt) {
            double sum_val = 0.0;
            double sum_err = 0.0;
            for (int it = 0; it < data_ndimx_; ++it) {
                int idx = jt * data_ndimx_ + it;
                sum_val += tmps_.at(itmp).at(idx)[0];
                sum_err += tmps_.at(itmp).at(idx)[1] * tmps_.at(itmp).at(idx)[1];
            }
            tmps_y_.back().at(jt)[0] = sum_val;
            tmps_y_.back().at(jt)[1] = std::sqrt(sum_err);
        }
    }
    
    sum_y_ = std::vector<std::array<double, 2>>(data_ndimy_, { 0.0, 0.0 });
    for (int jt = 0; jt < data_ndimy_; ++jt) {
        double sum_val = 0.0;
        double sum_err = 0.0;
        for (int it = 0; it < data_ndimx_; ++it) {
            int idx = jt * data_ndimx_ + it;
            sum_val += sum_.at(idx)[0];
            sum_err += sum_.at(idx)[1] * sum_.at(idx)[1];
        }
        sum_y_.at(jt)[0] = sum_val;
        sum_y_.at(jt)[1] = std::sqrt(sum_err);
    }

    smp_y_ = std::vector<std::array<double, 2>>(data_ndimy_, { 0.0, 0.0 });
        for (int jt = 0; jt < data_ndimy_; ++jt) {
        double sum_val = 0.0;
        double sum_err = 0.0;
            for (int it = 0; it < data_ndimx_; ++it) {
            int idx = jt * data_ndimx_ + it;
            sum_val += smp_.at(idx)[0];
            sum_err += smp_.at(idx)[1] * smp_.at(idx)[1];
        }
        smp_y_.at(jt)[0] = sum_val;
        smp_y_.at(jt)[1] = std::sqrt(sum_err);
    }
    
    if (num_ref_ > 0.0) {
        ref_y_ = std::vector<std::array<double, 2>>(data_ndimy_, { 0.0, 0.0 });
        for (int jt = 0; jt < data_ndimy_; ++jt) {
            double sum_val = 0.0;
            double sum_err = 0.0;
            for (int it = 0; it < data_ndimx_; ++it) {
                int idx = jt * data_ndimx_ + it;
                sum_val += ref_.at(idx)[0];
                sum_err += ref_.at(idx)[1] * ref_.at(idx)[1];
            }
            ref_y_.at(jt)[0] = sum_val;
            ref_y_.at(jt)[1] = std::sqrt(sum_err);
        }
    }
}
        

void ResultFit2D::set_histogram_result(bool hist_build, TAxis* hist_xaxis, TAxis* hist_yaxis, const std::string& hist_xtitle, const std::string& hist_ytitle, const std::string& hist_ztitle) {
    hist_tmps_.clear();
    if (hist_build && hist_xaxis != nullptr && hist_yaxis != nullptr) {
        hist_build_ = hist_build;
        hist_xaxis_ = (*hist_xaxis);
        hist_yaxis_ = (*hist_yaxis);
    }
    if (!hist_build_) return;
    
    for (int itmp = 0; itmp < data_ntmp_; ++itmp) {
        std::shared_ptr<TH2D> hist_tmp = std::make_shared<TH2D>(Form("hFit2D_tmp%d_%s", itmp, data_name_.c_str()), "", hist_xaxis_.GetNbins(), hist_xaxis_.GetXbins()->GetArray(), hist_yaxis_.GetNbins(), hist_yaxis_.GetXbins()->GetArray());
        for (int jt = 0; jt < data_ndimy_; ++jt) {
        for (int it = 0; it < data_ndimx_; ++it) {
            int idx = jt * data_ndimx_ + it;
            hist_tmp->SetBinContent(it+1, jt+1, tmps_.at(itmp).at(idx)[0]);
            hist_tmp->SetBinError  (it+1, jt+1, tmps_.at(itmp).at(idx)[1]);
        }}
        hist_tmp->GetXaxis()->SetTitle(hist_xtitle.c_str());
        hist_tmp->GetYaxis()->SetTitle(hist_ytitle.c_str());
        hist_tmp->GetZaxis()->SetTitle(hist_ztitle.c_str());
        hist_tmps_.push_back(hist_tmp);
    }
    
    hist_sum_ = std::make_shared<TH2D>(Form("hFit2D_sum_%s", data_name_.c_str()), "", hist_xaxis_.GetNbins(), hist_xaxis_.GetXbins()->GetArray(), hist_yaxis_.GetNbins(), hist_yaxis_.GetXbins()->GetArray());
    for (int jt = 0; jt < data_ndimy_; ++jt) {
    for (int it = 0; it < data_ndimx_; ++it) {
        int idx = jt * data_ndimx_ + it;
        hist_sum_->SetBinContent(it+1, jt+1, sum_.at(idx)[0]);
        hist_sum_->SetBinError  (it+1, jt+1, sum_.at(idx)[1]);
    }}
    hist_sum_->GetXaxis()->SetTitle(hist_xtitle.c_str());
    hist_sum_->GetYaxis()->SetTitle(hist_ytitle.c_str());
    hist_sum_->GetZaxis()->SetTitle(hist_ztitle.c_str());
    
    hist_smp_ = std::make_shared<TH2D>(Form("hFit2D_smp_%s", data_name_.c_str()), "", hist_xaxis_.GetNbins(), hist_xaxis_.GetXbins()->GetArray(), hist_yaxis_.GetNbins(), hist_yaxis_.GetXbins()->GetArray());
    for (int jt = 0; jt < data_ndimy_; ++jt) {
    for (int it = 0; it < data_ndimx_; ++it) {
        int idx = jt * data_ndimx_ + it;
        hist_smp_->SetBinContent(it+1, jt+1, smp_.at(idx)[0]);
        hist_smp_->SetBinError  (it+1, jt+1, smp_.at(idx)[1]);
    }}
    hist_smp_->GetXaxis()->SetTitle(hist_xtitle.c_str());
    hist_smp_->GetYaxis()->SetTitle(hist_ytitle.c_str());
    hist_smp_->GetZaxis()->SetTitle(hist_ztitle.c_str());
    
    if (num_ref_ > 0.0) {
        hist_ref_ = std::make_shared<TH2D>(Form("hFit2D_ref_%s", data_name_.c_str()), "", hist_xaxis_.GetNbins(), hist_xaxis_.GetXbins()->GetArray(), hist_yaxis_.GetNbins(), hist_yaxis_.GetXbins()->GetArray());
        for (int jt = 0; jt < data_ndimy_; ++jt) {
        for (int it = 0; it < data_ndimx_; ++it) {
            int idx = jt * data_ndimx_ + it;
            hist_ref_->SetBinContent(it+1, jt+1, ref_.at(idx)[0]);
            hist_ref_->SetBinError  (it+1, jt+1, ref_.at(idx)[1]);
        }}
        hist_ref_->GetXaxis()->SetTitle(hist_xtitle.c_str());
        hist_ref_->GetYaxis()->SetTitle(hist_ytitle.c_str());
        hist_ref_->GetZaxis()->SetTitle(hist_ztitle.c_str());
    }
    
    for (int itmp = 0; itmp < data_ntmp_; ++itmp) {
        std::shared_ptr<TH1D> hist_tmp_x = std::make_shared<TH1D>(Form("hFit2D_tmp%d_x_%s", itmp, data_name_.c_str()), "", hist_xaxis_.GetNbins(), hist_xaxis_.GetXbins()->GetArray());
        for (int it = 0; it < data_ndimx_; ++it) {
            hist_tmp_x->SetBinContent(it+1, tmps_x_.at(itmp).at(it)[0]);
            hist_tmp_x->SetBinError  (it+1, tmps_x_.at(itmp).at(it)[1]);
        }
        hist_tmp_x->GetXaxis()->SetTitle(hist_xtitle.c_str());
        hist_tmp_x->GetYaxis()->SetTitle(hist_ztitle.c_str());
        hist_tmps_x_.push_back(hist_tmp_x);
    }
    
    hist_sum_x_ = std::make_shared<TH1D>(Form("hFit2D_sum_x_%s", data_name_.c_str()), "", hist_xaxis_.GetNbins(), hist_xaxis_.GetXbins()->GetArray());
    for (int it = 0; it < data_ndimx_; ++it) {
        hist_sum_x_->SetBinContent(it+1, sum_x_.at(it)[0]);
        hist_sum_x_->SetBinError  (it+1, sum_x_.at(it)[1]);
    }
    hist_sum_x_->GetXaxis()->SetTitle(hist_xtitle.c_str());
    hist_sum_x_->GetYaxis()->SetTitle(hist_ztitle.c_str());
    
    hist_smp_x_ = std::make_shared<TH1D>(Form("hFit2D_smp_x_%s", data_name_.c_str()), "", hist_xaxis_.GetNbins(), hist_xaxis_.GetXbins()->GetArray());
    for (int it = 0; it < data_ndimx_; ++it) {
        hist_smp_x_->SetBinContent(it+1, smp_x_.at(it)[0]);
        hist_smp_x_->SetBinError  (it+1, smp_x_.at(it)[1]);
    }
    hist_smp_x_->GetXaxis()->SetTitle(hist_xtitle.c_str());
    hist_smp_x_->GetYaxis()->SetTitle(hist_ztitle.c_str());
    
    if (num_ref_ > 0.0) {
        hist_ref_x_ = std::make_shared<TH1D>(Form("hFit2D_ref_x_%s", data_name_.c_str()), "", hist_xaxis_.GetNbins(), hist_xaxis_.GetXbins()->GetArray());
        for (int it = 0; it < data_ndimx_; ++it) {
            hist_ref_x_->SetBinContent(it+1, ref_x_.at(it)[0]);
            hist_ref_x_->SetBinError  (it+1, ref_x_.at(it)[1]);
        }
        hist_ref_x_->GetXaxis()->SetTitle(hist_xtitle.c_str());
        hist_ref_x_->GetYaxis()->SetTitle(hist_ztitle.c_str());
    }

    for (int itmp = 0; itmp < data_ntmp_; ++itmp) {
        std::shared_ptr<TH1D> hist_tmp_y = std::make_shared<TH1D>(Form("hFit2D_tmp%d_y_%s", itmp, data_name_.c_str()), "", hist_yaxis_.GetNbins(), hist_yaxis_.GetXbins()->GetArray());
        for (int jt = 0; jt < data_ndimy_; ++jt) {
            hist_tmp_y->SetBinContent(jt+1, tmps_y_.at(itmp).at(jt)[0]);
            hist_tmp_y->SetBinError  (jt+1, tmps_y_.at(itmp).at(jt)[1]);
        }
        hist_tmp_y->GetXaxis()->SetTitle(hist_ytitle.c_str());
        hist_tmp_y->GetYaxis()->SetTitle(hist_ztitle.c_str());
        hist_tmps_y_.push_back(hist_tmp_y);
    }
    
    hist_sum_y_ = std::make_shared<TH1D>(Form("hFit2D_sum_y_%s", data_name_.c_str()), "", hist_yaxis_.GetNbins(), hist_yaxis_.GetXbins()->GetArray());
    for (int jt = 0; jt < data_ndimy_; ++jt) {
        hist_sum_y_->SetBinContent(jt+1, sum_y_.at(jt)[0]);
        hist_sum_y_->SetBinError  (jt+1, sum_y_.at(jt)[1]);
    }
    hist_sum_y_->GetXaxis()->SetTitle(hist_ytitle.c_str());
    hist_sum_y_->GetYaxis()->SetTitle(hist_ztitle.c_str());
    
    hist_smp_y_ = std::make_shared<TH1D>(Form("hFit2D_smp_y_%s", data_name_.c_str()), "", hist_yaxis_.GetNbins(), hist_yaxis_.GetXbins()->GetArray());
    for (int jt = 0; jt < data_ndimy_; ++jt) {
        hist_smp_y_->SetBinContent(jt+1, smp_y_.at(jt)[0]);
        hist_smp_y_->SetBinError  (jt+1, smp_y_.at(jt)[1]);
    }
    hist_smp_y_->GetXaxis()->SetTitle(hist_ytitle.c_str());
    hist_smp_y_->GetYaxis()->SetTitle(hist_ztitle.c_str());
    
    if (num_ref_ > 0.0) {
        hist_ref_y_ = std::make_shared<TH1D>(Form("hFit2D_ref_y_%s", data_name_.c_str()), "", hist_yaxis_.GetNbins(), hist_yaxis_.GetXbins()->GetArray());
        for (int jt = 0; jt < data_ndimy_; ++jt) {
            hist_ref_y_->SetBinContent(jt+1, ref_y_.at(jt)[0]);
            hist_ref_y_->SetBinError  (jt+1, ref_y_.at(jt)[1]);
        }
        hist_ref_y_->GetXaxis()->SetTitle(hist_ytitle.c_str());
        hist_ref_y_->GetYaxis()->SetTitle(hist_ztitle.c_str());
    }
}

DataFit2D::DataFit2D(const std::string& name, std::vector<TH2D*> htmps, TH2D* hsmp, TH2D* href, bool hist_build, const std::string& hist_xtitle, const std::string& hist_ytitle, const std::string& hist_ztitle) {
    clear();

    if (htmps.size() == 0) return;
    for (auto&& htmp : htmps) { if (htmp == nullptr) return; }
    if (hsmp == nullptr) return;

    hist_xtitle_ = hist_xtitle;
    hist_ytitle_ = hist_ytitle;
    hist_ztitle_ = hist_ztitle;

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
        hist_xaxis_ = (*(hsmp_->GetXaxis()));
        hist_yaxis_ = (*(hsmp_->GetYaxis()));
    }
    
    data_name_ = name;
    data_ndimx_ = hsmp_->GetXaxis()->GetNbins();
    data_ndimy_ = hsmp_->GetYaxis()->GetNbins();
    data_ntmp_ = data_tmps_.size();

    if (data_ndimx_ == 0) { clear(); return; }
    if (data_ndimy_ == 0) { clear(); return; }
    if (data_ntmp_ == 0) { clear(); return; }

    // Check
    int count_not_empty_tmps = 0;
    for (int it = 0; it < data_ntmp_; ++it) {
        if (data_num_tmps_.at(it) <= 0.0) { clear(); return; } 
        if (data_tmps_.at(it).size() != data_ndimy_*data_ndimx_) { clear(); return; }
        count_not_empty_tmps++;
    }
    if (count_not_empty_tmps == 0) { clear(); return; }
    
    if (data_num_smp_ <= 0.0) { clear(); return; }
    if (data_smp_.size() != data_ndimy_*data_ndimx_) { clear(); return; }
    
    if (href_ != nullptr && data_num_ref_ <= 0.0) { clear(); return; }
    if (href_ != nullptr && data_ref_.size() != data_ndimy_*data_ndimx_) { clear(); return; }

    if (!fit()) clear();
}

void DataFit2D::clear() {
    hist_build_ = false;
    hist_xaxis_ = TAxis();
    hist_yaxis_ = TAxis();
    hist_xtitle_ = "";
    hist_ytitle_ = "";
    hist_ztitle_ = "";

    data_name_ = "";
    data_ndimx_ = 0;
    data_ndimy_ = 0;
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

    result_ = ResultFit2D();
}

std::tuple<double, std::vector<std::array<double, 2>>> DataFit2D::build_data(TH2D* hist) {
    std::tuple<double, std::vector<std::array<double, 2>>> datas;
    std::get<0>(datas) = 0.0;

    if (hist == nullptr) return datas;
    std::vector<std::array<double, 2>> data(hist->GetYaxis()->GetNbins()*hist->GetXaxis()->GetNbins(), {0.0, 0.0});

    double sum_val = 0.0;
    for (int iby = 1; iby <= hist->GetYaxis()->GetNbins(); ++iby) {
    for (int ibx = 1; ibx <= hist->GetXaxis()->GetNbins(); ++ibx) {
        if (hist->GetBinContent(ibx, iby) <= 0.0) continue;
        int idx = (iby-1) * hist->GetXaxis()->GetNbins() + (ibx-1);
        data.at(idx)[0] = hist->GetBinContent(ibx, iby);
        data.at(idx)[1] = hist->GetBinError  (ibx, iby);
        sum_val += data.at(idx)[0];
    }}
   
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

bool DataFit2D::fit() {
    std::vector<double> params(data_ntmp_, 1.0/static_cast<double>(data_ntmp_));

    // CeresSolver: Cost Function
    ceres::CostFunction* cost_function = new VirtualDataFit2D(data_ndimx_, data_ndimy_, data_ntmp_, data_tmps_, data_smp_);

    // CeresSolver: Loss Function
    ceres::LossFunction* loss_function = new ceres::HuberLoss(3.0);

    // CeresSolver: Problem
    ceres::Problem problem;
    problem.AddResidualBlock(cost_function, loss_function, params.data());

    for (int it = 0; it < params.size(); ++it) {
        problem.SetParameterLowerBound(params.data(), it, 0.0);
        problem.SetParameterUpperBound(params.data(), it, 1.5);
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
   
    for (auto&& val : params) {
        if (!std::isfinite(val)) return false;
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
        if (!std::isfinite(errors[it])) return false;
    }

    int    ndof = static_cast<VirtualDataFit2D*>(cost_function)->ndof();
    double nchi = 2.0 * summary.final_cost / static_cast<double>(ndof);

    ResultFit2D result(data_name_, data_ndimx_, data_ndimy_, data_ntmp_, data_tmps_, data_num_smp_, data_smp_, data_num_ref_, data_ref_);
    result.set_fitting_result(ndof, nchi, params, errors);
    result.set_histogram_result(hist_build_, &hist_xaxis_, &hist_yaxis_, hist_xtitle_, hist_ytitle_, hist_ztitle_);

    result_ = result;
    return true;
}
        
bool VirtualDataFit2D::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const {
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
       
        double sum_val = smp_[it][0];
        double sum_err = smp_[it][1] * smp_[it][1];
        for (int itmp = 0; itmp < ntmp_; ++itmp) {
            if (!has_tmps_[itmp][it]) continue;
            sum_val -= wgt_tmps[itmp] * tmps_.at(itmp).at(it)[0];
            //sum_err += (wgt_tmps[itmp] * tmps_.at(itmp).at(it)[1]) * (wgt_tmps[itmp] * tmps_.at(itmp).at(it)[1]);
        }
        sum_err = std::sqrt(sum_err);

        res[it] += (sum_val / sum_err);
        if (hasJacb) {
            for (int itmp = 0; itmp < ntmp_; ++itmp) {
                if (!has_tmps_[itmp][it]) continue;
                jacb[it*npar_+itmp] += (-tmps_[itmp][it][0] / sum_err);
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

#endif // __DataFit2D_C__
