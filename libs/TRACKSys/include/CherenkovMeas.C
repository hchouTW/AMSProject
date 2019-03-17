#ifndef __TRACKLibs_CherenkovMeas_C__
#define __TRACKLibs_CherenkovMeas_C__


#include "Sys.h"
#include "Math.h"
#include "CherenkovMeas.h"


namespace TrackSys {


const double& CherenkovHit::search_closest_beta(double bta) {
    mode_ = -1; beta_ = -1.0;
    if (type_ <= 0) return beta_;

    if      (type_ == 1) { mode_ = 0; beta_ = dbta_; }
    else if (type_ == 2) { mode_ = 1; beta_ = rbtaA_; }
    else if (type_ == 4) { mode_ = 2; beta_ = rbtaB_; }
    else if (type_ == 3) {
        mode_ = (std::fabs(dbta_ - bta) <= std::fabs(rbtaA_ - bta)) ? 0 : 1;
        beta_ = (mode_ == 0 ? dbta_ : rbtaA_);
    }
    else if (type_ == 5) {
        mode_ = (std::fabs(dbta_ - bta) <= std::fabs(rbtaB_ - bta)) ? 0 : 2;
        beta_ = (mode_ == 0 ? dbta_ : rbtaB_);
    }
    else if (type_ == 6) {
        mode_ = (std::fabs(rbtaA_ - bta) <= std::fabs(rbtaB_ - bta)) ? 1 : 2;
        beta_ = (mode_ == 1 ? rbtaA_ : rbtaB_);
    }
    else if (type_ == 7) {
        mode_ = (std::fabs(rbtaA_ - bta) <= std::fabs(rbtaB_ - bta)) ? 1 : 2;
        beta_ = (mode_ == 1 ? rbtaA_ : rbtaB_);
        mode_ = (std::fabs(dbta_ - bta) <= std::fabs(beta_ - bta)) ? 0 : mode_;
        beta_ = (mode_ == 0 ? dbta_ : beta_);
    }
    return beta_;
}
        

void CherenkovHit::set_lmt_bta(double lmtl_bta, double lmtu_bta) {
    bool is_switch = (Numc::Compare(lmtl_bta, lmtu_bta) <= 0) && (lmtl_bta > 0.0 || lmtu_bta > 0.0);
    bool sw_lmtl   = (is_switch && lmtl_bta > 0.0);
    bool sw_lmtu   = (is_switch && lmtu_bta > 0.0);

    type_ = 0;
    if (dbta_  > 0.0 && (!sw_lmtl || dbta_  > lmtl_bta) && (!sw_lmtu || dbta_  < lmtu_bta)) type_ += 1;
    if (rbtaA_ > 0.0 && (!sw_lmtl || rbtaA_ > lmtl_bta) && (!sw_lmtu || rbtaA_ < lmtu_bta)) type_ += 2;
    if (rbtaB_ > 0.0 && (!sw_lmtl || rbtaB_ > lmtl_bta) && (!sw_lmtu || rbtaB_ < lmtu_bta)) type_ += 4;
}

        
CherenkovFit::CherenkovFit(const std::vector<CherenkovHit>& args_hits, double rfr_index, double width_bta, const MultiGaus& pdf_bta, double bta_crr, bool weighted) {
    rfr_index_ = rfr_index;
    width_bta_ = width_bta;
    pdf_bta_   = pdf_bta;
    bta_crr_   = bta_crr;
    weighted_  = weighted;
    if (!check()) { clear(); return; }
    if (args_hits.size() == 0) { clear(); return; }

    timer_.start();

    // (hits) group table
    auto&& gtable = make_group_table(args_hits);
    if (gtable.size() == 0) { clear(); return; }
    
    // build (stone cloud tumor)
    auto&& rlt = build(gtable, weighted_);
    stns_ = std::get<0>(rlt);
    clds_ = std::get<1>(rlt);
    tmrstns_ = std::get<2>(rlt);
    tmrclds_ = std::get<3>(rlt);
            
    // build (npe)
    auto&& npe = build_npe(args_hits, stns_, clds_, tmrstns_, tmrclds_);
    npe_total_ = npe[0];
    npe_stone_ = npe[1];
    npe_cloud_ = npe[2];
    npe_tumor_ = npe[3];
    npe_other_ = npe[4];

    // result
    succ_ = (stns_.size() != 0 || clds_.size() != 0 || tmrstns_.size() != 0 || tmrclds_.size() != 0);
    
    timer_.stop();

    if (!succ_) { clear(); return; }
}


void CherenkovFit::clear() {
    succ_ = false;
    weighted_ = false;

    stns_.clear();
    clds_.clear();
    tmrstns_.clear();
    tmrclds_.clear();
    
    npe_total_ = Numc::ZERO<>;
    npe_stone_ = Numc::ZERO<>;
    npe_cloud_ = Numc::ZERO<>;
    npe_tumor_ = Numc::ZERO<>;
    npe_other_ = Numc::ZERO<>;

    rfr_index_ = Numc::ONE<>;
    width_bta_ = Numc::ZERO<>; 
    scan_bta_  = std::move(MultiGaus());
    pdf_bta_   = std::move(MultiGaus());
    bta_crr_   = Numc::ONE<>; 
    lmtl_bta_  = Numc::ONE<>;
    lmtu_bta_  = Numc::ONE<>;

    timer_.clear();
}


bool CherenkovFit::check() {
    if (Numc::Compare(rfr_index_, Numc::ONE<>) <= 0) return false;
    if (Numc::Compare(width_bta_) <= 0) return false;
    if (Numc::Compare(bta_crr_)   <= 0) return false;
    scan_bta_ = std::move(MultiGaus(Robust::Option(Robust::Opt::ON, 3.0L, 0.5L), width_bta_));
    lmtl_bta_ = (Numc::ONE<> + WIDTH_CORE_COS) / rfr_index_ + width_bta_ * (Numc::THREE<> * Numc::ONE_TO_FOUR);
    lmtu_bta_ = Numc::ONE<> + width_bta_ * Numc::SIX<>;
    if (Numc::Compare(lmtl_bta_) <= 0 || Numc::Compare(lmtu_bta_) <= 0) return false;
    return true;
}


std::array<double, 5> CherenkovFit::build_npe(const std::vector<CherenkovHit>& args_hits, const std::vector<CherenkovStone>& stns, const std::vector<CherenkovCloud>& clds, const std::vector<CherenkovTumor>& tmrstns, const std::vector<CherenkovTumor>& tmrclds) {
    std::array<double, 5> npes({0.0, 0.0, 0.0, 0.0, 0.0});

    double npe_ttl = 0.0;
    for (auto&& hit : args_hits) npe_ttl += hit.npe();
    
    double npe_stn = 0.0;
    std::set<short> chann_stn;
    for (auto&& stn : stns) { 
        for (auto&& hit : stn.hits()) {
            if (chann_stn.find(hit.chann()) != chann_stn.end()) continue;
            chann_stn.insert(hit.chann());
            npe_stn += hit.npe();
        }
    }
    
    double npe_cld = 0.0;
    std::set<short> chann_cld;
    for (auto&& cld : clds) { 
        for (auto&& hit : cld.hits()) {
            if (chann_cld.find(hit.chann()) != chann_cld.end()) continue;
            chann_cld.insert(hit.chann());
            npe_cld += hit.npe();
        }
    }
    
    double npe_tmr = 0.0;
    std::set<short> chann_tmr;
    for (auto&& tmrstn : tmrstns) { 
        for (auto&& hit : tmrstn.hits()) {
            if (chann_tmr.find(hit.chann()) != chann_tmr.end()) continue;
            chann_tmr.insert(hit.chann());
            npe_tmr += hit.npe();
        }
    }
    for (auto&& tmrcld : tmrclds) { 
        for (auto&& hit : tmrcld.hits()) {
            if (chann_tmr.find(hit.chann()) != chann_tmr.end()) continue;
            chann_tmr.insert(hit.chann());
            npe_tmr += hit.npe();
        }
    }
    
    double npe_oth = npe_ttl - (npe_stn + npe_cld + npe_tmr);

    npes = std::move(std::array<double, 5>({ npe_ttl, npe_stn, npe_cld, npe_tmr, npe_oth }));
    return npes;
}


std::vector<std::vector<CherenkovHit>> CherenkovFit::make_group_table(const std::vector<CherenkovHit>& args_hits, bool opt_stone, bool opt_cloud) {
    std::vector<std::vector<CherenkovHit>> gtable;
    if (!opt_stone && !opt_cloud) return gtable;
    if (args_hits.size() == 0) return gtable;

    if (args_hits.size() == 1) {
        gtable.push_back(std::vector<CherenkovHit>({ args_hits.at(0) }));
        return gtable;
    }

    // discontact table
    std::vector<std::vector<short>> discontact_table(args_hits.size(), std::vector<short>(args_hits.size(), 0)); // set to contact
    for (int it = 0; it < args_hits.size(); ++it) discontact_table.at(it).at(it) = 1; // set self-self to discontact

    for (int it = 0; it < args_hits.size() - 1; ++it) {
        for (int jt = it + 1; jt < args_hits.size(); ++jt) {
            const CherenkovHit& ihit = args_hits.at(it);
            const CherenkovHit& jhit = args_hits.at(jt);

            double dist = Numc::INV_SQRT_TWO * std::hypot(ihit.cx() - jhit.cx(), ihit.cy() - jhit.cy()) / WIDTH_CORE_COO;
            if (opt_stone && dist < LMTMAX_GROUP_SGM) continue;

            if (opt_cloud && (ihit.hasDb()  && jhit.hasDb() ) && (std::fabs(ihit.dbta()  - jhit.dbta() ) / width_bta_) < LMTMAX_GROUP_SGM) continue;
            if (opt_cloud && (ihit.hasDb()  && jhit.hasRbA()) && (std::fabs(ihit.dbta()  - jhit.rbtaA()) / width_bta_) < LMTMAX_GROUP_SGM) continue;
            if (opt_cloud && (ihit.hasDb()  && jhit.hasRbB()) && (std::fabs(ihit.dbta()  - jhit.rbtaB()) / width_bta_) < LMTMAX_GROUP_SGM) continue;
            if (opt_cloud && (ihit.hasRbA() && jhit.hasDb() ) && (std::fabs(ihit.rbtaA() - jhit.dbta() ) / width_bta_) < LMTMAX_GROUP_SGM) continue;
            if (opt_cloud && (ihit.hasRbA() && jhit.hasRbA()) && (std::fabs(ihit.rbtaA() - jhit.rbtaA()) / width_bta_) < LMTMAX_GROUP_SGM) continue;
            if (opt_cloud && (ihit.hasRbA() && jhit.hasRbB()) && (std::fabs(ihit.rbtaA() - jhit.rbtaB()) / width_bta_) < LMTMAX_GROUP_SGM) continue;
            if (opt_cloud && (ihit.hasRbB() && jhit.hasDb() ) && (std::fabs(ihit.rbtaB() - jhit.dbta() ) / width_bta_) < LMTMAX_GROUP_SGM) continue;
            if (opt_cloud && (ihit.hasRbB() && jhit.hasRbA()) && (std::fabs(ihit.rbtaB() - jhit.rbtaA()) / width_bta_) < LMTMAX_GROUP_SGM) continue;
            if (opt_cloud && (ihit.hasRbB() && jhit.hasRbB()) && (std::fabs(ihit.rbtaB() - jhit.rbtaB()) / width_bta_) < LMTMAX_GROUP_SGM) continue;

            discontact_table.at(it).at(jt) = 1; // set to discontact
            discontact_table.at(jt).at(it) = 1; // set to discontact
        }
    }
   
    // alone table
    std::vector<short> alone_table(args_hits.size(), 1); // set to alone
    for (int it = 0; it < args_hits.size(); ++it) {
        for (int jt = 0; jt < args_hits.size(); ++jt) {
            if (discontact_table.at(it).at(jt) != 0) continue; // skip if discontact
            alone_table.at(it) = 0; // has contact
            break;
        }
    }
    
    // group table
    std::vector<std::set<int>> group_table;
    for (int it = 0; it < args_hits.size(); ++it) {
        if (alone_table.at(it) != 0) { // alone
            group_table.push_back(std::set<int>({it}));
            continue;
        }
        
        // check
        bool exist = false;
        for (auto&& group : group_table) {
            if (group.find(it) != group.end()) { // found
                exist = true;
                break;
            }
        }
        if (exist) continue;

        // find group
        int current_size = 0;
        std::set<int> group({it});
        while (current_size != group.size()) {
            current_size = group.size();
            std::set<int> tmpgr;
            for (auto&& idx : group) {
                for (int jdx = 0; jdx < args_hits.size(); ++jdx) {
                    if (discontact_table.at(idx).at(jdx) != 0) continue; // skip if discontact
                    if (group.find(jdx) != group.end()) continue; // already exist
                    tmpgr.insert(jdx);
                }
            }
            for (auto&& idx : tmpgr) group.insert(idx);
        }
        group_table.push_back(group);
    }
    
    // (hits) group table
    for (auto&& group : group_table) {
        std::vector<CherenkovHit> hits;
        for (auto&& idx : group) {
            hits.push_back(args_hits.at(idx));
        }
        if (hits.size() == 0) continue;
        gtable.push_back(hits);
    }
    
    return gtable;
}


std::tuple<std::vector<CherenkovStone>, std::vector<CherenkovCloud>, std::vector<CherenkovTumor>, std::vector<CherenkovTumor>> CherenkovFit::build(const std::vector<std::vector<CherenkovHit>>& args_gtable, bool weighted) {
    std::tuple<std::vector<CherenkovStone>, std::vector<CherenkovCloud>, std::vector<CherenkovTumor>, std::vector<CherenkovTumor>> result;
    std::vector<CherenkovStone>& rltstns = std::get<0>(result);
    std::vector<CherenkovCloud>& rltclds = std::get<1>(result);
    std::vector<CherenkovTumor>& rlttmrs_stone = std::get<2>(result);
    std::vector<CherenkovTumor>& rlttmrs_cloud = std::get<3>(result);
    if (args_gtable.size() == 0) return result;

    std::vector<CherenkovHit> other_hits;
    for (auto&& group_hits : args_gtable) {
        if (group_hits.size() < LMTMIN_GROUP_HIT) {
            for (auto&& hit : group_hits) other_hits.push_back(hit);
            continue;
        }
        auto&& sub_rlt = build(group_hits, weighted);
        std::vector<CherenkovStone>& sub_rltstns = sub_rlt.first;
        std::vector<CherenkovCloud>& sub_rltclds = sub_rlt.second;
        for (auto&& stn : sub_rltstns) rltstns.push_back(stn);
        for (auto&& cld : sub_rltclds) rltclds.push_back(cld);

        // find other hits
        for (auto&& hit : group_hits) {
            bool is_used_hit = false;
            
            for (auto&& stn : sub_rltstns) {
                for (auto&& stn_hit : stn.hits()) {
                    if (stn_hit.chann() != hit.chann()) continue;
                    is_used_hit = true;
                    break;
                }
                if (is_used_hit) break;
            }
            if (is_used_hit) continue;
            
            for (auto&& cld : sub_rltclds) {
                for (auto&& cld_hit : cld.hits()) {
                    if (cld_hit.chann() != hit.chann()) continue;
                    is_used_hit = true;
                    break;
                }
                if (is_used_hit) break;
            }
            if (is_used_hit) continue;

            other_hits.push_back(hit);
        }
    }
    if (rltstns.size() > 1) std::sort(rltstns.begin(), rltstns.end(), CherenkovStone_sort());
    if (rltclds.size() > 1) std::sort(rltclds.begin(), rltclds.end(), CherenkovCloud_sort());
    if (other_hits.size() > 1) std::sort(other_hits.begin(), other_hits.end(), CherenkovHit_sort());

    auto&& rlttmrs = build_tumor(other_hits);
    rlttmrs_stone = rlttmrs.first;
    rlttmrs_cloud = rlttmrs.second;

    return result;
}


std::pair<std::vector<CherenkovStone>, std::vector<CherenkovCloud>> CherenkovFit::build(const std::vector<CherenkovHit>& args_hits, bool weighted) {
    std::pair<std::vector<CherenkovStone>, std::vector<CherenkovCloud>> result;
    std::vector<CherenkovStone>& rltstns = result.first;
    std::vector<CherenkovCloud>& rltclds = result.second;
    if (args_hits.size() == 0) return result; 

    // stone
    std::vector<CherenkovHit> all_hits = args_hits;
    
    for (auto&& hit : all_hits) { hit.set_wgt(Numc::ONE<>); hit.set_cnt(Numc::ONE<>); }
    if (weighted) {
        double sumnpe = Numc::ZERO<>;
        for (auto&& hit : all_hits) sumnpe += hit.npe();
        for (auto&& hit : all_hits) hit.set_wgt(static_cast<double>(all_hits.size()) * (hit.npe() / sumnpe));
    }

    rltstns = std::move(fit_stone(all_hits));
    std::vector<CherenkovHit> remainder_hits = all_hits;
    
    while (rltstns.size() != 0) {
        remainder_hits.clear();
        for (auto&& hit : all_hits) {
            bool is_stn = false;
            for (auto&& stn : rltstns) {
                for (auto&& stn_hit : stn.hits()) {
                    if (stn_hit.chann() != hit.chann()) continue;
                    is_stn = true;
                    break;
                }
                if (is_stn) break;
            }
            if (is_stn) continue;
            remainder_hits.push_back(hit);
        }
        std::vector<CherenkovStone>&& other_stns = fit_stone(remainder_hits);
        if (other_stns.size() == 0) break;

        for (auto&& stn : other_stns) rltstns.push_back(stn);
        if (rltstns.size() > 1) std::sort(rltstns.begin(), rltstns.end(), CherenkovStone_sort());
    }
   
    // cloud
    std::vector<CherenkovHit> cld_hits;
    for (auto&& hit : remainder_hits) {
        hit.set_lmt_bta(lmtl_bta_, lmtu_bta_);
        if (hit.type() == 0) continue;
        
        hit.search_closest_beta(Numc::ONE<>);
        if (hit.mode() < 0) continue;

        cld_hits.push_back(hit);
    }

    std::vector<CherenkovCloud>&& cand_clds = fit_cloud(cld_hits);
    for (auto&& cand_cld : cand_clds) {
        CherenkovCloud&& cld = refit_cloud(cand_cld);
        if (!cld.status()) continue;
        rltclds.push_back(cld);
    }
    if (rltclds.size() > 1) std::sort(rltclds.begin(), rltclds.end(), CherenkovCloud_sort());

    return result;
}
        

std::pair<std::vector<CherenkovTumor>, std::vector<CherenkovTumor>> CherenkovFit::build_tumor(const std::vector<CherenkovHit>& args_hits, bool weighted) {
    std::pair<std::vector<CherenkovTumor>, std::vector<CherenkovTumor>> tmrs;
    std::vector<CherenkovTumor>& tmrs_stone = tmrs.first;
    std::vector<CherenkovTumor>& tmrs_cloud = tmrs.second;
    if (args_hits.size() == 0) return tmrs;

    std::vector<CherenkovHit> hits = args_hits;
    for (auto&& hit : hits) hit.set_lmt_bta(lmtl_bta_, lmtu_bta_);
    
    for (auto&& hit : hits) { hit.set_wgt(Numc::ONE<>); hit.set_cnt(Numc::ONE<>); }
    if (weighted) {
        double sumnpe = Numc::ZERO<>;
        for (auto&& hit : hits) sumnpe += hit.npe();
        for (auto&& hit : hits) hit.set_wgt(static_cast<double>(hits.size()) * (hit.npe() / sumnpe));
    }

    // tumor (stone)
    auto&& otable_stone = make_group_table(hits, true, false);
    for (auto otable : otable_stone) {
        if (otable.size() < LMTMIN_TUMOR_STONE_HITS) continue;
        bool is_exist_core = false;
        std::map<int, int> pmt_maps;
        for (auto&& hit : otable) {
            if (pmt_maps.count(hit.pmtid()) == 0) pmt_maps[hit.pmtid()] = 1;
            else pmt_maps[hit.pmtid()]++;
        }
        for (auto&& pmt : pmt_maps) {
            if (pmt.second >= LMTMIN_TUMOR_STONE_PMT_HITS) { is_exist_core = true; }
        }
        if (!is_exist_core) continue;
        
        CherenkovTumor&& tmr = fit_tumor(otable, 0);
        if (!tmr.status()) continue;
        tmrs_stone.push_back(tmr);
    }
    if (tmrs_stone.size() > 1) std::sort(tmrs_stone.begin(), tmrs_stone.end(), CherenkovTumor_sort());

    // tumor (cloud)
    auto&& otable_cloud = make_group_table(hits, false, true);
    for (auto otable : otable_cloud) {
        if (otable.size() < LMTMIN_TUMOR_CLOUD_HITS) continue;
        std::map<int, int> pmt_maps;
        for (auto&& hit : otable) {
            if (pmt_maps.count(hit.pmtid()) == 0) pmt_maps[hit.pmtid()] = 1;
            else pmt_maps[hit.pmtid()]++;
        }
        if (pmt_maps.size() < LMTMIN_TUMOR_CLOUD_PMTS) continue;
        
        CherenkovTumor&& tmr = fit_tumor(otable, 1);
        if (!tmr.status()) continue;
        tmrs_cloud.push_back(tmr);
    }
    if (tmrs_cloud.size() > 1) std::sort(tmrs_cloud.begin(), tmrs_cloud.end(), CherenkovTumor_sort());

    return tmrs;
}
        

CherenkovTumor CherenkovFit::fit_tumor(std::vector<CherenkovHit>& hits, short mode) {
    if (hits.size() == 0 || mode < 0 || mode > 1) return CherenkovTumor();
        
    std::map<int, int> pmt_maps;
    for (auto&& hit : hits) {
        if (pmt_maps.count(hit.pmtid()) == 0) pmt_maps[hit.pmtid()] = 1;
        else pmt_maps[hit.pmtid()]++;
    }
    short nhit = hits.size();
    short npmt = pmt_maps.size();
    if (nhit <= 0 || npmt <= 0) return CherenkovTumor();
        
    double weight = Numc::ZERO<>;
    for (auto&& hit : hits) weight += hit.wgt();
    if (Numc::EqualToZero(weight)) return CherenkovTumor();
    for (auto&& hit : hits) hit.set_cnt(static_cast<double>(hits.size()) * (hit.wgt() / weight));

    double sum_cx = Numc::ZERO<>;
    double sum_cy = Numc::ZERO<>;
    double sum_cw = Numc::ZERO<>;
    double npe    = Numc::ZERO<>;
    for (auto&& hit : hits) {
        sum_cx += hit.cnt() * hit.cx();
        sum_cy += hit.cnt() * hit.cy();
        sum_cw += hit.cnt();
        npe    += hit.npe();
    }
    double cx = sum_cx / sum_cw;
    double cy = sum_cy / sum_cw;

    double init_beta = Numc::ZERO<>;
    double count_beta = Numc::ZERO<>;
    for (auto&& hit : hits) {
        if (hit.type() == 0) continue;
        double multi = hit.cnt() / static_cast<double>(hit.hasDb() + hit.hasRbA() + hit.hasRbB());
        if (hit.hasDb() ) { init_beta += multi * hit.dbta() ; count_beta += multi; }
        if (hit.hasRbA()) { init_beta += multi * hit.rbtaA(); count_beta += multi; }
        if (hit.hasRbB()) { init_beta += multi * hit.rbtaB(); count_beta += multi; }
    }
    init_beta = (!Numc::EqualToZero(count_beta) ? (init_beta / count_beta) : Numc::ZERO<>);

    double beta = init_beta;
    if (Numc::Compare(beta) > 0) {
        bool succ = false;
        short iter = 1;
        while (!succ && iter <= LMTMAX_ITER) {
            double sum_bta = Numc::ZERO<>;
            double sum_cnt = Numc::ZERO<>;
            for (auto&& hit : hits) {
                if (hit.type() == 0) continue;
                sum_bta += hit.cnt() * hit.search_closest_beta(beta);
                sum_cnt += hit.cnt();
            }
            if (Numc::EqualToZero(sum_cnt)) break;
            double new_beta = (sum_bta / sum_cnt);
            succ = (iter >= LMTL_ITER && Numc::Compare(std::fabs(new_beta - beta), CONVG_TOLERANCE) <= 0);
            if (!succ) beta = new_beta;
            iter++;
        }
        if (!succ) beta = Numc::ZERO<>;
    }

    // refit
    bool succ = false;
    double nchi = Numc::ZERO<>;
    for (int iter = 0; iter <= LMTMAX_ITER && (beta > 0.0); ++iter) {
        double count = Numc::ZERO<>;
        double chisq = Numc::ZERO<>;
        double grdB  = Numc::ZERO<>;
        double cvBB  = Numc::ZERO<>;

        for (auto&& hit : hits) {
            double dlt = hit.search_closest_beta(beta) - beta;
            std::array<long double, 3>&& minib = pdf_bta_.minimizer(dlt);
            
            count += hit.cnt();
            chisq += hit.cnt() * (minib[0] * minib[0]);
            grdB  += hit.cnt() * (Numc::NEG<> * minib[2] * minib[1]);
            cvBB  += hit.cnt() * (minib[2] * minib[2]);
        }
        if (!Numc::Valid(count) || Numc::Compare(count, Numc::ONE<>) <= 0) break;

        double new_nchi = chisq / (count - Numc::ONE<>);
        if (!Numc::Valid(new_nchi)) break;

        double resB = grdB / cvBB;
        double new_beta = beta - resB;
        if (!Numc::Valid(resB) || Numc::Compare(new_beta) <= 0) break;
       
        if (iter >= LMTL_ITER) {
            succ = (Numc::Compare(std::fabs(new_beta - beta), CONVG_TOLERANCE) <= 0) && 
                   (Numc::Compare(std::fabs(new_nchi - nchi), CONVG_TOLERANCE) <= 0);
            if (succ) break;
        }

        beta = new_beta;
        nchi = new_nchi;
    }
    double cbta = bta_crr_ * beta;
    if (!succ) return CherenkovTumor();
    
    return CherenkovTumor(hits, mode, nhit, npmt, cx, cy, beta, cbta, npe);
}
        

std::vector<CherenkovStone> CherenkovFit::fit_stone(std::vector<CherenkovHit>& hits) {
    std::vector<CherenkovStone> stns;
    if (hits.size() == 0) return stns;
    
    std::vector<std::array<double, 3>>&& cand_cstns = clustering_stone(hits); // (cx cy wgt)
    if (cand_cstns.size() == 0) return stns;

    std::vector<std::vector<CherenkovHit>> cand_chits(cand_cstns.size(), std::vector<CherenkovHit>());
    for (auto&& hit : hits) {
        std::vector<double> elms;
        double sumprb = Numc::ZERO<>;
        for (auto&& stn : cand_cstns) {
            double dcx = (hit.cx() - stn[0]) / WIDTH_CORE_COO;
            double dcy = (hit.cy() - stn[1]) / WIDTH_CORE_COO;
            double nrm = Numc::INV_SQRT_TWO * std::hypot(dcx,  dcy);
            double pdf = std::exp(Numc::NEG<> * nrm * nrm);
            double prb = stn[2] * pdf;
            
            if (Numc::Compare(prb * hits.size(), CONVG_PROB_SGM70) < 0) {
                elms.push_back(Numc::ZERO<>);
                continue; 
            }
            elms.push_back(prb);
            sumprb += prb;
        }
        if (Numc::Compare(sumprb) <= 0) continue;
        for (auto&& elm : elms) elm = (elm / sumprb);
        
        for (int it = 0; it < cand_chits.size(); ++it) {
            if (Numc::Compare(elms.at(it), CONVG_PROB_SGM30) < 0) continue;
            std::vector<CherenkovHit>& chit = cand_chits.at(it);
            chit.push_back(hit);
        }
    }

    for (int it = 0; it < cand_cstns.size(); ++it) {
        std::array<double, 3>&     cand_cstn = cand_cstns.at(it);
        std::vector<CherenkovHit>& cand_chit = cand_chits.at(it);

        if (cand_chit.size() < LMTMIN_STONE_PMT_HITS_L) continue;
        std::sort(cand_chit.begin(), cand_chit.end(), CherenkovHit_sort());

        bool is_exist_core_l = false;
        bool is_exist_core_h = false;
        std::map<int, int> pmt_maps;
        for (auto&& hit : cand_chit) {
            if (pmt_maps.count(hit.pmtid()) == 0) pmt_maps[hit.pmtid()] = 1;
            else pmt_maps[hit.pmtid()]++;
        }
        for (auto&& pmt : pmt_maps) {
            if (pmt.second >= LMTMIN_STONE_PMT_HITS_L) { is_exist_core_l = true; }
            if (pmt.second >= LMTMIN_STONE_PMT_HITS_H) { is_exist_core_h = true; }
        }
        if (!is_exist_core_l) continue;
        if (!is_exist_core_h && cand_chit.size() < LMTMIN_STONE_HITS) continue;
        
        double weight = Numc::ZERO<>;
        for (auto&& hit : cand_chit) weight += hit.wgt();
        if (Numc::EqualToZero(weight)) continue;
        for (auto&& hit : cand_chit) hit.set_cnt(static_cast<double>(cand_chit.size()) * (hit.wgt() / weight));
        
        short nhit = cand_chit.size();
        short npmt = pmt_maps.size();

        double cnt = 0.;
        double npe = 0.;
        double chisq = 0.;
        double chicb = 0.;
        double chiqd = 0.;
        for (auto&& hit : cand_chit) {
            cnt += hit.cnt();
            npe += hit.npe();
            double dcx = (hit.cx() - cand_cstn[0]) / WIDTH_CORE_COO;
            double dcy = (hit.cy() - cand_cstn[1]) / WIDTH_CORE_COO;
            double nrm = Numc::INV_SQRT_TWO * std::hypot(dcx,  dcy);
            chisq += hit.cnt() * (nrm * nrm); 
            chicb += hit.cnt() * (nrm * nrm * nrm); 
            chiqd += hit.cnt() * (nrm * nrm * nrm * nrm); 
        }

        double ndof = (cnt - Numc::ONE<>);
        double nchi = chisq / ndof;
        double qlt  = Numc::NormQuality(nchi, ndof);
        if (!Numc::Valid(ndof) || Numc::Compare(ndof) <= 0) continue;
        if (!Numc::Valid(nchi)) continue;
        if (!Numc::Valid(qlt)) continue;

        double skewness = std::cbrt(chicb / cnt);
        if (!Numc::Valid(skewness)) continue;
        
        double kurtosis = std::sqrt(std::sqrt(chiqd / cnt));
        if (!Numc::Valid(kurtosis)) continue;
        
        stns.push_back(CherenkovStone(cand_chit, nhit, npmt, cand_cstn[0], cand_cstn[1], npe, cnt, nchi, qlt, skewness, kurtosis));
    }
    if (stns.size() > 1) std::sort(stns.begin(), stns.end(), CherenkovStone_sort());

    return stns;
}


std::vector<std::array<double, 3>> CherenkovFit::clustering_stone(std::vector<CherenkovHit>& hits) {
   std::vector<std::array<double, 3>> cstns; // (cx cy wgt)
    if (hits.size() == 0) return cstns;
    
    double weight = Numc::ZERO<>;
    for (auto&& hit : hits) weight += hit.wgt();
    if (Numc::EqualToZero(weight)) return cstns;
    for (auto&& hit : hits) hit.set_cnt(static_cast<double>(hits.size()) * (hit.wgt() / weight));
    for (auto&& hit : hits) cstns.push_back(std::array<double, 3>({ hit.cx(), hit.cy(), (hit.wgt() / weight) }));
    if (cstns.size() == 0) return cstns;

    bool succ = false;
    for (int iter = 1; iter <= LMTMAX_ITER; ++iter) {
        std::vector<std::array<double, 3>> mgaus(cstns.size(), std::array<double, 3>({ 0.0, 0.0, 0.0 })); // (cx cy wgt)
        double count = Numc::ZERO<>;

        for (auto&& hit : hits) {
            double sumprb = Numc::ZERO<>;
            std::vector<std::array<double, 3>> elms; // (cx cy wgt)
            for (auto&& stn : cstns) {
                if (Numc::Compare(stn[2]) <= 0) {
                    elms.push_back(std::array<double, 3>({ 0.0, 0.0, 0.0 }));
                    continue;
                }
                double nrmx = (hit.cx() - stn[0]) / WIDTH_CORE_COO;
                double nrmy = (hit.cy() - stn[1]) / WIDTH_CORE_COO;
                double nrm  = Numc::INV_SQRT_TWO * std::hypot(nrmx, nrmy);
                double pdf  = std::exp(Numc::NEG<> * nrm * nrm);
                double prb  = stn[2] * pdf;
                
                elms.push_back(std::array<double, 3>({ hit.cx(), hit.cy(), prb }));
                sumprb += prb;
            }
            if (!Numc::Valid(sumprb) || Numc::Compare(sumprb) <= 0) continue;
            for (auto&& elm : elms) elm.at(2) = hit.cnt() * (elm.at(2) / sumprb);
            count += hit.cnt();

            for (int it = 0; it < elms.size(); ++it) {
                std::array<double, 3>& elm  = elms.at(it);
                std::array<double, 3>& gaus = mgaus.at(it);
                if (Numc::Compare(elm[2]) <= 0) continue;
                gaus[0] += elm[2] * elm[0];
                gaus[1] += elm[2] * elm[1];
                gaus[2] += elm[2];
            }
        }
        if (!Numc::Valid(count) || Numc::Compare(count, Numc::ONE<>) <= 0) break;
        
        for (auto&& gaus : mgaus) {
            if (Numc::Compare(gaus[2]) <= 0) continue;
            gaus[0] = gaus[0] / gaus[2];
            gaus[1] = gaus[1] / gaus[2];
            gaus[2] = gaus[2] / count;
        }

        // remove empty
        if (iter >= LMTL_ITER) {
            for (int it = mgaus.size()-1; it >= 0; --it) {
                if (Numc::Compare(mgaus.at(it)[2], CONVG_EPSILON) < 0)
                    mgaus.erase(mgaus.begin() + it);
            }
            if (mgaus.size() == 0) break;
        }
        bool has_remove_empty = (mgaus.size() != cstns.size());
        
        // merge
        bool has_merge = false;
        if (iter >= LMTM_ITER && !has_remove_empty) {
            std::vector<std::vector<int>> sets;
            for (int it = 0; it < mgaus.size(); ++it) {
                std::vector<int>* set_ptr = nullptr;
                for (auto&& set : sets) {
                    for (auto&& idx : set) {
                        if (idx != it) continue;
                        else { set_ptr = &set; break; }
                    }
                    if (set_ptr != nullptr) break;
                }
                if (set_ptr == nullptr) {
                    sets.push_back(std::vector<int>({it}));
                    set_ptr = &(sets.back());
                }
               
                for (int jt = it+1; jt < mgaus.size(); ++jt) {
                    bool is_inside = false;
                    for (auto&& idx : (*set_ptr)) {
                        if (idx == jt) { is_inside = true; break; }
                    }
                    if (is_inside) continue;

                    auto&& gausA = mgaus.at(it);
                    auto&& gausB = mgaus.at(jt);
                    double prbA = gausA[2];
                    double prbB = gausB[2];
                    double newx = (gausA[0] * gausA[2] + gausB[0] * gausB[2]) / (gausA[2] + gausB[2]);
                    double newy = (gausA[1] * gausA[2] + gausB[1] * gausB[2]) / (gausA[2] + gausB[2]);
                    double dnxA = (newx - gausA[0]) / WIDTH_CORE_COO;
                    double dnyA = (newy - gausA[1]) / WIDTH_CORE_COO;
                    double dnxB = (newx - gausB[0]) / WIDTH_CORE_COO;
                    double dnyB = (newy - gausB[1]) / WIDTH_CORE_COO;
                    double dnA  = Numc::INV_SQRT_TWO * std::sqrt(dnxA * dnxA + dnyA * dnyA);
                    double dnB  = Numc::INV_SQRT_TWO * std::sqrt(dnxB * dnxB + dnyB * dnyB);
                    double prb  = gausA[2] * std::exp(Numc::NEG<> * dnA * dnA) + 
                                  gausB[2] * std::exp(Numc::NEG<> * dnB * dnB);
                    bool merge = (prb > CONVG_CLOSED * prbA || prb > CONVG_CLOSED * prbB);
                    if (merge) set_ptr->push_back(jt);
                }
            }

            has_merge = (mgaus.size() != sets.size());
            if (has_merge) {
                std::vector<std::array<double, 3>> newmg;
                for (auto&& set : sets) {
                    if (set.size() == 1) {
                        newmg.push_back(mgaus.at(set.at(0))); 
                        continue; 
                    }
                    double newcx  = Numc::ZERO<>;
                    double newcy  = Numc::ZERO<>;
                    double newwgt = Numc::ZERO<>;
                    for (auto&& idx : set) {
                        auto&& gaus = mgaus.at(idx);
                        newcx  += gaus[2] * gaus[0];
                        newcy  += gaus[2] * gaus[1];
                        newwgt += gaus[2];
                    }
                    newcx = (newcx / newwgt);
                    newcy = (newcy / newwgt);
                    newmg.push_back(std::array<double, 3>({ newcx, newcy, newwgt }));
                }
                mgaus = std::move(newmg);
            }
        }

        double normalized = Numc::ZERO<>;
        for (auto&& gaus : mgaus) normalized += gaus[2];
        for (auto&& gaus : mgaus) gaus[2] = gaus[2] / normalized;
        if (Numc::Compare(normalized) <= 0) break;

        if (iter >= LMTU_ITER && !has_remove_empty && !has_merge) {
            bool is_convg = true;
            for (int it = 0; it < mgaus.size(); ++it) {
                auto&& gaus = mgaus.at(it);
                auto&& stn  = cstns.at(it);
                bool match = (Numc::Compare(std::fabs(gaus[0] - stn[0]), CONVG_TOLERANCE) <= 0 && 
                              Numc::Compare(std::fabs(gaus[1] - stn[1]), CONVG_TOLERANCE) <= 0 && 
                              Numc::Compare(std::fabs(gaus[2] - stn[2]), CONVG_TOLERANCE) <= 0);
                if (!match) { is_convg = false; break; }
            }
            if (is_convg) { succ = true; break; }
        }
        cstns = mgaus;
    }
    
    if (!succ) cstns.clear();
    return cstns;
}
    

CherenkovCloud CherenkovFit::refit_cloud(const CherenkovCloud& cand_cld) {
    if (!cand_cld.status()) return CherenkovCloud();
    std::vector<CherenkovHit> cand_chit = cand_cld.hits();
    if (cand_chit.size() == 0) return CherenkovCloud();

    double weight = Numc::ZERO<>;
    for (auto&& hit : cand_chit) weight += hit.wgt();
    if (Numc::EqualToZero(weight)) return CherenkovCloud();
    for (auto&& hit : cand_chit) hit.set_cnt(static_cast<double>(cand_chit.size()) * (hit.wgt() / weight));

    double cand_beta = cand_cld.beta();
    double cand_cnt  = Numc::ZERO<>;
    double cand_nchi = Numc::ZERO<>;

    double cand_chisq = Numc::ZERO<>;
    double cand_chicb = Numc::ZERO<>;
    double cand_chiqd = Numc::ZERO<>;
    
    // refit cloud
    bool succ = false;
    for (int iter = 0; iter <= LMTMAX_ITER; ++iter) {
        double count = Numc::ZERO<>;
        double chisq = Numc::ZERO<>;
        double chicb = Numc::ZERO<>;
        double chiqd = Numc::ZERO<>;
        double grdB  = Numc::ZERO<>;
        double cvBB  = Numc::ZERO<>;

        for (auto&& hit : cand_chit) {
            double dlt = hit.search_closest_beta(cand_beta) - cand_beta;
            double nrm = dlt / width_bta_;
            
            if (Numc::Compare(std::fabs(nrm), Numc::FIVE<>) > 0) continue;
            std::array<long double, 3>&& minib = pdf_bta_.minimizer(dlt);
            
            count += hit.cnt();
            chisq += hit.cnt() * (minib[0] * minib[0]);
            chicb += hit.cnt() * (minib[0] * minib[0] * minib[0]);
            chiqd += hit.cnt() * (minib[0] * minib[0] * minib[0] * minib[0]);
            grdB  += hit.cnt() * (Numc::NEG<> * minib[2] * minib[1]);
            cvBB  += hit.cnt() * (minib[2] * minib[2]);
        }
        if (!Numc::Valid(count) || Numc::Compare(count, Numc::ONE<>) <= 0) break;

        double nchi = chisq / (count - Numc::ONE<>);
        if (!Numc::Valid(nchi)) break;

        double resB = grdB / cvBB;
        double beta = cand_beta - resB;
        if (!Numc::Valid(resB) || Numc::Compare(beta) <= 0) break;
       
        if (iter >= LMTL_ITER) {
            succ = (Numc::Compare(std::fabs(cand_beta - beta), CONVG_TOLERANCE) <= 0) && 
                   (Numc::Compare(std::fabs(cand_cnt - count), CONVG_TOLERANCE) <= 0) &&
                   (Numc::Compare(std::fabs(cand_nchi - nchi), CONVG_TOLERANCE) <= 0);
            if (succ) break;
        }

        cand_beta = beta;
        cand_cnt  = count;
        cand_nchi = nchi;

        cand_chisq = chisq;
        cand_chicb = chicb;
        cand_chiqd = chiqd;
    }
    if (!succ) return CherenkovCloud();
  
    // check result
    int count_within_three_sigma = 0;
    std::vector<CherenkovHit> cand_reduce_chit;
    for (auto&& hit : cand_chit) {
        double absnrm = std::fabs(hit.search_closest_beta(cand_beta) - cand_beta) / width_bta_;
        if (Numc::Compare(absnrm, Numc::FIVE<>) > 0) continue;
        cand_reduce_chit.push_back(hit);

        if (Numc::Compare(absnrm, Numc::THREE<>) < 0) count_within_three_sigma++;
    } 
    if (cand_reduce_chit.size() < LMTMIN_CLOUD_HITS) return CherenkovCloud();
    std::sort(cand_reduce_chit.begin(), cand_reduce_chit.end(), CherenkovHit_sort());
    
    if (count_within_three_sigma  < LMTMIN_CLOUD_HITS) return CherenkovCloud();
    
    std::map<int, int> pmt_maps;
    int count_pmt_over_two_hits = 0;
    for (auto&& hit : cand_reduce_chit) {
        if (pmt_maps.count(hit.pmtid()) == 0) pmt_maps[hit.pmtid()] = 1;
        else pmt_maps[hit.pmtid()]++;
    }
    for (auto&& pmt : pmt_maps) {
        if (pmt.second < LMTMIN_CLOUD_SETS) continue;
        count_pmt_over_two_hits++;
    }
    if (pmt_maps.size() < LMTMIN_CLOUD_PMTS) return CherenkovCloud();
  
    short cand_nhit = cand_reduce_chit.size();
    short cand_npmt = pmt_maps.size();
    
    double cand_npe = 0.;
    for (auto&& hit : cand_reduce_chit) cand_npe += hit.npe();

    double cand_cbta = cand_beta * bta_crr_;
    double cand_ndof = cand_cnt - Numc::ONE<>;
    double cand_qlt  = Numc::NormQuality(cand_nchi, cand_ndof);
    if (!Numc::Valid(cand_ndof) || Numc::Compare(cand_ndof) <= 0) return CherenkovCloud();
    if (!Numc::Valid(cand_nchi)) return CherenkovCloud();
    if (!Numc::Valid(cand_qlt)) return CherenkovCloud();
        
    double cand_skewness = std::cbrt(cand_chicb / cand_cnt);
    if (!Numc::Valid(cand_skewness)) return CherenkovCloud();
    
    double cand_kurtosis = std::sqrt(std::sqrt(cand_chiqd / cand_cnt));
    if (!Numc::Valid(cand_kurtosis)) return CherenkovCloud();
        
    return CherenkovCloud(cand_reduce_chit, cand_nhit, cand_npmt, cand_beta, cand_cbta, cand_npe, cand_cnt, cand_nchi, cand_qlt, cand_skewness, cand_kurtosis);
}


std::vector<CherenkovCloud> CherenkovFit::fit_cloud(std::vector<CherenkovHit>& hits) {
    std::vector<CherenkovCloud> cclds;
    if (hits.size() == 0) return cclds;
    
    std::vector<std::array<double, 2>>&& cand_cclds = clustering_cloud(hits); // (bta wgt)
    if (cand_cclds.size() == 0) return cclds;
    
    std::vector<std::vector<CherenkovHit>> cand_chits(cand_cclds.size(), std::vector<CherenkovHit>());
    for (auto&& hit : hits) {
        std::vector<double> elms;
        double sumprb = Numc::ZERO<>;
        for (auto&& cld : cand_cclds) {
            double nrm = (hit.search_closest_beta(cld[0]) - cld[0]) / width_bta_;
            double pdf = std::exp(-Numc::ONE_TO_TWO * nrm * nrm);
            double prb = cld[1] * pdf;
            
            if (Numc::Compare(prb * hits.size(), CONVG_PROB_SGM70) < 0) {
                elms.push_back(Numc::ZERO<>);
                continue; 
            }
            elms.push_back(prb);
            sumprb += prb;
        }
        if (Numc::Compare(sumprb) <= 0) continue;
        for (auto&& elm : elms) elm = (elm / sumprb);
        
        for (int it = 0; it < cand_chits.size(); ++it) {
            if (Numc::Compare(elms.at(it), CONVG_PROB_SGM30) < 0) continue;
            std::vector<CherenkovHit>& chit = cand_chits.at(it); 
            chit.push_back(hit);
        }
    }

    for (int it = 0; it < cand_cclds.size(); ++it) {
        std::array<double, 2>&     cand_ccld = cand_cclds.at(it);
        std::vector<CherenkovHit>& cand_chit = cand_chits.at(it);

        if (cand_chit.size() < LMTMIN_CLOUD_HITS) continue;
        std::sort(cand_chit.begin(), cand_chit.end(), CherenkovHit_sort());

        std::map<int, int> pmt_maps;
        int count_pmt_over_two_hits = 0;
        for (auto&& hit : cand_chit) {
            if (pmt_maps.count(hit.pmtid()) == 0) pmt_maps[hit.pmtid()] = 1;
            else pmt_maps[hit.pmtid()]++;
        }
        for (auto&& pmt : pmt_maps) {
            if (pmt.second < LMTMIN_CLOUD_SETS) continue;
            count_pmt_over_two_hits++;
        }
        if (pmt_maps.size() < LMTMIN_CLOUD_PMTS) continue;
        
        double weight = Numc::ZERO<>;
        for (auto&& hit : cand_chit) weight += hit.wgt();
        if (Numc::EqualToZero(weight)) continue;
        for (auto&& hit : cand_chit) hit.set_cnt(static_cast<double>(cand_chit.size()) * (hit.wgt() / weight));
        
        short nhit = cand_chit.size();
        short npmt = pmt_maps.size();

        double beta = cand_ccld[0];
        double cbta = cand_ccld[0] * bta_crr_;

        double cnt = 0.;
        double npe = 0.;
        double chisq = 0.;
        double chicb = 0.;
        double chiqd = 0.;
        for (auto&& hit : cand_chit) {
            cnt += hit.cnt();
            npe += hit.npe();
            
            double dlt = hit.search_closest_beta(cand_ccld[0]) - cand_ccld[0];
            std::array<long double, 3>&& minib = scan_bta_.minimizer(dlt);
            chisq += hit.cnt() * (minib[0] * minib[0]);
            chicb += hit.cnt() * (minib[0] * minib[0] * minib[0]);
            chiqd += hit.cnt() * (minib[0] * minib[0] * minib[0] * minib[0]);
        }

        double ndof = (cnt - Numc::ONE<>);
        double nchi = chisq / ndof;
        double qlt  = Numc::NormQuality(nchi, ndof);
        if (!Numc::Valid(ndof) || Numc::Compare(ndof) <= 0) continue;
        if (!Numc::Valid(nchi)) continue;
        if (!Numc::Valid(qlt)) continue;
    
        double skewness = std::cbrt(chicb / cnt);
        if (!Numc::Valid(skewness)) continue;
        
        double kurtosis = std::sqrt(std::sqrt(chiqd / cnt));
        if (!Numc::Valid(kurtosis)) continue;
        
        cclds.push_back(CherenkovCloud(cand_chit, nhit, npmt, beta, cbta, npe, cnt, nchi, qlt, skewness, kurtosis));
    }
    if (cclds.size() > 1) std::sort(cclds.begin(), cclds.end(), CherenkovCloud_sort());

    return cclds;
}


std::vector<std::array<double, 2>> CherenkovFit::clustering_cloud(std::vector<CherenkovHit>& hits) {
    // bayesian information criterion (BIC)
    std::vector<std::array<double, 2>> cclds; // (bta wgt)
    if (hits.size() == 0) return cclds;
    
    double weight_hits = Numc::ZERO<>;
    for (auto&& hit : hits) weight_hits += hit.wgt();
    if (Numc::EqualToZero(weight_hits)) return cclds;
    for (auto&& hit : hits) hit.set_cnt(static_cast<double>(hits.size()) * (hit.wgt() / weight_hits));

    for (auto&& hit : hits) {
        if (hit.type() == 0) continue;
        double multi = hit.cnt() / static_cast<double>(hit.hasDb() + hit.hasRbA() + hit.hasRbB());
        if (hit.hasDb() ) cclds.push_back(std::array<double, 2>({ hit.dbta() , multi }));
        if (hit.hasRbA()) cclds.push_back(std::array<double, 2>({ hit.rbtaA(), multi }));
        if (hit.hasRbB()) cclds.push_back(std::array<double, 2>({ hit.rbtaB(), multi }));
    }
    if (cclds.size() > 1) std::sort(cclds.begin(), cclds.end());
    if (cclds.size() == 0) return cclds;

    double weight = Numc::ZERO<>;
    for (auto&& cld : cclds) weight += cld[1];
    for (auto&& cld : cclds) cld[1] = (cld[1] / weight);
    if (Numc::EqualToZero(weight)) { cclds.clear(); return cclds; }

    bool succ = false;
    for (int iter = 1; iter <= LMTMAX_ITER; ++iter) {
        std::vector<std::array<double, 4>> mgaus(cclds.size(), std::array<double, 4>({ 0.0, 0.0, 0.0, 0.0 })); // (wgt grd cv res)
        double count = Numc::ZERO<>;

        for (auto&& hit : hits) {
            double sumprb = Numc::ZERO<>;
            std::vector<std::array<double, 3>> elms; // (wgt grd cv)
            for (auto&& cld : cclds) {
                if (Numc::Compare(cld[1]) <= 0) {
                    elms.push_back(std::array<double, 3>({ 0.0, 0.0, 0.0 }));
                    continue;
                }
                double dlt = hit.search_closest_beta(cld[0]) - cld[0];
                double nrm = dlt / width_bta_;
                double pdf = std::exp(-Numc::ONE_TO_TWO * nrm * nrm);
                double prb = cld[1] * pdf;
        
                std::array<long double, 3>&& minib = scan_bta_.minimizer(dlt);
                double grdB = (Numc::NEG<> * minib[2] * minib[1]);
                double cvBB = (minib[2] * minib[2]);
                
                elms.push_back(std::array<double, 3>({ prb, grdB, cvBB }));
                sumprb += prb;
            }
            if (!Numc::Valid(sumprb) || Numc::Compare(sumprb) <= 0) continue;
            for (auto&& elm : elms) elm.at(0) = hit.cnt() * (elm.at(0) / sumprb);
            count += hit.cnt();

            for (int it = 0; it < elms.size(); ++it) {
                std::array<double, 3>& elm  = elms.at(it);
                std::array<double, 4>& gaus = mgaus.at(it);
                if (Numc::Compare(elm[0]) <= 0) continue;
                gaus[0] += elm[0];
                gaus[1] += elm[0] * elm[1];
                gaus[2] += elm[0] * elm[2];
            }
        }
        if (!Numc::Valid(count) || Numc::Compare(count, Numc::ONE<>) <= 0) break;

        for (auto&& gaus : mgaus) {
            if (Numc::Compare(gaus[0]) <= 0) continue;
            gaus[3] = gaus[1] / gaus[2];
            gaus[0] = gaus[0] / count;
        }
       
        // new cclds
        std::vector<std::array<double, 2>> new_cclds; // (bta wgt)
        for (int it = 0; it < cclds.size(); ++it) {
            double bta = cclds.at(it)[0] - mgaus.at(it)[3];
            double wgt = mgaus.at(it)[0];
            new_cclds.push_back(std::array<double, 2>({ bta, wgt }));
        }

        // remove empty
        if (iter >= LMTL_ITER) {
            for (int it = new_cclds.size()-1; it >= 0; --it) {
                if (Numc::Compare(new_cclds.at(it)[1], CONVG_EPSILON) < 0)
                    new_cclds.erase(new_cclds.begin() + it);
            }
            if (new_cclds.size() == 0) break;
        }
        bool has_remove_empty = (new_cclds.size() != cclds.size());
        if (new_cclds.size() > 1) std::sort(new_cclds.begin(), new_cclds.end());
        
        // merge
        bool has_merge = false;
        if (iter >= LMTM_ITER && !has_remove_empty) {
            std::vector<std::array<int, 2>> sets;
            for (int it = 0; it < new_cclds.size(); ++it) {
                auto&& gausA = (it == 0) ? new_cclds.at(it) : new_cclds.at(it-1);
                auto&& gausB = new_cclds.at(it);
                double prbA = gausA[1];
                double prbB = gausB[1];

                double newx = (gausA[0] * gausA[1] + gausB[0] * gausB[1]) / (gausA[1] + gausB[1]);
                double dnxA = (newx - gausA[0]) / width_bta_;
                double dnxB = (newx - gausB[0]) / width_bta_;
                double prbx = gausA[1] * std::exp(-Numc::ONE_TO_TWO * dnxA * dnxA) + 
                              gausB[1] * std::exp(-Numc::ONE_TO_TWO * dnxB * dnxB);
                bool merge = (prbx > CONVG_CLOSED * prbA || prbx > CONVG_CLOSED * prbB);
        
                if (it != 0 && merge) { sets.back()[1] = it; continue; }
                sets.push_back(std::array<int, 2>({it, it}));
            }
            has_merge = (new_cclds.size() != sets.size());
            if (has_merge) {
                std::vector<std::array<double, 2>> newcclds;
                for (auto&& set : sets) {
                    if (set.at(0) == set.at(1)) {
                        newcclds.push_back(new_cclds.at(set.at(0))); 
                        continue; 
                    }
                    double newbta = Numc::ZERO<>;
                    double newwgt = Numc::ZERO<>;
                    for (int it = set.at(0); it <= set.at(1); ++it) {
                        auto&& gaus = new_cclds.at(it);
                        newbta += gaus[1] * gaus[0];
                        newwgt += gaus[1];
                    }
                    newbta = (newbta / newwgt);
                    newcclds.push_back(std::array<double, 2>({ newbta, newwgt }));
                }
                new_cclds = std::move(newcclds);
            }
        }

        double normalized = Numc::ZERO<>;
        for (auto&& cld : new_cclds) normalized += cld[1];
        for (auto&& cld : new_cclds) cld[1] = (cld[1] / normalized);
        if (Numc::Compare(normalized) <= 0) break;

        if (iter >= LMTU_ITER && !has_remove_empty && !has_merge) {
            bool is_convg = true;
            for (int it = 0; it < new_cclds.size(); ++it) {
                auto&& gaus = new_cclds.at(it);
                auto&& cld  = cclds.at(it);
                bool match = (Numc::Compare(std::fabs(gaus[0] - cld[0]), CONVG_TOLERANCE) <= 0 && 
                              Numc::Compare(std::fabs(gaus[1] - cld[1]), CONVG_TOLERANCE) <= 0);
                if (!match) { is_convg = false; break; }
            }
            if (is_convg) { succ = true; break; }
        }
        cclds = new_cclds;
    }

    if (!succ) cclds.clear();
    return cclds;
}


std::array<long double, 3> CherenkovMeas::minimizer(long double x, long double ibta) const {
    if (Numc::Compare(ibta, Numc::ONE<long double>) < 0)
        return std::array<long double, 3>({ Numc::ZERO<long double>, Numc::ZERO<long double>, Numc::ZERO<long double> });

    bool rescl = false;
    long double eftsgm = Numc::ONE<long double>;
    if (Numc::Compare(rfr_, Numc::ONE<long double>) > 0) {
        long double thres = rfr_ - Numc::TWO<long double> * sgm_;
        if (Numc::Compare(ibta, thres) > 0) {
            long double extsgm = (ibta - thres) / sgm_;
            eftsgm = Numc::ONE<long double> + extsgm * extsgm;
            rescl = true;
        }
    }

    long double sclx = (rescl ? (x / eftsgm) : x);
    return mgs_.minimizer(sclx);
}


} // namesapce TrackSys


#endif // __TRACKLibs_CherenkovMeas_C__
