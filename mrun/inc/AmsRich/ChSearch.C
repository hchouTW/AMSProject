#ifndef __ChSearch_C__
#define __ChSearch_C__


#include "ChSearch.h"


const double& ChHit::search_closest_beta(double bta) {
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
        

void ChHit::set_lmt_bta(double lmtl_bta, double lmtu_bta) {
    bool is_switch = (lmtl_bta <= lmtu_bta) && (lmtl_bta > 0.0 || lmtu_bta > 0.0);
    bool sw_lmtl   = (is_switch && lmtl_bta > 0.0);
    bool sw_lmtu   = (is_switch && lmtu_bta > 0.0);

    type_ = 0;
    if (dbta_  > 0.0 && (!sw_lmtl || dbta_  > lmtl_bta) && (!sw_lmtu || dbta_  < lmtu_bta)) type_ += 1;
    if (rbtaA_ > 0.0 && (!sw_lmtl || rbtaA_ > lmtl_bta) && (!sw_lmtu || rbtaA_ < lmtu_bta)) type_ += 2;
    if (rbtaB_ > 0.0 && (!sw_lmtl || rbtaB_ > lmtl_bta) && (!sw_lmtu || rbtaB_ < lmtu_bta)) type_ += 4;
    
    search_closest_beta(1.0);
}

   
ChFit::ChFit(const std::vector<ChHit>& args_hits, const std::array<double, 2>& pmtc, double rfr_index, double width_bta, double bta_crr) : ChFit() {
    pmtc_      = pmtc;
    rfr_index_ = rfr_index;
    width_bta_ = width_bta;
    bta_crr_   = bta_crr;
    if (!check()) { clear(); return; }
    if (args_hits.size() == 0) { succ_ = true; return; }

    // build (stone cloud)
    auto&& rlt = build(args_hits);
    stns_ = std::get<0>(rlt);
    clds_ = std::get<1>(rlt);
    ghts_ = std::get<2>(rlt);
    hits_ = std::get<3>(rlt);
   
    for (auto&& hit : hits_) {
        if (hit.cluster() == ChHit::Cluster::stone) { nhit_stone_++; npe_stone_ += hit.npe(); }
        if (hit.cluster() == ChHit::Cluster::cloud) { nhit_cloud_++; npe_cloud_ += hit.npe(); }
        if (hit.cluster() == ChHit::Cluster::ghost) { nhit_ghost_++; npe_ghost_ += hit.npe(); }
        if (hit.cluster() == ChHit::Cluster::other_inn) { nhit_other_inn_++; npe_other_inn_ += hit.npe(); }
        if (hit.cluster() == ChHit::Cluster::other_out) { nhit_other_out_++; npe_other_out_ += hit.npe(); }
    }
    nhit_total_ = nhit_stone_ + nhit_cloud_ + nhit_ghost_ + nhit_other_inn_ + nhit_other_out_;
    npe_total_  =  npe_stone_ +  npe_cloud_ +  npe_ghost_ +  npe_other_inn_ +  npe_other_out_;

    succ_ = true;
    if (!succ_) { clear(); return; }
}


void ChFit::clear() {
    succ_ = false;

    stns_.clear();
    clds_.clear();
    ghts_.clear();
    hits_.clear();
    
    nhit_total_ = 0;
    nhit_stone_ = 0;
    nhit_cloud_ = 0;
    nhit_ghost_ = 0;
    nhit_other_inn_ = 0;
    nhit_other_out_ = 0;

    npe_total_ = 0.0;
    npe_stone_ = 0.0;
    npe_cloud_ = 0.0;
    npe_ghost_ = 0.0;
    npe_other_inn_ = 0.0;
    npe_other_out_ = 0.0;

    pmtc_      = std::array<double, 2>({ 0.0, 0.0 });
    rfr_index_ = 1.0;
    width_bta_ = 0.0; 
    bta_crr_   = 1.0; 
    
    lmtl_bta_  = 1.0;
    lmtu_bta_  = 1.0;
}


bool ChFit::check() {
    if (rfr_index_ <= 1.0) return false;
    if (width_bta_ <= 0.0) return false;
    if (bta_crr_   <= 0.0) return false;
    lmtl_bta_ = (1.0 + 0.75 * WIDTH_CORE_COS) / rfr_index_ + width_bta_;
    lmtu_bta_ = 1.0 + width_bta_ * 6.0;
    if (lmtl_bta_ <= 0.0 || lmtu_bta_ <= 0.0) return false;
    return true;
}

std::tuple<std::vector<ChStone>, std::vector<ChCloud>, std::vector<ChCloud>, std::vector<ChHit>> ChFit::build(const std::vector<ChHit>& args_hits) {
    // stone
    std::vector<ChHit> candiate_stone_hits = args_hits;
    std::vector<ChStone>&& stones = build_stone(candiate_stone_hits);

    std::vector<ChHit> stones_hits;
    for (auto&& stone : stones) {
    for (auto&& hit : stone.hits()) {
        bool is_found_same_hit = false;
        for (auto&& oldhit : stones_hits) {
            if (hit.chann() != oldhit.chann()) continue;
            is_found_same_hit = true;
            break;
        }
        if (is_found_same_hit) continue;
        stones_hits.push_back(hit);
    }}
    if (stones_hits.size() != 0) std::sort(stones_hits.begin(), stones_hits.end(), ChHit_sort());

    // candiate cloud
    std::vector<ChHit> candiate_cloud_hits;
    for (auto&& hit : args_hits) {
        if (hit.type() == 0) continue;

        bool is_stone_hit = false;
        for (auto&& stone : stones) {
            for (auto&& stn_hit : stone.hits()) {
                if (hit.chann() != stn_hit.chann()) continue;
                is_stone_hit = true;
                break;
            }
            if (is_stone_hit) break;
        }
        if (is_stone_hit) continue;
        
        if (is_within_pmtc(hit.lx(), hit.ly())) continue;
        if (!is_within_detectable(hit.lx(), hit.ly())) continue;

        ChHit chhit = hit;
        chhit.set_lmt_bta(lmtl_bta_, lmtu_bta_);
        if (chhit.type() == 0) continue;
        
        candiate_cloud_hits.push_back(chhit);
    }
    std::vector<ChCloud>&& cand_clouds = build_cloud(candiate_cloud_hits);

    // cloud & ghost
    std::tuple<std::vector<ChCloud>, std::vector<ChCloud>>&& cgpair = build_cloud_and_ghost(cand_clouds);
    std::vector<ChCloud>& clouds = std::get<0>(cgpair);
    std::vector<ChCloud>& ghosts = std::get<1>(cgpair);
    
    std::vector<ChHit> clouds_hits;
    for (auto&& cloud : clouds) {
    for (auto&& hit : cloud.hits()) {
        bool is_found_same_hit = false;
        for (auto&& oldhit : clouds_hits) {
            if (hit.chann() != oldhit.chann()) continue;
            is_found_same_hit = true;
            break;
        }
        if (is_found_same_hit) continue;
        clouds_hits.push_back(hit);
    }}
    if (clouds_hits.size() != 0) std::sort(clouds_hits.begin(), clouds_hits.end(), ChHit_sort());
    
    std::vector<ChHit> ghosts_hits;
    for (auto&& ghost : ghosts) {
    for (auto&& hit : ghost.hits()) {
        bool is_exist_in_clouds = false;
        for (auto&& cld_hit : clouds_hits) {
            if (hit.chann() != cld_hit.chann()) continue;
            is_exist_in_clouds = true;
            break;
        }
        if (is_exist_in_clouds) continue;
        
        bool is_found_same_hit = false;
        for (auto&& oldhit : ghosts_hits) {
            if (hit.chann() != oldhit.chann()) continue;
            is_found_same_hit = true;
            break;
        }
        if (is_found_same_hit) continue;
        ghosts_hits.push_back(hit);
    }}
    if (ghosts_hits.size() != 0) std::sort(ghosts_hits.begin(), ghosts_hits.end(), ChHit_sort());

    // other
    std::vector<ChHit> others_inn_hits;
    std::vector<ChHit> others_out_hits;
    for (auto&& hit : args_hits) {
        bool is_stone_hit = false;
        for (auto&& stone_hit : stones_hits) {
            if (hit.chann() != stone_hit.chann()) continue;
            is_stone_hit = true;
            break;
        }
        if (is_stone_hit) continue;
       
        bool is_cloud_hit = false;
        for (auto&& cloud_hit : clouds_hits) {
            if (hit.chann() != cloud_hit.chann()) continue;
            is_cloud_hit = true;
            break;
        }
        if (is_cloud_hit) continue;

        bool is_ghost_hit = false;
        for (auto&& ghost_hit : ghosts_hits) {
            if (hit.chann() != ghost_hit.chann()) continue;
            is_ghost_hit = true;
            break;
        }
        if (is_ghost_hit) continue;

        if (is_within_pmtc(hit.lx(), hit.ly())) continue;
        if (!is_within_detectable(hit.lx(), hit.ly())) continue;
        
        ChHit chhit = hit;
        chhit.set_lmt_bta(lmtl_bta_, lmtu_bta_);
        if (chhit.type() > 0) chhit.set_cluster(ChHit::Cluster::other_inn);
        else                  chhit.set_cluster(ChHit::Cluster::other_out);
       
        if (chhit.cluster() == ChHit::Cluster::other_inn) others_inn_hits.push_back(chhit);
        if (chhit.cluster() == ChHit::Cluster::other_out) others_out_hits.push_back(chhit);
    }
    if (others_inn_hits.size() != 0) std::sort(others_inn_hits.begin(), others_inn_hits.end(), ChHit_sort());
    if (others_out_hits.size() != 0) std::sort(others_out_hits.begin(), others_out_hits.end(), ChHit_sort());

    std::vector<ChHit> chhits;
    for (auto&& hit : stones_hits) chhits.push_back(hit);
    for (auto&& hit : clouds_hits) chhits.push_back(hit);
    for (auto&& hit : ghosts_hits) chhits.push_back(hit);
    for (auto&& hit : others_inn_hits) chhits.push_back(hit);
    for (auto&& hit : others_out_hits) chhits.push_back(hit);

    return std::make_tuple(stones, clouds, ghosts, chhits);
}

std::vector<ChStone> ChFit::build_stone(const std::vector<ChHit>& args_hits) {
    if (args_hits.size() == 0) return std::vector<ChStone>();
    std::map<int, std::set<int>> group_pmts;
    for (int it = 0; it < args_hits.size(); ++it) {
        const ChHit& hit = args_hits.at(it);
        group_pmts[hit.pmtid()].insert(it);
    }

    const int MIN_CONTACT_IN_PMT = 2;
    const double MIN_CONTACT_WIDTH = 1.73205e+00 * WIDTH_CELL;

    // search group in self-pmt
    std::set<std::set<int>> groups;
    for (auto&& pmt : group_pmts) {
        std::vector<int>           gr_cnt;
        std::vector<std::set<int>> gr_idx;
        for (auto&& it : pmt.second) {
        for (auto&& jt : pmt.second) {
            if (it >= jt) continue;
            const ChHit& ihit = args_hits.at(it);
            const ChHit& jhit = args_hits.at(jt);
            int iloc = (ihit.chann() - ihit.pmtid() * 16);
            int jloc = (jhit.chann() - jhit.pmtid() * 16);
            int dist = std::abs((iloc / 4) - (jloc / 4)) + std::abs((iloc % 4) - (jloc % 4));
            if (dist > 1) continue;

            int gr_found = -1;
            for (int igr = 0; igr < gr_idx.size(); ++igr) {
                bool is_found = false;
                for (auto&& idx : gr_idx.at(igr)) {
                    if (idx != it && idx != jt) continue;
                    is_found = true;
                    break;
                }
                if (!is_found) continue;
                gr_found = igr;
                break;
            }

            if (gr_found < 0) {
                gr_cnt.push_back(1);
                gr_idx.push_back({ it, jt });
            }
            else {
                gr_cnt.at(gr_found)++;
                gr_idx.at(gr_found).insert(it);
                gr_idx.at(gr_found).insert(jt);
            }
        }}
        
        for (int igr = 0; igr < gr_idx.size(); ++igr) {
            if (gr_cnt.at(igr) < MIN_CONTACT_IN_PMT) continue;
            groups.insert(gr_idx.at(igr));
        }
    }
    if (groups.size() == 0) return std::vector<ChStone>();

    // search group in near hits
    short search_iter = 0;
    bool is_reached_stop_condition = false;
    while (!is_reached_stop_condition) {
        bool is_updated = false;
        std::set<std::set<int>> new_groups;
        for (auto&& group : groups) {
            // central location
            double sum_lx = 0.0, sum_ly = 0.0;
            for (auto&& ii : group) {
                sum_lx += args_hits.at(ii).lx();
                sum_ly += args_hits.at(ii).ly();
            }
            double lx = sum_lx / static_cast<double>(group.size());
            double ly = sum_ly / static_cast<double>(group.size());
            double limit = WIDTH_CELL * std::sqrt(2.0 * static_cast<double>(group.size()));

            std::set<int> new_group = group;
            for (int it = 0; it < args_hits.size(); ++it) {
                const ChHit& hit = args_hits.at(it);
                if (group.find(it) != group.end()) continue;
                double dist = std::hypot(hit.lx() - lx, hit.ly() - ly);
                if (dist > limit) continue;
                for (auto&& idx : group) {
                    const ChHit& gr_hit = args_hits.at(idx);
                    double min_dist = std::hypot(hit.lx() - gr_hit.lx(), hit.ly() - gr_hit.ly());
                    if (min_dist > MIN_CONTACT_WIDTH) continue;
                    new_group.insert(it);
                    is_updated = true;
                    break;
                }
            }
            new_groups.insert(new_group);
        }
        groups = new_groups;

        search_iter++;
        is_reached_stop_condition = (search_iter > LMTMAX_ITER || !is_updated);
    }

    // check stone condition
    std::vector<ChStone> stones;
    for (auto&& group : groups) {
        std::vector<ChHit> hits;
        for (auto&& idx : group) {
            hits.push_back(args_hits.at(idx));
            hits.back().set_cluster(ChHit::Cluster::stone);
        }
        
        std::map<int, int> pmts;
        for (auto&& hit : hits) {
            if (pmts.find(hit.pmtid()) == pmts.end()) pmts[hit.pmtid()] = 1;
            else                                      pmts[hit.pmtid()]++;
        }

        int max_nhit_in_pmt = 0;
        for (auto&& pmt : pmts) max_nhit_in_pmt = std::max(max_nhit_in_pmt, pmt.second);

        // central location
        double sum_lx = 0.0, sum_ly = 0.0;
        for (auto&& hit : hits) {
            sum_lx += hit.lx();
            sum_ly += hit.ly();
        }
        double lx = sum_lx / static_cast<double>(hits.size());
        double ly = sum_ly / static_cast<double>(hits.size());
        double dist_to_pmtc = cal_dist_to_pmtc(lx, ly);
        
        bool is_reached_stone_condition = 
            (max_nhit_in_pmt >= LMTMIN_STONE_HITS_L) && 
            ((is_within_pmtc(lx, ly) && hits.size() >= LMTMIN_STONE_HITS_L) || (hits.size() >= LMTMIN_STONE_HITS_H));
        if (!is_reached_stone_condition) continue;

        double npe = 0.0;
        for (auto&& hit : hits) npe += hit.npe();

        double chi = 0.0;
        for (auto&& hit : hits) {
            double res = std::hypot(hit.lx() - lx, hit.ly() - ly) / WIDTH_CELL;
            chi += (res * res);
        }
        double nchi = chi / static_cast<double>(hits.size() - 1);

        ChStone stone(hits, hits.size(), pmts.size(), lx, ly, npe, nchi, dist_to_pmtc);
        stones.push_back(stone);
    }
    if (stones.size() == 0) return std::vector<ChStone>();
    std::sort(stones.begin(), stones.end(), ChStone_sort());

    return stones;
}

std::vector<ChCloud> ChFit::build_cloud(const std::vector<ChHit>& args_hits) {
    if (args_hits.size() == 0) return std::vector<ChCloud>();
    std::vector<ChHit> all_hits = args_hits;
    const double MIN_CONTACT_DIST = 2.8284280; // 2 sigma of two hit = 2*sqrt(2)
   
    // search group by seeds
    std::set<int> seeds;
    for (int it = 0; it < all_hits.size(); ++it) {
        ChHit& hit = all_hits.at(it);
        hit.set_lmt_bta(lmtl_bta_, lmtu_bta_);
        if (hit.type() == 0) continue;
        if (hit.hasDb() ) seeds.insert(it*3+0);
        if (hit.hasRbA()) seeds.insert(it*3+1);
        if (hit.hasRbB()) seeds.insert(it*3+2);
    }
    
    std::vector<std::pair<std::set<int>, double>> groups;
    for (auto&& seed : seeds) {
        int seed_idx  = (seed / 3);
        int seed_type = (seed % 3);
        ChHit& seed_hit = all_hits.at(seed_idx);

        double seed_bta = 1.0;
        if      (seed_type == 0) seed_bta = seed_hit.dbta();
        else if (seed_type == 1) seed_bta = seed_hit.rbtaA();
        else if (seed_type == 2) seed_bta = seed_hit.rbtaB();

        double cand_bta = seed_bta;
        std::set<int> cand_group;
        cand_group.insert(seed_idx);
        double bdl_bta = seed_bta;
        double bdu_bta = seed_bta;

        short iter_of_seed = 0;
        bool is_reached_stop_condition_of_seed = false;
        while (!is_reached_stop_condition_of_seed) {
            bool is_updated = false;
            for (int idx = 0; idx < all_hits.size(); ++idx) {
                if (cand_group.find(idx) != cand_group.end()) continue;

                double hit_bta = -1;
                ChHit& hit = all_hits.at(idx);
                std::array<bool, 3> hit_triggers({ false, false, false });
                if (hit.hasDb()  && ((hit.dbta()  - bdl_bta) / width_bta_) > -MIN_CONTACT_DIST && ((hit.dbta()  - bdu_bta) / width_bta_) < MIN_CONTACT_DIST) hit_triggers[0] = true;
                if (hit.hasRbA() && ((hit.rbtaA() - bdl_bta) / width_bta_) > -MIN_CONTACT_DIST && ((hit.rbtaA() - bdu_bta) / width_bta_) < MIN_CONTACT_DIST) hit_triggers[1] = true;
                if (hit.hasRbB() && ((hit.rbtaB() - bdl_bta) / width_bta_) > -MIN_CONTACT_DIST && ((hit.rbtaB() - bdu_bta) / width_bta_) < MIN_CONTACT_DIST) hit_triggers[2] = true;
                bool hit_trigger = (hit_triggers[0] || hit_triggers[1] || hit_triggers[2]);
                if (!hit_trigger) continue;
                is_updated = true;

                cand_group.insert(idx);
                if (hit_triggers[0]) { bdl_bta = std::min(bdl_bta, hit.dbta() ); bdu_bta = std::max(bdu_bta, hit.dbta() ); }
                if (hit_triggers[1]) { bdl_bta = std::min(bdl_bta, hit.rbtaA()); bdu_bta = std::max(bdu_bta, hit.rbtaA()); }
                if (hit_triggers[2]) { bdl_bta = std::min(bdl_bta, hit.rbtaB()); bdu_bta = std::max(bdu_bta, hit.rbtaB()); }
            }
            iter_of_seed++;
            is_reached_stop_condition_of_seed = (iter_of_seed > LMTMAX_ITER || !is_updated);
        }
        if (cand_group.size() < LMTMIN_CLOUD_HITS) continue;

        for (int iter_of_smooth = 1; iter_of_smooth <= SMOOTH_NUM; ++iter_of_smooth) {
            short iter_of_seed_bta = 0;
            bool is_reached_stop_condition_of_seed_bta = false;
            while (!is_reached_stop_condition_of_seed_bta) {
                double sum_wgt = 0.0;
                double sum_bta = 0.0;
                for (auto&& idx : cand_group) {
                    ChHit& hit = all_hits.at(idx);
                    double hit_bta = hit.search_closest_beta(cand_bta);
                    double dlt_bta = std::abs(hit_bta - cand_bta) / width_bta_;
                    double hit_wgt = (dlt_bta <= SMOOTH_WIDTH[iter_of_smooth]) ? 1.0 : std::exp(-0.5 * (dlt_bta - SMOOTH_WIDTH[iter_of_smooth]) * (dlt_bta - SMOOTH_WIDTH[iter_of_smooth]));

                    sum_wgt += hit_wgt;
                    sum_bta += hit_wgt * hit_bta;
                }
                double bta = sum_bta / sum_wgt;
                double dlt = std::abs(bta - cand_bta) / width_bta_;
                cand_bta = bta;

                iter_of_seed_bta++;
                is_reached_stop_condition_of_seed_bta = (iter_of_seed_bta > LMTMAX_ITER || (dlt < 1.0e-3));
            }
        }

        std::set<int> rebuild_cand_group;
        for (auto&& idx : cand_group) {
            ChHit& hit = all_hits.at(idx);
            double hit_bta = hit.search_closest_beta(cand_bta);
            double dlt_bta = std::abs(hit_bta - cand_bta) / width_bta_;
            if (dlt_bta > (SMOOTH_BOUND + 2.0)) continue;
            rebuild_cand_group.insert(idx);
        }
        cand_group = rebuild_cand_group;
        
        if (cand_group.size() < LMTMIN_CLOUD_HITS) continue;
        if (!std::isfinite(cand_bta)) continue;

        std::pair<std::set<int>, double> group = std::make_pair(cand_group, cand_bta);

        bool is_searched_same_group = false;
        for (auto&& test_group : groups) {
            if (group.first.size() != test_group.first.size()) continue;
            double dlt_bta = std::abs((group.second - test_group.second) / width_bta_);
            if (dlt_bta > 5.0e-3) continue;
            is_searched_same_group = true;
            break;
        }
        if (is_searched_same_group) continue;
        
        bool is_searched_sub_group = false;
        for (auto&& test_group : groups) {
            if (group.first.size() > test_group.first.size()) continue;
            std::vector<int> overlap_sec(group.first.size());
            auto&& overlap_sec_iter  = std::set_intersection(
                group.first.begin(), group.first.end(), 
                test_group.first.begin(), test_group.first.end(), 
                overlap_sec.begin());
            overlap_sec.resize(overlap_sec_iter - overlap_sec.begin());
            if (overlap_sec.size() != group.first.size()) continue;
            is_searched_sub_group = true;
            break;
        }
        if (is_searched_sub_group) continue;
        
        groups.push_back(group);
    }
    if (groups.size() == 0) return std::vector<ChCloud>();

    // check cloud condition
    std::vector<ChCloud> clouds;
    for (auto&& group : groups) {
        std::vector<ChHit> hits;
        for (auto&& idx : group.first) {
            ChHit& hit = all_hits.at(idx);
            if (hit.type() == 0) continue;
            hit.search_closest_beta(group.second);
            hit.set_cluster(ChHit::Cluster::cloud);
            hits.push_back(hit);
        }
        
        std::map<int, int> pmts;
        for (auto&& hit : hits) {
            if (pmts.find(hit.pmtid()) == pmts.end()) pmts[hit.pmtid()] = 1;
            else                                      pmts[hit.pmtid()]++;
        }

        double beta = group.second;
        double cbta = group.second * bta_crr_;
        for (int iter = 1; iter <= 3; ++iter) {
            double sum_wgt = 0.0;
            double sum_bta = 0.0;
            for (auto&& hit : hits) {
                double hit_bta = hit.search_closest_beta(beta);
                double dlt_bta = std::abs(hit_bta - beta) / width_bta_;
                double hit_wgt = (dlt_bta <= SMOOTH_BOUND) ? 1.0 : std::exp(-0.5 * (dlt_bta - SMOOTH_BOUND) * (dlt_bta - SMOOTH_BOUND));
                sum_wgt += hit_wgt;
                sum_bta += hit_wgt * hit_bta;
            }
            beta = sum_bta / sum_wgt;
            cbta = beta * bta_crr_;
        }
        if (!std::isfinite(beta)) continue;
        if (!std::isfinite(cbta)) continue;
        
        bool is_reached_cloud_condition = (hits.size() >= LMTMIN_CLOUD_HITS && pmts.size() >= LMTMIN_CLOUD_PMTS);
        if (!is_reached_cloud_condition) continue;
        
        double npe = 0.0;
        for (auto&& hit : hits) npe += hit.npe();

        double chi = 0.0;
        for (auto&& hit : hits) {
            hit.search_closest_beta(beta);
            double res = (hit.beta() - beta) / width_bta_;
            chi += (res * res);
        }
        double nchi = chi / static_cast<double>(hits.size() - 1);

        int nhit_dir = 0;
        int nhit_rfl = 0;
        int nhit_ght = 0;
        for (auto&& hit : hits) {
            if (hit.mode() == 0) nhit_dir++;
            else                 nhit_rfl++;
        }

        ChCloud cloud(hits, hits.size(), pmts.size(), nhit_dir, nhit_rfl, nhit_ght, beta, cbta, npe, nchi);
        clouds.push_back(cloud);
    }
    if (clouds.size() == 0) return std::vector<ChCloud>();
    std::sort(clouds.begin(), clouds.end(), ChCloud_sort());
    
    return clouds;
}
        
std::tuple<std::vector<ChCloud>, std::vector<ChCloud>> ChFit::build_cloud_and_ghost(const std::vector<ChCloud>& args_clds) {
    std::vector<ChCloud> clouds;
    std::vector<ChCloud> ghosts;
    if ( args_clds.size() == 0) return make_tuple(clouds, ghosts);

    std::vector<std::set<int>> chann_set_all(args_clds.size());
    for (int ic = 0; ic < args_clds.size(); ++ic) {
        for (auto&& hit : args_clds.at(ic).hits()) chann_set_all.at(ic).insert(hit.chann());
    }

    std::vector<std::set<int>> chann_set_ght(args_clds.size());
    for (int ic = 0; ic < args_clds.size(); ++ic) {
    for (int jc = ic + 1; jc < args_clds.size(); ++jc) {
        std::vector<int> overlap_sec(std::min(chann_set_all.at(ic).size(), chann_set_all.at(jc).size()));
        auto&& overlap_iter = std::set_intersection(
                chann_set_all.at(ic).begin(), chann_set_all.at(ic).end(), 
                chann_set_all.at(jc).begin(), chann_set_all.at(jc).end(), 
                overlap_sec.begin());
        overlap_sec.resize(overlap_iter - overlap_sec.begin());
        for (auto&& chann : overlap_sec) {
            chann_set_ght.at(ic).insert(chann);
            chann_set_ght.at(jc).insert(chann);
        }
    }}

    for (int ic = 0; ic < args_clds.size(); ++ic) {
        bool is_cloud = ((chann_set_all.at(ic).size() - chann_set_ght.at(ic).size()) >= LMTMIN_CLOUD_ALONE_HITS);
        
        const ChCloud& cand_cld = args_clds.at(ic);
        std::vector<ChHit> hits = cand_cld.hits();
        for (auto&& hit : hits) hit.set_cluster( (is_cloud ? ChHit::Cluster::cloud : ChHit::Cluster::ghost) );
        int nhit_ght = chann_set_ght.at(ic).size();

        ChCloud cld(hits, cand_cld.nhit(), cand_cld.npmt(), cand_cld.nhit_dir(), cand_cld.nhit_rfl(), nhit_ght, cand_cld.beta(), cand_cld.cbta(), cand_cld.npe(), cand_cld.nchi());
        if (is_cloud) clouds.push_back(cld);
        else          ghosts.push_back(cld);
    }
    if (clouds.size() > 0) std::sort(clouds.begin(), clouds.end(), ChCloud_sort());
    if (ghosts.size() > 0) std::sort(ghosts.begin(), ghosts.end(), ChCloud_sort());

    return make_tuple(clouds, ghosts);
}


#endif // __ChSearch_C__
