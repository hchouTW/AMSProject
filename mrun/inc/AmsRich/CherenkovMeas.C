#ifndef __CherenkovMeas_C__
#define __CherenkovMeas_C__


#include "CherenkovMeas.h"


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
    bool is_switch = (lmtl_bta <= lmtu_bta) && (lmtl_bta > 0.0 || lmtu_bta > 0.0);
    bool sw_lmtl   = (is_switch && lmtl_bta > 0.0);
    bool sw_lmtu   = (is_switch && lmtu_bta > 0.0);

    type_ = 0;
    if (dbta_  > 0.0 && (!sw_lmtl || dbta_  > lmtl_bta) && (!sw_lmtu || dbta_  < lmtu_bta)) type_ += 1;
    if (rbtaA_ > 0.0 && (!sw_lmtl || rbtaA_ > lmtl_bta) && (!sw_lmtu || rbtaA_ < lmtu_bta)) type_ += 2;
    if (rbtaB_ > 0.0 && (!sw_lmtl || rbtaB_ > lmtl_bta) && (!sw_lmtu || rbtaB_ < lmtu_bta)) type_ += 4;
    
    search_closest_beta(1.0);
}

        
CherenkovFit::CherenkovFit(const std::vector<CherenkovHit>& args_hits, const std::array<double, 2>& pmtc, double rfr_index, double width_bta, double bta_crr) : CherenkovFit() {
    pmtc_      = pmtc;
    rfr_index_ = rfr_index;
    width_bta_ = width_bta;
    bta_crr_   = bta_crr;
    if (!check()) { clear(); return; }
    if (args_hits.size() == 0) { succ_ = true; return; }

    // build (stone cloud tumor)
    auto&& rlt = build(args_hits);
    stns_ = std::get<0>(rlt);
    clds_ = std::get<1>(rlt);
    tmrs_ = std::get<2>(rlt);
    gsts_ = std::get<3>(rlt);
    hits_ = std::get<4>(rlt);

    // build (npe)
    for (auto&& hit : hits_) {
        if (hit.cluster() == CherenkovHit::Cluster::stone) { nhit_stone_++; npe_stone_ += hit.npe(); }
        if (hit.cluster() == CherenkovHit::Cluster::cloud) { nhit_cloud_++; npe_cloud_ += hit.npe(); }
        if (hit.cluster() == CherenkovHit::Cluster::tumor) { nhit_tumor_++; npe_tumor_ += hit.npe(); }
        if (hit.cluster() == CherenkovHit::Cluster::ghost) { nhit_ghost_++; npe_ghost_ += hit.npe(); }
        if (hit.cluster() == CherenkovHit::Cluster::emery) { nhit_emery_++; npe_emery_ += hit.npe(); }
        if (hit.cluster() == CherenkovHit::Cluster::other) { nhit_other_++; npe_other_ += hit.npe(); }
        if (hit.cluster() == CherenkovHit::Cluster::other && hit.type() != 0) { nhit_other_inn_++; npe_other_inn_ += hit.npe(); }
        if (hit.cluster() == CherenkovHit::Cluster::other && hit.type() == 0) { nhit_other_out_++; npe_other_out_ += hit.npe(); }
        nhit_total_++; npe_total_ += hit.npe();
    }
    
    succ_ = true;
    if (!succ_) { clear(); return; }
}


void CherenkovFit::clear() {
    succ_ = false;

    stns_.clear();
    clds_.clear();
    tmrs_.clear();
    gsts_.clear();
    hits_.clear();

    nhit_total_ = 0;
    nhit_stone_ = 0;
    nhit_cloud_ = 0;
    nhit_tumor_ = 0;
    nhit_ghost_ = 0;
    nhit_emery_ = 0;
    nhit_other_ = 0;
    nhit_other_inn_ = 0;
    nhit_other_out_ = 0;

    npe_total_ = 0.0;
    npe_stone_ = 0.0;
    npe_cloud_ = 0.0;
    npe_tumor_ = 0.0;
    npe_ghost_ = 0.0;
    npe_emery_ = 0.0;
    npe_other_ = 0.0;
    npe_other_inn_ = 0.0;
    npe_other_out_ = 0.0;

    pmtc_      = std::array<double, 2>({ 0.0, 0.0 });
    rfr_index_ = 1.0;
    width_bta_ = 0.0; 
    bta_crr_   = 1.0; 
    lmtl_bta_  = 1.0;
    lmtu_bta_  = 1.0;
}


bool CherenkovFit::check() {
    if (rfr_index_ <= 1.0) return false;
    if (width_bta_ <= 0.0) return false;
    if (bta_crr_   <= 0.0) return false;
    lmtl_bta_ = (1.0 + 0.75 * WIDTH_CORE_COS) / rfr_index_ + width_bta_;
    lmtu_bta_ = 1.0 + width_bta_ * 6.0;
    if (lmtl_bta_ <= 0.0 || lmtu_bta_ <= 0.0) return false;
    return true;
}


std::vector<std::vector<CherenkovHit>> CherenkovFit::make_group_table(const std::vector<CherenkovHit>& args_hits, bool opt_stone, double width_stone, bool opt_cloud, double width_cloud) {
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

            double dist = INV_SQRT_TWO * std::hypot(ihit.cx() - jhit.cx(), ihit.cy() - jhit.cy()) / width_stone;
            if (opt_stone && dist < LMTMAX_GROUP_SGM) continue;

            if (opt_cloud && (ihit.hasDb()  && jhit.hasDb() ) && (std::fabs(ihit.dbta()  - jhit.dbta() ) / width_cloud) < LMTMAX_GROUP_SGM) continue;
            if (opt_cloud && (ihit.hasDb()  && jhit.hasRbA()) && (std::fabs(ihit.dbta()  - jhit.rbtaA()) / width_cloud) < LMTMAX_GROUP_SGM) continue;
            if (opt_cloud && (ihit.hasDb()  && jhit.hasRbB()) && (std::fabs(ihit.dbta()  - jhit.rbtaB()) / width_cloud) < LMTMAX_GROUP_SGM) continue;
            if (opt_cloud && (ihit.hasRbA() && jhit.hasDb() ) && (std::fabs(ihit.rbtaA() - jhit.dbta() ) / width_cloud) < LMTMAX_GROUP_SGM) continue;
            if (opt_cloud && (ihit.hasRbA() && jhit.hasRbA()) && (std::fabs(ihit.rbtaA() - jhit.rbtaA()) / width_cloud) < LMTMAX_GROUP_SGM) continue;
            if (opt_cloud && (ihit.hasRbA() && jhit.hasRbB()) && (std::fabs(ihit.rbtaA() - jhit.rbtaB()) / width_cloud) < LMTMAX_GROUP_SGM) continue;
            if (opt_cloud && (ihit.hasRbB() && jhit.hasDb() ) && (std::fabs(ihit.rbtaB() - jhit.dbta() ) / width_cloud) < LMTMAX_GROUP_SGM) continue;
            if (opt_cloud && (ihit.hasRbB() && jhit.hasRbA()) && (std::fabs(ihit.rbtaB() - jhit.rbtaA()) / width_cloud) < LMTMAX_GROUP_SGM) continue;
            if (opt_cloud && (ihit.hasRbB() && jhit.hasRbB()) && (std::fabs(ihit.rbtaB() - jhit.rbtaB()) / width_cloud) < LMTMAX_GROUP_SGM) continue;

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


std::tuple<std::vector<CherenkovStone>, std::vector<CherenkovCloud>, std::vector<CherenkovTumor>, std::vector<CherenkovGhost>, std::vector<CherenkovHit>> CherenkovFit::build(const std::vector<CherenkovHit>& args_hits) {
    std::tuple<std::vector<CherenkovStone>, std::vector<CherenkovCloud>, std::vector<CherenkovTumor>, std::vector<CherenkovGhost>, std::vector<CherenkovHit>> result;
    std::vector<CherenkovStone>& rltstns = std::get<0>(result);
    std::vector<CherenkovCloud>& rltclds = std::get<1>(result);
    std::vector<CherenkovTumor>& rlttmrs = std::get<2>(result);
    std::vector<CherenkovGhost>& rltgsts = std::get<3>(result);
    std::vector<CherenkovHit>&   rlthits = std::get<4>(result);
    if (args_hits.size() == 0) return result; 

    // stone
    rltstns = std::move(build_stone(args_hits, true));
    std::vector<CherenkovHit> remainder_hits = args_hits;
    
    while (rltstns.size() != 0) {
        remainder_hits.clear();
        for (auto&& hit : args_hits) {
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
        
        hit.search_closest_beta(1.0);
        if (hit.mode() < 0) continue;

        if (is_within_pmtc(hit.cx(), hit.cy())) continue;
        if (!is_within_detectable(hit.cx(), hit.cy())) continue;

        cld_hits.push_back(hit);
    }

    rltclds = std::move(build_cloud(cld_hits, false));

    // find tumor hits
    std::vector<CherenkovHit> tumor_hits;
    for (auto&& hit : args_hits) {
        bool is_used_hit = false;
        
        for (auto&& stn : rltstns) {
            for (auto&& stn_hit : stn.hits()) {
                if (stn_hit.chann() != hit.chann()) continue;
                is_used_hit = true;
                break;
            }
            if (is_used_hit) break;
        }
        if (is_used_hit) continue;
        
        for (auto&& cld : rltclds) {
            for (auto&& cld_hit : cld.hits()) {
                if (cld_hit.chann() != hit.chann()) continue;
                is_used_hit = true;
                break;
            }
            if (is_used_hit) break;
        }
        if (is_used_hit) continue;
        
        if (is_within_pmtc(hit.cx(), hit.cy())) continue;

        tumor_hits.push_back(hit);
    }
    if (tumor_hits.size() > 1) std::sort(tumor_hits.begin(), tumor_hits.end(), CherenkovHit_sort());
    
    rlttmrs = std::move(build_tumor(tumor_hits, false));
   
    // find ghost hits
    std::vector<CherenkovHit> ghost_hits;
    for (auto&& hit : args_hits) {
        bool is_used_hit = false;
        
        for (auto&& stn : rltstns) {
            for (auto&& stn_hit : stn.hits()) {
                if (stn_hit.chann() != hit.chann()) continue;
                is_used_hit = true;
                break;
            }
            if (is_used_hit) break;
        }
        if (is_used_hit) continue;
        
        for (auto&& cld : rltclds) {
            for (auto&& cld_hit : cld.hits()) {
                if (cld_hit.chann() != hit.chann()) continue;
                is_used_hit = true;
                break;
            }
            if (is_used_hit) break;
        }
        if (is_used_hit) continue;
        
        for (auto&& tmr : rlttmrs) {
            for (auto&& tmr_hit : tmr.hits()) {
                if (tmr_hit.chann() != hit.chann()) continue;
                is_used_hit = true;
                break;
            }
            if (is_used_hit) break;
        }
        if (is_used_hit) continue;
        
        if (is_within_pmtc(hit.cx(), hit.cy())) continue;

        ghost_hits.push_back(hit);
    }
    if (ghost_hits.size() > 1) std::sort(ghost_hits.begin(), ghost_hits.end(), CherenkovHit_sort());
    
    rltgsts = std::move(build_ghost(ghost_hits, false));
   
    // result
    rlthits.clear();
    for (auto&& stn : rltstns) {
        for (auto&& hit : stn.hits()) {
            rlthits.push_back(hit);
        }
    }
    for (auto&& cld : rltclds) {
        for (auto&& hit : cld.hits()) {
            rlthits.push_back(hit);
        }
    }
    for (auto&& tmr : rlttmrs) {
        for (auto&& hit : tmr.hits()) {
            rlthits.push_back(hit);
        }
    }
    for (auto&& gst : rltgsts) {
        for (auto&& hit : gst.hits()) {
            rlthits.push_back(hit);
        }
    }

    // other
    std::vector<CherenkovHit> othhits;
    for (auto&& hit : args_hits) {
        bool is_found = false;
        for (auto&& rlthit : rlthits) {
            if (hit.chann() != rlthit.chann()) continue;
            is_found = true;
            break;
        }
        if (is_found) continue;
        
        othhits.push_back(hit);
    }
    for (auto&& hit : othhits) {
        hit.set_wgt(1.0);
        hit.set_cnt(1.0);
        hit.set_lmt_bta(lmtl_bta_, lmtu_bta_);
        
        if (is_within_pmtc(hit.cx(), hit.cy())) { hit.set_cluster(CherenkovHit::Cluster::emery); }
        else                                    { hit.set_cluster(CherenkovHit::Cluster::other); }
        
        if (rltclds.size() == 0) continue;
        if (hit.type() == 0) continue;
        double beta = 1.0;

        std::vector<std::array<double, 2>> dltset;
        for (auto&& cld : rltclds) {
            double dlt = std::fabs(hit.search_closest_beta(cld.beta()) - cld.beta());
            dltset.push_back(std::array<double, 2>({ dlt, cld.beta() }));
        }
        if (dltset.size() > 1) std::sort(dltset.begin(), dltset.end());
        if (dltset.size() == 0) continue;

        hit.search_closest_beta(dltset.at(0)[1]);
    }
    for (auto&& hit : othhits) rlthits.push_back(hit);

    return result;
}


std::vector<CherenkovStone> CherenkovFit::build_stone(const std::vector<CherenkovHit>& args_hits, bool weighted) {
    std::vector<CherenkovStone> stns;
    if (args_hits.size() == 0) return stns;
    
    std::vector<CherenkovHit> hits = args_hits;
    for (auto&& hit : hits) { hit.set_wgt(1.0); hit.set_cnt(1.0); }
    if (weighted) {
        double sumnpe = 0.0;
        for (auto&& hit : hits) sumnpe += hit.npe();
        for (auto&& hit : hits) hit.set_wgt(static_cast<double>(hits.size()) * (hit.npe() / sumnpe));
    }

    auto&& gtable = make_group_table(hits, true, WIDTH_CORE_COO, false, width_bta_);
    for (auto table : gtable) {
        auto&& cstns = fit_stone(table);
        for (auto&& cstn : cstns) {
            if (!cstn.status()) continue;
            stns.push_back(cstn);
        }
    }
    if (stns.size() > 1) std::sort(stns.begin(), stns.end(), CherenkovStone_sort());

    return stns;
}
        

std::vector<CherenkovStone> CherenkovFit::fit_stone(std::vector<CherenkovHit>& hits) {
    std::vector<CherenkovStone> stns;
    if (hits.size() == 0) return stns;
    
    std::vector<std::array<double, 3>>&& cand_cstns = clustering_stone(hits); // (cx cy wgt)
    if (cand_cstns.size() == 0) return stns;

    std::vector<std::vector<CherenkovHit>> cand_chits(cand_cstns.size(), std::vector<CherenkovHit>());
    for (auto&& hit : hits) {
        std::vector<double> elms;
        double sumprb = 0.0;
        for (auto&& stn : cand_cstns) {
            double dcx = (hit.cx() - stn[0]) / WIDTH_COO;
            double dcy = (hit.cy() - stn[1]) / WIDTH_COO;
            double nrm = INV_SQRT_TWO * std::hypot(dcx,  dcy);
            double pdf = std::exp(-1.0 * nrm * nrm);
            double prb = stn[2] * pdf;
            
            if (prb * hits.size() < CONVG_PROB_SGM70) {
                elms.push_back(0.0);
                continue; 
            }
            elms.push_back(prb);
            sumprb += prb;
        }
        if (sumprb <= 0.0) continue;
        for (auto&& elm : elms) elm = (elm / sumprb);
        
        for (int it = 0; it < cand_chits.size(); ++it) {
            if (elms.at(it) < CONVG_PROB_SGM50) continue;
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
        
        double weight = 0.0;
        for (auto&& hit : cand_chit) weight += hit.wgt();
        if (weight == 0.0) continue;

        for (auto&& hit : cand_chit) hit.set_cnt(static_cast<double>(cand_chit.size()) * (hit.wgt() / weight));
        
        short nhit = cand_chit.size();
        short npmt = pmt_maps.size();

        double cnt = 0.;
        double npe = 0.;
        double chisq = 0.;
        double chisqc = 0.;
        for (auto&& hit : cand_chit) {
            cnt += hit.cnt();
            npe += hit.npe();
            double dcx = (hit.cx() - cand_cstn[0]) / WIDTH_COO;
            double dcy = (hit.cy() - cand_cstn[1]) / WIDTH_COO;
            double nrm = INV_SQRT_TWO * std::hypot(dcx,  dcy);
            chisq += hit.cnt() * (nrm * nrm); 
            
            double dcxc = (hit.cx() - cand_cstn[0]) / WIDTH_CELL;
            double dcyc = (hit.cy() - cand_cstn[1]) / WIDTH_CELL;
            double nrmc = INV_SQRT_TWO * std::hypot(dcxc,  dcyc);
            chisqc += (nrmc * nrmc);
        }

        double ndof = (cnt - 1.0);
        double nchi = chisq / ndof;
        double chic = chisqc / static_cast<double>(cand_chit.size() - 1);
        if (!std::isfinite(ndof) || ndof <= 0.0) continue;
        if (!std::isfinite(nchi)) continue;
        if (!std::isfinite(chic)) continue;

        double dist = cal_dist_to_pmtc(cand_cstn[0], cand_cstn[1]);

        for (auto&& hit : cand_chit) hit.set_cluster(CherenkovHit::Cluster::stone);

        stns.push_back(CherenkovStone(cand_chit, nhit, npmt, cand_cstn[0], cand_cstn[1], npe, cnt, nchi, chic, dist));
    }
    if (stns.size() > 1) std::sort(stns.begin(), stns.end(), CherenkovStone_sort());

    return stns;
}


std::vector<std::array<double, 3>> CherenkovFit::clustering_stone(std::vector<CherenkovHit>& hits) {
   std::vector<std::array<double, 3>> cstns; // (cx cy wgt)
    if (hits.size() == 0) return cstns;
    
    double weight = 0.0;
    for (auto&& hit : hits) weight += hit.wgt();
    if (weight == 0.0) return cstns;
    
    for (auto&& hit : hits) hit.set_cnt(static_cast<double>(hits.size()) * (hit.wgt() / weight));
    for (auto&& hit : hits) cstns.push_back(std::array<double, 3>({ hit.cx(), hit.cy(), (hit.wgt() / weight) }));
    if (cstns.size() == 0) return cstns;

    bool succ = false;
    for (int iter = 1; iter <= LMTMAX_ITER; ++iter) {
        std::vector<std::array<double, 3>> mgaus(cstns.size(), std::array<double, 3>({ 0.0, 0.0, 0.0 })); // (cx cy wgt)
        double count = 0.0;

        for (auto&& hit : hits) {
            double sumprb = 0.0;
            std::vector<std::array<double, 3>> elms; // (cx cy wgt)
            for (auto&& stn : cstns) {
                if (stn[2] <= 0.0) {
                    elms.push_back(std::array<double, 3>({ 0.0, 0.0, 0.0 }));
                    continue;
                }
                double nrmx = (hit.cx() - stn[0]) / WIDTH_COO;
                double nrmy = (hit.cy() - stn[1]) / WIDTH_COO;
                double nrm  = INV_SQRT_TWO * std::hypot(nrmx, nrmy);
                double pdf  = std::exp(-1.0 * nrm * nrm);
                double prb  = stn[2] * pdf;
                
                elms.push_back(std::array<double, 3>({ hit.cx(), hit.cy(), prb }));
                sumprb += prb;
            }
            if (!std::isfinite(sumprb) || sumprb <= 0.0) continue;
            for (auto&& elm : elms) elm.at(2) = hit.cnt() * (elm.at(2) / sumprb);
            count += hit.cnt();

            for (int it = 0; it < elms.size(); ++it) {
                std::array<double, 3>& elm  = elms.at(it);
                std::array<double, 3>& gaus = mgaus.at(it);
                if (elm[2] <= 0.0) continue;
                gaus[0] += elm[2] * elm[0];
                gaus[1] += elm[2] * elm[1];
                gaus[2] += elm[2];
            }
        }
        if (!std::isfinite(count) || count <= 1.0) break;
        
        for (auto&& gaus : mgaus) {
            if (gaus[2] <= 0.0) continue;
            gaus[0] = gaus[0] / gaus[2];
            gaus[1] = gaus[1] / gaus[2];
            gaus[2] = gaus[2] / count;
        }

        // remove empty
        if (iter >= LMTL_ITER) {
            for (int it = mgaus.size()-1; it >= 0; --it) {
                if (mgaus.at(it)[2] < CONVG_EPSILON)
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
                    double dnxA = (newx - gausA[0]) / WIDTH_COO;
                    double dnyA = (newy - gausA[1]) / WIDTH_COO;
                    double dnxB = (newx - gausB[0]) / WIDTH_COO;
                    double dnyB = (newy - gausB[1]) / WIDTH_COO;
                    double dnA  = INV_SQRT_TWO * std::sqrt(dnxA * dnxA + dnyA * dnyA);
                    double dnB  = INV_SQRT_TWO * std::sqrt(dnxB * dnxB + dnyB * dnyB);
                    double prb  = gausA[2] * std::exp(-1.0 * dnA * dnA) + 
                                  gausB[2] * std::exp(-1.0 * dnB * dnB);
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
                    double newcx  = 0.0;
                    double newcy  = 0.0;
                    double newwgt = 0.0;
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

        double normalized = 0.0;
        for (auto&& gaus : mgaus) normalized += gaus[2];
        for (auto&& gaus : mgaus) gaus[2] = gaus[2] / normalized;
        if (normalized <= 0.0) break;

        if (iter >= LMTU_ITER && !has_remove_empty && !has_merge) {
            bool is_convg = true;
            for (int it = 0; it < mgaus.size(); ++it) {
                auto&& gaus = mgaus.at(it);
                auto&& stn  = cstns.at(it);
                bool match = ((std::abs(gaus[0] - stn[0]) <= CONVG_TOLERANCE) && 
                              (std::abs(gaus[1] - stn[1]) <= CONVG_TOLERANCE) && 
                              (std::abs(gaus[2] - stn[2]) <= CONVG_TOLERANCE));
                if (!match) { is_convg = false; break; }
            }
            if (is_convg) { succ = true; break; }
        }
        cstns = mgaus;
    }
    
    if (!succ) cstns.clear();
    return cstns;
}


std::vector<CherenkovCloud> CherenkovFit::build_cloud(const std::vector<CherenkovHit>& args_hits, bool weighted) {
    std::vector<CherenkovCloud> clds;
    if (args_hits.size() == 0) return clds;
    
    std::vector<CherenkovHit> hits = args_hits;
    for (auto&& hit : hits) hit.set_lmt_bta(lmtl_bta_, lmtu_bta_);

    for (auto&& hit : hits) { hit.set_wgt(1.0); hit.set_cnt(1.0); }
    if (weighted) {
        double sumnpe = 0.0;
        for (auto&& hit : hits) sumnpe += hit.npe();
        for (auto&& hit : hits) hit.set_wgt(static_cast<double>(hits.size()) * (hit.npe() / sumnpe));
    }

    auto&& gtable = make_group_table(hits, false, WIDTH_CORE_COO, true, width_bta_);
    for (auto table : gtable) {
        auto&& cclds = fit_cloud(table);
        for (auto&& ccld : cclds) {
            if (!ccld.status()) continue;
            CherenkovCloud&& cld = refit_cloud(ccld);
            if (!cld.status()) continue;
            clds.push_back(cld);
        }
    }
    if (clds.size() > 1) std::sort(clds.begin(), clds.end(), CherenkovCloud_sort());

    // umbra (misjudge)
    for (auto&& cld : clds) {
        std::vector<CherenkovHit> ubr_hits;
        for (auto&& hit : cld.hits()) {
            bool hasDb  = (hit.hasDb()  && hit.mode() != 0);
            bool hasRbA = (hit.hasRbA() && hit.mode() != 1);
            bool hasRbB = (hit.hasRbB() && hit.mode() != 2);
            if (!hasDb && !hasRbA && !hasRbB) continue;
            ubr_hits.push_back(
                CherenkovHit(hit.chann(), hit.pmtid(),
                             (hasDb  ? hit.dbta()  : -1.0),
                             (hasRbA ? hit.rbtaA() : -1.0),
                             (hasRbB ? hit.rbtaB() : -1.0), 
                             hit.npe(), hit.cx(), hit.cy())
            );
        }
        if (ubr_hits.size() == 0) continue;
        
        for (auto&& hit : ubr_hits) hit.set_lmt_bta(lmtl_bta_, lmtu_bta_);
        
        for (auto&& hit : ubr_hits) { hit.set_wgt(1.0); hit.set_cnt(1.0); }
        if (weighted) {
            double sumnpe = 0.0;
            for (auto&& hit : ubr_hits) sumnpe += hit.npe();
            for (auto&& hit : ubr_hits) hit.set_wgt(static_cast<double>(ubr_hits.size()) * (hit.npe() / sumnpe));
        }
        
        auto&& ubr_gtable = make_group_table(ubr_hits, false, WIDTH_CORE_COO, true, width_bta_);
        for (auto table : ubr_gtable) {
            auto&& cubrs = fit_cloud(table);
            for (auto&& cubr : cubrs) {
                if (!cubr.status()) continue;
                CherenkovCloud&& ubr = refit_cloud(cubr);
                if (!ubr.status()) continue;

                double wgtgbl = static_cast<double>(ubr.hits().size()) / static_cast<double>(cld.hits().size());
                double wgtloc = static_cast<double>(ubr.hits().size()) / static_cast<double>(ubr_hits.size());
                if (wgtgbl <= (2.0 / 3.0)) continue;
                if (wgtloc <= (2.0 / 3.0)) continue;
                
                double misjudge = std::fabs(cld.cbta() - ubr.cbta()) / width_bta_;
                if (misjudge < cld.misjudge()) continue;
                cld.set_misjudge(misjudge);
            }
        }
    }

    return clds;
}
    

CherenkovCloud CherenkovFit::refit_cloud(const CherenkovCloud& cand_cld) {
    if (!cand_cld.status()) return CherenkovCloud();
    std::vector<CherenkovHit> cand_chit = cand_cld.hits();
    if (cand_chit.size() == 0) return CherenkovCloud();

    double weight = 0.0;
    for (auto&& hit : cand_chit) weight += hit.wgt();
    if (weight == 0.0) return CherenkovCloud();
    for (auto&& hit : cand_chit) hit.set_cnt(static_cast<double>(cand_chit.size()) * (hit.wgt() / weight));

    double cand_beta = cand_cld.beta();
    double cand_cnt  = 0.0;
    double cand_nchi = 0.0;

    // refit cloud
    bool succ = false;
    for (int iter = 0; iter <= LMTMAX_ITER; ++iter) {
        double count = 0.0;
        double chisq = 0.0;
        double grdB  = 0.0;
        double cvBB  = 0.0;
        
        for (auto&& hit : cand_chit) {
            double dlt = hit.search_closest_beta(cand_beta) - cand_beta;
            double nrm = dlt / width_bta_;

            double smooth = SMOOTH_BOUND[1] - static_cast<double>(iter) * SMOOTH_RATE;
            if (smooth < SMOOTH_BOUND[0]) smooth = SMOOTH_BOUND[0];
            
            double smnrm = std::fabs(nrm) - smooth;
            double smprb = (smnrm <= 0.0) ? 1.0 : std::exp(-0.5 * smnrm * smnrm);

            count += hit.cnt();
            chisq += hit.cnt() * (dlt / width_bta_) * (dlt / width_bta_);
            grdB  += hit.cnt() * smprb * (-1.0 * dlt / width_bta_ / width_bta_);
            cvBB  += hit.cnt() * smprb * (1.0 / width_bta_ / width_bta_);
        }
        if (!std::isfinite(count) || count <= 1.0) break;

        double nchi = chisq / (count - 1.0);
        if (!std::isfinite(nchi)) break;

        double resB = grdB / cvBB;
        double beta = cand_beta - resB;
        if (!std::isfinite(resB) || beta <= 0.0) break;
       
        if (iter >= LMTM_ITER) {
            succ = (std::abs(cand_beta - beta) <= CONVG_TOLERANCE) && 
                   (std::abs(cand_cnt - count) <= CONVG_TOLERANCE) &&
                   (std::abs(cand_nchi - nchi) <= CONVG_TOLERANCE);
            if (succ) break;
        }

        cand_beta = beta;
        cand_cnt  = count;
        cand_nchi = nchi;
    }
    if (!succ) return CherenkovCloud();
  
    // check result
    std::vector<CherenkovHit> cand_reduce_chit;
    for (auto&& hit : cand_chit) {
        double absnrm = std::fabs(hit.search_closest_beta(cand_beta) - cand_beta) / width_bta_;
        if (absnrm > 5.0) continue;
        cand_reduce_chit.push_back(hit);
    } 
    if (cand_reduce_chit.size() < LMTMIN_CLOUD_HITS) return CherenkovCloud();
    std::sort(cand_reduce_chit.begin(), cand_reduce_chit.end(), CherenkovHit_sort());
    
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
    double cand_nhit_npmt = static_cast<double>(cand_nhit) / static_cast<double>(cand_npmt);

    short cand_nhit_dir = 0;
    short cand_nhit_rfl = 0;
    for (auto&& hit : cand_reduce_chit) {
        if (hit.mode() == 0) cand_nhit_dir++;
        else                 cand_nhit_rfl++;
    }

    double cand_npe = 0.;
    for (auto&& hit : cand_reduce_chit) cand_npe += hit.npe();

    double cand_cbta = cand_beta * bta_crr_;
    double cand_ndof = cand_cnt - 1.0;
    if (!std::isfinite(cand_ndof) || cand_ndof <= 0.0) return CherenkovCloud();
    if (!std::isfinite(cand_nchi)) return CherenkovCloud();

    double cand_misjudge = 0.0; 

    for (auto&& hit : cand_reduce_chit) hit.set_cluster(CherenkovHit::Cluster::cloud);
    
    return CherenkovCloud(cand_reduce_chit, cand_nhit, cand_npmt, cand_nhit_dir, cand_nhit_rfl, cand_beta, cand_cbta, cand_npe, cand_cnt, cand_nchi, cand_misjudge);
}


std::vector<CherenkovCloud> CherenkovFit::fit_cloud(std::vector<CherenkovHit>& hits) {
    std::vector<CherenkovCloud> cclds;
    if (hits.size() == 0) return cclds;
    
    std::vector<std::array<double, 2>>&& cand_cclds = clustering_cloud(hits); // (bta wgt)
    if (cand_cclds.size() == 0) return cclds;
    
    std::vector<std::vector<CherenkovHit>> cand_chits(cand_cclds.size(), std::vector<CherenkovHit>());
    for (auto&& hit : hits) {
        std::vector<double> elms;
        double sumprb = 0.0;
        for (auto&& cld : cand_cclds) {
            double nrm = (hit.search_closest_beta(cld[0]) - cld[0]) / width_bta_;
            double pdf = std::exp(-0.5 * nrm * nrm);
            double prb = cld[1] * pdf;
            
            if (prb * hits.size() < CONVG_PROB_SGM70) {
                elms.push_back(0.0);
                continue; 
            }
            elms.push_back(prb);
            sumprb += prb;
        }
        if (sumprb <= 0.0) continue;
        for (auto&& elm : elms) elm = (elm / sumprb);
        
        for (int it = 0; it < cand_chits.size(); ++it) {
            if (elms.at(it) < CONVG_PROB_SGM30) continue;
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
        
        double weight = 0.0;
        for (auto&& hit : cand_chit) weight += hit.wgt();
        if (weight == 0.0) continue;
        for (auto&& hit : cand_chit) hit.set_cnt(static_cast<double>(cand_chit.size()) * (hit.wgt() / weight));
        
        short nhit = cand_chit.size();
        short npmt = pmt_maps.size();
        double nhit_npmt = static_cast<double>(nhit) / static_cast<double>(npmt);
        
        short nhit_dir = 0;
        short nhit_rfl = 0;
        for (auto&& hit : cand_chit) {
            if (hit.mode() == 0) nhit_dir++;
            else                 nhit_rfl++;
        }

        double beta = cand_ccld[0];
        double cbta = cand_ccld[0] * bta_crr_;

        double cnt = 0.;
        double npe = 0.;
        double chisq = 0.;
        for (auto&& hit : cand_chit) {
            cnt += hit.cnt();
            npe += hit.npe();
            
            double dlt = hit.search_closest_beta(cand_ccld[0]) - cand_ccld[0];
            chisq += hit.cnt() * (dlt / width_bta_) * (dlt / width_bta_);
        }

        double ndof = (cnt - 1.0);
        double nchi = chisq / ndof;
        if (!std::isfinite(ndof) || ndof <= 0.0) continue;
        if (!std::isfinite(nchi)) continue;
    
        double misjudge = 0.0;
    
        for (auto&& hit : cand_chit) hit.set_cluster(CherenkovHit::Cluster::cloud);
        
        cclds.push_back(CherenkovCloud(cand_chit, nhit, npmt, nhit_dir, nhit_rfl, beta, cbta, npe, cnt, nchi, misjudge));
    }
    if (cclds.size() > 1) std::sort(cclds.begin(), cclds.end(), CherenkovCloud_sort());

    return cclds;
}


std::vector<std::array<double, 2>> CherenkovFit::clustering_cloud(std::vector<CherenkovHit>& hits) {
    // bayesian information criterion (BIC)
    std::vector<std::array<double, 2>> cclds; // (bta wgt)
    if (hits.size() == 0) return cclds;
    
    double weight_hits = 0.0;
    for (auto&& hit : hits) weight_hits += hit.wgt();
    if (weight_hits == 0.0) return cclds;
    for (auto&& hit : hits) hit.set_cnt(static_cast<double>(hits.size()) * (hit.wgt() / weight_hits));

    for (auto&& hit : hits) {
        if (hit.type() == 0) continue;
        double multi = hit.cnt() / static_cast<double>(hit.hasDb() * 2.0 + hit.hasRbA() + hit.hasRbB());
        if (hit.hasDb() ) cclds.push_back(std::array<double, 2>({ hit.dbta() , 2.0 * multi }));
        if (hit.hasRbA()) cclds.push_back(std::array<double, 2>({ hit.rbtaA(), 1.0 * multi }));
        if (hit.hasRbB()) cclds.push_back(std::array<double, 2>({ hit.rbtaB(), 1.0 * multi }));
    }
    if (cclds.size() > 1) std::sort(cclds.begin(), cclds.end());
    if (cclds.size() == 0) return cclds;

    double weight = 0.0;
    for (auto&& cld : cclds) weight += cld[1];
    for (auto&& cld : cclds) cld[1] = (cld[1] / weight);
    if (weight == 0.0) { cclds.clear(); return cclds; }

    bool succ = false;
    for (int iter = 1; iter <= LMTMAX_ITER; ++iter) {
        std::vector<std::array<double, 4>> mgaus(cclds.size(), std::array<double, 4>({ 0.0, 0.0, 0.0, 0.0 })); // (wgt grd cv res)
        double count = 0.0;

        for (auto&& hit : hits) {
            double sumprb = 0.0;
            std::vector<std::array<double, 3>> elms; // (wgt grd cv)
            for (auto&& cld : cclds) {
                if (cld[1] <= 0.0) {
                    elms.push_back(std::array<double, 3>({ 0.0, 0.0, 0.0 }));
                    continue;
                }
                double dlt = hit.search_closest_beta(cld[0]) - cld[0];
                double nrm = dlt / width_bta_;
                double pdf = std::exp(-0.5 * nrm * nrm);
                double prb = cld[1] * pdf;
        
                double grdB = (-1.0 * dlt / width_bta_ / width_bta_);
                double cvBB = (1.0 / width_bta_ / width_bta_);
                
                elms.push_back(std::array<double, 3>({ prb, grdB, cvBB }));
                sumprb += prb;
            }
            if (!std::isfinite(sumprb) || sumprb <= 0.0) continue;
            for (auto&& elm : elms) elm.at(0) = hit.cnt() * (elm.at(0) / sumprb);
            count += hit.cnt();

            for (int it = 0; it < elms.size(); ++it) {
                std::array<double, 3>& elm  = elms.at(it);
                std::array<double, 4>& gaus = mgaus.at(it);
                if (elm[0] <= 0.0) continue;
                gaus[0] += elm[0];
                gaus[1] += elm[0] * elm[1];
                gaus[2] += elm[0] * elm[2];
            }
        }
        if (!std::isfinite(count) || count <= 1.0) break;

        for (auto&& gaus : mgaus) {
            if (gaus[0] <= 0.0) continue;
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
                if (new_cclds.at(it)[1] < CONVG_EPSILON)
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
                double prbx = gausA[1] * std::exp(-0.5 * dnxA * dnxA) + 
                              gausB[1] * std::exp(-0.5 * dnxB * dnxB);
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
                    double newbta = 0.0;
                    double newwgt = 0.0;
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

        double normalized = 0.0;
        for (auto&& cld : new_cclds) normalized += cld[1];
        for (auto&& cld : new_cclds) cld[1] = (cld[1] / normalized);
        if (normalized <= 0.0) break;

        if (iter >= LMTU_ITER && !has_remove_empty && !has_merge) {
            bool is_convg = true;
            for (int it = 0; it < new_cclds.size(); ++it) {
                auto&& gaus = new_cclds.at(it);
                auto&& cld  = cclds.at(it);
                bool match = ((std::abs(gaus[0] - cld[0]) <= CONVG_TOLERANCE) && 
                              (std::abs(gaus[1] - cld[1]) <= CONVG_TOLERANCE));
                if (!match) { is_convg = false; break; }
            }
            if (is_convg) { succ = true; break; }
        }
        cclds = new_cclds;
    }

    if (!succ) cclds.clear();
    return cclds;
}


std::vector<CherenkovTumor> CherenkovFit::build_tumor(const std::vector<CherenkovHit>& args_hits, bool weighted) {
    std::vector<CherenkovTumor> tmrs;
    if (args_hits.size() == 0) return tmrs;

    std::vector<CherenkovHit> hits = args_hits;
    for (auto&& hit : hits) hit.set_lmt_bta(lmtl_bta_, lmtu_bta_ + width_bta_ * 4.0);
    
    for (auto&& hit : hits) { hit.set_wgt(1.0); hit.set_cnt(1.0); }
    if (weighted) {
        double sumnpe = 0.0;
        for (auto&& hit : hits) sumnpe += hit.npe();
        for (auto&& hit : hits) hit.set_wgt(static_cast<double>(hits.size()) * (hit.npe() / sumnpe));
    }

    // tumor (cloud)
    auto&& gtable = make_group_table(hits, false, WIDTH_CORE_COO, true, width_bta_);
    for (auto&& table : gtable) {
        CherenkovTumor&& tmr = fit_tumor(table);
        if (!tmr.status()) continue;
        tmrs.push_back(tmr);
    }
    if (tmrs.size() > 1) std::sort(tmrs.begin(), tmrs.end(), CherenkovTumor_sort());

    return tmrs;
}
        

CherenkovTumor CherenkovFit::fit_tumor(std::vector<CherenkovHit>& hits) {
    if (hits.size() == 0) return CherenkovTumor();

    double weight = 0.0;
    for (auto&& hit : hits) weight += hit.wgt();
    if (weight == 0.0) return CherenkovTumor();
    for (auto&& hit : hits) hit.set_cnt(static_cast<double>(hits.size()) * (hit.wgt() / weight));

    double init_beta = 0.0;
    double count_beta = 0.0;
    for (auto&& hit : hits) {
        if (hit.type() == 0) continue;
        double multi = hit.cnt() / static_cast<double>(hit.hasDb() + hit.hasRbA() + hit.hasRbB());
        if (hit.hasDb() ) { init_beta += multi * hit.dbta() ; count_beta += multi; }
        if (hit.hasRbA()) { init_beta += multi * hit.rbtaA(); count_beta += multi; }
        if (hit.hasRbB()) { init_beta += multi * hit.rbtaB(); count_beta += multi; }
    }
    init_beta = (count_beta != 0.0 ? (init_beta / count_beta) : 0.0);

    double beta = init_beta;
    if (beta > 0.0) {
        short iter = 1;
        bool is_ok = false;
        while (!is_ok && iter <= LMTMAX_ITER) {
            double sum_bta = 0.0;
            double sum_cnt = 0.0;
            for (auto&& hit : hits) {
                if (hit.type() == 0) continue;
                sum_bta += hit.cnt() * hit.search_closest_beta(beta);
                sum_cnt += hit.cnt();
            }
            if (sum_cnt == 0.0) break;
            double new_beta = (sum_bta / sum_cnt);
            is_ok = (iter >= LMTL_ITER && std::abs(new_beta - beta) <= CONVG_TOLERANCE);
            if (!is_ok) beta = new_beta;
            iter++;
        }
        if (!is_ok) beta = 0.0;
    }

    // refit
    bool succ = false;
    double nchi = 0.0;
    for (int iter = 0; iter <= LMTMAX_ITER && (beta > 0.0); ++iter) {
        double count = 0.0;
        double chisq = 0.0;
        double grdB  = 0.0;
        double cvBB  = 0.0;

        for (auto&& hit : hits) {
            double dlt = hit.search_closest_beta(beta) - beta;
            double nrm = dlt / width_bta_;
            
            double smooth = SMOOTH_BOUND[1] - static_cast<double>(iter) * SMOOTH_RATE;
            if (smooth < SMOOTH_BOUND[0]) smooth = SMOOTH_BOUND[0];
            
            double smnrm = std::fabs(nrm) - smooth;
            double smprb = (smnrm <= 0.0) ? 1.0 : std::exp(-0.5 * smnrm * smnrm);
            
            count += hit.cnt();
            chisq += hit.cnt() * (dlt / width_bta_) * (dlt / width_bta_);
            grdB  += hit.cnt() * smprb * (-1.0 * dlt / width_bta_ / width_bta_);
            cvBB  += hit.cnt() * smprb * (1.0 / width_bta_ / width_bta_);
        }
        if (!std::isfinite(count) || count <= 1.0) break;

        double new_nchi = chisq / (count - 1.0);
        if (!std::isfinite(new_nchi)) break;

        double resB = grdB / cvBB;
        double new_beta = beta - resB;
        if (!std::isfinite(resB) || new_beta <= 0.0) break;
       
        if (iter >= LMTL_ITER) {
            succ = (std::abs(new_beta - beta) <= CONVG_TOLERANCE) && 
                   (std::abs(new_nchi - nchi) <= CONVG_TOLERANCE);
            if (succ) break;
        }

        beta = new_beta;
        nchi = new_nchi;
    }
    double cbta = bta_crr_ * beta;
    if (!succ) return CherenkovTumor();
    
    // check result
    std::vector<CherenkovHit> reduce_hits;
    for (auto&& hit : hits) {
        double absnrm = std::fabs(hit.search_closest_beta(beta) - beta) / width_bta_;
        if (absnrm > 5.0) continue;
        reduce_hits.push_back(hit);
    } 
    if (reduce_hits.size() < LMTMIN_TUMOR_HITS) return CherenkovTumor();
    std::sort(reduce_hits.begin(), reduce_hits.end(), CherenkovHit_sort());
    
    std::map<int, int> pmt_maps;
    for (auto&& hit : reduce_hits) {
        if (pmt_maps.count(hit.pmtid()) == 0) pmt_maps[hit.pmtid()] = 1;
        else pmt_maps[hit.pmtid()]++;
    }
    short nhit = reduce_hits.size();
    short npmt = pmt_maps.size();
    if (nhit < LMTMIN_TUMOR_HITS) return CherenkovTumor();
    if (npmt < LMTMIN_TUMOR_PMTS) return CherenkovTumor();
        
    double npe = 0.0;
    for (auto&& hit : reduce_hits) npe += hit.npe();

    for (auto&& hit : reduce_hits) hit.set_cluster(CherenkovHit::Cluster::tumor);

    return CherenkovTumor(reduce_hits, nhit, npmt, beta, cbta, npe);
}
        

std::vector<CherenkovGhost> CherenkovFit::build_ghost(const std::vector<CherenkovHit>& args_hits, bool weighted) {
    std::vector<CherenkovGhost> gsts;
    if (args_hits.size() == 0) return gsts;

    std::vector<CherenkovHit> hits = args_hits;
    for (auto&& hit : hits) hit.set_lmt_bta(lmtl_bta_, lmtu_bta_ + width_bta_ * 6.0);
    
    for (auto&& hit : hits) { hit.set_wgt(1.0); hit.set_cnt(1.0); }
    if (weighted) {
        double sumnpe = 0.0;
        for (auto&& hit : hits) sumnpe += hit.npe();
        for (auto&& hit : hits) hit.set_wgt(static_cast<double>(hits.size()) * (hit.npe() / sumnpe));
    }
    
    // ghost (cloud)
    auto&& gtable = make_group_table(hits, false, (1.0/3.0) * WIDTH_PMT, true, (1.5) * width_bta_);
    for (auto&& table : gtable) {
        if (table.size() < LMTMIN_GHOST_HITS) continue;
        
        short maxnhits = 0;
        auto&& gstones = make_group_table(table, true, (1.0/3.0) * WIDTH_PMT, false, (1.5) * width_bta_);
        for (auto&& gstone : gstones) maxnhits = std::max(maxnhits, static_cast<short>(gstone.size()));
        if (gstones.size() < LMTMIN_GHOST_PMTS && maxnhits <= LMTMIN_GHOST_HITS) continue;

        CherenkovGhost&& gst = fit_ghost(table);
        if (!gst.status()) continue;
        
        gsts.push_back(gst);
    }
    if (gsts.size() > 1) std::sort(gsts.begin(), gsts.end(), CherenkovGhost_sort());

    return gsts;
}


CherenkovGhost CherenkovFit::fit_ghost(std::vector<CherenkovHit>& hits) {
    if (hits.size() == 0) return CherenkovGhost();

    double weight = 0.0;
    for (auto&& hit : hits) weight += hit.wgt();
    if (weight == 0.0) return CherenkovGhost();
    for (auto&& hit : hits) hit.set_cnt(static_cast<double>(hits.size()) * (hit.wgt() / weight));

    double init_beta = 0.0;
    double count_beta = 0.0;
    for (auto&& hit : hits) {
        if (hit.type() == 0) continue;
        double multi = hit.cnt() / static_cast<double>(hit.hasDb() + hit.hasRbA() + hit.hasRbB());
        if (hit.hasDb() ) { init_beta += multi * hit.dbta() ; count_beta += multi; }
        if (hit.hasRbA()) { init_beta += multi * hit.rbtaA(); count_beta += multi; }
        if (hit.hasRbB()) { init_beta += multi * hit.rbtaB(); count_beta += multi; }
    }
    init_beta = (count_beta != 0.0 ? (init_beta / count_beta) : 0.0);

    double beta = init_beta;
    if (beta > 0.0) {
        short iter = 1;
        bool is_ok = false;
        while (!is_ok && iter <= LMTMAX_ITER) {
            double sum_bta = 0.0;
            double sum_cnt = 0.0;
            for (auto&& hit : hits) {
                if (hit.type() == 0) continue;
                sum_bta += hit.cnt() * hit.search_closest_beta(beta);
                sum_cnt += hit.cnt();
            }
            if (sum_cnt == 0.0) break;
            double new_beta = (sum_bta / sum_cnt);
            is_ok = (iter >= LMTL_ITER && std::abs(new_beta - beta) <= CONVG_TOLERANCE);
            if (!is_ok) beta = new_beta;
            iter++;
        }
        if (!is_ok) beta = 0.0;
    }
    if (beta <= 0.0) return CherenkovGhost();
    double cbta = bta_crr_ * beta;
    
    std::map<int, int> pmt_maps;
    for (auto&& hit : hits) {
        if (pmt_maps.count(hit.pmtid()) == 0) pmt_maps[hit.pmtid()] = 1;
        else pmt_maps[hit.pmtid()]++;
    }
    short nhit = hits.size();
    short npmt = pmt_maps.size();
    if (nhit < LMTMIN_GHOST_HITS) return CherenkovGhost();
    
    double npe = 0.0;
    for (auto&& hit : hits) npe += hit.npe();

    for (auto&& hit : hits) hit.set_cluster(CherenkovHit::Cluster::ghost);

    return CherenkovGhost(hits, nhit, npmt, beta, cbta, npe);
}







///////////////////////////////////////////////////////////////////////////////////
/*
CherenkovRayTrace::CherenkovRayTrace(double kind, double index, double cbta, const CherenkovCloud* cloud) : CherenkovRayTrace() {
    if (kind == 0 || index 
    bool has_cloud = (cloud != nullptr && cloud->status());
    if (!status_ || kind_ == KIND_EMPTY || (has_cloud ? cloud->cbta() : cbta) <= (1.0 / index_)) return rlt_trace;

    double cosv = 1.0 / (index_ * (has_cloud ? cloud->cbta() : cbta));
    double sinv = std::sqrt((1.0 - cosv) * (1.0 + cosv));
    if (!std::isfinite(cosv) || !std::isfinite(sinv)) return rlt_trace;
    int kind = (kind_ == KIND_AGL) ? 0 : 1;

    double lmtrfr = -std::cos(std::asin(1.0 / index_));

    SVecD<3> radp(radp_[0], radp_[1], radp_[2]);
    SVecD<3> radd(radd_[0], radd_[1], radd_[2]);
    OrthCoord orth(radd);
    
    double topx = radp[0] + (0.5 * RAD_HEIGHT[kind]) * (radd[0] / radd[2]);
    double topy = radp[1] + (0.5 * RAD_HEIGHT[kind]) * (radd[1] / radd[2]);
    double topz = radp[2] + (0.5 * RAD_HEIGHT[kind]);

    std::array<std::array<short,  TRACE_NSET>, TRACE_NPHI> border;  // 0 no,  1 pass
    std::array<std::array<short,  TRACE_NSET>, TRACE_NPHI> simphs;  // 0 no,  1 dir,  2 rfl
    std::array<std::array<double, TRACE_NSET>, TRACE_NPHI> simphsx;
    std::array<std::array<double, TRACE_NSET>, TRACE_NPHI> simphsy;
    for(auto&& phs : border) phs.fill(0);
    for(auto&& phs : simphs) phs.fill(0);
    for(auto&& phsx : simphsx) phsx.fill(0.0);
    for(auto&& phsy : simphsy) phsy.fill(0.0);

    for (int it = 0; it < TRACE_NPHI; ++it) {
    for (int is = 0; is < TRACE_NSET; ++is) {
        double locwgt = (static_cast<double>(is) + 0.5) / static_cast<double>(TRACE_NSET);
        
        double phi = TWO_PI * static_cast<double>(it) / static_cast<double>(TRACE_NPHI);
        SVecD<3>&& phivec = orth.tau() * std::cos(phi) + orth.rho() * std::sin(phi);
        SVecD<3>&& srcvec = radd * cosv + phivec * sinv;
        if (srcvec[2] >= lmtrfr) continue;

        double planex = topx + (-RAD_HEIGHT[kind]) * (locwgt * (radd[0] / radd[2]) + (1.0 - locwgt) * (srcvec[0] / srcvec[2]));
        double planey = topy + (-RAD_HEIGHT[kind]) * (locwgt * (radd[1] / radd[2]) + (1.0 - locwgt) * (srcvec[1] / srcvec[2]));
        double planez = topz + (-RAD_HEIGHT[kind]);
        if (std::hypot(planex, planey) > MIRROR_TOP_RADIUS) continue;
        
        short tile = RichRingR::getTileIndex(planex, planey);
        if (tile < 0 || tile != tile_) continue;
        border[it][is] = 1; 

        double cosvac = -1.0 * std::sqrt(1.0 - (index_ * index_) * (1.0 - srcvec[2]) * (1.0 + srcvec[2]));
        double sinvac = std::sqrt((1.0 - cosvac) * (1.0 + cosvac));
        if (!std::isfinite(cosvac) || !std::isfinite(sinvac)) continue;
        SVecD<2>&& vacr = LA::Unit(SVecD<2>(srcvec[0], srcvec[1])) * sinvac;
        if (!std::isfinite(vacr[0]) || !std::isfinite(vacr[1])) vacr = std::move(SVecD<2>());

        SVecD<3> vacp(planex, planey, planez);
        SVecD<3> vacd(vacr[0], vacr[1], cosvac);

        double dheight = (PMT_CZ - vacp[2]);
        double pmtx = vacp[0] + dheight * (vacd[0] / vacd[2]);
        double pmty = vacp[1] + dheight * (vacd[1] / vacd[2]);
        if (std::fabs(pmtx) < PMT_HOLE[0] && std::fabs(pmty) < PMT_HOLE[1]) continue;
        if (std::hypot(pmtx, pmty) < MIRROR_BTM_RADIUS) { 
            if (RichOffline::RichPMTsManager::FindChannel(pmtx, pmty) < 0) continue;
            simphs[it][is]  = 1; 
            simphsx[it][is] = pmtx; 
            simphsy[it][is] = pmty; 
            continue; 
        }

        double vactx = (vacd[0] / vacd[2]);
        double vacty = (vacd[1] / vacd[2]);
        double mirtr = ((MIRROR_TOP_RADIUS - MIRROR_BTM_RADIUS) / MIRROR_HEIGHT);
        double poly0 = (mirtr * mirtr - (vactx * vactx + vacty * vacty));
        double poly1 = 2.0 * (MIRROR_TOP_RADIUS * mirtr - (vacp[0] * vactx + vacp[1] * vacty));
        double poly2 = (MIRROR_TOP_RADIUS * MIRROR_TOP_RADIUS - (vacp[0] * vacp[0] + vacp[1] * vacp[1]));

        double det  = std::sqrt(poly1 * poly1 - 4.0 * poly0 * poly2);
        double sol1 = 0.5 * (-poly1 + det) / poly0;
        double sol2 = 0.5 * (-poly1 - det) / poly0;
        bool is_ok1 = (std::isfinite(sol1) && sol1 < 0.0 && sol1 > dheight);
        bool is_ok2 = (std::isfinite(sol2) && sol2 < 0.0 && sol2 > dheight);

        if (!is_ok1 && !is_ok2) continue;

        double sol = 0;
        if (is_ok1) sol = sol1;
        if (is_ok2) sol = sol2;
        vacp[0] = vacp[0] + sol * vactx;
        vacp[1] = vacp[1] + sol * vacty;
        vacp[2] = vacp[2] + sol;
        
        double cosmir = -1.0 * std::sqrt((mirtr * mirtr) / (1.0 + mirtr * mirtr));
        double sinmir = std::sqrt(1.0 - cosmir * cosmir);
        SVecD<2>&& mirr = LA::Unit(SVecD<2>(-vacp[0], -vacp[1])) * sinmir;
        SVecD<3> mird(mirr[0], mirr[1], cosmir);

        double cosrfl = std::fabs(LA::Dot(vacd, mird));
        if (!std::isfinite(cosrfl)) continue;
        
        vacd = std::move(LA::Unit(vacd + 2.0 * cosrfl * mird));
        if (vacd[2] >= 0.0) continue;

        dheight = (PMT_CZ - vacp[2]);
        pmtx = vacp[0] + dheight * (vacd[0] / vacd[2]);
        pmty = vacp[1] + dheight * (vacd[1] / vacd[2]);

        if (std::fabs(pmtx) < PMT_HOLE[0] && std::fabs(pmty) < PMT_HOLE[1]) continue;
        if (std::hypot(pmtx, pmty) < MIRROR_BTM_RADIUS) {
            if (RichOffline::RichPMTsManager::FindChannel(pmtx, pmty) < 0) continue;
            simphs[it][is]  = 2;
            simphsx[it][is] = pmtx;
            simphsy[it][is] = pmty;
        }
    }}
    
    short count_bd = 0;
    for(auto&& phs : border) {
    for (auto&& ph : phs) {
        count_bd += (ph != 0);
    }}
    double trace_bd = static_cast<double>(count_bd) / static_cast<double>(TRACE_NPHI * TRACE_NSET);
    
    short count = 0;
    for(auto&& phs : simphs) {
    for (auto&& ph : phs) {
        count += (ph != 0);
    }}
    double trace = static_cast<double>(count) / static_cast<double>(TRACE_NPHI * TRACE_NSET);
    double crrch = (trace > 0.0) ? (0.1 * ((index_ * index_ - 1.0) / (sinv * sinv)) * std::fabs(radd[2]) / trace) : 0.0;

    if (!std::isfinite(trace_bd) || !std::isfinite(trace) || !std::isfinite(crrch)) return rlt_trace;

    if (!has_cloud) {
        rlt_trace = std::move(std::array<double, 5>({ trace_bd, trace, 1.0, 1.0, crrch }));
        return rlt_trace;
    }

    const double thres = 3.4; // pmt width
    const std::vector<CherenkovHit>& hits = cloud->hits();
    std::vector<short>  hphi(hits.size(), -1);
    std::vector<double> hdxy(hits.size(), -1.0);
    std::vector<short>  hpmt(hits.size(), -1);
    for (int is = 0; is < TRACE_NSET; ++is) {
    for (int iphi = 0; iphi < TRACE_NPHI; ++iphi) {
        if (simphs[iphi][is] == 0) continue;
        for (int ih = 0; ih < hphi.size(); ++ih) {
            const CherenkovHit& hit = hits.at(ih);
            if (simphs[iphi][is] == 1 && (hit.mode() != 0)) continue;
            if (simphs[iphi][is] == 2 && (hit.mode() != 1 && hit.mode() != 2)) continue;
            
            double dist = std::hypot(hit.cx() - simphsx[iphi][is], hit.cy() - simphsy[iphi][is]);
            if (dist > thres || (hdxy.at(ih) >= 0.0 && hdxy.at(ih) < dist)) continue;
            
            hphi.at(ih) = iphi;
            hdxy.at(ih) = dist;
            hpmt.at(ih) = hit.pmtid();
        }
    }}

    for (int ih = hphi.size() - 1; ih >= 0; --ih) {
        if (hphi.at(ih) >= 0) continue;
        hphi.erase(hphi.begin() + ih);
        hpmt.erase(hpmt.begin() + ih);
    }
    std::set<short> setpmt;
    for (int ih = 0; ih < hpmt.size(); ++ih) setpmt.insert(hpmt.at(ih));
    
    double accuracy = static_cast<double>(hphi.size()) / static_cast<double>(hits.size());

    if (hphi.size() < 2 || setpmt.size() < 2) {
        rlt_trace = std::move(std::array<double, 5>({ trace_bd, trace, accuracy, 0.0, crrch }));
        return rlt_trace;
    }

    int uniformity_width = static_cast<int>(0.5 * static_cast<double>(count) / static_cast<double>(setpmt.size()));
    if (uniformity_width <= 0) uniformity_width = 1;

    std::set<int> withphi;
    for (int ih = 0; ih < hphi.size(); ++ih) {
        int init_phi = hphi.at(ih);
        int cnt_lw = 0;
        for (int iphi = init_phi; iphi > (init_phi - TRACE_NPHI); --iphi) {
            if (cnt_lw > uniformity_width) break;
            int idx_phi = (iphi >= 0) ? iphi : (iphi + TRACE_NPHI);
            for (int is = 0; is < TRACE_NSET; ++is) {
                if (simphs[idx_phi][is] == 0) continue;
                int idx = idx_phi * TRACE_NSET + is;
                withphi.insert(idx);
                cnt_lw++;
            }
        }
        int cnt_up = 0;
        for (int iphi = init_phi; iphi < (init_phi + TRACE_NPHI); ++iphi) {
            if (cnt_up > uniformity_width) break;
            int idx_phi = (iphi < TRACE_NPHI) ? iphi : (iphi - TRACE_NPHI);
            for (int is = 0; is < TRACE_NSET; ++is) {
                if (simphs[idx_phi][is] == 0) continue;
                int idx = idx_phi * TRACE_NSET + is;
                withphi.insert(idx);
                cnt_up++;
            }
        }
    }
    
    double uniformity = static_cast<double>(withphi.size()) / static_cast<double>(count);
    if (!std::isfinite(uniformity)) {
        rlt_trace = std::move(std::array<double, 5>({ trace_bd, trace, accuracy, 1.0, crrch }));
        return rlt_trace;
    }

    rlt_trace = std::move(std::array<double, 5>({ trace_bd, trace, accuracy, uniformity, crrch }));
    return rlt_trace;
}
*/



CherenkovRayTrace::CherenkovRayTrace(const std::array<double, 6>& part, double cbta, double index, short kind, short tile) : CherenkovRayTrace() {
    if (cbta <= 0.0 || index <= 1.0 || kind <= 0 || kind >= 3 || tile < 0) return;
    if (cbta <= (1.0 / index)) return;
    double height = RAD_HEIGHT[kind];
    
    double cosv = 1.0 / (index * cbta);
    double sinv = std::sqrt((1.0 - cosv) * (1.0 + cosv));
    if (!std::isfinite(cosv) || !std::isfinite(sinv)) return;

    double lmtrfr = -std::cos(std::asin(1.0 / index));

    SVecD<3> partp(part[0], part[1], part[2]);
    SVecD<3> partd(part[3], part[4], part[5]);
    OrthCoord orth(partd);
    
    double topx = partp[0] + (0.5 * height) * (partd[0] / partd[2]);
    double topy = partp[1] + (0.5 * height) * (partd[1] / partd[2]);
    double topz = partp[2] + (0.5 * height);

    std::vector<std::array<double, 3>> trace;
    std::array<std::array<short,  TRACE_NSET>, TRACE_NPHI> types; // 0 no, 1 rad, 2 dir, 3 rfl
    for (int it = 0; it < TRACE_NPHI; ++it) {
    for (int is = 0; is < TRACE_NSET; ++is) {
        types[it][is] = 0;
        double locwgt = (static_cast<double>(is) + 0.5) / static_cast<double>(TRACE_NSET);
        
        double phi = TWO_PI * static_cast<double>(it) / static_cast<double>(TRACE_NPHI);
        SVecD<3>&& phivec = orth.tau() * std::cos(phi) + orth.rho() * std::sin(phi);
        SVecD<3>&& srcvec = partd * cosv + phivec * sinv;
        if (srcvec[2] >= lmtrfr) continue;

        double planex = topx + (-height) * (locwgt * (partd[0] / partd[2]) + (1.0 - locwgt) * (srcvec[0] / srcvec[2]));
        double planey = topy + (-height) * (locwgt * (partd[1] / partd[2]) + (1.0 - locwgt) * (srcvec[1] / srcvec[2]));
        double planez = topz + (-height);
        if (std::hypot(planex, planey) > MIRROR_TOP_RADIUS) continue;
        
        short otile = RichRingR::getTileIndex(planex, planey);
        if (otile < 0 || otile != tile) continue;
        types[it][is] = 1;

        double cosvac = -1.0 * std::sqrt(1.0 - (index * index) * (1.0 - srcvec[2]) * (1.0 + srcvec[2]));
        double sinvac = std::sqrt((1.0 - cosvac) * (1.0 + cosvac));
        if (!std::isfinite(cosvac) || !std::isfinite(sinvac)) continue;
        
        SVecD<2>&& vacr = ROOT::Math::Unit(SVecD<2>(srcvec[0], srcvec[1])) * sinvac;
        if (!std::isfinite(vacr[0]) || !std::isfinite(vacr[1])) vacr = std::move(SVecD<2>());

        SVecD<3> vacp(planex, planey, planez);
        SVecD<3> vacd(vacr[0], vacr[1], cosvac);

        double dheight = (PMT_CZ - vacp[2]);
        double pmtx = vacp[0] + dheight * (vacd[0] / vacd[2]);
        double pmty = vacp[1] + dheight * (vacd[1] / vacd[2]);
        if (std::fabs(pmtx) < PMT_HOLE[0] && std::fabs(pmty) < PMT_HOLE[1]) continue;
        if (std::hypot(pmtx, pmty) < MIRROR_BTM_RADIUS) { 
            types[it][is] = 2;
            trace.push_back(std::array<double, 3>({ pmtx, pmty, phi })); 
            continue; 
        }

        double vactx = (vacd[0] / vacd[2]);
        double vacty = (vacd[1] / vacd[2]);
        double mirtr = ((MIRROR_TOP_RADIUS - MIRROR_BTM_RADIUS) / MIRROR_HEIGHT);
        double poly0 = (mirtr * mirtr - (vactx * vactx + vacty * vacty));
        double poly1 = 2.0 * (MIRROR_TOP_RADIUS * mirtr - (vacp[0] * vactx + vacp[1] * vacty));
        double poly2 = (MIRROR_TOP_RADIUS * MIRROR_TOP_RADIUS - (vacp[0] * vacp[0] + vacp[1] * vacp[1]));

        double det  = std::sqrt(poly1 * poly1 - 4.0 * poly0 * poly2);
        double sol1 = 0.5 * (-poly1 + det) / poly0;
        double sol2 = 0.5 * (-poly1 - det) / poly0;
        bool is_ok1 = (std::isfinite(sol1) && sol1 < 0.0 && sol1 > dheight);
        bool is_ok2 = (std::isfinite(sol2) && sol2 < 0.0 && sol2 > dheight);

        if (!is_ok1 && !is_ok2) continue;

        double sol = 0;
        if (is_ok1) sol = sol1;
        if (is_ok2) sol = sol2;
        vacp[0] = vacp[0] + sol * vactx;
        vacp[1] = vacp[1] + sol * vacty;
        vacp[2] = vacp[2] + sol;
        
        double cosmir = -1.0 * std::sqrt((mirtr * mirtr) / (1.0 + mirtr * mirtr));
        double sinmir = std::sqrt(1.0 - cosmir * cosmir);
        SVecD<2>&& mirr = ROOT::Math::Unit(SVecD<2>(-vacp[0], -vacp[1])) * sinmir;
        SVecD<3> mird(mirr[0], mirr[1], cosmir);

        double cosrfl = std::fabs(ROOT::Math::Dot(vacd, mird));
        if (!std::isfinite(cosrfl)) continue;
        
        vacd = std::move(ROOT::Math::Unit(vacd + 2.0 * cosrfl * mird));
        if (vacd[2] >= 0.0) continue;

        dheight = (PMT_CZ - vacp[2]);
        pmtx = vacp[0] + dheight * (vacd[0] / vacd[2]);
        pmty = vacp[1] + dheight * (vacd[1] / vacd[2]);

        if (std::fabs(pmtx) < PMT_HOLE[0] && std::fabs(pmty) < PMT_HOLE[1]) continue;
        if (std::hypot(pmtx, pmty) < MIRROR_BTM_RADIUS) {
            types[it][is] = 3;
            trace.push_back(std::array<double, 3>({pmtx, pmty, phi })); 
        }
    }}

    photon_trace_ = trace;
    photon_types_ = types;
    
    short count_border = 0;
    short count_trace  = 0;
    for (auto&& phis : photon_types_) {
        for (auto&& sets : phis) {
            count_border += (sets >= 1);
            count_trace  += (sets >= 2);
        }
    }

    border_ = static_cast<double>(count_border) / static_cast<double>(TRACE_NSET * TRACE_NPHI);
    trace_  = static_cast<double>(count_trace ) / static_cast<double>(TRACE_NSET * TRACE_NPHI);
}

void CherenkovRayTrace::cal(const CherenkovCloud* cloud) {
    accuracy_ = 0.0;
    phis_.clear();

    if (cloud == nullptr) return;
    if (cloud->hits().size() == 0) return;
    if (photon_trace_.size() == 0) return;
   
    short count_acc = 0;
    std::vector<double> phis;

    const double thres = 2.55; // pmt 0.75*width
    const std::vector<CherenkovHit>& hits = cloud->hits();
    for (auto&& hit : hits) {
        double phi = -1.;
        double mini_dist = 1000.;
        for (auto&& ph : photon_trace_) {
            double dist = std::hypot(hit.cx() - ph[0], hit.cy() - ph[1]);
            if (dist > mini_dist) continue;
            mini_dist = dist;
            phi = ph[2];
        }
        phis.push_back(phi);
        count_acc += (mini_dist < thres);
    }

    accuracy_ = static_cast<double>(count_acc) / static_cast<double>(hits.size());
    phis_ = phis;
}


const CherenkovRayTrace::SVecD<3> CherenkovRayTrace::OrthCoord::AXIS_X(1, 0, 0);
const CherenkovRayTrace::SVecD<3> CherenkovRayTrace::OrthCoord::AXIS_Y(0, 1, 0);
const CherenkovRayTrace::SVecD<3> CherenkovRayTrace::OrthCoord::AXIS_Z(0, 0, 1);

void CherenkovRayTrace::OrthCoord::reset(const SVecD<3>& org, const SVecD<3>& seed) {
    Double_t org_mag = ROOT::Math::Mag(org);
    if (org_mag == 0.0) return;
    SVecD<3>&& uorg = org / org_mag;

    SVecD<3> tag = seed; 
    Double_t tag_mag = ROOT::Math::Mag(tag);
    if (tag_mag != 0.0) tag /= tag_mag;
    else {
        tag = std::move(AXIS_X);
        Double_t dotx = std::fabs(ROOT::Math::Dot(uorg, tag));
        if (dotx == 1.0) {
            tag = std::move(AXIS_Y);
            Double_t doty = std::fabs(ROOT::Math::Dot(uorg, tag));
            if (doty == 1.0) {
                tag = std::move(AXIS_Z);
            }
        }
    }

    Double_t&&   dot = ROOT::Math::Dot(uorg, tag);
    SVecD<3>&& cross = ROOT::Math::Cross(uorg, tag);
    Double_t&&  norm = ROOT::Math::Mag(cross);

    org_ = std::move(uorg);
    tau_ = std::move((tag - dot * uorg) / norm);
    rho_ = std::move(cross / norm);
}

//////////////////////////////////////////////////////////////////////////////////////




#endif // __CherenkovMeas_C__
