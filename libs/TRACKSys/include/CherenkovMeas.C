#ifndef __TRACKLibs_CherenkovMeas_C__
#define __TRACKLibs_CherenkovMeas_C__


//#include "Sys.h"
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
    
    search_closest_beta(Numc::ONE<>);
}

        
CherenkovFit::CherenkovFit(const std::vector<CherenkovHit>& args_hits, const std::array<double, 2>& pmtc, double rfr_index, double width_bta, const MultiGaus& pdf_bta, double bta_crr) : CherenkovFit() {
    pmtc_      = pmtc;
    rfr_index_ = rfr_index;
    width_bta_ = width_bta;
    pdf_bta_   = pdf_bta;
    bta_crr_   = bta_crr;
    if (!check()) { clear(); return; }
    if (args_hits.size() == 0) { succ_ = true; return; }

    timer_.start();

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
    
    timer_.stop();

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

    npe_total_ = Numc::ZERO<>;
    npe_stone_ = Numc::ZERO<>;
    npe_cloud_ = Numc::ZERO<>;
    npe_tumor_ = Numc::ZERO<>;
    npe_ghost_ = Numc::ZERO<>;
    npe_emery_ = Numc::ZERO<>;
    npe_other_ = Numc::ZERO<>;
    npe_other_inn_ = Numc::ZERO<>;
    npe_other_out_ = Numc::ZERO<>;

    pmtc_      = std::array<double, 2>({ 0.0, 0.0 });
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
    lmtl_bta_ = (Numc::ONE<> + WIDTH_CORE_COS) / rfr_index_ + width_bta_;
    lmtu_bta_ = Numc::ONE<> + width_bta_ * Numc::SIX<>;
    if (Numc::Compare(lmtl_bta_) <= 0 || Numc::Compare(lmtu_bta_) <= 0) return false;
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

            double dist = Numc::INV_SQRT_TWO * std::hypot(ihit.cx() - jhit.cx(), ihit.cy() - jhit.cy()) / width_stone;
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
        
        hit.search_closest_beta(Numc::ONE<>);
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
        hit.set_wgt(Numc::ONE<>);
        hit.set_cnt(Numc::ONE<>);
        hit.set_lmt_bta(lmtl_bta_, lmtu_bta_);
        
        if (is_within_pmtc(hit.cx(), hit.cy())) { hit.set_cluster(CherenkovHit::Cluster::emery); }
        else                                    { hit.set_cluster(CherenkovHit::Cluster::other); }
        
        if (rltclds.size() == 0) continue;
        if (hit.type() == 0) continue;
        double beta = Numc::ONE<>;

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
    for (auto&& hit : hits) { hit.set_wgt(Numc::ONE<>); hit.set_cnt(Numc::ONE<>); }
    if (weighted) {
        double sumnpe = Numc::ZERO<>;
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
        double sumprb = Numc::ZERO<>;
        for (auto&& stn : cand_cstns) {
            double dcx = (hit.cx() - stn[0]) / WIDTH_COO;
            double dcy = (hit.cy() - stn[1]) / WIDTH_COO;
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
            if (Numc::Compare(elms.at(it), CONVG_PROB_SGM50) < 0) continue;
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
        double chisqc = 0.;
        for (auto&& hit : cand_chit) {
            cnt += hit.cnt();
            npe += hit.npe();
            double dcx = (hit.cx() - cand_cstn[0]) / WIDTH_COO;
            double dcy = (hit.cy() - cand_cstn[1]) / WIDTH_COO;
            double nrm = Numc::INV_SQRT_TWO * std::hypot(dcx,  dcy);
            chisq += hit.cnt() * (nrm * nrm); 
            
            double dcxc = (hit.cx() - cand_cstn[0]) / WIDTH_CELL;
            double dcyc = (hit.cy() - cand_cstn[1]) / WIDTH_CELL;
            double nrmc = Numc::INV_SQRT_TWO * std::hypot(dcxc,  dcyc);
            chisqc += (nrmc * nrmc);
        }

        double ndof = (cnt - Numc::ONE<>);
        double nchi = chisq / ndof;
        double chic = chisqc / static_cast<double>(cand_chit.size() - 1);
        if (!Numc::Valid(ndof) || Numc::Compare(ndof) <= 0) continue;
        if (!Numc::Valid(nchi)) continue;
        if (!Numc::Valid(chic)) continue;

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
                double nrmx = (hit.cx() - stn[0]) / WIDTH_COO;
                double nrmy = (hit.cy() - stn[1]) / WIDTH_COO;
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
                    double dnxA = (newx - gausA[0]) / WIDTH_COO;
                    double dnyA = (newy - gausA[1]) / WIDTH_COO;
                    double dnxB = (newx - gausB[0]) / WIDTH_COO;
                    double dnyB = (newy - gausB[1]) / WIDTH_COO;
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


std::vector<CherenkovCloud> CherenkovFit::build_cloud(const std::vector<CherenkovHit>& args_hits, bool weighted) {
    std::vector<CherenkovCloud> clds;
    if (args_hits.size() == 0) return clds;
    
    std::vector<CherenkovHit> hits = args_hits;
    for (auto&& hit : hits) hit.set_lmt_bta(lmtl_bta_, lmtu_bta_);

    for (auto&& hit : hits) { hit.set_wgt(Numc::ONE<>); hit.set_cnt(Numc::ONE<>); }
    if (weighted) {
        double sumnpe = Numc::ZERO<>;
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
        
        for (auto&& hit : ubr_hits) { hit.set_wgt(Numc::ONE<>); hit.set_cnt(Numc::ONE<>); }
        if (weighted) {
            double sumnpe = Numc::ZERO<>;
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
                if (Numc::Compare(wgtgbl, (Numc::TWO<> / Numc::THREE<>)) <= 0) continue;
                if (Numc::Compare(wgtloc, (Numc::TWO<> / Numc::THREE<>)) <= 0) continue;
                
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

    double weight = Numc::ZERO<>;
    for (auto&& hit : cand_chit) weight += hit.wgt();
    if (Numc::EqualToZero(weight)) return CherenkovCloud();
    for (auto&& hit : cand_chit) hit.set_cnt(static_cast<double>(cand_chit.size()) * (hit.wgt() / weight));

    double cand_beta = cand_cld.beta();
    double cand_cnt  = Numc::ZERO<>;
    double cand_nchi = Numc::ZERO<>;

    // refit cloud
    bool succ = false;
    for (int iter = 0; iter <= LMTMAX_ITER; ++iter) {
        double count = Numc::ZERO<>;
        double chisq = Numc::ZERO<>;
        double grdB  = Numc::ZERO<>;
        double cvBB  = Numc::ZERO<>;
        
        for (auto&& hit : cand_chit) {
            double dlt = hit.search_closest_beta(cand_beta) - cand_beta;
            double nrm = dlt / width_bta_;

            double smooth = SMOOTH_BOUND[1] - static_cast<double>(iter) * SMOOTH_RATE;
            if (smooth < SMOOTH_BOUND[0]) smooth = SMOOTH_BOUND[0];
            
            double smnrm = std::fabs(nrm) - smooth;
            double smprb = (Numc::Compare(smnrm) <= 0) ? Numc::ONE<> : std::exp(-Numc::ONE_TO_TWO * smnrm * smnrm);

            std::array<long double, 3>&& minib = pdf_bta_.minimizer(dlt);
            
            count += hit.cnt();
            chisq += hit.cnt() * (minib[0] * minib[0]);
            grdB  += hit.cnt() * smprb * (Numc::NEG<> * minib[2] * minib[1]);
            cvBB  += hit.cnt() * smprb * (minib[2] * minib[2]);
        }
        if (!Numc::Valid(count) || Numc::Compare(count, Numc::ONE<>) <= 0) break;

        double nchi = chisq / (count - Numc::ONE<>);
        if (!Numc::Valid(nchi)) break;

        double resB = grdB / cvBB;
        double beta = cand_beta - resB;
        if (!Numc::Valid(resB) || Numc::Compare(beta) <= 0) break;
       
        if (iter >= LMTM_ITER) {
            succ = (Numc::Compare(std::fabs(cand_beta - beta), CONVG_TOLERANCE) <= 0) && 
                   (Numc::Compare(std::fabs(cand_cnt - count), CONVG_TOLERANCE) <= 0) &&
                   (Numc::Compare(std::fabs(cand_nchi - nchi), CONVG_TOLERANCE) <= 0);
            if (succ) break;
        }

        cand_beta = beta;
        cand_cnt  = count;
        cand_nchi = nchi;
    }
    if (!succ) return CherenkovCloud();
  
    // check result
    int count_within_three_sigma = 0;
    int count_within_one_half_sigma = 0;
    std::vector<CherenkovHit> cand_reduce_chit;
    for (auto&& hit : cand_chit) {
        double absnrm = std::fabs(hit.search_closest_beta(cand_beta) - cand_beta) / width_bta_;
        if (Numc::Compare(absnrm, Numc::FIVE<>) > 0) continue;
        cand_reduce_chit.push_back(hit);

        if (Numc::Compare(absnrm, Numc::THREE<>) < 0) count_within_three_sigma++;
        if (Numc::Compare(absnrm, Numc::ONE<> + Numc::ONE_TO_TWO) < 0) count_within_one_half_sigma++;
    } 
    if (cand_reduce_chit.size() < LMTMIN_CLOUD_HITS) return CherenkovCloud();
    std::sort(cand_reduce_chit.begin(), cand_reduce_chit.end(), CherenkovHit_sort());
    
    if (count_within_three_sigma < LMTMIN_CLOUD_HITS) return CherenkovCloud();
    if (count_within_one_half_sigma <= Numc::ONE<short>) return CherenkovCloud();
    
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
    double cand_ndof = cand_cnt - Numc::ONE<>;
    if (!Numc::Valid(cand_ndof) || Numc::Compare(cand_ndof) <= 0) return CherenkovCloud();
    if (!Numc::Valid(cand_nchi)) return CherenkovCloud();

    double cand_misjudge = Numc::ZERO<>; 

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
            std::array<long double, 3>&& minib = scan_bta_.minimizer(dlt);
            chisq += hit.cnt() * (minib[0] * minib[0]);
        }

        double ndof = (cnt - Numc::ONE<>);
        double nchi = chisq / ndof;
        if (!Numc::Valid(ndof) || Numc::Compare(ndof) <= 0) continue;
        if (!Numc::Valid(nchi)) continue;
    
        double misjudge = Numc::ZERO<>;
    
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
    
    double weight_hits = Numc::ZERO<>;
    for (auto&& hit : hits) weight_hits += hit.wgt();
    if (Numc::EqualToZero(weight_hits)) return cclds;
    for (auto&& hit : hits) hit.set_cnt(static_cast<double>(hits.size()) * (hit.wgt() / weight_hits));

    for (auto&& hit : hits) {
        if (hit.type() == 0) continue;
        double multi = hit.cnt() / static_cast<double>(hit.hasDb() * Numc::TWO<> + hit.hasRbA() + hit.hasRbB());
        if (hit.hasDb() ) cclds.push_back(std::array<double, 2>({ hit.dbta() , Numc::TWO<> * multi }));
        if (hit.hasRbA()) cclds.push_back(std::array<double, 2>({ hit.rbtaA(), Numc::ONE<> * multi }));
        if (hit.hasRbB()) cclds.push_back(std::array<double, 2>({ hit.rbtaB(), Numc::ONE<> * multi }));
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


std::vector<CherenkovTumor> CherenkovFit::build_tumor(const std::vector<CherenkovHit>& args_hits, bool weighted) {
    std::vector<CherenkovTumor> tmrs;
    if (args_hits.size() == 0) return tmrs;

    std::vector<CherenkovHit> hits = args_hits;
    for (auto&& hit : hits) hit.set_lmt_bta(lmtl_bta_, lmtu_bta_ + width_bta_ * Numc::FOUR<>);
    
    for (auto&& hit : hits) { hit.set_wgt(Numc::ONE<>); hit.set_cnt(Numc::ONE<>); }
    if (weighted) {
        double sumnpe = Numc::ZERO<>;
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

    double weight = Numc::ZERO<>;
    for (auto&& hit : hits) weight += hit.wgt();
    if (Numc::EqualToZero(weight)) return CherenkovTumor();
    for (auto&& hit : hits) hit.set_cnt(static_cast<double>(hits.size()) * (hit.wgt() / weight));

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
        short iter = 1;
        bool is_ok = false;
        while (!is_ok && iter <= LMTMAX_ITER) {
            double sum_bta = Numc::ZERO<>;
            double sum_cnt = Numc::ZERO<>;
            for (auto&& hit : hits) {
                if (hit.type() == 0) continue;
                sum_bta += hit.cnt() * hit.search_closest_beta(beta);
                sum_cnt += hit.cnt();
            }
            if (Numc::EqualToZero(sum_cnt)) break;
            double new_beta = (sum_bta / sum_cnt);
            is_ok = (iter >= LMTL_ITER && Numc::Compare(std::fabs(new_beta - beta), CONVG_TOLERANCE) <= 0);
            if (!is_ok) beta = new_beta;
            iter++;
        }
        if (!is_ok) beta = Numc::ZERO<>;
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
            double nrm = dlt / width_bta_;
            
            double smooth = SMOOTH_BOUND[1] - static_cast<double>(iter) * SMOOTH_RATE;
            if (smooth < SMOOTH_BOUND[0]) smooth = SMOOTH_BOUND[0];
            
            double smnrm = std::fabs(nrm) - smooth;
            double smprb = (Numc::Compare(smnrm) <= 0) ? Numc::ONE<> : std::exp(-Numc::ONE_TO_TWO * smnrm * smnrm);
            
            std::array<long double, 3>&& minib = pdf_bta_.minimizer(dlt);
            
            count += hit.cnt();
            chisq += hit.cnt() * (minib[0] * minib[0]);
            grdB  += hit.cnt() * smprb * (Numc::NEG<> * minib[2] * minib[1]);
            cvBB  += hit.cnt() * smprb * (minib[2] * minib[2]);
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
    
    // check result
    int count_within_three_sigma = 0;
    std::vector<CherenkovHit> reduce_hits;
    for (auto&& hit : hits) {
        double absnrm = std::fabs(hit.search_closest_beta(beta) - beta) / width_bta_;
        if (Numc::Compare(absnrm, Numc::FIVE<>) > 0) continue;
        reduce_hits.push_back(hit);

        if (Numc::Compare(absnrm, Numc::THREE<>) < 0) count_within_three_sigma++;
    } 
    if (reduce_hits.size() < LMTMIN_TUMOR_HITS) return CherenkovTumor();
    std::sort(reduce_hits.begin(), reduce_hits.end(), CherenkovHit_sort());
    
    if (count_within_three_sigma  < LMTMIN_TUMOR_HITS) return CherenkovTumor();
    
    std::map<int, int> pmt_maps;
    for (auto&& hit : reduce_hits) {
        if (pmt_maps.count(hit.pmtid()) == 0) pmt_maps[hit.pmtid()] = 1;
        else pmt_maps[hit.pmtid()]++;
    }
    short nhit = reduce_hits.size();
    short npmt = pmt_maps.size();
    if (nhit < LMTMIN_TUMOR_HITS) return CherenkovTumor();
    if (npmt < LMTMIN_TUMOR_PMTS) return CherenkovTumor();
        
    double npe = Numc::ZERO<>;
    for (auto&& hit : reduce_hits) npe += hit.npe();

    for (auto&& hit : reduce_hits) hit.set_cluster(CherenkovHit::Cluster::tumor);

    return CherenkovTumor(reduce_hits, nhit, npmt, beta, cbta, npe);
}
        

std::vector<CherenkovGhost> CherenkovFit::build_ghost(const std::vector<CherenkovHit>& args_hits, bool weighted) {
    std::vector<CherenkovGhost> gsts;
    if (args_hits.size() == 0) return gsts;

    std::vector<CherenkovHit> hits = args_hits;
    for (auto&& hit : hits) hit.set_lmt_bta(lmtl_bta_, lmtu_bta_ + width_bta_ * Numc::SIX<>);
    
    for (auto&& hit : hits) { hit.set_wgt(Numc::ONE<>); hit.set_cnt(Numc::ONE<>); }
    if (weighted) {
        double sumnpe = Numc::ZERO<>;
        for (auto&& hit : hits) sumnpe += hit.npe();
        for (auto&& hit : hits) hit.set_wgt(static_cast<double>(hits.size()) * (hit.npe() / sumnpe));
    }
    
    // ghost (cloud)
    auto&& gtable = make_group_table(hits, false, Numc::ONE_TO_THREE * WIDTH_PMT, true, (Numc::ONE<> + Numc::ONE_TO_TWO) * width_bta_);
    for (auto&& table : gtable) {
        if (table.size() < LMTMIN_GHOST_HITS) continue;
        
        short maxnhits = 0;
        auto&& gstones = make_group_table(table, true, Numc::ONE_TO_THREE * WIDTH_PMT, false, (Numc::ONE<> + Numc::ONE_TO_TWO) * width_bta_);
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

    double weight = Numc::ZERO<>;
    for (auto&& hit : hits) weight += hit.wgt();
    if (Numc::EqualToZero(weight)) return CherenkovGhost();
    for (auto&& hit : hits) hit.set_cnt(static_cast<double>(hits.size()) * (hit.wgt() / weight));

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
        short iter = 1;
        bool is_ok = false;
        while (!is_ok && iter <= LMTMAX_ITER) {
            double sum_bta = Numc::ZERO<>;
            double sum_cnt = Numc::ZERO<>;
            for (auto&& hit : hits) {
                if (hit.type() == 0) continue;
                sum_bta += hit.cnt() * hit.search_closest_beta(beta);
                sum_cnt += hit.cnt();
            }
            if (Numc::EqualToZero(sum_cnt)) break;
            double new_beta = (sum_bta / sum_cnt);
            is_ok = (iter >= LMTL_ITER && Numc::Compare(std::fabs(new_beta - beta), CONVG_TOLERANCE) <= 0);
            if (!is_ok) beta = new_beta;
            iter++;
        }
        if (!is_ok) beta = Numc::ZERO<>;
    }
    if (Numc::Compare(beta) <= 0) return CherenkovGhost();
    double cbta = bta_crr_ * beta;
    
    std::map<int, int> pmt_maps;
    for (auto&& hit : hits) {
        if (pmt_maps.count(hit.pmtid()) == 0) pmt_maps[hit.pmtid()] = 1;
        else pmt_maps[hit.pmtid()]++;
    }
    short nhit = hits.size();
    short npmt = pmt_maps.size();
    if (nhit < LMTMIN_GHOST_HITS) return CherenkovGhost();
    
    double npe = Numc::ZERO<>;
    for (auto&& hit : hits) npe += hit.npe();

    for (auto&& hit : hits) hit.set_cluster(CherenkovHit::Cluster::ghost);

    return CherenkovGhost(hits, nhit, npmt, beta, cbta, npe);
}


std::array<long double, 3> CherenkovMeas::minimizer(long double x, long double ibta) const {
    //if (Numc::Compare(ibta, Numc::ONE<long double>) < 0)
    //    return std::array<long double, 3>({ Numc::ZERO<long double>, Numc::ZERO<long double>, Numc::ZERO<long double> });

    //bool rescl = false;
    //long double eftsgm = Numc::ONE<long double>;
    //if (Numc::Compare(rfr_, Numc::ONE<long double>) > 0) {
    //    long double thres = rfr_ - Numc::TWO<long double> * sgm_;
    //    if (Numc::Compare(ibta, thres) > 0) {
    //        long double extsgm = (ibta - thres) / sgm_;
    //        eftsgm = Numc::ONE<long double> + extsgm * extsgm;
    //        rescl = true;
    //    }
    //}

    //long double sclx = (rescl ? (x / eftsgm) : x);
    //return mgs_.minimizer(sclx);
    
    return mgs_.minimizer(x); // testcode
}


} // namesapce TrackSys


#endif // __TRACKLibs_CherenkovMeas_C__
