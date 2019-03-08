#ifndef __TRACKLibs_CherenkovMeas_C__
#define __TRACKLibs_CherenkovMeas_C__


#include "Sys.h"
#include "Math.h"
#include "CherenkovMeas.h"


namespace TrackSys {

        
CherenkovFit::CherenkovFit(const std::vector<CherenkovHit>& args_hits, const std::array<double, 2>& args_core, double width_bta, const std::array<double, 4>& args_bta, const std::array<double, 4>& args_npe, double bta_crr) {
    args_core_    = std::array<double, 2>({ args_core[0], args_core[1] });
    args_bta_     = std::array<double, 4>({ args_bta[0], args_bta[1], args_bta[2], args_bta[3] });
    args_npe_sig_ = std::array<double, 3>({ args_npe[0], args_npe[1], args_npe[2] });
    args_npe_nos_ = args_npe[3];
    width_bta_    = width_bta;
    bta_crr_      = bta_crr;
    if (!check()) { clear(); return; }
    if (args_hits.size() == 0) { clear(); return; }
    timer_.start();
            
    std::vector<CherenkovHit> org_hits;
    for (auto&& args_hit : args_hits) {
        if (!args_hit.status()) continue;
        double prbnpe = LandauGaus::Func(args_hit.npe(), args_npe_sig_[0], args_npe_sig_[1], args_npe_sig_[2]);
        double wgtnpe = prbnpe * (Numc::ONE<> + args_npe_nos_) / (prbnpe + args_npe_nos_);
        //if (Numc::Compare(wgtnpe) <= 0) continue;
        org_hits.push_back(args_hit);
        //org_hits.back().set_wgt(wgtnpe);
    }
    if (org_hits.size() == 0) { clear(); return; }
   
    // stone
    stns_ = std::move(fit_stone(org_hits));
    std::vector<CherenkovHit> cld_hits = org_hits;
    while (stns_.size() != 0) {
        cld_hits.clear();
        for (auto&& hit : org_hits) {
            bool is_stn = false;
            for (auto&& stn : stns_) {
                for (auto&& stn_hit : stn.hits()) {
                    if (stn_hit.chann() != hit.chann()) continue;
                    is_stn = true;
                    break;
                }
                if (is_stn) break;
            }
            if (is_stn) continue;
            cld_hits.push_back(hit);
        }
        std::vector<CherenkovStone>&& other_stns = fit_stone(cld_hits);
        if (other_stns.size() == 0) break;

        for (auto&& stn : other_stns) stns_.push_back(stn);
        if (stns_.size() > 1) std::sort(stns_.begin(), stns_.end(), CherenkovStone_sort());
    }
    
    // cloud
    std::vector<CherenkovCloud>&& cand_clds = fit_cloud(cld_hits);
    for (auto&& cand_cld : cand_clds) {
        CherenkovCloud&& cld = refit_cloud(cand_cld);
        if (!cld.status()) continue;
        clds_.push_back(cld);
    }
    if (clds_.size() > 1) std::sort(clds_.begin(), clds_.end(), CherenkovCloud_sort());
    
    succ_ = (stns_.size() != 0 || clds_.size() != 0);
    
    timer_.stop();

    if (!succ_) { clear(); return; }
}


void CherenkovFit::clear() {
    succ_ = false;
    stns_.clear();
    clds_.clear();

    args_bta_     = std::array<double, 4>({0, 0, 0, 0});
    args_npe_sig_ = std::array<double, 3>({0, 0, 0});
    args_npe_nos_ = Numc::ZERO<>;
    width_bta_    = Numc::ZERO<>; 
    bta_crr_      = Numc::ONE<>; 
    pdf_bta_      = std::move(MultiGaus());
    pdf_bta_acc_  = std::move(MultiGaus());

    timer_.clear();
}


bool CherenkovFit::check() {
    if (Numc::Compare(args_bta_[0]) <= 0 || 
        Numc::Compare(args_bta_[1]) <= 0 || 
        Numc::Compare(args_bta_[2]) <= 0 || 
        Numc::Compare(args_bta_[3]) <= 0) return false; 
    if (Numc::Compare(args_npe_sig_[0]) <= 0 || 
        Numc::Compare(args_npe_sig_[1]) <= 0 || 
        Numc::Compare(args_npe_sig_[2]) <= 0 || 
        Numc::Compare(args_npe_nos_)    <  0) return false;
    if (Numc::Compare(width_bta_) <= 0) return false;
    if (Numc::Compare(bta_crr_) < 0) return false;
    pdf_bta_     = std::move(MultiGaus(Robust::Option(Robust::Opt::ON, 3.0L, 0.5L), width_bta_));
    pdf_bta_acc_ = std::move(MultiGaus(Robust::Option(Robust::Opt::ON, 3.0L, 0.5L), args_bta_[0], args_bta_[1], args_bta_[2], args_bta_[3]));
    return true;
}
        

std::vector<CherenkovStone> CherenkovFit::fit_stone(std::vector<CherenkovHit>& hits) {
    std::vector<CherenkovStone> stns;
    if (hits.size() < 2) return stns;
    
    std::vector<std::array<double, 3>>&& cand_cstns = clustering_stone(hits); // (cx cy wgt)
    if (cand_cstns.size() == 0) return stns;

    std::vector<std::vector<CherenkovHit>> cand_chits(cand_cstns.size(), std::vector<CherenkovHit>());
    for (auto&& hit : hits) {
        std::vector<double> cnts;
        for (auto&& stn : cand_cstns) {
            double dcx = (hit.cx() - stn[0]) / WIDTH_CELL;
            double dcy = (hit.cy() - stn[1]) / WIDTH_CELL;
            double nrm = Numc::INV_SQRT_TWO * std::hypot(dcx,  dcy);
            double pdf = std::exp(Numc::NEG<> * nrm * nrm);
            double prb = stn[2] * pdf;
            
            if (Numc::Compare(prb * hits.size(), CONVG_PROB_SGM70) > 0) cnts.push_back(prb);
            else                                                        cnts.push_back(Numc::ZERO<>);
        }
        double sumprb = std::accumulate(cnts.begin(), cnts.end(), 0.0);
        if (Numc::Compare(sumprb) <= 0) continue;
        for (auto&& cnt : cnts) {
            cnt = (cnt / sumprb);
            if (Numc::Compare(cnt, CONVG_PROB_SGM30) < 0) cnt = 0.0;
        }
        sumprb = std::accumulate(cnts.begin(), cnts.end(), 0.0);
        for (auto&& cnt : cnts) cnt = (cnt / sumprb);
        
        for (int it = 0; it < cand_chits.size(); ++it) {
            if (Numc::EqualToZero(cnts.at(it))) continue;
            std::vector<CherenkovHit>& chit = cand_chits.at(it);
            chit.push_back(hit);
            chit.back().set_cnt(cnts.at(it));
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

        double dist = std::hypot((cand_cstn[0] - args_core_[0]), (cand_cstn[1] - args_core_[1]));

        double cnt = 0.;
        double npe = 0.;
        double chisq = 0.;
        for (auto&& hit : cand_chit) {
            cnt += hit.cnt();
            npe += hit.npe();
            double dcx = (hit.cx() - cand_cstn[0]) / WIDTH_CELL;
            double dcy = (hit.cy() - cand_cstn[1]) / WIDTH_CELL;
            double nrm = Numc::INV_SQRT_TWO * std::hypot(dcx,  dcy);
            chisq += nrm * nrm; 
        }

        double ndof = (cnt - Numc::ONE<>);
        double nchi = chisq / ndof;
        double qlt  = Numc::NormQuality(nchi, ndof);
        if (!Numc::Valid(ndof) || Numc::Compare(ndof) <= 0) continue;
        if (!Numc::Valid(nchi)) continue;
        if (!Numc::Valid(qlt)) continue;
        
        stns.push_back(CherenkovStone(cand_chit, cand_cstn[0], cand_cstn[1], dist, npe, cnt, nchi, qlt));
    }
    if (stns.size() > 1) std::sort(stns.begin(), stns.end(), CherenkovStone_sort());

    return stns;
}


std::vector<std::array<double, 3>> CherenkovFit::clustering_stone(std::vector<CherenkovHit>& hits) {
   std::vector<std::array<double, 3>> cstns; // (cx cy wgt)
    if (hits.size() == 0) return cstns;

    cstns = std::move(clustering_evolve_stone(hits));
    if (cstns.size() == 0) return cstns;

    bool succ = false;
    for (int iter = 1; iter <= LMTMAX_ITER; ++iter) {
        std::vector<std::array<double, 3>> mgaus(cstns.size(), std::array<double, 3>({ 0.0, 0.0, 0.0 })); // (cx cy wgt)
        double weight = Numc::ZERO<>;

        for (auto&& hit : hits) {
            double sumprb = Numc::ZERO<>;
            double sumchi = Numc::ZERO<>;
            std::vector<std::array<double, 3>> elms; // (cx cy wgt)
            for (auto&& cstn : cstns) {
                if (Numc::Compare(cstn[2]) <= 0) {
                    elms.push_back(std::array<double, 3>({ 0.0, 0.0, 0.0 }));
                    continue;
                }
                double nrmx = (hit.cx() - cstn[0]) / WIDTH_CELL;
                double nrmy = (hit.cy() - cstn[1]) / WIDTH_CELL;
                double nrm  = Numc::INV_SQRT_TWO * std::hypot(nrmx, nrmy);
                double pdf  = std::exp(Numc::NEG<> * nrm * nrm);
                double prb  = cstn[2] * pdf;
                
                elms.push_back(std::array<double, 3>({ hit.cx(), hit.cy(), prb }));
                sumprb += prb;
                sumchi += prb * nrm * nrm;
            }
            if (!Numc::Valid(sumprb) || Numc::Compare(sumprb) <= 0) continue;
            for (auto&& elm : elms) elm.at(2) = (elm.at(2) / sumprb);
            weight++;

            for (int it = 0; it < elms.size(); ++it) {
                std::array<double, 3>& elm  = elms.at(it);
                std::array<double, 3>& gaus = mgaus.at(it);
                if (Numc::Compare(elm[2]) <= 0) continue;
                gaus[0] += elm[2] * elm[0];
                gaus[1] += elm[2] * elm[1];
                gaus[2] += elm[2];
            }
        }
        if (!Numc::Valid(weight) || Numc::Compare(weight, Numc::ONE<>) <= 0) break;
        
        for (auto&& gaus : mgaus) {
            if (Numc::Compare(gaus[2]) <= 0) continue;
            gaus[0] = gaus[0] / gaus[2];
            gaus[1] = gaus[1] / gaus[2];
            gaus[2] = gaus[2] / weight;
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
                    double dnxA = (newx - gausA[0]) / WIDTH_CELL;
                    double dnyA = (newy - gausA[1]) / WIDTH_CELL;
                    double dnxB = (newx - gausB[0]) / WIDTH_CELL;
                    double dnyB = (newy - gausB[1]) / WIDTH_CELL;
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
                auto&& cstn = cstns.at(it);
                bool match = (Numc::Compare(std::fabs(gaus[0] - cstn[0]), CONVG_TOLERANCE) <= 0 && 
                              Numc::Compare(std::fabs(gaus[1] - cstn[1]), CONVG_TOLERANCE) <= 0 && 
                              Numc::Compare(std::fabs(gaus[2] - cstn[2]), CONVG_TOLERANCE) <= 0);
                if (!match) { is_convg = false; break; }
            }
            if (is_convg) { succ = true; break; }
        }
        cstns = mgaus;
    }
    
    if (!succ) cstns.clear();
    return cstns;
}
    

std::vector<std::array<double, 3>> CherenkovFit::clustering_evolve_stone(std::vector<CherenkovHit>& hits) {
    std::vector<std::array<double, 3>> seeds; // (cx cy wgt)
    if (hits.size() == 0) return seeds;
    
    std::map<int, std::vector<CherenkovHit*>> pmt_maps;
    std::sort(hits.begin(), hits.end(), CherenkovHit_sort());
    for (auto&& hit : hits) pmt_maps[hit.pmtid()].push_back(&hit);
    
    for (auto&& pmt_map : pmt_maps) {
        double cx = Numc::ZERO<>;
        double cy = Numc::ZERO<>;
        double wgt = static_cast<double>((pmt_map.second).size());
        for (auto&& hit : pmt_map.second) {
            cx += hit->cx();
            cy += hit->cy();
        }
        seeds.push_back(std::array<double, 3>({ (cx / wgt), (cy / wgt), wgt }));
    }

    if (seeds.size() != 0) {
        double weight = Numc::ZERO<>;
        for (auto&& seed : seeds) weight += seed[2];
        for (auto&& seed : seeds) seed[2] /= weight;
        if (!Numc::Valid(weight)) seeds.clear();
    }
    if (seeds.size() == 0) return seeds;

    bool succ = false;
    for (int iter = 1; iter <= LMTMAX_ITER; ++iter) {
        std::vector<std::array<double, 3>> mgaus(seeds.size(), std::array<double, 3>({ 0.0, 0.0, 0.0 })); // (cx cy wgt)
        double weight = Numc::ZERO<>;
        for (auto&& hit : hits) {
            int index = 0;
            for (auto&& seed : seeds) {
                double nrmx = (hit.cx() - seed[0]) / WIDTH_CELL;
                double nrmy = (hit.cy() - seed[1]) / WIDTH_CELL;
                double nrm  = Numc::INV_SQRT_TWO * std::hypot(nrmx, nrmy);
                if (Numc::Compare(nrm, Numc::THREE<>) > 0) { index++; continue; }

                std::array<double, 3>& gaus = mgaus.at(index);
                gaus[0] += hit.cx();
                gaus[1] += hit.cy();
                gaus[2] += Numc::ONE<>;
                index++;
            }
            weight++;
        }
        for (auto&& gaus : mgaus) {
            if (Numc::EqualToZero(gaus[2])) continue;
            gaus[0] = gaus[0] / gaus[2];
            gaus[1] = gaus[1] / gaus[2];
            gaus[2] = gaus[2] / weight;
        }

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
                auto&& gausA = mgaus.at(it);
                auto&& gausB = mgaus.at(jt);
                double dnx   = (gausA[0] - gausB[0]) / WIDTH_CELL;
                double dny   = (gausA[1] - gausB[1]) / WIDTH_CELL;
                double dnc   = Numc::INV_SQRT_TWO * std::hypot(dnx, dny);
                if (Numc::Compare(dnc, Numc::ONE<>) < 0) set_ptr->push_back(jt);
            }
        }
        
        bool has_merge = (mgaus.size() != sets.size());
        if (has_merge) {
            std::vector<std::array<double, 3>> newmg; // (cx cy wgt)
            for (auto&& set : sets) {
                if (set.size() == 1) {
                    newmg.push_back(mgaus.at(set.at(0))); 
                    continue; 
                }
                double newcx  = Numc::ZERO<>;
                double newcy  = Numc::ZERO<>;
                for (auto&& idx : set) {
                    auto&& gaus = mgaus.at(idx);
                    newcx += gaus[0];
                    newcy += gaus[1];
                }
                newcx = (newcx / static_cast<double>(set.size()));
                newcy = (newcy / static_cast<double>(set.size()));
                newmg.push_back(std::array<double, 3>({ newcx, newcy, 0.0 }));
            }
            mgaus = newmg;
        }

        if (iter >= LMTL_ITER && !has_merge) {
            bool is_stable = true;
            for (int it = 0; it < mgaus.size(); ++it) {
                auto&& gaus = mgaus.at(it);
                auto&& seed = seeds.at(it);
                bool match = (Numc::Compare(std::fabs(gaus[0] - seed[0]), CONVG_TOLERANCE) <= 0 && 
                              Numc::Compare(std::fabs(gaus[1] - seed[1]), CONVG_TOLERANCE) <= 0 && 
                              Numc::Compare(std::fabs(gaus[2] - seed[2]), CONVG_TOLERANCE) <= 0);
                if (!match) { is_stable = false; break; }
            }
            if (is_stable) { succ = true; break; }
        }
        seeds = mgaus;
    }
    
    if (!succ) seeds.clear();
    return seeds;
}


CherenkovCloud CherenkovFit::refit_cloud(const CherenkovCloud& cand_cld) {
    if (!cand_cld.status()) return CherenkovCloud();
    std::vector<CherenkovHit> cand_hits = cand_cld.hits();
    double cand_beta = cand_cld.beta();
    double cand_nchi = Numc::ZERO<>;
    double cand_weight = Numc::ZERO<>;

    // remove stone-like
    std::vector<CherenkovStone>&& stonelike = fit_stone(cand_hits);
    if (stonelike.size() != 0) { 
        for (auto&& stn : stonelike) stns_.push_back(stn);
        if (stns_.size() > 1) std::sort(stns_.begin(), stns_.end(), CherenkovStone_sort());
        
        std::vector<CherenkovHit> cld_hits;
        for (auto&& hit : cand_hits) {
            bool is_stn = false;
            for (auto&& stn : stonelike) {
                for (auto&& stn_hit : stn.hits()) {
                    if (stn_hit.chann() != hit.chann()) continue;
                    is_stn = true;
                    break;
                }
                if (is_stn) break;
            }
            if (is_stn) continue;
            cld_hits.push_back(hit);
        }
        cand_hits = cld_hits;
    }
    if (cand_hits.size() < LMTMIN_CLOUD_HITS) return CherenkovCloud();

    // refit cloud
    bool succ = false;
    for (int iter = 0; iter <= LMTMAX_ITER; ++iter) {
        double weight = Numc::ZERO<>;
        double chisq  = Numc::ZERO<>;
        double grdB   = Numc::ZERO<>;
        double cvBB   = Numc::ZERO<>;

        for (auto&& hit : cand_hits) {
            double dlt = hit.search_closed_beta(cand_beta) - cand_beta;
            double nrm = dlt / width_bta_;
            
            if (Numc::Compare(std::fabs(nrm), Numc::FIVE<>) > 0) continue;
            std::array<long double, 3>&& minib = pdf_bta_acc_.minimizer(dlt);
            
            weight += hit.wgt();
            chisq  += hit.wgt() * (minib[0] * minib[0]);
            grdB   += hit.wgt() * (Numc::NEG<> * minib[2] * minib[1]);
            cvBB   += hit.wgt() * (minib[2] * minib[2]);
        }
        if (!Numc::Valid(weight) || Numc::Compare(weight, Numc::ONE<>) <= 0) break;

        double nchi = chisq / (weight - Numc::ONE<>);
        if (!Numc::Valid(nchi)) break;

        double resB = grdB / cvBB;
        double beta = cand_beta - resB;
        if (!Numc::Valid(resB) || Numc::Compare(beta) <= 0) break;
       
        if (iter >= LMTL_ITER) {
            succ = (Numc::Compare(std::fabs(cand_beta - beta), CONVG_TOLERANCE) <= 0) && 
                   (Numc::Compare(std::fabs(cand_nchi - nchi), CONVG_TOLERANCE) <= 0) &&
                   (Numc::Compare(std::fabs(cand_weight - weight), CONVG_TOLERANCE) <= 0);
            if (succ) break;
        }

        cand_beta = beta;
        cand_nchi = nchi;
        cand_weight = weight;
    }
    if (!succ) return CherenkovCloud();
  
    // check result
    std::vector<CherenkovHit> cand_reduce_hits;
    for (auto&& hit : cand_hits) {
        double nrm = (hit.search_closed_beta(cand_beta) - cand_beta) / width_bta_;
        if (Numc::Compare(std::fabs(nrm), Numc::FIVE<>) > 0) continue;
        cand_reduce_hits.push_back(hit);
    } 
    if (cand_reduce_hits.size() < LMTMIN_CLOUD_HITS) return CherenkovCloud();
    std::sort(cand_reduce_hits.begin(), cand_reduce_hits.end(), CherenkovHit_sort());
    
    std::map<int, int> pmt_maps;
    for (auto&& hit : cand_reduce_hits) {
        if (pmt_maps.count(hit.pmtid()) == 0) pmt_maps[hit.pmtid()] = 1;
        else pmt_maps[hit.pmtid()]++;
    }
    if (pmt_maps.size() < LMTMIN_CLOUD_PMTS) return CherenkovCloud();
   
    short cand_nhit = cand_reduce_hits.size();
    short cand_npmt = pmt_maps.size();

    double cand_cbta = cand_beta * bta_crr_;
    double cand_ndof = cand_weight - Numc::ONE<>;
    double cand_qlt  = Numc::NormQuality(cand_nchi, cand_ndof);
    if (!Numc::Valid(cand_qlt)) return CherenkovCloud();

    double cand_npe = 0.;
    double cand_cnt = 0.;
    double cand_wgt = 0.;
    double cand_all = 0.;
    for (auto&& hit : cand_reduce_hits) {
        hit.set_cnt(Numc::ONE<>);
        cand_npe += hit.npe();
        cand_cnt += hit.cnt();
        cand_wgt += hit.wgt();
        cand_all += hit.wgt() * hit.cnt();
    }
        
    return CherenkovCloud(cand_reduce_hits, cand_beta, cand_cbta, cand_npe, cand_nhit, cand_npmt, cand_cnt, cand_wgt, cand_all, cand_nchi, cand_qlt);
}


std::vector<CherenkovCloud> CherenkovFit::fit_cloud(std::vector<CherenkovHit>& hits) {
    std::vector<CherenkovCloud> cclds;
    if (hits.size() == 0) return cclds;
    
    std::vector<std::array<double, 2>>&& cand_cclds = clustering_cloud(hits); // (bta wgt)
    if (cand_cclds.size() == 0) return cclds;
    
    std::vector<std::vector<CherenkovHit>> cand_chits(cand_cclds.size(), std::vector<CherenkovHit>());
    for (auto&& hit : hits) {
        std::vector<double> cnts;
        for (auto&& ccld : cand_cclds) {
            double nrm = (hit.search_closed_beta(ccld[0]) - ccld[0]) / width_bta_;
            double pdf = std::exp(-Numc::ONE_TO_TWO * nrm * nrm);
            double prb = ccld[1] * pdf;
            
            if (Numc::Compare(prb * hits.size(), CONVG_PROB_SGM50) > 0) cnts.push_back(prb);
            else                                                        cnts.push_back(Numc::ZERO<>);
        }
        double sumprb = std::accumulate(cnts.begin(), cnts.end(), 0.0);
        if (Numc::Compare(sumprb) <= 0) continue;
        for (auto&& cnt : cnts) {
            cnt = (cnt / sumprb);
            if (Numc::Compare(cnt, CONVG_PROB_SGM30) < 0) cnt = 0.0;
        }
        sumprb = std::accumulate(cnts.begin(), cnts.end(), 0.0);
        for (auto&& cnt : cnts) cnt = (cnt / sumprb);
        
        for (int it = 0; it < cand_chits.size(); ++it) {
            if (Numc::EqualToZero(cnts.at(it))) continue;
            std::array<double, 2>&     ccld = cand_cclds.at(it);
            std::vector<CherenkovHit>& chit = cand_chits.at(it); 
            chit.push_back(hit);
            chit.back().search_closed_beta(ccld[0]);
            chit.back().set_cnt(cnts.at(it));
        }
    }

    for (int it = 0; it < cand_cclds.size(); ++it) {
        std::array<double, 2>&     cand_ccld = cand_cclds.at(it);
        std::vector<CherenkovHit>& cand_chit = cand_chits.at(it);

        if (cand_chit.size() < LMTMIN_CLOUD_HITS) continue;
        std::sort(cand_chit.begin(), cand_chit.end(), CherenkovHit_sort());

        std::map<int, int> pmt_maps;
        for (auto&& hit : cand_chit) {
            if (pmt_maps.count(hit.pmtid()) == 0) pmt_maps[hit.pmtid()] = 1;
            else pmt_maps[hit.pmtid()]++;
        }
        if (pmt_maps.size() < LMTMIN_CLOUD_PMTS) continue;
        
        short nhit = cand_chit.size();
        short npmt = pmt_maps.size();

        double beta = cand_ccld[0];
        double cbta = cand_ccld[0] * bta_crr_;

        double npe = 0.;
        double cnt = 0.;
        double wgt = 0.;
        double all = 0.;
        double chisq = 0.;
        for (auto&& hit : cand_chit) {
            npe += hit.npe();
            cnt += hit.cnt();
            wgt += hit.wgt();
            all += hit.wgt() * hit.cnt();
            
            double dlt = hit.search_closed_beta(cand_ccld[0]) - cand_ccld[0];
            std::array<long double, 3>&& minib = pdf_bta_.minimizer(dlt);
            chisq += (hit.wgt() * hit.cnt()) * (minib[0] * minib[0]);
        }

        double ndof = (all - Numc::ONE<>);
        double nchi = chisq / ndof;
        double qlt  = Numc::NormQuality(nchi, ndof);
        if (!Numc::Valid(ndof) || Numc::Compare(ndof) <= 0) continue;
        if (!Numc::Valid(nchi)) continue;
        if (!Numc::Valid(qlt)) continue;
        
        cclds.push_back(CherenkovCloud(cand_chit, beta, cbta, npe, nhit, npmt, cnt, wgt, all, nchi, qlt));
    }
    if (cclds.size() > 1) std::sort(cclds.begin(), cclds.end(), CherenkovCloud_sort());

    return cclds;
}


std::vector<std::array<double, 2>> CherenkovFit::clustering_cloud(std::vector<CherenkovHit>& hits) {
    std::vector<std::array<double, 2>> cclds; // (bta wgt)
    if (hits.size() == 0) return cclds;

    cclds = std::move(clustering_evolve_cloud(hits)); // (bta wgt)
    if (cclds.size() == 0) return cclds;

    bool succ = false;
    for (int iter = 0; iter <= LMTMAX_ITER; ++iter) {
        std::vector<std::array<double, 4>> mgaus(cclds.size(), std::array<double, 4>({ 0.0, 0.0, 0.0, 0.0 })); // (wgt grd cv res)
        double weight = Numc::ZERO<>;

        for (auto&& hit : hits) {
            double sumprb = Numc::ZERO<>;
            std::vector<std::array<double, 3>> elms; // (wgt grd cv)
            for (auto&& ccld : cclds) {
                if (Numc::Compare(ccld[1]) <= 0) {
                    elms.push_back(std::array<double, 3>({ 0.0, 0.0, 0.0 }));
                    continue;
                }
                double dlt = hit.search_closed_beta(ccld[0]) - ccld[0];
                double nrm = dlt / width_bta_;
                double pdf = std::exp(-Numc::ONE_TO_TWO * nrm * nrm);
                double prb = ccld[1] * pdf;
        
                std::array<long double, 3>&& minib = pdf_bta_.minimizer(dlt);
                double grdB = (Numc::NEG<> * minib[2] * minib[1]);
                double cvBB = (minib[2] * minib[2]);
                
                elms.push_back(std::array<double, 3>({ prb, grdB, cvBB }));
                sumprb += prb;
            }
            if (!Numc::Valid(sumprb) || Numc::Compare(sumprb) <= 0) continue;
            for (auto&& elm : elms) elm.at(0) = hit.wgt() * (elm.at(0) / sumprb);
            weight += hit.wgt();

            for (int it = 0; it < elms.size(); ++it) {
                std::array<double, 3>& elm  = elms.at(it);
                std::array<double, 4>& gaus = mgaus.at(it);
                if (Numc::Compare(elm[0]) <= 0) continue;
                gaus[0] += elm[0];
                gaus[1] += elm[0] * elm[1];
                gaus[2] += elm[0] * elm[2];
            }
        }
        if (!Numc::Valid(weight) || Numc::Compare(weight, Numc::ONE<>) <= 0) break;
        
        for (auto&& gaus : mgaus) {
            if (Numc::Compare(gaus[0]) <= 0) continue;
            gaus[3] = gaus[1] / gaus[2];
            gaus[0] = gaus[0] / weight;
        }

        std::vector<std::array<double, 2>> new_cclds; // (bta wgt)
        for (int it = 0; it < cclds.size(); ++it) {
            double bta = cclds.at(it)[0] - mgaus.at(it)[3];
            double wgt = mgaus.at(it)[0];
            new_cclds.push_back(std::array<double, 2>({ bta, wgt }));
        }

        if (iter >= LMTM_ITER) {
            bool is_convg = true;
            for (int it = 0; it < cclds.size(); ++it) {
                auto&& gaus = new_cclds.at(it);
                auto&& ccld = cclds.at(it);
                bool match = (Numc::Compare(std::fabs(gaus[0] - ccld[0]), CONVG_TOLERANCE) <= 0 && 
                              Numc::Compare(std::fabs(gaus[1] - ccld[1]), CONVG_TOLERANCE) <= 0);
                if (!match) { is_convg = false; break; }
            }
            if (is_convg) { succ = true; break; }
        }
        cclds = new_cclds;
    }
    
    if (!succ) cclds.clear();
    return cclds;
}


std::vector<std::array<double, 2>> CherenkovFit::clustering_evolve_cloud(std::vector<CherenkovHit>& hits) {
    // bayesian information criterion (BIC)
    std::vector<std::array<double, 2>> cclds; // (bta wgt)
    if (hits.size() == 0) return cclds;

    std::vector<double> seed_beta;
    for (auto&& hit : hits) {
        if (hit.type()%2 == 1) seed_beta.push_back(hit.dbta());
        if (hit.type()/2 == 1) seed_beta.push_back(hit.rbta());
    }
    std::sort(seed_beta.begin(), seed_beta.end());

    for (int it = 0; it < seed_beta.size(); ++it) {
        double init_bta = seed_beta.at(it);
        double init_wgt = Numc::ONE<> / static_cast<double>(seed_beta.size());
        cclds.push_back(std::array<double, 2>({ init_bta, init_wgt })); // (bta wgt)
    }

    bool succ = false;
    for (int iter = 1; iter <= LMTMAX_ITER; ++iter) {
        std::vector<std::array<double, 2>> mgaus(cclds.size(), std::array<double, 2>({ 0.0, 0.0 })); // (bta wgt)
        double weight = Numc::ZERO<>;

        for (auto&& hit : hits) {
            double sumprb = Numc::ZERO<>;
            std::vector<std::array<double, 2>> elms; // (bta wgt)
            for (auto&& ccld : cclds) {
                if (Numc::Compare(ccld[1]) <= 0) {
                    elms.push_back(std::array<double, 2>({ 0.0, 0.0 }));
                    continue;
                }
                double bta = hit.search_closed_beta(ccld[0]);
                double nrm = (bta - ccld[0]) / width_bta_;
                double pdf = std::exp(-Numc::ONE_TO_TWO * nrm * nrm);
                double prb = ccld[1] * pdf;
                
                elms.push_back(std::array<double, 2>({ bta, prb }));
                sumprb += prb;
            }
            if (!Numc::Valid(sumprb) || Numc::Compare(sumprb) <= 0) continue;
            for (auto&& elm : elms) elm.at(1) = hit.wgt() * (elm.at(1) / sumprb);
            weight += hit.wgt();

            for (int it = 0; it < elms.size(); ++it) {
                std::array<double, 2>& elm  = elms.at(it);
                std::array<double, 2>& gaus = mgaus.at(it);
                if (Numc::Compare(elm.at(1)) <= 0) continue;
                
                gaus[0] += elm[1] * elm[0];
                gaus[1] += elm[1];
            }
        }
        if (!Numc::Valid(weight) || Numc::Compare(weight, Numc::ONE<>) <= 0) break;

        for (auto&& gaus : mgaus) {
            if (Numc::Compare(gaus[1]) <= 0) continue;
            gaus[0] = gaus[0] / gaus[1];
            gaus[1] = gaus[1] / weight;
        }

        // remove empty
        if (iter >= LMTL_ITER) {
            for (int it = mgaus.size()-1; it >= 0; --it) {
                if (Numc::Compare(mgaus.at(it)[1], CONVG_EPSILON) < 0)
                    mgaus.erase(mgaus.begin() + it);
            }
            if (mgaus.size() == 0) break;
        }
        bool has_remove_empty = (mgaus.size() != cclds.size());
        if (mgaus.size() >= 2) std::sort(mgaus.begin(), mgaus.end());
        
        // merge
        bool has_merge = false;
        if (iter >= LMTM_ITER && !has_remove_empty) {
            std::vector<std::array<int, 2>> sets;
            for (int it = 0; it < mgaus.size(); ++it) {
                auto&& gausA = (it == 0) ? mgaus.at(it) : mgaus.at(it-1);
                auto&& gausB = mgaus.at(it);
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
            has_merge = (mgaus.size() != sets.size());
            if (has_merge) {
                std::vector<std::array<double, 2>> newmg;
                for (auto&& set : sets) {
                    if (set.at(0) == set.at(1)) {
                        newmg.push_back(mgaus.at(set.at(0))); 
                        continue; 
                    }
                    double newbta = Numc::ZERO<>;
                    double newwgt = Numc::ZERO<>;
                    for (int it = set.at(0); it <= set.at(1); ++it) {
                        auto&& gaus = mgaus.at(it);
                        newbta += gaus[1] * gaus[0];
                        newwgt += gaus[1];
                    }
                    newbta = (newbta / newwgt);
                    newmg.push_back(std::array<double, 2>({ newbta, newwgt }));
                }
                mgaus = std::move(newmg);
            }
        }
        double normalized = Numc::ZERO<>;
        for (auto&& gaus : mgaus) normalized += gaus[1];
        for (auto&& gaus : mgaus) gaus[1] = (gaus[1] / normalized);
        if (Numc::Compare(normalized) <= 0) break;

        if (iter >= LMTU_ITER && !has_remove_empty && !has_merge) {
            bool is_stable = true;
            for (int it = 0; it < mgaus.size(); ++it) {
                auto&& gaus = mgaus.at(it);
                auto&& ccld = cclds.at(it);
                bool match = (Numc::Compare(std::fabs(gaus[0] - ccld[0]), CONVG_TOLERANCE) <= 0 && 
                              Numc::Compare(std::fabs(gaus[1] - ccld[1]), CONVG_TOLERANCE) <= 0);
                if (!match) { is_stable = false; break; }
            }
            if (is_stable) { succ = true; break; }
        }
        cclds = mgaus;
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
