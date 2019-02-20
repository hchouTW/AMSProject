#ifndef __TRACKLibs_CherenkovMeas_C__
#define __TRACKLibs_CherenkovMeas_C__


#include "Sys.h"
#include "Math.h"
#include "CherenkovMeas.h"


namespace TrackSys {

        
CherenkovFit::CherenkovFit(const std::vector<CherenkovHit>& args_hits, const std::array<double, 2>& scan_bta, const std::array<double, 4>& scan_npe, const std::array<double, 5>& args_bta) : scan_bta_(scan_bta), scan_npe_(scan_npe), args_bta_(args_bta) {
    std::vector<CherenkovHit> hits;
    for (auto&& hit : args_hits) { if (hit.status() && (hit.dist() > LMTL_DIST)) hits.push_back(hit); }
    if (hits.size() <= 2 || !check()) { clear(); return; }
 
    timer_.start();
    
    std::vector<std::pair<double, std::vector<CherenkovHit>>>&& clusters = clustering(hits);
    if (clusters.size() == 0) { clear(); return; }
    
    for (auto&& cluster : clusters) {
        CherenkovCls&& chcls = physicalFit(cluster);
        if (!chcls.status()) continue;
        clss_.push_back(chcls);
    }
    if (clss_.size() == 0) { clear(); return; }
    if (clss_.size() > 1) std::sort(clss_.begin(), clss_.end(), CherenkovCls_sort());

    for (auto&& cls : clss_) {
    for (auto&& hit : cls.hits()) {
        hits_.push_back(hit);
    }}
    if (hits_.size() <= 2) { clear(); return; }
    
    timer_.stop();
    succ_ = true;
}


void CherenkovFit::clear() {
    succ_ = false;
    hits_.clear();
    clss_.clear();

    scan_bta_ = std::array<double, 2>({0, 0});
    scan_npe_ = std::array<double, 4>({0, 0, 0, 0});
    args_bta_ = std::array<double, 5>({0, 0, 0, 0, 0});
    
    mgscan_ = std::move(MultiGaus());
    mgfit_  = std::move(MultiGaus());

    timer_.clear();
}


bool CherenkovFit::check() {
    if (Numc::Compare(scan_bta_[0]) <= 0 || 
        Numc::Compare(scan_bta_[1]) < 0) return false;
    if (Numc::Compare(scan_npe_[0]) <= 0 || 
        Numc::Compare(scan_npe_[1]) <= 0 || 
        Numc::Compare(scan_npe_[2]) <= 0 || 
        Numc::Compare(scan_npe_[3]) < 0) return false;
    if (Numc::Compare(args_bta_[0]) <= 0 || 
        Numc::Compare(args_bta_[1]) <= 0 || 
        Numc::Compare(args_bta_[2]) <= 0 || 
        Numc::Compare(args_bta_[3]) <= 0 || 
        Numc::Compare(args_bta_[4]) < 0) return false;
    mgscan_ = std::move(MultiGaus(Robust(Robust::Opt::OFF), scan_bta_[0]));
    mgfit_  = std::move(MultiGaus(Robust(Robust::Opt::OFF), args_bta_[0], args_bta_[1], args_bta_[2], args_bta_[3]));
    return true;
}


CherenkovCls CherenkovFit::physicalFit(const std::pair<double, std::vector<CherenkovHit>>& param) {
    double beta = param.first;
    std::vector<CherenkovHit> hits = param.second;
    if (Numc::Compare(beta) <= 0 || hits.size() == 0) return CherenkovCls();

    bool succ = false;
    for (short iter = 1; iter <= LMTU_ITER; ++iter) {
        double grdB = Numc::ZERO<>;
        double cvBB = Numc::ZERO<>;

        short cnt_nhit = 0;
        for (auto&& hit : hits) {
            double dlt_b = hit.search_closed_beta(beta) - beta;
            double bta_s = MultiGaus::Func(dlt_b, args_bta_[0], args_bta_[1], args_bta_[2], args_bta_[3]);
            double npe_s = LandauGaus::Func(hit.npe(), scan_npe_[0], scan_npe_[1], scan_npe_[2]);
            double bta_w = bta_s * (Numc::ONE<> + args_bta_[4]) / (bta_s + args_bta_[4]);
            double npe_w = npe_s * (Numc::ONE<> + scan_npe_[3]) / (npe_s + scan_npe_[3]);
            double wgt   = (bta_w * npe_w);
        
            std::array<long double, 3>&& minib = mgfit_.minimizer(dlt_b);
            grdB += wgt * (Numc::NEG<> * minib[2] * minib[1]);
            cvBB += wgt * (minib[2] * minib[2]);
            
            cnt_nhit++;
        }
        if (cnt_nhit != hits.size()) break;
        
        double rslB = grdB / cvBB;
        if (!Numc::Valid(rslB)) break;

        double newbta = beta - rslB;
        if (!Numc::Valid(newbta) || Numc::Compare(newbta) <= 0) break; 

        long double dlt = std::fabs(newbta / beta - Numc::ONE<long double>);
        long double rat = std::fabs(newbta - beta) / (newbta + beta);
        if (!Numc::Valid(dlt) || !Numc::Valid(rat)) break;

        if (iter >= Numc::TWO<short> * LMTL_ITER && dlt < CONVG_EPSILON && rat < CONVG_TOLERANCE) { 
            succ = true;
            beta = newbta; 
            break; 
        }
        else beta = newbta;
    }
    if (!succ) return CherenkovCls();
    
    short  ndof = hits.size() - Numc::ONE<short>;
    double nchi = Numc::ZERO<>;
   
    double cnt = Numc::ZERO<>;
    double nos = Numc::ZERO<>;
    double eft_ndof = -Numc::ONE<double>;
    double eft_nchi = Numc::ZERO<>;

    double sum_bw  = Numc::ZERO<>;
    double sum_sw  = Numc::ZERO<>;
    double sum_bsw = Numc::ZERO<>;
    for (auto&& hit : hits) {
        double dlt_b = hit.search_closed_beta(beta) - beta;
        double bta_s = MultiGaus::Func(dlt_b, args_bta_[0], args_bta_[1], args_bta_[2], args_bta_[3]);
        double npe_s = LandauGaus::Func(hit.npe(), scan_npe_[0], scan_npe_[1], scan_npe_[2]);
        double bta_w = bta_s * (Numc::ONE<> + args_bta_[4]) / (bta_s + args_bta_[4]);
        double npe_w = npe_s * (Numc::ONE<> + scan_npe_[3]) / (npe_s + scan_npe_[3]);
        double wgt   = (bta_w * npe_w);
    
        sum_bw  += bta_w;
        sum_sw  += npe_w;
        sum_bsw += wgt;

        std::array<long double, 3>&& minib = mgfit_.minimizer(dlt_b);
        nchi += minib[0] * minib[0];

        cnt += wgt;
        nos += (Numc::ONE<double> - wgt);

        eft_ndof += wgt;
        eft_nchi += wgt * minib[0] * minib[0];
    }
    if (!Numc::Valid(sum_bw)  || Numc::Compare(sum_bw)  <= 0) return CherenkovCls();
    if (!Numc::Valid(sum_sw)  || Numc::Compare(sum_sw)  <= 0) return CherenkovCls();
    if (!Numc::Valid(sum_bsw) || Numc::Compare(sum_bsw) <= 0) return CherenkovCls();
    
    nchi = (nchi / static_cast<double>(ndof));
    double quality = Numc::NormQuality(nchi, ndof); 
    if (!Numc::Valid(nchi) || !Numc::Valid(quality)) return CherenkovCls();
    
    double eta = (eft_nchi / eft_ndof);
    if (!Numc::Valid(eta) || Numc::Compare(eta) <= 0) return CherenkovCls();

    double compact_c = -Numc::TWO<> * std::log(sum_bsw / static_cast<double>(hits.size())); // compact of cluster
    double compact_b = -Numc::TWO<> * std::log(sum_bsw / sum_sw); // compact of beta
    double compact_s = -Numc::TWO<> * std::log(sum_bsw / sum_bw); // compact of signal(npe)

    if (!Numc::Valid(compact_c) || Numc::Compare(compact_c) <= 0) compact_c = Numc::ZERO<>;
    if (!Numc::Valid(compact_b) || Numc::Compare(compact_b) <= 0) compact_b = Numc::ZERO<>;
    if (!Numc::Valid(compact_s) || Numc::Compare(compact_s) <= 0) compact_s = Numc::ZERO<>;

    return CherenkovCls(hits, beta, cnt, nos, eta, ndof, nchi, quality, compact_c, compact_b, compact_s);
}


std::vector<std::pair<double, std::vector<CherenkovHit>>> CherenkovFit::clustering(std::vector<CherenkovHit>& hits) {
    std::vector<std::pair<double, std::vector<CherenkovHit>>> clusters;
    if (hits.size() == 0) return clusters;
    
    std::vector<std::pair<std::array<double, 2>, CherenkovHit*>> rawcls;
    for (auto&& hit : hits) {
        std::array<double, 2> candcls({ 0., 0. });
        for (short mode = 0; mode <= 1; ++mode) {
            if (mode == 0 && hit.type() != 1 && hit.type() != 3) continue;
            if (mode == 1 && hit.type() != 2 && hit.type() != 3) continue;
            double sbta = (mode == 0) ? hit.dbta() : hit.rbta();
            std::array<double, 2>&& rsl = clustering_evolve(hits, sbta);
            if (Numc::Compare(rsl[0]) <= 0) continue;
            if (Numc::Compare(candcls[0]) <= 0 || rsl[1] < candcls[1]) candcls = rsl;
        }
        if (Numc::Compare(candcls[0]) > 0) rawcls.push_back(std::make_pair(candcls, &hit));
    }
    if (rawcls.size() == 0) return clusters;
    std::sort(rawcls.begin(), rawcls.end());
    std::reverse(rawcls.begin(), rawcls.end());

    std::vector<std::array<short, 2>> presets;
    for (short it = 0; it < rawcls.size(); ++it) {
        double tmpb = (it == 0) ? 0.0 : (rawcls.at(it-1).first)[0];
        double tmpa = (rawcls.at(it).first)[0];
        if (it != 0 && std::fabs(tmpa - tmpb) < CONVG_CLUSTERING) { presets.back().at(1) = it; continue; }
        presets.push_back(std::array<short, 2>({it, it}));
    }
    std::vector<std::array<short, 2>> sets;
    for (auto&& set : presets) {
        int size = set.at(1) - set.at(0) + 1;
        if (size <= 2) continue;
        sets.push_back(set);
    }
    if (sets.size() == 0) return clusters;

    for (auto&& set : sets) {
        double beta = 0.0;
        std::vector<CherenkovHit> hits;
        for (short it = set.at(0); it <= set.at(1); ++it) {
            std::array<double, 2>& cls = (rawcls.at(it).first);
            CherenkovHit* hit = (rawcls.at(it).second);
            hits.push_back(*hit);
            beta += cls[0];
        }
        if (hits.size() == 0) continue;
        if (!Numc::Valid(beta)) continue;
        beta = (beta / static_cast<double>(hits.size()));
        clusters.push_back(std::make_pair(beta, hits));
    }

    return clusters;
}
       

std::array<double, 2> CherenkovFit::clustering_evolve(std::vector<CherenkovHit>& hits, double sbta) {
    std::array<double, 2> result({ 0., 0. });
    if (Numc::Compare(sbta) <= 0) return result;
    if (hits.size() == 0) return result;

    bool succ = false;
    double beta = sbta;
    double compact = Numc::ZERO<>;
    for (short iter = 1; iter <= LMTU_ITER; ++iter) {
        double chic = Numc::ZERO<>;
        double grdB = Numc::ZERO<>;
        double cvBB = Numc::ZERO<>;

        short cnt_nhit = 0;
        for (auto&& hit : hits) {
            double dlt_b = hit.search_closed_beta(beta) - beta;
            double bta_s = MultiGaus::Func(dlt_b, scan_bta_[0]);
            double npe_s = LandauGaus::Func(hit.npe(), scan_npe_[0], scan_npe_[1], scan_npe_[2]);
            double bta_w = bta_s * (Numc::ONE<> + scan_bta_[1]) / (bta_s + scan_bta_[1]);
            double npe_w = npe_s * (Numc::ONE<> + scan_npe_[3]) / (npe_s + scan_npe_[3]);
            double wgt   = (bta_w * npe_w);
            chic += wgt;

            std::array<long double, 3>&& minib = mgscan_.minimizer(dlt_b);
            grdB += wgt * (Numc::NEG<> * minib[2] * minib[1]);
            cvBB += wgt * (minib[2] * minib[2]);
            
            cnt_nhit++;
        }
        if (cnt_nhit != hits.size()) return result;
        
        chic = -Numc::TWO<> * std::log(chic / static_cast<double>(hits.size()));
        if (!Numc::Valid(chic) || Numc::Compare(chic) <= 0) chic = Numc::ZERO<>;

        double rslB = grdB / cvBB;
        if (!Numc::Valid(rslB)) return result;

        double newbta = beta - rslB;
        if (!Numc::Valid(newbta) || Numc::Compare(newbta) <= 0) return result;

        long double dlt = std::fabs(newbta / beta - Numc::ONE<long double>);
        long double rat = std::fabs(newbta - beta) / (newbta + beta);
        if (!Numc::Valid(dlt) || !Numc::Valid(rat)) return result;

        if (iter >= Numc::TWO<short> * LMTL_ITER && dlt < CONVG_EPSILON && rat < CONVG_TOLERANCE) { 
            succ = true;
            beta = newbta;
            compact = chic;
            break; 
        }
        else beta = newbta;
    }
    if (!succ) return result;

    result[0] = beta;
    result[1] = compact;
    return result; 
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
