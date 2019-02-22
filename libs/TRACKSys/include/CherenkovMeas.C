#ifndef __TRACKLibs_CherenkovMeas_C__
#define __TRACKLibs_CherenkovMeas_C__


#include "Sys.h"
#include "Math.h"
#include "CherenkovMeas.h"


namespace TrackSys {

        
CherenkovFit::CherenkovFit(const std::vector<CherenkovHit>& args_hits, const std::array<double, 4>& args_bta, const std::array<double, 2>& scan_bta, const std::array<double, 4>& scan_npe) {
    args_bta_ = std::array<double, 4>({ args_bta[0], args_bta[1], args_bta[2], args_bta[3] });
    scan_bta_sig_ = scan_bta[0];
    scan_bta_nos_ = scan_bta[1];
    scan_npe_sig_ = std::array<double, 3>({ scan_npe[0], scan_npe[1], scan_npe[2] });
    scan_npe_nos_ = scan_npe[3];
    
    timer_.start();
    if (args_hits.size() < 2 || !check()) { clear(); return; }
    succ_ = fit(args_hits);    
    timer_.stop();

    if (!succ_) { clear(); return; }
}


void CherenkovFit::clear() {
    succ_ = false;
    hits_.clear();
    clss_.clear();

    args_bta_ = std::array<double, 4>({0, 0, 0, 0});
    pdf_bta_  = std::move(MultiGaus());

    scan_bta_sig_ = 0;
    scan_bta_nos_ = 0;
    pdf_scan_bta_ = std::move(MultiGaus());
    
    scan_npe_sig_ = std::array<double, 3>({0, 0, 0});
    scan_npe_nos_ = Numc::ZERO<>;

    convg_epsilon_   = 0;
    convg_tolerance_ = 0;
    convg_closed_    = 0;

    timer_.clear();
}


bool CherenkovFit::check() {
    if (Numc::Compare(args_bta_[0]) <= 0 || 
        Numc::Compare(args_bta_[1]) <= 0 || 
        Numc::Compare(args_bta_[2]) <= 0 || 
        Numc::Compare(args_bta_[3]) <= 0) return false; 
    if (Numc::Compare(scan_bta_sig_) <= 0 || 
        Numc::Compare(scan_bta_nos_) <  0) return false;
    if (Numc::Compare(scan_npe_sig_[0]) <= 0 || 
        Numc::Compare(scan_npe_sig_[1]) <= 0 || 
        Numc::Compare(scan_npe_sig_[2]) <= 0 || 
        Numc::Compare(scan_npe_nos_)     <  0) return false;
    pdf_bta_      = std::move(MultiGaus(Robust::Option(Robust::Opt::ON, 4.0L, 0.5L), args_bta_[0], args_bta_[1], args_bta_[2], args_bta_[3]));
    pdf_scan_bta_ = std::move(MultiGaus(Robust::Option(Robust::Opt::ON, 3.0L, 0.5L), scan_bta_sig_));
    convg_epsilon_   = scan_bta_sig_ * CONVG_EPSILON;
    convg_tolerance_ = scan_bta_sig_ * CONVG_TOLERANCE;
    convg_closed_    = scan_bta_sig_ * CONVG_CLOSED;
    return true;
}

        
bool CherenkovFit::fit(const std::vector<CherenkovHit>& hits) {
    std::vector<short> hitsID;
    clss_.clear();
    hits_.clear();

    int iter = 0;
    const int max_iter = 100;
    while(true) {
        iter++;
        std::vector<CherenkovHit> cand_hits;
        for (auto&& hit : hits) {
            if (!hit.status()) continue;
            bool has = false;
            for (auto&& idx : hitsID) {
                if (idx == hit.index()) { has = true; break; }
            }
            if (has) continue;
            cand_hits.push_back(hit);
        }
        if (cand_hits.size() == 0) break;

        std::vector<std::pair<double, std::vector<CherenkovHit>>>&& clusters = clustering(cand_hits);
        if (clusters.size() == 0) break;
        
        std::vector<CherenkovCls> clss;
        for (auto&& cluster : clusters) {
            CherenkovCls&& chcls = physicalFit(cluster);
            if (!chcls.status()) continue;
            clss.push_back(chcls);
        }
        if (clss.size() == 0) break;
        if (clss.size() > 1) std::sort(clss.begin(), clss.end(), CherenkovCls_sort());

        CherenkovCls& main_cls = clss.at(0);
        clss_.push_back(main_cls);
        for (auto&& hit : main_cls.hits()) {
            hits_.push_back(hit);
            hitsID.push_back(hit.index());
        }

        if (iter >= max_iter) break;
    }
    if (clss_.size() > 1) std::sort(clss_.begin(), clss_.end(), CherenkovCls_sort());

    bool succ = (clss_.size() >= 1 && hits_.size() >= 1);
    return succ;
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
            double dlt_b  = hit.search_closed_beta(beta) - beta;
            double npe_s  = LandauGaus::Func(hit.npe(), scan_npe_sig_[0], scan_npe_sig_[1], scan_npe_sig_[2]);
            double wgtnpe = npe_s * (Numc::ONE<> + scan_npe_nos_) / (npe_s + scan_npe_nos_);
            
            std::array<long double, 3>&& minib = pdf_bta_.minimizer(dlt_b);
            grdB += wgtnpe * (Numc::NEG<> * minib[2] * minib[1]);
            cvBB += wgtnpe * (minib[2] * minib[2]);
            
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

        if (iter >= Numc::TWO<short> * LMTL_ITER && dlt < convg_epsilon_ && rat < convg_tolerance_) { 
            succ = true;
            beta = newbta; 
            break; 
        }
        else beta = newbta;
    }
    if (!succ) return CherenkovCls();
    
    short  ndof = hits.size() - Numc::ONE<short>;
    double nchi = Numc::ZERO<>;
    
    double sig_s2 = Numc::ZERO<>;
    double sig_s3 = Numc::ZERO<>;
    double sig_s4 = Numc::ZERO<>;
    double sig_s5 = Numc::ZERO<>;
    
    double nos_s2 = Numc::ZERO<>;
    double nos_s3 = Numc::ZERO<>;
    double nos_s4 = Numc::ZERO<>;
    double nos_s5 = Numc::ZERO<>;
    
    double scan_s2 = Numc::ZERO<>;
    double scan_s3 = Numc::ZERO<>;
    double scan_s4 = Numc::ZERO<>;
    double scan_s5 = Numc::ZERO<>;

    for (auto&& hit : hits) {
        double dlt_b  = hit.search_closed_beta(beta) - beta;
        double npe_s  = LandauGaus::Func(hit.npe(), scan_npe_sig_[0], scan_npe_sig_[1], scan_npe_sig_[2]);
        double wgtnpe = npe_s * (Numc::ONE<> + scan_npe_nos_) / (npe_s + scan_npe_nos_);
        
        std::array<long double, 3>&& minib = pdf_bta_.minimizer(dlt_b);
        nchi += minib[0] * minib[0];
        
        double bta_s     = MultiGaus::Func(dlt_b, args_bta_[0], args_bta_[1], args_bta_[2], args_bta_[3]);
        double wgtbta_s2 = bta_s * (Numc::ONE<> + PROB_SIGMA2) / (bta_s + PROB_SIGMA2);
        double wgtbta_s3 = bta_s * (Numc::ONE<> + PROB_SIGMA3) / (bta_s + PROB_SIGMA3);
        double wgtbta_s4 = bta_s * (Numc::ONE<> + PROB_SIGMA4) / (bta_s + PROB_SIGMA4);
        double wgtbta_s5 = bta_s * (Numc::ONE<> + PROB_SIGMA5) / (bta_s + PROB_SIGMA5);
        double wgt_s2    = (wgtbta_s2 * wgtnpe);
        double wgt_s3    = (wgtbta_s3 * wgtnpe);
        double wgt_s4    = (wgtbta_s4 * wgtnpe);
        double wgt_s5    = (wgtbta_s5 * wgtnpe);

        sig_s2 += wgt_s2;
        sig_s3 += wgt_s3;
        sig_s4 += wgt_s4;
        sig_s5 += wgt_s5;

        nos_s2 += Numc::ONE<> - wgt_s2;
        nos_s3 += Numc::ONE<> - wgt_s3;
        nos_s4 += Numc::ONE<> - wgt_s4;
        nos_s5 += Numc::ONE<> - wgt_s5;
        
        scan_s2 += wgt_s2;
        scan_s3 += wgt_s3;
        scan_s4 += wgt_s4;
        scan_s5 += wgt_s5;
    }
    
    nchi = (nchi / static_cast<double>(ndof));
    double quality = Numc::NormQuality(nchi, ndof); 
    if (!Numc::Valid(nchi) || !Numc::Valid(quality)) return CherenkovCls();
    
    if (!Numc::Valid(sig_s2) || !Numc::Valid(sig_s3) || !Numc::Valid(sig_s4) || !Numc::Valid(sig_s5)) return CherenkovCls(); 
    if (!Numc::Valid(nos_s2) || !Numc::Valid(nos_s3) || !Numc::Valid(nos_s4) || !Numc::Valid(nos_s5)) return CherenkovCls(); 

    scan_s2 = (scan_s2 / static_cast<double>(hits.size()));
    scan_s3 = (scan_s3 / static_cast<double>(hits.size()));
    scan_s4 = (scan_s4 / static_cast<double>(hits.size()));
    scan_s5 = (scan_s5 / static_cast<double>(hits.size()));
    if (!Numc::Valid(scan_s2) || !Numc::Valid(scan_s3) || !Numc::Valid(scan_s4) || !Numc::Valid(scan_s5)) return CherenkovCls(); 

    double compact_s2 = -Numc::TWO<> * std::log(scan_s2); // compact of cluster in sigma2
    double compact_s3 = -Numc::TWO<> * std::log(scan_s3); // compact of cluster in sigma3
    double compact_s4 = -Numc::TWO<> * std::log(scan_s4); // compact of cluster in sigma4
    double compact_s5 = -Numc::TWO<> * std::log(scan_s5); // compact of cluster in sigma5
    if (!Numc::Valid(compact_s2) || Numc::Compare(compact_s2) <= 0) compact_s2 = Numc::ZERO<>;
    if (!Numc::Valid(compact_s3) || Numc::Compare(compact_s3) <= 0) compact_s3 = Numc::ZERO<>;
    if (!Numc::Valid(compact_s4) || Numc::Compare(compact_s4) <= 0) compact_s4 = Numc::ZERO<>;
    if (!Numc::Valid(compact_s5) || Numc::Compare(compact_s5) <= 0) compact_s5 = Numc::ZERO<>;

    std::array<double, 4> sig({sig_s2, sig_s3, sig_s4, sig_s5});
    std::array<double, 4> nos({nos_s2, nos_s3, nos_s4, nos_s5});
    std::array<double, 4> compact({compact_s2, compact_s3, compact_s4, compact_s5});

    double signal_s2s4 = (Numc::Compare(sig_s4) > 0) ? (sig_s2 / sig_s4) : Numc::ZERO<>;
    double signal_s3s4 = (Numc::Compare(sig_s4) > 0) ? (sig_s3 / sig_s4) : Numc::ZERO<>;
    double noise_s2s4 = (Numc::Compare(nos_s2) > 0) ? (nos_s4 / nos_s2) : Numc::ZERO<>;
    double noise_s3s4 = (Numc::Compare(nos_s3) > 0) ? (nos_s4 / nos_s3) : Numc::ZERO<>;
    if (!Numc::Valid(signal_s2s4) || !Numc::Valid(signal_s3s4)) return CherenkovCls(); 
    if (!Numc::Valid(noise_s2s4) || !Numc::Valid(noise_s3s4)) return CherenkovCls(); 

    double sn_s2s4 = signal_s2s4 / (signal_s2s4 + noise_s2s4);
    double sn_s3s4 = signal_s3s4 / (signal_s3s4 + noise_s3s4);
    if (!Numc::Valid(sn_s2s4) || !Numc::Valid(sn_s3s4)) return CherenkovCls(); 

    return CherenkovCls(hits, beta, 
                        ndof, nchi, quality, 
                        sig, nos, compact, 
                        signal_s2s4, signal_s3s4, 
                        noise_s2s4, noise_s3s4, 
                        sn_s2s4, sn_s3s4);
}


std::vector<std::pair<double, std::vector<CherenkovHit>>> CherenkovFit::clustering(std::vector<CherenkovHit>& hits) {
    std::vector<std::pair<double, std::vector<CherenkovHit>>> clusters;
    if (hits.size() == 0) return clusters;
    
    std::vector<std::pair<std::array<double, 2>, CherenkovHit*>> rawcls;
    for (auto&& hit : hits) {
        for (short mode = 0; mode <= 1; ++mode) {
            if (mode == 0 && hit.type() != 1 && hit.type() != 3) continue;
            if (mode == 1 && hit.type() != 2 && hit.type() != 3) continue;
            double sbta = (mode == 0) ? hit.dbta() : hit.rbta();
            std::array<double, 2>&& rsl = clustering_evolve(hits, sbta);
            if (Numc::Compare(rsl[0]) <= 0) continue;
            rawcls.push_back(std::make_pair(rsl, &hit));
        }
    }
    if (rawcls.size() == 0) return clusters;
    std::sort(rawcls.begin(), rawcls.end());
    std::reverse(rawcls.begin(), rawcls.end());

    std::vector<std::array<short, 2>> presets;
    for (short it = 0; it < rawcls.size(); ++it) {
        double tmpb = (it == 0) ? 0.0 : (rawcls.at(it-1).first)[0];
        double tmpa = (rawcls.at(it).first)[0];
        if (it != 0 && std::fabs(tmpa - tmpb) < convg_closed_) { presets.back().at(1) = it; continue; }
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
        double cnt  = Numc::ZERO<>;
        double grdB = Numc::ZERO<>;
        double cvBB = Numc::ZERO<>;

        short cnt_nhit = 0;
        for (auto&& hit : hits) {
            double dlt_b  = hit.search_closed_beta(beta) - beta;
            double bta_s  = MultiGaus::Func(dlt_b, scan_bta_sig_);
            double npe_s  = LandauGaus::Func(hit.npe(), scan_npe_sig_[0], scan_npe_sig_[1], scan_npe_sig_[2]);
            double wgtbta = bta_s * (Numc::ONE<> + scan_bta_nos_) / (bta_s + scan_bta_nos_);
            double wgtnpe = npe_s * (Numc::ONE<> + scan_npe_nos_) / (npe_s + scan_npe_nos_);
            double wgt    = (wgtbta * wgtnpe);

            std::array<long double, 3>&& minib = pdf_scan_bta_.minimizer(dlt_b);
            grdB += wgt * (Numc::NEG<> * minib[2] * minib[1]);
            cvBB += wgt * (minib[2] * minib[2]);
            
            cnt += wgt;
            cnt_nhit++;
        }
        if (cnt_nhit != hits.size()) return result;
        
        double chic = -Numc::TWO<> * std::log(cnt / static_cast<double>(hits.size()));
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
