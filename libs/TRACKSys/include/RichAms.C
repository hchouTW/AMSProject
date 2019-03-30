#if defined(_PGTRACK_) || defined(__ROOTSHAREDLIBRARY__)
#ifndef __TRACKLibs_RichAms_C__
#define __TRACKLibs_RichAms_C__


#include "Sys.h"
#include "Math.h"
#include "CherenkovMeas.h"
#include "RichAms.h"


namespace TrackSys {


RichHitAms::RichHitAms(RichHitR* hit, double dbta, double rbtaA, double rbtaB, double dist) : RichHitAms() {
    if (hit == nullptr || (dbta <= 0 && rbtaA <= 0 && rbtaB <= 0)) return;
    status_ = true;
    chann_  = hit->Channel;
    pmtid_  = chann_ / NUM_CHANN_IN_2D_PMT;
    pixel_  = chann_ - pmtid_ * NUM_CHANN_IN_2D_PMT;
    locid_[0] = (pixel_ / NUM_CHANN_IN_1D_PMT);
    locid_[1] = (pixel_ % NUM_CHANN_IN_1D_PMT);

    type_  = (dbta > 0) + (rbtaA > 0) * 2 + (rbtaB > 0) * 4;
    dbta_  = (type_&1 == 1) ? dbta  : -1.0;
    rbtaA_ = (type_&2 == 2) ? rbtaA : -1.0;
    rbtaB_ = (type_&4 == 4) ? rbtaB : -1.0;
    npe_   = hit->Npe;
    cx_    = hit->Coo[0];
    cy_    = hit->Coo[1];
    cz_    = hit->Coo[2];
    dist_  = dist;
    hit_   = hit;
    
    cross_ = hit->IsCrossed();
}


CherenkovHit RichHitAms::trans() const {
    if (!status_) return CherenkovHit();
    return CherenkovHit(chann_, pmtid_, dbta_, rbtaA_, rbtaB_, npe_, cx_, cy_);
}


void RichHitAms::clear() {
    status_   = false;
    chann_    = -1; 
    pmtid_    = -1;
    pixel_    = -1;
    locid_[0] = -1;
    locid_[1] = -1;
    type_     = -1;
    dbta_     = -1;
    rbtaA_    = -1;
    rbtaB_    = -1;
    npe_      = -1;
    cx_       = 0;
    cy_       = 0;
    cz_       = 0;
    dist_     = 0;
    hit_      = nullptr;
    cross_    = false;
}


RichAms::RichAms(AMSEventR* event, TrTrackR* trtk) : RichAms() {
    if (event == nullptr || event->NRichHit() == 0) return;
    timer_.start();

    TrTrackR* cand_trtk = nullptr;
    if (trtk == nullptr && event->NParticle() != 0 && event->pParticle(0) != nullptr) {
        cand_trtk = event->pParticle(0)->pTrTrack();
    }
    else { cand_trtk = trtk; }
    if (cand_trtk == nullptr) return;
    
    AMSPoint pmtp;
    AMSDir   pmtd;
    cand_trtk->Interpolate(PMT_CZ, pmtp, pmtd);
    pmtp_ = pmtp;

    AMSPoint pnt;
    AMSDir   dir;
    cand_trtk->Interpolate(RichOffline::RICHDB::RICradpos(), pnt, dir);
    RichOffline::TrTrack track(pnt, dir);

    event_ = event;
    trtk_  = cand_trtk;
    
    status_ = build(track);

    timer_.stop();
    if (!status_) clear();
}


void RichAms::clear() {
    status_  = false;
    kind_    = KIND_EMPTY;
    tile_    = -1;
    index_   = 0;
    dist_    = 0;
    dirp_    = AMSPoint();
    dird_    = AMSDir();
    refp_    = AMSPoint();
    refd_    = AMSDir();
    
    radp_    = AMSPoint();
    radd_    = AMSDir();
    pmtp_    = AMSPoint();
    npe_col_ = 0.0;
    bta_crr_ = 1.0;
    
    is_good_in_geometry_ = false;
    is_bad_tile_ = false;

    event_   = nullptr;
    trtk_    = nullptr;
    rhhits_.clear();
    chhits_.clear();
    timer_.clear();
}


CherenkovFit RichAms::fit() const {
    if (!status_ || kind_ == KIND_EMPTY || chhits_.size() == 0) return CherenkovFit();
    std::array<double, 2> pmtc({ pmtp_[0], pmtp_[1] });
    if (kind_ == KIND_AGL) {
        return CherenkovFit(chhits_, pmtc, index_, AGL_BETA_WIDTH, AGL_BETA_PDF, bta_crr_);
    }
    else if (kind_ == KIND_NAF) {
        return CherenkovFit(chhits_, pmtc, index_, NAF_BETA_WIDTH, NAF_BETA_PDF, bta_crr_);
    }
    return CherenkovFit();
}
        

std::array<double, 4> RichAms::cal_trace(double cbta, const CherenkovCloud* cloud) const {
    std::array<double, 2> RAD_HEIGHT { 2.5, 0.5 }; // AGL, NAF
    std::array<double, 4> rlt_trace({ 0.0, 0.0, 1.0, 0.0 });
    bool has_cloud = (cloud != nullptr && cloud->status());
    if (!status_ || kind_ == KIND_EMPTY || (has_cloud ? cloud->cbta() : cbta) <= (Numc::ONE<> / index_)) return rlt_trace;

    double cosv = Numc::ONE<> / (index_ * (has_cloud ? cloud->cbta() : cbta));
    double sinv = std::sqrt((Numc::ONE<> - cosv) * (Numc::ONE<> + cosv));
    if (!Numc::Valid(cosv) || !Numc::Valid(sinv)) return rlt_trace;
    int kind = (kind_ == KIND_AGL) ? 0 : 1;

    double lmtrfr = -std::cos(std::asin(Numc::ONE<> / index_));

    SVecD<3> radp(radp_[0], radp_[1], radp_[2]);
    SVecD<3> radd(radd_[0], radd_[1], radd_[2]);
    OrthCoord orth(radd);
    
    double topx = radp[0] + (Numc::ONE_TO_TWO * RAD_HEIGHT[kind]) * (radd[0] / radd[2]);
    double topy = radp[1] + (Numc::ONE_TO_TWO * RAD_HEIGHT[kind]) * (radd[1] / radd[2]);
    double topz = radp[2] + (Numc::ONE_TO_TWO * RAD_HEIGHT[kind]);

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
        double locwgt = (static_cast<double>(is) + Numc::ONE_TO_TWO) / static_cast<double>(TRACE_NSET);
        
        double phi = Numc::TWO_PI * static_cast<double>(it) / static_cast<double>(TRACE_NPHI);
        SVecD<3>&& phivec = orth.tau() * std::cos(phi) + orth.rho() * std::sin(phi);
        SVecD<3>&& srcvec = radd * cosv + phivec * sinv;
        if (srcvec[2] >= lmtrfr) continue;

        double planex = topx + (-RAD_HEIGHT[kind]) * (locwgt * (radd[0] / radd[2]) + (Numc::ONE<> - locwgt) * (srcvec[0] / srcvec[2]));
        double planey = topy + (-RAD_HEIGHT[kind]) * (locwgt * (radd[1] / radd[2]) + (Numc::ONE<> - locwgt) * (srcvec[1] / srcvec[2]));
        double planez = topz + (-RAD_HEIGHT[kind]);
        if (std::hypot(planex, planey) > MIRROR_TOP_RADIUS) continue;
        
        short tile = RichRingR::getTileIndex(planex, planey);
        if (tile < 0 || tile != tile_) continue;
        border[it][is] = 1; 

        double cosvac = Numc::NEG<> * std::sqrt(Numc::ONE<> - (index_ * index_) * (Numc::ONE<> - srcvec[2]) * (Numc::ONE<> + srcvec[2]));
        double sinvac = std::sqrt((Numc::ONE<> - cosvac) * (Numc::ONE<> + cosvac));
        if (!Numc::Valid(cosvac) || !Numc::Valid(sinvac)) continue;
        SVecD<2>&& vacr = LA::Unit(SVecD<2>(srcvec[0], srcvec[1])) * sinvac;
        if (!Numc::Valid(vacr[0]) || !Numc::Valid(vacr[1])) vacr = std::move(SVecD<2>());

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
        double poly1 = Numc::TWO<> * (MIRROR_TOP_RADIUS * mirtr - (vacp[0] * vactx + vacp[1] * vacty));
        double poly2 = (MIRROR_TOP_RADIUS * MIRROR_TOP_RADIUS - (vacp[0] * vacp[0] + vacp[1] * vacp[1]));

        double det  = std::sqrt(poly1 * poly1 - Numc::FOUR<> * poly0 * poly2);
        double sol1 = Numc::ONE_TO_TWO * (-poly1 + det) / poly0;
        double sol2 = Numc::ONE_TO_TWO * (-poly1 - det) / poly0;
        bool is_ok1 = (Numc::Valid(sol1) && sol1 < 0.0 && sol1 > dheight);
        bool is_ok2 = (Numc::Valid(sol2) && sol2 < 0.0 && sol2 > dheight);

        if (!is_ok1 && !is_ok2) continue;

        double sol = 0;
        if (is_ok1) sol = sol1;
        if (is_ok2) sol = sol2;
        vacp[0] = vacp[0] + sol * vactx;
        vacp[1] = vacp[1] + sol * vacty;
        vacp[2] = vacp[2] + sol;
        
        double cosmir = Numc::NEG<> * std::sqrt((mirtr * mirtr) / (Numc::ONE<> + mirtr * mirtr));
        double sinmir = std::sqrt(Numc::ONE<> - cosmir * cosmir);
        SVecD<2>&& mirr = LA::Unit(SVecD<2>(-vacp[0], -vacp[1])) * sinmir;
        SVecD<3> mird(mirr[0], mirr[1], cosmir);

        double cosrfl = std::fabs(LA::Dot(vacd, mird));
        if (!Numc::Valid(cosrfl)) continue;
        
        vacd = std::move(LA::Unit(vacd + Numc::TWO<> * cosrfl * mird));
        if (vacd[2] >= Numc::ZERO<>) continue;

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
    
    if (!Numc::Valid(trace_bd) || !Numc::Valid(trace)) return rlt_trace;

    if (!has_cloud) {
        rlt_trace = std::move(std::array<double, 4>({ trace_bd, trace, 1.0, 0.0 }));
        return rlt_trace;
    }

    const double thres = 3.4; // pmt width
    const std::vector<CherenkovHit>& hits = cloud->hits();
    std::vector<short>  hphi(hits.size(), -1);
    std::vector<double> hdxy(hits.size(), -1.0);
    for (int is = 0; is < TRACE_NSET; ++is) {
    for (int it = 0; it < TRACE_NPHI; ++it) {
        if (simphs[it][is] == 0) continue;
        for (int ih = 0; ih < hphi.size(); ++ih) {
            const CherenkovHit& hit = hits.at(ih);
            if (simphs[it][is] == 1 && (hit.mode() != 0)) continue;
            if (simphs[it][is] == 2 && (hit.mode() != 1 && hit.mode() != 2)) continue;
            
            double dist = std::hypot(hit.cx() - simphsx[it][is], hit.cy() - simphsy[it][is]);
            if (dist > thres || (hdxy.at(ih) >= 0.0 && hdxy.at(ih) < dist)) continue;
            
            hphi.at(ih) = it;
            hdxy.at(ih) = dist;
        }
    }}

    for (int ih = hphi.size() - 1; ih >= 0; --ih) {
        if (hphi.at(ih) >= 0) continue;
        hphi.erase(hphi.begin() + ih);
    }
    double accuracy = static_cast<double>(hphi.size()) / static_cast<double>(hits.size());

    if (hphi.size() < 2) {
        rlt_trace = std::move(std::array<double, 4>({ trace_bd, trace, accuracy, 0.0 }));
        return rlt_trace;
    }

    std::vector<double> hnpas(hphi.size(), 0.0);
    for (int iphi = 0; iphi < TRACE_NPHI; ++iphi) {
        std::vector<std::pair<short, int>> hpair(hphi.size(), std::pair<short, int>({ -1, -1 }));
        for (int ih = 0; ih < hphi.size(); ++ih) {
            short dphi = std::abs(hphi.at(ih) - iphi);
            if (dphi > TRACE_NPHI / 2) dphi = (TRACE_NPHI - dphi);
            hpair.at(ih).first  = dphi;
            hpair.at(ih).second = ih;
        }
        std::sort(hpair.begin(), hpair.end());
        for (int jt = 0; jt < hpair.size() - 1; ++jt) {
            if (hpair.at(jt+1).first == hpair.at(jt).first) continue;
            hpair.resize(jt + 1);
            break;
        }

        double wgt = Numc::ONE<> / static_cast<double>(hpair.size());
        for (auto&& ipair : hpair) {
            for (int is = 0; is < TRACE_NSET; ++is) {
                hnpas.at(ipair.second) += wgt * (simphs[iphi][is] != 0);
            }
        }
    }
    std::vector<double> hprb(hphi.size(), 0.0);
    for (int ih = 0; ih < hphi.size(); ++ih) hprb.at(ih) = (hnpas.at(ih) / static_cast<double>(count));
    double avgprb = Numc::ONE<> / static_cast<double>(hphi.size());
    double nrmsgm = avgprb * (Numc::ONE<> - avgprb);
    
    double uniformity = Numc::ZERO<>;
    for (int it = 0; it < hphi.size(); ++it) {
        double res = (hprb.at(it) - avgprb) / nrmsgm;
        uniformity += res * res;
    }
    if (!Numc::Valid(uniformity)) {
        rlt_trace = std::move(std::array<double, 4>({ trace_bd, trace, accuracy, 0.0 }));
        return rlt_trace;
    }

    //TVectorD res(hphi.size());
    //TMatrixDSym cov(hphi.size());
    //for (int it = 0; it < hphi.size(); ++it) {
    //    res(it) = hprb.at(it) - avgprb;
    //    for (int jt = it; jt < hphi.size(); ++jt) {
    //        if (it == jt) cov(it, jt) = avgprb * (Numc::ONE<> - avgprb);
    //        else          cov(it, jt) = Numc::NEG<> * avgprb * avgprb;
    //    }
    //}
    //
    //double det = 0;
    //cov.Invert(&det);
    //double uniformity = cov.Similarity(res);
    //if (!Numc::Valid(det) || !Numc::Valid(uniformity)) {
    //    rlt_trace = std::move(std::array<double, 4>({ trace_bd, trace, accuracy, 0.0 }));
    //    return rlt_trace;
    //}

    rlt_trace = std::move(std::array<double, 4>({ trace_bd, trace, accuracy, uniformity }));
    return rlt_trace;
}
        

bool RichAms::build(RichOffline::TrTrack track) {
    RichOffline::RichRadiatorTileManager crossed_tile(&track);
    if(crossed_tile.getkind() == KIND_EMPTY) return false;
    
    kind_   = static_cast<short>(crossed_tile.getkind());
    tile_   = static_cast<short>(crossed_tile.getcurrenttile());
    index_  = static_cast<double>(crossed_tile.getindex());
    dist_   = static_cast<double>(crossed_tile.getdistance());
    dirp_   = crossed_tile.getemissionpoint();
    dird_   = crossed_tile.getemissiondir();
    refp_   = crossed_tile.getemissionpoint(1);
    refd_   = crossed_tile.getemissiondir(1);

    // Charge to AMS coordinate
    AMSPoint amsp = RichOffline::RichAlignment::RichToAMS(dirp_);
    AMSDir   amsd = RichOffline::RichAlignment::RichToAMS(dird_);
    if (std::cos(amsd.gettheta()) > 0) amsd = AMSDir(-amsd[0], -amsd[1], -amsd[2]);
    radp_ = amsp;
    radd_ = amsd;

    if (RichBetaUniformityCorrection::getHead() != nullptr) {
        RichBetaUniformityCorrection* ptr = RichBetaUniformityCorrection::getHead();
        double bta_crr = ptr->getCorrection(amsp[0], amsp[1], amsd[0], amsd[1]);
        if (bta_crr > 0) bta_crr_ = bta_crr;
    }

    double radius = std::hypot(radp_[0], radp_[1]);
    is_good_in_geometry_ = (radius < EXTERNAL_RAD_RADIUS);
    if (is_good_in_geometry_ && kind_ == KIND_AGL) {
        is_good_in_geometry_ = (std::max(std::fabs(radp_[0]), std::fabs(radp_[1])) > RAD_BOUNDARY[0]);
    }
    if (is_good_in_geometry_ && kind_ == KIND_NAF) {
        is_good_in_geometry_ = (std::max(std::fabs(radp_[0]), std::fabs(radp_[1])) < RAD_BOUNDARY[1]);
    }
    tile_ = RichRingR::getTileIndex(radp_[0], radp_[1]);
    if (tile_ < 0) return false;
        
    bool is_bad_tile = false;
    std::array<short, 7>  BAD_TILE_INDEX { 3, 7, 12, 20, 87, 100, 108 };
    for (int idx = 0; idx < BAD_TILE_INDEX.size(); ++idx) {
        if (tile_ == BAD_TILE_INDEX[idx]) is_bad_tile = true;
    }
    is_bad_tile_ = is_bad_tile;

    int bit = event_->nRichRing();

    const double betamin = 0.01;
    const double betamax = 5.00;
    for(RichOffline::RichRawEvent* rawhit = new RichOffline::RichRawEvent(event_); rawhit != nullptr; rawhit = rawhit->next()) {
        if(!rawhit->getbit(RichOffline::ok_status_bit)) continue;
        RichHitR* phit = rawhit->getpointer();

        float recs[3] = { 0., 0., 0. };
        rawhit->reconstruct(dirp_, refp_, dird_, refd_, betamin, betamax, recs, index_, kind_);

        if (recs[0] <= 0 && recs[1] <= 0 && recs[2] <= 0) continue;
        double dbta  = (recs[0] > 0) ? recs[0] : -1.0;
        double rbtaA = (recs[1] > 0) ? recs[1] : -1.0;
        double rbtaB = (recs[2] > 0) ? recs[2] : -1.0;
        
        if (phit == nullptr || phit->Npe <= 0) continue;
        double dist = std::hypot(phit->Coo[0] - pmtp_[0], phit->Coo[1] - pmtp_[1]);

        double pmt_radius = std::hypot(phit->Coo[0], phit->Coo[1]);
        if (pmt_radius > (MIRROR_BTM_RADIUS - WIDTH_CELL)) continue;

        RichHitAms rhhit(phit, dbta, rbtaA, rbtaB, dist);
        if (!rhhit.status()) continue;
        
        CherenkovHit&& chhit = rhhit.trans();
        if (!chhit.status()) continue;
        
        rhhits_.push_back(rhhit);
        chhits_.push_back(chhit);
    }
    if (rhhits_.size() > 1) std::sort(rhhits_.begin(), rhhits_.end(), RichHitAms_sort());
    if (chhits_.size() > 1) std::sort(chhits_.begin(), chhits_.end(), CherenkovHit_sort());

    for (auto&& hit : rhhits_) {
        if (hit.cross()) continue;
        npe_col_ += hit.npe();
    }

    return true;
}


} // namespace TrackSys


#endif // __TRACKLibs_RichAms_C__
#endif // _PGTRACK_ __ROOTSHAREDLIBRARY__ 
