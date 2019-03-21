#if defined(_PGTRACK_) || defined(__ROOTSHAREDLIBRARY__)
#ifndef __TRACKLibs_RichAms_C__
#define __TRACKLibs_RichAms_C__


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


std::array<double, 3> RichAms::cal_trace(double cbta) const {
    std::array<double, 3> rlt_trace({0., 0., 0.});
    if (!status_ || kind_ == KIND_EMPTY || cbta <= (Numc::ONE<> / index_)) return rlt_trace;

    double cosv = Numc::ONE<> / (index_ * cbta);
    double sinv = std::sqrt((Numc::ONE<> - cosv) * (Numc::ONE<> + cosv));
    if (!Numc::Valid(cosv) || !Numc::Valid(sinv)) return rlt_trace;
    int kind = (kind_ == KIND_AGL) ? 0 : 1;

    double lmtrfr = -std::cos(std::asin(Numc::ONE<> / index_));

    SVecD<3> radp(radp_[0], radp_[1], radp_[2]);
    SVecD<3> radd(radd_[0], radd_[1], radd_[2]);
    OrthCoord orth(radd);

    int count = 0;
    int count_dir = 0;
    int count_rfl = 0;
    const int nphi = 780;
    std::vector<std::array<double, 2>> simhits(nphi, std::array<double, 2>({0., 0.}));
    for (int it = 0; it < nphi; ++it) {
        double phi = Numc::TWO_PI * static_cast<double>(it) / static_cast<double>(nphi);
        SVecD<3>&& phivec = orth.tau() * std::cos(phi) + orth.rho() * std::sin(phi);
        SVecD<3>&& srcvec = radd * cosv + phivec * sinv;
        if (srcvec[2] >= lmtrfr) continue;

        double planex = radp[0] + (-Numc::ONE_TO_TWO * RAD_HEIGHT[kind]) * srcvec[0] / srcvec[2];
        double planey = radp[1] + (-Numc::ONE_TO_TWO * RAD_HEIGHT[kind]) * srcvec[1] / srcvec[2];
        double planez = radp[2] + (-Numc::ONE_TO_TWO * RAD_HEIGHT[kind]);
        if (std::hypot(planex, planey) > MIRROR_TOP_RADIUS) continue;
        
        short tile = RichRingR::getTileIndex(planex, planey);
        if (tile < 0 || tile != tile_) continue;

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
        if (RichOffline::RichPMTsManager::FindChannel(pmtx, pmty) < 0) continue;
        if (std::fabs(pmtx) < PMT_HOLE[0] && std::fabs(pmty) < PMT_HOLE[1]) continue;
        if (std::hypot(pmtx, pmty) < MIRROR_BTM_RADIUS) {
            simhits.at(count)[0] = pmtx;
            simhits.at(count)[1] = pmty;
            count_dir++;
            count++; 
            continue; 
        }

        double vactx = (vacd[0] / vacd[2]);
        double vacty = (vacd[1] / vacd[2]);
        double mirtr = ((MIRROR_BTM_RADIUS - MIRROR_TOP_RADIUS) / dheight);
        double poly0 = (vacp[0] * vacp[0] + vacp[1] * vacp[1] - MIRROR_TOP_RADIUS * MIRROR_TOP_RADIUS);
        double poly1 = Numc::TWO<> * (vacp[0] * vactx + vacp[1] * vacty - MIRROR_TOP_RADIUS * mirtr);
        double poly2 = (vactx * vactx + vacty * vacty - mirtr * mirtr);

        double det  = std::sqrt(poly1 * poly1 - Numc::FOUR<> * poly2 * poly0);
        double sol1 = Numc::ONE_TO_TWO * (-poly1 + det) / poly2;
        double sol2 = Numc::ONE_TO_TWO * (-poly1 - det) / poly2;
        bool is_ok1 = (Numc::Valid(sol1) && sol1 < 0.0 && sol1 > dheight);
        bool is_ok2 = (Numc::Valid(sol2) && sol2 < 0.0 && sol2 > dheight);
        if (!is_ok1 && !is_ok2) continue;

        double sol = 0;
        if (is_ok1) sol = sol1;
        if (is_ok2) sol = sol2;
        vacp[0] = vacp[0] + sol * vactx;
        vacp[1] = vacp[1] + sol * vacty;
        vacp[2] = vacp[2] + sol;
        
        double mirr = std::hypot(vacp[0], vacp[1]);
        double cosr = Numc::NEG<> * std::sin(std::atan(std::fabs(mirtr)));
        double sinr = std::sqrt((Numc::ONE<> - cosr) * (Numc::ONE<> + cosr));
        if (!Numc::Valid(cosr) || !Numc::Valid(sinr)) continue;
        SVecD<3> mirn(-vacp[0] * sinr / mirr, -vacp[1] * sinr / mirr, cosr);

        double cosrfl = std::fabs(LA::Dot(vacd, mirn));
        if (!Numc::Valid(cosrfl)) continue;
        
        vacd = std::move(LA::Unit(vacd + Numc::TWO<> * cosrfl * mirn));
        if (vacd[2] >= Numc::ZERO<>) continue;

        dheight = (PMT_CZ - vacp[2]);
        pmtx = vacp[0] + dheight * (vacd[0] / vacd[2]);
        pmty = vacp[1] + dheight * (vacd[1] / vacd[2]);
        if (RichOffline::RichPMTsManager::FindChannel(pmtx, pmty) < 0) continue;
        if (std::fabs(pmtx) < PMT_HOLE[0] && std::fabs(pmty) < PMT_HOLE[1]) continue;
        if (std::hypot(pmtx, pmty) > MIRROR_BTM_RADIUS) continue;
 
        simhits.at(count)[0] = pmtx;
        simhits.at(count)[1] = pmty;
        count_rfl++;
        count++;
    }
    simhits.resize(count);
    double trace     = static_cast<double>(count    ) / static_cast<double>(nphi);
    double trace_dir = static_cast<double>(count_dir) / static_cast<double>(nphi);
    double trace_rfl = static_cast<double>(count_rfl) / static_cast<double>(nphi);
    if (!Numc::Valid(trace) || !Numc::Valid(trace_dir) || !Numc::Valid(trace_rfl) || simhits.size() == 0) return rlt_trace;

    rlt_trace = std::move(std::array<double, 3>({ trace, trace_dir, trace_rfl }));
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
