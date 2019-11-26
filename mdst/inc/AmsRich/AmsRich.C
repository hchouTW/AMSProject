#ifndef __AmsRich_C__
#define __AmsRich_C__

#include "AmsRich.h"


AmsRichHit::AmsRichHit(RichHitR* hit, double dbeta, double rbetaA, double rbetaB) : AmsRichHit() {
    if (hit == nullptr || (dbeta <= 0 && rbetaA <= 0 && rbetaB <= 0)) return;
    status_ = true;
    chann_  = hit->Channel;
    pmtid_  = chann_ / NUM_CHANN_IN_2D_PMT;
    pixel_  = chann_ - pmtid_ * NUM_CHANN_IN_2D_PMT;

    type_   = (dbeta > 0) + (rbetaA > 0) * 2 + (rbetaB > 0) * 4;
    dbeta_  = (type_&1 == 1) ? dbeta  : -1.0;
    rbetaA_ = (type_&2 == 2) ? rbetaA : -1.0;
    rbetaB_ = (type_&4 == 4) ? rbetaB : -1.0;
    npe_    = hit->Npe;
    cx_     = hit->Coo[0];
    cy_     = hit->Coo[1];
    cz_     = hit->Coo[2];
    hit_    = hit;
}


void AmsRichHit::clear() {
    status_ = false;
    chann_  = -1; 
    pmtid_  = -1;
    pixel_  = -1;
    type_   = -1;
    dbeta_  = -1;
    rbetaA_ = -1;
    rbetaB_ = -1;
    npe_    = -1;
    cx_     = 0;
    cy_     = 0;
    cz_     = 0;
    hit_    = nullptr;
}


AmsRich::AmsRich(AMSEventR* event) : AmsRich() {
    if (event == nullptr || event->NRichHit() == 0) return;
    if (event->NParticle() == 0 || event->pParticle(0) == nullptr) return;
    TrTrackR* trtk = event->pParticle(0)->pTrTrack();
    if (trtk == nullptr) return;

    AMSPoint pmtp;
    AMSDir   pmtd;
    trtk->Interpolate(PMT_CZ, pmtp, pmtd);
    pmtp_ = pmtp;

    AMSPoint pnt;
    AMSDir   dir;
    trtk->Interpolate(RichOffline::RICHDB::RICradpos(), pnt, dir);
    RichOffline::TrTrack track(pnt, dir);

    event_ = event;
    trtk_  = trtk;
    
    status_ = build(track);

    if (!status_) clear();
}
        

void AmsRich::clear() {
    status_ = false;
    kind_   = KIND_EMPTY;
    tile_   = -1;
    index_  = 0;
    dist_   = 0;
    locx_   = 0;
    locy_   = 0;
    loctha_ = 0;
    locphi_ = 0;

    dirp_   = AMSPoint();
    dird_   = AMSDir();
    refp_   = AMSPoint();
    refd_   = AMSDir();
    
    radp_     = AMSPoint();
    radd_     = AMSDir();
    pmtp_     = AMSPoint();
    beta_crr_ = 1.0;
    
    is_good_geom_ = false;
    is_bad_tile_  = false;

    event_ = nullptr;
    trtk_  = nullptr;
    hits_.clear();
}


bool AmsRich::build(RichOffline::TrTrack track) {
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

    double tile_x = crossed_tile.get_tile_x(tile_);
    double tile_y = crossed_tile.get_tile_y(tile_);
    double tile_locx = dirp_[0] - tile_x;
    double tile_locy = dirp_[1] - tile_y;
    locx_ = tile_locx;
    locy_ = tile_locy;
    loctha_ = track._d.gettheta();
    locphi_ = track._d.getphi();

    // Charge to AMS coordinate
    AMSPoint amsp = RichOffline::RichAlignment::RichToAMS(dirp_);
    AMSDir   amsd = RichOffline::RichAlignment::RichToAMS(dird_);
    if (std::cos(amsd.gettheta()) > 0) amsd = AMSDir(-amsd[0], -amsd[1], -amsd[2]);
    radp_ = amsp;
    radd_ = amsd;

    if (RichBetaUniformityCorrection::getHead() != nullptr) {
        RichBetaUniformityCorrection* ptr = RichBetaUniformityCorrection::getHead();
        double beta_crr = ptr->getCorrection(amsp[0], amsp[1], amsd[0], amsd[1]);
        if (beta_crr > 0) beta_crr_ = beta_crr;
    }

    double radius = std::hypot(radp_[0], radp_[1]);
    is_good_geom_ = (radius < EXTERNAL_RAD_RADIUS);
    if (is_good_geom_ && kind_ == KIND_AGL) {
        is_good_geom_ = (std::max(std::fabs(radp_[0]), std::fabs(radp_[1])) > RAD_BOUNDARY[0]);
    }
    if (is_good_geom_ && kind_ == KIND_NAF) {
        is_good_geom_ = (std::max(std::fabs(radp_[0]), std::fabs(radp_[1])) < RAD_BOUNDARY[1]);
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
        double dbeta  = (recs[0] > 0) ? recs[0] : -1.0;
        double rbetaA = (recs[1] > 0) ? recs[1] : -1.0;
        double rbetaB = (recs[2] > 0) ? recs[2] : -1.0;
        
        if (phit == nullptr || phit->Npe <= 0) continue;

        AmsRichHit rhhit(phit, dbeta, rbetaA, rbetaB);
        if (!rhhit.status()) continue;
        hits_.push_back(rhhit);
    }
    if (hits_.size() > 1) std::sort(hits_.begin(), hits_.end(), AmsRichHit_sort());

    return true;
}
        

#endif // __AmsRich_C__
