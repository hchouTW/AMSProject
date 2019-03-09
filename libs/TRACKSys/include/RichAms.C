#if defined(_PGTRACK_) || defined(__ROOTSHAREDLIBRARY__)
#ifndef __TRACKLibs_RichAms_C__
#define __TRACKLibs_RichAms_C__


namespace TrackSys {


RichHitAms::RichHitAms(RichHitR* hit, double dbta, double rbta, double dist) : RichHitAms() {
    if (hit == nullptr || (dbta <= 0 && rbta <= 0)) return;
    status_ = true;
    chann_  = hit->Channel;
    pmtid_  = chann_ / 16;
    type_   = (dbta > 0) + (rbta > 0) * 2;
    dbta_   = (type_%2==1) ? dbta : -1.0;
    rbta_   = (type_/2==1) ? rbta : -1.0;
    npe_    = hit->Npe;
    cx_     = hit->Coo[0];
    cy_     = hit->Coo[1];
    cz_     = hit->Coo[2];
    dist_   = dist;
    hit_    = hit;
}


CherenkovHit RichHitAms::trans() const {
    if (!status_) return CherenkovHit();
    return CherenkovHit(chann_, pmtid_, dbta_, rbta_, npe_, cx_, cy_);
}


void RichHitAms::clear() {
    status_ = false;
    chann_  = -1; 
    pmtid_  = -1;
    type_   = -1;
    dbta_   = -1;
    rbta_   = -1;
    npe_    = -1;
    cx_     = 0;
    cy_     = 0;
    cz_     = 0;
    dist_   = 0;
    hit_    = nullptr;
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
    tile_    = 0;
    index_   = 0;
    dist_    = 0;
    dirp_    = AMSPoint();
    dird_    = AMSDir();
    refp_    = AMSPoint();
    refd_    = AMSDir();
    radp_    = AMSPoint();
    radd_    = AMSDir();
    pmtp_    = AMSPoint();
    bta_crr_ = 1.0;
    event_   = nullptr;
    trtk_    = nullptr;
    rhhits_.clear();
    chhits_.clear();
    timer_.clear();
}


bool RichAms::build(RichOffline::TrTrack track) {
    RichOffline::RichRadiatorTileManager crossed_tile(&track);
    if(crossed_tile.getkind() == KIND_EMPTY) return false;
    
    kind_  = static_cast<short>(crossed_tile.getkind());
    tile_  = static_cast<short>(crossed_tile.getcurrenttile());
    index_ = static_cast<double>(crossed_tile.getindex());
    dist_  = static_cast<double>(crossed_tile.getdistance());
    dirp_  = crossed_tile.getemissionpoint();
    dird_  = crossed_tile.getemissiondir();
    refp_  = crossed_tile.getemissionpoint(1);
    refd_  = crossed_tile.getemissiondir(1);

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

    int bit = event_->nRichRing();

    const double betamin = 0.01;
    const double betamax = 5.00;
    for(RichOffline::RichRawEvent* rawhit = new RichOffline::RichRawEvent(event_); rawhit != nullptr; rawhit = rawhit->next()) {
        float recs[3] = { 0., 0., 0. };
        if(!rawhit->getbit(RichOffline::ok_status_bit)) continue;
        rawhit->reconstruct(dirp_, refp_, dird_, refd_, betamin, betamax, recs, index_, kind_);
        rawhit->unsetbit(bit);
        
        if (recs[0] <= 0 && recs[1] <= 0 && recs[2] <= 0) continue;
        RichHitR* phit = rawhit->getpointer();
        double dbta = (recs[0] > 0) ? recs[0] : -1.0;
        double rbta = -1.0;
        if ((recs[1] > 0 || recs[2] > 0) && !(recs[1] > 0 && recs[2] > 0)) {
            rbta = (recs[1] > recs[2]) ? recs[1] : recs[2];
        }
        if (phit == nullptr || phit->Npe <= 0) continue;
        double dist = std::hypot(phit->Coo[0] - pmtp_[0], phit->Coo[1] - pmtp_[1]);

        RichHitAms rhhit(phit, dbta, rbta, dist);
        if (!rhhit.status()) continue;
        
        CherenkovHit&& chhit = rhhit.trans();
        if (!chhit.status()) continue;
        
        rhhits_.push_back(rhhit);
        chhits_.push_back(chhit);
    }
    if (rhhits_.size() > 1) std::sort(rhhits_.begin(), rhhits_.end(), RichHitAms_sort());
    if (chhits_.size() > 1) std::sort(chhits_.begin(), chhits_.end(), CherenkovHit_sort());

    if (chhits_.size() > 1) {
        double sumnpe = 0.0;
        for (auto&& hit : chhits_) sumnpe += hit.npe();
        for (auto&& hit : chhits_) hit.set_wgt(static_cast<double>(chhits_.size()) * hit.npe() / sumnpe);
    }
   
    int cnt = 0;
    for (auto&& hit : chhits_) {
        CERR("[P] PMT %3d WGT %14.8f DIST %14.8f NPE %14.8f [D] BTA %14.8f [R] BTA %14.8f\n", hit.pmtid(), hit.wgt(), rhhits_.at(cnt).dist(), hit.npe(), hit.dbta(), hit.rbta());
        cnt++;
    }


    return true;
}


} // namespace TrackSys


#endif // __TRACKLibs_RichAms_C__
#endif // _PGTRACK_ __ROOTSHAREDLIBRARY__ 
